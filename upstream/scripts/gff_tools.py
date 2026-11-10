#!/usr/bin/env python3
"""GFF3 validation + GFF3->GTF conversion for gene-level 3' counting.

validate: seqids vs FASTA (after an optional verified seqid map), unique IDs, parent->child
integrity (multi-parent exons allowed), feature-type counts, per-transcript inferred UTRs
(exon span beyond CDS; UTR features are not invented as GFF features, only measured),
genes without any 3' UTR, coordinate bounds.
to-gtf: writes exon records (gene_id = gene ID, transcript_id) for every transcript, splitting
multi-parent exons into one record per transcript, renaming seqids with the map; seqids absent
from the map/FASTA are dropped and reported.

Usage:
  gff_tools.py validate --gff X.gff3 --fasta-fai genome.fa.fai [--seqmap map.tsv] --out-prefix P
  gff_tools.py to-gtf   --gff X.gff3 [--seqmap map.tsv] [--keep-seqids-from genome.fa.fai] --out X.gtf
seqmap: TSV with columns gff_seqid, fasta_seqid (only identity-verified pairs).
"""
from __future__ import annotations

import argparse
import collections
import sys
from dataclasses import dataclass, field
from urllib.parse import unquote

import pandas as pd

TRANSCRIPT_TYPES = {"mRNA", "lncRNA", "transcript", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA",
                    "miRNA", "pseudogenic_transcript"}


def parse_attrs(s: str) -> dict[str, str]:
    out = {}
    for kv in s.strip().strip(";").split(";"):
        if not kv:
            continue
        k, _, v = kv.partition("=")
        out[k.strip()] = unquote(v.strip())
    return out


@dataclass
class Model:
    genes: dict = field(default_factory=dict)          # gene_id -> (seqid, start, end, strand, attrs)
    transcripts: dict = field(default_factory=dict)    # tx_id -> (gene_id, type, seqid, start, end, strand)
    exons: dict = field(default_factory=lambda: collections.defaultdict(list))  # tx -> [(s,e)]
    cds: dict = field(default_factory=lambda: collections.defaultdict(list))
    types: collections.Counter = field(default_factory=collections.Counter)
    problems: list = field(default_factory=list)
    dup_ids: collections.Counter = field(default_factory=collections.Counter)


def read_gff(path: str) -> Model:
    m = Model()
    seen_ids = collections.Counter()
    pending_children = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                if line.startswith("##FASTA"):
                    break
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                m.problems.append(("malformed_line", ln, line[:80]))
                continue
            seqid, src, ftype, start, end, score, strand, phase, attr = f
            start, end = int(start), int(end)
            if start > end:
                m.problems.append(("start_gt_end", ln, attr[:80]))
            a = parse_attrs(attr)
            m.types[ftype] += 1
            fid = a.get("ID")
            if fid and ftype not in ("CDS",):   # CDS IDs are legitimately repeated per segment
                seen_ids[fid] += 1
            parents = a.get("Parent", "").split(",") if a.get("Parent") else []
            if ftype in ("gene", "pseudogene", "ncRNA_gene"):
                m.genes[fid] = (seqid, start, end, strand, a)
            elif ftype in TRANSCRIPT_TYPES:
                if len(parents) != 1:
                    m.problems.append(("transcript_parent_count", ln, fid))
                for p in parents[:1]:
                    m.transcripts[fid] = (p, ftype, seqid, start, end, strand)
            elif ftype == "exon":
                pending_children.append(("exon", parents, start, end, ln))
            elif ftype == "CDS":
                pending_children.append(("CDS", parents, start, end, ln))
    for kind, parents, s, e, ln in pending_children:
        if not parents:
            m.problems.append((f"{kind}_without_parent", ln, ""))
        for p in parents:
            if p not in m.transcripts:
                m.problems.append((f"{kind}_orphan_parent", ln, p))
                continue
            (m.exons if kind == "exon" else m.cds)[p].append((s, e))
    for tx, (g, *_rest) in m.transcripts.items():
        if g not in m.genes:
            m.problems.append(("transcript_orphan_gene", 0, tx))
    m.dup_ids = collections.Counter({k: v for k, v in seen_ids.items() if v > 1})
    return m


def read_seqmap(path: str | None) -> dict[str, str] | None:
    if not path:
        return None
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df.gff_seqid, df.fasta_seqid))


def fai_lengths(path: str) -> dict[str, int]:
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(path)}


def tx_utr(m: Model, tx: str) -> tuple[int, int]:
    """(5' UTR bp, 3' UTR bp) inferred from exons minus CDS span, strand-aware."""
    ex = sorted(m.exons.get(tx, []))
    cd = m.cds.get(tx, [])
    if not ex or not cd:
        return (0, 0)
    cmin = min(s for s, _ in cd)
    cmax = max(e for _, e in cd)
    left = sum(max(0, min(e, cmin - 1) - s + 1) for s, e in ex if s < cmin)
    right = sum(max(0, e - max(s, cmax + 1) + 1) for s, e in ex if e > cmax)
    strand = m.transcripts[tx][5]
    return (left, right) if strand == "+" else (right, left)


def cmd_validate(a) -> None:
    m = read_gff(a.gff)
    smap = read_seqmap(a.seqmap)
    lens = fai_lengths(a.fasta_fai)
    rows = []

    def add(check, status, detail):
        rows.append({"check": check, "status": status, "detail": detail})

    add("feature_types", "info", "; ".join(f"{k}={v}" for k, v in sorted(m.types.items())))
    add("n_genes", "info", len(m.genes))
    add("n_transcripts", "info", len(m.transcripts))
    add("duplicate_ids", "pass" if not m.dup_ids else "fail",
        f"{len(m.dup_ids)} duplicated; e.g. {list(m.dup_ids)[:5]}")
    pc = collections.Counter(p[0] for p in m.problems)
    add("structure_problems", "pass" if not pc else "fail", dict(pc))
    multi_parent = sum(1 for t in m.exons for _ in [0])  # placeholder count of tx with exons
    add("transcripts_with_exons", "info", multi_parent)
    no_exon = [t for t in m.transcripts if t not in m.exons]
    add("transcripts_without_exons", "pass" if not no_exon else "warn", f"{len(no_exon)} e.g. {no_exon[:3]}")
    gseq = collections.Counter(v[0] for v in m.genes.values())
    unmapped_seq = {}
    oob = 0
    for s, n in gseq.items():
        fs = smap.get(s, None) if smap is not None else s
        if fs is None or fs not in lens:
            unmapped_seq[s] = n
    for g, (s, st, en, *_r) in m.genes.items():
        fs = smap.get(s) if smap is not None else s
        if fs in lens and en > lens[fs]:
            oob += 1
    add("gene_seqids_without_fasta_sequence", "pass" if not unmapped_seq else "warn",
        f"{sum(unmapped_seq.values())} genes on {unmapped_seq}")
    add("genes_out_of_bounds", "pass" if oob == 0 else "fail", oob)

    # per-gene UTR / 3' end summary (longest-3'UTR transcript)
    grows = []
    by_gene = collections.defaultdict(list)
    for tx, (g, t, *_r) in m.transcripts.items():
        by_gene[g].append(tx)
    for g, (s, st, en, strand, attrs) in m.genes.items():
        txs = by_gene.get(g, [])
        utr = [tx_utr(m, t) for t in txs]
        coding = any(m.cds.get(t) for t in txs)
        grows.append({"gene_id": g, "seqid": s, "start": st, "end": en, "strand": strand,
                      "n_transcripts": len(txs),
                      "transcript_types": ",".join(sorted({m.transcripts[t][1] for t in txs})),
                      "coding": coding,
                      "max_utr5_bp": max((u[0] for u in utr), default=0),
                      "max_utr3_bp": max((u[1] for u in utr), default=0),
                      "pseudogene": attrs.get("pseudogene", ""),
                      "transposable": attrs.get("transposable", ""),
                      "has_fasta_sequence": s not in unmapped_seq})
    gdf = pd.DataFrame(grows)
    gdf.to_csv(a.out_prefix + ".genes.tsv", sep="\t", index=False)
    cod = gdf[gdf.coding]
    add("coding_genes_with_3utr", "info",
        f"{(cod.max_utr3_bp > 0).sum()}/{len(cod)} coding genes have exon sequence 3' of CDS "
        f"(median {cod.max_utr3_bp[cod.max_utr3_bp > 0].median():.0f} bp); no explicit UTR features")
    out = pd.DataFrame(rows)
    out.to_csv(a.out_prefix + ".validation.tsv", sep="\t", index=False)
    pd.DataFrame(m.problems, columns=["problem", "line", "detail"]).to_csv(
        a.out_prefix + ".problems.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


def cmd_to_gtf(a) -> None:
    m = read_gff(a.gff)
    smap = read_seqmap(a.seqmap)
    keep = set(fai_lengths(a.keep_seqids_from)) if a.keep_seqids_from else None
    dropped = collections.Counter()
    n = 0
    with open(a.out + ".tmp", "w") as out:
        for tx, (g, ttype, s, st, en, strand) in sorted(m.transcripts.items(),
                                                        key=lambda kv: (kv[1][2], kv[1][3])):
            fs = smap.get(s) if smap is not None else s
            if fs is None or (keep is not None and fs not in keep):
                dropped[s] += 1
                continue
            gattr = f'gene_id "{a.gene_prefix}{g}"; transcript_id "{a.gene_prefix}{tx}"; transcript_type "{ttype}";'
            out.write(f"{fs}\t{a.source}\ttranscript\t{st}\t{en}\t.\t{strand}\t.\t{gattr}\n")
            for es, ee in sorted(m.exons.get(tx, [])):
                out.write(f"{fs}\t{a.source}\texon\t{es}\t{ee}\t.\t{strand}\t.\t{gattr}\n")
                n += 1
    import os
    os.replace(a.out + ".tmp", a.out)
    print(f"wrote {n} exon records; dropped transcripts by seqid: {dict(dropped)}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--gff", required=True)
    v.add_argument("--fasta-fai", required=True)
    v.add_argument("--seqmap")
    v.add_argument("--out-prefix", required=True)
    t = sub.add_parser("to-gtf")
    t.add_argument("--gff", required=True)
    t.add_argument("--seqmap")
    t.add_argument("--keep-seqids-from")
    t.add_argument("--gene-prefix", default="")
    t.add_argument("--source", default="gff_tools")
    t.add_argument("--out", required=True)
    a = ap.parse_args()
    {"validate": cmd_validate, "to-gtf": cmd_to_gtf}[a.cmd](a)


if __name__ == "__main__":
    main()
