#!/usr/bin/env python3
"""Validate the identity of the PI's focal gene before any expression interpretation.

Evidence collected (each row of results/target_gene/identity_validation.tsv is one piece):
  1. Official ID tables: Grapedia equivalences (no-copies and with-copies) and Zenodo
     correelation.list.txt, searched for every alias (and every ID they lead to).
  2. Sequence: blastn of the 228-nt query (both strands, megablast + blastn-short-ish settings)
     against v5.1 transcripts and CDS (extracted from the primary genome + GTF/GFF), Zenodo
     PN40024.cds.fa, Ensembl r62 cDNA, and the whole primary genome; spliced minimap2 to genome.
  3. Locus: v5.1 model(s) and Ensembl r62 (vendor) gene(s) overlapping the genomic hit
     (coordinates comparable because sequences are identical; CP<->chr map verified).
Usage: focal_gene_identity.py --root <project_root>
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd

QUERY = ("ATGAGCAGCAGAGAAAGCAACACTGGTAGTGGGAAGAATGAAGCCAAAGTTGGTGGCGATGGAAGCACTGGAAAGGGAAGTGGAACTTCTT"
         "GCAGTGGCGCTGGGATGATGGCGGCCCCAGGTCAGGCAGGTGAATGCATCTCCAGGGAGGCCTTTGAGAGCAACCCCAAGGCATACTTTGCT"
         "GATCTTCATGCTAAGGAAAAGGGCTCTCCTAAACCCAAGAAATAG")
ALIASES = ["Vitvi01g01965", "VIT_01s0011g04570"]
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovhsp"


def sh(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    E = root / "envs/rnaseq/bin"
    out = root / "results/target_gene"
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    ref = root / "references"
    B = ref / "pn40024_v5.1/build"
    (work / "query.fa").write_text(f">focal_query_228nt\n{QUERY}\n")
    assert len(QUERY) == 228
    rows = []

    # ---- 1. ID tables ----
    seen, frontier = set(), list(ALIASES)
    gp = {"no_copies": ref / "grapedia/all_on_T2T_ref_equivalences_filtered_no_copies_symbols.csv",
          "with_copies": ref / "grapedia/all_on_T2T_ref_equivalences_filtered_symbols.csv"}
    tabs = {k: pd.read_csv(v, sep="\t", dtype=str) for k, v in gp.items()}
    while frontier:
        q = frontier.pop()
        if q in seen:
            continue
        seen.add(q)
        for name, t in tabs.items():
            hits = t[(t.query_id == q) | (t.target_id == q)]
            for r in hits.itertuples():
                rows.append({"evidence_class": "id_table", "source": f"grapedia_{name}", "query": q,
                             "hit": r.target_id if r.query_id == q else r.query_id,
                             "detail": f"{r.query_id}({r.query_origin})->{r.target_id}({r.target_origin}); "
                                       f"min_gene%={r.min_gene_percent}; min_exon%={r.min_exon_percent}; "
                                       f"min_CDS%={r.min_CDS_percent}; same_strand={r.same_strand}; "
                                       f"symbol={r._12}; descriptor={r._13}"})
                nxt = r.target_id if r.query_id == q else r.query_id
                if nxt.startswith("Vitvi05_01") and nxt not in seen:
                    frontier.append(nxt)
    cl = pd.read_csv(ref / "pn40024_v5.1/correelation.list.txt", sep="\t", skiprows=1, dtype=str)
    for q in sorted(seen):
        m = cl[cl.apply(lambda r: any(q in str(v).split(",") for v in r.values), axis=1)]
        if m.empty:
            rows.append({"evidence_class": "id_table", "source": "zenodo_correelation.list", "query": q,
                         "hit": "NOT_FOUND", "detail": "alias absent from every column"})
        for r in m.itertuples(index=False):
            rows.append({"evidence_class": "id_table", "source": "zenodo_correelation.list", "query": q,
                         "hit": r[0], "detail": "; ".join(f"{c}={v}" for c, v in zip(cl.columns, r))})

    # ---- 2. sequence databases ----
    if not (work / "v51_transcripts.fa").exists():
        sh(f"{E}/gffread -w {work}/v51_transcripts.fa -g {B}/genome_nuc_org.fa {B}/pn_v51.gtf")
    gff_renamed = work / "v51_renamed.gff3"
    if not gff_renamed.exists():
        smap = dict(pd.read_csv(B / "seqmap_gff_to_fasta.tsv", sep="\t").values)
        with open(gff_renamed, "w") as o:
            for line in open(ref / "pn40024_v5.1/PN_T2T.v5.1.gff3"):
                f = line.split("\t")
                if line.startswith("#") or f[0] in smap:
                    o.write(line)
        sh(f"{E}/gffread -x {work}/v51_cds.fa -g {B}/genome_nuc_org.fa {gff_renamed}")
    if not (work / "ensembl_cdna.fa").exists():
        sh(f"zcat {ref}/ensembl_r62/Vitis_vinifera.ASM3070453v1.cdna.all.fa.gz > {work}/ensembl_cdna.fa")
    dbs = {"v5.1_transcripts": work / "v51_transcripts.fa", "v5.1_CDS": work / "v51_cds.fa",
           "zenodo_PN40024.cds(original_T2T)": ref / "pn40024_v5.1/PN40024.cds.fa",
           "ensembl_r62_cDNA(vendor)": work / "ensembl_cdna.fa", "primary_genome": B / "genome_nuc_org.fa"}
    for name, fa in dbs.items():
        db = work / ("db_" + name.split("(")[0].replace(".", "_"))
        if not Path(str(db) + ".nsq").exists() and not Path(str(db) + ".00.nsq").exists():
            sh(f"{E}/makeblastdb -dbtype nucl -in {fa} -out {db} -parse_seqids")
        for task in ("megablast", "blastn"):
            res = sh(f"{E}/blastn -task {task} -query {work}/query.fa -db {db} -evalue 1e-5 -max_target_seqs 50 "
                     f"-max_hsps 10 -outfmt '{OUTFMT}'")
            for l in res.strip().splitlines():
                f = l.split("\t")
                rows.append({"evidence_class": "sequence", "source": f"blastn_{task}:{name}", "query": "focal_query_228nt",
                             "hit": f[1], "pident": float(f[2]), "aln_len": int(f[3]), "mismatch": int(f[4]),
                             "gapopen": int(f[5]), "qstart": int(f[6]), "qend": int(f[7]), "sstart": int(f[8]),
                             "send": int(f[9]), "evalue": float(f[10]), "bitscore": float(f[11]),
                             "query_cov": int(f[3]) / 228, "subject_len": int(f[13]), "strand": f[14]})
    paf = sh(f"{E}/minimap2 -c -x splice --secondary=yes -N 10 {B}/genome_nuc_org.fa {work}/query.fa 2>/dev/null")
    for l in paf.strip().splitlines():
        f = l.split("\t")
        rows.append({"evidence_class": "sequence", "source": "minimap2_splice:primary_genome", "query": f[0],
                     "hit": f"{f[5]}:{int(f[7]) + 1}-{f[8]}", "strand": f[4], "qstart": int(f[2]) + 1, "qend": int(f[3]),
                     "pident": 100 * int(f[9]) / int(f[10]), "aln_len": int(f[10]), "query_cov": (int(f[3]) - int(f[2])) / 228,
                     "detail": " ".join(x for x in f[12:] if x.split(":")[0] in ("tp", "NM", "cg", "ts"))})
    seq = pd.DataFrame([r for r in rows if r["evidence_class"] == "sequence"])

    # ---- 3. locus: overlapping v5.1 and Ensembl genes at every genome hit (>=90% id, >=50% cov) ----
    ens_map = dict(pd.read_csv(B / "seqmap_ensembl_to_fasta.tsv", sep="\t")[["fasta_seqid", "ensembl_seqid"]].values)
    v51 = pd.read_csv(B / "v51.genes.tsv", sep="\t")
    smap = dict(pd.read_csv(B / "seqmap_gff_to_fasta.tsv", sep="\t").values)
    v51["fasta_seqid"] = v51.seqid.map(smap)
    egtf = pd.read_csv(ref / "ensembl_r62/build/ensembl_r62.gtf", sep="\t", comment="#", header=None,
                       names=["seq", "src", "type", "start", "end", "score", "strand", "phase", "attr"])
    egenes = egtf[egtf.type == "gene"].copy()
    egenes["gene_id"] = egenes.attr.str.extract(r'gene_id "([^"]+)"')
    egenes["name"] = egenes.attr.str.extract(r'gene_name "([^"]+)"')
    egenes["biotype"] = egenes.attr.str.extract(r'gene_biotype "([^"]+)"')
    ge = seq[(seq.source.str.startswith("blastn")) & (seq.source.str.endswith("primary_genome"))
             & (seq.pident >= 90) & (seq.query_cov >= 0.5)]
    for r in ge.itertuples():
        lo, hi = sorted((r.sstart, r.send))
        strand = "+" if r.strand == "plus" else "-"
        for g in v51[(v51.fasta_seqid == r.hit) & (v51.start <= hi) & (v51.end >= lo)].itertuples():
            rows.append({"evidence_class": "locus", "source": "v5.1_gene_overlapping_genome_hit", "query": f"{r.hit}:{lo}-{hi}({strand})",
                         "hit": g.gene_id, "strand": g.strand,
                         "detail": f"gene {g.seqid}:{g.start}-{g.end}({g.strand}); same_strand_as_query={g.strand == strand}; "
                                   f"n_tx={g.n_transcripts}; max_3UTR_bp={g.max_utr3_bp}; max_5UTR_bp={g.max_utr5_bp}; coding={g.coding}"})
        ecp = ens_map.get(r.hit)
        for g in egenes[(egenes.seq == ecp) & (egenes.start <= hi) & (egenes.end >= lo)].itertuples():
            rows.append({"evidence_class": "locus", "source": "ensembl_r62_gene_overlapping_genome_hit(vendor_ID)",
                         "query": f"{r.hit}:{lo}-{hi}({strand})", "hit": g.gene_id, "strand": g.strand,
                         "detail": f"gene {ecp}:{g.start}-{g.end}({g.strand}); name={g.name}; biotype={g.biotype}; "
                                   f"same_strand_as_query={g.strand == strand}"})
    df = pd.DataFrame(rows)
    df.to_csv(out / "identity_validation.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 140)
    print(df[df.evidence_class != "sequence"][["source", "query", "hit", "detail"]].to_string(index=False))
    print(seq.sort_values("bitscore", ascending=False).groupby("source").head(4)[
        ["source", "hit", "pident", "aln_len", "qstart", "qend", "sstart", "send", "strand", "evalue"]].to_string(index=False))


if __name__ == "__main__":
    main()
