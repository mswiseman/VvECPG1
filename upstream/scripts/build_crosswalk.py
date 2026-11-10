#!/usr/bin/env python3
"""Bidirectional gene-ID crosswalk centred on PN40024 T2T v5.1 gene IDs.

One row of the long table = one supported mapping edge between a v5.1 gene and another ID.
Evidence sources (all recorded with provenance; nothing positional/nearest-neighbour):
  grapedia_equivalence      official Grapedia v5.1<->{v4.3,v3,v2,v1,v0,8x,12XNCBI} (no-copies file)
  grapedia_liftoff_copy     additional edges present only in the with-copies file
  zenodo_correelation       Zenodo correelation.list.txt (original T2T Vitvi0000NN <-> v5.1)
  exon_overlap              Ensembl r62 (vendor `Vitis..g.....`) vs v5.1 on the identical assembly:
                            same-strand exonic bp overlap, fractions of each model
  cds_identity              identical CDS sequence (MD5) between vendor and v5.1 / original T2T
  chardonnay_*              appended by crosswalk_chardonnay.py when available
mapping_status (per non-v5.1 ID, from 'strong' edges): unique | one_to_many | many_to_one |
  ambiguous (one_to_many AND many_to_one) | unmapped. Strong vendor edge: same strand and
  >=50% of the exonic bp of either model overlapped, or identical CDS.
Usage: build_crosswalk.py --root <project_root>
"""
from __future__ import annotations

import argparse
import bisect
import collections
import hashlib
import os
import subprocess
from pathlib import Path

import pandas as pd
import pysam

NS = {"5.1": "pn_v5_1", "v4.3": "vitvi_v4", "v3": "vitvi_v3", "v2": "VIT_v2", "v1": "VIT_v1",
      "v0": "v0_12Xv0", "8x": "8x", "12XNCBI": "ncbi_12X_LOC"}


def exon_unions(gtf: str, seqmap: dict | None = None) -> dict:
    ex = collections.defaultdict(list)
    meta = {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) != 9 or f[2] != "exon":
            continue
        gid = f[8].split('gene_id "')[1].split('"')[0]
        chrom = seqmap.get(f[0], f[0]) if seqmap else f[0]
        ex[gid].append((int(f[3]), int(f[4])))
        meta[gid] = (chrom, f[6])
    out = {}
    for g, iv in ex.items():
        iv.sort()
        merged = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[g] = (meta[g][0], meta[g][1], [tuple(x) for x in merged])
    return out


def overlap_bp(a, b) -> int:
    i = j = tot = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s <= e:
            tot += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def cds_md5(fasta: str, tx_to_gene=None) -> pd.DataFrame:
    rows = []
    with pysam.FastxFile(fasta) as fh:
        for r in fh:
            s = r.sequence.upper()
            rows.append({"transcript_id": r.name, "cds_md5": hashlib.md5(s.encode()).hexdigest(), "cds_len": len(s)})
    return pd.DataFrame(rows)


def classify(edges: pd.DataFrame, other_col: str) -> pd.Series:
    """mapping_status per other-ID from strong edges (v5.1 gene <-> other)."""
    st = {}
    n_v = edges.groupby(other_col).pn_v5_1_gene_id.nunique()
    n_o = edges.groupby("pn_v5_1_gene_id")[other_col].nunique()
    for o, k in n_v.items():
        targets = edges[edges[other_col] == o].pn_v5_1_gene_id.unique()
        many_to_one = any(n_o[t] > 1 for t in targets)
        st[o] = ("ambiguous" if k > 1 and many_to_one else "one_to_many" if k > 1
                 else "many_to_one" if many_to_one else "unique")
    return pd.Series(st)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    E = root / "envs/rnaseq/bin"
    ref = root / "references"
    B = ref / "pn40024_v5.1/build"
    out = root / "results/crosswalk"
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)

    v51 = pd.read_csv(B / "v51.genes.tsv", sep="\t")
    smap = dict(pd.read_csv(B / "seqmap_gff_to_fasta.tsv", sep="\t").values)
    v51["chr"] = v51.seqid.map(smap).fillna(v51.seqid)
    vinfo = v51.set_index("gene_id")
    edges = []

    # ---- Grapedia ----
    gd = {}
    for kind, fn in (("no_copies", "all_on_T2T_ref_equivalences_filtered_no_copies_symbols.csv"),
                     ("with_copies", "all_on_T2T_ref_equivalences_filtered_symbols.csv")):
        t = pd.read_csv(ref / "grapedia" / fn, sep="\t", dtype=str)
        t = t[(t.query_origin == "5.1") | (t.target_origin == "5.1")]
        t = t[t.target_id.notna() & (t.target_id != "")]
        fw = t[t.query_origin == "5.1"].rename(columns={"query_id": "v51", "target_id": "other", "target_origin": "ns"})
        rv = t[t.target_origin == "5.1"].rename(columns={"target_id": "v51", "query_id": "other", "query_origin": "ns"})
        fw["direction"], rv["direction"] = "5.1->other", "other->5.1"
        gd[kind] = pd.concat([fw, rv])[["v51", "other", "ns", "direction", "min_gene_percent", "min_exon_percent",
                                        "min_CDS_percent", "same_strand", "overlap_score",
                                        "T2T.v5.1_gene_symbol", "T2T.v5.1_descriptor"]]
    nc = gd["no_copies"]
    keys_nc = set(zip(nc.v51, nc.other, nc.ns))
    for kind, t in gd.items():
        agg = t.groupby(["v51", "other", "ns"]).agg(
            directions=("direction", lambda x: ",".join(sorted(set(x)))),
            min_gene_percent=("min_gene_percent", "first"), min_exon_percent=("min_exon_percent", "first"),
            min_CDS_percent=("min_CDS_percent", "first"), same_strand=("same_strand", "first"),
            overlap_score=("overlap_score", "first"), symbol=("T2T.v5.1_gene_symbol", "first"),
            descriptor=("T2T.v5.1_descriptor", "first")).reset_index()
        if kind == "with_copies":
            agg = agg[[k not in keys_nc for k in zip(agg.v51, agg.other, agg.ns)]]
        for r in agg.itertuples():
            edges.append({"pn_v5_1_gene_id": r.v51, "other_id": r.other, "other_namespace": NS.get(r.ns, r.ns),
                          "source": "grapedia_equivalence" if kind == "no_copies" else "grapedia_liftoff_copy",
                          "evidence_type": "annotation_transfer_overlap(Grapedia)",
                          "identity": None, "coverage": float(r.min_exon_percent) / 100 if pd.notna(r.min_exon_percent) else None,
                          "reciprocal": r.directions == "5.1->other,other->5.1",
                          "notes": f"directions={r.directions}; min_gene%={r.min_gene_percent}; min_exon%={r.min_exon_percent}; "
                                   f"min_CDS%={r.min_CDS_percent}; same_strand={r.same_strand}; overlap_score={r.overlap_score}",
                          "symbol": r.symbol, "descriptor": r.descriptor})

    # ---- Zenodo correelation.list (original T2T Vitvi0000NN <-> v5.1) ----
    cl = pd.read_csv(ref / "pn40024_v5.1/correelation.list.txt", sep="\t", skiprows=1, dtype=str)
    for r in cl.itertuples(index=False):
        t2t, v = r[0], r[1]
        if isinstance(v, str) and v != "NA":
            for vv in v.split(","):
                edges.append({"pn_v5_1_gene_id": vv, "other_id": t2t, "other_namespace": "zenodo_original_T2T",
                              "source": "zenodo_correelation", "evidence_type": "published_correspondence_table",
                              "reciprocal": None, "notes": "Zenodo correelation.list.txt col PN40024_T2T -> PN40024_T2T.v5.1"})

    # ---- Ensembl r62 (vendor IDs) vs v5.1: exon overlap on identical assembly ----
    ens_map = dict(pd.read_csv(B / "seqmap_ensembl_to_fasta.tsv", sep="\t")[["ensembl_seqid", "fasta_seqid"]].values)
    eg = exon_unions(str(ref / "ensembl_r62/build/ensembl_r62.gtf"), ens_map)
    vg = exon_unions(str(B / "pn_v51.gtf"))
    idx = collections.defaultdict(list)
    for g, (c, s, ex) in vg.items():
        idx[(c, s)].append((ex[0][0], ex[-1][1], g))
    for k in idx:
        idx[k].sort()
    starts = {k: [x[0] for x in v] for k, v in idx.items()}
    maxspan = {k: max(x[1] - x[0] for x in v) for k, v in idx.items()}
    ov_rows = []
    for g, (c, s, ex) in eg.items():
        lo, hi = ex[0][0], ex[-1][1]
        L = idx.get((c, s), [])
        i0 = bisect.bisect_left(starts.get((c, s), []), lo - maxspan.get((c, s), 0))
        elen = sum(e - b + 1 for b, e in ex)
        for vs, ve, vgid in L[i0:]:
            if vs > hi:
                break
            if ve < lo:
                continue
            bp = overlap_bp(ex, vg[vgid][2])
            if bp > 0:
                vlen = sum(e - b + 1 for b, e in vg[vgid][2])
                ov_rows.append({"vendor_id": g, "v51": vgid, "overlap_bp": bp, "frac_vendor": bp / elen, "frac_v51": bp / vlen})
    ov = pd.DataFrame(ov_rows)

    # CDS identity (vendor vs v5.1 and vendor vs Zenodo original T2T CDS)
    eg_chr = work / "ensembl_r62_chrnames.gtf"
    if not eg_chr.exists():
        with open(eg_chr, "w") as o:
            for line in open(ref / "ensembl_r62/build/ensembl_r62.gtf"):
                f = line.split("\t")
                if not line.startswith("#") and f[0] in ens_map:
                    f[0] = ens_map[f[0]]
                    o.write("\t".join(f))
    if not (work / "ensembl_cds.fa").exists():
        subprocess.run(f"{E}/gffread -x {work}/ensembl_cds.fa -g {B}/genome_nuc_org.fa {eg_chr}", shell=True, check=True)
    v51_cds = root / "results/target_gene/work/v51_cds.fa"
    if not v51_cds.exists():
        raise SystemExit("run focal_gene_identity.py first (creates v5.1 CDS FASTA)")
    ecds = cds_md5(str(work / "ensembl_cds.fa"))
    etx2g = {}
    for line in open(eg_chr):
        f = line.split("\t")
        if len(f) == 9 and f[2] == "transcript":
            etx2g[f[8].split('transcript_id "')[1].split('"')[0]] = f[8].split('gene_id "')[1].split('"')[0]
    ecds["gene"] = ecds.transcript_id.map(etx2g)
    vcds = cds_md5(str(v51_cds))
    vcds["gene"] = vcds.transcript_id.str.replace(r"_t\d+$", "", regex=True)
    zcds = cds_md5(str(ref / "pn40024_v5.1/PN40024.cds.fa"))
    zcds["gene"] = zcds.transcript_id.str.replace(r"\.t\d+$", "", regex=True)
    ev = ecds.merge(vcds, on="cds_md5", suffixes=("_vendor", "_v51"))
    cds_pairs = set(zip(ev.gene_vendor, ev.gene_v51))
    ez = ecds.merge(zcds, on="cds_md5", suffixes=("_vendor", "_t2t"))

    for r in ov.itertuples():
        same_cds = (r.vendor_id, r.v51) in cds_pairs
        edges.append({"pn_v5_1_gene_id": r.v51, "other_id": r.vendor_id, "other_namespace": "vendor_ensembl_r62",
                      "source": "exon_overlap", "evidence_type": "coordinate_exon_overlap(identical_assembly)",
                      "identity": 1.0 if same_cds else None, "coverage": min(r.frac_vendor, r.frac_v51),
                      "reciprocal": None, "strong": (max(r.frac_vendor, r.frac_v51) >= 0.5) or same_cds,
                      "notes": f"overlap_bp={r.overlap_bp}; frac_of_vendor_exons={r.frac_vendor:.3f}; "
                               f"frac_of_v51_exons={r.frac_v51:.3f}; identical_CDS={same_cds}"})
    ovp = set(zip(ov.vendor_id, ov.v51)) if len(ov) else set()
    for (vend, v) in cds_pairs - ovp:
        edges.append({"pn_v5_1_gene_id": v, "other_id": vend, "other_namespace": "vendor_ensembl_r62",
                      "source": "cds_identity", "evidence_type": "identical_CDS_sequence_md5", "identity": 1.0,
                      "coverage": 1.0, "strong": True, "notes": "identical CDS but no same-strand exon overlap (copy elsewhere)"})
    # vendor <-> Zenodo original T2T (no v5.1 node): kept in a side table
    ez[["gene_vendor", "gene_t2t", "transcript_id_vendor", "transcript_id_t2t", "cds_len_vendor"]].drop_duplicates() \
        .to_csv(out / "vendor_to_zenodo_original_T2T_by_identical_CDS.tsv", sep="\t", index=False)

    long = pd.DataFrame(edges)
    long["strong"] = long["strong"].fillna(True).astype(bool)
    long = long.merge(v51[["gene_id", "chr", "start", "end", "strand"]].rename(columns={"gene_id": "pn_v5_1_gene_id"}),
                      on="pn_v5_1_gene_id", how="left")
    long["reference"] = "PN40024_T2T_v5.1"
    long["pn_v5_1_in_annotation"] = long.pn_v5_1_gene_id.isin(vinfo.index)
    long["pn_v5_1_quantifiable"] = long.pn_v5_1_gene_id.map(vinfo.has_fasta_sequence).fillna(False)
    # namespace-specific columns requested in the spec
    for col, nss in {"vitis_original_t2t_id": ["vendor_ensembl_r62"],
                     "vitvi_v4_or_v3_id": ["vitvi_v4", "vitvi_v3"],
                     "VIT_legacy_id": ["VIT_v1", "VIT_v2"],
                     "chardonnay_gene_id": ["chardonnay_hap1", "chardonnay_hap2"]}.items():
        long[col] = long.other_id.where(long.other_namespace.isin(nss))
    long["pn_v5_1_transcript_id"] = None
    long["reciprocal_best_hit"] = long.pop("reciprocal")

    # mapping status per (namespace, other_id) from strong edges
    long["mapping_status"] = None
    for ns_, sub in long[long.strong].groupby("other_namespace"):
        st = classify(sub, "other_id")
        m = long.other_namespace == ns_
        long.loc[m, "mapping_status"] = long.loc[m, "other_id"].map(st)
    long.loc[long.mapping_status.isna(), "mapping_status"] = "weak_only"

    cols = ["pn_v5_1_gene_id", "pn_v5_1_transcript_id", "other_id", "other_namespace", "vitis_original_t2t_id",
            "vitvi_v4_or_v3_id", "VIT_legacy_id", "chardonnay_gene_id", "reference", "chr", "start", "end", "strand",
            "source", "evidence_type", "identity", "coverage", "reciprocal_best_hit", "strong", "mapping_status",
            "pn_v5_1_in_annotation", "pn_v5_1_quantifiable", "symbol", "descriptor", "notes"]
    long = long[cols].sort_values(["pn_v5_1_gene_id", "other_namespace", "other_id", "source"])
    long.to_csv(out / "gene_id_crosswalk_long.tsv.tmp", sep="\t", index=False)
    os.replace(out / "gene_id_crosswalk_long.tsv.tmp", out / "gene_id_crosswalk_long.tsv")

    # ---- wide lookup (one row per v5.1 gene) ----
    strong = long[long.strong]
    def joined(ns_list):
        s = strong[strong.other_namespace.isin(ns_list)]
        return s.groupby("pn_v5_1_gene_id").other_id.apply(lambda x: ";".join(sorted(set(x))))
    wide = v51[["gene_id", "chr", "start", "end", "strand", "coding", "max_utr3_bp", "has_fasta_sequence"]].rename(
        columns={"gene_id": "pn_v5_1_gene_id", "has_fasta_sequence": "quantifiable_on_primary_genome"}).set_index("pn_v5_1_gene_id")
    for col, nss in {"vendor_ensembl_r62_ids": ["vendor_ensembl_r62"], "zenodo_original_T2T_ids": ["zenodo_original_T2T"],
                     "vitvi_v4_ids": ["vitvi_v4"], "vitvi_v3_ids": ["vitvi_v3"], "VIT_v1_ids": ["VIT_v1"],
                     "VIT_v2_ids": ["VIT_v2"], "ncbi_12X_ids": ["ncbi_12X_LOC"],
                     "chardonnay_hap1_ids": ["chardonnay_hap1"], "chardonnay_hap2_ids": ["chardonnay_hap2"]}.items():
        wide[col] = joined(nss)
    sym = long.dropna(subset=["symbol"]).groupby("pn_v5_1_gene_id").symbol.first()
    desc = long.dropna(subset=["descriptor"]).groupby("pn_v5_1_gene_id").descriptor.first()
    wide["symbol"], wide["descriptor"] = sym, desc
    vstat = strong[strong.other_namespace == "vendor_ensembl_r62"].groupby("pn_v5_1_gene_id").mapping_status \
        .apply(lambda x: ";".join(sorted(set(x))))
    wide["vendor_mapping_status"] = vstat.reindex(wide.index).fillna("no_vendor_gene")
    wide.reset_index().to_csv(out / "gene_id_lookup_wide.tsv", sep="\t", index=False)

    # ---- vendor-centric reports ----
    vend = pd.read_csv(root / "PPPJ8P_results/PPPJ8P-expression-matrix.tsv", sep="\t", usecols=["gene_id", "gene_name", "gene_biotype"])
    vs = long[long.other_namespace == "vendor_ensembl_r62"]
    best = vs[vs.strong].groupby("other_id").mapping_status.first()
    vend["mapping_status"] = vend.gene_id.map(best)
    weak = set(vs[~vs.strong].other_id)
    vend.loc[vend.mapping_status.isna() & vend.gene_id.isin(weak), "mapping_status"] = "weak_overlap_only"
    vend["mapping_status"] = vend.mapping_status.fillna("unmapped")
    vend["v51_targets"] = vend.gene_id.map(vs[vs.strong].groupby("other_id").pn_v5_1_gene_id.apply(lambda x: ";".join(sorted(set(x)))))
    vend.to_csv(out / "vendor_id_mapping_status.tsv", sep="\t", index=False)
    vend[vend.mapping_status.isin(["unmapped", "weak_overlap_only"])].to_csv(out / "unmapped_vendor_ids.tsv", sep="\t", index=False)
    long[long.mapping_status.isin(["one_to_many", "many_to_one", "ambiguous"])].to_csv(
        out / "ambiguous_mappings.tsv", sep="\t", index=False)
    stats = []
    for ns_, sub in long.groupby("other_namespace"):
        ids = sub.drop_duplicates("other_id").mapping_status.value_counts()
        for k, v in ids.items():
            stats.append({"namespace": ns_, "mapping_status": k, "n_ids": v})
    vc = vend.mapping_status.value_counts()
    for k, v in vc.items():
        stats.append({"namespace": "vendor_matrix_genes(37480)", "mapping_status": k, "n_ids": v, "fraction": v / len(vend)})
    n_v51 = len(v51)
    stats.append({"namespace": "v5.1_genes", "mapping_status": "with_any_vendor_gene", "n_ids": int(wide.vendor_ensembl_r62_ids.notna().sum()),
                  "fraction": wide.vendor_ensembl_r62_ids.notna().mean()})
    stats.append({"namespace": "v5.1_genes", "mapping_status": "total", "n_ids": n_v51})
    pd.DataFrame(stats).to_csv(out / "crosswalk_coverage_stats.tsv", sep="\t", index=False)
    print(pd.DataFrame(stats).to_string(index=False))


if __name__ == "__main__":
    main()
