#!/usr/bin/env python3
"""Collect per-sample QC into one row: fastp, STAR, UMIcollapse, idxstats (organelles),
featureCounts assignment (primary unique, antisense check), 3' profile summary.
Usage: sample_qc.py --sample S --genome G --fastp J --star-log L --dedup-log U --idxstats I
                    --fc-primary P.summary --fc-antisense A.summary --three-prime T --out O
"""
from __future__ import annotations

import argparse
import json
import re

import pandas as pd


def star_log(path: str) -> dict:
    d = {}
    for line in open(path):
        if "|" in line:
            k, v = [x.strip() for x in line.split("|", 1)]
            d[k] = v
    pct = lambda k: float(d[k].rstrip("%")) / 100
    return {"star_input_reads": int(d["Number of input reads"]),
            "star_unique": int(d["Uniquely mapped reads number"]),
            "star_unique_frac": pct("Uniquely mapped reads %"),
            "star_multi": int(d["Number of reads mapped to multiple loci"]),
            "star_multi_frac": pct("% of reads mapped to multiple loci"),
            "star_toomany_frac": pct("% of reads mapped to too many loci"),
            "star_unmapped_short_frac": pct("% of reads unmapped: too short"),
            "star_unmapped_other_frac": pct("% of reads unmapped: other"),
            "star_splices_per_read": int(d["Number of splices: Total"]) / max(1, int(d["Uniquely mapped reads number"]))}


def fc_summary(path: str, prefix: str) -> dict:
    s = pd.read_csv(path, sep="\t", index_col=0).iloc[:, 0]
    tot = s.sum()
    return {f"{prefix}_assigned": int(s["Assigned"]),
            f"{prefix}_assigned_frac": s["Assigned"] / tot if tot else 0,
            f"{prefix}_ambiguity": int(s.get("Unassigned_Ambiguity", 0)),
            f"{prefix}_multimapping": int(s.get("Unassigned_MultiMapping", 0)),
            f"{prefix}_nofeatures": int(s.get("Unassigned_NoFeatures", 0))}


def main() -> None:
    ap = argparse.ArgumentParser()
    for k in ("sample", "genome", "fastp", "star_log", "dedup_log", "idxstats", "fc_primary",
              "fc_antisense", "three_prime", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args()
    row = {"sample_uid": a.sample, "genome": a.genome}
    fp = json.load(open(a.fastp))
    row.update({"raw_reads": fp["summary"]["before_filtering"]["total_reads"],
                "fastp_pass_reads": fp["summary"]["after_filtering"]["total_reads"],
                "fastp_pass_frac": fp["summary"]["after_filtering"]["total_reads"] /
                fp["summary"]["before_filtering"]["total_reads"],
                "fastp_mean_len_after": fp["summary"]["after_filtering"]["read1_mean_length"]})
    row.update(star_log(a.star_log))
    txt = open(a.dedup_log).read()
    m_in = re.search(r"Number of input reads\s+(\d+)", txt)
    m_out = re.search(r"Number of reads after deduplicating\s+(\d+)", txt)
    if m_in and m_out:
        row.update({"dedup_input_alignments": int(m_in.group(1)), "dedup_output_alignments": int(m_out.group(1)),
                    "dedup_retained_frac": int(m_out.group(1)) / int(m_in.group(1))})
    idx = pd.read_csv(a.idxstats, sep="\t", header=None, names=["seq", "len", "mapped", "unmapped"])
    tot = idx.mapped.sum()
    row.update({"dedup_mapped_alignments": int(tot),
                "frac_MT": idx[idx.seq == "MT"].mapped.sum() / tot if tot else 0,
                "frac_chloroplast": idx[idx.seq.str.startswith("chl_")].mapped.sum() / tot if tot else 0})
    row.update(fc_summary(a.fc_primary, "fc_sense"))
    anti = fc_summary(a.fc_antisense, "fc_antisense")
    row["fc_antisense_assigned"] = anti["fc_antisense_assigned"]
    row["strand_sense_frac"] = row["fc_sense_assigned"] / max(1, row["fc_sense_assigned"] + anti["fc_antisense_assigned"])
    tp = pd.read_csv(a.three_prime, sep="\t").iloc[0].to_dict()
    row.update({k: v for k, v in tp.items() if k.startswith("frac_")})
    pd.DataFrame([row]).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
