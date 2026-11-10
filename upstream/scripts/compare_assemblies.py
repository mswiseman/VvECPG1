#!/usr/bin/env python3
"""Compare assemblies sequence-by-sequence (length + MD5 of uppercase, unwrapped sequence).

Builds 1:1 contig maps ONLY from exact normalized-sequence identity (never from names or
approximate lengths). Also reports soft-mask (lowercase) fraction and N content.
Usage:
  compare_assemblies.py --set NAME=fasta[,fasta...] --set NAME2=... --out-prefix reports/assembly_compare
Outputs <prefix>.sequences.tsv (one row per sequence) and <prefix>.pairs.tsv (exact matches
between sets, plus unmatched sequences).
"""
from __future__ import annotations

import argparse
import hashlib
import itertools

import pandas as pd
import pysam


def digest(fasta: str, set_name: str) -> list[dict]:
    rows = []
    with pysam.FastxFile(fasta) as fh:
        for rec in fh:
            s = rec.sequence
            up = s.upper()
            rows.append({"set": set_name, "file": fasta, "seqid": rec.name,
                         "description": rec.comment or "", "length": len(s),
                         "md5_upper": hashlib.md5(up.encode()).hexdigest(),
                         "frac_lowercase": sum(c.islower() for c in s) / len(s) if s else 0,
                         "n_count": up.count("N")})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", action="append", required=True, help="NAME=fasta[,fasta]")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    rows = []
    sets = []
    for spec in a.set:
        name, files = spec.split("=", 1)
        sets.append(name)
        for f in files.split(","):
            rows += digest(f, name)
    seqs = pd.DataFrame(rows)
    seqs.to_csv(a.out_prefix + ".sequences.tsv", sep="\t", index=False)

    pairs = []
    for s1, s2 in itertools.combinations(sets, 2):
        A = seqs[seqs.set == s1]
        B = seqs[seqs.set == s2]
        m = A.merge(B, on="md5_upper", suffixes=("_a", "_b"))
        for r in m.itertuples():
            pairs.append({"set_a": s1, "seqid_a": r.seqid_a, "set_b": s2, "seqid_b": r.seqid_b,
                          "length": r.length_a, "status": "identical_sequence"})
        for r in A[~A.md5_upper.isin(B.md5_upper)].itertuples():
            same_len = B[B.length == r.length].seqid.tolist()
            pairs.append({"set_a": s1, "seqid_a": r.seqid, "set_b": s2, "seqid_b": ",".join(same_len),
                          "length": r.length,
                          "status": "no_identical_sequence" + ("(same_length_differs)" if same_len else "")})
        for r in B[~B.md5_upper.isin(A.md5_upper)].itertuples():
            pairs.append({"set_a": s1, "seqid_a": "", "set_b": s2, "seqid_b": r.seqid, "length": r.length,
                          "status": "only_in_b"})
    pd.DataFrame(pairs).to_csv(a.out_prefix + ".pairs.tsv", sep="\t", index=False)
    print(pd.DataFrame(pairs).groupby(["set_a", "set_b", "status"]).size())


if __name__ == "__main__":
    main()
