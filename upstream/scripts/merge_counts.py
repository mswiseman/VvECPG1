#!/usr/bin/env python3
"""Merge per-sample featureCounts outputs into one gene x sample matrix (+ assignment summary).
Column names are sample_uids parsed from file names (<sid>.<mode>.txt). Integer modes are
verified to be integral. Usage: merge_counts.py --out M.tsv --summary-out S.tsv f1 f2 ...
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    cols, summ, meta = {}, {}, None
    for f in a.files:
        sid, mode = Path(f).name.split(".")[:2]
        d = pd.read_csv(f, sep="\t", comment="#", index_col=0)
        if meta is None:
            meta = d[["Chr", "Start", "End", "Strand", "Length"]].copy()
            meta["Chr"] = meta.Chr.str.split(";").str[0]
            meta["Start"] = d.Start.astype(str).str.split(";").map(lambda x: min(map(int, x)))
            meta["End"] = d.End.astype(str).str.split(";").map(lambda x: max(map(int, x)))
            meta["Strand"] = meta.Strand.str.split(";").str[0]
        cols[sid] = d.iloc[:, -1]
        summ[sid] = pd.read_csv(f + ".summary", sep="\t", index_col=0).iloc[:, 0]
    m = pd.DataFrame(cols)
    m = m[sorted(m.columns, key=lambda s: int(s.split("_")[1]))]
    if mode in ("primary", "antisense"):
        assert np.allclose(m.values, np.round(m.values)), "non-integer counts in integer mode"
        m = m.round().astype(int)
    out = meta.join(m)
    out.index.name = "gene_id"
    out.to_csv(a.out + ".tmp", sep="\t")
    os.replace(a.out + ".tmp", a.out)
    s = pd.DataFrame(summ)[m.columns]
    s.to_csv(a.summary_out, sep="\t")


if __name__ == "__main__":
    main()
