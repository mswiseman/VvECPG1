#!/usr/bin/env python3
"""3' bias / gene-end profile for a UMI-deduplicated BAM.

For each primary, uniquely mapped (NH==1) read, find same-strand (sense, forward-stranded
library) genes it overlaps or lies within <window> bp downstream of. Reports:
  * histogram of read 3'-end position relative to the annotated gene 3' end (bp; negative =
    inside gene upstream of the end, positive = downstream of the annotated end)
  * per-gene counts of sense reads that fall 1..window bp DOWNSTREAM of the annotated end
    (i.e. outside the model; lost to exon-based counting) vs inside the gene
Usage: three_prime_profile.py --bam X.bam --gtf G.gtf --sample S --out-prefix P [--window 1000]
"""
from __future__ import annotations

import argparse
import bisect
import collections

import numpy as np
import pandas as pd
import pysam


def gene_ends(gtf: str) -> dict[tuple[str, bool], tuple[list, list, list]]:
    spans = {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) != 9 or f[2] != "exon":
            continue
        gid = f[8].split('gene_id "')[1].split('"')[0]
        s, e = int(f[3]), int(f[4])
        k = (f[0], f[6] == "-", gid)
        cur = spans.get(k)
        spans[k] = (min(s, cur[0]), max(e, cur[1])) if cur else (s, e)
    idx = collections.defaultdict(list)
    for (chrom, rev, gid), (s, e) in spans.items():
        idx[(chrom, rev)].append((s, e, gid))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = ([x[0] for x in v], [x[1] for x in v], [x[2] for x in v])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--window", type=int, default=1000)
    a = ap.parse_args()
    W = a.window
    idx = gene_ends(a.gtf)
    maxlen = {k: max(e - s for s, e in zip(v[0], v[1])) for k, v in idx.items()}
    hist = collections.Counter()
    inside = collections.Counter()
    downstream = collections.Counter()
    n = n_unique = n_none = 0
    with pysam.AlignmentFile(a.bam) as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            n += 1
            if r.get_tag("NH") != 1:
                continue
            n_unique += 1
            rev = r.is_reverse
            key = (r.reference_name, rev)
            if key not in idx:
                n_none += 1
                continue
            starts, ends, gids = idx[key]
            # read 3' end in transcript orientation (forward-stranded: read strand == gene strand)
            p = r.reference_end if not rev else r.reference_start + 1
            lo = bisect.bisect_left(starts, p - maxlen[key] - W)
            hi = bisect.bisect_right(starts, p + W)
            best = None
            for i in range(lo, hi):
                s, e = starts[i], ends[i]
                end = s if rev else e
                d = (end - p) if rev else (p - end)   # >0 downstream of annotated end
                if s - (W if rev else 0) <= p <= e + (0 if rev else W):
                    if best is None or abs(d) < abs(best[0]):
                        best = (d, gids[i])
            if best is None:
                n_none += 1
                continue
            d, gid = best
            hist[max(-5000, min(W, d)) // 50 * 50] += 1
            if d > 0:
                downstream[gid] += 1
            else:
                inside[gid] += 1
    h = pd.Series(hist).sort_index().rename("reads").rename_axis("bin_bp_rel_gene_end")
    h.to_csv(a.out_prefix + ".three_prime_hist.tsv", sep="\t")
    g = pd.DataFrame({"inside": pd.Series(inside), "downstream_1_to_W": pd.Series(downstream)}).fillna(0).astype(int)
    g.index.name = "gene_id"
    g.to_csv(a.out_prefix + ".gene_end_reads.tsv", sep="\t")
    tot = sum(hist.values())
    within = sum(v for k, v in hist.items() if -500 <= k <= 0)
    summ = {"sample_uid": a.sample, "primary_reads": n, "unique_reads": n_unique,
            "unique_near_sense_gene": tot, "unique_no_sense_gene": n_none,
            "frac_within_500bp_upstream_of_end": within / tot if tot else np.nan,
            "frac_downstream_of_annotated_end": sum(v for k, v in hist.items() if k > 0) / tot if tot else np.nan}
    pd.DataFrame([summ]).to_csv(a.out_prefix + ".three_prime_summary.tsv", sep="\t", index=False)
    print(summ)


if __name__ == "__main__":
    main()
