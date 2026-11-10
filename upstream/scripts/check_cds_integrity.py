#!/usr/bin/env python3
"""Verify that GFF CDS coordinates fit the genome: translate every mRNA's CDS and report
start codon, terminal stop, internal stops and length%3. A high pass rate is direct evidence
that the annotation belongs to this exact sequence (wrong builds give ~random codons).
Usage: check_cds_integrity.py --gff X.gff3 --fasta genome.fa --seqmap map.tsv --out out.tsv
"""
from __future__ import annotations

import argparse
import collections

import pandas as pd
import pysam

CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(bases):
    for j, b in enumerate(bases):
        for k, c in enumerate(bases):
            CODON[a + b + c] = aas[16 * i + 4 * j + k]
COMP = str.maketrans("ACGTN", "TGCAN")


def translate(s: str) -> str:
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--seqmap", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    smap = dict(pd.read_csv(a.seqmap, sep="\t").values)
    fa = pysam.FastaFile(a.fasta)
    cds = collections.defaultdict(list)
    info = {}
    for line in open(a.gff):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) != 9 or f[2] != "CDS":
            continue
        for p in [kv.split("=", 1)[1] for kv in f[8].split(";") if kv.startswith("Parent=")][0].split(","):
            cds[p].append((int(f[3]), int(f[4]), int(f[7]) if f[7] != "." else 0))
            info[p] = (f[0], f[6])
    rows = []
    for tx, segs in cds.items():
        seqid, strand = info[tx]
        if seqid not in smap:
            continue
        segs.sort()
        s = "".join(fa.fetch(smap[seqid], st - 1, en).upper() for st, en, _ in segs)
        if strand == "-":
            s = s.translate(COMP)[::-1]
        first_phase = (segs[0] if strand == "+" else segs[-1])[2]
        s = s[first_phase:]
        p = translate(s)
        rows.append({"transcript_id": tx, "len_mod3": len(s) % 3, "starts_ATG": s[:3] == "ATG",
                     "ends_stop": p.endswith("*"), "internal_stops": p[:-1].count("*")})
    df = pd.DataFrame(rows)
    df.to_csv(a.out, sep="\t", index=False)
    n = len(df)
    print(f"transcripts={n} ATG_start={df.starts_ATG.mean():.4f} terminal_stop={df.ends_stop.mean():.4f} "
          f"no_internal_stop={(df.internal_stops == 0).mean():.4f} len_mod3_0={(df.len_mod3 == 0).mean():.4f}")


if __name__ == "__main__":
    main()
