#!/usr/bin/env python3
"""Arabidopsis homologs for PN40024 v5.1 genes (DIAMOND blastp, longest protein per gene).

Reference: Ensembl Plants r62 TAIR10 peptides (references/arabidopsis, CHECKSUMS-verified at download).
Output: results/annotation/v51_arabidopsis_homologs.tsv with, per v5.1 gene, the best Arabidopsis hit
(gene, symbol, description, identity, coverage, e-value) and whether it is a reciprocal best hit (RBH).
Best hits are family-level homology, not proven orthology; RBH is the stricter flag.
Usage: arabidopsis_homologs.py --root <project_root> --threads N
"""
from __future__ import annotations

import argparse
import gzip
import re
import subprocess
from pathlib import Path

import pandas as pd

COLS = "qseqid sseqid pident length qlen slen evalue bitscore"


def sh(cmd: str) -> None:
    print("+", cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    root = Path(a.root)
    E = root / "envs/rnaseq/bin"
    out = root / "results/annotation"
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    final = out / "v51_arabidopsis_homologs.tsv"
    if final.exists():
        print("exists:", final)
        return
    v51 = root / "results/crosswalk/work_chardonnay/v51_longest_protein.fa"
    if not v51.exists():
        raise SystemExit(f"missing {v51}; run make chardonnay-analysis first")

    # longest Arabidopsis protein per gene + gene annotation
    at_fa = work / "at_longest.fa"
    ann = {}
    best: dict[str, str] = {}
    name = None
    seqs: dict[str, list[str]] = {}
    with gzip.open(root / "references/arabidopsis/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                g = re.search(r"gene:(\S+)", line).group(1)
                sym = re.search(r"gene_symbol:(\S+)", line)
                desc = re.search(r"description:(.*?)(\s\[Source.*)?$", line.strip())
                ann[g] = (sym.group(1) if sym else "", desc.group(1) if desc else "")
                name = (g, line[1:].split()[0])
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    for (g, _), parts in seqs.items():
        s = "".join(parts).rstrip("*")
        if g not in best or len(s) > len(best[g]):
            best[g] = s
    if not at_fa.exists():
        with open(at_fa, "w") as o:
            for g, s in best.items():
                o.write(f">{g}\n{s}\n")
    sh(f"{E}/diamond makedb --in {at_fa} -d {work}/at --quiet")
    sh(f"{E}/diamond makedb --in {v51} -d {work}/v51 --quiet")
    fw, rv = work / "v51_to_at.tsv", work / "at_to_v51.tsv"
    sh(f"{E}/diamond blastp --sensitive -q {v51} -d {work}/at -o {fw} -p {a.threads} --max-target-seqs 5 -e 1e-5 --quiet --outfmt 6 {COLS}")
    sh(f"{E}/diamond blastp --sensitive -q {at_fa} -d {work}/v51 -o {rv} -p {a.threads} --max-target-seqs 5 -e 1e-5 --quiet --outfmt 6 {COLS}")
    cols = COLS.split()
    f = pd.read_csv(fw, sep="\t", names=cols)
    r = pd.read_csv(rv, sep="\t", names=cols)
    fb = f.sort_values(["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid")
    rb = r.sort_values(["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid")
    rbest = dict(zip(rb.qseqid, rb.sseqid))
    fb = fb.assign(
        query_cov=lambda d: (d.length / d.qlen).round(3),
        subject_cov=lambda d: (d.length / d.slen).round(3),
        rbh=lambda d: [rbest.get(s) == q for q, s in zip(d.qseqid, d.sseqid)],
        at_symbol=lambda d: [ann.get(s, ("", ""))[0] for s in d.sseqid],
        at_description=lambda d: [ann.get(s, ("", ""))[1] for s in d.sseqid])
    fb = fb.rename(columns={"qseqid": "pn_v5_1_gene_id", "sseqid": "at_gene", "pident": "pident_best"})
    fb = fb[["pn_v5_1_gene_id", "at_gene", "at_symbol", "at_description", "pident_best", "query_cov", "subject_cov",
             "evalue", "bitscore", "rbh"]]
    allg = pd.read_csv(v51, sep="\t", header=None, comment=None)[0]
    allg = allg[allg.str.startswith(">")].str[1:]
    fb = pd.DataFrame({"pn_v5_1_gene_id": allg}).merge(fb, how="left")
    tmp = final.with_suffix(".tmp")
    fb.to_csv(tmp, sep="\t", index=False)
    tmp.rename(final)
    print(f"{fb.at_gene.notna().sum()}/{len(fb)} v5.1 genes with an Arabidopsis hit (e<=1e-5); RBH: {int(fb.rbh.fillna(False).sum())}")


if __name__ == "__main__":
    main()
