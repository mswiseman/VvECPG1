#!/usr/bin/env python3
"""Strand-aware coverage at the focal locus (chr01:4,083,500-4,088,500) from UMI-deduplicated,
uniquely mapped (NH==1) primary alignments, normalized to molecules per million assigned (DESeq2
size factors are used for group summaries elsewhere). Tracks: v5.1 models, vendor Ensembl r62
models, then coverage for the 6 controls + the samples with highest/lowest focal counts.
Writes PDF/SVG/PNG figure, per-sample stranded bedGraphs and an IGV session for manual review.
Usage: focal_locus_plot.py --root <project_root>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
import pysam

CHR, START, END = "chr01", 4083500, 4088500
FOCAL = "Vitvi05_01chr01g05750"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
SENSE, ANTI = "#2a78d6", "#eb6834"   # categorical slots 1-2 (validated reference palette)


def models(gtf: str, chrom: str, lo: int, hi: int) -> dict:
    out = {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) != 9 or f[2] != "exon" or f[0] != chrom:
            continue
        s, e = int(f[3]), int(f[4])
        if e < lo or s > hi:
            continue
        gid = f[8].split('gene_id "')[1].split('"')[0]
        out.setdefault(gid, {"strand": f[6], "exons": []})["exons"].append((s, e))
    return out


def stranded_cov(bam: str, chrom: str, lo: int, hi: int):
    plus = np.zeros(hi - lo + 1)
    minus = np.zeros(hi - lo + 1)
    with pysam.AlignmentFile(bam) as b:
        for r in b.fetch(chrom, lo, hi):
            if r.is_secondary or r.is_supplementary or r.get_tag("NH") != 1:
                continue
            arr = minus if r.is_reverse else plus
            for bs, be in r.get_blocks():
                a, z = max(bs + 1, lo), min(be, hi)
                if a <= z:
                    arr[a - lo:z - lo + 1] += 1
    return plus, minus


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    out = root / "results/pn40024_v5.1/target_gene"
    (out / "bedgraph").mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(root / "metadata/samples.tsv", sep="\t").set_index("sample_uid")
    cm = pd.read_csv(root / "results/pn_v51/counts/pn_v51.primary.counts.tsv", sep="\t", index_col=0)
    cnt = cm.loc[FOCAL, [c for c in cm.columns if c.startswith("PPPJ8P_")]].astype(float)
    lib = cm[[c for c in cm.columns if c.startswith("PPPJ8P_")]].sum()
    cpm = cnt / lib * 1e6
    ctrls = list(s[s.is_control].index)
    treated = cpm.drop(ctrls).sort_values()
    show = ctrls + list(treated.index[-3:][::-1]) + list(treated.index[:3])
    v51 = models(str(root / "references/pn40024_v5.1/build/pn_v51.gtf"), CHR, START, END)
    ens_map = dict(pd.read_csv(root / "references/pn40024_v5.1/build/seqmap_ensembl_to_fasta.tsv", sep="\t")[["fasta_seqid", "ensembl_seqid"]].values)
    ens = models(str(root / "references/ensembl_r62/build/ensembl_r62.gtf"), ens_map[CHR], START, END)
    x = np.arange(START, END + 1)
    fig, axes = plt.subplots(len(show) + 2, 1, figsize=(10, 1.0 * (len(show) + 2) + 1), sharex=True,
                             gridspec_kw={"height_ratios": [0.8, 0.8] + [1] * len(show)})
    fig.patch.set_facecolor(SURF)
    for ax, (title, mdl) in zip(axes[:2], (("PN40024 v5.1", v51), ("Vendor Ensembl r62", ens))):
        for i, (gid, m) in enumerate(sorted(mdl.items(), key=lambda kv: min(e[0] for e in kv[1]["exons"]))):
            y = i % 2
            lo_, hi_ = min(e[0] for e in m["exons"]), max(e[1] for e in m["exons"])
            ax.plot([max(lo_, START), min(hi_, END)], [y, y], color=INK2, lw=0.8)
            for es, ee in m["exons"]:
                ax.add_patch(plt.Rectangle((es, y - 0.3), ee - es, 0.6, color=SENSE if gid == FOCAL else INK2, lw=0))
            ax.text(max(lo_, START), y + 0.42, f"{gid} ({m['strand']})", fontsize=6, color=INK)
        ax.set_ylim(-0.6, 1.9)
        ax.set_yticks([])
        ax.set_ylabel(title, rotation=0, ha="right", va="center", fontsize=7, color=INK)
        for sp in ax.spines.values():
            sp.set_visible(False)
    bams = {sid: root / f"results/pn_v51/bam/{sid}.dedup.bam" for sid in show}
    ymax = 0
    covs = {}
    for sid in show:
        p, m = stranded_cov(str(bams[sid]), CHR, START, END)
        scale = 1e6 / lib[sid]
        covs[sid] = (p * scale, m * scale)
        ymax = max(ymax, covs[sid][0].max(), covs[sid][1].max())
        with open(out / "bedgraph" / f"{sid}.focal_locus.plus.bedGraph", "w") as fp, \
                open(out / "bedgraph" / f"{sid}.focal_locus.minus.bedGraph", "w") as fm:
            for arr, fh in ((covs[sid][0], fp), (covs[sid][1], fm)):
                fh.write(f'track type=bedGraph name="{sid} {s.loc[sid, "vendor_sample_name"]}"\n')
                for i in np.flatnonzero(arr):
                    fh.write(f"{CHR}\t{START + i - 1}\t{START + i}\t{arr[i]:.4f}\n")
    for ax, sid in zip(axes[2:], show):
        p, m = covs[sid]
        ax.fill_between(x, 0, p, color=SENSE, lw=0, step="mid")
        ax.fill_between(x, 0, -m, color=ANTI, lw=0, step="mid")
        ax.axhline(0, color=GRID, lw=0.5)
        ax.set_ylim(-ymax * 0.15, ymax * 1.05)
        tag = "control" if sid in ctrls else ("high" if sid in treated.index[-3:] else "low")
        ax.set_ylabel(f"{s.loc[sid, 'vendor_sample_name']}\n{sid.split('_')[1]} b{s.loc[sid, 'batch']} [{tag}]\n"
                      f"{cnt[sid]:.0f} mol; {cpm[sid]:.0f} CPM", rotation=0, ha="right", va="center", fontsize=6, color=INK)
        ax.tick_params(labelsize=6, colors=INK2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[-1].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    axes[-1].set_xlabel(f"{CHR} position (bp); blue = + strand (sense for {FOCAL}), orange = − strand (plotted below 0); "
                        "molecules per million assigned", fontsize=7, color=INK)
    fig.suptitle(f"Focal locus {FOCAL} (aliases Vitvi01g01965 / VIT_01s0011g04570): UMI-deduplicated unique coverage",
                 fontsize=9, color=INK)
    fig.tight_layout(rect=(0.0, 0, 1, 0.97))
    for ext in ("pdf", "svg"):
        fig.savefig(out / f"focal_locus_coverage.{ext}")
    fig.savefig(out / "focal_locus_coverage.png", dpi=300)
    # IGV session
    tracks = "\n".join(
        f'    <Resource path="bedgraph/{sid}.focal_locus.{st}.bedGraph"/>' for sid in show for st in ("plus", "minus"))
    (out / "focal_locus_igv_session.xml").write_text(
        f'<?xml version="1.0" encoding="UTF-8"?>\n<Session genome="../../../references/pn40024_v5.1/build/genome_nuc_org.fa" '
        f'locus="{CHR}:{START}-{END}" version="8">\n  <Resources>\n'
        f'    <Resource path="../../../references/pn40024_v5.1/build/pn_v51.gtf"/>\n{tracks}\n  </Resources>\n</Session>\n')
    pd.DataFrame({"sample_uid": cnt.index, "focal_molecules": cnt.values, "focal_CPM": cpm.values}).to_csv(
        out / "focal_raw_counts_cpm.tsv", sep="\t", index=False)
    print("focal locus plot done;", "shown:", show)


if __name__ == "__main__":
    main()
