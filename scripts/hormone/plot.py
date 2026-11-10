#!/usr/bin/env python3
"""Focal gene (Vitvi05_01chr01g05750) at 24 h: each treatment vs untreated controls.

Controls (PPPJ8P_94-96, batch 1) were sampled at 24 h (user, 2026-11-06), so the
24 h contrasts are time-matched. Bars = mean DESeq2-normalised counts (batch-1 size
factors) +/- SE, n = 3; dots = replicates. Significance = DESeq2 Wald test from the
existing within-batch-1 model (de/tables/*_24h_vs_control_b1.all_genes.tsv), BH-adjusted
across the 10 treatment comparisons (the focal gene is pre-specified). The genome-wide
BH padj is reported in the TSV too.
"""
import glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def holm(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.maximum.accumulate(p[o] * (n - np.arange(n)))
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

from pathlib import Path
PACKAGE = Path(__file__).resolve().parents[2]
ROOT = str(PACKAGE / "source_addenda/C_de_tables")
FOCAL = "Vitvi05_01chr01g05750"
OUTDIR = PACKAGE / "results/hormone"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT = str(OUTDIR / "focal_24h_vs_control")

# --- stats -----------------------------------------------------------------
rows = []
for f in sorted(glob.glob(f"{ROOT}/de_tables/*_24h_vs_control_b1.all_genes.tsv.gz")):
    d = pd.read_csv(f, sep="\t")
    r = d.loc[d.gene_id == FOCAL].iloc[0]
    rows.append(dict(condition=r.condition, label=r.label, log2FC=r.log2FC, lfcSE=r.lfcSE,
                     log2FC_shrunk=r.log2FC_shrunk_apeglm, pvalue=r.pvalue,
                     padj_genomewide=r.padj_BH, control_samples=r.control_samples))
st = pd.DataFrame(rows)
assert len(st) == 10
st["padj_10_tests_BH"] = bh(st.pvalue)
st["padj_10_tests_holm"] = holm(st.pvalue)

# --- per-sample normalised counts -----------------------------------------
ex = pd.read_csv(f"{ROOT}/focal/focal_normalized_expression_per_sample.tsv", sep="\t")
ctrl_ids = st.control_samples.iloc[0].split(",")
ctrl = ex[ex.sample_uid.isin(ctrl_ids)].assign(condition="control")
trt = ex[(ex.batch == 1) & (ex.time_hours.astype(str) == "24") & (~ex.is_control)]
assert len(ctrl) == 3 and len(trt) == 30

hormones = ["ABA", "ACC", "GA3", "MeJA", "NAA"]
# order: control, then hormone x increasing dose
order = ["control"] + list(
    trt.drop_duplicates("condition").sort_values(["hormone", "dose_value"])
       .loc[lambda x: x.hormone.isin(hormones)].condition)
lab = {"control": "Untreated\ncontrol"}
for c, h, dv, du in trt.drop_duplicates("condition")[["condition", "hormone", "dose_value", "dose_unit"]].values:
    lab[c] = f"{h}\n{dv:g} {du.replace('uM', 'µM')}"
hmap = dict(trt.drop_duplicates("condition")[["condition", "hormone"]].values)

allx = pd.concat([ctrl, trt])
summ = allx.groupby("condition").norm.agg(["mean", "std", "count"]).loc[order]
summ["se"] = summ["std"] / np.sqrt(summ["count"])

out = summ.reset_index().merge(st.drop(columns="control_samples"), on="condition", how="left")
out.to_csv(f"{OUT}.tsv", sep="\t", index=False)

# --- plot ------------------------------------------------------------------
pal = {"ABA": "#2a78d6", "ACC": "#eb6834", "GA3": "#1baf7a", "MeJA": "#eda100", "NAA": "#e87ba4"}
CTRL = "#8a8986"
TXT, TXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=200)
x = np.arange(len(order))
for i, c in enumerate(order):
    col = CTRL if c == "control" else pal[hmap[c]]
    m, se = summ.loc[c, "mean"], summ.loc[c, "se"]
    ax.bar(i, m, width=0.72, color=col, edgecolor="white", linewidth=1.5, zorder=2)
    ax.errorbar(i, m, yerr=se, color=TXT, capsize=4, lw=1.2, capthick=1.2, zorder=3)
    pts = allx.loc[allx.condition == c, "norm"].values
    jit = np.linspace(-0.16, 0.16, len(pts))
    ax.scatter(i + jit, pts, s=16, facecolor="white", edgecolor=TXT, lw=0.8, zorder=4)
    if c != "control":
        r = st.set_index("condition").loc[c]
        top = max(m + se, pts.max())
        s = stars(r.padj_10_tests_BH)
        ax.text(i, top * 1.03, s, ha="center", va="bottom",
                fontsize=11 if s != "ns" else 8, color=TXT if s != "ns" else TXT2,
                fontweight="bold" if s != "ns" else "normal")

cm = summ.loc["control", "mean"]
ax.axhline(cm, color=CTRL, lw=1, ls=(0, (4, 3)), zorder=1)
ax.set_xticks(x, [lab[c] for c in order], fontsize=8, color=TXT)
ax.tick_params(axis="x", length=0, pad=4)
ax.set_ylabel("Normalised counts (DESeq2)", color=TXT, fontsize=9)
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1000:g}k"))
ax.tick_params(axis="y", labelsize=8, colors=TXT2)
ax.set_ylim(0, allx.norm.max() * 1.15)
ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
for s in ["top", "right", "left"]:
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(TXT2)
# No title/footnote: caption carries the description (user, 2026-11-06).
fig.tight_layout()
for ext in ["png", "pdf", "svg"]:
    fig.savefig(f"{OUT}_bar.{ext}", facecolor="white")
print(out[["condition", "mean", "se", "log2FC", "pvalue", "padj_10_tests_BH",
           "padj_10_tests_holm", "padj_genomewide"]].to_string(index=False))
