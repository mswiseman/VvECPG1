#!/usr/bin/env python3
"""Figures for the 24 h focal-gene vs embryogenic-programme follow-up (scripts/focal_24h_programme.R, *_boot.R).

Outputs (results/pn40024_v5.1/focal_24h_programme/):
  focal_24h_programme_scatter.{png,pdf,svg}   ten 24 h contrasts: VvECPG1 shrunken log2FC vs programme shift
                                              (error bars = bootstrap 2.5-97.5 percentiles if the bootstrap was run)
  within_group_residuals.{png,pdf,svg}        replicate deviations from group means (technically adjusted matrix):
                                              A treated libraries, all times; B 24 h treated libraries
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.path.join(os.environ.get("VVECPG1_UPSTREAM_ROOT", "."), "results/pn40024_v5.1/focal_24h_programme")
PAL = {"ABA": "#2a78d6", "ACC": "#eb6834", "GA3": "#1baf7a", "MeJA": "#eda100", "NAA": "#e87ba4"}
TXT, TXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": TXT2, "axes.labelcolor": TXT,
                     "xtick.color": TXT2, "ytick.color": TXT2})

ct = pd.read_csv(f"{D}/focal_24h_contrasts.tsv", sep="\t")
S = pd.read_csv(f"{D}/focal_24h_summary.tsv", sep="\t")
row = lambda a, m: S[(S.analysis == a) & S.metric.str.startswith(m)].iloc[0]
bootp = f"{D}/focal_24h_bootstrap_contrasts.tsv"
boot = pd.read_csv(bootp, sep="\t").set_index("contrast") if os.path.exists(bootp) else None
bsum = pd.read_csv(f"{D}/focal_24h_bootstrap_summary.tsv", sep="\t").set_index("metric") if os.path.exists(f"{D}/focal_24h_bootstrap_summary.tsv") else None

# ---------------- 1. 24 h scatter ----------------
fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=200)
fig.subplots_adjust(left=0.13, right=0.97, top=0.93, bottom=0.25)
ax.axhline(0, color=TXT2, lw=0.8, ls=(0, (4, 3)), zorder=1)
ax.grid(color=GRID, lw=0.6, zorder=0)
lowdose = ct.groupby("hormone").concentration.transform("min") == ct.concentration
for r, lo in zip(ct.itertuples(), lowdose):
    x, y = r.focal_log2FC_shrunk_apeglm, r.programme_shift
    if boot is not None:
        b = boot.loc[r.contrast]
        ax.errorbar(x, y, xerr=[[x - b.focal_q025], [b.focal_q975 - x]], yerr=[[y - b.shift_q025], [b.shift_q975 - y]],
                    fmt="none", ecolor=PAL[r.hormone], alpha=0.35, lw=1, zorder=2)
    ax.scatter(x, y, s=58, marker="o" if lo else "s", color=PAL[r.hormone], edgecolor="white", lw=1.2, zorder=3)
# direct labels (hand-tuned offsets to avoid collisions)
off = {"ABA 5 mg/L": (0, -12, "center"), "ABA 7.5 mg/L": (6, -11, "left"), "ACC 1 mg/L": (5, 8, "left"), "ACC 2.5 mg/L": (-4, 9, "right"),
       "GA3 1 mg/L": (-8, 6, "right"), "GA3 2 mg/L": (6, -9, "left"), "MeJA 100 µM": (6, 4, "left"), "MeJA 75 µM": (6, -3, "left"),
       "NAA 10 mg/L": (6, -9, "left"), "NAA 5 mg/L": (6, 5, "left")}
for r in ct.itertuples():
    lab = r.label.replace(", 24 h", "")
    dx, dy, ha = off.get(lab, (6, 4, "left"))
    ax.annotate(lab, (r.focal_log2FC_shrunk_apeglm, r.programme_shift), xytext=(dx, dy), textcoords="offset points",
                fontsize=7.8, color=TXT, ha=ha, va="center", zorder=4)
p, sp = row("24 h contrasts", "Pearson r, focal shrunken LFC vs programme shift"), row("24 h contrasts", "Spearman")
lo_ = S[S.analysis == "24 h contrasts, leave one out"].value
txt = (f"Pearson r = {p.value:.2f}; Spearman ρ = {sp.value:.2f} (10 contrasts)\n"
       f"Random gene sets, same contrasts: r ≥ {p.value:.2f} in {100 * (p.p_upper * 1001 - 1) / 1000:.0f}% (median {p.null_median:.2f})\n"
       f"Leave one contrast out: r = {lo_.min():.2f} to {lo_.max():.2f}")
if bsum is not None:
    b = bsum.loc["pearson_r_shift"]
    txt += f"\nBootstrap ({int(b.n_boot)} refits) 95% interval for r: {b.q025:.2f} to {b.q975:.2f}"
ax.text(0.015, 0.985, txt, transform=ax.transAxes, fontsize=7.6, color=TXT2, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=GRID, lw=0.8))
ax.set_xlabel("VvECPG1 (Vitvi05_01chr01g05750) log2 fold change vs time-matched untreated control\n(DESeq2, apeglm-shrunken; 24 h, batch 1)", fontsize=8.6)
ax.set_ylabel("Embryogenic programme shift (log2)\nsignature mean shrunken log2FC − random-set median", fontsize=8.6)
ax.set_title("24 h hormone treatments vs untreated callus", loc="left", fontsize=10, color=TXT, fontweight="bold")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
hz = [plt.Line2D([], [], marker="o", lw=0, color=v, markersize=7, label=k) for k, v in PAL.items()]
hz += [plt.Line2D([], [], marker="o", lw=0, color=TXT2, markersize=6, label="lower dose"),
       plt.Line2D([], [], marker="s", lw=0, color=TXT2, markersize=6, label="higher dose")]
ax.legend(handles=hz, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=7, frameon=False, fontsize=7.6, handletextpad=0.2, columnspacing=1)
cap = ("Programme = 880 of the 964 atlas embryogenic-signature genes (higher in embryogenic vs non-embryogenic callus and in early vs late "
       "somatic embryos; VvECPG1 excluded). Random sets: 1,000 sets matched on expression decile. "
       + ("Bars: bootstrap 2.5–97.5 percentiles (controls resampled once per refit, shared by all ten contrasts)." if boot is not None else ""))
fig.text(0.13, 0.012, cap, fontsize=6.6, color=TXT2, wrap=True)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{D}/focal_24h_programme_scatter.{ext}", facecolor="white")
plt.close(fig)

# ---------------- 2. within-group residuals ----------------
R = pd.read_csv(f"{D}/within_group_residuals.tsv", sep="\t")
R = R[R.matrix == "adjusted"]
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4), dpi=200)
fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.2, wspace=0.25)
for ax, (sub, ttl) in zip(axes, [("treated only (controls excluded)", "A. Treated libraries, all times"),
                                 ("24 h treated groups only", "B. Treated libraries, 24 h")]):
    d = R[R.subset == sub]
    ax.axhline(0, color=GRID, lw=0.8, zorder=0); ax.axvline(0, color=GRID, lw=0.8, zorder=0)
    for h, g in d.groupby("hormone"):
        ax.scatter(g.focal_resid, g.programme_resid, s=22, color=PAL[h], edgecolor="white", lw=0.6, zorder=3, label=h)
    s = row("within-group, adjusted", "Pearson r of residuals from group means: " + sub)
    ax.text(0.02, 0.97, f"r = {s.value:.2f}  ({int(s.n_libraries)} libraries, {int(s.n_groups)} groups, residual df {int(s.residual_df)})\n"
            f"random gene sets: median {s.null_median:.2f}, 97.5th percentile {s.null_q975:.2f}", transform=ax.transAxes,
            fontsize=7.4, color=TXT2, va="top")
    ax.set_title(ttl, loc="left", fontsize=9.5, fontweight="bold", color=TXT)
    ax.set_xlabel("VvECPG1, deviation from group mean (z, technically adjusted)", fontsize=8.2)
    ax.set_ylabel("Programme score, deviation from group mean\n(mean z of signature genes)", fontsize=8.2)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
axes[0].legend(loc="lower right", frameon=False, fontsize=7.4, handletextpad=0.2)
fig.text(0.08, 0.02, "Each point is one library. Residuals are taken within each hormone × dose × time group; the correlation is pooled over all residuals "
         "(not an average of per-group correlations). Untreated controls excluded. Covariation among replicate cultures, not a hormone effect.",
         fontsize=6.6, color=TXT2)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{D}/within_group_residuals.{ext}", facecolor="white")
print("figures written")
