"""Rung 7 (issue #7) figures from the CSV files of data-raw/rung7-power-summary.R.

    python data-raw/rung7-power-figure.py [resultsDir]

Writes fig-rung7-beta.png, fig-rung7-recovery.png and fig-rung7-calibration.png
into resultsDir (default data-raw/rung7-results).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

res = sys.argv[1] if len(sys.argv) > 1 else "data-raw/rung7-results"
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlelocation": "left", "savefig.dpi": 300,
})
# One colour per fossil-density x size combination, threaded through all figures
combos = [("low", "small"), ("high", "small"), ("low", "large"), ("high", "large")]
colours = {("low", "small"): "#9ecae1", ("high", "small"): "#3182bd",
           ("low", "large"): "#fdae6b", ("high", "large"): "#e6550d"}
labels = {c: f"{c[1]} tree, {c[0]} \u03c8" for c in combos}
nChars = [10, 30, 100]


def letter(ax, s):
    ax.text(-0.18, 1.08, s, transform=ax.transAxes, fontweight="bold",
            fontsize=10, va="bottom")


def rate_panel(ax, d, y, lo, hi, title):
    for k, c in enumerate(combos):
        x = d[(d.psiLevel == c[0]) & (d["size"] == c[1])].sort_values("nChar")
        if not len(x):
            continue
        off = 1 + (k - 1.5) * 0.04
        ax.errorbar(x.nChar * off, x[y], yerr=[x[y] - x[lo], x[hi] - x[y]],
                    color=colours[c], marker="o", ms=3.5, lw=1.2, capsize=0,
                    elinewidth=0.8, label=labels[c])
    ax.set_xscale("log")
    ax.set_xticks(nChars)
    ax.set_xticklabels([str(n) for n in nChars])
    ax.minorticks_off()
    ax.set_xlim(8, 125)
    ax.set_ylim(-0.03, 1.03)
    ax.set_title(title)
    ax.set_xlabel("characters")


# ---- Figure 1: the (1,1) vs (0,1) contrast ---------------------------------
con = pd.read_csv(os.path.join(res, "rung7-contrasts.csv"))
runs = pd.read_csv(os.path.join(res, "rung7-runs.csv"))
b1 = con[con.contrast.str.startswith("beta | gamma = 1")]
fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.5), gridspec_kw={"wspace": 0.38})
for ax, truth, ttl in [(axes[0], "(1,1)", "True (1,1): log BF > 0"),
                       (axes[1], "(0,1)", "True (0,1): log BF < 0")]:
    rate_panel(ax, b1[b1.truth == truth], "pRight", "lo", "hi", ttl)
axes[0].set_ylabel("fraction of datasets")
axes[1].set_yticklabels([])
axes[0].legend(frameon=False, loc="lower right", handlelength=1.2)
ax = axes[2]
d = runs[runs.truth.isin(["(1,1)", "(0,1)"]) & (runs.nChar == 100)]
for truth, col, mk in [("(1,1)", "#3182bd", "o"), ("(0,1)", "#636363", "^")]:
    x = d[d.truth == truth]
    ax.scatter(x.s, x.logBFbeta1, s=9, color=col, marker=mk, alpha=0.75,
               lw=0, label=f"true {truth}")
ax.axhline(0, color="#999999", lw=0.7)
ax.set_xscale("log")
ax.set_xlabel("true s (stationary sd of log-tempo)")
ax.set_ylabel("log BF, (1,1) vs (0,1)")
ax.set_title("Evidence for \u03b2 grows with s (100 characters)")
ax.legend(frameon=False, loc="upper left", handletextpad=0.2)
for a, s in zip(axes, "abc"):
    letter(a, s)
fig.savefig(os.path.join(res, "fig-rung7-beta.png"), bbox_inches="tight")
plt.close(fig)

# ---- Figure 2: recovery of the true cell -------------------------------------
rec = pd.read_csv(os.path.join(res, "rung7-recovery.csv"))
fig, axes = plt.subplots(1, 4, figsize=(7.1, 2.2), sharey=True,
                         gridspec_kw={"wspace": 0.08})
for ax, truth in zip(axes, ["(0,0)", "(0,1)", "(1,0)", "(1,1)"]):
    rate_panel(ax, rec[rec.truth == truth], "pCorrect", "lo", "hi",
               f"True cell {truth}")
    ax.axhline(0.25, color="#999999", lw=0.7, ls=":")
axes[0].set_ylabel("fraction selecting the true cell")
axes[0].text(9, 0.27, "chance", fontsize=6, color="#777777")
axes[3].legend(frameon=False, loc="lower right", handlelength=1.2)
fig.savefig(os.path.join(res, "fig-rung7-recovery.png"), bbox_inches="tight")
plt.close(fig)

# ---- Figure 3: calibration -----------------------------------------------------
cal = pd.read_csv(os.path.join(res, "rung7-calibration.csv"))
avg = pd.read_csv(os.path.join(res, "rung7-averaged.csv"))
fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.5), gridspec_kw={"wspace": 0.4})
ax = axes[0]
ax.plot([0, 1], [0, 1], color="#999999", lw=0.7)
for lab, col, mk in [("size-corrected", "#3182bd", "o"), ("uniform prior", "#636363", "s")]:
    x = cal[cal.probabilities == lab]
    ax.errorbar(x.meanP, x.observed, yerr=[x.observed - x.lo, x.hi - x.observed],
                color=col, marker=mk, ms=3.5, lw=0, elinewidth=0.8, label=lab)
ax.set_xlabel("posterior model probability (bin mean)")
ax.set_ylabel("fraction that is the true cell")
ax.set_title("Posterior model probabilities are calibrated")
ax.legend(frameon=False, loc="upper left")
ax.set_xlim(-0.03, 1.03)
ax.set_ylim(-0.03, 1.03)
ax = axes[1]
cellNames = ["(0,0)", "(0,1)", "(1,0)", "(1,1)"]
for k, n in enumerate(nChars):
    x = avg[avg.nChar == n].set_index("cell").loc[cellNames]
    pos = np.arange(4) + (k - 1) * 0.22
    ax.errorbar(pos, x.meanPCorrected, yerr=1.96 * x.se, fmt="o", ms=3.5,
                color=["#c6dbef", "#6baed6", "#08519c"][k], elinewidth=0.8,
                label=f"{n} characters")
ax.axhline(0.25, color="#999999", lw=0.7)
ax.set_xticks(range(4))
ax.set_xticklabels(cellNames)
ax.set_xlabel("cell")
ax.set_ylabel("mean posterior probability")
ax.set_title("Data-averaged posterior = 1/4")
ax.legend(frameon=False, loc="upper right", fontsize=6)
for a, s in zip(axes, "ab"):
    letter(a, s)
fig.savefig(os.path.join(res, "fig-rung7-calibration.png"), bbox_inches="tight")
plt.close(fig)
print("figures written to", res)
