"""Figures for the full-size rung 6 SBC (issue #21).

Reads the CSVs written by data-raw/rung6-sbc-full-summary.R and writes
fig-rung6-full-ecdf.png (ECDF differences with simultaneous 95% bands) and
fig-rung6-full-ranks.png (rank histograms with pointwise 99% binomial bands).
Usage: python data-raw/rung6-sbc-full-figure.py [dir]
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

d = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "data-raw/rung6-full-results"
try:
    apply_figure_style()  # noqa: F821 (provided by the figure-style kernel helpers)
except NameError:
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 6, "ytick.labelsize": 6, "savefig.dpi": 300,
                         "axes.spines.top": False, "axes.spines.right": False})

ranks = pd.read_csv(f"{d}/rung6-full-ranks.csv")
ecdf = pd.read_csv(f"{d}/rung6-full-ecdf.csv")
band = pd.read_csv(f"{d}/rung6-full-ecdf-band.csv")
tests = pd.read_csv(f"{d}/rung6-full-tests.csv")
order = ["lambda", "mu", "psi", "theta", "s", "rate", "rootBin", "clock", "treeLength",
         "meanNodeAge", "nAncestor", "oldestFossil", "youngestFossil"]
label = {"lambda": "λ", "mu": "μ", "psi": "ψ", "theta": "θ", "s": "s", "rate": "rate c",
         "rootBin": "root bin", "clock": "clock ∑∫eˣdt", "treeLength": "tree length",
         "meanNodeAge": "mean node age", "nAncestor": "sampled ancestors",
         "oldestFossil": "oldest fossil age", "youngestFossil": "youngest fossil age"}
short = dict(label, clock="clock", treeLength="tree length", meanNodeAge="node age",
             nAncestor="SAs", oldestFossil="oldest fossil", youngestFossil="youngest fossil")
cells = [(0, 0), (0, 1), (1, 0), (1, 1)]
# Okabe-Ito: CVD-safe, one hue per (beta, gamma) cell, used in both figures
colour = {1: "#0072B2", 2: "#E69F00", 3: "#009E73", 4: "#CC79A7"}
cellName = {i + 1: f"β={b}, γ={g}" for i, (b, g) in enumerate(cells)}
nPerCell = ranks.groupby("cell")["g"].nunique()
nFail = int((~tests["pass"]).sum())

# ECDF differences
fig, axes = plt.subplots(3, 5, figsize=(7.2, 4.6), sharex=True, sharey=True)
axes = axes.ravel()
for ax, q in zip(axes, order):
    sub = ecdf[ecdf.quantity == q]
    nMin = sub.n.min()
    b = band[band.n == nMin]
    ax.fill_between(b.z, b.lo, b.hi, color="0.88", lw=0, zorder=0)
    ax.axhline(0, color="0.5", lw=0.5, zorder=1)
    for cell in range(1, 5):
        s = sub[sub.cell == cell]
        ax.step(s.z, s["diff"], where="post", color=colour[cell], lw=0.9, zorder=2)
        t = tests[(tests.quantity == q) & (tests.cell == cell)]
        if len(t) and not bool(t["pass"].iloc[0]):
            ax.plot(s.z.iloc[s["diff"].abs().argmax()], s["diff"].abs().max() *
                    np.sign(s["diff"].iloc[s["diff"].abs().argmax()]), "o", ms=3,
                    mfc="none", mec="k", zorder=3)
    ax.set_title(label[q], loc="left")
    ax.margins(x=0.02)
for ax in axes[len(order):]:
    ax.axis("off")
handles = [plt.Line2D([], [], color=colour[c], lw=1.2) for c in range(1, 5)]
handles.append(plt.Rectangle((0, 0), 1, 1, color="0.88"))
axes[-1].legend(handles, [cellName[c] for c in range(1, 5)] +
                ["simultaneous 95% band\n(smallest n of the panel)"],
                loc="center", frameon=False)
for ax in axes[10:13]:
    ax.set_xlabel("fractional rank")
for ax in axes[[0, 5, 10]]:
    ax.set_ylabel("ECDF − uniform")
verdict = ("every quantity stays within its band" if nFail == 0 else
           f"{nFail} of {len(tests)} cell × quantity tests fail after Holm correction")
nText = (f"{nPerCell.min()}" if nPerCell.min() == nPerCell.max() else
         f"{nPerCell.min()}–{nPerCell.max()}")
fig.suptitle(f"Rank ECDF differences ({nText} replicates per cell): "
             f"{verdict}", x=0.01, ha="left", fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(f"{d}/fig-rung6-full-ecdf.png")

# Rank histograms, 10 bins, pointwise 99% binomial band
nRank = 100  # ranks 0..99
nBin = 10
fig, axes = plt.subplots(4, len(order), figsize=(12, 4.8), sharex=True, sharey="row")
for i in range(4):
    for j, q in enumerate(order):
        ax = axes[i, j]
        rk = ranks[(ranks.cell == i + 1) & (ranks.quantity == q)]["rank"].dropna()
        n = len(rk)
        counts = np.bincount((rk // (nRank // nBin)).astype(int), minlength=nBin)
        lo, hi = stats.binom.ppf([0.005, 0.995], n, 1 / nBin)
        ax.axhspan(lo, hi, color="0.88", lw=0, zorder=0)
        ax.axhline(n / nBin, color="0.5", lw=0.5, zorder=1)
        ax.bar(np.arange(nBin), counts, width=0.85, color=colour[i + 1], zorder=2)
        out = (counts < lo) | (counts > hi)
        if out.any():
            ax.bar(np.arange(nBin)[out], counts[out], width=0.85, color="none",
                   edgecolor="k", lw=0.8, zorder=3)
        ax.set_xticks([])
        if i == 0:
            ax.set_title(short[q], loc="left", fontsize=7)
        if j == 0:
            ax.set_ylabel(f"{cellName[i + 1]}\ncount")
nOut = 0
for i in range(4):
    for q in order:
        rk = ranks[(ranks.cell == i + 1) & (ranks.quantity == q)]["rank"].dropna()
        counts = np.bincount((rk // (nRank // nBin)).astype(int), minlength=nBin)
        lo, hi = stats.binom.ppf([0.005, 0.995], len(rk), 1 / nBin)
        nOut += int(((counts < lo) | (counts > hi)).sum())
fig.suptitle(f"Rank histograms (10 bins of 10 ranks; grey: pointwise 99% binomial band): "
             f"{nOut} of {4 * len(order) * nBin} bins outside (≈{0.01 * 4 * len(order) * nBin:.0f} expected)",
             x=0.01, ha="left", fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(f"{d}/fig-rung6-full-ranks.png")
print(f"bins outside 99% band: {nOut}; failing tests: {nFail}")
