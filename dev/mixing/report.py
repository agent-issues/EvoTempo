"""Tables and figure for dev/mixing/README.md (issue #19).

Usage: python dev/mixing/report.py <resultsDir>
Reads chains.csv, rhat.csv, probes.csv and oracles.csv written by
dev/mixing/summary.R; writes tables.md into <resultsDir>.
"""
import sys
import numpy as np
import pandas as pd

d = sys.argv[1] if len(sys.argv) > 1 else "dev/mixing/results"
ch = pd.read_csv(f"{d}/chains.csv")
probes = pd.read_csv(f"{d}/probes.csv")
oracles = pd.read_csv(f"{d}/oracles.csv")
rhat = pd.read_csv(f"{d}/rhat.csv")
treeName = {1: "small", 2: "large"}
qs = ["theta", "s", "rate", "lambda", "rootBin", "clock", "charLL"]
ch["hours"] = ch["secs"] / 3600
ch["ess.clockMed"] = ch["essClockMed"]
for q in qs + ["clockMed"]:
    ch[f"eph.{q}"] = ch[f"ess.{q}"] / ch["hours"]
    ch[f"epk.{q}"] = ch[f"ess.{q}"] / ch["nIter"] * 1000

# Projected engine cost of anchored-joint proposals: in the prototype they
# run R loops over segments; an engine version costs about one `joint`
# proposal (same ODE solve, messages, whole-tree draw and character
# likelihood).  Replace their measured seconds by tries x the mean seconds
# per joint proposal of the base runs on the same tree and character count.
jointCost = (ch[ch.config == "base"].assign(
    c=lambda x: x["share.joint"] * x["secPerIter"] / 6)
    .groupby(["tree", "nChar"])["c"].mean())
ajCols = [c for c in ch.columns if c.startswith("share.anchored")]
proj = []
for _, r in ch.iterrows():
    shareAj = np.nansum([r[c] for c in ajCols])
    if shareAj > 0 and (r.tree, r.nChar) in jointCost.index:
        # anchored proposals per iteration: one per parameter in the move
        nAj = sum(len(c.split("[")[1].split("]")[0].split(",")) for c in ajCols
                  if np.isfinite(r[c]) and r[c] > 0)
        perIt = r.secPerIter * (1 - shareAj) + nAj * jointCost[(r.tree, r.nChar)]
        proj.append(r.secs * perIt / r.secPerIter)
    else:
        proj.append(r.secs)
ch["secsProj"] = proj
for q in qs + ["clockMed"]:
    ch[f"ephp.{q}"] = ch[f"ess.{q}"] / (ch["secsProj"] / 3600)

out = []
cell = ch.groupby(["tree", "nChar", "config"])
mean = cell.mean(numeric_only=True)
nrep = cell.size()

def Fmt(x, digits=0):
    if not np.isfinite(x):
        return "-"
    return f"{x:.{digits}f}" if abs(x) < 1e5 else f"{x:.2g}"

# Table 1: probe summary (medians over states)
pt = probes.groupby(["tree", "nChar", "proposal"]).median(numeric_only=True)
out.append("### Fresh path proposals at fixed parameters (median over 10 states)\n")
out.append("| tree | chars | proposal | sd charLL | sd log weight | acceptance | sd log clock |")
out.append("|---|---|---|---|---|---|---|")
order = ["tree", "anchoredInf", "anchored4", "anchored2", "anchored1"]
for (t, n), _ in probes.groupby(["tree", "nChar"]):
    for p in order:
        if (t, n, p) in pt.index:
            r = pt.loc[(t, n, p)]
            out.append(f"| {treeName[t]} | {n} | {p} | {r.sdChar:.2f} | {r.sdLw:.3f} | "
                       f"{r.accept:.2g} | {r.sdLogClock:.3f} |")
out.append("")
out.append("### Conditional spread of log parameters given the paths (median over states)\n")
out.append("| tree | chars | nJumps | theta | s | rate | lambda | mu | psi |")
out.append("|---|---|---|---|---|---|---|---|---|")
for (t, n), g in probes[probes.proposal == "tree"].groupby(["tree", "nChar"]):
    m = g.median(numeric_only=True)
    out.append(f"| {treeName[t]} | {n} | {m.nJumps:.0f} | {m['condSdLog.theta']:.3f} | "
               f"{m['condSdLog.s']:.3f} | {m['condSdLog.rate']:.3f} | "
               f"{m['condSdLog.lambda']:.3f} | {m['condSdLog.mu']:.3f} | "
               f"{m['condSdLog.psi']:.3f} |")
out.append("")

# Table 2: time shares and acceptance of the current schedule
out.append("### Current schedule (`base`): acceptance and share of CPU time\n")
out.append("| tree | chars | s / iteration | tree acc | tree share | joint acc | joint share | "
           "local acc (node/tip/branch) | local share | sweep share |")
out.append("|---|---|---|---|---|---|---|---|---|---|")
for (t, n, c), r in mean.iterrows():
    if c != "base":
        continue
    out.append(f"| {treeName[t]} | {n} | {r.secPerIter:.3f} | {r['acc.tree']:.2g} | "
               f"{r['share.tree']:.2f} | {r['acc.joint']:.2g} | {r['share.joint']:.2f} | "
               f"{r['acc.node']:.2f}/{r['acc.tip']:.2f}/{r['acc.branch']:.2f} | "
               f"{r['share.node'] + r['share.tip'] + r['share.branch']:.2f} | "
               f"{r['share.parameters']:.2f} |")
out.append("")

# Table 3: posterior spread, a convergence check across samplers
out.append("### Posterior sd of log parameters by sampler (mean of replicates)\n")
cfgs = ["base", "shift", "L25", "L25S", "L25SA", "L25SN"]
out.append("| tree | chars | quantity | " + " | ".join(cfgs) + " |")
out.append("|---|---|---|" + "---|" * len(cfgs))
for (t, n), _ in ch[ch.nChar > 0].groupby(["tree", "nChar"]):
    for q in ["theta", "s", "rate"]:
        vals = [Fmt(mean.loc[(t, n, c)][f"postSdLog.{q}"], 3) if (t, n, c) in mean.index
                else "-" for c in cfgs]
        out.append(f"| {treeName[t]} | {n} | {q} | " + " | ".join(vals) + " |")
out.append("")

# Table 4: gains
show = ["base", "shift", "L25", "L25S", "L25SA", "L25SN"]
for label, pre in [("ESS per hour of CPU, as measured", "eph"),
                   ("ESS per 1000 iterations", "epk")]:
    out.append(f"### {label} (mean of replicates; `n` = replicates)\n")
    out.append("| tree | chars | config | n | " + " | ".join(qs + ["branch clock (median)"]) +
               " |")
    out.append("|---|---|---|---|" + "---|" * (len(qs) + 1))
    for (t, n), _ in ch[ch.nChar > 0].groupby(["tree", "nChar"]):
        for c in show:
            if (t, n, c) not in mean.index:
                continue
            r = mean.loc[(t, n, c)]
            vals = [Fmt(r[f"{pre}.{q}"], 1 if pre == "epk" else 0) for q in qs + ["clockMed"]]
            out.append(f"| {treeName[t]} | {n} | {c} | {nrep[(t, n, c)]} | " +
                       " | ".join(vals) + " |")
    out.append("")
out.append("### Gain of `L25S` over `base`, ESS per CPU hour (ratio of means)\n")
out.append("| tree | chars | " + " | ".join(qs + ["branch clock"]) + " |")
out.append("|---|---|" + "---|" * (len(qs) + 1))
for (t, n), _ in ch[ch.nChar > 0].groupby(["tree", "nChar"]):
    if (t, n, "L25S") in mean.index and (t, n, "base") in mean.index:
        a = mean.loc[(t, n, "L25S")]
        b = mean.loc[(t, n, "base")]
        out.append(f"| {treeName[t]} | {n} | " + " | ".join(
            f"{a[f'eph.{q}'] / b[f'eph.{q}']:.1f}" for q in qs + ["clockMed"]) + " |")
out.append("")
out.append("### Projected engine cost: anchored proposals priced as one `joint` proposal\n")
out.append("| tree | chars | config | " + " | ".join(qs) + " |")
out.append("|---|---|---|" + "---|" * len(qs))
for (t, n), _ in ch[ch.nChar > 0].groupby(["tree", "nChar"]):
    for c in ["L25S", "L25SA", "aj1", "ajInf"]:
        if (t, n, c) in mean.index:
            r = mean.loc[(t, n, c)]
            out.append(f"| {treeName[t]} | {n} | {c} | " +
                       " | ".join(Fmt(r[f"ephp.{q}"]) for q in qs) + " |")
out.append("")

# Table 5: R-hat
out.append("### Split R-hat over replicate chains\n")
rq = ["theta", "s", "rate", "lambda", "rootBin", "clock", "charLL"]
out.append("| cell | " + " | ".join(rq) + " |")
out.append("|---|" + "---|" * len(rq))
for _, r in rhat.sort_values("cell").iterrows():
    if any(k in r.cell for k in ["-base", "-L25S", "-shift", "-L25SA", "-L25SN"]):
        out.append(f"| {r.cell} | " + " | ".join(Fmt(r[f'rhat.{q}'], 2) for q in rq) + " |")
out.append("")

# Table 6: oracles
out.append("### Invariance oracles (rung 4b tree, one character)\n")
out.append("| variant | parameter | reference | z at q05 / q25 / q50 / q75 / q95 | max abs z |")
out.append("|---|---|---|---|---|")
for (v, p), g in oracles.groupby(["variant", "parameter"]):
    z = g.sort_values("p").z.values
    out.append(f"| {v} | {p} | {g.reference.iloc[0]} | " +
               " / ".join(f"{x:+.2f}" for x in z) + f" | {np.max(np.abs(z)):.2f} |")
out.append("")
open(f"{d}/tables.md", "w").write("\n".join(out))

print("wrote", f"{d}/tables.md")
