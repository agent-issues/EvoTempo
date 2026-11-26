#!/usr/bin/env python3
"""Summaries of a BEAST2 SA run: sampling proportion, SA count and per-fossil
P(sampled ancestor) (zero-length leaf branches) after burn-in.
    python3 beast-sa-summary.py <log> <trees> [burn-in fraction]"""
import math, re, sys
logf, treef = sys.argv[1], sys.argv[2]; burn = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
rows = [l.rstrip("\n").split("\t") for l in open(logf) if l.strip() and not l.startswith("#")]
head, rows = rows[0], rows[1:]
rows = rows[int(len(rows) * burn):]
ix = {h: i for i, h in enumerate(head)}
def hpd(x, p=0.95):
    x = sorted(x); n = len(x); k = int(p * n); w = [x[i + k] - x[i] for i in range(n - k)]; i = w.index(min(w)); return x[i], x[i + k]
def col(n): return [float(r[ix[n]]) for r in rows]
out = {}
for n in ["samplingProportion", "turnover", "diversificationRate", "SACount", "origin", "ucldMean.c:penguins_morph", "ucldMean.c:penguins_dna", "ucldStdev.c:penguins_morph", "ucldStdev.c:penguins_dna"]:
    if n in ix: out[n] = col(n)
if all(k in out for k in ("samplingProportion", "turnover", "diversificationRate")):
    lam = [d / (1 - v) for d, v in zip(out["diversificationRate"], out["turnover"])]
    mu = [v * l for v, l in zip(out["turnover"], lam)]
    psi = [s * m / (1 - s) for s, m in zip(out["samplingProportion"], mu)]
    out["lambda_minus_mu_minus_psi"] = [l - m - p for l, m, p in zip(lam, mu, psi)]
for n in ["rate.c:penguins_morph.mean", "rate.c:penguins_morph.coefficientOfVariation", "rate.c:penguins_dna.mean", "rate.c:penguins_dna.coefficientOfVariation"]:
    if n in ix: out[n] = col(n)
print("n samples", len(rows))
for n, x in out.items():
    lo, hi = hpd(x); print("%s,%.5g,%.5g,%.5g" % (n, sum(x) / len(x), lo, hi))
if "SACount" in out and "samplingProportion" in out:
    by = {}
    for a, s in zip(out["SACount"], out["samplingProportion"]): by.setdefault(int(a), []).append(s)
    print("SA count: n, mean s:", "; ".join("%d:%d,%.3f" % (k, len(v), sum(v) / len(v)) for k, v in sorted(by.items())))
# trees
txt = open(treef).read()
tr = re.search(r"Translate(.*?);", txt, re.S | re.I).group(1)
trans = dict(re.findall(r"(\d+)\s+([^,\s]+)", tr))
trees = re.findall(r"^tree\s+\S+\s*=\s*(?:\[[^\]]*\]\s*)?(.*);\s*$", txt, re.M)
trees = trees[int(len(trees) * burn):]
cnt = {}
for t in trees:
    t = re.sub(r"\[[^\]]*\]", "", t)
    for lab, bl in re.findall(r"[(,](\d+):([0-9.eE+-]+)", t):
        if float(bl) == 0.0: cnt[trans[lab]] = cnt.get(trans[lab], 0) + 1
print("trees", len(trees))
for k in sorted(trans.values()):
    if k in cnt or True:
        p = cnt.get(k, 0) / len(trees)
        if p > 0: print("PSA,%s,%.3f" % (k, p))
