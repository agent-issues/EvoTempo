#!/usr/bin/env python3
"""Check every sampled tree: fossil ages within their taxa.tsv ranges, and how
often a fossil is a sampled ancestor at the root (zero-length child of the root).
    python3 check-fossil-ages.py <taxa.tsv> <trees file> [...]"""
import re, sys
rng = {}
for i, l in enumerate(open(sys.argv[1])):
    f = l.rstrip("\n").split("\t")
    if i == 0: ix = {h: j for j, h in enumerate(f)}; continue
    rng[f[ix["taxon"]]] = (float(f[ix["min_age"]]), float(f[ix["max_age"]]))
def parse(nwk):
    nwk = re.sub(r"\[[^\]]*\]", "", nwk)
    par, lab, bl, stack, cur, i = [-1], [None], [0.0], [], 0, 0
    while i < len(nwk):
        ch = nwk[i]
        if ch == "(":
            par.append(cur); lab.append(None); bl.append(0.0); stack.append(cur); cur = len(par) - 1; i += 1
        elif ch == ",":
            cur = stack[-1]; par.append(cur); lab.append(None); bl.append(0.0); cur = len(par) - 1; i += 1
        elif ch == ")": cur = stack.pop(); i += 1
        elif ch == ":":
            j = i + 1
            while j < len(nwk) and nwk[j] not in ",();": j += 1
            bl[cur] = float(nwk[i + 1:j]); i = j
        elif ch == ";": break
        else:
            j = i
            while j < len(nwk) and nwk[j] not in ":,();": j += 1
            lab[cur] = nwk[i:j]; i = j
    return par, lab, bl
for tf in sys.argv[2:]:
    lines = open(tf).read().split("\n")[1:]
    n = bad = rootsa = 0; worst = 0.0
    for l in lines:
        if not l.strip().endswith(";"): continue
        par, lab, bl = parse(l.split("\t")[-1])
        m = len(par); depth = [0.0] * m
        for v in range(1, m): depth[v] = depth[par[v]] + bl[v]
        kids = [0] * m
        for v in range(1, m): kids[par[v]] += 1
        tips = {lab[v]: v for v in range(m) if kids[v] == 0}
        dmax = max(depth[v] for v in tips.values())
        n += 1
        for nm, v in tips.items():
            lo, hi = rng[nm]; age = dmax - depth[v]
            if hi > 0 and (age < lo - 1e-6 or age > hi + 1e-6):
                bad += 1; worst = max(worst, lo - age, age - hi)
            if hi > 0 and bl[v] == 0.0 and par[v] == 0: rootsa += 1
    print("%s: trees %d, fossil ages outside range %d (worst %.3g Myr), root sampled ancestors %d" % (tf.split("/")[-1], n, bad, worst, rootsa))
