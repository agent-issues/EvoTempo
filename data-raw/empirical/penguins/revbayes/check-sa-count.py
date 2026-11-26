#!/usr/bin/env python3
"""Compare the logged num_samp_anc with the number of zero-length fossil tip
branches in the tree sampled at the same iteration (root children included).
    python3 check-sa-count.py <taxa.tsv> <run output dir> [...]"""
import re, sys
fos = set()
for i, l in enumerate(open(sys.argv[1])):
    f = l.rstrip().split("\t")
    if i and float(f[2]) > 0: fos.add(f[0])
for od in sys.argv[2:]:
    import glob
    lf = [x for x in glob.glob(od + "/penguins_*_run*.log")][0]; tf = lf[:-4] + ".trees"
    L = [l.rstrip("\n").split("\t") for l in open(lf) if l.strip()]
    ix = {h: i for i, h in enumerate(L[0])}
    logged = {r[0]: int(float(r[ix["num_samp_anc"]])) for r in L[1:] if len(r) == len(L[0])}
    n = same = minus1 = other = 0; rootsa = 0
    for l in open(tf).read().split("\n")[1:]:
        if not l.strip().endswith(";"): continue
        it = l.split("\t")[0]; nwk = re.sub(r"\[[^\]]*\]", "", l.split("\t")[-1])
        zl = sum(1 for nm, bl in re.findall(r"([A-Za-z_]+):([0-9.]+)", nwk) if nm in fos and float(bl) == 0.0)
        if it not in logged: continue
        n += 1; d = zl - logged[it]
        same += d == 0; minus1 += d == 1; other += d not in (0, 1)
    print("%s: %d trees; zero-length fossil tips == logged %d, == logged+1 %d, other %d" % (od.split("/")[-2], n, same, minus1, other))
