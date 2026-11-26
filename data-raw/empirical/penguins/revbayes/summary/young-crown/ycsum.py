#!/usr/bin/env python3
# Block summaries of the PR #99 young-crown test chains (RevBayes fd807d10).
#   python3 ycsum.py [block iterations, default 1000]
import os, re, sys, glob
BLK = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
YC = "/nobackup/pjjg18/evotempo-verify-99/youngcrown"
src = open("/nobackup/pjjg18/evotempo-penguins-beast2/beast-posterior-summary.py").read()
ns = {}
exec(src.split("burn = float")[0].replace("sys.argv[1]", "0"), ns)
exec("def parse" + src.split("def parse")[1].split("par_cols")[0], ns)
NODES, parse = ns["NODES"], ns["parse"]
EXT = NODES["crown_age"]
CF = ["Spheniscus_urbinai", "Spheniscus_megaramphus", "Spheniscus_muizoni",
      "Pygoscelis_grandis", "Madrynornis_mirandus"]
def ncf(nwk):
    par, lab, bl = parse(nwk)
    n = len(par); kids = [[] for _ in range(n)]
    for v in range(1, n): kids[par[v]].append(v)
    tip = {lab[v]: v for v in range(n) if not kids[v]}
    def anc(v):
        o = []
        while v >= 0: o.append(v); v = par[v]
        return o
    common = None
    for t in EXT:
        a = anc(tip[t]); common = a if common is None else [u for u in common if u in set(a)]
    cr = common[0]
    return sum(1 for f in CF if cr in set(anc(tip[f])))
print("chain block iters n Posterior Likelihood Prior crown min_crown cf_in_crown SA clock_dna")
for d in sorted(glob.glob(YC + "/*-*-*/")):
    lg = glob.glob(d + "output/*[0-9].log")
    if not lg: continue
    rows = [l.rstrip("\n").split("\t") for l in open(lg[0])]
    head, rows = rows[0], [r for r in rows[1:] if len(r) == len(rows[0])]
    ix = {h: i for i, h in enumerate(head)}
    trf = lg[0][:-4] + ".trees"
    cfin = {}
    for l in open(trf):
        p = l.rstrip("\n").split("\t")
        if len(p) > 1 and p[0].isdigit() and p[-1].endswith(";"):
            cfin[p[0]] = ncf(re.sub(r"\[[^\]]*\]", "", p[-1]))
    blocks = {}
    MIN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    rows = [r for r in rows if int(r[0]) >= MIN]
    for r in rows:
        b = int(r[0]) // BLK
        blocks.setdefault(b, []).append(r)
    for b in sorted(blocks):
        R = blocks[b]; m = lambda k: sum(float(r[ix[k]]) for r in R) / len(R)
        cf = [cfin[r[0]] for r in R if r[0] in cfin]
        fy = sum(1 for x in cf if x <= 2) / max(len(cf), 1)
        print("%s %d %d-%d %d %.1f %.1f %.1f %.2f %.2f %.2f %.2f %.4f fyoung=%.3f maxPost=%.1f" % (
            os.path.basename(d.rstrip("/")), b, int(R[0][0]), int(R[-1][0]), len(R), m("Posterior"), m("Likelihood"), m("Prior"),
            m("crown_age"), min(float(r[ix["crown_age"]]) for r in R), sum(cf) / max(len(cf), 1), m("num_samp_anc"), m("clock_dna"),
            fy, max(float(r[ix["Posterior"]]) for r in R)))
