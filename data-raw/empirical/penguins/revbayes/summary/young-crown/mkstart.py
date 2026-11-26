#!/usr/bin/env python3
# Convert one BEAST2 post2 sample (tree + parameters) into a RevBayes start
# tree (Newick, taxon names, extant tips at age 0, sampled ancestors as
# zero-length tips) and a Rev file of setValue() calls.
#   python3 mkstart.py <sample> <taxa.tsv> <loci.txt> <outstem>
import json, sys, re
sys.setrecursionlimit(10000)
smp, taxaf, locif, stem = sys.argv[1:5]
T = "/nobackup/pjjg18/evotempo-verify-99/trace/"
src = open("/nobackup/pjjg18/evotempo-penguins-beast2/beast-posterior-summary.py").read()
ns = {}
exec("import re\n" + "def parse" + src.split("def parse")[1].split("par_cols")[0], ns)
parse = ns["parse"]
J = json.load(open(T + "post2_rows.json"))
head, row, trans = J["head"], J["rows"][smp], J["trans"]
ix = {h: i for i, h in enumerate(head)}
val = lambda k: float(row[ix[k]])
tree = [l.split("\t", 1)[1].strip() for l in open(T + "post2_trees.tsv") if l.split("\t", 1)[0] == smp][0]
par, lab, bl = parse(tree)
n = len(par)
kids = [[] for _ in range(n)]
for v in range(1, n): kids[par[v]].append(v)
depth = [0.0] * n
for v in range(1, n): depth[v] = depth[par[v]] + bl[v]
rng = {}
for l in open(taxaf).read().splitlines()[1:]:
    t, a, b = l.strip().split("\t"); rng[t] = (float(a), float(b))
name = {v: trans.get(lab[v], lab[v]) for v in range(n) if not kids[v]}
dmax = max(depth[v] for v in name if rng[name[v]][1] == 0)
age = [dmax - d for d in depth]
nclamp = 0
for v, nm in name.items():
    a, b = rng[nm]
    if b == 0: age[v] = 0.0
    if bl[v] == 0.0: age[v] = age[par[v]]   # sampled ancestor
    if age[v] < a or age[v] > b:
        nclamp += 1
        print("WARN %s age %.12f outside [%g,%g]" % (nm, age[v], a, b))
def nwk(v):
    if kids[v]:
        s = "(" + ",".join(nwk(k) for k in kids[v]) + ")"
    else:
        s = name[v]
    if par[v] >= 0:
        s += ":%.12f" % (age[par[v]] - age[v])
    return s
open(stem + ".nwk", "w").write(nwk(0) + ";\n")
loci = [l.strip() for l in open(locif) if l.strip()]
L = ["diversification.setValue(%r)" % val("diversificationRate"),
     "turnover.setValue(%r)" % val("turnover"),
     "sampling_prop.setValue(%r)" % val("samplingProportion"),
     "origin_time.setValue(%r)" % val("origin"),
     "clock_morph.setValue(%r)" % val("ucldMean.morph"),
     "clock_dna.setValue(%r)" % val("ucldMean.dna"),
     "alpha_morph.setValue(%r)" % val("gammaShape.morph")]
L.append('if (clock_model == "shared") { ucld_sd.setValue(%r) } else { ucld_sd_morph.setValue(%r); ucld_sd_dna.setValue(%r) }'
         % (0.7, val("ucldStdev.morph"), max(val("ucldStdev.dna"), 0.05)))
mr = [val("mutationRate." + l) for l in loci]
L.append("locus_rates.setValue(simplex(v(%s)))" % ", ".join(repr(x) for x in mr))
for i, l in enumerate(loci, 1):
    L.append("alpha_dna[%d].setValue(%r)" % (i, val("gammaShape." + l)))
    for r in ["ac", "at", "cg", "ct", "gt"]:
        L.append("r_%s[%d].setValue(%r)" % (r, i, val("rate%s_%s" % (r.upper(), l))))
    L.append("pi[%d].setValue(simplex(v(%s)))" % (i, ", ".join(repr(val("freqParameter_%s%d" % (l, k))) for k in range(1, 5))))
open(stem + ".Rev", "w").write("\n".join(L) + "\n")
ext = [v for v in name if rng[name[v]][1] == 0]
print("sample", smp, "posterior", row[ix["posterior"]], "root age %.4f" % age[0], "origin", val("origin"),
      "n SA", sum(1 for v in name if bl[v] == 0.0), "clamp warnings", nclamp, "mutationRate sum %.6f" % sum(mr))
