#!/usr/bin/env python3
# Per-sample traces of the six BEAST2 posterior chains (post1..post6):
# log columns + crown age + crown-fossil placement, computed from the trees.
# Writes trace99.csv, blocks99.csv and post2 young-crown candidate trees.
import os, re, sys, json
B = "/nobackup/pjjg18/evotempo-penguins-beast2"
sys.argv = ["x", "0"]
src = open(os.path.join(B, "beast-posterior-summary.py")).read()
# reuse NODES and parse() only
ns = {}
exec(src.split("burn = float")[0], ns)
exec("def parse" + src.split("def parse")[1].split("par_cols")[0], ns)
NODES, parse = ns["NODES"], ns["parse"]
EXT = NODES["crown_age"]
CF = ["Spheniscus_urbinai", "Spheniscus_megaramphus", "Spheniscus_muizoni",
      "Pygoscelis_grandis", "Madrynornis_mirandus"]
cols = ["posterior", "likelihood", "prior", "birthDeath", "SACount", "origin", "TreeHeight",
        "ucldMean.dna", "ucldMean.morph", "ucldStdev.dna", "ucldStdev.morph", "samplingProportion",
        "diversificationRate", "turnover"]
out = open("trace99.csv", "w")
out.write("chain,sample," + ",".join(cols) + ",morphLik,dnaLik,crown_age,node2_is_crown,n_cf_in_crown\n")
cand = {}
for c in range(1, 7):
    rd = os.path.join(B, "post%d" % c)
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(rd, "Mkv_morph+dna.log")) if l.strip() and not l.startswith("#")]
    head, rows = rows[0], rows[1:]
    rows = [r for r in rows if len(r) == len(head)]
    ix = {h: i for i, h in enumerate(head)}
    mcols = [i for h, i in ix.items() if h.startswith("morphTreeLikelihood")]
    dcols = [i for h, i in ix.items() if h.startswith("treeLikelihood.")]
    txt = open(os.path.join(rd, "Mkv_morph+dna.trees")).read()
    trans = dict(re.findall(r"(\d+)\s+([^,\s;]+)", re.search(r"Translate(.*?);", txt, re.S | re.I).group(1)))
    trees = re.findall(r"^tree\s+STATE_(\d+)\s*=\s*(?:\[[^\]]*\]\s*)?(.*;)\s*$", txt, re.M)
    tdict = {int(s): t for s, t in trees}
    for r in rows:
        smp = int(r[ix["Sample"]])
        t = tdict.get(smp)
        if t is None:
            continue
        par, lab, bl = parse(t)
        n = len(par)
        depth = [0.0] * n
        for v in range(1, n): depth[v] = depth[par[v]] + bl[v]
        kids = [[] for _ in range(n)]
        for v in range(1, n): kids[par[v]].append(v)
        name = {v: trans.get(lab[v], lab[v]) for v in range(n) if not kids[v]}
        tip = {nm: v for v, nm in name.items()}
        dmax = max(depth[tip[e]] for e in EXT)
        def anc(v):
            o = []
            while v >= 0: o.append(v); v = par[v]
            return o
        def mrca(taxa):
            common = None
            for tx in taxa:
                a = anc(tip[tx]); common = a if common is None else [u for u in common if u in set(a)]
            return common[0]
        cr = mrca(EXT); n2 = mrca(NODES["node_2"])
        crown = dmax - depth[cr]
        ca = set(anc(tip[CF[0]]))
        ncf = sum(1 for f in CF if cr in set(anc(tip[f])))
        ml = sum(float(r[i]) for i in mcols); dl = sum(float(r[i]) for i in dcols)
        vals = [r[ix[k]] for k in cols]
        out.write("%d,%d,%s,%.4f,%.4f,%.5f,%d,%d\n" % (c, smp, ",".join(vals), ml, dl, crown, int(n2 == cr), ncf))
        if c == 2:
            cand.setdefault("rows", []).append((smp, float(r[ix["posterior"]]), crown, ncf))
            if smp % 5000000 == 0 or True:
                pass
    if c == 2:
        # store trees + full parameter rows for later selection
        json.dump({"head": head, "rows": {r[ix["Sample"]]: r for r in rows},
                   "trans": trans}, open("post2_rows.json", "w"))
        with open("post2_trees.tsv", "w") as f:
            for s, t in trees: f.write("%s\t%s\n" % (s, t))
out.close()
print("done")
