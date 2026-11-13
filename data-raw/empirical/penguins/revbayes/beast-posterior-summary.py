#!/usr/bin/env python3
"""Pooled posterior summaries of BEAST2 SA runs of the penguin XML.

    python3 beast-posterior-summary.py <burn-in fraction> <run dir> [<run dir> ...]

Each run dir holds Mkv_morph+dna.log and Mkv_morph+dna.trees.  Writes to
./summary/: beast2_parameters.csv (pooled mean, median, 95 pc HPD; per-chain
min ESS and summed ESS; includes the common-ancestor ages of the extant nodes
of Gavryushkina et al. 2017 Fig. 5, defined as in penguins_fbd.Rev) and
beast2_fossils.csv (P(sampled ancestor), P(in crown)) and beast2_basal_split.txt
(how often each node is the crown root itself).  A sampled ancestor is
a zero-length leaf branch; a fossil is in the crown if it descends from the
MRCA of the 19 extant species.
"""
import math, os, re, sys, json
NODES = json.loads('{"crown_age": ["Aptenodytes_forsteri", "Aptenodytes_patagonicus", "Eudyptes_chrysocome", "Eudyptes_chrysolophus", "Eudyptes_filholi", "Eudyptes_moseleyi", "Eudyptes_pachyrhynchus", "Eudyptes_robustus", "Eudyptes_schlegeli", "Eudyptes_sclateri", "Eudyptula_minor", "Megadyptes_antipodes", "Pygoscelis_adeliae", "Pygoscelis_antarctica", "Pygoscelis_papua", "Spheniscus_demersus", "Spheniscus_humboldti", "Spheniscus_magellanicus", "Spheniscus_mendiculus"], "node_2": ["Megadyptes_antipodes", "Eudyptes_chrysocome", "Eudyptes_chrysolophus", "Eudyptes_filholi", "Eudyptes_moseleyi", "Eudyptes_pachyrhynchus", "Eudyptes_robustus", "Eudyptes_schlegeli", "Eudyptes_sclateri", "Eudyptula_minor", "Spheniscus_demersus", "Spheniscus_humboldti", "Spheniscus_magellanicus", "Spheniscus_mendiculus"], "node_3": ["Eudyptula_minor", "Spheniscus_demersus", "Spheniscus_humboldti", "Spheniscus_magellanicus", "Spheniscus_mendiculus"], "node_4": ["Aptenodytes_forsteri", "Aptenodytes_patagonicus", "Pygoscelis_adeliae", "Pygoscelis_antarctica", "Pygoscelis_papua"], "node_5": ["Pygoscelis_adeliae", "Pygoscelis_antarctica", "Pygoscelis_papua"], "node_6": ["Megadyptes_antipodes", "Eudyptes_chrysocome", "Eudyptes_chrysolophus", "Eudyptes_filholi", "Eudyptes_moseleyi", "Eudyptes_pachyrhynchus", "Eudyptes_robustus", "Eudyptes_schlegeli", "Eudyptes_sclateri"], "node_7": ["Pygoscelis_antarctica", "Pygoscelis_papua"], "node_8": ["Eudyptes_chrysocome", "Eudyptes_chrysolophus", "Eudyptes_filholi", "Eudyptes_moseleyi", "Eudyptes_pachyrhynchus", "Eudyptes_robustus", "Eudyptes_schlegeli", "Eudyptes_sclateri"], "node_9": ["Eudyptes_chrysocome", "Eudyptes_filholi", "Eudyptes_moseleyi", "Eudyptes_chrysolophus", "Eudyptes_schlegeli"], "node_10": ["Eudyptes_pachyrhynchus", "Eudyptes_robustus", "Eudyptes_sclateri"], "node_11": ["Spheniscus_demersus", "Spheniscus_humboldti", "Spheniscus_magellanicus", "Spheniscus_mendiculus"], "node_12": ["Aptenodytes_forsteri", "Aptenodytes_patagonicus"], "node_13": ["Spheniscus_humboldti", "Spheniscus_mendiculus"], "node_14": ["Spheniscus_demersus", "Spheniscus_magellanicus"], "node_15": ["Eudyptes_chrysocome", "Eudyptes_filholi", "Eudyptes_moseleyi"], "node_16": ["Eudyptes_pachyrhynchus", "Eudyptes_robustus"], "node_17": ["Eudyptes_chrysolophus", "Eudyptes_schlegeli"], "node_18": ["Eudyptes_chrysocome", "Eudyptes_filholi"]}')
burn = float(sys.argv[1]); runs = sys.argv[2:]

def ess(x):
    n = len(x); m = sum(x) / n; v = sum((a - m) ** 2 for a in x) / n
    if v == 0: return float("nan")
    acf = lambda k: sum((x[i] - m) * (x[i + k] - m) for i in range(n - k)) / (n * v)
    s, k = 0.0, 1
    while k + 1 < n:
        p = acf(k) + acf(k + 1)
        if p <= 0: break
        s += p; k += 2
    return n / (1 + 2 * s)

def hpd(x, p=0.95):
    x = sorted(x); n = len(x); k = int(p * n); w = [x[i + k] - x[i] for i in range(n - k)]; i = w.index(min(w)); return x[i], x[i + k]

def parse(nwk):
    nwk = re.sub(r"\[[^\]]*\]", "", nwk)
    par, lab, bl, stack, cur, i = [-1], [None], [0.0], [], 0, 0
    while i < len(nwk):
        ch = nwk[i]
        if ch == "(":
            par.append(cur); lab.append(None); bl.append(0.0); stack.append(cur); cur = len(par) - 1; i += 1
        elif ch == ",":
            cur = stack[-1]; par.append(cur); lab.append(None); bl.append(0.0); cur = len(par) - 1; i += 1
        elif ch == ")":
            cur = stack.pop(); i += 1
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

par_cols = ["samplingProportion", "turnover", "diversificationRate", "SACount", "origin", "TreeHeight",
            "ucldMean.morph", "ucldMean.dna", "ucldStdev.morph", "ucldStdev.dna",
            "rate.morph.mean", "rate.morph.coefficientOfVariation", "rate.dna.mean", "rate.dna.coefficientOfVariation"]
SPH = ["Spheniscus_demersus", "Spheniscus_humboldti", "Spheniscus_magellanicus", "Spheniscus_mendiculus"]
chains = {}
fossil_sa, fossil_crown, ntree, muiz = {}, {}, 0, 0
for rd in runs:
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(rd, "Mkv_morph+dna.log")) if l.strip() and not l.startswith("#")]
    head, rows = rows[0], rows[1:]
    rows = [r for r in rows if len(r) == len(head)]
    rows = rows[int(len(rows) * burn):]
    ix = {h: i for i, h in enumerate(head)}
    ch = {}
    for n in par_cols:
        if n in ix: ch[n] = [float(r[ix[n]]) for r in rows]
    lam = [d / (1 - v) for d, v in zip(ch["diversificationRate"], ch["turnover"])]
    mu = [v * l for v, l in zip(ch["turnover"], lam)]
    psi = [s * m / (1 - s) for s, m in zip(ch["samplingProportion"], mu)]
    ch["lambda_minus_mu_minus_psi"] = [l - m - p for l, m, p in zip(lam, mu, psi)]
    txt = open(os.path.join(rd, "Mkv_morph+dna.trees")).read()
    trans = dict(re.findall(r"(\d+)\s+([^,\s;]+)", re.search(r"Translate(.*?);", txt, re.S | re.I).group(1)))
    trees = re.findall(r"^tree\s+\S+\s*=\s*(?:\[[^\]]*\]\s*)?(.*;)\s*$", txt, re.M)
    trees = trees[int(len(trees) * burn):]
    ages = {k: [] for k in NODES}
    for t in trees:
        par, lab, bl = parse(t)
        n = len(par)
        depth = [0.0] * n
        for v in range(1, n): depth[v] = depth[par[v]] + bl[v]
        kids = [[] for _ in range(n)]
        for v in range(1, n): kids[par[v]].append(v)
        name = {v: trans.get(lab[v], lab[v]) for v in range(n) if not kids[v]}
        tip = {nm: v for v, nm in name.items()}
        dmax = max(depth[tip[e]] for e in NODES["crown_age"])
        def anc(v):
            out = []
            while v >= 0: out.append(v); v = par[v]
            return out
        def mrca(taxa):
            common = None
            for tx in taxa:
                a = anc(tip[tx])
                common = a if common is None else [u for u in common if u in set(a)]
            return common[0]
        def subtree(u):
            out, st = set(), [u]
            while st:
                w = st.pop(); out.add(w); st.extend(kids[w])
            return out
        for k, taxa in NODES.items():
            ages[k].append(dmax - depth[mrca(taxa)])
        below = subtree(mrca(NODES["crown_age"]))
        for v, nm in name.items():
            if nm in NODES["crown_age"]: continue
            fossil_sa[nm] = fossil_sa.get(nm, 0) + (bl[v] == 0.0)
            fossil_crown[nm] = fossil_crown.get(nm, 0) + (v in below)
        mv = tip["Spheniscus_muizoni"]
        if bl[mv] == 0.0:
            sib = [u for u in kids[par[mv]] if u != mv][0]
            sub = subtree(sib)
            if all(tip[e] in sub for e in SPH): muiz += 1
        ntree += 1
    ch.update(ages)
    chains[rd] = ch
os.makedirs("summary", exist_ok=True)
keys = list(next(iter(chains.values())).keys())
with open("summary/beast2_parameters.csv", "w") as f:
    f.write("parameter,mean,median,hpd95_lower,hpd95_upper,ess_min_chain,ess_total,n_chains,samples_total\n")
    for k in keys:
        x = [v for ch in chains.values() for v in ch[k]]
        s = sorted(x); med = s[len(s) // 2]; lo, hi = hpd(x)
        e = [ess(ch[k]) for ch in chains.values()]
        f.write("%s,%.6g,%.6g,%.6g,%.6g,%.0f,%.0f,%d,%d\n" % (k, sum(x) / len(x), med, lo, hi, min(e), sum(e), len(chains), len(x)))
with open("summary/beast2_fossils.csv", "w") as f:
    f.write("taxon,p_sampled_ancestor,p_in_crown\n")
    for nm in sorted(fossil_sa, key=lambda z: -fossil_sa[z]):
        f.write("%s,%.4f,%.4f\n" % (nm, fossil_sa[nm] / ntree, fossil_crown[nm] / ntree))
with open("summary/beast2_basal_split.txt", "w") as f:
    for k in NODES:
        if k == "crown_age": continue
        x = [a for ch in chains.values() for a in ch[k]]; cr = [a for ch in chains.values() for a in ch["crown_age"]]
        f.write("P(%s common-ancestor age = crown age) = %.3f\n" % (k, sum(1 for a, b in zip(x, cr) if a >= b - 1e-9) / len(x)))
with open("summary/beast2_headline.txt", "w") as f:
    f.write("chains: %s\ntrees after burn-in: %d\nP(S. muizoni is a sampled ancestor of the four extant Spheniscus): %.3f\n" % (" ".join(runs), ntree, muiz / ntree))
print(open("summary/beast2_headline.txt").read())
