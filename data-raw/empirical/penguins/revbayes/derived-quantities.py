#!/usr/bin/env python3
"""Derived FBD and clock quantities from the penguin RevBayes chains.

    python3 derived-quantities.py <run root> <clock> [burn-in fraction] [tree thinning]

Writes derived_<clock>.csv (mean, median, 95% HPD of each quantity, pooled over
chains after burn-in) and sa_vs_s_<clock>.csv (sampling proportion by number
of sampled ancestors) to ./summary/.

* psi = s mu / (1 - s) and lambda - mu - psi (Gavryushkina et al. 2017 report
  the latter as 0.019 [-0.061, 0.099]).
* Clock: the parametric coefficient of variation of a lognormal with log-SD
  sigma, sqrt(exp(sigma^2) - 1), and, from the _rates.trees files, the
  coefficient of variation of the sampled branch rates across branches
  (unweighted over all branches; unweighted over branches of positive length;
  and weighted by branch length), plus the length-weighted mean branch rate.
"""
import glob, math, os, re, sys

root, clock = sys.argv[1], sys.argv[2]
burnin = float(sys.argv[3]) if len(sys.argv) > 3 else 0.25
thin = int(sys.argv[4]) if len(sys.argv) > 4 else 5

def complete_lines(path):
    with open(path) as f:
        lines = f.read().split("\n")
    return lines[0], [l for l in lines[1:-1] if l]   # drop a possibly partial last line

def hpd(x, prob=0.95):
    x = sorted(x); n = len(x); k = int(math.floor(prob * n))
    if k >= n: return x[0], x[-1]
    w = [x[i + k] - x[i] for i in range(n - k)]
    i = w.index(min(w)); return x[i], x[i + k]

def summ(x):
    x = [v for v in x if v == v]
    s = sorted(x); n = len(s)
    med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    lo, hi = hpd(x)
    return sum(x) / n, med, lo, hi, n

def cv(v, w=None):
    if w is None:
        m = sum(v) / len(v); var = sum((a - m) ** 2 for a in v) / (len(v) - 1)
    else:
        W = sum(w); m = sum(a * b for a, b in zip(v, w)) / W
        var = sum(b * (a - m) ** 2 for a, b in zip(v, w)) / W
    return math.sqrt(var) / m, m

out = {}
def add(k, v): out.setdefault(k, []).append(v)
sa_s = []
logs = sorted(f for f in glob.glob(os.path.join(root, clock + "_run*", "output", "penguins_%s_run*.log" % clock)))
for lf in logs:
    head, rows = complete_lines(lf)
    names = head.split("\t"); ix = {n: i for i, n in enumerate(names)}
    rows = rows[int(len(rows) * burnin):]
    for r in rows:
        v = r.split("\t")
        g = lambda n: float(v[ix[n]])
        s, mu, lam = g("sampling_prop"), g("mu"), g("lambda")
        psi = s * mu / (1 - s)
        add("sampling_prop", s); add("lambda", lam); add("mu", mu); add("psi", psi)
        add("lambda_minus_mu", lam - mu); add("lambda_minus_mu_minus_psi", lam - mu - psi)
        add("num_samp_anc", g("num_samp_anc"))
        add("clock_morph", g("clock_morph")); add("clock_dna", g("clock_dna"))
        for sd in ("ucld_sd", "ucld_sd_morph", "ucld_sd_dna"):
            if sd in ix:
                add(sd, g(sd)); add("cv_param_" + sd, math.sqrt(math.exp(g(sd) ** 2) - 1))
        sa_s.append((int(g("num_samp_anc")), s))
        # how often a node's common ancestor is the crown root itself (the
        # published basal split, Aptenodytes + Pygoscelis vs the rest, makes
        # node_2 and node_4 younger than the crown)
        for nd in ("node_2", "node_4"):
            if nd in ix: add("p_%s_is_crown" % nd, 1.0 if g(nd) >= g("crown_age") - 1e-9 else 0.0)

    tf = lf[:-4] + "_rates.trees"
    if os.path.exists(tf):
        _, trows = complete_lines(tf)
        trows = trows[int(len(trows) * burnin)::thin]
        pat = re.compile(r":([0-9.eE+-]+)\[&dna_branch_rates=([0-9.eE+-]+),morph_branch_rates=([0-9.eE+-]+)\]")
        for t in trows:
            m = pat.findall(t.split("\t")[-1])
            bl = [float(a) for a, _, _ in m]; rd = [float(b) for _, b, _ in m]; rm = [float(c) for _, _, c in m]
            pos = [i for i, b in enumerate(bl) if b > 0]
            for lab, rr in (("morph", rm), ("dna", rd)):
                add("cv_branch_all_" + lab, cv(rr)[0])
                add("cv_branch_poslen_" + lab, cv([rr[i] for i in pos])[0])
                c_w, m_w = cv(rr, bl)
                add("cv_branch_lenweighted_" + lab, c_w); add("mean_rate_lenweighted_" + lab, m_w)
            add("n_branches_zero_length", len(bl) - len(pos))

os.makedirs("summary", exist_ok=True)
with open("summary/derived_%s.csv" % clock, "w") as f:
    f.write("quantity,mean,median,hpd95_lower,hpd95_upper,n\n")
    for k, x in out.items():
        f.write("%s,%.6g,%.6g,%.6g,%.6g,%d\n" % ((k,) + summ(x)))
with open("summary/sa_vs_s_%s.csv" % clock, "w") as f:
    f.write("num_samp_anc,n_samples,mean_sampling_prop\n")
    for k in sorted(set(a for a, _ in sa_s)):
        xs = [b for a, b in sa_s if a == k]
        f.write("%d,%d,%.4f\n" % (k, len(xs), sum(xs) / len(xs)))
    n = len(sa_s); ma = sum(a for a, _ in sa_s) / n; mb = sum(b for _, b in sa_s) / n
    cov = sum((a - ma) * (b - mb) for a, b in sa_s) / n
    va = sum((a - ma) ** 2 for a, _ in sa_s) / n; vb = sum((b - mb) ** 2 for _, b in sa_s) / n
    f.write("# correlation(num_samp_anc, sampling_prop) = %.3f over %d samples\n" % (cov / math.sqrt(va * vb), n))
print("chains:", len(logs))
