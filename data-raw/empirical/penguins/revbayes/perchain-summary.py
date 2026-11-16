#!/usr/bin/env python3
"""Per-chain posterior summaries of the penguin RevBayes FBD chains.

    python3 perchain-summary.py <taxa.tsv> <burn-in fraction> <run dir> [<run dir> ...]

Each run dir holds output/penguins_<clock>_run<k>.log and .trees.  For every
chain, after discarding the first <burn-in fraction> of its samples (all
samples kept, the last one included), writes to stdout a CSV with one row per
chain and quantity:
  num_samp_anc, sampling_prop (means, from the log),
  p_node_2_is_crown (node_2 common-ancestor age equals the crown age, as in
  derived-quantities.py),
  p_sa:<fossil>, p_crown:<fossil> (zero-length fossil branch as printed;
  descent from the MRCA of the extant taxa).
Fossils are the taxa of <taxa.tsv> with max_age > 0.
"""
import csv, re, sys

def parse(nwk):
    nwk = re.sub(r"\[[^\]]*\]", "", nwk).strip().rstrip(";")
    par, lab, bl = [-1], [None], [0.0]
    cur, i, n = 0, 0, len(nwk)
    while i < n:
        ch = nwk[i]
        if ch == "(":
            par.append(cur); lab.append(None); bl.append(0.0); cur = len(par) - 1; i += 1
        elif ch == ",":
            cur = par[cur]
            par.append(cur); lab.append(None); bl.append(0.0); cur = len(par) - 1; i += 1
        elif ch == ")":
            cur = par[cur]; i += 1
        elif ch == ":":
            j = i + 1
            while j < n and nwk[j] not in ",():;": j += 1
            bl[cur] = float(nwk[i + 1:j]); i = j
        else:
            j = i
            while j < n and nwk[j] not in ",():;": j += 1
            lab[cur] = nwk[i:j]; i = j
    return par, lab, bl

def main():
    taxa, burn, dirs = sys.argv[1], float(sys.argv[2]), sys.argv[3:]
    rows = list(csv.DictReader(open(taxa), delimiter="\t"))
    extant = {r["taxon"] for r in rows if float(r["max_age"]) == 0}
    fossils = sorted(r["taxon"] for r in rows if float(r["max_age"]) > 0)
    w = csv.writer(sys.stdout); w.writerow(["chain", "quantity", "value", "n"])
    for d in dirs:
        name = d.rstrip("/").split("/")[-1]
        stem = "%s/output/penguins_%s" % (d, name)
        lg = list(csv.DictReader(open(stem + ".log"), delimiter="\t"))
        lg = lg[int(burn * len(lg)):]
        for q in ["num_samp_anc", "sampling_prop"]:
            w.writerow([name, q, "%.6g" % (sum(float(r[q]) for r in lg) / len(lg)), len(lg)])
        p2 = sum(abs(float(r["node_2"]) - float(r["crown_age"])) < 1e-9 for r in lg) / len(lg)
        w.writerow([name, "p_node_2_is_crown", "%.4f" % p2, len(lg)])
        tr = [l.rstrip("\n").split("\t") for l in open(stem + ".trees")][1:]
        tr = tr[int(burn * len(tr)):]
        sa = dict.fromkeys(fossils, 0); cr = dict.fromkeys(fossils, 0)
        for t in tr:
            par, lab, bl = parse(t[-1])
            kids = [[] for _ in par]
            for v, p in enumerate(par):
                if p >= 0: kids[p].append(v)
            tip = {l: v for v, l in enumerate(lab) if l is not None and not kids[v]}
            anc = lambda v: [v] + (anc(par[v]) if par[v] >= 0 else [])
            common = None
            for e in extant:
                a = anc(tip[e]); common = a if common is None else [u for u in common if u in set(a)]
            m = common[0]
            for f in fossils:
                v = tip[f]; sa[f] += bl[v] == 0.0; cr[f] += m in anc(v)
        for f in fossils:
            w.writerow([name, "p_sa:" + f, "%.4f" % (sa[f] / len(tr)), len(tr)])
            w.writerow([name, "p_crown:" + f, "%.4f" % (cr[f] / len(tr)), len(tr)])

if __name__ == "__main__":
    main()
