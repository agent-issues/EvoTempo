# Data-free check of the RevBayes comparator's tree moves (#66): number of
# sampled ancestors without data, RevBayes a2a8e2db and v1.4.1 against BEAST2.
#
#   python3 data-raw/empirical/penguins/cet/revbayes/prior-check/prior-check-summary.py
#
# Reads the archived summaries in this folder (sa_vs_s_separate.csv: histogram
# of num_samp_anc; convergence_separate.csv: its ESS over both chains;
# derived_separate.csv: sampling proportion; beast2/beast2_prior_summary.txt:
# the BEAST2 histogram) and writes prior-check-summary.csv.  MCSE = sd /
# sqrt(ESS); the BEAST2 summary records no ESS, so its MCSE is not given.
import csv, math, os
here = os.path.dirname(os.path.abspath(__file__))

def stats(pairs):
    n = sum(c for _, c in pairs)
    m = sum(k * c for k, c in pairs) / n
    sd = math.sqrt(sum((k - m) ** 2 * c for k, c in pairs) / n)
    return m, sd, n

rows = []
for build in ("a2a8e2db", "v1.4.1"):
    d = os.path.join(here, build)
    with open(os.path.join(d, "sa_vs_s_separate.csv")) as f:
        pairs = [(int(r["num_samp_anc"]), int(r["n_samples"]))
                 for r in csv.DictReader(l for l in f if not l.startswith("#"))]
    with open(os.path.join(d, "convergence_separate.csv")) as f:
        ess = [float(r["ess_total"]) for r in csv.DictReader(f) if r["parameter"] == "num_samp_anc"][0]
    with open(os.path.join(d, "derived_separate.csv")) as f:
        sp = [float(r["mean"]) for r in csv.DictReader(f) if r["quantity"] == "sampling_prop"][0]
    m, sd, n = stats(pairs)
    rows.append(dict(program="RevBayes " + build, samples=n, meanSA=m, sdSA=sd, essSA=ess,
                     mcseSA=sd / math.sqrt(ess), meanSamplingProportion=sp))
with open(os.path.join(here, "beast2", "beast2_prior_summary.txt")) as f:
    txt = f.read().splitlines()
hist = [l for l in txt if l.startswith("SA count")][0].split("mean s: ", 1)[1]
pairs = [(int(x.split(":")[0]), int(x.split(":")[1].split(",")[0])) for x in hist.split("; ")]
sp = [float(l.split(",")[1]) for l in txt if l.startswith("samplingProportion")][0]
m, sd, n = stats(pairs)
rows.append(dict(program="BEAST2 2.4.8 + SA 1.1.7", samples=n, meanSA=m, sdSA=sd, essSA="",
                 mcseSA="", meanSamplingProportion=sp))
with open(os.path.join(here, "prior-check-summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r)
