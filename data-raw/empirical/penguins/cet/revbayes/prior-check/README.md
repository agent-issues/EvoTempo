# Data-free check of the comparator's tree moves (#66)

The morphology-only comparator (`../penguins_morph_fbd.Rev`) keeps
`mvCollapseExpandFossilBranch`, which in RevBayes v1.4.1 is biased when fossil
ages are uncertain (fixed upstream in c389174b). Every penguin fossil has an age
range, so the comparator runs on the project's development build a2a8e2db
(`../run-hamilton.sh`). This folder records the data-free run behind that choice.

## What was run (Hamilton, 2026-11-02)

- Model: the total-evidence baseline script `penguins_fbd.Rev` (separate morphological
  and DNA clocks) with `under_prior = true`, which calls `ignoreAllData()`. Two chains,
  3600 samples each after burn-in. The baseline is not yet on GitHub (it will live in
  `../../../revbayes/`). The copy committed here comes from the baseline's
  work-in-progress archive, saved on 2026-11-02 after this run. Its move block (lines
  104-114) is the one the verifier read in the run directory; it was not compared byte
  for byte with the copy that ran (see below).
- Tree moves: exactly the comparator's six, with the same weights: `mvFNPR` 15,
  `mvNarrow` 15, `mvCollapseExpandFossilBranch` 10, `mvNodeTimeSlideUniform` 40,
  `mvRootTimeSlideUniform` 5, and `mvFossilTipTimeUniform` 1 per ranged tip. These
  were read in the run directory by the #43 verifier on 2026-11-02. The run has
  sampled ancestors, an origin and ranged tip ages, and differs from the comparator
  only in the FBD parameterisation (d, turnover, sampling proportion), which these
  moves do not touch.
- RevBayes a2a8e2db (`/nobackup/pjjg18/revbayes-dev-a2a8e2db`, `rb --version`
  1.4.2-preview), and the same on v1.4.1 (`/nobackup/pjjg18/revbayes-1.4.1`).
- Reference: BEAST2 2.4.8 with SA 1.1.7, the source analysis's XML with the data
  removed (`/nobackup/pjjg18/evotempo-penguins-summaries/beast2-prior/penguins_prior.xml`,
  not committed: 0.6 MB), 9001 samples.

The summaries were made by the baseline's `summarise-hamilton.sh` and archived to
`/nobackup/pjjg18/evotempo-penguins-summaries/prior-{dev-a2a8e2db,v1.4.1}/`. The
run directory (`/nobackup/pjjg18/evotempo-penguins-rbdev-prior`, with its driver,
logs and trees) was deleted on 2026-11-03 before this commit, in another session's
scratch cleanup. So the exact driver is not committed, and the numbers below come
from the archived summaries only.

## Files

| File | Content |
|---|---|
| `a2a8e2db/`, `v1.4.1/` | The archived summaries: `headline_separate.txt`, `parameters_separate.csv`, `convergence_separate.csv` (ESS, PSRF), `derived_separate.csv`, `sa_vs_s_separate.csv` (histogram of the sampled-ancestor count), `fossils_separate.csv` |
| `beast2/beast2_prior_summary.txt` | BEAST2 data-free summary, with the sampled-ancestor histogram |
| `penguins_fbd.Rev` | The model script (copy from the baseline's archive; see above) |
| `prior-check-summary.py`, `prior-check-summary.csv` | Mean, sd, ESS and MCSE of the sampled-ancestor count, and the mean sampling proportion, per program |

## Result

| Program | Samples | Sampled ancestors: mean (sd) | ESS | MCSE | Sampling proportion |
|---|---|---|---|---|---|
| RevBayes a2a8e2db | 7202 | 17.18 (6.83) | 873 | 0.23 | 0.266 |
| RevBayes v1.4.1 | 7202 | 31.97 (3.00) | 2013 | 0.07 | 0.735 |
| BEAST2 | 9001 | 16.57 (6.85) | not recorded | - | 0.254 |

a2a8e2db agrees with BEAST2: the difference of 0.61 is 2.6 a2a8e2db MCSE, or
about 1.9 combined MCSE if BEAST2's MCSE is similar, and the histograms have the
same spread. v1.4.1 is far off (32 against 16.6). Not checked: node ages against
an exact oracle, and the v1.4.1 run's move set (only the a2a8e2db run's was read).
