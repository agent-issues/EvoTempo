# Rung 5 with a free start age (#18)

Validation of the origin start (`start = "origin"`: stem path, origin bin, single-lineage
survival) and of the start-age moves (`rootAge`, `originAge`, and `nodeAge` on the root under an
origin) added to `TempoMCMC()` in this PR. Everything here was run on Hamilton with the C++
engine, at branch `origin-root-age-18` (package code of commit 8b4a55f; later commits touch
only `data-raw/`).

## 1. Exact one-dimensional oracle (`data-raw/rung5-age-oracle.R`)

The tree is fixed: `((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7)`, with root age 4 and a stem of
0.6; a and c are fossils and f is a sampled ancestor. The model has K = 16, z in [-3, 3],
rho = 0.5 and extant conditioning. Parameters are fixed at the rung-5 prior medians. Only the
start age and the paths move (`tree` 1, `node` 2, `tip` 2, `branch` 2 and the age move 2 per
iteration), over 32 chains x 20000 iterations.

The reference is prior x `TreeLikelihood()` on a 321-point grid of the start age:

- crown start, root age with a uniform prior on [3.4, 7];
- origin start, origin age with a uniform prior on [4.2, 8].

z is the excess of the chains' CDF at each exact quantile over its nominal value, divided by
its standard error across chains. For the mean, z is the chain-mean difference over its
between-chain standard error.

| case | p | exact | mcmc | cdfAtExact | z |
|---|---|---|---|---|---|
| origin | q05 | 4.274 | 4.273 | 0.0505 | 0.92 |
| origin | q25 | 4.616 | 4.616 | 0.2498 | -0.16 |
| origin | q50 | 5.187 | 5.19 | 0.4988 | -0.76 |
| origin | q75 | 6.076 | 6.076 | 0.75 | 0.01 |
| origin | q95 | 7.418 | 7.422 | 0.9497 | -0.55 |
| origin | mean | 5.432 | 5.437 |  | 2.06 |
| root | q05 | 3.446 | 3.446 | 0.0503 | 0.55 |
| root | q25 | 3.658 | 3.658 | 0.2501 | 0.1 |
| root | q50 | 4.024 | 4.023 | 0.5008 | 0.65 |
| root | q75 | 4.647 | 4.639 | 0.7521 | 1.9 |
| root | q95 | 5.932 | 5.919 | 0.9509 | 1.42 |
| root | mean | 4.257 | 4.259 |  | 1.42 |

Both moves recover the exact posterior: the largest quantile |z| is 1.90. The origin mean
(z = 2.06) uses a Riemann sum over the grid; the CDF at all five quantiles agrees to
within |z| <= 0.92.

## 2. Prior-only MCMC vs forward simulation (`data-raw/rung5-age-validation.R`)

This is rung 5 (`data-raw/rung5-validation.R`): exactly 3 sampled extant tips and 2 fossil
samples, K = 16, z in [-3, 3], rho = 0.5, the rung-5 priors, and every tree move (`nodeAge` 2,
`fossilAge` 2, `slide` 2, `saJump` 2) with path moves. Fossil ranges are wide enough never to
bind. The three cases are:

- `root`: crown start with a free root age, prior U(3, 5), and the `rootAge` move;
- `origin`: origin fixed at 5, with the root moved by `nodeAge`;
- `originFree`: origin age with prior U(4.5, 6) and the `originAge` move.

In `fixed` runs the parameters are held at the prior medians. In `free` runs they are
sampled too. Each run uses 16 chains x 40000 iterations (thin 10) and 4000 (seed 0) or
8000 (later seeds) simulated trees.

The simulator draws the start age (and the parameters, when free) from the prior. It
resimulates at those values until the conditioning holds, then accepts the tree when the
tip and fossil counts match and, under an origin, the root is not a sampled ancestor.
Rejection may redraw both the start age and the parameters.

The tests are, per statistic, the mean and the CDF at five simulated quantiles
(q05-q95), plus the probabilities of each value of `nAncestor` and `cherries`. z uses the
spread over the 16 independent chains together with the simulation's binomial or sampling
error. Statistics: the free parameters, `treeLength`, `meanNodeAge`, `rootAge`,
`originAge`, `oldestFossil`, `youngestFossil`, `clock`, `inverseClock` (below the root),
`nAncestor` and `cherries`.

### Per run

| run | tests | max |z| |
|---|---|---|
| origin-fixed | 50 | 2.46 |
| origin-fixed-seed1000 | 50 | 2.08 |
| origin-fixed-seed2000 | 50 | 2.35 |
| origin-free | 80 | 2.61 |
| origin-free-seed1000 | 80 | 2.27 |
| originFree-fixed | 56 | 3.04 |
| originFree-fixed-seed1000 | 56 | 2.55 |
| originFree-fixed-seed2000 | 56 | 3.14 |
| originFree-fixed-seed3000 | 56 | 1.69 |
| originFree-fixed-seed4000 | 56 | 2.65 |
| originFree-fixed-seed5000 | 56 | 3.61 |
| originFree-free | 86 | 1.94 |
| originFree-free-seed1000 | 86 | 2.05 |
| root-fixed | 50 | 1.72 |
| root-fixed-seed1000 | 50 | 2.9 |
| root-fixed-seed2000 | 50 | 2.05 |
| root-free | 80 | 2.26 |
| root-free-seed1000 | 80 | 2.0 |

### Pooled over independent replicates (seed offsets 0, 1000, ...; z summed / sqrt(n))

| case | replicates | tests | max |z| pooled | at | Holm min p |
|---|---|---|---|---|---|
| origin-fixed | 3 | 50 | 2.25 | oldestFossil q75 | 0.978 |
| origin-free | 2 | 80 | 2.36 | youngestFossil q75 | 0.997 |
| originFree-fixed | 6 | 56 | 2.93 | rootAge q75 | 0.191 |
| originFree-free | 2 | 86 | 2.19 | inverseClock q95 | 0.983 |
| root-fixed | 3 | 50 | 2.25 | inverseClock q95 | 0.998 |
| root-free | 2 | 80 | 2.51 | meanNodeAge q75 | 0.972 |

Pooled z of the means of selected statistics:

| case | clock | lambda | meanNodeAge | oldestFossil | originAge | rootAge | theta | treeLength |
|---|---|---|---|---|---|---|---|---|
| origin-fixed | -0.66 |  | -1.08 | 1.57 |  | 0.03 |  | -0.92 |
| origin-free | 0.31 | -0.39 | -0.65 | 1.95 |  | 1.88 | 0.05 | 0.84 |
| originFree-fixed | -2.62 |  | -1.68 | -1.26 | -0.95 | -2.64 |  | -2.0 |
| originFree-free | 0.42 | -1.23 | 0.96 | -0.19 | 0.61 | 0.16 | -0.7 | 0.21 |
| root-fixed | 0.8 |  | -0.73 | 0.58 |  | -0.23 |  | -0.79 |
| root-free | 0.11 | 0.97 | -0.97 | 0.36 |  | 0.77 | 0.36 | -0.23 |

Reading these results:

- No case shows a significant discrepancy after a Holm correction over its tests. The
  smallest Holm-adjusted p is 0.19 (`originFree`, fixed parameters). Single runs reach
  |z| of up to 3.6. This is consistent with 50-86 correlated tests per run whose error comes
  from 16 chains.
- `originFree` with fixed parameters had the largest pooled values, all in one direction:
  `rootAge` mean -2.64, `clock` -2.62 and `treeLength` -2.0, which would mean slightly younger
  trees in the MCMC (about 1% of the root age). The signs varied between replicates. Section 3
  settles this flag: it was **noise**. The six replicates above are not independent of the
  flag, because the seed-0 run is what raised it.

## 3. The `originFree` flag: 240 further replicates (`data-raw/rung5-age-pool.R`)

The fixed-parameter `originFree` case was rerun with 240 new seed offsets (6000 to 245000 in
steps of 1000), in three batches decided in advance:

- A: 40 replicates, the test of the flag;
- B: 80, a confirmation of the largest value in A (`rootAge` q75);
- C: 120, a confirmation of the largest value in B (`clock` mean).

Controls were 10 replicates each of `origin` and `root` with fixed parameters (offsets 3000 to
12000), plus 10 of `originFree` with free parameters (2000 to 11000). Settings are as in
section 2 (16 chains x 40000 iterations, 8000 simulated trees), on Hamilton with the C++ engine
at commit b99e810, whose package code is that of 8b4a55f.

All chains of a case target the same distribution, and all simulated trees of a case come from
the same simulator. `rung5-age-pool.R` therefore pools them: each mean or CDF estimate from all
chains is compared with the estimate from all simulated trees. Its standard error combines the
spread of the independent chains with the simulation's sampling error. Results are in
`rung5-age-pooled.csv` and, per replicate, in `rung5-age-pooled-per-run.csv` (z of the means).

| case | replicates | chains | simulated trees | tests | max \|z\| pooled | at | Holm min p |
|---|---|---|---|---|---|---|---|
| originFree-fixed | 240 | 3840 | 1920000 | 56 | 1.44 | originAge q50 | 1.00 |
| originFree-free | 10 | 160 | 80000 | 86 | 2.43 | originAge q75 | 1.00 |
| origin-fixed | 10 | 160 | 80000 | 50 | 2.31 | oldestFossil q25 | 1.00 |
| root-fixed | 10 | 160 | 80000 | 50 | 1.37 | treeLength q50 | 1.00 |

`originFree-fixed`, pooled over all 240 replicates. The relative difference is
(MCMC - simulation) / simulation, given with its standard error:

| statistic | simulation | MCMC | relative difference | z |
|---|---|---|---|---|
| rootAge | 3.6939 | 3.6940 | +0.003% +- 0.024% | 0.13 |
| originAge | 5.2124 | 5.2121 | -0.006% +- 0.006% | -0.91 |
| treeLength | 9.1507 | 9.1503 | -0.004% +- 0.026% | -0.16 |
| meanNodeAge | 1.6617 | 1.6615 | -0.008% +- 0.044% | -0.18 |
| clock | 9.0304 | 9.0286 | -0.020% +- 0.035% | -0.56 |

Per-replicate z of the means, by batch (mean, sd; a calibrated test gives about 0 and 1.07):

| batch | replicates | rootAge | clock | treeLength |
|---|---|---|---|---|
| A | 40 | 0.09, 1.08 | 0.07, 0.95 | 0.07, 1.07 |
| B | 80 | -0.01, 1.01 | -0.26, 1.00 | -0.09, 0.98 |
| C | 120 | -0.00, 0.99 | 0.08, 1.03 | 0.02, 0.93 |

Reading these results:

- **Verdict: noise, not bias.** Each value that stood out failed to replicate in the next
  independent batch:
  - `rootAge` q75 had z = -3.10 over A (Holm p 0.11; the six earlier runs gave +2.93, the
    opposite sign), and -0.04 over B;
  - `clock` mean had z = -2.36 over B, and +0.87 over C.
- Pooled over all 240 replicates, no test exceeds |z| = 1.44.
- The 95% interval for the relative bias of the mean root age is (-0.04%, +0.05%). A 1% bias
  would have given z of about 40.
- The per-replicate z behave as calibrated tests should: their spread matches the t law
  expected from 16 chains. In batch A, the scatter of the simulated mean root age between
  replicates is 1.02 times its standard error. The wider scatter seen in the first six
  replicates (3.682 to 3.743) was chance.
- No change to the moves or the simulator was needed, so none was made.

## Conditioning choice

The start age is a parameter with its own prior, and the survival conditioning is
evaluated at the current start age: target = prior(start age) x p(tree, paths | Theta,
start age) / C(start age). The crown start uses extant conditioning of both root daughters;
the origin start uses that of the single stem lineage. This is documented in
`dev/spec/dCET-spec.md` section 6.1.

## Deviations and limits

- Section 3 commits only the pooled summaries, not the 270 per-replicate result files of
  `rung5-age-validation.R` (about 1.3 GB of `.rds`). Each one can be regenerated with
  `EVOTEMPO_SEED=<offset> Rscript data-raw/rung5-age-validation.R <case> fixed 40000 16 8000`.

- A sampled ancestor at the root (a fossil on the stem ancestral to all other samples) is
  not supported. The simulator rejects such trees, and `TempoMCMC()` stops on them.
- In the seed-0 runs of the fixed-origin case, `rung5-age-origin-*.csv` contain degenerate
  rows for `originAge`, which is constant there (NA or very large z). The script now skips a
  start age held fixed, and those rows are excluded from the tables above.
- All tests use one small tree size (3 extant + 2 fossil samples) and the rung-5 grid
  (K = 16). The default grid (K = 26) was not run.
- Start-age moves were not run with characters in these checks (rung 5 is prior-only).
  The engine tests cover that combination for R/C++ agreement only.
