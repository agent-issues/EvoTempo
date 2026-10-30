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
- `originFree` with fixed parameters has the largest pooled values, all in one direction:
  `rootAge` mean -2.64, `clock` -2.62 and `treeLength` -2.0, meaning slightly younger trees
  in the MCMC. That case was replicated six times for this reason. The signs vary between
  replicates (seed 1000 is positive), and the exact oracle above, which isolates
  `originAge`, shows no bias. This is recorded as an **unresolved flag**, not a verified
  discrepancy. If it is real, the effect is about 1% of the root age.

## Conditioning choice

The start age is a parameter with its own prior, and the survival conditioning is
evaluated at the current start age: target = prior(start age) x p(tree, paths | Theta,
start age) / C(start age). The crown start uses extant conditioning of both root daughters;
the origin start uses that of the single stem lineage. This is documented in
`dev/spec/dCET-spec.md` section 6.1.

## Deviations and limits

- A sampled ancestor at the root (a fossil on the stem ancestral to all other samples) is
  not supported. The simulator rejects such trees, and `TempoMCMC()` stops on them.
- In the seed-0 runs of the fixed-origin case, `rung5-age-origin-*.csv` contain degenerate
  rows for `originAge`, which is constant there (NA or very large z). The script now skips a
  start age held fixed, and those rows are excluded from the tables above.
- All tests use one small tree size (3 extant + 2 fossil samples) and the rung-5 grid
  (K = 16). The default grid (K = 26) was not run.
- Start-age moves were not run with characters in these checks (rung 5 is prior-only).
  The engine tests cover that combination for R/C++ agreement only.
