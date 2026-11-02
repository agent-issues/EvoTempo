# Rung 5 with a free start age (#18)

Validation of the origin start (`start = "origin"`: stem path, origin bin, single-lineage
survival) and of the start-age moves (`rootAge`, `originAge`, and `nodeAge` on the root under an
origin) added to `TempoMCMC()` in this PR. Everything here was run on Hamilton with the C++
engine, at branch `origin-root-age-18` (package code of commit 8b4a55f, unchanged up to
fbfc5bc; the later start-up guards of #55 and sample columns of #56 change neither the moves
nor the target).

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
| origin | mean | 5.436 | 5.437 |  | 0.51 |
| root | q05 | 3.446 | 3.446 | 0.0503 | 0.55 |
| root | q25 | 3.658 | 3.658 | 0.2501 | 0.1 |
| root | q50 | 4.024 | 4.023 | 0.5008 | 0.65 |
| root | q75 | 4.647 | 4.639 | 0.7521 | 1.9 |
| root | q95 | 5.932 | 5.919 | 0.9509 | 1.42 |
| root | mean | 4.262 | 4.259 |  | -1.65 |

Both moves recover the exact posterior: the largest |z| is 1.90 (root q75). The exact means
use the trapezoid rule over the grid, as the CDF does. The first version of this table used
an equal-weight sum, which overweights the ends of the grid and gave exact means of 5.432
(origin, z = 2.06) and 4.257 (root, z = 1.42); the mean rows of the CSVs were recomputed
from the same 321-point grid (`TreeLikelihood()` is deterministic, and the recomputed
equal-weight means reproduce the stored ones exactly), with the chains' standard errors
unchanged (#57).

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
| origin-fixed | 3 | 50 | 2.25 | oldestFossil q75 | 1.00 |
| origin-free | 2 | 80 | 2.36 | youngestFossil q75 | 1.00 |
| originFree-fixed | 6 | 56 | 2.93 | rootAge q75 | 0.191 |
| originFree-free | 2 | 86 | 2.19 | inverseClock q95 | 1.00 |
| root-fixed | 3 | 50 | 2.25 | inverseClock q95 | 1.00 |
| root-free | 2 | 80 | 2.51 | meanNodeAge q75 | 0.972 |

Holm min p is the Holm-adjusted two-sided p of the largest pooled |z|, min(1, m p) over the
case's m tests, from the per-run CSVs (pooled z = sum of the per-run z / sqrt(n); tests
present in every replicate). Four values in the first version of this table (0.978, 0.997,
0.983, 0.998) did not reproduce from the CSVs; m p is 1.22, 1.48, 2.43 and 1.22 there (#57).

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
  `rootAge` mean -2.64, `clock` -2.62 and `treeLength` -2.0. The gap in the mean root age
  averaged over the six replicates is -0.45% of the root age (not about 1%, as first stated),
  and it lies on the **simulation** side: the six replicates' MCMC means (3.690 to 3.701) all
  agree with the 240-replicate pooled mean of section 3 (3.6940), while their simulated means
  scatter from 3.682 to 3.743 around the pooled simulated mean (3.6939). Section 3 settles
  this flag: it was **noise**, likeliest chance plus selection (the seed-0 run is what raised
  it, so the six replicates are not independent of the flag).

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
- Detectable bias (half-width of the 95% interval for the relative bias of the mean, 1.96
  standard errors): with all 240 replicates, 0.05% for `rootAge` and 0.01% for `originAge`;
  with the 10-replicate controls, 0.11-0.23% (`rootAge`) and 0.06% (`originAge`,
  `originFree`-free); with the 2-6 replicates of section 2, 0.22-0.62% (`rootAge`) and
  0.09-0.16% (`originAge`). Section 2 alone therefore could not resolve a bias of the size
  first suspected; section 3 can.

The batch-level values quoted above for `rootAge` q75 (A, B) and the `clock` mean (B, C) were
computed by `rung5-age-pool.R` over each batch's per-replicate `.rds` files, which are no
longer kept (see below), so they cannot be regenerated from committed files.
`rung5-age-batches.R` writes what the committed per-replicate z of the means support,
`rung5-age-pooled-per-batch.csv`: per batch and statistic, the mean and sd of the z and their
Stouffer sum (z sum / sqrt(n)). For the `clock` mean these are -2.31 (B) and +0.90 (C), close
to the pooled -2.36 and +0.87; quantile tests per batch are not available (#57).

## 4. Start-age moves under joint and survival conditioning (#54)

Sections 1-3 use extant conditioning, marginal over the start bin. The oracle of section 1
was rerun with the other conditionings that `TempoMCMC()` permits with a free start age:
joint conditioning on the origin start, and `condition = "survival"` on both starts
(`Rscript data-raw/rung5-age-oracle.R <root|origin> 20000 32 <joint|survival>`; 32 chains x
20000 iterations each, seeds 19900 + i and 20900 + i, otherwise as in section 1). Results are
in `rung5-age-oracle-<case>-<variant>.csv`.

| variant | p | exact | mcmc | cdfAtExact | z |
|---|---|---|---|---|---|
| origin-joint | q05 | 4.275 | 4.275 | 0.0500 | 0.01 |
| origin-joint | q25 | 4.617 | 4.620 | 0.2485 | -1.11 |
| origin-joint | q50 | 5.188 | 5.194 | 0.4981 | -1.14 |
| origin-joint | q75 | 6.078 | 6.080 | 0.7497 | -0.27 |
| origin-joint | q95 | 7.420 | 7.424 | 0.9495 | -0.70 |
| origin-joint | mean | 5.437 | 5.439 |  | 0.65 |
| root-survival | q05 | 3.445 | 3.445 | 0.0499 | -0.15 |
| root-survival | q25 | 3.655 | 3.657 | 0.2490 | -0.79 |
| root-survival | q50 | 4.019 | 4.022 | 0.4985 | -0.99 |
| root-survival | q75 | 4.639 | 4.640 | 0.7499 | -0.07 |
| root-survival | q95 | 5.925 | 5.924 | 0.9501 | 0.16 |
| root-survival | mean | 4.257 | 4.258 |  | 0.23 |
| origin-survival | q05 | 4.274 | 4.273 | 0.0505 | 0.77 |
| origin-survival | q25 | 4.615 | 4.617 | 0.2490 | -0.78 |
| origin-survival | q50 | 5.186 | 5.188 | 0.4989 | -0.61 |
| origin-survival | q75 | 6.075 | 6.070 | 0.7512 | 0.81 |
| origin-survival | q95 | 7.419 | 7.421 | 0.9498 | -0.35 |
| origin-survival | mean | 5.435 | 5.434 |  | -0.38 |

The largest of the 18 |z| is 1.14. An earlier, independent run of the same comparison by the
reviewers' probe script (24 chains x 15000 iterations, seeds 43001-43024, at the merge of
`main` into this branch) gave a largest |z| of 1.85 (origin-joint q05). In the same probe,
seeded chains with characters, every tree move and an age move agreed between the R and
C++ engines under joint (origin) and survival (crown) conditioning, to 3e-11 relative, with
identical acceptance; `tests/testthat/test-Origin.R` now includes both cases.

## Conditioning choice

The start age is a parameter with its own prior, and the survival conditioning is
evaluated at the current start age: target = prior(start age) x p(tree, paths | Theta,
start age) / C(start age). The crown start uses extant conditioning of both root daughters;
the origin start uses that of the single stem lineage. This is documented in
`dev/spec/dCET-spec.md` section 6.1.

## Deviations and limits

- Section 3 commits only the pooled summaries, not the 270 per-replicate result files of
  `rung5-age-validation.R` (about 1.3 GB of `.rds`), which were deleted with the Hamilton
  scratch directory. Each one can be regenerated with
  `EVOTEMPO_SEED=<offset> Rscript data-raw/rung5-age-validation.R <case> fixed 40000 16 8000`.
- Provenance of the seed-0 runs of section 2: they cannot have used the script as committed
  in 8b4a55f, where `perWorker` is not exported to the PSOCK workers (fixed in c0f3b36), so
  they must have run a modified copy. The package code (`R/`, `src/`, `tests/`) is identical
  from 8b4a55f to fbfc5bc, the head of this PR before the #57 corrections. The `sha` recorded in their `.rds` files cannot now be read, because those
  files were deleted with the scratch directory (#57).
- `.AgeMove()` rebuilds the frozen model at the current tree whenever a tree move has changed
  a node age since the model was last built (the key of `.SameSet()` includes node ages,
  which are points of the E grid). This costs one rebuild per age move after an accepted tree move;
  it is needed for the reverse proposal density and was left unchanged (#57, F11).

- A sampled ancestor at the root (a fossil on the stem ancestral to all other samples) is
  not supported. The simulator rejects such trees, and `TempoMCMC()` stops on them.
- In the seed-0 runs of the fixed-origin case, `rung5-age-origin-*.csv` contain degenerate
  rows for `originAge`, which is constant there (NA or very large z). The script now skips a
  start age held fixed, and those rows are excluded from the tables above.
- All tests use one small tree size (3 extant + 2 fossil samples) and the rung-5 grid
  (K = 16). The default grid (K = 26) was not run.
- Start-age moves were not run with characters in these checks (rung 5 is prior-only).
  The engine tests cover that combination for R/C++ agreement only.
