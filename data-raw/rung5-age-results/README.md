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
was rerun with three other conditionings with a free start age: joint conditioning on the
origin start, and `condition = "survival"` on both starts
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
identical acceptance; `tests/testthat/test-Origin.R` now includes both cases. The probe's
output is `pr43-probe-output.txt` (the probe script itself was not committed).

**These runs have no power to tell the conditionings apart (#69).** At the fixture
parameters the survival probability barely depends on the start age. The largest shift
of an exact quantile from the extant-conditioned one is 0.0024 for origin-joint, 0.0010
for origin-survival and 0.0078 for root-survival, which is 0.53, 0.16 and 1.61 SE. So a
chain that ignored `joint` or `condition = "survival"` would pass this table as well.
Section 5 repeats the oracle at parameters where the conditionings differ, with a
negative control.

## 5. High-power start-age oracle with negative controls (#69)

Section 4 cannot tell the conditionings apart. This section reruns the oracle of section 1 at
parameters where they differ, for every conditioning that `TempoMCMC()` accepts with a free
start age, and scores each run against every other conditioning's exact posterior as a
negative control.

**Design** (`rung5-age-oracle.R <case> 60000 32 <variant> power`). Same tree, grid (K = 16,
z in [-3, 3]), moves and age priors as section 1, with lambda = mu = 1, psi = 0.15,
theta = 0.4, s = 1 and rho = 0.35. The process is critical in every bin, so the survival
probability falls steadily with the start age, and the bins spread over e^-3 to e^3, so
conditioning given the start bin (joint) differs from marginal conditioning. The exact
reference uses a 641-point grid. The variants are `extant`, `joint` (extant, joint),
`survival`, `survival-joint` and `none`, on both starts: 10 runs of 32 chains x 60000
iterations, seeds 118900 + i, 119900 + i, 120900 + i, 121900 + i and 122900 + i.
Each run took 12-20 min on 16 cores (Hamilton, EvoTempo 589f866 with this commit's
`data-raw/`).

**Power, before running** (`rung5-age-oracle.R <case> 60000 32 design power`;
`rung5-age-oracle-<case>-power-design.csv`, and the exact posteriors in `-design-exact.csv`).
For chains targeting A, scored against B's exact posterior, the expected |z| is
|F_A(q_B(p)) - p| / sqrt(p(1 - p) / ESS) at B's quantiles, and the difference of means over
sd_A / sqrt(ESS). The assumed ESS is 0.078 per chain iteration, half the rate observed in
section 4. Over the one-flag pairs (one conditioning flag changed: joint and marginal,
extant and survival, extant and none), the smallest expected max |z| was 14.5 (origin) and
30.3 (root). At the section-4 parameters it is 0.20 and 0.90
(`rung5-age-oracle-<case>-design.csv`).

**Criteria, and when they were set** (times UTC, 3 November 2026; Slurm accounting, file
times and the scripts that ran are in `power-oracle-provenance.txt`).

- *Negative controls, fixed before the runs.* The header of `rung5-age-oracle.R` in the copy
  that ran (the Hamilton overlay of 08:10:55, on 589f866; the design ran at 08:22-08:23 and
  the ten chain runs were submitted at 08:54:22) requires every run to reject each one-flag
  alternative with max |z| >= 5. That is the only bar written down before any chain ran.
- *Own target, chosen after the results.* No bar for a run against its own exact posterior
  was written down beforehand. The bar used below, |z| < 3.5 for all six statistics, was
  chosen after all ten runs had finished. Root-none (largest |z| 3.30) finished at 09:08:50
  and root-extant (3.01) at 09:14:04, and the bar first appears in this README in commit
  dc47aed (10:40:39). A bar of 3.0 would have failed both runs. The own-target tables are
  therefore a post hoc reading, not a pre-registered test.
- *False alarms.* Each z divides by the between-chain standard error of 32 chains, so its
  reference distribution is close to t with 31 degrees of freedom, not N(0, 1). Under t31,
  P(|z| >= 3.5) = 0.0014 per statistic. If every run were exact, the chance that at least
  one of the 60 statistics of the ten runs reaches 3.5 would be about 8.2%: 12.9% over the
  96 statistics that include the six root-extant replicates below, and 13.6% with the
  root-none replicate as well. (An earlier version of this section gave 3%, the
  normal-theory figure.) Statistics within a run are positively correlated, so the true
  rates are somewhat lower than these independence calculations.

**Own target** (`rung5-age-oracle-<case>-power-<variant>.csv`): z of the chains against their
own conditioning's exact posterior.

Origin start:

| conditioning | q05 | q25 | q50 | q75 | q95 | mean |
|---|---|---|---|---|---|---|
| extant | 0.06 | 0.17 | -0.29 | 0.90 | 1.74 | -0.64 |
| joint | 0.61 | -0.69 | -0.92 | -0.20 | -1.66 | 1.91 |
| survival | 0.92 | 1.12 | 1.97 | 1.32 | 1.15 | -1.79 |
| survival-joint | -1.96 | -1.96 | -1.30 | 0.07 | 0.98 | 0.96 |
| none | 1.36 | 0.76 | -0.64 | -1.45 | -2.10 | 1.04 |

Crown (root) start:

| conditioning | q05 | q25 | q50 | q75 | q95 | mean |
|---|---|---|---|---|---|---|
| extant | 1.25 | 0.32 | 1.01 | 2.38 | 3.01 | -2.46 |
| joint | 2.03 | 0.91 | 0.42 | -0.28 | -0.86 | -0.40 |
| survival | -1.12 | -1.12 | -1.46 | -0.69 | -0.10 | 1.16 |
| survival-joint | 1.02 | 0.81 | 0.54 | -0.46 | 0.07 | -0.19 |
| none | -0.08 | -0.24 | -3.30 | -1.00 | 0.80 | 1.61 |

The largest |z| is 3.30 (root, none, q50), so every run is inside the post hoc bar of 3.5.

**Replicates: sequential, and decided on the results.** Root-extant had a coherent pattern
(q75 2.38, q95 3.01, mean -2.46; the section-1 run at the old design has the same signs).
Root-none had the largest single |z| (3.30), so it was rerun as well. None of these reruns
was planned before the first results. They ran in two batches, each decided after seeing
the results before it (`EVOTEMPO_RUNG5_SEED_OFFSET`; files `-power-<variant>-seed<offset>.csv`):

1. `confirm.sh`, submitted at 09:17:20, after the root-none (09:08:50) and root-extant
   (09:14:04) results: root-extant at seed offsets 50000 and 60000, and root-none at 50000.
2. `confirm2.sh`, submitted at 09:38:56, after the first two root-extant replicates had
   finished (09:34:47 and 09:37:18): root-extant at 70000, 80000, 90000 and 100000.

That makes six root-extant replicates (192 chains) and one root-none replicate (32 chains).

| replicate (seed offset) | q05 | q25 | q50 | q75 | q95 | mean |
|---|---|---|---|---|---|---|
| 0 | 1.25 | 0.32 | 1.01 | 2.38 | 3.01 | -2.46 |
| 50000 | -1.40 | 0.67 | -0.25 | 0.59 | 1.07 | -0.44 |
| 60000 | 1.41 | 1.75 | 1.17 | 1.18 | 0.74 | -1.89 |
| 70000 | 1.57 | 1.37 | 0.75 | 1.10 | 1.21 | -1.70 |
| 80000 | 1.34 | -0.54 | -0.80 | 0.72 | -0.42 | 0.13 |
| 90000 | -0.20 | 1.13 | 0.53 | 0.23 | 0.51 | -0.75 |
| 100000 | 0.60 | 0.09 | 0.51 | -1.30 | -2.82 | 1.10 |
| Stouffer, the 6 new replicates | 1.36 | 1.83 | 0.78 | 1.03 | 0.12 | -1.45 |
| Stouffer, all 7 | 1.73 | 1.81 | 1.11 | 1.85 | 1.25 | -2.27 |

The root-none replicate gave max |z| 0.98.

The six new root-extant replicates are weaker than the first run (largest Stouffer |z|
1.83), but they keep its direction in every statistic: positive z at the quantiles and a
negative z for the mean, so the chains sit slightly young. The inverse-variance pooled mean
bias (chain mean minus exact mean) is:

- over the 6 unselected replicates: -0.00104 +- 0.00065 (z -1.60; 95% interval -0.0023 to
  +0.0002);
- over all 7 runs: -0.00144 +- 0.00060 (z -2.41).

The 7-run figure includes the run that triggered the replication, so selection inflates
it; the 6-run figure is the one to quote. Either way the bias is small, 0.12-0.16% of the
posterior sd (0.89), but it is not excluded. It is tracked in #72.

**Negative controls** (`-cross.csv`): max |z| over the six statistics when each run's chains
are scored against each conditioning's exact posterior. The required one-flag alternatives
are in bold and the run's own target is in brackets.

Origin start:

| chains \ scored against | extant | joint | survival | survival-joint | none |
|---|---|---|---|---|---|
| extant | (1.7) | **22.3** | **48.4** | 69.3 | **41.8** |
| joint | **27.3** | (1.9) | 35.3 | **62.0** | 26.9 |
| survival | **52.2** | 28.6 | (2.0) | **19.7** | 8.6 |
| survival-joint | 84.3 | **56.4** | **24.5** | (2.0) | 32.6 |
| none | **52.1** | 23.4 | 9.3 | 35.4 | (2.1) |

Crown (root) start:

| chains \ scored against | extant | joint | survival | survival-joint | none |
|---|---|---|---|---|---|
| extant | (3.0) | **47.2** | **85.5** | 125.6 | **68.7** |
| joint | **53.1** | (2.0) | 40.1 | **82.6** | 22.3 |
| survival | **108.4** | 46.5 | (1.5) | **51.2** | 19.8 |
| survival-joint | 147.3 | **90.2** | **46.3** | (1.0) | 65.6 |
| none | **91.6** | 26.5 | 23.6 | 76.3 | (3.3) |

Every one-flag alternative is rejected, the weakest at 19.7 (origin survival chains against
survival-joint). The realised max |z| are 1.36-1.93 times the design's predictions on the
one-flag pairs and 1.40-2.04 times on the two-flag pairs, because the design assumed half
the observed ESS rate. Even pairs that differ in two flags, which the
design did not require, are rejected at 8.6 or more.

**What this shows, and what it does not.** With a free start age, for crown and origin
starts under extant, joint, survival, survival + joint and `condition = "none"`:

- every run rejects each one-flag alternative, at the bar fixed before the runs;
- every run matches its own exact posterior within the bar chosen after the results.

This covers the combinations that section 4 left untested (crown + joint, joint +
survival, none). The exact reference is the package's own conditioned `TreeLikelihood()`,
so the oracle shows that `TempoMCMC()`'s start-age moves sample that density and that each
conditioning flag reaches the sampler. It does not show that the conditioning formulas
themselves are right. Those are checked independently only in part:

- extant conditioning, against forward simulation, in section 2;
- `condition = "none"`, against simulation, in rung 3 (`data-raw/rung3-validation.R`,
  Monte Carlo on tiny trees);
- at beta = 0, `none`, `extant` and `survival` reduce to the closed-form FBD likelihood
  (`tests/testthat/test-TreeLikelihood.R`).

Survival conditioning and joint conditioning with rates that vary across bins have no
check against simulation. Further limits: fixed parameters, one small tree, no characters,
K = 16, and the small root-extant bias above (#72).

## Conditioning choice

The start age is a parameter with its own prior, and the survival conditioning is
evaluated at the current start age: target = prior(start age) x p(tree, paths | Theta,
start age) / C(start age). The crown start uses extant conditioning of both root daughters;
the origin start uses that of the single stem lineage. This is documented in
`dev/spec/dCET-spec.md` section 6.1.

## Deviations and limits

- Section 5: the power design was chosen with the help of a Python port of the E/D
  equations. The port and its comparison with `TreeLikelihood()` were not committed, so no
  agreement is claimed for them here. An earlier version of section 5 claimed agreement to
  1e-10; that claim is withdrawn.

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
