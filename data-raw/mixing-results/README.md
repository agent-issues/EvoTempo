# Mixing results for #19: shift move and move schedule

Validation of the `shift` move and the `TempoMoves()` default schedule added to
`TempoMCMC()` in this PR. The specification and the prototype study are in
`dev/mixing/README.md` (#41). Everything here was run on Hamilton with the C++ engine,
at the branch head `shift-move-19` (package code of commit 5f099ad).

Files:

- `mixing-schedule-{chains,gains,rhat}.csv`: ESS per CPU hour, old vs new schedule
  (`data-raw/mixing-schedule.R`);
- `rung4b-shift-fixed.csv`: shift move at fixed parameters vs the exact root-bin posterior
  (`data-raw/rung4b-shift-fixed.R`);
- `shift-oracles-between-chain.csv`: one-parameter oracles with the shift move
  (`data-raw/rung4b-one-parameter.R ... shift`);
- `sbc/`: the short rung-6 SBC with the new schedule (`data-raw/rung6-sbc-full.R` with
  `EVOTEMPO_SBC_SCHEDULE=default`; summary by `data-raw/rung6-sbc-full-summary.R`).

The `.rds` outputs (chains, oracle draws) are not committed.

## 1. The shift move is exact

### One free parameter at a time (rung 4b tree, one character)

Each run frees one parameter, with path moves (`tree` 1, `node` 2, `tip` 2, `branch` 2), the
random-walk sweep and one `shift` per iteration. The shift then compensates only that
parameter. The reference is the exact one-dimensional posterior from a 161-point
`TreeLikelihood(character = )` grid. z is the excess of the chains' CDF at each exact
quantile over its nominal value, divided by its standard error across independent chains.

| run | parameter | nChain | nIter | z q05 | z q25 | z q50 | z q75 | z q95 |
|---|---|---|---|---|---|---|---|---|
| four | lambda | 16 | 12000 | -2.86 | -2.18 | -0.81 | -0.85 | -0.45 |
| four | mu | 16 | 12000 | -0.75 | -0.6 | 0.22 | -0.05 | -1.32 |
| four | rate | 16 | 12000 | 1.21 | 0.95 | -2.07 | -1.65 | -1.18 |
| four | theta | 16 | 12000 | -0.23 | -0.75 | -1.06 | -0.67 | -0.34 |
| lambdaControl | lambda | 16 | 24000 | 0.44 | -1.17 | -0.59 | -0.99 | 0.17 |
| lambdaShift | lambda | 16 | 24000 | -0.23 | 1.03 | 0.58 | 0.75 | 2.2 |
| lambdaShift64 | lambda | 64 | 24000 | -0.79 | -0.88 | -0.22 | -0.92 | 0.41 |

- `four` is one run with `rate`, `lambda`, `mu` and `theta` freed in turn (16 chains each).
  Its `lambda` lower tail (z = -2.86 at q05) was followed up with two independent
  replicates. A 16-chain replicate gave |z| <= 2.20, with its largest value in the upper tail.
  A 64-chain run gave max |z| = 0.92. A control without the shift move (`lambdaControl`,
  16 chains) gave max |z| = 1.17. Neither replicate confirms the q05 value, so it is read
  as a fluctuation.
- The batch-means z-scores that the script prints (within-chain) are in the logs; they
  agree except for the same `lambda` q05 value (3.48 in `four`, 0.23 in the 16-chain
  replicate).

### Fixed parameters: root-bin posterior

The tree and character are those of rung 4b, with parameters at the truth. The `shift` chains
use `node` 2, `tip` 2, `branch` 2 and `shift` 2 per iteration; with no parameter free, the
shift moves bins only. The `reference` chains use `tree` 1, `node` 2, `tip` 2 and `branch` 2.
Both run 16 chains x 20000 iterations. The exact root-bin posterior is
p(k | T, c) proportional to pi_k S_k exp(L_k), where L_k is `TreeLikelihood(character = )`
with a point-mass start prior on bin k. The MCSE is from batch means.

| config | quantity | exact | mcmc | mcse | z |
|---|---|---|---|---|---|
| shift | P(rootBin=2) | 0.001449 | 0.001309 | 0.0001287 | -1.09 |
| shift | P(rootBin=3) | 0.01712 | 0.01719 | 0.0008792 | 0.08 |
| shift | P(rootBin=4) | 0.08713 | 0.08732 | 0.002481 | 0.08 |
| shift | P(rootBin=5) | 0.2101 | 0.2121 | 0.00353 | 0.58 |
| shift | P(rootBin=6) | 0.274 | 0.2771 | 0.00262 | 1.2 |
| shift | P(rootBin=7) | 0.2201 | 0.219 | 0.002854 | -0.4 |
| shift | P(rootBin=8) | 0.121 | 0.1198 | 0.002832 | -0.43 |
| shift | P(rootBin=9) | 0.04915 | 0.04795 | 0.001875 | -0.64 |
| shift | P(rootBin=10) | 0.01563 | 0.01439 | 0.0009233 | -1.35 |
| shift | P(rootBin=11) | 0.004237 | 0.003719 | 0.0003698 | -1.4 |
| shift | E(rootBin) | 6.252 | 6.235 | 0.02433 | -0.69 |
| reference | P(rootBin=2) | 0.001449 | 0.001494 | 8.858e-05 | 0.5 |
| reference | P(rootBin=3) | 0.01712 | 0.0175 | 0.0002773 | 1.39 |
| reference | P(rootBin=4) | 0.08713 | 0.08768 | 0.0006097 | 0.91 |
| reference | P(rootBin=5) | 0.2101 | 0.2103 | 0.0008841 | 0.25 |
| reference | P(rootBin=6) | 0.274 | 0.2751 | 0.001004 | 1.15 |
| reference | P(rootBin=7) | 0.2201 | 0.2191 | 0.0008506 | -1.19 |
| reference | P(rootBin=8) | 0.121 | 0.1197 | 0.0006136 | -2.21 |
| reference | P(rootBin=9) | 0.04915 | 0.04894 | 0.0004264 | -0.48 |
| reference | P(rootBin=10) | 0.01563 | 0.01605 | 0.0002357 | 1.78 |
| reference | P(rootBin=11) | 0.004237 | 0.004059 | 0.0001245 | -1.43 |
| reference | E(rootBin) | 6.252 | 6.245 | 0.003031 | -2.01 |

Shift acceptance was 0.25 (burn-in included, as for every acceptance rate and `acc.*` column
in this directory). The mean log-tempo of the two samplers agrees
(z = -0.82).

## 2. Calibration is unchanged: short rung-6 SBC

These are the gamma = 1 cells of the full rung-6 SBC (#21): the same replicate seeds
(g = 2, 4, ..., 200, i.e. the first 50 replicates of cells 2 and 4), the same data, tree
moves and ESS-targeted chain extension. The only change is the path schedule, now
`TempoMoves()` (no `tree`, no `joint`; `node`, `tip` and `branch` at 0.25 x the
nodes, tips and branches; one `shift`) plus the same tree moves. There are 30 characters,
K = 16 and z in [-3, 3].

| label | nReplicate | nPass | nQuantity | pass | minPChisqHolm | minPKSHolm | fracReplicatesAllEss50 | medianIter | fracExtended | coreHours |
|---|---|---|---|---|---|---|---|---|---|---|
| beta=0, gamma=1 | 50 | 13 | 13 | True | 0.8 | 0.15 | 1 | 20000 | 0.2 | 27.13 |
| beta=1, gamma=1 | 50 | 13 | 13 | True | 1.0 | 1.0 | 1 | 40000 | 0.54 | 39.64 |

All 13 ranked quantities pass in both cells (chi-square and KS, Holm-adjusted). Every
replicate reached ESS >= 50 on every quantity. The two cells took 66.76 core-hours in total
(`runtimeMin` in `sbc/rung6-full-runs.csv`).

The same 100 replicates under the #21 schedule, from `data-raw/rung6-full-results/`:

| schedule | cell | replicates | median minutes | total core-hours | extended | median min ESS |
|---|---|---|---|---|---|---|
| #21 (tree, node 2, tip 2, branch 2, joint 1) | beta=0, gamma=1 | 50 | 39.2 | 40.9 | 0.48 | 67.8 |
| #21 (tree, node 2, tip 2, branch 2, joint 1) | beta=1, gamma=1 | 50 | 62.6 | 60.4 | 0.76 | 64.1 |
| TempoMoves() (this PR) | beta=0, gamma=1 | 50 | 28.0 | 27.1 | 0.20 | 87.6 |
| TempoMoves() (this PR) | beta=1, gamma=1 | 50 | 48.1 | 39.6 | 0.54 | 74.9 |

Median ESS per CPU hour by quantity, from the same replicates (ESS of the stored draws per
hour of chain run time):

| cell | quantity | oldESSperH | newESSperH | oldESSper1k | newESSper1k | gain |
|---|---|---|---|---|---|---|
| beta=0, gamma=1 | clock | 272 | 579 | 6.21 | 13.42 | 2.13 |
| beta=0, gamma=1 | lambda | 1307 | 1432 | 30.63 | 29.13 | 1.1 |
| beta=0, gamma=1 | meanNodeAge | 366 | 409 | 9.02 | 8.43 | 1.12 |
| beta=0, gamma=1 | rate | 283 | 612 | 6.76 | 13.32 | 2.16 |
| beta=0, gamma=1 | rootBin | 471 | 587 | 10.68 | 12.82 | 1.25 |
| beta=0, gamma=1 | s | 703 | 726 | 16.56 | 15.98 | 1.03 |
| beta=0, gamma=1 | theta | 104 | 177 | 2.65 | 4.22 | 1.71 |
| beta=0, gamma=1 | treeLength | 341 | 373 | 8.28 | 8.04 | 1.09 |
| beta=1, gamma=1 | clock | 242 | 560 | 6.12 | 14.02 | 2.31 |
| beta=1, gamma=1 | lambda | 549 | 1229 | 14.35 | 30.26 | 2.24 |
| beta=1, gamma=1 | meanNodeAge | 143 | 170 | 3.82 | 4.47 | 1.19 |
| beta=1, gamma=1 | rate | 304 | 585 | 7.86 | 14.17 | 1.93 |
| beta=1, gamma=1 | rootBin | 515 | 765 | 13.0 | 17.1 | 1.49 |
| beta=1, gamma=1 | s | 312 | 406 | 8.52 | 9.69 | 1.3 |
| beta=1, gamma=1 | theta | 68 | 125 | 1.65 | 3.08 | 1.85 |
| beta=1, gamma=1 | treeLength | 133 | 167 | 3.58 | 4.23 | 1.26 |

The #21 runs used package commit a2c3aa1. Between a2c3aa1 and 5f099ad (#42), `src/` is
unchanged and the only other `R/` change is #39's solver-failure handling, so run times are
confounded only by node; they are therefore only indicative. The ESS-per-1000-iterations
columns do not depend on run time.

## 3. ESS per CPU hour: old schedule vs `TempoMoves()`

The fixtures are the two fixed true trees of `dev/mixing/` (small: 17 extant tips; large: 75),
cell beta = gamma = 1, K = 16, 30 / 100 / 300 nested characters, with all six parameters free.
Chains start at the truth and run 2000 burn-in plus 5000 sampling iterations, 1000 stored,
with 3 replicates per cell. `old` is the schedule used before this PR (`tree` 1, `node` 2,
`tip` 2, `branch` 2, `joint` 1, `minAcceptance = 0`). `new` is the default `TempoMoves()`.
CPU time is the whole `TempoMCMC()` call (user + system), scaled to the sampling phase.
ESS is Geyer's initial monotone sequence estimator. Unlike `dev/mixing/`, these chains run
the engine's own `TempoMCMC()`, not the prototype driver.

### Gain in ESS per CPU hour (new / old; ratio of replicate means)

| tree | nChar | theta | s | c (rate) | lambda | rootBin | clock | logPost |
|---|---|---|---|---|---|---|---|---|
| small | 30 | 3.1 | 1.6 | 4.8 | 5.5 | 1.2 | 5 | 3.2 |
| small | 100 | 2.5 | 1.1 | 13 | 4.1 | 3.5 | 16 | 1.4 |
| small | 300 | 4.2 | 1.9 | 32 | 2 | 14 | 44 | 1.7 |
| large | 30 | 0.79 | 0.42 | 1.5 | 0.39 | 1.8 | 2.9 | 1.7 |
| large | 100 | 1 | 0.42 | 1.8 | 0.37 | 1.1 | 2.5 | 0.6 |
| large | 300 | 1.3 | 0.31 | 1.5 | 0.25 | 1.4 | 4.5 | 1.5 |

### ESS per CPU hour (old / new)

| tree | nChar | theta | s | rate | lambda | rootBin | clock | logPost |
|---|---|---|---|---|---|---|---|---|
| small | 30 | 308 / 958 | 894 / 1441 | 1205 / 5780 | 2679 / 14795 | 1531 / 1775 | 1005 / 5017 | 377 / 1223 |
| small | 100 | 141 / 352 | 364 / 385 | 211 / 2698 | 2742 / 11305 | 372 / 1320 | 151 / 2408 | 267 / 381 |
| small | 300 | 149 / 621 | 107 / 207 | 107 / 3387 | 5050 / 10307 | 123 / 1765 | 85 / 3793 | 594 / 1007 |
| large | 30 | 48 / 38 | 301 / 125 | 199 / 301 | 1320 / 513 | 106 / 186 | 69 / 199 | 50 / 86 |
| large | 100 | 48 / 48 | 321 / 136 | 196 / 351 | 3168 / 1187 | 187 / 212 | 127 / 312 | 165 / 100 |
| large | 300 | 30 / 41 | 156 / 48 | 178 / 270 | 2249 / 560 | 174 / 244 | 57 / 257 | 59 / 86 |

### CPU seconds per iteration

| tree | nChar | secPerIterOld | secPerIterNew | ratio |
|---|---|---|---|---|
| small | 30 | 0.046 | 0.035 | 0.77 |
| small | 100 | 0.065 | 0.048 | 0.74 |
| small | 300 | 0.065 | 0.052 | 0.79 |
| large | 30 | 0.111 | 0.179 | 1.62 |
| large | 100 | 0.098 | 0.189 | 1.92 |
| large | 300 | 0.113 | 0.283 | 2.51 |

Per iteration, the new schedule runs 137 local moves on the large tree (`node` 42, `tip` 27,
`branch` 68) and 36 on the small tree (11, 7, 18), each with a full `MkLikelihood()` re-prune,
against 6 under `old`. On the large tree it costs 1.6-2.5x more per iteration than `old`. On the
small tree it is cheaper, because the dropped `joint` sweep (one FFBS draw per parameter) costs
more than the extra local moves. At 300 characters the acceptance rates (burn-in included, as
`TempoMCMC()` then reported them) were `node` 0.76-0.80,
`tip` 0.78-0.85, `branch` 0.90-0.92 and `shift` 0.33-0.41.

### Posterior sd of log parameters (old / new; mean of replicates)

| tree | nChar | theta | s | rate | lambda | clock |
|---|---|---|---|---|---|---|
| small | 30 | 0.410 / 0.399 | 0.293 / 0.268 | 0.198 / 0.193 | 0.187 / 0.184 | 0.180 / 0.174 |
| small | 100 | 0.362 / 0.378 | 0.184 / 0.206 | 0.115 / 0.178 | 0.159 / 0.180 | 0.100 / 0.171 |
| small | 300 | 0.295 / 0.353 | 0.166 / 0.180 | 0.082 / 0.170 | 0.164 / 0.181 | 0.077 / 0.168 |
| large | 30 | 0.191 / 0.408 | 0.183 / 0.228 | 0.075 / 0.106 | 0.104 / 0.120 | 0.061 / 0.098 |
| large | 100 | 0.178 / 0.378 | 0.105 / 0.128 | 0.032 / 0.088 | 0.095 / 0.112 | 0.025 / 0.085 |
| large | 300 | 0.208 / 0.341 | 0.069 / 0.103 | 0.027 / 0.098 | 0.096 / 0.124 | 0.023 / 0.097 |

### Split R-hat over the 3 replicate chains

| tree | nChar | config | theta | s | rate | lambda | rootBin | clock | logPost |
|---|---|---|---|---|---|---|---|---|---|
| small | 30 | new | 1.07 | 1.02 | 1.0 | 1.0 | 1.02 | 1.0 | 1.08 |
| small | 30 | old | 1.05 | 1.02 | 1.01 | 1.0 | 1.01 | 1.01 | 1.06 |
| small | 100 | new | 1.02 | 1.01 | 1.01 | 1.0 | 1.01 | 1.01 | 1.03 |
| small | 100 | old | 1.48 | 1.05 | 1.38 | 1.07 | 1.1 | 1.5 | 1.2 |
| small | 300 | new | 1.05 | 1.28 | 1.01 | 1.0 | 1.02 | 1.01 | 1.04 |
| small | 300 | old | 1.17 | 1.22 | 2.3 | 1.1 | 1.32 | 2.5 | 1.11 |
| large | 30 | new | 1.22 | 1.13 | 1.03 | 1.03 | 1.06 | 1.04 | 1.05 |
| large | 30 | old | 2.63 | 1.06 | 1.07 | 1.03 | 1.06 | 1.11 | 1.13 |
| large | 100 | new | 1.05 | 1.13 | 1.05 | 1.02 | 1.15 | 1.05 | 1.1 |
| large | 100 | old | 1.32 | 1.63 | 1.48 | 1.03 | 1.36 | 1.65 | 1.3 |
| large | 300 | new | 1.19 | 1.21 | 1.03 | 1.02 | 1.01 | 1.03 | 1.03 |
| large | 300 | old | 2.84 | 2.49 | 2.88 | 1.11 | 1.09 | 3.32 | 2.2 |

Reading these together:

- **Small tree.** The new schedule is better for every quantity. The gains for `c` and
  the tree clock are 5-44x, and for theta 2.5-4.2x.
- **Large tree.** `old` chains do not converge in `c`, the clock and theta. At 300
  characters their R-hat is 2.9 for `c`, 3.3 for the clock and 2.8 for theta. They also
  under-disperse: at 100-300 characters the posterior sd of log c is 0.03 against
  0.09-0.10, and that of
  log theta 0.18-0.21 against 0.34-0.41. Their ESS for these quantities is therefore
  overstated, and the gains below 1 for `lambda` and `s` (0.25-0.42) compare against chains
  that have not converged. The sd of log lambda is 0.10 under `old` against 0.11-0.12 under
  `new`.
- **Large tree, new schedule.** theta and s are still slow (R-hat 1.05-1.22 for 3 x 5000
  iterations), as the #19 study found. All large-tree theta and s values are
  **provisional**.

## Deviations and limits

- The NCMC theta move (optional in #19) is not implemented. The engine has no incremental
  character pruning, and in the prototype NCMC did not beat the schedule per CPU hour;
  it is left for #26.
- The ESS comparison uses 3 replicates per cell, on two fixed trees and one model cell.
  Gains are **provisional**, with roughly +/-30-50% error where the ESS per chain is below 20.
- The short SBC covers the two gamma = 1 cells only, with 50 replicates each (task: "e.g.
  50 replicates per cell at gamma = 1").
- `charLL` is not a column of `TempoMCMC()` samples, so section 3 reports `logPost` where
  `dev/mixing/` reported the character log-likelihood.
