# `wfnpr` validation results (#110, move 4)

Derivation: `dev/spec/weighted-regraft.md`.

## Topology and clade frequencies at fixed ages and parameters

`data-raw/wfnpr-validation.R`, move code as of 0ad095c. With gamma = 0 the topology has
probability proportional to `TreeLikelihood()` times the Mk likelihood, over the topologies that
fixed-height regrafts reach. The statistics are every clade's probability and every topology of
probability above 0.02. Each case ran once with seed 0, 32 chains of 20 000 iterations per move
(the first tenth discarded); z uses the between-chain standard error. `data-raw/wfnpr-pool.R`
pools every chain of a case and writes `pooled.csv` (against the exact values) and
`pooled-fnpr.csv` (`wfnpr` against `fnpr` chains).

Trees:

- **sa5**: five tips, two fossils and a sampled ancestor; 18 topologies, 8 of probability 0 (no
  extant tip below a root daughter); never more than two targets.
- **extant5**: five extant tips; 180 topologies; up to three targets, so that a cap of 3
  (`WFNPR_MAX=3`) draws subsets of two candidates and a cap of 2 subsets of one. The default cap,
  8, never truncates on either tree.

Against the exact values (`pooled.csv`; `pMax` is the probability that the largest of as many
independent |t| values, on 31 degrees of freedom, exceeds the one observed):

| Tree | Case | Cap | Statistics | `wfnpr` max abs z | rms z | pMax | `fnpr` max abs z |
|---|---|---|---|---|---|---|---|
| sa5 | prior only | 8 | 17 | 1.51 | 0.78 | 0.92 | 1.98 |
| sa5 | 4 characters, gamma = 0 | 8 | 17 | 1.33 | 0.79 | 0.97 | 2.15 |
| extant5 | prior only | 3 | 25 | 2.78 | 1.20 | 0.21 | 1.83 |
| extant5 | 5 characters, gamma = 0 | 8 | 39 | 2.11 | 1.06 | 0.82 | 2.25 |
| extant5 | 5 characters, gamma = 0 | 3 | 39 | 2.10 | 0.86 | 0.83 | 2.25 |
| extant5 | 5 characters, gamma = 0 | 2 | 39 | 2.16 | 1.07 | 0.78 | 2.25 |

The three extant5 character cases share one set of `fnpr` control chains (same seeds).

With gamma = 1 the characters depend on the tempo paths, and the exact probabilities come from
pruning on the product of the tempo bins and the character patterns (`ProductLogLik()` in the
script, checked there against the factorized likelihood at gamma = 0; `<tree>-<case>-exact.csv`).
Against them (`pooled.csv`):

| Tree | Case | Cap | Chains per move | Statistics | `wfnpr` max abs z | rms z | pMax | `fnpr` max abs z | rms z | pMax |
|---|---|---|---|---|---|---|---|---|---|---|
| extant5 | 5 characters, gamma = 1 | 3 | 32 | 39 | 2.80 | 1.12 | 0.29 | 2.92 | 1.55 | 0.23 |
| sa5 | 4 characters, gamma = 1 | 8 | 320 | 17 | 2.27 | 1.44 | 0.34 | 2.82 | 1.39 | 0.084 |

The sa5 case has three runs, because the first two were not clean when `wfnpr` was compared with
`fnpr`, the only reference then available:

| Run | Chains per move | topology 002, z `wfnpr` - `fnpr` | clade a+c | largest of 17 |
|---|---|---|---|---|
| seed 0 (local, 0ad095c) | 32 | 3.32 | 2.48 | 3.32 |
| seed 1000 (local, 0ad095c) | 32 | 0.74 | 1.26 | 1.49 |
| seed 2000 (Hamilton Slurm 19572313, 3d4eda7) | 256 | 0.36 | 0.69 | 3.26 (clade c+d) |

The third run was a test fixed before it was submitted (in the header of its job script): those
two statistics, 256 new chains, verdict from the new chains alone, |z| >= 3 for either as evidence
of a bias and |z| < 2 for both as noise. By that rule the first run's outlier was noise. The same
run then put a different statistic, clade c+d, at -3.26 between the moves (3.12 over all 320
chains), which is why the exact values were added. Against them the two moves sit on opposite
sides of the exact value for that clade, 0.4724: `wfnpr` at z = -2.27 and `fnpr` at +2.16 over 320
chains each. Neither move exceeds 3 on any statistic, and `wfnpr` is no further from the exact
values than the `fnpr` control. The rms z is about 1.4 for both moves, against 0.8 to 1.2 in the
gamma = 0 cases; the statistics are correlated, so this is not a test, and whatever it reflects it
is common to both moves.

For the gamma = 0 cases, `pooled-fnpr.csv` also gives `wfnpr` against `fnpr`: max abs z 1.81 to
2.70.

Acceptance, `wfnpr` / `fnpr`: sa5 0.336 / 0.300 (prior), 0.267 / 0.233 (characters), 0.255 / 0.216
(gamma = 1); extant5 with characters 0.280 / 0.233 at the default cap, 0.271 at cap 3 and 0.233 at
cap 2 (one candidate: the `fnpr` kernel).

### First runs

Before the script wrote per-chain output it was run on sa5 with four random characters and
rate 0.3 (`wfnpr-validation*.csv`; topologies only, probabilities 0.077 to 0.133):

| File | Case | Chains x iterations | max abs zExact | max abs zFnpr |
|---|---|---|---|---|
| `wfnpr-validation.csv` | prior only | 8 x 10 000 | 1.58 | 1.33 |
| `wfnpr-validation.csv` | 4 characters, gamma = 0 | 8 x 10 000 | 3.06 | 1.52 |
| `wfnpr-validation-characters-50000.csv` | the same, other seeds | 16 x 20 000 | 1.71 | 1.90 |
| `wfnpr-validation.csv` | 4 characters, gamma = 1 | 8 x 10 000 | - | 1.86 |
| `wfnpr-validation-characters-70000-max2.csv` | 4 characters, gamma = 0, cap 2 | 12 x 15 000 | 3.53 | 4.73 |
| `wfnpr-validation-characters-90000-max2.csv` | the same, other seeds | 32 x 20 000 | 2.69 | 2.97 |

Two of these exceeded 3 and were each followed by a larger run with new seeds. Pooling the runs of
a case (means and standard errors of the two files, runs weighted by chains x iterations, since
there is no per-chain output): max abs z 1.07 at the default cap (24 chains) and 2.47 at cap 2
(44 chains), where the `fnpr` control gives 2.32 and 1.83. With 8 chains z has 7 degrees of
freedom: |z| = 3.06 has probability 0.018 for one topology and 0.17 for the largest of ten.

## Rung 5 at fixed parameters, prior only (`data-raw/rung5-fixed-parameters.R`)

`EVOTEMPO_RUNG5_TOPOLOGY=wfnpr`: `wfnpr` is the only topology move; 8 chains of 4000 iterations
(the script's design) against 2000 trees from rejection simulation. Hamilton Slurm 19567926 (coarse
grid) and 19567927 (K = 26), source 6c0bc9c;
`data-raw/rung5-results/rung5-fixed-parameters[-K26]-wfnpr.rds`.

| statistic | z, coarse grid | z, K = 26 |
|---|---|---|
| treeLength | 0.23 | 0.48 |
| meanNodeAge | 0.36 | -0.30 |
| oldestFossil | -0.41 | 0.63 |
| youngestFossil | -0.60 | 1.19 |
| clock | -1.04 | 1.12 |
| inverseClock | 0.81 | 0.15 |
| nAncestor | -0.30 | -0.89 |
| cherries | 0.44 | 0.83 |

The same script with `fnpr` (Slurm 19567924, 19567925; not committed here, they belong with
#118) gave max |z| 2.40 (cherries, coarse grid) and 0.77 (K = 26).

## Penguin b0g0 pilot

Production schedule plus `CET_FNPR=1 CET_WFNPR=1` (6-7 of each per iteration), 4 replicates from
the dispersed start trees, 1000 burn-in + 8000 iterations; Hamilton Slurm 19567923, source 6c0bc9c,
output `/nobackup/pjjg18/evotempo-wfnpr/runs/`. Acceptance over the seven segments without burn-in
(7000 iterations per replicate), range over replicates:

| move | acceptance |
|---|---|
| `slide` | 0.039-0.044 |
| `fnpr` | 0.105-0.110 |
| `wfnpr` | 0.172-0.175 |
| `saJump` | 0.297-0.307 |
| `nodeAge` | 0.502-0.508 |

0.456-0.463 s per iteration. At one state of replicate 1 (400 proposals each, Slurm 19568085): a
`slide` proposal takes 5.9 ms, `fnpr` 5.1 ms and `wfnpr` 14.9 ms; one character likelihood takes
3.6 ms and `.LogAugmented()` 0.14 ms; `wfnpr` scores 3.25 candidates per proposal on average (17% of
proposals have no target, 22% reach the cap of 7). These timings predate the change that skips
the likelihood of a candidate of zero density.

Measured again at the same state with the build at 3d4eda7 (Slurm 19573107, a faster node:
character likelihood 2.65 ms): `slide` 4.5 ms, `fnpr` 3.9 ms, `wfnpr` 11.4 ms per proposal, 3.04
candidates scored per proposal. The checkpoints written at 6c0bc9c are read by this build and give
the same character likelihood to the six decimals printed.

Infeasible draws (the numerical guard that rejects a move when a bridge or a weight underflows):
none in 8000 proposals of each of `slide`, `fnpr` and `wfnpr`, 2000 from the last state of each pilot
replicate with only that move running (so parameters and the other paths stay fixed).
