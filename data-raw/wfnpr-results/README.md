# `wfnpr` validation results (#110, move 4)

Derivation: `dev/spec/weighted-regraft.md`. All runs use the move as committed in 185b8a5
(Hamilton source 6c0bc9c).

## Topology frequencies at fixed ages and parameters (`data-raw/wfnpr-validation.R`)

Five-tip tree with a sampled ancestor, fossil tips and extant tips; K = 7 bins. Fixed-height
regrafts reach 18 topologies, 8 of which have probability 0 (the conditioning on extant
descendants of both root daughters). `zExact` compares the mean frequency over chains with the
exact probability, using the between-chain standard error; `zFnpr` compares with `fnpr` chains.

| File | Case | Chains x iterations | max abs zExact | max abs zFnpr | acceptance `wfnpr` / `fnpr` |
|---|---|---|---|---|---|
| `wfnpr-validation.csv` | prior only (gamma = 0, no characters) | 8 x 10 000 | 1.58 | 1.33 | 0.336 / 0.299 |
| `wfnpr-validation.csv` | 4 characters, gamma = 0 | 8 x 10 000 | 3.06 | 1.52 | 0.332 / 0.291 |
| `wfnpr-validation-characters-50000.csv` | the same, other seeds | 16 x 20 000 | 1.71 | 1.90 | 0.332 / 0.291 |
| `wfnpr-validation.csv` | 4 characters, gamma = 1 (no exact oracle) | 8 x 10 000 | - | 1.86 | 0.307 / 0.267 |
| `wfnpr-validation-characters-70000-max2.csv` | 4 characters, gamma = 0, `WFNPR_MAX=2` | 12 x 15 000 | 3.53 | 4.73 | 0.292 / 0.290 |
| `wfnpr-validation-characters-90000-max2.csv` | the same, other seeds | 32 x 20 000 | 2.69 | 2.97 | 0.291 / 0.292 |

The one |z| above 3 (topology t02, 8 chains, 7 degrees of freedom) did not recur in the larger
replicate, where the same topology has z = 0.51. No chain visited a topology of probability 0.

`WFNPR_MAX=2` sets the candidate cap to 2, so that each proposal scores one target drawn uniformly
from `fnpr`'s: on this tree no proposal has more than two targets, and the default cap of 8 never
subsamples. In the first run topology t09 has z = 3.53 for `wfnpr` and, against the same exact
value, z = -3.25 for the `fnpr` control chains (hence 4.73 between them); t15 has z = -3.31 for
`fnpr`. In the 32-chain rerun t09 has z = -1.17 and the largest is t08 at 2.69 (0.12 in the first
run). The deviations do not recur in sign or place, and the validated `fnpr` control scatters as
widely, so they are read as noise in standard errors estimated from few chains. They are not as
clean as the other rows.

## Rung 5 at fixed parameters, prior only (`data-raw/rung5-fixed-parameters.R`)

`EVOTEMPO_RUNG5_TOPOLOGY=wfnpr`: `wfnpr` is the only topology move; 8 chains of 4000 iterations
against 2000 trees from rejection simulation. Hamilton Slurm 19567926 (coarse grid) and 19567927
(K = 26); `data-raw/rung5-results/rung5-fixed-parameters[-K26]-wfnpr.rds`.

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
the dispersed start trees, 1000 burn-in + 8000 iterations; Hamilton Slurm 19567923, output
`/nobackup/pjjg18/evotempo-wfnpr/runs/`. Acceptance over the seven segments without burn-in (7000
iterations per replicate), range over replicates:

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
proposals have no target, 22% reach the cap of 7).
