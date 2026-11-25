# `wfnpr` validation results (#110, move 4)

Derivation: `dev/spec/weighted-regraft.md`. All runs at commit 6c0bc9c.

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

The one |z| above 3 (topology t02, 8 chains, 7 degrees of freedom) did not recur in the larger
replicate, where the same topology has z = 0.51. No chain visited a topology of probability 0.

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
