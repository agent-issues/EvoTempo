# `sfnpr` validation results (#110 move 3)

Plan and derivation: `dev/plans/spliced-fnpr.md`. Every run below is at commit cbe1d7a, on Hamilton, on 26 Nov 2026.

## Exact topology oracle at fixed ages and parameters (`data-raw/sfnpr-validation.R`)

- **Exact values:** with one character, `TreeLikelihood(character = )` is the exact likelihood of the time tree and the character, with tempo integrated out by product-space pruning. P(time tree) is proportional to it over every time tree that fixed-height regrafts reach.
- **Chains:** 32 per move, 20 000 iterations each (the first tenth discarded). The schedule is `node`, `tip`, `branch` and one topology move: `sfnpr`, or `fnpr` as the control.
- **Grid and parameters:** K = 7 on z in [-2, 2], with lambda 0.5, mu 0.2, psi 0.1, theta 0.3, s 1, rate 0.4 and rho 0.7.
- **Statistics:** each time tree of probability above 0.02 and each clade. z uses the between-chain standard error.

| Tree | Case | Time trees (P > 0) | Statistics | Slurm | `sfnpr` max \|z\| | rms z | `fnpr` max \|z\| | rms z |
|---|---|---|---|---|---|---|---|---|
| sa5 (2 fossils, a sampled ancestor) | prior only | 18 (10) | 17 | 19573170 | 2.62 | 1.17 | 2.25 | 1.20 |
| sa5 | 1 character, gamma = 1 | 18 (10) | 17 | 19573171 | 1.30 | 0.77 | 1.96 | 1.14 |
| extant5 | prior only | 180 (180) | 25 | 19573172 | 2.14 | 0.88 | 2.24 | 0.94 |
| extant5 | 1 character, gamma = 1 | 180 (180) | 25 | 19573173 | 2.52 | 0.92 | 2.10 | 1.10 |

For 17-25 statistics (t on 31 degrees of freedom), a max \|z\| near 2.5 and an rms z near 1 are what an exact sampler gives. Per-statistic values are in `<tree>-<case>-g1-0.csv`.

## Rung 5 at fixed parameters, prior only (`data-raw/rung5-fixed-parameters.R`)

`EVOTEMPO_RUNG5_TOPOLOGY=sfnpr`: `sfnpr` is the only topology move, with 8 chains of 4000 iterations against 2000 trees from rejection simulation. Results: `data-raw/rung5-results/rung5-fixed-parameters[-K26]-sfnpr.rds`.

| statistic | z, coarse grid (Slurm 19573174) | z, K = 26 (Slurm 19573175) |
|---|---|---|
| treeLength | 0.53 | 0.45 |
| meanNodeAge | 0.48 | -0.15 |
| oldestFossil | 0.19 | -0.37 |
| youngestFossil | -0.92 | 0.54 |
| clock | -1.21 | 1.22 |
| inverseClock | 1.29 | -0.17 |
| nAncestor | -0.58 | -0.50 |
| cherries | 0.84 | 0.35 |
