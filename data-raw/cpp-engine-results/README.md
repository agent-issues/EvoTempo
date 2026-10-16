# C++ engine vs R reference

Output of `data-raw/cpp-engine-agreement.R` (commit f267821), run on one core of
a Hamilton `shared` node (R 4.5.1, g++ 11.2, `-O2`).
Grids: K16 = 16 bins on z in [-3, 3]; K26 = 26 bins on z in [-5, 5] (the default).

Fixtures:

- rung 3 trees (`data-raw/rung3-validation.R`): `sampledAncestor` (5 tips, a sampled
  ancestor and a stem; survival conditioning, origin start, removal 0.4),
  `ultrametric` (12 extant tips; crown conditioning, joint), `fossilTips` (10 tips;
  no conditioning, removal 1), each at the 20 rung-3 parameter points with psi = 0.15;
- `rung4`: the rung-4 data tree at its true parameters and 19 draws from `TempoPriors()`;
- `rung4b`: the same tree and parameters with the rung-4b binary character
  (product-space likelihood, `TreeLikelihood(character = )`).

## Agreement of the log-likelihood (`TreeLikelihood()`)

"Tight R" is the R engine at `rtol = 1e-13`, `atol = 1e-18` (first 3 points of each
fixture without a character): the C++ solver is about ten times nearer it than lsode is
at the default tolerances.

| Fixture | Grid | n | max rel. diff C++ vs R | C++ vs tight R | R vs tight R |
|---|---|---|---|---|---|
| fossilTips | K16 | 20 | 8.2e-11 | 5.3e-12 | 5.7e-11 |
| rung4 | K16 | 20 | 1.1e-10 | 1.2e-11 | 1.1e-10 |
| rung4b | K16 | 20 | 9.1e-11 | NA | NA |
| sampledAncestor | K16 | 20 | 8.1e-11 | 5.0e-12 | 6.2e-11 |
| ultrametric | K16 | 20 | 3.1e-10 | 4.2e-11 | 2.6e-10 |
| fossilTips | K26 | 20 | 7.1e-11 | 5.7e-12 | 5.6e-11 |
| rung4 | K26 | 20 | 1.8e-10 | 1.8e-11 | 2.0e-10 |
| rung4b | K26 | 20 | 1.5e-10 | NA | NA |
| sampledAncestor | K26 | 20 | 6.4e-11 | 5.9e-12 | 5.0e-11 |
| ultrametric | K26 | 20 | 1.9e-10 | 3.5e-11 | 1.6e-10 |

## Augmented density, weights and characters (rung-4 tree, 5 points x 10 FFBS draws)

Draws use the same seed in both engines.  `logWeight` values are near zero
(|w| < 0.005), so their relative differences are large although the absolute
differences are below 7e-10; they inherit the R engine's E-grid error.

| Quantity | Grid | n | max rel. diff |
|---|---|---|---|
| charLogLik | K16 | 50 | 3.0e-13 |
| charLogLik | K26 | 50 | 3.0e-13 |
| ffbsDraw | K16 | 5 | bins identical in 5/5 draws |
| ffbsDraw | K26 | 5 | bins identical in 5/5 draws |
| frozenLogLik | K16 | 5 | 7.0e-12 |
| frozenLogLik | K26 | 5 | 8.5e-12 |
| logAugmented | K16 | 50 | 3.6e-11 |
| logAugmented | K26 | 50 | 6.0e-11 |
| logWeight | K16 | 50 | 3.1e-07 |
| logWeight | K26 | 50 | 3.2e-07 |

## Run times (median per call)

`frozenLogLik` times `.Precompute()`, `.AddSlices()` and `.Messages()`, whose
eigendecompositions stay in R.

| Fixture | Quantity | Grid | R (ms) | C++ (ms) | Speed-up |
|---|---|---|---|---|---|
| rung4 | charLogLik | K16 | 1.12 | 0.095 | 11.8x |
| rung4 | charLogLik | K26 | 1.13 | 0.095 | 11.9x |
| rung4 | ffbsDraw | K16 | 36.2 | 0.5 | 72.4x |
| rung4 | ffbsDraw | K26 | 46.3 | 0.9 | 51.4x |
| rung4 | frozenLogLik | K16 | 29 | 8 | 3.6x |
| rung4 | frozenLogLik | K26 | 37 | 11 | 3.4x |
| rung4 | logAugmented | K16 | 0.685 | 0.06 | 11.4x |
| rung4 | logAugmented | K26 | 0.73 | 0.07 | 10.4x |
| fossilTips | logLik | K16 | 202 | 14.5 | 14.0x |
| rung4 | logLik | K16 | 107 | 8.5 | 12.6x |
| rung4b | logLik | K16 | 185 | 24 | 7.7x |
| sampledAncestor | logLik | K16 | 105 | 8 | 13.1x |
| ultrametric | logLik | K16 | 192 | 14.5 | 13.2x |
| fossilTips | logLik | K26 | 444 | 34.5 | 12.9x |
| rung4 | logLik | K26 | 186 | 20 | 9.3x |
| rung4b | logLik | K26 | 411 | 51.5 | 8.0x |
| sampledAncestor | logLik | K26 | 233 | 17.5 | 13.3x |
| ultrametric | logLik | K26 | 481 | 38.5 | 12.5x |

## TempoMCMC, rung-4 tree with the rung-4b character

20 iterations of `tree`, `node`, `tip`, `branch` and `joint` moves, five free
parameters plus `rate`; both engines finish at the same `logPost` (difference
< 4e-10).

| Grid | R (s/iteration) | C++ (s/iteration) | Speed-up |
|---|---|---|---|
| K16 | 0.697 | 0.084 | 8.3x |
| K26 | 0.807 | 0.132 | 6.1x |
