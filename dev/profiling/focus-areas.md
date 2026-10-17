# Profiling focus areas

Ranked by (cost per call) x (calls per `TempoMCMC()` iteration) on the rung-6 SBC
workload (`drivers/fixture-rung6.R`), the next large compute job (#21).
Statuses: `NEW`, `PROFILED`, `OPTIMISED`, `AT-LIMIT`, `SKIPPED`.
Round records are in `log.md`; current timings in `baselines.md`.

| # | Area | Files | Why hot | Baseline cost | Last profiled | Status |
|---|------|-------|---------|---------------|---------------|--------|
| 1 | E grid | `src/tree_likelihood.cpp` (`EGridCpp`), `src/radau5.h` | One stiff solve per `.Precompute()`, i.e. per parameter proposal (10 per iteration with `joint`) | 2.13 -> 0.77 ms (Linux, K = 16) | 2026-10-16 | OPTIMISED |
| 2 | Path draws (FFBS bridges) | `src/paths.cpp` (`Bridge`, `SampleBranch`) | Every path move and every `joint` proposal; Brent root-finding per jump, 2K transcendentals per evaluation | 2.19 -> 1.91 ms per tree draw (Linux); Brent end points fixed | 2026-10-16 | PROFILED |
| 3 | Slice eigendecompositions | `R/Paths.R` (`.AddSlices`), `src/slices.cpp` | Per new parameter set: `eigen()` on nSlice symmetric tridiagonal K x K matrices | 1.47 -> 0.48 ms (Linux) | 2026-10-16 | OPTIMISED |
| 4 | Parameter-move recomputation | `R/MCMC.R` (`Target`, `.JointSweep`) | `rate` proposals rebuilt the frozen model although `rate` does not enter it | 2 x (E grid + slices + messages) per iteration | 2026-10-16 | OPTIMISED |
| 5 | Frozen messages | `src/paths.cpp` (`MessagesCpp`) | Per new parameter set and per accepted tree move; ~7 % of an iteration | 1.8 ms (Win), 0.8 ms (Linux) | 2026-10-16 | PROFILED |
| 6 | `.Call` marshalling | `src/frozen.h`, `src/tree.h`, `src/paths.cpp` | `FrozenSet`, `Msg` and `TreeView` are rebuilt from R lists on every call (name lookups, copies of the E grid and every transfer matrix): `strcmp`, `Rf_translateChar`, `malloc` ~5-6 % | -- | 2026-10-16 | PROFILED (#26 item 1) |
| 7 | Tree moves | `R/TreeMoves.R` | Pure-R transfer rows/columns, path densities and bin weights per tree move | 2.6-3.8 ms per `nodeAge` / `slide` / `saJump` (Linux), ~17 % of an `all` iteration | -- | NEW (#26 item 4) |
| 8 | Survival conditioning solve | `R/Paths.R` (`.Precompute`), `ExtinctionCpp` | Second E solve (psi = 0) per `.Precompute()` when `condition = "extant"` and psi > 0 | ~3 % of an iteration | 2026-10-16 | PROFILED |
| 9 | Character likelihood | `src/characters.cpp`, `R/Characters.R` | Whole-tree Mk pruning per path move | 0.2 ms per evaluation | 2026-10-16 | PROFILED (#26 item 3) |
| 10 | Marginal likelihood | `src/tree_likelihood.cpp` (`PruneTreeCpp`) | `MarginalMCMC()`, rung 3-4 oracles; joint (E, D) Radau per branch | -- | -- | NEW |

SKIPPED unless asked: `TempoTree()`, `DCETModel()`, forward simulators (validation scripts only),
`.PathSummary()` / `.TreeStats()` (once per stored sample; `thin` is ~48 in the SBC).
