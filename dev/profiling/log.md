# Profiling log

One entry per round. Findings that need someone else's attention are GitHub
issues labelled `profiling` (the issue number is the finding's id); fixes made in
the round itself are recorded here and in the fixing PR only.

## Round 1 — 2026-10-16 — `TempoMCMC()` on the rung-6 SBC workload

**Target.** The full `TempoMCMC()` iteration, because the queued compute (#21:
full-size rung-6 SBC on Hamilton, with `joint` and tree moves) is a long run of
exactly this loop. Drivers: `drivers/fixture-rung6.R` (one rung-6 replicate:
25 tips, 5 fossils, 3 sampled ancestors, 30 binary characters, K = 16,
z in [-3, 3], nSlice = 16, six free parameters), `drivers/tempo-mcmc.R` (the
chain; move sets `all` = #21 design, `rung6` = current `data-raw/rung6-sbc.R`),
`drivers/tempo-components.R` (each building block timed alone; runs on Hamilton
too), `drivers/egrid-accuracy.R` (E-grid error against a tight reference and
lsode).

**Environment notes.**
- The profile skill's `build-symboled-lib.ps1` installed with
  `--no-byte-compile`; the R JIT then compiled package closures during the run
  and appeared as ~14 % `cmpfun` time. Fixed in the skill (byte-compiles now);
  all numbers below are from byte-compiled builds.
- The installed MkPrime (May 2026) predates `inst/include`; local builds use
  headers from `agent-issues/MkPrime@main` staged in `dev/profiling/.deps-lib`
  (git-ignored). CI installs MkPrime from `Remotes:`, so nothing changes there.
- Windows `Rprof` ticks at ~16 ms whatever `interval` says: shares are valid,
  absolute times are not. `.PathSummary()` / `.TreeStats()` run once per stored
  sample, so `thin = 1` in the driver inflates them ~50x relative to the SBC
  (`thin` ~48); they are not targets.

**Where the time went (before).** Windows, `all`, Rprof of the chain alone:
C++ ~75 % of the iteration — `.SampleTreeCpp` 27 %, `.EGridCpp` 22 %,
`.MessagesCpp` 7 %, `.LogWeightCpp` 5 %, `.ExtinctionCpp` 3 %; R `eigen()` in
`.AddSlices()` 4 %. By move: `joint` 56 % (six parameters x [E grid, slices,
messages, FFBS draw]), parameter sweep 17 %, tree moves 10 %, path moves 10 %.
VTune (EvoTempo.dll): `.SampleTreeCpp` is ~45 % MinGW's statically linked x87
`exp` (`_expl_internal`) inside `Bridge`'s Brent root-finder; `.EGridCpp` is
Radau steps (`_divdc3` complex division, banded LU, right-hand sides);
`.Call` marshalling (`strcmp`, `Rf_translateChar`, `malloc` from rebuilding
`FrozenSet` / `Msg` / `TreeView` per call) ~5-6 %.

**Linux cross-check (Hamilton, gcc 14.2, R 4.5.1, `test` partition).** 85 ms
per `all` iteration, 30 ms per `rung6` iteration. `.SampleTree` is 3.3x faster
than on Windows against ~1.7-1.9x for everything else, confirming the x87 `exp`
cost as a MinGW artefact; on Linux the E grid is the top target (~30 % of `all`,
~42 % of `rung6`), then path draws (~18 %) and slices (~12 %).

**Fixed this round** (PR from branch `profile-tempo-mcmc`; no issues filed):

1. *`rate` proposals rebuilt the frozen model.* `rate` does not enter E or the
   generator, yet each `rate` proposal in the parameter sweep and in `joint`
   recomputed `.Precompute()` (and in `joint` also `.AddSlices()` and
   `.Messages()`). Sets now carry a key (E parameters and node ages,
   `.SetKey()`); equal keys reuse the set. Bit-identical chains (checked on 40
   `all` iterations and by a test that forces recomputation).
   Windows `all`: 206 -> 187 ms/iter (-9.5 %).
2. *E grid stopped at every output age.* `.EGridCpp()` integrated to each of
   the ~490 grid ages (504 Radau steps; 66 natural steps cover the same span).
   It now steps freely and interpolates with Radau IIA's continuous output (the
   collocation polynomial, as RADAU5's CONTR5). That output has stage order 3,
   so the solver runs at a 100x tighter tolerance; at rtol = 1e-10 the grid's
   error against a 1e-13 reference is then 2-4e-10 in E (forced steps:
   0.9-18e-10; lsode, the R engine: 3.4-14e-10) across the fixture, the
   default grid, theta = 2, and lambda = 2 over 30 time units
   (`drivers/egrid-accuracy.R`). Steps 504 -> 134 (fixture). Micro-bench
   (bench::mark, Windows): 3.49 -> 1.20 ms (K = 16), 5.4 -> 2.25 ms (K = 26).
3. *`eigen()` per slice.* `.AddSlices()` now computes E at the slice midpoints
   and the tridiagonal eigendecompositions in C++ (`src/slices.cpp`, LAPACK
   `dstev`) when the engine is `cpp`; the R engine keeps `eigen()`.
   Eigenvalues agree to 1e-12, rebuilt matrices to 1e-12 (test-Engine.R).
   Linux: 1469 -> 477 us per call.
4. *Brent's end points re-evaluated.* `Bridge()` passed F(0) and F(remain) to
   the root-finder by evaluating J again (2K transcendentals each), although
   J(0) = 0 and J(remain) is the vector that gave the pick weight. Now passed
   directly: bit-identical draws (checked on a 40-iteration chain).
   `.SampleTree()` 7.1 -> 5.25 ms (Windows), 2.19 -> 1.91 ms (Linux).

Cumulative, interleaved A/B x3:

| Platform | `all` ms/iter | `rung6` ms/iter |
|----------|---------------|-----------------|
| Windows (MinGW) | 129 -> 89 (-31 %) | 41.6 -> 29 (-30 %) |
| Linux (Hamilton, gcc) | 83.8 -> 54.3 (-35 %) | 29.0 -> 18.8 (-35 %) |

No OS trade-off: every change is faster on both. Against the old build, a
60-iteration `all` chain makes identical accept/reject decisions; `logPost`
differs by <= 3e-12 relative. Full test suite passes (4236 expectations);
changed files at 100 % line coverage except pre-existing gaps (Radau
Newton-failure and step-rejection branches; R `.Bridge()` infeasible branch).
`data-raw/rung4b-likelihood-check.R` rerun on both grids: importance-sampled
and exact likelihoods move by <= 5e-10, z-scores by < 1e-7; result files
updated. Not rerun (multi-hour cluster jobs): rungs 3, 4a, 5, 6 and the
default-grid validation; the changes move log densities by ~1e-10, far below
their Monte Carlo error.

**Where the time goes now (Linux, `all`, 54 ms).** `joint` ~44 % (per
parameter: E grid 0.77, slices 0.48, messages 0.73, tree draw 1.9 ms), tree
moves ~17 % (`nodeAge` 2.6, `slide` 3.8, `saJump` 2.6 ms, pure R), parameter
sweep ~13 %, path moves ~13 %. The tree draw (focus area 2) and the R tree
moves (focus area 7) are the next targets.

**Measured but not changed.**
- *Bridge root-finding* (focus area 2): Windows-heavy because of x87 `exp`.
  Cheaper evaluations by hoisting `exp(lambda remain)` out of J are unsafe
  (0 x inf at the default grid). About 5-6 Brent evaluations remain per jump;
  a bracketed Newton step using the integrand (a by-product of the same
  transcendentals) should cut them to ~3, at the price of jump times that
  differ from the R engine's `uniroot()` at ~1e-12 (within test tolerance).
- *Marshalling* (focus area 6): belongs to #26 item 1 (state held in C++).
- *Survival solve* (focus area 8): a second E solve with psi = 0 when
  `condition = "extant"`; cannot reuse the E grid.

last_focus: 1
