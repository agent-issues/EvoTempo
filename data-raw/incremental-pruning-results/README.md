# Incremental character pruning: benchmark (#26, item 3)

Benchmark of the persistent character cache of `.CharacterLikelihood()`
(`dev/plans/incremental-pruning.md`) against the package before the change, i.e. the head of #43,
`origin-root-age-18` at b99e810. Everything ran on Hamilton with the C++ engine.

Files:

- `bench-chains.csv`, `bench-gains.csv`: per chain and per cell. Each Slurm task ran the same
  seeded chain under `base` (b99e810) and then under `inc` (this branch, package code of a1bd997),
  one after the other on the same node;
- `conv-chains.csv`, `conv-rhat.csv`: four longer large-tree chains per character count, under
  `inc`;
- the summary script is `data-raw/incremental-pruning-bench.R`; the chains come from
  `data-raw/mixing-schedule.R chain <tree> <nChar> new <rep> <nIter> <burnin>`.

The fixtures are those of `dev/mixing/` and `data-raw/mixing-results/`: two fixed true trees
(small: 17 extant tips, 69 path branches; large: 75 extant tips, 272 path branches), cell
beta = gamma = 1, K = 16, z in [-3, 3], 30 / 100 / 300 nested characters, all six parameters
free, start at the truth, default schedule `TempoMoves()` (`node`, `tip` and `branch` at
0.25 x the nodes, tips and branches, one `shift`). The `.rds` chains are not committed.

## 1. Identical chains, lower cost

All 18 pairs of chains (2 trees x 3 character counts x 3 replicates; 2000 burn-in plus 5000
sampling iterations) gave **identical** samples and acceptance rates under `base` and `inc`.
The incremental likelihood is bitwise equal to full pruning, and random-number use is unchanged.
So the ESS is the same, and the gain in ESS per CPU hour is the ratio of CPU times, for
every quantity.

| tree | nChar | s/iter base | s/iter inc | gain | theta ESS/h base | theta ESS/h inc | s ESS/h base | s ESS/h inc | c ESS/h base | c ESS/h inc |
|---|---|---|---|---|---|---|---|---|---|---|
| small | 30 | 0.037 | 0.032 | 1.18 | 907 | 1068 | 1384 | 1630 | 5362 | 6314 |
| small | 100 | 0.043 | 0.033 | 1.31 | 395 | 518 | 433 | 567 | 3031 | 3973 |
| small | 300 | 0.049 | 0.034 | 1.45 | 608 | 884 | 227 | 330 | 3478 | 5060 |
| large | 30 | 0.191 | 0.137 | 1.39 | 38 | 53 | 112 | 157 | 273 | 381 |
| large | 100 | 0.205 | 0.127 | 1.61 | 44 | 71 | 128 | 206 | 319 | 514 |
| large | 300 | 0.298 | 0.134 | 2.23 | 39 | 86 | 46 | 102 | 257 | 572 |

- The cost per iteration under `inc` hardly depends on the number of characters any more
  (0.127-0.137 s on the large tree at 30-300 characters, against 0.191-0.298 s before).
  What remains is the path-move machinery: `.SampleClade()`, `.LogWeight()`,
  `.SampleBranch()` and the R loop (items 1 and 4 of #26).
- A profile of 300 iterations on the large tree with 300 characters, before the length memo was
  added, put `MkLikelihood()` at 53% of `base`'s time (42% in `.MkLikelihoodCpp()` itself, the rest in
  converting the character matrix at every call). After the pruning became incremental,
  recomputing every branch length took 21% of the time, and the pruning itself 4%. That is why
  the cache also memoizes branch lengths. Profiles on different nodes are not comparable in
  absolute time; the paired chains above are.

## 2. Large-tree convergence of theta and s

Four chains per character count, 20000 sampling iterations after 2000 of burn-in, 1000
stored draws each, under `inc` (seeds of replicates 11-14). Split R-hat over the four chains
uses the first quarter, the first half and all of each chain's draws:

| nChar | iterations | CPU h / chain | theta | s | c (rate) | lambda | rootBin | clock | logPost |
|---|---|---|---|---|---|---|---|---|---|
| 30 | 5000 | 0.22 | 1.14 | 1.23 | 1.04 | 1.01 | 1.09 | 1.05 | 1.02 |
| 30 | 10000 | 0.43 | 1.09 | 1.14 | 1.02 | 1.00 | 1.06 | 1.02 | 1.01 |
| 30 | 20000 | 0.86 | 1.06 | 1.03 | 1.02 | 1.01 | 1.04 | 1.02 | 1.05 |
| 100 | 5000 | 0.22 | 1.15 | 1.15 | 1.05 | 1.02 | 1.05 | 1.06 | 1.13 |
| 100 | 10000 | 0.44 | 1.11 | 1.05 | 1.02 | 1.01 | 1.04 | 1.02 | 1.09 |
| 100 | 20000 | 0.89 | 1.02 | 1.06 | 1.01 | 1.00 | 1.03 | 1.01 | 1.01 |
| 300 | 5000 | 0.21 | 1.31 | 1.87 | 1.04 | 1.03 | 1.00 | 1.04 | 1.20 |
| 300 | 10000 | 0.43 | 1.22 | 1.44 | 1.03 | 1.01 | 1.02 | 1.03 | 1.10 |
| 300 | 20000 | 0.85 | 1.15 | 1.23 | 1.02 | 1.01 | 1.01 | 1.02 | 1.08 |

- **theta and s are still not converged at 300 characters**, even at 20000 iterations
  (0.85 CPU hours per chain): theta 1.15, s 1.23. At 30 and 100 characters they reach
  1.02-1.06 at 20000 iterations, but not at 5000.
- Per chain, theta ESS is 10-43 and s ESS 6-73 per 20000 iterations
  (`conv-chains.csv`), about 1 per 1000 iterations. The cache makes each iteration cheaper
  (1.4-2.2x on this tree). It does not change how slowly the local path moves refresh the
  paths that pin theta and s (`dev/mixing/README.md`, diagnosis 3). The sampler itself needs
  a better theta/s move.
- For comparison at equal CPU: a `base` chain of 5000 iterations costs 0.27-0.41 CPU hours on
  this tree. In that budget `inc` runs about 7000-11000 iterations. The R-hat values at
  10000 iterations are then the relevant ones: theta 1.09-1.22, s 1.05-1.44.

## 3. NCMC for theta and s, with the cache

The NCMC prototype of `dev/mixing/` (`L25SN`: `L25S` plus NCMC moves on theta and s, 8 stages,
relaxation 0.25 x branches) was rerun on the large tree against `L25S`. Its character
likelihood was the cached closure (`dev/mixing/run.R` with `CharLLFun()` replaced by
`EvoTempo:::.CharacterLikelihood()`), built on this branch. There were 2 replicates per cell:
`L25S` ran 2000 + 5000 iterations and `L25SN` 1000 + 3000. Time is the job's wall time
(single-threaded), scaled to the sampling phase.

| nChar | config | s / iteration | theta ESS / 1000 it | s ESS / 1000 it | theta ESS / h | s ESS / h |
|---|---|---|---|---|---|---|
| 100 | L25S | 0.137 | 4.4 | 3.0 | 114 | 79 |
| 100 | L25SN | 1.590 | 9.2 | 7.9 | 21 | 18 |
| 300 | L25S | 0.167 | 2.9 | 2.3 | 63 | 50 |
| 300 | L25SN | 1.702 | 7.5 | 5.8 | 16 | 12 |

- NCMC still gives 2.1-2.6x more theta and s ESS per iteration. Each of its iterations now
  costs 10-12x one `L25S` iteration, against 6-7.5x in the prototype study. The cache made
  the `L25S` iteration cheaper, but NCMC's cost lies elsewhere: per iteration, 2 x 8
  frozen-model rebuilds (`.Precompute()`, slices, messages) and about 950 relaxation moves.
  Those moves are dominated by path proposals, not by characters.
- Per CPU hour, NCMC is therefore 4-5.5x **worse** than `L25S` for theta and s on this tree.
  It was **not** added to the engine. The prototype study's open question ("an engine with
  incremental pruning would change this verdict") is answered no, for this implementation.
- Both samplers are far from converged at these lengths. ESS values of 6-32 per chain make
  every ratio here provisional (+/-30-50%).

## Deviations and limits

- `base` is b99e810 (the #43 head at the time). `inc` is the package code of a1bd997. That
  commit's RcppExports were regenerated on Hamilton; the committed b7f7f2a has the same
  generated files.
- Two fixed trees only, and one model cell, as in `dev/mixing/`. Three replicates per cell for
  the paired benchmark and four per cell for convergence. ESS values of 10-40 carry roughly
  +/-30% error. The CPU ratios are tight (within-cell spread under 3%), because each pair ran
  on one node.
- CPU time is the whole `TempoMCMC()` call (user + system), scaled to the sampling phase.
- NCMC was measured with the `dev/mixing/` prototype, not as an engine move. It ran 2 replicates
  of 3000 iterations, on the large tree only, at 100 and 300 characters.

## Reproduce (Hamilton; `lib-base` = b99e810, `lib-inc` = this branch)

```sh
# per task: the same seeded chain under each library, on one node
R_LIBS=lib-base EVOTEMPO_MIXING_OUT=bench-base Rscript data-raw/mixing-schedule.R chain <tree> <nChar> new <rep> 5000 2000
R_LIBS=lib-inc  EVOTEMPO_MIXING_OUT=bench-inc  Rscript data-raw/mixing-schedule.R chain <tree> <nChar> new <rep> 5000 2000
# large-tree convergence, replicates 11-14
R_LIBS=lib-inc  EVOTEMPO_MIXING_OUT=conv-inc   Rscript data-raw/mixing-schedule.R chain 2 <nChar> new <rep> 20000 2000
Rscript data-raw/incremental-pruning-bench.R bench-base bench-inc conv-inc
```
