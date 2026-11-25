# Cached character likelihood and C++ Mkv correction: time per proposal (#110 move 2)

`data-raw/mk-cache-benchmark.R` on the penguin data: 55 tips; 201 characters in
six K-partitions; 132 Mkv constant patterns. Parameters are fixed, each move type
runs alone, `nIter = 8` with 5 proposals per iteration (40 proposals), and the
parameter sweep makes 6 proposals per iteration. PR #119 makes two changes, and
the three modes keep them apart:

| mode | pruning | Mkv correction |
|---|---|---|
| `main` (the code path before #119) | one-shot | R (`ConstantPatternLogLik()`) |
| `cppMkv` | one-shot | C++ (`PatternLogLik()` → `.MkSiteLogLik()`) |
| `cached` | cached | C++ |

These are laptop smoke tests (Windows, R-devel, other jobs running) of about a
minute per cell, so each figure is good to perhaps ±20-30%. `same` means that
the samples of all three modes agree to 1e-10. CSVs: `smoke-b1g1.csv`,
`smoke-b0g0.csv`, and `-call.csv` for the time per call.

## Per iteration (production schedule, with the sweep)

| cell | main | cppMkv | cached | Mkv in C++ | cache | total |
|---|---|---|---|---|---|---|
| b1g1 | 550 ms | 334 ms | 292 ms | 1.64x | 1.14x | 1.88x |
| b0g0 | 474 ms | 323 ms | 283 ms | 1.47x | 1.14x | 1.67x |

The independent review of #119 timed 40 iterations of the b1g1 production
schedule on a quieter machine, with alternating runs. It measured 754 ms for
`main` and 440 ms after (1.7x). In that run the cache alone was worth **about
5%** of an iteration (3-7%: `cppMkv` 408 and 399 ms, cached 395 and 370 ms).

Either way, **most of the gain comes from moving the Mkv correction into C++.
The cache adds 5-15%.**

## Per call of the character likelihood

One `PartitionedMk()` call makes 12 `.MkCpp` calls: 6 partitions and 6
constant-pattern matrices. In the cached mode, one tip edge changes per call.

| | ms per call |
|---|---|
| main | 3.9-4.0 |
| cppMkv | 1.1-1.3 |
| cached | 1.0-1.1 |

The review measured 5.8, 1.66 and 1.27 ms; the cached figure falls to 0.89 ms
when nothing has changed. Caching does little per call on this tree, for two
reasons:

- **The changes reach a large share of the tree.** One changed tip edge dirties
  about 30% of the internal nodes, and in the production schedule a call
  recomputes 25% of them on average.
- **Most of a cached call is fixed cost.** About 70% of it does no pruning: the
  R wrappers, the 12 calls into C++, the key lookup and the checkpoint.

## Per proposal (ms)

| move | b1g1 main | cppMkv | cached | b0g0 main | cppMkv | cached |
|---|---|---|---|---|---|---|
| slide | 8.23 | 6.21 | 5.90 | 8.52 | 4.54 | 4.01 |
| saJump | 6.07 | 4.43 | 3.43 | 5.30 | 3.33 | 2.28 |
| nodeAge | 8.40 | 4.18 | 3.92 | 7.82 | 4.59 | 4.62 |
| fossilAge | 9.19 | 6.52 | 5.51 | 8.20 | 4.43 | 5.28 |
| originAge | 14.80 | 14.73 | 12.12 | 13.96 | 12.42 | 11.74 |
| node | 4.78 | 4.00 | 4.09 | 3.95 | 3.22 | 2.88 |
| tip | 5.60 | 4.93 | 3.30 | 5.11 | 3.07 | 2.44 |
| branch | 3.42 | 2.38 | 1.37 | 4.01 | 3.09 | 0.75 |
| shift | 4.91 | 2.53 | 1.93 | - | - | - |
| tree | - | - | - | 8.96 | 6.01 | 3.21 |
| parameters | 9.16 | 7.61 | 5.74 | 9.35 | 6.85 | 5.64 |

`shift` calls the character likelihood only because the benchmark fixes
`rate`. With `rate` free, as in production, the move is character-invariant and
is unaffected.

## What dominates now

The review's Rprof of 80 iterations of the b1g1 production schedule, after
#119:

| item | share of the run |
|---|---|
| character likelihood (C++ call 17%, its R wrappers 12.5%, `.CharacterLengths()` 4%) | 34% |
| `.Messages()` | 15% |
| `.Precompute()` | 13.5% |
| bridges and path draws | 11.5% |
| `.LogAugmented()` | 5% |
| `.LogWeight()` | 5% |

About 70% of the run's self time is in C++ kernels. **The character likelihood
is still the largest single item, but no single item holds the ~20x gap to
RevBayes estimated in #110.** Removing the likelihood entirely would give about
1.5x.

## Fixed cost of an evaluation: one C++ call per character set (after #119)

As above, about 70% of a cached evaluation did no pruning: twelve `.MkCpp()`
calls, each with a content-keyed lookup and a copy of the tip states, behind
three layers of R wrappers.  `PartitionedMk()` now
registers the six partitions and six constant-pattern matrices once
(`.MkSet()`) and scores them in one call (`.MkSetLogLik()`), with one diff of
the tree for all twelve pruners.  The values and the order in which they are
combined are unchanged, so chains are bitwise those of the earlier build.

`data-raw/mk-set-benchmark.R micro 3000` (ms per `PartitionedMk()` call; laptop,
other jobs running, so read the ratios; two runs each, less one disturbed
measurement in each of the first two columns of the last row):

| since the last call | before | thin R path, twelve C++ calls | one C++ call | of which C++ |
|---|---|---|---|---|
| nothing changed (no pruning) | 0.43-0.48 | 0.21-0.28 | 0.04-0.05 | 0.013-0.017 |
| one tip edge changed | 0.69-0.74 | 0.44-0.50 | 0.26-0.27 | 0.23-0.24 |
| every edge changed | 0.92 | 0.57-0.60 | 0.37-0.39 | 0.35 |

So the R wrappers cost about 0.2-0.3 ms per evaluation and the eleven extra C++
calls, with their lookups and copies, another 0.2 ms; about 0.03 ms of R and
0.015 ms of glue remain.

`mk-set-benchmark.R iter 40` (production b1g1 schedule, about 84 evaluations per
iteration, same seed; ms per iteration, four runs each, alternating):

| | ms per iteration | relative |
|---|---|---|
| before (the cached build above) | 202-209 | 1 |
| thin R path, twelve C++ calls (`twelve`) | 187-192 | 0.92 |
| one C++ call | 164-196 | 0.86 |

The two effects are about 8% (R wrappers) and 6% (single entry point) of an
iteration.  The character likelihood falls from 33% to 16% of the profile
(`Rprof`, 80 iterations), behind `.Messages()` (20%) and `.Precompute()` (18%).

What is left of an evaluation is mostly `Pruner::set_tree()`, not pruning:
`mk-set-benchmark.R phases 20000` gives 0.13-0.14 ms of a 0.20-0.23 ms
evaluation with one tip edge changed (`checkpoint()`: 0.007 ms).
`dev/plans/mk-edge-updates.md` proposes how to avoid it.

## Longer run

The longer run is recommended on Hamilton, in a new library (not
`lib-48e6bf7`):

    Rscript data-raw/mk-cache-benchmark.R 1 1 200 10 mk-cache-b1g1.csv
