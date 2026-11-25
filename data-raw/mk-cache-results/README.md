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

## Longer run

The longer run is recommended on Hamilton, in a new library (not
`lib-48e6bf7`):

    Rscript data-raw/mk-cache-benchmark.R 1 1 200 10 mk-cache-b1g1.csv
