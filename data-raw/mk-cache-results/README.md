# Cached character likelihood: time per proposal (#110 move 2)

`data-raw/mk-cache-benchmark.R` on the penguin data (55 tips, 201 characters in
six K-partitions, 132 Mkv constant patterns), fixed parameters, one move type per
run, `nIter = 8`, 5 proposals per iteration (40 proposals; the parameter sweep:
6 proposals per iteration).  "before": one-shot pruning and the R Mkv correction
(`ConstantPatternLogLik()`); "after": cached pruner for both.  Laptop smoke test
(Windows, R-devel), about 10 s per cell, so each figure is good to perhaps
+/-20%; `same` = samples of the two modes agree to 1e-10.

| move | b1g1 before | after | speed-up | b0g0 before | after | speed-up |
|---|---|---|---|---|---|---|
| slide | 5.92 | 3.13 | 1.9 | 5.47 | 2.84 | 1.9 |
| saJump | 3.48 | 2.24 | 1.6 | 3.65 | 2.44 | 1.5 |
| nodeAge | 5.58 | 2.94 | 1.9 | 6.60 | 3.73 | 1.8 |
| fossilAge | 6.31 | 4.72 | 1.3 | 6.10 | 3.14 | 1.9 |
| originAge | 11.93 | 10.15 | 1.2 | 10.27 | 9.09 | 1.1 |
| node | 4.68 | 3.22 | 1.5 | 2.41 | 1.88 | 1.3 |
| tip | 4.27 | 1.98 | 2.2 | 2.63 | 1.13 | 2.3 |
| branch | 3.07 | 1.15 | 2.7 | 3.86 | 0.92 | 4.2 |
| shift | 2.80 | 1.48 | 1.9 | - | - | - |
| tree | - | - | - | 6.81 | 3.11 | 2.2 |
| parameters | 8.15 | 5.57 | 1.5 | 5.88 | 5.05 | 1.2 |

(ms per proposal; CSVs `smoke-b1g1.csv`, `smoke-b0g0.csv`.)  `shift` calls the
character likelihood only because the benchmark fixes `rate`; with `rate` free, as
in production, the move is character-invariant and unaffected.

Weighted by the production b1g1 schedule (`node` 14, `tip` 11, `branch` 24,
`shift` 1, `nodeAge` 10, `fossilAge` 9, `slide` 6, `saJump` 4, `originAge` 1, and
a 6-parameter sweep), with `shift` unchanged: 412 ms -> 240 ms of move work per
iteration, 1.7x.

The character likelihood itself fell from 4.4 ms to 0.4-0.7 ms per call. Most of
that came from moving the Mkv correction (3.4 ms in R) into the C++ pruner; the
caching of the 201 data columns, which cost 0.6 ms one-shot, saves less on a
tree this small.  What remains of a proposal is R-level move work (proposal,
bridges, `.LogAugmented()`, `.Messages()` on acceptance), which now dominates:
the 20x per-proposal gap to RevBayes estimated in #110 is not in the character
likelihood.

## Fixed cost of an evaluation: one C++ call per character set

A review of the cache (PR #119) found that about 70% of a cached evaluation did
no pruning: twelve `.MkCpp()` calls, each with a content-keyed lookup and a copy
of the tip states, behind three layers of R wrappers.  `PartitionedMk()` now
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
| before (PR #119) | 202-209 | 1 |
| thin R path, twelve C++ calls (`twelve`) | 187-192 | 0.92 |
| one C++ call | 164-196 | 0.86 |

The two effects are about 8% (R wrappers) and 6% (single entry point) of an
iteration.  The character likelihood falls from 33% to 16% of the profile
(`Rprof`, 80 iterations), behind `.Messages()` (20%) and `.Precompute()` (18%).

What is left of an evaluation is mostly `Pruner::set_tree()`, not pruning:
`mk-set-benchmark.R phases 20000` gives 0.13-0.14 ms of a 0.20-0.23 ms
evaluation with one tip edge changed (`checkpoint()`: 0.007 ms).
`dev/plans/mk-edge-updates.md` proposes how to avoid it.

Longer run (recommended on Hamilton, in a new library, not `lib-48e6bf7`):

    Rscript data-raw/mk-cache-benchmark.R 1 1 200 10 mk-cache-b1g1.csv
