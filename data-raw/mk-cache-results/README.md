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

Longer run (recommended on Hamilton, in a new library, not `lib-48e6bf7`):

    Rscript data-raw/mk-cache-benchmark.R 1 1 200 10 mk-cache-b1g1.csv
