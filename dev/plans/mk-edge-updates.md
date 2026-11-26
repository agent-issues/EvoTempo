# Edge-length updates for the cached character likelihood (proposal, #110)

Not implemented. A proposal for the step after the single entry point
(`.MkSetLogLik()`), with the measurements behind it.

## What an evaluation costs now

Penguin data (55 tips, six K-partitions and six constant-pattern matrices,
twelve `MkPrime::Pruner`s), one tip edge changed per proposal, timed inside the
Pruner by `data-raw/mk-set-benchmark.R phases 20000` (laptop; ms per evaluation
of all twelve):

| phase | one tip edge | every edge | one tip edge, `set_edge_length()` |
|---|---|---|---|
| `accept()` / `reject()` | 0.007 | 0.02-0.03 | 0.007 |
| `checkpoint()` | 0.007 | 0.008-0.011 | 0.007 |
| `set_tree()` | 0.13-0.14 | 0.14-0.19 | 0.0005 |
| `site_likelihoods()` (the pruning) | 0.06 | 0.18-0.25 | 0.06 |
| total | 0.20-0.23 | 0.35-0.48 | 0.06-0.08 |

`checkpoint()` does not matter. `set_tree()` does: it is about 60% of an
evaluation in which one edge changed, more than the pruning it saves, because
for every Pruner it validates the whole tree, builds a vector of child vectors,
diffs it against the last tree and rebuilds the postorder. With the R and Rcpp
glue now at about 0.03 ms, `set_tree()` is the largest fixed cost left: about
0.13 ms of an evaluation in which one length changed (0.29 ms as a whole
`PartitionedMk()` call, by `micro`, where each call is an accepted proposal).

## Proposal

`Pruner::set_edge_length(edge, length)` already exists and marks the same nodes
dirty as `set_tree()` does for a changed length (the parent of the edge and its
ancestors), so the result stays bitwise that of a one-shot evaluation.

1. **No move needs to report anything.** `Evaluate()` in `src/characters.cpp`
   already compares the new tree with the last one (`x == s.last`, `NDiff()`).
   When `parent` and `child` are equal, the edges whose lengths differ are known
   from that comparison, and `edge` is their index in the arrays last given to
   `set_tree()`, which are those same arrays. So: equal topology ->
   `set_edge_length()` for each differing edge; otherwise `set_tree()` as now.
   This keeps the property that the cache never holds chain state and that every
   new move is cached without hooks.
2. The comparison must be against the tree the Pruner actually holds, which
   after `reject()` is `s.prev` swapped back into `s.last`; `Evaluate()` already
   maintains that.
3. Moves that change the topology (`slide`, `saJump`) still pay for
   `set_tree()`. They are about 10 of 84 proposals in the production schedule.
   Making those cheap is a change to MkPrime (a `set_tree()` that reuses its
   buffers, or a regraft call), not to EvoTempo.
4. A proposal that changes every length (the `rate` update of the parameter
   sweep, and the `s` update when gamma != 0) gains little: there `set_tree()` is a third or so of the cost and
   the pruning the rest.

Expected gain: at most 0.13 ms on the roughly 70 evaluations per iteration that
change lengths but not topology (the path and age moves, and the `rate` and,
with gamma != 0, `s` updates, which change every length; the `lambda`, `mu`,
`psi` and `theta` updates change none and already return before `set_tree()`),
so up to 9 ms of 175 ms (5%). Cost: a dozen lines in `Evaluate()`, no change to MkPrime, R or any
move.

## Tests it would need

- The existing fuzz in `test-MkCache.R` (set against one uncached call per
  matrix, `identical()`, with `EVOTEMPO_MK_CACHE_CHECK = 1`), which already mixes
  length-only and topology proposals with accepts, rejects and returns to old
  states. It must reach both branches: assert on a counter of
  `set_edge_length()` calls, as `.MkCache()$recomputed` does for nodes.
- A length-only proposal straight after a rejected topology proposal, and the
  reverse, since those are the cases where the Pruner's tree and the last tree
  given could come apart.
- Zero-length edges and an unchanged tree (no call at all).
