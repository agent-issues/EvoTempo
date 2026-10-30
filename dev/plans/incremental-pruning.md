# Incremental character pruning in `TempoMCMC()` (#26, item 3)

## Problem

With characters, every move of `TempoMCMC()` scores its proposal with `MkLikelihood()`. Under the
C++ engine that called `MkPrime::log_likelihood()`, which builds a fresh pruner and prunes every
node for every character. A local path move (`node`, `tip`, `branch`) changes the tempo-weighted
length of one to three branches, and a tree move (`nodeAge`, `fossilAge`, `slide`, `saJump`)
changes a few branches and at most two parent pointers. Yet each one paid for a whole-tree
pruning. On the large `dev/mixing/` tree (272 path branches), `TempoMoves()` runs about 137
local moves per iteration, and most of the time goes on this pruning.

## Design

1. **Owner.** `.CharacterLikelihood(characters, nState, gamma)` (R/Characters.R) is called once
   per `TempoMCMC()` run. It returns the function `CharLL(tree, paths, par, set)` that the moves
   already take. Under the C++ engine, the closure holds an external pointer
   (`Rcpp::XPtr`, deleted by its finaliser) to a `CharacterCache` (src/characters.cpp). The
   cache owns a `MkPrime::Pruner<JCModel>`.
2. **Owned storage.** The MkPrime lesson (agent-issues/MkPrime#460) is that nothing kept across
   `.Call()`s may point into R memory, because ALTREP wrappers can swap their payload. The pruner
   copies the tip states at construction, and copies the parent vector and branch lengths into
   `std::vector`s at every call. No R object is held between calls. A serialised pointer comes
   back null; using it is an error, not a crash.
3. **Evaluation.** Each call passes the tree's `parent` vector and every branch length
   (`.CharacterLengths()`, unchanged). If the parent vector equals the pruner's own topology,
   the cache calls `set_edge_length()` on every edge. This is a no-op for unchanged lengths and
   otherwise marks the branch's parent and its ancestors dirty. If the topology differs (`slide`),
   `set_tree()` diffs the children of every node and marks only the nodes whose children or child
   lengths changed, and their ancestors. `log_likelihood()` then recomputes the dirty nodes only.
   Edges are the non-root nodes in increasing order, as in `.MkLikelihoodCpp()`. Tree moves
   never renumber nodes, so an edge keeps its index.
4. **Proposals.** The first evaluation after a resolution opens a `checkpoint()`. The moves report
   each Metropolis decision with `.ResolveCharacters(CharLL, accepted)`, which calls `accept()`
   or `reject()`. `reject()` restores the journaled CLs, the topology and the lengths in place,
   without recomputation. The calls are in `.PathMove()`, `.TreeMove()`, `.AgeMove()`,
   `.ShiftMove()`, `.JointSweep()`, the parameter sweep (`.ParameterSweep(Resolve = )`) and after
   each assignment of `st$charLL` in `TempoMCMC()`. For any other function, such as the plain
   closures that tests pass to the moves, `.ResolveCharacters()` is a no-op.
5. **Correctness does not depend on the bookkeeping.** The pruner's own topology and lengths
   are the reference for the next diff. A checkpoint is always a consistent snapshot
   (CLs, dirty flags, tree). So a proposal left unresolved, a missed `reject()` or a
   discarded R state costs only recomputation. It never gives a wrong value.
6. **Identity with full pruning.** A node's CLs are a deterministic function of its children's
   CLs and edge lengths, combined in the same order (decreasing edge index). By induction from
   the tips, an incrementally maintained pruner holds bitwise the same CLs as a fresh one. So
   `CharLL` returns exactly `.MkLikelihoodCpp()`, which itself agrees with the R reference to
   ~1e-15. Random-number use is unchanged, so seeded chains are identical with and without the
   cache, and between engines as before.
7. **R engine.** It keeps full re-pruning in R (`MkLikelihood(engine = "R")`) as the reference.
   Option `EvoTempo.incremental = FALSE` also restores full C++ pruning. It exists for tests and
   benchmarks only and is not documented for users.

## What is not incremental

- `.CharacterLengths()` still sums every branch's path at each call. That is O(path pieces),
  small beside pruning `nChar x nState` values per node.
- Moves that change every length recompute the whole tree, as before: `rate` in the sweep,
  `joint`, `tree`, and `shift` when `rate` is not compensated. Rejections there are still
  undone from the journal.

## Consumers that rebind the character likelihood

`data-raw/rung7-power.R` tempers the character likelihood by rebinding a namespace function.
`TempoMCMC()` no longer calls `MkLikelihood()` directly, so the script now rebinds
`.CharacterLikelihood()` and keeps the resolve hook. Its existing per-chain check (the final
`charLL / power` against a recomputation) would catch a tempering that silently stopped working.

## Tests

`tests/testthat/test-CharacterCache.R`:

- every evaluation during 120 random accepted and rejected path and tree moves, plus shifts,
  is identical to `.MkLikelihoodCpp()` and within 1e-10 of the R engine, and so is the stored
  `charLL` after every move;
- on a 16-tip balanced tree: a tip-branch change recomputes exactly its 4 ancestors; a rejection
  then recomputes nothing; a regraft recomputes the union of the two parents' ancestor paths;
  unresolved and repeated resolutions stay correct;
- full-pruning fallbacks (R engine, the option off, zero characters), error paths;
- seeded `TempoMCMC()` chains with every move type are identical with and without the cache.
