# Cached character likelihood ("move 2" of #110)

Goal: each `TempoMCMC()` proposal recomputes only the dirty part of the Mk
character likelihood, instead of pruning all characters on the whole tree.

## Baseline (penguins, 55 tips, 201 characters, b1g1, laptop)

| piece | ms per call |
|---|---|
| rebound `MkLikelihood()` (`PartitionedMk()`), total | 4.4 |
| C++ pruning of the six K-partitions (one-shot, the then `.MkLikelihoodCpp`) | 0.6 |
| R Mkv correction, 132 constant patterns (`ConstantPatternLogLik()`) | 3.4 |

Whole proposals cost 3-12 ms (slide 5.9, nodeAge 5.6, fossilAge 6.3, path moves
3-5, originAge 12; `data-raw/mk-cache-results/`), so the character likelihood is
a third to most of a proposal, and **most of it is the
R-level Mkv correction, not the C++ pruning**. Caching the C++ partials alone saves
at most 0.6 ms per proposal. The plan therefore caches both, by running the
constant-pattern columns through the same cached pruner. Even with a free
character likelihood, the other R-level work of a proposal (proposal, bridges,
`.LogAugmented()`, `.Messages()` on acceptance) bounds the speed-up per proposal
at about 1.5-3x, not the 20x estimated in #110.

## Where the cache lives

Not in `TempoMCMC()` and not in `st`. The penguin and trilobite pipelines rebind
`MkLikelihood()` in the namespace and call the original once per K-partition,
with different matrices and `nState`, and they save `st` to disk and reload it
(checkpoint/resume), which would leave an external pointer null. So:

- `src/characters.cpp` keeps a small registry of `MkPrime::Pruner<JCModel>`
  objects (at most 32, least recently used evicted), keyed by `nTip`, `nState`
  and the full contents of the tip-state matrix (`memcmp`).
- The cpp branch of `MkLikelihood()` looks up (or creates) the entry and asks it
  for the likelihood of the given tree and lengths. A new internal
  `.MkSiteLogLik()` does the same and returns one log-likelihood per character,
  for the Mkv correction.
- Calls from `.mkOriginal()` in the pipelines hit the cache with no change to
  `TempoMCMC()` or to any move. The FNPR move (#110 move 1) gets caching for free.

## Never trust the cache to hold the chain state

Every call passes the whole tree (parent of each node) and the whole vector of
character lengths. The entry compares them with the tree and lengths of its last
evaluation and lets `Pruner::set_tree()` mark the changed nodes and their
ancestors dirty (`set_tree()` diffs both topology and lengths, so a separate
`set_edge_length()` with edge-index bookkeeping is not needed). A missed
accept/reject signal therefore costs time, never correctness: state injected
behind the engine's back (resume), other trees scored in between
(`CetStartLogPost()`, `RecordingStats()`), or a new move type all just give a
larger diff.

`checkpoint()` / `accept()` / `reject()` are a speed layer that infers the MCMC
decision from the next call. Each entry keeps the inputs and site likelihoods of
its last evaluation (`last`) and, while a proposal is open, of the state before it
(`prev`). On a new call with input `x`:

1. `x == last`: return the stored result.
2. Proposal open and `x == prev` (the proposal was rejected and the current
   state is being rescored): `reject()`, which restores the journalled partials,
   and return the stored `prev` result.
3. Proposal open and `x` differs from `prev` in fewer nodes than from `last`
   (the typical rejected proposal followed by a new proposal from the current
   state): `reject()`, then continue from `prev`. Otherwise `accept()`.
4. `checkpoint()`, `set_tree(x)`, recompute the dirty nodes, store `x` as `last`.

The first call on an entry builds the tree without a checkpoint.

Edge index = order of the child node among non-root nodes, as in the one-shot
path, so it is stable across tree moves (node numbers never change; `slide`
changes parents, `saJump` changes a tip's length to or from zero).

Cached and one-shot results are **bitwise identical**: a node's partials depend
only on its children's partials and edge lengths, a clean node has exactly the
inputs it had when computed, and the site likelihoods are summed in the same
order with the same `-Inf` rule as `Pruner::log_likelihood()`. So seeded R and
C++ chains still match (`test-Engine.R`), and random-number use is unchanged.

## Proposal types

None needs code. What each changes, and so what the diff marks dirty:

| proposal | changes | dirty |
|---|---|---|
| `slide` | parent of pruned node and its parent, lengths of p, sib, v, target | the old and new paths to the root |
| `saJump` | tip length to/from 0, lengths of p and sib | those nodes' ancestors |
| `nodeAge`, `fossilAge` | lengths of the node's or tip's edges | their ancestors |
| `originAge` / `rootAge` | stem (ignored) and root-adjacent lengths; with gamma = 1 the redrawn start-branch paths | root region |
| `node`, `tip`, `branch` | with gamma != 0, lengths of the redrawn branches; with gamma = 0 lengths change only by round-off of the dwell sums | their ancestors (often nothing) |
| `tree`, `joint` | every path, so (gamma != 0) every length | everything |
| `shift` | with `rate` free (or gamma = 0) the character likelihood is invariant and `.ShiftMove()` does not call it; with `rate` fixed and gamma != 0, every length changes | - / everything |
| parameter sweep | only `rate` changes lengths (all of them) | everything; other parameters leave lengths unchanged, so the call returns the stored result |
| FNPR (move 1) | like `slide` | like `slide` |

Gamma-coupled lengths need nothing special: `.CharacterLengths()` integrates
`rate * exp(gamma x)` along each path, and its output is the input that is diffed.
Sampled ancestors are tips on zero-length edges; a zero length is an ordinary
length. The model (`JCModel(nState)`) and rate categories (`{1}`) are fixed for an
entry, because `nState` is part of the key; `set_model()` / `set_rates()` are never
called, and a change of `rate` reaches the pruner as a change of every length.

## Pipelines that rebind `MkLikelihood()`

- **Penguins**: partition calls go through `.mkOriginal()` and are cached. The
  Mkv correction moves to `.MkSiteLogLik()` on one constant matrix per K (the
  pattern's observed tips in state 0, the rest unknown; its site likelihood is
  exactly `P(all observed tips in state 0)`). `ConstantPatternLogLik()` stays as
  the R reference; check-mkv.R compares the two. The cached route is used only
  with the cpp engine and when the installed EvoTempo has `.MkSiteLogLik()`, so the
  scripts keep working with older builds (Hamilton `lib-48e6bf7`).
- **Trilobites**: the `cm$cpp` partitions are cached; the `CetPrune()` columns
  (ordered and ambiguous characters) stay one-shot. Out of scope here.

## Debug check and escape hatch

- `EVOTEMPO_MK_CACHE_CHECK=1` (environment variable, read on each call): after
  every cached evaluation, including the restore of a rejected proposal, the
  entry's partials are re-read (no recomputation) and compared with a fresh
  one-shot evaluation; any difference is an error. Tests run chains with it set
  in all four (beta, gamma) cells, with `slide`, `saJump`, `rate` changes and a
  two-partition rebinding wrapper.
- `options(EvoTempo.mkCache = FALSE)` turns the cache off (one-shot path), for
  benchmarks and as a fallback. `.MkCache(clear = TRUE)` empties the registry; `.MkCache()` reports its size and
  the number of partials recomputed.

## Correctness criterion

Cached == one-shot, bitwise, at every evaluation (debug check); identical seeded
R/C++ chains (`test-Engine.R`); 100% coverage of new lines.

## Outcome (after review)

Measured after implementation (`data-raw/mk-cache-results/README.md`): the
cache itself saves 5-15% of an iteration on the penguin data; moving the Mkv
correction into C++ gives most of the 1.7-1.9x.  One changed edge dirties about
30% of the internal nodes of this tree, and most of a cached call is fixed
cost.  Larger trees, or fewer, larger character matrices, would gain more from
the cache.
