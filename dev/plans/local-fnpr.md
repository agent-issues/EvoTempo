# Local fixed-node-height prune and regraft (`fnpr`): plan (#110)

## Why

Penguin stage-1 topology does not converge (#110). The only topology move, `slide`, draws the
regraft height uniformly on [age(v), root age]. Shallow crowns (Eudyptes, Spheniscus, about 1 Myr
under a 61 Myr root) are therefore almost never proposed into. `slide` acceptance is 4.0% even at
beta = gamma = 0, where the path moves are exact. RevBayes reaches 11-15% with `mvFNPR` and `mvNarrow`,
which keep node heights and regraft locally.

## Proposal

State: tree (ages, parents, tip types), node bins `bin`, branch paths `paths`.

1. Choose v uniformly from `.SlideNodes(tree)`: the path-carrying non-root nodes whose parent p is a
   non-root bifurcation. Let pp = parent(p) and sib = the other child of p.
2. The **pruned tree** is the tree with v's clade and p removed and sib hanging from pp. Its
   lineages crossing age(p) are `C = .Crossing(tree, age(p), p, clade(v))`, which includes sib.
3. Draw an auxiliary threshold H = age(p) + (R - age(p)) * B^(-U), where U ~ U(-1, 1), R is the
   root age and B = 100 (`.fnprSpan`).
   - With probability 1/2 (U < 0), H > R, so every crossing lineage is a target: a global
     `mvFNPR`.
   - Otherwise H - age(p) is log-uniform over two orders of magnitude of the excess above age(p).
   - A first draft used U ~ U(0, 1). It reached the root with probability 0, so regrafts across
     the root never happened; the unit test on the 5-tip fixture caught this.
4. Targets: T = {y in C, y != sib : MRCA_pruned(sib, y) has age <= H}. MRCA_pruned is the most
   recent common ancestor in the pruned tree; node ages are those of the current tree, which the
   pruned tree shares. If T is empty, return NULL (rejection).
5. Choose y uniformly from T. Regraft: p becomes the parent of v and y, and a child of
   x = parent(y), at its unchanged age(p). In the pruned tree sib's parent is pp.
6. Tempo:
   - draw the bin of p from `.NodeLogWeights(newTree, p, lambda, bin[x], kids = (y, v), ends = bin[c(y, v)])`, exactly as `slide` does;
   - redraw frozen bridges on {p, y, v, sib} between their (new) end bins;
   - all other bins and paths are kept.

## Reversibility and Hastings ratio

The reverse move from the proposed state picks the same v (so the same p), the same U (so the same
H, since age(p) and R are unchanged) and target sib.

- **Same pruned tree.** Removing v's clade and p from the proposed state gives the same pruned
  tree, with the same node ages. Hence the same crossing set C and the same MRCA ages.
- **Equal target counts.** "MRCA_pruned(a, b) has age <= H" is an equivalence relation on C, by the
  three-point inequality of ancestor ages in a rooted tree:
  age(MRCA(a, c)) <= max(age(MRCA(a, b)), age(MRCA(b, c))). Both MRCAs lie on b's ancestral path,
  so the older of them is a common ancestor of a and c. MRCA(a, c) is that node or a younger one,
  for instance when the two MRCAs are the same node. So if both ages are <= H, so is age(MRCA(a, c)). sib and y are in the same class, so T_fwd = class \ {sib} and T_rev = class \ {y} have
  the same size. The 1/|T| terms cancel.
- **Height.** No new height is drawn. The density of H is the same in both directions, since H is a
  function of U, age(p) and R only, so it cancels; the move is a mixture over U of reversible kernels.
- **Choice of v.** The reverse picks v from `.SlideNodes(newTree)`. Its size is twice the number of
  non-root bifurcations that are not sampled-ancestor nodes. The move changes neither the root nor
  any `sampled` flag, so the size is the same in both directions and the term cancels. (`slide`'s
  log|cand| ratio is likewise always 0.)
- **Tempo.** The forward proposal density is q_bin(new bin of p) * q_bridges(new paths on p, y, v, sib).
  The reverse is q_bin(old bin of p, computed in the old geometry with kTop = bin[pp], kids
  (sib, v)) * q_bridges(old paths on p, sib, v, y in the old tree).

Putting these together:

    log H = [log q_bridges_old(p, sib, v, y) + log q_bin_old(bin[p])]
          - [log q_bridges_new(p, y, v, sib) + log q_bin_new(bin'[p])]
          (|T_fwd| = |T_rev|, so the count terms are omitted; the code asserts the equality
           and that the old sib is a reverse target with `stopifnot`)

This is `slide`'s Hastings ratio without its height draw: in `slide`, `at` is uniform on
[age(v), R] in both directions and cancels. `.TreeMove()` adds the exact augmented density, the
character likelihood and the 2^|V| factor (unchanged, since the number of bifurcations is fixed).

**Every redrawn branch has its reverse density** (AGENTS.md, Tree-move proposals): p (new top x), y
(new top p), v (new top bin of p) and sib (new top pp). All four are redrawn and scored in both
directions.

## Edge cases

- **Root.** p is never the root (`.SlideNodes`). x may be the root, and pp may be the root. Targets
  never include the root branch (`.Crossing` excludes it), so with an origin the stem never changes.
- **Sampled ancestors.**
  - p is never sampled (`.SlideNodes`).
  - Sampled-ancestor tips have zero length, so they never cross age(p) and are never targets.
  - y may be a sampled-ancestor node's branch; regrafting onto it is valid.
- **`condition = "extant"` with a root start.** A move that leaves a root daughter without extant
  tips gets an exact augmented density of -Inf and is rejected, as for `slide`.
- **Trilobite pipeline.** Its ingroup constraint is enforced by rebinding `.ProposeTreeMove`. A new
  `switch` branch there is covered automatically.

## Engines

The proposal is R code that calls the same kernels as `slide`: `.NodeLogWeights`, `.DrawBridges`,
`.LogBridges` and `.Rebuild`, which dispatch to C++ under `engine = "cpp"`. The new R-level random
draws (`.Sample1`, `runif`) are the same in both engines, so seeded chains must stay identical
(`test-Engine.R`). No new C++ kernel is needed. The MRCA computation is O(n) per proposal and runs
once.

## Files

- `R/TreeMoves.R`: `.PrunedParent()`, `.MRCAAge()`, `.FnprTargets()`, and the `fnpr` branch of `.ProposeTreeMove`.
- `R/MCMC.R`: `fnpr` in `moveNames` and `treeMoves`; documentation.
- `tests/testthat/test-TreeMoves.R`, `test-Engine.R`: validity, symmetry of target counts, and
  engine identity.
- `data-raw/rung5-fixed-parameters.R`: optional move set (`EVOTEMPO_RUNG5_TOPOLOGY=fnpr`) to check
  `fnpr` against the simulator without `slide`.
- `data-raw/empirical/penguins/cet/common.R`: `CET_FNPR` count scale, default 0 so that production is
  unchanged.

## Correctness criteria

1. Prior-only rung 5 at fixed parameters (`rung5-fixed-parameters.R`) with `fnpr` as the only
   topology move: tree statistics, including `cherries`, agree with rejection simulation
   (|z| < 3 across statistics).
2. Seeded R and C++ chains with `fnpr` are identical.
3. Unit test (`test-TreeMoves.R`, "fnpr Hastings ratio matches the proposal densities in both
   directions"): over a random walk of states, with crown and origin starts, `logHastings` equals
   log q(reverse) - log q(forward), recomputed from the definition of the move. The test also checks
   that the target counts are equal, that the old sib is a reverse target, that ages are unchanged
   and that the topology changes. It fails under each of the review's Hastings mutants (#118):
   - reverse bin term dropped;
   - sib's reverse bridge dropped;
   - a constant added;
   - forward bridges dropped.
4. Afterwards, on Hamilton: SBC as in #21, and the penguin payoff (clade spread per CPU-hour).

## Review of the derivation

A `math-prover` subagent checked points 1-5 above against the code `slide` reuses, before
implementation. It verified:
- the pruned tree is the same in both directions;
- the MRCA relation is an equivalence;
- the auxiliary threshold H is valid;
- the set of redrawn branches is complete, including when x or pp is the root, with an origin start,
  and with sampled-ancestor or fossil endpoints;
- excluding sib and returning NULL keep detailed balance.

It found one error: the |cand| ratio is identically 0, not variable as first claimed. Two
implementation notes from that review are applied:
- the MRCA is computed on the pruned parent vector;
- x is read before `parent` changes.
