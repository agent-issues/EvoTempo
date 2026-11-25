# Likelihood-weighted fixed-height regraft (`wfnpr`): derivation (#110, move 4)

`wfnpr` prunes a subtree, scores a bounded set of reattachment positions at the unchanged height of
its parent by their full target density (tree prior, tempo paths and characters), and chooses among
them in proportion to those scores. It generalizes `fnpr` (`dev/plans/local-fnpr.md`), and with one
alternative position it *is* `fnpr`.

## Notation

State: tree (parents, ages, tip types), node bins `bin`, branch paths `paths`. pi is the exact
target: `.LogAugmented()` plus the character log-likelihood (plus priors of quantities that do not
change here). As in `fnpr`:

- v is the pruned node, p its parent (a non-root, non-sampled bifurcation), pp the parent of p and
  sib the other child of p. p keeps its age a = age(p).
- The **pruned tree** is the tree without v's clade and p, with sib hanging from pp. A *lineage* is
  a branch of the pruned tree that crosses a; it is named by its lower node y, and x(y) is its
  parent in the pruned tree. The lineage of sib runs from pp to sib.
- Attaching at lineage y means: p becomes the child of x(y) and the parent of y and v.

For a lineage y, two tempo configurations are defined on it:

- **intact**, I_y: one path from x(y) to y, between the fixed bins bin[x(y)] and bin[y];
- **split**, S_y = (b_y, P_y, L_y, V_y): the bin b_y of p, the path P_y on p's branch
  (x(y) to p, bins bin[x(y)] to b_y), the path L_y on y's branch (p to y, bins b_y to bin[y]) and
  the path V_y on v's branch (p to v, bins b_y to bin[v]).

The state attached at y, written s_y, has S_y on lineage y and the intact path on every other
lineage. Everything else (all other bins and paths, all ages, the pruned topology) is the *base*,
shared by all s_y.

The frozen proposal densities are those `fnpr` uses:

- qI_y(I_y): the frozen bridge density of the intact path (`.LogBridges` on branch y hanging from
  x(y));
- qS_y(S_y) = qbin_y(b_y) * bridge(P_y) * bridge(L_y) * bridge(V_y), where qbin_y is the normalized
  `.NodeLogWeights(tree_y, p, lambda, bin[x(y)], kids = (y, v), ends = bin[c(y, v)])`.

**Key property.** qI_y and qS_y depend only on the base: on the geometry of the pruned tree, the
ages a and age(v), the fixed bins bin[x(y)], bin[y], bin[v], and the frozen model `set`. They do not
depend on where the subtree is currently attached. (The start age does not change, so `set` is the
same for every candidate.)

## The move

1. Choose v uniformly from `.SlideNodes(tree)` (as `slide` and `fnpr`). Current attachment: c = sib.
2. Draw the threshold H exactly as `fnpr`: H = a + (R - a) B^(-U), U ~ U(-1, 1), R the root age,
   B = `.fnprSpan`. The **class** of c is C = {c} + `.FnprTargets()` = the lineages crossing a
   whose MRCA with c in the pruned tree is no older than H.
3. **Bound.** If |C| - 1 > M - 1, choose a subset A of C \ {c} of size M - 1 uniformly at random;
   otherwise A = C \ {c}. The candidate set is S = A + {c}, so |S| <= M (`.wfnprMax`, 8). If A is
   empty, return NULL. Candidates are processed in increasing node number.
4. **Ensemble.** Draw the intact path I_c on c's lineage (pp to sib) once from qI_c. For each
   y in A, draw S_y from qS_y (bin of p, then bridges on p, y and v, in the tree attached at y).
   The intact paths I_y, y in A, and the split S_c are those of the current state. This defines
   s_y for every y in S: s_y carries S_y, and I_z for all z in S \ {y}.
5. **Weights.** w_y = pi(s_y) * qI_y(I_y) / qS_y(S_y) for every y in S.
6. **Choice (Metropolized Gibbs; Liu 1996).** Choose y in A with probability w_y / (W - w_c),
   W = sum of w over S, and accept with probability min(1, (W - w_c) / (W - w_y)).

Cost: |A| <= M - 1 evaluations of the character likelihood and of `.LogAugmented()` (the current
state's are known), with no reference draws.

## Why it is correct

Let q(S-ensemble | s_c) be the density of generating, from the current state s_c, the threshold, the
subset and all the draws of step 4.

**Claim 1: the same ensemble is reachable from every member, with density**

    q(E | s_y) = g(H) * u(S) * qI_y(I_y) * prod_{z in S, z != y} qS_z(S_z)

where g is the density of H and u(S) the probability of the subset.

- From s_y the move picks the same v with probability 1 / |`.SlideNodes`|. This count is the same in
  every s_y: it is twice the number of non-root, non-sampled bifurcations, and the move changes
  neither the root nor any `sampled` flag (as in `fnpr`).
- g(H) depends only on U, a and R: identical from every member.
- The pruned tree, and so the crossing lineages and their MRCA ages, are the same from every member.
  "MRCA no older than H" is an equivalence relation on them (`dev/plans/local-fnpr.md`), so the
  class C is the same set from every member.
- Subset: from s_y, S arises when the subset S \ {y} of C \ {y} is drawn. If |C| <= M this has
  probability 1; otherwise 1 / choose(|C| - 1, M - 1). Either way u(S) is the same for every y in S.
- Draws: from s_y the move draws I_y and S_z for z in S \ {y}; by the key property their densities
  are qI_y(I_y) and qS_z(S_z) whatever y is. The states it builds are exactly the s_z of the
  ensemble: s_z has S_z on z, the *new* I_y on y, and the *current* intact paths elsewhere.

**Claim 2: detailed balance.** Write E for the full ensemble ({I_z}, {S_z}, z in S, plus the base).
By Claim 1,

    pi(s_y) q(E | s_y) = g u [prod_z qS_z(S_z)] * pi(s_y) qI_y(I_y) / qS_y(S_y) = k(E) * w_y

with k(E) free of y. Hence, conditional on E, the current member is y with probability w_y / W, and
any kernel K_E on S that is reversible with respect to w gives

    pi(s_c) q(E | s_c) K_E(c, y) = k(E) w_c K_E(c, y) = k(E) w_y K_E(y, c) = pi(s_y) q(E | s_y) K_E(y, c).

Integrating over the parts of E that are not determined by (s_c, s_y) gives detailed balance of the
marginal kernel with respect to pi. (This is the multiple-proposal construction of Tjelmeland 2004,
with the weights pi(s_y) q(rest | s_y).) The dimension of the auxiliary draws is the same in both
directions, because every ensemble holds one intact and one split configuration per member lineage.

The Metropolized Gibbs kernel K_E(c, y) = w_y / (W - w_c) * min(1, (W - w_c) / (W - w_y)) for
y != c is w-reversible: w_c K_E(c, y) = w_c w_y min(1 / (W - w_c), 1 / (W - w_y)), symmetric in
(c, y). It moves at least as often as plain Gibbs (Liu 1996; Peskun ordering), which is why it is
used instead of always-accepted Gibbs.

**Gibbs or multiple-try Metropolis?** Multiple-try Metropolis (Liu, Liang & Wong 2000) with
independent draws per try needs reference draws from the selected state, about 2|A| - 1 likelihoods,
because its weights differ between the forward and reverse sets. Here the single shared draw I_c
makes the ensemble closed (Claim 1), so the weights are the same from every member and no references
are needed: |A| likelihoods. The shared draw is essential. If each candidate drew its own path on
c's lineage, s_z could not be reached from s_y (their paths on c's lineage would differ), and
Claim 1 would fail.

## Relation to `fnpr`, and the weights in terms of the code

For a candidate y, with tree_y the tree attached at y and `tree` the current one:

    log w_y - log w_c = [log pi(s_y) - log pi(s_c)]
                      + [log qS_c(S_c) + log qI_y(I_y)] - [log qS_y(S_y) + log qI_c(I_c)]

The first bracket of proposal terms is `fnpr`'s reverse density `OldQ(p, sib, v, y) + qbin_old`,
evaluated on the current tree; the second is its forward density `fwd$logQ + draw$logQ`. So
log w_y - log w_c is `fnpr`'s complete log acceptance ratio for target y (without the target-count
term, which is 0), and with |A| = 1 the move accepts with min(1, w_y / w_c): it has `fnpr`'s
kernel. (It consumes random numbers differently, so seeded chains differ.)

In the code, log w_c = 0 and

    log w_y = .LogAugmented(tree_y, ...) - .LogAugmented(tree, ...) + charLL_y - charLL_c
            + .LogBridges(tree, set, c(p, sib, v, y), paths) + log qbin_old(bin[p])
            - log qI_c(I_c) - log qS_y(S_y).

The 2^(bifurcations) factor of `.TreeMove()` is the same for every candidate and is omitted.

## Tempo paths, sampled ancestors and FBD terms

- **Paths.** Exactly the four branches `fnpr` redraws change between s_c and s_y: p, y, v and sib;
  each has its forward and reverse bridge density in w (AGENTS.md, Tree-move proposals).
- **Sampled ancestors.** p is never a sampled-ancestor node and v's parent is never the root
  (`.SlideNodes`). Sampled-ancestor tips have zero length and never cross a, so they are never
  targets. A target lineage may end at a sampled-ancestor node, and v may be a fossil tip or a
  sampled-ancestor node; neither changes any `sampled` flag, so psi factors are unchanged and enter
  only through `.LogAugmented()`.
- **FBD terms.** No age changes, so the only prior terms that differ between candidates are the
  speciation factor at p (lambda at bin b_y), the path terms on the four branches, and the
  conditioning (`condition = "extant"`: a candidate that leaves a root daughter without extant tips
  has pi = 0, weight 0, and is never chosen). All are inside `.LogAugmented()`, which is evaluated in
  full for every candidate.
- **Origin.** The root branch is never a target (`.Crossing`), so the stem is untouched.

## Infeasible draws and constraints

- A bridge that underflows signals `infeasibleProposal`; the whole move is then rejected, as in
  every other move. In exact arithmetic every bridge is feasible, so this is a numerical guard and
  not part of the kernel. It is no more exact than in the other moves (which draws are needed
  differs between members), and it fires more often per proposal: 2|A| + 1 bridges and |A| bin
  draws, against four bridges and one bin draw for `fnpr`.
- A candidate with a non-finite weight gets weight 0. If no candidate other than the current one
  has positive weight, the move returns NULL.
- W - w_y is computed as a log-sum-exp over the other members, never by subtraction.
- A hard topological constraint (the trilobite ingroup, applied by a wrapper of
  `.ProposeTreeMove()` that vetoes the returned tree) is a factor 1[constraint] of the target. The
  kernel above satisfies detailed balance pairwise with respect to pi, so vetoing a chosen state
  outside the constraint keeps detailed balance with respect to pi * 1[constraint]. Giving such
  candidates weight 0 before the choice would also be exact and waste fewer proposals; it needs
  the constraint inside the package, which it is not.

## Scoring interface

Every candidate's character log-likelihood comes from one function, `.CandidateLogLik()`, which
today calls the chain's `CharLL` closure (`MkLikelihood()`, whichever engine is active, and whatever
a pipeline has rebound it to): a full pruning pass per candidate. With a persistent
`MkPrime::Pruner` it becomes checkpoint, four `set_edge_length` calls and a re-attachment,
`log_likelihood`, `reject`; nothing else in the move changes.

## Review

A `math-prover` subagent tried to refute Claims 1 and 2, the subset bound, the code-level weights,
the constraint veto and the infeasibility guard against the code, before the move was written, and
confirmed each. Its two corrections (same kernel as `fnpr` but not the same seeded chain; the
infeasibility guard fires more often) are in the text above.

## Correctness criteria

1. For a fixed ensemble, the weights computed from two different members agree up to a common
   factor (unit test; catches a missing or unshared proposal term).
2. On a tiny tree with a sampled ancestor, at fixed ages and parameters, the chain's topology
   frequencies match the exact probabilities, proportional to `TreeLikelihood()` times the
   character likelihood (gamma = 0), over all topologies reachable at those ages.
3. Prior only: tree statistics agree with and without the move, and with rejection simulation
   (rung 5 at fixed parameters).
4. Seeded R and C++ chains with `wfnpr` are identical.
