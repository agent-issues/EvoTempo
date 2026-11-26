# Spliced fixed-height regraft (`sfnpr`): plan (#110, move 3)

## Why

`fnpr` (#118) and `slide` redraw frozen-model bridges on the four branches whose top changes
(p, y, v and sib) and redraw p's bin. The frozen model ignores the characters. When gamma != 0 the
redrawn paths change the tempo-weighted character lengths of four branches at random, so a
proposal that would improve the topology is often rejected because of the paths.

The cost is measured in the penguin chains: `slide` acceptance falls from 0.040 at gamma = 0 to
0.029 at gamma = 1, and falls smoothly along the gamma ladder (#110). `sfnpr` keeps the tempo
history of every lineage. Only the subtree's own branch, and only sometimes, gets a new path.

## State and notation

- **State:** the tree (ages, parents, tip types), node bins `bin`, and branch paths
  `paths[[u]] = list(t, k)`. `t` holds the jump ages in decreasing order; `k` holds the bins, one
  more than `t`. The path of branch u runs from `.TopAge(tree, u)` to `age[u]`, with
  `k[1] = .TopBin(tree, bin, u)` and `k[length(k)] = bin[u]`.
- **Notation:** v is the pruned node and p its parent (a non-root bifurcation). pp is p's parent,
  sib is v's sibling, and a = age(p). The target is y, with x = parent(y). The threshold H and the
  target sets T are exactly those of `fnpr` (`.FnprTargets`, `.fnprSpan`).

## Proposal

1. **Choose v, H and y** exactly as in `fnpr`: v uniform on `.SlideNodes(tree)`; H from U ~ U(-1, 1);
   y uniform on `T_fwd = .FnprTargets(tree, p, sib, clade(v), H)`.
2. **Rewire** as in `fnpr`: sib hangs from pp, p from x, y from p. Ages are unchanged.
3. **Splice and split the paths deterministically:**
   - `sib' = concat(paths[[p]], paths[[sib]])`, that is `t = c(t_p, t_sib)` and
     `k = c(k_p, k_sib[-1])`. This is continuous, because `k_p[last] = bin[p] = k_sib[1]`. It is
     the lineage's own tempo history from pp to sib.
   - Split `paths[[y]]` at a. With `j = sum(t_y > a)`, the new path of p is
     `p' = list(t = t_y[1:j], k = k_y[1:(j + 1)])`, the new path of y is
     `y' = list(t = t_y[-(1:j)], k = k_y[(j + 1):end])`, and the new bin of p is `k' = k_y[j + 1]`,
     the bin y occupies at age a. A jump exactly at a has probability 0.
4. **The pruned branch v.** Its interval [a, age(v)] is unchanged, but its top bin changes from
   k = bin[p] to k'.
   - If k' = k, keep `paths[[v]]` (no draw).
   - If k' != k, draw a frozen bridge on v from k' to `bin[v]` (`.DrawBridges`).
5. All other bins and paths are kept.

## Reversibility

- **The reverse move** from the proposed state chooses the same v, U (so the same H) and target sib.
  - The pruned tree is unchanged (as for `fnpr`), so the target counts are equal and the old sib is
    a reverse target. Shown for `fnpr` and checked by its `stopifnot`.
  - The reverse splices p' and y', which gives back y's old path exactly. It splits sib' at a, which
    gives back p's and sib's old paths, and returns `bin[p] = k_sib'` at a = the old k.
- **The pruned branch** is decided by the same condition in both directions, k' = k:
  - **k' = k:** both directions keep v's path. Each step is a deterministic bijection on the paths;
    it only relabels the jump times of the spliced and split paths. Its Jacobian is 1, so the paths
    contribute no proposal density.
  - **k' != k:** the forward move draws v's path from the frozen bridge k' -> bin[v] in the new
    tree. The reverse draws from k -> bin[v] in the old tree. Both densities enter the ratio.
- **The splice and split map** (paths of p, sib, y) <-> (paths of p', sib', y') is a bijection:
  - it is defined almost everywhere, since jumps exactly at a have probability 0;
  - it preserves the total number of jumps and copies every jump time, so it preserves the
    reference measure.
  - The number of jumps per branch changes, but the map pairs states of equal total dimension.

## Hastings ratio

- **Counts:** the choice of v and the target counts cancel (`fnpr`), and so does the density of H.
- **Bin of p:** set deterministically by the split, so there is no bin proposal term.
- **Ratio:**

      log H = 0                                                          if k' = k
      log H = log q_bridge(old v path; old tree, k -> bin[v])
            - log q_bridge(new v path; new tree, k' -> bin[v])           if k' != k

- **Target ratio:** `.TreeMove()` evaluates the exact augmented density, the character likelihood
  and the 2^|V| factor (unchanged) in full.

The exact augmented density (`.LogAugmented`) is a sum of three kinds of terms:
- terms over (bin, time interval) segments, `(Q_kk - sigma_k) dt + 2 lambda_k (G_k(hi) - G_k(lo))`,
  with G a function of (age, bin) only;
- `log Q` per jump;
- node factors (`lambda[bin]` at each bifurcation, `psi` at sampled-ancestor nodes), tip values,
  and the start terms.

So for k' = k the path terms of p, sib and y are unchanged by splicing, and only the character
likelihood can change. Otherwise the density changes by the v-branch path terms and by
`lambda[k'] / lambda[k]` at p (that is e^{beta (x_k' - x_k)}). This is the efficiency argument, not
a correctness requirement: the acceptance uses the full density.

**Character lengths.** For each lineage the tempo integral over time is preserved:
- sib' carries p's and sib's integrals;
- p' and y' share y's integral;
- v keeps its own when k' = k.

The character likelihood therefore changes through the topology, plus v's length when it is
redrawn. With gamma = 0, sfnpr and fnpr have the same target ratio; only their proposals differ.

## Edge cases

- **Root.** x or pp may be the root.
  - With `start = "root"` the root carries no path, but the root's children do. The new p' then
    runs from the root, and so does sib'.
  - With an origin the stem never changes, because y is never the root branch (`.Crossing`).
- **Sampled ancestors.**
  - p is never sampled.
  - y may end at a sampled-ancestor node; its lower part y' keeps that node's bin.
  - Sampled-ancestor tips have no path and are never v, sib or y.
- **Fossil tips:** sib, y or v may be a fossil tip; its bin is never changed.
- **`condition = "extant"`:** handled as for `fnpr`; the exact density rejects.
- **Trilobite constraint:** covered by the rebinding of `.ProposeTreeMove`.
- **Partial-likelihood cache (#119, #122):** `sfnpr` changes the character lengths of p, sib and y,
  and of v only when it redraws. The cache diffs lengths and topology, so it recomputes only those
  branches and their ancestors, which is fewer than `fnpr`.

## Engines

As for `fnpr`: the proposal is R code. The splice and split are pure list operations, and the
v bridge uses the dispatching kernels `.DrawBridges` and `.LogBridges`. Random-number use is
identical in both engines, so seeded chains agree. No new C++.

## Files

- `R/TreeMoves.R`: `.SplicePaths()`, `.SplitPath()`, and the `sfnpr` branch of `.ProposeTreeMove`.
- `R/MCMC.R`: `sfnpr` in `moveNames` and `treeMoves`; documentation.
- `tests/testthat/test-TreeMoves.R`:
  - splice/split round trip;
  - continuity of every path with its bins;
  - the Hastings identity recomputed from the definition over a random walk, with crown and origin
    starts and both the k' = k and k' != k branches;
  - preserved per-lineage tempo integrals.
- `tests/testthat/test-Engine.R`: `sfnpr` in the tree-move identity chain.
- `data-raw/rung5-fixed-parameters.R`: `EVOTEMPO_RUNG5_TOPOLOGY=sfnpr`.
- `data-raw/sfnpr-validation.R` (new): an exact topology oracle at fixed ages and parameters with
  **gamma = 1 and one character**.
  - With one character, `TreeLikelihood(character = )` gives the exact product-space likelihood, so
    P(topology) is proportional to its exponential over the topologies reachable by fixed-height
    regrafts.
  - Chains with `sfnpr` (plus path moves) are compared with the exact topology and clade
    probabilities. `fnpr` chains run alongside as a control. A prior-only case (no character) is
    included.
- `data-raw/empirical/penguins/cet/common.R`: `CET_SFNPR` count scale, default 0.

## Correctness criteria

1. **Unit test:** `logHastings` equals the ratio recomputed from the definition, over a random walk
   that hits both branches, with crown and origin starts. The test must fail under mutants:
   - the bridge ratio dropped;
   - the bridge ratio's sign flipped;
   - a constant added;
   - v kept when k' != k.
2. **Seeded R and C++ chains with `sfnpr`** agree (`test-Engine.R`, 1e-7 as for the other chains;
   1e-8 relative on log-likelihoods).
3. **Rung 5 at fixed parameters, prior only,** with `sfnpr` as the only topology move: |z| < 3 on
   both grids, with seed replicates for any |z| > 2.
4. **The gamma = 1 one-character topology oracle** (`sfnpr-validation.R`): topology and clade
   frequencies against the exact values, with max |z| consistent with the number of statistics.
