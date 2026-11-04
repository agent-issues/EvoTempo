# Hymenoptera: conventional tip-dating baselines in RevBayes 1.4.1

UCLN + FBD comparator for the CET analyses of the Hymenoptera data (Ronquist et al. 2012; fossil
age intervals and FBD settings of Zhang et al. 2016; inputs in `..`). Three configurations, each
run as six independent chains:

| Config | Data | Clock |
|---|---|---|
| `te-shared` | 343 morphological characters + 4748 sites in 6 DNA partitions | one UCLN clock for all partitions (the EvoTempo comparator; MrBayes links its clock across partitions in the same way) |
| `te-separate` | same | independent UCLN clocks for morphology and DNA |
| `morph` | morphology only | one UCLN clock (like-for-like with CET runs that score morphology only) |

The reported runs use the RevBayes development build at commit `a2a8e2db`
(`rb --version`: `1.4.2-preview`, `/nobackup/pjjg18/revbayes-dev-a2a8e2db` on Hamilton). They do not
use the v1.4.1 release. In v1.4.1, `mvCollapseExpandFossilBranch` is biased when fossil ages are
ranges: expand draws the new parent age above the top of the fossil's age range, but collapse is
allowed from any parent age above the fossil's sampled age. This was found on the penguin data and
fixed upstream after v1.4.1 (`c389174b`). Every fossil here has a stratum interval, so the
maintainer chose the pinned development build. A first set of 22 v1.4.1 chains (about 2.5 days)
was superseded; the development-build chains started from their last states. The cluster's
`revbayes/1.2.5` module is not used.

## Files

| File | Content |
|---|---|
| `prepare-data.R` | Writes `data/`: one NEXUS matrix per morphological partition (unordered K = 2-7, ordered K = 3-5; K = number of observed states, contiguous from 0 in every analysed character) and per DNA partition; `taxa.tsv` (tip-age intervals); `partitions.Rev` (partition names and sizes, taxon sets) |
| `start-tree.R` | Starting tree (see below): `data/start_tree.nwk`, `data/start_tree_placement.csv` |
| `hymenoptera_fbd.Rev` | The model. A driver sets `config`, `rep`, `n_gen`, `sample_every`, `burnin_gen` |
| `run-hamilton.sh` | One chain (Slurm). `RESUME=true` continues it from its checkpoint |
| `launch-hamilton.sh` | One directory per chain under a run root. On `shared` (3-day limit), each chain is followed by resume jobs (`--dependency=afterany`) |
| `summarise-hamilton.sh`, `summarise.R`, `mcc.Rev` | Pooled post-burn-in trees, MCC trees and the CSVs in `summary/` |
| `check-diagnostics.R` | Checks that the coda-free ESS, PSRF and HPD in `summarise.R` reproduce coda |
| `sa-scale-check.Rev` | Prior-only check of RevBayes's tree-scaling move (see **RevBayes up-down move**) |
| `published.csv`, `published-fossils.csv` | Zhang et al. 2016 Table 5 and Table S1 values used for comparison |
| `summary/` | Results (below). Pooled tree files and traces are not committed |

Run from this directory: `Rscript prepare-data.R; Rscript start-tree.R`; then, on Hamilton,
`bash launch-hamilton.sh <run root> <config> 1 6 shared 2` per configuration and, once the chains
have run, `bash summarise-hamilton.sh <run root> <burn-in fraction>`.

## Model and priors

| Component | This analysis (RevBayes) | Zhang et al. 2016 (MrBayes 3.2.6; constant-rate FBD, random sampling, IGR clock) |
|---|---|---|
| Tree prior | `dnFBDP(rootAge, lambda, mu, psi, rho)`, sampled ancestors allowed, `condition = "survival"` | `brlenspr = clock:fossilization`, `samplestrat = random` |
| Net diversification d = λ − μ | Exponential(rate 100) | Exponential(100) (the "intermediate" prior of their Table 4) |
| Turnover r = μ/λ; fossil-sampling proportion s = ψ/(μ + ψ) | Beta(1, 1) each | Beta(1, 1) each |
| Extant sampling ρ | 0.0005, random sampling | 0.0005 |
| Root age | offset exponential, minimum 315 Ma, mean 396 Ma: `dnExponential(1/81, offset = 315)` | `treeagepr = offsetexp(315, 396)` |
| Holometabola (with all fossils) | monophyletic (`dnConstrainedTopology`); age calibration offset exponential, minimum 302 Ma, mean 396 Ma | same constraint and calibration |
| Fossil ages | uniform on the stratum intervals of `../hymenoptera_tips.csv` (`mvFossilTipTimeUniform`) | same intervals |
| Clock | UCLN: base rate ~ LogNormal(−7.1, 0.5) changes per site (character) per Myr, averaged over partitions by site; log-SD ~ Exponential(mean 0.587405); non-centred branch multipliers exp(σz − σ²/2), z ~ N(0, 1) | IGR, `igrvarpr = exp(37)`; `clockratepr = lognorm(−7.1, 0.5)` |
| Partition rates | flat Dirichlet on (multiplier × proportion of sites), so that the site-weighted mean multiplier is 1 (`ratepr = variable`) | same |
| Morphology | Mkv (`coding = "variable"`) + Γ4, one shape ~ Exponential(1); the 43 analysed ordered characters as adjacent-state, equal-rate matrices; 10 characters excluded | Mkv + Γ; same ordered and excluded sets |
| DNA | 12S+16S, 18S, CO1 positions 1-2, EF1α positions 1-2, EF1α position 3: GTR + Γ4; 28S: SYM + Γ4 (equal frequencies). Flat Dirichlet exchangeabilities and frequencies; shapes ~ Exponential(1); CO1 position 3 excluded | same partitions and models (GTR + Γ, SYM + Γ; no +I) |
| MCMC | six independent chains per configuration, no heating | 4 runs × 4 Metropolis-coupled chains, 5 × 10⁷ generations |

### Differences from the source analysis

1. **UCLN instead of IGR.** The CET comparator is the uncorrelated lognormal clock. MrBayes's IGR
   variance prior (`igrvarpr = exp(37)`, a variance per unit time) has no UCLN counterpart, so the
   log-SD prior is the RevBayes clock tutorial's Exponential(mean 0.587405). Under that prior, the
   95 % range of a branch rate spans about one order of magnitude.
2. **The separate-clock and morphology-only runs** have no counterpart in the source. Both base
   rates of `te-separate`, and the base rate of `morph`, take the source's LogNormal(−7.1, 0.5).
3. **Calibration densities.** The root-age prior multiplies the FBD density conditioned on the root
   age, as `treeagepr` does in MrBayes. The Holometabola calibration is a clamped exponential whose
   offset carries the node age, contributing exp(−(t − 302)/94)/94 for t > 302 Ma, as a MrBayes
   `calibrate` statement does.
4. **Partition rate multipliers scale the Γ site rates**, not the branch rates. The product is the
   same. The point is that branch rates built element by element let RevBayes recompute the
   likelihood along one path after a single-branch move; a vector expression marks every branch
   as changed, which made the chains about 4× slower.
5. **No tree-scaling moves** (see below). Node ages move by `mvNodeTimeSlideUniform`, the root by
   a slide on the root age, and fossils by `mvFossilTipTimeUniform`. Topology moves are `mvFNPR`,
   `mvNarrow` and `mvCollapseExpandFossilBranch`. RevBayes has no Metropolis coupling here, so
   mixing rests on independent chains.
6. **Random extant sampling only.** `dnFBDP` has no diversified-sampling option, so the published
   comparison is Zhang et al.'s FBD (Rnd) analysis.

## Starting tree

On the penguin data, chains started from random FBD trees stalled with crown fossils outside a
crown that was too young. Here every chain starts from a tree built from the source study
(`start-tree.R`):

- **Backbone.** `tree_1` of `../raw/MrBayes_Commands.nex`: the 68 extant taxa, in the topology of
  Ronquist et al. 2012's fixed-topology analyses. A node with n extant descendants gets age
  400 × ((n − 1)/67)^0.5 Ma (root 400 Ma, Holometabola 394 Ma, crown Hymenoptera 375 Ma).
- **Fossils.** Each fossil is placed at the midpoint of its interval, oldest first. It joins the
  lineage of its most similar extant hymenopteran (proportion of differing states over the
  characters scored in both) on the first branch, walking rootwards, whose parent is older than
  the fossil.

All 45 fossils therefore start inside an old crown, the opposite of the penguin trap. The chains
move fossils and node ages freely from there.

## RevBayes up-down move

The models use no tree-scaling move. In RevBayes 1.4.1, `TreeUtilities::rescaleSubtree()`
(`src/core/utils/TreeUtilities.cpp:1608-1631`) scales a node and the interior nodes below it. It
stops at a node with a sampled-ancestor child, leaving that node and its subtree unscaled. The
Hastings ratios of the moves that call it count nodes in other ways:

- `mvUpDownScale` counts `Tree::getNumberOfInteriorNodes()`, which leaves out the root
  (`UpDownScaleProposal.cpp:208-256`, `Tree.cpp:836-853`).
- `mvSubtreeScale` and `mvLayeredScaleProposal` use `getNumberOfNodesInSubtree(false) - 1`.
- `mvTreeScale` counts all interior nodes.

`sa-scale-check.Rev` adds one move at a time to a reference move set and compares the
stationary distributions. The tree prior is `dnFBDP` with λ = 1, μ = 0.8, ψ = 0.6, ρ = 1, and there
are no data. x ~ Exponential(1) is independent of the tree and is paired with it in the up-down
variants. Each variant ran 4 chains × 400 000 iterations, sampled every 10, with 10 % burn-in. Cells
give the pooled mean ± MCSE; full table in `summary/scale-move-check.csv`:

| Variant | Tree | Added move | x (expected 1) | root age | tree length |
|---|---|---|---|---|---|
| A | 4 extant + 6 fossils, origin | none (reference) | 1.013 ± 0.011 | 3.489 ± 0.002 | 8.221 ± 0.013 |
| B | as A | up-down (x up, tree down) | stuck at 0 | 3.724 ± 0.002 | 8.881 ± 0.013 |
| C | 10 extant, origin | none (reference) | 0.993 ± 0.011 | 1.644 ± 0.002 | 6.666 ± 0.009 |
| D | as C | up-down (x up, tree down) | 1.584 ± 0.008 | 1.574 ± 0.002 | 6.370 ± 0.008 |
| E | as A | `mvSubtreeScale` | 0.999 ± 0.011 | 3.607 ± 0.002 | 9.319 ± 0.011 |
| F | as C | `mvSubtreeScale` | 0.990 ± 0.011 | 1.639 ± 0.002 | 6.649 ± 0.007 |
| G | as A | `mvLayeredScaleProposal` | 0.991 ± 0.011 | 3.725 ± 0.002 | 9.851 ± 0.013 |
| H | as C | `mvLayeredScaleProposal` | 1.000 ± 0.011 | 1.702 ± 0.002 | 7.188 ± 0.007 |
| I | 10 extant, stochastic root age | none (reference) | 1.016 ± 0.011 | 1.119 ± 0.002 | 5.065 ± 0.008 |
| J | as I | up-down (x up, tree and root age down) | 1.004 ± 0.009 | 1.118 ± 0.002 | 5.060 ± 0.006 |
| K | as C | `mvTreeScale` | 1.006 ± 0.011 | 1.642 ± 0.002 | 6.658 ± 0.006 |

What the variants show:

- **Up-down with a root-age variable is unbiased (J).** Leaving the root out of the count is
  therefore a convention: the move expects the root age to be scaled as a variable of its own.
- **That convention breaks on origin-conditioned trees.** There the root is a node of the tree
  and cannot be added to the move, and the move is biased (D; root age −21 MCSE).
- **Sampled ancestors bias up-down (B) and `mvSubtreeScale` (E)**, because nodes are counted that
  were not scaled. Without them `mvSubtreeScale` (F) and `mvTreeScale` (K) agree with the reference.
- **`mvLayeredScaleProposal` is biased even without sampled ancestors (H).**
- **With sampled ancestors, up-down is biased even when the root age is a variable of the move.**
  Variants L and M repeat I and J on 4 extant taxa and 6 fossils. The root is constrained to
  split e1 from the rest, and `root_time ~ dnExponential(1, offset = 3)`. Values are M minus L,
  with z in MCSE:

  | Build | x | root age | tree length | sampled ancestors |
  |---|---|---|---|---|
  | 1.4.1 | −0.830 (z −60) | +0.210 (z +123) | +0.547 (z +32) | −0.088 (z −16) |
  | development `a2a8e2db` | −0.863 (z −64) | +0.214 (z +126) | +0.562 (z +33) | −0.081 (z −15) |

  Up-down is therefore not usable on any FBD tree with sampled ancestors in either build, and the
  `te-*` chains run without it.
- **Fossils must be kept off the root when the root age is a variable.** In an earlier L/M design
  the topology was unconstrained, and the root ran away to about 13 700 in both builds. The
  cause is a separate fault. With `dnFBDP(rootAge = root_time)` and
  `mvCollapseExpandFossilBranch(tree, root_time)`, collapsing a fossil that is a child of the root
  sets the root node's age to the fossil's age, while `root_time` stays unchanged
  (`CollapseExpandFossilBranchProposal.cpp`, lines 169–192 in `a2a8e2db`). Each later move on
  `root_time` then carries the fossil, now a sampled ancestor at the root, with the root. No
  age-range check stops this, and a root-age prior bounded above the oldest fossil does not
  prevent it. Origin-conditioned trees are not affected.

  The analysis here is safe for two reasons. All 45 fossils lie inside the constrained
  Holometabola clade, and the root lies outside it. We also checked every fifth sampled tree of
  every chain: no tree has a fossil outside its age range or a sampled ancestor at the root.

The penguin baseline's round-1 and round-2 model (`../../penguins/revbayes/penguins_fbd.Rev`)
used up-down on an origin-conditioned FBD tree with sampled ancestors.

## Build check

This is a data-free check of the `morph` model with `mymodel.ignoreAllData()`: 2 chains per build,
200 000 iterations sampled every 100, and 10 % burn-in. Cells give mean ± MCSE, with the two chain
means in brackets. Logs are on Hamilton under
`/nobackup/pjjg18/evotempo-baselines-rb141/prior-check/`.

| Quantity | v1.4.1 | development `a2a8e2db` |
|---|---|---|
| sampled ancestors | 0.039 ± 0.008 [0.015, 0.064] | 0.060 ± 0.004 [0.002, 0.118] |
| sampling proportion s | 0.0010 ± 0.0007 | 0.0029 ± 0.0011 |
| root age (Ma) | 559 ± 163 [329, 788] | 797 ± 8 [328, 1265] |

Neither build mixes under this vague prior: within each build the two chains sit far apart, so
the MCSEs understate the uncertainty and the comparison is inconclusive. Sampled ancestors are rare
in both builds (< 0.12 per tree), because ρ = 0.0005 and d ~ Exponential(100). The move to the
development build follows from the known `mvCollapseExpandFossilBranch` bug with interval ages,
not from this check.

## Results

{{RESULTS}}

## Compute

{{COMPUTE}}
