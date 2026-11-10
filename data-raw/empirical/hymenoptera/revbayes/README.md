# Hymenoptera: conventional tip-dating baselines in RevBayes

UCLN + FBD comparator for the CET analyses of the Hymenoptera data (Ronquist et al. 2012; fossil
age intervals and FBD settings of Zhang et al. 2016; inputs in `..`). Three configurations, run as
independent chains (8 for each `te-*` configuration, 6 for `morph`):

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
`bash launch-hamilton.sh <run root> <config> 1 <n> shared 2` per configuration (n = 8 for `te-*`, 6 for `morph`) and, once the chains
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
| MCMC | 8 (`te-*`) or 6 (`morph`) independent chains, no heating | 4 runs × 4 Metropolis-coupled chains, 5 × 10⁷ generations |

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

All chains ran on Hamilton under the development build `a2a8e2db`, on `shared` in 3-day jobs, each
resumed from its checkpoint. They started from the last states of the superseded v1.4.1 chains.

| Config | Chains | Iterations per chain | Sampled every | Samples per chain after burn-in (shortest chain) |
|---|---|---|---|---|
| `te-shared` | 8 | 247 150–319 950 | 50 | 3708 |
| `te-separate` | 8 | 273 700–341 350 | 50 | 4106 |
| `morph` | 6 | 250 000–284 000 | 100 | 1875 |

`summarise.R` discards the first 25 % of each chain. The `te-*` chains were stopped once all three
criteria held. All files are in `summary/`: `headline_<config>.txt`, `convergence_<config>.csv`,
`parameters_<config>.csv`, `published_<config>.csv`, `fossils_<config>.csv` (per fossil: P(sampled
ancestor), P(in crown Hymenoptera), and Zhang et al.'s Table S1 values), `clades_<config>.csv` and
`mcc_<config>.nex`. MCC node ages are in Ma, because the extant tips are at 0.

### Convergence

| Config | ASDSF (max SDSF) | Lowest pooled ESS | Highest PSRF | Lowest per-chain ESS |
|---|---|---|---|---|
| `te-shared` | 0.0119 (0.037), 182 clades | 251 (crown Hymenoptera age) | 1.030 (Hymenoptera total-group age) | 20 (crown age) |
| `te-separate` | 0.0139 (0.067), 179 clades | 289 (crown age) | 1.017 (crown age) | 21 (crown age) |
| `morph` | 0.0109 (0.033), 145 clades | 1096 (Prior) | 1.008 (Holometabola age) | 154 (Prior) |

All three configurations meet the criteria: ASDSF < 0.05, PSRF < 1.05 and pooled ESS ≥ 200 for
every logged parameter. In the `te-*` runs the node ages and the clock rate mix slowly along
the clock-rate–age ridge, because no tree-scaling move is available (see **RevBayes up-down move**).
Their per-chain ESS is only about 20–60 for the root, Holometabola and Hymenoptera ages and the
clock rates. The ages and clock rates therefore rest on agreement among the 8 chains (PSRF ≤ 1.03),
not on long single chains. Their HPD bounds are less precise than the pooled ESS suggests. The
other parameters mix better: λ and μ have pooled ESS > 1000, and the rest > 2500.

### Estimates

Median (95 % HPD).

| Parameter | `te-shared` | `te-separate` | `morph` |
|---|---|---|---|
| root (Ma) | 412.1 (340.5–497.1) | 388.4 (324.7–467.5) | 357.1 (315.1–427.9) |
| Holometabola (Ma) | 397.4 (330.3–473.1) | 374.2 (313.0–444.2) | 341.1 (302.1–402.7) |
| crown Hymenoptera (Ma) | 351.4 (299.8–414.3) | 325.3 (278.4–383.0) | 312.4 (264.3–380.5) |
| net diversification λ − μ (per Myr) | 0.0076 (0.0045–0.0109) | 0.0078 (0.0045–0.0112) | 0.0055 (0.0015–0.0099) |
| turnover μ/λ | 0.9993 (0.9986–0.9998) | 0.9994 (0.9988–0.9998) | 0.9999 (0.9997–1.0000) |
| sampling proportion s | 2.6e−7 (0.8e−7–5.5e−7) | 2.4e−7 (0.8e−7–5.0e−7) | 5.1e−8 (1.5e−8–1.1e−7) |
| clock rate (per site per Myr) | 9.55e−4 (7.52e−4–1.17e−3) | DNA 9.95e−4 (7.97e−4–1.21e−3); morphology 1.23e−3 (0.88e−3–1.64e−3) | 1.92e−3 (1.32e−3–2.64e−3) |
| UCLN log-SD | 0.566 (0.462–0.675) | DNA 0.525 (0.415–0.651); morphology 1.063 (0.881–1.253) | 1.051 (0.862–1.270) |
| sampled ancestors (mean) | 0.033 | 0.030 | 0.013 |

The sampled-ancestor count has median 0 and 95 % HPD 0–0 in every configuration. No fossil has
P(sampled ancestor) above 0.0032. Across fossils, these probabilities rank like Zhang et al.'s
Table S1 values for FBD (Rnd): Spearman ρ = 0.95 (`te-shared`), 0.90 (`te-separate`) and 0.75
(`morph`). Their probabilities are larger (maximum 0.016, *Eoxyela*). Every fossil falls in crown
Hymenoptera with probability ≥ 0.97 (`te-shared`, `morph`) or ≥ 0.85 (`te-separate`).

### Comparison with Zhang et al. 2016

Zhang et al. 2016, Table 5 (MrBayes, IGR clock, total evidence), median (95 % HPD) in Ma:

| Analysis | Root | Crown Hymenoptera | Proportion of fossils that are ancestors |
|---|---|---|---|
| **This analysis, `te-shared`** (FBD, random sampling, UCLN) | **412.1 (340.5–497.1)** | **351.4 (299.8–414.3)** | **0.00 (0.00–0.00)** |
| FBD (Rnd): the matching analysis | 406.4 (328.1–500.3) | 346.6 (291.9–426.5) | 0.00 (0.00–0.02) |
| FBD (Div) | 322.7 (315.1–342.6) | 279.4 (255.1–304.3) | 0.02 (0.00–0.07) |
| pcFBD (Rnd) | 369.6 (315.6–433.3) | 325.6 (282.5–374.9) | 0.33 (0.20–0.47) |
| pcFBD (Div) | 328.2 (315.0–358.9) | 251.7 (237.9–310.9) | 0.02 (0.00–0.07) |
| Uniform (U) | 345.1 (315.0–402.2) | 306.0 (289.3–341.3) | – |

`te-shared` matches the analysis it reproduces, FBD (Rnd). The medians differ by about 5 Myr for both
ages, the HPDs nearly coincide, and the proportion of ancestral fossils is 0 in both. The UCLN
clock therefore gives the same answer as their IGR clock. The other two configurations have no
published counterpart. With a separate morphological clock, the ages are about 25 Myr younger
than `te-shared`. With morphology alone, they are younger again, and the root's HPD
reaches the root-age prior's minimum (315 Ma).


## Compute

Core-hours on Hamilton (Slurm `CPUTimeRAW`, one core per chain job):

| Jobs | Core-hours |
|---|---|
| Reported chains, development build: `te-shared` 8 chains (`evotempo-hyme-dev-te-shared-*`) | 1519.9 |
| Reported chains, development build: `te-separate` 8 chains | 1519.9 |
| Reported chains, development build: `morph` 6 chains | 432.0 |
| Superseded v1.4.1 chains: `te-shared` 632.2, `te-separate` 632.2, `morph` 432.0 (`evotempo-hyme-rb141-*`) | 1696.4 |
| Build check, both builds (`evotempo-hyme-prior-check`) | 32.6 |
| Scale-move checks A–M, both builds (`evotempo-revbayes-updown-check*`), and fossil-age check of all chains | 3.6 |
| Smoke tests | 0.9 |
| **Total** | **5205.5** |

The summary jobs ran on 4 cores, about 2 core-hours in all, shared with the trilobite baseline.
