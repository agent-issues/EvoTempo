# Hymenoptera: conventional tip-dating baselines in RevBayes 1.4.1

UCLN + FBD comparator for the CET analyses of the Hymenoptera data (Ronquist et al. 2012; fossil
age intervals and FBD settings of Zhang et al. 2016; inputs in `..`). Three configurations, each
run as six independent chains:

| Config | Data | Clock |
|---|---|---|
| `te-shared` | 343 morphological characters + 4748 sites in 6 DNA partitions | one UCLN clock for all partitions (the EvoTempo comparator; MrBayes links its clock across partitions in the same way) |
| `te-separate` | same | independent UCLN clocks for morphology and DNA |
| `morph` | morphology only | one UCLN clock (like-for-like with CET runs that score morphology only) |

All runs use RevBayes v1.4.1 (tag `v1.4.1`, commit `3a504eb`), built on Hamilton
(`/nobackup/pjjg18/revbayes-1.4.1`), not the cluster's `revbayes/1.2.5` module.

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

`mvUpDownScale` with a tree (and `mvTreeScale`) is not used. In RevBayes 1.4.1,
`TreeUtilities::rescaleSubtree()` scales the root, but it does not rescale a node with a
sampled-ancestor child, or anything below it. The Hastings ratio counts
`Tree::getNumberOfInteriorNodes()`, which is every interior node except the root.
`sa-scale-check.Rev` compares prior-only chains with and without
`mvUpDownScale(x up, tree down)` for an independent x ~ Exponential(1). There are no data,
400 000 iterations and 10 % burn-in:

| Tree | Moves | x mean (expected 1) | root age | tree length |
|---|---|---|---|---|
| 4 extant + 6 fossils (sampled ancestors possible) | reference | 0.999 | 3.497 | 8.27 |
| | + up-down | 0.000 (collapses) | 3.729 | 8.90 |
| 10 extant, no fossils | reference | 1.018 | 1.634 | 6.64 |
| | + up-down | 1.571 | 1.573 | 6.36 |

Both biases have the sign the miscount predicts. With no sampled ancestors, one node too few is
counted (the root). With sampled ancestors, nodes are counted that the move did not scale. The
penguin baseline (`../../penguins/revbayes/penguins_fbd.Rev`) uses this move.

## Results

{{RESULTS}}

## Compute

{{COMPUTE}}
