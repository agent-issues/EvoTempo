# Trilobites: CET analysis pipeline (morphology only)

Everything needed to run the four (beta, gamma) cells of the covariant evolutionary tempo model
on the Cambrian trilobite data of Paterson, Edgecombe & Lee (2019) (`../README.md`), with the
merged RevBayes FBD + UCLN baseline (`../revbayes/`, #78) as the comparator.  It follows the
penguin pipeline (`../../penguins/cet/`); the differences are listed in each section and in
**Deviations**.  The data have no DNA, so there is a single (morphology-only) stage.

**Production is not launched by this PR.** It starts after the independent review.

## Files

| File | Role |
|---|---|
| `prepare.R` | Converts `../trilobites_morph.nex` and `../trilobites_characters.csv` into `data/characters.csv` and `data/partitions.csv`, and draws one starting state per production chain from the RevBayes posterior (`data/start-tree-01.nwk` ... `16.nwk`, `data/start-values.csv`, `data/start-tree-log.txt`) |
| `common.R` | Time origin and conditioning, priors, moves, data loading, and the three namespace rebindings: character likelihood, per-sample recording and exact resume, ingroup constraint |
| `prune.cpp` | C++ pruning (via `Rcpp::sourceCpp()`) of the columns `MkLikelihood()` cannot score: ambiguous cells, ordered characters, and the constant columns of the Mkv corrections |
| `check-mkv.R`, `revbayes/check-oracle.Rev`, `revbayes/check-oracle.sh` | Fixed-tree checks against RevBayes: the character likelihood against `dnPhyloCTMC(coding = "variable")`, and the conditioning and time origin against `dnFBDP(rho = 0, condition = "sampling")` |
| `check-window.R`, `check-window.sh` | Forward-simulation check of the all-fossil set-up (rho = 0, origin start, survival, beta = 1) |
| `check-priors.R` | Monte Carlo check of the priors matched to the comparator |
| `check-sbc.R`, `check-sbc.sh` | Free-parameter simulation-based calibration at rho = 0, beta = 1 (SBC-lite) |
| `run-chain.R` | One chain of one cell, in checkpointed segments |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R`, `check-pilot.R` | Smoke run on the `test` partition and its checks (timing, acceptance, exact resume, `bayes-factors.R eval` on every smoke chain including two replicate-2 chains); pilot summary for the cost projection |
| `bayes-factors.R` | Cross-model evaluation of recorded states (`eval`) and bridge-sampling Bayes factors (`combine`), also reweighted to the comparator's joint FBD prior |
| `check-labels.R` | Label and `combine` checks on synthetic chains, no MCMC (#100, #101): plain-decimal labels, rung columns, the expected chain set and its opt-out |
| `check-combine.sh`, `check-combine.R` | `eval` and `combine` on tiny real chains with a ladder rung (b0.5g0) present, `combine` stopping on a missing `cross.csv`, and the md5 resume guard (#100-#102) |
| `summarise.R` | Convergence (MkPrime R-hat and ESS), parameters, ages in Ma, clade and sampled-ancestor frequencies; CET and RevBayes side by side |
| `smoke-results/`, `window-check/`, `sbc-check/` | Outputs quoted below |

## Data (`prepare.R`)

- **Characters.** 106 characters: the 107 of the matrix less character 17 (one state in all taxa), as
  the comparator.  K = highest observed state + 1 and the partitions are those of
  `../revbayes/prepare-data.R`: unordered u2 (79), u3 (13), u4 (5), u5 (2), u6 (1), u7 (1), u9 (2),
  u16 (1); ordered o4 (character 47) and o5 (character 45).  `?` and `-` are unknown.  The 84
  polymorphic cells (in 17 characters: 3, 5, 6, 7, 10, 12, 14, 42, 45, 47, 72, 76, 78, 81, 86, 94, 107)
  are kept as **ambiguity** among the listed states, as RevBayes reads them.  (The penguin
  pipeline set its 32 ambiguous cells to unknown, because its comparator was built from the same
  recoded matrix; here the comparator already exists and treats them as ambiguous.)
- **Tip ages.** Point ages: age = Ma - 485.4 (0.6-33.6 Myr).  Each tip's range is its age, so no
  move changes a tip age (`fossilAge` is not used).
- **Constraint.** Ingroup (all but the three fallotaspidoids) monophyletic, as the comparator and
  the source analyses (**Ingroup constraint**).

## Time origin and conditioning

All 107 tips are Cambrian fossils, and the comparator runs the FBD process from the origin to the
end of the Cambrian sampling window (485.4 Ma, GTS 2012) with no sampling there (rho = 0),
conditioned on at least one sample (`condition = "sampling"`, origin mode).  EvoTempo has no offset
argument: age 0 is the end of the process.  Shifting every age by 485.4 Myr therefore *is* the
offset, and the matching model is

    DCETModel(rho = 0, condition = "survival", start = "origin")

`condition = "survival"` divides by the probability that the stem lineage leaves at least one
sample (spec section 6.1, row "Origin"), evaluated at the current origin age; with rho = 0 that
sample can only be a fossil.  `condition = "extant"` (the penguin default) is impossible here (no
extant tip: likelihood zero), and `condition = "none"` would drop the survival term that RevBayes
includes.  Marginal (not joint) conditioning over the origin bin, as the penguin pipeline.

Rung 5 (`data-raw/rung5-age-results/README.md`) validated free origin ages under survival
conditioning only against the package's own `TreeLikelihood()`, with rho > 0 and extant tips; it
has no check of survival conditioning with bin-varying rates against simulation.  So the set-up
was checked here before being relied on:

1. **Exact check at beta = 0, on the real start tree** (`check-mkv.R` with
   `revbayes/check-oracle.Rev`, RevBayes development fd807d10; start tree 01, which has 7 sampled
   ancestors).  With beta = 0 the dCET tree density equals the FBD density for any generator (spec
   section 7.2).  Over six (lambda, mu, psi, origin) points, RevBayes's `dnFBDP` lnProbability minus
   `TreeLikelihood()` under `CetModel(0, 0)` equals the known offset **-ln(107!) + nBif ln 2**, here
   -327.627246176 with nBif = 99 bifurcations (internal nodes without a sampled ancestor, the root
   included under an origin): maximum residual 3.4e-8 (C++ engine) and 6.7e-7 (R engine)
   (`smoke-results/fbd-check.csv`).  The offset is therefore a constant only for a given tree: it
   changes by ln 2 with each sampled ancestor (review of #91, F2; the independent verifier found it on
   all six earlier start trees, with nBif 101-106).  `TempoMCMC()` adds nBif ln 2 to every tree-move
   ratio (`.TreeMove()`), so the tree density it samples is RevBayes's + ln(107!), the same constant
   for every tree.  The samples' `logPost` does **not** include nBif ln 2, so it is not the log of the
   sampled target, and it is not comparable with RevBayes's `Posterior` (F10).  Negative controls fail
   as they should: `condition = "none"` gives a spread of 2.76 nats over the grid, and putting the
   present at the youngest tip (486.0 Ma) instead of 485.4 Ma gives 3.16 nats.
2. **Forward simulation at beta = 1** (`check-window.R`; design and criteria fixed in its header
   before any run): K = 16, z in [-3, 3], rho = 0, origin age ~ U(3, 10), 4 fossil samples, fixed
   parameters, every tree move.  Prior-only `TempoMCMC()` under `survival` and under `none`, each
   scored against forward simulation under both.

   Replicate 1 (Slurm 19429515; 32 chains x 40 000 iterations per conditioning, 8000 simulated
   trees each, 9.9 core-hours; `window-check/window-verdict.csv`, `window-check.csv`): both
   own-target runs pass (max |z| 1.89 for survival, 2.92 for none, bar 3.5).  The negative controls
   show that the test separates the conditionings, but one falls short of the pre-registered power
   bar of 5: `none` chains against `survival` simulations give max |z| 6.16 (originAge q25), and
   `survival` chains against `none` simulations 4.86 (originAge q50).  By the rule fixed beforehand,
   replicate 1 therefore shows no discrepancy for the all-fossil, rho = 0 tree density and moves, but
   lacks full power for the conditioning.  The stem's survival varies only from 0.339 to 0.481 over
   the origin prior, which limits the power.  **Replicate 2** (decided after replicate 1; same design
   and criteria, new seeds (`CET_WINDOW_SEED=100000`), 64 chains and 16 000 simulated trees per
   conditioning; Slurm 19430059, 21.0 allocated core-hours; `window-check/window-verdict-seed100000.csv`,
   `window-check-seed100000.csv`, `window-log-seed100000.txt`): both own-target runs pass (max |z|
   2.76 for survival at meanNodeAge q50, 2.02 for none at treeLength q95; bar 3.5), and both negative
   controls clear the power bar (`none` chains against `survival` simulations 8.67, originAge mean;
   `survival` chains against `none` simulations 7.76, originAge q25).  Under the rule below, which
   was fixed before the result existed, the window check **passes** with full power for the
   conditioning, and it does not block production.  The `.rds` outputs (41 MB of raw draws) stay
   in the Hamilton work directory.  The exact beta = 0 check (item 1) does not depend on this.

   **Decision rule for replicate 2** (written at 15:25 GMT on 11 November (commit 8c6959a), while Slurm 19430059
   was still PENDING and before any of its output existed; review of #91, F5):
   - replicate 2 is judged **on its own**, against the same bars as replicate 1 (own target: every
     |z| < 3.5; power: each negative control max |z| >= 5); it is not pooled with replicate 1;
   - if either own-target run has a statistic with |z| >= 3.5, the window check **fails**, and
     production does not start until the cause is found;
   - if both own-target runs pass but a negative control misses the power bar, the result is
     reported as **limited power** for the survival conditioning, which does not block production.

3. **Free-parameter SBC-lite at rho = 0, beta = 1** (`check-sbc.R`; review of #91, F4; design and
   pass rule fixed in its header before any run).  The window check holds the parameters fixed, so
   the parameter sweep, joint sweep and shift move were untested against simulation in this set-up.
   200 replicates: parameters and origin age from narrow priors, forward simulation with rho = 0 and
   survival, trees of 3-12 fossil samples (a data-only selection, which keeps ranks uniform), then
   `TempoMCMC()` with free parameters and origin from an independent prior draw (2000 burn-in +
   19 800 iterations thinned to 99 draws); rank uniformity for lambda, mu, psi, theta, s and the
   origin age by chi-square over 10 bins, Holm-adjusted p >= 0.01 to pass.  **Result: PASS** (Slurm 19431496; 200
   replicates, 0 failed; 48.7 allocated core-hours; `sbc-check/sbc-summary.csv`, `sbc-ranks.csv`,
   `sbc-log.txt`).  Raw chi-square p: lambda 0.941, mu 0.994, psi 0.485, theta 0.350, s 0.177,
   origin age 0.392 (Holm-adjusted 1 for all; bar 0.01).  Mean normalised ranks 0.494-0.509;
   median ESS of the 99 thinned draws 94.3-98.0 per quantity.  Simulated trees had 5.4 samples on
   average (3-12) and 1.83 sampled ancestors, so this checks the samplers and the free-parameter
   moves at rho = 0, not behaviour on a 107-tip tree; 200 replicates over 10 bins detect only
   gross miscalibration.  Its parameters (lambda about 0.5, mu about 0.4) are also far from the
   trilobite regime (the comparator's lambda about 3-13 with turnover about 0.99), so it says
   nothing about mixing or calibration at high turnover (#91 delta review; #102).  A trial of 8 replicates x 990 iterations (Slurm 19431323) ran end to end.

**Engine findings** (not used by `TempoMCMC()`, whose target is `.LogAugmented()` and whose frozen
messages use eigendecomposed slice exponentials): `FBDTreeLikelihood()` returns NaN with rho = 0
and no extant tip (it adds 0 x log 0); `TreeLikelihood()`'s default `atol = 1e-14` on the rescaled D
is too loose on long branches where D decays far (grid case 6, psi 0.3 and a 100-Myr stem: 9.3
nats too high, the engines disagree), so `check-mkv.R` uses `atol = 1e-30`, at which both engines
agree with RevBayes.  An issue draft is with the lead.

## Character model (`common.R`, `prune.cpp`)

`PartitionedMk()` (rebinding `MkLikelihood()`) scores the comparator's character model without its
Gamma rate variation:

- unordered characters without ambiguous cells: `MkLikelihood()` per K (C++ engine);
- characters with ambiguous cells (15 unordered) and the two ordered characters: `CetPrune()`
  (C++), with an R reference (`SymLogLik()`, `OrdLogLik()`) used by the R engine and by the checks;
- ordered characters: adjacent-state changes at equal rates, rescaled to mean rate 1 (rate
  K / (2(K - 1)), as RevBayes's `fnOrderedRateMatrix(maxState = K)`, whose matrices the oracle
  prints), uniform root;
- the Mkv correction once per missing-data pattern and partition, as RevBayes
  `coding = "variable"` (#65): 81 unordered patterns plus one per ordered character.  A cell is
  observed for the masks unless it is `?`, `-` or ambiguous over **all** K states: RevBayes treats
  such a cell (character 12 has one `{01}` in a binary character) as missing, and counting it as
  observed moved partition u2 by 1.4e-3 nats at rate 0.005 before this was found.

`check-mkv.R` on start tree 01 with a strict clock (`smoke-results/mkv-check.csv`, both engines):

| Rate | CET (per pattern) | RevBayes `coding = "variable"` | Difference | Mkv correction |
|---|---|---|---|---|
| 0.005 | -2805.185423 | -2805.185423 | -4.7e-10 | 31.70 |
| 0.0236 | -2667.526174 | -2667.526174 | 1.1e-10 | 4.36 |
| 0.05 | -3123.850355 | -3123.850355 | -7.1e-11 | 1.98 |

Every partition agrees to 5.1e-10 (10 partitions x 3 rates x 2 engines); the route through the R pruners alone (no C++) agrees
with the C++ route to 1e-8, the engines to 1.5e-11, and power t scales the corrected lnL exactly.
RevBayes prints 14 significant digits (`setOption("outputPrecision", "14")`).

Cost of one evaluation on start tree 01 (`smoke-results/mkv-check-timing.csv`, mean of 200):


| Part | ms |
|---|---|
| `PartitionedMk()`, C++ engine (all of it) | 1.40 |
| of which `MkLikelihood()` groups | 0.34 |
| of which `CetPrune()` (ambiguous, ordered, Mkv columns) | 1.01 |
| `PartitionedMk()`, R engine (R pruners) | 45.0 |

The R pruners alone took 18.4 ms per evaluation (about 3.4 s per MCMC iteration), which is why
`prune.cpp` exists; the R versions remain the reference.

## Priors (`common.R`)

| Parameter | CET prior | RevBayes comparator (`../revbayes/trilobites_fbd.Rev`) |
|---|---|---|
| lambda, mu, psi | independent: the **marginals** of the comparator's prior, truncated to the central 1 - 2 pnorm(-4) of their mass | joint: d = lambda - mu ~ LogNormal(-3.5, 1.5), turnover ~ U(0, 1), sampling proportion ~ U(0, 1) |
| s (tempo SD) | the comparator's UCLN log-SD prior: CV ~ Exp(rate 3), sigma = sqrt(ln(1 + CV^2)), on (0, 3] | ucld_sd from ucld_cv ~ Exp(3) |
| rate (clock at tempo 0) | U(0, 100) | clock_morph (mean branch rate) ~ U(0, 100) |
| theta | LogNormal(log 0.05, 1), truncated at 4 sdlog (as penguins; reversion time median 20 Myr) | - |
| origin | U(0.6, 1000.6) Myr before 485.4 Ma (U(0, 1000) above the youngest tip) | the same |
| rho; conditioning | 0; survival of the stem (at least one sample) | 0; `condition = "sampling"` |
| among-character rate variation | none | Gamma(4), shape ~ Exp(mean 0.5) |

- **lambda, mu, psi.** `TempoMCMC()` takes one independent prior per parameter, evaluated
  parameter by parameter in the parameter sweep, the joint sweep and the shift move, so the
  comparator's joint prior cannot be expressed.  (Injecting the joint-prior factor through
  `.LogAugmented()` would miss the joint sweep, which accepts on `.FrozenLikelihood()`.)  The CET
  priors are therefore the exact marginals that the comparator's prior induces
  (`FbdMarginalPriors()`: log lambda = log d + Exp(1), log mu = log d + Logistic, log psi = log d +
  Logistic + Logistic, by discrete convolution).  They match RevBayes's marginal by marginal
  (`check-priors.R`: |CDF difference| <= 6.4e-4 at nine quantiles, `smoke-results/prior-check.csv`)
  but not jointly: the joint prior has mu < lambda and a dependence between all three, and the
  comparator's posterior turnover (0.992) equals its data-free value (#85, F12), so these rates are
  prior-dominated and the difference can matter.  `bayes-factors.R` and `summarise.R` therefore also
  report every estimate and Bayes factor **importance-reweighted** to the joint prior
  (`JointPriorLogWeight()`; zero weight where mu >= lambda), with the weights' effective sample size.
  The weight ESS is reported twice: MkPrime's bulk ESS of the normalised weights over chains
  (`WeightEss()`, `MkPrime:::.EssBulk()`; autocorrelation-aware) and the iid Kish ESS (review of
  #91, F3).  In the pilot, lambda's own ESS was only 4-12 per 2000 iterations, so whether the
  reweighting has enough ESS can only be seen in production.  **Decision for the maintainer**
  (options under **Open decisions**).
- **"Comparator-matched" holds at unit tempo only** (review of #91, F8).  The priors are on lambda,
  mu and psi at x = 0.  In beta = 1 cells the per-bin rates are lambda e^x and mu e^x while psi is
  not tempo-scaled (delta = 0), so the comparator's (d, turnover, sampling proportion) correspond to
  CET's only at x = 0; in gamma = 1 cells `rate` is the clock at x = 0, not the mean branch rate that
  the comparator's `clock_morph` prior describes.  The matching is exact for cell (0, 0) only.
- **s.** As the penguin pipeline: the tempo SD takes the comparator's UCLN log-SD prior, so tempo
  variation and UCLN rate variation start from the same prior.  The bound 3 is numerical (the
  penguin pilot's s = 86.8 overflow); the mass above it is exp(-270).
- **Starting values** (`CetStart()`): the RevBayes sample whose tree starts chain k (draw k of
  `prepare.R`; **Starting trees**): lambda, mu, psi, rate = clock_morph, s = ucld_sd; theta 0.05.

## Ingroup constraint

EvoTempo's tree moves know no topology constraint.  `common.R` rebinds `.ProposeTreeMove()`: a
proposed tree in which no node has exactly the 104 ingroup tips below it (sampled ancestors count
as tips) is returned as `NULL`, which `.TreeMove()` treats as a rejection.  The target is then the
posterior restricted to constrained trees; the proposal kernel is unchanged, so every allowed move
keeps its Hastings ratio.  Rejections are counted per move (`acceptance.csv`,
`constraintRejections`): in the pilot 14-17% of `slide` proposals; no other move breaks the
constraint.  The wrapper does not change the random-number stream of allowed proposals.

## Starting trees

107 tips exceed ConsTree `Quartet()`'s limit of 100, and the project does not use MCC trees, so each
chain starts from one RevBayes posterior sample.  Rule (`prepare.R`, adapted from the penguin
pipeline's `draw-start-trees.R`, branch penguin-start-tree-r4, 511cd4e; review of #91, F6): pool the
post-burn-in samples of the six RevBayes chains (each without its first 25%, as
`../revbayes/summarise.R`): 18 438 samples.  Draw 16 uniformly at random without replacement
(`set.seed(20261111); sample.int()`).  Draw k starts production chain k = 4 (rep - 1) + cell
(Slurm task k of `run-array.sh`), with that sample's tree and scalar values; a ladder rung or power
chain uses draw <rep>.  Tip ages are reset to the exact point ages (the RevBayes Newick has 6
decimals; the largest tip-implied age spread was 4e-6 Myr), zero-length tips become sampled
ancestors at their exact age, and `origin_time` gives the stem; no age needed raising
(`data/start-tree-log.txt`).  The RevBayes outputs were copied (read-only source,
`/nobackup/pjjg18/evotempo-baselines-rb141/runs-dev/trilobites`, log md5s checked) to
`/nobackup/pjjg18/evotempo-trilobite-cet/rbposterior/`.

The first version of this pipeline started each replicate from its RevBayes chain's
maximum-Posterior sample.  The independent review showed those starts were concentrated (3.0-4.5 sd
above their chain's mean posterior, five of six in the youngest 8% of root ages, clock at the
85th-99th percentile), which weakens R-hat as a test of slow topology mixing.  The random draws
span the posterior (percentiles in the pooled posterior):

| Draw | Chain | RevBayes chain | Iteration | Posterior pct. | Root (Ma) | Root pct. | Sampled ancestors |
|---|---|---|---|---|---|---|---|
| 01 | b0g0 rep 1 | run5 | 584200 | 57.3 | 529.55 | 66.9 | 7 |
| 02 | b0g1 rep 1 | run3 | 567400 | 82.7 | 526.40 | 15.1 | 4 |
| 03 | b1g0 rep 1 | run3 | 644200 | 19.8 | 528.59 | 51.2 | 3 |
| 04 | b1g1 rep 1 | run3 | 319000 | 43.8 | 526.26 | 13.3 | 7 |
| 05 | b0g0 rep 2 | run1 | 219200 | 2.1 | 530.97 | 83.2 | 4 |
| 06 | b0g1 rep 2 | run6 | 760600 | 32.1 | 531.52 | 87.8 | 6 |
| 07 | b1g0 rep 2 | run6 | 686800 | 25.2 | 531.06 | 83.9 | 1 |
| 08 | b1g1 rep 2 | run5 | 268000 | 40.8 | 525.28 | 4.6 | 7 |
| 09 | b0g0 rep 3 | run6 | 563000 | 33.4 | 527.49 | 31.4 | 2 |
| 10 | b0g1 rep 3 | run1 | 257400 | 12.9 | 530.59 | 79.6 | 1 |
| 11 | b1g0 rep 3 | run6 | 491400 | 66.6 | 530.65 | 80.3 | 0 |
| 12 | b1g1 rep 3 | run5 | 678400 | 41.9 | 532.12 | 91.3 | 2 |
| 13 | b0g0 rep 4 | run3 | 558800 | 19.3 | 529.79 | 70.4 | 4 |
| 14 | b0g1 rep 4 | run3 | 461000 | 45.0 | 526.44 | 15.7 | 2 |
| 15 | b1g0 rep 4 | run4 | 249200 | 93.2 | 529.65 | 68.5 | 1 |
| 16 | b1g1 rep 4 | run1 | 278000 | 80.3 | 533.07 | 95.1 | 1 |

Every chain now has its own start tree, and so its own tip order.  Every consumer of the character
matrix takes the tip order of the tree it scores (`LoadCharacters(tipLabel)`), and `run-chain.R`
refuses to resume a checkpoint whose tree has another tip order (**Moves and chain mechanics**).

## Moves and chain mechanics

`CetMoves()` gives every chain the same move counts whatever its start tree (review of #91, F7;
the first version took them from each start tree, so they differed by replicate: nodeAge 26/27, tip
26/27, branch 52/53).  The counts are those for a 107-tip tree without sampled ancestors: with
gamma = 1, 27 `node`, 27 `tip`, 53 `branch`, 1 `shift`; with gamma = 0, `TempoMoves()`'s fixed
schedule (1 `tree`, 2 `node`, 2 `tip`, 2 `branch`, 1 `joint`, 1 `shift`); plus `nodeAge` 27 (0.25 x
106 bifurcations), `slide` 16, `saJump` 11 (0.1 x tips) and `originAge` 1; no `fossilAge` (point ages).
Per-sample trees, states and exact resume work as in the penguin pipeline (rebound `.TreeStats()`).
The checkpoint records the grid (`CET_NBIN`, `CET_ZMAX`), thinning, segment length, beta, gamma,
power, Mkv variant, start-draw index and the start tree (its md5, its path relative to `cet/`,
and the full path); a run with any of them different, from a start tree with another md5, or whose
tree has another tip order than the loaded characters, stops (#76; review of #91, F1).  The resume
guard compares the md5, not the absolute path, so resubmitting from another checkout of the same
commit continues the chain (#102); checkpoints written before #102 have no md5 and are compared by
absolute path.  `bayes-factors.R eval` (review of #91, F1 and A1-A4):
- scores the characters in the chain's own tip order (`states[[1]]$tree$tipLabel`, the same for
  every state of a chain, asserted), not in start tree 01's;
- refuses a chain without a checkpoint, without a recorded variant or with another grid;
- stops when a recomputed log posterior, or a recomputed character log-likelihood, differs from the
  chain's by more than `CET_BF_TOL` (1e-6); the second check also covers power-0 ladder chains,
  whose log posterior does not contain the characters (A4);
- writes `cross.csv` only after the checks pass (temporary file, then rename; a failed eval removes
  any old `cross.csv`), and `combine` refuses any `cross.csv` whose recorded checks exceed the
  tolerance or that has no `checkChar` column (A2);
- names its columns with `FormatValue()` (plain decimal, one number at a time): `format()` on the
  vector c(0, 0.5, 1) gave `aug_b0.0`, `aug_b0.5`, `aug_b1.0`, which `combine` could not find, so
  any ladder rung made it crash (#100; the penguin pipeline's #75 fix).  Chain labels use the same
  function, and `ParseLabel()` / `Chains()` (now in `common.R`) stop on a label they cannot read.

`combine` checks the chains it finds against the expected set (#101): the launch record
`<run root>/launch.csv` (one row per launched group of chains: `cell`, `reps`, optional `power`
and `mkv`; see `ExpectedChains()` in `common.R`), or, if there is none, the production launch of
4 cells x 4 replicates at power 1 (`CET_BF_REPS` changes the 4).  It stops, and reports no Bayes
factor, if an expected chain has no complete segment, if an expected power-1 chain has no
`cross.csv` (its eval failed or was not run), or if it finds a chain that the record does not
list.  Before, a chain without `cross.csv` dropped out silently (#91 delta review, N2: b0g0 -> b0g1
moved from 7.13 to 7.54 nats, exit 0).  The only opt-out is explicit and per chain:
`CET_BF_ALLOW_MISSING=b0g0_run3,...`, which `combine` echoes as "Left out on request".  The
character-power stepping-stone now uses only chains of the `CET_MKV` variant (as penguin #75).
Before these fixes, eval on a replicate-2 chain scored tree-01 character rows against another tip
order (3 of 107 positions shared): the self-check stopped it, but only after it had written a
`cross.csv` that `combine` would have used (gamma edges wrong by up to 114 nats per state).
Log names carry the variant and power (#77).  `CET_MKV` takes `pattern` (default) or `none`; the
penguin pipeline's `allTips` variant is not carried over.

## Smoke run, pilot and cost projection

Hamilton, EvoTempo main 48e6bf7 (`lib-48e6bf7`), one core per chain.

**Smoke** (`smoke.sh`, `test` partition, 40 burn-in + 2 x 100 iterations per chain; rerun for the
fixes of the #91 review as Slurm 19431320, checked by 19431322; `smoke-results/timing.csv`,
`acceptance.csv`, `resume.txt`, `eval.txt`, `check-smoke.txt`): every cell ran from its own draw
(01-04), plus cells (1, 1) and (0, 0) of replicate 2 (draws 08 and 05); 0.47-0.62 s per iteration
after burn-in.  **Exact resume:** cell (1, 1) in two segments of 100 and in one of 200 gave identical
samples (maximum absolute difference 0 over 10 samples x 23 columns) and identical trees.
**`eval` on every smoke chain**, the two replicate-2 chains included: exit 0 and `cross.csv` written
for all six; maximum |recomputed - logPost| 4.6e-12 to 1.6e-9, maximum |recomputed - charLL|
<= 5e-12 (review of #91, A3).  **Guard tests** (Slurm 19431540, `smoke-results/guard-tests.txt`):
`eval` with `CET_BF_TOL = 1e-13` on a replicate-2 chain (residual 1.6e-9) exits 1 and leaves no
`cross.csv`; `combine` runs on the smoke chains and refuses them once one `cross.csv`'s check column
is corrupted.  Peak memory in the first smoke 303-304 MB per chain.

**Pilot** (`smoke.sh` with 500 burn-in + 3000 iterations, thin 10, one chain per cell, `shared`;
Slurm 19429527; `smoke-results/pilot-*`).  The pilot predates the #91 review: each chain started
from the old maximum-Posterior start tree 1, with that tree's move counts (which differ from the
fixed counts by at most one per move), so its rates are indicative only.  ESS (MkPrime) of the 2000
iterations after the first segment:

| Cell | s/iteration (segments 2-3) | Slowest quantity: ESS / 2000 it. | lambda | s | theta | rootAge | charLL |
|---|---|---|---|---|---|---|---|
| (0, 0) | 0.606 | charLL 3.6 | 8.6 | 219.2 | 191.0 | 10.2 | 3.6 |
| (0, 1) | 0.751 | logPost 3.9 | 4.5 | 4.7 | 8.6 | 8.3 | 14.9 |
| (1, 0) | 0.556 | treeLength 3.6 | 8.7 | 159.4 | 159.0 | 4.6 | 4.1 |
| (1, 1) | 0.721 | rootAge 5.0 | 11.8 | 7.0 | 5.3 | 5.0 | 8.5 |

Pilot acceptance (last 1000 iterations): `slide` 1.8-2.1%, `saJump` 2.7-5.5%, `nodeAge` 42-56%,
`originAge` 22-32%, local path moves 60-100%, `shift` 11-39%, parameter sweep 25-31%; with
gamma = 0, whole-tree draws 100% and `joint` 26-28%.  As in the penguin pilot, `slide` is the only
topology move and the likely bottleneck.  Pilot estimates are single short chains from one start
and are **not results**; for orientation only: root (first split) means 532.1 / 527.3 / 531.0 /
522.4 Ma in cells (0,0) / (0,1) / (1,0) / (1,1), against 528.5 Ma in the comparator.

**Projection** (provisional: single 3000-iteration chains understate slow drifts).  The slowest
quantities reach 1.8-2.5 ESS per 1000 iterations, as in the penguin pilot, so the penguin length is
kept: 200 000 iterations per chain (thin 100, segments of 10 000), about 25 000 iterations for a
pooled ESS of 200 with an 8x margin.  The precision target should come from the comparison to be
resolved (user decision, 2026-11-11), so the length is to be revisited after the first segments.

| Item | Core-hours |
|---|---|
| Production: 4 cells x 4 chains x 200 000 iterations at 0.556-0.751 s (cell means after burn-in; mean 0.658 s) | ~585 (31-42 h per chain, one 3-day job each) |
| `bayes-factors.R eval` | < 5 |
| Ladder rungs, if a Bayes-factor edge needs them: 4 chains x ~37 h per rung, up to 3 per edge, 4 edges | 0-1760 |
| Character-power stepping-stone, best cell only (16 powers x 2 chains x 50 000 iterations) | ~290 |
| Checks run for this PR (oracle, window checks, smoke, pilot, prior check) | 12.8 in the first round (Slurm CPU time of every job); the fix round adds about 1 (smoke, oracle, guards) plus the SBC-lite (48.7) and window replicate 2 (21.0), both as allocated core-hours (elapsed x 16 cores) |

Disk: about 20 KB per stored state (`seg*_states.rds`), so about 40 MB per chain and 0.65 GB for 16
chains; remove them once `eval` has run.

## Running on Hamilton

    W=/nobackup/pjjg18/evotempo-trilobite-cet     # src/ = checkout, lib-<sha>/ = EvoTempo build
    cd $W/src
    Rscript data-raw/empirical/trilobites/cet/prepare.R $W/rbposterior
    sbatch -J evotempo-trilobite-cet-oracle data-raw/empirical/trilobites/cet/revbayes/check-oracle.sh
    sbatch -J evotempo-trilobite-cet-smoke --array=1-7 data-raw/empirical/trilobites/cet/smoke.sh
    Rscript data-raw/empirical/trilobites/cet/check-smoke.R $W     # also runs eval on every smoke chain
    sbatch -J evotempo-trilobite-cet --array=1-16 \
      --export=ALL,CET_NITER=200000,CET_THIN=100,CET_SEG=10000 data-raw/empirical/trilobites/cet/run-array.sh
    # launch record (only needed beyond the 4 x 4 production chains; one row per launched group):
    printf 'cell,reps\nb0g0,4\nb0g1,4\nb1g0,4\nb1g1,4\n' > $W/runs/launch.csv   # + e.g. b0.5g0,4 per rung
    # per chain, then once:
    Rscript data-raw/empirical/trilobites/cet/bayes-factors.R eval $W/runs/b1g1_run1   # must exit 0
    Rscript data-raw/empirical/trilobites/cet/bayes-factors.R combine $W/runs        # refuses failed evals and missing chains
    Rscript data-raw/empirical/trilobites/cet/summarise.R $W/runs $W/rbposterior $W/runs/summary

`common.R` compiles `prune.cpp` at start (about 20 s per job).  The Bayes-factor procedure is the
penguin pipeline's (its README, "Bayes factors"); it ran end to end on the smoke chains
(`smoke-results/eval.txt`, `guard-tests.txt`); those numbers carry no information and are not
reported.  `summarise.R` ran on the pilot
and the comparator (`smoke-results/pilot-summary-*`): for the comparator it gives root 528.76 Ma
(mean; 95% interval 524.82-534.09), against the median 528.53 (HPD 524.53-533.57) in
`../revbayes/README.md`, with R-hat <= 1.001 and ESS >= 4370 from MkPrime.  (Its `headline` printed
R-hat as -Inf for the single-chain pilot cells; fixed since, display only.).

## Deviations from the comparator and the penguin pipeline

1. **No among-character rate variation** in CET; the comparator has Gamma(4) with posterior shape
   0.52 (0.37-0.70), i.e. strong rate heterogeneity.  CET's tempo varies along branches, not among
   characters, so CET-versus-RevBayes differences confound clock model with rate heterogeneity
   (engine gap #47).
2. **Independent marginal priors** on lambda, mu, psi instead of the joint (d, turnover, sampling
   proportion) prior, with importance reweighting to the joint prior (**Priors**).
3. **Ingroup constraint** enforced by rejecting proposals (an equivalent target, not a change of
   model), with one edge case not compared with RevBayes: a fallotaspidoid sampled ancestor on the
   ingroup's stem leaves the ingroup monophyletic below it, so it is allowed here; how RevBayes's
   clade check treats that state was not examined.  None of the three fallotaspidoids is a sampled
   ancestor in the 16 starting samples (`data/start-tree-log.txt`).
4. **A sampled ancestor at the root** (a fossil on the stem ancestral to all other samples) is
   excluded by EvoTempo's samplers (rung 5 README) but allowed by RevBayes in origin mode; the
   comparator's root HPD (lower bound 524.5 Ma, older than the oldest tip, 519.0 Ma) shows it carries
   no posterior mass there (#85, F10).
5. Starting trees: one random RevBayes posterior sample per chain (rule above), not an MCC tree.
6. Polymorphic cells kept ambiguous (the penguin pipeline set them to unknown).
7. No summary trees in `summarise.R` (ConsTree not installed on Hamilton; Quartet <= 100 tips).

## Open decisions for the maintainer

- **Priors on lambda, mu, psi**: (a) run as is (marginals) and report the joint-prior reweighting
  where its ESS allows; (b) add joint or transformed priors to `TempoMCMC()` (#48), then run with
  the comparator's prior exactly; or (c) run a matched RevBayes comparator with CET's independent
  priors (and without Gamma), as the penguin pipeline did, at about 4 x 7 days of `long`
  (all-fossil trees cannot be resumed from RevBayes checkpoints).
- **Gamma rate variation**: compare CET with the existing comparator (confounded), or add a
  matched no-Gamma comparator (as (c)), or wait for #47.
- **Run length**: 200 000 iterations per chain as projected, or a target derived from the
  comparison to be resolved.
- Production also waits for Gamma rate variation in the package (#97), per the maintainer.
