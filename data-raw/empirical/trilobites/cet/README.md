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
| `prepare.R` | Converts `../trilobites_morph.nex` and `../trilobites_characters.csv` into `data/characters.csv` and `data/partitions.csv`, and the RevBayes posterior into six starting states (`data/start-tree-<k>.nwk`, `data/start-values.csv`, `data/start-tree-log.txt`) |
| `common.R` | Time origin and conditioning, priors, moves, data loading, and the three namespace rebindings: character likelihood, per-sample recording and exact resume, ingroup constraint |
| `prune.cpp` | C++ pruning (via `Rcpp::sourceCpp()`) of the columns `MkLikelihood()` cannot score: ambiguous cells, ordered characters, and the constant columns of the Mkv corrections |
| `check-mkv.R`, `revbayes/check-oracle.Rev`, `revbayes/check-oracle.sh` | Fixed-tree checks against RevBayes: the character likelihood against `dnPhyloCTMC(coding = "variable")`, and the conditioning and time origin against `dnFBDP(rho = 0, condition = "sampling")` |
| `check-window.R`, `check-window.sh` | Forward-simulation check of the all-fossil set-up (rho = 0, origin start, survival, beta = 1) |
| `check-priors.R` | Monte Carlo check of the priors matched to the comparator |
| `run-chain.R` | One chain of one cell, in checkpointed segments |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R`, `check-pilot.R` | Smoke run on the `test` partition and its checks (timing, acceptance, exact resume); pilot summary for the cost projection |
| `bayes-factors.R` | Cross-model evaluation of recorded states (`eval`) and bridge-sampling Bayes factors (`combine`), also reweighted to the comparator's joint FBD prior |
| `summarise.R` | Convergence (MkPrime R-hat and ESS), parameters, ages in Ma, clade and sampled-ancestor frequencies; CET and RevBayes side by side |
| `smoke-results/`, `window-check/` | Outputs quoted below |

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
   `revbayes/check-oracle.Rev`, RevBayes development fd807d10).  With beta = 0 the dCET tree density
   equals the FBD density for any generator (spec section 7.2).  Over six (lambda, mu, psi, origin)
   points, `TreeLikelihood()` under `CetModel(0, 0)` differs from RevBayes's `dnFBDP` lnProbability by
   a constant (-324.854657465; the two use different labelling conventions), constant to 1.8e-8
   (C++ engine) and 3.4e-7 (R engine) (`smoke-results/fbd-check.csv`).  Negative controls fail as
   they should: `condition = "none"` gives a spread of 2.76 nats, and putting the present at the
   youngest tip (486.0 Ma) instead of 485.4 Ma gives 3.25 nats.
2. **Forward simulation at beta = 1** (`check-window.R`; design and criteria fixed in its header
   before any run): K = 16, z in [-3, 3], rho = 0, origin age ~ U(3, 10), 4 fossil samples, fixed
   parameters, every tree move.  Prior-only `TempoMCMC()` under `survival` and under `none`, each
   scored against forward simulation under both.  Replicate 1 (Slurm 19429515; 32 chains x 40 000 iterations per conditioning, 8000 simulated
   trees each, 9.9 core-hours; `window-check/window-verdict.csv`, `window-check.csv`): both
   own-target runs pass (max |z| 1.89 for survival, 2.92 for none, bar 3.5).  The negative controls
   show that the test separates the conditionings, but one falls short of the pre-registered power
   bar of 5: `none` chains against `survival` simulations give max |z| 6.16 (originAge q25), and
   `survival` chains against `none` simulations 4.86 (originAge q50).  By the rule fixed beforehand,
   replicate 1 therefore shows no discrepancy for the all-fossil, rho = 0 tree density and moves, but
   lacks full power for the conditioning.  The stem's survival varies only from 0.339 to 0.481 over
   the origin prior, which limits the power.  **Replicate 2** (decided after replicate 1; same design
   and criteria, new seeds, 64 chains and 16 000 simulated trees per conditioning; Slurm 19430059)
   is queued: the cluster is full until about 21:00 on 11 November, and its result will be added
   to this PR as a separate commit.  The exact beta = 0 check (item 1) does not depend on it.

   **Decision rule for replicate 2** (written at 15:40 GMT on 11 November, while Slurm 19430059
   was still PENDING and before any of its output existed; review of #91, F5):
   - replicate 2 is judged **on its own**, against the same bars as replicate 1 (own target: every
     |z| < 3.5; power: each negative control max |z| >= 5); it is not pooled with replicate 1;
   - if either own-target run has a statistic with |z| >= 3.5, the window check **fails**, and
     production does not start until the cause is found;
   - if both own-target runs pass but a negative control misses the power bar, the result is
     reported as **limited power** for the survival conditioning, which does not block production.

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

`check-mkv.R` on start tree 1 with a strict clock (`smoke-results/mkv-check.csv`, both engines):

| Rate | CET (per pattern) | RevBayes `coding = "variable"` | Difference | Mkv correction |
|---|---|---|---|---|
| 0.005 | -2936.746640 | -2936.746640 | 4.9e-10 | 41.50 |
| 0.0236 | -2666.195244 | -2666.195244 | 3.6e-11 | 6.11 |
| 0.05 | -2960.208140 | -2960.208140 | 7.2e-11 | 2.66 |

Every partition agrees to 1.1e-9 (10 partitions x 3 rates x 2 engines); the route through the R pruners alone (no C++) agrees
with the C++ route to 1e-8, the engines to 2.2e-11, and power t scales the corrected lnL exactly.
RevBayes prints 14 significant digits (`setOption("outputPrecision", "14")`).

Cost of one evaluation on start tree 1 (`smoke-results/mkv-check-timing.csv`, mean of 200):


| Part | ms |
|---|---|
| `PartitionedMk()`, C++ engine (all of it) | 1.40 |
| of which `MkLikelihood()` groups | 0.35 |
| of which `CetPrune()` (ambiguous, ordered, Mkv columns) | 1.01 |
| `PartitionedMk()`, R engine (R pruners) | 45.6 |

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
  Whether that reweighting has enough ESS can only be seen in production.  **Decision for the
  maintainer** (options under **Open decisions**).
- **s.** As the penguin pipeline: the tempo SD takes the comparator's UCLN log-SD prior, so tempo
  variation and UCLN rate variation start from the same prior.  The bound 3 is numerical (the
  penguin pilot's s = 86.8 overflow); the mass above it is exp(-270).
- **Starting values** (`CetStart()`): the RevBayes sample whose tree starts the replicate:
  lambda, mu, psi, rate = clock_morph, s = ucld_sd; theta 0.05.

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
replicate starts from one RevBayes posterior sample.  Rule (`prepare.R`): in RevBayes chain k
(k = 1..6), drop the first 25% of samples (as `../revbayes/summarise.R`) and take the sample with
the highest logged Posterior; CET replicate k starts from it (replicates above 6 reuse them).  Tip
ages are reset to the exact point ages (the RevBayes Newick has 6 decimals; the largest
tip-implied age spread was 4e-6 Myr), zero-length tips become sampled ancestors at their exact age,
and `origin_time` gives the stem.  No age needed raising.  The RevBayes outputs were copied
(read-only source, `/nobackup/pjjg18/evotempo-baselines-rb141/runs-dev/trilobites`, log md5s
checked) to `/nobackup/pjjg18/evotempo-trilobite-cet/rbposterior/`.

| k | Iteration | RevBayes Posterior | Root (Ma) | Origin (Ma) | Sampled ancestors |
|---|---|---|---|---|---|
| 1 | 484800 | -3477.770 | 525.77 | 527.70 | 3 |
| 2 | 427200 | -3474.515 | 529.52 | 529.78 | 2 |
| 3 | 639400 | -3484.564 | 523.96 | 529.03 | 0 |
| 4 | 752600 | -3481.456 | 525.66 | 525.69 | 1 |
| 5 | 221200 | -3450.545 | 523.51 | 523.61 | 5 |
| 6 | 646000 | -3474.525 | 525.41 | 528.79 | 0 |

Using a different RevBayes chain per replicate (rather than one tree for all, as the penguin
pipeline's MCC tree) gives dispersed starts for R-hat.

## Moves and chain mechanics

`CetMoves()` = `TempoMoves()` (gamma = 1: 27 `node`, 26 `tip`, 53 `branch`, 1 `shift`; gamma = 0:
1 `tree`, 2 `node`, 2 `tip`, 2 `branch`, 1 `joint`, 1 `shift`) plus `nodeAge` 26 (0.25 x
bifurcations), `slide` 16, `saJump` 11 (0.1 x tips), `originAge` 1; no `fossilAge` (point ages).
Per-sample trees, states and exact resume work as in the penguin pipeline (rebound `.TreeStats()`).
The checkpoint now records the grid (`CET_NBIN`, `CET_ZMAX`), thinning and segment length as well
as beta, gamma, power and Mkv variant, and a run with any of them different stops (#76); `eval`
refuses a chain without a checkpoint, without a recorded variant or with another grid, and stops
when a recomputed log posterior differs from the chain's by more than `CET_BF_TOL` (1e-6) (#76).
Log names carry the variant and power (#77).  `CET_MKV` takes `pattern` (default) or `none`; the
penguin pipeline's `allTips` variant is not carried over.

## Smoke run, pilot and cost projection

Hamilton, EvoTempo main 48e6bf7 (`lib-48e6bf7`), one core per chain.

**Smoke** (`smoke.sh`, `test` partition, 40 burn-in + 2 x 100 iterations per cell; Slurm 19429516;
`smoke-results/timing.csv`, `acceptance.csv`, `resume.txt`): every cell ran; 0.45-0.63 s per
iteration after burn-in.  **Exact resume:** cell (1, 1) in two segments of 100 and in one of 200
gave identical samples (maximum absolute difference 0 over 10 samples x 23 columns) and identical
trees.  Peak memory 303-304 MB per chain.

**Pilot** (`smoke.sh` with 500 burn-in + 3000 iterations, thin 10, one chain per cell, `shared`;
Slurm 19429527; `smoke-results/pilot-*`).  ESS (MkPrime) of the 2000 iterations after the first
segment:

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
| Checks run for this PR (oracle, window checks, smoke, pilot, prior check) | 12.8 (Slurm CPU time of every job of this PR; replicate 2 of the window check adds about 20) |

Disk: about 20 KB per stored state (`seg*_states.rds`), so about 40 MB per chain and 0.65 GB for 16
chains; remove them once `eval` has run.

## Running on Hamilton

    W=/nobackup/pjjg18/evotempo-trilobite-cet     # src/ = checkout, lib-<sha>/ = EvoTempo build
    cd $W/src
    Rscript data-raw/empirical/trilobites/cet/prepare.R $W/rbposterior
    sbatch -J evotempo-trilobite-cet-oracle data-raw/empirical/trilobites/cet/revbayes/check-oracle.sh
    sbatch -J evotempo-trilobite-cet-smoke --array=1-5 data-raw/empirical/trilobites/cet/smoke.sh
    Rscript data-raw/empirical/trilobites/cet/check-smoke.R $W
    sbatch -J evotempo-trilobite-cet --array=1-16 \
      --export=ALL,CET_NITER=200000,CET_THIN=100,CET_SEG=10000 data-raw/empirical/trilobites/cet/run-array.sh
    # per chain, then once:
    Rscript data-raw/empirical/trilobites/cet/bayes-factors.R eval $W/runs/b1g1_run1
    Rscript data-raw/empirical/trilobites/cet/bayes-factors.R combine $W/runs
    Rscript data-raw/empirical/trilobites/cet/summarise.R $W/runs $W/rbposterior $W/runs/summary

`common.R` compiles `prune.cpp` at start (about 20 s per job).  The Bayes-factor procedure is the
penguin pipeline's (its README, "Bayes factors"); it ran end to end on the smoke chains
(`smoke-results/bf-eval.txt`, `bf-combine.txt`: recomputed log posteriors match the chains' to <= 2.1e-9); those numbers carry no information and are not reported.  `summarise.R` ran on the pilot
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
   ancestor in the six starting samples (`data/start-tree-log.txt`).
4. **A sampled ancestor at the root** (a fossil on the stem ancestral to all other samples) is
   excluded by EvoTempo's samplers (rung 5 README) but allowed by RevBayes in origin mode; the
   comparator's root HPD (lower bound 524.5 Ma, older than the oldest tip, 519.0 Ma) shows it carries
   no posterior mass there (#85, F10).
5. Starting trees: one RevBayes sample per replicate (rule above), not an MCC tree.
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
