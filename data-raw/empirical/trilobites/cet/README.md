# Trilobites: CET analysis pipeline (morphology only)

Everything needed to run the four (beta, gamma) cells of the covariant evolutionary tempo model
on the Cambrian trilobite data of Paterson, Edgecombe & Lee (2019) (`../README.md`), with the
merged RevBayes FBD + UCLN baseline (`../revbayes/`, #78) as the comparator.  It follows the
penguin pipeline (`../../penguins/cet/`); the differences are listed in each section and in
**Deviations**.  The data have no DNA, so there is a single (morphology-only) stage.

**Production is not launched by this PR.** It starts after the independent review.

The characters are scored by the package's `MkCharacters()` (#97, PR #98): the comparator's
character model **with** its Gamma(4) rate variation among characters and a free shape alpha
(**Character model**).  The first version of the pipeline (PR #91) scored them with its own code
(`PartitionedMk()`, `prune.cpp`) and had no Gamma; that code is gone, and its no-Gamma values are
reproduced by the package (`check-mkv.R`).

## Files

| File | Role |
|---|---|
| `prepare.R` | Converts `../trilobites_morph.nex` and `../trilobites_characters.csv` into `data/characters.csv` and `data/partitions.csv`, and draws one starting state per production chain from the RevBayes posterior (`data/start-tree-01.nwk` ... `16.nwk`, `data/start-values.csv`, `data/start-tree-log.txt`) |
| `common.R` | Time origin and conditioning, priors, moves, data loading (characters as `MkCharacters()`), and the three namespace rebindings: a power on the character likelihood, per-sample recording and exact resume, ingroup constraint |
| `check-mkv.R`, `revbayes/check-oracle.Rev`, `revbayes/check-oracle.sh` | Fixed-tree checks against RevBayes: the character likelihood, with and without Gamma rates, against `dnPhyloCTMC(coding = "variable")`, and the conditioning and time origin against `dnFBDP(rho = 0, condition = "sampling")` |
| `check-identity.R`, `check-identity.sh` | Whether the character-free path of `TempoMCMC()` (that of the window check and SBC-lite) gives identical chains under the build that ran those checks and the build the chains use |
| `check-window.R`, `check-window.sh` | Forward-simulation check of the all-fossil set-up (rho = 0, origin start, survival, beta = 1) |
| `check-priors.R` | Monte Carlo check of the priors matched to the comparator |
| `check-sbc.R`, `check-sbc.sh` | Free-parameter simulation-based calibration at rho = 0, beta = 1 (SBC-lite) |
| `run-chain.R` | One chain of one cell, in checkpointed segments |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R`, `check-pilot.R` | Smoke run on the `test` partition and its checks (timing, acceptance, peak memory, exact resume in one process and in a second process, the predictive stop under a fake wall limit, a kill between `acceptance.csv` and the checkpoint, `bayes-factors.R eval` on every smoke chain including two replicate-2 chains); pilot summary for the cost projection |
| `check-guards.sh` | The checkpoint's build, data and prior guards: a resume under another build, changed characters, a changed alpha prior or a pre-Gamma checkpoint stops, naming the field; under an identical copy of the install it continues |
| `bayes-factors.R` | Cross-model evaluation of recorded states (`eval`) and bridge-sampling Bayes factors (`combine`), also reweighted to the comparator's joint FBD prior |
| `check-labels.R` | Label and `combine` checks on synthetic chains, no MCMC (#100, #101, #103 review): plain-decimal labels, rung columns, the expected chain set and its opt-outs, power ladders, stale outputs, freshness (`smoke-results/check-labels.txt`; `smoke-results/mutation-check.txt`: mutants M4 and M5 of the #103 review fail it) |
| `check-combine.sh`, `check-combine.R` | `eval` and `combine` on tiny real chains with a ladder rung (b0.5g0) present, `combine` stopping on a missing `cross.csv`, the old (pre-#100) script reproducing both faults, and the md5 resume guard (#100-#102) (`smoke-results/check-combine.txt`) |
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

**Items 2 and 3 after the switch to `MkCharacters()`.**  They ran at 48e6bf7 and use no characters.
`check-identity.sh` (Slurm 19473189, `smoke-results/identity.txt`) runs, under 48e6bf7 and under
the build the chains now use (bf8ce1c), a free-parameter prior-only chain with every path, tree and
origin move and the parameter, `joint` and `shift` moves, under survival and none, and a small
`check-window.R` (300 iterations, 2 chains, 32 simulated trees per conditioning): samples,
acceptance and simulations are `identical()` between the builds.  The character-free path of
`TempoMCMC()` is unchanged, so both results stand and the window check needs no rerun.  Item 1 was
rerun on bf8ce1c: the same residuals (3.4e-8 C++, 6.7e-7 R; spreads of the negative controls 2.76
and 3.16 nats).

**Engine findings** (not used by `TempoMCMC()`, whose target is `.LogAugmented()` and whose frozen
messages use eigendecomposed slice exponentials): `FBDTreeLikelihood()` returns NaN with rho = 0
and no extant tip (it adds 0 x log 0); `TreeLikelihood()`'s default `atol = 1e-14` on the rescaled D
is too loose on long branches where D decays far (grid case 6, psi 0.3 and a 100-Myr stem: 9.3
nats too high, the engines disagree), so `check-mkv.R` uses `atol = 1e-30`, at which both engines
agree with RevBayes.  An issue draft is with the lead.

## Character model (`common.R`, `MkCharacters()`)

`LoadCharacters()` builds the characters with the package's `MkCharacters()` (#97, PR #98), with
row names equal to the tip labels of the tree they are scored on (`MkLikelihood()` checks them):

    MkCharacters(cells, nState = parts$k, ordered = parts$ordered, coding = "variable", nRateCat = 4)

- K per character from `data/partitions.csv` (the comparator's K, from the matrix symbols);
  unordered characters Mk (RevBayes `fnJC(K)`); characters 45 and 47 ordered: adjacent-state
  changes at equal rates, rescaled to mean rate 1 (`fnOrderedRateMatrix(maxState = K)`); uniform
  root; character 17 excluded (`prepare.R`);
- cells as `prepare.R` writes them: one state (`"0"`, `"12"` in the 16-state character 3), a set
  (`"0/1"`, read as the ambiguous set {0, 1}), or `""` (unknown).  A set over all K states is
  unknown, as RevBayes (character 12 has one `{01}` in a binary character);
- `coding = "variable"`: one Mkv correction per missing-data mask, number of states and ordering,
  as RevBayes `dnPhyloCTMC(coding = "variable")` (#65), averaged over the rate categories as RevBayes
  does (one correction per mask, not per category).  `CET_MKV=none` gives `coding = "all"`;
- **Gamma rates among characters**: 4 equiprobable categories with the mean rate of each
  (`GammaRates()`, RevBayes `fnDiscretizeGamma(alpha, alpha, 4)`), one shape `alpha` for all
  characters, as the comparator.  `alpha` is a parameter of `TempoMCMC()` (`start$alpha`, prior
  `charPriors$alpha`; `CetCharPriors()`: Exponential(rate 2), the comparator's), updated by the
  package's scale move once per iteration (step adapted during burn-in);
- the power of the stepping-stone chains multiplies the package's log-likelihood
  (`PoweredMk()`, the only rebinding of `MkLikelihood()`; at power 1 it returns the package's value
  unchanged).

`check-mkv.R` on start tree 01 with a strict clock (`smoke-results/mkv-check.csv`, both engines),
against RevBayes with the comparator's matrices (`../revbayes/prepare-data.R`):

| Rate | Gamma shape | CET (`MkLikelihood()`) | RevBayes `coding = "variable"` | Mkv correction |
|---|---|---|---|---|
| 0.005 | none | -2805.185423 | -2805.185423 | 31.70 |
| 0.005 | 0.25 | -2493.735803 | -2493.735803 | 96.48 |
| 0.005 | 0.6 | -2538.950564 | -2538.950564 | 66.80 |
| 0.005 | 2 | -2627.377324 | -2627.377324 | 43.70 |
| 0.0236 | none | -2667.526174 | -2667.526174 | 4.36 |
| 0.0236 | 0.25 | -2433.332262 | -2433.332262 | 56.76 |
| 0.0236 | 0.6 | -2414.623148 | -2414.623148 | 25.52 |
| 0.0236 | 2 | -2448.941189 | -2448.941189 | 7.90 |
| 0.05 | none | -3123.850355 | -3123.850355 | 1.98 |
| 0.05 | 0.25 | -2460.105471 | -2460.105471 | 43.57 |
| 0.05 | 0.6 | -2438.300303 | -2438.300303 | 14.26 |
| 0.05 | 2 | -2555.296907 | -2555.296907 | 3.28 |

Gamma off (nRateCat = 1): every partition agrees with RevBayes to 5.1e-10 (10 partitions x 3 rates x
2 engines), and the totals equal those of PR #91's own code (`PartitionedMk()`) to 1e-5, their printed
precision.  Gamma on (nRateCat = 4, alpha 0.25, 0.6, 2): every partition agrees to 5.7e-8
(10 x 9 x 2) and the totals to 1.1e-7, within the bar of 1e-6 set beforehand.  The whole matrix equals the sum of its
partitions, the engines agree to 2.3e-12, and power t gives t x lnL.  RevBayes prints 14 significant
digits (`setOption("outputPrecision", "14")`).  The calibration of the alpha move is checked in
PR #98 (`data-raw/gamma-acrv-results/`: alpha oracle and SBC with Gamma), not here.

Cost of one evaluation on start tree 01 (`smoke-results/mkv-check-timing.csv`; C++ mean of 200):

| Character likelihood | ms |
|---|---|
| `MkLikelihood()`, C++, Gamma 4 categories, `coding = "variable"` (the chains) | 4.50 |
| `MkLikelihood()`, C++, Gamma 4 categories, `coding = "all"` | 3.92 |
| `MkLikelihood()`, C++, no Gamma, `coding = "variable"` | 1.27 |
| `MkLikelihood()`, R engine, Gamma 4 categories (mean of 10) | 409 |
| `MkLikelihood()`, R engine, no Gamma (mean of 10) | 103 |

PR #91's own code took 1.40 ms (C++, no Gamma).  Four rate categories cost about 3.5 times as
much as one; the ordered characters now use uniformization (#113: +50 % per evaluation on the
Hymenoptera matrix), but here they are 2 of 106 characters.

## Priors (`common.R`)

| Parameter | CET prior | RevBayes comparator (`../revbayes/trilobites_fbd.Rev`) |
|---|---|---|
| lambda, mu, psi | independent: the **marginals** of the comparator's prior, truncated to the central 1 - 2 pnorm(-4) of their mass | joint: d = lambda - mu ~ LogNormal(-3.5, 1.5), turnover ~ U(0, 1), sampling proportion ~ U(0, 1) |
| s (tempo SD) | the comparator's UCLN log-SD prior: CV ~ Exp(rate 3), sigma = sqrt(ln(1 + CV^2)), on (0, 3] | ucld_sd from ucld_cv ~ Exp(3) |
| rate (clock at tempo 0) | U(0, 100) | clock_morph (mean branch rate) ~ U(0, 100) |
| theta | LogNormal(log 0.05, 1), truncated at 4 sdlog (as penguins; reversion time median 20 Myr) | - |
| origin | U(0.6, 1000.6) Myr before 485.4 Ma (U(0, 1000) above the youngest tip) | the same |
| rho; conditioning | 0; survival of the stem (at least one sample) | 0; `condition = "sampling"` |
| among-character rate variation | Gamma(4), shape alpha ~ Exponential(rate 2) (`CetCharPriors()`) | the same: Gamma(4), `alpha_morph ~ dnExponential(2)` (mean 0.5) |

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
  `prepare.R`; **Starting trees**): lambda, mu, psi, rate = clock_morph, s = ucld_sd,
  alpha = alpha_morph; theta 0.05.

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
The Gamma shape alpha has the package's scale move, once per iteration in every cell (the
comparator: `mvScale`, weight 3).
Per-sample trees, states and exact resume work as in the penguin pipeline (rebound `.TreeStats()`);
the resume also restores alpha's adapted step (`alphaStep`).
The checkpoint records the grid (`CET_NBIN`, `CET_ZMAX`), thinning, segment length, beta, gamma,
power, Mkv variant, number of Gamma rate categories (so a checkpoint from before Gamma rates
stops: `smoke-results/guard-gamma.txt`), the alpha prior (a tag and an md5 of the code of
`CetCharPriors()` and `PoweredMk()`), md5s of `data/characters.csv` and `data/partitions.csv`, and
the EvoTempo build: `packageVersion()`, the commit SHA that the build step writes into the install
(`<lib>/EvoTempo/BUILD_SHA`; **Running on Hamilton**) and an md5 over the installed code (review of
#114, R114-3, R114-6).  A resubmission therefore resumes only under the same build, data and
prior; a stop names the fields that differ (`smoke-results/guards.txt`).  A reinstall, even of the
same commit, changes the md5 (`check-guards.sh` case F), so resume a chain with the library it
started with; `CET_ALLOW_REINSTALL=1` resumes under a reinstall of the same commit, comparing the
version and `BUILD_SHA` only, and notes it in `progress.txt` (#114 delta review, N2).  The
checkpoint also records the start-draw index and the start tree (its md5, its path relative to `cet/`,
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
#103 review fixes to `combine` (and `eval`):
- **Outputs.** `combine` deletes `bayes-factors.csv`, `model-probabilities.csv`, `closure.txt` and
  `chains-used.csv` before anything else, so a stop, or a run with fewer edges, leaves no output of
  an earlier run (R103-3).  `chains-used.csv` lists every chain of the variant with its role
  (bridge, stepping-stone, or left out) and eval's record of it (R103-12).
- **Ladders.** As the penguin pipeline (#88), `combine` stops if a cell's character-power ladder
  lacks t = 0 or t = 1 or has fewer than three rungs, both for the chains found and once the chains
  left out with `CET_BF_ALLOW_MISSING` are removed: a chain left out is left out everywhere, and
  leaving out a cell's power-1 chains (its t = 1 rung) stops `combine` (R103-4; before, it reported
  the marginal of a truncated ladder with exit 0).
- **Freshness.** `eval` records in `cross.csv` how many sample rows it read and an md5 over the
  `seg*_samples.csv` files (`evalSamples`, `evalSegMd5`; not the states files, which may be deleted
  after eval), and its tolerance (`evalTol`).  `combine` refuses a `cross.csv` whose chain has
  changed since, or that records nothing (written before this fix: rerun eval), applies the
  stricter of `evalTol` and its own `CET_BF_TOL`, and requires every chain it uses to be complete
  (the last line of `progress.txt` is "... complete"); `CET_BF_ALLOW_INCOMPLETE=<labels>` names
  exceptions, which it echoes and records in `chains-used.csv` (R103-2).
- **Launch record.** `combine` also stops on two directories for one chain (`b0g0_run1` and
  `b0g0_run01`), on a `reps` that is not a whole number, on a launched chain directory with no
  segment yet that the record does not list, and on a record with no chain of the `CET_MKV`
  variant (R103-7 to R103-10).
- **Resume guard.** A checkpoint written before #102 has no md5; a resume from another path then
  stops saying so, and a resume from the original absolute path continues (R103-11).
`check-labels.R` sections 6-11 and `check-combine.R` test each of these.

Before these fixes, eval on a replicate-2 chain scored tree-01 character rows against another tip
order (3 of 107 positions shared): the self-check stopped it, but only after it had written a
`cross.csv` that `combine` would have used (gamma edges wrong by up to 114 nats per state).
Log names carry the variant and power (#77).  `CET_MKV` takes `pattern` (default) or `none`; the
penguin pipeline's `allTips` variant is not carried over.

## Smoke run, pilot and cost projection

Hamilton, one core per chain.  Since the switch to `MkCharacters()`: EvoTempo bf8ce1c (PR #98's
gamma-acrv-97 at da4dbd0 with #103 merged; `lib-bf8ce1c`), work directory
`/nobackup/pjjg18/evotempo-trilobite-gamma`.  Before it (PR #91): main 48e6bf7.

**Smoke with Gamma rates** (`smoke.sh`, `test` partition, 40 burn-in + 2 x 100 iterations per
chain; tasks 1-10, Slurm 19569340, checked by 19569341; `smoke-results/timing.csv`,
`acceptance.csv`, `resume.txt`, `stop-resume.txt`, `eval.txt`, `peak-memory.csv`,
`check-smoke.txt`): every cell ran from its own draw (01-04), plus cells (1, 1) and (0, 0) of
replicate 2 (draws 08 and 05).  After burn-in: 0.74 s per iteration in (0, 0), 0.76 in (1, 0), 1.10
in (0, 1) and 1.10 in (1, 1), with ten tasks sharing one node; the first smoke with Gamma (Slurm
19473061, seven tasks) gave 0.709, 0.724, 1.052 and 1.133 s, against 0.471, 0.487, 0.581 and 0.590 s
before Gamma (same draws and seeds): x 1.49-1.51 with gamma = 0 and x 1.81-1.92 with gamma = 1,
whose many local path moves each rescore the characters.  Acceptance of the alpha move 29-40%; the other moves
as before (`slide` 2.2-2.7%, `saJump` 3.7-7.3%, `nodeAge` 53-57%).  **Exact resume:** cell (1, 1)
in two segments of 100 in one process (task 4) and in one segment of 200 (task 5) gave identical
samples (maximum absolute difference 0 over 10 samples x 25 columns, alpha included) and identical
trees.  **Resume in a second process** (task 8, review of #114, R114-2): stopped by `CET_HOURS=0`
after segment 1 and resumed by a new `Rscript` from the checkpoint on disk, it equals task 4 and
task 5 (0 over 20 samples x 26 columns; trees identical).

**Predictive stop** (task 9, R114-1), in four processes, with fake wall limits set from times
measured in the same task, and `check-smoke.R` checking each inequality below from the logged
times.  The first process stopped after segment 1 (`CET_HOURS=0`): segment time 62.0 s (`segSecs`,
in the checkpoint), process wall time 64.7 s.  The second started the same chain afresh in its own
directory under `CET_WALL_HOURS` = (64.670 + 1.25 x 61.995) / 3600 = 0.0395 h, reran segment 1
(62.2 s; samples and trees equal the first process's) and stopped before segment 2: 0.018 h
elapsed + 1.5 x 0.017 h passes the limit, while 1.5 x 0.017 h alone and elapsed + 1.0 x 0.017 h do
not.  This process tests the elapsed term and the in-process update of the last segment's time:
without either, or with a factor of 1.0, it would not have stopped (#114 delta review, D114-9a,
9b).  The third resumed the first chain under `CET_WALL_HOURS` = (1.25 x 61.995 + 5) / 3600 =
0.0229 h and stopped before its first segment (0.001 h elapsed + 1.5 x 0.017 h > 0.0229 h), as a
resubmitted job must, using the segment time carried over in the checkpoint; elapsed + 1.0 x it
would not have stopped it.  This process tests only the carried-over time and the 1.5 factor.  The
fourth resumed and completed the chain.

The two stops hold only within bounds, and outside them the check fails rather than passes.  The
fresh process needs its segment 1 to take 0.9-1.12 x the first process's (here 1.003 x).  The
resumed process needs its start-up before the check to take less than 0.25 x `segSecs` + 5 s,
about 20 s (here under 4 s).  (#114 delta and follow-up checks, N1: fixed fake limits of 0.02 h
and 0.043 h were below the elapsed time, or failed unmodified code on slower nodes.)

**Kill before a checkpoint** (task 10):
killed after writing segment 2's `acceptance.csv` rows and before its checkpoint, then resumed:
`acceptance.csv` holds each segment once, and the samples and trees equal task 9's (0 over 20 x 26).
**`eval` on every smoke chain**, the two replicate-2 chains included: exit 0 and `cross.csv` written
for all six; maximum |recomputed - logPost| 4.1e-12 to 4.7e-9 (the recomputation includes alpha's
prior), maximum |recomputed - charLL| <= 5e-12.  Peak memory 158-184 MB per process
(`peak-memory.csv`); the lowest is task 9's third process, which stops before running a segment.
**Guards** (`check-guards.sh`, Slurm 19476662, `smoke-results/guards.txt`;
R114-3, R114-6): a resume under an install with another `BUILD_SHA`, from changed `characters.csv`,
with the alpha prior's rate changed, under a reinstall of the same commit (other code md5), or from
PR #91's pre-Gamma checkpoint stops, naming exactly the fields that differ; with
`CET_ALLOW_REINSTALL=1` the reinstall resumes and says so in `progress.txt`, and an identical copy
of the install continues to completion (7 of 7).
**`check-combine.sh`** (#103's tiny real chains: `eval` and `combine` with a ladder rung, the
expected chain set and the resume guard; Slurm 19476412-14) passes on this code
(`smoke-results/check-combine.txt`), and so does `summarise.R` on the smoke chains and the
comparator, with the new bulk-ESS column (Slurm 19476470).

**Before Gamma (PR #91).**  Smoke Slurm 19431320, checked by 19431322: 0.47-0.62 s per iteration;
exact resume over 10 samples x 23 columns; `eval` on six chains, |recomputed - logPost| 4.6e-12 to
1.6e-9.  **Guard tests** (Slurm 19431540, `smoke-results/guard-tests.txt`): `eval` with
`CET_BF_TOL = 1e-13` on a replicate-2 chain (residual 1.6e-9) exits 1 and leaves no `cross.csv`;
`combine` runs on the smoke chains and refuses them once one `cross.csv`'s check column is
corrupted.  Peak memory in the first smoke 303-304 MB per chain.

**Pilot with Gamma rates** (`smoke.sh` with 500 burn-in + 3000 iterations, thin 10, one chain per
cell from draws 01-04, `shared`; Slurm 19473080, summarised by 19473317 (`check-pilot.R`);
`smoke-results/pilot-gamma-*`).  ESS (MkPrime) of the 2000 iterations after the first segment:

| Cell | s/iteration (segments 2-3) | Slowest quantity: ESS / 2000 it. | lambda | s | theta | rootAge | charLL | alpha |
|---|---|---|---|---|---|---|---|---|
| (0, 0) | 0.927 | charLL 10.0 | 10.4 | 222.7 | 144.6 | 15.0 | 10.0 | 86.9 |
| (0, 1) | 1.248 | theta 3.7 | 4.5 | 44.7 | 3.7 | 14.7 | 18.6 | 172.1 |
| (1, 0) | 0.975 | charLL 4.2 | 7.1 | 196.7 | 136.3 | 19.9 | 4.2 | 159.5 |
| (1, 1) | 1.587 | logPost 3.7 | 3.8 | 5.2 | 3.8 | 9.7 | 23.3 | 166.7 |

The slowest quantities reach 1.9-5.0 ESS per 1000 iterations, as before Gamma (1.8-2.5), and the
alpha move mixes well (ESS 87-172 per 2000 iterations; acceptance 20-32%).  Not results, for
orientation only: alpha's means 0.50-0.57 per cell, against the comparator's posterior 0.52
(0.37-0.70).

**Pilot before Gamma** (`smoke.sh` with 500 burn-in + 3000 iterations, thin 10, one chain per cell,
`shared`; Slurm 19429527; `smoke-results/pilot-*`, without `gamma` in the name).  This pilot predates the #91 review: each chain started
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
quantities reach 1.9-5.0 ESS per 1000 iterations in the pilot with Gamma (1.8-2.5 before it), as in
the penguin pilot, so the penguin length is kept as the starting plan: 200 000 iterations per chain
(thin 100, segments of 10 000), about 25 000 iterations for a pooled ESS of 200 with an 8x margin.
Pilot ESS from 200 draws is a weak guide, and lambda's pilot means moved a long way between pilots
of a cell (b0g0 5.26 -> 1.69, b1g1 8.53 -> 11.23; within the Gamma pilot, b1g1's 95% range 3.4-28.1
at ESS 3.8), which suggests slow drift along the lambda ~ mu ridge (review of #114, R114-8), so the
length is set by the **run-length rule** below, not by this projection.  Seconds per iteration with Gamma rates, from the pilot
with Gamma (segments 2-3): 0.93, 1.25, 0.97 and 1.59 s (mean 1.18 s).  Before that pilot ran, the
estimate from the smoke ratios (the pre-Gamma pilot's rates times each cell's smoke ratio) was
0.83-1.38 s (mean 1.12 s, about 1000 core-hours).

| Item | Core-hours |
|---|---|
| Production: 4 cells x 4 chains x 200 000 iterations at 0.93-1.59 s (pilot with Gamma; mean 1.18 s) | **~1050** (52-88 h per chain; see **Wall clock** below: the chains of cells (0, 1) and (1, 1) need a second job) |
| `bayes-factors.R eval` | < 5 |
| Ladder rungs, if a Bayes-factor edge needs them: 4 chains x ~66 h per rung, up to 3 per edge, 4 edges | 0-3200 |
| Character-power stepping-stone, best cell only (16 powers x 2 chains x 50 000 iterations) | ~520 |
| Checks run for PR #91 (oracle, window checks, smoke, pilot, prior check, SBC-lite) | 12.8 in the first round (Slurm CPU time of every job); the fix round about 1 (smoke, oracle, guards) plus the SBC-lite (48.7) and window replicate 2 (21.0), both as allocated core-hours (elapsed x 16 cores) |
| Checks for the switch to `MkCharacters()` (oracle, identity, smoke, check-combine, guard) | < 1, plus about 5 for the pilot with Gamma (4 chains, 54-108 min each); fix round of the #114 review (smoke, guards, check-combine, summarise) < 1 |

**Wall clock** (review of #114, R114-1).  `shared` kills a job at 72 h, and a 10 000-iteration
segment takes 2.6-4.4 h, so a stop checked only after each segment (`CET_HOURS`, then 70) did not
guarantee a clean stop: at its pilot rate b1g1 would have stopped at 71.3 h, and 2% slower it would
have been killed inside segment 16; b0g1 would have finished at 70.5-71.0 h inside its first job, and
10% slower it would have been killed in segment 19 (the reviewer and verifier agree to 0.1 h).
`run-chain.R` now stops cleanly **before** a segment if the time elapsed plus 1.5 x the last
segment's would pass `CET_WALL_HOURS - CET_WALL_MARGIN` (72 - 0.25 h), and keeps `CET_HOURS` as a
hard cap, now 65 h (the launch sets it explicitly).  `CET_WALL_HOURS` is not read from Slurm: it
must match the job's `--time` (`run-array.sh`; #114 delta review, N3).  At the Gamma pilot's rates (burn-in 0.91,
2.06, 0.85 and 1.47 s per iteration; segments 2.58, 3.47, 2.71 and 4.41 h) the chains of cells
(0, 0) and (1, 0) complete in about 52 and 55 h; the chains of (0, 1) stop after segment 19 (about
67 h) and those of (1, 1) after segment 15 (about 67 h), and each finishes in a second job of about
3.5 and 22 h.  A kill can no longer leave a segment half-recorded: `acceptance.csv` is rewritten
with any earlier rows of the segment replaced (a kill between it and the checkpoint used to append
the segment's rows twice on resubmission).  Smoke tasks 8-10 test the predictive stop with a fake
limit, a resume in a second process, and a kill between `acceptance.csv` and the checkpoint.

**Run-length rule** (review of #114, R114-8; fixed before production).  After the first 2 segments
(20 000 iterations) of every chain of a cell, and again at the planned end, `summarise.R` gives per
quantity the MkPrime rank-normalised split R-hat over the cell's 4 chains (`rhat`) and the bulk ESS
over chains (`bulkEss`, `MkPrime:::.EssBulk()`), after discarding the first 25% of each chain.  The
quantities: lambda, mu, psi, theta, s, rate, alpha, logPost, charLL, meanX, treeLength, rootMa,
originMa.  A cell's run is long enough when **every** quantity has R-hat <= 1.01 and bulk ESS >= 400.
- At 20 000 iterations: project the length that reaches bulk ESS 400 from the ESS per iteration.
  If the projection exceeds 400 000 iterations per chain, report to the maintainer before
  continuing (the moves would need tuning rather than more iterations).
- At 200 000 iterations: if the rule is not met, extend all 4 chains of the cell by 100 000
  iterations (resubmit with `CET_NITER=300000`; the checkpoint continues) and check again;
  repeat.  Estimates and Bayes factors are reported only for cells that meet the rule, and a cell
  that does not is reported as such, with its R-hat and ESS.

Disk: about 20 KB per stored state (`seg*_states.rds`), so about 40 MB per chain and 0.65 GB for 16
chains; remove them once `eval` has run.

## Running on Hamilton

    W=/nobackup/pjjg18/evotempo-trilobite-gamma   # src/ = checkout, lib-<sha>/ = EvoTempo build
    cd $W/src
    # build: record the commit in the install, which every checkpoint records (R114-3)
    sha=$(git rev-parse --short HEAD); mkdir -p $W/lib-$sha
    R CMD INSTALL --no-docs -l $W/lib-$sha . && git rev-parse HEAD > $W/lib-$sha/EvoTempo/BUILD_SHA
    echo $sha > $W/out/build-sha.txt
    Rscript data-raw/empirical/trilobites/cet/prepare.R $W/rbposterior
    sbatch -J evotempo-trilobite-cet-oracle data-raw/empirical/trilobites/cet/revbayes/check-oracle.sh
    sbatch -J evotempo-trilobite-cet-identity data-raw/empirical/trilobites/cet/check-identity.sh
    sbatch -J evotempo-trilobite-cet-smoke --array=1-10 data-raw/empirical/trilobites/cet/smoke.sh
    sbatch -J evotempo-trilobite-cet-guards data-raw/empirical/trilobites/cet/check-guards.sh
    Rscript data-raw/empirical/trilobites/cet/check-smoke.R $W     # also runs eval on every smoke chain
    sbatch -J evotempo-trilobite-cet --array=1-16 \
      --export=ALL,CET_NITER=200000,CET_THIN=100,CET_SEG=10000,CET_HOURS=65 \
      data-raw/empirical/trilobites/cet/run-array.sh
    # resubmit the same command for the chains that stopped before the limit (cells (0, 1), (1, 1));
    # do NOT rebuild lib-<sha> while chains await resubmission: a reinstall, even of the same
    # commit, changes the installed code's md5 and the resume stops; if it had to be rebuilt,
    # resubmit with --export=ALL,...,CET_ALLOW_REINSTALL=1 (compares version and BUILD_SHA only).
    # The predictive stop assumes --time=3-00:00:00 (CET_WALL_HOURS=72): with a shorter --time,
    # export CET_WALL_HOURS to match.
    # launch record (only needed beyond the 4 x 4 production chains; one row per launched group):
    printf 'cell,reps\nb0g0,4\nb0g1,4\nb1g0,4\nb1g1,4\n' > $W/runs/launch.csv   # + e.g. b0.5g0,4 per rung
    # per chain, then once:
    Rscript data-raw/empirical/trilobites/cet/bayes-factors.R eval $W/runs/b1g1_run1   # must exit 0
    Rscript data-raw/empirical/trilobites/cet/bayes-factors.R combine $W/runs        # refuses failed evals and missing chains
    Rscript data-raw/empirical/trilobites/cet/summarise.R $W/runs $W/rbposterior $W/runs/summary

The scripts' defaults (`W`, the `#SBATCH -o` log paths) name the first work directory,
`/nobackup/pjjg18/evotempo-trilobite-cet`; export `W` (and pass `-o`) to use another.  The
Bayes-factor procedure is the
penguin pipeline's (its README, "Bayes factors"); it ran end to end on the smoke chains
(`smoke-results/eval.txt`, `guard-tests.txt`); those numbers carry no information and are not
reported.  `summarise.R` ran on the pilot
and the comparator (`smoke-results/pilot-summary-*`): for the comparator it gives root 528.76 Ma
(mean; 95% interval 524.82-534.09), against the median 528.53 (HPD 524.53-533.57) in
`../revbayes/README.md`, with R-hat <= 1.001 and ESS >= 4370 from MkPrime.  (Its `headline` printed
R-hat as -Inf for the single-chain pilot cells; fixed since, display only.).

## Deviations from the comparator and the penguin pipeline

1. **Among-character rate variation now matches** the comparator (Gamma(4), one shape with the
   comparator's prior; `MkCharacters()`, #97), so CET-versus-RevBayes differences no longer confound
   the clock model with rate heterogeneity (the comparator's posterior shape is 0.52, 0.37-0.70).
   PR #91 had no Gamma; this was its first deviation.
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
  priors, as the penguin pipeline did, at about 4 x 7 days of `long`
  (all-fossil trees cannot be resumed from RevBayes checkpoints).
- **Run length**: set by the run-length rule (**Smoke run, pilot and cost projection**); the
  precision target of the comparison to be resolved may raise it.
- Production waits for PR #98 (Gamma rate variation, #97) to merge, and for the independent review
  of this PR; it builds EvoTempo from main once both are in.  Every check here ran on `lib-bf8ce1c`
  (#98's head da4dbd0 plus #103); before launch, `revbayes/check-oracle.sh` and `check-identity.sh`
  are rerun on the production build (with #113's fixes, if merged, which touch the ordered
  characters' P(t) and cell parsing; the matrix has no cell that #113 N2 would misread).
