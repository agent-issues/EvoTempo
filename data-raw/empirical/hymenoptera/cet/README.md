# Hymenoptera: CET analysis pipeline (stage 1, morphology only)

Everything needed to run the four (beta, gamma) cells of the covariant
evolutionary tempo model on the Hymenoptera data of Ronquist et al. (2012),
with the fossil intervals of Zhang et al. (2016) (`../README.md`), matched to
the morphology-only RevBayes comparator of PR #79
(`../revbayes/hymenoptera_fbd.Rev`, config `morph`).  It follows the penguin
pipeline (`../../penguins/cet/`); this README lists what differs and why.
Issue #86.

**Status of this branch.**  Stacked on #98 (`gamma-acrv-97`): the
characters now use the package's model (`MkCharacters()`: ordered
characters, ambiguity sets, Mkv per missing-data mask, Gamma(4) rates among
characters with a free shape alpha ~ Exp(1)), which is the comparator's
character model.  The pipeline-local character kernel (`mk-columns.cpp`,
`HymMk()`) is gone.  Production chains start from 16 random draws from the
comparator's posterior (`draw-start-trees.R`), and the Bayes-factor scripts
follow the penguin pipeline after #103, #104, #105 and the #107 review
(labels, expected chain set, gated `cross.csv`, freshness record, ladder
checks, segment caps).  Checked on Hamilton: the fixed-tree character
likelihood against RevBayes with and without Gamma (`check-mkv.R`), the
label, chain-set and segment-cap checks (`check-labels.R`,
`check-segcap.R`), eval's gate on a smoke chain (`check-eval.sh`), and the
smoke run of the four cells (timing, exact resume), and a pilot that found
two problems to resolve before production (**Smoke run and cost
projection**), now diagnosed (**Diagnosis of the pilot's two problems**):
the runaway of cell (1, 0) is a weakly identified ridge in (lambda, s,
theta), to be removed by reparameterising the tempo process in the package
(maintainer decision 2026-11-25; the pipeline-level fix C piloted here is
superseded); the old roots come from the strict (and tempo) clock's fit to
the characters, not from the tree prior: RevBayes with a strict clock
(`morph-strict`) gives the same old roots (576 Ma against CET's 565 Ma).  The joint-prior check
(`check-prior.R`) is from the previous version: it uses no characters, so the
change of character model does not touch it.  Production is not launched.

## Files

| File | Role |
|---|---|
| `prepare.R` | Converts `../hymenoptera_morph.nex`, `../hymenoptera_characters.csv`, `../hymenoptera_tips.csv` and the RevBayes `morph` posterior into `data/` (characters, partitions, reference start tree) |
| `rb-trees.R` | Reads the comparator's posterior (`ReadRbPosterior()`) and converts its trees to CET start trees (`ConvertRbTree()`); used by `prepare.R` and `draw-start-trees.R` |
| `draw-start-trees.R` | 16 start trees drawn at random from the comparator's posterior (seed 20261111), one per production chain |
| `data/characters.csv`, `data/partitions.csv` | 113 taxa x 343 characters (cells as states, ambiguity sets `{01}`, or empty for unknown); per character K, ordered, partition `u<K>` / `o<K>` |
| `data/start-tree.nwk`, `data/start-tree-log.txt`, `data/start-tree-sample.csv` | Reference start tree (highest posterior sample; crown, no root edge); the selection and every change; its comparator sample's parameters |
| `data/start-trees/` | `start-tree-<kk>.nwk` for production chains k = 1..16, `draws.csv` (comparator sample, its parameters, conversion record, md5 per tree), `draws-log.txt` |
| `common.R` | Model, priors, moves, data loading, chain labels and the expected chain set, and the namespace rebindings (power and Holometabola term around the package's character likelihood; joint FBD-rate prior; FBD block moves; per-sample recording and exact resume) |
| `check-mkv.R`, `revbayes/check-mkv.Rev`, `revbayes/run-check-mkv.sh` | Fixed-tree check of the character likelihood against RevBayes, with and without Gamma(4), `coding = "variable"` / `"all"`, on two trees |
| `check-prior.R` | Check of the joint FBD-rate prior and the moves against an independent sampler |
| `run-chain.R` | One chain of one cell, in checkpointed segments (as the penguin pipeline); writes the FBD block moves' acceptance to `acceptance.csv` |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R` | Smoke run (five tasks; `test` partition) and its checks (timing, acceptance, start trees, step sizes, exact resume) |
| `pilot-summary.R` | ESS per 1000 iterations of a pilot, for the cost projection |
| `bayes-factors.R` | Cross-model evaluation of recorded states (`eval`) and bridge-sampling Bayes factors (`combine`), as the penguin pipeline |
| `check-labels.R`, `check-segcap.R`, `check-eval.sh` | Checks of labels, the expected chain set and `combine`'s refusals (synthetic roots), of the segment caps, and of `eval`'s gate on a real chain (ports of the penguin scripts) |
| `diag-ridge.R`, `diag-atol.R`, `diag-prior-scale.R`, `diag-rootage.R`, `diag-clock.R`, `check-prior-scale.R`, `summarise-rbstrict.R` | Diagnosis of the pilot's two problems (**Diagnosis of the pilot's two problems**); outputs in `smoke-results/diagnostics/` |
| `smoke-results/` | `prior-check*.csv`; `mkv-check-*.csv`, `revbayes-check-mkv-fd807d10-*.txt`, `revbayes-data-check.txt`; smoke and pilot outputs; check outputs |

## Data conversion (`prepare.R`)

- **Characters.** The 343 characters that the source analysis (MrBayes
  `EXSET`, `../hymenoptera_characters.csv`) and the comparator analyse,
  43 of them ordered (`TYPESET ORD`).  K is the number of observed states;
  every character uses the contiguous states 0..K-1 (checked, as in
  `../revbayes/prepare-data.R`).  Partitions: u2 246, u3 41, u4 7, u5 4,
  u6 1, u7 1, o3 33, o4 9, o5 1 characters.  `?` and `-` are unknown
  (20 729 cells).  The 22 polymorphic or ambiguous cells (14 characters) are
  **kept as state sets**, as RevBayes reads them (the penguin pipeline set
  them to unknown, because `MkLikelihood()` then took one state per cell;
  `MkCharacters()` takes sets).  No character is dropped: each varies
  among its unambiguous cells.  `LoadHymenoptera()` builds
  `MkCharacters(cells, nState = K, ordered = ordered, coding = "variable",
  nRateCat = 4)` with the tip labels as row names, in the start tree's tip
  order (the package refuses rows that are not).
- **Tip ages.** Fossil intervals from `../hymenoptera_tips.csv` (45 fossils,
  80-242 Ma) as `TempoTree(ranges = )`; `fossilAge` redraws within them.
- **Constraint and calibration.** See **Model**.

## Model, priors and the comparator (`common.R`)

`DCETModel(beta, gamma, rho = 0.0005, condition = "extant", start = "root")`
at the default grid (K = 26 bins, z in [-5, 5], central scheme), alpha = 1,
delta = 0, removal 0, equilibrium start-bin prior.  Time in Myr.

| Component | CET (this pipeline) | RevBayes `morph` comparator (#79) |
|---|---|---|
| Start, conditioning | crown (`start = "root"`); each root daughter has a sampled extant descendant (`condition = "extant"`), marginal over the root bin | `dnFBDP(rootAge = root_time, condition = "survival")`: the same event |
| rho | 0.0005 | 0.0005 |
| lambda, mu, psi | d = lambda - mu ~ Exp(100), turnover ~ U(0, 1), sampling proportion ~ U(0, 1), as a joint prior on (lambda, mu, psi) (`FbdLogPrior()`) | the same, on (d, turnover, sampling proportion) |
| Root age | offset exponential, min 315, mean 396 Ma (`agePriors$root`, `rootAge` move) | `root_time ~ dnExponential(1/81, offset = 315)` |
| Holometabola (111 tips, all fossils) | monophyletic, age ~ offset exponential min 302, mean 396 Ma, unnormalised (`TreeTerm()`) | `dnConstrainedTopology(clade_holo)` and the clamped-exponential calibration |
| Characters | 343, Mk per (ordering, K), ordered = adjacent changes at equal rates, mean rate 1; ambiguity sets summed; Mkv per missing-data mask | the same (`fnJC(K)`, `fnOrderedRateMatrix(K)`, `coding = "variable"`) |
| Among-character rates | Gamma(4), alpha ~ Exp(1) (`MkCharacters(nRateCat = 4)`, `charPriors$alpha`, `.AlphaMove()`) | Gamma(4), alpha ~ Exp(1) (posterior mean 1.10) |
| Clock | tempo-weighted: rate ~ LogNormal(-7.1, 0.5), truncated at +-4 sdlog | UCLN: clock_rate ~ LogNormal(-7.1, 0.5) |
| Rate variation | s ~ Exp(mean 0.587405) on (0, 3] (0.61% of the mass above 3 removed); theta ~ LogNormal(log(1/120), 1), truncated (reversion time median 120 Myr) | ucld_sd ~ Exp(mean 0.587405) |
| Tree moves | `nodeAge` 28, `fossilAge` 12, `slide` 17, `saJump` 5, `rootAge` 1 per iteration | FNPR, narrow, collapse/expand, node-time slide, fossil tip time, root slide |

Notes.

- **Start type and conditioning.**  The data support a crown start with
  extant conditioning on both root daughters.  Only Orthoptera and
  Paraneoptera (both extant) lie outside the constrained Holometabola, and
  all 45 fossils lie inside it, so in every tree with non-zero prior both
  root daughters have extant descendants, and no fossil can be a root child
  or a sampled ancestor at the root (which TempoMCMC() does not allow with a
  crown start).  The comparator conditions on the root age with a prior on
  it, which is what `start = "root"` with `agePriors$root` does; an origin
  start would need an origin prior that the comparator does not have.
  Marginal conditioning over the root bin (`joint = FALSE`) as the penguin
  pipeline; both variants passed the #69 high-power oracle with root start
  and extant conditioning.
- **Among-character rate variation** (#97, #98): Gamma(4) rates with shape
  alpha, free, prior Exp(1), as the comparator.  It matters: on the
  reference start tree with a strict clock, Gamma(4) at alpha 1.10 raises
  the morphology log-likelihood by 384.9, 192.1 and 1159.2 nats at rates
  0.0007, 0.002 and 0.006 (`check-mkv.R`).  Without it CET and the
  comparator would differ in the character model by far more than the clock.
- **theta** has no comparator counterpart.  Its prior scales the penguin
  one (median reversion time 20 Myr, a third of the penguin root age) by
  the tree depth: a third of the comparator's posterior root age (~360 Ma).
- **s** takes the comparator's UCLN log-sd prior, as the penguin pipeline
  (user decision 2026-11-02 for penguins); the upper bound 3 is the penguin
  pipeline's numerical bound.
- **Starting values** (`CetStart()`): lambda, mu, psi, rate (clock_rate)
  and alpha (alpha_morph) of the comparator sample that gave the chain's
  start tree (**Starting trees**), so each production chain starts at its own
  joint draw from the comparator's posterior; theta and s, which the
  comparator lacks, are drawn from their priors from a stream of their own
  (`set.seed(9000 + k)`), and redrawn while the start log posterior is not
  finite (none was: one draw each, `progress.txt`).

### Character likelihood against RevBayes (`check-mkv.R`)

`revbayes/check-mkv.Rev` (RevBayes development build fd807d10,
1.5.0-preview) scores each of the nine partitions with a strict clock at rates
0.0007, 0.002 and 0.006, from the comparator's data regenerated by
`../revbayes/prepare-data.R` (identical to the chains' `data/`, 11 files:
`smoke-results/revbayes-data-check.txt`), in five settings: `coding =
"variable"` and `"all"` without rate variation, `"variable"` with Gamma(4) at
alpha 1.10 and at 0.35, and `"all"` with Gamma(4) at 1.10.  Each value is
printed as floor(lnL) and remainder, so RevBayes is read to about 1e-11.
`check-mkv.R` builds each partition with `MkCharacters()` and scores it with
`MkLikelihood()` in both engines.  Two trees: the reference start tree and
production chain 6's start tree (`draw06`, one of the two draws whose
conversion moved an age).  Hamilton, Slurm 19473116 (RevBayes) and 19473122
(CET), EvoTempo `gamma-acrv-97` @ da4dbd0.

Totals over the nine partitions, reference tree
(`smoke-results/mkv-check-start-tree-totals.csv`):

| Rate | Setting | CET | RevBayes | Difference |
|---|---|---|---|---|
| 0.0007 | variable, no Gamma | -7165.14272496 | -7165.14272496 | 2.1e-10 |
| 0.0007 | variable, Gamma 1.10 | -6780.28057500 | -6780.28057178 | -3.2e-6 |
| 0.0007 | variable, Gamma 0.35 | -6666.04384834 | -6666.04384834 | 1.2e-9 |
| 0.0007 | all, no Gamma | -7205.49420280 | -7205.49420280 | 2.2e-10 |
| 0.0007 | all, Gamma 1.10 | -6892.31276921 | -6892.31276651 | -2.7e-6 |
| 0.002 | variable, no Gamma | -6762.85187518 | -6762.85187518 | -1.3e-11 |
| 0.002 | variable, Gamma 1.10 | -6570.76939443 | -6570.76939350 | -9.3e-7 |
| 0.002 | variable, Gamma 0.35 | -6600.63075082 | -6600.63075081 | -1.3e-9 |
| 0.006 | variable, no Gamma | -7855.82208368 | -7855.82208368 | -5.5e-12 |
| 0.006 | variable, Gamma 1.10 | -6696.58849120 | -6696.58849052 | -6.8e-7 |
| 0.006 | variable, Gamma 0.35 | -6681.25941560 | -6681.25941560 | -2.6e-9 |

- Largest difference per partition: Gamma off 1.9e-10 (reference) and
  4.1e-9 (draw06); Gamma on 8.5e-7 and 8.4e-7 (relative 1.4e-9 and
  1.2e-9), against a tolerance of 1e-6 relative.  The Gamma-on differences
  come from RevBayes's discretisation, whose AS 91 chi-square quantiles put
  the third category at alpha 1.10 at 1.014777101 against `GammaRates()`'s
  1.014777123 (2.2e-8; `mkv-check-*-gamma-rates.csv`); at alpha 0.35 the
  rates agree to 1e-11 and so do the lnL (1e-9), as #98's validation found
  on this matrix.
- C++ and R engines agree to 1.5e-15 (relative); the nine partitions sum
  to the whole matrix as `LoadHymenoptera()` builds it; the earlier no-Gamma
  totals of the pipeline's own kernel are reproduced to all printed digits.
- Power 0.5 halves the corrected lnL (Gamma mixture included) exactly, and
  the Holometabola term is added unpowered; a swapped-label tree violating
  the constraint scores -Inf; characters whose rows are not in the tree's
  tip order are refused by the package.
- Cost per evaluation of the whole matrix (C++, rate 0.002): 2.4 ms without
  Gamma, 8.6 ms with Gamma(4); 0.1 ms for the tree term; 753 ms in the R
  engine with Gamma.  (The pipeline's own kernel took 3.4 ms without Gamma.)

### The comparator's posterior is at extreme FBD rates

The comparator's posterior has lambda ~ mu ~ 44 per Myr (95% HPD 22-68),
turnover 0.99987, psi 2.2e-6, sampling proportion 5.7e-8
(`../revbayes/summary/parameters_morph.csv`): with rho = 0.0005 the tree
prior favours vast numbers of unsampled lineages.  Two consequences:

1. Numerics.  At these rates the dCET likelihood is finite and accurate: on
   the start tree, `TreeLikelihood()` with beta = 0 equals
   `FBDTreeLikelihood()` to 2.3e-6 (lambda 42.48) and 2.9e-6 (lambda 69.5)
   nats, and one evaluation takes 0.05 s (beta = 0) to 0.3 s (beta = 1,
   s = 1) (feasibility job 19428283).  At penguin-like rates (lambda 0.2,
   mu 0.15, psi 0.05) with rho = 0.0005 `TreeLikelihood()` is 37 nats
   above FBD.  Diagnosed (2026-11-11; issue draft for the maintainer): an
   accuracy failure of `TreeLikelihood()` on long branches, where D decays
   below the absolute tolerance within one branch (atol 1e-60 closes it;
   grid, conditioning and nSlice do not matter).  `TempoMCMC()`'s target is
   not affected: its implied marginal equals FBD to 7.5e-10 in the same
   case, and changes of `.LogAugmented()` with lambda, mu, psi equal those
   of FBD to < 1e-6.  This pipeline uses `TreeLikelihood()` only in this
   feasibility check.
2. Mixing.  lambda and mu lie on a ridge lambda - mu = d << lambda, along
   which `TempoMCMC()`'s one-parameter moves cannot travel (the shift move
   moves along it only when beta != 0).  `common.R` therefore adds three
   **FBD block moves** (below), as the comparator's moves on d, turnover and
   sampling proportion.

### Joint FBD-rate prior

`TempoMCMC()` evaluates priors one parameter at a time, as
`priors[[p]](value)`, at four call sites: `.ParameterSweep()` and
`.JointSweep()` (one parameter changes), `.ShiftMove()` (lambda and mu
rescaled together when beta != 0) and `.LogPrior()` (logPost).  The prior
density of (lambda, mu, psi) is

    p = 100 exp(-100 d) (lambda - mu) mu / (lambda^2 (mu + psi)^2),   lambda > mu > 0, psi > 0.

- `FbdConditional(p)` gives the priors of lambda, mu, psi: the joint
  density at the caller's current parameters (`par` in `.ParameterSweep()`,
  `st$par` in `.JointSweep()`) with `p` replaced.  They identify their caller
  by function identity and stop at any other call site, so an engine change
  that adds one fails loudly.
- `ShiftJoint()` wraps `.ShiftMove()`: the first rescaled FBD rate carries
  the joint density (old state for the old value, rescaled state otherwise)
  and the others return 0, so the move scores
  FbdLogPrior(new) - FbdLogPrior(old).
- `LogPriorJoint()` replaces `.LogPrior()`; `CetLogPrior()` is the same for
  outputs.

### FBD block moves

`FbdBlockMoves()` runs in the `shift` slot after the shift move (once per
iteration in every cell; `.DropMoves()` never drops `shift`), `CET_FBD_BLOCK`
rounds (default 1).  Three random-walk Metropolis moves on
u = (log d, logit turnover, logit sampling proportion), with
lambda = d (1 + e^u2), mu = d e^u2, psi = mu e^u3, each changing one
coordinate; the target is the paths-conditional density of the engine's
parameter sweep (augmented density at the current paths; the characters do
not depend on these rates) times the prior in u (Exp(100) on d, flat
turnover and sampling proportion, log/logit Jacobians).  That prior is
written independently of `FbdLogPrior()`, so the two kinds of move target
the same posterior only if `FbdLogPrior()`'s Jacobian is right.  Step sizes
adapt during burn-in as the engine's (`.Adapt()`, target 0.3, every 20
calls) and are checkpointed.

### Check of the prior and the moves (`check-prior.R`)

An 11-tip subtree of the start tree (Orthoptera, Paraneoptera, six extant
Hymenoptera, three fossils at fixed ages), no characters, theta and s
fixed, lambda, mu, psi free, rho 0.0005, 16 bins.  Oracle: adaptive
random-walk Metropolis on u with `TreeLikelihood()` (paths integrated out)
and the prior written in u, using neither `FbdLogPrior()` nor its Jacobian.
CET: `TempoMCMC()` with the default moves, 4 chains x 20 000 iterations.
Negative control: the same with the Jacobian left out of `FbdLogPrior()`.
Pass: |z| < 3.5 for every statistic (means and 0.1/0.5/0.9 quantiles of
log d, logit turnover, logit sampling proportion and log lambda;
z = difference / combined MCSE); the control must be rejected.

Three runs, all reported (`smoke-results/`):

| Run | Moves | Statistics | beta = 0: max \|z\| CET / control | beta = 1: max \|z\| CET / control |
|---|---|---|---|---|
| 1 (`prior-check-noblock.csv`, Slurm 19429501) | engine only | raw scale | 69.2 / 57.1 | 4.25 / 332 |
| 2 (`prior-check-rawscale.csv`, Slurm 19429537) | with FBD block moves | raw scale | 3.88 / 57.6 | 40.2 / 264 |
| 3 (`prior-check.csv`, statistics of run 2's draws, Slurm 19429687) | with FBD block moves | smooth scales | **1.36** / 57.6 | **2.30** / 56.2 |

- Run 1 found the mixing failure: at beta = 0 the chains did not reach the
  posterior (mean log lambda -1.13 against 2.51, ESS 26), while at beta = 1,
  where the shift move travels along the ridge, CET agreed with the oracle.
  This is why the block moves exist.
- Run 2 agreed on every statistic except the 0.9 quantile of the raw
  sampling proportion (z 3.88 and -40.2), whose MCSE comes from a kernel
  density at the quantile; on that long-tailed scale (mean 9e-7, 0.9
  quantile 1.6e-6) the estimate gave the oracle an SE of 0.2%, far too
  small.  The means of the same variable agreed (z 0.40 and -1.02).
- Run 3 recomputes the statistics of run 2's draws on smooth scales
  (means and 0.1/0.5/0.9 quantiles of log d, logit turnover, logit sampling
  proportion and log lambda; 16 per cell; the scales were fixed before the
  run-3 numbers were seen, but after run 2's failure, so run 3 is not a
  pre-registered test).  CET ESS 1 832-10 912 (beta = 0) and 1 747-11 703
  (beta = 1) per statistic, oracle ESS >= 9 611.  The control is rejected
  in both cells.

So `TempoMCMC()` with `FbdConditional()`, `ShiftJoint()`,
`LogPriorJoint()` and the block moves samples the intended posterior on
this target, at beta = 0 and 1, and the check can detect a missing Jacobian.
It does not exercise characters, tree moves or the Holometabola term.

## Starting trees

`prepare.R <rbmorph root>` reads the six chains of the comparator's `morph`
posterior (copied from
`/nobackup/pjjg18/evotempo-baselines-rb141/runs-dev/hymenoptera/morph_run{1..6}/output/`
to `/nobackup/pjjg18/evotempo-hym-cet/rbmorph/`, md5-checked, because the
source may be deleted).  ConsTree's Quartet consensus supports at most 100
tips (113 here), and the project does not use MCC trees, so the start is
one posterior sample.  Rule: drop the first 25% of each chain (as the
comparator's summaries), and take the sample with the highest Posterior of
the pooled remainder (ties: lowest chain, then earliest iteration).

This is the **reference start tree** (ladder rungs, power chains and checks;
`prepare.R` reproduces it byte for byte through `rb-trees.R`).

Result (`data/start-tree-log.txt`): run 5, iteration 132 100, Posterior
-7590.739 (pooled median about -7672), Likelihood -6113.807; root 323.84 Ma,
Holometabola 322.82 Ma, crown Hymenoptera 272.17 Ma; lambda 69.50,
mu 69.50, psi 1.64e-6, clock rate 0.00308, alpha 1.17; no sampled ancestors,
every fossil inside its interval, so `prepare.R` changed no age.  The sample
is a high-posterior state of the comparator's model, which includes UCLN
branch rates that CET does not have; it is a start, not an estimate.

**Production chains** start from their own trees (#89 review, R89-5, as the
penguin pipeline): `draw-start-trees.R <rbmorph root>` draws 16 samples
uniformly at random without replacement (seed 20261111) from the 11 934
pooled post-burn-in samples, and chain k = 4 (rep - 1) + cell starts from
draw k (`data/start-trees/start-tree-<kk>.nwk`) and from that sample's
lambda, mu, psi, clock rate and alpha (`CetStart()`).  The draws come from
all six comparator chains (runs 1, 2, 3, 4, 5, 6), with root ages 322.7-595.3 Ma,
Holometabola 316.4-510.9 Ma, Posterior -7723.2 to -7657.4, clock rate
0.00147-0.00239 and alpha 0.78-1.51; none has a sampled ancestor, and the
conversion moved one age in two of them (draw 6, 0.0012 Myr; draw 16,
0.0099 Myr).  The start tree's path and md5 go into the checkpoint; a resumed
chain must load the same file, and `eval` rebuilds the characters in the
chain's own tip order and stops if any recorded state has another.

## Moves and chain mechanics

As the penguin pipeline (`CetMoves()` = `TempoMoves()` plus tree moves in
proportion to the tree, with `rootAge` in place of `originAge`; segments,
checkpoints and exact resume through the `.TreeStats()` hook), plus the FBD
block moves and `TempoMCMC()`'s own move on alpha (`.AlphaMove()`, a scale
move after the parameter sweep; acceptance row `alpha`).  The checkpoint
holds the step sizes of the parameter sweep, joint sweep, root-age move,
alpha move and FBD block moves, the move counts after burn-in, the state
and the random-number state; resuming checks the cell, power, Mkv variant,
number of Gamma categories and the start tree (path and md5).  The samples
gain `alpha`, `treeTerm`, `logRootPrior` (instead of `logOriginPrior`), and
the ages of Holometabola, crown and total-group Hymenoptera for comparison
with the comparator's monitors; `logPrior` includes the joint FBD-rate prior
and the alpha prior.  `RecordingStats()` checks at every sample that the
stored character value equals power x charLL + treeTerm.

## Bayes factors

`bayes-factors.R` is the penguin pipeline's after #103, #104, #105 and the
#107 review: chain labels in plain decimal notation (`FormatValue()`,
`ParseLabel()`, which stops on an unreadable label); the expected chain set
from `<run root>/launch.csv` (or 4 cells x `CET_BF_REPS` replicates), with
`CET_BF_ALLOW_MISSING` / `CET_BF_ALLOW_INCOMPLETE` as the only opt-outs and
`chains-used.csv` listing every chain's role; `eval` writes `cross.csv`
(through a temporary file) only when its checks pass (`check`, the
recomputed log posterior minus the chain's, and `checkChar`, the
recomputed charLL minus the chain's, within `CET_BF_TOL`, default 1e-4), and
records what it read (`evalSamples`, `evalSegMd5`, `evalTol`); `combine`
refuses a stale or failing `cross.csv`, deletes its earlier outputs first,
and stops on incomplete power ladders, also after leaving chains out;
`CET_MAX_SEG` / `CET_MAX_SEG_RUNG` cap the segments read.  `eval` loads the
characters in the chain's own start tree's tip order (from its checkpoint,
md5-checked) and stops if any recorded state has another order.
Hymenoptera: the tree term, the root age prior and the alpha prior are the
same at every model point and cancel; `eval` adds them to its check, and
evaluates each gamma's character likelihood at the state's alpha.

## Convergence summaries

Not part of this PR (no chains yet).  The project uses MkPrime's
rank-normalised split R-hat and ESS (Vehtari et al. 2021), not PSRF:
MkPrime exports only `ConvergenceDiagnostics()`, which takes an
`MkPosterior` object, so a summary script will call the internal
`MkPrime:::.Rhat()` (matrix of draws x chains) and `MkPrime:::.EssMatrix()`
on the CET samples, or build an `MkPosterior`.  Summary trees: ConsTree
(Quartet consensus has a 100-tip limit, so a majority-rule or another
ConsTree method at 113 tips).

## Smoke run and cost projection

Hamilton, EvoTempo `gamma-acrv-97` @ da4dbd0 (`lib-da4dbd0`), this branch's
files, `test` partition, one core per chain.  `smoke.sh` (Slurm 19473176,
five tasks): 40 burn-in + 40 iterations per cell, thin 2, in two segments of
20 (the burn-in adapts the step sizes twice); task 5 runs cell (1, 1) in one
segment.  Each cell is replicate 1, so cells 1-4 start from draws 1-4.
`check-smoke.R` (`smoke-results/check-smoke.txt`):

| Cell (beta, gamma) | s/iteration, segment 1 (incl. burn-in) | s/iteration, segment 2 | Previous kernel, no Gamma (segment 2) |
|---|---|---|---|
| (0, 0) | 1.378 | 1.281 | 0.806 |
| (0, 1) | 2.391 | 2.349 | 1.093 |
| (1, 0) | 1.378 | 1.350 | 0.819 |
| (1, 1) | 2.396 | 2.328 | 1.116 |

- An evaluation of the characters costs 8.6 ms with Gamma(4) against 3.4 ms
  for the previous kernel without it, so iterations cost 1.6x (gamma = 0) and
  2.1x (gamma = 1) as much; with gamma = 1 the path moves (28 `node`, 29
  `tip`, 56 `branch` per iteration) rescore the characters, with gamma = 0
  they do not.  The ordered partitions' uniformization (#113) is part of this
  cost.
- Exact resume: cell (1, 1) in two segments equals one segment (max
  |difference| 0 over 10 samples x 26 columns; trees and states identical),
  with the adapted step sizes restored from the checkpoint (alpha step
  1.42-1.57 after burn-in).  The four cells started from four different
  trees; every sample has alpha.
- Peak memory 187 MB per chain (`/usr/bin/time`).
- Acceptance (segment 2, `smoke-results/smoke-acceptance.csv`): alpha
  0.20-0.35; FBD block moves d 0.30-0.50, turnover 0.20-0.35, sampling
  proportion 0.30-0.40; `nodeAge` 0.32-0.52, `fossilAge` 0.81-0.90, `rootAge`
  0.45-0.85, `slide` 0-0.021, `saJump` 0, parameter sweep 0.11-0.18.
- `eval` and `combine` on copies of the four smoke chains with a
  `launch.csv` (`smoke-results/smoke-eval.txt`, `smoke-combine.txt`): eval
  0.04 s per state for both betas and both gammas; |check| <= 1.45e-7 and
  |checkChar| <= 4.6e-12; combine runs end to end (15 states per chain, so
  its numbers mean nothing).  `check-eval.sh` on a copy of the (0, 1) chain:
  with CET_BF_TOL = 1e-15 eval exits 1 and leaves no `cross.csv`; with the
  default it writes one with both check columns and no temporary file.
- `check-labels.R` (57 checks) and `check-segcap.R` pass
  (`smoke-results/check-labels.txt`, `check-segcap.txt`).

**Pilot** (for the cost projection): one chain per cell (chains 1-4, draws
1-4), 200 burn-in + 3000 iterations, thin 10, segments of 100, as chained
15-minute `test` jobs (`run-array.sh` with `CET_OUT=$W/pilot-g`,
`CET_HOURS=0.12`; Slurm 19473251-19473267), so every chain resumed from its
checkpoint in a new process at each job boundary.  `pilot-summary.R`
(`smoke-results/pilot-summary.txt`, `pilot-ess.csv`, `pilot-speed.csv`;
rank-normalised bulk and tail ESS after dropping 25% of each chain):

| Cell | Segments done | s/iteration | Slowest quantity | Bulk ESS | Tail ESS | ESS per 1000 it. | Root age, mean (Ma) |
|---|---|---|---|---|---|---|---|
| (0, 0) | 30 / 30 | 1.33 | mean node age | 10.8 | 17.1 | 4.79 | 573 |
| (0, 1) | 30 / 30 | 2.30 | s | 1.4 | 10.6 | 0.63 | 577 |
| (1, 0) | **9 / 30** | 1.71 | charLL | 1.6 | 18.0 | 2.28 | 493 |
| (1, 1) | 30 / 30 | 2.36 | Holometabola age | 1.3 | 1.3 | 0.58 | 603 |

Two findings block production as configured (diagnosed in **Diagnosis of
the pilot's two problems** below):

1. **Cell (1, 0) runs away along a lambda-tempo ridge.**  lambda (= mu)
   rose from about 90 to a median of 1650 (maximum 3430) per Myr by segment
   9, while the mean tempo fell to x = -3.35 and s rose to 2.1-2.5, so that
   lambda e^x stayed near the comparator's ~50; the cost of a segment rose
   from 128 s to 419 s, and segment 10 did not finish in 14 consecutive
   15-minute jobs (the resumed chain replays the same iterations, so it
   stalls at the same place).  With beta = 1 and gamma = 0 only
   lambda e^x enters the likelihood: the characters do not see x, and the
   comparator's FBD-rate prior (d ~ Exp(100), turnover ~ U(0, 1)) does not
   bound lambda, so only the tempo prior (width proportional to s <= 3)
   limits the ridge.  In (1, 1) the characters see e^x and the clock prior
   on `rate` pins x, and lambda stays at 20-35; the penguin pipeline's
   (1, 0) cell is held by its truncated log-normal prior on lambda.
   Options: a bounded prior on lambda (or on turnover), a tighter bound on
   s, or a location constraint on the tempo.
2. **CET trees are much older than the comparator's.**  Every cell's root
   age drifts to 490-600 Ma (comparator: median 357, 90% 321-428 Ma), with
   root-age ESS of 2-17 in 2250 iterations, so the chains have not settled.
   The previous pilot (no Gamma, one start tree) had the same in its gamma = 0
   cells (root medians 429 and 515 Ma) but not in its gamma = 1 cells (326,
   341 Ma).  This needs a diagnosis (which terms favour deep roots: tree
   prior at turnover near 1, strict or tempo clock against UCLN, the
   calibrations) before production.

Projection (**provisional**; the slowest quantity of (0, 1) and (1, 1) gives
0.58 ESS per 1000 iterations, so a pooled ESS of 200 per quantity (50 per
chain, 4 chains) needs about 86 000 iterations per chain, and with the
penguin pipeline's 4x margin 350 000, as the previous pilot; (1, 0) is
excluded until finding 1 is resolved, and finding 2 may change everything):

| Item | Core-hours |
|---|---|
| Production: 4 cells x 4 chains x 350 000 iterations at 1.38 s (gamma = 0) and 2.43 s (gamma = 1) per iteration (8 x 134 h + 8 x 236 h) | ~2 960 |
| Wall time per chain | 134 h (gamma = 0) / 236 h (gamma = 1): 2 and 4 consecutive 70-h jobs per chain (checkpointed; resume shown above) |
| `bayes-factors.R eval` | 0.04 s per state (both betas, both gammas): ~1 core-hour in all |
| Ladder rungs, if needed (4 chains x up to 3 per edge x 4 edges, at the gamma = 1 speed) | up to ~11 400 |

## Diagnosis of the pilot's two problems

Scripts `diag-ridge.R`, `diag-atol.R`, `diag-prior-scale.R`,
`diag-rootage.R`, `diag-clock.R`; outputs in `smoke-results/diagnostics/`.
Hamilton, EvoTempo da4dbd0, Slurm 19476083-19476389.

### Cell (1, 0): a ridge along s, not a bug

- **Where the tempo process puts x.**  The generator is
  Q = diag(e^(alpha x)) A, with A the Ornstein-Uhlenbeck operator in
  z = x / s (mean 0, stationary sd 1) and alpha = 1 (`TempoGenerator()`).
  The time change tilts the lineage-stationary law to
  pi(z) proportional to phi(z) e^(-s z): mean x = -s^2 and
  E_pi[e^x] = e^(-s^2 / 2) (on the grid, s = 1: -0.987 and 0.609; s = 2:
  -3.80 and 0.142; `ridge-location.csv`).  The root-bin prior
  (`rootPrior = "equilibrium"`, Perron vector of Q^T + R) keeps
  E[e^x] at 0.93-1.05 for s = 0.25-3.  So the location is pinned, at -s^2,
  and lambda, the rate at x = 0, is not the rate lineages have: with
  beta = 1 a lineage's long-run birth rate is about lambda e^(-s^2 / 2).
- **The runaway follows that law.**  Along the pilot chain
  (`ridge-pilot.csv`), lambda rose from 94 to 1651 while s rose from 1.1 to
  2.2 and the mean path x fell from -0.29 to -3.35; the tree-weighted rate
  lambda <e^x> (`clock` / `treeLength`) rose only from 85 to 145, and
  theta <e^x> stayed at about 0.003.  With beta = 1 and gamma = 0 only the
  tree sees the tempo; it favours more rate heterogeneity among lineages
  (larger s), which moves the lineages' tempo down and makes lambda grow
  like e^(s^2 / 2).  The comparator's prior on (d, turnover, sampling
  proportion) puts a tail of only 1/lambda on log lambda, and the prior on s
  is weak, so nothing holds the chain.
- **The posterior is flat along it; the cost explodes.**  On one fixed tree
  (the chain's last), along the straight line in log-parameter space from
  the segment-1 state to the last state (`ridge-profile.csv`): the marginal
  log-likelihood (`TreeLikelihood()`, paths integrated out, atol 1e-60)
  rose -830.87, -827.97, -825.63 at t = 0, 0.2, 0.4 while the log prior fell
  -11.25, -13.40, -15.88, so likelihood + prior stayed at -842.12, -841.37,
  -841.51.  One evaluation took 0.5, 6.5 and 53 s; at t = 0.6 (lambda 495,
  s 1.69) the atol = 1e-60 solve failed after 400 s (the default atol gave
  -823.995), and at t = 0.8 and 1 (s 2.04, 2.46) both failed.
- **Penguins.**  The penguin (1, 0) production chains at 600 000 iterations
  (`ridge-penguin.csv`, 4 chains x 6 windows of 100 000) show no drift:
  lambda median 0.50-0.53, s 0.057-0.071, mean x -0.05 to -0.07.  With
  s about 0.065 the stationary law's location (-s^2 = -0.004) is at 0 and
  the ridge does not exist; their lambda prior is also a truncated
  log-normal.
- **#90.**  `check-prior.R`'s oracle uses `TreeLikelihood()` at the default
  atol.  At 200 of its posterior draws per beta (`diag-atol.R`), default
  minus atol 1e-60 lies in [-1.8e-7, 7.5e-9] (sd 1.3e-8 at beta 0, 3.3e-9 at
  beta 1), and atol 1e-60 against `FBDTreeLikelihood()` has sd 5.6e-7: the
  oracle is not affected.  `diag-ridge.R`'s profile reports both atols.

**Fix options, tested** (`common.R` switches, defaults unchanged:
`CET_FBD_SCALE = "lineage"` puts the comparator's prior on the
lineage-averaged rates k lambda, k mu with k = E_pi[e^(beta x)]
(`LineageMean()`; k = 1 at beta = 0, so cells with beta = 0 are unchanged;
s then enters the joint prior and the FBD block moves use the general
target); `CET_S_MAX` bounds the prior on s, here at 2, where the stationary
law, centred at z = -s, still lies within the grid z in [-5, 5]: its
lowest-bin mass is 9.4e-5 at s = 2 against 1.4e-3 at s = 3.)  Three pilots
of cell (1, 0) from the original pilot's start (chain 3: draw 3, same seeds),
200 burn-in + 3000 iterations (`fix-arms.csv`; quantiles 5/50/95% after 25%
burn-in, bulk ESS):

| Arm | Segments | Seconds per segment, median / max | lambda | lambda <e^x> | s | ESS s / root / tree length |
|---|---|---|---|---|---|---|
| original pilot | 9 (stalled) | 135 / 419 | 77 / 228 / 2340 | 66 / 109 / 211 | 1.03 / 1.48 / 2.33 | 14.6 / 1.9 / 1.8 |
| A: lineage-averaged prior | 20 (stalls) | 133 / 656 | 93 / 419 / 1360 | 63 / 121 / 207 | 1.00 / 1.63 / 2.02 | 16.8 / 22.8 / 3.6 |
| B: s <= 2 | 30 | 128 / 140 | 78 / 298 / 1240 | 59 / 103 / 163 | 1.07 / 1.68 / 1.98 | 21.1 / 1.9 / 1.4 |
| C: both | 30 | 128 / 135 | 46 / 119 / 592 | 50 / 86 / 143 | 0.84 / 1.28 / 1.87 | 25.4 / 25.0 / 11.3 |

- A alone makes it worse: it removes the 1/lambda brake, and the chain
  stalls as the original did.
- B alone stops the stall, but s piles up against the bound (95% quantile
  1.98), so the bound does all the work, and lambda still reaches 1240.
- C keeps lambda near the comparator's range (median 119; lambda <e^x> 86
  against the comparator's lambda of about 44), keeps s away from the bound,
  and mixes best; segment cost stays at 128-135 s.
- **Decision (maintainer, 2026-11-25): fix C is superseded and must not
  become a default.**  The tempo process is to be reparameterised in the
  package so that its location is identified (for example a stationary mean
  of 0, or lambda defined at the tempo mean).  The switches
  `CET_FBD_SCALE` and `CET_S_MAX` stay as diagnostic options, off by default.
- Longer pilot of C (for the record; Slurm 19477859, `$W/fixC-long`): cell
  (1, 0), 2 chains (draws 3 and 7), 1000 burn-in + 20 000 iterations, thin
  10, last 75% (`fixC-long.txt`, `fixC-long.csv`).  Both completed; 1.3-1.6
  s per iteration (median segment 660-813 s per 500 iterations, slowest
  2481 s).  Quantiles 5/50/95%: log lambda 3.71 / 5.20 / 6.82 (lambda
  41 / 181 / 914), s 0.83 / 1.35 / 1.88 (P(s > 1.9) = 4.1%), lambda <e^x>
  39.5 / 73.9 / 145, root age 550 / 649 / 770 Ma, clock rate
  0.000988 / 0.001148 / 0.001333; R-hat 1.006-1.065 for parameters and
  ages, 1.155 for charLL (bulk ESS 9.5); window medians of lambda 125-255
  and s 1.21-1.52, no runaway.
- Prior check of C with s free (`check-prior-scale.R`; oracle in
  (log dbar, logit turnover, logit sampling, log s) with `TreeLikelihood()`,
  negative control without the 2 log k Jacobian): **not completed**.  The
  first run (40 000 iterations per chain) exceeded its 10-h limit; the second
  (Slurm 19502537) finished 11 of 12 chains and failed on CET chain 803, where
  a proposal of s beyond the bound reached `LineageMean()`, whose generator
  overflowed (`eigen()`: "infinite or missing values"); `CetPriors()` now
  rejects s outside its support before computing k.  The rerun for chain 803
  (Slurm 19566813) was cancelled when C was superseded.  The 11 chains are in
  `$W/prior-check-scale-v2/chain-*.rds`; no statistics were computed.
- `diag-prior-scale.R` (unit checks): the default prior is unchanged
  (beta 0 against beta 1: 0; against the hand sum: 1.8e-15);
  `LineageMean(s, 1)` / e^(-s^2 / 2) = 1.0000, 1.0002, 1.0033, 1.0166 at
  s = 0.25, 0.5, 1, 1.5 (grid against continuum); the lineage prior equals
  its hand computation to 3.6e-15; the block-move Jacobian identity holds
  to 5.3e-15.  With the defaults, smoke chains of cells (1, 0) and (1, 1)
  are bit-identical to those before the change (`identity2.txt`), and
  `check-labels.R` and `check-eval.sh` pass.

### Root ages: the strict clock, not the prior

- **Priors and conditioning match the comparator term by term**
  (`../revbayes/hymenoptera_fbd.Rev`): d ~ Exp(100), turnover and sampling
  proportion ~ Beta(1, 1), rho 0.0005; root ~ Exp(1/81) offset 315 Ma
  (`CetAgePriors()`); `condition = "survival"` in RevBayes is the event that
  both root lineages leave sampled extant descendants, CET's
  `condition = "extant"`; Holometabola constraint as
  `dnConstrainedTopology()`; the calibration exp(-(t - 302) / 94) / 94 as
  RevBayes's clamped `holo_cal`.
- **Prior-only runs.**  CET: cell (0, 0) at power 0 (`CET_POWER = 0`: the
  tree term only), 2 chains x 20 000 iterations from draws 1 and 5, with
  and without the calibration (`CET_HOLO_CAL`).  RevBayes: the comparator
  (`config = "morph"`, a2a8e2db) with `under_prior = true`, 2 chains x
  150 000 generations (`$W/rbprior`, Slurm 19476129-19476130).  RevBayes's
  `ignoreAllData()` also ignores the clamped calibration node: 40% of its
  prior samples have Holometabola younger than 302 Ma, impossible under the
  calibration; reweighting by the calibration density restores it.
  Quantiles 5/50/95% after 25% burn-in (`diag-rootage.R`, `rootage.csv`):

| Source | Samples (bulk ESS of root age) | Root age (Ma) | Holometabola age (Ma) |
|---|---|---|---|
| CET prior, calibration | 1500 (748) | 316 / 328 / 362 | 303 / 313 / 335 |
| RevBayes prior, reweighted by the calibration | 4502 (1155; Kish 911) | 316 / 328 / 368 | 303 / 313 / 339 |
| CET prior, no calibration | 1500 (724) | 316 / 325 / 355 | 251 / 288 / 324 |
| RevBayes prior, no calibration | 4502 (1155) | 316 / 339 / 1393 | 260 / 314 / 1351 |
| RevBayes posterior | 11 934 (1518) | 321 / 357 / 428 | 310 / 341 / 403 |
| CET pilot (0, 0), Gamma(4) | 225 (15) | 506 / 565 / 659 | 495 / 553 / 636 |
| CET pilot (0, 0), no Gamma (`CET_NRATECAT = 1`, same start) | 225 (10) | 551 / 622 / 694 | 535 / 594 / 661 |
| CET pilot (0, 1), Gamma(4) | 225 (17) | 527 / 575 / 618 | 508 / 562 / 562 |
| CET pilot (1, 1), Gamma(4) | 225 (2) | 537 / 605 / 679 | 533 / 587 / 640 |
| previous CET pilot (0, 0), no Gamma, reference start | 225 (29) | 386 / 429 / 487 | 370 / 420 / 464 |

- **With the calibration, the two priors agree** (root median 328 Ma in
  both; 95% quantile 362 against 368 Ma), and the comparator's posterior (357 Ma) stays near it.  The CET
  posteriors lie 230-290 Myr above the prior, so the old roots come from
  CET's likelihood, not from its tree prior, root prior, calibration or
  conditioning.
- **Without the calibration they differ, in mode-finding, not in density**
  (`prior-parameters.txt`): RevBayes's prior has a second mode with low
  turnover (0.21-0.84 among samples older than 600 Ma, which are 39% of its
  samples; correlation of root age with log lambda -0.974), which CET's prior
  chains, started at turnover about 1, never visit (turnover 1 throughout,
  root at most 425 Ma): CET's moves act on logit turnover, where the two
  modes are about 10 units apart, while RevBayes slides turnover on (0, 1).
  The calibration suppresses that mode (Holometabola ages of 600+ Ma get
  weight e^-3 to e^-12), which is why the calibrated priors agree; both
  programs' posteriors have turnover about 1.  A CET sampler limitation
  worth noting, not the cause of the old roots.
- **Not Gamma-specific.**  Without Gamma, from the same start, the (0, 0)
  root is older still (median 622 against 565 Ma).  (The previous pilot's
  429 Ma came from a chain started on the 324 Ma reference tree with root-age
  ESS 29: it had not finished drifting.)
- **The strict clock prefers CET's trees by about 237 nats**
  (`diag-clock.R`, `clock-summary.csv`): on 20 trees from each source, with
  branch lengths = rate x duration and the rate maximised, the character
  log-likelihood (Gamma(4), alpha 1.1) has median -6540.5 on the
  comparator's posterior trees (root 346 Ma, tree length 7204 Myr, best rate
  0.00142) and -6303.0 on CET (0, 0) trees (root 582 Ma, 9226 Myr, best rate
  0.00115); the expected number of changes, rate x tree length, is about the
  same (10.40 and 10.52).
- **RevBayes with a strict clock reproduces CET's old roots.**  The
  comparator with every branch rate equal to `clock_rate` (config
  `morph-strict` of `../revbayes/hymenoptera_fbd.Rev`; otherwise unchanged:
  FBD prior, root prior, constraint, calibration, Mkv + Gamma(4) with alpha
  ~ Exp(1), clock_rate ~ LogNormal(-7.1, 0.5)), RevBayes fd807d10, 4 chains
  x 500 000 generations, sampled every 100, 25% burn-in
  (Slurm 19479952-19479955, `evotempo-hym-rbstrict-<k>`, 54 h each;
  `summarise-rbstrict.R`, `rbstrict-summary.csv`, `rb-morph-lnl.txt`):

| | RevBayes strict | CET (0, 0), Gamma(4) (pilot, 1 chain) | RevBayes UCLN (the comparator) |
|---|---|---|---|
| Root age, median (95%) | 575.9 (476.0-707.3) Ma | 565.2 (498.9-675.7) | 357.3 (318.4-445.3) |
| Holometabola age, median | 558.2 Ma | 552.5 | 341.4 |
| Clock rate, median (95%) | 0.001142 (0.000930-0.001389) | 0.001135 (0.000943-0.001331) | 0.001919 (0.001383-0.002738) |
| Morphology lnL, median (95%) | -6292.3 (-6313.2 to -6274.2) | -6304.6 (-6328.5 to -6290.6) | -6138.9 (-6166.5 to -6114.8) |
| alpha, median | 1.178 | 1.157 | 1.090 |
| lambda, median | 33.3 | 33.6 | 42.4 |
| R-hat (root, clock, lnL) | 1.0006, 1.0011, 1.0002 | (one chain) | 1.0061, 1.0041, 1.0018 |
| Bulk ESS (root, clock, lnL) | 3109, 3370, 8529 | 15, 14, 12 | 1128, 1163, 3232 |

  RevBayes's `Likelihood` column includes the clamped calibration node
  `holo_cal` (median -7.27 strict, -4.96 UCLN); the morphology lnL above
  subtracts it (CET's `charLL` excludes the calibration).  Every R-hat of
  the strict chains is <= 1.0011 and every bulk ESS >= 3039.
  **So the old roots are what this model gives with a strict clock**:
  CET (0, 0), which has the comparator's tree prior and character model
  with a strict clock, agrees with RevBayes's strict clock (root 565 against
  576 Ma, clock rate 0.001135 against 0.001142; CET's single pilot chain has
  ESS 12-15, so its 12-nat lower lnL is within its uncertainty).  The
  comparator's young roots come from UCLN: its branch rates gain about 154
  nats of morphology log-likelihood (-6138.9 against -6292.3) and let the
  tree stay near the prior's ages.  Neither CET clock (strict, gamma = 0;
  tempo, gamma = 1) does this on these data.  The like-for-like comparator
  of CET's gamma = 0 cells is therefore `morph-strict`, and of its gamma = 1
  cells a question for the model, not a bug.
- **#90.**  Nothing in this diagnosis depends on `TreeLikelihood()` at the
  default atol: `diag-ridge.R` reports both atols, `diag-atol.R` clears
  `check-prior.R`, and the prior-only runs use `TempoMCMC()`'s target.

## Engine gaps

- Character model: done by #97 / #98 (ordered characters, ambiguity sets,
  Mkv per mask, Gamma rates among characters); the power on the character
  likelihood and tree-only terms (constraint, calibration) are still added
  here.
- Joint or transformed priors, a move on FBD rates along the lambda-mu
  ridge, output and resume: #48.

## Deviations from the comparator

- Tempo-weighted clock in place of UCLN (the point of the comparison); s
  bounded at 3.
- Gamma(4) rates computed exactly (`GammaRates()`); RevBayes's
  `fnDiscretizeGamma()` approximates the quantiles (lnL differences up to
  about 1e-6).
- CET tree moves differ (no FNPR or narrow exchange; `slide` is the only
  topology move).
- Each production chain starts from one comparator posterior sample (tree
  and the shared parameters).
