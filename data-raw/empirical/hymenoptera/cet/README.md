# Hymenoptera: CET analysis pipeline (stage 1, morphology only)

Everything needed to run the four (beta, gamma) cells of the covariant
evolutionary tempo model on the Hymenoptera data of Ronquist et al. (2012),
with the fossil intervals of Zhang et al. (2016) (`../README.md`), matched to
the morphology-only RevBayes comparator of PR #79
(`../revbayes/hymenoptera_fbd.Rev`, config `morph`).  It follows the penguin
pipeline (`../../penguins/cet/`); this README lists what differs and why.
Issue #86.

**Status of this branch.**  Data conversion, starting tree, character
kernel, joint FBD-rate prior and FBD block moves are in place and checked on
Hamilton: the joint-prior check (`check-prior.R`), the fixed-tree check of
the character likelihood against RevBayes (`check-mkv.R`, which also
measures the Gamma(4) difference), and the smoke run of the four cells
(timing, exact resume).  The cost projection rests on a pilot (see
**Smoke run and cost projection**).  Production is not launched.

## Files

| File | Role |
|---|---|
| `prepare.R` | Converts `../hymenoptera_morph.nex`, `../hymenoptera_characters.csv`, `../hymenoptera_tips.csv` and a RevBayes `morph` posterior sample into `data/` |
| `data/characters.csv`, `data/partitions.csv` | 113 taxa x 343 characters (cells as states, ambiguity sets `{01}`, or empty for unknown); per character K, ordered, partition `u<K>` / `o<K>` |
| `data/start-tree.nwk`, `data/start-tree-log.txt` | Starting time tree (crown, no root edge); the selection and every change |
| `common.R` | Model, priors, moves, data loading, and the namespace rebindings (character likelihood with ordered characters, ambiguity sets and per-mask Mkv; Holometabola constraint and calibration; joint FBD-rate prior; FBD block moves; per-sample recording and exact resume) |
| `mk-columns.cpp` | C++ kernel (compiled by `Rcpp::sourceCpp()`): per-column Mk log-likelihood for one partition with a symmetric rate matrix |
| `check-mkv.R`, `revbayes/check-mkv.Rev`, `revbayes/run-check-mkv.sh` | Fixed-tree check of the character likelihood against RevBayes `coding = "variable"` / `"all"`, and the Gamma(4) difference |
| `check-prior.R` | Check of the joint FBD-rate prior and the moves against an independent sampler |
| `run-chain.R` | One chain of one cell, in checkpointed segments (as the penguin pipeline); writes the FBD block moves' acceptance to `acceptance.csv` |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R` | Smoke run (five tasks; `test` partition) and its checks (timing, acceptance, exact resume) |
| `bayes-factors.R` | Cross-model evaluation of recorded states and bridge-sampling Bayes factors (as the penguin pipeline) |
| `smoke-results/` | `prior-check*.csv`; `mkv-check*.csv`, `revbayes-check-mkv-fd807d10.txt`, `revbayes-data-check.txt`; `smoke-timing.csv`, `smoke-acceptance.csv`; pilot outputs |

## Data conversion (`prepare.R`)

- **Characters.** The 343 characters that the source analysis (MrBayes
  `EXSET`, `../hymenoptera_characters.csv`) and the comparator analyse,
  43 of them ordered (`TYPESET ORD`).  K is the number of observed states;
  every character uses the contiguous states 0..K-1 (checked, as in
  `../revbayes/prepare-data.R`).  Partitions: u2 246, u3 41, u4 7, u5 4,
  u6 1, u7 1, o3 33, o4 9, o5 1 characters.  `?` and `-` are unknown
  (20 729 cells).  The 22 polymorphic or ambiguous cells (14 characters) are
  **kept as state sets**, as RevBayes reads them (the penguin pipeline set
  them to unknown, because `MkLikelihood()` takes one state per cell; the
  kernel here takes any tip partials).  No character is dropped: each varies
  among its unambiguous cells.
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
| Among-character rates | **none** | **Gamma(4), alpha ~ Exp(1)** (posterior mean 1.10) |
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
- **Among-character rate variation** is the one part of the comparator's
  model that CET lacks: `TempoMCMC()` has no free parameter for a site-rate
  shape, and the character kernel would need one.  It matters: on the start
  tree with a strict clock, RevBayes' Gamma(4) model at alpha = 1.10 (the
  comparator's posterior mean) raises the morphology log-likelihood by
  384.9, 192.1 and 1159.2 nats at rates 0.0007, 0.002 and 0.006
  (`check-mkv.R`, `smoke-results/mkv-check-totals.csv`).  CET's
  tempo-dependent rates (gamma = 1) vary among branches and are shared by
  all characters, so they do not model this (the numbers are at fixed
  strict-clock rates; the difference at each model's own posterior is not
  measured).  CET and the comparator
  therefore differ in the character model by far more than the clock;
  a RevBayes `morph` run without Gamma would be the like-for-like
  comparator (**decision for the user**), and adding a shape parameter to
  the kernel (one free parameter, 4x the kernel cost) is the alternative.
- **theta** has no comparator counterpart.  Its prior scales the penguin
  one (median reversion time 20 Myr, a third of the penguin root age) by
  the tree depth: a third of the comparator's posterior root age (~360 Ma).
- **s** takes the comparator's UCLN log-sd prior, as the penguin pipeline
  (user decision 2026-11-02 for penguins); the upper bound 3 is the penguin
  pipeline's numerical bound.
- **Starting values** (`CetStart()`): the comparator's posterior medians
  (lambda 42.48, d 0.0055, psi 2.18e-6, rate 0.0019), theta at its prior
  median, s 0.3; replicates > 1 jitter d, lambda, psi, theta, s and rate by
  exp(N(0, 0.2^2)).

### Character likelihood against RevBayes (`check-mkv.R`)

On the start tree with a strict clock at rates 0.0007, 0.002, 0.006,
`revbayes/check-mkv.Rev` scores each of the nine partitions on the RevBayes
development builds fd807d10 (1.5.0-preview) and a2a8e2db (1.4.2-preview,
the comparator's build; byte-identical output), from the comparator's data
regenerated by `../revbayes/prepare-data.R` (identical to the chains'
`data/`, 11 files: `smoke-results/revbayes-data-check.txt`).  Each value is
printed as floor(lnL) and remainder, so RevBayes is read to ~1e-11.

| Rate | CET per mask | RevBayes `variable` | Difference | CET none = RevBayes `all` | RevBayes `variable` + Gamma(4), alpha 1.10 |
|---|---|---|---|---|---|
| 0.0007 | -7165.142725 | -7165.142725 | 1.5e-10 | -7205.494203 | -6780.280572 |
| 0.002 | -6762.851875 | -6762.851875 | -5.5e-12 | -6765.118297 | -6570.769393 |
| 0.006 | -7855.822084 | -7855.822084 | -2.2e-11 | -7855.964459 | -6696.588491 |

Per partition, every difference is <= 1.5e-10 (`smoke-results/mkv-check.csv`);
the C++ kernel equals the pure-R reference to <= 1.5e-10 and EvoTempo's
`MkLikelihood()` on unambiguous unordered columns to <= 1.4e-11; power 0.5
halves the corrected lnL exactly and the Holometabola term is added
unpowered; a swapped-label tree violating the constraint scores -Inf.
Mask convention: the first run matched every partition but u2, which was
off by up to 6.8e-3.  In a binary character the ambiguity {01} covers every
state, and RevBayes treats it as a gap in the Mkv masks; `BuildPartitions()`
now does the same (a set covering all K states is unscored).  Cost:
3.4 ms per `HymMk()` evaluation (C++ kernel), 0.1 ms for the tree term,
201 ms for the R reference (`smoke-results/mkv-check-timing.csv`).

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

## Starting tree

`prepare.R <rbmorph root>` reads the six chains of the comparator's `morph`
posterior (copied from
`/nobackup/pjjg18/evotempo-baselines-rb141/runs-dev/hymenoptera/morph_run{1..6}/output/`
to `/nobackup/pjjg18/evotempo-hym-cet/rbmorph/`, md5-checked, because the
source may be deleted).  ConsTree's Quartet consensus supports at most 100
tips (113 here), and the project does not use MCC trees, so the start is
one posterior sample.  Rule: drop the first 25% of each chain (as the
comparator's summaries), and take the sample with the highest Posterior of
the pooled remainder (ties: lowest chain, then earliest iteration).

Result (`data/start-tree-log.txt`): run 5, iteration 132 100, Posterior
-7590.739 (pooled median about -7672), Likelihood -6113.807; root 323.84 Ma,
Holometabola 322.82 Ma, crown Hymenoptera 272.17 Ma; lambda 69.50,
mu 69.50, psi 1.64e-6, clock rate 0.00308, alpha 1.17; no sampled ancestors,
every fossil inside its interval, so `prepare.R` changed no age.  The sample
is a high-posterior state of the comparator's model, which includes UCLN
branch rates that CET does not have; it is a start, not an estimate.

## Moves and chain mechanics

As the penguin pipeline (`CetMoves()` = `TempoMoves()` plus tree moves in
proportion to the tree, with `rootAge` in place of `originAge`; segments,
checkpoints and exact resume through the `.TreeStats()` hook), plus the FBD
block moves.  The samples gain `treeTerm`, `logRootPrior` (instead of
`logOriginPrior`), and the ages of Holometabola, crown and total-group
Hymenoptera for comparison with the comparator's monitors.  `RecordingStats()`
checks at every sample that the stored character value equals
power x charLL + treeTerm.

## Bayes factors

As the penguin pipeline (`bayes-factors.R`); the tree term is the same at
every model point and cancels, and is kept in `cross.csv` for the logPost
check.

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

Hamilton, EvoTempo main @ 48e6bf7 (`lib-48e6bf7`), branch files as committed,
`test` partition, one core per chain.  `smoke.sh` (Slurm 19430092, five
tasks): 10 burn-in + 40 iterations per cell, thin 2, in two segments of 20;
task 5 runs cell (1, 1) in one segment.

| Cell (beta, gamma) | s/iteration, segment 1 (incl. burn-in) | s/iteration, segment 2 |
|---|---|---|
| (0, 0) | 0.800 | 0.806 |
| (0, 1) | 1.166 | 1.093 |
| (1, 0) | 0.798 | 0.819 |
| (1, 1) | 1.148 | 1.116 |

- Exact resume: cell (1, 1) in two segments equals one segment
  (max |difference| 0 over 10 samples x 25 columns; trees identical).  The
  pilot below also resumes across processes at every job boundary.
- Peak memory 276-277 MB per chain (`/usr/bin/time`, 282 408-283 944 kB).
- Acceptance (segment 2, `smoke-results/smoke-acceptance.csv`): FBD block
  moves d 0.15-0.40, turnover 0.25-0.35, sampling proportion 0.45-0.55;
  `nodeAge` 0.39-0.51, `fossilAge` 0.78-0.85, `rootAge` 0.50-0.85, `slide`
  0.009-0.029, `saJump` 0 (20 iterations; the comparator's posterior has
  almost no sampled ancestors, 0.013 on average), parameter sweep
  0.15-0.21.  With gamma = 1, `TempoMoves()` schedules no `tree` or `joint`
  move (28 `node`, 29 `tip`, 56 `branch`, 1 `shift`), as in the penguin
  pipeline; with gamma = 0, `tree` accepted 0.65-1.00 and `joint`
  0.24-0.39.

**Pilot** (one chain per cell, 200 burn-in + 3000 iterations, thin 10,
segments of 200; `run-array.sh` with `CET_OUT=$W/pilot`, `CET_HOURS=0.12`,
as nine chained 15-minute `test` jobs, Slurm 19430101-19430112, so every
chain resumed from its checkpoint in a new process at each job boundary).  Time per
sampling iteration in its segments: 0.78-0.82 s with gamma = 0 and
1.16-1.23 s with gamma = 1.  `pilot-summary.R` (rank-normalised bulk and
tail ESS, MkPrime's `.EssBulk()` / `.EssTail()`, after dropping 25% of each
chain; `smoke-results/pilot-ess.csv`): slowest quantity per cell, as
min(bulk, tail) ESS per 1000 iterations:

| Cell | Slowest quantity | Bulk ESS | Tail ESS | ESS per 1000 it. |
|---|---|---|---|---|
| (0, 0) | crown Hymenoptera age | 2.7 | 20.5 | 1.22 |
| (0, 1) | s | 1.5 | 29.9 | 0.68 |
| (1, 0) | tree length | 1.8 | 13.3 | 0.78 |
| (1, 1) | Holometabola age | 1.9 | 1.3 | 0.57 |

These are single 2250-iteration windows with ESS below 3, so they measure
little more than a drift; read them as an order of magnitude
(**provisional**).  The slow quantities are node ages and tree length, i.e.
tree mixing (`slide` accepted 1-3% in the smoke), and s with gamma = 1.

Projection (provisional).  A pooled ESS of 200 per quantity (50 per chain,
4 chains) at 0.57 ESS per 1000 iterations needs about 88 000 iterations per
chain; with the penguin pipeline's 4x margin, 350 000.  At 0.80 s
(gamma = 0) and 1.2 s (gamma = 1) per iteration:

| Item | Core-hours |
|---|---|
| Production: 4 cells x 4 chains x 350 000 iterations (8 x 78 h + 8 x 117 h) | ~1 560 |
| Wall time per chain | 78 h (gamma = 0) / 117 h (gamma = 1): 2-3 consecutive 70-h jobs per chain (checkpointed; resume shown above) |
| `bayes-factors.R eval` | not measured here (penguin: 0.05 s per state for all four model points) |
| Ladder rungs, if needed (4 chains x up to 3 per edge x 4 edges) | up to ~5 600 |

Before production, a longer pilot (e.g. 20 000 iterations, 2 chains per
cell, ~25-45 h each on `shared`) would give usable ESS rates; and more tree
moves (`CET_TREE_MOVES`), or the FNPR / narrow moves RevBayes has, may pay
for themselves, since tree mixing is the bottleneck.

## Engine gaps

- Character model (ordered characters, ambiguity sets, Mkv, among-character
  rate variation): #47.
- Joint or transformed priors, a move on FBD rates along the lambda-mu
  ridge, output and resume: #48.

## Deviations from the comparator

- No among-character rate variation (comparator: Gamma(4)).
- Tempo-weighted clock in place of UCLN (the point of the comparison); s
  bounded at 3.
- CET tree moves differ (no FNPR or narrow exchange; `slide` is the only
  topology move).
- The starting tree is one comparator posterior sample.
