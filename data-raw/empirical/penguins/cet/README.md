# Penguins: CET analysis pipeline (stage 1, morphology only)

Everything needed to run the four (beta, gamma) cells of the covariant
evolutionary tempo model on the penguin data of Gavryushkina et al. (2017)
(`../README.md`), and the morphology-only RevBayes comparator, on Hamilton.

Staging (user decision, 2026-10-30):

1. **Stage 1 (this folder): morphology only.** Four CET cells
   (beta, gamma) in {0, 1}^2, compared with each other by Bayes factors and with
   a morphology-only RevBayes FBD + UCLN run.
2. Stage 2: DNA as a tempo-independent partition with its own UCLN clock.
3. Stage 3: DNA under tempo, with its own gamma.

Stages 2-3 need engine work; see the design issue linked under **Engine gaps**.

The production runs (stage 1) were launched on 2026-11-11 from main 48e6bf7,
starting from a sample of the round-4 total-evidence RevBayes posterior (see
**Starting tree** and **Production runs**).

## Files

| File | Role |
|---|---|
| `select-start-tree.R` | Picks the starting tree: the post-burn-in RevBayes sample closest, by quartet distance, to the quartet consensus of the pooled samples |
| `prepare.R` | Converts `../penguins_morph.nex`, `../penguins_tips.csv`, `../penguins_constraints.csv` (empty) and a dated RevBayes tree (the sample chosen by `select-start-tree.R`) into `data/` and `revbayes/data/` |
| `data/characters.csv`, `data/partitions.csv` | 55 taxa x 201 characters (0-based states, empty = unknown); number of states K per character |
| `data/start-tree.nwk`, `data/start-tree-log.txt` | Starting time tree with the origin as root edge; every age `prepare.R` changed |
| `data/start-tree-selection.txt`, `data/quartet-consensus.nwk` | Record of the starting-tree choice (chain, iteration, distances) and the quartet consensus it was measured against |
| `data/start-tree-round1.nwk` | The earlier starting tree (from the round-1 MCC tree); `check-mkv.R`'s RevBayes oracle was computed on it |
| `common.R` | Model, priors, moves, data loading, and the two namespace rebindings (partitioned Mkv likelihood with one correction per missing-data pattern; per-sample recording and exact resume) |
| `check-mkv.R` | Fixed-tree check of the partitioned Mkv likelihood against RevBayes `coding = "variable"` on `data/start-tree-round1.nwk`, and its cost (#65) |
| `run-chain.R` | One chain of one cell, in checkpointed segments |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R` | Smoke run on the `test` partition and its checks (timing, acceptance, exact resume) |
| `bayes-factors.R` | Cross-model evaluation of recorded states (`eval`) and bridge-sampling Bayes factors (`combine`) |
| `check-labels.R` | Checks of chain labels (plain decimal notation, unparseable labels stop), of the Mkv-variant filter of `combine`'s stepping-stone (#75), and of `combine` on incomplete ladders and edges (#88) |
| `summarise.R` | Convergence (ESS, PSRF, ASDSF), clade and sampled-ancestor frequencies, crown/root/origin ages, CET and RevBayes side by side |
| `revbayes/penguins_morph_fbd.Rev`, `revbayes/run-hamilton.sh`, `revbayes/data/` | Morphology-only RevBayes comparator (UCLN, and strict as a pipeline check) on the project's RevBayes development build a2a8e2db (#66) |
| `revbayes/prior-check/` | Data-free check of the comparator's tree moves on a2a8e2db: sampled-ancestor count against BEAST2 and v1.4.1 (#66) |
| `smoke-results/` | Smoke and pilot outputs quoted below |

## Data conversion (`prepare.R`)

- **Characters.** The 202 characters of `../penguins_morph.nex`, partitioned
  by the number of states K of its `CHARSET statesK` (from the
  outgroup-inclusive matrix, as in the source analysis).  `?` and `-` are
  unknown.  The 32 polymorphic or ambiguous cells (`{01}` 23, `{12}` 8,
  `{45}` 1; in 21 characters) are set to unknown, because `MkLikelihood()`
  takes one state or `NA` per cell.  Character 81 (state 1 everywhere bar two
  `{01}`) is then constant and is dropped, as RevBayes's `coding = "variable"`
  also drops it in the total-evidence baseline.  Result: 201 characters,
  K = 2 (129), 3 (50), 4 (13), 5 (4), 6 (4), 7 (1).  The RevBayes comparator
  reads the same recoded matrix (`revbayes/data/morph_K.nex`).
- **Tip ages.** `TempoTree(ranges = )` holds each fossil's [min, max] interval
  from `../penguins_tips.csv`, and the `fossilAge` move redraws a fossil's age
  uniformly within it, so the tip-age prior is uniform on the interval, as in
  BEAST2 and the RevBayes runs.  No binning is needed: rung 6 used 1-Myr bins
  only to generate synthetic ranges for simulated fossils.
- **Constraints.** None (the source analysis has none; `prepare.R` checks the
  file is empty).

## Model and priors (`common.R`)

`DCETModel(beta, gamma, rho = 1, condition = "extant", start = "origin")` at
the default grid (K = 26 bins, z in [-5, 5], central scheme), alpha = 1,
delta = 0, removal 0, equilibrium start-bin prior.  Characters: Mk per
K-partition with the Mkv correction per missing-data pattern (as RevBayes
`coding = "variable"`; see **Mkv correction** below), no among-character rate
variation.  Time in Myr.

| Parameter | CET prior | RevBayes morphology-only comparator | Total-evidence baseline (source analysis) |
|---|---|---|---|
| lambda, mu | LogNormal(log 0.1, 1.5), truncated at +-4 sdlog | same, untruncated | via d ~ LogNormal(-3.5, 1.5), turnover ~ U(0, 1) |
| psi | LogNormal(log 0.05, 1.5), truncated | same, untruncated | via sampling proportion ~ U(0, 1) |
| theta | LogNormal(log 0.05, 1), truncated (reversion time median 20 Myr) | - | - |
| s | Gamma(0.5396, scale 0.3819) on (0, 3] | ucld_sd: the same Gamma, untruncated | ucld_sd: the same Gamma |
| rate (morphological clock) | LogNormal(-5.5, 2), truncated | clock_morph: LogNormal(-5.5, 2) | clock_morph: LogNormal(-5.5, 2) |
| origin | U(0, 160) Ma (`agePriors$origin`, `originAge` move) | U(0, 160) | U(0, 160) |
| rho; conditioning | 1; a sampled extant descendant of the stem | 1; `condition = "survival"` (same event) | same |
| among-character rates | none | none | Gamma(4), shape ~ U(0, 10) |

Notes.

- `TempoMCMC()` takes independent per-parameter priors, so the source
  analysis's prior on (d, turnover, sampling proportion) cannot be used.  The
  comparator uses the CET priors on lambda, mu, psi instead, so that CET and
  RevBayes differ only in the clock.  The truncation at +-4 sdlog (as rung 7)
  removes 6.3e-5 of prior mass per parameter.
- s takes the UCLN log-sd prior, so that the CET tempo variation and the UCLN
  rate variation start from the same prior.  This prior has 16% of its mass
  below s = 0.01, where the four cells are indistinguishable, which pulls
  Bayes factors towards 0; a log-normal prior is the obvious alternative and
  can be checked by reweighting (Z'/Z = E_post[pi'(s) / pi(s)]).  **Decision
  for the user.**  The upper bound 3 (8.7e-5 of the mass) is numerical: in the
  12-minute pilot, cell (0, 0), where s enters only the path density, adapted
  its step for s to 12.8 and proposed s = 86.8; the self-rate theta e^x then
  overflows at x = 5 s, and the `eigen()` failure in `TempoEquilibrium()`
  stops `TempoMCMC()` with an error instead of rejecting the proposal.
- Starting values (`CetStart()`): near the round-1 total-evidence posterior
  means (lambda 0.2, mu 0.15, psi 0.05, rate 0.017; theta 0.05, s 0.2),
  multiplied by exp(N(0, 0.2^2)) per replicate > 1.

### Mkv correction (#65)

Every one of the 201 characters has unscored tips (median 34 of 55), and the
matrix holds only characters that vary among their **scored** tips (character
81 was dropped for being constant there).  The correct conditional is
therefore the one RevBayes `dnPhyloCTMC(coding = "variable")` uses, with `?`
treated as a gap: one correction per distinct missing-data pattern,

    log L = sum_K [ log Mk(D_K) - sum_{patterns m in K} n_m log(1 - K P_m) ],

where P_m is the probability that every tip observed in pattern m is in state
0 (the same for every state under symmetric Mk with a uniform root), from
pruning one constant column with the other tips unknown.  The 201 characters
fall into 132 patterns within the six K-partitions.  `ConstantPatternLogLik()`
prunes all 132 columns at once: in a constant column every partial
likelihood has the form (L0, L1, ..., L1), so two numbers per node and pattern
suffice.  The correction sits inside the term raised to the character power.

Until #65 the pipeline used one constant column per partition with **every**
tip observed.  That under-corrects: on the start tree with a strict clock it
moved the fixed-tree clock MLE from 0.0181 (per pattern) to 0.0211, against
0.0222 with no correction (verifier, #65).  `.cet$mkv = "allTips"` keeps
that variant for comparison.

`check-mkv.R` checks the likelihood on the round-1 start tree
(`data/start-tree-round1.nwk`, the start tree when the oracle was computed) at three
strict-clock rates (`smoke-results/mkv-check.csv`, both engines):

| Rate | Per pattern (this pipeline) | RevBayes `coding = "variable"` | Difference | All tips observed (before #65) | Correction, per pattern | Correction, all tips |
|---|---|---|---|---|---|---|
| 0.005 | -1778.517905 | -1778.51784 | -6.5e-05 | -1861.343283 | 185.59 | 102.77 |
| 0.017 | -1656.300536 | -1656.30055 | 1.4e-05 | -1693.677249 | 47.18 | 9.80 |
| 0.05 | -1820.971062 | -1820.97102 | -4.2e-05 | -1827.535176 | 6.61 | 0.048 |

The differences from RevBayes are within its print precision (7 significant
digits per partition, summed over six partitions); the all-tips column
reproduces the pipeline's earlier values to 5e-7.  `ConstantPatternLogLik()`
equals `MkLikelihood()` on each of the 132 constant columns to 2e-14, the two
engines agree to 7e-13, and power t scales the corrected log-likelihood exactly.

Cost of one `PartitionedMk()` evaluation on the start tree
(`smoke-results/mkv-check-timing.csv`, Hamilton cn120, AMD EPYC 7702, one core; mean of 500 evaluations):

| Engine | No correction | All tips (before #65) | Per pattern |
|---|---|---|---|
| C++ | 0.41 ms | 0.66 ms | 2.34 ms |
| R | 22.47 ms | 40.79 ms | 24.26 ms |

The per-pattern correction is computed in R in both engines (one vectorised
pass over the 132 patterns), so with the C++ engine it costs 1.7 ms more than
the all-tips correction per evaluation, 3.6x the time of the all-tips
evaluation.  In the smoke run (`smoke.sh` at this commit, `test` partition,
the same node for both variants, `CET_MKV=allTips` against the default) this
raises the time per MCMC iteration from 0.18-0.20 s to 0.28-0.33 s, 1.41-1.85x
depending on the cell and segment, 1.603x on average
(`smoke-results/mkv-smoke-timing.csv`; the per-segment times, acceptance rates
and the exact-resume check of the two variants are in
`smoke-results/check-smoke-pattern.txt` and `check-smoke-allTips.txt`).  The
cost projection below therefore takes the pre-#65 pilot rate of 0.19 s per
iteration times the mean ratio 1.603, i.e. 0.305 s per iteration.  A C++
implementation in the package (#47) would remove most of the extra cost.

## Starting tree

The total-evidence baseline found that random FBD starting trees trap the
chains (crown 3-7 Ma, crown fossils outside the crown), so every chain starts
from a tree of the RevBayes total-evidence posterior.  No MCC tree is used:

    Rscript data-raw/empirical/penguins/cet/select-start-tree.R <run root> <out dir> [burn-in] [clock]
    Rscript data-raw/empirical/penguins/cet/prepare.R <out dir>/start-tree-sample.nex $(cat <out dir>/start-tree-origin.txt)

`select-start-tree.R` pools the post-burn-in samples of every chain, computes
their quartet consensus with `ConsTree::Quartet()` (default settings; at most
100 tips), and takes the sample whose topology is closest to it by the
symmetric quartet distance 2d + r1 + r2 (`Quartet::QuartetStatus()`; the
quantity `Quartet()` minimises), together with that sample's origin age.

The committed tree (`data/start-tree-selection.txt`) comes from the converged
**round-4 shared-clock** baseline (RevBayes development a2a8e2db,
`/nobackup/pjjg18/evotempo-penguins-rbdev-r4/shared_run13..16`; ASDSF 0.010):
burn-in 25% (iteration >= 37 500 of 150 000), 7504 pooled samples (1876 per
chain), 55 tips, Q = 341 055 quartets per tree.  The consensus resolves 43 of
52 splits.  Distances to it range from 10 347 to 95 871 (median 31 857).  The
closest sample, the only one at the minimum, is **shared_run14, iteration
97 140**, at distance 10 347 (0.0152 of 2Q); the next closest are at 11 307
and 11 395.  It is typical of the posterior: crown 11.76 Ma (pooled median
11.95), origin 63.142 Ma (median 64.09), 16 sampled ancestors (median 15).
ConsTree 1.0.1 (ms609/ConsTree 6124244), Quartet 1.4.0.

`prepare.R` strips annotations, puts extant tips at 0, reads zero-length
fossil tips as sampled ancestors, clamps fossil (and sampled-ancestor node)
ages into their intervals, makes parents at least 0.01 Myr older than their
children, turns a sampled ancestor that cannot stay one (e.g. at the root,
which the samplers do not allow) into a fossil tip, and adds the origin as
root edge.  On the round-4 sample it changed **nothing**
(`data/start-tree-log.txt`): root 63.0588 Ma, origin 63.1419 Ma (root edge
0.0831 Myr), 16 sampled ancestors, and no age outside its interval or closer
than 0.01 Myr to a child.  It rewrote `data/start-tree.nwk`,
`data/start-tree-log.txt`, `revbayes/data/start_tree.nwk` and
`revbayes/data/start_origin.txt` (63.14191).  The character files and
`revbayes/data/morph_K.nex` / `taxa.tsv` it regenerates were byte-identical.
Ages are as RevBayes printed them (6 decimals per branch).

The earlier starting tree, from the **round-1** MCC tree (v1.4.1, shared
clock; origin = posterior median 64.39 Ma; root 60.99 Ma, 27 sampled
ancestors; *Waimanu manneringi* turned from a root sampled ancestor into a
fossil tip), is kept as `data/start-tree-round1.nwk` for `check-mkv.R`.  The
smoke and pilot results quoted below ran from it unless stated otherwise.

## Moves and chain mechanics

`CetMoves()` = `TempoMoves()` (with gamma = 1: 14 `node`, 7 `tip`,
21 `branch`, 1 `shift`; with gamma = 0: 1 `tree`, 2 `node`, 2 `tip`,
2 `branch`, 1 `joint`, 1 `shift`) plus tree moves in proportion to the tree:
`nodeAge` 7 (0.25 x bifurcations), `fossilAge` 9 (0.25 x fossils), `slide` 5,
`saJump` 4, `originAge` 1; `CET_TREE_MOVES` scales them.  These counts are
for the round-1 tree (27 sampled ancestors).  On the round-4 tree (16 sampled
ancestors, so more bifurcations) they are, with gamma = 1, 14 `node`, 10 `tip`,
23 `branch`, 1 `shift`, and, for every cell, `nodeAge` 10, `fossilAge` 9,
`slide` 6, `saJump` 4, `originAge` 1 (gamma = 0 path moves unchanged;
`progress.txt` of each chain).

`TempoMCMC()` returns neither per-sample trees nor a way to continue a chain,
so `common.R` rebinds `.TreeStats()` (called by `TempoMCMC()` at every stored
sample, and once before its loop).  At each sample the hook reads the chain
state from the caller's frame and records the Newick tree, the state (tree,
parameters, paths, bins) and the character log-likelihood; it also captures
everything the sampling phase carries between iterations (`st`, the step
sizes `step`, `jointStep`, `ageStep`, the move counts, the RNG state).  On
the first call of a later segment it injects the captured state.  A chain
therefore runs as segments of `CET_SEG` iterations with a checkpoint after
each, and continues **exactly**: in the smoke run, cell (1, 1) in two
segments of 100 and in one segment of 200 gave identical samples (maximum
absolute difference 0 over 10 samples x 21 columns) and identical trees.  A
job stops cleanly after `CET_HOURS` (70) and resubmitting the array continues
every chain.  The hook depends on `TempoMCMC()` internals (local variable
names, the pre-loop `.TreeStats()` call), so it must be rechecked
(`smoke.sh` task 5) after any engine change; the engine issue below asks for
this as API.

Chain labels and sample columns (#74).  A chain's directory is
`<cell>[_mkv<variant>][_t<power>]_run<rep>`: the Mkv variant appears only when
it is not the default (`CET_MKV=allTips` or `none`), and the power only when
it is not 1.  The checkpoint records the variant, power, beta and gamma, and
`run-chain.R` stops rather than resume a checkpoint of another target.
`bayes-factors.R eval` takes the variant from the chain's label (and checks it
against the checkpoint), and `combine` bridges only chains of the variant
given by `CET_MKV`. Its character-power stepping-stone also uses only chains of
that variant (#75; before, it pooled the charLL of every variant's power chains).
Numbers in labels are written in plain decimal notation, with the fewest
digits that read back exactly (`FormatValue()` in `common.R`; e.g. `_t0.0001`,
not `_t1e-04`), and a directory holding chain output whose label
`ParseLabel()` cannot read stops `eval` and `combine` instead of being
skipped (#75).  `combine` also stops, naming the missing rung, when a cell's
power ladder does not run from t = 0 to t = 1 (the t = 1 rung is the cell's
power-1 chains, which must be under the same root), and leaves out a bridge
edge that has no chains at one of its end points, saying which (#88).
`check-labels.R` checks all of these.
In `seg*_samples.csv`, `logPrior` is the log prior of the
parameters (as in `cross.csv`), and `logOriginPrior` is that of the origin
age; the chain's `logPost` includes both.  Before #74 the samples' `logPrior`
included the origin term (-log 160), unlike `cross.csv`.

## Running on Hamilton

    W=/nobackup/pjjg18/evotempo-penguin-cet     # src/ = checkout, lib-<sha>/ = EvoTempo build
    cd $W/src
    Rscript data-raw/empirical/penguins/cet/select-start-tree.R \
      /nobackup/pjjg18/evotempo-penguins-rbdev-r4 $W/start-tree 0.25 shared
    Rscript data-raw/empirical/penguins/cet/prepare.R $W/start-tree/start-tree-sample.nex \
      $(cat $W/start-tree/start-tree-origin.txt)
    sbatch -J evotempo-penguin-cet-smoke --array=1-5 data-raw/empirical/penguins/cet/smoke.sh
    Rscript data-raw/empirical/penguins/cet/check-smoke.R $W
    sbatch -J evotempo-penguin-cet --array=1-16 \
      --export=ALL,CET_NITER=200000,CET_THIN=100,CET_SEG=10000 data-raw/empirical/penguins/cet/run-array.sh
    # per chain, then once:
    Rscript data-raw/empirical/penguins/cet/bayes-factors.R eval $W/runs/b1g1_run1
    Rscript data-raw/empirical/penguins/cet/bayes-factors.R combine $W/runs
    Rscript data-raw/empirical/penguins/cet/summarise.R $W/runs $W/rbmorph $W/runs/summary

Jobs are named `evotempo-penguin-cet-*` / `evotempo-penguin-rbmorph-*`.
Remove `$W/runs/*/seg*_states.rds` once `bayes-factors.R eval` has run (they
are the bulk of the output).

RevBayes comparator, one chain (from a copy of `revbayes/` with `data/`):

    sbatch -J evotempo-penguin-rbmorph-ucln-1 --export=ALL,CLOCK=ucln,REP=1 run-hamilton.sh

## Bayes factors (`bayes-factors.R`)

Rung 7 (`data-raw/rung7-power.R`, PR #40) compared cells on a fixed true tree:
tree-only evidence by quadrature (beta = 0) or bridge sampling on
`MarginalMCMC()` draws (beta = 1), the gamma = 0 character evidence by 1-D
quadrature, the gamma = 1 character evidence by a character-power
stepping-stone ladder, and the beta contrast at gamma = 1 by a two-sample
bridge between the (0, 1) and (1, 1) posteriors.  With tree uncertainty the
quadratures and the tree-only evidence have no counterpart (the latter would
integrate over all labelled time trees on these tips).  The bridge does
generalise:

- All four cells share one state space x = (time tree with fossil ages and
  sampled ancestors, tempo paths and bins, lambda, mu, psi, theta, s, rate)
  and one prior on it, and for cells A, B
  log q_B(x) - log q_A(x) = [log p(T, paths | Theta, beta_B) - log p(T, paths | Theta, beta_A)]
  + [log Mkv(D | T, paths, rate, gamma_B) - log Mkv(D | T, paths, rate, gamma_A)],
  computable at every recorded state (`eval`; the chain's own value, with the
  origin age prior that `TempoMCMC()` has added to logPost since #56, is
  checked against its logPost: maximum residuals 4.6e-12 to 9.3e-8 over the 20 smoke states of each cell, `smoke-results/bf-eval.txt`, against a constant 5.075 = log 160 before #67).
- Two-sample bridge sampling (Meng & Wong 1996; `BridgeTwo()` from rung 7)
  on the draws of both posteriors estimates log Z_B / Z_A: the Bayes factor
  for all the data (characters, tip ages, sampled ancestors), with the tree
  integrated by the MCMC.
- Four edges: beta at gamma = 0 and at gamma = 1; gamma at beta = 0 and at
  beta = 1.  The square must close:
  e(00->10) + e(10->11) - e(00->01) - e(01->11) = 0 up to error, a check no
  single edge gives.  Posterior model probabilities (equal prior odds) take
  log Z(1, 1) as the mean of its two paths.
- Errors: the bridge relative-MSE SE (Fruhwirth-Schnatter 2004; rung 7 found it
  ~1.6x too small) and the spread over replicate pairs (rep k of A with rep k
  of B).
- Overlap.  The gamma edges are the risky ones: under gamma = 0 the paths
  ignore the characters, so the (beta, 0) draws may carry paths that score
  poorly under gamma = 1.  If an edge has bridge-weight ESS < 100, one-sided
  importance estimates more than 2 nats apart, or the square fails to close by
  more than 3 SE, add ladder rungs at intermediate model values along that
  edge (cell labels such as `b1g0.5`; `DCETModel()` takes numeric beta and
  gamma), and `combine` bridges consecutive rungs and sums them.
- Comparison with RevBayes.  Full-data evidences are not comparable across
  programs (the tree densities differ in conventions and in the lambda/mu/psi
  priors' truncation), but log p(characters | tip ages, model) is: the
  RevBayes comparator's `ML=true` runs `powerPosterior()` and the
  stepping-stone sampler, and CET chains at character powers t (`CET_POWER`,
  labels `<cell>_t<t>`, t = 0 samples the CET tree prior) give the same
  quantity by stepping-stone in `combine`.  Optional, and only for the best
  CET cell.

The procedure ran end to end on the smoke chains (15 states per cell); those
numbers carry no information and are not reported.

## Smoke run and cost projection

Hamilton, EvoTempo `origin-root-age-18` @ b99e810, `test` partition, one core.
**These smoke and pilot runs predate #65.** They used the all-tips Mkv
correction (now `CET_MKV=allTips`), so the times, acceptance rates, ESS rates
and crown ages in this section are those of the earlier target.  Under the
per-pattern correction each iteration costs 1.603x as much on average (see
**Mkv correction**); acceptance rates in a 240-iteration smoke of both
variants are in `smoke-results/check-smoke-pattern.txt` and
`check-smoke-allTips.txt`.  The ESS per iteration has not been measured under
the per-pattern target, and it may differ.

| Cell (beta, gamma) | s/iteration, smoke (200 it.) | s/iteration, pilot (3000 it.) | Pilot ESS of slowest quantity / 3000 it. |
|---|---|---|---|
| (0, 0) | 0.193 | 0.188 | 14 (treeLength) |
| (0, 1) | 0.172-0.178 | 0.171-0.173 | 6 (theta); s 8 |
| (1, 0) | 0.191-0.200 | 0.192-0.195 | 15 (treeLength) |
| (1, 1) | 0.177 | 0.180-0.195 | 7 (s); treeLength 7, charLL 8 |

The pilot (`smoke-results/pilot-*`) is one 12-minute chain per cell (500
burn-in + 3000 iterations, thin 10) from the round-1 MCC tree, with the final
priors; the 200-iteration smoke (`smoke-results/timing.csv`,
`acceptance.csv`) ran before the s prior was bounded.  Pilot acceptance
(last 1000 iterations): `slide` 2.3-4.3%, `saJump` 25-32%, `nodeAge` 46-52%,
`fossilAge` 61-67%, `originAge` 24-35%, local path moves 86-100%, `shift`
54-65%, parameter sweep 28-36%; with gamma = 0, whole-tree draws 87-100% and
`joint` 27-34%.  `slide` is the only topology move, and at 2-4% it is the
likely bottleneck for tree mixing (RevBayes adds FNPR and narrow exchange).
Pilot crown ages (provisional, single short chains, pre-#65 target):
12.2-13.5 Ma, against 12.5 Ma in the round-1 total-evidence baseline.

Projection (provisional: ESS from single 3000-iteration pre-#65 chains
understate slow drifts, and may not carry over to the per-pattern target).
The slowest quantities reach about 2 ESS per 1000 iterations.  A pooled ESS
of 200 per quantity (50 per chain) then needs about 25 000 iterations per
chain; with a 4x margin, 200 000 iterations per chain.  At the pre-#65 rate
of 0.19 s that was 10.6 h per chain; at 0.19 s x 1.603 = 0.305 s it is
16.9 h.  The table uses 0.305 s per iteration.

| Item | Core-hours |
|---|---|
| Production: 4 cells x 4 chains x 200 000 iterations | ~270 (170 before #65) |
| `bayes-factors.R eval` (about 0.05 s per state for all four model points: 20 states took 0.8 s in `smoke-results/bf-eval.txt`, 1.1 s in a rerun) | < 2 |
| Ladder rungs, if needed: 4 chains x 16.9 h per rung, up to 3 per edge, 4 edges | 0-810 (0-510 before #65) |
| Character-power stepping-stone, best cell only (16 powers x 2 chains x 50 000 it.) | ~135 (85 before #65) |
| RevBayes comparator (first runs biased by `mvUpDownScale` and by v1.4.1's `mvCollapseExpandFossilBranch`; rerun on development fd807d10, Slurm 19429797-19429804): UCLN 4 x 300 000 it. at ~0.43 s (first hour; 0.23 s in the smoke), strict 4 x 300 000 at ~0.05 s | ~145 + ~17 |
| RevBayes `powerPosterior` (UCLN, 32 powers x 10 000 it.) | ~20 per run |

Wall time: one 3-day job slot per chain suffices (about 17 h per chain at
0.305 s per iteration); peak memory in the pre-#65 pilot was 178-181 MB per
chain.

## Production runs (stage 1, launched 2026-11-11)

- Code: EvoTempo main **48e6bf7** (checkout `$W/src`, built into `$W/lib-48e6bf7`;
  `$W/out/build-sha.txt`), with this branch's start tree (`prepare.R` output as
  committed here); `$W = /nobackup/pjjg18/evotempo-penguin-cet`.
- Checks at 48e6bf7: `check-mkv.R` passed (round-1 tree; per-pattern lnL
  minus RevBayes `variable` -6.5e-5, 1.4e-5, -4.2e-5 at rates 0.005, 0.017,
  0.05, both engines).  `smoke.sh` (Slurm 19428287, `test`, from the round-4
  tree) and `check-smoke.R` (`smoke-results/check-smoke-48e6bf7.txt`): 0.30-0.36
  s per iteration, exact resume (maximum difference 0 over 10 samples x 23
  columns, trees identical), no non-finite value in any sample column, peak
  memory 169-171 MB.
- Array: Slurm **19428297** (`evotempo-penguin-cet`, tasks 1-16 = 4 cells x 4
  replicates), `CET_NITER=200000, CET_THIN=100, CET_SEG=10000`, burn-in 2000,
  `CET_MKV` default (`pattern`), s ~ Gamma(0.5396, scale 0.3819) on (0, 3],
  origin U(0, 160).  Submitted through `$W/jobs/run-array-wrap.sh`, which adds
  the batch header Hamilton jobs need (`set +e; set +o pipefail; source
  /etc/profile || true; module purge || true; module load r/4.5.1`) and then
  sources `run-array.sh` with the same `#SBATCH` settings.  To continue the
  chains, resubmit the same command:

      sbatch -J evotempo-penguin-cet --array=1-16 \
        --export=ALL,CET_NITER=200000,CET_THIN=100,CET_SEG=10000 $W/jobs/run-array-wrap.sh

- Segment 1 (2000 burn-in + 10 000 iterations; all 16 chains done at about
  11:20 GMT): 0.289-0.414 s per iteration (means: (0, 0) 0.310, (0, 1) 0.369,
  (1, 0) 0.303, (1, 1) 0.363), i.e. about 302 core-hours and 16-23 h per chain
  for 202 000 iterations.  No non-finite values.  Acceptance ranges over the four
  replicates: `slide` 2.8-4.2%, `saJump` 26-32%, `nodeAge` 41-51%,
  `fossilAge` 59-66%, `originAge` 22-39%, `shift` 35-61%, parameter sweep
  25-36%, `joint` (gamma = 0) 28-36%, `tree` 100% in (0, 0) and 79-82% in (1, 0).
  In b1g1_run4, logPost reached +82 at the end of segment 1 (others about
  -2000), at theta = 1.22, with an ordinary charLL (-1607) and logPrior
  (-10.5).  The #89 review traced this to the tempo-path terms of
  `.LogAugmented()` (`R/Paths.R`).  Write N_kj for the number of k -> j jumps
  and T_k for the time in bin k.  The jump and holding terms are
  A = sum_kj N_kj log Q_kj + sum_k T_k Q_kk, and under the CTMC
  (E[N_kj] = Q_kj E[T_k]) they have expectation
  E[A] = L sum_k p_k sum_{j != k} Q_kj (log Q_kj - 1), with L the total path
  length and p the occupation law.  Q is proportional to theta, so A grows
  roughly like N (log theta + 0.81) ~ theta log theta.  It changes sign near
  theta = 0.44.  Going from theta of about 0.16 (typical in b1g1) to 1.22
  predicts about +2160 nats; the chain shows +2160 to +2240 in logPost, and
  +2179 in `.LogAugmented()`.  These magnitudes use the reviewers' jump
  counts and are provisional: the state-level recomputation was held.
  Further evidence that the target is right:
  - the C++ engine, the R engine and the logged CSV agree at that state;
  - in b0g0, theta and s recover their priors (KS p = 0.49 and 0.73);
  - the chain is not stuck (theta 0.095 -> 1.632 -> 0.30).
  **logPost is therefore not comparable across theta, or between cells.**  It
  is not a fit statistic, and convergence summaries should not use it (#94).
  Theta mixes slowly in the gamma = 1 cells, which have no `tree` or `joint`
  move (#95).
- Restart decided (#89 review, R89-5): all 16 chains above start from the same
  tree, with parameter jitter exp(N(0, 0.2^2)) and none for replicate 1,
  which weakens every between-chain diagnostic.  Production is to be
  restarted with a distinct round-4 posterior sample per chain and wider
  jitter.  **Not yet done:** the step that draws the dispersed start trees
  was held by a safety-approval check, so array 19428297 is still running.

## Morphology-only RevBayes comparator

No morphology-only baseline existed (the total-evidence baseline runs shared
and separate clocks with DNA), so `revbayes/` adds one, matched to the CET runs
in data, character model, tree prior, priors on lambda, mu, psi, tip ages and
starting tree, and differing only in the clock: UCLN (non-centred), and strict
as the RevBayes counterpart of cell (0, 0), a check of the pipeline (the two
should agree on trees and FBD rates).  Smoke (600 iterations incl. 100
tuning): strict 25 s, UCLN 137 s.  Four chains of each, 300 000 iterations
(sample every 100, checkpoints every 5000), were launched on 2026-10-31 under
`/nobackup/pjjg18/evotempo-penguin-cet/rbmorph/<clock>_run<k>/` (Slurm jobs
19299263-19299270).

**Those first comparator runs are biased and must be rerun**, for two
reasons.
1. Their script had an `mvUpDownScale` move (`clock_morph` up, `fbd_tree`
   down). That move's Hastings ratio is wrong when the tree has sampled
   ancestors and on origin-conditioned trees (shown by prior-only runs), and
   the penguin tree has both. `penguins_morph_fbd.Rev` no longer has the move.
2. They ran on RevBayes v1.4.1 (`time_*.txt`), whose
   `mvCollapseExpandFossilBranch` is biased when fossil ages are uncertain
   (fixed upstream in c389174b). All 36 penguin fossils have age ranges, and
   the comparator keeps the move (`penguins_morph_fbd.Rev:77`). Without data,
   v1.4.1 gives about 32 sampled ancestors, against 16.6 in BEAST2 and 17.2
   on a2a8e2db (below). `run-hamilton.sh` therefore defaults `RB` to the
   project's development build a2a8e2db
   (`/nobackup/pjjg18/revbayes-dev-a2a8e2db/projects/cmake/build/rb`, the
   build the Hymenoptera and trilobite baselines use) and writes the binary
   and its `rb --version` to `out/rb_version_<CLOCK>_<REP>*.txt` (#66).

**Rerun (launched 2026-11-11).**  The comparator is rerunning on RevBayes
development **fd807d10** (1.5.0-preview; `/nobackup/pjjg18/revbayes-dev-fd807d10`),
the project's current build, rather than a2a8e2db.  Outputs are in
`/nobackup/pjjg18/evotempo-penguin-cet/rbmorph-fd807d10/{ucln,strict}_run1..4`,
from Slurm jobs 19429797-19429800 (`evotempo-penguin-rbmorph-ucln-1..4`) and
19429801-19429804 (`evotempo-penguin-rbmorph-strict-1..4`).  At 14:28 GMT
ucln-1 was running and the other seven were pending.  No comparator results
exist yet.  **CET cell (0, 0) is not interpreted until the strict comparator
agrees with it.**  The strict clock is the RevBayes counterpart of (0, 0), and
CET's b0g0 currently gives lambda of about 0.49 < mu of about 0.55
(provisional, early segments).  A disagreement there would be a tree-prior
question affecting all four cells.

**Data-free check of the comparator's tree moves** (`revbayes/prior-check/`).
On a2a8e2db, the total-evidence baseline script
(`../revbayes/penguins_fbd.Rev`, `under_prior = true`, which calls
`ignoreAllData()`) was run with exactly the comparator's six tree moves and
weights (`mvFNPR` 15, `mvNarrow` 15, `mvCollapseExpandFossilBranch` 10,
`mvNodeTimeSlideUniform` 40, `mvRootTimeSlideUniform` 5,
`mvFossilTipTimeUniform` 1 per ranged tip), with sampled ancestors, an origin
and ranged tip ages. It differs from the comparator only in the FBD
parameterisation, which these moves do not touch, and its result can be
compared with a data-free BEAST2 run of the source XML:

| Program | Sampled ancestors, mean (sd) | MCSE | Sampling proportion |
|---|---|---|---|
| RevBayes a2a8e2db | 17.18 (6.83) | 0.23 (ESS 873) | 0.266 |
| RevBayes v1.4.1 | 31.97 (3.00) | 0.07 (ESS 2013) | 0.735 |
| BEAST2 2.4.8 + SA 1.1.7 | 16.57 (6.85) | ESS not recorded | 0.254 |

a2a8e2db is close to BEAST2 but not shown to agree with it: the gap of 0.61
sampled ancestors is 2.6 a2a8e2db MCSE (0.23), or about 1.9 combined MCSE if
BEAST2's MCSE were the same; BEAST2's ESS was not recorded, so its MCSE is
unknown.  A residual difference of about 3.5% is therefore not excluded.  The
check does refute v1.4.1 (32 against 16.6 sampled ancestors), which is what
the choice of build needs. Details and caveats (the run
directory was deleted before this commit, so only the archived summaries are
kept) are in `revbayes/prior-check/README.md`.

What this does not check: node ages against an exact oracle.

## Engine gaps

Workarounds live in `common.R`; the issues describe the engine support that
would replace them.

- Character model: per-character number of states, Mkv correction, ambiguous
  states, among-character rate variation: #47.
- `TempoMCMC()` output and control: per-sample trees and states, character
  log-likelihood column, character power, checkpoint/resume, bounded step
  adaptation, errors other than `infeasibleProposal` in proposals, topology
  mixing, joint or transformed priors: #48.
- DNA (stages 2-3): #49.

## Deviations from the source analysis

- 32 ambiguous cells as unknown; character 81 dropped (201 characters).
- No among-character rate variation (source: Gamma(4)).
- Independent priors on lambda, mu, psi instead of the (d, turnover, sampling
  proportion) prior; s and the clocks keep the source priors (s bounded at 3).
- No DNA (stage 1 by design).
- Starting tree: one sample of the round-4 total-evidence posterior (shared clock), not a summary tree.
