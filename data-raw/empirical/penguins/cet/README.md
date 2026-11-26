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

The production runs (stage 1) were launched on 2026-11-11.  The package code
is main 48e6bf7, and the scripts are this branch at 921253b.  Each of the 16
chains starts from its own random sample of the round-4 total-evidence RevBayes
posterior (see **Starting tree** and **Production runs**).

## Files

| File | Role |
|---|---|
| `select-start-tree.R` | Picks the starting tree: the post-burn-in RevBayes sample closest, by quartet distance, to the quartet consensus of the pooled samples |
| `draw-start-trees.R` | Draws the 16 production start trees: seeded draws without replacement from the same post-burn-in pool, draw k for chain k = 4 (rep - 1) + cell |
| `prepare.R` | Converts `../penguins_morph.nex`, `../penguins_tips.csv`, `../penguins_constraints.csv` (empty) and a dated RevBayes tree (the sample chosen by `select-start-tree.R`) into `data/` and `revbayes/data/` |
| `data/characters.csv`, `data/partitions.csv` | 55 taxa x 201 characters (0-based states, empty = unknown); number of states K per character |
| `data/start-tree.nwk`, `data/start-tree-log.txt` | Starting time tree with the origin as root edge; every age `prepare.R` changed |
| `data/start-tree-selection.txt`, `data/quartet-consensus.nwk` | Record of the starting-tree choice (chain, iteration, distances) and the quartet consensus it was measured against |
| `data/start-trees/` | The 16 production start trees: `start-tree-<kk>.nex` (RevBayes samples), `start-tree-<kk>.nwk` and `-log.txt` (after `prepare.R`), `draws.csv` (chain, iteration and ages of each draw) |
| `data/start-tree-round1.nwk` | The earlier starting tree (from the round-1 MCC tree); `check-mkv.R`'s RevBayes oracle was computed on it |
| `common.R` | Model, priors, moves, data loading, and the two namespace rebindings (partitioned Mkv likelihood with one correction per missing-data pattern; per-sample recording and exact resume) |
| `check-mkv.R` | Fixed-tree check of the partitioned Mkv likelihood against RevBayes `coding = "variable"` on `data/start-tree-round1.nwk`, and its cost (#65) |
| `run-chain.R` | One chain of one cell, in checkpointed segments |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R` | Smoke run on the `test` partition and its checks (timing, acceptance, exact resume) |
| `bayes-factors.R` | Cross-model evaluation of recorded states (`eval`) and bridge-sampling Bayes factors (`combine`) |
| `check-labels.R` | Checks of chain labels (plain decimal notation, unparseable labels stop), of the Mkv-variant filter of `combine`'s stepping-stone (#75), and of `combine` on incomplete ladders and edges (#88), of the expected chain set (launch record, #101; #103 review), and of `combine` refusing a `cross.csv` that failed `eval`'s checks (#104); `smoke-results/check-labels.txt` |
| `check-eval.sh` | `eval`'s side of #104 on a copy of a real chain: no `cross.csv` unless the checks pass; `smoke-results/check-eval.txt` |
| `summarise.R` | Convergence (ESS, PSRF, ASDSF), clade and sampled-ancestor frequencies, crown/root/origin ages, CET and RevBayes side by side |
| `revbayes/penguins_morph_fbd.Rev`, `revbayes/run-hamilton.sh`, `revbayes/data/` | Morphology-only RevBayes comparator (UCLN, and strict as a pipeline check) on the project's RevBayes development build a2a8e2db (#66); `ML=true ML_MODE=stone` runs its stepping-stone one stone per Slurm task |
| `revbayes/prior-check/` | Data-free check of the comparator's tree moves on a2a8e2db: sampled-ancestor count against BEAST2 and v1.4.1 (#66) |
| `ml-ladder.sh`, `ml-root.sh`, `ml-summary.R`, `ml-helpers.R` | CET character-power ladders, their run root, and the common-footing marginal likelihoods against the RevBayes comparator ("Common-footing marginal likelihoods"); `ml-helpers.R` holds the summary's tested helpers (`tests/testthat/test-PenguinMlHelpers.R`) |
| `stage1-results.R`, `stage1-figures.R` | Stage-1 summaries and figures ("Stage-1 results") |
| `stage1-results/` | Stage-1 result tables, consensus trees, figures and Bayes-factor outputs |
| `s-prior-reweight.R` | Sensitivity of the Bayes factors to the prior on s ("Stage-1 results") |
| `rung-array.sh`, `data/start-trees/rungs/` | Ladder-rung arrays and the start trees of the dispersed rungs |
| `jobs/` | Slurm wrappers of the stage-1 jobs |
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
  can be checked by reweighting (Z'/Z = E_post[pi'(s) / pi(s)]).  The check
  (`s-prior-reweight.R`; "Stage-1 results") moves the gamma Bayes factor by
  -1.5 to +1.9 nats under four log-normal alternatives.  **Decision
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

**Production start trees (dispersed).** The quartet-nearest tree above is the
reference start tree, used by the smoke run and `check-smoke.R`.  The
production chains do not start from it.  Starting 16 chains from one tree
weakens every between-chain diagnostic (#89 review, R89-5), so each chain
starts from its own posterior sample:

    Rscript data-raw/empirical/penguins/cet/draw-start-trees.R <run root> <out dir> [n = 16] [seed = 20261111]
    Rscript data-raw/empirical/penguins/cet/prepare.R <out dir>/start-tree-<kk>.nex $(cat <out dir>/start-tree-<kk>-origin.txt) <prepared out>

`draw-start-trees.R` pools the same 7504 round-4 samples as
`select-start-tree.R` (25% burn-in, `shared_run13..16`).  It then draws 16 of
them uniformly without replacement (`set.seed(20261111); sample.int`).  Draw k
goes to chain k = 4 (rep - 1) + cell, with cells in the order b0g0, b0g1, b1g0,
b1g1, which is task k of `run-array.sh`.

`prepare.R` turned *W. manneringi* from a root sampled ancestor into a fossil
tip 0.01 Myr below the root in 8 draws (k = 1, 2, 3, 4, 8, 9, 10, 14).  In
k = 6 it raised one node by 0.0046 Myr; the other 7 draws were unchanged.

Ages are in Ma.  The start logPost is after `CetStart()` (below).

| k | cell | rep | chain | iteration | origin | root | crown | SAs (RevBayes) | `prepare.R` change | start logPost |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | b0g0 | 1 | shared_run16 | 56820 | 66.884 | 60.759 | 11.498 | 12 | W. manneringi root SA -> fossil tip | -2679.1 |
| 2 | b0g1 | 1 | shared_run16 | 50280 | 63.313 | 61.054 | 11.308 | 15 | W. manneringi root SA -> fossil tip | -2289.2 |
| 3 | b1g0 | 1 | shared_run15 | 134760 | 61.627 | 60.751 | 11.630 | 15 | W. manneringi root SA -> fossil tip | -2180.2 |
| 4 | b1g1 | 1 | shared_run13 | 46980 | 62.298 | 61.113 | 11.515 | 15 | W. manneringi root SA -> fossil tip | -2286.0 |
| 5 | b0g0 | 2 | shared_run15 | 141960 | 64.165 | 63.207 | 11.731 | 10 | none | -2724.1 |
| 6 | b0g1 | 2 | shared_run15 | 92280 | 63.951 | 61.790 | 11.512 | 21 | node 108: 39.6504 -> 39.6550 | -2214.4 |
| 7 | b1g0 | 2 | shared_run16 | 78060 | 64.848 | 63.490 | 12.097 | 21 | none | -2413.4 |
| 8 | b1g1 | 2 | shared_run16 | 103260 | 61.525 | 61.112 | 11.650 | 11 | W. manneringi root SA -> fossil tip | -2365.1 |
| 9 | b0g0 | 3 | shared_run13 | 42060 | 71.654 | 61.173 | 11.380 | 20 | W. manneringi root SA -> fossil tip | -2175.3 |
| 10 | b0g1 | 3 | shared_run13 | 134040 | 73.624 | 60.925 | 13.292 | 20 | W. manneringi root SA -> fossil tip | -2290.5 |
| 11 | b1g0 | 3 | shared_run14 | 65340 | 72.174 | 66.177 | 11.150 | 15 | none | -2845.2 |
| 12 | b1g1 | 3 | shared_run13 | 111900 | 84.975 | 61.263 | 12.346 | 9 | none | -2425.6 |
| 13 | b0g0 | 4 | shared_run13 | 137460 | 64.936 | 63.538 | 12.313 | 15 | none | -1989.7 |
| 14 | b0g1 | 4 | shared_run15 | 84240 | 63.487 | 60.698 | 12.188 | 12 | W. manneringi root SA -> fossil tip | -2525.9 |
| 15 | b1g0 | 4 | shared_run15 | 67980 | 66.367 | 62.145 | 12.319 | 11 | none | -2007.2 |
| 16 | b1g1 | 4 | shared_run13 | 74760 | 69.613 | 69.607 | 12.941 | 15 | none | -2243.3 |

The draws span origin 61.5-85.0 Ma, root 60.7-69.6 Ma, crown 11.15-13.29 Ma
and 9-21 sampled ancestors, against a single tree before.  `data/start-trees/`
holds every file; their md5s are recorded in each chain's checkpoint.

**Start parameters (`CetStart()`).**
- Each chain jitters the RevBayes-derived starting parameters by
  exp(N(0, 1)), including replicate 1, with seed 9000 + k.  It redraws while the
  start logPost is not finite.
- All 16 chains were finite at the first draw (logPost -2845 to -1990).
- `run-chain.R` records the start-tree file, its md5 and the start parameters
  in the checkpoint.  A resume stops if the file or its tip order differs.
- `bayes-factors.R eval` loads each chain's own start tree, checks its md5, and
  stops unless every state's tip labels equal the character rows.

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
`progress.txt` of each chain).  `CET_FNPR` (default 0, as in the stage-1
chains) adds `fnpr` moves, `ceiling(CET_FNPR x 0.15 x bifurcations)` per
iteration (#110).

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
Exact continuation also needs the same EvoTempo build for every segment.
From #110 (PR #119), builds that provide `.MkSiteLogLik()` compute the Mkv
correction in C++ (`PatternLogLik()`), which differs from the R
`ConstantPatternLogLik()` of older builds by round-off (about 1e-13 in the log
likelihood); the checkpoint does not record the build, so a chain resumed across
that change continues a valid chain, but not bitwise.

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
`combine` checks the chains it finds against the expected set (#101): the
launch record `<run root>/launch.csv` (one row per launched group of chains:
`cell`, `reps`, optional `power` and `mkv`; `ExpectedChains()` in `common.R`),
or, without one, the production launch of 4 cells x 4 replicates at power 1
(`CET_BF_REPS` changes the 4).  It stops, and reports no Bayes factor, if an
expected chain has no complete segment, if an expected power-1 chain has no
`cross.csv` (its eval failed or was not run), or if it finds a chain the record
does not list; before, a chain without `cross.csv` dropped out of the bridges
without a message.  The only opt-out is explicit and per chain:
`CET_BF_ALLOW_MISSING=b0g0_run3,...`, echoed as "Left out on request".
Since the #103 review, `combine` also deletes its outputs (`bayes-factors.csv`,
`model-probabilities.csv`, `closure.txt`, `chains-used.csv`) before anything
else, so a stop leaves none from an earlier run (R103-3); leaves a chain out of
the stepping-stone too, stopping if a ladder then lacks t = 1 (R103-4);
requires every chain it uses to be complete (last `progress.txt` line;
`CET_BF_ALLOW_INCOMPLETE=<labels>` names exceptions) and refuses a `cross.csv`
whose chain has changed since eval read it, or that does not record what eval
read (`evalSamples`, `evalSegMd5`, an md5 of the `seg*_samples.csv` files;
R103-2); stops on two directories for one chain, a non-integer `reps`, a
launched chain directory missing from the record, and a record with no chain
of the variant (R103-7 to R103-10); writes `chains-used.csv` (R103-12); and
names a column missing from `cross.csv` (R103-13).  A `cross.csv` written
before these fixes is refused: rerun `eval`.
With #104, `eval` writes `cross.csv` only after its checks pass, records its
tolerances (`evalTol`, `evalTolChar`), and `combine` refuses a `cross.csv` that
fails the stricter of these and its own `CET_BF_TOL`, `CET_BF_TOL_CHAR` (R103-6
is resolved by #104).  `check-labels.R` checks combine's side of all of these;
eval's side (the gated write, the record, an early stop) is checked by
`check-eval.sh` on a real chain (#105 review, R105-10).
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
    # launch record, one row per launched group (needed once rungs or power chains exist):
    printf 'cell,reps\nb0g0,4\nb0g1,4\nb1g0,4\nb1g1,4\n' > $W/runs/launch.csv   # + e.g. b0g0.5,4 per rung
    # per chain, then once:
    Rscript data-raw/empirical/penguins/cet/bayes-factors.R eval $W/runs/b1g1_run1
    Rscript data-raw/empirical/penguins/cet/bayes-factors.R combine $W/runs    # stops on missing chains or cross.csv
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
  Since #104, `eval` also recomputes the character log-likelihood
  (`checkChar`, against the chain's charLL).
  - It first deletes any earlier `cross.csv`, so no stop leaves one behind
    (#105 review, R105-1).
  - It writes `cross.csv` (temporary file, removed on exit, then rename)
    only when every state's |`check`| is within `CET_BF_TOL` (default 1e-4)
    and every |`checkChar`| within `CET_BF_TOL_CHAR` (default 1e-8).  The
    file records both tolerances.
  - `combine` refuses a `cross.csv` that fails these checks.  It also refuses
    one with a non-finite check, with no `checkChar`, `evalTol` or
    `evalTolChar` column, or with a recorded tolerance that is not a positive
    number, naming the chain (`check-labels.R` section 8, `check-eval.sh`).
  - Tolerances: the recomputed log posterior differs from the chain's by up to
    1.86e-5 (a 200 000-iteration rung) and 1.11e-5 (the 400 000-iteration
    production cells), about 1e-8 relative.  That is drift in the running
    value, hence a tolerance wider than the trilobite 1e-6.  charLL is
    recomputed exactly (max |`checkChar`| 5.0e-12 over the 40 chains), hence
    1e-8.
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

The first launch is recorded first, then the dispersed-start restart that
replaced it.

- Code (first launch): EvoTempo main **48e6bf7** (checkout `$W/src`, built into `$W/lib-48e6bf7`;
  `$W/out/build-sha.txt`), with this branch's start tree (`prepare.R` output as
  committed here); `$W = /nobackup/pjjg18/evotempo-penguin-cet`.
- Checks at 48e6bf7: `check-mkv.R` passed (round-1 tree; per-pattern lnL
  minus RevBayes `variable` -6.5e-5, 1.4e-5, -4.2e-5 at rates 0.005, 0.017,
  0.05, both engines).  `smoke.sh` (Slurm 19428287, `test`, from the round-4
  tree) and `check-smoke.R` (`smoke-results/check-smoke-48e6bf7.txt`): 0.30-0.36
  s per iteration, exact resume (maximum difference 0 over 10 samples x 23
  columns, trees identical), no non-finite value in any sample column, peak
  memory 169-171 MB.
- First-launch array: Slurm **19428297** (`evotempo-penguin-cet`, tasks 1-16 = 4 cells x 4
  replicates), `CET_NITER=200000, CET_THIN=100, CET_SEG=10000`, burn-in 2000,
  `CET_MKV` default (`pattern`), s ~ Gamma(0.5396, scale 0.3819) on (0, 3],
  origin U(0, 160).  Submitted through `$W/jobs/run-array-wrap.sh`, which adds
  the batch header Hamilton jobs need (`set +e; set +o pipefail; source
  /etc/profile || true; module purge || true; module load r/4.5.1`) and then
  sources `run-array.sh` with the same `#SBATCH` settings.

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
- **First launch, cancelled.** Array 19428297 started all 16 chains from the
  same tree, with jitter exp(N(0, 0.2^2)) and none for replicate 1 (R89-5).  It
  was cancelled at 4-7 of 20 segments.  Its runs and logs are kept in
  `$W/archive-19428297-common-start` (39 MB, without the saved-state files).
  The segment-1 figures above are from this launch.
- **Production: Slurm 19431480** (`evotempo-penguin-cet-dispersed`, tasks 1-16,
  launched about 15:37 GMT on 2026-11-11), with the dispersed start trees and
  `CetStart()` jitter described under "Starting tree".
  - Code: checkout of this branch at **921253b**, i.e. main fb24ed6's scripts
    (merged in 5ce1e31) plus this branch's fixes.  The package is unchanged from
    48e6bf7, so the build in `$W/lib-48e6bf7` is used.
  - Settings as above: `CET_NITER=200000, CET_THIN=100, CET_SEG=10000`,
    burn-in 2000, `CET_MKV=pattern`.
- **Smoke check of the restart** (cell (1, 1), chain k = 4; Slurm 19431233 tasks
  4-5 and test job 19431304):
  - 0.338-0.339 s per iteration; peak memory about 174 MB.
  - Segmented and uninterrupted runs are identical (maximum difference 0 over
    10 samples x 23 columns, trees identical), with no non-finite values.
  - A resume in a new process completed segment 3 at 0.327 s per iteration.
  - Resuming with a different start tree stops with an error.
  - `eval` on this chain: all 30 states pass the tip-label check, and the
    maximum |recomputed - logged logPost| is 5.1e-9.
- **Segment 1 of 19431480:** all 16 chains done, no non-finite values.
  - Speed: 0.288-0.417 s per iteration (cell means (0, 0) 0.306, (0, 1) 0.351,
    (1, 0) 0.308, (1, 1) 0.389).  That is about 304 core-hours, and 16-23 h per
    chain.
  - Acceptance: `slide` 2.8-4.2%, `saJump` 25-31%, `nodeAge` 42.5-51.5%, parameter
    sweep 18.5-40.8%.
- To continue the chains, resubmit:

      sbatch -J evotempo-penguin-cet-dispersed --array=1-16 \
        --export=ALL,CET_NITER=200000,CET_THIN=100,CET_SEG=10000 $W/jobs/run-array-wrap.sh


## Stage-1 results (400 000 iterations; provisional)

**Status: provisional.**
- **Code.** The Bayes factors use the code of #105 (`penguin-eval-check-104`;
  the 400 000-iteration numbers below were recomputed with its head at the
  time, ee4b364).  #103 and #105 are not merged.
- **Chains still running.** These chains are being extended:
  - production to 600 000 iterations;
  - rungs to 600 000;
  - four new rungs added.
  Every number here will be recomputed then (same scripts, no segment cap).
- **The 400 000-iteration state is reproducible** with
  `CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20` (stage1-results.R, bayes-factors.R).
  The run root must hold only the 40 chains of that state:
  `bash jobs/make-bf400k-root.sh $W/runs $W/bf400k-root` builds it.  The
  script makes links to the 16 production and 24 rung chains (not the
  gamma = 0.125 / 0.375 rungs), copies `stage1-results/launch.csv` (10
  groups x 4) as its launch record, and checks that every chain has its
  capped segments.  Then:
  - eval each chain as `jobs/bf400k-eval.sh` does;
  - combine on that root;
  - run `stage1-results.R $W/bf400k-root ...`.
- **Acceptance gate for the gamma Bayes factor** (#107 review):
  - every rung passes R-hat <= 1.01 and bulk ESS >= 400 (meanX included);
  - each bridge step has Kish ESS >= 100 in both directions;
  - the chain-half SE is < 0.5;
  - the small squares close within 3 SE.
  At 400 000 iterations the first two conditions fail (below).

### Runs

- **Extension to 400 000 iterations: Slurm 19448311** (`evotempo-penguin-cet-dispersed-ext400k`,
  tasks 1-16; `jobs/run-array-wrap.sh`).  The 16 production chains of
  19431480 (200 000 iterations) were continued from segment 21 with
  `CET_NITER=400000`.  Each task took about 24 h; peak memory was 169 MB.
  - Reason: at 200 000 iterations cells (0, 1) and (1, 1) failed the flags:
    (0, 1) treeLength R-hat 1.015, bulk ESS 262, and charLL bulk ESS 289;
    (1, 1) theta bulk ESS 358.  ASDSF was about 0.03 in every cell, against
    0.005 for RevBayes.
- **gamma rungs: Slurm 19448246-19448251** (`evotempo-penguin-cet-rung-b<beta>g<gamma>`;
  `jobs/rung-wrap.sh`, now `rung-array.sh` with `RUNG_START=production`).
  - Cells: gamma = 0.25, 0.5 and 0.75 at beta = 0 and at beta = 1, 4 replicates
    each, 200 000 iterations at power 1.  `CET_SEED` = 510000 + rep (b0g0.25),
    520000 + rep, 530000 + rep, and 610000 + rep, 620000 + rep, 630000 + rep at
    beta = 1.  About 23-27 h per chain.
  - **Shared starts (R107-15).**  Rung replicate r starts from:
    - the start tree of the gamma = 1 production chain with the same beta and
      replicate (k = 4 (r - 1) + 2 for beta = 0, + 4 for beta = 1);
    - **and** from `CetStart(r)`, i.e. the start parameters of production
      chain k = r.  `run-chain.R` uses index = rep for labels other than the
      four cells.

    The rungs' replicates are therefore not dispersed independently of the
    production chains.
  - Why rungs: at 200 000 iterations the gamma edges had one-sided gaps of 19.7
    (beta = 0) and 16.7 (beta = 1) nats, and the README overlap rule requires
    rungs.
- **Launched in the #107 fix round (2026-11-16):**
  - **Production extension to 600 000 iterations:** Slurm 19466805
    (`evotempo-penguin-cet-dispersed-600k`).
  - **Rung extension to 600 000 iterations** (resumed from segment 20, same
    code 921253b): Slurm 19472375 (b0g0.25), 19472376 (b0g0.5), 19472377
    (b0g0.75), 19472378 (b1g0.25), 19472379 (b1g0.5), 19472380 (b1g0.75)
    (`evotempo-penguin-cet-rung-<cell>-600k`).
  - **New rungs gamma = 0.125 and 0.375 at beta = 0 and 1**, 4 replicates,
    600 000 iterations, with dispersed starts:
    - Start trees: draws 17-32 of `draw-start-trees.R` (seed 20261111,
      n = 32; draws 1-16 are identical to the production draws, so no tree is
      shared), converted by `prepare.R` (`data/start-trees/rungs/`,
      `draws.csv`).
    - Start parameters: `CetStart(k)`, k = 17-32 (`CET_START_INDEX`, new in
      `run-chain.R`).
    - Jobs: Slurm 19472381 (b0g0.125, seed 540000 + rep), 19472382
      (b0g0.375, 550000 + rep), 19472383 (b1g0.125, 640000 + rep), 19472384
      (b1g0.375, 650000 + rep).
    - `rung-array.sh` with `RUNG_START=dispersed`, from this branch at
      4638d62.
    - `prepare.R` changes: six draws needed one.  In three
      (k = 18, 19, 22) *W. manneringi* went from a root sampled ancestor to a
      fossil tip; in three (k = 20, 25, 26) a node was raised by at most
      0.0074 Myr.
  - **Launch record:** `$W/runs/launch.csv` lists the 14 groups (production,
    6 rungs, 4 new rungs).  The 400 000-iteration record is
    `stage1-results/launch.csv`.
- **Stale outputs.**  `cross.csv` files and Bayes-factor outputs written on
  2026-11-13 while chains were still running were deleted before evaluation
  (`stage1-results/gate.txt`).  Every `cross.csv` written before the #105 fix
  round lacks `evalTolChar` and is refused by `combine`, so the 600 000
  re-evaluation must run on all chains.

### Convergence

`stage1-results.R` (Slurm 19472732, `CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20`, 66 min, 1.9 GB).
- **Burn-in and diagnostics.**  Burn-in is the first 25% of each chain's
  recorded samples.  R-hat and bulk/tail ESS use MkPrime's rank-normalised
  split versions (`versions.txt` records the md5 of the deparsed functions,
  whose source is in `mkprime-functions.R`).
- **Quantities diagnosed:** lambda, mu, psi, theta, s, rate, origin, crown and
  root age, nSA, charLL, treeLength, ucldSd, and the tempo-path state meanX
  and rootBin (R107-14).  logPost is not used.
- **Flags:** R-hat > 1.01 or bulk ESS < 400.
- **ASDSF** is over rooted clades with frequency >= 0.1.  The last column is
  the largest replicate range of the probability that one of the 18 baseline
  clades is monophyletic among extant taxa.

| analysis | chains x draws | max R-hat (parameter) | min bulk ESS (parameter) | min tail ESS | flagged (R-hat / bulk ESS) | ASDSF (max SDSF) | max replicate range of a baseline clade |
|---|---|---|---|---|---|---|---|
| (0, 0) | 4 x 3000 | 1.0090 (charLL) | 400.0 (charLL) | 1820 | charLL 1.009 / 400 | 0.0231 (0.098) | 0.20 |
| (0, 1) | 4 x 3000 | 1.0045 (treeLength) | 735.3 (treeLength) | 1756 | none | 0.0252 (0.133) | 0.30 |
| (1, 0) | 4 x 3000 | 1.0086 (charLL) | 1000.5 (treeLength) | 2781 | none | 0.0219 (0.090) | 0.20 |
| (1, 1) | 4 x 3000 | 1.0102 (charLL) | 703.1 (theta) | 1284 | charLL 1.010 / 735 | 0.0236 (0.104) | 0.24 |
| RevBayes strict | 4 x 4501 | 1.0003 (nSA) | 11144.1 (treeLength) | 14427 | none | 0.0056 (0.014) | 0.02 |
| RevBayes UCLN | 4 x 4501 | 1.0002 (rate) | 9715.6 (treeLength) | 13374 | none | 0.0052 (0.014) | 0.02 |
| b0g0.25 (rung) | 4 x 1500 | 1.0360 (s) | 97.9 (s) | 273 | theta 1.005 / 389; s 1.036 / 98; rate 1.017 / 239; charLL 1.030 / 151; meanX 1.014 / 193 | 0.0377 (0.177) | – |
| b0g0.5 (rung) | 4 x 1500 | 1.0132 (s) | 374.9 (treeLength) | 1201 | s 1.013 / 611; charLL 1.012 / 530; treeLength 1.010 / 375 | 0.0333 (0.103) | – |
| b0g0.75 (rung) | 4 x 1500 | 1.0076 (charLL) | 443.0 (treeLength) | 1096 | none | 0.0328 (0.159) | – |
| b1g0.25 (rung) | 4 x 1500 | 1.0211 (s) | 233.1 (theta) | 225 | theta 1.019 / 233; s 1.021 / 299; rate 1.010 / 392; charLL 1.017 / 344; treeLength 1.011 / 463; meanX 1.014 / 453 | 0.0348 (0.172) | – |
| b1g0.5 (rung) | 4 x 1500 | 1.0124 (treeLength) | 311.6 (treeLength) | 704 | theta 1.005 / 327; crownAge 1.012 / 338; treeLength 1.012 / 312 | 0.0325 (0.140) | – |
| b1g0.75 (rung) | 4 x 1500 | 1.0077 (charLL) | 363.4 (theta) | 808 | theta 1.006 / 363 | 0.0330 (0.123) | – |

- **Production cells.**  The scalar parameters pass, apart from two charLL
  values at the threshold: (0, 0) bulk ESS 399.99 and (1, 1) R-hat 1.0102.
  meanX and rootBin pass in every cell, with bulk ESS >= 2069.
- **Rungs (R107-1): 5 of the 6 rung cells fail the flags.**  The worst are
  b0g0.25 (s R-hat 1.036, bulk ESS 98; charLL 1.030 / 151; meanX 1.014 / 193)
  and b1g0.25 (theta 1.019 / 233; s 1.021; meanX 1.014).  Only b0g0.75
  passes.  The rungs mix slowly at gamma = 0.25-0.5, so they were extended
  and new rungs added (above).
- **Topology is not converged (#110).**  ASDSF is 0.022-0.025 in the cells and
  0.032-0.038 in the rungs, against 0.005 for RevBayes.
  - Within a cell, the probability of one of the 18 baseline clades differs
    between replicates by up to 0.20, 0.30, 0.20 and 0.24, against 0.02 for
    RevBayes.  Clades whose range exceeds 0.15: 3, 1, 2 and 3.
  - The largest between-replicate SD of such a probability is 0.098, 0.133,
    0.090 and 0.104, so its 4-replicate SE reaches 0.05-0.07 (R107-12).
  - Crown-age medians agree within 0.17 Myr between replicates.
- **Consistency** (`consistency-checks.csv`): for every CET chain, root and
  origin ages read from the trees match the logged values to 6e-6.
  Tree-based nSA is used for both programs, because RevBayes's
  `num_samp_anc` omits a sampled ancestor that is a child of the root.

### Posterior estimates

Pooled post-burn-in draws: median [2.5%, 97.5%]; P(lambda > mu) is the
posterior probability.  Ages are in Myr and rates per Myr.  `rate` is CET's
clock-rate parameter and RevBayes's `clock_morph`; under rate variation
(gamma = 1, UCLN) it is the location of a rate distribution.

| parameter | (0, 0) | (0, 1) | (1, 0) | (1, 1) | RevBayes strict | RevBayes UCLN |
|---|---|---|---|---|---|---|
| crownAge | 12.32 [10.70, 14.49] | 13.17 [11.21, 15.79] | 12.32 [10.68, 14.50] | 13.30 [11.21, 16.23] | 12.27 [10.64, 14.34] | 13.00 [11.09, 15.92] |
| rootAge | 62.13 [60.78, 67.04] | 62.04 [60.75, 66.95] | 62.18 [60.77, 67.16] | 62.21 [60.77, 69.05] | 61.58 [60.57, 66.44] | 61.39 [60.55, 65.35] |
| originAge | 64.89 [61.45, 75.38] | 64.71 [61.41, 75.62] | 64.96 [61.51, 75.86] | 65.19 [61.43, 79.19] | 64.41 [61.13, 75.01] | 64.01 [61.05, 74.13] |
| lambda | 0.502 [0.263, 1.01] | 0.443 [0.228, 0.914] | 0.517 [0.259, 1.16] | 0.562 [0.212, 1.38] | 0.488 [0.252, 0.991] | 0.476 [0.244, 1.01] |
| mu | 0.564 [0.266, 1.12] | 0.494 [0.215, 1.04] | 0.575 [0.254, 1.30] | 0.561 [0.185, 1.44] | 0.546 [0.245, 1.10] | 0.535 [0.234, 1.13] |
| turnover | 1.11 [0.832, 1.39] | 1.11 [0.786, 1.41] | 1.10 [0.829, 1.38] | 0.988 [0.790, 1.20] | 1.11 [0.815, 1.40] | 1.11 [0.802, 1.41] |
| diversification | -0.0600 [-0.200, 0.0692] | -0.0496 [-0.193, 0.0777] | -0.0563 [-0.213, 0.0725] | 0.0064 [-0.143, 0.102] | -0.0578 [-0.197, 0.0733] | -0.0567 [-0.203, 0.0757] |
| P(lambda > mu) | 0.188 | 0.244 | 0.208 | 0.552 | 0.206 | 0.212 |
| psi | 0.0648 [0.0327, 0.117] | 0.0712 [0.0363, 0.128] | 0.0647 [0.0326, 0.118] | 0.0668 [0.0341, 0.123] | 0.0671 [0.0339, 0.125] | 0.0702 [0.0346, 0.130] |
| samplingProportion | 0.103 [0.0327, 0.274] | 0.126 [0.0396, 0.334] | 0.101 [0.0289, 0.285] | 0.107 [0.0288, 0.337] | 0.110 [0.0341, 0.302] | 0.116 [0.0343, 0.323] |
| rate | 0.0200 [0.0162, 0.0248] | 0.0192 [0.0089, 0.0385] | 0.0200 [0.0161, 0.0246] | 0.0222 [0.0094, 0.0405] | 0.0202 [0.0163, 0.0250] | 0.0225 [0.0166, 0.0321] |
| theta | 0.0497 [0.0073, 0.338] | 0.0976 [0.0213, 0.383] | 0.0538 [0.0077, 0.374] | 0.145 [0.0175, 0.669] | – | – |
| s | 0.104 [0.0004, 0.999] | 0.830 [0.461, 1.43] | 0.0661 [0.0002, 0.474] | 0.741 [0.434, 1.33] | – | – |
| ucldSd | – | – | – | – | – | 0.649 [0.409, 0.966] |
| nSA | 13.00 [7.00, 20.00] | 14.00 [8.00, 21.00] | 13.00 [7.00, 20.00] | 13.00 [7.00, 20.00] | 14.00 [7.00, 21.00] | 14.00 [7.00, 21.00] |
| charLL | -1636 [-1650, -1625] | -1611 [-1627, -1597] | -1636 [-1649, -1625] | -1611 [-1627, -1598] | -1636 [-1650, -1625] | -1606 [-1621, -1592] |

**CET (0, 0) against RevBayes strict (R107-7).**  These are the same model up
to one difference: RevBayes allows a sampled ancestor at the root and
EvoTempo does not (#112).  RevBayes puts one there in 31.5% (strict) and
41.6% (UCLN) of posterior draws.
- So the comparison is made both unconditionally and on the RevBayes draws
  with no root sampled ancestor (`comparison.csv`, computed by
  `stage1-results.R`).
- The values are means of the 4 chain means, with z from chain-mean SEs; each
  side has 3 degrees of freedom, so |z| near 2 is weak evidence.

| quantity | CET (0, 0) | RevBayes strict, all draws (z) | RevBayes strict, no root SA (z) |
|---|---|---|---|
| lambda | 0.537 | 0.522 (+8.52) | 0.537 (+0.02) |
| mu | 0.598 | 0.580 (+6.45) | 0.597 (+0.46) |
| psi | 0.0673 | 0.0702 (-8.02) | 0.0668 (+1.24) |
| rate | 0.0201 | 0.0203 (-2.37) | 0.0201 (+0.84) |
| crownAge | 12.39 | 12.33 (+2.97) | 12.40 (-0.65) |
| rootAge | 62.60 | 62.14 (+23.89) | 62.65 (-2.44) |
| originAge | 65.79 | 65.34 (+8.80) | 65.90 (-2.02) |
| nSA | 13.20 | 13.71 (-11.69) | 13.11 (+1.61) |
| charLL | -1637 | -1636 (-1.03) | -1636 (-1.45) |
| P(lambda > mu) | 0.188 | 0.206 (-2.97) | 0.193 (-0.90) |

- **Unconditionally, the two differ:** root age by 0.46 Myr (z 24), nSA by
  0.5, and lambda, mu, psi by 2-5%.  The earlier statement that they
  "agree within Monte Carlo error" was wrong.
- **Given no root sampled ancestor, they agree:** |z| <= 1.6 for lambda, mu,
  psi, rate, crown age, nSA and charLL.  Root and origin age differ by 0.05
  and 0.11 Myr (z -2.4 and -2.0), which is borderline on these SEs.
- **lambda < mu is shared by the comparator.**  P(lambda > mu) is 0.188 in
  (0, 0) against 0.206 in RevBayes strict (0.193 given no root sampled
  ancestor).  It comes from the data and the FBD priors, not from CET.  The
  exception is cell (1, 1), at 0.552 (0.5525).
- **gamma = 1 cells against UCLN.**  These are different models (CET's
  autocorrelated tempo against uncorrelated lognormal rates), and their
  differences are expected:
  - crown age 13.17 / 13.30 against 13.00;
  - charLL -1611 against -1606;
  - P(*S. muizoni* is a sampled ancestor of extant *Spheniscus*) 0.59 against
    0.57.
- **theta and s.**  In the gamma = 0 cells, theta and s stay near their
  priors (s median 0.10 / 0.07).  With gamma = 1, s is 0.83 / 0.74 and theta
  0.098 / 0.145.
- **Quartet consensus trees** (`ConsTree::Quartet`, pooled draws) are in
  `stage1-results/consensus-*.nwk`, with split support in
  `consensus-splits.csv`.

### Bayes factors (provisional)

**Evaluation.**
- eval ran on the 40 chains through the link root `bf400k-root` (the 40
  chains plus the 400 000-iteration launch record), with
  `CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20` and #105 code (Slurm 19472716).
- All 40 tasks exited 0.  Max |check| was 1.86e-5 (b1g0.25_run1) and max
  |checkChar| 5.0e-12 (`eval-checks.csv`).
- combine (Slurm 19472717) found the 40 expected chains and reproduced the
  2026-11-16 08:21 Bayes factors exactly.

log BF of `to` over `from` (`bayes-factors.csv`, `bf-steps.csv`,
`small-squares.csv`).  The SE columns are:
- bridge: Fruhwirth-Schnatter;
- replicates: spread over replicate pairs;
- chain-half: spread over the 8 (replicate, half) units, which also sees slow
  mixing within chains (R107-2).

| from | to | log BF | SE bridge / replicates / chain-half | per replicate | min Kish ESS | one-sided bracket |
|---|---|---|---|---|---|---|
| b0g0 | b1g0 | -0.217 | 0.007 / 0.004 / 0.005 | -0.21556 -0.22141 -0.22593 -0.20547 | 941.9 | -0.26 to -0.05 |
| b0g1 | b1g1 | -0.584 | 0.052 / 0.058 / 0.060 | -0.44047 -0.54546 -0.63287 -0.71044 | 23.0 | -1.17 to 0.07 |
| b0g0 | b0g1 | 9.905 | 0.224 / 0.293 / 0.210 | 10.1857 9.1354 9.9518 10.5059 | 6.7 | 8.20 to 12.27 |
| b1g0 | b1g1 | 9.800 | 0.144 / 0.057 / 0.193 | 9.7350 9.7128 9.9650 9.8027 | 4.5 | 9.30 to 10.22 |

Per bridge step:

| step | log BF | SE | Kish ESS forward / reverse (of n0 / n1) | one-sided estimates | gap |
|---|---|---|---|---|---|
| b0g0 -> b1g0 | -0.217 | 0.007 | 941.9 / 1703.0 (of 12000 / 12000) | -0.26, -0.05 | 0.21 |
| b0g1 -> b1g1 | -0.584 | 0.052 | 34.9 / 23.0 (of 12000 / 12000) | -1.17, 0.07 | 1.25 |
| b0g0 -> b0g0.25 | 1.702 | 0.183 | 18.8 / 516.2 (of 12000 / 6000) | 0.64, 1.83 | 1.19 |
| b0g0.25 -> b0g0.5 | 5.227 | 0.109 | 9.7 / 9.1 (of 6000 / 6000) | 4.58, 6.30 | 1.72 |
| b0g0.5 -> b0g0.75 | 2.132 | 0.059 | 9.7 / 6.7 (of 6000 / 6000) | 2.15, 2.78 | 0.64 |
| b0g0.75 -> b0g1 | 0.845 | 0.036 | 24.2 / 87.4 (of 6000 / 12000) | 0.84, 1.35 | 0.51 |
| b1g0 -> b1g0.25 | 0.347 | 0.037 | 257.0 / 1023.9 (of 12000 / 6000) | 0.23, 0.37 | 0.14 |
| b1g0.25 -> b1g0.5 | 3.828 | 0.130 | 4.5 / 90.7 (of 6000 / 6000) | 3.94, 4.06 | 0.12 |
| b1g0.5 -> b1g0.75 | 3.712 | 0.041 | 89.6 / 9.9 (of 6000 / 6000) | 3.38, 3.75 | 0.37 |
| b1g0.75 -> b1g1 | 1.913 | 0.028 | 137.3 / 121.0 (of 6000 / 12000) | 1.75, 2.04 | 0.29 |

- **Closures.**  The square closes to 0.262 nats (SE 0.272).  The four small
  squares close to 0.061, 0.158, 0.064 and -0.021, all within 3 SE.
  - The closures cannot see a bias common to both gamma edges, because those
    edges enter with opposite signs (R107-4).
  - The overlap diagnostics are therefore the check that matters for gamma.
- **gamma.**  The data favour gamma = 1 over gamma = 0 by **about 10 nats**
  (9.90 at beta = 0, 9.80 at beta = 1).
  - Error: the chain-half SE is about 0.2 (0.21, 0.19).  Since the rung
    chains fail the flags, an SE of **about 0.3** is the honest quote; the
    bridge SE (0.14 at beta = 1) understates it.
  - **Provisional**, because overlap is poor.  On 6 of the 8 gamma steps the
    smaller Kish ESS is 4.5-24, against the gate's 100.
  - The summed one-sided estimates (8.20 to 12.27 at beta = 0; 9.30 to 10.22
    at beta = 1) are a diagnostic of that overlap, not an interval.  Halves of
    each chain show no sample-size bias (#107 review, R107-3).
- **What gamma measures (R107-5).**  In cell (0, 1), beta = 0, so the tempo
  process drives only the character rates (spec section 7.2: an
  autocorrelated relaxed clock whose log-rate follows a mean-reverting
  process, on a constant-rate FBD tree).
  - The gamma Bayes factor is therefore evidence for **tempo-modulated,
    autocorrelated rate variation among lineages**, not for covariation.
  - Whether that beats *uncorrelated* rate variation (the RevBayes UCLN
    comparator) is the common-footing question of #111.
- **Covariation is the beta contrast.**  Given gamma = 1 it is **-0.58**
  (SE about 0.06; chain halves show no drift), and given gamma = 0 it is -0.22.
  - Coupling diversification to the tempo process is not supported.  The
    gamma edges are the same at both beta (9.90, 9.80).
  - At intermediate gamma the beta edges are -1.63, -3.19 and -1.67.
  - Model probabilities (equal prior weight): (0, 1) 0.61, (1, 1) 0.39, the
    gamma = 0 cells < 1e-4.  This split is the least robust number (about
    +-0.1).
- **Sensitivity to the prior on s (R107-16; the open decision under "Model and
  priors").**  The gamma Bayes factor depends directly on the s prior,
  because (0, 0) does not depend on s.  (The (1, 0) evidence does: there s
  enters the tree prior through the tempo coupling, with log E[w] up to
  -0.64 under LN(log 0.5, 1).)
  - Method: reweighting the posterior draws, Z'/Z = E_post[pi'(s)/pi(s)]
    (`s-prior-reweight.R`; Kish ESS of the weights >= 2252 of 12000).  The
    per-chain values agree within 0.09 nats, except under LN(log 0.2, 0.5) in
    (1, 1), where they spread by 0.26.  It gives the shift in the gamma Bayes factor
    under log-normal alternatives on (0, 3].

  | alternative prior on s | shift, beta = 0 | shift, beta = 1 |
  |---|---|---|
  | log-normal with the Gamma's median (0.10) and 97.5% quantile (1.0) | -0.39 | -0.44 |
  | LN(log 0.2, 1) | +0.36 | +0.47 |
  | LN(log 0.2, 0.5) | -1.53 | -1.20 |
  | LN(log 0.5, 1) | +1.33 | +1.87 |

  The sign and order of magnitude of the gamma Bayes factor are robust.  Its
  value moves by up to about 2 nats with the s prior, which is more than its
  Monte Carlo error.  The beta Bayes factor at gamma = 1 moves by at most 0.3.

### Files

| File | Role |
|---|---|
| `stage1-results.R` | Convergence (cells and rungs, meanX, rootBin), ASDSF, per-chain and pooled summaries with P(lambda > mu), the CET-RevBayes comparison with and without root sampled ancestors, quartet consensus, consistency checks, versions; `CET_MAX_SEG` / `CET_MAX_SEG_RUNG` caps |
| `stage1-figures.R` | `fig-posteriors.png`, `fig-b0g0-strict.png`, `fig-traces.png` |
| `s-prior-reweight.R` | Sensitivity of the Bayes factors to the prior on s, by reweighting |
| `rung-array.sh` | Rung arrays (`RUNG_START=production` or `dispersed`) |
| `data/start-trees/rungs/` | Start trees of the gamma = 0.125 and 0.375 rungs (draws 17-32) |
| `jobs/` | The Slurm wrappers used for every job quoted here (R107-9) |
| `stage1-results/` | The outputs quoted above (400 000 iterations; `versions.txt`), the Bayes-factor tables, `bf-steps.csv`, `small-squares.csv`, eval checks, gate and launch records, `s-prior-reweight.csv`.  The 9 MB `samples.csv.gz` stays on Hamilton (`$W/stage1-results-400k-v3/`) |

## Morphology-only RevBayes comparator

No morphology-only baseline existed (the total-evidence baseline runs shared
and separate clocks with DNA), so `revbayes/` adds one, matched to the CET runs
in data, character model, tree prior, priors on lambda, mu, psi, tip ages and
starting tree, and differing in the clock: UCLN (non-centred), and strict
as the RevBayes counterpart of cell (0, 0), a check of the pipeline (the two
should agree on trees and FBD rates).  One further difference is structural:
RevBayes allows a sampled ancestor at the root, and EvoTempo does not (#112;
"Deviations").  Smoke (600 iterations incl. 100
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

**Rerun (launched 2026-11-11).**
- **Build and jobs.**  The comparator was rerun on RevBayes development
  **fd807d10** (1.5.0-preview; `/nobackup/pjjg18/revbayes-dev-fd807d10`),
  the project's current build, rather than a2a8e2db.  Outputs are in
  `/nobackup/pjjg18/evotempo-penguin-cet/rbmorph-fd807d10/{ucln,strict}_run1..4`,
  from Slurm 19429797-19429800 (`evotempo-penguin-rbmorph-ucln-1..4`) and
  19429801-19429804 (`evotempo-penguin-rbmorph-strict-1..4`).
- **Length.** Each chain ran 600 000 generations.  Strict run 1 took 4.3 h;
  UCLN runs 1 and 2 took 44 and 51 h.
- **Results** are under "Stage-1 results".  Strict agrees with CET (0, 0)
  given no sampled ancestor at the root, the one difference between the two
  (#112).

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

## Common-footing marginal likelihoods (stage 1)

**Question.**
- **What the gamma Bayes factor measures.**  The data favour gamma = 1 over
  gamma = 0 by about 10 nats (Bayes factors in "Stage-1 results", #107).  In
  cell (0, 1), beta = 0, so the tempo process drives only the character rates.
  That Bayes factor is therefore evidence for among-lineage rate variation
  that is autocorrelated along the tree (tempo is a CTMC on the tree).  It
  is not evidence for covariation between diversification and character change.
- **Covariation is the beta contrast.**  Within CET, beta = 1 against
  beta = 0 tests covariation, and it has no support (log BF -0.22 at
  gamma = 0, -0.58 at gamma = 1).
- **What this comparison asks.**  Does CET's autocorrelated rate variation
  fit the characters better or worse than the uncorrelated lognormal rate
  variation of the RevBayes UCLN comparator?  Equivalently, is the gamma
  Bayes factor specific to autocorrelated tempo, or would any rate variation
  earn it?  CET (0, 1) against UCLN answers this.  CET (1, 1) against UCLN adds
  the covariant tree prior, whose contribution beta has already measured.
  This needs marginal likelihoods on a common footing.

**What each program estimates.**
- Both programs anneal only the character likelihood:
  - CET's character-power chains (`CET_POWER`, `ml-ladder.sh`) raise
    `PartitionedMk()` to the power t;
  - RevBayes's `powerPosterior()` raises only the clamped `phy_morph` nodes.
- The FBD tree (with its tip-age intervals), the clock and the FBD parameters
  stay at their prior.  Each program therefore estimates
  log p(characters | tip-age intervals, model) = log E_prior[L(characters)],
  under its own normalised prior.
- The two tree densities differ by -log n! and the 1/width of the tip-age
  intervals, which are constants and cancel in that normalisation.  The
  orientation factor 2^|V| (`TempoMCMC()` targets the unoriented labelled
  tree) is not constant: |V| = n - nSA - 1 varies with the number of sampled
  ancestors.  It does not cancel, but both programs carry it identically
  (`R/MCMC.R`; RevBayes `lnProbTreeShape`), so it is shared.
- This is why the full-data evidences of "Bayes factors" are not comparable
  across programs, but these are.

**The remaining difference: a sampled ancestor at the root.**  EvoTempo's
samplers do not allow a sampled ancestor at the root; RevBayes does.  At the
stage-1 posterior, RevBayes puts one there in 31.5% (strict) and 41.6% (UCLN)
of draws.  Suppose the two priors agree on every other tree; the stage-1
posteriors support this (root age given no root sampled ancestor: 62.60 in
CET (0, 0) against 62.65 in RevBayes strict).  Then CET's prior is RevBayes's
conditioned on A = {no sampled ancestor at the root}, and

    log Z_CET = log Z_RB + log P_RB(A | characters) - log P_RB(A).

Here P_RB(A) comes from the trees of RevBayes's power-0 stones (prior), and
P_RB(A | characters) from the comparator's posterior trees.  The other
known difference is CET's truncation of its priors: lambda, mu, psi and rate at
+-4 sdlog (6e-5 of prior mass each), s on (0, 3]; theta has no RevBayes
counterpart.  Together these are below 1e-3 nats.

**Anchor.**  CET (0, 0) and RevBayes strict are the same model on that support.
Their log Z must therefore differ by the predicted root-SA term, up to
Monte Carlo error.  If they do, the comparison below is on a common footing.
If they do not, the gap is reported with its likely causes.  Those include
the root sampled-ancestor convention, the epsilon conventions for zero-length
branches and tip ages, and the FBD conditioning.

**Comparison.**
- CET (0, 1) = ladder (0, 0) + bridge (0, 0) -> (0, 1).  At fixed beta, the
  prior on (tree, paths, parameters) does not depend on gamma.  The bridge's
  full-data Bayes factor along a gamma edge is therefore a ratio of character
  marginals given tip ages.
- CET (1, 1) = ladder (1, 0) + bridge (1, 0) -> (1, 1).  A beta edge changes
  the tree prior, so the (0, 0) ladder cannot reach beta = 1; (1, 0) needs
  its own ladder.
- Both are compared with RevBayes UCLN restricted to CET's support (UCLN
  log Z plus its own root-SA term).  The unrestricted values are also reported.
- How to read the result:
  - CET (0, 1) - UCLN > 0: autocorrelated tempo fits better than uncorrelated
    lognormal rates.
  - About 0 (within SE): the gamma Bayes factor reflects rate variation in
    general.
  - < 0: the uncorrelated model fits better.
  - Neither sign speaks to covariation, which is the beta contrast.

**Design.**
- **Powers.** t_k = (k/64)^5, k = 0..64, i.e. RevBayes's `cats = 65`,
  `alpha = 0.2` schedule ((65 - k)/64)^5.  The two programs agree closely,
  not exactly: CET's powers are exact dyadic rationals, while RevBayes
  computes its powers by bisection in `Beta::quantile`
  (`DistributionBeta.cpp`), which matches to about 1e-6 relative and gives
  1.2e-302 for power 0.  The effect on log Z is below 1e-3 nats.
  The estimates on every second power show the discretisation error.
- **RevBayes.** fd807d10, strict and UCLN, 4 replicates, one Slurm task per
  stone (`runOneStone()`, then `summarize()`).
  - Generations per stone: strict 40 000, UCLN 20 000.  Each stone also has
    a 10 000-generation pre-burn-in with tuning, and 25% of its samples are
    discarded.  Samples every 10 generations.
  - Parameter and tree monitors on stones 1 (power 1) and 65 (power 0),
    every 100 generations.
  - Seed rep x 7919 + 202 + 100003 x stone.
  - Directories: `$W/rbml-fd807d10/<clock>_run<rep>/`.  Jobs:
    `evotempo-penguin-rbml-<clock>-<rep>` (arrays 1-65) and
    `...-summarize` (afterok).
- **CET.** Cells (0, 0) and (1, 0), the 64 powers t < 1, 4 replicates each
  (Slurm arrays 1-128 for replicates 1-2, then 129-256 for replicates 3-4).
  - 100 000 iterations each, plus 2000 burn-in; thinning 100; 25% of
    samples discarded.  The t = 1 rung is the cell's production chains,
    linked into the ladder root by `ml-root.sh`.
  - Start tree: the production start tree of the same cell and replicate.
    Seed `ML_SEED` + 100 j + rep (700000 for (0, 0), 800000 for (1, 0)).
  - Start parameters: `CetStart(rep)`, not production's `CetStart(k)`, so
    both cells share them at a given replicate.  Every power of a replicate
    starts from the same tree and parameters.
  - Replicates 5 and above need start trees 17 and above, which are not yet
    drawn (`ml-ladder.sh` header).
  - Code: the production scripts and build (921253b, `lib-48e6bf7`).
  - Root `$W/ml-runs`.  Jobs: `evotempo-penguin-cet-ml-<cell>`.
- **Estimation.** `ml-summary.R` applies the same estimators to both programs:
  - stepping-stone with `bayes-factors.R`'s `SteppingStone()` variance,
    per replicate and pooled.  The replicate is the unit: the SE of the
    pooled CET estimate is max(estimator SE, sd/sqrt(n)) over the n
    replicates (`CetPooledSe()`), as for the RevBayes mean and the bridge
    edges;
  - every replicate's ladder ends in all of the cell's production chains at
    t = 1.  Stepping-stone does not use those draws, but path sampling does,
    so the per-replicate path-sampling estimates are correlated;
  - it stops if a ladder lacks the t = 1 rung or a power, or if a chain is
    missing (dangling link) or incomplete (`ML_ALLOW_PARTIAL=1` warns instead);
  - path sampling (trapezoid) as a second estimator.
  - It also computes P(root SA) under the prior and the posterior, the anchor,
    and the comparison.
  - It checks the two programs' priors directly: CET's t = 0 chains against
    RevBayes's power-0 stones given no root sampled ancestor, on lambda, mu,
    psi, origin, root and crown age and nSA (`prior-check.csv`).
- **Target:** SE < 0.5 nats on each log Z, judged by both the estimator SE
  and the replicate spread.
- **Smoke tests** (Hamilton, 16 Nov):
  - The per-stone RevBayes path ran with 5 stones of 200 generations.  The
    powers are 1, (3/4)^5, (2/4)^5, (1/4)^5 and 0.  The power-0 and power-1
    stones wrote their logs and trees, and `summarize()` merged the stone files.
  - The CET ladder wrapper gave labels `b0g0_t0_run1` and
    `b0g0_t0.000000029802322387695312_run1`.  At t = 0 the logged charLL is
    the unpowered value (-3200 to -4500 under the prior).
  - `ml-summary.R`'s stepping-stone on the RevBayes smoke stones agrees with
    RevBayes's own `steppingStoneSampler()` to 1e-4 nats.

Results: pending (this PR is updated when the chains finish).

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
- Starting trees: each production chain starts from its own random sample of
  the round-4 total-evidence posterior (shared clock), not from a summary tree.
- No sampled ancestor at the root: EvoTempo's samplers exclude it
  (`.NoRootAncestor`), while the RevBayes comparator allows it (#112).  At the stage-1 posterior, RevBayes puts one
  there in 31.5% (strict) and 41.6% (UCLN) of draws.  Excluding it shifts
  root and origin age (CET's are 0.46 Myr older unconditionally) and also
  lambda, mu, psi and nSA ("Stage-1 results").  Given no root sampled
  ancestor, CET (0, 0) and RevBayes strict agree.
