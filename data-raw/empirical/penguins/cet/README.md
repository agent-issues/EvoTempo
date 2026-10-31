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

The full CET runs are **not** launched by this PR.  They need PRs #42 and #43
merged (this branch is stacked on #43) and, for the headline, the starting
tree from the round-2 total-evidence RevBayes baseline (see **Starting tree**).

## Files

| File | Role |
|---|---|
| `prepare.R` | Converts `../penguins_morph.nex`, `../penguins_tips.csv`, `../penguins_constraints.csv` (empty) and a RevBayes MCC tree into `data/` and `revbayes/data/` |
| `data/characters.csv`, `data/partitions.csv` | 55 taxa x 201 characters (0-based states, empty = unknown); number of states K per character |
| `data/start-tree.nwk`, `data/start-tree-log.txt` | Starting time tree with the origin as root edge; every age `prepare.R` changed |
| `common.R` | Model, priors, moves, data loading, and the two namespace rebindings (partitioned Mkv likelihood; per-sample recording and exact resume) |
| `run-chain.R` | One chain of one cell, in checkpointed segments |
| `run-array.sh` | Slurm array for the production chains (16 tasks = 4 cells x 4 replicates) |
| `smoke.sh`, `check-smoke.R` | Smoke run on the `test` partition and its checks (timing, acceptance, exact resume) |
| `bayes-factors.R` | Cross-model evaluation of recorded states (`eval`) and bridge-sampling Bayes factors (`combine`) |
| `summarise.R` | Convergence (ESS, PSRF, ASDSF), clade and sampled-ancestor frequencies, crown/root/origin ages, CET and RevBayes side by side |
| `revbayes/penguins_morph_fbd.Rev`, `revbayes/run-hamilton.sh`, `revbayes/data/` | Morphology-only RevBayes v1.4.1 comparator (UCLN, and strict as a pipeline check) |
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
delta = 0, removal 0, equilibrium start-bin prior.  Characters: Mkv per
K-partition, no among-character rate variation.  Time in Myr.

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

## Starting tree

The total-evidence baseline found that random FBD starting trees trap the
chains (crown 3-7 Ma, crown fossils outside the crown), so every chain starts
from the RevBayes total-evidence MCC tree:

    Rscript data-raw/empirical/penguins/cet/prepare.R <mcc.nex> <origin age>

`prepare.R` strips annotations, puts extant tips at 0, reads zero-length
fossil tips as sampled ancestors, clamps fossil (and sampled-ancestor node)
ages into their intervals, makes parents at least 0.01 Myr older than their
children, turns a sampled ancestor that cannot stay one (e.g. at the root,
which the samplers do not allow) into a fossil tip, and adds the origin as
root edge.  The committed `data/start-tree.nwk` comes from the **round-1**
MCC tree (`/nobackup/pjjg18/evotempo-penguins-rb1.4.1/summary_round1/mcc_shared.nex`,
shared clock; origin = posterior median 64.39 Ma): root 60.99 Ma,
27 sampled ancestors, one change (*Waimanu manneringi*, a sampled ancestor at
the root, became a fossil tip at 60.98 Ma).  Before the production runs,
rerun `prepare.R` on the round-2 MCC tree once the baseline sub-agent has
summarised round 2.  Fallback for testing without a baseline: the BEAST2 XML
starting tree of the baseline (`../revbayes/data/start_tree.nwk` on branch
`penguins-revbayes`), or any dated tree with these 55 tips, through the same
`prepare.R` path.

## Moves and chain mechanics

`CetMoves()` = `TempoMoves()` (with gamma = 1: 14 `node`, 7 `tip`,
21 `branch`, 1 `shift`; with gamma = 0: 1 `tree`, 2 `node`, 2 `tip`,
2 `branch`, 1 `joint`, 1 `shift`) plus tree moves in proportion to the tree:
`nodeAge` 7 (0.25 x bifurcations), `fossilAge` 9 (0.25 x fossils), `slide` 5,
`saJump` 4, `originAge` 1; `CET_TREE_MOVES` scales them.

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

## Running on Hamilton

    W=/nobackup/pjjg18/evotempo-penguin-cet     # src/ = checkout, lib-<sha>/ = EvoTempo build
    cd $W/src
    Rscript data-raw/empirical/penguins/cet/prepare.R <round-2 mcc.nex> <origin>
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
  computable at every recorded state (`eval`; the chain's own value is
  checked against its logPost: residuals 5e-12 to 6e-9 in the smoke run).
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
Pilot crown ages (provisional, single short chains): 12.2-13.5 Ma, against
12.5 Ma in the round-1 total-evidence baseline.

Projection (provisional: ESS from single 3000-iteration chains understate
slow drifts).  The slowest quantities reach about 2 ESS per 1000 iterations.
A pooled ESS of 200 per quantity (50 per chain) then needs about 25 000
iterations per chain; with a 4x margin, 200 000 iterations per chain at
0.19 s = 10.6 h.

| Item | Core-hours |
|---|---|
| Production: 4 cells x 4 chains x 200 000 iterations | 170 |
| `bayes-factors.R eval` (~0.02 s per state for all four model points) | < 2 |
| Ladder rungs, if needed: 4 chains x 10.6 h per rung, up to 3 per edge, 4 edges | 0-510 |
| Character-power stepping-stone, best cell only (16 powers x 2 chains x 50 000 it.) | ~85 |
| RevBayes comparator (running): UCLN 4 x 300 000 it. at ~0.43 s (first hour; 0.23 s in the smoke), strict 4 x 300 000 at ~0.05 s | ~145 + ~17 |
| RevBayes `powerPosterior` (UCLN, 32 powers x 10 000 it.) | ~20 per run |

Wall time: one 3-day job slot per chain suffices (10.6 h); peak memory in the
pilot was 178-181 MB per chain.

## Morphology-only RevBayes comparator

No morphology-only baseline existed (the total-evidence baseline runs shared
and separate clocks with DNA), so `revbayes/` adds one, matched to the CET runs
in data, character model, tree prior, priors on lambda, mu, psi, tip ages and
starting tree, and differing only in the clock: UCLN (non-centred), and strict
as the RevBayes counterpart of cell (0, 0), a check of the pipeline (the two
should agree on trees and FBD rates).  Smoke (600 iterations incl. 100
tuning): strict 25 s, UCLN 137 s.  Four chains of each, 300 000 iterations
(sample every 100, checkpoints every 5000), were launched on 2026-10-31 under
`/nobackup/pjjg18/evotempo-penguin-cet/rbmorph/<clock>_run<k>/`.

## Engine gaps

Workarounds live in `common.R`; the issues describe the engine support that
would replace them.

- Character model: per-character number of states, Mkv correction, ambiguous
  states, among-character rate variation: <ISSUE_CHARS>.
- `TempoMCMC()` output and control: per-sample trees and states, character
  log-likelihood column, character power, checkpoint/resume, bounded step
  adaptation, errors other than `infeasibleProposal` in proposals, topology
  mixing, joint or transformed priors: <ISSUE_MCMC>.
- DNA (stages 2-3): <ISSUE_DNA>.

## Deviations from the source analysis

- 32 ambiguous cells as unknown; character 81 dropped (201 characters).
- No among-character rate variation (source: Gamma(4)).
- Independent priors on lambda, mu, psi instead of the (d, turnover, sampling
  proportion) prior; s and the clocks keep the source priors (s bounded at 3).
- No DNA (stage 1 by design).
- Starting tree from the round-1 MCC tree, pending round 2.
