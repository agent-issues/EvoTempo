# Mixing of `TempoMCMC()` with many characters (issue #19)

Scripts in this folder diagnose why `TempoMCMC()` mixes slowly once characters are present,
and prototype new moves outside the engine. Every prototype calls engine internals through
`:::`. No file under `R/` or `src/` is changed. This note specifies the moves that worked,
so that the engine owner can implement them.

## Summary

1. **Most of the CPU goes on moves that are almost never accepted.** With 30 or more
   characters, a fresh whole-tree FFBS draw (`tree`, and the path part of `joint`) changes
   the character log-likelihood by 5.6-106 nats (sd over draws at fixed parameters). The
   frozen-E importance weight barely varies (sd 0.02-0.04): the frozen proposal is
   accurate, but it ignores the characters. On the large tree `tree` and `joint` are
   accepted at rates of order 1e-4 at every character count, yet `joint` takes 61-67% of
   the CPU time and `tree` another 6-7%.
2. **The character rate `c` sits on a ridge with the overall tempo level.** The
   characters fix `c * clock_e` on each branch, and the tree fixes `lambda * e^(beta x)`.
   Given the paths, log `c` has sd 0.013-0.04 at 100-300 characters, while its posterior
   sd is 0.07-0.17. Chains without a move along this ridge never cross it. On the large
   tree, `base` chains have a posterior sd of log `c` of 0.02-0.03, against 0.07-0.12 with
   the shift move: the current sampler under-disperses `c` and overstates its ESS.
3. **`theta` (and to a lesser degree `s`) are tied to the fine structure of the paths.**
   Given the paths, log `theta` has sd 0.054-0.090 (the jump count, N_J = 126-320, fixes
   `theta`), while its posterior sd is about 0.3-0.4. Paths are refreshed only by local
   moves (6 per iteration, against 60-210 branches) because the global moves are
   rejected. `theta` is therefore Gibbs-coupled to a slowly refreshed latent field.
4. **What to implement.** This list is in order of value per unit of engine work:
   1. the **shift move** (below), which is exact, costs one parameter proposal, and
      leaves the character likelihood unchanged. It is accepted 24-42% of the time;
   2. a **move schedule that adapts to the characters**: drop `tree` and `joint` (or run
      them with a probability fixed at the end of burn-in, set from their burn-in
      acceptance), and spend the time on local path moves scaled to the tree, about
      0.25 x (number of branches) per iteration;
   3. optionally, for `theta`: the **anchored joint move** (every-slice anchors) or the
      **NCMC move**. On the small tree these give 2-5x more `theta` ESS per iteration
      than (1)+(2). Per unit of CPU, they do not beat (1)+(2) in this prototype.
5. **Gains of (1)+(2) (`L25S`) over the current schedule (`base`)**, in ESS per CPU hour
   (tables below): on the small tree, `c` gains 10-66x, `lambda` 1-20x, the total clock
   8-78x, the character log-likelihood 1-4x and the root bin 1.6-18x. On the large tree,
   `c` gains 3-9x, the clock 2-17x and the root bin 0.4-10x. `theta` gains 1.4-5x, but on
   the large tree `theta` is not converged under any sampler tested in 3000-5000
   iterations (split R-hat and sampler-dependent posterior sd below), so those `theta`
   numbers are **provisional**.

## Set-up

- `common.R`: `Driver()` copies the loop of `TempoMCMC()` on a fixed tree, without tree
  moves. It runs path moves, `joint`, any extra moves, then the parameter sweep, and
  adapts every 20 burn-in iterations. With no extra moves it reproduces `TempoMCMC()`
  draw for draw: a seeded comparison is identical under both engines
  (`run.R check`; checked on Hamilton with the C++ engine at a2c3aa1, and locally with the
  R engine). The candidate moves are `ShiftMove()`, `AnchoredJoint()` and `NcmcMove()`.
- `run.R`: the trees (`sim`), chains (`chain`), probes (`probe`) and oracles (`oracle`).
  Configurations are listed in `configs`. `summary.R` writes the CSVs in `results/`, and
  `report.py` writes `results/tables.md`.
- Model and priors are those of rung 7 (#7): cell beta = gamma = 1, K = 16, z in [-3, 3]
  (bin width h = 0.4), rho = 0.5, `nSlice` = 16. Priors are log-normal, truncated at
  +/-4 sdlog; psi uses the "high" level.
- Two fixed true trees were simulated from the rung-7 simulator with seed 19019:
  - **small**: crown 8 Myr, 17 extant tips, 9 fossil tips and 19 sampled ancestors;
    truth lambda 0.60, mu 0.15, psi 0.26, theta 0.25, s 0.41, c 0.073;
  - **large**: crown 11 Myr, 75 extant tips, 31 fossil tips and 62 sampled ancestors;
    truth lambda 0.54, mu 0.19, psi 0.21, theta 0.080, s 0.40, c 0.135.
- 300 characters were simulated on each tree. Analyses use the first 30, 100 or 300, so
  the character sets are nested. All six parameters are free.
- Chains start at the true parameters with an FFBS draw of the paths. Each chain runs
  2000 burn-in and 5000 sampling iterations (NCMC chains: 1000 + 3000), with 1000 stored
  draws, and 2 replicates per cell. ESS is Geyer's initial monotone sequence estimator,
  as in rungs 6-7. CPU time is the time spent inside the moves, scaled to the sampling
  phase.

## Diagnosis

{{DIAGNOSIS}}

## Candidate moves

Notation: the state is (Theta, P), with Theta = (lambda, mu, psi, theta, s, c) and P the
node bins and branch paths. The target is
pi(Theta, P) = prior(Theta) x A(P; Theta) x L(D | l(P, c, s)), where A is the exact
augmented density (`.LogAugmented()`) and l_e = c * int_e e^(gamma s z(t)) dt. The tempo
exponents are alpha (generator), beta (lambda, mu), delta (psi) and gamma (characters).

### 1. Shift move (recommended)

- **Proposal.** Draw m in {-1, +1} with probability 1/2 each. Replace every bin k (all
  node bins, including sampled-ancestor and tip bins, and every piece of every path) by
  k + m, and keep all jump times. Then x moves by d = m h s everywhere, with
  h = (z_max - z_min) / (K - 1). Multiply each compensated free parameter by
  exp(-e_p d), where e_p is its tempo exponent: c by exp(-gamma d), lambda and mu by
  exp(-beta d), psi by exp(-delta d) (unchanged at the default delta = 0), and theta by
  exp(-alpha d).
- **Grid edge.** If any bin would leave 1..K, reject.
- **Acceptance.**
  log r = log pi(Theta', P') - log pi(Theta, P) + sum_{p compensated} log(p' / p).
- **Invariance.** T_m is a deterministic bijection of the extended state with
  T_m^-1 = T_-m, and m is chosen symmetrically. In log coordinates, T_m translates the
  compensated parameters by -e_p d, where d depends only on s, which T_m leaves unchanged,
  so its Jacobian there is 1. In natural coordinates that is prod p' / p (Green 1995).
  A rejected edge proposal is a self-transition, and its reverse would also leave the grid.
- **Character likelihood.** It is unchanged exactly when c is compensated:
  l_e' = c e^(-gamma d) int e^(gamma s (z + m h)) = l_e. No pruning is needed; the
  prototype recomputes it only as a check.
- **Cost.** One `.Precompute()` (E ODE) and one `.LogAugmented()`; on acceptance, slices
  and messages. That is the same as one parameter-sweep proposal.
- **Why it helps.** The characters fix c x clock and the tree fixes the per-bin
  speciation rates. The overall tempo level is identified only by the OU pull towards
  z = 0 and the root prior. The move travels along this ridge without disturbing either.
- **Oracles (rung 4b tree, one character).** Each compensation has a one-parameter oracle
  against the exact posterior (`rateShift`, `thetaShift`, `lambdaShift`, `muShift`), and
  the four-parameter version has one against `MarginalMCMC(character = )` (`multi`).
  Results are in the oracle table below. {{SHIFT_ORACLE_VERDICT}}

### 2. Move schedule (recommended)

- With characters, the `tree` move and the path part of `joint` are independence
  proposals from a law that ignores the characters. Their acceptance falls roughly as
  exp(-var(charLL) / 2) (probe table). They should not take most of the CPU.
- Proposal for the engine: during burn-in, record the acceptance of `tree` and `joint`.
  At the end of burn-in, fix each one's frequency, e.g. drop it when its acceptance is
  below 1%. Changing move frequencies only during burn-in keeps the sampling phase a
  fixed, valid mixture.
- Without characters, or with gamma = 0, `joint` remains the move that makes theta mix
  (rung 4b, #20). Keep it there.
- Local moves: run about 0.25 x (number of path branches) `branch` moves per iteration,
  0.25 x (internal nodes) `node` moves and 0.25 x (tips) `tip` moves. In this prototype
  each local move costs about 2-3 ms on the large tree, almost all of it in
  `MkLikelihood()`. An incremental (dirty-node) pruning in the engine would cut this
  cost, and the case for more local moves grows with it.

### 3. Anchored joint move (optional, for theta and s)

- **Anchors.** On branch v, the anchors are its end bins plus the path's bin at each
  interior slice bound b_i with i mod k = phase. The phase is drawn uniformly from
  0..k-1 on every proposal, independently of the state. With k = Inf, the anchors are the
  node bins only.
- **Proposal.** Theta'_p = Theta_p f with log f ~ U(-step/2, step/2). Then each
  anchor-to-anchor segment is redrawn from the frozen endpoint-conditioned law at Theta':
  `.SampleBranch()` on the sub-pieces, with both ends fixed.
- **Proposal density.**
  q(P' | anchors, Theta') = prod_v frozen_v(P'_v) / prod_seg T_seg[a, b], where frozen_v
  is `.LogPathExact() - .PathWeight()` and T_seg is the product of the frozen piece
  transfers between anchors a and b (`.LogTransferRow()` on the sub-pieces).
- **Acceptance.**
  log r = log pi(Theta', P') - log q(P' | a, Theta') - log pi(Theta, P) + log q(P | a, Theta) + log f.
- **Invariance.** The anchors are a function of the state that the proposal leaves fixed,
  so q(P | a(P'), Theta) and q(P' | a(P), Theta') are densities on the same fibre.
  The ratio is then the ordinary MH ratio for a joint (parameter, latent) proposal. With
  k = Inf this is `joint` with the node bins held, and `.LogBridges()` gives the same
  density, checked to 1e-14.
- **Compensated s.** `AnchoredJoint(comp = )`, configuration `L25SAc`, also multiplies
  c, lambda, mu and theta by A(s) / A(s'), where A is the trapezoid exposure
  sum int e^(e_p s z) dt over the anchors. A is a function of the anchors only, so the map
  is an involution in log coordinates, and its Jacobian there is 1.
- **Cost.** In an engine, about one `joint` proposal. The prototype loops over segments
  in R and is 2-3x slower. The projected-cost table prices it as one `joint` proposal.

### 4. NCMC parameter move (optional, for theta and s)

This is the non-equilibrium candidate MC move of Nilmeier et al. (2011).

- **Protocol.** Theta_p moves to Theta_p f in J geometric stages. After each stage
  j < J, the paths relax under n random-scan local moves (`node`, `tip`, `branch` with
  probabilities 0.3, 0.2, 0.5) at Theta_j.
- **Work.** W = sum_j [log pi(Theta_j, P_{j-1}) - log pi(Theta_{j-1}, P_{j-1})].
- **Acceptance.** min(1, exp(W) f).
- **Invariance.** Each local move is an MH kernel reversible for pi(. | Theta_j), and so
  are their random mixture and its powers. The reverse protocol (factor 1/f) visits the
  stages in reverse order, which gives detailed balance on the space of the intermediate
  paths.
- **Cost.** J x (`.Precompute()` + slices + messages + n local moves). Prototype:
  J = 8, n = 0.25 x (branches).

## Validation of the prototypes

Every oracle runs on the rung 4b tree with one character. The new move is the only
update of the parameter, except in the shift oracles, where a shift alone is not
irreducible and the random-walk sweep is kept. Results are compared with the exact 1-D
posterior on a 161-point `TreeLikelihood(character = )` grid, or with
`MarginalMCMC(character = )` for several free parameters. z is (MCMC - exact) / MCSE of
the quantile, using batch means as in `data-raw/rung4b-one-parameter.R`.

{{ORACLES}}

## Gains

{{GAINS}}

## Deviations and limits

{{DEVIATIONS}}

## Reproduce

```sh
Rscript dev/mixing/run.R sim                 # trees.rds
Rscript dev/mixing/run.R check               # Driver() == TempoMCMC(); density identities
Rscript dev/mixing/run.R probe <tree> <nChar>
Rscript dev/mixing/run.R chain <tree> <nChar> <config> <rep> [nIter] [burnin]
Rscript dev/mixing/run.R oracle <variant> [nIter] [nChain] [block]
Rscript dev/mixing/summary.R dev/mixing/results && python dev/mixing/report.py dev/mixing/results
```

The task lists of the Hamilton runs are in `results/tasks*.txt`: one line per Slurm array
task, run as `Rscript dev/mixing/run.R <line>`, with
`EVOTEMPO_MIXING_OUT` pointing at the results folder.
