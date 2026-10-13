# EvoTempo master plan

Five phases. Each validation rung gates the next and has a GitHub issue (label `validation`) that is closed
by the PR that meets it. This file is the plan; issue state is the status.

## Phase 0: Groundwork (done before this repository existed)

Scoping memo, dCET mathematical specification, precedents review, Budd & Mann numerics notes.

## Phase 1: Simulator and numerics

1. Repository scaffold, AGENTS.md, BRIEF.md, CI.
2. dCET numerics in R: generator (central / upwind / SqRA), stiff E-ODE with rho and psi, mean diversity,
   second moment, survival conditioning, tempo distribution (replicator-mutation), equilibrium Perron vector,
   branch-duration moments, clade-size distribution; K-convergence study. **Rung 1** identities as tests.
3. Reproduction of Budd & Mann (Figs 1-5, eq. A.26) with (a) their exact discretization and (b) converged dCET;
   table of paper vs paper-scheme vs converged; discretization error vs K and x-range.
4. Forward simulators: exact Gillespie dCET; continuous CET exact in tempo-time; rho, psi, sampled ancestors;
   characters on tempo-weighted branch lengths; reconstructed trees. **Rung 2**: Monte Carlo vs numerics.

## Phase 2: Reference inference in R

1. Marginal dCET tree likelihood (D-equation pruning, FBD tips, sampled ancestors, conditioning).
   **Rung 3**: agreement with RevBayes `dnCDBDP`; beta = 0 collapse to analytic FBD.
2. Data-augmented MCMC: FFBS of tempo paths, endpoint-conditioned bridges, node-bin Gibbs, parameter moves,
   characters via pruning on l_e. **Rung 4**: augmentation vs marginal; product-space oracle with one character.
3. **Rung 5**: prior-only MCMC vs forward simulation.

## Phase 3: C++ engine, SBC and power (Hamilton)

1. MkPrime exports likelihood kernel headers to `inst/include`; EvoTempo `LinkingTo` MkPrime.
2. C++ hot paths: E/D solvers, path sampling, pruning.
3. **Rung 6**: simulation-based calibration per (beta, gamma) cell, 20-40-tip trees.
4. **Rung 7**: power study (characters x fossil density x tree size), including (1,1) vs (0,1).

## Phase 4: Empirical tests

Tip-dated datasets (penguins pilot; Hymenoptera; trilobites; placentals second wave). Fit the 2 x 2 family
plus a UCLN+FBD comparator in RevBayes; stepping-stone marginal likelihoods; posterior-predictive checks.
