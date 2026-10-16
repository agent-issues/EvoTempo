# EvoTempo 0.0.0.9000

- Package scaffold.
- dCET numerics: `TempoGrid()`, `TempoGenerator()` (central, upwind, square-root
  approximation), `BuddMannGenerator()`, `BinRates()`, `ExtinctionProbability()`,
  `FBDExtinction()`, `DiversityMoments()`, `SurvivorDiversity()`,
  `TempoDistribution()`, `TempoEquilibrium()`, `AncestralTempo()`,
  `BranchExpectation()`, `BranchDurationDensity()` and `CladeSizeDistribution()`
  ([#1](https://github.com/agent-issues/EvoTempo/issues/1)).
- Forward simulators `SimulateDCET()` (exact) and `SimulateCET()` (exact in
  tempo-time), with fossil sampling, removal and extant sampling;
  `ReconstructedTree()`, `BranchClock()` and `SimulateCharacters()`
  ([#2](https://github.com/agent-issues/EvoTempo/issues/2)).
- Marginal dCET tree likelihood: `TempoTree()`, `DCETModel()` and
  `TreeLikelihood()` (extant sampling, fossil tips, sampled ancestors, fossil
  sampling, removal, root or origin start, crown/survival conditioning), and the
  closed-form `FBDTreeLikelihood()`
  ([#3](https://github.com/agent-issues/EvoTempo/issues/3)).
- Data-augmented inference on fixed time trees: `TempoMCMC()` (FFBS on a
  frozen-extinction proposal, endpoint-conditioned bridges by direct sampling,
  node, tip and branch path moves, exact Metropolis-Hastings correction,
  parameter random walks), `MarginalMCMC()` on the exact marginal likelihood,
  `SampleTempoPaths()` (importance-weighted stochastic mapping) and
  `TempoPriors()` ([#4](https://github.com/agent-issues/EvoTempo/issues/4)).
