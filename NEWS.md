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
