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
- Characters and time-tree moves: `MkLikelihood()` (pure-R Mk pruning, a
  stand-in for the MkPrime kernels); `TempoMCMC()` gains characters on
  tempo-weighted branch lengths and moves on node ages, fossil ages within
  stratigraphic ranges (`TempoTree(ranges = )`), subtree position and
  sampled ancestor <-> tip; `TreeLikelihood(character = )` gives the exact
  product-space likelihood of one character, and `MarginalMCMC()` accepts it
  ([#5](https://github.com/agent-issues/EvoTempo/issues/5)).
- C++ engine (`engine = "cpp"`, the default; `TempoEngine()`): the marginal
  tree likelihood, extinction probabilities and the E grid by a stiff Radau
  IIA solver; frozen messages, FFBS path draws, bridges, the augmented
  density and exact/frozen weights; and `MkLikelihood()` by the MkPrime
  header kernels (`LinkingTo: MkPrime`).  `engine = "R"` keeps the pure-R
  reference; seeded samplers draw identical paths in either engine
  ([#16](https://github.com/agent-issues/EvoTempo/issues/16)).
- `TempoMCMC()` gains `joint` moves, which propose a parameter together with a
  fresh FFBS draw of every path; without them `theta` mixes too slowly for
  reliable posterior tails, especially with characters
  ([#20](https://github.com/agent-issues/EvoTempo/issues/20)).
- Faster `TempoMCMC()` with the C++ engine, about a third less time per
  iteration on the rung-6 SBC workload: the E grid is interpolated from the
  stiff solver's continuous output instead of stopping at every grid age; slice
  eigendecompositions use LAPACK's tridiagonal solver; `rate` proposals reuse
  the frozen model; bridges skip two redundant root-finder evaluations
  ([#26](https://github.com/agent-issues/EvoTempo/issues/26)).
- Faster tempo-path draws with the C++ engine: bridge jump times come from a
  safeguarded Newton-Halley iteration instead of Brent's method (about 2.5
  instead of 7.9 evaluations per jump), so the FFBS tree draw takes ~40 % less
  time ([#26](https://github.com/agent-issues/EvoTempo/issues/26)).
- `TempoMCMC()` tree moves run in C++: `nodeAge`, `slide` and `saJump` are
  two to two-and-a-half times faster
  ([#26](https://github.com/agent-issues/EvoTempo/issues/26)).
