# EvoTempo: brief

## Concept

Budd & Mann (2025, *Syst. Biol.* 74:917-934) propose a **covariant evolutionary tempo** (CET): each lineage
carries a heritable log-tempo *x* that scales speciation (lambda e^x), extinction (mu e^x) and molecular or
morphological change (c e^x). Log-tempo follows a modified Ornstein-Uhlenbeck process,
dx = -theta x e^x dt + sqrt(2 theta s^2 e^x) dW. The paper gives only forward PDEs (clade-size distributions,
mean diversity, tempo distributions), with no tree likelihood and no inference. **EvoTempo** builds that
inference.

## Approach

1. **Discretized tempo (dCET).** Log-tempo is a K-state CTMC on a grid equally spaced in z = x/s, with
   generator Q = diag(e^(alpha x)) A, where A discretizes the OU operator (central differences by default;
   upwind and square-root approximation available). dCET is a state-dependent fossilized-birth-death-SSE model with
   a hidden tempo state. It has exact Gillespie simulation, exact path densities, K-dimensional
   extinction-probability ODEs, and converges to CET as K grows.
2. **Data augmentation.** Tempo paths on the branches of the reconstructed tree are sampled; hidden lineages
   are marginalized through E_k(t); characters are scored by pruning on branch lengths
   l_e = c int e^(gamma x) dt.
3. **A 2 x 2 nested model family** over tempo->diversification (beta) and tempo->character rate (gamma):
   (0,0) constant FBD + strict clock; (0,1) autocorrelated clock; (1,0) tempo-dependent FBD-SSE; (1,1) CET.
4. **Tip ages are essential.** Branch change is Exp(lambda + mu) whatever the tempo (B&M eq. 12), so the CET
   signal lies in the mismatch between character change and real time; empirical tests need tip-dated data.

## Defaults (decided 2026-10-12)

alpha = 1 in every model cell; start-bin prior = equilibrium (Perron) vector of Q^T + R; crown/extant survival
conditioning; constant fossil sampling psi (psi e^x optional); central-difference generator with h <= 2/z_max;
inference z-range +-5; OU mean fixed at 0; removal-after-sampling fixed at 0; rho fixed externally.

## Architecture

Standalone R package. A pure-R reference implementation comes first (validation rungs 1-5), with C++ hot paths later.
Character likelihoods will come from MkPrime kernels via `LinkingTo` once MkPrime exports its headers to
`inst/include`. RevBayes is the rung-3 oracle (`dnCDBDP`) and hosts the UCLN+FBD comparator.

## Expected outputs

- dCET mathematical specification (generator, E-ODE, augmented density, identifiability).
- Forward simulators (dCET exact; continuous CET exact in tempo-time) and a reproduction of Budd & Mann's figures.
- R reference inference passing validation rungs 1-6; C++ engine with SBC and power study.
- Empirical comparison of the 2 x 2 (beta, gamma) family plus a UCLN comparator on tip-dated datasets.

Target journal: *Systematic Biology* (placeholder; to be confirmed).
