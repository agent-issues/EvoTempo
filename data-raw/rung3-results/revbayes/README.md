# Rung 3: TreeLikelihood() vs RevBayes dnCDBDP (issue #15)

RevBayes development branch (revbayes-ms checkout, commit 361c92a, built 2026-10-16), `rb.exe` called by `data-raw/rung3-revbayes.R`. The summary tables (`table.csv`, `convergence.csv`, `sa-diagnosis.csv`) come from `data-raw/rung3-revbayes-summary.R`. `engines.csv` and `engines-summary.csv` come from `data-raw/rung3-revbayes-engines.R`.

## Design
- **Trees.** Four trees, with node ages rounded to multiples of 2^-10 so that tip ages are exact:
  - `extant`: rung-3 ultrametric tree, 12 tips;
  - `fossil`: rung-3 non-ultrametric tree, 10 tips, of which 9 are fossils;
  - `sa`: rung-3a collapse tree in crown form, 5 tips, 1 sampled ancestor;
  - `sa2`: the fossil tree with two fossil tips made sampled ancestors.
- **Parameters.** The 20 castor points of `rung3-validation.R` (lambda, mu, theta, s), with psi = 0.15 in every bin, rho = 0.7 and removal 0. dnCDBDP has no removal parameter: sampled lineages continue.
- **Grids.**
  - K = 16, z in [-3, 3];
  - K = 26, z in [-5, 5] (the default);
  - K = 16 with beta = 0, for the sampled-ancestor diagnosis.
- **dnCDBDP inputs.** Q (`fnFreeK`, `rescaled = false`), per-bin lambda, mu and psi, and `pi` set to EvoTempo's equilibrium root prior. All are taken from EvoTempo's `.Setup()`, so both programs see the same model.
- **Conditions.** `condition = "time"` and `"survival"`.
- **Engine.** TreeLikelihood() used the R engine (deSolve lsode, rtol 1e-10). The C++ engine could not be built on the machine used (no toolchain available). `rung3-revbayes-engines.R` evaluates both engines wherever the C++ engine is compiled.

## Reconciliation of conventions (StateDependentSpeciationExtinctionProcess.cpp)
1. **Root.** dnCDBDP multiplies lambda_i into the root and weights root states by pi: L = sum_i pi_i lambda_i D_L,i D_R,i. EvoTempo omits lambda at the root, so it is given root weights w = pi lambda / sum(pi lambda), and log sum(pi lambda) is added under "time".
2. **Survival.** dnCDBDP divides by sum_i pi_i lambda_i (1 - E_i)^2, where E includes fossil sampling ("any sampled descendant"). This equals EvoTempo `condition = "survival"`, `joint = FALSE`, with the same weights. EvoTempo's crown "extant" conditioning has no dnCDBDP counterpart and is not compared.
3. **Tree shape.** dnCDBDP adds `lnProbTreeShape()` = (nTip - nSA - 1) log 2 - log(nExtant!). Here log(n!) comes from `RbMath::lnFactorial()`, a truncated Stirling series with no exact small-n branch: it is off by -5.06e-4 at n = 1 and -2.1e-5 at n = 2. The comparison reproduces that approximation; the error itself is recorded in `info$stirlingError`.
4. **Printing.** Rev rounds a Real to 6 significant digits in string concatenation (`"x " + lnL`), while `print(lnL)` honours `outputPrecision`. Values are therefore printed on their own line.
5. **Integration.** dnCDBDP integrates with fixed-step Dormand-Prince (dt = rootAge / nTimeSlices; odeint `integrate_adaptive` with an uncontrolled stepper). Stiff points need many slices:
   - at 500 slices many K = 26 cases return NaN;
   - at 1000 one K = 16 case is still off by 0.13;
   - the stiffest K = 26 point (s = 1.49) needs 64000 slices.

## Results
Max |Delta logL| over the 20 points, at 16000 slices (64000 for K = 26 point 5):

| Tree | Grid | time | survival | Verdict |
|---|---|---|---|---|
| extant | K = 16 | 4.7e-9 | 4.6e-9 | Pass |
| extant | K = 26 | 3.8e-9 | 3.8e-9 | Pass |
| fossil | K = 16 | 5.4e-9 | 5.7e-9 | Pass |
| fossil | K = 26 | 7.4e-7 | 7.7e-7 | Pass |
| sa (1 SA) | K = 16, 26 | 1.87 | 1.87 | Fail: dnCDBDP sampled-ancestor term |
| sa2 (2 SAs) | K = 16, 26 | 6.44 | 6.44 | Fail: dnCDBDP sampled-ancestor term |

**Convergence in nTimeSlices** (6 points spanning s, survival; max |Delta|):

| Tree, grid | 1000 | 4000 | 16000 |
|---|---|---|---|
| extant, K = 16 | 2.9e-9 | 3.4e-9 | 3.4e-9 |
| extant, K = 26 | 1.5e-6 | 1.2e-6 | 2.4e-9 |
| fossil, K = 16 | 0.13 | 4.7e-9 | 4.7e-9 |
| fossil, K = 26 | 4.4e-9 | 1.8e-5 | 2.3e-7 |

For K = 26 point 5 with the "time" condition, the fossil tree gives 7.4e-5 at 16000 slices, then 5.1e-9 at both 64000 and 256000.

## Sampled-ancestor diagnosis
dnCDBDP handles a sampled-ancestor tip like any fossil tip, so it contributes psi * E(t) (`computeNodeProbability`). Its parent then omits lambda. The FBD term for a sampled ancestor with removal 0 is psi alone, because the lineage continues through the other child. dnCDBDP therefore multiplies the likelihood by an extra E_k(t) at every sampled ancestor, inside the sum over states.

When the rates do not depend on the bin (beta = 0), E is state-independent. The predicted offset is then sum over sampled ancestors of log E(t_SA), with E from the closed-form FBD (`FBDExtinction()`). After subtracting it, the sampled-ancestor trees agree to 1.6e-8 (`sa`) and 1.3e-8 (`sa2`), and the trees without sampled ancestors to 1.5e-8 (`sa-diagnosis.csv`). With beta = 1 the factor is state-dependent and cannot be removed outside the pruning.

EvoTempo's sampled-ancestor term agrees with the closed-form FBD (rung 3a, 7e-9) and with Monte Carlo simulation (rung 3c). The two failing rows are therefore a dnCDBDP convention, or bug, and not an EvoTempo error.
