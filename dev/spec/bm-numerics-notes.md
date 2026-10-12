<!-- Copied from the project artifact store on 2026-10-12. {{artifact:...}} links refer to that store and do not resolve on GitHub. -->

# Reference numerics for validation rung 1: Budd & Mann's PDE solutions

*2026-10-12, version 3. Written from the actual Dryad code. Companion to [dCET-spec.md]({{artifact:dde49858-dfd2-4474-9b36-c7ef3ceb7be0}}).*

**Source.** [BuddMannTempoSI.R]({{artifact:c05255f0-a1cc-4b52-a523-f30188cc1779}}) (Dryad doi:10.5061/dryad.q573n5ts2). The file was supplied by the user and its
sha-256 matches Dryad's metadata (`229618841af84dd29c510acc9aef79bc6e9cc55efd509d41cc3e2129a478f977`, 17,629 bytes).

**Re-runs.**
- [bm_reproduce.R]({{artifact:65cc1388-dd0f-4aa9-9097-d6d0f3d42abe}}) runs B&M's own `pdeG`, `pdeM`, `pdeMM2`, `pdeC` and `evalM` **verbatim** with deSolve 1.42 (conda) and
  R 4.5.3. ReacTran 1.4.3.2's R functions `tran.1D`, `setup.grid.1D` and `setup.prop.1D` are sourced from the CRAN source tarball,
  because the compiled package could not be loaded in this sandbox. These functions are pure R, so the numerics are unchanged.
  Output: [bm_reproduce_output.txt]({{artifact:27000047-2dbe-426e-b8b3-9eccbaee2542}}).
- [bm_scheme_replica.py]({{artifact:5d0a72a5-a98b-4ca4-b955-69183cc78aab}}) is an independent Python replica of the same discrete operator, plus grid refinement.
  Output: [bm_scheme_replica_output.txt]({{artifact:cb32d2e0-313a-42a4-875b-5f344416dd5d}}). It reproduces the R values (mean 57 836 in both).
- **Not re-run:** `evalQ`, the PGF/FFT inversion (MM = 2²¹ complex ODE solves; "days" per the code comment). Hence not reproduced:
  Fig 1A's distribution, P(n < 1000), P(n > 50 000), Fig 1B and Fig 4B.

## 1. The scheme, block by block

### 1.1 Grid and spatial discretisation

| Item | What the code does |
|---|---|
| x-grid | `setup.grid.1D(x.up=-10, x.down=10, N=lx)` with `lx = 101`. This is a **finite-volume** grid of 101 cells, width h = 20/101 = 0.19802, centres x_j = −10 + (j − ½)h ∈ [−9.901, 9.901]. Cell 51 is centred at x = 0 (5×10⁻¹⁵). `x0idx = 51` |
| e^x factor | `setup.prop.1D(func=exp)`; `ex$mid` = e^{x_j} at cell centres, multiplying the whole right-hand side (diffusion, drift and reaction) pointwise. Backward (non-conservative) form throughout |
| Diffusion θs²∂²_x | `tran.1D(C, D=1, v=0, dx=Grid)` × θs². Standard 3-point stencil (C_{j+1} − 2C_j + C_{j−1})/h² |
| Drift −θx∂_x | `tran.1D(C, D=0, v=1, dx=Grid)` × θ·x_j. The default `AFDW = 1` (full upstream weighting) with v = +1 sets interface concentrations to the upstream (left) cell, giving tran2$dC_j = −(C_j − C_{j−1})/h. **The direction of differencing is fixed by v = +1, not by the sign of the drift.** |
| Boundaries | `tran.1D` defaults `C.up = C[1]`, `C.down = C[N]`: zero gradient, so the diffusive flux is 0 at both ends (B&M's "no flow"). The advective term vanishes in cell 1 (interface value C.up = C₁). Cell N keeps −θx_N(C_N − C_{N−1})/h |

**How the fixed v = 1 interacts with the drift sign.** The backward drift coefficient is b(x) = −θx. A monotone (upwind) scheme
differences towards the side the drift comes from: backward for b < 0, i.e. x > 0, and forward for b > 0, i.e. x < 0. The code
always uses the backward difference, so it is **upwind for x > 0 and downwind for x < 0**. As a CTMC the semi-discrete operator has
- up-rate u_j = θs²/h²;
- down-rate d_j = θs²/h² + θx_j/h.

Consequences:
1. **Drift is exact.** h(u_j − d_j) = −θx_j.
2. **Diffusion is asymmetric.** h²(u_j + d_j) = 2θ(s² + x_jh/2). The effective diffusion is θ(s² + 0.099x), i.e. +20% at x = 2 and
   −20% at x = −2. The symmetric upwind proxy used in v1 of these notes adds θ|x|h/2 instead.
3. **Negative rates.** d_j < 0 for x_j < −s²/h = −5.05. In 24 cells (x ≤ −5.149) the operator is not a valid generator. Effects are
   small because e^x < e^{−5} there, but positivity is not guaranteed.
4. **First-order convergence** to the same limit as the dCET central scheme (§3).

### 1.2 Time stepping

| Item | What the code does |
|---|---|
| Real-valued PDEs (G at z = 0, mean, second moment, M arrays, C/H, Fig 5) | `deSolve::vode` with `atol = 1e-15` and **default `rtol = 1e-6`**. vode's default method is stiff BDF with an internally generated full Jacobian. Output on the integer grid t = 0, 1, …, 500 (Fig 5: 0–20 by 0.01) |
| Complex PGF (`evalQ`) | `deSolve::zvode` on U = G − μ/λ (`pdeU`, reaction λU(U − 1 + μ/λ)), `atol = 1e-18`, default rtol, single output at t = 500. The shift makes U → 0 as G → μ/λ |
| Tolerance sensitivity | Re-running Fig 1A with `rtol = 1e-10` changes neither meanclade (57 835.9) nor meanexclade (425 699.0). Time-stepping error is negligible; all discrepancies below are spatial |

### 1.3 PGF inversion (`evalQ`, not re-run)

- Evaluation points: MM = 2²¹ = 2 097 152, z_m = exp(2πi m/MM) on the **unit** circle (no damping radius).
- Computation: chunks of 2¹⁰ points via `parallel::mcmapply` on 6 cores (`mc.cores > 1` is not supported on Windows).
- Inversion: `Q[i,] = abs(fft(z[i,]))/MM`. R's `fft` is the unnormalised forward DFT, so Q[i, n+1] = P_n for start cell i, aliased
  modulo MM. `abs()` hides small negative or complex numerical errors.
- **Caveats (not quantified):**
  1. Mass at n ≥ 2²¹ wraps onto small n. The conditioned experienced size is ≈ 4×10⁵, so the tail beyond 2×10⁶ is not negligible
     a priori.
  2. zvode's default rtol of 10⁻⁶ limits the absolute accuracy of the far-tail P_n.
- The code comment says MM "needs to be at least 2^20".

### 1.4 What each figure plots

| Figure | Quantity in code | Notes |
|---|---|---|
| 1A curve | `Q[x0idx, 2:(1+2*meanexclade)]`, the **unconditioned** P_n, n ≥ 1, start x = 0, T = 500 | Not renormalised by S. Sums to S₀(500) ≈ 0.021 over n ≥ 1 |
| 1A blue/red lines | `meanclade` = N/S = E(n \| n > 0); `meanexclade` = E[n²]/E[n], from `pdeMM2` (second raw moment) | Both from the 101-cell scheme |
| 1A dashed | `dgeom(n, 1/meanclade)`, a *normalised* geometric pmf | **Normalisation mismatch** with the unconditioned curve (factor 1/S ≈ 47, ≈ 1.7 decades on the log axis). From code inspection; the figure was not re-rendered (needs Q) |
| 1B | `p = Q[x0idx, 2:2^21]` renormalised (survivors); P(share of one sister < w) for two independent sisters from x = 0, w in steps of 0.05 | Uses the normalised p (correct) |
| 2A | `meanNfromx[, idx]` for x₀ ≈ −2, 0, 2 (eq 4), including extinct clades | |
| 2B solid | `(meanNfromx − meanNfromx_correction)/(1 − D[501,])`, A.13 with C from `pdeC` | See bug §2.1 |
| 2B dashed | `rowSums(M[,,i]·S_y(500 − t))/S_i(500)` (A.16–A.17) | Correct |
| 3A | `evalmx(M)`: mean ± sd of y under M[t, ·, x₀]/ΣM, r = θ, x₀ ∈ {0, −2, 2} | M[t, y, x] = N^y_x(t) obtained by solving the *backward* equation from a one-hot initial condition for each y (101 solves); equivalent to A.21 |
| 3B | the same from x₀ = 0 for r = 0 (`Mlowr`), θ (`M`) and 2θ (`Mhighr`) | Finite-t curves; t = 500 values are not the equilibria (§3) |
| 4A | Eq 7: Px · `meanNfromx[i,]` normalised, mean ± sd, reversed for plotting (label "100 Ma" = look-back 100 Myr) | `Px = M[501,,lx]` normalised: the r = θ distribution at t = 500 from the **top** cell (x = 9.90), used as an equilibrium proxy (mean 0.0005, sd 0.831) |
| 4B | Eq 8: Σ_{x₀} p(x_t \| x₀)·P(n_T = `cladesize` \| x₀)·Px(x₀)/P(n_T), with `cladesize = round(margexmean)` | Needs Q. `margmean = Σ Px·N_x/S_x` (joint conditioning per x₀) |
| 5A | `pdeM` with r → −(λ+μ) from initial condition e^x(λ+μ) (**positive** IC, confirming the eq 9 typo), x₀ ∈ {−2, 0, 2}, t ∈ [0, 20] | Each curve rescaled to peak 1.3 (shapes only) |
| 5B | `dexp(ts, λ+μ)`, **analytic**, not a numerical solution of eq 11 | Rescaled likewise |

### 1.5 Parameters

λ = 0.51, μ = 0.5 (r = 0.01), θ = 0.01, s = 1, T = 500 Myr, as in the paper. `Mhighr` uses r = 0.02 and `Mlowr` r = 0.

## 2. Bugs and fragilities found in the code

1. **`pdeC` reaction term.** The second state C2 is H − μ/λ: the initial condition is `D − mu/lambda`, and the C-equation
   coefficient ((2C2 − 1)λ + μ) equals λ(2H − 1) − μ only under that shift. But C2 is advanced with λC2(C2 − 1) instead of
   λC2(C2 − 1 + μ/λ) (cf. `pdeU`, which has it). Effect on the Fig 2B survivor-conditioned mean from x₀ = 0, comparing code with
   corrected code (verbatim re-run):
   - t = 100: 157.79 vs 161.04 (−2.0%);
   - t = 250: 1530.46 vs 1533.92 (−0.2%);
   - t = 400: 13548.45 vs 13551.00 (−0.02%).

   Not visible on the log-scale figure.
2. **Asymmetric fixed-direction differencing of the drift** (§1.1). It inflates high-tempo diffusion and creates negative rates for
   x < −5.05. This is the main source of the Fig 1A offsets (§3).
3. **`evalmx` misplaced parenthesis.** `which.min(abs(Grid$x.mid)-x0[1])` should be `which.min(abs(Grid$x.mid-x0[1]))`. It is
   **harmless as called**: x0[1] = 0 in all three calls, so both give cell 51. The other two lookups are written correctly.
4. **Hard-coded `102:202`** in `margexmean` selects the second-moment columns. It is correct only for lx = 101 (= lx+1 : 2lx).
5. **Fig 1A normalisation mismatch** between the unconditioned P_n curve and the normalised geometric reference (§1.4).
6. **Fig 3B / A.26.** At r = 0 the plotted curve from x₀ = 0 reaches only −0.78 by 500 Myr (§3). The paper's "−1" is the asymptote,
   which low-tempo dynamics approach slowly.
7. `calcQ = F` is the default, so the script errors on a fresh session (Q undefined in the Fig 1 block) unless Q is pre-loaded.
   Output paths are hard-coded to `~/Dropbox/...`.

## 3. Rung-1 reference values: paper text vs B&M code vs converged dCET

Core parameters, x₀ = 0 unless stated.
- **B&M code**: verbatim R re-run (101 cells), or the Python replica where marked †.
- **Converged dCET**: central scheme. Fig 1A uses K = 401 on z ∈ [−10, 10] ([dcet_scratch_output.txt]({{artifact:ff706844-85ce-4af4-bb6e-4b1ceae37fdc}})); other rows use
  K = 161 on [−8, 8], or K = 401 on [−10, 10] where marked ‡.

| Quantity | Paper text | B&M code | Converged dCET | Code vs converged |
|---|---|---|---|---|
| S₀(500) | — | 0.02138 | 0.02154 | −0.7% |
| Fig 1A mean clade size E(n \| n > 0) | c. 60 000 | **57 836** | **52 644** (B&M-scheme Richardson limit: 52 621†) | +9.9% |
| Fig 1A experienced size E[n²]/E[n] | c. 400 000 | **425 699** | **352 938** (Richardson: 352 652†) | +20.6% |
| ratio experienced/mean | "c. 8" | 7.36 | 6.70 | |
| Fig 1A P(n < 1000 \| n > 0); P(n > 50 000 \| n > 0) | ≃ 1/3; ≃ 1/4 | not re-run (needs `evalQ`) | not computed | open |
| Fig 2A N_x(100), x₀ = −2/0/2 | — | 1.19 / 3.75 / 11.76 | 1.19 / 3.70 / 11.47 | |
| Fig 2A N_x(500), x₀ = −2/0/2 | — | 52.38 / 1236.38 / 4373.25 | 53.51 / 1134.77 / 3971.16 | −2% / +9% / +10% |
| Fig 2A asymptotic slope | "fixed rate" | 0.014511 (log-slope 400–500) | 0.014332 (= Perron Λ) | +1.2% |
| Fig 3A, r = θ, x₀ = 0: mean, sd at t = 500 | mean 0 | 0.000, 0.831 | 0.000, 0.833 | |
| Fig 3A mean from x₀ = −2 at t = 10/50/100/200/500 | — | −1.953 / −1.834 / −1.658 / −1.138 / −0.038 | −1.973 / −1.854 / −1.676 / −1.134 / −0.037 | |
| Fig 3A mean from x₀ = +2 at t = 10/50/100/200/500 | — | 1.165 / 0.411 / 0.168 / 0.036 / 0.000 | 1.181 / 0.422 / 0.174 / 0.038 / 0.001 | |
| Fig 3B as plotted, x₀ = 0, t = 500: r = 0 / 2θ | "−1" (A.26 limit) | −0.777† / 0.711† | −0.804‡ / 0.635‡ (equilibria: −1 / 0.634) | |
| Fig 4A mean ancestral x, look-back 100/200/500 Myr | "plateau ≃ 0.6 at c. 100 Ma" | **0.462 / 0.548 / 0.572** | **0.443 / 0.528 / 0.552** | +0.02 |
| Fig 5A mean branch duration, x₀ = −2/0/2 | — | code integrates only to 20 Myr: 5.49 (94% of mass) / 0.9998 / 0.141 | 7.244 / 0.9998 / 0.138 (untruncated) | |
| Fig 5B w-law | Exp(λ+μ) | plotted analytically | exact (mean 0.990099 in every bin) | |

**B&M scheme under refinement** (Python replica†; first-order convergence to the dCET limit):

| N cells | h | mean \| survival | experienced | Fig 4A at 100 / 500 Myr |
|---|---|---|---|---|
| 101 | 0.198 | 57 836 | 425 692 | 0.462 / 0.572 |
| 201 | 0.0995 | 55 129 | 387 156 | 0.452 / 0.561 |
| 401 | 0.0499 | 53 857 | 369 526 | 0.448 / 0.556 |
| 801 | 0.0250 | 53 239 | 361 089 | 0.445 / 0.554 |

**Conclusions for rung 1.**
1. The paper's Fig 1A numbers are the code's 101-cell values rounded: 57 836 → "c. 60 000" and 425 699 → "c. 400 000". They
   overstate the converged PDE values by about 10% (mean) and 21% (experienced) because of the asymmetric first-order drift
   differencing.
2. The Fig 4A "≃ 0.6 at c. 100 Ma" is a loose verbal reading of the code's own curve (0.46 at 100 Ma, 0.57 at 500 Ma). It is not a
   different quantity.
3. **Acceptance criteria.**
   - The dCET implementation, run with B&M's operator (101 cells, fixed backward drift differencing, zero-gradient ends),
     must reproduce the "B&M code" column to solver tolerance.
   - Run with its own converged scheme, it must reproduce the "converged dCET" column.
   - The paper text is a qualitative target only.
