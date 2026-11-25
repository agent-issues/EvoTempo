# Tempo location: the centred parameterisation

*Spec change for issue #116, 2026-11-25. This note is design only; no code follows it yet. Equation and section numbers
without a prefix refer to [dCET-spec.md](dCET-spec.md). The Hymenoptera pilot files cited here are on branch
`hymenoptera-cet-pipeline`, in `data-raw/empirical/hymenoptera/cet/smoke-results/diagnostics/`.*

## L.1 Problem

The generator Q = diag(e^{αx})A, with A the OU block in z = x/s, gives the lineage-time stationary law
N(−αs², s²) (§5.2). The location of the tempo is therefore tied to its spread: when s grows, the centre of the law falls
by αs². The rate parameters (λ, μ, ψ, c, θ) are rates at x = 0, which lies αs standard deviations above that centre.

In Hymenoptera cell (1, 0) (β = 1, γ = 0) only the tree sees the tempo, and it sees it only through λe^{x} along the
paths. The data pin λ at the level where the paths sit; in `ridge-pilot.csv`, λe^{⟨x⟩} stays at 38–70 over nine
segments. The paths stay within about one standard deviation of the stationary centre. So when the likelihood
favours larger s, the centre drags the paths down and λ follows as e^{αβs²}: from 94 to 1651 as s went from 1.1 to 2.2.
Only the 1/λ tail of the prior on λ resisted (`ridge-profile.csv`: likelihood plus prior flat, −842.1, −841.4, −841.5).

On the fixed window z ∈ [−5, 5] the law sits at z = −αs. So the top bins carry rates about e^{αβs²} (birth) and
e^{α²s²} (Q) larger than the identified rates need, and the reflecting floor sits within 3 sd of the mass. One
likelihood evaluation took 0.5, 6.5, 53 and 400 s at s = 0.96, 1.16, 1.40 and 1.69, and then failed.

## L.2 The centred model

**Definition.** Replace (0.1) by

$$dx = -\theta\,(x - \alpha s^2)\,e^{\alpha x}dt + \sqrt{2\theta s^2 e^{\alpha x}}\,dW. \qquad (L.1)$$

This is the §9.1 family with OU mean m = αs², chosen so that the time change tilts the law back to the origin.

**Proposition L.1.** For every α, the lineage-time stationary law of (L.1) is N(0, s²).

*Proof.* The generator is e^{αx}𝓛, where 𝓛 is an OU generator with stationary density N(x; αs², s²). A time-changed
diffusion has stationary density proportional to π_𝓛(x)e^{−αx} (§5.2). Here
exp(−(x − αs²)²/(2s²) − αx) = exp(−x²/(2s²)) · exp(−α²s²/2), which is proportional to N(0, s²). ∎

At α = 0, m = 0 and (L.1) is the current model; nothing changes under the α ≡ γ convention (§7.3) when γ = 0.

**Parameters.** The per-bin rates keep their form: λ_k = λe^{βx_k}, μ_k = μe^{βx_k}, ψ_k = ψe^{δx_k},
ℓ_e = c Σ_k e^{γx_k} d_{e,k}. What changes is the reference point. x = 0 is now the mean and median log-tempo of a lineage
in tempo equilibrium, so:

- λ, μ, ψ (when δ = 1), c (when γ = 1) and θ (when α = 1) are **rates at the median lineage tempo**;
- the lineage-time mean of λe^{βx} is λe^{β²s²/2};
- in the (0, 1) cell, c is the median of a lognormal relaxed clock with log-sd s;
- β, γ, δ and α keep their meaning as couplings;
- with β = 0, λ and μ are unchanged, and cell (0, 0) is unchanged entirely.

**Grid.** The grid keeps its formula: z_k equally spaced on [z_min, z_max] (default [−5, 5]), x_k = s z_k. The window is
now centred on the lineage-stationary law (N(0, 1) in z) for every s. In z the tempo-time OU has mean αs, so the
OU block depends on s through its drift centre:

$$(\mathcal{L}_z f)(z) = -\theta (z - \alpha s) f'(z) + \theta f''(z).$$

Consequence 1 of §1.1 (A independent of s) no longer holds; consequences 2 and 3 still do. In particular, paths stored as
bin indices stay valid when s moves.

**Discretisation.** Positivity of the central scheme now needs h · max_k |z_k − αs| ≤ 2. At K = 26 on [−5, 5]
(h = 0.4) this fails for every s > 0, so **SqRA becomes the default scheme** for the centred model. SqRA applied to
A with potential U(z) = (z − αs)²/2, then multiplied by e^{αx_k}, gives

$$Q_{k,k\pm1} = \frac{\theta}{h^2}\exp\!\Big(\frac{\alpha(x_k + x_{k\pm1})}{2} - \frac{z_{k\pm1}^2 - z_k^2}{4}\Big). \qquad (L.2)$$

The rates are always positive and the rows sum to zero. Detailed balance gives Q_{k,k+1}/Q_{k+1,k} = φ(z_{k+1})/φ(z_k), so
**the stationary law of Q is exactly the grid Gaussian** π_k ∝ φ(z_k), with mean x exactly 0 on a symmetric window.

Expanding (L.2) to first order in h, the jump moments are h(Q_{k,k+1} − Q_{k,k−1}) = θe^{αx_k}(αs − z_k) + O(h²) and
h²(Q_{k,k+1} + Q_{k,k−1}) = 2θe^{αx_k} + O(h²). Both are consistent with (L.1) in z.

The other schemes stay available:
- *central*, with drift −θ(z_k − αs). It stops with a message when the positivity bound fails; the minimum K is
  1 + (z_max − z_min)·max(|z_min − αs|, |z_max − αs|)/2.
- *upwind*, which is always positive and O(h).

**Start priors.**
- `equilibrium` keeps its definition (the Perron vector of Qᵀ + R).
- `stationary` is now the grid Gaussian.
- `uniform` is unchanged.
- A numeric per-bin `rootPrior` keeps x_k = s z_k but sits in a different law, so existing numeric priors must be
  reviewed.

## L.3 Relation to the current model

**Exact map (continuum).** By §9.1 with a = αs², the current model (m = 0) at (θ₀, λ₀, μ₀, ψ₀, c₀, s) and the centred
model at (θ, λ, μ, ψ, c, s) define the same law of trees, samples and characters when

$$\theta = \theta_0 e^{-\alpha^2 s^2},\quad \lambda = \lambda_0 e^{-\alpha\beta s^2},\quad \mu = \mu_0 e^{-\alpha\beta s^2},\quad
\psi = \psi_0 e^{-\alpha\delta s^2},\quad c = c_0 e^{-\alpha\gamma s^2}, \qquad (L.3)$$

with x = x₀ + αs². s, ρ and the couplings are unchanged.

**Exact map (grid).** The centred model on [z_min, z_max] equals the current model on the window
[z_min − αs, z_max − αs] at the parameters (L.3), for every scheme. Both generators share the same jump rates in
x-space and the same per-bin rates. It is not equal to the current model on the unshifted window; the difference is the
truncation. This identity is the main unit test (L.6 V1). It also means the centred model can be evaluated **with the
current code**, by passing `zRange = c(-5, 5) - alpha * s` (L.6 V0).

**Jacobian.** In coordinates (log θ, log λ, log μ, log ψ, log c, s), (L.3) is a shear with unit determinant. In natural
coordinates, |∂(new)/∂(old)| = exp(−αs² · Σ_p e_p) over the free rate parameters p. Here e_p is α for θ, β for λ and μ,
δ for ψ and γ for each clock rate. A run of the current model under prior p₀ is therefore reproduced by the centred
model under

$$p(\theta,\lambda,\mu,\psi,c,s) = p_0\big(\text{(L.3)}^{-1}\big)\cdot\exp\!\big(\alpha s^2 \textstyle\sum_p e_p\big)$$

on the shifted window. This prior couples s to the rates; independent priors on the new parameters define a
different, and intended, joint prior.

**Bayes factors.** At equivalent parameter points the likelihoods agree; the evidence changes only through the prior.
With priors stated on the new parameters, rate priors refer to the same quantity (the rate at the median lineage
tempo) in every (β, γ) cell. Under the current model a prior on λ₀ in a β = 1 cell sits e^{αs²} above the cells'
common reference. Cross-cell Bayes factors under the two parameterisations differ by this prior-mapping term.

## L.4 What the change does and does not do

**Does.**
- λ now tracks what the tree identifies. In the pilot, the rate at the stationary centre (λ₀e^{−s²}) changed by
  ×0.46 over segments 1–9, against ×17.6 for λ₀. The residual is the paths' offset from the centre, O(s) in x, not
  O(s²).
- At fixed identified rates the top-bin rates fall by e^{αβs²} (birth) and e^{α²s²} (Q): ×17 at s = 1.69, ×127 at
  s = 2.2.
- The low-end truncation goes.

**Does not.**
- It does not bound s. The prior on λ no longer penalises growth in s; that brake was an artefact of λ's meaning. So s
  is held only by its own prior and the likelihood, and the s prior must be stated on its merits. s is now the
  log-sd of lineage rates around their median, the role of σ in a lognormal relaxed clock. The trilobite and penguin
  pipelines already match it to the comparators' UCLN σ prior.
- It does not remove the intrinsic stiffness. The top-bin rates still grow as e^{α s z_max}. If V0 shows the solve
  still failing near s ≈ 2 at plausible identified rates, the window (z_max) or K needs revisiting before (1, 0) is
  feasible, and reparameterisation alone will not rescue it.

## L.5 Candidates considered

| | Reference point | Verdict |
|---|---|---|
| **A (adopted)** | Lineage-stationary median, law N(0, s²); OU mean αs²; window centred | Tracks the identified rate best in the pilot (×0.46); cuts the e^{αs²} stiffness excess; closed-form map |
| B | As A, relabelled only (old drift and old window) | Same target as fix C's arm A, which stalled; numerics unchanged. Rejected |
| C | Unit lineage-mean rate, law N(−s²/2, s²); OU mean (α − ½)s² | Its reference point still moves with s (λ₀e^{−s²/2} changed ×3.0 in the pilot). Second choice |
| D | Species equilibrium p\* | Map depends on (λ, μ, θ); no closed-form Jacobian; paths do not sit at p\*. Rejected |
| E | Current model plus priors (fix C) | Superseded by the maintainer's decision of 25 Nov |

## L.6 Implementation and validation plan

### Interface

- One argument, `location = c("centred", "bm")`, on `DCETModel()`, `TempoGenerator()`, `SimulateCET()` and any other
  constructor of Q or of the continuous process. The default is `"centred"`; `"bm"` is the current model (m = 0).
- `DCETModel()`'s default scheme becomes `"sqra"` for `"centred"` and stays `"central"` for `"bm"`.
- `TempoGrid()` is unchanged.

### R engine (the reference)

- `.OUBlock()` takes the drift centre ζ = αs (0 under `"bm"`): U = (z − ζ)²/2 for SqRA, drift −θ(z − ζ) for central
  and upwind, and the centre-aware positivity check.
- `.Setup()` passes ζ.
- `SimulateCET()`'s exact OU in tempo time (`.OUSegment()`) gets mean αs². The Gillespie dCET simulator uses Q and needs
  no change.
- `TempoEquilibrium()`, `.RootPrior()`, `BinRates()`, the shift move (`.Shifted()`, unchanged map and Jacobian) and
  the tree moves take Q and rates from `.Setup()` and need no change.
- Any code that rebuilds Q outside `.Setup()` must pass the centre. Audit this with
  `grep -n "TempoGenerator(\|\.OUBlock(" R/`.
- Documentation: `DCETModel()`, `TempoPriors()` (what λ means), the AGENTS.md model summary, and the empirical
  READMEs' "comparator-matched at unit tempo" notes, which become "at the median lineage tempo".

### C++ engine

- The C++ engine receives Q and the per-bin rates from R. Confirm that no kernel rebuilds them: `src/tree_moves.cpp`
  `.Rebuild()` and the slice set-up.
- C++ must agree with R to ≤ 1e−8 on the rung 3/4 fixtures, rerun at the new default.
- Seeded `TempoMCMC()` chains must be identical in both engines (`test-Engine.R`).

### Validation

**V0 (go/no-go, before any code).** Rerun the `diag-ridge.R` profile (the six points of `ridge-profile.csv`, on the
chain's last tree) with `DCETModel(scheme = "sqra", zRange = c(-5, 5) - s)` at each point's parameters. By L.3 this is
the centred model at the mapped parameters. Report log-likelihood and seconds per point beside the current grid. If the
points at s ≥ 1.69 still fail or exceed about 60 s, revise the window or K in this spec before implementing.

**V1 (unit tests, both engines).**
1. The L.3 grid identity, centred against `"bm"` on the shifted window, for TreeLikelihood at 1e−10 relative. Cover
   α, β, γ ∈ {0, 1}, δ ∈ {0, 1}, `start` root and origin, and every scheme.
2. The SqRA stationary law equals φ(z_k) normalised to 1e−12, with mean x = 0.
3. The β = 0 collapse to the analytic FBD (§2) still holds.
4. `"bm"` reproduces current results exactly: seeded chains are bit-identical to those before the change.
5. Coverage: 100% of the new lines.

**V2 (rungs; rerun at the default grid, `EVOTEMPO_GRID=default`).**
- Rung 1 is unchanged under `"bm"`.
- Rung 2: the continuous simulator under `"centred"` against the numerics (lineage-stationary mean 0, survival, N_x),
  plus the negative control.
- Rung 3: castor, Monte Carlo and RevBayes `dnCDBDP` with the new Q passed as a custom matrix.
- Rung 4: 4a–4d, `rung4b-likelihood-check.R`, `rung4b-theta-identity.R`, `rung4b-one-parameter.R` (every free
  parameter, and the `shift` arm with 16 or more chains) and `rung4b-shift-fixed.R`.
- Rung 5: `rung5-validation.R`, `rung5-fixed-parameters.R`, `rung5-age-oracle.R` and `rung5-age-validation.R`.
- Update their result files.

**V3 (SBC, rung 6).** Rerun at full size on Hamilton for every (β, γ) cell, with priors on the new parameters.

**V4 (power, rung 7).** Optional. Rerun only if model-selection conclusions are to be cited under the new model.

**V5 (Hymenoptera (1, 0) re-pilot).**
- Design: the fix-arms design (chain 3, draw 3, same seeds, 200 burn-in + 3000 iterations), then 2 chains × 20 000 as
  `fixC-long`.
- Switches: `CET_FBD_SCALE` and `CET_S_MAX` off.
- Priors: the comparator's FBD-block prior on the new (λ, μ, ψ); for s, the comparator's UCLN σ prior, as in the
  trilobite pipeline.
- Report: λ, s, λe^{⟨x⟩}, seconds per segment and R-hat.
- Pass when: window medians of λ and s show no monotone drift; segment cost stays within twice that of cell (0, 0);
  and s does not pile against its support edge.
- Recheck that cell (1, 1) is qualitatively unchanged.

The penguin and trilobite pipelines are not part of this plan until the maintainer answers the question in #116.
