# Tempo location: the mean-rate parameterisation

*Spec change for issue #116, 2026-11-25. This note is design only; no code follows it yet. Equation and section numbers
without a prefix refer to [dCET-spec.md](dCET-spec.md). The Hymenoptera pilot files cited here are on branch
`hymenoptera-cet-pipeline`, in `data-raw/empirical/hymenoptera/cet/smoke-results/diagnostics/`.*

## L.1 Problem

The generator Q = diag(e^{αx})A, with A the OU block in z = x/s, gives the lineage-time stationary law
N(−αs², s²) (§5.2). The location of the tempo is therefore tied to its spread: when s grows, the centre of the law falls
by αs². The rate parameters (λ, μ, ψ, c, θ) are rates at x = 0, which lies αs standard deviations above that centre,
and the rate a lineage actually has on average is λe^{−αβs² + β²s²/2} (λe^{−s²/2} at α = β = 1).

**The Hymenoptera (1, 0) runaway.** In this cell (β = 1, γ = 0) only the tree sees the tempo, and it sees it only
through λe^{x} along the paths. Across the nine segments of `ridge-pilot.csv` the only stable quantities are the rates at
the paths' own level:

| Quantity | log-sd over segments |
|---|---|
| λe^{⟨x⟩} | 0.17 |
| λ⟨e^x⟩ | 0.22 |
| θ⟨e^x⟩ | 0.27 |
| λ (the current parameter) | 0.99 |

The path level itself fell from x = −0.29 to −3.35 as s went from 1.1 to 2.2, following the falling centre, and λ
followed it (94 → 1651). Only the 1/λ tail of the prior on λ resisted. In `ridge-profile.csv` likelihood plus prior stay
flat (−842.1, −841.4, −841.5) while the likelihood alone rises.

**The numerical cost.** On the fixed window z ∈ [−5, 5] the lineage law sits at z = −αs. One likelihood evaluation
took 0.5, 6.5, 53 and 400 s at s = 0.96, 1.16, 1.40 and 1.69, and then failed. At those points the top-bin birth rate
λe^{5s} was 6×10³, 4×10⁴, 3×10⁵ and 2×10⁶ per Myr.

## L.2 The mean-rate model

**Definition.** Replace (0.1) by

$$dx = -\theta\,\big(x - (\alpha - \tfrac12)s^2\big)\,e^{\alpha x}dt + \sqrt{2\theta s^2 e^{\alpha x}}\,dW. \qquad (L.1)$$

This is the §9.1 family with OU mean m(s) = (α − ½)s².

**Proposition L.1.** For every α, the lineage-time stationary law of (L.1) is N(−s²/2, s²). Hence E_π[e^{x}] = 1, and
E_π[e^{ex}] = 1 for every exponent e ∈ {0, 1}.

*Proof.* The generator is e^{αx}𝓛, with 𝓛 an OU generator with stationary density N(x; m, s²). A time-changed diffusion
has stationary density proportional to π_𝓛(x)e^{−αx} (§5.2). Here
exp(−(x − m)²/(2s²) − αx) is proportional to exp(−(x − m + αs²)²/(2s²)), that is N(m − αs², s²) = N(−s²/2, s²).
The lognormal mean is then e^{−s²/2 + s²/2} = 1. ∎

**Parameters.** The per-bin rates keep their form: λ_k = λe^{βx_k}, μ_k = μe^{βx_k}, ψ_k = ψe^{δx_k},
ℓ_e = c Σ_k e^{γx_k} d_{e,k}. Because every exponent is 0 or 1, **each rate parameter is its own lineage-time mean**.
That is the ergodic average along one non-branching lineage in tempo equilibrium. It is not the mean over
contemporaneous species (p\*), nor over the lineages of a reconstructed tree.
- λ and μ are the mean speciation and extinction rates of a lineage in tempo equilibrium;
- c is the mean clock rate. In cell (0, 1), c and s are the mean and log-sd of the per-unit-time marginal rate law.
  Those are the parameters the comparators' UCLN priors describe (`ucld_mean`, `ucld_sd`). Rates here are
  autocorrelated along lineages, so branch-average rates are not lognormal with sd s; the match is to the marginal law,
  not to UCLN's branch rates;
- θ is the mean real-time reversion rate when α = 1;
- ψ is the mean fossil-sampling rate when δ = 1.

The couplings β, γ, δ and α keep their meaning. With β = 0, λ and μ are unchanged, and cell (0, 0) is unchanged
entirely. The median rate is λe^{−s²/2}. If a coupling γ ∈ (0, 1) is added later, E_π[e^{γx}] = e^{−γ(1−γ)s²/2}, and c
is no longer exactly the mean.

**Grid.** The z-window [z_min, z_max] (default [−5, 5]) is centred on the lineage law, and

$$x_k = s z_k - s^2/2, \qquad z_k = z_{\min} + (k-1)h.$$

`TempoGrid()` gains this offset. In z, the tempo-time OU has mean αs, so the generator depends on s through its drift
centre, (𝓛_z f)(z) = −θ(z − αs)f′ + θf″. Consequences 1 and 2 of §1.1 no longer hold: A depends on s through αs, and x
through −s²/2 as well as e^{sz_k}. Consequence 3 still holds: paths stored as bin indices stay valid when s moves.

**Discretisation.** The central scheme's check is now h · max_k |z_k − αs| ≤ 2. With α > 0, at K = 26 on [−5, 5]
(h = 0.4), the check fails for every s > 0. Positivity itself, with the binding rates at the reflecting ends discarded,
holds for s ≤ 0.4; at α = 0 the bound is met for every s. Under `"mean"` the bound depends on s, so a chain with free s
can leave the valid region mid-run. **SqRA therefore becomes the default scheme.** SqRA applied to A with potential
U(z) = (z − αs)²/2, then multiplied by e^{αx_k}, gives

$$Q_{k,k\pm1} = \frac{\theta}{h^2}\exp\!\Big(\frac{\alpha(x_k + x_{k\pm1})}{2} - \frac{z_{k\pm1}^2 - z_k^2}{4}\Big). \qquad (L.2)$$

The rates are always positive and the rows sum to zero. Detailed balance gives Q_{k,k+1}/Q_{k+1,k} = φ(z_{k+1})/φ(z_k), so
**the stationary law of Q is exactly the grid Gaussian** π_k ∝ φ(z_k).

Expanding (L.2) to first order in h, the jump moments are h(Q_{k,k+1} − Q_{k,k−1}) = θe^{αx_k}(αs − z_k) + O(h²) and
h²(Q_{k,k+1} + Q_{k,k−1}) = 2θe^{αx_k} + O(h²). Both are consistent with (L.1).

On the grid, E_π[e^{x}] = Σ_k φ(z_k)e^{sz_k − s²/2} / Σ_k φ(z_k). This differs from 1 by the truncation of the
e^{x}-tilted law N(s, 1) above the top of the grid. At K = 26 on [−5, 5] the grid values are −1.2×10⁻⁵ at s = 1,
−6.4×10⁻⁴ at s = 2 and −1.3×10⁻² at s = 3. The continuum tail 1 − Φ(5 − s) gives −3×10⁻⁵, −1.3×10⁻³ and −2.3×10⁻²: the sum
behaves as the integral to z_max + h/2, so the endpoint term is the same order as the tail, not exponentially small. A grid-exact offset, −log of that sum, is a possible variant. The closed form is preferred for a
simple map and Jacobian.

The other schemes stay available:
- *central*, with drift −θ(z_k − αs). When the bound fails it signals an infeasible proposal, so MCMC moves reject
  it, and it is an error in a direct call. With free s under `"mean"`, this truncates the support of s to the bound.
- *upwind*, which is always positive and O(h).

**Start priors.**
- `equilibrium` keeps its definition (the Perron vector of Qᵀ + R).
- `stationary` is now the grid Gaussian, exactly so only under SqRA.
- `uniform` is unchanged.
- A numeric per-bin `rootPrior` sits in a different law, so existing numeric priors must be reviewed.

## L.3 Relation to the current model

**Exact map (continuum).** By §9.1 with a = m(s) = (α − ½)s², the current model (m = 0) at (θ₀, λ₀, μ₀, ψ₀, c₀, s) and
the mean-rate model at (θ, λ, μ, ψ, c, s) define the same law of trees, samples and characters when

$$\theta = \theta_0 e^{-\alpha m},\quad \lambda = \lambda_0 e^{-\beta m},\quad \mu = \mu_0 e^{-\beta m},\quad
\psi = \psi_0 e^{-\delta m},\quad c = c_0 e^{-\gamma m}, \qquad (L.3)$$

with x = x₀ + m. At α = 1, m = s²/2 and every coupled rate is multiplied by e^{−s²/2}. At α = 0 (the α ≡ γ convention
with γ = 0), m = −s²/2. s, ρ and the couplings are unchanged.

**Exact map (grid).** The mean-rate model on [z_min, z_max] equals the current model on the window
[z_min − αs, z_max − αs] at the parameters (L.3), for every scheme. Both generators share the same jump rates in
x-space and the same per-bin rates and start priors. It is not equal to the current model on the unshifted window; the
difference is the truncation. This identity is the main unit test (L.6 V1). It also means the new model can be
evaluated **with the current code**, by passing `zRange = c(-5, 5) - alpha * s` (L.6 V0).

**Jacobian.** In coordinates (log θ, log λ, log μ, log ψ, log c, s), (L.3) is a shear with unit determinant. In natural
coordinates, |∂(new)/∂(old)| = exp(−m(s) Σ_p e_p) over the free rate parameters p. Here e_p is α for θ, β for λ and μ,
δ for ψ and γ for each clock rate. A run of the current model under prior p₀ is therefore reproduced, on the shifted
window, by the new model under

$$p(\theta,\lambda,\mu,\psi,c,s) = p_0\big(\text{(L.3)}^{-1}\big)\cdot\exp\!\big(m(s) \textstyle\sum_p e_p\big).$$

This prior couples s to the rates. Independent priors on the new parameters define a different, and intended, joint
prior. Reproducing an old run also needs s fixed, or every coupled parameter free: a fixed θ₀ becomes the s-dependent
constraint θ = θ₀e^{−αm(s)}. Posteriors and Bayes factors convert between the parameterisations only under this mapped prior.

**Bayes factors.** At equivalent parameter points the likelihoods agree; the evidence changes only through the prior.
With priors stated on the new parameters, every rate prior refers to the same quantity (the lineage-mean rate) in
every (β, γ) cell, and matches the comparators' mean-rate priors. Under the current model a prior on c₀ in a γ = 1 cell
sits e^{γm} above the mean clock (e^{s²/2} at α = γ = 1) that the comparator's prior describes; this is the "comparator-matched at unit tempo
only" caveat of the trilobite README (F8).

## L.4 What the change does and does not do

**Does.**
- The rate parameters become lineage means, whose meaning does not depend on s. The current parameters refer to x = 0,
  which lies αs lineage-sd above the law's centre.
- At fixed identified rates, the top-bin rates fall. The window top moves from x = 5s to x = 5s − s²/2, and the
  parameters are rescaled by (L.3). The net factor at α = 1 is e^{−βs²} on the top-bin birth rate (1/17 at s = 1.69,
  1/127 at s = 2.2) and e^{−s²} on the θe^{αx_K} factor of the top row of Q. The A row changes too, through the drift
  centre.
- The low-end truncation goes.

**Does not.**
- **It does not make the pilot's drift vanish.** The data pin rates at the paths' level (log-sd 0.17–0.27 above). Over
  the nine pilot segments, the candidate parameters still vary:

  | Candidate | log-sd of λ | log-sd of θ |
  |---|---|---|
  | Current | 0.99 | 0.69 |
  | Mean-rate (this spec) | 0.64 | 0.56 |
  | Median-centred (L.5 A) | 0.69 | 0.89 |

  The paths sit at an O(s) offset from any stationary reference point, set by the root prior, the conditioning and
  the slow reversion (θ⟨e^x⟩ ≈ 0.003/Myr against a tree depth of about 500 Myr). The pilot does not separate the
  mean-rate and median-centred choices. The choice rests on interpretation and comparator matching.
- **It does not bound s.** The 1/λ brake was an artefact of λ's meaning. Now s is held only by its own prior and the
  likelihood, and the s prior must be stated on its merits. Under this parameterisation s is the log-sd of the marginal rate
  law, as UCLN's σ is (see L.2 on the difference). The trilobite and penguin pipelines already use the comparators'
  UCLN σ prior.
- **It does not remove intrinsic stiffness, and its speed-up has a cost.** The top-bin rates still grow as
  e^{αsz_max}. The gain comes from no longer representing log-tempo more than 5 lineage-sd above the lineage law. But the
  species equilibrium p\* (the default root prior) and the e^{βx}-tilted lineage law also sit above it. At the pilot's
  parameters p\* is centred near z = +1.5 at s = 2 and z = +2.4 at s = 3 (from `ridge-location.csv`, shifted by s). The
  tilted law is N(s, 1). The window must still cover both; if it does not, z_max must grow (an asymmetric window is
  allowed, §1.6), and part of the speed-up is lost. V0 checks this.

## L.5 Candidates considered

| | Reference point | Verdict |
|---|---|---|
| **Mean rate (adopted)** | Lineage law N(−s²/2, s²); OU mean (α − ½)s²; window centred | Every rate is a lineage mean, matching the comparators' conventions; best of the stationary choices on θ in the pilot, tied on λ |
| A, median-centred | Lineage law N(0, s²); OU mean αs²; window centred | The same model on the same window; differs only in which label carries the prior (medians). Equally defensible; the pilot does not separate it from the mean-rate choice. Second choice |
| B, relabel only | Mean or median on the old drift and old window | Same target as fix C's arm A, which stalled; numerics unchanged. Rejected |
| D, species equilibrium p\* | The map depends on (λ, μ, θ) | No closed-form Jacobian; the paths do not sit at p\* either. Rejected |
| E, current model plus priors (fix C) | — | Superseded by the maintainer's decision of 25 Nov |

## L.6 Implementation and validation plan

### Interface

- One argument, `location = c("mean", "bm")`, on `DCETModel()`, `TempoGrid()`, `TempoGenerator()`, `SimulateCET()` and
  any other constructor of the grid, Q or the continuous process. The default is `"mean"`; `"bm"` is the current model
  (m = 0, x = sz).
- `DCETModel()`'s default scheme becomes `"sqra"` for `"mean"` and stays `"central"` for `"bm"`.

### R engine (the reference)

- `TempoGrid()`: x_k = sz_k − s²/2 under `"mean"`.
- `.OUBlock()` takes the drift centre ζ = αs (0 under `"bm"`): U = (z − ζ)²/2 for SqRA, drift −θ(z − ζ) for central
  and upwind, and the centre-aware positivity check.
- `.Setup()` passes ζ.
- `SimulateCET()`'s exact OU in tempo time (`.OUSegment()`) gets mean (α − ½)s². The Gillespie dCET simulator uses Q
  and the grid, and needs no change.
- `TempoEquilibrium()`, `.RootPrior()`, `BinRates()`, the tree moves and the shift move (`.Shifted()`: bin step
  d = hs, unchanged map and Jacobian) take Q, rates and the grid from `.Setup()`.
- Any code that rebuilds the grid or Q outside `.Setup()` must pass the location. Audit with
  `grep -n "TempoGrid(\|TempoGenerator(\|\.OUBlock(" R/`.
- Documentation: `DCETModel()`, `TempoPriors()` (what the rates mean), the AGENTS.md model summary, and the
  empirical READMEs' "comparator-matched at unit tempo" notes.

### C++ engine

- The C++ engine receives Q, the grid x and the per-bin rates from R. Confirm that no kernel rebuilds them:
  `src/tree_moves.cpp` `.Rebuild()` and the slice set-up.
- C++ must agree with R to ≤ 1e−8 on the rung 3/4 fixtures, rerun at the new default.
- Seeded `TempoMCMC()` chains must be identical in both engines (`test-Engine.R`).

### Validation

**V0 (go/no-go, before any code).** Rerun the `diag-ridge.R` profile (the six points of `ridge-profile.csv`, on the
chain's last tree) with `DCETModel(scheme = "sqra", zRange = c(-5, 5) - s)` at each point's parameters. By L.3 this is
the new model at the mapped parameters. At each point, report:
1. the log-likelihood and seconds per evaluation, beside the current grid;
2. the log-likelihood on a reference window and grid, for example K = 51 on [−7, 7] − s, with the difference from (1);
3. the masses in the top bin and the top three bins of the root prior p\*, of the e^{x}-tilted lineage law and, where
   available, of the root-bin posterior.

Go when the points up to s = 2.2 evaluate in under about 60 s, agree with the reference to < 0.05 log units, and have
top-three-bin masses below 10⁻³. Otherwise, revise the window (z_max) or K in this spec before implementing.

**V0 outcome (2026-11-26, #123; decision pending with the maintainer).** V0 was run and #123 implemented without
revising the spec. That is a departure from the rule above, and it is recorded here.

| Criterion | Result |
|---|---|
| Time | Met to s = 2.04 (46-76 s); s = 2.46 takes 117-187 s |
| Accuracy | Not met. K = 26 against K = 51 on the same window differs by 0.125, 0.130, 0.189 and 0.235 log units at s = 0.96, 1.40, 2.04, 2.46. The old default's error at the same points is 0.200 and 0.182 |
| Top-three mass, root prior | Met (≤ 7.3×10⁻⁴) |
| Top-three mass, tilted law | Not met, and **unmeetable as written**: under SqRA the tilted law is exactly φ(z_k)e^{sz_k}, whose top-three mass at K = 26 on [−5, 5] is about 1 − Φ(4 − s). That exceeds 10⁻³ for every s > 0.93 |

Two faults in the gate as written:
- **The reference.** K = 51 on [−7, 7] − s has h = 0.28, is itself 0.04-0.05 from the h = 0.2 value, and does not
  solve at s ≥ 2. A usable reference is specified by spacing: h ≤ 0.2 on the same window, plus a wider window at the
  same h.
- **The tilted-law criterion.** It is a condition on z_max (about z_max ≥ s + 4.1 at h = 0.4), not on the model.

The options (raise K, widen or skew the window, or relax the thresholds with the reason stated) are set out on #117 for
the maintainer. Whichever is chosen applies to every dataset.

**V1 (unit tests, both engines).**
1. The L.3 grid identity, mean-rate against `"bm"` on the shifted window, for TreeLikelihood at 1e−10 relative. Cover
   α, β, γ, δ ∈ {0, 1}, `start` root and origin, and every scheme.
2. The SqRA stationary law equals φ(z_k) normalised to 1e−12, and E_π[e^{x}] matches its grid formula.
3. The β = 0 collapse to the analytic FBD (§2) still holds.
4. `"bm"` reproduces current results exactly: seeded chains are bit-identical to those before the change.
5. Coverage: 100% of the new lines.

**V2 (rungs; rerun at the default grid, `EVOTEMPO_GRID=default`).**
- Rung 1 is unchanged under `"bm"`.
- Rung 2: the continuous simulator under `"mean"` against the numerics (lineage mean of e^{x} equal to 1, survival,
  N_x), plus the negative control.
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
- Priors: the comparator's FBD-block prior on the new (λ, μ, ψ), now lineage-mean rates as in the comparator; for s,
  the comparator's UCLN σ prior, as in the trilobite pipeline.
- Report: λ, s, λe^{⟨x⟩}, seconds per segment and R-hat.
- Pass when: window medians of λ and s show no monotone drift; segment cost stays within twice that of cell (0, 0);
  and s does not pile against its support edge.
- Recheck that cell (1, 1) is qualitatively unchanged.

**Scope (maintainer, 2026-11-25).** The aim is one model setup, with the same parameterisation and the same moves, for
every dataset; there are to be no per-dataset variants. The mean-rate model therefore becomes the model for penguins,
trilobites and Hymenoptera alike:
- penguins are rerun under it, together with the rerun that the new tree moves (#110) require anyway;
- trilobite production waits for it;
- `location = "bm"` stays only as the Budd & Mann reproduction fixture for rung 1, never as an option for an empirical
  analysis.

**Legacy objects.** A model or grid without a `location` (saved before #116) is read as `"bm"` with its stored scheme,
and raises a `legacyLocation` warning. The empirical pipelines' checkpoints record `location` and `scheme`, so:
- resuming a chain of another location is refused;
- scoring one in a Bayes-factor eval is refused;
- a checkpoint from before #116 is refused, because it records no location.

If V0 or V5 shows that the mean-rate model does not work for one dataset, this spec is revised for all of them. No
dataset-specific switch is added. The `CET_FBD_SCALE` and `CET_S_MAX` switches on `hymenoptera-cet-pipeline` are removed
when that pipeline moves to the new model.
