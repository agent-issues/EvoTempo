<!-- Copied from the project artifact store on 2026-10-12. {{artifact:...}} links refer to that store and do not resolve on GitHub. -->

# Discretised Covariant Evolutionary Tempo (dCET): mathematical specification

*Groundwork document, 2026-10-12; version 2 (B&M code inspected). Specification only; no package code. Companions: [precedents-review.md]({{artifact:b04de770-aebf-4f00-8545-6bbe2a444857}})
(literature) and [dryad-numerics-notes.md]({{artifact:c1c85348-fc5e-40ef-a611-6d5954eb2b55}}) (reference numerics written from B&M's Dryad code, which was inspected and partly re-run; see §11 and the notes).
Every number labelled "scratch" was computed for this document with a throw-away script
([dcet_scratch_checks.py]({{artifact:ec1e79d0-9b19-419f-9654-43660a095567}}); output in [dcet_scratch_output.txt]({{artifact:ff706844-85ce-4af4-bb6e-4b1ceae37fdc}})) that implements only the formulas below. These numbers are provisional reference
values, not validated outputs.*

**Notation.** B&M = Budd & Mann 2025, Syst Biol 74:917–934; equation numbers "B&M (n)" / "A.n" refer to it.
Real time runs forward as t. **Age** (time before present) is τ, so τ = 0 is the present. On a branch e the parent
(older) end is at age τ_e^p and the child end at age τ_e^c < τ_e^p. Vectors over tempo bins are bold or carry index k.
diag(v) is the diagonal matrix of v and **1** the all-ones vector.

---

## 0. Model in one paragraph

Each lineage carries a log-tempo x(t). Its rates are speciation λe^{βx}, extinction μe^{βx}, fossil sampling ψe^{δx}, and
character change c·e^{γx}. Log-tempo evolves by B&M's modified OU process, generalised with a self-rate exponent α:

$$dx = -\theta x\,e^{\alpha x}dt + \sqrt{2\theta s^2 e^{\alpha x}}\,dW. \qquad (0.1)$$

B&M's CET is (α, β, γ, δ) = (1, 1, 1, ·). They do not model fossils, so δ is new. dCET replaces x by a K-state CTMC on a grid.
Inference augments the tempo path on the branches of the reconstructed (sampled) tree, marginalises hidden lineages
through extinction probabilities E_k(τ), and scores characters by pruning on branch lengths
ℓ_e = c Σ_k e^{γx_k}·(dwell of e in bin k).

**Fixed choices** (from memo §2.7 and this document): the OU mean is fixed at 0 (§9.1); the removal-after-sampling probability
r_rm is fixed, default 0 (Truman et al.: unidentifiable if free); ρ is fixed from external richness. α = 1 by default; see §7.3 and open issue 1.

---

## 1. The generator (spec item a)

### 1.1 Grid

Let z = x/s. In z the B&M process (α = 1) is

$$dz = -\theta z\,e^{s z}dt + \sqrt{2\theta e^{s z}}\,dW,$$

and its infinitesimal generator acting on smooth f is

$$(\mathcal{G}f)(z) = e^{s z}\left[-\theta z f'(z) + \theta f''(z)\right] = e^{s z}(\mathcal{L}f)(z). \qquad (1.1)$$

The OU part 𝓛 does not involve s. **Bins are equally spaced in z:** z_k = z_min + (k−1)h, k = 1..K, h = (z_max − z_min)/(K−1),
with x_k = s z_k. Consequences:
1. The OU block A (below) is independent of s.
2. s enters only through the state rates e^{s z_k}·(…) and e^{βs z_k}, e^{γs z_k}.
3. A move on s keeps a sampled path valid, because paths are stored as bin indices.

z_min and z_max need not be symmetric (§1.6).

### 1.2 Three admissible discretisations of 𝓛

Write the tridiagonal generator A with up-rate u_k = A_{k,k+1}, down-rate d_k = A_{k,k−1}, A_kk = −(u_k + d_k), and all other
entries 0. The **reflecting ends** are d_1 = 0 and u_K = 0.

| Scheme | u_k | d_k | Off-diagonals ≥ 0 | Order | Exact discrete identity |
|---|---|---|---|---|---|
| **Upwind** (memo baseline) | θ/h² + θ·max(−z_k,0)/h | θ/h² + θ·max(z_k,0)/h | always | O(h) | (A z)_k = −θz_k (interior) ⇒ A.28 exact |
| **Central** | θ/h² − θz_k/(2h) | θ/h² + θz_k/(2h) | iff h·max_k\|z_k\| ≤ 2 | O(h²) | (A z)_k = −θz_k (interior) ⇒ A.28 exact |
| **SqRA** (square-root approximation) | (θ/h²)·exp(−(U_{k+1}−U_k)/2) | (θ/h²)·exp(−(U_{k−1}−U_k)/2) | always | O(h²) | stationary law ∝ e^{−U_k}, U = z²/2 ⇒ A.26 exact on the grid |

**Proof that off-diagonals are non-negative.**
1. *Upwind.* Each rate is θ/h² > 0 plus a term θ·max(±z_k, 0)/h ≥ 0.
2. *SqRA.* Each rate is a positive constant times an exponential, hence > 0.
3. *Central.* u_k ≥ 0 ⇔ θ/h² ≥ θz_k/(2h) ⇔ z_k h ≤ 2, and d_k ≥ 0 ⇔ −z_k h ≤ 2. Both hold for all k iff h·max_k|z_k| ≤ 2.

In every case the row sums are zero by construction, so A is a valid CTMC generator. The full generator is

$$Q = \mathrm{diag}(e^{\alpha s z_k})\,A, \qquad Q_{kj} = e^{\alpha x_k}A_{kj}. \qquad (1.2)$$

Left-multiplying by a positive diagonal preserves the signs of the off-diagonals and zero row sums, so Q is a CTMC generator. ∎

**Consistency with B&M (1) as K → ∞ (h → 0, fixed z-range).** For smooth f and an interior bin, the first two jump moments are:
1. *Drift.* h(u_k − d_k) = −θz_k **exactly** for central and upwind. For SqRA, expand ΔU_± = ±z_k h + h²/2, so
   u_k − d_k = −(2θ/h²)e^{−h²/4} sinh(z_k h/2), and h(u_k − d_k) = −θz_k(1 + O(h²)).
2. *Diffusion.* h²(u_k + d_k) = 2θ exactly for central, 2θe^{−h²/4}cosh(z_k h/2) = 2θ(1 + O(h²)) for SqRA, and 2θ + θ|z_k|h for
   upwind. The extra θ|z_k|h is first-order numerical diffusion, the source of the upwind bias in §10.
3. *Higher jump moments* are O(h).

Hence (A f)(z_k) = −θz_k f′ + θf″ + O(h^p), with p = 1 for upwind and 2 otherwise. Then (Qf)(z_k) = e^{αsz_k}(𝓛f)(z_k) + O(h^p), so Q
is consistent with 𝓖, and with B&M (1) at α = 1 in x-units (dx = s dz). Non-negative off-diagonals with zero row sums make
e^{Qt} a stochastic matrix, so the semi-discrete linear equations (4), (9), (11), (A.15) are max-norm stable, and the
reaction-term Riccati equations (2), (2.1) preserve [0, 1]. Consistency plus stability then gives convergence on the truncated
domain (Lax-equivalence argument for the linear equations; a monotone-scheme comparison argument for the Riccati ones). The truncation error from z_max is separate (§1.6).

**Reflecting boundaries correspond to B&M's no-flux condition.** In the backward (Kolmogorov) equations the end rows read
(A f)_1 = u_1(f_2 − f_1) and (A f)_K = d_K(f_{K−1} − f_K). These are the discrete homogeneous Neumann condition ∂f/∂x = 0, which
B&M impose on G at x_min, x_max. In the forward (Fokker–Planck) equations, 1ᵀQᵀ = 0 conserves mass: zero flux.

**Reversibility.** Any irreducible tridiagonal generator is a birth–death chain and satisfies detailed balance. A has
stationary law π^A_k ∝ Π_{j<k} u_j/d_{j+1}, and Q has π^Q_k ∝ π^A_k e^{−αx_k} (§5.2). So Q, and the tilted generators
Q + diag(g) of §3, are similar to symmetric matrices via diag(√π^Q). This guarantees real spectra and stable
eigendecompositions (direct endpoint-conditioned sampling, §8.2).

### 1.3 Recommended default

- **Central**, with K chosen so that h ≤ 2/max(|z_min|, |z_max|). It is second order, and it satisfies the discrete A.28 identity
  exactly (§5.3), which makes the r = θ checks sharp.
- **SqRA** wherever the positivity bound would force very large K, for example wide asymmetric ranges.
- **Upwind** only as a reference: the scratch runs (§10) show its O(h) numerical diffusion is large enough to bias headline
  quantities by about 13% at h = 0.2.
- B&M's own scheme is *not* this symmetric upwind scheme. It is a fixed backward difference: upwind for x > 0, downwind for x < 0, with
  effective diffusion θ(s² + xh/2) and negative rates for x < −s²/h. See §11 item 7 and the numerics notes. dCET must be able to reproduce it
  exactly (as an optional 'B&M operator') for rung 1.

### 1.4 Scratch verification of the generator identities

At θ = 0.01, s = 1, z ∈ [−8, 8] ([dcet_scratch_checks.py]({{artifact:ec1e79d0-9b19-419f-9654-43660a095567}})):

| Scheme, K | mean x at r = 0 (theory −s² = −1) | mean x at r = θ (theory 0) | A.28 residual at r = θ |
|---|---|---|---|
| upwind 81 / 161 / 321 | −1.1242 / −1.0602 / −1.0296 | 0.0000 | ≤ 1e−11 |
| central 81 / 161 / 321 | −0.9967 / −0.9992 / −0.9998 | 0.0000 | ≤ 1e−13 |
| SqRA 81 / 161 / 321 | −1.0000 / −1.0000 / −1.0000 | +0.0037 / +0.0009 / +0.0002 | −5.4e−5 / −1.3e−5 / −3.3e−6 |

### 1.5 Stiffness

The largest exit rate is q_max = max_k e^{αx_k}(u_k + d_k) ≈ e^{αs z_max}(2θ/h² + θz_max/h).
- At B&M's range (z_max = 10, s = 1), θ = 0.01 and h = 0.1: q_max ≈ 2.2×10⁴ × 3 ≈ 7×10⁴ per Myr.
- At an inference range of z_max = 5 with h = 0.2: q_max ≈ 148 × 0.75 ≈ 110 per Myr.

E- and D-solvers must be implicit (Radau/BDF/TR-BDF2/Rosenbrock), or exponential integrators on a fixed grid.

### 1.6 Range

- The truncation caps rate multipliers at e^{s z_max}.
- *B&M's x ∈ (−10, 10) range.* Scratch (§10): for clades started at x_0 = 0, mean and experienced clade sizes at T = 500 Myr are identical
  to 4 significant figures for z_max ∈ {5, 8, 10}. The upper tail is therefore not driving Fig 1A.
- *Inference range.* [z_min, z_max] ⊇ [m_π − 5sd_π, m_π + 5sd_π] for the relevant equilibrium π (§5), plus root-state
  coverage. Asymmetric ranges are allowed. Ancestral tempos are shifted upwards (B&M Fig 4A), so z_max matters more than z_min.

---

## 2. Extinction probabilities with extant and fossil sampling (spec item b)

**Rates in bin k:** λ_k = λe^{βx_k}, μ_k = μe^{βx_k}, ψ_k = ψe^{δx_k}, with δ ∈ {0, 1} (δ = 0: constant ψ; δ = 1: ψe^x), and
Q from (1.2). Extant lineages are sampled independently with probability ρ at τ = 0. A sampled fossil lineage is removed with
probability r_rm and otherwise continues.

**Definition.** E_k(τ) = P(a lineage in bin k at age τ leaves no sampled descendant: no fossil sample on it or on any
descendant, and no sampled extant descendant).

**Derivation.** Condition on what happens in the first dτ of real time after age τ + dτ, i.e. in the age interval [τ, τ+dτ]:
1. extinction, probability μ_k dτ → 1;
2. fossil sampling, ψ_k dτ → 0 (the sample is itself a sampled descendant, whatever r_rm);
3. speciation, λ_k dτ → E_k(τ)²;
4. tempo jump to j, Q_kj dτ → E_j(τ);
5. nothing, 1 − (λ_k+μ_k+ψ_k+q_k)dτ → E_k(τ), where q_k = −Q_kk.

Hence

$$\frac{dE_k}{d\tau} = \mu_k - (\lambda_k+\mu_k+\psi_k)E_k + \lambda_k E_k^2 + \sum_j Q_{kj}E_j,\qquad E_k(0) = 1-\rho. \qquad (2.1)$$

Since Σ_j Q_kj = 0, the last term equals Σ_{j≠k} Q_kj(E_j − E_k), which is the memo's form. r_rm does not enter E.

**Positivity.** [0,1]^K is forward-invariant:
- at E_k = 0, dE_k/dτ = μ_k + Σ_{j≠k} Q_kj E_j ≥ 0;
- at E_k = 1, dE_k/dτ = −ψ_k + Σ_{j≠k} Q_kj(E_j − 1) ≤ 0.

Solvers should preserve this.

**Relation to B&M (2).** With ρ = 1 (so E(0) = 0), ψ = 0 and β = α = 1, (2.1) is the semi-discretisation of B&M (2) at z = 0:
E_k(τ) ≈ G_{x_k}(τ, 0). Also:
- E = μ/λ·**1** is a fixed point whenever β is constant across bins and ψ = 0, because μ − (λ+μ)μ/λ + λ(μ/λ)² = 0 and Q**1** = 0.
  This is the tempo-invariant long-run extinction probability μ/λ (B&M p. 919). §9.3 explains why it is tempo-invariant.
- *Extant-descendant variant.* Ẽ_k(τ) = P(no sampled *extant* descendant) solves (2.1) with ψ_k set to 0 (fossil sampling does
  not affect extant descendants when r_rm = 0). It is needed for crown conditioning (§6). If r_rm > 0, a removing sample acts as an
  extinction for this purpose. Replace μ_k by μ_k + r_rm ψ_k and set ψ_k to 0 in the Ẽ equation.

**Special case β = 0, δ = 0.** All λ_k, μ_k, ψ_k are equal. The constant vector E(τ)**1** with E solving the scalar FBD
Riccati equation satisfies (2.1) exactly, because Q**1** = 0. This is Stadler's (2010) p₀ in closed form, an exact test of the solver.

---

## 3. D-equation, augmented path density, and the sum-over-paths identity (spec item c)

### 3.1 D-equation

D_{e,k}(τ) is the probability density of the observed subtree descending from the point at age τ on branch e (topology, sampling
events, sampling times, and later the characters), given the lineage is in bin k at τ. Over [τ, τ+dτ]:
1. no event → D_k(τ);
2. a jump k → j → D_j(τ);
3. a hidden speciation whose daughter subtree is unsampled → λ_k dτ · 2E_k(τ) D_k(τ). The factor 2 counts which (labelled)
   daughter carries the observed subtree;
4. extinction or fossil sampling at an unobserved point → 0.

Hence

$$\frac{dD_k}{d\tau} = \big(-(\lambda_k+\mu_k+\psi_k) + 2\lambda_kE_k(\tau)\big)D_k + \sum_j Q_{kj}D_j. \qquad (3.1)$$

In vector form, d**D**/dτ = (Q + diag(**g**(τ)))**D** with g_k(τ) := −(λ_k+μ_k+ψ_k) + 2λ_kE_k(τ). The equation is linear in **D**, so

$$\mathbf{D}_e(\tau^p_e) = \mathcal{M}_e\,\mathbf{D}_e(\tau^c_e),\qquad \partial_\tau \mathcal{M}(\tau,\tau^c_e) = (Q+\mathrm{diag}\,\mathbf{g}(\tau))\,\mathcal{M}(\tau,\tau^c_e),\ \mathcal{M}(\tau^c_e,\tau^c_e)=I. \qquad (3.2)$$

This is MTFBD's q-equation (precedents §7) with the uniform type-change kernel replaced by Q. It is the CET analogue of the
QuaSSE D-equation with state-dependent diffusion.

### 3.2 Augmented path density on one branch

A tempo path on branch e is piecewise constant. It starts in k_0 at the parent end and jumps k_{i−1} → k_i at ages
τ^p_e > a_1 > … > a_m > τ^c_e (forward in real time from parent to child). Writing X(τ) for the bin at age τ, define

$$p_e(\text{path}) = \prod_{i=1}^{m} Q_{k_{i-1}k_i}\;\exp\!\Big(-\!\int_{\tau^c_e}^{\tau^p_e} q_{X(\tau)}\,d\tau\Big)\;\exp\!\Big(\int_{\tau^c_e}^{\tau^p_e} g_{X(\tau)}(\tau)\,d\tau\Big). \qquad (3.3)$$

The first two factors are the CTMC path density under Q. The third is the probability that the lineage experiences no
*observed* event and that every hidden speciation's side subtree is unsampled. It equals
exp(−∫(λ+μ+ψ)_X) · Σ_n (∫2λ_XE_X)^n/n!. This matches memo §2.5.

### 3.3 Proposition (sum over paths recovers D)

For every pair of bins k, j,

$$\mathcal{M}_{kj}(\tau^p_e,\tau^c_e) = \sum_{m\ge0}\ \sum_{k_1,\dots,k_{m-1}}\ \int\!\!\cdots\!\!\int_{\tau^p_e>a_1>\dots>a_m>\tau^c_e} p_e(\text{path with } k_0=k,\ k_m=j)\,da_1\cdots da_m.$$

Consequently Σ_j (sum over paths ending in j) × D_{e,j}(τ^c_e) = D_{e,k}(τ^p_e).

*Proof.* Split Q + diag(**g**(τ)) = J + Λ(τ). Here J is the off-diagonal part of Q (entries Q_kj, k ≠ j, all ≥ 0) and Λ(τ) is
diagonal with Λ_kk(τ) = −q_k + g_k(τ). Let U(a, b) be the diagonal propagator U_kk(a, b) = exp(∫_b^a Λ_kk(τ)dτ). It solves
(3.2) with J = 0. The Dyson (interaction-picture) expansion of the time-ordered exponential is

$$\mathcal{M}(\tau^p,\tau^c) = \sum_{m\ge0}\int_{\tau^p>a_1>\dots>a_m>\tau^c} U(\tau^p,a_1)\,J\,U(a_1,a_2)\,J\cdots J\,U(a_m,\tau^c)\,da_{1:m}.$$

It converges absolutely: ‖J‖ < ∞ for finite K, and the m-th term is bounded by (‖J‖Δ)^m e^{‖Λ‖Δ}/m!. Writing out the (k, j)
entry of the m-th term as a sum over intermediate indices gives

Π_i U_{k_{i−1}k_{i−1}}(a_{i−1}, a_i) · Π_i J_{k_{i−1}k_i}, with a_0 = τ^p and a_{m+1} = τ^c.

This is exactly (3.3) integrated over jump times. The left factor carries the oldest state, so the index order matches forward-time
jumps k_{i−1} → k_i at rate Q_{k_{i−1}k_i}. Multiplying by **D**(τ^c) and using (3.2) proves the claim. ∎

*Scratch check.* For K = 4, a random Q, constant g and Δ = 0.7, a Monte Carlo average of exp(∫g_X) over 5×10⁴ prior paths per
start state reproduced expm((Q+diag g)Δ) entrywise to within 1.5×10⁻³ (maximum entry 0.29). This is consistent with Monte Carlo error.

### 3.4 Full augmented density

Let 𝒯 be the reconstructed, sampled, time-calibrated tree, with:
- bifurcation nodes 𝒱 (including the root if root-conditioned);
- extant tips (n_0 of them);
- fossil tips 𝓕;
- sampled ancestors 𝒮 (degree-2 nodes with a zero-length tip; the path passes through them continuously);
- branches 𝓔 (an SA splits its lineage into two branches).

The path state at a node is shared by all incident branches, since daughters inherit the parent's tempo.

$$p(\mathcal{T},\text{paths}\mid\Theta) = \frac{\pi(k_{\text{start}})\ \prod_{e\in\mathcal E} p_e(\text{path}_e)\ \prod_{v\in\mathcal V}\lambda_{k(v)}\ \rho^{\,n_0}\ \prod_{f\in\mathcal F}\psi_{k(f)}\big[r_{rm} + (1-r_{rm})E_{k(f)}(\tau_f)\big]\ \prod_{a\in\mathcal S}\psi_{k(a)}(1-r_{rm})}{\text{conditioning (§6)}}. \qquad (3.4)$$

Here k_start is the bin at the origin (stem start) or at the root. Then:
- **Marginal tree likelihood.** Summing (3.4) over all paths and applying Proposition 3.3 branch by branch, postorder, gives the standard SSE pruning
  recursion with tip values D = ρ (extant tips) or ψ_k[r_rm + (1−r_rm)E_k] (fossil tips), node combination
  D_k ← λ_k D_{left,k} D_{right,k}, SA factor ψ_k(1−r_rm) times continuation, and root value Σ_k π_k D_k.
  This is the quantity that must agree with RevBayes `dnCDBDP` in rung 3 (§10).
- **Characters.** p(data | 𝒯, paths, Θ) is ordinary pruning with ℓ_e from §7.

**Labelling convention.** (3.4) is the density of an *oriented* tree with labelled tips, with factor λ (not 2λ) per observed
bifurcation and 2λE inside branches. Other conventions differ by parameter-free combinatorial constants, such as 2^{|𝒱|}
or 1/n!. These do not affect inference but must be matched in rung 3: compare log-likelihood *differences* across parameter values,
or fix the constant once.

---

## 4. Tip, node and root terms (spec item d)

| Element | Factor | Notes |
|---|---|---|
| Extant tip | ρ (= D_k(0) for all k) | Its bin at τ = 0 is unobserved; the path still runs to τ = 0 and enters ℓ_e |
| Fossil tip at τ_f | ψ_{k}[r_rm + (1−r_rm)E_k(τ_f)] | k = bin at τ_f. With r_rm = 0: ψ_kE_k(τ_f). Stratigraphic-range uncertainty: integrate τ_f over its range (MCMC move) |
| Sampled ancestor at τ_a | ψ_k(1−r_rm) | Zero-length tip; ℓ = 0 on the tip edge, so characters are observed at the node |
| Bifurcation | λ_{k(v)} | Both daughters start in k(v) |
| Root / origin start | π(k) | Options in §5.4 |

**Origin vs root.**
1. *Origin at age τ_o.* One stem branch from τ_o to the root, carrying its own path, plus a bifurcation factor λ_{k(root)} at the root.
2. *Root at age τ_r.* No stem. The root bin k_r ~ π is shared by both daughter branches. The root factor λ_{k_r} is
   convention-dependent: it is usually omitted when conditioning on the root (Stadler 2010) and cancels against the conditioning in some
   formulations. Fix one convention and state it in code.

---

## 5. Equilibrium tempo distribution and root prior (spec item e)

### 5.1 Discrete replicator–mutation equation

Let n_k(t) be the expected number of extant species in bin k (B&M's N^y). Forward in time each species in k jumps at rates Q_kj
and grows at net rate r_k := λ_k − μ_k (fossil sampling does not change counts if r_rm = 0; otherwise use r_k − r_rm ψ_k). So

$$\dot{\mathbf n} = (Q^\top + R)\mathbf n,\qquad R = \mathrm{diag}(r_k). \qquad (5.1)$$

Since Qᵀ = Aᵀdiag(e^{αx}), (Qᵀ**n**)_k is the discrete forward operator applied to e^{αx}n. This is the semi-discretisation of B&M A.23
(with α = β = 1, r_k = re^{x_k}). The normalised **p** = **n**/**1**ᵀ**n** obeys

$$\dot{\mathbf p} = (Q^\top+R)\mathbf p - (\mathbf 1^\top R\mathbf p)\,\mathbf p, \qquad (5.2)$$

because **1**ᵀQᵀ = (Q**1**)ᵀ = 0. This is the discrete A.24: the term **1**ᵀR**p** = Σ_k r_kp_k is r⟨e^y⟩, with **no** minus sign (§11).

**Proposition (equilibrium = Perron vector).** Qᵀ + R is Metzler (off-diagonals ≥ 0) and irreducible, because A is tridiagonal with
positive off-diagonals. By Perron–Frobenius it has a simple real dominant eigenvalue Λ with a strictly positive eigenvector **p**\*.
The equilibria of (5.2) are exactly the eigenvectors with **1**ᵀR**p** equal to their eigenvalue, and only **p**\* is non-negative.
Moreover

$$\Lambda = \mathbf 1^\top(Q^\top+R)\mathbf p^* = \sum_k r_k p^*_k = \langle r e^{\beta x}\rangle_{p^*}, \qquad (5.3)$$

and Λ is the long-run exponential growth rate of mean diversity, the common slope of B&M Fig 2A. Proof: premultiply
(Qᵀ+R)**p**\* = Λ**p**\* by **1**ᵀ and use **1**ᵀQᵀ = 0 and **1**ᵀ**p**\* = 1. ∎ Global convergence of (5.2) from any non-negative start
follows from the linear flow (5.1) and the spectral gap. This is B&M's "regardless of starting tempo" (Fig 3A).

**Scratch** (core parameters, central, K = 161, z ∈ [−8, 8]): Λ = 0.014332 /Myr. This equals ⟨re^x⟩_{p*} and the log-slope of
N_0(t) over t ∈ [400, 500] (0.0143323). The baseline r = 0.01 is *not* the clade growth rate: tempo selection raises it by 43%.

### 5.2 r = 0 (B&M A.26)

With R = 0, Qᵀ**p** = 0 ⇔ Aᵀ(diag(e^{αx})**p**) = 0 ⇔ diag(e^{αx})**p** ∝ π^A ⇔ p_k ∝ π^A_k e^{−αx_k}.

In the continuum, π^A is N(0, 1) in z, i.e. N(0, s²) in x. So p(y) ∝ e^{−αy}e^{−y²/(2s²)} = N(−αs², s²), with

$$p(y) = \frac{e^{-\alpha^2 s^2/2}}{\sqrt{2\pi s^2}}\,e^{-\alpha y}e^{-y^2/(2s^2)},\qquad \langle y\rangle = -\alpha s^2.$$

At α = 1 this gives mean **−s²**, normaliser e^{−s²/2}. B&M print exp(−s/2) and ⟨y⟩ = −1. Both are correct only for s = 1 (§11).
Interpretation: a single lineage spends time ∝ 1/tempo in each state, so the lineage-time-stationary law is tilted by e^{−x}.

Under SqRA, π^A is the exact grid-Gaussian, so the discrete mean is −s² up to exponentially small grid error. Under upwind it is
−s²(1 + O(h)) (§1.4).

### 5.3 r = θ (B&M A.28) and the r = θ mean dynamics

Multiply (5.1) at equilibrium by x_k and sum. Use **x**ᵀQᵀ**p** = (Q**x**)ᵀ**p**, and the exact interior identity (A**z**)_k = −θz_k
for the central and upwind schemes, so (Q**x**)_k = −θx_ke^{αx_k}. Then, with α = β = 1,

$$(r-\theta)\langle x e^{x}\rangle - r\langle x\rangle\langle e^{x}\rangle = \text{(boundary terms, exponentially small for a wide range)}. \qquad (5.4)$$

The continuum version, from A.27 by parts, is (r−θ)⟨ye^y⟩ − **r**⟨y⟩⟨e^y⟩ = 0. B&M's printed A.28 omits the factor r on the
second term (§11). With the correct form:
- at r = θ, ⟨y⟩ = 0, as B&M state;
- at r = 0, it gives ⟨ye^y⟩ = 0, consistent with N(−s², s²).

The printed form at r = 0 would force ⟨y⟩ = 0, contradicting A.26.

**New exact identity (useful for Fig 3A).** Along the transient (5.2), with α = β = 1 and central or upwind discretisation,
d⟨x⟩/dt = (Q**x**)ᵀ**p** + Σ_k x_k r_k p_k − ⟨x⟩**1**ᵀR**p** = (r−θ)⟨xe^x⟩ − r⟨x⟩⟨e^x⟩ (+ boundary terms). At r = θ:

$$\frac{d\langle x\rangle}{dt} = -r\langle e^x\rangle_t\,\langle x\rangle \quad\Rightarrow\quad \langle x\rangle_t = \langle x\rangle_0\exp\Big(-r\!\int_0^t\langle e^x\rangle_u\,du\Big).$$

So a clade started at x_0 = 0 has mean log-tempo **exactly 0 for all t**: the flat black line in B&M Fig 3A. Clades started at
x_0 = ±2 relax monotonically, faster when ⟨e^x⟩ is large, which is the asymmetry B&M describe. Scratch: from x_0 = 0 the mean is
0.000 at t = 10, 50, 100, 200 and 500. From x_0 = +2 it is 1.181, 0.422, 0.174, 0.038, 0.001; from x_0 = −2 it is −1.973, −1.854,
−1.676, −1.134, −0.037.

### 5.4 Root/origin prior options

1. **Equilibrium of contemporaneous species, p\* (default).** The start lineage is "a random species in an old clade at
   equilibrium", the assumption behind B&M Fig 4 ("the clade originates with an ancestor drawn from the equilibrium").
   Compute it as the Perron vector of Qᵀ + R, with R built from (β, λ, μ) at the current parameters. It changes with every
   parameter move but is cheap (sparse tridiagonal eigenproblem; shift-invert or a few power iterations warm-started).
2. **Lineage-stationary law** π^Q ∝ π^A e^{−αx} (the r = 0 case). This is the natural prior when β = 0, and there it coincides with option 1.
3. **Free start bin** (uniform over bins, or a parameter x_0). This is needed for B&M-style "clade starting at x = 0"
   reproduction, and is the most robust to misspecified dynamics. It is weakly identified.

With survival conditioning, option 1 combined with marginal conditioning (§6) is the generative model "draw x_0 ~ p\*, simulate,
redraw x_0 if extinct". Joint conditioning with a fixed x_0 is B&M's A.11/A.17 setting.

---

## 6. Survival conditioning (spec item d, continued)

Let S_k(τ) = 1 − E_k(τ) (any sampled descendant) and S̃_k = 1 − Ẽ_k (a sampled *extant* descendant).

| Setting | Divide (3.4) by | Remarks |
|---|---|---|
| Origin, at least one sampled descendant | Σ_k π_k S_k(τ_o) (marginal) or S_{k_o}(τ_o) (joint) | Standard FBD origin conditioning |
| Root, both daughters with sampled descendants | Σ_k π_k S_k(τ_r)² (marginal; = BDS eq 6) or S_{k_r}(τ_r)² (joint) | |
| Crown (root), both daughters with **extant** sampled descendants | Σ_k π_k S̃_k(τ_r)² or S̃_{k_r}(τ_r)² | The right choice when the root is the crown node of an extant clade with fossils |
| None | 1 | Only when τ_o carries an informative prior |

- **Marginal vs joint.** Marginal reweights the root-bin posterior by S_k²/Σπ S², while joint does not. They define different
  generative models (BDS; Capobianco & Höhna 2025). B&M condition on survival given x (joint). Recommendation: joint when the start bin
  is a parameter or B&M reproduction is the goal; marginal with option 5.4(1) otherwise. Open issue 3.
- **Why conditioning matters.** Survival is strongly tempo-dependent, and in the direction opposite to naive expectation.
  Scratch values of S_k(T) for x_0 = −2 / 0 / +2:

  | T | x_0 = −2 | x_0 = 0 | x_0 = +2 |
  |---|---|---|---|
  | 10 Myr | 0.594 | 0.177 | 0.055 |
  | 100 Myr | 0.127 | 0.037 | 0.025 |
  | 500 Myr | 0.033 | 0.0215 | 0.0203 |

  All approach 1 − μ/λ = 0.0196. A near-critical process at high tempo reaches its asymptotic extinction risk sooner, while at low
  tempo less has happened. Survival conditioning therefore up-weights **low**-tempo starts, most strongly for short horizons. B&M's
  "push of the past" and the high-tempo ancestry of large clades act through the *size-biased* terms (N_x, D with many sampled
  tips), not through S. The two effects pull in opposite directions, which makes the marginal/joint choice consequential.

---

## 7. Characters and the (β, γ) family (spec item f)

### 7.1 Branch lengths

Given the path on e with dwell d_{e,k} = total time in bin k,

$$\ell_e = c\sum_k e^{\gamma x_k}\,d_{e,k} = c\int_e e^{\gamma x(t)}dt\quad\text{(exact for piecewise-constant paths).}$$

Characters (Mk, Mk′, MkN, JC(k), F81, DNA models) are conditionally independent given {ℓ_e} and are scored by pruning. For a
reversible model on a rooted tree the two root edges may be merged (memo §6.1), so unrooted MkPrime kernels apply unchanged. With
γ = α = 1, ℓ_e = c·w_e where w_e = ∫_e e^x dt is B&M's "myrs-equivalent" molecular change (B&M (10)).

### 7.2 The nested family

The self-rate exponent α is held fixed across the family (default α = 1; §7.3).

| (β, γ) | Diversification | Character clock | Exact reduction |
|---|---|---|---|
| (0, 0) | λ, μ constant | strict, ℓ_e = c·t_e | **Constant-rate FBD + strict clock, exactly** (proof below). Tempo is a free-floating nuisance with no likelihood contribution; θ and s are prior-only and should be removed or fixed |
| (0, 1) | constant FBD | ℓ_e = c∫e^{x}dt with x the (α-)OU tempo, start bin ~ π^Q | **Autocorrelated relaxed clock with log-rate following a mean-reverting process** on a constant-FBD tree. α = 1: rate-modulated OU (log-rate stationary N(−s², s²), mean rate e^{−s²/2}). α = 0: Aris-Brosou & Yang-type OU log-rate in real time. The tree likelihood is independent of the paths |
| (1, 0) | λe^x, μe^x | strict, ℓ_e = c·t_e | **Tempo-dependent hidden-state FBD-SSE** (BDS/QuaSSE-type) with a strict clock. Characters carry no tempo information beyond the tree; tempo is informed only by tree shape and branching times |
| (1, 1) | λe^x, μe^x | ℓ_e = c∫e^x dt | **CET** (B&M) + FBD sampling |

**Proof for (0, 0).** With β = 0 and δ = 0, g_k(τ) = g(τ) is the same for all k, because E_k = E for all k (§2). Then
Q + diag(**g**) = Q + g(τ)I, so 𝓜 = e^{∫g}·e^{QΔ}. Every tip and node factor is k-free, and Σ_j(e^{QΔ})_kj = 1. Summing (3.4) over
paths therefore gives the constant-rate FBD density with branch factors e^{∫g}, i.e. Stadler's q-function. The path marginal is the
CTMC prior, which integrates to 1. With γ = 0, ℓ_e = c·t_e. ∎ In (0, 1) the same argument applies to the tree factor. The path
then enters only through ℓ_e, so its conditional prior given the tree is the plain Q-CTMC started from π^Q.

**Generalisations (not in the 2×2 core):** a free γ ∈ [0, 1] (partial coupling; confounded with s, §9); δ = 1 (fossilisation tracks
tempo); β applied separately to λ and μ.

### 7.3 The self-rate exponent α

B&M motivate the e^x self-rate by tempo being genetically encoded and mutating at the molecular rate (Appendix "schematic"). Logically
that ties α to the molecular coupling γ: if characters, including the tempo-encoding loci, evolve at a constant rate (γ = 0), the
tempo-encoding genome does too. Two defensible conventions:
1. **α ≡ 1 in every cell.** Only the couplings β and γ are toggled, and the tempo dynamics are identical across cells, so Bayes factors
   isolate the couplings. (0, 1) is then "CET clock without CET diversification". *Recommended default.*
2. **α ≡ γ.** This is more faithful to the genetic-encoding argument. (1, 0) then has plain OU tempo in real time and (0, 0) is unchanged.

Open issue 1.

---

## 8. Path sampling (inference-facing consequences of §3)

### 8.1 Exact FFBS of the tree factor

Fix a segment grid on each branch (BDS default: L = 5000 segments per root height). In the postorder pass store **D** at segment
ends. In the preorder pass:
1. sample the start bin ∝ π_k D_k at the root or origin;
2. along each segment [b_i, b_{i+1}] (ages decreasing), sample the end bin j ∝ 𝓜_{kj}(b_i, b_{i+1}) D_j(b_{i+1}).
   𝓜 over the segment is computed by solving (3.2) forward with a one-hot start, as in BDS eq 8, or by expm with **g** frozen.

By Proposition 3.3 this is an exact draw of the segment-endpoint bins from p(paths | 𝒯, Θ) without characters. The same recursion
with **D** including character partials is unavailable (§0 of memo: the joint state space is exponential). Hence: propose with the
tree-only FFBS, and accept by MH on the character likelihood ratio. Because the proposal is the exact tree-conditional, the MH ratio
reduces to the **character likelihood ratio alone**.

### 8.2 Within-segment bridges (exact dwell times)

ℓ_e needs dwell times, so painting segments constant (Freyman & Höhna; BDS) is not enough. Given endpoint bins (k, j) on a segment
of length Δ, sample the path from the tilted endpoint-conditioned law ∝ (3.3). With **g** frozen at its segment midpoint value, this
is an endpoint-conditioned CTMC with sub-stochastic generator Q + diag(**g**):
- *Uniformisation.* Rate Ω ≥ max_k(q_k − g_k); cost ∝ ΩΔ.
- *Direct sampling.* Eigendecomposition of the symmetrised generator; stable by reversibility.

Hobolth & Stone (2009): neither dominates. Uniformisation suits short segments and low-tempo bins, direct sampling long segments
and high-tempo bins. The freeze-**g** error is O(Δ·max|ġ|). Either include it in the convergence study, or remove it with an MH
correction using the exact (3.3) weight.

### 8.3 Local proposals

- Node-bin Gibbs given the three incident branches: weights π-or-λ_k × Π_{3 branches}(segment bridge densities). Follow it with bridge
  redraws (Quintero-style triad move).
- Partial-branch FFBS on a random sub-interval (cheaper; better acceptance when characters are informative).
- Parameter moves keep paths in index space. When θ or s change, 𝓜, E and π\* are recomputed, and dwell times are unchanged. When
  β, λ, μ or ψ change, the path density (3.3) and E change.

---

## 9. Identifiability (spec item g)

### 9.1 The OU-mean ridge

Allow an OU mean m: dx = −θ(x − m)e^{αx}dt + √(2θs²e^{αx})dW. For any a, set x′ = x + a and m′ = m + a. Then

dx′ = −θe^{−αa}(x′ − m′)e^{αx′}dt + √(2θe^{−αa}s²e^{αx′})dW,

λe^{βx} = (λe^{−βa})e^{βx′}, ce^{γx} = (ce^{−γa})e^{γx′}, ψe^{δx} = (ψe^{−δa})e^{δx′}. So

$$(m,\theta,\lambda,\mu,\psi,c) \sim (m+a,\ \theta e^{-\alpha a},\ \lambda e^{-\beta a},\ \mu e^{-\beta a},\ \psi e^{-\delta a},\ c\,e^{-\gamma a})$$

defines the same law of trees, sampling and characters: an exact one-dimensional ridge. **Fix m = 0** (B&M). Then λ means "speciation at
τ = 1", and θ, c are likewise rates at unit tempo.

With m fixed at 0 the ridge is no longer exact, but it remains a weak direction: the characters fix c·∫e^{γx}, the tree fixes
λe^{βx}, and only the OU pull towards z = 0 and the start prior pin the tempo level. The `shift` move of `TempoMCMC()` (#19)
follows it on the grid: every bin k → k ± 1 (x → x + d, d = ±hs) with θ, λ, μ, ψ, c multiplied by e^{−αd}, e^{−βd}, e^{−βd},
e^{−δd}, e^{−γd}. This is the reparameterisation above with a = d applied to the free parameters, so the character
lengths are unchanged exactly; the acceptance ratio carries the change in the augmented density (generator drift, start
prior and conditioning, which depend on z) and the Jacobian Π p′/p.

### 9.2 Time scale

dCET is time-homogeneous. Scaling all ages by κ and all rates (λ, μ, ψ, θ, c) by 1/κ leaves the law of 𝒯 in rescaled time and of the
characters unchanged; s, ρ, β, γ are invariant. Without tip ages or calibrations only rate×time products are identified, as in any
clock model. Fossil ages in Myr fix κ.

### 9.3 Tempo-time reformulation, exact in dCET

Let W(t) = ∫_0^t e^{αx(u)}du along a lineage. **Lemma (random time change of a CTMC).** If X has generator diag(v)A with v > 0, then
Y(w) := X(W^{−1}(w)) is a CTMC with generator A. Proof: a sojourn in k of real duration T ~ Exp(v_k a_k), with a_k = −A_kk, has
w-length v_kT ~ Exp(a_k), and the jump probabilities A_kj/a_k are unchanged. Event processes with real-time intensity η_k have w-time
intensity η_k/v_k. ∎

Apply it with v_k = e^{αx_k}. In w-time each lineage carries an OU-CTMC A (independent of s apart from the grid), speciation
λe^{(β−α)x}, extinction μe^{(β−α)x}, fossil sampling ψe^{(δ−α)x}, and character rate ce^{(γ−α)x}. **At α = β = γ = 1** (CET) the
branching and the character clock are homogeneous in w. The tempo path then enters only through:
1. fossil sampling if δ = 0;
2. the real-time stopping rule: a lineage is stopped (and sampled with probability ρ) when its real-time clock
   ∫e^{−x}dw reaches the present.

**Corollary 1 (B&M eq 12, exact for every K and every scheme).** In the complete process (no reconstruction) with α = β = γ = 1, the
w-length of a branch (birth to its next speciation or extinction) is Exp(λ+μ), whatever the starting bin.

*Proof 1.* In w-time the branch ends at the first event of a Poisson process of constant rate λ+μ, independent of Y. ∎

*Proof 2 (Laplace transform).* φ_k(u) = E[e^{−uw} | k] solves 0 = Qφ − (λ+μ)diag(e^x)φ − u·diag(e^x)φ + (λ+μ)e^x. Since Q = diag(e^x)A,
this becomes [(λ+μ+u)I − A]φ = (λ+μ)**1**. A**1** = 0 then gives φ = (λ+μ)/(λ+μ+u)·**1**. ∎

Scratch: mean w = [(λ+μ)diag(e^x) − Q]⁻¹e^x = 0.990099 = 1/(λ+μ) in every bin, to printed precision. Eq 12 is exact in dCET, so it
tests the **simulator and path bookkeeping**, not the discretisation. If ψe^{δx} has δ = 1 and fossil sampling ends a branch, the law
is Exp(λ+μ+ψ). If β ≠ γ or α ≠ γ the result fails, which gives an equally sharp negative control.

**Corollary 2.** The tempo invariance of the long-run extinction probability μ/λ follows. Survival to w = ∞ in a homogeneous BD is
1 − μ/λ, and real time t → ∞ implies w → ∞ almost surely for the ergodic tempo chain.

**Implications for inference (1,1).**
- *Characters see the w-tree.* Its branching is homogeneous BD in w, but its tips are not ultrametric in w. Extant-only data
  therefore constrain tempo only through how the w-tree's non-ultrametricity and branch-length/topology pattern can be produced by
  OU paths that equalise real-time depths. This is the relaxed-clock signal, and an autocorrelated clock (0, 1) can mimic much of it.
  The discriminating feature is that in (1,1) **branching density per unit w is constant** (fast lineages branch more per Myr but not
  per unit change), whereas in (0,1) branching density per Myr is constant. Discrimination should rest on that.
- *Eq 12 is a complete-tree statement.* In the *reconstructed* tree, hidden speciations are removed and branches are survival-conditioned,
  so reconstructed branch w-lengths are **not** Exp(λ+μ). They are longer and tempo-dependent through E. Do not use eq 12 as a check
  on reconstructed trees.
- *Tip ages are direct evidence.* A fossil tip pins a real-time depth while characters give its w-depth, so tip dates test the
  w ↔ t mapping directly (memo recommendation 4).

### 9.4 Other weak directions to expect

- **(s, c) under γ = 1.** Mean rate ∝ c·E[e^x], so c and s trade off through the e^{s²/2}-type factor. θ controls autocorrelation time.
- **(θ, s, λ, μ) under β = 1** with extant-only trees: Quintero et al. found extinction diffusion unrecoverable on trees of ≤ 200 tips.
  Expect similar for μ/λ separately. Fossils help (Truman et al.; MTFBD), but only if fossils are placed by characters.
- **Free γ vs s.** Both change rate variance: γs is the sd of the log character rate, so they are confounded unless β pins s through tree
  shape.

---

## 10. Checkable identities and numerical tests by validation rung (spec item h)

Paper values are quoted from B&M. Scratch values are from [dcet_scratch_checks.py]({{artifact:ec1e79d0-9b19-419f-9654-43660a095567}}): central scheme unless stated, core parameters
λ = 0.51, μ = 0.5, θ = 0.01, s = 1, α = β = γ = 1, ρ = 1, ψ = 0. They are **provisional** and must be reproduced by the
reference implementation.

### Rung 1 — numerics vs B&M and internal identities

"B&M code" values are from a verbatim re-run of the Dryad functions ([bm_reproduce.R]({{artifact:65cc1388-dd0f-4aa9-9097-d6d0f3d42abe}}) →
[bm_reproduce_output.txt]({{artifact:27000047-2dbe-426e-b8b3-9eccbaee2542}})), cross-checked by an independent replica ([bm_scheme_replica.py]({{artifact:5d0a72a5-a98b-4ca4-b955-69183cc78aab}}) →
[bm_scheme_replica_output.txt]({{artifact:cb32d2e0-313a-42a4-875b-5f344416dd5d}})).

| # | Quantity | Paper | dCET scratch | Tolerance / note |
|---|---|---|---|---|
| 1.1 | E-solver at β = 0 vs Stadler p₀ (FBD, ρ < 1, ψ > 0) | — | exact by §2 | 1e−10 relative |
| 1.2 | E_k(τ) → μ/λ = 0.980392 as τ → ∞, all k | μ/λ, tempo-invariant | fixed point | 1e−8 |
| 1.3 | S_0(500) (survival from x_0 = 0) | not reported | 0.02154 | K-converged |
| 1.4 | Fig 1A mean clade size, surviving, T = 500, x_0 = 0 | "c. 60 000" | **52 644** converged (K = 401, z ∈ ±10); **B&M code: 57 836** (verbatim re-run) | Paper value = code value rounded; code is +9.9% vs converged |
| 1.5 | Fig 1A experienced clade size = 1 + M₂/N, with M₂ = E[n(n−1)] from ∂_tM₂ = (Q+R)M₂ + 2λ diag(e^x)N², M₂(0) = 0 | "c. 400 000" | **352 938** converged; **B&M code: 425 699** | Ratio: paper "c. 8", code 7.36, converged 6.70; code is +20.6% |
| 1.6 | Fig 1A P(n < 1000 \| n > 0) ≈ 1/3, P(n > 50 000 \| n > 0) ≈ 1/4 | as stated | **not computed** | Needs PGF solve on \|z\| = 1 + FFT (A.2); rung-1 implementation task |
| 1.7 | Fig 1B P(larger sister ≥ 20× smaller) ≈ 0.5 | ≈ 0.5 | not computed | Needs P_n |
| 1.8 | Fig 2A N_x(500), x_0 = −2/0/2 | figure only | converged 53.5 / 1134.8 / 3971; B&M code 52.4 / 1236.4 / 4373.3 | N_x(100): converged 1.19 / 3.70 / 11.47; code 1.19 / 3.75 / 11.76 |
| 1.9 | Fig 2A common asymptotic slope = Λ (Perron) = ⟨re^x⟩_{p*} | "fixed rate" | **0.014332 /Myr** converged; B&M code 0.014511 (log-slope 400–500) | Eigenvalue vs fitted log-slope agree to 1e−6 |
| 1.10 | Fig 3A equilibrium (r = θ): mean, sd | mean 0 | 0 (exact, central/upwind), sd 0.833 | Mean residual ≤ 1e−12 |
| 1.11 | Fig 3A transients; x_0 = 0 mean ≡ 0 | flat black line | 0 at all t (exact, §5.3) | From ±2: values in §5.3 |
| 1.12 | Fig 3B equilibrium mean vs r (θ = 0.01) | increases with r | equilibria: r = 0: −1; 0.005: −0.409; 0.01: 0; 0.02: 0.634. **As plotted** (x_0 = 0, t = 500): r = 0: −0.804 converged / −0.777 code; r = 2θ: 0.635 / 0.711 | The r = 0 curve has not reached −1 by 500 Myr |
| 1.13 | A.26 mean at r = 0 | −1 | −s² (= −1 at s = 1); SqRA exact, central O(h²), upwind O(h) | Test also at s = 0.5, 2 (should give −0.25, −4) to catch the A.26 typo |
| 1.14 | A.28 residual (corrected form) | — | ≤ 1e−11 (central/upwind) | Printed form fails at r = 0 |
| 1.15 | Fig 4A mean ancestral log-tempo (eq 7, start law Px = M[501,,lx] normalised), look-back 100 / 200 / 500 Myr | "plateau ≈ 0.6 at c. 100 Ma" | converged **0.443 / 0.528 / 0.552**; B&M code **0.462 / 0.548 / 0.572** | Same definition; the paper text is a loose reading of the code's curve (§11 item 8) |
| 1.16 | Fig 5A mean branch duration, x_0 = −2/0/2 | "lower tempo → longer" | 7.244 / 0.9998 / 0.138 Myr | Solve [(λ+μ)diag(e^x) − Q]τ = **1**. The B&M code integrates only to 20 Myr (x_0 = −2: 94% of mass); Fig 5 curves are rescaled to unit peak and Fig 5B is analytic |
| 1.17 | Fig 5B branch w ~ Exp(λ+μ) | eq 12 | mean 0.990099 in all bins (exact) | Machine precision |
| 1.18 | K-convergence of 1.3–1.16 at h ∈ {0.4, 0.2, 0.1, 0.05} | — | central O(h²), upwind O(h) | Report observed order |

### Rung 2 — forward simulator vs numerics

1. Branch w-lengths from complete simulated trees: KS test against Exp(λ+μ) for starts in several bins, n ≥ 10⁴ branches (exact for any K).
   Negative control: β ≠ γ must fail.
2. Branch real durations: means 7.244 / 1.000 / 0.138 for x_0 = −2/0/2.
3. Survival at T vs S_k(T), within binomial CI. Mean diversity at moderate T, e.g. N_0(100) = 3.70, avoiding T = 500 where the variance
   is extreme.
4. Tempo distribution among extant species at T vs (5.1) normalised.
5. Continuous-CET simulator (exact OU in w; real-time integral for stopping) vs dCET: differences shrink at the rate of §1.2.
6. Sampling layer: fossil count ~ Poisson(ψ × total lineage time) when δ = 0, with E-based reconstructed-tree statistics vs the ODEs.

### Rung 3 — tree likelihood vs RevBayes `dnCDBDP`

- Same Q (passed as a custom anagenetic rate matrix), λ_k, μ_k, root frequencies and conditioning. Log-likelihood differences across
  ≥ 20 parameter points must agree to 1e−8. The constant offset is the labelling convention, §3.4.
- Check whether `dnCDBDP` accepts ψ. If not, test ψ via the β = 0 collapse to the analytic FBD, and the extant-only part via dnCDBDP.
- β = 0 collapse to the constant-rate (F)BD density for arbitrary Q (§7.2 proof).

### Rung 4 — augmentation vs marginal

- **4a.** Fixed tree, no characters. DA posterior of (λ, μ, θ, s), and per-branch posterior mean dwell in each bin, must equal
  (i) MCMC on the marginal likelihood of rung 3, and (ii) FFBS stochastic-mapping summaries.
- **4b.** One binary (or JC(k)) character. The joint (bin, state) process is a covarion-type CTMC on K×k states with generator
  Q ⊗ I + diag(ce^{γx}) ⊗ Q_char. The exact marginal likelihood is a MuSSE-type pruning on that product space. The DA posterior must match.
  This is the first test that includes characters.
- **4c.** Path-sampler unit test (Lepage et al. 2006 mgf). For endpoints (i, j) over Δ, with Q̂ = Q + diag(**g**) frozen:

  E[e^{−uℓ} | i → j] = [e^{(Q̂ − u·diag(ce^{γx}))Δ}]_{ij} / [e^{Q̂Δ}]_{ij}.

  Compare the empirical Laplace transform of sampled ℓ at several u.
- **4d.** Sum-over-paths identity (Proposition 3.3) by Monte Carlo on small K (done in scratch for K = 4).

### Rung 5 — prior-only MCMC

With the likelihood switched off but conditioning kept: tree statistics (tip/fossil counts, root age if origin-dated, γ-statistic,
imbalance), path statistics (occupancy, jumps per branch), and ℓ_e distribution must match forward simulation under identical
conditioning.

### Rungs 6–7 — SBC and power

SBC: 20–40-tip trees, per (β, γ) cell, rank histograms for all parameters and for selected path functionals. Use parameter priors that
keep clade sizes simulable (memo §2.7: shorter T or n-conditioning). Power: Bayes factors across characters × fossil density × tree
size. Include the (1,1)-vs-(0,1) contrast, which §9.3 predicts to be the hardest.

---

## 11. Issues found in Budd & Mann 2025

Items 1–5 were checked against rendered crops of the typeset PDF ([bm-equation-crops.png]({{artifact:cb01110d-f871-4ec9-9127-86a880c08016}})).

1. **Eq 9 initial condition sign.** Main text: f_x(t = 0) = **−**e^x(λ+μ). A.31 text: f_x(t = 0) = e^x(λ+μ). A density must be
   positive, and f = ∂F/∂t with F(0) = 0 and A.30 gives e^x(λ+μ)(1 − 0). **Main-text minus is a typo.** Eq 11's f_x(w = 0) = λ+μ is
   consistent with the positive form.
2. **Stray minus in ⟨e^y⟩ under A.24.** Printed as ⟨e^y⟩ = −∫e^{y′}p(y′)dy′. With the minus the replicator term
   rp(e^y − ⟨e^y⟩) would not conserve probability. Main-text eq 5 has no minus. **Typo; ⟨e^y⟩ = +∫e^y p.**
3. **A.26 normaliser and mean.** Printed normaliser exp(−s/2)/√(2πs²) and ⟨y⟩ = −1. The correct values are e^{−s²/2}/√(2πs²) and
   ⟨y⟩ = −s² (§5.2). Both printed forms are correct only at s = 1, the paper's value. Harmless for the paper's figures, but wrong if
   re-used with other s.
4. **A.28 missing r.** Printed (r−θ)⟨ye^y⟩ − ⟨y⟩⟨e^y⟩ = 0. Integration by parts of A.27 gives (r−θ)⟨ye^y⟩ − **r**⟨y⟩⟨e^y⟩ = 0. The
   paper's conclusion (⟨y⟩ = 0 at r = θ) holds for either form when r ≠ 0. At r = 0, however, the printed form implies ⟨y⟩ = 0,
   contradicting A.26. The discrete check (§1.4) confirms the corrected form to 1e−11.
5. **A.13 left-hand side.** Printed Σ n_t P(n_t | n_T **=** 0) = [Σ n_tP(n_t) − Σ n_tP(n_t, n_T = 0)]/P(n_T > 0). The left side should
   be conditioned on n_T **>** 0, as the next line E(n_t | n_T > 0) shows. Typo.
6. **A.14/A.15 initial conditions. Resolved: the code solves per t (correct), but `pdeC` has a reaction-term bug.** The code loops over
   the 501 output times. For each t it runs a fresh `vode` solve over [500 − t, 500] from C = H = G(T − t, 0), with H carried as
   C2 = H − μ/λ. This is the correct per-t procedure, so the text is ambiguous but the code is right. However, C2 is advanced with
   λC2(C2 − 1) instead of λC2(C2 − 1 + μ/λ) (the latter is what `pdeU` uses). Effect on the Fig 2B survivor-conditioned mean from
   x₀ = 0 (verbatim re-run vs corrected): −2.0% at t = 100, −0.2% at t = 250, −0.02% at t = 400. Invisible on the figure.
7. **Fig 1A numbers. Confirmed, with the mechanism refined.** The scheme is ReacTran `tran.1D` on 101 finite-volume cells over
   (−10, 10), h = 0.19802. The drift term uses `v = 1`, `AFDW = 1`, i.e. a backward difference fixed in direction. It is upwind for
   x > 0 and downwind for x < 0, so the effective diffusion is θ(s² + xh/2) (+20% at x = 2) and down-rates are negative for
   x < −s²/h ≈ −5.05 (24 cells). The v1 hypothesis was a symmetric upwind scheme (59 717 / 439 940). The actual scheme gives
   **57 836 / 425 699** in a verbatim R re-run of B&M's functions. These are the paper's "c. 60 000" and "c. 400 000", and vode
   rtol = 1e−10 leaves them unchanged. Under refinement the B&M scheme converges at first order to the dCET limit
   (Richardson 52 621 / 352 652 vs central 52 644 / 352 938). The paper's values therefore overstate the PDE solution by 9.9% (mean)
   and 20.6% (experienced). The qualitative conclusions stand.
8. **Fig 4A. Refuted as a definitional discrepancy.** The code computes eq 7 with start law Px = M[501,,lx] normalised: the r = θ
   tempo distribution after 500 Myr starting from the top cell, a proxy for the equilibrium (mean 0.0005, sd 0.831). That matches
   §5.4 option 1. The code's own curve is 0.462 / 0.548 / 0.572 at look-back 100 / 200 / 500 Myr, against converged 0.443 / 0.528 /
   0.552. The paper's "plateauing at x ≃ 0.6 at c. 100 Ma" is a loose verbal description; discretisation adds about +0.02.
9. **Notation in A.2.** The Cauchy/Fourier inversion prefactor is 1/(2π). The extracted text shows "1/2∫", probably a rendering
   artefact (not checked on the page image).
10. **Sign conventions otherwise consistent.** Eq 1, A.4 and eq 5 (forward operator +θ∂(xe^xp)/∂x + s²θ∂²(e^xp)/∂x²) are mutually
    consistent with the backward operator e^x(−θx∂_x + θs²∂²_x) in eqs 2, 4, 9 and 11.
11. **Code: Fig 1A normalisation mismatch** (code inspection; figure not re-rendered). The plotted curve is the *unconditioned* P_n
    (it sums to S ≈ 0.021 over n ≥ 1), while the dashed geometric comparator `dgeom(n, 1/meanclade)` is normalised. Relative to the
    comparator the CET curve is drawn a factor ≈ 1/S ≈ 47 too low on the log axis.
12. **Code: `evalmx` misplaced parenthesis.** `which.min(abs(x) − x0[1])` should be `which.min(abs(x − x0[1]))`. It is harmless as
    called, because x0[1] = 0 in every call.
13. **Code: hard-coded `102:202`** in `margexmean`. Correct only for lx = 101.
14. **Code: Fig 5.** Eq 9 is solved from the positive initial condition e^x(λ+μ), which confirms item 1 as a typo. Fig 5B plots
    `dexp(λ+μ)` analytically rather than solving eq 11, and all Fig 5 curves are rescaled to a common peak height.
15. **Code: PGF inversion** (`evalQ`, not re-run). MM = 2²¹ points on the unit circle (no damping radius), zvode default rtol 10⁻⁶,
    `abs(fft)`. Aliasing of mass at n ≥ 2²¹ and rtol-limited tail accuracy are not quantified, so P(n < 1000) ≃ 1/3 and
    P(n > 50 000) ≃ 1/4 remain unverified.

---

## 12. Parameter summary

| Symbol | Meaning | Status |
|---|---|---|
| λ, μ | speciation/extinction at τ = 1 | free |
| θ | mean-reversion rate (per unit tempo-time) | free |
| s | stationary sd of log-tempo (OU in tempo-time) | free |
| c | character rate at τ = 1 | free (one per partition) |
| ψ | fossil sampling rate (at τ = 1 if δ = 1) | free |
| ρ | extant sampling probability | fixed (external); random vs diversified sampling as sensitivity |
| r_rm | removal after sampling | fixed 0 |
| α | tempo self-rate exponent | fixed 1 (open issue 1) |
| β, γ | couplings | {0, 1} per model cell; γ ∈ [0, 1] stretch |
| δ | fossil-sampling coupling | {0, 1} sensitivity axis |
| m | OU mean | fixed 0 (§9.1) |
| K, z_min, z_max, scheme | discretisation | not parameters; convergence-checked (BDS protocol + §1.3 bound) |
| π | start-bin prior | §5.4 option (open issue 2) |


---

## 13. Open issues needing a decision

1. **Self-rate exponent α.** α ≡ 1 in every (β, γ) cell (recommended: isolates the couplings) or α ≡ γ (closer to B&M's genetic-encoding argument). §7.3.
2. **Start-bin prior.** Equilibrium Perron vector p\* (recommended default), lineage-stationary π^Q, or a free start bin. §5.4.
3. **Survival conditioning.** Marginal or joint, and origin vs root vs crown (extant-descendant) conditioning. The choice matters
   because survival up-weights low-tempo starts (§6). Recommended: crown/extant conditioning for clades with an extant crown,
   joint when the start bin is a parameter.
4. **Fossil-sampling coupling δ.** δ = 0 (constant ψ) as default with δ = 1 as a sensitivity axis, or the reverse.
5. **Discretisation scheme.** Central (recommended, h ≤ 2/z_max) vs SqRA vs upwind, and the inference z-range (default ±5 in z,
   asymmetric allowed).
6. **Reference values for rung 1.** Proposed rule: implement B&M's exact operator as a test fixture; it must reproduce the code column of the
   numerics notes to solver tolerance. The converged dCET column is the validation reference, and the paper text is qualitative only. Still
   unverified: P(n < 1000), P(n > 50 000), Fig 1B and Fig 4B, all of which need the PGF/FFT inversion (§11 item 15). Decide whether rung 1
   should include a full `evalQ`-scale inversion (a 2²¹-point FFT; heavy; Hamilton) or a damped-contour inversion for n ≤ 10⁴ only.
7. **RevBayes oracle.** Whether `dnCDBDP` supports fossil sampling ψ, and where the local RevBayes checkout lives (memo decision 4).
