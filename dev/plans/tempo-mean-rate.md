# Plan: the mean-rate tempo location (#116)

The specification is `dev/spec/tempo-location.md` (PR #117), which this implements. Decisions are recorded on #116:
one model for every dataset; penguins switch and rerun everything; trilobites wait for the new model.

## Design decisions

1. **One argument, `location = c("mean", "bm")`.** It lives on the grid object: `TempoGrid(nBin, zRange, s,
   location)`.
   - `"mean"` sets x_k = s z_k − s²/2, and the generator's OU drift centre is ζ = αs in z.
   - `"bm"` sets x_k = s z_k with ζ = 0, the current model.
   - `TempoGenerator()` reads the location from the grid, so a grid and its generator cannot disagree.
   - `DCETModel(location = "mean")` passes it to `.Setup()`.
   - The default is `"mean"` everywhere. `"bm"` is used only by rung 1 and rung 2 code that reproduces Budd & Mann,
     and by tests of B&M identities (A.26, A.28, eq 12).
2. **Scheme default.** `scheme = NULL` means `"sqra"` for `"mean"` and `"central"` for `"bm"`. This applies in
   `TempoGenerator()` and `DCETModel()`. The central and upwind schemes take ζ: drift −θ(z − ζ). The central positivity
   check uses max|z − ζ|.
3. **SqRA with centre.** The potential is U = (z − ζ)²/2, so the stationary law of Q is exactly the grid Gaussian (spec
   L.2).
4. **Old saved models.** A `DCETModel` without a `location` element (checkpoints written before this change) is read
   as `"bm"` with its stored scheme. Resuming an old chain therefore reproduces it, and the pipelines' own build check
   still refuses a new lib.
5. **`SimulateCET()` gains `location`.** The OU in tempo time has mean m = s²/2 under `"mean"` (α = 1 in that
   simulator) and 0 under `"bm"`. `SimulateDCET()` takes Q and x, so it needs no change.
6. **No C++ change.** The C++ engine reads `set$Q`, the per-bin rates, the start prior and x from R. It never builds the
   grid or Q. This keeps the change disjoint from the FNPR proposal and the incremental-likelihood work (#110), which
   touch `src/` and the move code. Shared files touched here: `R/Generator.R`, `R/TreeLikelihood.R` (`DCETModel`,
   `.Setup`), `R/Simulate.R`, and the tests.
7. **The shift move is unchanged.** x moves by hs at fixed s under either location (spec L.6).

## Correctness criteria

1. **Grid identity (spec L.3).** `"mean"` on [a, b] at parameters (θ, λ, μ, ψ, c) equals `"bm"` on
   [a − αs, b − αs] at parameters θ₀ = θe^{αm}, λ₀ = λe^{βm}, μ₀ likewise, ψ₀ = ψe^{δm}, c₀ = ce^{γm}, with
   m = (α − ½)s². Check `TreeLikelihood()` to 1e−10 relative, in both engines, for every scheme, α, β, γ, δ ∈ {0, 1},
   and start at root and origin.
2. **Stationary law.** Under SqRA, `TempoEquilibrium(Q)$p` equals the normalised φ(z_k) to 1e−12, and
   Σ p e^{x} = Σφ e^{sz−s²/2}/Σφ.
3. **The `"bm"` path is unchanged.** Seeded chains and likelihoods under `location = "bm"` match `origin/main`
   exactly.
4. **Engines agree.** R and C++ agree to ≤ 1e−8, and seeded `TempoMCMC()` chains are identical in both engines, at the
   new default.
5. **Coverage.** 100% of new or changed lines.

## Validation (spec L.6)

- **V0.** The ridge-profile numerics check on Hamilton.
- **V2 rungs touched.**
  - Rung 2: the continuous simulator under `"mean"`.
  - Rung 3: the β = 0 collapse and the castor and Monte Carlo comparisons at the new default.
  - Rungs 4 and 5: internal oracles at the new default grid.
- **V3.** Full SBC on Hamilton.
- **V5.** The Hymenoptera (1, 0) re-pilot.
- **Penguin full rerun.** The launch plan is in the PR body; it is not submitted.
