# Tempo location (#116): V0 ridge profile and rung reruns under the mean-rate model

The spec is `dev/spec/tempo-location.md`; the implementation plan is `dev/plans/tempo-mean-rate.md`.

Validation used the package built at 2e7836d, with the rung scripts at b8ab46a. The runs were on Hamilton, in
`/nobackup/pjjg18/evotempo-meanrate/`, at the default grid (`EVOTEMPO_GRID=default`: K = 26, z in [-5, 5], SqRA) with
`location = "mean"`. V0 ran locally.

## V0: the Hymenoptera (1, 0) ridge profile (`v0/`)

**Script and inputs.** `data-raw/tempo-location-v0.R` evaluates `TreeLikelihood()` at six points along the ridge, on the
last tree of the pilot chain (`pilot-g`). The points are geometric interpolations t = 0, 0.2, ..., 1 between the chain's
first and last parameters. The rates are mapped to lineage-mean rates by multiplying by exp(-s^2/2).

**Decomposition.** `v0/decomp.R` and `v0/decomp2.R` split each point's error into its window and grid-spacing parts.
Their output is in `v0/decomposition-t*.txt`. Both scripts read the pilot chain from a local path; it is not in the repo.

| t | s | secs | logLik | logLik, reference (K = 51, z in [-7, 7]) | difference | top-3 mass, root prior | top-3 mass, e^x-tilted law |
|---|---|---|---|---|---|---|---|
| 0.0 | 0.96 | 0.31 | -830.972 | -831.055 | 0.083 | 7.3e-4 | 1.1e-3 |
| 0.2 | 1.16 | 0.85 | -828.056 | -828.150 | 0.095 | 3.4e-4 | 2.1e-3 |
| 0.4 | 1.40 | 6.25 | -825.693 | -825.773 | 0.080 | 2.0e-4 | 4.4e-3 |
| 0.6 | 1.69 | 26.6 | -824.078 | -824.086 | 0.008 | 1.5e-4 | 9.8e-3 |
| 0.8 | 2.04 | 46.3 | -823.188 | did not solve | | 1.2e-4 | 0.024 |
| 1.0 | 2.46 | 116.8 | -823.624 | did not solve | | 1.3e-4 | 0.058 |

Under the old model (`bm`, the current grid) the same points took 0.5, 6.5, 53 and 400 s, and the evaluation failed
beyond them. The mean-rate model equals `bm` on the shifted window to <= 1.4e-5. The full table is in `v0/v0-*.csv`.

**Against the spec's go criteria (§V0):**

- **Time: met.** Every point evaluates, and the points up to s = 2.2 take under 60 s.
- **Agreement with the reference to < 0.05 log units: not met** at s <= 1.4 (0.08-0.095). The reference does not
  solve at s >= 2. The decomposition puts almost all of the gap in grid spacing (K = 26 against K = 51 on the same
  window). On the same points, the old default's spacing error is as large or larger (0.18-0.20). Widening the window
  to [-7, 7] at the same h changes the log-likelihood by <= 0.001 at s <= 1.4.
- **Top-three-bin mass < 1e-3: met for the root prior, not met for the e^x-tilted law.** The tilted law exceeds 1e-3
  at every point, rising to 6% at s = 2.5.

So the mean-rate default removes the failure on the ridge and is at least as accurate as the old default where both solve. It
does not reach the spec's absolute thresholds at large s. Whether to raise K, or to accept the K = 26 spacing error, is
for the maintainer; see the PR.

## Rungs at the default grid (`rungs/`)

| Rung | Check | Result |
|---|---|---|
| 1 | Budd & Mann reproduction | Unchanged: the scripts are pinned to `location = "bm"` |
| 2 | Simulators vs numerics, `EVOTEMPO_LOCATION=mean` (`rung2-mean.rds`) | Agree. The x0 = 2 dCET branch-duration z = -2.6 is the same as in the `bm` run (a known discretisation offset) |
| 3 | beta = 0 collapse (`rung3-collapse.csv`) | <= 4.8e-10 relative |
| 3 | Monte Carlo (`rung3-montecarlo*.csv`) | chi-squared p = 0.32, max abs(z) 2.80 |
| 4a(i) | Augmented vs marginal posterior, 20 000 iterations (`rung4a-quantiles.csv`, `rung4a-summary.csv`) | max abs(z) 2.45 over 25 quantiles; KS p >= 0.085 |
| 4a(ii) | Local moves vs weighted FFBS, 16 chains, iterations 10 001-20 000 (`rung4a-ffbs-pooled16.csv`) | max abs(z) 0.84 (branch means), 1.67 (occupancy, bins with a non-zero reference) |
| 4b | Likelihood check (`rung4b-likelihood-check.csv`) | abs(z) <= 1.34 |
| 4b | One-parameter exact posterior: lambda, mu, theta (8 chains); s; psi (16 chains x 6000) | All abs(z) <= 1.65 (`log-onepar2.txt`, `log-ones.txt`) |
| 4b | One-parameter, `shift` arm on theta (16 chains) | abs(z) <= 1.92 (`log-shift.txt`) |
| 4c | Bridge Laplace transform | max abs(z) 1.95 |
| 4d | Sum over paths | max abs(z) 2.19 |
| 5 | Prior-only MCMC vs simulator (`rung5-*.csv`) | 55 quantiles, max abs(z) 2.00; KS p >= 0.10 |
| 5 | Tree moves at fixed parameters (`rung5-fixed-parameters.csv`) | abs(z) <= 1.91 |
| 4b | Theta identity (`rung4b-theta-identity.R`), 5 thetas | abs(z) <= 0.82 without the character; <= 2.27 with it |
| 4b | Shift move at fixed parameters (`rung4b-shift-fixed.R`) | Root-bin posterior: shift arm abs(z) <= 1.73, reference arm <= 1.77 |
| 5 | Start-age oracle (`rung5-age-oracle.R`), root / origin | max abs(z) 1.16 / 1.51 |
| 5 | Start-age validation (`rung5-age-validation.R`) | root: 2.59 over 80 tests; origin: 2.88 over 80; originFree: 2.02 over 86 |
| 6 | Full SBC, 200 replicates per cell (`sbc/`) | All four cells pass; 52 tests, smallest raw p 0.040, Holm-adjusted p = 1 |

**Notes on the rungs.**

- **4a(i).** A first run at 5000 iterations gave psi q50 z = 5.7 and s KS p = 0.001. At 20 000 iterations both pass, so
  the 5000-iteration result was mixing error. The old (`bm`) run at its own length flagged mu (z = 3.58) and psi
  (KS p = 0.01) in the same way.
- **4a(ii).** The report's own table (`rung4a-ffbs.csv`: 8 chains, seeds 501-508, 20 000 iterations, all after
  `burnin = 100`) shows branch means 1.5-3 sigma above FFBS, all of one sign. A second set of 8 chains (seeds 601-608)
  was run to check this. Splitting the 16 chains in half shows a start-up transient, not a sampler bias.
  - Iterations 1-10 000: branch max abs(z) 3.38, all positive.
  - Iterations 10 001-20 000: branch max abs(z) 0.84, 10 of 16 positive; occupancy max abs(z) 1.67.

  The cause is the short burn-in. The local chains start from paths up to +/-0.45 from the reference, and the
  lag-100 autocorrelation of the branch mean is 0.53, so `burnin = 100` does not remove the transient. The committed
  table uses the second half. `rung4-local-chains.R` keeps its 100-iteration burn-in; a longer burn-in there would
  make the report's own table usable.

  Occupancy bin 25 has an undefined z. The FFBS reference has no draws there, so its SE is 0, while the chains visit
  the bin with frequency 1e-5.
- **4b one-parameter oracle.** `rung4b-one-parameter.R` now rebuilds stored pre-#116 models under `EVOTEMPO_LOCATION`,
  and its grid includes s.
- **Not run.**
  - Rung 3 castor MuSSE: castor is not installed on Hamilton.
  - Rung 3 RevBayes `dnCDBDP` with the new Q.
  - Rung 7, which is optional (V4).
  - V5, the Hymenoptera re-pilot, which belongs to the Hymenoptera pipeline.
- **Theta, shift and start-age scripts.** These use their own fixtures, not the default grid: the rung 4b fixture,
  and K = 16 on [-3, 3] for the start-age scripts. They ran under `location = "mean"` (Hamilton jobs 19572305-11;
  `log-theta-shift-age.txt`).
- **Logs.** `log-*.txt` are the Hamilton job logs. Each records the lib SHA and the source SHA.
