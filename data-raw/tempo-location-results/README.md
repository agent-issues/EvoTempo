# Tempo location (#116): V0 ridge profile and rung reruns under the mean-rate model

The spec is `dev/spec/tempo-location.md`; the implementation plan is `dev/plans/tempo-mean-rate.md`.

**Provenance.** The package was built at 2e7836d, on Hamilton in `/nobackup/pjjg18/evotempo-meanrate/`. V0 ran
locally. The `src=` SHA printed in the job logs is **not reliable** for runs before the review of #123. Those runs used
a source copy into which individual scripts had been copied after it was made:
- `log-ones.txt` names b8ab46a, but its s arm exists only from a46b4bc;
- `log-theta-shift-age.txt` names b8ab46a, but its two rung-4b scripts ran with the location rebuild added in b0b7861.

Their outputs show the `"mean"` model (the deterministic `exact` columns differ from the `"bm"` values). Runs after
the review use a clean `git archive` of the named commit.

**Grids.** The rung 4 and rung 5 scripts that honour `EVOTEMPO_GRID=default` ran at the default grid (K = 26, z in
[-5, 5], SqRA, `location = "mean"`). The others used their own grids, all under `"mean"`:
- rung 6 SBC: K = 16 on [-3, 3] (see below);
- the start-age scripts: K = 16 on [-3, 3];
- the theta-identity and shift-fixed scripts: the rung 4b fixture.

## V0: the Hymenoptera (1, 0) ridge profile (`v0/`)

**Script and inputs.** `data-raw/tempo-location-v0.R` evaluates `TreeLikelihood()` at six points along the ridge, on the
last tree of the pilot chain (`pilot-g`). The points are geometric interpolations t = 0, 0.2, ..., 1 between the chain's
first and last parameters. The rates are mapped to lineage-mean rates by multiplying by exp(-s^2/2).

**Decomposition.** `v0/decomp.R` and `v0/decomp2.R` split each point's error into its window and grid-spacing parts.
Their output is in `v0/decomposition-t*.txt`. Both scripts read the pilot chain from a local path; it is not in the repo.

**Accuracy against grid spacing, same window [-5, 5].** The columns compare K = 26 (h = 0.4) with K = 51 (h = 0.2),
from `v0/decomposition-t*.txt`. There is no decomposition at t = 0.2 or 0.6.

| t | s | secs, K = 26 | logLik K = 26 | logLik K = 51 | difference | old default (`bm`, central), K = 26 vs 51 |
|---|---|---|---|---|---|---|
| 0.0 | 0.96 | 0.31 | -830.972 | -831.097 | **0.125** | 0.200 |
| 0.4 | 1.40 | 6.25 | -825.693 | -825.822 | **0.130** | 0.182 |
| 0.8 | 2.04 | 46-76 | -823.188 | -823.377 | **0.189** | failed |
| 1.0 | 2.46 | 117-187 | -823.624 | -823.859 | **0.235** | failed |

**Notes on the table.**
- **The K = 51, [-7, 7] reference** in `v0/v0-*.csv` (`logLikRef`) has h = 0.28. It is itself 0.04-0.05 from the
  h = 0.2 values, and it does not solve at s >= 2, so it is not a usable reference. The 0.008 at s = 1.69 is off the
  trend: a poor reference solve.
- **Window truncation** ([-7, 7] or an upper bound of 7 at h = 0.4) changes the log-likelihood by <= 0.001 at
  s <= 1.4. It is untested at s >= 2, where the larger window does not solve.
- **Times.** The CSV and the decomposition time the same point differently: at s = 2.04, 46 s and 76 s.
- **Old-default times:** 0.5, 6.5, 53 and 400 s at the first four points, then failure. The mean-rate model equals
  `"bm"` on the shifted window to <= 1.4e-5.

**Top-three-bin mass.**
- **Root prior:** 7.3e-4 at s = 0.96, falling to 1.3e-4 (`top3Prior`).
- **e^x-tilted lineage law** (`top3Tilt`): 1.1e-3, 2.1e-3, 4.4e-3, 9.8e-3, 0.024 and 0.058. Under SqRA this law is
  exactly phi(z_k) e^{s z_k} normalised, so its top-three mass at K = 26 on [-5, 5] is about 1 - Phi(4 - s). That
  exceeds 1e-3 for every s > 0.93, so the criterion was unmeetable at this window from the outset.

**Against the spec's go criteria (§V0).**
- **Time: met** up to s = 2.04; s = 2.46 takes 117-187 s.
- **Accuracy < 0.05 log units: not met at any point.** The miss is 0.125-0.235, all grid spacing. It is smaller than
  the old default's (0.18-0.20).
- **Top-three mass < 1e-3: met for the root prior; not met, and not meetable at this window, for the tilted law.**

The spec says to revise K or the window before implementing when V0 fails; #123 implemented without doing so. The
decision (raise K, widen the window, or relax the thresholds) is the maintainer's, on #117.

## Rungs (`rungs/`)

| Rung | Check | Result |
|---|---|---|
| 1 | Budd & Mann reproduction | Unchanged: the scripts are pinned to `location = "bm"`. `"bm"` chains are `identical()` to `origin/main`'s (`test-Location.R`) |
| 2 | Simulators vs numerics, `EVOTEMPO_LOCATION=mean` (`rung2-mean.rds`) | Agree. The x0 = 2 dCET branch-duration z = -2.6 is the same as in the `bm` run (a known discretisation offset) |
| 3 | beta = 0 collapse (`rung3-collapse.csv`) | Absolute difference <= 4.8e-10. Uses random matrices: says nothing about #116 |
| 3 | Monte Carlo (`rung3-montecarlo*.csv`) | chi-squared p = 0.32, max abs(z) 2.80 |
| 4a(i) | Augmented vs marginal posterior, 20 000 iterations (`rung4a-quantiles.csv`, `rung4a-summary.csv`) | max abs(z) 2.45 over 25 quantiles; KS p >= 0.085 (see notes: chain length chosen after a failure) |
| 4a(ii) | Local moves vs weighted FFBS, first run: 16 chains x 20 000, all iterations (`rung4a-ffbs-pooled16.csv`) | 16 of 16 branch means above FFBS, z 0.94-2.73; occupancy max abs(z) 2.36 where the reference mass is >= 1e-3 |
| 4a(ii) | Pre-registered rerun, 32 chains x 50 000, 64 000 reference draws (`rung4a-rerun/`) | **Pass under both locations.** `"mean"`: local arm max abs(z) 2.25 over 30 quantities, branch-average offset z 1.26 (12 of 16 branches above); diagnostic arm 3.05, offset 2.00. `"bm"`: local 1.54, offset 0.12 (8 of 16); diagnostic 2.05, offset 1.41 |
| — | C++ vs R engine at the new default (`cpp-engine/`) | Max relative difference in log-likelihood 6.7e-10 (<= 1e-8); seeded chains give the same final log posterior in both engines |
| 4b | Likelihood check (`rung4b-likelihood-check.csv`) | abs(z) <= 1.34 |
| 4b | One-parameter exact posterior, 8 chains (`log-onepar2.txt`) | lambda q25 z = -2.02; **psi q50 z = 2.73**; others abs(z) < 2 |
| 4b | One-parameter, s (`log-ones.txt`); psi again at 16 chains x 6000 (`log-onepsi.txt`) | s: abs(z) <= 0.58; psi: abs(z) <= 1.65 (q50 1.46) |
| 4b | One-parameter, `shift` arm on theta, 16 chains (`log-shift.txt`) | max abs(z) 1.985; all five DA z positive (0.44-1.92) |
| 4b | Theta identity (`rung4b-theta-identity.R`), 5 values of theta | abs(z) <= 0.82 without the character; <= 2.27 with it |
| 4b | Shift move at fixed parameters (`rung4b-shift-fixed.R`) | Root-bin posterior: shift arm abs(z) <= 1.73, reference arm <= 1.77 |
| 4c | Bridge Laplace transform | max abs(z) 1.95 |
| 4d | Sum over paths | max abs(z) 2.19. Random matrices: says nothing about #116 |
| 5 | Prior-only MCMC vs simulator (`rung5-*.csv`) | 55 quantiles, max abs(z) 2.00, KS p >= 0.10; discrete `rootBin` = 7: z = 2.40 |
| 5 | Tree moves at fixed parameters (`rung5-fixed-parameters.csv`) | abs(z) <= 1.91 |
| 5 | Start-age oracle (`rung5-age-oracle.R`), root / origin | max abs(z) 1.16 / 1.51 |
| 5 | Start-age validation (`rung5-age-validation.R`) | root: 2.59 over 80 tests; origin: 2.88 over 80; originFree: 2.02 over 86 |
| 6 | Full SBC, 200 replicates per cell, **K = 16 on [-3, 3]** (`sbc/`) | All four cells pass; 52 tests, smallest raw p 0.040, Holm-adjusted p = 1 |

**Notes on the rungs.**

- **Optional stopping.** Three checks were lengthened or repeated after a first failure, and the passing run is the one
  reported:
  - **4a(i):** psi q50 z = 5.7 and s KS p = 0.001 at 5000 iterations; passes at 20 000.
  - **psi in the one-parameter oracle:** q50 z = 2.73 at 8 chains; passes at 16 chains x 6000 (q50 z = 1.46).
  - **4a(ii):** see below.

  psi's median has come out high three times (5.7, 2.73, 1.46). Read the psi results as open, not settled.
- **4a(ii).**
  - **What #123 reported.** It first replaced the full 16 x 20 000 table with its second half and called the excess a
    start-up transient. That was wrong: `TempoMCMC()` starts from a frozen-FFBS draw, close to an exact target draw,
    and the split point was chosen after the failure. The full table is restored here.
  - **What is committed.** Seeds 501-508 and 601-608 (`rung4-local-chains.R 20000 500` and `... 20000 600`); 50 block
    means per chain in `rung4a-local-blockmeans.rds`; the pooling code `rung4a-pool.R`, which reproduces
    `rung4a-ffbs-pooled16.csv`.
  - **The pre-`#116` record.** The old `"bm"` tables show the same sign on several branches, with SEs 3-4 times larger.
  - **The rerun.** Pre-registered in `dev/plans/rung4a-local-rerun.md` and committed (f507e1b) before running, under
    `"mean"` and `"bm"`, with a whole-tree-proposal diagnostic arm. It ran from a clean archive of f507e1b; the logs
    and the outputs record location, scheme, SHA and script md5.
    - **Result.** Every arm passes the pre-registered rule (max abs(z) < 3.3, offset abs(z) < 3) under both locations.
      The offset the first run showed (+0.02 in branch-mean log-tempo, all 16 branches positive) is not reproduced at
      twice the chains and 2.5 times the length.
    - **What remains.** Under `"mean"` the local arm still leans positive: 12 of 16 branches, offset z 1.26. The
      diagnostic arm's largest z is 3.05, within the rule but near it. Read the first run as a fluctuation of the
      size these z-scores allow when 16 correlated branch means share one overall level. The branch-average offset
      is 0.0066 +/- 0.0053 under `"mean"` and -0.0006 +/- 0.0049 under `"bm"`. The first run's 0.019 +/- 0.010 is 1.2
      SE from the former. A bias up to about 0.017 under `"mean"` is not excluded at this size.
  - **Bin 25.** Occupancy bin 25 has reference mass 0 and SE 0, so its z (3.76 for the first 8 chains, 2.89 here) is
    not defined. Excluding it is part of the stated rule (reference mass >= 1e-3).
- **Rung 6 SBC grid and provenance.** The run used the lib built at 2e7836d, whose default location is `"mean"`
  (SqRA). It did not set `EVOTEMPO_GRID`, so `rung6-sbc-full.R` used K = 16 on [-3, 3]. Every one of the 800
  replicate files records `grid` = 16 bins, [-3, 3] and `sha` = 2e7836d (`sbc/provenance.csv`).
  - **Launch.** `sbatch -J sbc --array=1-800 -c 1 --mem=3G -t 08:00:00 --export=ALL,RUN=sbc,RUNDIR=keep,EVOTEMPO_SBC_OUT=<out>,CMD=Rscript data-raw/rung6-sbc-full.R task $SLURM_ARRAY_TASK_ID 1 run.sh`.
    `run.sh` is the wrapper in `rungs/run.sh`. The committed `rung6-sbc-full-array.sh` was not used.
  - **What SBC can show.** Simulation and inference share `.Setup()`, so SBC cannot detect a wrong offset or drift
    centre; the L.3 identity and rung 2 do that.
  - **Not done.** A K = 26 rerun (about 1650 core-hours by the rung 6 README) has not been done.
- **4b one-parameter oracle and the two rung 4b scripts.** These rebuild stored pre-#116 fixtures under
  `EVOTEMPO_LOCATION`.
- **Stale fixtures.** The committed `data-raw/rung4-results/rung4-K26.rds` and `rung5-results/rung4b.rds` are still
  `"bm"` models. They now warn when read (`legacyLocation`). The `"mean"` RDS files from these runs are on Hamilton,
  not committed, so the rung tables here cannot be regenerated from the branch alone.
  `validation-report-default-grid/` still describes the `"bm"` results.
- **Not run:**
  - rung 3 castor MuSSE (not installed on Hamilton);
  - rung 3 RevBayes `dnCDBDP` with the new Q;
  - rung 7 (optional, V4);
  - V5, the Hymenoptera re-pilot, which belongs to the Hymenoptera pipeline.
- **Logs.** `log-*.txt` are the Hamilton job logs; see Provenance above.
