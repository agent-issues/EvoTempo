# Pre-registration: verify43c root-age bias test (PR #43 head dc47aed)

Written 2026-11-03 11:17 UTC, before any chain of this test was submitted. Author: verifier (red-team, peer tier).
Nothing below may be changed after submission; any departure is reported as a deviation.

## Question
Is the small 'young' root-age bias reported for rung 5 (root start, extant conditioning, power
design; pooled mean bias -0.00144 +- 0.00060 over 7 runs, -0.00104 +- 0.00065 over the 6
unselected replicates) real? If it is, does it localise to the prior bound (root-age prior
U(3.4, 7), whose lower bound carries the posterior mode) or to the target?

## Code
- EvoTempo at dc47aedffd36ee4e7dbba64da6fcfff487db31ee (package files identical to 589f866),
  installed from the GitHub tarball into /nobackup/pjjg18/evotempo-verify-pr43c/lib
  (module r/4.5.1, gcc 11.2).
- Arms A and C: data-raw/rung5-age-oracle.R exactly as committed at dc47aed (md5 c353ce92d64d98f81988613f734882f0).
- Arm B: verify/oracle-gamma.R (md5 9d06d64c8126acb6eacb12d702448587), the committed script with three edits only
  (unified diff in the report): with EVOTEMPO_VERIFY_PRIOR=gamma the root-age prior is
  3.4 + Gamma(shape 3, scale 0.6) truncated to (3.4, 12); the exact grid is 1601 points on
  (3.4, 12) (spacing 0.0054, against 0.0056 in arm A); output suffix '-gamma'. Everything else
  (tree, parameters, moves, burn-in 1000, thin 5, nSlice 16, K 16, scoring) is unchanged.

## Design (all power design: lambda = mu = 1, psi 0.15, theta 0.4, s 1, rho 0.35)
Every run = 32 chains x 60000 iterations on 16 cores, `<case> 60000 32 extant power`.
| arm | case | prior | runs | chains | EVOTEMPO_RUNG5_SEED_OFFSET | chain seeds |
|---|---|---|---|---|---|---|
| A (primary) | root | U(3.4, 7) (committed) | 40 | 1280 | 200000 + 1000 j, j = 1..40 | 118900 + offset + i |
| B (bound away from mass) | root | 3.4 + Gamma(3, 0.6) on (3.4, 12) | 24 | 768 | 400000 + 1000 j, j = 1..24 | 118900 + offset + i |
| C (origin control) | origin | U(4.2, 8) (committed) | 16 | 512 | 300000 + 1000 j, j = 1..16 | 118900 + offset + i |
i = 1..32. No seed coincides with any earlier rung-5 run (those used offsets 0-100000, so seeds <= 222932).
Before the chains, one `design` call of arm B's script (exact posteriors only, no chains) checks that
the exact grid on (3.4, 12) is computable; it gives no information about chain behaviour.
A 2561-point Simpson-rule exact posterior mean (extant conditioning) is computed for each arm, as a
check of the 641/1601-point trapezoid reference.

## Statistics (computed by verify/analyse.R, over ALL chains of an arm pooled)
- Primary: mean root (arm C: origin) age bias b = mean(chain means) - exact mean;
  SE = sd(chain means) / sqrt(N chains); z = b / SE. Reference N(0, 1) (N >= 512).
  Relative bias r = b / exact posterior sd.
- Secondary: for p in (0.05, 0.25, 0.5, 0.75, 0.95), z_p = (mean_i F_i(q_p) - p) / (sd_i F_i / sqrt(N)),
  F_i the chain-i empirical CDF at the exact quantile q_p.
- Calibration diagnostic: KS test of the per-run mean z (32 chains each) against t31.
- Exact reference: if the Simpson 2561-point mean differs from the script's trapezoid mean by more
  than 0.0002, the Simpson value is used for b (and this is reported).

## Decision rules (fixed now)
1. Arm A, primary. Bias REAL if |z_A| >= 3 (two-sided p 0.0027).
   NO BIAS OF THE REPORTED SIZE if |z_A| < 3 and the 99% interval b_A +- 2.576 SE_A excludes -0.0010
   (the 6-replicate estimate). Otherwise INCONCLUSIVE (report the interval).
   Power: SE_A ~ 0.0016 / sqrt(40) = 0.00025, so a true bias of -0.0010 gives E[z] ~ -4.0 (power ~84%),
   and -0.00144 gives ~ -5.7.
2. Only if 1 says REAL, localisation:
   - bound: |z_B| < 2 and the relative-bias difference (r_A - r_B) / sqrt(SE_rA^2 + SE_rB^2) has
     |z| >= 2.5 -> 'localises to the bound region'; |z_B| >= 3 with the sign of z_A -> 'not the bound
     (target or move in general)'; otherwise 'localisation inconclusive'.
   - start type: |z_C| >= 3 with the sign of z_A -> 'also at the origin start (not specific to the
     rootAge move)'; |z_C| < 2 -> 'root-start specific'; otherwise inconclusive.
3. Secondary quantile statistics are reported with a flag at |z| >= 3.5; they do not change the verdict.
4. Arms B and C are also each read on their own: |z| >= 3 on the primary statistic is a new finding.
5. All completed runs are used. A run that fails is reported and not replaced; no further runs are
   added after any result is seen. If fewer than 75% of an arm's runs complete, that arm is
   reported as underpowered.

## Budget
About 400 core-hours on `shared` (root run ~5.2 core-h, origin ~3.7). Scratch
/nobackup/pjjg18/evotempo-verify-pr43c/, deleted afterwards except a small results CSV.
