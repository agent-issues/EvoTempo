# Rung 4a(ii) rerun: local path moves vs weighted FFBS (pre-registration, review of #123, B1)

This file is committed before any of the runs it describes. The result is reported whole: no subsetting by iteration,
chain or seed after the fact.

## Why

At 16 chains × 20 000 iterations under `location = "mean"`, every one of the 16 branch means came out above the
weighted-FFBS reference (z 0.94-2.73). #123 then reported only iterations 10 001-20 000 and called the excess a
start-up transient. That is wrong: `TempoMCMC()` starts from a frozen-FFBS draw (`.SampleTree()`), which is close to
an exact target draw, and the split point was chosen after seeing the failure.

The old `"bm"` tables have the same sign on several branches, but with standard errors 3-4 times larger. This rerun
asks whether a real offset exists. If so, it asks whether it predates #116, and whether it lies in the local moves or
in the reference.

## Design

**Script.** `data-raw/rung4a-local-rerun.R <location> [nChain] [nIter]`. Run it once per location, `mean` and `bm`.

**Data and model**
- **Data tree and true parameters:** `data-raw/rung4-results/rung4a-rerun-fixture.rds`. This is the tree and truth
  of the `"mean"` rung-4 fixture that the flagged run used, with λ = 0.5, μ = 0.2, ψ = 0.1, θ = 0.1, s = 0.8. Both
  locations use the same tree.
- **Model:** `DCETModel(rho = 0.5, location = <location>)` at the default grid (K = 26, z in [−5, 5]). The scheme is
  the location's default: SqRA for `mean`, central for `bm`. Parameters are fixed at the truth. `nSlice` = 20.

**Arms (seeds fixed in the script)**

| Arm | Runs | Moves | Seeds |
|---|---|---|---|
| Local | 32 chains × 50 000 iterations, `burnin = 100` | `node = 2, tip = 2, branch = 2` | 41 001-41 032 |
| Reference | 16 jobs × 4 000 importance-weighted frozen-FFBS draws (`SampleTempoPaths()`), 64 000 in all | — | 42 001-42 016 |
| Diagnostic | 16 chains × 20 000 iterations, `burnin = 100` | `tree = 1`: whole-tree FFBS proposals with exact MH | 43 001-43 016 |

The diagnostic arm is a third estimator. It shows whether a discrepancy sits in the local moves or in the reference.

**Saved output**
- Per chain: the means of `branchX` (16 branches) and `occupancy` (26 bins) over 50 equal blocks, and the acceptance
  rates.
- Per reference job: the weighted sums and sums of squares needed for the IS mean and its SE.

## Statistics

- **Per quantity** q (16 branch means; occupancy bins whose reference mass is ≥ 1e-3):
  z = (mean over chains − IS reference) / sqrt(sd over chains² / nChain + SE_IS²).
  SE_IS is computed as in `validation-report-phase2.R`.
- **Overall offset:** the branch-average offset (the mean of the 16 branch means, minus the reference's) and its z,
  with the same between-chain SE.

## Decision rule, fixed in advance

**Per location and arm, the arm passes when both hold:**
- max |z| < 3.3 over the quantities above (about 0.05 family-wise over 30-40 correlated tests);
- the overall-offset |z| < 3.

**Reading the outcome:**

| Local arm, `mean` | Local arm, `bm` | Conclusion |
|---|---|---|
| passes | passes | 4a(ii) is reported as passing at this size |
| fails | fails | Pre-existing; a separate issue, not a blocker for #123 |
| fails | passes | A defect introduced by #116; blocks #123 |

The diagnostic arm, set against the reference, then locates any failure:
- the diagnostic arm agrees with the local arm but not with the reference: the reference is suspect;
- the diagnostic arm agrees with the reference: the local moves are suspect.

Hamilton, the `shared` partition, 16 cores per location. Results are committed to
`data-raw/tempo-location-results/rung4a-rerun/`.
