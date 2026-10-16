# Validation at the default grid (K = 26, z in [-5, 5]; issue #17)

All scripts take `EVOTEMPO_GRID=default` and write `-K26` result files. `validation-report-phase2.R` run with the same variable writes the tables in this directory. Each rung reuses the K = 16 design, data tree and character, so the two grids are compared on the same data.

## Deviations
1. **castor (rung 3b) was not run.** At K = 26, one 12-tip likelihood takes more than 15 min at z in [-3, 3], and castor 1.8.7 exits without a message at z in [-5, 5].
2. **The rung 3 Monte Carlo caps daughter clades at 20000 branches** and counts capped clades as "other". Such clades cannot be among the tiny configurations counted, so the cap is safe.
3. **The rung-4b posterior comparison was not run.** The product space has 52 states, which makes the marginal chains impractically slow. Rung 4b is represented by its fixed-parameter likelihood check.
4. **Rung 6 (SBC) was not rerun.** It was outside the scope of #17.
5. **Rung 4a(i) was run twice.** The first run used the original moves; the second added `joint` moves (see #20). Both runs share the same marginal chains.

## Results
| Rung | Check | Result | Verdict |
|---|---|---|---|
| 3a | beta = 0 collapse (grid-independent) | max abs(delta log L) = 7.1e-9 | Pass |
| 3c | Monte Carlo, 2.8 million crown pairs, 19 cells | chi2 = 17.2 on 19 df (p = 0.58), max abs(z) = 2.33 | Pass |
| 4a(i) | augmented vs marginal, 8 + 8 chains x 5000, original moves | 23/25 quantiles abs(z) <= 2; theta q75 z = -2.25; theta ESS 326 vs 6818 | Slow theta mixing |
| 4a(i) | same, with joint moves | 24/25 abs(z) <= 2; mu q25 z = 3.58; KS p >= 0.067 except psi (0.010); theta ESS 7564 | Pass with one flag (below) |
| 4a(i) | one free parameter vs exact 1-D posterior (mu, psi) | augmented abs(z) <= 2.07, marginal abs(z) <= 1.63 | Pass |
| 4a(ii) | local moves vs FFBS maps, 36 summaries | max abs(z) = 1.12 | Pass |
| 4c | bridge Laplace transform, 9 cases | max abs(z) = 1.57 | Pass |
| 4d | Proposition 3.3 (grid-independent) | 16 entries, max abs(z) = 2.19 | Pass |
| 4b | fixed-parameter likelihood vs product space, 6 points | abs(z) <= 1.43 for 5 points; one at -2.29 | Pass |
| 5 | prior-only MCMC vs 3200 simulated trees (105468 prior draws; 3483 exploded) | 51/55 quantiles abs(z) <= 2 (max 2.58); KS p >= 0.095 except theta (0.033); cherry counts abs(z) <= 2.21 | Pass |
| 5 | fixed parameters, tree moves only | 8 statistics, abs(z) <= 1.30 | Pass |

**The flag in 4a(i).** The single mu q25 value (z = 3.58) rests on between-chain standard errors from 8 chains. Two further checks argue against a real discrepancy:
- with one free parameter, the augmented and marginal samplers both match the exact posterior for mu and psi at this grid;
- the psi KS value comes from the same marginal chains in both runs, so it is not independent evidence.

Neither rung fails. The one-parameter checks are in `data-raw/rung5-results/rung4b-one-parameter-mu-psi-3000-nochar-joint-rung4-K26.rds`.
