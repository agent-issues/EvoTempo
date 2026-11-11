# PR #78 F3 oracle: pre-registered comparison criteria (written BEFORE any oracle run)

Verifier of review78 (agent-issues/EvoTempo PR #78, head a3ad184). Date: 2026-11-10.

## Question
Under the data-free trilobite FBD model, RevBayes v1.4.1 gives about 100.6 sampled ancestors (SA) and
s = 0.953; dev a2a8e2db gives 6.49 SA and s = 0.0111 (PR prior-check, 2 chains x 200k each). Which is
the correct target, if either?

## Target
This is the FBD tree density with 107 fixed point tip ages (trilobites_tips.csv) and the ingroup
constraint (104 taxa monophyletic, as an indicator). It is integrated over (d, nu, s, origin) under the
priors d ~ LogNormal(-3.5, 1.5), nu ~ U(0,1), s ~ U(0,1) and origin ~ U(0,1000) above the youngest tip,
with rho = 0, removal probability r = 0, origin conditioning and condition on >= 1 sampled lineage.

## Runs
- A (oracle): BEAST 2.4.8 + SA 1.1.7, sampleFromPrior, SABirthDeathModel(origin estimated,
  diversificationRate/turnover/samplingProportion, rho = 0, r = 0, conditionOnSampling = true), with
  MRCAPrior monophyly on the ingroup and the standard SA operators. 4 chains.
  KNOWN DIFFERENCE: BEAST2 measures time from the youngest tip (486.0 Ma); its "present" with rho = 0
  is at the youngest tip, so it cannot put the end of the window at 485.4 Ma. A therefore uses
  offset 0 (window end = youngest tip).
- B (system under test, matched to A): RevBayes dev a2a8e2db. It has the PR's tree block and its tree
  and parameter moves verbatim (trilobites_fbd.Rev lines 41-90), with clock and characters removed
  (data-free, so they do not affect the tree target). It uses offset 0: tip ages - 0.6, youngest tip 0,
  origin ~ U(0, 1000). 4 chains.
- B' (link to the PR): as B, but with the PR's own offset of 0.6 (window end 485.4 Ma). 4 chains.
- V (contrast): RevBayes v1.4.1, same as B (offset 0). 2 chains.
- S (second, fully independent oracle; small exact case): a forward simulation of the FBD process with
  fixed lambda = 1.0, mu = 0.5, psi = 0.4 and origin x0 = 2.5, rho = 0, r = 0, conditioned on exactly
  4 sampled individuals. Statistics: P(#SA = m) for m = 0..3, E[sorted sample ages], E[first-split
  (root) age]. Compared with RevBayes dev (and v1.4.1 as a positive control) on 4 taxa with ages free
  on [0, 2.5], fixed parameters and origin, using the same move types plus mvFossilTipTimeUniform.
  This tests the SA weighting of the density and the moves, but with free rather than point ages.

## Estimators
Per chain, discard the first 10 % as burn-in. ESS per chain uses Geyer's initial monotone sequence;
pooled MCSE = pooled sd / sqrt(sum of chain ESS). The between-chain check is the split R-hat on the
primary statistics.

## Criteria (fixed now)
Primary statistics: the mean SA count and the mean s. Secondary: mean turnover, log d, origin, root
age and ingroup age.
1. Precision gate: MCSE(SA) <= 0.25 and MCSE(s) <= 0.002 in A and in B. That is more than 300x below
   the v1.4.1-dev gap (94 SA, 0.94 s). R-hat <= 1.05. If the gate fails, extend the chains once
   (doubling them) and report the result either way.
2. Agreement A vs B: z = (mean_A - mean_B) / sqrt(MCSE_A^2 + MCSE_B^2). Agreement means |z| <= 4 for
   BOTH primaries. Secondaries with |z| > 4 are reported as flags.
3. Verdict for F3:
   (a) A agrees with B, and V differs from A by more than 10 SA: the dev build is correct for this
       configuration and v1.4.1 is wrong. F3's "nothing independent shows dev is right" is settled.
   (b) A is within 10 SA of V and far (> 10 SA) from B: the dev build is wrong.
   (c) Otherwise: unresolved. Report the size of the discrepancy.
4. B' vs the PR prior-check dev (6.4925 +- 0.128 SA; s 0.01109 +- 0.00037): |z| <= 4 expected, since
   removing data-free components does not change the target. B vs B' measures the offset effect and
   has no criterion.
5. S: for each statistic, |z| <= 4 for dev vs the simulation (z uses MCSE from MCMC and the binomial or
   iid SE from the simulation). The chi-square test on the SA-count distribution has an effective
   sample size (ESS)-adjusted p >= 0.001. v1.4.1 is a positive control; it is expected to fail,
   because of the known age-range collapse/expand bug. If it passes, the test's power is in doubt.


## Amendment 1 (2026-11-10, after the smoke test, BEFORE any production oracle run)
- RevBayes cannot express offset 0. With rho = 0, a tip at age exactly 0 is scored as a present-day
  rho-sample and has probability 0: on the start tree, lnP is -inf at offset 0 and -1096.44 at offset
  0.001 (login-node check, dev build). B therefore uses offset epsilon = 0.001 Myr: tip ages + 0.001
  and origin ~ U(0.001, 1000.001), which is still U(0, 1000) above the youngest tip. The density is
  continuous in the offset, and B' measures the full 0.6 effect, so this is not expected to matter at
  the precision gate. If B and B' differ by more than 4 sigma, the size of the 0.001 offset effect is
  bounded by scaling (0.001/0.6 of the B-B' difference); it is reported rather than assumed.
- BEAST2: sampleFromPrior requires a 'likelihood' node. The XML has none, so the flag is removed;
  the posterior is then the prior only, which is the same target.
- S: dnFBDP needs taxon max age < origin, so the ages are free on [0, 2.49]. The simulation
  likewise rejects realisations with any sample age > 2.49. This conditions both on the same event.
