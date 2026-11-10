# #72 / verify43c pre-registered test: deviations and implementation choices

Written 2026-11-03 11:49 UTC. The 80-run array (Slurm 19337298) was submitted at about 11:49 UTC. No chain result had been seen when this was written.
Pre-registration: prereg-verify43c.md (artifact 77a52c3e). Variant script: oracle-gamma.R (artifact f36f958f).

## Unchanged from the pre-registration
- **Code.** EvoTempo dc47aedffd36ee4e7dbba64da6fcfff487db31ee.
  - data-raw/rung5-age-oracle.R has md5 c353ce92d64d98f81988613f734882f0, verified locally and on Hamilton.
  - oracle-gamma.R has md5 9d06d64c8126acb6eacb12d702448587, verified. Diffed against the committed script: only the three declared edits.
- **Arms, runs, chains, arguments and offsets.** `<case> 60000 32 extant power`.
  - A: 40 runs, offsets 201000..240000.
  - B: 24 runs, offsets 401000..424000, with EVOTEMPO_VERIFY_PRIOR=gamma.
  - C: 16 runs, offsets 301000..316000.
  - Chain seed = 118900 + offset + i, from the script's own seedBase.
- **The arm-B `design` call** ran first and succeeded. The exact grid on (3.4, 12) is computable, with an extant exact mean of 4.64027. It gives no information about chain behaviour.
- **Statistics and decision rules.** As pre-registered, implemented in analyse.R. The rules are written into the code before any data were seen.

## Deviations
1. **Scratch path.** /nobackup/pjjg18/evotempo-bias72/ (with lib/ beneath it), in place of /nobackup/pjjg18/evotempo-verify-pr43c/. This was a task instruction and has no effect on the results.
2. **Install source.**
   - EvoTempo was installed from a package-only subset of the dc47aed GitHub REST tarball: DESCRIPTION, NAMESPACE, R/, src/, man/, inst/, data/, LICENSE.md and .Rbuildignore. data-raw/, dev/ and manuscript/ were dropped. The installed code is identical.
   - Dependencies, including MkPrime for LinkingTo headers, came from the existing project library /nobackup/pjjg18/evotempo-dev/lib.
   - Module: r/4.5.1. The compiler is the default gcc 11.2.
3. **Scripts the pre-registration names did not exist.** verify/analyse.R and the fine-grid script were in the verifier's cell that the safety check held. I wrote them before seeing any chain data, as analyse.R, fine-exact.R and extract.R. Where the pre-registration text left a choice open, I resolved it as follows:
   - **(a) r uses the trapezoid exact sd.** This is the posterior sd that the script's own Exact() would give on its grid (641 points; 1601 for B). fine-exact.R recomputes it, and checks its trapezoid mean against the chains' exact mean.
   - **(b) The per-run z for the KS diagnostic** uses the same exact mean as b: the script's trapezoid mean, unless the Simpson substitution applies. The script's own per-run z is also reported.
   - **(c) The secondary z_p** use the script's trapezoid exact quantiles q_p.
   - **(d) Localisation and start-type quantities** are computed for every outcome and shown as descriptive values. They count as verdicts only if rule 1 says REAL.
   - **(e) SE_r** = SE / exact sd.
4. **First fine-exact attempt failed.** My script had a closure-export bug ('tree' not found on the workers). I fixed it and reran it as a separate job. It is chain-independent.
5. **Order of submission.** The array was submitted before the Simpson check finished. That check does not depend on the chains, so the order has no effect.
6. **Slurm resources.** -c 16, --mem 12G, -t 4 h per run, on `shared`. These choices have no effect on the results.

## Post-run notes
Added after the results were seen. None of these notes changes the analysis.

- **Completion.** All 80 runs finished (exit 0) and every chain kept 12000 samples. No run failed, and no runs were added or replaced. The pooled chains were A 1280, B 768 and C 512.
- **Simpson check.** The Simpson-minus-trapezoid differences in exact mean were A 4.5e-6, B 1.1e-8 and C 3.2e-6. All are below the 0.0002 threshold, so the trapezoid reference was used for b.
- **Transcription check.** fine-exact.R reproduced each chain-side trapezoid mean exactly (difference 0).
- **Acceptance-column bug in extract.R.**
  - The diagnostic column `ageMoveAcceptance` took the first acceptance name containing "Age". That was `nodeAge`, which has 0 acceptance because it is not in the schedule.
  - I re-extracted it in a separate 1-core job, by name (rootAge/originAge), from the same .rds files. Acceptance is 0.284-0.318 in every run.
  - The corrected column is in bias72-perrun.csv v2. It is diagnostic only and not part of any pre-registered statistic.
- **SSH outage.** SSH from the agent to the Hamilton login node failed for a while during the array (connection closed, exit 255). The jobs were unaffected.
- **Hamilton cleanup.** Scratch (225 MB) and the four job workdirs were deleted after the CSVs were saved as artifacts.
