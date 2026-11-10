# Build check and independent oracle (trilobite baseline)

## Build check (data-free, both RevBayes builds)

| File | Content |
|---|---|
| `prior-check.sh` | The Slurm script as run on Hamilton on 2 November 2026 (job 19329128): 2 data-free chains per build, 200 000 iterations sampled every 100. The header says how to repeat it with `run-hamilton.sh` (`UNDER_PRIOR=true`) |
| `prior-summary.R` | Mean, MCSE and ESS per build, after 10 % burn-in. It drops the last log line, which may be incomplete |
| `prior-summary.csv` | Its output; `README.md` (**Build check**) gives the table |
| `rb-versions.txt` | `rb --version` of the two builds |

## Independent oracle (`oracle/`)

The independent reviewer of PR #78 ran this oracle on 10 November 2026; their report is
verify78-report.md. The comparison criteria were pre-registered (`prereg-verify78.md`) before any
production run.

- **Oracle A:** BEAST 2.4.8 with SA 1.1.7, prior only. It uses the same 107 fixed tip ages,
  ingroup constraint and tree priors as this baseline, with ρ = 0, r = 0 and conditioning on
  sampling.
- **Compared with:**
  - the development build `a2a8e2db` (B): the baseline's tree block and moves, without clock and
    characters;
  - B′: as B, at this baseline's own 0.6-Myr window offset;
  - v1.4.1 (V).
- **Second, small exact case (S):** compared against forward simulation of the FBD process.

| Statistic | BEAST2 (A) | `a2a8e2db` (B) | z (A − B) | v1.4.1 (V) |
|---|---|---|---|---|
| sampled ancestors | 6.247 ± 0.052 | 6.282 ± 0.048 | −0.49 | 100.63 ± 0.21 |
| sampling proportion s | 0.01032 ± 0.00014 | 0.01037 ± 0.00012 | −0.28 | 0.9551 ± 0.0018 |
| turnover | 0.99131 | 0.99126 | +0.25 | 0.663 |
| root (Ma) | 519.625 | 519.627 | −0.28 | 519.367 |

All seven statistics agree between A and B (|z| ≤ 0.49; R̂ ≤ 1.003), against a pre-registered
criterion of |z| ≤ 4. In the small exact case S, the development build matches the forward
simulation (|z| ≤ 1.28; χ² p = 0.66), and v1.4.1, the positive control, fails.

The pre-registered verdict is that the development build is correct for this model and v1.4.1
is not. The v1.4.1 bias appears with point ages in a tree-only model. It is therefore a second
sampled-ancestor fault, separate from the interval-age `mvCollapseExpandFossilBranch` bug, and
its cause has not been identified.

Files:
- `prereg-verify78.md`: the pre-registration and its amendment. The 0.001-Myr offset in B was
  set before any production run.
- `oracle-comparison.csv`: every comparison, with z.
- `oracle-summary.csv`: per run and statistic — chains, samples, mean, sd, ESS, MCSE, R̂ and
  chain means.

The oracle's inputs and traces are kept by the review (artifact `verify78-oracle-inputs.tar.gz`,
chain logs `A_chain*.log` and `B_rb_*.log`). They are not committed here.
