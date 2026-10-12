# EvoTempo: agent development notes

You MUST read the `r-conventions` skill before writing any code, and `coverage-gate` before committing.
Read `BRIEF.md` for scope and `dev/plans/MASTERPLAN.md` for the phase plan.

## GitHub: the `agent-issues` mirror

Two remotes, which are not interchangeable:

| Remote | Repo | Role |
|--------|------|------|
| `origin` | `agent-issues/EvoTempo` (fork) | **All agent work.** Branches, PRs, issues, Discussions, CI. |
| `upstream` | `ms609/EvoTempo` | Release repo. Fetch only; push URL is `no-push-use-gha`. |

The fork's default branch is `main`, so `Fixes #N` in a PR body closes the issue on merge.

**Never push to `upstream`, under any identity.** Work reaches `ms609/EvoTempo` only when the human
syncs the fork.

**Everything you create on GitHub must be authored by `ms609-agent`.** The token is `CLAUDE_GH_TOKEN`.
Apply it per command, never globally:

```bash
GH_TOKEN=$CLAUDE_GH_TOKEN gh pr create --base main --head <branch> --reviewer ms609 ...
git -c user.name=ms609-agent -c user.email=313734811+ms609-agent@users.noreply.github.com commit ...
```

`git push` must go out as `ms609-agent`. Use a one-shot credential, with an empty `credential.helper=`
first so that the push cannot fall back to another identity. Where a shell credential helper is
unavailable (e.g. sandboxed Windows), pass a per-command header instead:
`git -c credential.helper= -c "http.https://github.com/.extraheader=AUTHORIZATION: basic <base64 ms609-agent:TOKEN>" push origin <branch>`.
Never print, log or commit the token. No force-pushes; no branch deletions.

**Never commit to `main`.** Every change, documentation included, goes on a branch and reaches `main`
through a PR with `ms609` as reviewer. Agents do not merge.

## Where things live

| Record | Home |
|--------|------|
| Tasks, bugs, validation rungs | GitHub issues (labels `validation`, `chore`, `red-team`, `sev:*`) |
| Work in progress | `in-progress` label on the issue; PR linked with `Fixes #N` |
| Blockers needing a human | `needs-escalation` label |
| Design discussion, round records | GitHub Discussions |
| Phase plan | `dev/plans/MASTERPLAN.md` (plan, not a status ledger) |
| Mathematical specification | `dev/spec/` |

**No committed status ledgers.** No file carries a status column, a to-do list or a completed-tasks log.
GitHub issues are the queue; their open/closed state is the status.

## Build and test

Never install to the default library. Use a per-agent library:

```bash
R CMD INSTALL --library=.agent-X .
Rscript -e ".libPaths(c('.agent-X', .libPaths())); testthat::test_dir('tests/testthat')"
```

Regenerate documentation with a plain `roxygen2::roxygenise()`. Coverage target: 100% of reachable lines
(`covr::package_coverage()`); `# nocov` only where truly unreachable. Long-running validation
(Monte Carlo, convergence studies) lives in `data-raw/`, not in `tests/`; tests use small, seeded runs.

## Package layout

| Path | Content |
|------|---------|
| `R/Generator.R` | Grid, generator schemes, Budd & Mann's discretization, per-bin rates |
| `R/Extinction.R` | E_k(t) with rho, psi; analytic FBD p0 |
| `R/Moments.R` | Mean diversity, second moment, survivor means, tempo distribution, equilibrium, branch moments |
| `R/CladeSize.R` | Clade-size distribution by inversion of the probability generating function |
| `R/Simulate*.R` | Forward simulators (dCET Gillespie; continuous CET exact in tempo-time); sampling; characters |
| `data-raw/` | Budd & Mann reproduction, convergence and Monte Carlo validation scripts and reports |

## Naming

See `r-conventions`: functions BigCamelCase, variables camelCase, private functions `.DotPrefix`.
Mathematical symbols keep their spec names in documentation (lambda, mu, psi, rho, theta, s, alpha, beta,
gamma, delta).
