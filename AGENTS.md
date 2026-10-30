# EvoTempo: agent development notes

This file orients any agent arriving at the codebase. Read it in full before touching code.

## Start here

1. `BRIEF.md`: scope. `dev/plans/MASTERPLAN.md`: the phase plan (plan, not status).
2. `dev/spec/dCET-spec.md`: the mathematics; every symbol and equation number used in code comments
   refers to it. `dev/spec/bm-numerics-notes.md`: how Budd & Mann's own code discretises the model,
   and the bugs and typos found in it.
3. The validation reports: `data-raw/validation-report/` (rungs 1-2), `data-raw/validation-report-phase2/`
   (rungs 3-6), `data-raw/validation-report-default-grid/` (all rungs at the default grid) and
   `data-raw/cpp-engine-results/` (C++ vs R agreement and run times).
4. Open GitHub issues on `agent-issues/EvoTempo`: the work queue. Read the issue you are picking up
   *and* any issue it links before starting.
5. If your environment provides the `r-conventions` and `coverage-gate` skills, read them before
   writing code and before committing. If it does not, the essentials are under **Naming** and
   **Build and test** below.

## The model in brief

EvoTempo implements the covariant evolutionary tempo (CET) model of Budd & Mann (2025,
Syst. Biol. 74:917-934): every lineage carries a heritable log-tempo x that follows an
Ornstein-Uhlenbeck process (reversion theta, stationary sd s), and tempo e^x scales speciation
(lambda e^(beta x)), extinction (mu e^(beta x)), fossil sampling (psi) and character change
(rate e^(gamma x)). beta and gamma in {0, 1} give the 2x2 nested model family that the empirical
work compares. The package works with a **discretised** model (dCET): x lives on K bins of a grid in
z = x / s (default K = 26, z in [-5, 5], central-difference generator), which turns the model into a
K-state birth-death-sampling process with an OU-like rate matrix `Q`.

Inference has two layers:

- **Marginal likelihood** (`TreeLikelihood()`): tempo integrated out exactly by solving the
  backward E/D equations along each branch of a dated tree with fossil tips and sampled ancestors,
  conditioned on crown survival. `MarginalMCMC()` samples parameters on this likelihood. It is the
  oracle that every augmented result is tested against.
- **Data augmentation** (`TempoMCMC()`): the chain carries an explicit tempo path on every branch
  (bins and jump times), so that characters can be scored on tempo-weighted branch lengths and tree
  moves can be made. Paths are drawn by forward filtering-backward sampling (FFBS) on "frozen"
  time slices of the extinction probability E, with endpoint-conditioned bridges
  (Hobolth & Stone 2009); importance weights correct from the frozen to the exact model.
  Moves (argument `moves`): `tree`, `node`, `tip`, `branch` (path redraws); `joint` (a parameter
  together with a fresh draw of every path; makes theta mix without characters); `shift` (every bin
  +/-1, with `rate`, `lambda`, `mu`, `psi`, `theta` rescaled by their tempo exponents so that the
  character likelihood is unchanged; #19); `nodeAge`, `fossilAge`, `slide`, `saJump` (tree moves); `rootAge`, `originAge` (start-age moves,
  with priors in `agePriors`; #18); plus a scale sweep over parameters. The default schedule is `TempoMoves()`: with characters (and
  gamma != 0) it omits `tree` and `joint`, whose whole-tree draws ignore the characters and are
  almost never accepted, and runs `node`, `tip` and `branch` at 0.25 x the nodes, tips and branches
  per iteration, plus one `shift`. Fractional counts are allowed. `minAcceptance` drops `tree` and
  `joint` after burn-in if their burn-in acceptance is below it (default 1%).

## Code path of one `TempoMCMC()` run

`TempoTree()` (ages, tip types: extant / fossil / ancestor) and `DCETModel()` (grid, scheme, beta,
gamma, conditioning, start) build the inputs; `.ForModel()` marks the tree with `origin = TRUE` for
`start = "origin"`. `.Precompute()` builds per-bin rates and the E grid
(`.EGrid()`, Hermite-interpolated by `.Hermite()`); `.AddSlices()` adds the frozen-E slices and their
eigendecompositions; `.Messages()` runs the backward FFBS pass; `.SampleTree()` / `.SampleClade()` /
`.SampleBranch()` / `.Bridge()` draw paths; `.LogAugmented()` is the exact log augmented density;
`.FrozenLikelihood()` and `.LogWeight()` give the frozen likelihood and importance weights.
`.PathMove()`, `.TreeMove()` (proposals in `R/TreeMoves.R`), `.JointSweep()`, `.ShiftMove()` (map in
`.Shifted()`) and `.ParameterSweep()` are the moves; `.MoveCount()` resolves fractional counts. `MkLikelihood()` scores characters on lengths from `.CharacterLengths()`;
inside `TempoMCMC()` the moves call the closure from `.CharacterLikelihood()`, which under the C++ engine keeps a
persistent MkPrime pruner and recomputes only changed clades (#26; `dev/plans/incremental-pruning.md`). A move that
scores characters must report its decision with `.ResolveCharacters(CharLL, accepted)`.

Infeasible proposals (a bridge or a weight that underflows) signal a condition of class
`infeasibleProposal` via `.Infeasible()`; every move catches it and rejects. Never let one escape
as an error, and never treat it as acceptance.


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
The only exception is a direct push that the maintainer explicitly asks for, for that change only.

### Stacked PRs

When PRs depend on each other, base each on its predecessor's branch; only the bottom of the stack
targets `main`. Merging a PR whose base is not `main` merges it into that base branch, **not** into
`main`, and agents never delete branches, so GitHub does not retarget the next PR automatically.
This has twice left finished work stranded on a feature branch. So:

- after the bottom PR merges, retarget the next PR to `main` (`gh pr edit <n> --base main`) before
  it is merged, and say so in the PR thread;
- when told that work is merged, confirm it is on `main` (`git log origin/main`, or the compare API)
  before reporting it as such, and before stacking new work on it;
- when a branch already merged into `main` has received further merges, bring them over with a
  fresh PR from that branch to `main`.

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
| `R/Simulate.R` | Forward simulators (dCET Gillespie; continuous CET exact in tempo-time); fossils, rho; `BranchClock()` |
| `R/Reconstruct.R` | Sampled and complete trees from a history (sampled ancestors as zero-length tips); characters |
| `R/TreeLikelihood.R` | `TempoTree()` (ages, tip types), `DCETModel()`, marginal likelihood by joint E/D integration, closed-form FBD |
| `R/Paths.R` | Tempo paths: exact augmented density (E and its integral on a Hermite grid), frozen-E slices, FFBS messages, direct-sampling bridges, importance weights |
| `R/MCMC.R` | `TempoMCMC()` (data augmentation), `MarginalMCMC()`, `SampleTempoPaths()`, priors and parameter moves |
| `R/Characters.R` | `MkLikelihood()` (Mk pruning on tempo-weighted lengths; the C++ engine uses MkPrime kernels) |
| `R/Engine.R` | `TempoEngine()`: `engine = "cpp"` (default) / `"R"` switch, option `EvoTempo.engine`; R is the reference oracle |
| `R/Integrate.R` | Thin wrappers over `deSolve` (`lsode`) for the R engine |
| `R/RcppExports.R`, `src/RcppExports.cpp` | Generated by `Rcpp::compileAttributes()`; never edit by hand |
| `src/radau5.h`, `src/band.h`, `src/dcet_ode.h` | Radau IIA (order 5) stiff integrator, banded LU, dCET backward equations (E; E,D; E,G) |
| `src/tree_likelihood.cpp` | C++ tree pruning, `ExtinctionProbability()`, E grid |
| `src/paths.cpp`, `src/paths.h`, `src/frozen.h`, `src/tree.h` | C++ frozen messages, FFBS, bridges, augmented density, weights; read-only views of `set` and `TempoTree` |
| `src/tree_moves.cpp` | C++ kernels of the tree moves: transfer rows, bridge densities and draws, node and tip bin weights, `.Rebuild()` |
| `src/characters.cpp` | `MkLikelihood()` via `MkPrime::Pruner<JCModel>`; the persistent incremental `CharacterCache` behind `.CharacterLikelihood()` |
| `R/TreeMoves.R` | Node-age, fossil-age, subtree-slide and sampled-ancestor moves with frozen-model bin and path proposals |
| `data-raw/bm-*.R` | Budd & Mann reproduction (`bm/` holds their verbatim code), convergence study, clade-size inversion |
| `data-raw/rung2-validation.R`, `data-raw/validation-report.R` | Monte Carlo validation and the rung 1-2 report |
| `data-raw/rung3-validation.R` | Rung 3: beta = 0 collapse, castor MuSSE comparison, Monte Carlo on tiny trees |
| `data-raw/rung3-revbayes.R`, `-summary.R`, `-engines.R` | Rung 3 vs RevBayes `dnCDBDP` (`REVBAYES` = path to `rb`); conventions and results in `data-raw/rung3-results/revbayes/README.md` |
| `data-raw/rung4-validation.R` | Rung 4: augmented vs marginal posterior, FFBS summaries, bridge Laplace transform, sum over paths |
| `data-raw/rung5-validation.R` | Rung 4b (one character, product-space oracle) and rung 5 (prior-only MCMC vs simulator) |
| `data-raw/cpp-engine-agreement.R` | C++ vs R agreement and run times on the rung 3/4 fixtures (K = 16, 26); results in `data-raw/cpp-engine-results/` |
| `data-raw/empirical/` | Curated tip-dated datasets (penguins, Hymenoptera, trilobites): raw downloads with checksums, NEXUS matrices, tip-age tables, constraints, per-dataset README; `build.R` regenerates everything from `raw/` |
| `manuscript/` | Quarto manuscript (build-ignored; rendered by `render-manuscript.yaml`). Not in `vignettes/`: R CMD check treats a `.qmd` there as a vignette without an engine |
| `tests/testthat/test-Engine.R` | C++ vs R agreement for every ported function, and identical seeded chains in both engines |
| `data-raw/rung4b-likelihood-check.R` | Rung 4b: importance-sampled likelihood from FFBS maps vs the exact product-space likelihood at fixed parameters |
| `data-raw/rung5-fixed-parameters.R` | Rung 5 at fixed parameters (tree moves only) vs rejection simulation; regression check for tree-move Hastings ratios |
| `data-raw/rung4b-one-parameter.R` | One free parameter at a time vs the exact 1-D posterior from a likelihood grid (rung 4b tree, or the rung-4 tree at the default grid) |
| `data-raw/rung4b-shift-fixed.R` | Shift move at fixed parameters: root-bin posterior given the rung 4b tree and character vs the exact posterior from point-mass start priors (#19) |
| `data-raw/mixing-schedule.R` | ESS per CPU hour of `TempoMCMC()` under the pre-#19 schedule and `TempoMoves()` on the `dev/mixing/` fixed trees; results and the shift-move oracles in `data-raw/mixing-results/README.md` |
| `data-raw/rung5-age-oracle.R` | Start-age moves on a fixed tree at fixed parameters vs the exact 1-D posterior of the root or origin age (prior x `TreeLikelihood()` on a grid) (#18) |
| `data-raw/rung5-age-validation.R` | Rung 5 with a free start age (#18): prior-only MCMC vs rejection simulation for a free crown age, a fixed origin and a free origin, at fixed or free parameters; results in `data-raw/rung5-age-results/README.md` |
| `data-raw/rung4b-theta-identity.R` | Importance-sampling identity (augmented to marginal likelihood) across theta, with or without the character |
| `data-raw/rung6-sbc.R` | Rung 6: simulation-based calibration per (beta, gamma) cell on fixed true time trees with characters |
| `data-raw/rung6-sbc-full*.{R,sh,py}` | Rung 6 full size (#21): SBC with tree moves, 200 replicates per cell as Slurm array tasks on Hamilton; results and verdict in `data-raw/rung6-full-results/` |
| `data-raw/rung7-power.R`, `-array.sh`, `-summary.R`, `-figure.py` | Rung 7: power of (beta, gamma) model selection by marginal likelihoods on fixed true trees (quadrature, bridge sampling, character-power stepping-stone, posterior-to-posterior bridge); results and caveats in `data-raw/rung7-results/README.md` |
| `data-raw/validation-report-phase2.R` | Summary CSVs for the rungs 3-6 report (`data-raw/validation-report-phase2/README.md`) |

## Moves that change parameters and paths together

A move that maps (parameters, paths) deterministically, like `shift`, needs the Jacobian of the map in
the Hastings ratio (`prod p' / p` for multiplicative parameter changes), and must reject, not clamp,
proposals that leave the grid. The exact one-parameter oracle (`data-raw/rung4b-one-parameter.R`,
seventh argument `shift`) checks each compensation separately; batch-means z-scores there scatter
by about +/-2 at 8 chains, so use 16 or more.

## Origin start and start-age moves

With `start = "origin"` (#18) the root edge (`tree$stem`) carries a path: `paths[[root]]` runs from
the origin to the root, the origin bin is element `nTip + nNode + 1` of the bin vector, `.PathNodes()`
includes the root, `.Messages()` returns the message at the origin as `top`, and the root is an
ordinary bifurcation (lambda factor). Read the top of any branch with `.TopAge()` / `.TopBin()` (R) or
`TreeView::TopAge()` / `TopIndex()` (C++), never `age[parent[v]]`, which is `NA` at the root.
Survival is that of the single stem lineage (`set$survive` unsquared); a sampled ancestor at the root
is not supported. `nodeAge` moves the root below a fixed origin; `originAge` moves the origin; with
`start = "root"`, `rootAge` moves the root. A start-age move rebuilds the frozen model (its E grid
and slices end at the start age), so its reverse bridge density uses a frozen model built afresh at
the current tree, not a stale `st$set`. The conditioning is evaluated at the current start age, so
the target is the age prior times the survival-conditioned likelihood (spec section 6.1).

## Tree-move proposals

The reverse proposal density of a tree move must include every path the forward move
redraws, including branches whose endpoints are unchanged but whose parent changes
(the regraft target of `slide`). Omitting one is invisible when paths rarely jump
(coarse grid, small theta) and shows up as old node ages on finer grids;
`data-raw/rung5-fixed-parameters.R` detects it.

## Default grid

Validation scripts run at the package default grid (K = 26, z in [-5, 5]) when
`EVOTEMPO_GRID=default` is set, writing `-K26` result files; the summary is
`data-raw/validation-report-default-grid/README.md`. Run such jobs one at a time:
two 15-worker clusters at once have lost all their workers on the development
machine.

## Validation simulators

The likelihood is conditioned on survival, so a simulator used as a reference must
resimulate at fixed parameters until the crown survives; redrawing parameters on
failure weights them by P(survival) (seen in rung 5 as biased lambda and mu).
Rejection on other data features (tip counts, tree size) may redraw parameters.

## RevBayes comparisons

`dnCDBDP` scores a sampled ancestor as psi * E(t) rather than psi, and its tree-shape term uses a
Stirling approximation to log(n!) that is off by 5e-4 at n = 1.  Rev rounds Reals to 6 significant
digits in string concatenation; print values on their own line.  Its fixed-step integrator needs
`nTimeSlices` of 16000-64000 at the default grid.

## Naming

See `r-conventions`: functions BigCamelCase, variables camelCase, private functions `.DotPrefix`.
Mathematical symbols keep their spec names in documentation (lambda, mu, psi, rho, theta, s, alpha, beta,
gamma, delta).

## C++ engine

`engine = "cpp"` is the default; `"R"` runs the original R code, which stays as the reference
implementation. Rules for any change to either side:

- Both engines must give the same log-likelihoods to <= 1e-8 relative on the rung 3/4 fixtures
  (`data-raw/cpp-engine-agreement.R`) and **identical** samples from seeded `TempoMCMC()` chains
  (`test-Engine.R`). A change that alters random-number consumption on one side breaks the second
  check: change both sides together.
- New C++ is exported through `// [[Rcpp::export]]`; regenerate the glue with
  `Rcpp::compileAttributes()` and the docs with `roxygen2::roxygenise()`.
- `LinkingTo: MkPrime` uses the header-only API in MkPrime's `inst/include/MkPrime/`
  (`Pruner<Model>`, edge-aware transitions). MkPrime is installed from `agent-issues/MkPrime`
  (`Remotes:`). Engine-agnostic MkPrime improvements belong in MkPrime, not here.
- Profile before optimising. As first measured, the C++ kernels ran 8-70x faster than R while a full
  `TempoMCMC()` iteration gained only 6-8x, because the loop, the tree moves and the slice
  eigendecompositions still run in R (issue #26).

## CI

`.github/workflows/`: `R-CMD-check.yaml` (ubuntu, macOS, windows), `test-coverage.yaml` (covr,
uploads to Codecov) and `render-manuscript.yaml` (on changes under `manuscript/`). A failure
that appears only in `test-coverage` usually means a test depends on the random-number state or on
round-off: covr runs the tests in a different order and state from R CMD check.

## Compute

- **Local Windows sandbox** (Claude Science): no usable C++ compiler, so the package cannot be
  installed from source there. Pure-R work runs with `pkgload::load_all()`. `covr::package_coverage()`
  needs `sh` and fails; use `covr::file_coverage()` on the R and test files. `R CMD build/check` need
  `TMPDIR` set to a writable folder. Run heavy validation jobs one at a time, with at most four cores:
  two 15-worker clusters at once have killed every worker on the development machine.
- **Hamilton HPC** (Durham; SSH over the university VPN; Slurm): compiling and longer runs. Use
  partition `shared` (`test` for jobs under 15 min), `module load r/4.5.1 gcc/14.2`,
  `OMP_NUM_THREADS=1`, work directories under `/nobackup/pjjg18/`, and the project library
  `.libPaths(c('/nobackup/pjjg18/evotempo-dev/lib', .libPaths()))`. Ask before installing new packages.
  `Permission denied (publickey)` means the VPN has lapsed; do not touch keys.
- **RevBayes** for cross-checks: on the maintainer's machine,
  `C:\Users\pjjg18\GitHub\revbayes-ms\build\rb.exe` (branch `development`; no `dnMkPrime`, which
  lives on `feature/dnmkprime-pr1-native-integration`). Its source lists are generated: after pulling,
  rerun `projects/cmake/regenerate.sh` before building.

## Validation ladder

Issues #1-#7 are the validation rungs: (1) reproduce Budd & Mann; (2) simulators vs analytic moments;
(3) `TreeLikelihood()` vs closed forms, castor MuSSE, Monte Carlo and RevBayes `dnCDBDP`;
(4) augmented vs marginal posteriors and path-sampler identities; (5) prior-only MCMC vs the
simulator, and tree moves at fixed parameters; (6) simulation-based calibration per (beta, gamma)
cell; (7) power study. Any change to the likelihood, the sampler or a move must keep the relevant
rung scripts passing; rerun them (in `data-raw/`, not `tests/`) and update their result files.

## Testing pitfalls

- Tests use small, seeded runs. A test that needs a rare branch of the code (an infeasible bridge,
  an underflow) must force it, for instance with `testthat::local_mocked_bindings()`, rather than
  hope that a random draw or round-off produces it.
- Conditioning on survival affects references: see **Validation simulators** above.
- The exact one-parameter posterior on a likelihood grid (`data-raw/rung4b-one-parameter.R`) is the
  sharpest check that an MCMC discrepancy is mixing error rather than a wrong target.
