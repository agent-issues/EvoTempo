# CI workflows awaiting installation

The `ms609-agent` token cannot write `.github/workflows/` (GitHub requires the *Workflows* permission
for that path). A maintainer should move these files to `.github/workflows/`:

- `R-CMD-check.yaml`: R CMD check on Ubuntu, Windows and macOS (from the new-methods-project template).
- `render-manuscript.yaml`: renders `vignettes/manuscript.qmd` (template).
- `test-coverage.yaml`: `covr::codecov()`, used with `codecov.yml`'s 100% patch target. Needs the
  `CODECOV_TOKEN` secret.
