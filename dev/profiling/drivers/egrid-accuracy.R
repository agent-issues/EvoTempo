# Accuracy and cost of the C++ E grid (.EGridCpp: free steps + continuous
# output at 100-fold tighter tolerance) against a tight-tolerance reference
# and the R engine (lsode), on the rung-6 fixture and on the default grid
# (K = 26, z in [-5, 5]) at stress parameters.  Round 1 of dev/profiling/log.md
# records the same table for the earlier stop-at-every-age solver.
# Usage: Rscript dev/profiling/drivers/egrid-accuracy.R   (from the package root)
lib <- Sys.getenv("EVOTEMPO_LIB", "")
if (nzchar(lib)) library(EvoTempo, lib.loc = lib) else library(EvoTempo)
source("dev/profiling/drivers/fixture-rung6.R")
ns <- asNamespace("EvoTempo")
Case <- function(label, par, mdl, maxAge, extra) {
  set <- ns$.Setup(par, mdl)
  ages <- sort(unique(c(seq(0, maxAge, length.out = 401),
                        maxAge * 10^seq(-6, 0, length.out = 61), extra)))
  Run <- function(rtol) {
    c0 <- ns$.SolverCountsCpp()[["steps"]]
    eg <- ns$.EGridCpp(ages, set$Q, set$lambda, set$mu, set$psi, mdl$rho, rtol,
                       1e-4 * rtol)
    c(eg, steps = ns$.SolverCountsCpp()[["steps"]] - c0)
  }
  ref <- Run(1e-13)
  rEng <- local({
    old <- options(EvoTempo.engine = "R")
    on.exit(options(old))
    ns$.EGrid(set$Q, set$lambda, set$mu, set$psi, mdl$rho, maxAge, extra)
  })
  cpp <- Run(1e-10)
  Err <- function(a, v) signif(max(abs(a[[v]] - ref[[v]])), 2)
  cat("\n==", label, "\n")
  print(data.frame(method = c("C++ 1e-10", "R lsode 1e-10"),
                   steps = c(cpp$steps, NA),
                   errE = c(Err(cpp, "E"), Err(rEng, "E")),
                   errG = c(Err(cpp, "G"), Err(rEng, "G")),
                   errDE = c(Err(cpp, "dE"), Err(rEng, "dE"))), row.names = FALSE)
}
Case("rung-6 fixture (K = 16)", truth, model, tt$age[tt$root], tt$age)
def <- DCETModel(rho = 0.5)
Case("default grid, fixture parameters", truth, def, tt$age[tt$root], tt$age)
Case("default grid, theta = 2, s = 1.2", modifyList(truth, list(theta = 2, s = 1.2)),
     def, tt$age[tt$root], tt$age)
Case("default grid, lambda = 2, mu = 1.8, age 30",
     modifyList(truth, list(lambda = 2, mu = 1.8)), def, 30, numeric(0))
