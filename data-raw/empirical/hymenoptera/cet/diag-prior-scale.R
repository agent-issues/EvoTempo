# Hymenoptera CET pipeline: unit checks of the FBD-rate prior options
# (common.R: CET_FBD_SCALE, CET_S_MAX).
#
#   Rscript data-raw/empirical/hymenoptera/cet/diag-prior-scale.R
#
# 1. Default (CET_FBD_SCALE = "none"): CetPriors(0) and CetPriors(1) give the
#    same log prior, equal to FbdLogPrior() plus the theta, s and rate priors.
# 2. LineageMean(s, 1) against e^(-s^2 / 2) (continuous stationary law) and
#    LineageMean(s, 0) = 1.
# 3. "lineage": the joint prior at beta = 1 equals the comparator's density
#    at (k lambda, k mu, psi) written out here, + 2 log k + the s prior, with
#    k = LineageMean(s, 1); at beta = 0 it equals the default.
# 4. The block-move target identity: FbdUPrior(u) - FbdLogPrior(par(u)) is
#    log |d(lambda, mu, psi) / du| = log(psi mu d) (lambda = d (1 + e^u2),
#    mu = d e^u2, psi = mu e^u3, d = e^u1: the determinant is psi (lambda mu -
#    mu^2)), and FbdFromU() inverts FbdToU().
suppressMessages(library(EvoTempo))
Check <- function(ok, what) {
  cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n")
  if (!isTRUE(ok)) assign("failed", TRUE, globalenv())
}
failed <- FALSE
Sys.setenv(CET_FBD_SCALE = "none")
source("data-raw/empirical/hymenoptera/cet/common.R")
set.seed(1)
pars <- lapply(1:20, function(i) list(lambda = exp(stats::runif(1, 2, 8)), mu = NA,
                                      psi = exp(stats::runif(1, -15, -10)),
                                      theta = exp(stats::runif(1, -6, -3)),
                                      s = stats::runif(1, 0.1, 2.9), rate = exp(stats::rnorm(1, -7.1, 0.3))))
pars <- lapply(pars, function(p) { p$mu <- p$lambda - exp(stats::runif(1, -7, -4)); p })
p0 <- CetPriors(0); p1 <- CetPriors(1)
d01 <- max(abs(vapply(pars, function(p) LogPriorJoint(p, p0) - LogPriorJoint(p, p1), 0)))
dman <- max(abs(vapply(pars, function(p) LogPriorJoint(p, p1) - (FbdLogPrior(p) + p1$theta(p$theta) +
  ExpTrunc(1 / 0.587405, 3)(p$s) + p1$rate(p$rate)), 0)))
Check(d01 == 0 && dman < 1e-12, sprintf("1. default: beta 0 vs 1 differ by %g; vs hand sum %g", d01, dman))
lm <- vapply(c(0.25, 0.5, 1, 1.5), function(s) LineageMean(s, 1), 0)
Check(all(abs(lm / exp(-c(0.25, 0.5, 1, 1.5)^2 / 2) - 1) < 0.02) &&
        abs(LineageMean(1.3, 0) - 1) < 1e-12,
      sprintf("2. LineageMean(s, 1) / exp(-s^2 / 2) = %s; LineageMean(1.3, 0) = %.15g",
              paste(signif(lm / exp(-c(0.25, 0.5, 1, 1.5)^2 / 2), 5), collapse = " "), LineageMean(1.3, 0)))
.cet$fbdScale <- "lineage"
q1 <- CetPriors(1); q0 <- CetPriors(0)
Comparator <- function(l, m, p) {
  d <- l - m; r <- m / l; sp <- p / (m + p)
  # density of (d, r, sp) times |d(d, r, sp) / d(l, m, p)| = (l - m) m / (l^2 (m + p)^2)
  stats::dexp(d, 100, log = TRUE) + log(d) + log(m) - 2 * log(l) - 2 * log(m + p)
}
dl <- max(abs(vapply(pars, function(p) {
  k <- LineageMean(p$s, 1)
  LogPriorJoint(p, q1) - (Comparator(k * p$lambda, k * p$mu, p$psi) + 2 * log(k) +
                            ExpTrunc(1 / 0.587405, 3)(p$s) + q1$theta(p$theta) + q1$rate(p$rate))
}, 0)))
d0 <- max(abs(vapply(pars, function(p) LogPriorJoint(p, q0) - LogPriorJoint(p, p0), 0)))
Check(dl < 1e-10 && d0 == 0 && identical(sort(JointNames(q1)), sort(c(fbdNames, "s"))) &&
        identical(sort(JointNames(q0)), sort(fbdNames)),
      sprintf("3. lineage: beta 1 vs hand %g; beta 0 vs default %g; joint names %s / %s", dl, d0,
              paste(JointNames(q1), collapse = ","), paste(JointNames(q0), collapse = ",")))
jac <- max(abs(vapply(pars, function(p) {
  u <- FbdToU(p)
  (FbdUPrior(u) - FbdLogPrior(p)) - (log(p$psi) + log(p$mu) + log(p$lambda - p$mu))
}, 0)))
inv <- max(abs(vapply(pars, function(p) max(abs(unlist(FbdFromU(FbdToU(p))) /
                                                  unlist(p[c("lambda", "mu", "psi")]) - 1)), 0)))
Check(jac < 1e-9 && inv < 1e-9,
      sprintf("4. FbdUPrior(u) - FbdLogPrior(par(u)) - log(psi mu d) = %g; FbdFromU(FbdToU()) relative error %g", jac, inv))
if (failed) stop("diag-prior-scale.R FAILED")
cat("diag-prior-scale.R: all checks passed\n")
