# Trilobite CET pipeline: Monte Carlo check of the priors matched to the
# RevBayes comparator (common.R, CetPriors()).
#
#   Rscript data-raw/empirical/trilobites/cet/check-priors.R [out dir]
#
# Draws 2e6 values of (d, nu, s) and of the UCLN CV from the comparator's
# priors, maps them to (lambda, mu, psi) and sigma, and compares their
# quantiles with the CDFs of the CET priors (numerical integration of
# exp(prior)).  Checks, stopping on failure:
# 1. lambda, mu, psi: |empirical CDF - CET CDF| <= 0.002 at the empirical
#    0.1%, 1%, 5%, 25%, 50%, 75%, 95%, 99%, 99.9% quantiles (that inside the
#    truncation, renormalised; the Monte Carlo SE is <= 3.5e-4);
# 2. s: the same for the UCLN log-SD prior;
# 3. JointFbdLogPrior() integrates to 1 over (lambda, mu, psi): checked as
#    E[exp(JointFbdLogPrior - log q)] = 1 under a broad proposal q (Student-t(3) on the log scale),
#    to within 4 Monte Carlo SE.
# Writes <out dir>/prior-check.csv (default <cetDir>/smoke-results).
suppressMessages(library(EvoTempo))
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
outDir <- if (length(args)) args[[1]] else file.path(cetDir, "smoke-results")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
set.seed(4711)
n <- 2e6
d <- stats::rlnorm(n, -3.5, 1.5)
nu <- stats::runif(n)
s <- stats::runif(n)
lambda <- d / (1 - nu)
mu <- nu * lambda
psi <- s * mu / (1 - s)
cv <- stats::rexp(n, 3)
sigma <- sqrt(log1p(cv^2))
pr <- CetPriors()
probs <- c(0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999)
fail <- character()
rows <- list()
for (p in c("lambda", "mu", "psi", "s")) {
  x <- switch(p, lambda = lambda, mu = mu, psi = psi, s = sigma)
  f <- pr[[p]]
  if (p == "s") {
    lo <- 1e-9
    hi <- 3
    inside <- x <= hi
  } else {
    lr <- attr(f, "logRange")
    lo <- exp(lr[1])
    hi <- exp(lr[2])
    inside <- x >= lo & x <= hi
  }
  xi <- x[inside]
  # CDF of the CET prior on a fine log grid
  y <- seq(log(lo), log(hi), length.out = 40001)
  dens <- exp(vapply(exp(y), f, 0)) * exp(y)
  cdf <- cumsum(c(0, (dens[-1] + dens[-length(dens)]) / 2 * diff(y)))
  qs <- stats::quantile(xi, probs, names = FALSE)
  cet <- stats::approx(y, cdf / cdf[length(cdf)], log(qs))$y
  rows[[p]] <- data.frame(parameter = p, prob = probs, quantile = qs, cetCdf = cet,
                          difference = cet - probs, massOutside = mean(!inside),
                          integral = cdf[length(cdf)])
  if (any(abs(cet - probs) > 0.002)) fail <- c(fail, paste(p, "marginal differs from the comparator's"))
}
res <- do.call(rbind, rows)
print(res, digits = 4, row.names = FALSE)
# Joint density normalisation by importance sampling
m <- 1e6
# Student-t(3) proposal on the log scale: heavier tails than the target's
# exponential tails in log psi and log lambda, so the weights have finite variance
LogT <- function(x, loc, sc) stats::dt((x - loc) / sc, 3, log = TRUE) - log(sc)
ql <- -1 + 2.5 * stats::rt(m, 3)
qm <- -1.5 + 2.5 * stats::rt(m, 3)
qp <- -3 + 3 * stats::rt(m, 3)
lq <- LogT(ql, -1, 2.5) - ql + LogT(qm, -1.5, 2.5) - qm + LogT(qp, -3, 3) - qp
lj <- mapply(JointFbdLogPrior, exp(ql), exp(qm), exp(qp))
w <- exp(lj - lq)
est <- mean(w)
se <- stats::sd(w) / sqrt(m)
cat(sprintf("Joint prior normalisation: %.4f (SE %.4f)\n", est, se))
if (abs(est - 1) > 4 * se) fail <- c(fail, "JointFbdLogPrior() does not integrate to 1")
res <- rbind(res, data.frame(parameter = "joint", prob = NA, quantile = NA, cetCdf = est,
                             difference = est - 1, massOutside = NA, integral = se))
write.csv(res, file.path(outDir, "prior-check.csv"), row.names = FALSE)
if (length(fail)) stop("check-priors.R FAILED: ", paste(fail, collapse = "; "))
cat("check-priors.R: all checks passed\n")
