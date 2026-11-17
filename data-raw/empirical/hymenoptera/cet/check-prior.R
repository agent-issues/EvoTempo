# Hymenoptera CET pipeline: check that TempoMCMC(), with common.R's joint
# FBD-rate prior (FbdConditional(), ShiftJoint(), LogPriorJoint()), samples
# the intended posterior of (lambda, mu, psi).
#
#   Rscript data-raw/empirical/hymenoptera/cet/check-prior.R <out dir> [cores]
#
# Target: a fixed 11-tip subtree of the start tree (Orthoptera, Paraneoptera,
# six extant Hymenoptera, three fossils at fixed ages), no characters, theta
# and s fixed, lambda, mu, psi free, rho = 0.0005, crown start, 16 bins
# z in [-3, 3].  Posterior p(lambda, mu, psi | T) ~ TreeLikelihood(T) x
# prior, the prior being the comparator's d ~ Exp(100), turnover ~ U(0, 1),
# sampling proportion ~ U(0, 1).
# - Oracle: random-walk Metropolis directly on u = (log d, logit turnover,
#   logit sampling proportion), whose prior is the product of the three
#   priors (times the log/logit Jacobians), with TreeLikelihood() (paths
#   integrated out) as the likelihood.  It uses neither FbdLogPrior() nor
#   its Jacobian, so it checks the algebra too.
# - CET: TempoMCMC() with the default moves (tree, node, tip, branch, joint,
#   shift, parameter sweep), on the augmented target.  beta = 0: the shift
#   move does not touch lambda, mu (only sweeps and joint moves see the
#   prior); beta = 1: the shift move rescales lambda and mu together
#   (ShiftJoint()).
# - Negative control: the same TempoMCMC() runs with the Jacobian
#   |d(d, r, s) / d(lambda, mu, psi)| left out of FbdLogPrior(), a plausible
#   implementation error.  The check must reject it.
# Compared: posterior means and 0.1, 0.5, 0.9 quantiles of log d, logit
# turnover, logit sampling proportion and log lambda;
# z = difference / sqrt(MCSE_CET^2 + MCSE_oracle^2) (MCSE of a quantile by
# the indicator-ESS and delta method).  CET_PC_REUSE=1 recomputes the
# statistics from <out dir>/prior-check-draws.rds.  Pass: every |z| < 3.5 for beta = 0
# and 1; the control must have some |z| >= 3.5.  Writes <out dir>/prior-check.csv
# and stops with an error if a check fails.
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/hymenoptera/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
outDir <- args[[1]]
cores <- if (length(args) >= 2) as.integer(args[[2]]) else 1L
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
nChain <- as.integer(Sys.getenv("CET_PC_CHAINS", "4"))
nIterCet <- as.integer(Sys.getenv("CET_PC_NITER", "20000"))
nIterOracle <- as.integer(Sys.getenv("CET_PC_ORACLE", "40000"))
burnFrac <- 0.2

# ---- Target -----------------------------------------------------------------------
full <- read.tree(file.path(cetDir, "data", "start-tree.nwk"))
tips <- read.csv(file.path(hymDir, "hymenoptera_tips.csv"))
info <- tips[match(full$tip.label, tips$taxon), ]
hymExtant <- full$tip.label[info$status == "extant" & info$group == "Hymenoptera"]
fossils <- full$tip.label[info$status == "fossil"]
keep <- c("Orthoptera", "Paraneoptera", hymExtant[1:6], fossils[1:3])
sub <- drop.tip(full, setdiff(full$tip.label, keep))
si <- tips[match(sub$tip.label, tips$taxon), ]
depth <- node.depth.edgelength(sub)
extant <- si$status == "extant"
age <- max(depth[which(extant)]) - depth
age[which(extant)] <- 0
tipType <- ifelse(extant, "extant", "fossil")
tree <- TempoTree(sub, ages = age, tipType = tipType)
MakeModel <- function(beta) {
  DCETModel(nBin = 16, zRange = c(-3, 3), beta = beta, gamma = 0, rho = 5e-4,
            condition = "extant", start = "root")
}
fixed <- list(theta = 1 / 120, s = 0.5)
cat(sprintf("Subtree: %d tips (%d fossils), root %.3f Ma\n", tree$nTip,
            sum(tipType == "fossil"), max(tree$age)))

FromU <- function(u) {
  d <- exp(u[1])
  r <- stats::plogis(u[2])
  sp <- stats::plogis(u[3])
  lambda <- d / (1 - r)
  mu <- r * lambda
  list(lambda = lambda, mu = mu, psi = sp * mu / (1 - sp))
}
Derived <- function(lambda, mu, psi) {
  data.frame(logD = log(lambda - mu), logitTurnover = log(mu) - log(lambda - mu),
             logitSampling = log(psi) - log(mu), logLambda = log(lambda),
             turnover = mu / lambda, sampling = psi / (mu + psi))
}
Oracle <- function(beta, seed) {
  set.seed(seed)
  m <- MakeModel(beta)
  LogPost <- function(u) {
    d <- exp(u[1])
    r <- stats::plogis(u[2])
    sp <- stats::plogis(u[3])
    prior <- log(100) - 100 * d + u[1] + log(r) + log1p(-r) + log(sp) + log1p(-sp)
    ll <- tryCatch(TreeLikelihood(tree, c(FromU(u), fixed), m), error = function(e) -Inf)
    prior + ll
  }
  u <- c(log(0.01), stats::qlogis(0.9), stats::qlogis(0.1))
  lp <- LogPost(u)
  S <- diag(c(0.5, 0.5, 0.5)^2)
  out <- matrix(NA_real_, nIterOracle, 3)
  acc <- 0
  for (i in seq_len(nIterOracle)) {
    prop <- u + drop(t(chol(S)) %*% stats::rnorm(3))
    lq <- LogPost(prop)
    if (log(stats::runif(1)) < lq - lp) {
      u <- prop
      lp <- lq
      acc <- acc + 1
    }
    out[i, ] <- u
    # Adapt the proposal covariance during the first burn-in fraction only
    if (i <= burnFrac * nIterOracle && i %% 500 == 0 && i >= 1000) {
      S <- 2.38^2 / 3 * stats::cov(out[(i %/% 2):i, ]) + diag(1e-6, 3)
    }
  }
  keepI <- seq_len(nIterOracle) > burnFrac * nIterOracle
  p <- lapply(seq_len(sum(keepI)), function(i) FromU(out[keepI, ][i, ]))
  res <- Derived(vapply(p, `[[`, 0, "lambda"), vapply(p, `[[`, 0, "mu"), vapply(p, `[[`, 0, "psi"))
  attr(res, "acceptance") <- acc / nIterOracle
  res
}
Cet <- function(beta, seed, control = FALSE) {
  set.seed(seed)
  if (control) {
    # Negative control: no Jacobian (FbdConditional() and ShiftJoint() look
    # FbdLogPrior() up in the global environment)
    assign("FbdLogPrior", function(par) {
      l <- par[["lambda"]]; m <- par[["mu"]]; p <- par[["psi"]]
      if (!(m > 0 && l > m && p > 0)) return(-Inf)
      log(100) - 100 * (l - m)
    }, envir = globalenv())
  }
  m <- MakeModel(beta)
  pri <- CetPriors()[fbdNames]
  start <- c(list(lambda = 0.1, mu = 0.09, psi = 0.01), fixed)
  run <- TempoMCMC(tree, m, start, pri, nIter = nIterCet, burnin = 1000, thin = 5)
  smp <- run$samples
  smp <- smp[seq_len(nrow(smp)) > burnFrac * nrow(smp), ]
  res <- Derived(smp$lambda, smp$mu, smp$psi)
  attr(res, "acceptance") <- run$acceptance
  res
}

jobs <- c(lapply(seq_len(nChain), function(k) list(kind = "oracle", beta = 0, seed = 100 + k)),
          lapply(seq_len(nChain), function(k) list(kind = "oracle", beta = 1, seed = 200 + k)),
          lapply(seq_len(nChain), function(k) list(kind = "cet", beta = 0, seed = 300 + k)),
          lapply(seq_len(nChain), function(k) list(kind = "cet", beta = 1, seed = 400 + k)),
          lapply(seq_len(nChain), function(k) list(kind = "control", beta = 0, seed = 500 + k)),
          lapply(seq_len(nChain), function(k) list(kind = "control", beta = 1, seed = 600 + k)))
rdsFile <- file.path(outDir, "prior-check-draws.rds")
if (nzchar(Sys.getenv("CET_PC_REUSE")) && file.exists(rdsFile)) {
  # Recompute the statistics from saved draws (draws on the raw scale are
  # converted to the coordinates of Derived())
  saved <- readRDS(rdsFile)
  stopifnot(identical(saved$jobs, jobs))
  draws <- lapply(saved$draws, function(d) {
    lambda <- exp(d$logLambda)
    mu <- d$turnover * lambda
    psi <- d$sampling * mu / (1 - d$sampling)
    Derived(lambda, mu, psi)
  })
  cat("Statistics from saved draws:", rdsFile, "\n")
} else {
  t0 <- proc.time()[["elapsed"]]
  draws <- parallel::mclapply(jobs, function(j) {
    tryCatch(if (j$kind == "oracle") Oracle(j$beta, j$seed) else Cet(j$beta, j$seed, j$kind == "control"),
             error = function(e) structure(list(), error = conditionMessage(e)))
  }, mc.cores = cores, mc.preschedule = FALSE)
  cat(sprintf("%d chains in %.1f min\n", length(jobs), (proc.time()[["elapsed"]] - t0) / 60))
  errs <- vapply(draws, function(d) if (is.null(attr(d, "error"))) "" else attr(d, "error"), "")
  if (any(nzchar(errs))) stop("chains failed: ", paste(unique(errs[nzchar(errs)]), collapse = "; "))
  saveRDS(list(jobs = jobs, draws = draws), rdsFile)
}

# Mean and quantiles with MCSE, pooling chains (MCSE from per-chain ESS)
Stat <- function(chains, col, what) {
  vals <- lapply(chains, `[[`, col)
  x <- unlist(vals)
  if (what == "mean") {
    est <- mean(x)
    ess <- sum(vapply(vals, Ess, 0))
    se <- stats::sd(x) / sqrt(ess)
  } else {
    q <- as.numeric(what)
    est <- stats::quantile(x, q, names = FALSE)
    ind <- lapply(vals, function(v) as.numeric(v <= est))
    ess <- sum(vapply(ind, Ess, 0))
    # delta method: se(quantile) = sqrt(q (1 - q) / ess) / density at the quantile
    dens <- stats::density(x, n = 512)
    f <- stats::approx(dens$x, dens$y, est)$y
    se <- sqrt(q * (1 - q) / ess) / f
  }
  c(est = est, se = se, ess = ess)
}
# On smooth scales (log d, logit turnover, logit sampling proportion,
# log lambda): the delta-method MCSE of a quantile needs the density at the
# quantile, which a kernel estimate gives badly on the raw sampling
# proportion (mean 9e-7, 0.9 quantile 1.6e-6, long right tail); a first
# version on the raw scale reported an oracle SE of 0.2% for that quantile
cols <- c("logD", "logitTurnover", "logitSampling", "logLambda")
stats <- rbind(data.frame(col = cols, what = "mean"),
               data.frame(col = rep(cols, each = 3), what = rep(c("0.1", "0.5", "0.9"), 4)))
kinds <- vapply(jobs, `[[`, "", "kind")
betas <- vapply(jobs, `[[`, 0, "beta")
rows <- list()
for (b in c(0, 1)) for (k in c("cet", "control")) {
  o <- draws[kinds == "oracle" & betas == b]
  x <- draws[kinds == k & betas == b]
  for (i in seq_len(nrow(stats))) {
    so <- Stat(o, stats$col[i], stats$what[i])
    sx <- Stat(x, stats$col[i], stats$what[i])
    rows[[length(rows) + 1L]] <- data.frame(
      beta = b, sampler = k, statistic = paste(stats$what[i], stats$col[i]),
      oracle = so[["est"]], oracleSE = so[["se"]], cet = sx[["est"]], cetSE = sx[["se"]],
      z = (sx[["est"]] - so[["est"]]) / sqrt(sx[["se"]]^2 + so[["se"]]^2),
      essOracle = so[["ess"]], essCet = sx[["ess"]])
  }
}
res <- do.call(rbind, rows)
print(res, digits = 4, row.names = FALSE)
utils::write.csv(res, file.path(outDir, "prior-check.csv"), row.names = FALSE)
acc <- unlist(lapply(draws[kinds == "oracle"], function(d) attr(d, "acceptance")))
if (length(acc)) cat(sprintf("oracle acceptance %.3f-%.3f\n", min(acc), max(acc)))
fail <- character()
for (b in c(0, 1)) {
  zc <- res$z[res$beta == b & res$sampler == "cet"]
  zx <- res$z[res$beta == b & res$sampler == "control"]
  cat(sprintf("beta %d: max |z| CET %.2f, control %.2f\n", b, max(abs(zc)), max(abs(zx))))
  if (any(abs(zc) >= 3.5)) fail <- c(fail, sprintf("beta %d: CET differs from the oracle", b))
  if (!any(abs(zx) >= 3.5)) fail <- c(fail, sprintf("beta %d: the control is not rejected", b))
}
if (length(fail)) stop("check-prior.R FAILED: ", paste(fail, collapse = "; "))
cat("check-prior.R: all checks passed\n")
