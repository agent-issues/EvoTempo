# Hymenoptera CET pipeline: check of the lineage-averaged FBD-rate prior with
# s free (common.R, CET_FBD_SCALE = "lineage", CET_S_MAX; fix C of README,
# "Cell (1, 0)"), as check-prior.R checks the default prior with s fixed.
#
#   Rscript data-raw/empirical/hymenoptera/cet/check-prior-scale.R <out dir> [cores]
#
# Target: check-prior.R's 11-tip subtree, no characters, beta = 1, gamma = 0,
# the production grid (26 bins, z in [-5, 5]), theta fixed at 1/120;
# lambda, mu, psi and s free; prior: the comparator's (d ~ Exp(100), turnover
# and sampling proportion ~ U(0, 1)) on the lineage-averaged rates
# k lambda, k mu (k = LineageMean(s, 1)), and s ~ Exp(mean 0.587405) on
# (0, 2].
# - Oracle: random-walk Metropolis on v = (log dbar, logit turnover, logit
#   sampling proportion, log s), dbar = k (lambda - mu), sampling proportion
#   psi / (k mu + psi); prior written directly in v (the product of the four
#   priors and the log / logit Jacobians); likelihood TreeLikelihood() (paths
#   integrated out; atol 1e-60, #90).  It does not use FbdLogPrior(), so it
#   checks the 2 log k Jacobian and the joint prior in s.
# - CET: TempoMCMC() with CetPriors(1) (joint prior on lambda, mu, psi, s;
#   parameter and joint sweeps, shift move, FBD block moves with the general
#   target).
# - Negative control: CET with the Jacobian 2 log k left out.
# Compared: mean and 0.1, 0.5, 0.9 quantiles of log dbar, logit turnover,
# logit sampling proportion, log s and log lambda; z as check-prior.R.  Pass:
# every |z| < 3.5 for CET; the control must have some |z| >= 3.5.  Writes
# <out dir>/prior-check-scale.csv and -draws.rds.
Sys.setenv(CET_FBD_SCALE = "lineage", CET_S_MAX = "2", CET_NBIN = "26", CET_ZMAX = "5")
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
nIterCet <- as.integer(Sys.getenv("CET_PC_NITER", "30000"))
nIterOracle <- as.integer(Sys.getenv("CET_PC_ORACLE", "12000"))
burnFrac <- 0.2
sMax <- 2

# ---- Target (as check-prior.R) ---------------------------------------------------
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
tree <- TempoTree(sub, ages = age, tipType = ifelse(extant, "extant", "fossil"))
model <- DCETModel(nBin = 26, zRange = c(-5, 5), beta = 1, gamma = 0, rho = 5e-4,
                   condition = "extant", start = "root")
fixed <- list(theta = 1 / 120)

FromV <- function(v) {
  s <- exp(v[4])
  k <- LineageMean(s, 1)
  dbar <- exp(v[1]); r <- stats::plogis(v[2]); sp <- stats::plogis(v[3])
  lbar <- dbar / (1 - r); mbar <- r * lbar
  list(lambda = lbar / k, mu = mbar / k, psi = sp * mbar / (1 - sp), s = s)
}
Derived <- function(lambda, mu, psi, s) {
  k <- vapply(s, LineageMean, 0, beta = 1)
  data.frame(logDbar = log(k * (lambda - mu)), logitTurnover = log(mu) - log(lambda - mu),
             logitSampling = log(psi) - log(k * mu), logS = log(s), logLambda = log(lambda))
}
Oracle <- function(seed) {
  set.seed(seed)
  LogPost <- function(v) {
    s <- exp(v[4])
    if (s > sMax) return(-Inf)
    dbar <- exp(v[1]); r <- stats::plogis(v[2]); sp <- stats::plogis(v[3])
    prior <- log(100) - 100 * dbar + v[1] + log(r) + log1p(-r) + log(sp) + log1p(-sp) +
      stats::dexp(s, 1 / 0.587405, log = TRUE) - stats::pexp(sMax, 1 / 0.587405, log.p = TRUE) + v[4]
    ll <- tryCatch(TreeLikelihood(tree, c(FromV(v), fixed), model, atol = 1e-60),
                   error = function(e) -Inf)
    prior + ll
  }
  v <- c(log(0.01), stats::qlogis(0.9), stats::qlogis(0.1), log(0.5))
  lp <- LogPost(v)
  S <- diag(c(0.5, 0.5, 0.5, 0.3)^2)
  out <- matrix(NA_real_, nIterOracle, 4)
  acc <- 0
  for (i in seq_len(nIterOracle)) {
    prop <- v + drop(t(chol(S)) %*% stats::rnorm(4))
    lq <- LogPost(prop)
    if (log(stats::runif(1)) < lq - lp) {
      v <- prop; lp <- lq; acc <- acc + 1
    }
    out[i, ] <- v
    if (i <= burnFrac * nIterOracle && i %% 500 == 0 && i >= 1000) {
      S <- 2.38^2 / 4 * stats::cov(out[(i %/% 2):i, ]) + diag(1e-6, 4)
    }
  }
  keepI <- seq_len(nIterOracle) > burnFrac * nIterOracle
  p <- lapply(which(keepI), function(i) FromV(out[i, ]))
  res <- Derived(vapply(p, `[[`, 0, "lambda"), vapply(p, `[[`, 0, "mu"), vapply(p, `[[`, 0, "psi"),
                 vapply(p, `[[`, 0, "s"))
  attr(res, "acceptance") <- acc / nIterOracle
  res
}
Cet <- function(seed, control = FALSE) {
  set.seed(seed)
  if (control) {
    # Negative control: the Jacobian 2 log k left out (CetPriors()'s joint
    # prior looks FbdLogPrior() up in the global environment)
    orig <- FbdLogPrior
    assign("FbdLogPrior", function(par, beta = 0) {
      k <- if (beta != 0) LineageMean(par[["s"]], beta) else 1
      orig(par, beta) - 2 * log(k)
    }, envir = globalenv())
  }
  pri <- CetPriors(1)[c(fbdNames, "s")]
  start <- c(list(lambda = 0.1, mu = 0.09, psi = 0.01, s = 0.5), fixed)
  run <- TempoMCMC(tree, model, start, pri, nIter = nIterCet, burnin = 1000, thin = 5)
  smp <- run$samples
  smp <- smp[seq_len(nrow(smp)) > burnFrac * nrow(smp), ]
  res <- Derived(smp$lambda, smp$mu, smp$psi, smp$s)
  attr(res, "acceptance") <- run$acceptance
  res
}
jobs <- c(lapply(seq_len(nChain), function(k) list(kind = "oracle", seed = 700 + k)),
          lapply(seq_len(nChain), function(k) list(kind = "cet", seed = 800 + k)),
          lapply(seq_len(nChain), function(k) list(kind = "control", seed = 900 + k)))
rdsFile <- file.path(outDir, "prior-check-scale-draws.rds")
if (nzchar(Sys.getenv("CET_PC_REUSE")) && file.exists(rdsFile)) {
  draws <- readRDS(rdsFile)$draws
} else {
  t0 <- proc.time()[["elapsed"]]
  # Each chain's draws are saved as it finishes (<out dir>/chain-<kind>-<seed>.rds)
  # and reused by a rerun, so a job that hits its time limit loses only the
  # chains still running
  draws <- parallel::mclapply(jobs, function(j) {
    f <- file.path(outDir, sprintf("chain-%s-%d.rds", j$kind, j$seed))
    if (file.exists(f)) return(readRDS(f))
    t1 <- proc.time()[["elapsed"]]
    d <- tryCatch(if (j$kind == "oracle") Oracle(j$seed) else Cet(j$seed, j$kind == "control"),
                  error = function(e) structure(list(), error = conditionMessage(e)))
    if (is.null(attr(d, "error"))) saveRDS(d, f)
    message(sprintf("%s %s seed %d: %.1f min", format(Sys.time()), j$kind, j$seed,
                    (proc.time()[["elapsed"]] - t1) / 60))
    d
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
cols <- c("logDbar", "logitTurnover", "logitSampling", "logS", "logLambda")
stats <- rbind(data.frame(col = cols, what = "mean"),
               data.frame(col = rep(cols, each = 3), what = rep(c("0.1", "0.5", "0.9"), length(cols))))
kinds <- vapply(jobs, `[[`, "", "kind")
rows <- list()
for (k in c("cet", "control")) {
  o <- draws[kinds == "oracle"]
  x <- draws[kinds == k]
  for (i in seq_len(nrow(stats))) {
    so <- Stat(o, stats$col[i], stats$what[i])
    sx <- Stat(x, stats$col[i], stats$what[i])
    rows[[length(rows) + 1L]] <- data.frame(
      sampler = k, statistic = paste(stats$what[i], stats$col[i]),
      oracle = so[["est"]], oracleSE = so[["se"]], cet = sx[["est"]], cetSE = sx[["se"]],
      z = (sx[["est"]] - so[["est"]]) / sqrt(sx[["se"]]^2 + so[["se"]]^2),
      essOracle = so[["ess"]], essCet = sx[["ess"]])
  }
}
res <- do.call(rbind, rows)
options(width = 200)
print(res, digits = 4, row.names = FALSE)
utils::write.csv(res, file.path(outDir, "prior-check-scale.csv"), row.names = FALSE)
acc <- unlist(lapply(draws[kinds == "oracle"], function(d) attr(d, "acceptance")))
cat(sprintf("oracle acceptance %.3f-%.3f\n", min(acc), max(acc)))
zc <- res$z[res$sampler == "cet"]
zx <- res$z[res$sampler == "control"]
cat(sprintf("max |z| CET %.2f, control %.2f\n", max(abs(zc)), max(abs(zx))))
fail <- c(if (any(abs(zc) >= 3.5)) "CET differs from the oracle",
          if (!any(abs(zx) >= 3.5)) "the control is not rejected")
if (length(fail)) stop("check-prior-scale.R FAILED: ", paste(fail, collapse = "; "))
cat("check-prior-scale.R: all checks passed\n")
