# Rung 4b localisation (#20): one free parameter at a time, against the exact
# one-dimensional posterior from TreeLikelihood(character = ) on a grid.
# Run from the package root:
#   Rscript data-raw/rung4b-one-parameter.R [nIter] [parameters] [nDA] [character]
# e.g. "12000 theta 8 no" frees theta only, 8 DA chains, without the character;
# a fifth argument "joint" adds joint parameter-path moves; a seventh argument
# "shift" adds one shift move per iteration (#19), which compensates the free
# parameter (`rate`, `lambda`, `mu`, `psi` or `theta`) for a shift of every bin.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args)) as.integer(args[[1]]) else 3000L
free <- if (length(args) > 1) strsplit(args[[2]], ",")[[1]] else c("lambda", "psi", "theta")
nDA <- if (length(args) > 2) as.integer(args[[3]]) else 4L
useCharacter <- !(length(args) > 3 && args[[4]] == "no")
useJoint <- length(args) > 4 && args[[5]] == "joint"
# Sixth argument "rung4-K26": the rung-4 tree and model at the default grid
# (no character; #17)
source4 <- length(args) > 5 && args[[6]] == "rung4-K26"
useShift <- length(args) > 6 && args[[7]] == "shift"
tag <- if (length(args)) paste0("-", paste(free, collapse = "-"), "-", nIter,
                                if (useCharacter) "" else "-nochar",
                                if (useJoint) "-joint" else "",
                                if (source4) "-rung4-K26" else "",
                                if (useShift) "-shift" else "") else ""
burnin <- 200L
r4b <- readRDS(if (source4) "data-raw/rung4-results/rung4-K26.rds" else "data-raw/rung5-results/rung4b.rds")
if (source4) useCharacter <- FALSE
tt <- r4b$tt; model <- r4b$model; truth <- r4b$truth
character <- if (useCharacter) r4b$character else NULL
priors <- c(TempoPriors(), rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
gridFor <- function(p) {
  m <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.1), theta = log(0.1),
         rate = log(0.3))[[p]]
  s <- c(lambda = 0.5, mu = 0.5, psi = 0.5, theta = 0.75, rate = 0.5)[[p]]
  exp(seq(m - 3.5 * s, m + 3.5 * s, length.out = 161))
}
jobs <- rbind(
  data.frame(kind = "grid", p = rep(free, each = 161), i = rep(1:161, length(free))),
  data.frame(kind = "DA", p = rep(free, each = nDA), i = rep(seq_len(nDA), length(free))),
  data.frame(kind = "marginal", p = free, i = 1))
cl <- makeCluster(min(15, nrow(jobs)))
invisible(clusterCall(cl, function() suppressMessages(pkgload::load_all(".", quiet = TRUE))))
clusterExport(cl, c("tt", "model", "truth", "character", "priors", "gridFor", "jobs",
                    "nIter", "burnin", "useJoint", "useShift"))
t0 <- Sys.time()
# Chains first (long), grid points fill the remaining workers
order <- c(which(jobs$kind != "grid"), which(jobs$kind == "grid"))
res <- clusterApplyLB(cl, order, function(j) {
  job <- jobs[j, ]
  set.seed(7000 + j)
  if (job$kind == "grid") {
    par <- truth; par[[job$p]] <- gridFor(job$p)[job$i]
    return(TreeLikelihood(tt, par, model, character = character))
  }
  if (job$kind == "DA") {
    return(TempoMCMC(tt, model, truth, priors[job$p], nIter = nIter, burnin = burnin,
                     nSlice = 20, moves = c(tree = 1, node = 2, tip = 2, branch = 2, joint = as.numeric(useJoint),
                                            shift = as.numeric(useShift)),
                     characters = if (!is.null(character)) matrix(character),
                     minAcceptance = 0)$samples[[job$p]])
  }
  MarginalMCMC(tt, model, truth, priors[job$p], nIter = nIter, burnin = burnin,
               character = character)[[job$p]]
})
stopCluster(cl)
res <- res[order(order)]
runtime <- difftime(Sys.time(), t0, units = "mins")
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
# Batch-means standard error of a quantile: SE of the indicator mean at the
# exact quantile, divided by the exact density there
BatchVar <- function(x, nBatch = 30) {
  m <- floor(length(x) / nBatch) * nBatch
  stats::var(colMeans(matrix(x[seq_len(m)], ncol = nBatch))) / nBatch
}
out <- do.call(rbind, lapply(free, function(p) {
  g <- gridFor(p)
  ll <- unlist(res[jobs$kind == "grid" & jobs$p == p])
  lp <- ll + priors[[p]](g)
  # Posterior density on log scale (uniform spacing), trapezoid CDF
  dens <- exp(lp - max(lp)) * g
  cdf <- c(0, cumsum((dens[-1] + dens[-length(dens)]) / 2)); cdf <- cdf / max(cdf)
  lx <- log(g)
  do.call(rbind, lapply(c("DA", "marginal"), function(kind) {
    chains <- res[jobs$kind == kind & jobs$p == p]
    do.call(rbind, lapply(probs, function(pr) {
      qx <- stats::approx(cdf, lx, pr)$y
      f <- stats::approx(lx, dens / sum((dens[-1] + dens[-length(dens)]) / 2 * diff(lx)), qx)$y
      ind <- lapply(chains, function(ch) as.numeric(log(ch) <= qx))
      seF <- sqrt(sum(vapply(ind, BatchVar, 0) * vapply(ind, length, 0)^2)) /
        sum(vapply(ind, length, 0))
      est <- stats::quantile(log(unlist(chains)), pr, names = FALSE)
      data.frame(parameter = p, kind = kind, p = pr, exact = exp(qx), mcmc = exp(est),
                 mcse = exp(qx) * seF / f, z = (est - qx) / (seF / f))
    }))
  }))
}))
saveRDS(list(jobs = jobs, res = res, out = out, useCharacter = useCharacter, nIter = nIter, burnin = burnin,
             useShift = useShift, runtime = runtime), paste0("data-raw/rung5-results/rung4b-one-parameter", tag, ".rds"))
print(out, digits = 3)
cat("runtime", format(runtime), "\n")
