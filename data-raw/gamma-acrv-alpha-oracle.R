# Validation 2 of #97: the alpha move of TempoMCMC() against the exact
# one-dimensional posterior of the Gamma shape on a fixed tree with fixed
# paths.
#
#   Rscript data-raw/gamma-acrv-alpha-oracle.R <data dir> <gamma> <out prefix> [nChain] [nIter]
#
# <data dir> is written by `gamma-acrv-revbayes.R prepare` (fixed time tree
# with fossils and a sampled ancestor; 74 unordered and ordered characters
# with missing cells and ambiguity sets), scored with coding = "variable" and
# four Gamma categories.  Each chain holds every parameter, the tree and the
# paths fixed (no moves, `priors = list()`), so only alpha moves, with prior
# log alpha ~ N(0, 1).  The target is then
#   p(alpha | data) propto prior(alpha) L(alpha), L = MkLikelihood() on the
# chain's own (fixed) tempo-weighted lengths, which with gamma = 1 depend on
# the chain's initial path draw.  The exact posterior CDF F_chain is computed
# on a grid of log alpha (trapezoid rule, 2001 points over +-8 posterior sd)
# for each chain; the probability integral transforms u = F_chain(alpha_i)
# are uniform under the correct target.  Pooled over chains, it reports the
# batch-means z-scores of mean(u) - 1/2 and of P(u < q) - q for q in
# {0.05, 0.25, 0.5, 0.75, 0.95}, and a KS test on thinned u.
suppressMessages({
  library(EvoTempo)
  library(ape)
})
ns <- asNamespace("EvoTempo")
args <- commandArgs(TRUE)
dir <- args[[1]]
gamma <- as.numeric(args[[2]])
out <- args[[3]]
nChain <- if (length(args) > 3) as.integer(args[[4]]) else 16L
nIter <- if (length(args) > 4) as.integer(args[[5]]) else 20000L
burnin <- 1000L

tr <- read.tree(file.path(dir, "tree.nwk"))
depth <- node.depth.edgelength(tr)
age <- max(depth) - depth
age[abs(age) < 1e-9] <- 0
len <- numeric(Ntip(tr) + tr$Nnode)
len[tr$edge[, 2]] <- tr$edge.length
tip <- seq_len(Ntip(tr))
tt <- TempoTree(tr, ages = age,
                tipType = ifelse(age[tip] == 0, "extant",
                                 ifelse(len[tip] < 1e-12, "ancestor", "fossil")))
raw <- read.csv(file.path(dir, "cells.csv"), colClasses = "character",
                check.names = FALSE, na.strings = character(0))
cells <- as.matrix(raw[match(tt$tipLabel, raw$taxon), -1])
parts <- read.csv(file.path(dir, "partitions.csv"))
partOf <- rep(seq_len(nrow(parts)), parts$nChar)
mk <- MkCharacters(cells, nState = parts$nState[partOf], ordered = parts$ordered[partOf],
                   coding = "variable", nRateCat = 4)
model <- DCETModel(nBin = 11, zRange = c(-2, 2), gamma = gamma)
par <- list(lambda = 3, mu = 1.5, psi = 0.8, theta = 0.5, s = 0.6, rate = 0.6, alpha = 1)
aPrior <- function(a) stats::dlnorm(a, 0, 1, log = TRUE)

ExactCdf <- function(lengths) {
  LogPost <- function(u) {
    vapply(u, function(ui) {
      aPrior(exp(ui)) + ui + MkLikelihood(tt, mk, lengths, alpha = exp(ui))
    }, 0)
  }
  # Mode and curvature, then a grid over +-8 sd
  opt <- stats::optimize(function(u) -LogPost(u), c(-7, 7))
  h <- 1e-3
  curv <- -(LogPost(opt$minimum + h) - 2 * LogPost(opt$minimum) +
              LogPost(opt$minimum - h)) / h^2
  sd <- 1 / sqrt(max(curv, 1e-6))
  u <- seq(opt$minimum - 8 * sd, opt$minimum + 8 * sd, length.out = 2001)
  lp <- LogPost(u)
  d <- exp(lp - max(lp))
  cdf <- c(0, cumsum((d[-1] + d[-length(d)]) / 2 * diff(u)))
  cdf <- cdf / cdf[length(cdf)]
  list(u = u, cdf = cdf, mean = sum(u * d) / sum(d), sd = sd,
       tailMass = d[c(1, length(d))])
}

BatchVar <- function(x, nBatch = 25) {
  m <- floor(length(x) / nBatch) * nBatch
  stats::var(colMeans(matrix(x[seq_len(m)], ncol = nBatch))) / nBatch
}
PooledZ <- function(chains, exact) {
  n <- vapply(chains, length, 0)
  est <- sum(vapply(chains, sum, 0)) / sum(n)
  se <- sqrt(sum(vapply(chains, BatchVar, 0) * n^2)) / sum(n)
  c(estimate = est, exact = exact, z = (est - exact) / se)
}

t0 <- Sys.time()
chains <- parallel::mclapply(seq_len(nChain), function(k) {
  set.seed(9700 + k)
  run <- TempoMCMC(tt, model, par, priors = list(), nIter = nIter, burnin = burnin,
                   nSlice = 8, characters = mk, charPriors = list(alpha = aPrior),
                   moves = c(node = 0), minAcceptance = 0)
  st <- run$state
  lengths <- ns$.CharacterLengths(st$tree, st$set$grid$x, st$paths, st$par$rate, gamma)
  ex <- ExactCdf(lengths)
  pit <- stats::approx(ex$u, ex$cdf, log(run$samples$alpha), rule = 2)$y
  list(pit = pit, logAlpha = log(run$samples$alpha), exactMean = ex$mean,
       exactSd = ex$sd, tailMass = ex$tailMass, acceptance = run$acceptance[["alpha"]],
       burnin = run$burninAcceptance[["alpha"]])
}, mc.cores = min(nChain, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))))
pits <- lapply(chains, `[[`, "pit")
rows <- rbind(
  data.frame(statistic = "mean(u)", t(PooledZ(pits, 0.5))),
  do.call(rbind, lapply(c(0.05, 0.25, 0.5, 0.75, 0.95), function(q) {
    data.frame(statistic = sprintf("P(u < %.2f)", q),
               t(PooledZ(lapply(pits, function(u) as.numeric(u < q)), q)))
  })),
  data.frame(statistic = "mean(log alpha - exact mean)",
             t(PooledZ(lapply(chains, function(ch) ch$logAlpha - ch$exactMean), 0))))
rows$gamma <- gamma
ks <- stats::ks.test(unlist(lapply(pits, function(u) u[seq(1, length(u), by = 50)])),
                     "punif")
summ <- data.frame(gamma = gamma, nChain = nChain, nIter = nIter,
                   exactMeanLogAlpha = vapply(chains, `[[`, 0, "exactMean"),
                   exactSdLogAlpha = vapply(chains, `[[`, 0, "exactSd"),
                   maxTailDensity = vapply(chains, function(ch) max(ch$tailMass), 0),
                   acceptance = vapply(chains, `[[`, 0, "acceptance"))
print(rows, digits = 4, row.names = FALSE)
print(summ, digits = 4, row.names = FALSE)
cat(sprintf("KS on thinned u (every 50th, n = %d): D = %.4f, p = %.3f\n",
            length(unlist(lapply(pits, function(u) u[seq(1, length(u), by = 50)]))),
            ks$statistic, ks$p.value))
cat(sprintf("max |z| = %.2f; runtime %.1f min\n", max(abs(rows$z)),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(rows, paste0(out, "-z.csv"), row.names = FALSE)
write.csv(summ, paste0(out, "-chains.csv"), row.names = FALSE)
writeLines(sprintf("gamma %s: KS D %.4f p %.4f; max|z| %.3f", gamma, ks$statistic,
                   ks$p.value, max(abs(rows$z))), paste0(out, "-ks.txt"))
