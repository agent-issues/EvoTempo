# Summary of the full-size rung 6 SBC (data-raw/rung6-sbc-full.R, issue #21).
# Reads <dir>/replicates/rep*.rds and writes, in <dir>:
#   rung6-full-ranks.csv  one row per replicate and quantity (rank, truth, ESS)
#   rung6-full-runs.csv   one row per replicate (size, iterations, run time, acceptance)
#   rung6-full-tests.csv  chi-square (10 bins) and KS uniformity per cell and
#                         quantity, with Holm and Benjamini-Hochberg adjustment
#                         over all cell x quantity tests of each kind
#   rung6-full-tests-ess50.csv  the same, on replicates with ESS >= 50 throughout
#   rung6-full-ess.csv    ESS per cell and quantity
#   rung6-full-cells.csv  per-cell verdict, ESS and core-hours
#   rung6-full-ecdf.csv   ECDF differences, and rung6-full-ecdf-band.csv the
#                         simultaneous 95% band (Sailynoja et al. 2022)
# Usage: Rscript data-raw/rung6-sbc-full-summary.R [dir]
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[[1]] else "data-raw/rung6-full-results"
files <- list.files(file.path(dir, "replicates"), "^rep[0-9]+\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)
quantities <- names(res[[1]]$ranks)
nRank <- res[[1]]$nRank
stopifnot(all(vapply(res, `[[`, 0L, "nRank") == nRank))
cellLabel <- function(b, g) sprintf("beta=%d, gamma=%d", b, g)

ranks <- do.call(rbind, lapply(res, function(r) {
  data.frame(g = r$g, cell = r$cell, beta = r$beta, gamma = r$gamma,
             quantity = quantities, rank = unname(r$ranks),
             truth = unname(r$truth[quantities]), ess = unname(r$ess[quantities]),
             essRanked = unname(r$essRanked[quantities]),
             posteriorMean = unname(r$posteriorMean[quantities]),
             posteriorSd = unname(r$posteriorSd[quantities]))
}))
write.csv(ranks, file.path(dir, "rung6-full-ranks.csv"), row.names = FALSE)

runs <- do.call(rbind, lapply(res, function(r) {
  acc <- colMeans(r$acceptance)
  data.frame(g = r$g, cell = r$cell, beta = r$beta, gamma = r$gamma, seed = r$seed,
             nTip = r$nTip, nFossil = r$nFossil, nAncestor = r$nAncestor,
             tries = r$tries, nIter = r$nIter, nExtend = r$nExtend,
             chainRetries = length(r$chainErrors),
             minEss = min(r$ess, na.rm = TRUE),
             minEssQuantity = names(which.min(r$ess)), runtimeMin = r$runtime,
             t(setNames(acc, paste0("acc.", names(acc)))), sha = r$sha)
}))
write.csv(runs, file.path(dir, "rung6-full-runs.csv"), row.names = FALSE)

nBin <- 10L
Tests <- function(ranks) {
  tests <- do.call(rbind, lapply(split(ranks, list(ranks$cell, ranks$quantity)), function(d) {
    rk <- d$rank[!is.na(d$rank)]
    n <- length(rk)
    counts <- tabulate(rk %/% ((nRank + 1L) / nBin) + 1L, nBin)
    u <- (rk + stats::runif(n)) / (nRank + 1L)
    data.frame(cell = d$cell[1], beta = d$beta[1], gamma = d$gamma[1],
               quantity = d$quantity[1], n = n,
               chisq = suppressWarnings(stats::chisq.test(counts)$statistic),
               pChisq = suppressWarnings(stats::chisq.test(counts)$p.value),
               ksD = stats::ks.test(u, "punif")$statistic,
               pKS = stats::ks.test(u, "punif")$p.value,
               minBin = min(counts), maxBin = max(counts))
  }))
  tests$pChisqHolm <- stats::p.adjust(tests$pChisq, "holm")
  tests$pKSHolm <- stats::p.adjust(tests$pKS, "holm")
  tests$pChisqBH <- stats::p.adjust(tests$pChisq, "BH")
  tests$pKSBH <- stats::p.adjust(tests$pKS, "BH")
  tests$pass <- tests$pChisqHolm >= 0.05 & tests$pKSHolm >= 0.05
  tests[order(tests$cell, match(tests$quantity, quantities)), ]
}
set.seed(6021)
tests <- Tests(ranks)
write.csv(tests, file.path(dir, "rung6-full-tests.csv"), row.names = FALSE)
# Sensitivity: only replicates whose every ranked quantity reached ESS 50
set.seed(6023)
okG <- runs$g[runs$minEss >= 50]
write.csv(Tests(ranks[ranks$g %in% okG, ]), file.path(dir, "rung6-full-tests-ess50.csv"),
          row.names = FALSE)

ess <- do.call(rbind, lapply(split(ranks, list(ranks$cell, ranks$quantity)), function(d) {
  e <- d$ess[!is.na(d$ess)]
  data.frame(cell = d$cell[1], beta = d$beta[1], gamma = d$gamma[1],
             quantity = d$quantity[1], n = length(e), min = min(e),
             q05 = unname(stats::quantile(e, 0.05)), median = stats::median(e),
             fracAtLeast50 = mean(e >= 50),
             medianRanked = stats::median(d$essRanked, na.rm = TRUE))
}))
ess <- ess[order(ess$cell, match(ess$quantity, quantities)), ]
write.csv(ess, file.path(dir, "rung6-full-ess.csv"), row.names = FALSE)

cellsOut <- do.call(rbind, lapply(split(runs, runs$cell), function(d) {
  t <- tests[tests$cell == d$cell[1], ]
  e <- ess[ess$cell == d$cell[1], ]
  data.frame(cell = d$cell[1], label = cellLabel(d$beta[1], d$gamma[1]),
             nReplicate = nrow(d), nPass = sum(t$pass), nQuantity = nrow(t),
             pass = all(t$pass), minPChisqHolm = min(t$pChisqHolm),
             minPKSHolm = min(t$pKSHolm),
             fracReplicatesAllEss50 = mean(d$minEss >= 50),
             minQuantityFracEss50 = min(e$fracAtLeast50),
             medianIter = stats::median(d$nIter), fracExtended = mean(d$nExtend > 0),
             coreHours = sum(d$runtimeMin) / 60, medianTip = stats::median(d$nTip),
             medianFossil = stats::median(d$nFossil))
}))
write.csv(cellsOut, file.path(dir, "rung6-full-cells.csv"), row.names = FALSE)

# ECDF differences of the fractional ranks, evaluated where the discrete
# uniform ECDF is exact: F(j / (L + 1)) = P(rank < j), j = 1..L
ecdf <- do.call(rbind, lapply(split(ranks, list(ranks$cell, ranks$quantity)), function(d) {
  rk <- d$rank[!is.na(d$rank)]
  z <- seq_len(nRank) / (nRank + 1)
  data.frame(cell = d$cell[1], beta = d$beta[1], gamma = d$gamma[1],
             quantity = d$quantity[1], n = length(rk), z = z,
             diff = vapply(seq_len(nRank), function(j) mean(rk < j), 0) - z)
}))
write.csv(ecdf, file.path(dir, "rung6-full-ecdf.csv"), row.names = FALSE)

# Simultaneous 95% band: pointwise binomial level calibrated by simulation so
# that a uniform sample of size n stays inside at every z with probability 0.95
Band <- function(n, nSim = 4000, conf = 0.95) {
  z <- seq_len(nRank) / (nRank + 1)
  sims <- replicate(nSim, {
    rk <- sample.int(nRank + 1L, n, replace = TRUE) - 1L
    cumsum(tabulate(rk + 1L, nRank + 1L))[seq_len(nRank)]
  })
  # Each simulated ECDF's most extreme pointwise binomial tail probability
  lower <- stats::pbinom(sims, n, z)
  upper <- stats::pbinom(sims - 1L, n, z, lower.tail = FALSE)
  extreme <- apply(pmin(lower, upper), 2, min)
  alpha <- stats::quantile(extreme, 1 - conf, type = 1)
  data.frame(n = n, z = z,
             lo = stats::qbinom(alpha, n, z) / n - z,
             hi = stats::qbinom(alpha, n, z, lower.tail = FALSE) / n - z,
             pointwiseAlpha = unname(alpha))
}
set.seed(6022)
band <- do.call(rbind, lapply(sort(unique(ecdf$n)), Band))
write.csv(band, file.path(dir, "rung6-full-ecdf-band.csv"), row.names = FALSE)
print(cellsOut)
