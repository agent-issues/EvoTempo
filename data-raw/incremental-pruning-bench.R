# Benchmark of incremental character pruning (#26, item 3; dev/plans/incremental-pruning.md).
# Chains are those of data-raw/mixing-schedule.R (`new` = TempoMoves() schedule, fixed
# trees of dev/mixing/, cell beta = gamma = 1, K = 16, 30 / 100 / 300 nested characters).
# Each task ran the same seeded chain under the package before this change (`base`) and
# after it (`inc`), one after the other on the same node; the samples must be identical,
# so the ESS is shared and the gain in ESS per CPU hour is the ratio of CPU times.
# Usage, from the package root:
#   Rscript data-raw/incremental-pruning-bench.R <baseDir> <incDir> <convDir> [outDir]
# where each *Dir holds the `chains/` folder written by mixing-schedule.R, and convDir the
# longer large-tree chains (20000 iterations) run under `inc`.
args <- commandArgs(trailingOnly = TRUE)
baseDir <- args[[1]]
incDir <- args[[2]]
convDir <- args[[3]]
outDir <- if (length(args) > 3) args[[4]] else "data-raw/incremental-pruning-results"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
logged <- c("theta", "s", "rate", "lambda", "clock")

files <- list.files(file.path(incDir, "chains"), "\\.rds$")
pairs <- do.call(rbind, lapply(files, function(f) {
  b <- readRDS(file.path(baseDir, "chains", f))
  n <- readRDS(file.path(incDir, "chains", f))
  data.frame(tree = c("small", "large")[n$tree], nChar = n$nChar, rep = n$rep,
             nIter = n$nIter, burnin = n$burnin,
             identical = identical(b$samples, n$samples) &&
               identical(b$acceptance, n$acceptance),
             cpuBase = b$cpu, cpuInc = n$cpu, speedup = b$cpu / n$cpu,
             secPerIterBase = b$cpu / (b$nIter + b$burnin),
             secPerIterInc = n$cpu / (n$nIter + n$burnin),
             t(n$ess), shaBase = b$sha, shaInc = n$sha)
}))
pairs <- pairs[order(match(pairs$tree, c("small", "large")), pairs$nChar, pairs$rep), ]
write.csv(pairs, file.path(outDir, "bench-chains.csv"), row.names = FALSE)

quantities <- c("theta", "s", "rate", "lambda", "rootBin", "clock", "logPost")
cells <- do.call(rbind, lapply(split(pairs, interaction(pairs$tree, pairs$nChar, drop = TRUE)),
                               function(d) {
  hoursBase <- d$cpuBase * d$nIter / (d$nIter + d$burnin) / 3600
  hoursInc <- d$cpuInc * d$nIter / (d$nIter + d$burnin) / 3600
  out <- data.frame(tree = d$tree[1], nChar = d$nChar[1], nRep = nrow(d),
                    allIdentical = all(d$identical),
                    secPerIterBase = mean(d$secPerIterBase),
                    secPerIterInc = mean(d$secPerIterInc),
                    gain = mean(hoursBase) / mean(hoursInc))
  for (q in quantities) {
    out[[paste0(q, "PerHourBase")]] <- mean(d[[q]]) / mean(hoursBase)
    out[[paste0(q, "PerHourInc")]] <- mean(d[[q]]) / mean(hoursInc)
  }
  out
}))
cells <- cells[order(match(cells$tree, c("small", "large")), cells$nChar), ]
write.csv(cells, file.path(outDir, "bench-gains.csv"), row.names = FALSE)

# Split R-hat over chains (each split in halves), as data-raw/mixing-schedule.R
SplitRhat <- function(xs) {
  halves <- unlist(lapply(xs, function(x) {
    n <- length(x) %/% 2
    list(x[seq_len(n)], x[n + seq_len(n)])
  }), recursive = FALSE)
  n <- length(halves[[1]])
  m <- vapply(halves, mean, 0)
  w <- mean(vapply(halves, stats::var, 0))
  b <- n * stats::var(m)
  sqrt(((n - 1) / n * w + b / n) / w)
}
Ess <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4 || stats::var(x) == 0) return(NA_real_)
  ac <- stats::acf(x, lag.max = n - 1, plot = FALSE, demean = TRUE)$acf[, 1, 1]
  nPair <- length(ac) %/% 2
  pairSum <- ac[2 * seq_len(nPair) - 1] + ac[2 * seq_len(nPair)]
  stop1 <- which(pairSum <= 0)[1]
  if (!is.na(stop1)) pairSum <- pairSum[seq_len(stop1 - 1)]
  tau <- if (length(pairSum)) -1 + 2 * sum(cummin(pairSum)) else 1
  min(n, n / max(tau, 1e-8))
}
Get <- function(r, q) if (q %in% logged) log(r$samples[[q]]) else r$samples[[q]]

convFiles <- list.files(file.path(convDir, "chains"), "\\.rds$", full.names = TRUE)
conv <- lapply(convFiles, readRDS)
byChar <- split(conv, vapply(conv, `[[`, 0, "nChar"))
# R-hat on the first `frac` of each chain's stored draws: 0.25 of 20000 iterations is
# 5000, the length of the mixing-schedule chains
rhat <- do.call(rbind, lapply(names(byChar), function(nc) {
  rs <- byChar[[nc]]
  do.call(rbind, lapply(c(0.25, 0.5, 1), function(frac) {
    keep <- function(x) x[seq_len(floor(length(x) * frac))]
    nIter <- rs[[1]]$nIter * frac
    data.frame(nChar = as.integer(nc), nChain = length(rs), nIter = nIter,
               cpuHoursPerChain = mean(vapply(rs, function(r) {
                 r$cpu * nIter / (r$nIter + r$burnin) / 3600
               }, 0)),
               t(vapply(quantities, function(q) {
                 SplitRhat(lapply(rs, function(r) keep(Get(r, q))))
               }, 0)))
  }))
}))
write.csv(rhat, file.path(outDir, "conv-rhat.csv"), row.names = FALSE)
convEss <- do.call(rbind, lapply(conv, function(r) {
  data.frame(nChar = r$nChar, rep = r$rep, nIter = r$nIter, cpu = r$cpu, t(r$ess))
}))
write.csv(convEss, file.path(outDir, "conv-chains.csv"), row.names = FALSE)
print(pairs[, c("tree", "nChar", "rep", "identical", "cpuBase", "cpuInc", "speedup")],
      digits = 3, row.names = FALSE)
print(cells[, c("tree", "nChar", "allIdentical", "secPerIterBase", "secPerIterInc", "gain")],
      digits = 3, row.names = FALSE)
print(rhat, digits = 3, row.names = FALSE)
