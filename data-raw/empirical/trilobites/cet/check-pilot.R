# Trilobite CET pipeline: summary of a pilot run (smoke.sh with SMOKE_* set,
# one chain per cell), for the cost projection.
#
#   Rscript data-raw/empirical/trilobites/cet/check-pilot.R <pilot dir> [out dir]
#
# Per cell: seconds per iteration after burn-in, and per quantity the MkPrime
# ESS (common.R, Ess()) of the samples after the first segment, its rate per
# 1000 iterations, and the mean and 95% range.  Writes pilot-convergence.csv
# and pilot-timing.csv to [out dir] (default <pilot dir>).
suppressMessages(library(EvoTempo))
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
pilot <- args[[1]]
outDir <- if (length(args) > 1) args[[2]] else pilot
dirs <- sort(list.dirs(pilot, recursive = FALSE))
dirs <- dirs[grepl("_run[0-9]+$", basename(dirs))]
quantities <- c("lambda", "mu", "psi", "theta", "s", "rate", "logPost", "charLL", "meanX",
                "treeLength", "rootAge", "originAge", "nAncestor", "clock")
conv <- list()
timing <- list()
for (d in dirs) {
  smp <- do.call(rbind, lapply(sort(Sys.glob(file.path(d, "seg*_samples.csv"))), read.csv))
  acc <- read.csv(file.path(d, "acceptance.csv"))
  seg <- unique(acc[, c("segment", "secs", "burnin")])
  iterPerSeg <- max(smp$iteration[smp$segment == 1])
  timing[[d]] <- data.frame(chain = basename(d), segment = seg$segment, secs = seg$secs,
                            secsPerIteration = seg$secs / (iterPerSeg + seg$burnin))
  keep <- smp$segment > 1
  nIter <- diff(range(smp$iteration[keep])) + smp$iteration[2] - smp$iteration[1]
  for (q in intersect(quantities, names(smp))) {
    x <- smp[[q]][keep]
    if (q %in% c("rootAge", "originAge")) x <- x + windowEnd
    e <- Ess(x)
    conv[[paste(d, q)]] <- data.frame(chain = basename(d), quantity = q, iterations = nIter,
                                      ess = e, essPer1000 = 1000 * e / nIter, mean = mean(x),
                                      q025 = unname(stats::quantile(x, 0.025)),
                                      q975 = unname(stats::quantile(x, 0.975)))
  }
}
conv <- do.call(rbind, conv)
timing <- do.call(rbind, timing)
print(timing, row.names = FALSE, digits = 4)
print(conv, row.names = FALSE, digits = 4)
slow <- do.call(rbind, lapply(split(conv, conv$chain), function(x) x[which.min(x$ess), ]))
cat("Slowest quantity per chain:\n")
print(slow[, c("chain", "quantity", "ess", "essPer1000")], row.names = FALSE, digits = 3)
write.csv(conv, file.path(outDir, "pilot-convergence.csv"), row.names = FALSE)
write.csv(timing, file.path(outDir, "pilot-timing.csv"), row.names = FALSE)
