# Hymenoptera CET pipeline: mixing summary of a pilot (one chain per cell,
# run-array.sh with CET_OUT=<pilot dir>), for the cost projection.
#
#   Rscript data-raw/empirical/hymenoptera/cet/pilot-summary.R <pilot dir> [burn fraction]
#
# Per cell: seconds per iteration (acceptance.csv, sampling segments only),
# rank-normalised bulk and tail ESS (MkPrime's internal .EssBulk(),
# .EssTail(); Vehtari et al. 2021) of each monitored quantity after dropping
# the first <burn fraction> (default 0.25) of the samples, and ESS per 1000
# iterations.  One chain per cell, so no R-hat.  Writes
# <pilot dir>/pilot-ess.csv and pilot-speed.csv.
args <- commandArgs(trailingOnly = TRUE)
dir <- args[[1]]
burn <- if (length(args) >= 2) as.numeric(args[[2]]) else 0.25
EssBulk <- get(".EssBulk", asNamespace("MkPrime"))
EssTail <- get(".EssTail", asNamespace("MkPrime"))
cells <- c("b0g0", "b0g1", "b1g0", "b1g1")
qs <- c("lambda", "mu", "psi", "theta", "s", "rate", "logD", "logitTurnover", "logPost",
        "charLL", "rootAge", "holoAge", "hymCrownAge", "treeLength", "meanNodeAge", "nAncestor")
ess <- list(); speed <- list()
for (cl in cells) {
  d <- file.path(dir, paste0(cl, "_run1"))
  smp <- do.call(rbind, lapply(sort(Sys.glob(file.path(d, "seg*_samples.csv"))), utils::read.csv))
  smp$logD <- log(smp$lambda - smp$mu)
  smp$logitTurnover <- log(smp$mu) - smp$logD
  keep <- smp[seq_len(nrow(smp)) > floor(burn * nrow(smp)), ]
  thin <- stats::median(diff(smp$iteration))
  nIt <- nrow(keep) * thin
  acc <- utils::read.csv(file.path(d, "acceptance.csv"))
  segs <- unique(acc[acc$burnin == 0, c("segment", "secs")])
  seg <- readRDS(file.path(d, "checkpoint.rds"))$seg
  speed[[cl]] <- data.frame(cell = cl, segments = max(smp$segment), samples = nrow(smp),
                            iterationsKept = nIt, secsPerIteration = sum(segs$secs) / (nrow(segs) * seg))
  for (q in intersect(qs, names(keep))) {
    x <- keep[[q]]
    b <- tryCatch(EssBulk(x), error = function(e) NA_real_)
    t <- tryCatch(EssTail(x), error = function(e) NA_real_)
    ess[[length(ess) + 1L]] <- data.frame(cell = cl, quantity = q, mean = mean(x), sd = stats::sd(x),
                                          bulkEss = b, tailEss = t,
                                          minEssPer1000 = 1000 * min(b, t, na.rm = TRUE) / nIt)
  }
}
ess <- do.call(rbind, ess)
speed <- do.call(rbind, speed)
utils::write.csv(ess, file.path(dir, "pilot-ess.csv"), row.names = FALSE)
utils::write.csv(speed, file.path(dir, "pilot-speed.csv"), row.names = FALSE)
options(width = 200)
print(speed, row.names = FALSE)
print(ess, digits = 4, row.names = FALSE)
slow <- do.call(rbind, lapply(split(ess, ess$cell), function(x) x[which.min(x$minEssPer1000), ]))
cat("Slowest quantity per cell:\n")
print(slow[, c("cell", "quantity", "bulkEss", "tailEss", "minEssPer1000")], digits = 4, row.names = FALSE)
