# Hymenoptera CET pipeline: diagnosis of the old root ages in the pilot
# (README, "Root ages").
#
#   Rscript data-raw/empirical/hymenoptera/cet/diag-rootage.R <out.csv> <label>=<source> ...
#
# Each <source> is a CET chain directory (seg*_samples.csv) or a RevBayes log
# (*.log, columns root_age and holo_age); several sources may share a label
# (pooled).  Drops the first 25% of each source.  For each label: number of
# samples, bulk ESS of the root age (MkPrime's .EssBulk(), summed over
# sources), the 5, 50 and 95% quantiles of the root and Holometabola ages, and
# the share of Holometabola ages below 302 Ma (impossible under the
# calibration).  Labels ending in "-nocal" are also reported reweighted by the
# calibration density exp(-(t - 302) / 94) / 94 (label "<label>-reweighted"):
# the prior with the calibration, if the run lacked it (RevBayes's
# ignoreAllData() ignores the clamped calibration node).
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[[1]]
spec <- args[-1]
labels <- sub("=.*$", "", spec)
paths <- sub("^[^=]*=", "", spec)
EssBulk <- get(".EssBulk", asNamespace("MkPrime"))
Read <- function(p) {
  x <- if (grepl("\\.log$", p)) {
    y <- utils::read.delim(p)
    data.frame(root = y$root_age, holo = y$holo_age)
  } else {
    y <- do.call(rbind, lapply(sort(Sys.glob(file.path(p, "seg*_samples.csv"))), utils::read.csv))
    data.frame(root = y$rootAge, holo = y$holoAge)
  }
  x[seq_len(nrow(x)) > floor(0.25 * nrow(x)), ]
}
WQ <- function(v, w, q) {
  o <- order(v); v <- v[o]; w <- cumsum(w[o]) / sum(w)
  vapply(q, function(qq) v[which(w >= qq)[1]], 0)
}
rows <- list()
for (lab in unique(labels)) {
  xs <- lapply(paths[labels == lab], Read)
  ess <- sum(vapply(xs, function(x) tryCatch(EssBulk(x$root), error = function(e) NA_real_), 0))
  x <- do.call(rbind, xs)
  Row <- function(name, w) {
    rq <- WQ(x$root, w, c(0.05, 0.5, 0.95)); hq <- WQ(x$holo, w, c(0.05, 0.5, 0.95))
    data.frame(label = name, n = nrow(x), essRoot = ess, kishEss = sum(w)^2 / sum(w^2),
               root05 = rq[1], root50 = rq[2], root95 = rq[3],
               holo05 = hq[1], holo50 = hq[2], holo95 = hq[3],
               holoBelow302 = sum(w[x$holo < 302]) / sum(w))
  }
  rows[[length(rows) + 1L]] <- Row(lab, rep(1, nrow(x)))
  if (grepl("-nocal$", lab)) {
    w <- ifelse(x$holo > 302, exp(-(x$holo - 302) / 94) / 94, 0)
    rows[[length(rows) + 1L]] <- Row(paste0(lab, "-reweighted"), w)
  }
}
res <- do.call(rbind, rows)
options(width = 200)
print(res, digits = 4, row.names = FALSE)
utils::write.csv(res, outFile, row.names = FALSE)
