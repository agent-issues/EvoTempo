# Hymenoptera: summary of the RevBayes strict-clock comparator
# (hymenoptera_fbd.Rev, config "morph-strict"; README, "Root ages") against the
# UCLN comparator ("morph") and CET cell (0, 0).
#
#   Rscript data-raw/empirical/hymenoptera/cet/summarise-rbstrict.R <out.csv> <burn-in> \
#     strict=<log> ... ucln=<log> ... [cet=<CET chain dir> ...]
#
# Each chain's first <burn-in> fraction is dropped and chains of one label are
# cut to a common length.  Per label and quantity: mean, 2.5/50/97.5%
# quantiles, and for labels with two or more chains MkPrime's rank-normalised
# split R-hat (.Rhat()) and bulk and tail ESS (.EssBulk(), .EssTail()) on the
# draws x chains matrix.  Quantities: root age, Holometabola age, clock rate,
# the character log-likelihood (RevBayes: Likelihood, the morphology lnL of
# dnPhyloCTMC(coding = "variable") with Gamma(4); CET: charLL, the same
# quantity, check-mkv.R), alpha_morph, turnover, lambda, sampled ancestors.
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[[1]]
burn <- as.numeric(args[[2]])
spec <- args[-(1:2)]
labels <- sub("=.*$", "", spec)
paths <- sub("^[^=]*=", "", spec)
mk <- asNamespace("MkPrime")
Rhat <- get(".Rhat", mk); EssBulk <- get(".EssBulk", mk); EssTail <- get(".EssTail", mk)
Read <- function(p) {
  if (grepl("\\.log$", p)) {
    y <- utils::read.delim(p, check.names = FALSE)
    x <- data.frame(rootAge = y$root_age, holoAge = y$holo_age, clockRate = y$clock_rate,
                    charLnL = y$Likelihood, alpha = y$alpha_morph, turnover = y$turnover,
                    lambda = y$lambda, nSA = y$num_samp_anc)
  } else {
    y <- do.call(rbind, lapply(sort(Sys.glob(file.path(p, "seg*_samples.csv"))), utils::read.csv))
    x <- data.frame(rootAge = y$rootAge, holoAge = y$holoAge, clockRate = y$rate,
                    charLnL = y$charLL, alpha = if (is.null(y$alpha)) NA else y$alpha,
                    turnover = y$mu / y$lambda, lambda = y$lambda, nSA = y$nAncestor)
  }
  x[seq_len(nrow(x)) > floor(burn * nrow(x)), ]
}
rows <- list()
for (lab in unique(labels)) {
  xs <- lapply(paths[labels == lab], Read)
  n <- min(vapply(xs, nrow, 0L))
  xs <- lapply(xs, function(x) utils::tail(x, n))
  for (q in names(xs[[1]])) {
    m <- vapply(xs, function(x) x[[q]], numeric(n))
    if (!is.matrix(m)) m <- matrix(m, ncol = length(xs))
    v <- as.vector(m)
    if (all(is.na(v))) next
    multi <- ncol(m) >= 2 && stats::var(v) > 0
    rows[[length(rows) + 1L]] <- data.frame(
      label = lab, quantity = q, chains = ncol(m), drawsPerChain = n, mean = mean(v),
      q025 = stats::quantile(v, 0.025, names = FALSE), q50 = stats::median(v),
      q975 = stats::quantile(v, 0.975, names = FALSE),
      rhat = if (multi) tryCatch(Rhat(m), error = function(e) NA_real_) else NA_real_,
      essBulk = if (stats::var(v) > 0) tryCatch(EssBulk(m), error = function(e) NA_real_) else NA_real_,
      essTail = if (stats::var(v) > 0) tryCatch(EssTail(m), error = function(e) NA_real_) else NA_real_)
  }
}
res <- do.call(rbind, rows)
options(width = 220)
print(res, digits = 5, row.names = FALSE)
utils::write.csv(res, outFile, row.names = FALSE)
