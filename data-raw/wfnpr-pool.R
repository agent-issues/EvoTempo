# Pools the per-chain results of data-raw/wfnpr-validation.R: every chain of
# every run with the same tree, case and candidate cap (whatever the seed)
# counts once, so that no run is set aside.  For each statistic, z is the mean
# over chains minus the exact value, over the between-chain standard error;
# with n chains it has n - 1 degrees of freedom, and `pMax` is the probability
# that the largest |t| of as many independent statistics exceeds the observed
# one (the statistics are correlated, so it is a guide, not a test).
# Run from the package root: Rscript data-raw/wfnpr-pool.R
dir <- "data-raw/wfnpr-results"
files <- list.files(dir, "-chains\\.csv$", full.names = TRUE)
chains <- do.call(rbind, lapply(files, utils::read.csv))
chains$group <- sub("-[0-9]+$", "", chains$label)
chains <- chains[chains$statistic != "acceptance", ]
pooled <- do.call(rbind, lapply(split(chains, chains[c("group", "move")], drop = TRUE), function(g) {
  n <- length(unique(paste(g$label, g$chain)))
  mean <- tapply(g$value, g$statistic, mean)
  se <- tapply(g$value, g$statistic, stats::sd) / sqrt(n)
  exact <- tapply(g$exact, g$statistic, function(x) x[[1]])
  # Statistics that never vary (clades of probability 0 or 1) carry no test
  z <- ((mean - exact) / se)[se > 0]
  top <- if (length(z) && !anyNA(z)) max(abs(z)) else NA
  data.frame(
    group = g$group[[1]], move = g$move[[1]], chains = n, statistics = length(z),
    maxAbsZ = signif(top, 3), worst = if (is.na(top)) NA else names(z)[[which.max(abs(z))]],
    rmsZ = signif(sqrt(mean(z^2)), 3),
    pMax = signif(1 - (1 - 2 * stats::pt(-top, n - 1))^length(z), 2)
  )
}))
rownames(pooled) <- NULL
utils::write.csv(pooled, file.path(dir, "pooled.csv"), row.names = FALSE)
print(pooled)
