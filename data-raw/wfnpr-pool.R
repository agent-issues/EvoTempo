# Pools the per-chain results of data-raw/wfnpr-validation.R: every chain of
# every run with the same tree, case and candidate cap (whatever the seed)
# counts once, so that no run is set aside.  For each statistic, and with the
# between-chain standard errors:
#   pooled.csv        z of each move's mean against the exact value (for runs
#                     made before the gamma = 1 oracle, taken from
#                     <tree>-<case>-exact.csv);
#   pooled-fnpr.csv   z of the `wfnpr` mean against the `fnpr` mean.
# `pMax` is the probability that the largest |t| of as many independent
# statistics exceeds the observed one (n - 1 degrees of freedom; for the
# comparison of moves, a normal reference).  The statistics are correlated, so
# it is a guide, not a test.  Statistics that never vary (clades of
# probability 0 or 1) carry no test and are left out.
# Run from the package root: Rscript data-raw/wfnpr-pool.R
dir <- "data-raw/wfnpr-results"
files <- list.files(dir, "-chains\\.csv$", full.names = TRUE)
chains <- do.call(rbind, lapply(files, utils::read.csv))
chains$group <- sub("-[0-9]+$", "", chains$label)
chains <- chains[chains$statistic != "acceptance", ]
chains$id <- paste(chains$label, chains$chain)
for (case in unique(sub("-max.*", "", chains$group[is.na(chains$exact)]))) {
  known <- utils::read.csv(file.path(dir, paste0(case, "-exact.csv")))
  fill <- is.na(chains$exact) & startsWith(chains$group, paste0(case, "-max"))
  chains$exact[fill] <- known$exact[match(chains$statistic[fill], known$statistic)]
}

# Mean and standard error over chains of each statistic that every run of the
# group reports (with gamma = 1 the topologies listed depend on the run)
Summarize <- function(g) {
  n <- length(unique(g$id))
  count <- tapply(g$value, g$statistic, length)
  g <- g[g$statistic %in% names(count)[count == n], ]
  list(n = n, mean = tapply(g$value, g$statistic, mean),
       se = tapply(g$value, g$statistic, stats::sd) / sqrt(n),
       exact = tapply(g$exact, g$statistic, function(x) x[[1]]))
}
Row <- function(group, z, p, ...) {
  top <- if (length(z) && !anyNA(z)) max(abs(z)) else NA
  data.frame(
    group = group, ..., statistics = length(z), maxAbsZ = signif(top, 3),
    worst = if (is.na(top)) NA else names(z)[[which.max(abs(z))]],
    rmsZ = signif(sqrt(mean(z^2)), 3),
    pMax = if (is.na(top)) NA else signif(1 - (1 - p(top))^length(z), 2)
  )
}
byGroup <- split(chains, chains$group)
exact <- do.call(rbind, lapply(names(byGroup), function(group) {
  do.call(rbind, lapply(c("wfnpr", "fnpr"), function(move) {
    s <- Summarize(byGroup[[group]][byGroup[[group]]$move == move, ])
    z <- ((s$mean - s$exact) / s$se)[s$se > 0]
    z <- z[!is.na(z)]  # topologies that an early gamma = 1 run listed by frequency
    Row(group, z, function(t) 2 * stats::pt(-t, s$n - 1), move = move, chains = s$n)
  }))
}))
moves <- do.call(rbind, lapply(names(byGroup), function(group) {
  s <- lapply(c(wfnpr = "wfnpr", fnpr = "fnpr"), function(move) {
    Summarize(byGroup[[group]][byGroup[[group]]$move == move, ])
  })
  both <- intersect(names(s$wfnpr$mean), names(s$fnpr$mean))
  se <- sqrt(s$wfnpr$se[both]^2 + s$fnpr$se[both]^2)
  z <- ((s$wfnpr$mean[both] - s$fnpr$mean[both]) / se)[se > 0]
  Row(group, z, function(t) 2 * stats::pnorm(-t), chainsPerMove = s$wfnpr$n)
}))
utils::write.csv(exact[!is.na(exact$maxAbsZ), ], file.path(dir, "pooled.csv"), row.names = FALSE)
utils::write.csv(moves, file.path(dir, "pooled-fnpr.csv"), row.names = FALSE)
print(exact[!is.na(exact$maxAbsZ), ], row.names = FALSE)
print(moves, row.names = FALSE)
