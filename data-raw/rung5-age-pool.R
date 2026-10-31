# Rung 5 with a free start age (#18): pooled analysis of independent
# replicates of data-raw/rung5-age-validation.R.  All chains of all replicates
# of a case target the same distribution, and all simulated trees of a case
# come from the same simulator, so they can be pooled: z is the difference of
# the pooled chain-mean (or chain-CDF) estimate and the pooled simulation
# estimate over its standard error, the spread of the independent chains
# together with the simulation's sampling error.  Per-replicate z of the means
# are summarised too; with 16 chains each they follow roughly a t15 law
# (sd 1.07).
# Usage: Rscript data-raw/rung5-age-pool.R <dir with rung5-age-*.rds> [outDir]
#   Replicates whose file names match EVOTEMPO_POOL_EXCLUDE (a regular
#   expression) are left out.
args <- commandArgs(trailingOnly = TRUE)
inDir <- args[[1]]
outDir <- if (length(args) > 1) args[[2]] else "data-raw/rung5-age-results"
files <- list.files(inDir, "^rung5-age-.*\\.rds$", full.names = TRUE)
exclude <- Sys.getenv("EVOTEMPO_POOL_EXCLUDE")
if (nzchar(exclude)) files <- files[!grepl(exclude, basename(files))]
runs <- lapply(files, function(f) {
  x <- readRDS(f)
  list(case = paste0(x$case, if (x$fixed) "-fixed" else "-free"), sims = x$sims,
       chains = lapply(x$chains, `[[`, "samples"), file = basename(f))
})
quantities <- c("lambda", "mu", "psi", "theta", "s", "treeLength", "meanNodeAge", "rootAge",
                "originAge", "oldestFossil", "youngestFossil", "clock", "inverseClock",
                "nAncestor", "cherries")
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
Finite <- function(v) v[is.finite(v)]
pooled <- list()
perRun <- list()
for (cs in unique(vapply(runs, `[[`, "", "case"))) {
  R <- Filter(function(x) x$case == cs, runs)
  free <- grepl("-free$", cs)
  for (q in intersect(quantities, names(R[[1]]$sims))) {
    if (!free && q %in% c("lambda", "mu", "psi", "theta", "s")) next
    s <- Finite(unlist(lapply(R, function(x) x$sims[[q]])))
    if (stats::var(s) == 0) next
    per <- unlist(lapply(R, function(x) lapply(x$chains, function(ch) Finite(ch[[q]]))),
                  recursive = FALSE)
    Row <- function(stat, simEst, chainEst, simVar) {
      data.frame(case = cs, quantity = q, stat = stat, nRep = length(R),
                 nChain = length(chainEst), nSim = length(s), sim = simEst,
                 mcmc = mean(chainEst),
                 z = (mean(chainEst) - simEst) /
                   sqrt(stats::var(chainEst) / length(chainEst) + simVar / length(s)))
    }
    cm <- vapply(per, mean, 0)
    pooled[[length(pooled) + 1L]] <- Row("mean", mean(s), cm, stats::var(s))
    if (q %in% c("nAncestor", "cherries")) {
      for (v in sort(unique(s))) {
        ps <- mean(s == v)
        pooled[[length(pooled) + 1L]] <- Row(paste0("P(=", v, ")"), ps,
                                             vapply(per, function(x) mean(x == v), 0),
                                             ps * (1 - ps))
      }
    } else {
      for (pr in probs) {
        qs <- stats::quantile(s, pr, names = FALSE)
        pooled[[length(pooled) + 1L]] <- Row(sprintf("q%02d", round(100 * pr)), pr,
                                             vapply(per, function(x) mean(x <= qs), 0),
                                             pr * (1 - pr))
      }
    }
    for (x in R) {
      ss <- Finite(x$sims[[q]])
      m <- vapply(x$chains, function(ch) mean(Finite(ch[[q]])), 0)
      perRun[[length(perRun) + 1L]] <- data.frame(
        case = cs, run = sub("\\.rds$", "", x$file), quantity = q, sim = mean(ss),
        mcmc = mean(m),
        z = (mean(m) - mean(ss)) / sqrt(stats::var(m) / length(m) + stats::var(ss) / length(ss)))
    }
  }
}
pooled <- do.call(rbind, pooled)
perRun <- do.call(rbind, perRun)
pooled$holm <- stats::ave(2 * stats::pnorm(-abs(pooled$z)), pooled$case,
                          FUN = function(p) stats::p.adjust(p, "holm"))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
write.csv(pooled, file.path(outDir, "rung5-age-pooled.csv"), row.names = FALSE)
write.csv(perRun, file.path(outDir, "rung5-age-pooled-per-run.csv"), row.names = FALSE)
for (cs in unique(pooled$case)) {
  p <- pooled[pooled$case == cs, ]
  r <- perRun[perRun$case == cs, ]
  cat(sprintf("%s: %d replicates, %d tests, max |z| %.2f (%s %s), Holm min p %.3f\n", cs,
              p$nRep[[1]], nrow(p), max(abs(p$z)), p$quantity[which.max(abs(p$z))],
              p$stat[which.max(abs(p$z))], min(p$holm)))
  print(aggregate(z ~ quantity, r, function(z) c(mean = mean(z), sd = stats::sd(z))))
}
