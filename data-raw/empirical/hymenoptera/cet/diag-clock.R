# Hymenoptera CET pipeline: do the characters under a strict clock prefer
# old trees?  (README, "Root ages".)
#
#   Rscript data-raw/empirical/hymenoptera/cet/diag-clock.R <out.csv> <rbmorph root> <label>=<CET chain dir> ...
#
# Draws 20 trees from each source: the comparator's posterior (RevBayes,
# UCLN; pooled post-burn-in samples of the six chains; converted with
# rb-trees.R) and each CET chain (last 75% of its seg*_trees.nwk).  On each
# tree, with branch lengths = rate x duration (strict clock, as CET cell
# (0, 0)), maximises the character log-likelihood over the rate (Gamma(4),
# alpha 1.1, coding "variable"; and without Gamma), and reports the root age,
# the tree length (Myr), the maximum and the maximising rate, and the
# log-likelihood at the rate prior's median, exp(-7.1).
suppressMessages(library(EvoTempo))
source("data-raw/empirical/hymenoptera/cet/common.R")
source(file.path(cetDir, "rb-trees.R"))
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[[1]]
rbRoot <- args[[2]]
spec <- args[-(1:2)]
set.seed(20261117)
tips <- read.csv(file.path(hymDir, "hymenoptera_tips.csv"))
tmp <- tempfile(fileext = ".nwk")
Trees <- list()
post <- ReadRbPosterior(rbRoot)
Trees[["rb-posterior"]] <- lapply(post$newick[sample.int(nrow(post), 20)], function(nw) {
  write.tree(ConvertRbTree(nw, tips), tmp, digits = 10)
  readLines(tmp)
})
for (sp in spec) {
  lab <- sub("=.*$", "", sp); d <- sub("^[^=]*=", "", sp)
  nw <- unlist(lapply(sort(Sys.glob(file.path(d, "seg*_trees.nwk"))), readLines))
  nw <- nw[seq_along(nw) > floor(0.25 * length(nw))]
  Trees[[lab]] <- as.list(nw[sample.int(length(nw), min(20, length(nw)))])
}
Profile <- function(tt, chars, alpha) {
  dur <- numeric(tt$nTip + tt$nNode)
  ch <- which(!is.na(tt$parent))
  dur[ch] <- tt$age[tt$parent[ch]] - tt$age[ch]
  f <- function(lr) MkLikelihood(tt, chars, exp(lr) * dur, alpha = alpha)
  o <- stats::optimize(f, c(log(1e-4), log(2e-2)), maximum = TRUE, tol = 1e-4)
  c(maxLnL = o$objective, rateHat = exp(o$maximum), lnLAtPriorMedian = f(-7.1),
    treeLength = sum(dur))
}
rows <- list()
for (lab in names(Trees)) {
  for (i in seq_along(Trees[[lab]])) {
    writeLines(Trees[[lab]][[i]], tmp)
    dat <- LoadHymenoptera(tmp)
    dat1 <- LoadHymenoptera(tmp, nRateCat = 1L)
    g <- Profile(dat$tree, dat$characters, 1.1)
    n <- Profile(dat1$tree, dat1$characters, NULL)
    rows[[length(rows) + 1L]] <- data.frame(source = lab, tree = i, rootAge = max(dat$tree$age),
      holoAge = CladeAge(dat$tree, .cet$holo), treeLength = g[["treeLength"]],
      maxLnLGamma = g[["maxLnL"]], rateHatGamma = g[["rateHat"]], lnLGammaAtPriorMedian = g[["lnLAtPriorMedian"]],
      maxLnLNoGamma = n[["maxLnL"]], rateHatNoGamma = n[["rateHat"]])
  }
  utils::write.csv(do.call(rbind, rows), outFile, row.names = FALSE)
}
res <- do.call(rbind, rows)
s <- do.call(rbind, lapply(split(res, res$source), function(x) data.frame(
  source = x$source[1], n = nrow(x), rootAge = stats::median(x$rootAge), treeLength = stats::median(x$treeLength),
  maxLnLGamma = stats::median(x$maxLnLGamma), rateHatGamma = stats::median(x$rateHatGamma),
  rateTimesTreeLength = stats::median(x$rateHatGamma * x$treeLength),
  maxLnLNoGamma = stats::median(x$maxLnLNoGamma), rateHatNoGamma = stats::median(x$rateHatNoGamma))))
options(width = 200)
print(s, digits = 6, row.names = FALSE)
utils::write.csv(s, sub("\\.csv$", "-summary.csv", outFile), row.names = FALSE)
