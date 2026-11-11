# Penguin CET pipeline: fixed-tree check of the partitioned Mkv likelihood
# (PartitionedMk() in common.R) against RevBayes (#65).
#
#   Rscript data-raw/empirical/penguins/cet/check-mkv.R [out.csv]
#
# Run from the package root with EvoTempo installed (or CET_LOAD_ALL=1 for
# pkgload::load_all()).  On the round-1 start tree (data/start-tree-round1.nwk,
# the start tree until the round-4 tree replaced it; its branch lengths equal
# those of the revbayes/data/start_tree.nwk of that time bar the origin stem,
# which MkLikelihood() ignores, and the oracle below was computed on it) with
# a strict clock, branch length = rate x duration, at rates 0.005, 0.017 and
# 0.05:
# 1. PartitionedMk() with the per-pattern correction (the default) must equal
#    RevBayes's dnPhyloCTMC(type = "Standard", coding = "variable") lnL summed
#    over the six K-partitions.  Oracle: an independent verifier's RevBayes
#    runs (v1.4.1 and a2a8e2db, identical), per partition to 7 significant
#    digits, on revbayes/data/start_tree.nwk (issue #65).  Tolerance 1e-3.
# 2. The earlier all-tips-observed variant (.cet$mkv = "allTips") must still
#    give the values that the verifier reproduced for it.
# 3. ConstantPatternLogLik() must equal MkLikelihood() on each pattern's
#    constant column (unobserved tips unknown), in both engines, to 1e-10.
# 4. The correction is inside the powered term: power t gives t x lnL.
# 5. Cost of one PartitionedMk() evaluation per variant and engine
#    (milliseconds, mean of CET_MKV_REPS evaluations, default 200).
# Writes [out.csv] (default <cetDir>/smoke-results/mkv-check.csv) and stops
# with an error if a check fails.
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/penguins/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
outFile <- if (length(args)) args[[1]] else file.path(cetDir, "smoke-results", "mkv-check.csv")
nRep <- as.integer(Sys.getenv("CET_MKV_REPS", "200"))

# The oracle values are for this tree, not for the current start tree
dat <- LoadPenguins(file.path(cetDir, "data", "start-tree-round1.nwk"))
tt <- dat$tree
chars <- dat$characters
pat <- attr(chars, "mkv")
dur <- numeric(tt$nTip + tt$nNode)
child <- which(!is.na(tt$parent))
dur[child] <- tt$age[tt$parent[child]] - tt$age[child]
cat(sprintf("%d characters, %d missing-data patterns within K-partitions\n",
            ncol(chars), length(pat$k)))

oracle <- data.frame(
  rate = c(0.005, 0.017, 0.05),
  revbayesVariable = c(-1778.51784, -1656.30055, -1820.97102),
  allTipsBefore = c(-1861.343283, -1693.677249, -1827.535176))
engines <- c("cpp", "R")
fail <- character()
rows <- list()
for (i in seq_len(nrow(oracle))) {
  len <- oracle$rate[i] * dur
  for (eng in engines) {
    v <- c(pattern = PartitionedMk(tt, chars, len, engine = eng, power = 1, mkv = "pattern"),
           allTips = PartitionedMk(tt, chars, len, engine = eng, power = 1, mkv = "allTips"),
           none = PartitionedMk(tt, chars, len, engine = eng, power = 1, mkv = "none"))
    half <- PartitionedMk(tt, chars, len, engine = eng, power = 0.5, mkv = "pattern")
    lc <- ConstantPatternLogLik(tt, pat, len)
    ref <- vapply(seq_along(pat$k), function(j) {
      col <- matrix(ifelse(pat$observed[, j], 0L, NA_integer_), ncol = 1)
      .mkOriginal(tt, col, len, pat$k[j], eng)
    }, 0)
    rows[[length(rows) + 1L]] <- data.frame(
      rate = oracle$rate[i], engine = eng, pattern = v[["pattern"]],
      revbayesVariable = oracle$revbayesVariable[i],
      diffPattern = v[["pattern"]] - oracle$revbayesVariable[i],
      allTips = v[["allTips"]], allTipsBefore = oracle$allTipsBefore[i],
      diffAllTips = v[["allTips"]] - oracle$allTipsBefore[i], none = v[["none"]],
      correctionPattern = v[["pattern"]] - v[["none"]],
      correctionAllTips = v[["allTips"]] - v[["none"]],
      maxDiffConstColumn = max(abs(lc - ref)), powerHalfDiff = half - 0.5 * v[["pattern"]])
  }
}
res <- do.call(rbind, rows)
print(res, digits = 10, row.names = FALSE)
if (any(abs(res$diffPattern) > 1e-3)) fail <- c(fail, "per-pattern lnL differs from RevBayes coding = \"variable\"")
if (any(abs(res$diffAllTips) > 1e-3)) fail <- c(fail, "all-tips variant differs from its earlier values")
if (any(res$maxDiffConstColumn > 1e-10)) fail <- c(fail, "ConstantPatternLogLik() differs from MkLikelihood()")
if (any(abs(res$powerHalfDiff) > 1e-9)) fail <- c(fail, "power does not scale the corrected lnL")
engDiff <- max(abs(res$pattern[res$engine == "cpp"] - res$pattern[res$engine == "R"]))
cat(sprintf("max |cpp - R| (per-pattern lnL) = %.3g\n", engDiff))
if (engDiff > 1e-8 * max(abs(res$pattern))) fail <- c(fail, "engines disagree")

# Cost per evaluation, at the rate of the start values (0.017)
len <- 0.017 * dur
Ms <- function(f) {
  f()
  1000 * system.time(for (r in seq_len(nRep)) f())[["elapsed"]] / nRep
}
timing <- do.call(rbind, lapply(engines, function(eng) {
  data.frame(engine = eng, variant = c("none", "allTips", "pattern"),
             msPerEvaluation = c(
               Ms(function() PartitionedMk(tt, chars, len, engine = eng, power = 1, mkv = "none")),
               Ms(function() PartitionedMk(tt, chars, len, engine = eng, power = 1, mkv = "allTips")),
               Ms(function() PartitionedMk(tt, chars, len, engine = eng, power = 1, mkv = "pattern"))))
}))
timing$ratioToAllTips <- timing$msPerEvaluation /
  rep(timing$msPerEvaluation[timing$variant == "allTips"], each = 3)
print(timing, digits = 4, row.names = FALSE)
dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, outFile, row.names = FALSE)
utils::write.csv(timing, sub("\\.csv$", "-timing.csv", outFile), row.names = FALSE)
if (length(fail)) stop("check-mkv.R FAILED: ", paste(fail, collapse = "; "))
cat("check-mkv.R: all checks passed\n")
