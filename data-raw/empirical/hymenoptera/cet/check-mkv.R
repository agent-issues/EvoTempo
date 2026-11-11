# Hymenoptera CET pipeline: fixed-tree check of the character likelihood
# (HymMk() in common.R) against RevBayes.
#
#   Rscript data-raw/empirical/hymenoptera/cet/check-mkv.R <rb output> [out.csv]
#
# Run from the package root with EvoTempo installed (or CET_LOAD_ALL=1).
# <rb output> is the stdout of revbayes/check-mkv.Rev (revbayes/run-check-mkv.sh
# writes smoke-results/revbayes-check-mkv-<build>.txt).  On the committed
# start tree with a strict clock (branch length = rate x duration) at rates
# 0.0007, 0.002 and 0.006, for each of the nine (ordering, K) partitions:
# 1. HymMk() per pattern (mkv = "pattern", the default) must equal RevBayes
#    dnPhyloCTMC(coding = "variable"), and mkv = "none" must equal
#    coding = "all".  Tolerance 1e-6 relative per partition (RevBayes values
#    are printed as floor(lnL) and the remainder to 7 significant digits).
# 2. The C++ kernel (mk-columns.cpp) must equal the pure-R reference
#    ColumnLogLikR() on every column, to 1e-9.
# 3. On unordered partitions, the columns without ambiguous cells must equal
#    EvoTempo's MkLikelihood() (C++ and R engines), to 1e-9.
# 4. The correction is inside the powered term (power t gives t x lnL), and
#    the tree term (Holometabola constraint and calibration) is outside it.
# 5. Reports the RevBayes Gamma(4) site-rate model (alpha 1.10) against no
#    site-rate variation: the part of the comparator's character model that
#    CET does not have.
# 6. Cost of one HymMk() evaluation (C++ kernel; mean of CET_MKV_REPS,
#    default 200) and of the R reference.
# Writes [out.csv] (default <cetDir>/smoke-results/mkv-check.csv) and
# <out>-timing.csv, and stops with an error if a check fails.
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/hymenoptera/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
rbFile <- args[[1]]
outFile <- if (length(args) >= 2) args[[2]] else file.path(cetDir, "smoke-results", "mkv-check.csv")
nRep <- as.integer(Sys.getenv("CET_MKV_REPS", "200"))

# ---- RevBayes values ----
rb <- gsub("\"", "", trimws(readLines(rbFile)))
at <- grep("^RESULT ", rb)
stopifnot(length(at) == 3 * 9 * 3)
num <- function(s) as.numeric(trimws(s))
rbv <- do.call(rbind, lapply(at, function(i) {
  f <- as.integer(strsplit(sub("^RESULT ", "", rb[i]), " ")[[1]])
  data.frame(rateIndex = f[1], part = f[2], coding = f[3], lnL = num(rb[i + 1]) + num(rb[i + 2]))
}))
rbRoot <- num(rb[grep("^ROOTAGE", rb) + 1])
rates <- c(0.0007, 0.002, 0.006)

dat <- LoadHymenoptera()
tt <- dat$tree
chars <- dat$characters
parts <- attr(chars, "hym")
cat(sprintf("Start tree: root %.6f Ma (RevBayes reads %.6f)\n", max(tt$age), rbRoot))
# RevBayes partition index -> name (data/partitions.Rev morph_names)
pl <- grep("^PART ", rb, value = TRUE)
rbNames <- setNames(sub("^PART [0-9]+ ", "", pl), as.integer(sub("^PART ([0-9]+) .*$", "\\1", pl)))
partOrder <- vapply(parts, `[[`, "", "label")
stopifnot(setequal(rbNames, partOrder))
rbv$part <- match(rbNames[as.character(rbv$part)], partOrder)
cat("CET partitions:", paste(partOrder, collapse = " "), "\n")
dur <- numeric(tt$nTip + tt$nNode)
child <- which(!is.na(tt$parent))
dur[child] <- tt$age[tt$parent[child]] - tt$age[child]

fail <- character()
rows <- list()
for (ri in seq_along(rates)) {
  len <- rates[ri] * dur
  for (pi in seq_along(parts)) {
    g <- parts[[pi]]
    cpp <- PartitionLogLik(tt, g, len, "cpp")
    llC <- MkColumnLogLik(as.integer(tt$postorder), tt$kids, tt$nTip, tt$root, len, g$V, g$ev,
                          g$tip, g$nCol)
    llR <- ColumnLogLikR(tt, g, len)
    # EvoTempo MkLikelihood() on unambiguous unordered columns
    mkDiff <- NA_real_
    if (!g$ordered) {
      cols <- g$columns
      x <- chars[, cols, drop = FALSE]
      clean <- which(colSums(is.na(x) & nzchar(dat$cells[, cols, drop = FALSE])) == 0)
      if (length(clean)) {
        mine <- sum(llC[clean])
        mkDiff <- max(abs(c(.mkOriginal(tt, x[, clean, drop = FALSE], len, g$K, "cpp") - mine,
                            .mkOriginal(tt, x[, clean, drop = FALSE], len, g$K, "R") - mine)))
      }
    }
    v1 <- rbv$lnL[rbv$rateIndex == ri & rbv$part == pi & rbv$coding == 1]
    v2 <- rbv$lnL[rbv$rateIndex == ri & rbv$part == pi & rbv$coding == 2]
    v3 <- rbv$lnL[rbv$rateIndex == ri & rbv$part == pi & rbv$coding == 3]
    rows[[length(rows) + 1L]] <- data.frame(
      rate = rates[ri], partition = g$label, nChar = g$nData, nMask = g$nMask,
      pattern = cpp[["data"]] + cpp[["correction"]], revbayesVariable = v1,
      diffPattern = cpp[["data"]] + cpp[["correction"]] - v1,
      none = cpp[["data"]], revbayesAll = v2, diffNone = cpp[["data"]] - v2,
      correction = cpp[["correction"]], revbayesGamma = v3, gammaMinusNoGamma = v3 - v1,
      maxDiffCppR = max(abs(llC - llR)), maxDiffMkLikelihood = mkDiff)
  }
}
res <- do.call(rbind, rows)
print(res, digits = 10, row.names = FALSE)
tolRel <- function(a, b) abs(a - b) > 1e-6 * abs(b) + 1e-6
if (any(tolRel(res$pattern, res$revbayesVariable))) {
  fail <- c(fail, "per-pattern lnL differs from RevBayes coding = \"variable\"")
}
if (any(tolRel(res$none, res$revbayesAll))) fail <- c(fail, "uncorrected lnL differs from RevBayes coding = \"all\"")
if (any(res$maxDiffCppR > 1e-9)) fail <- c(fail, "C++ kernel differs from the R reference")
if (any(res$maxDiffMkLikelihood > 1e-9, na.rm = TRUE)) fail <- c(fail, "kernel differs from MkLikelihood()")

tot <- aggregate(cbind(pattern, revbayesVariable, none, revbayesAll, revbayesGamma) ~ rate, res, sum)
tot$diffPattern <- tot$pattern - tot$revbayesVariable
tot$gammaMinusNoGamma <- tot$revbayesGamma - tot$revbayesVariable
print(tot, digits = 10, row.names = FALSE)

# Power and tree term
len <- 0.002 * dur
full <- HymMk(tt, chars, len, power = 1, treeTerm = FALSE)
half <- HymMk(tt, chars, len, power = 0.5, treeTerm = FALSE)
withTree <- HymMk(tt, chars, len, power = 0.5, treeTerm = TRUE)
holo <- CladeAge(tt, .cet$holo)
powerDiff <- half - 0.5 * full
treeDiff <- withTree - half - HoloCalibration(holo)
cat(sprintf("power 0.5: %.3g; tree term: %.3g (Holometabola age %.4f, calibration %.6f)\n",
            powerDiff, treeDiff, holo, HoloCalibration(holo)))
if (abs(powerDiff) > 1e-9) fail <- c(fail, "power does not scale the corrected lnL")
if (abs(treeDiff) > 1e-9) fail <- c(fail, "tree term is not added unpowered")
broken <- tt
broken$tipLabel[match(c("Orthoptera", .cet$holo[1]), broken$tipLabel)] <-
  broken$tipLabel[match(c(.cet$holo[1], "Orthoptera"), broken$tipLabel)]
if (is.finite(TreeTerm(broken))) fail <- c(fail, "constraint violation not detected")

# Cost per evaluation at rate 0.002
Ms <- function(f, n = nRep) {
  f()
  1000 * system.time(for (r in seq_len(n)) f())[["elapsed"]] / n
}
timing <- data.frame(
  what = c("HymMk, C++ kernel, per pattern", "HymMk, C++ kernel, no correction",
           "TreeTerm", "HymMk, R reference"),
  msPerEvaluation = c(Ms(function() HymMk(tt, chars, len)),
                      Ms(function() HymMk(tt, chars, len, mkv = "none")),
                      Ms(function() TreeTerm(tt)),
                      Ms(function() HymMk(tt, chars, len, impl = "R"), 5)))
print(timing, digits = 4, row.names = FALSE)
dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, outFile, row.names = FALSE)
utils::write.csv(tot, sub("\\.csv$", "-totals.csv", outFile), row.names = FALSE)
utils::write.csv(timing, sub("\\.csv$", "-timing.csv", outFile), row.names = FALSE)
if (length(fail)) stop("check-mkv.R FAILED: ", paste(fail, collapse = "; "))
cat("check-mkv.R: all checks passed\n")
