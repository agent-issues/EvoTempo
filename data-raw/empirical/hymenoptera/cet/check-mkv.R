# Hymenoptera CET pipeline: fixed-tree check of the character likelihood
# (the package's MkCharacters() / MkLikelihood(), #97, #98, as common.R uses
# them) against RevBayes.
#
#   Rscript data-raw/empirical/hymenoptera/cet/check-mkv.R <rb output> [tree] [out.csv]
#
# Run from the package root with EvoTempo installed (or CET_LOAD_ALL=1).
# <rb output> is the stdout of revbayes/check-mkv.Rev on [tree] (default the
# reference start tree data/start-tree.nwk; revbayes/run-check-mkv.sh writes
# smoke-results/revbayes-check-mkv-fd807d10-<tag>.txt).  With a strict clock
# (branch length = rate x duration) at rates 0.0007, 0.002 and 0.006, for each
# of the nine (ordering, K) partitions, both engines (C++ and R):
# 1. Gamma off (nRateCat = 1): coding = "variable" must equal RevBayes
#    dnPhyloCTMC(coding = "variable") and coding = "all" must equal
#    coding = "all".
# 2. Gamma on (nRateCat = 4): coding = "variable" at alpha 1.10 and 0.35, and
#    coding = "all" at alpha 1.10, must equal RevBayes with
#    siteRates = fnDiscretizeGamma(alpha, alpha, 4).
#    Tolerance for 1 and 2: 1e-6 relative per partition (RevBayes values are
#    printed as floor(lnL) and the remainder to 14 significant digits; its
#    Gamma quantiles are approximate (AS 91), which moves lnL by up to about
#    1e-6 at alpha 1.1: README of data-raw/gamma-acrv-results).  The RevBayes
#    rates are reported against GammaRates().
# 3. The C++ and R engines agree to 1e-9 (relative) per partition.
# 4. The nine partitions sum to the likelihood of the whole matrix as
#    LoadHymenoptera() builds it (one MkCharacters() object), to 1e-8.
# 5. The power multiplies the corrected lnL (Gamma mixture included), and the
#    tree term (Holometabola constraint and calibration) is added outside it;
#    a tree that breaks the constraint gives -Inf.
# 6. Characters whose rows are not in the tree's tip order are matched by row
#    name or refused, never scored against the wrong tips.
# 7. Cost of one evaluation of the whole matrix (C++; mean of CET_MKV_REPS,
#    default 200) with and without Gamma, and of the tree term.
# Writes [out.csv] (default <cetDir>/smoke-results/mkv-check-<tree>.csv),
# <out>-totals.csv and <out>-timing.csv, and stops with an error if a check
# fails.
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/hymenoptera/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
rbFile <- args[[1]]
treeFile <- if (length(args) >= 2) args[[2]] else referenceTree
treeTag <- sub("\\.nwk$", "", basename(treeFile))
outFile <- if (length(args) >= 3) args[[3]] else {
  file.path(cetDir, "smoke-results", sprintf("mkv-check-%s.csv", treeTag))
}
nRep <- as.integer(Sys.getenv("CET_MKV_REPS", "200"))

# ---- RevBayes values ----
rb <- gsub("\"", "", trimws(readLines(rbFile)))
num <- function(s) as.numeric(trimws(s))
at <- grep("^RESULT ", rb)
stopifnot(length(at) == 3 * 9 * 5)
rbv <- do.call(rbind, lapply(at, function(i) {
  f <- as.integer(strsplit(sub("^RESULT ", "", rb[i]), " ")[[1]])
  data.frame(rateIndex = f[1], part = f[2], coding = f[3], lnL = num(rb[i + 1]) + num(rb[i + 2]))
}))
atR <- grep("^RATES ", rb)
stopifnot(length(atR) == 8)
rbRates <- do.call(rbind, lapply(atR, function(i) {
  f <- as.integer(strsplit(sub("^RATES ", "", rb[i]), " ")[[1]])
  data.frame(alphaIndex = f[1], category = f[2], rate = num(rb[i + 1]) + num(rb[i + 2]))
}))
rbRoot <- num(rb[grep("^ROOTAGE", rb) + 1])
rates <- c(0.0007, 0.002, 0.006)
alphas <- c(1.10, 0.35)
# RevBayes codings (check-mkv.Rev): CET coding, number of categories, alpha
codings <- data.frame(cd = 1:5, coding = c("variable", "all", "variable", "variable", "all"),
                      nRateCat = c(1L, 1L, 4L, 4L, 4L), alpha = c(NA, NA, 1.10, 0.35, 1.10))

dat <- LoadHymenoptera(treeFile)
tt <- dat$tree
cat(sprintf("Tree %s: root %.6f Ma (RevBayes reads %.6f)\n", treeFile, max(tt$age), rbRoot))
fail <- character()
if (abs(max(tt$age) - rbRoot) > 1e-5) fail <- c(fail, "RevBayes read another tree")
gr <- do.call(rbind, lapply(seq_along(alphas), function(a) {
  x <- rbRates[rbRates$alphaIndex == a, ]
  data.frame(alpha = alphas[a], category = x$category, revbayes = x$rate,
             GammaRates = GammaRates(alphas[a], 4)[x$category])
}))
gr$difference <- gr$GammaRates - gr$revbayes
cat("Discretised Gamma rates, RevBayes (AS 91 quantiles) and GammaRates():\n")
print(gr, digits = 15, row.names = FALSE)

# RevBayes partition index -> name (data/partitions.Rev morph_names)
pl <- grep("^PART ", rb, value = TRUE)
rbNames <- setNames(sub("^PART [0-9]+ ", "", pl), as.integer(sub("^PART ([0-9]+) .*$", "\\1", pl)))
parts <- dat$partitions
labels <- unique(parts$partition)
stopifnot(setequal(rbNames, labels))
dur <- numeric(tt$nTip + tt$nNode)
child <- which(!is.na(tt$parent))
dur[child] <- tt$age[tt$parent[child]] - tt$age[child]
PartChars <- function(lab, coding, nRateCat) {
  cols <- which(parts$partition == lab)
  MkCharacters(dat$cells[, cols, drop = FALSE], nState = parts$nState[cols],
               ordered = parts$ordered[cols], coding = coding, nRateCat = nRateCat)
}
rows <- list()
for (ri in seq_along(rates)) {
  len <- rates[ri] * dur
  for (lab in labels) {
    pi <- as.integer(names(rbNames)[rbNames == lab])
    for (k in seq_len(nrow(codings))) {
      cd <- codings[k, ]
      x <- PartChars(lab, cd$coding, cd$nRateCat)
      a <- if (cd$nRateCat > 1L) cd$alpha else NULL
      cpp <- MkLikelihood(tt, x, len, engine = "cpp", alpha = a)
      r <- MkLikelihood(tt, x, len, engine = "R", alpha = a)
      v <- rbv$lnL[rbv$rateIndex == ri & rbv$part == pi & rbv$coding == cd$cd]
      rows[[length(rows) + 1L]] <- data.frame(
        rate = rates[ri], partition = lab, nChar = x$nChar, coding = cd$coding,
        nRateCat = cd$nRateCat, alpha = if (is.null(a)) NA else a, cet = cpp, revbayes = v,
        diff = cpp - v, relDiff = (cpp - v) / abs(v), diffR = r - v,
        relCppR = abs(cpp - r) / abs(r))
    }
  }
}
res <- do.call(rbind, rows)
print(res, digits = 10, row.names = FALSE)
Bad <- function(a, b) abs(a - b) > 1e-6 * abs(b) + 1e-6
if (any(Bad(res$cet[res$nRateCat == 1], res$revbayes[res$nRateCat == 1]))) {
  fail <- c(fail, "Gamma off: lnL differs from RevBayes")
}
if (any(Bad(res$cet[res$nRateCat > 1], res$revbayes[res$nRateCat > 1]))) {
  fail <- c(fail, "Gamma on: lnL differs from RevBayes")
}
if (any(Bad(res$diffR + res$revbayes, res$revbayes))) fail <- c(fail, "R engine differs from RevBayes")
if (any(res$relCppR > 1e-9)) fail <- c(fail, "C++ and R engines differ")
tot <- aggregate(cbind(cet, revbayes) ~ rate + coding + nRateCat + alpha, transform(res, alpha = ifelse(is.na(alpha), 0, alpha)), sum)
tot$alpha[tot$alpha == 0] <- NA
tot$diff <- tot$cet - tot$revbayes
tot <- tot[order(tot$rate, tot$coding, tot$nRateCat, -tot$alpha), ]
print(tot, digits = 12, row.names = FALSE)
cat(sprintf("Largest |CET - RevBayes| per partition: Gamma off %.3g, Gamma on %.3g (relative %.3g); C++ vs R %.3g (relative)\n",
            max(abs(res$diff[res$nRateCat == 1])), max(abs(res$diff[res$nRateCat > 1])),
            max(abs(res$relDiff[res$nRateCat > 1])), max(res$relCppR)))

# Partitions sum to the whole matrix
len <- 0.002 * dur
sumParts <- sum(vapply(labels, function(lab) MkLikelihood(tt, PartChars(lab, "variable", 4L), len,
                                                           alpha = 1.1), 0))
whole <- MkLikelihood(tt, dat$characters, len, alpha = 1.1)
cat(sprintf("Whole matrix (LoadHymenoptera(), alpha 1.1, rate 0.002): %.8f; sum of partitions %.8f\n",
            whole, sumParts))
if (abs(whole - sumParts) > 1e-8) fail <- c(fail, "partitions do not sum to the whole matrix")

# Power and tree term
full <- CetMk(tt, dat$characters, len, alpha = 1.1, power = 1, treeTerm = FALSE)
half <- CetMk(tt, dat$characters, len, alpha = 1.1, power = 0.5, treeTerm = FALSE)
withTree <- CetMk(tt, dat$characters, len, alpha = 1.1, power = 0.5, treeTerm = TRUE)
holo <- CladeAge(tt, .cet$holo)
powerDiff <- half - 0.5 * full
treeDiff <- withTree - half - HoloCalibration(holo)
cat(sprintf("power 0.5: %.3g; tree term: %.3g (Holometabola age %.4f, calibration %.6f)\n",
            powerDiff, treeDiff, holo, HoloCalibration(holo)))
if (abs(full - whole) > 1e-9) fail <- c(fail, "CetMk() at power 1 differs from MkLikelihood()")
if (abs(powerDiff) > 1e-9) fail <- c(fail, "power does not scale the corrected lnL")
if (abs(treeDiff) > 1e-9) fail <- c(fail, "tree term is not added unpowered")
broken <- tt
broken$tipLabel[match(c("Orthoptera", .cet$holo[1]), broken$tipLabel)] <-
  broken$tipLabel[match(c(.cet$holo[1], "Orthoptera"), broken$tipLabel)]
if (is.finite(TreeTerm(broken))) fail <- c(fail, "constraint violation not detected")

# Rows out of tip order: matched by name, or refused
perm <- c(2L, 1L, seq.int(3L, nrow(dat$cells)))
shuffled <- CetCharacters(dat$cells[perm, ], parts)
shuf <- tryCatch(MkLikelihood(tt, shuffled, len, alpha = 1.1), error = function(e) conditionMessage(e))
cat("Rows out of tip order:", if (is.numeric(shuf)) sprintf("lnL %.8f (in order %.8f)", shuf, whole) else
  paste("refused:", shuf), "\n")
if (is.numeric(shuf) && abs(shuf - whole) > 1e-8) fail <- c(fail, "rows out of tip order scored against the wrong tips")

# Cost per evaluation at rate 0.002
Ms <- function(f, n = nRep) {
  f()
  1000 * system.time(for (r in seq_len(n)) f())[["elapsed"]] / n
}
noGamma <- LoadHymenoptera(treeFile, nRateCat = 1L)$characters
allCoding <- LoadHymenoptera(treeFile, mkv = "none")$characters
timing <- data.frame(
  what = c("whole matrix, C++, variable, no Gamma", "whole matrix, C++, variable, Gamma(4)",
           "whole matrix, C++, all, Gamma(4)", "TreeTerm", "whole matrix, R, variable, Gamma(4)"),
  msPerEvaluation = c(Ms(function() MkLikelihood(tt, noGamma, len)),
                      Ms(function() MkLikelihood(tt, dat$characters, len, alpha = 1.1)),
                      Ms(function() MkLikelihood(tt, allCoding, len, alpha = 1.1)),
                      Ms(function() TreeTerm(tt)),
                      Ms(function() MkLikelihood(tt, dat$characters, len, engine = "R", alpha = 1.1), 5)))
print(timing, digits = 4, row.names = FALSE)
dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, outFile, row.names = FALSE)
utils::write.csv(tot, sub("\\.csv$", "-totals.csv", outFile), row.names = FALSE)
utils::write.csv(gr, sub("\\.csv$", "-gamma-rates.csv", outFile), row.names = FALSE)
utils::write.csv(timing, sub("\\.csv$", "-timing.csv", outFile), row.names = FALSE)
if (length(fail)) stop("check-mkv.R FAILED: ", paste(fail, collapse = "; "))
cat("check-mkv.R: all checks passed\n")
