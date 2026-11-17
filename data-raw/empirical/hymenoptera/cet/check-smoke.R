# Hymenoptera CET pipeline: checks on a smoke run (smoke.sh), as the penguin
# pipeline's check-smoke.R, without its fixed segment length.
#
#   [CET_MKV=<variant>] Rscript data-raw/empirical/hymenoptera/cet/check-smoke.R <W> [smoke dir]
#
# Reports seconds per iteration per cell and segment, acceptance rates
# (including the FBD block moves and the alpha move), each chain's start tree
# (progress.txt) and the step sizes its checkpoint holds, and checks that the
# four cells started from four different start trees, that every sample has
# an alpha, and that the resumed chain (two segments, <smoke dir>) equals the
# uninterrupted one (one segment, <smoke dir>-identity) exactly: samples,
# trees and states.  Writes <smoke dir>/timing.csv and acceptance.csv; exits
# with an error if a check fails.
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args)) args[[1]] else "/nobackup/pjjg18/evotempo-hym-cet"
smoke <- if (length(args) >= 2) args[[2]] else file.path(W, "smoke")
mkv <- Sys.getenv("CET_MKV", "pattern")
sfx <- if (mkv == "pattern") "" else paste0("_mkv", mkv)
cells <- c("b0g0", "b0g1", "b1g0", "b1g1")
labels <- paste0(cells, sfx)
Samples <- function(dir) {
  do.call(rbind, lapply(sort(Sys.glob(file.path(dir, "seg*_samples.csv"))), utils::read.csv))
}
timing <- do.call(rbind, lapply(seq_along(labels), function(i) {
  dir <- file.path(smoke, paste0(labels[i], "_run1"))
  acc <- utils::read.csv(file.path(dir, "acceptance.csv"))
  seg <- unique(acc[, c("segment", "secs", "burnin")])
  smp <- Samples(dir)
  it <- readRDS(file.path(dir, "checkpoint.rds"))$seg + seg$burnin
  data.frame(cell = cells[i], segment = seg$segment, secs = seg$secs, iterations = it,
             secsPerIteration = seg$secs / it,
             finalLogPost = tapply(smp$logPost, smp$segment, utils::tail, 1)[as.character(seg$segment)],
             finalCharLL = tapply(smp$charLL, smp$segment, utils::tail, 1)[as.character(seg$segment)],
             finalRootAge = tapply(smp$rootAge, smp$segment, utils::tail, 1)[as.character(seg$segment)],
             finalHoloAge = tapply(smp$holoAge, smp$segment, utils::tail, 1)[as.character(seg$segment)],
             finalLambda = tapply(smp$lambda, smp$segment, utils::tail, 1)[as.character(seg$segment)],
             finalAlpha = tapply(smp$alpha, smp$segment, utils::tail, 1)[as.character(seg$segment)])
}))
print(timing, row.names = FALSE)
utils::write.csv(timing, file.path(smoke, "timing.csv"), row.names = FALSE)
acc <- do.call(rbind, lapply(seq_along(labels), function(i) {
  a <- utils::read.csv(file.path(smoke, paste0(labels[i], "_run1"), "acceptance.csv"))
  cbind(cell = cells[i], a)
}))
last <- acc[acc$segment == max(acc$segment), c("cell", "move", "rate")]
wide <- stats::reshape(last, idvar = "move", timevar = "cell", direction = "wide")
print(wide, row.names = FALSE, digits = 3)
utils::write.csv(acc, file.path(smoke, "acceptance.csv"), row.names = FALSE)
b11 <- paste0("b1g1", sfx, "_run1")
a <- utils::read.csv(file.path(smoke, b11, "seg0002_samples.csv"))
b <- Samples(file.path(paste0(smoke, "-identity"), b11))
b <- b[b$iteration %in% a$iteration, ]
cols <- setdiff(intersect(names(a), names(b)), c("iteration", "segment"))
diff <- max(abs(as.matrix(a[, cols]) - as.matrix(b[, cols])), na.rm = TRUE)
ta <- readLines(file.path(smoke, b11, "seg0002_trees.nwk"))
tb <- utils::tail(readLines(Sys.glob(file.path(paste0(smoke, "-identity"), b11, "seg*_trees.nwk"))), length(ta))
sa <- readRDS(file.path(smoke, b11, "seg0002_states.rds"))
sb <- utils::tail(unlist(lapply(sort(Sys.glob(file.path(paste0(smoke, "-identity"), b11, "seg*_states.rds"))),
                                readRDS), recursive = FALSE), length(sa))
cat(sprintf("Resume check: %d samples x %d columns, max |difference| = %g; trees identical: %s; states identical: %s\n",
            nrow(a), length(cols), diff, identical(ta, tb), identical(sa, sb)))
fail <- character()
if (!isTRUE(diff == 0) || !identical(ta, tb) || !identical(sa, sb)) fail <- c(fail, "resumed chain differs")
# Start trees and step sizes
starts <- do.call(rbind, lapply(seq_along(labels), function(i) {
  dir <- file.path(smoke, paste0(labels[i], "_run1"))
  ckp <- readRDS(file.path(dir, "checkpoint.rds"))
  cap <- ckp$captured
  data.frame(cell = cells[i], startTree = basename(ckp$startTree), md5 = substr(ckp$startTreeMd5, 1, 8),
             nRateCat = ckp$nRateCat, alphaStep = cap$alphaStep,
             step = paste(names(cap$step), signif(cap$step, 3), sep = "=", collapse = " "),
             blockStep = paste(signif(cap$blockStep, 3), collapse = " "),
             alphaInSamples = all(is.finite(Samples(dir)$alpha)))
}))
print(starts, row.names = FALSE)
if (anyDuplicated(starts$md5)) fail <- c(fail, "cells share a start tree")
if (!all(starts$alphaInSamples)) fail <- c(fail, "alpha missing from the samples")
if (!any(acc$move == "alpha")) fail <- c(fail, "no alpha move in acceptance.csv")
if (length(fail)) stop("check-smoke.R FAILED: ", paste(fail, collapse = "; "))
cat("check-smoke.R: all checks passed\n")
