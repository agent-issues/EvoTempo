# Penguin CET pipeline: checks of the segment cap (CET_MAX_SEG,
# CET_MAX_SEG_RUNG; common.R, SegCap() and SegFiles(); #107 review, R107-9) on
# synthetic chain directories.
#
#   Rscript data-raw/empirical/penguins/cet/check-segcap.R
#
# Stops at the first failed check.
source("data-raw/empirical/penguins/cet/common.R")
root <- tempfile("segcap")
MakeChain <- function(label, nSeg, complete = TRUE) {
  d <- file.path(root, label)
  dir.create(d, recursive = TRUE)
  for (k in seq_len(nSeg)) {
    utils::write.csv(data.frame(logPost = rnorm(100), charLL = rnorm(100), iteration = (k - 1) * 100 + 1:100),
                     file.path(d, sprintf("seg%04d_samples.csv", k)), row.names = FALSE)
    saveRDS(as.list(seq_len(100)), file.path(d, sprintf("seg%04d_states.rds", k)))
  }
  writeLines(c(sprintf("x segment %d/%d done", nSeg, nSeg), if (complete) "x complete" else "x resuming"),
             file.path(d, "progress.txt"))
  d
}
prod <- MakeChain("b0g0_run1", 46, complete = FALSE)
rung <- MakeChain("b0g0.25_run1", 23, complete = FALSE)
ladder <- MakeChain("b0g0_t0.5_run1", 12, complete = TRUE)
Check <- function(x, msg) {
  if (!isTRUE(x)) stop("FAILED: ", msg)
  cat("ok  ", msg, "\n")
}
Sys.unsetenv(c("CET_MAX_SEG", "CET_MAX_SEG_RUNG"))
Check(length(SegFiles(prod)) == 46, "no cap: every segment")
Check(!ChainComplete(prod) && ChainComplete(ladder), "no cap: completeness from progress.txt")
Sys.setenv(CET_MAX_SEG = "40", CET_MAX_SEG_RUNG = "20")
Check(length(SegFiles(prod)) == 40 && length(SegFiles(prod, "states.rds")) == 40, "production cell capped at CET_MAX_SEG")
Check(length(SegFiles(rung)) == 20 && length(SegFiles(ladder)) == 12, "rung and power chains capped at CET_MAX_SEG_RUNG")
Check(identical(basename(SegFiles(rung)), sprintf("seg%04d_samples.csv", 1:20)), "the cap keeps the first segments")
Check(ChainComplete(prod) && ChainComplete(rung) && !ChainComplete(ladder),
      "capped: complete when the capped set is full (12 < 20 is not)")
Check(SegmentFingerprint(prod)$samples == 4000, "fingerprint reads the capped rows")
fp40 <- SegmentFingerprint(prod)$md5
Sys.setenv(CET_MAX_SEG = "41")
Check(SegmentFingerprint(prod)$md5 != fp40, "a different cap gives a different fingerprint")
Sys.unsetenv("CET_MAX_SEG_RUNG")
Check(length(SegFiles(rung)) == 23, "CET_MAX_SEG_RUNG defaults to CET_MAX_SEG")
Sys.setenv(CET_MAX_SEG = "x")
Check(inherits(try(SegCap("b0g0_run1"), silent = TRUE), "try-error"), "a non-numeric cap stops")
unlink(root, recursive = TRUE)
cat("check-segcap.R: all checks passed\n")
