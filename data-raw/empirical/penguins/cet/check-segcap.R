# Penguin CET pipeline: checks of the segment cap (CET_MAX_SEG,
# CET_MAX_SEG_RUNG; common.R, ParseCap(), SegCap(), SegFiles(),
# ChainComplete(); #107 review R107-9, #105 delta review C1, C2, C4) on
# synthetic chain directories.  C3 (the cap in eval's record, checked by
# combine) is in check-labels.R section 8.
#
#   Rscript data-raw/empirical/penguins/cet/check-segcap.R
#
# Stops at the first failed check.
source("data-raw/empirical/penguins/cet/common.R")
root <- tempfile("segcap")
MakeChain <- function(label, nSeg, ckpSegment = nSeg, complete = TRUE) {
  d <- file.path(root, label)
  dir.create(d, recursive = TRUE)
  for (k in seq_len(nSeg)) {
    utils::write.csv(data.frame(logPost = rnorm(100), charLL = rnorm(100), iteration = (k - 1) * 100 + 1:100),
                     file.path(d, sprintf("seg%04d_samples.csv", k)), row.names = FALSE)
    saveRDS(as.list(seq_len(100)), file.path(d, sprintf("seg%04d_states.rds", k)))
  }
  # run-chain.R writes the checkpoint atomically after a segment's files
  saveRDS(list(segment = ckpSegment), file.path(d, "checkpoint.rds"))
  writeLines(c(sprintf("x segment %d/%d done", nSeg, nSeg), if (complete) "x complete" else "x resuming"),
             file.path(d, "progress.txt"))
  d
}
prod <- MakeChain("b0g0_run1", 46, complete = FALSE)
rung <- MakeChain("b0g0.25_run1", 23, complete = FALSE)
ladder <- MakeChain("b0g0_t0.5_run1", 12, complete = TRUE)
# Segment 40 of this chain is still being written: its samples file exists,
# the checkpoint records segment 39
writing <- MakeChain("b1g0_run1", 40, ckpSegment = 39, complete = FALSE)
gap <- MakeChain("b1g1_run1", 46, complete = FALSE)
unlink(file.path(gap, "seg0001_samples.csv"))
Check <- function(x, msg) {
  if (!isTRUE(x)) stop("FAILED: ", msg)
  cat("ok  ", msg, "\n")
}
Stops <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
Sys.unsetenv(c("CET_MAX_SEG", "CET_MAX_SEG_RUNG"))
Check(length(SegFiles(prod)) == 46, "no cap: every segment")
Check(!ChainComplete(prod) && ChainComplete(ladder), "no cap: completeness from progress.txt")
Sys.setenv(CET_MAX_SEG = "40", CET_MAX_SEG_RUNG = "20")
Check(length(SegFiles(prod)) == 40 && length(SegFiles(prod, "states.rds")) == 40, "production cell capped at CET_MAX_SEG")
Check(length(SegFiles(rung)) == 20, "rung capped at CET_MAX_SEG_RUNG")
Check(identical(basename(SegFiles(rung)), sprintf("seg%04d_samples.csv", 1:20)), "the cap reads segments 1..m")
Check(ChainComplete(prod) && ChainComplete(rung), "capped: complete when the checkpoint has reached the cap")
Check(!ChainComplete(ladder) && Stops(SegFiles(ladder)),
      "capped: a power chain with 12 < 20 segments is not complete, and reading it stops")
Check(!ChainComplete(writing),
      "C1: segment 40 written but not checkpointed (checkpoint at 39): not complete under cap 40")
Check(Stops(SegFiles(gap)) && !ChainComplete(gap),
      "C2: seg0001 missing: reading under the cap stops, and the chain is not complete (before: segments 2..41)")
Check(Stops(SegmentFingerprint(gap)), "C2: the fingerprint stops too")
Check(SegmentFingerprint(prod)$samples == 4000, "fingerprint reads the capped rows")
fp40 <- SegmentFingerprint(prod)$md5
Sys.setenv(CET_MAX_SEG = "41")
Check(SegmentFingerprint(prod)$md5 != fp40, "a different cap gives a different fingerprint")
Sys.unsetenv("CET_MAX_SEG_RUNG")
Check(SegCap("b0g0.25_run1") == 41L && Stops(SegFiles(rung)),
      "CET_MAX_SEG_RUNG defaults to CET_MAX_SEG (41 > the rung's 23 segments: reading stops)")
for (bad in c("40.7", "-1", "0", "4e1", "x", "40 x")) {
  Sys.setenv(CET_MAX_SEG = bad)
  Check(Stops(SegCap("b0g0_run1")), sprintf("C4: CET_MAX_SEG = '%s' stops", bad))
}
Sys.setenv(CET_MAX_SEG = "40", CET_MAX_SEG_RUNG = "2.5")
Check(Stops(SegCap("b0g0.25_run1")) && SegCap("b0g0_run1") == 40L,
      "C4: CET_MAX_SEG_RUNG = '2.5' stops for a rung (a production cell uses CET_MAX_SEG)")
Sys.setenv(CET_MAX_SEG = "", CET_MAX_SEG_RUNG = "")
Check(SegCap("b0g0_run1") == 0L && SegCap("b0g0.25_run1") == 0L, "empty means no cap")
unlink(root, recursive = TRUE)
cat("check-segcap.R: all checks passed\n")
