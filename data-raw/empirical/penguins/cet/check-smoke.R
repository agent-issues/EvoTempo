# Penguin CET pipeline: checks on a smoke run (smoke.sh).
#
#   Rscript data-raw/empirical/penguins/cet/check-smoke.R <W>
#
# Reports seconds per iteration per cell, acceptance rates, and whether the
# resumed chain (two segments) equals the uninterrupted one (one segment).
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args)) args[[1]] else "/nobackup/pjjg18/evotempo-penguin-cet"
labels <- c("b0g0", "b0g1", "b1g0", "b1g1")
timing <- do.call(rbind, lapply(labels, function(l) {
  acc <- read.csv(file.path(W, "smoke", paste0(l, "_run1"), "acceptance.csv"))
  seg <- unique(acc[, c("segment", "secs", "burnin")])
  smp <- do.call(rbind, lapply(sort(Sys.glob(file.path(W, "smoke", paste0(l, "_run1"),
                                                        "seg*_samples.csv"))), read.csv))
  data.frame(cell = l, segment = seg$segment, secs = seg$secs,
             iterations = 100 + seg$burnin,
             secsPerIteration = seg$secs / (100 + seg$burnin),
             finalLogPost = tapply(smp$logPost, smp$segment, utils::tail, 1)[seg$segment],
             finalCharLL = tapply(smp$charLL, smp$segment, utils::tail, 1)[seg$segment])
}))
print(timing, row.names = FALSE)
write.csv(timing, file.path(W, "smoke", "timing.csv"), row.names = FALSE)
acc <- do.call(rbind, lapply(labels, function(l) {
  a <- read.csv(file.path(W, "smoke", paste0(l, "_run1"), "acceptance.csv"))
  cbind(cell = l, a)
}))
wide <- reshape(acc[acc$segment == 2, c("cell", "move", "rate")], idvar = "move",
                timevar = "cell", direction = "wide")
print(wide, row.names = FALSE, digits = 3)
write.csv(acc, file.path(W, "smoke", "acceptance.csv"), row.names = FALSE)
a <- read.csv(file.path(W, "smoke", "b1g1_run1", "seg0002_samples.csv"))
b <- read.csv(file.path(W, "smoke-identity", "b1g1_run1", "seg0001_samples.csv"))
b <- b[b$iteration > 100, ]
cols <- setdiff(intersect(names(a), names(b)), c("iteration", "segment"))
diff <- max(abs(as.matrix(a[, cols]) - as.matrix(b[, cols])), na.rm = TRUE)
ta <- readLines(file.path(W, "smoke", "b1g1_run1", "seg0002_trees.nwk"))
tb <- utils::tail(readLines(file.path(W, "smoke-identity", "b1g1_run1", "seg0001_trees.nwk")), 10)
cat(sprintf("Resume check: max |difference| over %d samples x %d columns = %g; trees identical: %s\n",
            nrow(a), length(cols), diff, identical(ta, tb)))
