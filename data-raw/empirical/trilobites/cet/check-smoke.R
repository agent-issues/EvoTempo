# Trilobite CET pipeline: checks on a smoke run (smoke.sh).
#
#   [CET_MKV=<variant>] Rscript data-raw/empirical/trilobites/cet/check-smoke.R <W> [out dir]
#
# Reports seconds per iteration per cell and segment, acceptance rates (and
# tree proposals rejected by the ingroup constraint), and whether the resumed
# chain (two segments) equals the uninterrupted one (one segment); stops if it
# does not.  Writes timing.csv, acceptance.csv and resume.txt to [out dir]
# (default <W>/smoke).
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args)) args[[1]] else "/nobackup/pjjg18/evotempo-trilobite-cet"
outDir <- if (length(args) > 1) args[[2]] else file.path(W, "smoke")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
mkv <- Sys.getenv("CET_MKV", "pattern")
sfx <- if (mkv == "pattern") "" else paste0("_mkv", mkv)
cells <- c("b0g0", "b0g1", "b1g0", "b1g1")
labels <- paste0(cells, sfx)
timing <- do.call(rbind, lapply(labels, function(l) {
  d <- file.path(W, "smoke", paste0(l, "_run1"))
  acc <- read.csv(file.path(d, "acceptance.csv"))
  seg <- unique(acc[, c("segment", "secs", "burnin")])
  smp <- do.call(rbind, lapply(sort(Sys.glob(file.path(d, "seg*_samples.csv"))), read.csv))
  data.frame(cell = cells[match(l, labels)], segment = seg$segment, secs = seg$secs,
             iterations = 100 + seg$burnin, secsPerIteration = seg$secs / (100 + seg$burnin),
             finalLogPost = tapply(smp$logPost, smp$segment, utils::tail, 1)[seg$segment],
             finalCharLL = tapply(smp$charLL, smp$segment, utils::tail, 1)[seg$segment],
             finalOriginMa = tapply(smp$originAge, smp$segment, utils::tail, 1)[seg$segment] + 485.4)
}))
print(timing, row.names = FALSE)
write.csv(timing, file.path(outDir, "timing.csv"), row.names = FALSE)
acc <- do.call(rbind, lapply(labels, function(l) {
  a <- read.csv(file.path(W, "smoke", paste0(l, "_run1"), "acceptance.csv"))
  cbind(cell = cells[match(l, labels)], a)
}))
wide <- reshape(acc[acc$segment == 2, c("cell", "move", "rate")], idvar = "move",
                timevar = "cell", direction = "wide")
print(wide, row.names = FALSE, digits = 3)
cr <- aggregate(constraintRejections ~ cell + move, acc, sum)
print(cr[cr$constraintRejections > 0, ], row.names = FALSE)
write.csv(acc, file.path(outDir, "acceptance.csv"), row.names = FALSE)
b11 <- paste0("b1g1", sfx, "_run1")
a <- read.csv(file.path(W, "smoke", b11, "seg0002_samples.csv"))
b <- read.csv(file.path(W, "smoke-identity", b11, "seg0001_samples.csv"))
b <- b[b$iteration > 100, ]
cols <- setdiff(intersect(names(a), names(b)), c("iteration", "segment"))
diff <- max(abs(as.matrix(a[, cols]) - as.matrix(b[, cols])), na.rm = TRUE)
ta <- readLines(file.path(W, "smoke", b11, "seg0002_trees.nwk"))
tb <- utils::tail(readLines(file.path(W, "smoke-identity", b11, "seg0001_trees.nwk")), 10)
msg <- sprintf("Resume check: max |difference| over %d samples x %d columns = %g; trees identical: %s",
               nrow(a), length(cols), diff, identical(ta, tb))
cat(msg, "\n")
writeLines(msg, file.path(outDir, "resume.txt"))
if (!(diff == 0 && identical(ta, tb))) stop("resumed chain differs from the uninterrupted one")
