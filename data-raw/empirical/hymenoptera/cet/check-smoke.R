# Hymenoptera CET pipeline: checks on a smoke run (smoke.sh), as the penguin
# pipeline's check-smoke.R, without its fixed segment length.
#
#   [CET_MKV=<variant>] Rscript data-raw/empirical/hymenoptera/cet/check-smoke.R <W> [smoke dir]
#
# Reports seconds per iteration per cell and segment, acceptance rates
# (including the FBD block moves), and whether the resumed chain (two
# segments, <smoke dir>) equals the uninterrupted one (one segment,
# <smoke dir>-identity).  Writes <smoke dir>/timing.csv and acceptance.csv.
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
             finalLambda = tapply(smp$lambda, smp$segment, utils::tail, 1)[as.character(seg$segment)])
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
cat(sprintf("Resume check: %d samples x %d columns, max |difference| = %g; trees identical: %s\n",
            nrow(a), length(cols), diff, identical(ta, tb)))
