# Trilobite CET pipeline: checks on a smoke run (smoke.sh).
#
#   [CET_MKV=<variant>] Rscript data-raw/empirical/trilobites/cet/check-smoke.R <W> [out dir]
#
# Reports seconds per iteration per cell and segment, acceptance rates (and
# tree proposals rejected by the ingroup constraint), and whether the resumed
# chain (two segments) equals the uninterrupted one (one segment); stops if it
# does not.  Then runs `bayes-factors.R eval` on every smoke chain, including
# the replicate-2 chains (smoke.sh tasks 6-7, whose start trees have other tip
# orders), and stops unless every eval exits 0 and writes cross.csv (review of
# #91, A3).  Run from the package root.  Writes timing.csv, acceptance.csv,
# resume.txt and eval.txt to [out dir] (default <W>/smoke).
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
             finalAlpha = tapply(smp$alpha, smp$segment, utils::tail, 1)[seg$segment],
             finalOriginMa = tapply(smp$originAge, smp$segment, utils::tail, 1)[seg$segment] + 485.4)
}))
print(timing, row.names = FALSE)
write.csv(timing, file.path(outDir, "timing.csv"), row.names = FALSE)
# Peak memory per smoke process (/usr/bin/time -v in smoke.sh's .err logs;
# review of #114, R114-5)
mem <- do.call(rbind, lapply(sort(Sys.glob(file.path(W, "out", "smoke*.err"))), function(f) {
  kb <- as.numeric(sub(".*: ", "", grep("Maximum resident set size", readLines(f), value = TRUE)))
  if (length(kb)) data.frame(log = basename(f), process = seq_along(kb), peakMB = round(kb / 1024, 1))
}))
print(mem, row.names = FALSE)
write.csv(mem, file.path(outDir, "peak-memory.csv"), row.names = FALSE)
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

# ---- Review of #114: resume in a second process, predictive stop, kill (tasks 8-10) ----
ReadChain <- function(d) {
  list(samples = do.call(rbind, lapply(sort(Sys.glob(file.path(d, "seg*_samples.csv"))), read.csv)),
       trees = unlist(lapply(sort(Sys.glob(file.path(d, "seg*_trees.nwk"))), readLines)),
       progress = readLines(file.path(d, "progress.txt")),
       acc = read.csv(file.path(d, "acceptance.csv")))
}
# Max |difference| over every column but `segment` (which differs only when
# the segment lengths do), and whether the trees are identical
Compare <- function(a, b) {
  cols <- setdiff(intersect(names(a$samples), names(b$samples)), "segment")
  ok <- nrow(a$samples) == nrow(b$samples)
  d <- if (ok) max(abs(as.matrix(a$samples[, cols]) - as.matrix(b$samples[, cols])), na.rm = TRUE) else Inf
  list(diff = d, trees = identical(a$trees, b$trees), n = nrow(a$samples), cols = length(cols))
}
Has <- function(ch, pattern) any(grepl(pattern, ch$progress))
Complete <- function(ch) grepl(" complete$", trimws(utils::tail(ch$progress, 1)))
AccOnce <- function(ch, nSeg) {
  !anyDuplicated(ch$acc[, c("segment", "move")]) && identical(sort(unique(ch$acc$segment)), seq_len(nSeg))
}
lab4 <- paste0("b1g1", sfx, "_run1")
lab1 <- paste0("b0g0", sfx, "_run1")
two <- ReadChain(file.path(W, "smoke-twoproc", lab4))
inProc <- ReadChain(file.path(W, "smoke", lab4))
oneSeg <- ReadChain(file.path(W, "smoke-identity", lab4))
pre <- ReadChain(file.path(W, "smoke-predict", lab1))
kil <- ReadChain(file.path(W, "smoke-kill", lab1))
c1 <- Compare(two, inProc)
c2 <- Compare(two, oneSeg)
c3 <- Compare(kil, pre)
checks <- c(
  "task 8: first process stopped at CET_HOURS = 0 after segment 1" = Has(two, "segment 1/2 done") &&
    Has(two, "stopping after .* \\(CET_HOURS\\)"),
  "task 8: second process resumed from the checkpoint on disk and completed" =
    Has(two, "resuming after segment 1") && Complete(two),
  "task 8: = task 4 (two segments, one process), every column and tree" = c1$diff == 0 && c1$trees,
  "task 8: = task 5 (one segment), every column but segment, and every tree" = c2$diff == 0 && c2$trees,
  "task 9: predictive stop before segment 2 under the fake wall limit" =
    isTRUE(grep("stopping before segment 2: .* would pass", pre$progress)[1] <
             grep("segment 2/4 done", pre$progress)[1]) && Has(pre, "segment 1/4 done"),
  "task 9: only the 1.5 factor stopped it (elapsed + 1.0 x last segment < limit)" = local({
    m <- regmatches(pre$progress, regexec(
      "stopping before segment 2: ([0-9.]+) h elapsed \\+ 1\\.5 x ([0-9.]+) h .* would pass ([0-9.]+) h",
      pre$progress))
    m <- Filter(length, m)
    length(m) == 1L && {
      v <- as.numeric(m[[1]][2:4])
      v[1] + v[2] < v[3] && v[1] + 1.5 * v[2] > v[3]
    }
  }),
  "task 9: resumed and completed (4 segments)" = Has(pre, "resuming after segment 1") && Complete(pre) &&
    max(pre$samples$segment) == 4,
  "task 10: killed after segment 2's acceptance.csv, before its checkpoint" =
    Has(kil, "test kill after segment 2") && Has(kil, "resuming after segment 1") && Complete(kil),
  "task 10: acceptance.csv holds each segment once" = AccOnce(kil, 4),
  "task 9: acceptance.csv holds each segment once" = AccOnce(pre, 4),
  "task 10: = task 9 (other interruption), every column and tree" = c3$diff == 0 && c3$trees)
stopMsg <- c(
  sprintf("task 8 vs task 4: max |difference| %g over %d samples x %d columns; trees identical %s",
          c1$diff, c1$n, c1$cols, c1$trees),
  sprintf("task 8 vs task 5: max |difference| %g over %d samples x %d columns; trees identical %s",
          c2$diff, c2$n, c2$cols, c2$trees),
  sprintf("task 10 vs task 9: max |difference| %g over %d samples x %d columns; trees identical %s",
          c3$diff, c3$n, c3$cols, c3$trees),
  "task 9 progress:", paste(" ", pre$progress[-1]),
  "task 10 progress:", paste(" ", kil$progress[-1]),
  "task 8 progress:", paste(" ", two$progress[-1]),
  paste(ifelse(checks, "ok  ", "FAIL"), names(checks)))
writeLines(stopMsg)
writeLines(stopMsg, file.path(outDir, "stop-resume.txt"))
if (!all(checks)) stop("stop/resume checks failed: ", paste(names(checks)[!checks], collapse = "; "))

# eval on every smoke chain (the start-tree path of bayes-factors.R)
dirs <- sort(Sys.glob(file.path(W, "smoke", paste0("b*", sfx, "_run*"))))
stopifnot(any(grepl("_run2$", dirs)))
ev <- do.call(rbind, lapply(dirs, function(d) {
  unlink(file.path(d, "cross.csv"))
  out <- suppressWarnings(system2("Rscript", c("data-raw/empirical/trilobites/cet/bayes-factors.R",
                                               "eval", d), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  data.frame(chain = basename(d), exit = if (is.null(status)) 0L else status,
             crossWritten = file.exists(file.path(d, "cross.csv")),
             startTree = sub("^.*start tree (\\S+) \\(draw ([0-9]+)\\).*$", "draw \\2",
                             grep("start tree", readLines(file.path(d, "progress.txt")), value = TRUE)[1]),
             message = utils::tail(grep("states", out, value = TRUE), 1))
}))
print(ev, row.names = FALSE)
writeLines(c(msg, utils::capture.output(print(ev, row.names = FALSE))), file.path(outDir, "eval.txt"))
if (any(ev$exit != 0 | !ev$crossWritten)) stop("bayes-factors.R eval failed on a smoke chain")
