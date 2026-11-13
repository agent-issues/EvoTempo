# Trilobite CET pipeline: checks of bayes-factors.R eval and combine on tiny
# real chains with a ladder rung present (#100), of combine's expected chain
# set (#101), and of run-chain.R's start-tree resume guard (#102).
#
#   Rscript data-raw/empirical/trilobites/cet/check-combine.R <check root> [<old bayes-factors.R>]
#
# Stage 3 of check-combine.sh, which ran the chains (<check root>/runs: cells
# b0g0, b0g1, b1g0, b1g1 and the rung b0.5g0, replicates 1-2) and eval on each.
# Checks:
# 1. eval exited 0 on all 10 chains and wrote cross.csv with the columns
#    aug_b0, aug_b0.5, aug_b1, chr_g0, chr_g1, its checks within CET_BF_TOL.
# 2. With launch.csv listing the 5 labels x 2 replicates, combine exits 0,
#    reports 4 edges, and bridges b0g0 -> b1g0 through the rung (3 points).
# 3. One cross.csv moved away: combine stops, naming the chain, and writes no
#    bayes-factors.csv; with CET_BF_ALLOW_MISSING naming it, exit 0 and the
#    chain is reported as left out.  Without launch.csv (production default 4
#    cells x 4 replicates) combine stops: run3/run4 missing, rung not listed.
# 4. If the old script is given (main before #100/#101): its eval wrote
#    aug_b0.0 / aug_b0.5 / aug_b1.0 and its combine crashes on them (#100 as
#    reproduced in the #91 delta review); on the 4 cells without the rung,
#    with one cross.csv removed, it exits 0 with a different log Bayes factor
#    (#101, silent drop), where the new combine stops.
# 5. Resume guard (#102): rerunning run-chain.R on b0g0_run1 with
#    CET_START_TREE a copy of its start tree at another absolute path
#    continues (exit 0, "resuming"); with another draw's tree it stops.
# Writes <check root>/check-combine.txt; stops if a check fails.
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
oldBF <- if (length(args) > 1 && file.exists(args[[2]])) args[[2]]
runs <- file.path(root, "runs")
tol <- as.numeric(Sys.getenv("CET_BF_TOL", "1e-6"))
lines <- character()
Say <- function(...) {
  s <- sprintf(...)
  cat(s, "\n", sep = "")
  lines <<- c(lines, s)
}
fail <- character()
Check <- function(ok, what) {
  Say("%s %s", if (isTRUE(ok)) "ok  " else "FAIL", what)
  if (!isTRUE(ok)) fail <<- c(fail, what)
}
rscript <- file.path(R.home("bin"), "Rscript")
Run <- function(script, args, env = character()) {
  out <- suppressWarnings(system2(rscript, c(script, args), env = env, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(out = as.character(out), status = if (is.null(status)) 0L else status)
}
newBF <- file.path(cetDir, "bayes-factors.R")
Combine <- function(r, script = newBF, env = character()) {
  unlink(file.path(r, "bayes-factors.csv"))
  res <- Run(script, c("combine", r), c("CET_MKV=pattern", env))
  res$edges <- if (file.exists(file.path(r, "bayes-factors.csv"))) {
    utils::read.csv(file.path(r, "bayes-factors.csv"))
  }
  res
}
Show <- function(res, pattern = "Error|Expected|Left out|no cross|not in|no complete|column|Edge|missing value") {
  paste(grep(pattern, res$out, value = TRUE), collapse = " / ")
}
EdgeText <- function(e) {
  if (is.null(e)) "none" else paste(sprintf("%s->%s %.4f (%d points)", e$from, e$to, e$logBF, e$points),
                                    collapse = ", ")
}
cells4 <- c("b0g0", "b0g1", "b1g0", "b1g1")
labs <- c(cells4, "b0.5g0")
chains <- as.vector(outer(labs, 1:2, function(l, r) sprintf("%s_run%d", l, r)))
LastExit <- function(f) {
  if (!file.exists(f)) return(NA_integer_)
  x <- grep("^exit [0-9]+$", readLines(f), value = TRUE)
  if (length(x)) as.integer(sub("exit ", "", utils::tail(x, 1))) else NA_integer_
}

# ---- 1. eval on every chain ---------------------------------------------------------
ev <- do.call(rbind, lapply(chains, function(ch) {
  f <- file.path(runs, ch, "cross.csv")
  x <- if (file.exists(f)) utils::read.csv(f, check.names = FALSE)
  data.frame(chain = ch, exit = LastExit(file.path(root, "logs", paste0("eval-", ch, ".txt"))),
             cross = !is.null(x), states = if (is.null(x)) NA else nrow(x),
             columns = if (is.null(x)) "" else paste(grep("^(aug|chr)_", names(x), value = TRUE),
                                                     collapse = " "),
             maxCheck = if (is.null(x)) NA else max(abs(as.matrix(x[, c("check", "checkChar")]))))
}))
Say("eval results:")
for (k in seq_len(nrow(ev))) {
  Say("  %-12s exit %s, cross.csv %s, %s states, max |check| %.3g, columns %s", ev$chain[k],
      format(ev$exit[k]), ev$cross[k], format(ev$states[k]), ev$maxCheck[k], ev$columns[k])
}
want <- "aug_b0 aug_b0.5 aug_b1 chr_g0 chr_g1"
Check(all(ev$exit %in% 0L) && all(ev$cross), "eval exited 0 and wrote cross.csv on all 10 chains")
Check(all(ev$columns == want), sprintf("every cross.csv has the columns %s (#100)", want))
Check(all(ev$maxCheck <= tol), sprintf("every cross.csv's checks are within %g", tol))

# ---- 2. combine with the rung ----------------------------------------------------
utils::write.csv(data.frame(cell = labs, reps = 2), file.path(runs, "launch.csv"),
                 row.names = FALSE, quote = FALSE)
full <- Combine(runs)
Say("combine, rung present, exit %d: %s", full$status, Show(full, "Expected|Error"))
Say("  edges: %s", EdgeText(full$edges))
e <- full$edges
Check(full$status == 0L && !is.null(e) && nrow(e) == 4L && all(is.finite(e$logBF)),
      "combine with the rung present exits 0 with 4 finite edges (#100; before: crash)")
Check(!is.null(e) && any(e$from == "b0g0" & e$to == "b1g0" & e$points == 3L),
      "b0g0 -> b1g0 is bridged through the rung b0.5g0 (3 points)")

# ---- 3. a missing cross.csv --------------------------------------------------------
gone <- file.path(runs, "b0g1_run2", "cross.csv")
invisible(file.rename(gone, paste0(gone, ".bak")))
res <- Combine(runs)
Say("b0g1_run2/cross.csv moved away, exit %d: %s", res$status, Show(res))
Check(res$status != 0L && any(grepl("no cross.csv.*b0g1_run2", res$out)) && is.null(res$edges),
      "combine stops on a missing cross.csv, names b0g1_run2, writes no bayes-factors.csv (#101)")
res <- Combine(runs, env = "CET_BF_ALLOW_MISSING=b0g1_run2")
Say("... with CET_BF_ALLOW_MISSING=b0g1_run2, exit %d: %s; edges %s", res$status,
    Show(res, "Left out|Error"), EdgeText(res$edges))
Check(res$status == 0L && any(grepl("Left out on request.*b0g1_run2", res$out)) &&
        !is.null(res$edges) && nrow(res$edges) == 4L,
      "explicit opt-out: exit 0, b0g1_run2 reported as left out")
invisible(file.rename(paste0(gone, ".bak"), gone))
invisible(file.rename(file.path(runs, "launch.csv"), file.path(runs, "launch.bak")))
res <- Combine(runs)
Say("no launch.csv (production default), exit %d: %s", res$status, Show(res))
Check(res$status != 0L && any(grepl("production default", res$out)) &&
        any(grepl("no complete segment.*b0g0_run3", res$out)) &&
        any(grepl("not in the launch record.*b0.5g0_run1", res$out)),
      "without launch.csv: stops, run3/run4 of the 4 x 4 default missing and the rung not listed")
invisible(file.rename(file.path(runs, "launch.bak"), file.path(runs, "launch.csv")))

# ---- 4. The old script, for comparison --------------------------------------------
if (!is.null(oldBF)) {
  old <- file.path(root, "old-runs")
  cols <- vapply(chains, function(ch) {
    f <- file.path(old, ch, "cross.csv")
    if (!file.exists(f)) return("(no cross.csv)")
    paste(grep("^aug_", names(utils::read.csv(f, check.names = FALSE)), value = TRUE), collapse = " ")
  }, "")
  Say("old eval columns: %s", paste(unique(cols), collapse = " | "))
  Check(all(cols == "aug_b0.0 aug_b0.5 aug_b1.0"),
        "old eval named the columns aug_b0.0 / aug_b0.5 / aug_b1.0 (the #100 fault)")
  res <- Combine(old, oldBF)
  Say("old combine, rung present, exit %d: %s", res$status, Show(res, "Error|missing value"))
  Check(res$status != 0L, "old combine fails with the rung present (reproduces #100)")

  # #101 before and after, on the 4 cells without the rung (new cross.csv files)
  r4 <- file.path(root, "norung")
  unlink(r4, recursive = TRUE)
  dir.create(r4)
  for (ch in chains[!grepl("^b0.5g0", chains)]) file.copy(file.path(runs, ch), r4, recursive = TRUE)
  before <- Combine(r4, oldBF)
  file.remove(file.path(r4, "b0g0_run1", "cross.csv"))
  after <- Combine(r4, oldBF)
  Edge01 <- function(e) if (is.null(e)) NA else e$logBF[e$from == "b0g0" & e$to == "b0g1"]
  Say("old combine, 4 cells x 2: complete exit %d, b0g0->b0g1 %.4f; without b0g0_run1/cross.csv exit %d, %.4f",
      before$status, Edge01(before$edges), after$status, Edge01(after$edges))
  Check(after$status == 0L && isTRUE(Edge01(before$edges) != Edge01(after$edges)),
        "old combine drops the chain silently: exit 0, b0g0 -> b0g1 changes (reproduces #101)")
  utils::write.csv(data.frame(cell = cells4, reps = 2), file.path(r4, "launch.csv"),
                   row.names = FALSE, quote = FALSE)
  res <- Combine(r4)
  Say("new combine on the same root, exit %d: %s", res$status, Show(res))
  Check(res$status != 0L && any(grepl("no cross.csv.*b0g0_run1", res$out)),
        "new combine stops on the same root and names b0g0_run1")
}

# ---- 5. Resume guard (#102) -------------------------------------------------------
ckp <- readRDS(file.path(runs, "b0g0_run1", "checkpoint.rds"))
Say("b0g0_run1 checkpoint: startTreeRel %s, md5 %s", format(ckp$startTreeRel), format(ckp$startTreeMd5))
Check(identical(ckp$startTreeRel, "data/start-tree-01.nwk") &&
        identical(ckp$startTreeMd5, unname(tools::md5sum(StartTreeFile(1L)))),
      "the checkpoint records the start tree relative to cet/ and its md5")
copyDir <- file.path(root, "tree-copy")
dir.create(copyDir, showWarnings = FALSE)
invisible(file.copy(StartTreeFile(1L), copyDir, overwrite = TRUE))
env <- c("CET_NITER=40", "CET_BURNIN=10", "CET_THIN=10", "CET_SEG=40", paste0("CET_OUT=", runs))
runChain <- file.path(cetDir, "run-chain.R")
res <- Run(runChain, c("1", "1"), c(env, paste0("CET_START_TREE=", file.path(copyDir, "start-tree-01.nwk"))))
prog <- readLines(file.path(runs, "b0g0_run1", "progress.txt"))
Say("resume with a copy of start-tree-01 at %s: exit %d; progress: %s", copyDir, res$status,
    utils::tail(prog, 2)[1])
Check(res$status == 0L && any(grepl("resuming after segment 1", prog)),
      "same tree (same md5) at another absolute path: the chain resumes (#102; before: stop)")
res <- Run(runChain, c("1", "1"), c(env, paste0("CET_START_TREE=", normalizePath(StartTreeFile(2L)))))
Say("resume with start-tree-02: exit %d: %s", res$status,
    paste(grep("Error|start tree", res$out, value = TRUE), collapse = " / "))
Check(res$status != 0L && any(grepl("a different start tree", res$out)),
      "another start tree (other md5): the resume stops")

writeLines(lines, file.path(root, "check-combine.txt"))
if (length(fail)) stop("check-combine.R FAILED: ", paste(fail, collapse = "; "))
cat("check-combine.R: all checks passed\n")
