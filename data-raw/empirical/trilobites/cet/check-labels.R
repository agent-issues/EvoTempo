# Trilobite CET pipeline: checks of the chain labels, of ladder-rung column
# names, and of the expected chain set in bayes-factors.R combine (#100, #101).
#
#   Rscript data-raw/empirical/trilobites/cet/check-labels.R [out.txt]
#
# Run from the package root with EvoTempo installed (or CET_LOAD_ALL=1 for
# pkgload::load_all()); needs no chain output (synthetic chain directories in
# tempdir()).  check-combine.sh runs eval and combine on real chains.  Checks:
# 1. FormatValue() writes rung values and powers in plain decimal notation,
#    one number at a time (c(0, 0.5, 1) -> "0", "0.5", "1", where format() on
#    the vector gave "0.0", "0.5", "1.0"; #100), ParseLabel() reads them back
#    exactly, and the production labels b0g0 ... b1g1 are unchanged.
# 2. ParseLabel() stops on a label it cannot read, and Chains() stops on a
#    chain directory with such a label, while it ignores summary/.
# 3. combine with a ladder rung b0.5g0 present (cross.csv columns aug_b0,
#    aug_b0.5, aug_b1 as eval now writes them) exits 0 and bridges
#    b0g0 -> b1g0 through 3 points (#100); with the old column names
#    (aug_b0.0, aug_b0.5, aug_b1.0) it stops and says which column is missing.
# 4. Expected chain set (#101): with launch.csv listing 4 cells x 2 replicates,
#    combine exits 0; after one cross.csv is removed it stops, names that
#    chain and writes no bayes-factors.csv; with CET_BF_ALLOW_MISSING naming
#    it, it exits 0 and says it left the chain out; a chain missing entirely,
#    a chain not in the record, an opt-out naming an unexpected chain, and
#    the production default (4 x 4) against 4 x 2 chains all stop.
# 5. combine stops on a run root that does not exist.
# Stops with an error if a check fails.
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
outFile <- if (length(args)) args[[1]] else NULL
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
Errors <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

# ---- 1. Labels in plain decimal notation, read back exactly -------------------
rungs <- c(0, 0.25, 0.5, 0.75, 1)
Check(identical(paste0("aug_b", format(c(0, 0.5, 1))), c("aug_b0.0", "aug_b0.5", "aug_b1.0")),
      "(the #100 fault) format(c(0, 0.5, 1)) pads to aug_b0.0 / aug_b0.5 / aug_b1.0")
Check(identical(paste0("aug_b", vapply(c(0, 0.5, 1), FormatValue, "")),
                c("aug_b0", "aug_b0.5", "aug_b1")),
      "FormatValue() one at a time gives aug_b0 / aug_b0.5 / aug_b1")
powers <- c((1:15 / 16)^(1 / 0.3), 1e-05, 1e-04, 1 / 3)
new <- vapply(powers, FormatValue, "")
back <- vapply(new, function(s) ParseLabel(sprintf("b1g1_mkvnone_t%s_run1", s))$power, 0)
Check(!any(grepl("e", new)) && all(back == powers),
      sprintf("FormatValue() writes %d powers without scientific notation; ParseLabel() reads them back exactly",
              length(powers)))
Check(identical(vapply(seq_len(4), function(i) CellLabel(cellTable$beta[i], cellTable$gamma[i]), ""),
                  c("b0g0", "b0g1", "b1g0", "b1g1")),
      "production labels unchanged (b0g0, b0g1, b1g0, b1g1)")
lab <- ParseLabel(paste0(CellLabel(0.5, 0), "_t0.0001_run3"))
Check(identical(CellLabel(0.5, 0), "b0.5g0") && lab$beta == 0.5 && lab$gamma == 0 &&
        lab$power == 1e-04 && lab$rep == 3L && lab$mkv == "pattern",
      "ParseLabel() reads a ladder-rung label (b0.5g0_t0.0001_run3)")

# ---- 2. Unreadable labels stop -----------------------------------------------------
Check(Errors(ParseLabel("b1g1_t1e-04_run1")), "ParseLabel() stops on b1g1_t1e-04_run1")
Check(Errors(ParseLabel("b0.5g0_mkvallTips_run1")), "ParseLabel() stops on an unknown Mkv variant")
root1 <- tempfile("labels")
dir.create(file.path(root1, "summary"), recursive = TRUE)
writeLines("x", file.path(root1, "summary", "clades.csv"))
dir.create(file.path(root1, "b1g1_run1"))
utils::write.csv(data.frame(logPost = 0, charLL = 0), file.path(root1, "b1g1_run1", "seg0001_samples.csv"))
ok <- Chains(root1)
Check(nrow(ok) == 1L && ok$label == "b1g1_run1", "Chains() ignores summary/ and keeps b1g1_run1")
dir.create(file.path(root1, "b1g1_t1e-04_run1"))
writeLines("start", file.path(root1, "b1g1_t1e-04_run1", "progress.txt"))
msg <- tryCatch({
  Chains(root1)
  ""
}, error = conditionMessage)
Check(grepl("b1g1_t1e-04_run1", msg), "Chains() stops on a chain directory it cannot parse")
unlink(root1, recursive = TRUE)

# ---- Synthetic chains -------------------------------------------------------------
# seg0001_samples.csv and a cross.csv as eval writes it (checks 0), with
# log q(beta, gamma) = -4 beta - 2 gamma plus a little noise, so the bridge log
# BF is about -4 per unit beta and -2 per unit gamma
nSmp <- 40L
set.seed(10101)
MakeChain <- function(root, label, rep, betas = c(0, 1), gammas = c(0, 1), oldNames = FALSE,
                      cross = TRUE) {
  d <- file.path(root, sprintf("%s_run%d", label, rep))
  dir.create(d, recursive = TRUE)
  ll <- -1500 + stats::rnorm(nSmp)
  utils::write.csv(data.frame(logPost = ll, charLL = ll, iteration = seq_len(nSmp)),
                   file.path(d, "seg0001_samples.csv"), row.names = FALSE)
  if (!cross) return(invisible(d))
  Name <- function(v) if (oldNames) format(v) else vapply(v, FormatValue, "")
  z <- stats::rnorm(nSmp)
  z2 <- stats::rnorm(nSmp)
  x <- data.frame(sample = seq_len(nSmp))
  for (k in seq_along(betas)) x[[paste0("aug_b", Name(betas)[k])]] <- -4 * betas[k] + 0.05 * betas[k] * z
  for (k in seq_along(gammas)) x[[paste0("chr_g", Name(gammas)[k])]] <- -2 * gammas[k] + 0.05 * gammas[k] * z2
  x$logPrior <- 0
  x$logJointW <- stats::rnorm(nSmp, sd = 0.1)
  x$check <- 0
  x$checkChar <- 0
  utils::write.csv(x, file.path(d, "cross.csv"), row.names = FALSE)
  invisible(d)
}
Launch <- function(root, cells, reps = 2) {
  utils::write.csv(data.frame(cell = cells, reps = reps), file.path(root, "launch.csv"),
                   row.names = FALSE, quote = FALSE)
}
rscript <- file.path(R.home("bin"), "Rscript")
Combine <- function(root, env = character()) {
  unlink(file.path(root, "bayes-factors.csv"))
  out <- suppressWarnings(system2(rscript,
                                  c(file.path(cetDir, "bayes-factors.R"), "combine", root),
                                  env = c("CET_MKV=pattern", "CET_BF_BURN=0.25", env),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(out = as.character(out), status = if (is.null(status)) 0L else status,
       edges = if (file.exists(file.path(root, "bayes-factors.csv"))) {
         utils::read.csv(file.path(root, "bayes-factors.csv"))
       })
}
Show <- function(res, pattern = "Error|Expected|Left out|no cross|not in|no complete|column|Edge") {
  paste(grep(pattern, res$out, value = TRUE), collapse = " / ")
}
cells4 <- c("b0g0", "b0g1", "b1g0", "b1g1")

# ---- 3. A ladder rung present (#100) ------------------------------------------------
root3 <- tempfile("rung")
for (lab in c(cells4, "b0.5g0")) for (rep in 1:2) MakeChain(root3, lab, rep, betas = c(0, 0.5, 1))
Launch(root3, c(cells4, "b0.5g0"))
res <- Combine(root3)
e <- res$edges
Say("rung b0.5g0 present, exit %d: %s; edges %s", res$status, Show(res, "Expected|Error"),
    if (is.null(e)) "none" else paste(sprintf("%s->%s %.3f (%d points)", e$from, e$to, e$logBF, e$points),
                                      collapse = ", "))
Check(res$status == 0L && !is.null(e) && nrow(e) == 4L,
      "combine with a rung present exits 0 and reports all 4 edges")
Check(!is.null(e) && any(e$from == "b0g0" & e$to == "b1g0" & e$points == 3L) &&
        isTRUE(abs(e$logBF[e$from == "b0g0" & e$to == "b1g0"] + 4) < 0.1),
      "b0g0 -> b1g0 runs through the rung (3 points), log BF about -4")
unlink(root3, recursive = TRUE)

root3b <- tempfile("rungOld")
for (lab in c(cells4, "b0.5g0")) for (rep in 1:2) {
  MakeChain(root3b, lab, rep, betas = c(0, 0.5, 1), oldNames = TRUE)
}
Launch(root3b, c(cells4, "b0.5g0"))
res <- Combine(root3b)
Say("rung present, old column names (aug_b0.0 ...), exit %d: %s", res$status, Show(res, "Error|column"))
Check(res$status != 0L && any(grepl("no column aug_b0", res$out)) && is.null(res$edges),
      "old-style columns: combine stops, names the missing column aug_b0, writes no bayes-factors.csv")
unlink(root3b, recursive = TRUE)

# ---- 4. Expected chain set (#101) --------------------------------------------------
root4 <- tempfile("set")
for (lab in cells4) for (rep in 1:2) MakeChain(root4, lab, rep)
Launch(root4, cells4)
full <- Combine(root4)
Say("4 cells x 2, launch.csv: exit %d: %s", full$status, Show(full, "Expected"))
Check(full$status == 0L && !is.null(full$edges) && nrow(full$edges) == 4L,
      "complete chain set: combine exits 0 with 4 edges")
Check(any(grepl("Expected chains \\(mkv pattern\\): 8 .*launch.csv.*found 8; missing 0", full$out)),
      "combine reports the expected set (8, from launch.csv) and that none is missing")

invisible(file.rename(file.path(root4, "b0g0_run1", "cross.csv"), file.path(root4, "b0g0_run1", "cross.bak")))
res <- Combine(root4)
Say("b0g0_run1/cross.csv removed: exit %d: %s", res$status, Show(res))
Check(res$status != 0L && any(grepl("no cross.csv.*b0g0_run1", res$out)) && is.null(res$edges),
      "missing cross.csv: combine stops, names b0g0_run1, writes no bayes-factors.csv (before: exit 0)")
res <- Combine(root4, "CET_BF_ALLOW_MISSING=b0g0_run1")
Say("... with CET_BF_ALLOW_MISSING=b0g0_run1: exit %d: %s", res$status, Show(res, "Left out|Error"))
Check(res$status == 0L && any(grepl("Left out on request \\(CET_BF_ALLOW_MISSING\\): b0g0_run1", res$out)) &&
        !is.null(res$edges) && nrow(res$edges) == 4L,
      "explicit opt-out: exit 0, says b0g0_run1 was left out")
res <- Combine(root4, "CET_BF_ALLOW_MISSING=b0g0_run9")
Check(res$status != 0L && any(grepl("not expected: b0g0_run9", res$out)),
      "opt-out naming an unexpected chain (b0g0_run9) stops")
invisible(file.rename(file.path(root4, "b0g0_run1", "cross.bak"), file.path(root4, "b0g0_run1", "cross.csv")))

unlink(file.path(root4, "b1g1_run2"), recursive = TRUE)
res <- Combine(root4)
Say("b1g1_run2 directory removed: exit %d: %s", res$status, Show(res))
Check(res$status != 0L && any(grepl("no complete segment.*b1g1_run2", res$out)),
      "an expected chain with no directory: combine stops and names it")
MakeChain(root4, "b1g1", 2)

MakeChain(root4, "b0.5g0", 1, betas = c(0, 0.5, 1))
res <- Combine(root4)
Say("extra chain b0.5g0_run1 not in launch.csv: exit %d: %s", res$status, Show(res))
Check(res$status != 0L && any(grepl("not in the launch record.*b0.5g0_run1", res$out)),
      "a chain not in the launch record: combine stops and names it")
unlink(file.path(root4, "b0.5g0_run1"), recursive = TRUE)

unlink(file.path(root4, "launch.csv"))
res <- Combine(root4)
Say("no launch.csv (default 4 x 4) with 4 x 2 chains: exit %d: %s", res$status, Show(res))
Check(res$status != 0L && any(grepl("production default", res$out)) &&
        any(grepl("no complete segment.*b0g0_run3", res$out)),
      "without launch.csv the production default (4 x 4) applies, and run3/run4 are missing")
res <- Combine(root4, "CET_BF_REPS=2")
Check(res$status == 0L, "CET_BF_REPS=2 sets the default to 4 x 2: exit 0")
unlink(root4, recursive = TRUE)

# ---- 5. Nonexistent root -------------------------------------------------------------
res <- Combine(file.path(tempdir(), "no-such-root"))
Check(res$status != 0L && any(grepl("does not exist", res$out)),
      "combine stops on a run root that does not exist")

if (!is.null(outFile)) writeLines(lines, outFile)
if (length(fail)) stop("check-labels.R FAILED: ", paste(fail, collapse = "; "))
cat("check-labels.R: all checks passed\n")
