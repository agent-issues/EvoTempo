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
# #103 review (synthetic chains are written as eval and run-chain.R write
# them: a complete progress.txt, and cross.csv with evalSamples, evalSegMd5 and
# evalTol; the Combine() helper no longer deletes outputs, so the checks see
# what combine leaves behind):
# 6. Power ladders (R103-4; verifier V3, V4): a complete pattern ladder of b1g1
#    (t = 0, 0.3, 0.6 and the power-1 chains) gives the stepping-stone computed
#    here; opting out b1g1's power-1 chains (CET_BF_ALLOW_MISSING) stops
#    combine with no marginal (before: -915 with 3 rungs, exit 0); a root with
#    only t = 0, 0.3, 0.6 stops ("lacks t = 1").
# 7. Whole-cell opt-out with a rung (V9; kills mutant M4, the Edge end-point
#    guard disabled): b0g0 left out with rung b0.5g0 present gives no edge from
#    b0g0, only b0g1 -> b1g1 and b1g0 -> b1g1, and says so.
# 8. Mkv filter of the stepping-stone (V13; kills mutant M5, the stepping-stone
#    reading Chains(root)): mkvnone power chains under the same root leave the
#    pattern marginal unchanged.
# 9. Stale outputs (R103-3; V1, V2): after a stop, no bayes-factors.csv,
#    model-probabilities.csv, closure.txt or chains-used.csv remain; a
#    whole-cell opt-out leaves a 2-edge bayes-factors.csv and no
#    model-probabilities.csv or closure.txt from the earlier 4-edge run.
# 10. Freshness (R103-2; V6): an unfinished chain stops combine unless named in
#    CET_BF_ALLOW_INCOMPLETE; a segment added after eval, a cross.csv that does
#    not record what eval read, and a check within CET_BF_TOL but beyond the
#    evalTol eval recorded (or the reverse) all stop combine.
# 11. Launch record (R103-7, -8, -9, -10; V5, V7, V11, V8): two directories for
#    one chain, reps = 2.7, launched chain directories with no segment that the
#    record does not list, and a record with no chain of the CET_MKV variant
#    all stop combine.
# 12. chains-used.csv (R103-12) lists every chain with its role, including the
#    chains left out.
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
                      cross = TRUE, charLL = -1500, progress = "2026-11-16 12:00:00 complete") {
  d <- file.path(root, sprintf("%s_run%d", label, rep))
  dir.create(d, recursive = TRUE)
  ll <- charLL + stats::rnorm(nSmp)
  utils::write.csv(data.frame(logPost = ll, charLL = ll, iteration = seq_len(nSmp)),
                   file.path(d, "seg0001_samples.csv"), row.names = FALSE)
  if (!is.null(progress)) writeLines(progress, file.path(d, "progress.txt"))
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
  # What eval read and its tolerance, as eval records them (R103-2)
  fp <- SegmentFingerprint(d)
  x$evalSamples <- fp$samples
  x$evalSegMd5 <- fp$md5
  x$evalTol <- 1e-6
  utils::write.csv(x, file.path(d, "cross.csv"), row.names = FALSE)
  invisible(d)
}
Launch <- function(root, cells, reps = 2, power = 1, mkv = "pattern") {
  utils::write.csv(data.frame(cell = cells, reps = reps, power = vapply(power, FormatValue, ""),
                              mkv = mkv), file.path(root, "launch.csv"), row.names = FALSE, quote = FALSE)
}
rscript <- file.path(R.home("bin"), "Rscript")
# combine itself must remove any earlier outputs (R103-3): this helper does not
Combine <- function(root, env = character()) {
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

# ---- 6. Power ladders (R103-4) ---------------------------------------------------
Marginal <- function(res) {
  hit <- grep("log p\\(characters \\| tip ages\\)", res$out, value = TRUE)
  if (length(hit) == 1L) as.numeric(sub(".*tip ages\\) = (-?[0-9.]+) .*", "\\1", hit)) else NA_real_
}
burn <- 0.25
Stone <- function(root, label, ladder) {
  logL <- lapply(ladder, function(t) unlist(lapply(1:2, function(k) {
    d <- file.path(root, sprintf("%s%s_run%d", label, if (t != 1) paste0("_t", FormatValue(t)) else "", k))
    x <- utils::read.csv(file.path(d, "seg0001_samples.csv"))$charLL
    x[seq_along(x) > floor(burn * length(x))]
  })))
  sum(vapply(seq_len(length(ladder) - 1L), function(k) LogMeanExp((ladder[k + 1] - ladder[k]) * logL[[k]]), 0))
}
Ladder <- function(root) {
  for (lab in cells4) for (k in 1:2) MakeChain(root, lab, k)
  for (t in c(0, 0.3, 0.6)) for (k in 1:2) {
    MakeChain(root, paste0("b1g1_t", FormatValue(t)), k, cross = FALSE, charLL = -1500 - 30 * (1 - t))
  }
  Launch(root, c(cells4, rep("b1g1", 3)), 2, c(1, 1, 1, 1, 0, 0.3, 0.6))
}
root6 <- tempfile("ladder")
Ladder(root6)
truth6 <- Stone(root6, "b1g1", c(0, 0.3, 0.6, 1))
res <- Combine(root6)
Say("V3a full pattern ladder of b1g1, exit %d: %s (expected %.3f)", res$status,
    Show(res, "log p|Error"), truth6)
Check(res$status == 0L && isTRUE(abs(Marginal(res) - truth6) < 1e-3) &&
        any(grepl("b1g1 \\(mkv pattern\\).*4 rungs, 8 chains", res$out)),
      "V3a: a complete ladder gives the stepping-stone computed here (4 rungs, 8 chains)")
for (k in 1:2) invisible(file.remove(file.path(root6, sprintf("b1g1_run%d", k), "cross.csv")))
res <- Combine(root6, "CET_BF_ALLOW_MISSING=b1g1_run1,b1g1_run2")
Say("V3b b1g1's power-1 chains opted out, exit %d: %s", res$status,
    Show(res, "log p|Error|lacks|left out"))
Check(res$status != 0L && is.na(Marginal(res)) &&
        any(grepl("cell b1g1 \\(mkv pattern\\): power ladder 0, 0.3, 0.6 lacks t = 1", res$out)) &&
        any(grepl("once b1g1_run1, b1g1_run2 \\(CET_BF_ALLOW_MISSING\\) are left out", res$out)),
      "V3b: opting out the t = 1 rung stops combine with no marginal (before: -915, 3 rungs, exit 0)")
Check(!file.exists(file.path(root6, "bayes-factors.csv")),
      "V3b: and leaves no bayes-factors.csv")
unlink(root6, recursive = TRUE)
root6b <- tempfile("ladderOnly")
for (t in c(0, 0.3, 0.6)) for (k in 1:2) {
  MakeChain(root6b, paste0("b1g1_t", FormatValue(t)), k, cross = FALSE, charLL = -1500 - 30 * (1 - t))
}
Launch(root6b, rep("b1g1", 3), 2, c(0, 0.3, 0.6))
res <- Combine(root6b)
Say("V4 ladder-only root (t = 0, 0.3, 0.6), exit %d: %s", res$status, Show(res, "log p|Error|lacks"))
Check(res$status != 0L && is.na(Marginal(res)) &&
        any(grepl("power ladder 0, 0.3, 0.6 lacks t = 1", res$out)),
      "V4: a ladder without t = 1 stops combine with no marginal (before: -915, exit 0)")
unlink(root6b, recursive = TRUE)

# ---- 7. Whole-cell opt-out with a rung (V9, mutant M4) ----------------------------
root7 <- tempfile("cellOut")
for (lab in c(cells4, "b0.5g0")) for (k in 1:2) MakeChain(root7, lab, k, betas = c(0, 0.5, 1))
Launch(root7, c(cells4, "b0.5g0"))
for (k in 1:2) invisible(file.remove(file.path(root7, sprintf("b0g0_run%d", k), "cross.csv")))
res <- Combine(root7, "CET_BF_ALLOW_MISSING=b0g0_run1,b0g0_run2")
e <- res$edges
Say("V9 b0g0 opted out, rung b0.5g0 present, exit %d: %s; edges %s", res$status, Show(res, "^Edge|Error"),
    if (is.null(e)) "none" else paste(sprintf("%s->%s %.4f (%d points)", e$from, e$to, e$logBF, e$points),
                                      collapse = ", "))
Check(res$status == 0L && !is.null(e) && !any(e$from == "b0g0") &&
        setequal(paste(e$from, e$to), c("b0g1 b1g1", "b1g0 b1g1")),
      "V9: no edge from the opted-out cell b0g0 (M4 reports b0g0->b1g0 from the b0.5g0 rung); 2 edges left")
Check(any(grepl("Edge b0g0 -> b1g0: no chains with cross.csv at an end point; edge left out", res$out)),
      "V9: combine says the b0g0 -> b1g0 edge is left out")
used <- utils::read.csv(file.path(root7, "chains-used.csv"))
Check(nrow(used) == 10L && sum(used$role == "left out (CET_BF_ALLOW_MISSING)") == 2L &&
        all(used$label[used$role != "bridge"] %in% c("b0g0_run1", "b0g0_run2")) &&
        all(!is.na(used$evalSegMd5[used$role == "bridge"])),
      "R103-12: chains-used.csv lists all 10 chains, the 2 left out, and eval's record of the 8 used")
unlink(root7, recursive = TRUE)

# ---- 8. Mkv filter of the stepping-stone (V13, mutant M5) ---------------------------
root8 <- tempfile("mkvFilter")
Ladder(root8)
res0 <- Combine(root8)
for (k in 1:2) MakeChain(root8, "b1g1_mkvnone_t0.3", k, cross = FALSE, charLL = -3000)
res1 <- Combine(root8)
Say("V13 pattern ladder %.4f; with mkvnone t = 0.3 chains (charLL -3000) under the root, exit %d: %.4f",
    Marginal(res0), res1$status, Marginal(res1))
Check(res1$status == 0L && isTRUE(abs(Marginal(res1) - Marginal(res0)) < 1e-9) &&
        any(grepl("4 rungs, 8 chains", res1$out)),
      "V13: mkvnone power chains do not join the pattern ladder (M5 and main pool them: -1520.70)")
unlink(root8, recursive = TRUE)

# ---- 9. Stale outputs (R103-3) -------------------------------------------------------
outs <- c("bayes-factors.csv", "model-probabilities.csv", "closure.txt", "chains-used.csv")
root9 <- tempfile("stale")
for (lab in cells4) for (k in 1:2) MakeChain(root9, lab, k)
Launch(root9, cells4)
res <- Combine(root9)
Check(res$status == 0L && all(file.exists(file.path(root9, outs))), "V1: a complete run writes all 4 outputs")
invisible(file.rename(file.path(root9, "b0g0_run1", "cross.csv"), file.path(root9, "b0g0_run1", "cross.bak")))
res <- Combine(root9)
Say("V1 after a stop (b0g0_run1/cross.csv moved away), exit %d; outputs left: %s", res$status,
    paste(outs[file.exists(file.path(root9, outs))], collapse = ", "))
Check(res$status != 0L && !any(file.exists(file.path(root9, outs))),
      "V1: a stop leaves none of the earlier outputs (before: all 3, byte-identical)")
invisible(file.rename(file.path(root9, "b0g0_run1", "cross.bak"), file.path(root9, "b0g0_run1", "cross.csv")))
res <- Combine(root9)
for (k in 1:2) invisible(file.remove(file.path(root9, sprintf("b1g1_run%d", k), "cross.csv")))
res <- Combine(root9, "CET_BF_ALLOW_MISSING=b1g1_run1,b1g1_run2")
Say("V2 whole-cell opt-out after a 4-edge run, exit %d; edges %d; outputs: %s", res$status,
    if (is.null(res$edges)) 0L else nrow(res$edges), paste(outs[file.exists(file.path(root9, outs))], collapse = ", "))
Check(res$status == 0L && nrow(res$edges) == 2L &&
        !any(file.exists(file.path(root9, c("model-probabilities.csv", "closure.txt")))),
      "V2: a 2-edge run leaves no model-probabilities.csv or closure.txt from the 4-edge run")
unlink(root9, recursive = TRUE)

# ---- 10. Freshness (R103-2) ----------------------------------------------------------
root10 <- tempfile("fresh")
for (lab in cells4) for (k in 1:2) MakeChain(root10, lab, k)
Launch(root10, cells4)
writeLines(c("2026-11-16 start", "2026-11-16 segment 1/20 done"), file.path(root10, "b0g1_run1", "progress.txt"))
res <- Combine(root10)
Say("V6a b0g1_run1 unfinished (segment 1/20), exit %d: %s", res$status, Show(res, "not complete|Error"))
Check(res$status != 0L && any(grepl("not complete.*b0g1_run1", res$out)),
      "V6a: an unfinished chain stops combine (before: exit 0)")
res <- Combine(root10, "CET_BF_ALLOW_INCOMPLETE=b0g1_run1")
used <- if (file.exists(file.path(root10, "chains-used.csv"))) utils::read.csv(file.path(root10, "chains-used.csv"))
Check(res$status == 0L && any(grepl("Used although not complete \\(CET_BF_ALLOW_INCOMPLETE\\): b0g1_run1", res$out)) &&
        isTRUE(used$incompleteAllowed[used$label == "b0g1_run1"]),
      "V6a: CET_BF_ALLOW_INCOMPLETE=b0g1_run1 uses it, says so and records it in chains-used.csv")
writeLines("2026-11-16 complete", file.path(root10, "b0g1_run1", "progress.txt"))
utils::write.csv(data.frame(logPost = -1500, charLL = -1500, iteration = seq_len(nSmp)),
                 file.path(root10, "b1g0_run2", "seg0002_samples.csv"), row.names = FALSE)
res <- Combine(root10)
Say("V6b segment added to b1g0_run2 after eval, exit %d: %s", res$status, Show(res, "stale|Error"))
Check(res$status != 0L && any(grepl("b1g0_run2: cross.csv is stale: eval read 40 samples .* now has 80", res$out)),
      "V6b: a segment added after eval stops combine (before: exit 0, logBF unchanged)")
invisible(file.remove(file.path(root10, "b1g0_run2", "seg0002_samples.csv")))
g <- file.path(root10, "b1g1_run1", "cross.csv")
x <- utils::read.csv(g, check.names = FALSE)
utils::write.csv(x[, setdiff(names(x), c("evalSamples", "evalSegMd5", "evalTol"))], g, row.names = FALSE)
res <- Combine(root10)
Check(res$status != 0L && any(grepl("b1g1_run1: cross.csv does not record what eval read", res$out)),
      "R103-2: a cross.csv without eval's record stops combine")
utils::write.csv(x, g, row.names = FALSE)
y <- x
y$check[3] <- 5e-7
y$evalTol <- 1e-7
utils::write.csv(y, g, row.names = FALSE)
res <- Combine(root10)
Check(res$status != 0L && any(grepl("b1g1_run1: cross.csv fails eval's checks .* > tolerance 1e-07", res$out)),
      "R103-2: |check| 5e-7 is within CET_BF_TOL 1e-6 but not the recorded evalTol 1e-7: stop")
y$evalTol <- 1e-6
utils::write.csv(y, g, row.names = FALSE)
res <- Combine(root10)
Check(res$status == 0L, "R103-2: the same check with evalTol 1e-6 passes")
res <- Combine(root10, "CET_BF_TOL=1e-7")
Check(res$status != 0L && any(grepl("tolerance 1e-07", res$out)),
      "R103-2: combine's own CET_BF_TOL 1e-7, stricter than evalTol 1e-6, stops")
utils::write.csv(x, g, row.names = FALSE)

# ---- 11. Launch record hardening (R103-7, -8, -9, -10) ---------------------------------
dir.create(file.path(root10, "b0g0_run01"))
invisible(file.copy(list.files(file.path(root10, "b0g0_run1"), full.names = TRUE), file.path(root10, "b0g0_run01")))
res <- Combine(root10)
Say("V5 b0g0_run01 beside b0g0_run1, exit %d: %s", res$status, Show(res, "same chain|Error"))
Check(res$status != 0L && any(grepl("two directories hold the same chain: .*b0g0_run01", res$out)),
      "V5: two directories for one chain stop combine (before: counted twice, exit 0)")
unlink(file.path(root10, "b0g0_run01"), recursive = TRUE)
Launch(root10, cells4, "2.7")
res <- Combine(root10)
Check(res$status != 0L && any(grepl("reps must be a whole number >= 1 \\(got '2.7'\\)", res$out)),
      "V7: reps = 2.7 stops combine (before: truncated to 2)")
Launch(root10, cells4)
dir.create(file.path(root10, "b0.5g0_run1"))
writeLines("2026-11-16 start", file.path(root10, "b0.5g0_run1", "progress.txt"))
dir.create(file.path(root10, "b1g1_run3"))
saveRDS(list(segment = 0), file.path(root10, "b1g1_run3", "checkpoint.rds"))
res <- Combine(root10)
Say("V11 launched dirs without a segment, not in the record, exit %d: %s", res$status, Show(res, "not in|Error"))
Check(res$status != 0L &&
        any(grepl("not in the launch record.*b0.5g0_run1 \\(no complete segment yet\\).*b1g1_run3 \\(no complete segment yet\\)",
                  res$out)),
      "V11: launched chain directories with no segment and not in the record stop combine (before: exit 0)")
unlink(file.path(root10, c("b0.5g0_run1", "b1g1_run3")), recursive = TRUE)
res <- Combine(root10)
Check(res$status == 0L, "the restored root combines again (exit 0)")
unlink(root10, recursive = TRUE)
root11 <- tempfile("empty")
for (lab in cells4) for (k in 1:2) MakeChain(root11, paste0(lab, "_mkvnone"), k)
Launch(root11, cells4, 2, 1, "none")
res <- Combine(root11)
Say("V8 only mkvnone chains and record, CET_MKV=pattern, exit %d: %s", res$status, Show(res, "Error|no expected"))
Check(res$status != 0L && any(grepl("no expected chain of Mkv variant pattern", res$out)),
      "V8: a record with no chain of the CET_MKV variant stops combine (before: exit 0, no Bayes factor)")
unlink(root11, recursive = TRUE)

# ---- 5. Nonexistent root -------------------------------------------------------------
res <- Combine(file.path(tempdir(), "no-such-root"))
Check(res$status != 0L && any(grepl("does not exist", res$out)),
      "combine stops on a run root that does not exist")

if (!is.null(outFile)) writeLines(lines, outFile)
if (length(fail)) stop("check-labels.R FAILED: ", paste(fail, collapse = "; "))
cat("check-labels.R: all checks passed\n")
