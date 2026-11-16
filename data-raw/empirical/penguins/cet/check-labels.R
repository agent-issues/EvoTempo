# Penguin CET pipeline: checks of the chain labels and of the Mkv-variant
# filter of the character-power stepping-stone in bayes-factors.R (#75).
#
#   Rscript data-raw/empirical/penguins/cet/check-labels.R [out.txt]
#
# Run from the package root with EvoTempo installed (or CET_LOAD_ALL=1 for
# pkgload::load_all()); needs no chain output.  Checks:
# 1. FormatValue() writes every power of the 33- and 17-point ladders
#    (k/K)^(1/0.3) and the powers named in #75 in plain decimal notation, and
#    ParseLabel() reads the power of the resulting label back exactly; the
#    same for ladder-rung cell labels such as CellLabel(1e-04, 1).
# 2. ParseLabel() stops on a label it cannot read (b1g1_t1e-04_run1), and
#    Chains() stops on a chain directory with such a label, while it ignores a
#    subdirectory without chain output (summary/).
# 3. `bayes-factors.R combine` on a synthetic run root that holds power
#    ladders of cell (1, 1) for both Mkv variants (4 powers x 2 replicates
#    each; charLL drawn so that the variants differ by 40 nats) reports, for
#    each CET_MKV, the stepping-stone estimate computed here from that
#    variant's chains only, and not the estimate from both variants pooled.
# 4. combine stops, naming the missing rung, and prints no marginal when a
#    ladder does not reach t = 1: a root holding only a power ladder
#    (0, 0.0001, 0.3; F1b of the #88 review), and an allTips ladder without
#    t = 1 beside a complete pattern ladder (F1a; pattern still works) (#88).
# 5. A bridge edge needs chains at both end points: with b0g0, b0g1, b1g0 and
#    a rung b0.5g1 but no b1g1 (F2a), and with rungs b0.25g1, b0.75g1 only
#    (F2b), combine leaves out the b0g1 -> b1g1 edge and says so, and still
#    reports the complete edges (#88).
# 6. combine stops on a run root that does not exist (#88).
# 7. Expected chain set (#101): with launch.csv listing 4 cells x 2 replicates,
#    combine exits 0; with one cross.csv removed it stops, names the chain and
#    writes no bayes-factors.csv (before: exit 0, chain silently dropped); with
#    CET_BF_ALLOW_MISSING naming it, exit 0 and the chain is reported as left
#    out; a chain not in the record, an expected chain with no directory, and
#    the production default (4 x 4) against 4 x 2 chains all stop.
# 8. eval's checks reach combine (#104; #105 review): with 4 cells x 2
#    replicates, combine exits 0; it stops, naming the chain and writing no
#    bayes-factors.csv, on a cross.csv with |check| 1e-3 in one state, with
#    |checkChar| 1e-3 (M4) or 1e-7 (> CET_BF_TOL_CHAR 1e-8, while |check| 1e-7
#    passes), with an NA check (M5), with no checkChar column (written before
#    #104), with no evalTol column, or with evalTol NA or "abc" (R105-7/8);
#    CET_BF_TOL = 0.01 does not loosen a file eval checked at 1e-4, and accepts
#    |check| 1e-3 in a file eval checked at 0.01.  eval's side (no cross.csv
#    unless the checks pass, the record it writes) needs a real chain:
#    check-eval.sh.
# 9. #103 review: a stop leaves none of the earlier outputs (R103-3); opting
#    out a cell's power-1 chains stops combine when that leaves the ladder
#    without t = 1 (R103-4); an unfinished chain stops it unless named in
#    CET_BF_ALLOW_INCOMPLETE, as do a segment added after eval and a cross.csv
#    that does not record what eval read (R103-2); two directories for one
#    chain (R103-7), reps = 2.7 (R103-8), a launched directory with no segment
#    missing from the record (R103-9) and a record with no chain of the
#    CET_MKV variant (R103-10) stop it; chains-used.csv lists the chains
#    (R103-12); a cross.csv without a rung's column names it (R103-13).
# Since #101 every synthetic root of checks 3-5 carries a launch.csv, and its
# power-1 chains a cross.csv, as combine now requires; since #104 a cross.csv
# has check and checkChar columns (and, since the #105 review, evalTol and
# evalTolChar); since the #103 review every synthetic chain
# has a complete progress.txt and its cross.csv records what eval read
# (evalSamples, evalSegMd5).  Check 9 also covers the stricter of CET_BF_TOL
# and the evalTol eval records (R103-2).
# Stops with an error if a check fails.
# The checks build short synthetic chains: a segment cap inherited from the
# environment (CET_MAX_SEG, #107 review) would make them incomplete
Sys.unsetenv(c("CET_MAX_SEG", "CET_MAX_SEG_RUNG"))
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/penguins/cet/common.R")
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
# Launch record of a synthetic root (#101): one row per group of chains
Launch <- function(root, cell, reps = 2, power = 1, mkv = "pattern") {
  utils::write.csv(data.frame(cell = cell, reps = reps, power = vapply(power, FormatValue, ""),
                              mkv = mkv), file.path(root, "launch.csv"), row.names = FALSE, quote = FALSE)
}

# ---- 1. Labels in plain decimal notation, read back exactly -------------------
powers <- list(`33-point ladder` = (0:32 / 32)^(1 / 0.3),
               `17-point ladder` = (0:16 / 16)^(1 / 0.3),
               `#75 examples and others` = c(1e-05, 5e-05, 1e-04, 2e-04, 5e-04, 1e-10, 1 / 3, 0.5))
for (nm in names(powers)) {
  p <- powers[[nm]]
  p <- p[p > 0 & p < 1]
  old <- vapply(p, format, "")
  new <- vapply(p, FormatValue, "")
  back <- vapply(new, function(s) {
    ParseLabel(sprintf("b1g1_mkvallTips_t%s_run1", s))$power
  }, 0)
  Say("%s: %d powers in (0, 1); format() wrote %d in scientific notation (%s)", nm, length(p),
      sum(grepl("e", old)), paste(old[grepl("e", old)], collapse = ", "))
  Check(!any(grepl("e", new)), sprintf("%s: FormatValue() never writes scientific notation", nm))
  Check(all(back == p), sprintf("%s: ParseLabel() reads every power back exactly", nm))
}
Check(identical(CellLabel(1e-04, 1), "b0.0001g1"), "CellLabel(1e-04, 1) is b0.0001g1")
lab <- ParseLabel(paste0(CellLabel(1e-04, 0.5), "_t0.0001_run3"))
Check(lab$beta == 1e-04 && lab$gamma == 0.5 && lab$power == 1e-04 && lab$rep == 3L &&
        lab$mkv == "pattern", "ParseLabel() reads a ladder-rung cell label")
Check(identical(vapply(c(0, 1, 0.5, 0.25), FormatValue, ""), c("0", "1", "0.5", "0.25")),
      "FormatValue() keeps the existing labels (b0g0 ... b1g1, 0.5)")
Check(identical(CellLabel(0, 1), "b0g1"), "production labels unchanged (CellLabel(0, 1) = b0g1)")

# ---- 2. Unreadable labels stop ---------------------------------------------------
Check(Errors(ParseLabel("b1g1_t1e-04_run1")), "ParseLabel() stops on b1g1_t1e-04_run1")
Check(Errors(ParseLabel("b1e-04g1_run1")), "ParseLabel() stops on b1e-04g1_run1")
root1 <- tempfile("labels")
dir.create(file.path(root1, "summary"), recursive = TRUE)
writeLines("x", file.path(root1, "summary", "clades.csv"))
dir.create(file.path(root1, "b1g1_run1"))
utils::write.csv(data.frame(logPost = 0, charLL = 0), file.path(root1, "b1g1_run1",
                                                                 "seg0001_samples.csv"))
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

# ---- 3. combine: the stepping-stone uses one Mkv variant ---------------------------
root2 <- tempfile("ladder")
ladder <- c(0, 1e-04, 0.3, 1)
nSmp <- 40L
burn <- 0.25
set.seed(7501)
truth <- list()
for (mkv in c("pattern", "allTips")) {
  for (t in ladder) {
    for (rep in 1:2) {
      label <- paste0(CellLabel(1, 1), MkvSuffix(mkv), if (t != 1) paste0("_t", FormatValue(t)))
      d <- file.path(root2, sprintf("%s_run%d", label, rep))
      dir.create(d, recursive = TRUE)
      charLL <- (if (mkv == "pattern") -1700 else -1660) - 30 * (1 - t) +
        stats::rnorm(nSmp, sd = 2 + 3 * (1 - t))
      utils::write.csv(data.frame(logPost = charLL, charLL = charLL, iteration = seq_len(nSmp)),
                       file.path(d, "seg0001_samples.csv"), row.names = FALSE)
      writeLines("2026-11-16 12:00:00 complete", file.path(d, "progress.txt"))
      # A power-1 chain must have a cross.csv (#101), with eval's checks (#104)
      # and its record of what it read (#103 review); no bridge uses this one
      fp <- SegmentFingerprint(d)
      if (t == 1) utils::write.csv(data.frame(sample = 1, check = 0, checkChar = 0,
                                              evalSamples = fp$samples, evalSegMd5 = fp$md5,
                                              evalTol = 1e-4, evalTolChar = 1e-8),
                                   file.path(d, "cross.csv"),
                                   row.names = FALSE)
      truth[[length(truth) + 1L]] <- data.frame(mkv = mkv, power = t, rep = rep,
                                                charLL = charLL[-seq_len(floor(burn * nSmp))])
    }
  }
}
truth <- do.call(rbind, truth)
grid2 <- expand.grid(power = ladder, mkv = c("pattern", "allTips"), stringsAsFactors = FALSE)
Launch(root2, "b1g1", 2, grid2$power, grid2$mkv)
Stone <- function(x) {
  sum(vapply(seq_len(length(ladder) - 1L), function(k) {
    LogMeanExp((ladder[k + 1] - ladder[k]) * x$charLL[x$power == ladder[k]])
  }, 0))
}
expected <- c(pattern = Stone(truth[truth$mkv == "pattern", ]),
              allTips = Stone(truth[truth$mkv == "allTips", ]))
pooled <- Stone(truth)
Say("synthetic ladder: expected pattern %.3f, allTips %.3f; both variants pooled %.3f",
    expected[["pattern"]], expected[["allTips"]], pooled)
rscript <- file.path(R.home("bin"), "Rscript")
Combine <- function(root, mkv = "pattern") {
  out <- suppressWarnings(system2(rscript,
                                  c(file.path(cetDir, "bayes-factors.R"), "combine", root),
                                  env = c(paste0("CET_MKV=", mkv), "CET_BF_BURN=0.25"),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(out = as.character(out), status = if (is.null(status)) 0L else status)
}
for (mkv in names(expected)) {
  out <- Combine(root2, mkv)$out
  hit <- grep("log p\\(characters \\| tip ages\\)", out, value = TRUE)
  Say("combine, CET_MKV=%s: %s", mkv, if (length(hit)) paste(hit, collapse = " / ") else
    paste(utils::tail(out, 5), collapse = " / "))
  got <- if (length(hit) == 1L) {
    as.numeric(sub(".*tip ages\\) = (-?[0-9.]+) .*", "\\1", hit))
  } else NA_real_
  Check(isTRUE(abs(got - expected[[mkv]]) < 1e-3),
        sprintf("combine (CET_MKV=%s) reports that variant's stepping-stone estimate", mkv))
  Check(isTRUE(abs(got - pooled) > 1),
        sprintf("combine (CET_MKV=%s) does not pool the two variants", mkv))
  Check(length(hit) == 1L && grepl(sprintf("mkv %s\\).*4 rungs, 8 chains", mkv), hit),
        sprintf("combine (CET_MKV=%s) names its variant and uses 4 rungs x 2 chains", mkv))
}
unlink(root2, recursive = TRUE)

# Synthetic chain directory: seg0001_samples.csv (charLL), and optionally a
# cross.csv as `eval` writes it, with log q(beta, gamma) = -4 beta - 2 gamma
# (plus a little draw-to-draw noise), so the bridge log BF is about -4 per
# unit beta and -2 per unit gamma
MakeChain <- function(root, label, rep, charLL = -1700, cross = FALSE,
                      betas = c(0, 1), gammas = c(0, 1), progress = "2026-11-16 12:00:00 complete") {
  d <- file.path(root, sprintf("%s_run%d", label, rep))
  dir.create(d, recursive = TRUE)
  ll <- charLL + stats::rnorm(nSmp)
  utils::write.csv(data.frame(logPost = ll, charLL = ll, iteration = seq_len(nSmp)),
                   file.path(d, "seg0001_samples.csv"), row.names = FALSE)
  if (!is.null(progress)) writeLines(progress, file.path(d, "progress.txt"))
  if (cross) {
    z <- stats::rnorm(nSmp)
    x <- data.frame(sample = seq_len(nSmp))
    for (b in betas) x[[paste0("aug_b", FormatValue(b))]] <- -4 * b + 0.05 * b * z
    z2 <- stats::rnorm(nSmp)
    for (g in gammas) x[[paste0("chr_g", FormatValue(g))]] <- -2 * g + 0.05 * g * z2
    x$logPrior <- 0
    x$check <- 0
    x$checkChar <- 0
    # What eval read, as eval records it (#103 review, R103-2)
    fp <- SegmentFingerprint(d)
    x$evalSamples <- fp$samples
    x$evalSegMd5 <- fp$md5
    x$evalTol <- 1e-4      # the tolerances eval applied (#105 review)
    x$evalTolChar <- 1e-8
    utils::write.csv(x, file.path(d, "cross.csv"), row.names = FALSE)
  }
}
HasMarginal <- function(out) any(grepl("log p\\(characters", out))

# ---- 4. combine stops on a ladder that does not reach t = 1 -------------------
root3 <- tempfile("noTop")
for (t in c(0, 1e-04, 0.3)) for (rep in 1:2) {
  MakeChain(root3, paste0("b1g1_t", FormatValue(t)), rep, charLL = -1700 - 30 * (1 - t))
}
res <- Combine(root3)
Say("F1b (power ladder 0, 0.0001, 0.3 only), exit %d: %s", res$status,
    paste(grep("ladder|lacks|Error", res$out, value = TRUE), collapse = " / "))
Check(res$status != 0L, "F1b: combine stops when the ladder lacks t = 1")
Check(!HasMarginal(res$out), "F1b: combine prints no marginal")
Check(any(grepl("cell b1g1 \\(mkv pattern\\): power ladder 0, 0.0001, 0.3 lacks t = 1", res$out)),
      "F1b: combine names the cell and the missing rung t = 1")
unlink(root3, recursive = TRUE)

root4 <- tempfile("mixedTop")
for (t in c(0, 1e-04, 0.3, 1)) for (rep in 1:2) {
  MakeChain(root4, paste0("b1g1", if (t != 1) paste0("_t", FormatValue(t))), rep,
            charLL = -1700 - 30 * (1 - t), cross = t == 1)
  if (t != 1) {
    MakeChain(root4, paste0("b1g1_mkvallTips_t", FormatValue(t)), rep,
              charLL = -1660 - 30 * (1 - t))
  }
}
Launch(root4, "b1g1", 2, c(0, 1e-04, 0.3, 1))
resA <- Combine(root4, "allTips")
resP <- Combine(root4, "pattern")
Say("F1a (allTips 0, 0.0001, 0.3 beside a complete pattern ladder): allTips exit %d (%s); pattern exit %d (%s)",
    resA$status, paste(grep("lacks", resA$out, value = TRUE), collapse = " / "), resP$status,
    paste(grep("log p\\(characters", resP$out, value = TRUE), collapse = " / "))
Check(resA$status != 0L && !HasMarginal(resA$out) &&
        any(grepl("mkv allTips\\): power ladder 0, 0.0001, 0.3 lacks t = 1", resA$out)),
      "F1a: CET_MKV=allTips stops, naming the missing t = 1, with no marginal")
Check(resP$status == 0L && sum(grepl("b1g1 \\(mkv pattern\\).*4 rungs, 8 chains", resP$out)) == 1L,
      "F1a: CET_MKV=pattern still reports its complete 4-rung ladder")
unlink(root4, recursive = TRUE)

# ---- 5. Edges need both end points ------------------------------------------------
EdgeRows <- function(root) {
  f <- file.path(root, "bayes-factors.csv")
  if (file.exists(f)) utils::read.csv(f) else NULL
}
root5 <- tempfile("F2a")
for (lab in c("b0g0", "b0g1", "b1g0", "b0.5g1")) for (rep in 1:2) {
  MakeChain(root5, lab, rep, cross = TRUE, betas = c(0, 0.5, 1))
}
Launch(root5, c("b0g0", "b0g1", "b1g0", "b0.5g1"))
res <- Combine(root5)
e <- EdgeRows(root5)
Say("F2a (b0g0, b0g1, b1g0, rung b0.5g1, no b1g1), exit %d: %s; edges %s", res$status,
    paste(grep("^Edge", res$out, value = TRUE), collapse = " / "),
    if (is.null(e)) "none" else paste(sprintf("%s->%s %.3f (%d points)", e$from, e$to, e$logBF,
                                              e$points), collapse = ", "))
Check(res$status == 0L && !is.null(e) && !any(e$from == "b0g1" & e$to == "b1g1"),
      "F2a: no b0g1 -> b1g1 edge without b1g1 chains (before: the b0g1 -> b0.5g1 bridge, -2)")
Check(any(grepl("Edge b0g1 -> b1g1: no chains with cross.csv at b1g1; edge left out", res$out)),
      "F2a: combine says which end point of b0g1 -> b1g1 is missing")
Check(!is.null(e) && nrow(e) == 2L && all(abs(e$logBF[e$from == "b0g0" & e$to == "b1g0"] + 4) < 0.1) &&
        all(abs(e$logBF[e$from == "b0g0" & e$to == "b0g1"] + 2) < 0.1),
      "F2a: the two complete edges are reported (about -4 and -2)")
unlink(root5, recursive = TRUE)

root6 <- tempfile("F2b")
for (lab in c("b0g0", "b1g0", "b0.25g1", "b0.75g1")) for (rep in 1:2) {
  MakeChain(root6, lab, rep, cross = TRUE, betas = c(0, 0.25, 0.75, 1))
}
Launch(root6, c("b0g0", "b1g0", "b0.25g1", "b0.75g1"))
res <- Combine(root6)
e <- EdgeRows(root6)
Say("F2b (rungs b0.25g1, b0.75g1 only), exit %d: %s", res$status,
    paste(grep("^Edge", res$out, value = TRUE), collapse = " / "))
Check(res$status == 0L && !is.null(e) && !any(e$from == "b0g1" & e$to == "b1g1"),
      "F2b: no b0g1 -> b1g1 edge from the interior rungs alone (before: -2.0)")
Check(any(grepl("Edge b0g1 -> b1g1: no chains with cross.csv at b0g1 and b1g1", res$out)),
      "F2b: combine names both missing end points")
unlink(root6, recursive = TRUE)

# ---- 6. Nonexistent root ------------------------------------------------------------
res <- Combine(file.path(tempdir(), "no-such-root"))
Say("nonexistent root, exit %d: %s", res$status,
    paste(grep("exist", res$out, value = TRUE), collapse = " / "))
Check(res$status != 0L && any(grepl("does not exist", res$out)),
      "combine stops on a run root that does not exist")

# ---- 7. Expected chain set (#101) ------------------------------------------------
# combine itself must remove any earlier outputs (#103 review, R103-3): this
# helper does not
CombineEnv <- function(root, env = character()) {
  out <- suppressWarnings(system2(rscript, c(file.path(cetDir, "bayes-factors.R"), "combine", root),
                                  env = c("CET_MKV=pattern", "CET_BF_BURN=0.25", env),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(out = as.character(out), status = if (is.null(status)) 0L else status,
       edges = EdgeRows(root))
}
Problems <- function(res) {
  paste(grep("Expected|Left out|no cross|not in|no complete|not expected", res$out, value = TRUE),
        collapse = " / ")
}
root7 <- tempfile("set")
cells4 <- c("b0g0", "b0g1", "b1g0", "b1g1")
for (lab in cells4) for (rep in 1:2) MakeChain(root7, lab, rep, cross = TRUE)
Launch(root7, cells4)
res <- CombineEnv(root7)
Say("#101 4 cells x 2 with launch.csv, exit %d: %s", res$status, Problems(res))
Check(res$status == 0L && !is.null(res$edges) && nrow(res$edges) == 4L,
      "#101: complete chain set, combine exits 0 with 4 edges")
invisible(file.rename(file.path(root7, "b1g0_run2", "cross.csv"), file.path(root7, "b1g0_run2", "cross.bak")))
res <- CombineEnv(root7)
Say("#101 b1g0_run2/cross.csv removed, exit %d: %s", res$status, Problems(res))
Check(res$status != 0L && any(grepl("no cross.csv.*b1g0_run2", res$out)) && is.null(res$edges),
      "#101: a missing cross.csv stops combine, naming b1g0_run2 (before: exit 0, dropped)")
res <- CombineEnv(root7, "CET_BF_ALLOW_MISSING=b1g0_run2")
Say("#101 ... with CET_BF_ALLOW_MISSING=b1g0_run2, exit %d: %s", res$status, Problems(res))
Check(res$status == 0L && any(grepl("Left out on request.*b1g0_run2", res$out)) &&
        !is.null(res$edges) && nrow(res$edges) == 4L,
      "#101: explicit opt-out, exit 0 and b1g0_run2 reported as left out")
res <- CombineEnv(root7, "CET_BF_ALLOW_MISSING=b1g0_run7")
Check(res$status != 0L && any(grepl("not expected: b1g0_run7", res$out)),
      "#101: an opt-out naming an unexpected chain stops")
invisible(file.rename(file.path(root7, "b1g0_run2", "cross.bak"), file.path(root7, "b1g0_run2", "cross.csv")))
MakeChain(root7, "b1g0.5", 1, cross = TRUE, gammas = c(0, 0.5, 1))
res <- CombineEnv(root7)
Check(res$status != 0L && any(grepl("not in the launch record.*b1g0.5_run1", res$out)),
      "#101: a chain not in the launch record stops combine")
unlink(file.path(root7, "b1g0.5_run1"), recursive = TRUE)
unlink(file.path(root7, "b0g1_run1"), recursive = TRUE)
res <- CombineEnv(root7)
Check(res$status != 0L && any(grepl("no complete segment.*b0g1_run1", res$out)),
      "#101: an expected chain with no directory stops combine")
MakeChain(root7, "b0g1", 1, cross = TRUE)
unlink(file.path(root7, "launch.csv"))
res <- CombineEnv(root7)
Say("#101 no launch.csv, 4 x 2 chains, exit %d: %s", res$status, Problems(res))
Check(res$status != 0L && any(grepl("production default", res$out)) &&
        any(grepl("no complete segment.*b0g0_run3", res$out)),
      "#101: without launch.csv the production default (4 x 4) applies; run3/run4 missing")
unlink(root7, recursive = TRUE)

# ---- 8. eval's checks reach combine (#104) ---------------------------------------
root8 <- tempfile("checks")
for (lab in cells4) for (rep in 1:2) MakeChain(root8, lab, rep, cross = TRUE)
Launch(root8, cells4)
res <- CombineEnv(root8)
Check(res$status == 0L && !is.null(res$edges) && nrow(res$edges) == 4L,
      "#104: cross.csv files that pass their checks, combine exits 0 with 4 edges")
f <- file.path(root8, "b0g1_run1", "cross.csv")
good <- utils::read.csv(f, check.names = FALSE)
Refused <- function(mod, msg, pattern, env = character()) {
  utils::write.csv(mod, f, row.names = FALSE)
  res <- CombineEnv(root8, env)
  Say("#104/#105 %s, exit %d: %s", msg, res$status,
      paste(grep("b0g1_run1:", res$out, value = TRUE), collapse = " / "))
  Check(res$status != 0L && any(grepl(paste0("b0g1_run1: ", pattern), res$out)) && is.null(res$edges),
        sprintf("#104/#105: %s stops combine, naming b0g1_run1", msg))
  utils::write.csv(good, f, row.names = FALSE)
}
bad <- good; bad$check[7] <- 1e-3
Refused(bad, "one state with |check| = 1e-3", "cross.csv fails eval's checks \\(max \\|check\\| 0.001")
Refused(bad, "|check| = 1e-3 under CET_BF_TOL = 0.01 but evalTol 1e-4 (no loosening)",
        "cross.csv fails eval's checks", "CET_BF_TOL=0.01")
bad <- good; bad$checkChar[5] <- 1e-3
Refused(bad, "one state with |checkChar| = 1e-3 (M4)", "cross.csv fails eval's checks.*max \\|checkChar\\| 0.001")
bad <- good; bad$checkChar[5] <- 1e-7
Refused(bad, "|checkChar| = 1e-7 > CET_BF_TOL_CHAR 1e-8", "cross.csv fails eval's checks.*checkChar\\| 1e-07")
bad <- good; bad$check[3] <- NA
Refused(bad, "an NA check (M5)", "cross.csv fails eval's checks: 1 non-finite")
bad <- good; bad$evalTol <- NULL
Refused(bad, "no evalTol column", "cross.csv records no evalTol")
bad <- good; bad$evalTol <- NA
Refused(bad, "evalTol NA", "cross.csv records evalTol NA")
bad <- good; bad$evalTol <- "abc"
Refused(bad, "evalTol \"abc\"", "cross.csv records evalTol abc")
bad <- good; bad$check[7] <- 1e-7
utils::write.csv(bad, f, row.names = FALSE)
res <- CombineEnv(root8)
Check(res$status == 0L && !is.null(res$edges) && nrow(res$edges) == 4L,
      "#105: |check| 1e-7 is within CET_BF_TOL (1e-4): accepted")
bad <- good; bad$check[7] <- 1e-3; bad$evalTol <- 0.01
utils::write.csv(bad, f, row.names = FALSE)
res <- CombineEnv(root8, "CET_BF_TOL=0.01")
Check(res$status == 0L && !is.null(res$edges) && nrow(res$edges) == 4L,
      "#104: CET_BF_TOL = 0.01 accepts |check| = 1e-3 in a file eval checked at 0.01")
utils::write.csv(good, f, row.names = FALSE)
g <- file.path(root8, "b1g1_run2", "cross.csv")
old <- utils::read.csv(g, check.names = FALSE)
old$checkChar <- NULL
utils::write.csv(old, g, row.names = FALSE)
res <- CombineEnv(root8)
Say("#104 b1g1_run2 cross.csv without checkChar, exit %d: %s", res$status,
    paste(grep("b1g1_run2:", res$out, value = TRUE), collapse = " / "))
Check(res$status != 0L && any(grepl("b1g1_run2: cross.csv has no check or checkChar column", res$out)) && is.null(res$edges),
      "#104: a cross.csv from before #104 (no checkChar) stops combine, naming b1g1_run2")
unlink(root8, recursive = TRUE)

# ---- 9. #103 review ---------------------------------------------------------------
outs <- c("bayes-factors.csv", "model-probabilities.csv", "closure.txt", "chains-used.csv")
root9 <- tempfile("review103")
for (lab in cells4) for (rep in 1:2) MakeChain(root9, lab, rep, cross = TRUE)
Launch(root9, cells4)
res <- CombineEnv(root9)
Check(res$status == 0L && all(file.exists(file.path(root9, outs))),
      "#103: a complete run writes bayes-factors.csv, model-probabilities.csv, closure.txt, chains-used.csv")
invisible(file.rename(file.path(root9, "b0g0_run1", "cross.csv"), file.path(root9, "b0g0_run1", "cross.bak")))
res <- CombineEnv(root9)
Check(res$status != 0L && !any(file.exists(file.path(root9, outs))),
      "R103-3: a stop leaves none of the earlier outputs (before: all three, unchanged)")
invisible(file.rename(file.path(root9, "b0g0_run1", "cross.bak"), file.path(root9, "b0g0_run1", "cross.csv")))
writeLines(c("2026-11-16 start", "2026-11-16 segment 1/20 done"), file.path(root9, "b0g1_run2", "progress.txt"))
res <- CombineEnv(root9)
Say("R103-2 b0g1_run2 unfinished, exit %d: %s", res$status, Problems(res))
Check(res$status != 0L && any(grepl("not complete.*b0g1_run2", res$out)),
      "R103-2: an unfinished chain stops combine")
res <- CombineEnv(root9, "CET_BF_ALLOW_INCOMPLETE=b0g1_run2")
Check(res$status == 0L && any(grepl("Used although not complete.*b0g1_run2", res$out)),
      "R103-2: CET_BF_ALLOW_INCOMPLETE=b0g1_run2 uses it and says so")
writeLines("2026-11-16 complete", file.path(root9, "b0g1_run2", "progress.txt"))
utils::write.csv(data.frame(logPost = -1700, charLL = -1700, iteration = seq_len(nSmp)),
                 file.path(root9, "b1g0_run1", "seg0002_samples.csv"), row.names = FALSE)
res <- CombineEnv(root9)
Check(res$status != 0L && any(grepl("b1g0_run1: cross.csv is stale", res$out)),
      "R103-2: a segment added after eval stops combine")
invisible(file.remove(file.path(root9, "b1g0_run1", "seg0002_samples.csv")))
g <- file.path(root9, "b1g1_run1", "cross.csv")
x <- utils::read.csv(g, check.names = FALSE)
utils::write.csv(x[, setdiff(names(x), c("evalSamples", "evalSegMd5"))], g, row.names = FALSE)
res <- CombineEnv(root9)
Check(res$status != 0L && any(grepl("b1g1_run1: cross.csv does not record what eval read", res$out)),
      "R103-2: a cross.csv without eval's record stops combine")
utils::write.csv(x, g, row.names = FALSE)
dir.create(file.path(root9, "b0g0_run01"))
invisible(file.copy(list.files(file.path(root9, "b0g0_run1"), full.names = TRUE), file.path(root9, "b0g0_run01")))
res <- CombineEnv(root9)
Check(res$status != 0L && any(grepl("two directories hold the same chain: .*b0g0_run01", res$out)),
      "R103-7: two directories for one chain stop combine")
unlink(file.path(root9, "b0g0_run01"), recursive = TRUE)
Launch(root9, cells4, "2.7")
res <- CombineEnv(root9)
Check(res$status != 0L && any(grepl("reps must be a whole number", res$out)), "R103-8: reps = 2.7 stops combine")
Launch(root9, cells4)
dir.create(file.path(root9, "b0g0.5_run1"))
writeLines("2026-11-16 start", file.path(root9, "b0g0.5_run1", "progress.txt"))
res <- CombineEnv(root9)
Check(res$status != 0L && any(grepl("not in the launch record.*b0g0.5_run1 \\(no complete segment yet\\)", res$out)),
      "R103-9: a launched directory with no segment, missing from the record, stops combine")
unlink(file.path(root9, "b0g0.5_run1"), recursive = TRUE)
res <- CombineEnv(root9)
used <- if (file.exists(file.path(root9, "chains-used.csv"))) utils::read.csv(file.path(root9, "chains-used.csv"))
Check(res$status == 0L && !is.null(used) && nrow(used) == 8L && all(used$role == "bridge") &&
        all(grepl("^md5:", used$evalSegMd5)),
      "R103-12: chains-used.csv lists the 8 chains, their role and eval's record")
unlink(root9, recursive = TRUE)
root9b <- tempfile("ladderOut")
for (lab in cells4) for (rep in 1:2) MakeChain(root9b, lab, rep, cross = TRUE)
for (t in c(0, 0.3)) for (rep in 1:2) MakeChain(root9b, paste0("b1g1_t", FormatValue(t)), rep)
Launch(root9b, c(cells4, "b1g1", "b1g1"), 2, c(1, 1, 1, 1, 0, 0.3))
for (rep in 1:2) invisible(file.remove(file.path(root9b, sprintf("b1g1_run%d", rep), "cross.csv")))
res <- CombineEnv(root9b, "CET_BF_ALLOW_MISSING=b1g1_run1,b1g1_run2")
Say("R103-4 b1g1's power-1 chains opted out, exit %d: %s", res$status,
    paste(grep("lacks|log p", res$out, value = TRUE), collapse = " / "))
Check(res$status != 0L && !HasMarginal(res$out) &&
        any(grepl("once b1g1_run1, b1g1_run2 \\(CET_BF_ALLOW_MISSING\\) are left out", res$out)) &&
        any(grepl("power ladder 0, 0.3 lacks t = 1", res$out)),
      "R103-4: opting out the t = 1 rung stops combine with no marginal (before: stopifnot)")
unlink(root9b, recursive = TRUE)
root9c <- tempfile("rungColumn")
for (lab in cells4) for (rep in 1:2) MakeChain(root9c, lab, rep, cross = TRUE)
for (rep in 1:2) MakeChain(root9c, "b0g0.5", rep, cross = TRUE, gammas = c(0, 0.5, 1))
Launch(root9c, c(cells4, "b0g0.5"))
res <- CombineEnv(root9c)
Check(res$status != 0L && any(grepl("cross.csv has no column chr_g0.5; rerun eval", res$out)),
      "R103-13: a cross.csv without the rung's column names it (before: 'missing value where TRUE/FALSE needed')")
unlink(root9c, recursive = TRUE)
root9d <- tempfile("empty")
for (lab in cells4) for (rep in 1:2) MakeChain(root9d, paste0(lab, "_mkvnone"), rep, cross = TRUE)
Launch(root9d, cells4, 2, 1, "none")
res <- CombineEnv(root9d)
Check(res$status != 0L && any(grepl("no expected chain of Mkv variant pattern", res$out)),
      "R103-10: a record with no chain of the CET_MKV variant stops combine")
unlink(root9d, recursive = TRUE)
# evalTol: the stricter of CET_BF_TOL (default 1e-4) and the tolerance eval
# recorded applies (R103-2)
root9e <- tempfile("evalTol")
for (lab in cells4) for (rep in 1:2) MakeChain(root9e, lab, rep, cross = TRUE)
Launch(root9e, cells4)
g <- file.path(root9e, "b1g0_run2", "cross.csv")
x <- utils::read.csv(g, check.names = FALSE)
x$check[2] <- 5e-5
x$evalTol <- 1e-5
utils::write.csv(x, g, row.names = FALSE)
res <- CombineEnv(root9e)
Check(res$status != 0L && any(grepl("b1g0_run2: cross.csv fails eval's checks", res$out)),
      "R103-2: |check| 5e-5 within CET_BF_TOL 1e-4 but beyond the recorded evalTol 1e-5 stops combine")
x$evalTol <- 1e-4
utils::write.csv(x, g, row.names = FALSE)
res <- CombineEnv(root9e)
Check(res$status == 0L, "R103-2: the same check with evalTol 1e-4 passes")
unlink(root9e, recursive = TRUE)

if (!is.null(outFile)) writeLines(lines, outFile)
if (length(fail)) stop("check-labels.R FAILED: ", paste(fail, collapse = "; "))
cat("check-labels.R: all checks passed\n")
