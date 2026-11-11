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
# Stops with an error if a check fails.
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
      truth[[length(truth) + 1L]] <- data.frame(mkv = mkv, power = t, rep = rep,
                                                charLL = charLL[-seq_len(floor(burn * nSmp))])
    }
  }
}
truth <- do.call(rbind, truth)
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
                      betas = c(0, 1), gammas = c(0, 1)) {
  d <- file.path(root, sprintf("%s_run%d", label, rep))
  dir.create(d, recursive = TRUE)
  ll <- charLL + stats::rnorm(nSmp)
  utils::write.csv(data.frame(logPost = ll, charLL = ll, iteration = seq_len(nSmp)),
                   file.path(d, "seg0001_samples.csv"), row.names = FALSE)
  if (cross) {
    z <- stats::rnorm(nSmp)
    x <- data.frame(sample = seq_len(nSmp))
    for (b in betas) x[[paste0("aug_b", FormatValue(b))]] <- -4 * b + 0.05 * b * z
    z2 <- stats::rnorm(nSmp)
    for (g in gammas) x[[paste0("chr_g", FormatValue(g))]] <- -2 * g + 0.05 * g * z2
    x$logPrior <- 0
    x$check <- 0
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
            charLL = -1700 - 30 * (1 - t))
  if (t != 1) {
    MakeChain(root4, paste0("b1g1_mkvallTips_t", FormatValue(t)), rep,
              charLL = -1660 - 30 * (1 - t))
  }
}
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

if (!is.null(outFile)) writeLines(lines, outFile)
if (length(fail)) stop("check-labels.R FAILED: ", paste(fail, collapse = "; "))
cat("check-labels.R: all checks passed\n")
