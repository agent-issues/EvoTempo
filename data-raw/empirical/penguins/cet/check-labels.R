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
for (mkv in names(expected)) {
  out <- suppressWarnings(system2(rscript,
                                  c(file.path(cetDir, "bayes-factors.R"), "combine", root2),
                                  env = c(paste0("CET_MKV=", mkv), "CET_BF_BURN=0.25"),
                                  stdout = TRUE, stderr = TRUE))
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

if (!is.null(outFile)) writeLines(lines, outFile)
if (length(fail)) stop("check-labels.R FAILED: ", paste(fail, collapse = "; "))
cat("check-labels.R: all checks passed\n")
