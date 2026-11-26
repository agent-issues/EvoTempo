# Trilobite CET pipeline: fixed-tree checks of the character likelihood and of
# the tree density's conditioning against RevBayes, on the start tree
# (data/start-tree-01.nwk).
#
#   Rscript data-raw/empirical/trilobites/cet/check-mkv.R write <oracle dir>
#   (run RevBayes: revbayes/check-oracle.sh does both steps and this one)
#   Rscript data-raw/empirical/trilobites/cet/check-mkv.R compare <oracle dir> [out dir]
#
# `write` writes the oracle inputs for revbayes/check-oracle.Rev: the start
# tree with branch lengths rate x duration at strict-clock rates 0.005, 0.0236
# (the comparator's posterior median clock) and 0.05, the same tree in Myr,
# the Gamma shapes, and a grid of (lambda, mu, psi, origin).
#
# `compare` checks, and stops with an error if a check fails (tolerances fixed
# before the run):
# 1. Characters, Gamma off.  The package's MkLikelihood() on MkCharacters()
#    (LoadCharacters(): coding = "variable", per-character K and ordering
#    from partitions.csv) with nRateCat = 1, per RevBayes partition (u2-u16,
#    o4, o5) and in total, must equal RevBayes dnPhyloCTMC(coding =
#    "variable") with the comparator's Q and no site-rate variation, read from
#    the comparator's own matrices (../revbayes/prepare-data.R, which treats
#    polymorphisms as ambiguous): |difference| <= 1e-6 per partition, both
#    engines.  The totals must also equal those of the pipeline's own
#    character code before the switch to MkCharacters() (PartitionedMk(), PR
#    #91, smoke-results/mkv-check.csv there) to 1e-5 (they were reported to 6
#    decimals).
# 2. Characters, Gamma on.  The same with nRateCat = 4 and the shape alpha
#    at 0.25, 0.6 and 2, against dnPhyloCTMC(siteRates =
#    fnDiscretizeGamma(alpha, alpha, 4), coding = "variable"), as the
#    comparator: |difference| <= 1e-6 per partition, both engines.
# 3. The whole matrix as one MkCharacters() object (the pipeline's) equals
#    the sum of its partitions to 1e-8; the engines agree to 1e-8 relative;
#    PoweredMk() at power 1 returns the package's value exactly and at power
#    t gives t x lnL.
# 4. Conditioning and time origin.  With beta = 0 the dCET tree density equals
#    the FBD density for any generator (spec section 7.2), so TreeLikelihood()
#    under CetModel(0, 0) (rho = 0, condition = "survival", start = "origin",
#    ages before 485.4 Ma) must differ from RevBayes dnFBDP(rho = 0,
#    condition = "sampling") by a constant over the grid of parameters and
#    origins: spread (max - min of the difference) <= 1e-6 (relative to
#    max |lnL| 1e-8 when larger), with both engines.  The offset is known
#    exactly (review of #91, F2): RevBayes - TreeLikelihood() =
#    -ln(n!) + nBif ln 2, with n = 107 tips and nBif the bifurcations (internal
#    nodes without a sampled ancestor, the root included under an origin), so
#    it is a constant only for a given tree (nBif changes with the number of
#    sampled ancestors); each point's residual from it must be <= the same
#    tolerance.  TempoMCMC() adds nBif ln 2 in every tree-move ratio
#    (.TreeMove()), so the sampled tree density is RevBayes's + ln(n!), the
#    same constant for every tree; the samples' logPost does not include
#    nBif ln 2.  (FBDTreeLikelihood(),
#    the closed form, gives NaN here: with rho = 0 and no extant tip it adds
#    0 x log(0); TreeLikelihood() is unaffected.)  Negative controls, which
#    must FAIL that
#    test by a wide margin (spread > 1e-2): condition = "none"; and the
#    present at the youngest tip (486.0 Ma) instead of 485.4 Ma.
# 5. Cost per MkLikelihood() evaluation (ms), by engine, rate categories and
#    coding.
# Writes <out dir>/mkv-check.csv, fbd-check.csv and mkv-check-timing.csv
# (default <cetDir>/smoke-results/).
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1]]
oracleDir <- args[[2]]
outDir <- if (length(args) > 2) args[[3]] else file.path(cetDir, "smoke-results")
rates <- c(0.005, 0.0236, 0.05)
alphas <- c(0.25, 0.6, 2)

dat <- LoadTrilobites(StartTreeFile(1))
tt <- dat$tree
dur <- numeric(tt$nTip + tt$nNode)
child <- which(!is.na(tt$parent))
dur[child] <- tt$age[tt$parent[child]] - tt$age[child]
rootAge <- tt$age[[tt$root]]
grid <- data.frame(lambda = c(3.9, 0.5, 1.0, 0.2, 2.0, 0.1),
                   mu = c(3.87, 0.3, 0.9, 0.05, 1.98, 0.09),
                   psi = c(0.01, 0.05, 0.02, 0.1, 0.005, 0.3),
                   origin = rootAge + c(2.5, 0.5, 5, 20, 50, 100))
AsPhylo <- function(tree, scale = 1) {
  s <- tree
  s$age <- s$age * scale
  ph <- read.tree(text = TempoNewick(s, digits = 12))
  ph$root.edge <- NULL
  ph
}

if (mode == "write") {
  dir.create(oracleDir, recursive = TRUE, showWarnings = FALSE)
  for (r in seq_along(rates)) {
    write.tree(AsPhylo(tt, rates[r]), file.path(oracleDir, sprintf("mkv-tree-%d.nwk", r)),
               digits = 12)
  }
  write.tree(AsPhylo(tt), file.path(oracleDir, "fbd-tree.nwk"), digits = 12)
  V <- function(name, x) sprintf("%s = v(%s)", name, paste(sprintf("%.12g", x), collapse = ", "))
  writeLines(c(V("lam_v", grid$lambda), V("mu_v", grid$mu), V("psi_v", grid$psi),
               V("ori_v", grid$origin)), file.path(oracleDir, "fbd-params.Rev"))
  writeLines(V("alpha_v", alphas), file.path(oracleDir, "gamma-params.Rev"))
  write.csv(grid, file.path(oracleDir, "fbd-grid.csv"), row.names = FALSE)
  cat("oracle inputs written to", oracleDir, "\n")
  quit(save = "no")
}

stopifnot(mode == "compare")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
# RevBayes write() ends no record with a newline, so the records run together:
# split them by pattern (kind, case, partition, lnL with 14 significant digits)
ReadOracle <- function(file) {
  txt <- paste(readLines(file, warn = FALSE), collapse = "")
  rec <- regmatches(txt, gregexpr("(mkvg|mkv|fbd)\t[0-9]+\t[A-Za-z0-9]+\t-?[0-9.]+(e[-+]?[0-9]+)?", txt))[[1]]
  f <- do.call(rbind, strsplit(rec, "\t"))
  data.frame(kind = f[, 1], case = as.integer(f[, 2]), partition = f[, 3], lnL = as.numeric(f[, 4]))
}
rb <- ReadOracle(file.path(oracleDir, "revbayes-oracle.tsv"))
cells <- dat$cells
parts <- dat$parts
nPart <- length(unique(parts$part))
stopifnot(sum(rb$kind == "mkv") == 3 * nPart, sum(rb$kind == "mkvg") == 9 * nPart,
          sum(rb$kind == "fbd") == nrow(grid))
fail <- character()
# One partition (or the whole matrix) as MkCharacters(), rows in tree order
Mk <- function(cols, nRateCat, coding = "variable") {
  MkCharacters(cells[, cols, drop = FALSE], nState = parts$k[cols], ordered = parts$ordered[cols],
               coding = coding, nRateCat = nRateCat)
}
partMk <- lapply(setNames(nm = unique(parts$part)), function(part) {
  cols <- which(parts$part == part)
  list(off = Mk(cols, 1L), on = Mk(cols, 4L))
})
whole <- list(off = LoadCharacters(tt$tipLabel, nRateCat = 1L)$characters,
              on = dat$characters)
stopifnot(whole$on$nRateCat == 4L, whole$on$coding == "variable",
          identical(whole$on$tipLabel, tt$tipLabel))
# Totals of PartitionedMk() on this tree (PR #91, smoke-results/mkv-check.csv)
before <- c(-2805.185423, -2667.526174, -3123.850355)
rows <- list()
for (r in seq_along(rates)) {
  len <- rates[r] * dur
  for (a in c(NA, seq_along(alphas))) {
    gam <- !is.na(a)
    alpha <- if (gam) alphas[a] else NULL
    kind <- if (gam) "mkvg" else "mkv"
    case <- if (gam) 3L * (r - 1L) + a else r
    for (eng in c("cpp", "R")) {
      for (part in names(partMk)) {
        ref <- rb$lnL[rb$kind == kind & rb$case == case & rb$partition == part]
        stopifnot(length(ref) == 1L)
        rows[[length(rows) + 1L]] <- data.frame(
          rate = rates[r], alpha = if (gam) alphas[a] else NA, engine = eng, partition = part,
          nChar = sum(parts$part == part),
          cet = MkLikelihood(tt, partMk[[part]][[if (gam) "on" else "off"]], len, engine = eng,
                             alpha = alpha),
          revbayes = ref)
      }
      w <- whole[[if (gam) "on" else "off"]]
      total <- MkLikelihood(tt, w, len, engine = eng, alpha = alpha)
      sumParts <- sum(vapply(rows[(length(rows) - nPart + 1L):length(rows)], `[[`, 0, "cet"))
      rows[[length(rows) + 1L]] <- data.frame(
        rate = rates[r], alpha = if (gam) alphas[a] else NA, engine = eng, partition = "total",
        nChar = nrow(parts), cet = total,
        revbayes = sum(rb$lnL[rb$kind == kind & rb$case == case]))
      if (abs(total - sumParts) > 1e-8) {
        fail <- c(fail, sprintf("whole matrix differs from the sum of its partitions (rate %g, %s)",
                                rates[r], eng))
      }
      if (!identical(PoweredMk(tt, w, len, engine = eng, alpha = alpha, power = 1), total) ||
          abs(PoweredMk(tt, w, len, engine = eng, alpha = alpha, power = 0.5) - 0.5 * total) > 1e-9) {
        fail <- c(fail, "PoweredMk() does not return t x the package's lnL")
      }
      if (!gam && abs(total - before[r]) > 1e-5) {
        fail <- c(fail, sprintf("Gamma-off total %.6f differs from PartitionedMk()'s %.6f (rate %g)",
                                total, before[r], rates[r]))
      }
      none <- MkLikelihood(tt, LoadCharacters(tt$tipLabel, mkv = "none",
                                              nRateCat = if (gam) 4L else 1L)$characters,
                           len, engine = eng, alpha = alpha)
      cat(sprintf("rate %g, alpha %s, %s: total %.8f, sum of partitions %.8f, Mkv correction %.4f\n",
                  rates[r], if (gam) format(alphas[a]) else "- (no Gamma)", eng, total, sumParts,
                  total - none))
    }
  }
}
res <- do.call(rbind, rows)
res$difference <- res$cet - res$revbayes
print(res, digits = 10, row.names = FALSE)
for (gam in c(FALSE, TRUE)) {
  d <- res$difference[res$partition != "total" & is.na(res$alpha) == !gam]
  cat(sprintf("Gamma %s: max |CET - RevBayes| per partition %.3g over %d values\n",
              if (gam) "on" else "off", max(abs(d)), length(d)))
  if (any(!is.finite(d)) || max(abs(d)) > 1e-6) {
    fail <- c(fail, sprintf("a partition's lnL differs from RevBayes coding = \"variable\" (Gamma %s)",
                            if (gam) "on" else "off"))
  }
}
tot <- res[res$partition == "total", ]
engDiff <- max(abs(tot$cet[tot$engine == "cpp"] - tot$cet[tot$engine == "R"]))
cat(sprintf("max |cpp - R| (total) = %.3g\n", engDiff))
if (engDiff > 1e-8 * max(abs(tot$cet))) fail <- c(fail, "engines disagree")
write.csv(res, file.path(outDir, "mkv-check.csv"), row.names = FALSE)

# ---- Conditioning and time origin (beta = 0 against dnFBDP) --------------------
rbF <- rb$lnL[rb$kind == "fbd"][order(rb$case[rb$kind == "fbd"])]
stopifnot(length(rbF) == nrow(grid))
WithOrigin <- function(tree, origin) {
  tree$stem <- origin - tree$age[[tree$root]]
  tree
}
Shifted <- function(tree, by) {
  tree$age <- tree$age - by
  tree$range <- tree$range - by
  tree
}
# TreeLikelihood() integrates D rescaled to 1 at the start of each branch,
# with absolute tolerance atol (default 1e-14).  On a long branch where D
# decays below atol (grid case 6: psi 0.3 and a 100-Myr stem) the default
# loses accuracy: there it is 9.3 nats too high, and the engines disagree
# (column cetDefaultTol).  TempoMCMC() does not use TreeLikelihood() (its
# target is .LogAugmented(), which integrates no D), so the check uses
# atol = 1e-30, at which both engines agree with RevBayes.
atol <- 1e-30
ModelNone <- function() {
  m <- CetModel(0, 0)
  m$condition <- "none"
  m
}
fbd <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
  g <- grid[j, ]
  par <- list(lambda = g$lambda, mu = g$mu, psi = g$psi, theta = 0.05, s = 0.6)
  tj <- WithOrigin(tt, g$origin)
  youngest <- min(tt$age[seq_len(tt$nTip)])
  data.frame(case = j, lambda = g$lambda, mu = g$mu, psi = g$psi, originMa = g$origin + windowEnd,
             revbayes = rbF[j],
             cet = TreeLikelihood(tj, par, CetModel(0, 0), atol = atol),
             cetR = TreeLikelihood(tj, par, CetModel(0, 0), atol = atol, engine = "R"),
             cetDefaultTol = TreeLikelihood(tj, par, CetModel(0, 0)),
             controlNone = TreeLikelihood(tj, par, ModelNone(), atol = atol),
             controlYoungestTip = TreeLikelihood(Shifted(tj, youngest - 1e-9), par, CetModel(0, 0),
                                                 atol = atol))
}))
Spread <- function(x) diff(range(x))
spread <- c(cet = Spread(fbd$revbayes - fbd$cet), cetR = Spread(fbd$revbayes - fbd$cetR),
            controlNone = Spread(fbd$revbayes - fbd$controlNone),
            controlYoungestTip = Spread(fbd$revbayes - fbd$controlYoungestTip))
nBif <- sum(!tt$sampled)
offset <- -lgamma(tt$nTip + 1) + nBif * log(2)
fbd$expectedOffset <- offset
fbd$residual <- fbd$revbayes - fbd$cet - offset
print(fbd, digits = 12, row.names = FALSE)
cat(sprintf("RevBayes - CET: mean %.12g; -ln(%d!) + %d ln 2 = %.12g; max |residual| %.3g\n",
            mean(fbd$revbayes - fbd$cet), tt$nTip, nBif, offset, max(abs(fbd$residual))))
print(spread, digits = 4)
tol <- max(1e-6, 1e-8 * max(abs(fbd$revbayes)))
if (max(abs(fbd$residual)) > tol || max(abs(fbd$revbayes - fbd$cetR - offset)) > tol) {
  fail <- c(fail, "RevBayes - CET differs from -ln(n!) + nBif ln 2")
}
if (spread[["cet"]] > tol || spread[["cetR"]] > tol) {
  fail <- c(fail, "beta = 0 tree density is not RevBayes's dnFBDP up to a constant")
}
if (spread[["controlNone"]] < 1e-2 || spread[["controlYoungestTip"]] < 1e-2) {
  fail <- c(fail, "a negative control passed: the conditioning check has no power")
}
fbd$spreadTolerance <- tol
write.csv(fbd, file.path(outDir, "fbd-check.csv"), row.names = FALSE)

# ---- Cost ------------------------------------------------------------------------
len <- 0.0236 * dur
nRep <- as.integer(Sys.getenv("CET_MKV_REPS", "200"))
# Without the partial-likelihood cache (#110), which would make every repetition
# on the same tree a cache hit
Ms <- function(f, n = nRep) {
  old <- options(EvoTempo.mkCache = FALSE)
  on.exit(options(old))
  f()
  1000 * system.time(for (i in seq_len(n)) f())[["elapsed"]] / n
}
noneOn <- LoadCharacters(tt$tipLabel, mkv = "none")$characters
timing <- data.frame(
  part = c("MkLikelihood, Gamma 4 categories (cpp; the chains)", "MkLikelihood, no Gamma (cpp)",
           "MkLikelihood, Gamma 4, coding all (cpp)", "MkLikelihood, Gamma 4 (R)",
           "MkLikelihood, no Gamma (R)"),
  msPerEvaluation = c(
    Ms(function() MkLikelihood(tt, whole$on, len, engine = "cpp", alpha = 0.6)),
    Ms(function() MkLikelihood(tt, whole$off, len, engine = "cpp")),
    Ms(function() MkLikelihood(tt, noneOn, len, engine = "cpp", alpha = 0.6)),
    Ms(function() MkLikelihood(tt, whole$on, len, engine = "R", alpha = 0.6), max(5L, nRep %/% 20L)),
    Ms(function() MkLikelihood(tt, whole$off, len, engine = "R"), max(5L, nRep %/% 20L))))
print(timing, digits = 4, row.names = FALSE)
write.csv(timing, file.path(outDir, "mkv-check-timing.csv"), row.names = FALSE)
if (length(fail)) stop("check-mkv.R FAILED: ", paste(fail, collapse = "; "))
cat("check-mkv.R: all checks passed\n")
