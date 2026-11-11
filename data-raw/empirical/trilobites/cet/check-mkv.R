# Trilobite CET pipeline: fixed-tree checks of the character likelihood and of
# the tree density's conditioning against RevBayes, on the start tree
# (data/start-tree-1.nwk).
#
#   Rscript data-raw/empirical/trilobites/cet/check-mkv.R write <oracle dir>
#   (run RevBayes: revbayes/check-oracle.sh does both steps and this one)
#   Rscript data-raw/empirical/trilobites/cet/check-mkv.R compare <oracle dir> [out dir]
#
# `write` writes the oracle inputs for revbayes/check-oracle.Rev: the start
# tree with branch lengths rate x duration at strict-clock rates 0.005, 0.0236
# (the comparator's posterior median clock) and 0.05, the same tree in Myr,
# and a grid of (lambda, mu, psi, origin).
#
# `compare` checks, and stops with an error if a check fails:
# 1. Characters.  PartitionedMk() per RevBayes partition (u2-u16, o4, o5) and
#    in total must equal RevBayes dnPhyloCTMC(coding = "variable") with the
#    comparator's Q and no site-rate variation, read from the comparator's own
#    matrices (../revbayes/prepare-data.R, which treats polymorphisms as
#    ambiguous): |difference| <= 1e-4 per partition, both engines.
# 2. The same total by an independent route (AllByR(): every unordered
#    character through the R pruner SymLogLik(), none through the C++ engine)
#    to 1e-8; the engines agree to 1e-8 relative; power t gives t x lnL.
# 3. Conditioning and time origin.  With beta = 0 the dCET tree density equals
#    the FBD density for any generator (spec section 7.2), so TreeLikelihood()
#    under CetModel(0, 0) (rho = 0, condition = "survival", start = "origin",
#    ages before 485.4 Ma) must differ from RevBayes dnFBDP(rho = 0,
#    condition = "sampling") by a constant over the grid of parameters and
#    origins (the two use different labelling conventions, which give a
#    constant): spread (max - min of the difference) <= 1e-6 (relative to
#    max |lnL| 1e-8 when larger), with both engines.  (FBDTreeLikelihood(),
#    the closed form, gives NaN here: with rho = 0 and no extant tip it adds
#    0 x log(0); TreeLikelihood() is unaffected.)  Negative controls, which
#    must FAIL that
#    test by a wide margin (spread > 1e-2): condition = "none"; and the
#    present at the youngest tip (486.0 Ma) instead of 485.4 Ma.
# 4. Cost per PartitionedMk() evaluation (ms), by part and engine.
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
  rec <- regmatches(txt, gregexpr("(mkv|fbd)\t[0-9]+\t[A-Za-z0-9]+\t-?[0-9.]+(e[-+]?[0-9]+)?", txt))[[1]]
  f <- do.call(rbind, strsplit(rec, "\t"))
  data.frame(kind = f[, 1], case = as.integer(f[, 2]), partition = f[, 3], lnL = as.numeric(f[, 4]))
}
rb <- ReadOracle(file.path(oracleDir, "revbayes-oracle.tsv"))
stopifnot(sum(rb$kind == "mkv") == 30, sum(rb$kind == "fbd") == nrow(grid))
fail <- character()
cells <- dat$cells
parts <- dat$parts
PartLL <- function(part, len, engine) {
  cols <- which(parts$part == part)
  x <- matrix(0L, nrow(cells), length(cols))
  attr(x, "cet") <- BuildCharacters(cells[, cols, drop = FALSE], parts[cols, ])
  PartitionedMk(tt, x, len, engine = engine, power = 1, mkv = "pattern")
}
rows <- list()
for (r in seq_along(rates)) {
  len <- rates[r] * dur
  for (eng in c("cpp", "R")) {
    for (part in unique(parts$part)) {
      ref <- rb$lnL[rb$kind == "mkv" & rb$case == r & rb$partition == part]
      rows[[length(rows) + 1L]] <- data.frame(rate = rates[r], engine = eng, partition = part,
                                              nChar = sum(parts$part == part),
                                              cet = PartLL(part, len, eng), revbayes = ref)
    }
    total <- PartitionedMk(tt, dat$characters, len, engine = eng, power = 1)
    rows[[length(rows) + 1L]] <- data.frame(
      rate = rates[r], engine = eng, partition = "total", nChar = nrow(parts), cet = total,
      revbayes = sum(rb$lnL[rb$kind == "mkv" & rb$case == r]))
    byR <- AllByR(tt, cells, parts, len)
    half <- PartitionedMk(tt, dat$characters, len, engine = eng, power = 0.5)
    none <- PartitionedMk(tt, dat$characters, len, engine = eng, power = 1, mkv = "none")
    if (abs(byR - total) > 1e-8) fail <- c(fail, sprintf("AllByR differs at rate %g (%s)", rates[r], eng))
    if (abs(half - 0.5 * total) > 1e-9) fail <- c(fail, "power does not scale the corrected lnL")
    cat(sprintf("rate %g %s: total %.8f, by R pruners %.8f, Mkv correction %.4f\n",
                rates[r], eng, total, byR, total - none))
  }
}
res <- do.call(rbind, rows)
res$difference <- res$cet - res$revbayes
print(res, digits = 10, row.names = FALSE)
if (any(abs(res$difference[res$partition != "total"]) > 1e-4)) {
  fail <- c(fail, "a partition's lnL differs from RevBayes coding = \"variable\"")
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
print(fbd, digits = 12, row.names = FALSE)
cat("constant RevBayes - CET:", format(mean(fbd$revbayes - fbd$cet), digits = 12), "\n")
print(spread, digits = 4)
tol <- max(1e-6, 1e-8 * max(abs(fbd$revbayes)))
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
cm <- attr(dat$characters, "cet")
nRep <- as.integer(Sys.getenv("CET_MKV_REPS", "200"))
Ms <- function(f) {
  f()
  1000 * system.time(for (i in seq_len(nRep)) f())[["elapsed"]] / nRep
}
timing <- data.frame(
  part = c("PartitionedMk (cpp)", "PartitionedMk (R)", "C++ groups (MkLikelihood)",
           "ambiguous + ordered + Mkv columns (CetPrune, C++)", "ambiguous (SymLogLik, R)",
           "ordered (OrdLogLik, R)", "Mkv patterns (ConstantPatternLogLik, R)", "no correction (cpp)"),
  msPerEvaluation = c(
    Ms(function() PartitionedMk(tt, dat$characters, len, engine = "cpp", power = 1)),
    Ms(function() PartitionedMk(tt, dat$characters, len, engine = "R", power = 1)),
    Ms(function() for (p in cm$cpp) .mkOriginal(tt, p$x, len, p$k, "cpp")),
    Ms(function() {
      kids <- tt$kids
      storage.mode(kids) <- "integer"
      CetPrune(as.integer(tt$postorder), kids, tt$nTip, tt$root, len, cm$pack$tip, cm$pack$k,
               cm$pack$ordId, cm$pack$ordV, cm$pack$ordLam)
    }),
    Ms(function() SymLogLik(tt, cm$sym, len)),
    Ms(function() lapply(cm$ord, function(o) OrdLogLik(tt, o, len))),
    Ms(function() ConstantPatternLogLik(tt, cm$mkv, len)),
    Ms(function() PartitionedMk(tt, dat$characters, len, engine = "cpp", power = 1, mkv = "none"))))
print(timing, digits = 4, row.names = FALSE)
write.csv(timing, file.path(outDir, "mkv-check-timing.csv"), row.names = FALSE)
if (length(fail)) stop("check-mkv.R FAILED: ", paste(fail, collapse = "; "))
cat("check-mkv.R: all checks passed\n")
