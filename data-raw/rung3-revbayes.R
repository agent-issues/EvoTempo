# Rung 3 (#15): TreeLikelihood() vs RevBayes dnCDBDP.
#
# Trees: the rung-3 castor trees (12 extant tips; 10 tips with fossils), the
# collapse tree with a sampled ancestor (crown form), and the fossil tree with
# two fossil tips turned into sampled ancestors.  Node ages are rounded to
# multiples of 2^-10 so that tip ages, and hence RevBayes' extant/fossil
# classification, are exact.  Parameters: the 20 castor points of
# rung3-validation.R, psi = 0.15 per bin, rho = 0.7, removal 0 (dnCDBDP has
# no removal), at K = 16 (z in [-3, 3]) and K = 26 (z in [-5, 5]).
#
# Reconciliation (see StateDependentSpeciationExtinctionProcess.cpp):
#   RevBayes multiplies lambda_i into the root and weights root states by pi;
#   EvoTempo omits lambda at the root, so it is given root weights
#   w = pi * lambda / sum(pi * lambda).  Then
#     condition "time":     RB = TL(condition "none", w) + log sum(pi lambda) + shape
#     condition "survival": RB = TL(condition "survival", w, joint = FALSE) + shape
#   where shape = (nTip - nSA - 1) log 2 - log(nExtant!) (lnProbTreeShape, with
#   RevBayes' Stirling approximation to log(nExtant!)),
#   and RevBayes' survival is "any sampled descendant" (E includes psi).
#   EvoTempo's crown "extant" conditioning has no dnCDBDP counterpart.
#
# Run from the package root: Rscript data-raw/rung3-revbayes.R [nTimeSlices...]
# Environment: REVBAYES (executable), REVBAYES_POINTS (e.g. "1,5,9"),
# REVBAYES_CONDITIONS ("time,survival"), REVBAYES_GRIDS ("K16,K26,K16b0"),
# REVBAYES_TAG (suffix of the result file), EVOTEMPO_NO_COMPILE (load without
# compiling and use the R engine).  The stored results were computed with the
# R engine; data-raw/rung3-revbayes-engines.R re-evaluates the EvoTempo side
# with each available engine.
# Grid K16b0 has beta = 0 (state-independent rates), where the sampled-ancestor
# treatment of dnCDBDP can be checked against the closed-form FBD E(t).
noCompile <- nzchar(Sys.getenv("EVOTEMPO_NO_COMPILE"))
suppressWarnings(suppressMessages(
  pkgload::load_all(".", quiet = TRUE, compile = if (noCompile) FALSE else NA)))
if (noCompile) options(EvoTempo.engine = "R")
library(ape)
args <- commandArgs(trailingOnly = TRUE)
slices <- if (length(args)) as.numeric(args) else c(500, 5000, 50000)
rb <- Sys.getenv("REVBAYES", "rb")
outDir <- "data-raw/rung3-results/revbayes"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Round node ages to multiples of 2^-10; edge lengths are then exact
Dyadic <- function(tree) {
  depth <- node.depth.edgelength(tree)
  age <- round((max(depth) - depth) * 1024) / 1024
  tree$edge.length <- age[tree$edge[, 1]] - age[tree$edge[, 2]]
  tree
}
# Turn fossil tip `tip` into a sampled ancestor: its parent moves down to it
ToAncestor <- function(tree, tip) {
  e <- which(tree$edge[, 2] == tip)
  p <- tree$edge[e, 1]
  sib <- which(tree$edge[, 1] == p & tree$edge[, 2] != tip)
  pe <- which(tree$edge[, 2] == p)
  shift <- tree$edge.length[e]
  stopifnot(tree$edge.length[sib] > shift, length(pe) == 1)
  tree$edge.length[e] <- 0
  tree$edge.length[sib] <- tree$edge.length[sib] - shift
  tree$edge.length[pe] <- tree$edge.length[pe] + shift
  tree
}
set.seed(3)
ultra <- rcoal(12)
ultra$edge.length <- ultra$edge.length * 3
set.seed(11)
nonUltra <- rtree(10)
nonUltra$edge.length <- nonUltra$edge.length * 2
trees <- list(extant = Dyadic(ultra), fossil = Dyadic(nonUltra),
              sa = Dyadic(read.tree(text = "((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7);")))
# Two sampled ancestors on the fossil tree: fossil tips whose sibling branch is
# long enough, excluding root daughters
ft <- trees$fossil
tt <- TempoTree(ft)
cand <- which(tt$tipType == "fossil")
done <- 0
for (tip in cand) {
  e <- which(ft$edge[, 2] == tip)
  p <- ft$edge[e, 1]
  if (p == Ntip(ft) + 1) next
  sib <- which(ft$edge[, 1] == p & ft$edge[, 2] != tip)
  if (ft$edge.length[sib] > ft$edge.length[e] + 1 / 64 && done < 2) {
    ft <- ToAncestor(ft, tip)
    done <- done + 1
  }
}
trees$sa2 <- ft
for (nm in names(trees)) {
  write.tree(trees[[nm]], file.path(outDir, paste0(nm, ".nwk")), digits = 16)
}
info <- do.call(rbind, lapply(names(trees), function(nm) {
  t <- TempoTree(trees[[nm]])
  data.frame(tree = nm, nTip = t$nTip, nSA = sum(t$tipType == "ancestor"),
             nExtant = sum(t$tipType == "extant"), rootAge = max(t$age))
}))
# RevBayes' RbMath::lnFactorial(): a truncated Stirling series with no exact
# small-n branch (error -5.0e-4 at n = 1, -2.1e-5 at n = 2), used as is in
# lnProbTreeShape(); the shape term reproduces it, and its error is recorded
RbLnFactorial <- function(n) {
  ifelse(n == 0, 0, (n + 0.5) * log(n) - n + 0.918938533204672722 +
           (1 / n) * (1 / 12 + (1 / n)^2 * (-1 / 360)))
}
info$shape <- (info$nTip - info$nSA - 1) * log(2) - RbLnFactorial(info$nExtant)
info$stirlingError <- RbLnFactorial(info$nExtant) - lfactorial(info$nExtant)

set.seed(5)
points <- data.frame(lambda = runif(20, 0.3, 1), mu = runif(20, 0, 0.3),
                     theta = runif(20, 0.05, 1), s = runif(20, 0.5, 1.5))
grids <- list(K16 = list(nBin = 16, zRange = c(-3, 3), beta = 1),
              K26 = list(nBin = 26, zRange = c(-5, 5), beta = 1),
              K16b0 = list(nBin = 16, zRange = c(-3, 3), beta = 0))
Env <- function(name, default) {
  v <- Sys.getenv(name)
  if (nzchar(v)) strsplit(v, ",")[[1]] else default
}
cases <- expand.grid(point = as.integer(Env("REVBAYES_POINTS", seq_len(nrow(points)))),
                     tree = names(trees),
                     condition = Env("REVBAYES_CONDITIONS", c("time", "survival")),
                     grid = Env("REVBAYES_GRIDS", c("K16", "K26")),
                     stringsAsFactors = FALSE)
tag <- Sys.getenv("REVBAYES_TAG")
cases$id <- seq_len(nrow(cases))

# EvoTempo side, and the Rev code for each case
Num <- function(x) formatC(x, digits = 17, format = "g")
Vec <- function(x) paste0("v(", paste(Num(x), collapse = ", "), ")")
revCase <- character(nrow(cases))
cases$evo <- NA_real_
cases$saOffset <- NA_real_
for (i in seq_len(nrow(cases))) {
  cs <- cases[i, ]
  g <- grids[[cs$grid]]
  par <- c(as.list(points[cs$point, ]), psi = 0.15)
  base <- DCETModel(nBin = g$nBin, zRange = g$zRange, rho = 0.7, beta = g$beta)
  set <- .Setup(par, base)
  w <- set$prior * set$lambda
  model <- DCETModel(nBin = g$nBin, zRange = g$zRange, rho = 0.7, removal = 0, beta = g$beta,
                     rootPrior = w / sum(w),
                     condition = if (cs$condition == "time") "none" else "survival")
  tl <- TreeLikelihood(trees[[cs$tree]], par, model)
  cases$evo[i] <- tl + info$shape[info$tree == cs$tree] +
    if (cs$condition == "time") log(sum(w)) else 0
  # Predicted dnCDBDP offset if each sampled ancestor carries psi * E(t)
  # instead of psi (exact when rates do not depend on the bin)
  tt <- TempoTree(trees[[cs$tree]])
  saAge <- tt$age[seq_len(tt$nTip)][tt$tipType == "ancestor"]
  cases$saOffset[i] <- if (length(saAge) && g$beta == 0) {
    sum(log(FBDExtinction(saAge, set$lambda[1], set$mu[1], 0.15, 0.7)))
  } else NA
  Q <- set$Q
  rows <- vapply(seq_len(nrow(Q)), function(r) {
    q <- Q[r, ]
    q[r] <- 0
    paste0("[", paste(Num(q), collapse = ", "), "]")
  }, "")
  revCase[i] <- paste0(
    "Q <- fnFreeK([", paste(rows, collapse = ", "), "], rescaled=false)\n",
    "x ~ dnCDBDP(rootAge=tree_", cs$tree, ".rootAge(), speciationRates=",
    Vec(set$lambda), ", extinctionRates=", Vec(set$mu), ", psi=", Vec(set$psi),
    ", Q=Q, delta=1.0, pi=simplex(", Vec(set$prior), "), rho=0.7, condition=\"",
    cs$condition, "\", nTimeSlices=NTS)\n",
    "x.clamp(tree_", cs$tree, ")\n",
    "print(\"CASE ", cs$id, "\")\nprint(x.lnProbability())\n")
}

# RevBayes side: one run per nTimeSlices
header <- paste0("setOption(\"outputPrecision\", \"16\")\n",
                 paste0("tree_", names(trees), " <- readTrees(\"",
                        file.path(outDir, paste0(names(trees), ".nwk")),
                        "\")[1]\n", collapse = ""))
runs <- do.call(rbind, lapply(slices, function(nts) {
  script <- file.path(outDir, sprintf("dnCDBDP%s-%d.Rev", tag, as.integer(nts)))
  writeLines(c(header, gsub("NTS", format(nts, scientific = FALSE), revCase), "q()"),
             script)
  t0 <- Sys.time()
  out <- system2(rb, script, stdout = TRUE, stderr = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # Values are printed on the line after their label: string concatenation
  # in Rev rounds to 6 significant digits, print() honours outputPrecision
  at <- grep("^ *CASE [0-9]+ *$", out)
  data.frame(id = as.integer(sub("^ *CASE ([0-9]+) *$", "\\1", out[at])),
             nTimeSlices = nts, rb = suppressWarnings(as.numeric(out[at + 1])),
             seconds = secs)
}))
res <- merge(cases, runs, by = "id")
res$delta <- res$rb - res$evo
saveRDS(list(cases = cases, runs = runs, res = res, info = info, points = points,
             trees = trees), file.path(outDir, paste0("rung3-revbayes", tag, ".rds")))
print(info)
print(aggregate(abs(delta) ~ grid + tree + condition + nTimeSlices, res, max), digits = 4)
