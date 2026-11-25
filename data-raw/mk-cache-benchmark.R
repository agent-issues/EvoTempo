# Time per proposal of TempoMCMC() on the penguin data, by move type, with and
# without the cached character likelihood (#110 move 2;
# dev/plans/mk-likelihood-cache.md).
#
#   Rscript data-raw/mk-cache-benchmark.R [beta gamma nIter count out.csv]
#
# Run from the package root with EvoTempo installed.  "before" is the one-shot
# path: options(EvoTempo.mkCache = FALSE) and the R Mkv correction
# (ConstantPatternLogLik()); "after" is the cached pruner for the partitions and
# the constant-pattern columns.  Each move type runs alone (`count` proposals
# per iteration, `nIter` iterations, fixed parameters); the cost of a run with
# no moves (set-up, per-iteration statistics) is subtracted.  "parameters" is
# the parameter sweep alone (one proposal per free parameter).  Both modes run
# from the same seed; `same` reports whether their samples agree to 1e-10 (the
# R and C++ Mkv corrections differ by round-off, so not bitwise).
suppressMessages(library("EvoTempo"))
source("data-raw/empirical/penguins/cet/common.R")
args <- suppressWarnings(as.numeric(commandArgs(trailingOnly = TRUE)[1:4]))
beta <- if (!is.na(args[[1]])) args[[1]] else 1
gamma <- if (!is.na(args[[2]])) args[[2]] else 1
nIter <- if (!is.na(args[[3]])) args[[3]] else 2
count <- if (!is.na(args[[4]])) args[[4]] else 5
outFile <- commandArgs(trailingOnly = TRUE)[5]

dat <- LoadPenguins()
model <- CetModel(beta, gamma)
priors <- CetPriors()
par <- list(lambda = 0.2, mu = 0.15, psi = 0.05, theta = 0.05, s = 0.2, rate = 0.017)
cachedPattern <- PatternLogLik
Mode <- function(mode) {
  options(EvoTempo.mkCache = mode == "after")
  PatternLogLik <<- if (mode == "after") {
    cachedPattern
  } else {
    function(tree, pat, lengths, engine) ConstantPatternLogLik(tree, pat, lengths)
  }
}
Run <- function(moves, update = FALSE) {
  set.seed(1)
  t0 <- Sys.time()
  run <- TempoMCMC(dat$tree, model, par, priors, nIter = nIter, burnin = 0, thin = 1,
                   nSlice = 20, moves = moves, characters = dat$characters,
                   agePriors = CetAgePriors(), updateParameters = update,
                   minAcceptance = 0)
  list(secs = as.numeric(Sys.time() - t0, units = "secs"), samples = run$samples)
}
moveTypes <- c("slide", "saJump", "nodeAge", "fossilAge", "originAge", "node", "tip",
               "branch", if (gamma != 0) "shift", if (gamma == 0) "tree",
               "parameters")
rows <- list()
for (mode in c("before", "after")) {
  Mode(mode)
  get(".MkCache", .ns)(clear = TRUE)
  Run(c(slide = 1, node = 1))  # warm-up (byte compilation, first allocations)
  none <- Run(c(node = 0))[["secs"]]
  for (mv in moveTypes) {
    res <- if (mv == "parameters") Run(c(node = 0), update = TRUE) else Run(setNames(count, mv))
    nProp <- nIter * if (mv == "parameters") length(priors) else count
    rows[[length(rows) + 1L]] <- data.frame(
      move = mv, mode = mode, msPerProposal = 1000 * (res$secs - none) / nProp,
      samples = I(list(res$samples)))
  }
}
Mode("after")
res <- do.call(rbind, rows)
before <- res[res$mode == "before", ]
after <- res[res$mode == "after", ]
out <- data.frame(
  beta = beta, gamma = gamma, move = before$move,
  msBefore = round(before$msPerProposal, 2), msAfter = round(after$msPerProposal, 2),
  speedUp = round(before$msPerProposal / after$msPerProposal, 2),
  same = mapply(function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-10)),
                before$samples, after$samples))
print(out, row.names = FALSE)
if (!is.na(outFile)) {
  write.csv(out, outFile, row.names = FALSE)
}
