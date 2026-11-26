# Time per proposal of TempoMCMC() on the penguin data, by move type, with and
# without the cached character likelihood (#110 move 2;
# dev/plans/mk-likelihood-cache.md).
#
#   Rscript data-raw/mk-cache-benchmark.R [beta gamma nIter count out.csv]
#
# Run from the package root with EvoTempo installed.  Three modes separate the
# two changes of #119:
# - main: one-shot pruning and the R Mkv correction (ConstantPatternLogLik()),
#   the code path of main before #119;
# - cppMkv: one-shot pruning, Mkv correction in C++ (PatternLogLik());
# - cached: as cppMkv, with the cache.
# Each move type runs alone (`count` proposals per iteration, `nIter`
# iterations, fixed parameters); the cost of a run with no moves (set-up,
# per-iteration statistics) is subtracted.  "parameters" is the parameter sweep
# alone (one proposal per free parameter); "schedule" is the production
# schedule (CetMoves()) with the sweep, in ms per iteration.  All modes run
# from the same seed; `same` reports whether their samples agree to 1e-10 (the
# R and C++ Mkv corrections differ by round-off, so not bitwise).  Also prints
# the time of one PartitionedMk() call per mode (cached: one tip edge changed
# per call).
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
cppPattern <- PatternLogLik
modes <- c("main", "cppMkv", "cached")
Mode <- function(mode) {
  options(EvoTempo.mkCache = mode == "cached")
  PatternLogLik <<- if (mode == "main") {
    function(tree, pat, lengths, engine) ConstantPatternLogLik(tree, pat, lengths)
  } else {
    cppPattern
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
               "parameters", "schedule")
schedule <- CetMoves(dat$tree, dat$characters, model)
tt <- get(".ForModel", .ns)(dat$tree, model)
set.seed(2)
lens <- stats::runif(tt$nTip + tt$nNode, 0.01, 0.3)
CallMs <- function(reps = 200) {
  t0 <- Sys.time()
  for (i in seq_len(reps)) {
    l <- lens
    l[1 + i %% tt$nTip] <- lens[1 + i %% tt$nTip] * (1 + 1e-3 * i)
    PartitionedMk(tt, dat$characters, l)
  }
  1000 * as.numeric(Sys.time() - t0, units = "secs") / reps
}
callMs <- numeric(0)
rows <- list()
for (mode in modes) {
  Mode(mode)
  get(".MkCache", .ns)(clear = TRUE)
  Run(c(slide = 1, node = 1))  # warm-up (byte compilation, first allocations)
  callMs[[mode]] <- CallMs()
  none <- Run(c(node = 0))[["secs"]]
  for (mv in moveTypes) {
    res <- switch(mv,
                  parameters = Run(c(node = 0), update = TRUE),
                  schedule = Run(schedule, update = TRUE),
                  Run(setNames(count, mv)))
    nProp <- nIter * switch(mv, parameters = length(priors), schedule = 1, count)
    rows[[length(rows) + 1L]] <- data.frame(
      move = mv, mode = mode, msPerProposal = 1000 * (res$secs - none) / nProp,
      samples = I(list(res$samples)))
  }
}
Mode("cached")
res <- do.call(rbind, rows)
By <- function(mode) res[res$mode == mode, ]
out <- data.frame(
  beta = beta, gamma = gamma, move = By("main")$move,
  msMain = round(By("main")$msPerProposal, 2),
  msCppMkv = round(By("cppMkv")$msPerProposal, 2),
  msCached = round(By("cached")$msPerProposal, 2),
  mkvGain = round(By("main")$msPerProposal / By("cppMkv")$msPerProposal, 2),
  cacheGain = round(By("cppMkv")$msPerProposal / By("cached")$msPerProposal, 2),
  total = round(By("main")$msPerProposal / By("cached")$msPerProposal, 2),
  same = mapply(function(a, b, c) {
    isTRUE(all.equal(a, b, tolerance = 1e-10)) && isTRUE(all.equal(a, c, tolerance = 1e-10))
  }, By("main")$samples, By("cppMkv")$samples, By("cached")$samples))
print(out, row.names = FALSE)
cat("ms per PartitionedMk() call:", paste(names(callMs), round(callMs, 2), sep = " ",
                                           collapse = ", "), "
")
if (!is.na(outFile)) {
  write.csv(out, outFile, row.names = FALSE)
  write.csv(data.frame(mode = names(callMs), msPerCall = callMs),
            sub("\\.csv$", "-call.csv", outFile), row.names = FALSE)
}
