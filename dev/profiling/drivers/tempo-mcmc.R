# Profiling driver: TempoMCMC() on a rung-6-style replicate (issue #21 design):
# beta = gamma = 1, K = 16, z in [-3, 3], crown age 8, 20-40 tips, 30 binary
# characters, nSlice = 16; path moves of the rung-6 script plus `joint` and the
# four tree moves that the full-size SBC adds.
#
# Usage: Rscript dev/profiling/drivers/tempo-mcmc.R [nIter] [moveSet]
#   moveSet: "all" (default), "rung6" (path moves only), or a single move name
# Library: EVOTEMPO_LIB (default: installed EvoTempo).
lib <- Sys.getenv("EVOTEMPO_LIB", "")
if (nzchar(lib)) {
  library(EvoTempo, lib.loc = lib)
} else {
  library(EvoTempo)
}
args <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args)) as.integer(args[[1]]) else 20L
moveSet <- if (length(args) > 1) args[[2]] else "all"

source("dev/profiling/drivers/fixture-rung6.R")  # run from the package root

moves <- switch(moveSet,
  all = c(tree = 1, node = 3, tip = 3, branch = 3, joint = 1, nodeAge = 1,
          fossilAge = 1, slide = 1, saJump = 1),
  rung6 = c(tree = 1, node = 3, tip = 3, branch = 3),
  setNames(1, moveSet))

rprof <- Sys.getenv("EVOTEMPO_RPROF", "")  # profile the chain alone
set.seed(5813)
if (nzchar(rprof)) Rprof(rprof, interval = 0.02)
t0 <- proc.time()
run <- TempoMCMC(tt, model, truth, priors, nIter = nIter, burnin = 0, thin = 1,
                 nSlice = 16, moves = moves, characters = chars)
elapsed <- (proc.time() - t0)[["elapsed"]]
if (nzchar(rprof)) Rprof(NULL)
cat("Elapsed:", round(elapsed, 2), "s;", round(1000 * elapsed / nIter, 1),
    "ms/iteration\n")
cat("logPost[last]:", format(run$samples$logPost[nIter], digits = 12), "\n")
print(round(run$acceptance, 3))
out <- Sys.getenv("EVOTEMPO_OUT", "")  # save samples, to compare two builds
if (nzchar(out)) saveRDS(run[c("samples", "acceptance")], out)
