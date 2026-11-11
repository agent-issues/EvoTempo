# Gamma off (#97): seeded TempoMCMC() chains on a character matrix must be
# identical() under main and under the Gamma branch, in both engines, with
# every move type.  Run once per installed build, then compare:
#
#   Rscript data-raw/gamma-acrv-identity.R run <out.rds>   # with the build first on .libPaths()
#   Rscript data-raw/gamma-acrv-identity.R compare <main.rds> <branch.rds>
args <- commandArgs(TRUE)
if (args[[1]] == "compare") {
  a <- readRDS(args[[2]])
  b <- readRDS(args[[3]])
  stopifnot(identical(names(a), names(b)))
  same <- vapply(names(a), function(k) identical(a[[k]], b[[k]]), TRUE)
  print(data.frame(run = names(a), identical = same,
                   nSample = vapply(a, function(r) nrow(r$samples), 0L)),
        row.names = FALSE)
  cat(if (all(same)) "ALL IDENTICAL\n" else "DIFFERENT\n")
  quit(status = if (all(same)) 0 else 1)
}
library(EvoTempo)
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils with age ranges;
# f a sampled ancestor (the fixture of tests/testthat/test-Engine.R)
phylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
  class = "phylo")
tree <- TempoTree(phylo, ranges = cbind(c(0.2, 0, 0.5, 0.5, 0), c(2.5, 0, 3, 3.5, 0)))
par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.4)
priors <- c(TempoPriors()[c("lambda", "mu", "psi", "theta", "s")],
            rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
set.seed(1)
chars <- matrix(sample(0:2, 5 * 8, replace = TRUE), 5)
chars[2, 3] <- NA
runs <- list()
for (engine in c("cpp", "R")) {
  for (gamma in c(0, 1)) {
    model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, gamma = gamma)
    moves <- c(tree = 1, node = 1, tip = 1, branch = 1, joint = 1, shift = 1,
               nodeAge = 1, fossilAge = 1, slide = 1, saJump = 1)
    set.seed(2)
    run <- TempoMCMC(tree, model, par, priors, nIter = 30, burnin = 20, nSlice = 4,
                     characters = chars, nState = 3, moves = moves,
                     minAcceptance = 0, engine = engine)
    runs[[paste0(engine, "_gamma", gamma)]] <-
      run[c("samples", "occupancy", "acceptance", "burninAcceptance", "moves")]
  }
}
saveRDS(runs, args[[2]])
cat("saved", length(runs), "runs to", args[[2]], "\n")
