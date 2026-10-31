# Start-age moves (#18) on a fixed tree at fixed parameters: the posterior of
# the root age (start = "root", move `rootAge`) or of the origin age
# (start = "origin", move `originAge`) from TempoMCMC() chains with path
# moves, against the exact one-dimensional posterior, prior x TreeLikelihood()
# on a grid of start ages (trapezoid CDF), as in data-raw/rung4b-one-parameter.R.
# Tree: ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7), root age 4, stem 0.6;
# a, c fossils; f a sampled ancestor (the fixture of tests/testthat/test-Origin.R).
# z-scores use the spread of the CDF at each exact quantile over independent
# chains.
# Usage: Rscript data-raw/rung5-age-oracle.R <root|origin> [nIter] [nChain]
# Environment: EVOTEMPO_RUNG5_OUT (default data-raw/rung5-age-results)
suppressMessages(library(EvoTempo))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]
stopifnot(case %in% c("root", "origin"))
nIter <- if (length(args) > 1) as.integer(args[[2]]) else 20000L
nChain <- if (length(args) > 2) as.integer(args[[3]]) else 16L
outDir <- Sys.getenv("EVOTEMPO_RUNG5_OUT", "data-raw/rung5-age-results")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
origin <- case == "origin"
phylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
  class = "phylo")
tree <- TempoTree(phylo)
model <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5,
                   start = if (origin) "origin" else "root")
par <- list(lambda = 0.5, mu = 0.2, psi = 0.15, theta = 0.2, s = 0.6)
ageRange <- if (origin) c(4.2, 8) else c(3.4, 7)
prior <- function(a) stats::dunif(a, ageRange[1], ageRange[2], log = TRUE)
agePriors <- if (origin) list(origin = prior) else list(root = prior)
root <- tree$root
WithAge <- function(a) {
  tr <- tree
  if (origin) tr$stem <- a - tr$age[[root]] else tr$age[root] <- a
  tr
}
grid <- seq(ageRange[1], ageRange[2], length.out = 321)
cl <- makeCluster(16)
invisible(clusterCall(cl, function() suppressMessages(library(EvoTempo))))
clusterExport(cl, c("tree", "model", "par", "WithAge", "origin", "root", "agePriors",
                    "prior", "ageRange", "nIter"))
logL <- unlist(parLapply(cl, grid, function(a) TreeLikelihood(WithAge(a), par, model)))
t0 <- Sys.time()
chains <- parLapply(cl, seq_len(nChain), function(i) {
  set.seed(18900 + i)
  run <- TempoMCMC(tree, model, par, list(), nIter = nIter, burnin = 1000, thin = 5,
                   nSlice = 16, agePriors = agePriors,
                   moves = c(tree = 1, node = 2, tip = 2, branch = 2,
                             if (origin) c(originAge = 2) else c(rootAge = 2)))
  list(age = run$samples[[if (origin) "originAge" else "rootAge"]],
       acceptance = run$acceptance)
})
runtime <- difftime(Sys.time(), t0, units = "mins")
stopCluster(cl)
lp <- logL + prior(grid)
dens <- exp(lp - max(lp))
cdf <- c(0, cumsum((dens[-1] + dens[-length(dens)]) / 2 * diff(grid)))
cdf <- cdf / max(cdf)
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
out <- do.call(rbind, lapply(probs, function(pr) {
  q <- stats::approx(cdf, grid, pr)$y
  f <- vapply(chains, function(ch) mean(ch$age <= q), 0)
  data.frame(case = case, p = pr, exact = q,
             mcmc = stats::quantile(unlist(lapply(chains, `[[`, "age")), pr, names = FALSE),
             cdfAtExact = mean(f), z = (mean(f) - pr) / (stats::sd(f) / sqrt(length(f))))
}))
exactMean <- sum(grid * dens) / sum(dens)
cm <- vapply(chains, function(ch) mean(ch$age), 0)
out <- rbind(out, data.frame(case = case, p = NA, exact = exactMean, mcmc = mean(cm),
                             cdfAtExact = NA, z = (mean(cm) - exactMean) / (stats::sd(cm) / sqrt(length(cm)))))
acceptance <- colMeans(do.call(rbind, lapply(chains, `[[`, "acceptance")))
saveRDS(list(case = case, grid = grid, logL = logL, chains = chains, out = out,
             acceptance = acceptance, nIter = nIter, runtime = runtime,
             sha = Sys.getenv("EVOTEMPO_SHA")),
        file.path(outDir, paste0("rung5-age-oracle-", case, ".rds")))
write.csv(out, file.path(outDir, paste0("rung5-age-oracle-", case, ".csv")), row.names = FALSE)
print(out, digits = 4, row.names = FALSE)
print(round(acceptance[acceptance > 0], 3))
cat("runtime", format(runtime), "\n")
