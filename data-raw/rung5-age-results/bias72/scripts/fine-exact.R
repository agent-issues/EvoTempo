# verify43c / #72: 2561-point Simpson-rule exact posterior mean (extant conditioning, power
# design) for each arm, and the script's own trapezoid mean/sd on its grid (641 or 1601 points)
# as a transcription check against the chains' exact mean. Definitions copied from
# data-raw/rung5-age-oracle.R at dc47aed (and verify/oracle-gamma.R for arm B).
suppressMessages(library(EvoTempo)); library(parallel)
phylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
  class = "phylo")
tree <- TempoTree(phylo); root <- tree$root
par <- list(lambda = 1, mu = 1, psi = 0.15, theta = 0.4, s = 1)
cl <- makeCluster(16)
invisible(clusterCall(cl, function() suppressMessages(library(EvoTempo))))
res <- list()
for (arm in c("A", "B", "C")) {
  origin <- arm == "C"; gammaPrior <- arm == "B"
  model <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.35,
                     start = if (origin) "origin" else "root", condition = "extant", joint = FALSE)
  ageRange <- if (origin) c(4.2, 8) else if (gammaPrior) c(3.4, 12) else c(3.4, 7)
  prior <- if (gammaPrior) {
    function(a) ifelse(a > 3.4 & a < 12, stats::dgamma(a - 3.4, 3, scale = 0.6, log = TRUE), -Inf)
  } else function(a) stats::dunif(a, ageRange[1], ageRange[2], log = TRUE)
  LL <- function(grid) unlist(parLapply(cl, grid, function(a, tree, root, origin, par, model) {
    tr <- tree; if (origin) tr$stem <- a - tr$age[[root]] else tr$age[root] <- a
    TreeLikelihood(tr, par, model) }, tree, root, origin, par, model))
  Dens <- function(grid, logL) { lp <- logL + prior(grid); exp(lp - max(lp)) }
  g1 <- seq(ageRange[1], ageRange[2], length.out = if (gammaPrior) 1601 else 641)
  d1 <- Dens(g1, LL(g1)); w1 <- c(0.5, rep(1, length(g1) - 2), 0.5)
  mT <- sum(w1 * g1 * d1) / sum(w1 * d1); sT <- sqrt(sum(w1 * (g1 - mT)^2 * d1) / sum(w1 * d1))
  n <- 2561; g2 <- seq(ageRange[1], ageRange[2], length.out = n)
  d2 <- Dens(g2, LL(g2)); w2 <- c(1, rep(c(4, 2), (n - 3) / 2), 4, 1)
  mS <- sum(w2 * g2 * d2) / sum(w2 * d2); sS <- sqrt(sum(w2 * (g2 - mS)^2 * d2) / sum(w2 * d2))
  res[[arm]] <- data.frame(arm = arm, nTrap = length(g1), trapMean = mT, trapSd = sT,
                           nSimpson = n, simpsonMean = mS, simpsonSd = sS, diff = mS - mT)
  print(res[[arm]], digits = 10)
}
stopCluster(cl)
write.csv(do.call(rbind, res), "fine-exact.csv", row.names = FALSE)
