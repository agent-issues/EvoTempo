# Start-age moves (#18) on a fixed tree at fixed parameters: the posterior of
# the root age (start = "root", move `rootAge`) or of the origin age
# (start = "origin", move `originAge`) from TempoMCMC() chains with path
# moves, against the exact one-dimensional posterior, prior x TreeLikelihood()
# on a grid of start ages (trapezoid CDF), as in data-raw/rung4b-one-parameter.R.
# Tree: ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7), root age 4, stem 0.6;
# a, c fossils; f a sampled ancestor (the fixture of tests/testthat/test-Origin.R).
# z-scores use the spread of the CDF at each exact quantile over independent
# chains.
#
# Usage: Rscript data-raw/rung5-age-oracle.R <root|origin> [nIter] [nChain] [variant] [design]
#   variant: the conditioning of the chains
#     extant          extant conditioning, marginal over the start bin (default)
#     joint           extant, joint conditioning (given the start bin)
#     survival        condition = "survival", marginal
#     survival-joint  condition = "survival", joint (#69)
#     none            condition = "none" (#69)
#     design          no chains: the exact posteriors under every conditioning
#                     and the power to tell them apart at nIter x nChain (#69)
#   design: default  the rung-5 fixture parameters (lambda 0.5, mu 0.2,
#                    psi 0.15, theta 0.2, s 0.6, rho 0.5; sections 1 and 4 of
#                    data-raw/rung5-age-results/README.md).  At these the
#                    conditionings differ by less than 2 SE (#69).
#           power    (#69) lambda = mu = 1, psi 0.15, theta 0.4, s 1, rho 0.35:
#                    critical in every bin, so the survival probability falls
#                    with the start age, and spread over the bins (s = 1 on
#                    z in [-3, 3]), so joint and marginal conditioning differ;
#                    finer exact grid (641 points).
#   Results of a variant other than extant go to files with a "-<variant>"
#   suffix (#54), and those of the power design to "-power-<variant>" (the
#   extant variant included).
# Every run also scores its chains against the exact posterior of every other
# conditioning (`-cross.csv`).  This is the negative control: the chains must
# reject each conditioning that differs from their own in one flag
# (joint/marginal, extant/survival, extant/none), with max |z| >= 5 over the
# six statistics; at the power design the exact differences are large enough
# for that (variant `design`).
# Environment: EVOTEMPO_RUNG5_OUT (default data-raw/rung5-age-results);
# EVOTEMPO_RUNG5_SEED_OFFSET (see seedBase);
# EVOTEMPO_RUNG5_ESS_RATE, the ESS per chain iteration assumed by `design`
# (default 0.078, half the rate of the 32 x 20000 runs of section 4).
suppressMessages(library(EvoTempo))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]
stopifnot(case %in% c("root", "origin"))
nIter <- if (length(args) > 1) as.integer(args[[2]]) else 20000L
nChain <- if (length(args) > 2) as.integer(args[[3]]) else 16L
variant <- if (length(args) > 3) args[[4]] else "extant"
design <- if (length(args) > 4) args[[5]] else "default"
conditionings <- c("extant", "joint", "survival", "survival-joint", "none")
stopifnot(variant %in% c(conditionings, "design"), design %in% c("default", "power"))
power <- design == "power"
suffix <- paste0(if (power) "-power", if (variant != "extant" || power) paste0("-", variant))
# EVOTEMPO_RUNG5_SEED_OFFSET (default 0): independent replicates, written to
# files with a "-seed<offset>" suffix
seedOffset <- as.integer(Sys.getenv("EVOTEMPO_RUNG5_SEED_OFFSET", "0"))
if (seedOffset != 0) suffix <- paste0(suffix, "-seed", seedOffset)
seedBase <- c(extant = 18900, joint = 19900, survival = 20900, `survival-joint` = 21900,
              none = 22900, design = 0)[[variant]] + (if (power) 100000 else 0) + seedOffset
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
rho <- if (power) 0.35 else 0.5
Model <- function(v) {
  DCETModel(nBin = 16, zRange = c(-3, 3), rho = rho,
            start = if (origin) "origin" else "root",
            condition = switch(v, survival = , `survival-joint` = "survival",
                               none = "none", "extant"),
            joint = v %in% c("joint", "survival-joint"))
}
models <- setNames(lapply(conditionings, Model), conditionings)
model <- if (variant == "design") NULL else models[[variant]]
par <- if (power) {
  list(lambda = 1, mu = 1, psi = 0.15, theta = 0.4, s = 1)
} else {
  list(lambda = 0.5, mu = 0.2, psi = 0.15, theta = 0.2, s = 0.6)
}
ageRange <- if (origin) c(4.2, 8) else c(3.4, 7)
prior <- function(a) stats::dunif(a, ageRange[1], ageRange[2], log = TRUE)
agePriors <- if (origin) list(origin = prior) else list(root = prior)
root <- tree$root
WithAge <- function(a) {
  tr <- tree
  if (origin) tr$stem <- a - tr$age[[root]] else tr$age[root] <- a
  tr
}
grid <- seq(ageRange[1], ageRange[2], length.out = if (power) 641 else 321)
cl <- makeCluster(16)
invisible(clusterCall(cl, function() suppressMessages(library(EvoTempo))))
clusterExport(cl, c("tree", "model", "models", "par", "WithAge", "origin", "root", "agePriors",
                    "prior", "ageRange", "nIter", "seedBase"))
# Exact log-likelihood on the grid under every conditioning
logLAll <- sapply(conditionings, function(v) {
  m <- models[[v]]
  unlist(parLapply(cl, grid, function(a, m) TreeLikelihood(WithAge(a), par, m), m))
})
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
wTrap <- c(0.5, rep(1, length(grid) - 2), 0.5)
Exact <- function(logL) {
  lp <- logL + prior(grid)
  dens <- exp(lp - max(lp))
  cdf <- c(0, cumsum((dens[-1] + dens[-length(dens)]) / 2 * diff(grid)))
  cdf <- cdf / max(cdf)
  # Trapezoid rule, as for the CDF (an equal-weight sum overweights the ends
  # of the grid; #57)
  m <- sum(wTrap * grid * dens) / sum(wTrap * dens)
  list(cdf = cdf, q = stats::approx(cdf, grid, probs)$y, mean = m,
       sd = sqrt(sum(wTrap * (grid - m)^2 * dens) / sum(wTrap * dens)))
}
exact <- setNames(lapply(conditionings, function(v) Exact(logLAll[, v])), conditionings)
# Conditionings that differ from each in one flag: the negative controls
oneFlag <- list(extant = c("joint", "survival", "none"),
                joint = c("extant", "survival-joint"),
                survival = c("survival-joint", "extant"),
                `survival-joint` = c("survival", "joint"),
                none = "extant")

if (variant == "design") {
  essRate <- as.numeric(Sys.getenv("EVOTEMPO_RUNG5_ESS_RATE", "0.078"))
  ess <- essRate * nIter * nChain
  # Chains targeting A, scored against B's exact posterior: expected |z| =
  # shift / SE with SE = sqrt(p (1 - p) / ESS) for the CDF at B's quantiles,
  # sd_A / sqrt(ESS) for the mean
  pw <- do.call(rbind, lapply(conditionings, function(a) {
    do.call(rbind, lapply(setdiff(conditionings, a), function(b) {
      fA <- stats::approx(grid, exact[[a]]$cdf, exact[[b]]$q)$y
      zq <- abs(fA - probs) / sqrt(probs * (1 - probs) / ess)
      zm <- abs(exact[[a]]$mean - exact[[b]]$mean) / (exact[[a]]$sd / sqrt(ess))
      data.frame(case = case, design = design, target = a, against = b,
                 oneFlag = b %in% oneFlag[[a]], ess = ess,
                 zq05 = zq[1], zq25 = zq[2], zq50 = zq[3], zq75 = zq[4], zq95 = zq[5],
                 zMean = zm, zMax = max(zq, zm))
    }))
  }))
  ex <- do.call(rbind, lapply(conditionings, function(v) {
    data.frame(case = case, design = design, conditioning = v,
               p = c(probs, NA), exact = c(exact[[v]]$q, exact[[v]]$mean))
  }))
  write.csv(ex, file.path(outDir, paste0("rung5-age-oracle-", case, suffix, "-exact.csv")),
            row.names = FALSE)
  write.csv(pw, file.path(outDir, paste0("rung5-age-oracle-", case, suffix, ".csv")),
            row.names = FALSE)
  print(ex, digits = 6, row.names = FALSE)
  print(pw[, c("target", "against", "oneFlag", "zMax")], digits = 3, row.names = FALSE)
  cat(sprintf("Assumed ESS %.0f; smallest expected max|z| over one-flag pairs: %.2f\n", ess,
              min(pw$zMax[pw$oneFlag])))
  stopCluster(cl)
  quit(save = "no")
}

t0 <- Sys.time()
chains <- parLapply(cl, seq_len(nChain), function(i) {
  set.seed(seedBase + i)
  run <- TempoMCMC(tree, model, par, list(), nIter = nIter, burnin = 1000, thin = 5,
                   nSlice = 16, agePriors = agePriors,
                   moves = c(tree = 1, node = 2, tip = 2, branch = 2,
                             if (origin) c(originAge = 2) else c(rootAge = 2)))
  list(age = run$samples[[if (origin) "originAge" else "rootAge"]],
       acceptance = run$acceptance)
})
runtime <- difftime(Sys.time(), t0, units = "mins")
stopCluster(cl)
# Chains scored against the exact posterior under conditioning `v`
Score <- function(v) {
  e <- exact[[v]]
  rows <- lapply(seq_along(probs), function(k) {
    f <- vapply(chains, function(ch) mean(ch$age <= e$q[k]), 0)
    data.frame(case = case, p = probs[k], exact = e$q[k],
               mcmc = stats::quantile(unlist(lapply(chains, `[[`, "age")), probs[k], names = FALSE),
               cdfAtExact = mean(f), z = (mean(f) - probs[k]) / (stats::sd(f) / sqrt(length(f))))
  })
  cm <- vapply(chains, function(ch) mean(ch$age), 0)
  rbind(do.call(rbind, rows),
        data.frame(case = case, p = NA, exact = e$mean, mcmc = mean(cm), cdfAtExact = NA,
                   z = (mean(cm) - e$mean) / (stats::sd(cm) / sqrt(length(cm)))))
}
out <- Score(variant)
logL <- logLAll[, variant]
cross <- do.call(rbind, lapply(conditionings, function(v) {
  cbind(target = variant, against = v, oneFlag = v %in% oneFlag[[variant]], Score(v))
}))
acceptance <- colMeans(do.call(rbind, lapply(chains, `[[`, "acceptance")))
saveRDS(list(case = case, variant = variant, design = design, par = par, rho = rho, grid = grid,
             logL = logL, logLAll = logLAll, chains = chains, out = out, cross = cross,
             acceptance = acceptance, nIter = nIter, runtime = runtime,
             sha = Sys.getenv("EVOTEMPO_SHA")),
        file.path(outDir, paste0("rung5-age-oracle-", case, suffix, ".rds")))
write.csv(out, file.path(outDir, paste0("rung5-age-oracle-", case, suffix, ".csv")),
          row.names = FALSE)
write.csv(cross, file.path(outDir, paste0("rung5-age-oracle-", case, suffix, "-cross.csv")),
          row.names = FALSE)
print(out, digits = 4, row.names = FALSE)
crossMax <- tapply(abs(cross$z), cross$against, max)
print(round(crossMax, 2))
cat(sprintf("own max|z| = %.2f; one-flag alternatives min max|z| = %.2f (%s)\n",
            max(abs(out$z)), min(crossMax[oneFlag[[variant]]]),
            paste(oneFlag[[variant]], collapse = ", ")))
print(round(acceptance[acceptance > 0], 3))
cat("runtime", format(runtime), "\n")
