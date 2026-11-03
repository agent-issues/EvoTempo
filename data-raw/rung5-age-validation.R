# Rung 5 with a free start age (#18): prior-only TempoMCMC() with every tree
# move against rejection sampling from the forward simulator, as
# data-raw/rung5-validation.R and data-raw/rung5-fixed-parameters.R, for
#   root        start = "root", crown age free with a uniform prior on
#               [crownLo, crownHi] (move `rootAge`);
#   origin      start = "origin", origin age fixed at originAge (the root
#               age moves with `nodeAge`);
#   originFree  start = "origin", origin age free with a uniform prior on
#               [originLo, originHi] (move `originAge`).
# A second argument "fixed" holds the parameters at the prior medians, which
# isolates the tree and age moves; "free" also samples the parameters.
# Simulator: draw the start age from its prior (and the parameters from theirs
# when free); resimulate at those values until the conditioning holds (both
# crown daughters, or the stem lineage, with a sampled extant descendant);
# accept if there are exactly nExtant extant tips and nFossil fossil samples
# and, with an origin, the root is not a sampled ancestor (TempoMCMC() does
# not support one).  Rejection on these data features may redraw the start age
# and the parameters (AGENTS.md, "Validation simulators").
# Usage, from the package root (EvoTempo installed, or EVOTEMPO_LOAD_ALL=1):
#   Rscript data-raw/rung5-age-validation.R <case> <fixed|free> [nIter] [nChain] [nSim]
# Environment: EVOTEMPO_RUNG5_OUT (default data-raw/rung5-age-results),
#   EVOTEMPO_CORES (default 16), EVOTEMPO_SEED (offset added to every seed, for
#   independent replicates; default 0, which also leaves the file names
#   unchanged)
if (nzchar(Sys.getenv("EVOTEMPO_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
  for (fn in c(".Setup", ".Sample1", ".TreeStats")) {
    assign(fn, get(fn, asNamespace("EvoTempo")))
  }
}
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]
stopifnot(case %in% c("root", "origin", "originFree"))
fixedPar <- length(args) > 1 && args[[2]] == "fixed"
nIter <- if (length(args) > 2) as.integer(args[[3]]) else 40000L
nChain <- if (length(args) > 3) as.integer(args[[4]]) else 16L
nSim <- if (length(args) > 4) as.integer(args[[5]]) else 4000L
nCore <- as.integer(Sys.getenv("EVOTEMPO_CORES", "16"))
outDir <- Sys.getenv("EVOTEMPO_RUNG5_OUT", "data-raw/rung5-age-results")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
seedOffset <- as.integer(Sys.getenv("EVOTEMPO_SEED", "0"))
tag <- paste0(case, if (fixedPar) "-fixed" else "-free",
              if (seedOffset != 0) paste0("-seed", seedOffset))
burnin <- 1000L
thin <- 10L
nSlice <- 16L

origin <- case != "root"
model <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5,
                   start = if (origin) "origin" else "root")
nExtant <- 3
nFossil <- 2
crownLo <- 3; crownHi <- 5
originAge <- 5
originLo <- 4.5; originHi <- 6
# Priors of rung 5 (narrow, to keep simulated clades small)
meanlog <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.15), theta = log(0.2),
             s = log(0.6))
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.25)
priors <- TempoPriors(meanlog, sdlog)
medians <- as.list(exp(meanlog))
PriorDraw <- function() {
  if (fixedPar) medians else as.list(setNames(stats::rlnorm(5, meanlog, sdlog), names(meanlog)))
}
StartDraw <- switch(case,
                    root = function() stats::runif(1, crownLo, crownHi),
                    origin = function() originAge,
                    originFree = function() stats::runif(1, originLo, originHi))
agePriors <- switch(case,
                    root = list(root = function(a) stats::dunif(a, crownLo, crownHi, log = TRUE)),
                    origin = list(),
                    originFree = list(origin = function(a) {
                      stats::dunif(a, originLo, originHi, log = TRUE)
                    }))
moves <- c(tree = 1, node = 1, tip = 1, branch = 1, nodeAge = 2, fossilAge = 2,
           slide = 2, saJump = 2,
           switch(case, root = c(rootAge = 2), origin = NULL, originFree = c(originAge = 2)))

# Join two daughter histories under a zero-length root edge (crown start)
Join <- function(h1, h2, x0) {
  n1 <- nrow(h1$edges)
  e1 <- h1$edges
  e2 <- h2$edges
  e1$parent <- ifelse(is.na(e1$parent), 1L, e1$parent + 1L)
  e2$parent <- ifelse(is.na(e2$parent), 1L, e2$parent + 1L + n1)
  f1 <- h1$fossils
  f2 <- h2$fossils
  f1$edge <- f1$edge + 1L
  f2$edge <- f2$edge + 1L + n1
  structure(list(edges = rbind(data.frame(parent = NA_integer_, start = 0, end = 0,
                                          type = "speciation", sampled = NA), e1, e2),
                 paths = c(list(list(t = 0, x = x0)), h1$paths, h2$paths),
                 fossils = rbind(f1, f2), time = h1$time),
            class = "TempoHistory", kind = "discrete")
}
NExtant <- function(h) sum(h$edges$sampled %in% TRUE)

# One accepted draw: list(tree, stats)
SimOne <- function() {
  repeat {
    par <- PriorDraw()
    T0 <- StartDraw()
    set <- .Setup(par, model)
    nLineage <- if (origin) 1 else 2
    repeat {
      k <- .Sample1(set$prior)
      h <- tryCatch(lapply(seq_len(nLineage), function(i) {
        SimulateDCET(T0, set$grid$x, set$Q, set$lambda, set$mu, set$psi, start = k,
                     rho = model$rho, maxEdges = 5000)
      }), error = function(e) NULL)
      if (is.null(h)) break
      ext <- vapply(h, NExtant, 0)
      # Conditioning: each lineage from the start has a sampled extant descendant
      if (all(ext > 0)) break
    }
    if (is.null(h)) next
    fos <- vapply(h, function(hi) nrow(hi$fossils), 0)
    if (sum(ext) != nExtant || sum(fos) != nFossil) next
    history <- if (origin) h[[1]] else Join(h[[1]], h[[2]], set$grid$x[k])
    tree <- ReconstructedTree(history, gamma = 1)
    inv <- ReconstructedTree(history, gamma = -1)
    if (!origin) tree$root.edge <- NULL
    tt <- TempoTree(tree, present = T0)
    if (origin && tt$sampled[tt$root - tt$nTip]) next
    stats <- c(unlist(par), rootBin = k,
               .TreeStats(tt, attr(tree, "characterLength"), attr(inv, "characterLength")),
               if (origin) c(originAge = T0))
    return(list(tree = tree, T0 = T0, par = par, stats = stats))
  }
}

cl <- makeCluster(nCore)
invisible(clusterCall(cl, function() {
  if (nzchar(Sys.getenv("EVOTEMPO_LOAD_ALL"))) {
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
  } else {
    suppressMessages(library(EvoTempo))
  }
  NULL
}))
clusterExport(cl, c(".Setup", ".Sample1", ".TreeStats", "SimOne", "PriorDraw", "StartDraw",
                    "Join", "NExtant", "model", "origin", "fixedPar", "medians", "meanlog",
                    "sdlog", "nExtant", "nFossil", "crownLo", "crownHi", "originAge",
                    "originLo", "originHi", "case"))
t0 <- Sys.time()
perWorker <- ceiling(nSim / nCore)
clusterExport(cl, c("perWorker", "seedOffset"))
sims <- parLapply(cl, seq_len(nCore), function(i) {
  set.seed(18500 + i + seedOffset)
  do.call(rbind, lapply(seq_len(perWorker), function(j) SimOne()$stats))
})
simTime <- difftime(Sys.time(), t0, units = "mins")
sims <- as.data.frame(do.call(rbind, sims))

# Starting tree: one simulated tree, fossil ranges wide enough never to bind
set.seed(18000 + seedOffset)
start <- SimOne()
tt0 <- TempoTree(start$tree, present = start$T0)
tt <- TempoTree(start$tree, present = start$T0,
                ranges = cbind(0, ifelse(tt0$tipType == "extant", 0, 1e3)))
clusterExport(cl, c("tt", "priors", "agePriors", "moves", "nIter", "burnin", "thin",
                    "nSlice", "start", "seedOffset"))
t0 <- Sys.time()
chains <- parLapply(cl, seq_len(nChain), function(i) {
  set.seed(18700 + i + seedOffset)
  run <- TempoMCMC(tt, model, if (fixedPar) medians else start$par,
                   if (fixedPar) list() else priors, nIter = nIter, burnin = burnin,
                   thin = thin, nSlice = nSlice, moves = moves, agePriors = agePriors)
  list(samples = run$samples, acceptance = run$acceptance)
})
chainTime <- difftime(Sys.time(), t0, units = "mins")
stopCluster(cl)

# z-scores of MCMC - simulation for the mean and five quantiles of each
# statistic; the MCMC error is the between-chain spread of chain estimates
quantities <- intersect(c(names(meanlog), "treeLength", "meanNodeAge", "rootAge",
                          "originAge", "oldestFossil", "youngestFossil", "clock",
                          "inverseClock", "nAncestor", "cherries"), names(sims))
if (fixedPar) quantities <- setdiff(quantities, names(meanlog))
# A start age that is held fixed is not tested
quantities <- quantities[vapply(quantities, function(q) stats::var(sims[[q]]) > 0, TRUE)]
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
rows <- list()
for (q in quantities) {
  s <- sims[[q]]
  s <- s[is.finite(s)]
  per <- lapply(chains, function(ch) ch$samples[[q]][is.finite(ch$samples[[q]])])
  # Mean
  cm <- vapply(per, mean, 0)
  rows[[length(rows) + 1L]] <- data.frame(
    quantity = q, stat = "mean", sim = mean(s), mcmc = mean(cm),
    z = (mean(cm) - mean(s)) / sqrt(stats::var(cm) / length(cm) + stats::var(s) / length(s)))
  if (q %in% c("nAncestor", "cherries")) {
    for (v in sort(unique(s))) {
      ps <- mean(s == v)
      pc <- vapply(per, function(x) mean(x == v), 0)
      rows[[length(rows) + 1L]] <- data.frame(
        quantity = q, stat = paste0("P(=", v, ")"), sim = ps, mcmc = mean(pc),
        z = (mean(pc) - ps) / sqrt(stats::var(pc) / length(pc) + ps * (1 - ps) / length(s)))
    }
    next
  }
  for (pr in probs) {
    qs <- stats::quantile(s, pr, names = FALSE)
    # Proportion of draws below the simulated quantile
    pc <- vapply(per, function(x) mean(x <= qs), 0)
    rows[[length(rows) + 1L]] <- data.frame(
      quantity = q, stat = sprintf("q%02d", round(100 * pr)), sim = pr, mcmc = mean(pc),
      z = (mean(pc) - pr) / sqrt(stats::var(pc) / length(pc) + pr * (1 - pr) / length(s)))
  }
}
out <- do.call(rbind, rows)
acceptance <- colMeans(do.call(rbind, lapply(chains, `[[`, "acceptance")))
saveRDS(list(case = case, fixed = fixedPar, model = model, moves = moves, sims = sims,
             chains = chains, out = out, acceptance = acceptance, nIter = nIter,
             burnin = burnin, thin = thin, simTime = simTime, chainTime = chainTime,
             sha = Sys.getenv("EVOTEMPO_SHA")),
        file.path(outDir, paste0("rung5-age-", tag, ".rds")))
write.csv(out, file.path(outDir, paste0("rung5-age-", tag, ".csv")), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
cat(sprintf("%s: max |z| %.2f over %d tests; sims %s, chains %s\n", tag, max(abs(out$z)),
            nrow(out), format(simTime), format(chainTime)))
print(round(acceptance[acceptance > 0], 3))
