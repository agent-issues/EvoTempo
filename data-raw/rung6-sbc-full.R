# Rung 6, full size (issue #21): simulation-based calibration (Talts et al.
# 2018) per (beta, gamma) cell, with time-tree uncertainty.
#
# Generalises data-raw/rung6-sbc.R (Phase 2; its results in
# data-raw/rung6-results/ are left untouched).  Differences:
# - tree moves (nodeAge, fossilAge, slide, saJump) as well as the path moves
#   and `joint` parameter-path moves, so topology, node ages, fossil ages and
#   sampled-ancestor status are all inferred;
# - fossil ages are data only through their 1-Myr stratigraphic bin
#   [floor(a), floor(a) + 1].  Given the true age the bin is deterministic, so
#   the likelihood of the range is flat within it and ranks stay uniform;
# - the chain starts from an independent prior draw of the parameters and from
#   a tree randomised by a prior-only run of the tree moves, never at the truth;
# - long chains with an ESS target: ESS (Geyer initial monotone sequence) of
#   every ranked quantity on the stored draws; a chain that misses the target
#   is extended (re-adapting, then continuing from its final state);
# - ranks on nRank = 99 equally spaced stored draws (ranks 0..99, 100 values);
# - ranked quantities: every free parameter, the root bin, the clock
#   (sum_e int_e e^x dt) and the tree statistics treeLength, meanNodeAge,
#   nAncestor, oldestFossil, youngestFossil.
# Survival is imposed by resimulating at fixed parameters (AGENTS.md,
# "Validation simulators"); rejection on tip count is a data feature and may
# redraw parameters.
#
# Replicate g (1, 2, ...) belongs to cell (g - 1) %% 4 + 1 and has seed
# 21000 + g, so any subset can be run in any order, and reruns skip done work.
#
# Usage, from a checkout with EvoTempo installed (or EVOTEMPO_LOAD_ALL=1):
#   Rscript data-raw/rung6-sbc-full.R task <taskId> <nPerTask>
#     runs replicates (taskId - 1) * nPerTask + 1 .. taskId * nPerTask
#   Rscript data-raw/rung6-sbc-full.R reps <g1,g2,...>
# Environment:
#   EVOTEMPO_SBC_OUT   output folder (default data-raw/rung6-full-results)
#   EVOTEMPO_GRID      "default" for K = 26, z in [-5, 5]; else K = 16, z in [-3, 3]
#   EVOTEMPO_SBC_ITER  base post-burn-in iterations, one value or one per cell
#   EVOTEMPO_SBC_MAXEXT  maximum number of chain extensions (default 2)
#   EVOTEMPO_SBC_MOVES   overrides of the move counts, e.g. "nodeAge=12,slide=8"
#   EVOTEMPO_LOAD_ALL  if set, pkgload::load_all(".") instead of library()
if (nzchar(Sys.getenv("EVOTEMPO_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
  for (fn in c(".Setup", ".Sample1", ".TreeStats")) {
    assign(fn, get(fn, asNamespace("EvoTempo")))
  }
}
args <- commandArgs(trailingOnly = TRUE)
EnvNum <- function(name, default) {
  v <- Sys.getenv(name)
  if (nzchar(v)) as.numeric(strsplit(v, ",")[[1]]) else default
}
outDir <- Sys.getenv("EVOTEMPO_SBC_OUT", "data-raw/rung6-full-results")
dir.create(file.path(outDir, "replicates"), recursive = TRUE, showWarnings = FALSE)

defaultGrid <- identical(Sys.getenv("EVOTEMPO_GRID"), "default")
model0 <- if (defaultGrid) {
  list(nBin = 26, zRange = c(-5, 5))
} else {
  list(nBin = 16, zRange = c(-3, 3))
}
cells <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
# Base chain length per cell, from the timing pilot (README): theta mixes
# slowly once characters depend on tempo (gamma = 1)
nIterCell <- rep_len(EnvNum("EVOTEMPO_SBC_ITER", c(10000, 20000, 10000, 20000)),
                     nrow(cells))
maxExtend <- EnvNum("EVOTEMPO_SBC_MAXEXT", 2)
burnin <- 1000L        # first segment: adaptation and burn-in
reburn <- 200L         # each extension re-adapts proposal scales first
nStore <- 1000L        # stored draws per segment
nRank <- 99L           # draws ranked (equally spaced over the stored draws)
essTarget <- 50
nPre <- 300L           # prior-only tree-move iterations to randomise the start
nSlice <- 16L
nChar <- 30L
crown <- 8
tipRange <- c(20, 40)  # tips excluding sampled ancestors
binWidth <- 1          # stratigraphic resolution of fossil ages (Myr)
mainMoves <- c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1,
               nodeAge = 12, fossilAge = 4, slide = 8, saJump = 3)
preMoves <- c(tree = 1, nodeAge = 4, fossilAge = 2, slide = 4, saJump = 2)
# EVOTEMPO_SBC_MOVES, e.g. "nodeAge=12,slide=8", overrides entries of mainMoves
for (kv in strsplit(Sys.getenv("EVOTEMPO_SBC_MOVES"), ",")[[1]]) {
  kv <- strsplit(kv, "=")[[1]]
  mainMoves[[kv[1]]] <- as.numeric(kv[2])
}

meanlog <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.05), theta = log(0.2),
             s = log(0.6), rate = log(0.1))
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.25, rate = 0.4)
priors <- TempoPriors(meanlog, sdlog)
PriorDraw <- function() {
  as.list(setNames(stats::rlnorm(6, meanlog, sdlog), names(meanlog)))
}
treeQ <- c("treeLength", "meanNodeAge", "nAncestor", "oldestFossil",
           "youngestFossil")
quantities <- c(names(meanlog), "rootBin", "clock", treeQ)

# Join two daughter histories under a zero-length root edge
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

# Effective sample size of one chain: Geyer's initial monotone sequence
# estimator (as in Stan), capped at the number of draws
Ess <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4 || stats::var(x) == 0) {
    return(NA_real_)
  }
  ac <- stats::acf(x, lag.max = n - 1, plot = FALSE, demean = TRUE)$acf[, 1, 1]
  nPair <- length(ac) %/% 2
  pairSum <- ac[2 * seq_len(nPair) - 1] + ac[2 * seq_len(nPair)]
  stop1 <- which(pairSum <= 0)[1]
  if (!is.na(stop1)) {
    pairSum <- pairSum[seq_len(stop1 - 1)]
  }
  tau <- if (length(pairSum)) -1 + 2 * sum(cummin(pairSum)) else 1
  min(n, n / max(tau, 1e-8))
}

RankOf <- function(draws, value) {
  if (!is.finite(value) || anyNA(draws)) {
    return(NA_integer_)
  }
  sum(draws < value) + sample.int(sum(draws == value) + 1L, 1L) - 1L
}

Simulate <- function(cell, model) {
  tries <- 0
  repeat {
    tries <- tries + 1
    truth <- PriorDraw()
    set <- .Setup(truth, model)
    # Conditioned on survival: resimulate at the same parameters until both
    # crown daughters have sampled extant descendants
    repeat {
      k <- .Sample1(set$prior)
      h <- tryCatch(lapply(1:2, function(i) {
        SimulateDCET(crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                     start = k, rho = model$rho, maxEdges = 5000)
      }), error = function(e) NULL)
      if (is.null(h) ||
          all(vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0) > 0)) {
        break
      }
    }
    # More than 5000 branches: survives but exceeds the tip range
    if (is.null(h)) next
    history <- Join(h[[1]], h[[2]], set$grid$x[k])
    tree <- ReconstructedTree(history, gamma = cell$gamma, rate = truth$rate)
    nTip <- sum(attr(tree, "tipType") != "ancestor")
    if (nTip >= tipRange[1] && nTip <= tipRange[2]) break
  }
  list(truth = truth, k = k, history = history, tree = tree, tries = tries)
}

Replicate <- function(g) {
  file <- file.path(outDir, "replicates", sprintf("rep%04d.rds", g))
  if (file.exists(file)) {
    return(invisible(NULL))
  }
  t0 <- Sys.time()
  cellId <- (g - 1L) %% nrow(cells) + 1L
  cell <- cells[cellId, ]
  model <- DCETModel(nBin = model0$nBin, zRange = model0$zRange, rho = 0.5,
                     beta = cell$beta, gamma = cell$gamma)
  set.seed(21000L + g)
  sim <- Simulate(cell, model)
  tree <- sim$tree
  chars <- SimulateCharacters(tree, nChar, 2)
  clk <- ReconstructedTree(sim$history, gamma = 1, rate = 1)
  inv <- ReconstructedTree(sim$history, gamma = -1, rate = 1)
  tree$root.edge <- NULL
  tt0 <- TempoTree(tree, present = crown)
  fossil <- tt0$tipType != "extant"
  tipAge <- tt0$age[seq_len(tt0$nTip)]
  lo <- ifelse(fossil, pmin(floor(tipAge / binWidth) * binWidth, crown - binWidth), 0)
  tt <- TempoTree(tree, present = crown,
                  ranges = cbind(lo, ifelse(fossil, lo + binWidth, 0)))
  truthStats <- .TreeStats(tt, attr(clk, "characterLength"),
                           attr(inv, "characterLength"))
  truth <- c(unlist(sim$truth), rootBin = sim$k, truthStats)

  # Independent start: prior draw of parameters; tree randomised by prior-only
  # tree moves (no characters, parameters fixed) from the true tree
  start <- PriorDraw()
  pre <- TempoMCMC(tt, model, start, priors, nIter = nPre, nSlice = nSlice,
                   moves = preMoves, updateParameters = FALSE)
  startTree <- pre$state$tree
  nIter <- nIterCell[[cellId]]
  Chain <- function(tr, par, burn) {
    TempoMCMC(tr, model, par, priors, nIter = nIter, burnin = burn,
              thin = nIter %/% nStore, nSlice = nSlice, moves = mainMoves,
              characters = chars)
  }
  run <- Chain(startTree, start, burnin)
  samples <- run$samples
  acceptance <- list(run$acceptance)
  ess <- vapply(samples[quantities], Ess, 0)
  nExtend <- 0
  while (any(ess < essTarget, na.rm = TRUE) && nExtend < maxExtend) {
    nExtend <- nExtend + 1
    run <- Chain(run$state$tree, run$state$par, reburn)
    samples <- rbind(samples, run$samples)
    acceptance[[length(acceptance) + 1L]] <- run$acceptance
    ess <- vapply(samples[quantities], Ess, 0)
  }
  keep <- unique(round(seq(nrow(samples) / nRank, nrow(samples), length.out = nRank)))
  ranked <- samples[keep, quantities]
  ranks <- vapply(quantities, function(q) RankOf(ranked[[q]], truth[[q]]), 0L)
  out <- list(g = g, cell = cellId, beta = cell$beta, gamma = cell$gamma,
              seed = 21000L + g, truth = truth, ranks = ranks, nRank = length(keep),
              ess = ess, essRanked = vapply(ranked, Ess, 0),
              nIter = nIter * (1 + nExtend), nExtend = nExtend,
              nTip = sum(tt$tipType != "ancestor"), nFossil = sum(fossil),
              nAncestor = sum(tt$tipType == "ancestor"), tries = sim$tries,
              acceptance = do.call(rbind, acceptance),
              posteriorMean = colMeans(samples[quantities], na.rm = TRUE),
              posteriorSd = vapply(samples[quantities], stats::sd, 0, na.rm = TRUE),
              grid = model0, moves = mainMoves, nStore = nStore,
              runtime = as.numeric(difftime(Sys.time(), t0, units = "mins")),
              packageVersion = as.character(utils::packageVersion("EvoTempo")),
              sha = Sys.getenv("EVOTEMPO_SHA"))
  tmp <- paste0(file, ".tmp")
  saveRDS(out, tmp)
  file.rename(tmp, file)
  invisible(out)
}

if (length(args) && args[[1]] == "task") {
  id <- as.integer(args[[2]])
  n <- as.integer(args[[3]])
  reps <- (id - 1L) * n + seq_len(n)
} else if (length(args) && args[[1]] == "reps") {
  reps <- as.integer(strsplit(args[[2]], ",")[[1]])
} else {
  stop("Usage: rung6-sbc-full.R task <id> <nPerTask> | reps <g1,g2,...>")
}
for (g in reps) {
  res <- Replicate(g)
  if (!is.null(res)) {
    cat(sprintf("rep %d cell %d: %.1f min, nIter %d, min ESS %.0f (%s)\n", g, res$cell,
                res$runtime, res$nIter, min(res$ess, na.rm = TRUE),
                names(which.min(res$ess))))
  }
}
