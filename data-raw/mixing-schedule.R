# ESS per CPU hour of TempoMCMC() under the move schedule before #19 and the
# new default (TempoMoves(): with characters, no `tree` or `joint`; local
# path moves at 0.25 x the nodes, tips and branches; one `shift` per
# iteration), on the fixed trees of the mixing study (dev/mixing/README.md:
# small and large, cell beta = gamma = 1, K = 16, z in [-3, 3], 30 / 100 /
# 300 nested characters, all six parameters free, start at the truth).
# Unlike dev/mixing/, chains run the engine's own TempoMCMC().
# Usage, from the package root with EvoTempo installed:
#   Rscript data-raw/mixing-schedule.R chain <tree> <nChar> <config> <rep> [nIter] [burnin]
#   Rscript data-raw/mixing-schedule.R summary
# Environment: EVOTEMPO_MIXING_OUT (default data-raw/mixing-results)
suppressMessages(library(EvoTempo))
args <- commandArgs(trailingOnly = TRUE)
outDir <- Sys.getenv("EVOTEMPO_MIXING_OUT", "data-raw/mixing-results")
dir.create(file.path(outDir, "chains"), recursive = TRUE, showWarnings = FALSE)
nSlice <- 16L
# Model and priors of the mixing study (rung 7, #7)
model <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5, beta = 1, gamma = 1)
meanlog <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.2), theta = log(0.2),
             s = log(0.6), rate = log(0.1))
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.5, rate = 0.4)
zTrunc <- 4
logMass <- log(stats::pnorm(zTrunc) - stats::pnorm(-zTrunc))
priors <- lapply(names(meanlog), function(p) {
  m <- meanlog[[p]]
  s <- sdlog[[p]]
  force(m)
  force(s)
  function(x) {
    if (!(x > 0) || abs(log(x) - m) > zTrunc * s) {
      return(-Inf)
    }
    stats::dlnorm(x, m, s, log = TRUE) - logMass
  }
})
names(priors) <- names(meanlog)
oldMoves <- c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1)
quantities <- c("theta", "s", "rate", "lambda", "rootBin", "clock", "logPost")

# Geyer initial monotone sequence ESS (as rungs 6-7 and dev/mixing)
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

mode <- if (length(args)) args[[1]] else stop("see header for usage")
if (mode == "chain") {
  treeId <- as.integer(args[[2]])
  nChar <- as.integer(args[[3]])
  config <- args[[4]]
  rep <- as.integer(args[[5]])
  nIter <- if (length(args) > 5) as.integer(args[[6]]) else 5000L
  burnin <- if (length(args) > 6) as.integer(args[[7]]) else 2000L
  file <- file.path(outDir, "chains", sprintf("t%d-c%03d-%s-r%d.rds", treeId, nChar,
                                              config, rep))
  if (file.exists(file)) quit(save = "no")
  tr <- readRDS("dev/mixing/results/trees.rds")[[treeId]]
  chars <- tr$chars[, seq_len(nChar), drop = FALSE]
  moves <- switch(config, old = oldMoves, new = NULL, stop("Unknown config"))
  set.seed(1900L + 100L * treeId + 10L * rep + nChar)
  t0 <- proc.time()
  run <- TempoMCMC(tr$tt, model, tr$truth, priors, nIter = nIter, burnin = burnin,
                   thin = nIter %/% 1000L, nSlice = nSlice, moves = moves,
                   characters = chars, minAcceptance = if (config == "old") 0 else 0.01)
  cpu <- unname(sum((proc.time() - t0)[c("user.self", "sys.self")]))
  out <- list(tree = treeId, nChar = nChar, config = config, rep = rep, nIter = nIter,
              burnin = burnin, cpu = cpu, samples = run$samples[quantities],
              acceptance = run$acceptance, moves = run$moves,
              ess = vapply(run$samples[quantities], Ess, 0),
              sha = Sys.getenv("EVOTEMPO_SHA"))
  saveRDS(out, file)
  cat(sprintf("%s: %.0f s CPU, min ESS %.1f\n", basename(file), cpu, min(out$ess)))
}

if (mode == "summary") {
  files <- list.files(file.path(outDir, "chains"), "\\.rds$", full.names = TRUE)
  res <- lapply(files, readRDS)
  chains <- do.call(rbind, lapply(res, function(r) {
    # CPU time scaled to the sampling phase
    hours <- r$cpu * r$nIter / (r$nIter + r$burnin) / 3600
    data.frame(tree = c("small", "large")[r$tree], nChar = r$nChar, config = r$config,
               rep = r$rep, cpuHours = hours, quantity = names(r$ess),
               ess = unname(r$ess), essPerHour = unname(r$ess) / hours,
               essPer1000 = unname(r$ess) / r$nIter * 1000,
               sd = vapply(names(r$ess), function(q) {
                 x <- r$samples[[q]]
                 stats::sd(if (q %in% c("theta", "s", "rate", "lambda", "clock")) log(x) else x)
               }, 0))
  }))
  write.csv(chains, file.path(outDir, "mixing-schedule-chains.csv"), row.names = FALSE)
  key <- interaction(chains$tree, chains$nChar, chains$quantity, drop = TRUE)
  cells <- do.call(rbind, lapply(split(chains, key), function(d) {
    o <- d[d$config == "old", ]
    n <- d[d$config == "new", ]
    data.frame(tree = d$tree[1], nChar = d$nChar[1], quantity = d$quantity[1],
               nRep = nrow(n), essPerHourOld = mean(o$essPerHour),
               essPerHourNew = mean(n$essPerHour), gain = mean(n$essPerHour) / mean(o$essPerHour),
               essPer1000Old = mean(o$essPer1000), essPer1000New = mean(n$essPer1000),
               sdOld = mean(o$sd), sdNew = mean(n$sd),
               secPerIterOld = mean(o$cpuHours) * 3600 / res[[1]]$nIter,
               secPerIterNew = mean(n$cpuHours) * 3600 / res[[1]]$nIter)
  }))
  cells <- cells[order(match(cells$tree, c("small", "large")), cells$nChar,
                       match(cells$quantity, quantities)), ]
  write.csv(cells, file.path(outDir, "mixing-schedule-gains.csv"), row.names = FALSE)
  # Split R-hat over the replicate chains of each cell
  SplitRhat <- function(xs) {
    halves <- unlist(lapply(xs, function(x) {
      n <- length(x) %/% 2
      list(x[seq_len(n)], x[n + seq_len(n)])
    }), recursive = FALSE)
    n <- length(halves[[1]])
    m <- vapply(halves, mean, 0)
    w <- mean(vapply(halves, stats::var, 0))
    b <- n * stats::var(m)
    sqrt(((n - 1) / n * w + b / n) / w)
  }
  rhat <- do.call(rbind, lapply(split(res, vapply(res, function(r) {
    paste(r$tree, r$nChar, r$config)
  }, "")), function(rs) {
    data.frame(tree = c("small", "large")[rs[[1]]$tree], nChar = rs[[1]]$nChar,
               config = rs[[1]]$config, t(vapply(quantities, function(q) {
                 xs <- lapply(rs, function(r) {
                   x <- r$samples[[q]]
                   if (q %in% c("theta", "s", "rate", "lambda", "clock")) log(x) else x
                 })
                 SplitRhat(xs)
               }, 0)))
  }))
  write.csv(rhat, file.path(outDir, "mixing-schedule-rhat.csv"), row.names = FALSE)
  print(cells[, c("tree", "nChar", "quantity", "gain", "essPerHourOld", "essPerHourNew")],
        digits = 3, row.names = FALSE)
  print(rhat, digits = 3, row.names = FALSE)
}
