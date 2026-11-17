# Trilobite CET pipeline: forward-simulation check of the all-fossil set-up
# before it is relied on (README, "Conditioning").
#
#   Rscript data-raw/empirical/trilobites/cet/check-window.R [nIter] [nChain] [nSim]
#
# The trilobite analysis needs, together, four things that no earlier
# validation run combined: no extant sample (rho = 0, every tip a fossil, the
# youngest tip older than the end of the process at age 0), an origin start
# with a free origin age, condition = "survival" (the stem lineage leaves at
# least one sample), and tempo-dependent diversification (beta = 1, rates
# that vary across bins).  Rung 5 (data-raw/rung5-age-results/README.md)
# checked survival conditioning with bin-varying rates only against the
# package's own TreeLikelihood(), never against simulation, and always with
# rho > 0 and extant tips.  check-mkv.R checks the beta = 0 case exactly
# against RevBayes; this script checks beta = 1 against simulation.
#
# Design (as data-raw/rung5-age-validation.R): K = 16, z in [-3, 3],
# beta = 1, delta = 0, rho = 0, start = "origin", marginal conditioning,
# parameters fixed at lambda 0.6, mu 0.5, psi 0.15, theta 0.2, s 0.6, origin
# age ~ U(3, 10), exactly 4 fossil samples (tips or sampled ancestors).
#   Simulator, "survival": draw the origin age from its prior; draw a start bin
#   and simulate forward to age 0 with rho = 0, repeating at the same origin
#   age (new start bin each time) until the stem leaves at least one fossil;
#   accept if it left exactly 4 and the root is not a sampled ancestor.
#   Simulator, "none": the same without the repeat (one forward run per
#   origin draw).
#   MCMC: prior-only TempoMCMC() with every tree move (nodeAge, fossilAge,
#   slide, saJump, originAge) and path moves, fossil ranges [0, 1000] (never
#   binding, so fossil ages are free with a flat prior, as the simulated ones
#   are), once with condition = "survival" and once with condition = "none".
# Tests (rung 5): for treeLength, meanNodeAge, rootAge, originAge,
# oldestFossil, youngestFossil, clock and inverseClock, the mean and the CDF
# at five simulated quantiles; for nAncestor and cherries, each value's
# probability; z uses the spread over independent chains and the
# simulation's sampling error.
# Criteria, fixed before any run:
#   (a) each MCMC run against the simulator of its own conditioning: every
#       |z| < 3.5 (z is close to t with nChain - 1 degrees of freedom);
#   (b) power: each MCMC run against the OTHER conditioning's simulator must
#       give max |z| >= 5; otherwise the check is reported as lacking power
#       for the conditioning (it still checks the all-fossil tree density).
# Environment: CET_WINDOW_OUT (default data-raw/empirical/trilobites/cet/window-check),
#   EVOTEMPO_CORES (default 16), CET_WINDOW_SEED (offset added to every seed,
#   for an independent replicate; default 0, which leaves the file names
#   unchanged; otherwise they end in -seed<offset>).
if (nzchar(Sys.getenv("EVOTEMPO_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
library(parallel)
ns <- asNamespace("EvoTempo")
args <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args) > 0) as.integer(args[[1]]) else 40000L
nChain <- if (length(args) > 1) as.integer(args[[2]]) else 32L
nSim <- if (length(args) > 2) as.integer(args[[3]]) else 8000L
nCore <- as.integer(Sys.getenv("EVOTEMPO_CORES", "16"))
outDir <- Sys.getenv("CET_WINDOW_OUT", "data-raw/empirical/trilobites/cet/window-check")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
seedOffset <- as.integer(Sys.getenv("CET_WINDOW_SEED", "0"))
tag <- if (seedOffset != 0) paste0("-seed", seedOffset) else ""
burnin <- 1000L
thin <- 10L
nSlice <- 16L
nFossil <- 4L
originLo <- 3
originHi <- 10
par <- list(lambda = 0.6, mu = 0.5, psi = 0.15, theta = 0.2, s = 0.6)
Model <- function(condition) {
  DCETModel(nBin = 16, zRange = c(-3, 3), beta = 1, rho = 0, condition = condition,
            start = "origin")
}
agePriors <- list(origin = function(a) stats::dunif(a, originLo, originHi, log = TRUE))
moves <- c(tree = 1, node = 1, tip = 1, branch = 1, nodeAge = 2, fossilAge = 2, slide = 2,
           saJump = 2, originAge = 2)

# Design check: survival of the stem over the origin prior's range
set <- get(".Setup", ns)(par, Model("survival"))
ages <- c(originLo, (originLo + originHi) / 2, originHi)
eAt <- ExtinctionProbability(ages, set$Q, set$lambda, set$mu, set$psi, 0)
survive <- 1 - as.vector(eAt %*% set$prior)
cat(sprintf("Survival of the stem at origin age %g / %g / %g: %.4f / %.4f / %.4f\n",
            ages[1], ages[2], ages[3], survive[1], survive[2], survive[3]))

SimOne <- function(condition) {
  set <- .Setup(par, Model(condition))
  repeat {
    T0 <- stats::runif(1, originLo, originHi)
    repeat {
      k <- .Sample1(set$prior)
      h <- tryCatch(SimulateDCET(T0, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                                 start = k, rho = 0, maxEdges = 5000),
                    error = function(e) NULL)
      if (is.null(h) || condition == "none" || nrow(h$fossils) > 0) break
    }
    if (is.null(h) || nrow(h$fossils) != nFossil) next
    tree <- ReconstructedTree(h, gamma = 1)
    if (is.null(tree)) next
    inv <- ReconstructedTree(h, gamma = -1)
    tt <- TempoTree(tree, present = T0)
    if (tt$sampled[tt$root - tt$nTip] || any(tt$tipType == "extant")) next
    stats <- c(.TreeStats(tt, attr(tree, "characterLength"), attr(inv, "characterLength")),
               originAge = T0)
    return(list(tree = tree, T0 = T0, stats = stats))
  }
}

cl <- makeCluster(nCore)
invisible(clusterCall(cl, function() {
  if (nzchar(Sys.getenv("EVOTEMPO_LOAD_ALL"))) {
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
  } else {
    suppressMessages(library(EvoTempo))
  }
  ns <- asNamespace("EvoTempo")
  for (fn in c(".Setup", ".Sample1", ".TreeStats")) assign(fn, get(fn, ns), envir = globalenv())
  NULL
}))
clusterExport(cl, c("SimOne", "Model", "par", "nFossil", "originLo", "originHi"))
perWorker <- ceiling(nSim / nCore)
clusterExport(cl, c("perWorker", "seedOffset"))
sims <- list()
simTime <- list()
for (cond in c("survival", "none")) {
  t0 <- Sys.time()
  clusterExport(cl, "cond")
  s <- parLapply(cl, seq_len(nCore), function(i) {
    set.seed(87000 + i + seedOffset + if (cond == "none") 500 else 0)
    do.call(rbind, lapply(seq_len(perWorker), function(j) SimOne(cond)$stats))
  })
  sims[[cond]] <- as.data.frame(do.call(rbind, s))
  simTime[[cond]] <- difftime(Sys.time(), t0, units = "mins")
  cat(sprintf("simulated %d trees (%s) in %s\n", nrow(sims[[cond]]), cond, format(simTime[[cond]])))
}

for (fn in c(".Setup", ".Sample1", ".TreeStats")) assign(fn, get(fn, ns))
set.seed(87900 + seedOffset)
start <- SimOne("survival")
tt0 <- TempoTree(start$tree, present = start$T0)
tt <- TempoTree(start$tree, present = start$T0, ranges = cbind(0, rep(1e3, tt0$nTip)))
clusterExport(cl, c("tt", "agePriors", "moves", "nIter", "burnin", "thin", "nSlice"))
chains <- list()
chainTime <- list()
for (cond in c("survival", "none")) {
  t0 <- Sys.time()
  clusterExport(cl, "cond")
  chains[[cond]] <- parLapply(cl, seq_len(nChain), function(i) {
    set.seed(87700 + i + seedOffset + if (cond == "none") 500 else 0)
    run <- TempoMCMC(tt, Model(cond), par, list(), nIter = nIter, burnin = burnin, thin = thin,
                     nSlice = nSlice, moves = moves, agePriors = agePriors)
    list(samples = run$samples, acceptance = run$acceptance)
  })
  chainTime[[cond]] <- difftime(Sys.time(), t0, units = "mins")
  cat(sprintf("ran %d chains (%s) in %s\n", nChain, cond, format(chainTime[[cond]])))
}
stopCluster(cl)
saveRDS(list(sims = sims, chains = chains), file.path(outDir, paste0("window-raw", tag, ".rds")))

Score <- function(sims, chains) {
  quantities <- intersect(c("treeLength", "meanNodeAge", "rootAge", "originAge", "oldestFossil",
                            "youngestFossil", "clock", "inverseClock", "nAncestor", "cherries"),
                          names(sims))
  # meanNodeAge is NA for a tree with no node below the root (two tips); such
  # draws are left out of that statistic's tests on both sides
  quantities <- quantities[vapply(quantities, function(q) {
    isTRUE(stats::var(sims[[q]], na.rm = TRUE) > 0)
  }, TRUE)]
  rows <- list()
  for (q in quantities) {
    s <- sims[[q]]
    s <- s[is.finite(s)]
    per <- lapply(chains, function(ch) ch$samples[[q]][is.finite(ch$samples[[q]])])
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
    for (pr in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
      qs <- stats::quantile(s, pr, names = FALSE)
      pc <- vapply(per, function(x) mean(x <= qs), 0)
      rows[[length(rows) + 1L]] <- data.frame(
        quantity = q, stat = sprintf("q%02d", round(100 * pr)), sim = pr, mcmc = mean(pc),
        z = (mean(pc) - pr) / sqrt(stats::var(pc) / length(pc) + pr * (1 - pr) / length(s)))
    }
  }
  do.call(rbind, rows)
}
out <- do.call(rbind, lapply(c("survival", "none"), function(mc) {
  do.call(rbind, lapply(c("survival", "none"), function(sm) {
    cbind(mcmc = mc, simulation = sm, Score(sims[[sm]], chains[[mc]]))
  }))
}))
write.csv(out, file.path(outDir, paste0("window-check", tag, ".csv")), row.names = FALSE)
verdict <- do.call(rbind, lapply(split(out, paste(out$mcmc, out$simulation)), function(d) {
  data.frame(mcmc = d$mcmc[1], simulation = d$simulation[1], tests = nrow(d),
             maxAbsZ = max(abs(d$z)), at = paste(d$quantity, d$stat)[which.max(abs(d$z))])
}))
verdict$role <- ifelse(verdict$mcmc == verdict$simulation, "own target", "negative control")
verdict$criterion <- ifelse(verdict$role == "own target", "max |z| < 3.5", "max |z| >= 5")
verdict$pass <- ifelse(verdict$role == "own target", verdict$maxAbsZ < 3.5, verdict$maxAbsZ >= 5)
print(verdict, row.names = FALSE, digits = 3)
write.csv(verdict, file.path(outDir, paste0("window-verdict", tag, ".csv")), row.names = FALSE)
acceptance <- lapply(chains, function(ch) colMeans(do.call(rbind, lapply(ch, `[[`, "acceptance"))))
print(lapply(acceptance, function(a) round(a[a > 0], 3)))
saveRDS(list(par = par, nFossil = nFossil, origin = c(originLo, originHi), survive = survive,
             sims = sims, out = out, verdict = verdict, acceptance = acceptance, nIter = nIter,
             nChain = nChain, simTime = simTime, chainTime = chainTime,
             chainMeans = lapply(chains, function(ch) t(vapply(ch, function(x) colMeans(x$samples), numeric(ncol(ch[[1]]$samples)))))),
        file.path(outDir, paste0("window-check", tag, ".rds")))
cat(sprintf("check-window.R: own-target max |z| %.2f / %.2f; negative-control min of max |z| %.2f\n",
            verdict$maxAbsZ[verdict$role == "own target"][1],
            verdict$maxAbsZ[verdict$role == "own target"][2],
            min(verdict$maxAbsZ[verdict$role != "own target"])))
