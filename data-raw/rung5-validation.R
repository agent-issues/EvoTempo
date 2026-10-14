# Rung 4b and rung 5.
#   4b  One binary character on a fixed tree: data-augmented posterior of
#       (lambda, mu, psi, theta, s, rate) vs MCMC on the exact product-space
#       (tempo bin x character state) likelihood.
#   5   Prior-only MCMC with every tree move (node ages, fossil ages within
#       [0, T], subtree slide, sampled ancestor <-> tip) and free parameters,
#       vs rejection sampling from the forward simulator SimulateDCET() under
#       identical conditioning: crown age T, both crown daughters with sampled
#       extant descendants, exactly nExtant sampled extant tips and nFossil
#       fossil samples.  Compared: parameter marginals, tree statistics and the
#       tempo-weighted lengths sum(int e^x dt), sum(int e^-x dt).
# Run from the package root: Rscript data-raw/rung5-validation.R [part] [nIter]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
part <- if (length(args)) args[[1]] else "both"
nIter <- if (length(args) > 1) as.integer(args[[2]]) else 6000L
burnin <- 500L
nChain <- 8L
outDir <- "data-raw/rung5-results"
dir.create(outDir, showWarnings = FALSE)
priors <- TempoPriors()
nCore <- 16L

# 4b ---------------------------------------------------------------------------
if (part %in% c("4b", "both")) {
  r4 <- readRDS("data-raw/rung4-results/rung4.rds")
  tt <- r4$tt
  # K = 11 on z in [-2, 2]: the product space (K x 2 states) makes the
  # marginal chains about twice as costly as in rung 4
  model4b <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.5)
  nIter4b <- 3000L
  burnin4b <- 300L
  truth4b <- c(r4$truth, rate = 0.3)
  set.seed(42)
  # One character simulated on the tempo-weighted lengths of a sampled path
  st0 <- .AddSlices(tt, .Precompute(tt, truth4b, model4b), 20)
  draw <- .SampleTree(tt, st0, .Messages(tt, st0))
  len <- .CharacterLengths(tt, st0$grid$x, draw$paths, truth4b$rate, 1)
  phylo4b <- r4$tree
  edgeLen <- len[phylo4b$edge[, 2]]
  character <- as.vector(SimulateCharacters(phylo4b, 1, 2, edgeLen))
  priors4b <- c(priors, rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  # Four chains of each kind (memory limits on the development machine)
  jobs4b <- expand.grid(i = seq_len(4), kind = c("DA", "marginal"),
                        stringsAsFactors = FALSE)
  cl <- makeCluster(nrow(jobs4b))
  invisible(clusterCall(cl, function() {
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
  }))
  clusterExport(cl, c("tt", "model4b", "truth4b", "priors4b", "character",
                      "nIter4b", "burnin4b", "jobs4b"))
  t0 <- Sys.time()
  chains4b <- parLapply(cl, seq_len(nrow(jobs4b)), function(j) {
    set.seed(2000 + j)
    if (jobs4b$kind[j] == "DA") {
      TempoMCMC(tt, model4b, truth4b, priors4b, nIter = nIter4b, burnin = burnin4b,
                nSlice = 20, moves = c(tree = 1, node = 2, tip = 2, branch = 2),
                characters = matrix(character))[c("samples", "acceptance")]
    } else {
      list(samples = MarginalMCMC(tt, model4b, truth4b, priors4b, nIter = nIter4b,
                                  burnin = burnin4b, character = character))
    }
  })
  stopCluster(cl)
  saveRDS(list(tt = tt, model = model4b, truth = truth4b, character = character,
               jobs = jobs4b, chains = chains4b, nIter = nIter4b, burnin = burnin4b,
               runtime = difftime(Sys.time(), t0, units = "mins")),
          file.path(outDir, "rung4b.rds"))
}

# 5 ----------------------------------------------------------------------------
if (part %in% c("5", "both")) {
  model5 <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5)
  crown <- 4
  nExtant <- 3
  nFossil <- 2
  # Narrower priors than the defaults keep simulated clades small
  meanlog5 <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.15),
                theta = log(0.2), s = log(0.6))
  sdlog5 <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.25)
  priors5 <- TempoPriors(meanlog5, sdlog5)
  PriorDraw <- function() {
    as.list(setNames(stats::rlnorm(5, meanlog5, sdlog5), names(meanlog5)))
  }
  # Join two daughter histories under a zero-length root edge
  Join <- function(h1, h2, x0) {
    n1 <- nrow(h1$edges)
    e1 <- h1$edges
    e2 <- h2$edges
    e1$parent <- ifelse(is.na(e1$parent), 1L, e1$parent + 1L)
    e2$parent <- ifelse(is.na(e2$parent), 1L, e2$parent + 1L + n1)
    root <- data.frame(parent = NA_integer_, start = 0, end = 0,
                       type = "speciation", sampled = NA)
    f1 <- h1$fossils
    f2 <- h2$fossils
    f1$edge <- f1$edge + 1L
    f2$edge <- f2$edge + 1L + n1
    structure(list(edges = rbind(root, e1, e2),
                   paths = c(list(list(t = 0, x = x0)), h1$paths, h2$paths),
                   fossils = rbind(f1, f2), time = h1$time),
              class = "TempoHistory", kind = "discrete")
  }
  Simulate <- function(seed, n) {
    set.seed(seed)
    out <- list()
    tries <- 0
    exploded <- 0
    while (length(out) < n) {
      tries <- tries + 1
      par <- PriorDraw()
      set <- .Setup(par, model5)
      # The likelihood is conditioned on survival: resimulate at the same
      # parameters until both crown daughters have sampled extant descendants
      repeat {
        k <- .Sample1(set$prior)
        h <- tryCatch(lapply(1:2, function(i) {
          SimulateDCET(crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                       start = k, rho = model5$rho, maxEdges = 5000)
        }), error = function(e) NULL)
        if (is.null(h)) break
        ext <- vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0)
        if (all(ext > 0)) break
      }
      if (is.null(h)) {
        # More than 5000 branches: survives, but fails the tip and fossil
        # counts in practice
        exploded <- exploded + 1
        next
      }
      fos <- vapply(h, function(hi) nrow(hi$fossils), 0)
      if (sum(ext) != nExtant || sum(fos) != nFossil) next
      joint <- Join(h[[1]], h[[2]], set$grid$x[k])
      tree <- ReconstructedTree(joint, gamma = 1)
      inv <- ReconstructedTree(joint, gamma = -1)
      tt <- TempoTree(tree, present = crown)
      out[[length(out) + 1L]] <- c(unlist(par), rootBin = k,
                                   .TreeStats(tt, attr(tree, "characterLength"),
                                              attr(inv, "characterLength")))
    }
    list(stats = do.call(rbind, out), tries = tries, exploded = exploded)
  }
  cl <- makeCluster(nCore)
  invisible(clusterCall(cl, function() {
    suppressMessages(pkgload::load_all(".", quiet = TRUE))
  }))
  clusterExport(cl, c("model5", "crown", "nExtant", "nFossil", "PriorDraw",
                      "Join", "Simulate", "meanlog5", "sdlog5"))
  simFile <- file.path(outDir, "rung5-simulations.rds")
  if (file.exists(simFile)) {
    cached <- readRDS(simFile)
    sims <- cached$sims
    simTime <- cached$simTime
  } else {
    t0 <- Sys.time()
    sims <- parLapply(cl, seq_len(nCore), function(i) Simulate(500 + i, 200))
    simTime <- difftime(Sys.time(), t0, units = "mins")
    saveRDS(list(sims = sims, simTime = simTime), simFile)
  }
  simStats <- as.data.frame(do.call(rbind, lapply(sims, `[[`, "stats")))
  # Starting tree: one simulated tree with the target counts
  set.seed(1)
  repeat {
    par <- PriorDraw()
    set <- .Setup(par, model5)
    k <- .Sample1(set$prior)
    h <- tryCatch(lapply(1:2, function(i) {
      SimulateDCET(crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                   start = k, rho = model5$rho, maxEdges = 5000)
    }), error = function(e) NULL)
    if (is.null(h)) next
    ext <- vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0)
    fos <- vapply(h, function(hi) nrow(hi$fossils), 0)
    if (all(ext > 0) && sum(ext) == nExtant && sum(fos) == nFossil) break
  }
  startTree <- ReconstructedTree(Join(h[[1]], h[[2]], set$grid$x[k]))
  startTree$root.edge <- NULL
  tt5 <- TempoTree(startTree, present = crown)
  tt5 <- TempoTree(startTree, present = crown,
                   ranges = cbind(0, ifelse(tt5$tipType == "extant", 0, crown)))
  start5 <- as.list(exp(meanlog5))
  clusterExport(cl, c("tt5", "start5", "priors5", "nIter", "burnin"))
  # Chains are cached separately from the simulations (delete to rerun)
  chainFile <- file.path(outDir, "rung5-chains.rds")
  if (file.exists(chainFile)) {
    cachedChains <- readRDS(chainFile)
    chains5 <- cachedChains$chains
    chainTime <- cachedChains$runtime
  } else {
    t0 <- Sys.time()
    chains5 <- parLapply(cl, seq_len(nChain), function(i) {
      set.seed(3000 + i)
      TempoMCMC(tt5, model5, start5, priors5, nIter = nIter, burnin = burnin,
                nSlice = 16,
                moves = c(tree = 1, node = 1, tip = 1, branch = 1, nodeAge = 2,
                          fossilAge = 2, slide = 2, saJump = 2))[c("samples", "acceptance")]
    })
    chainTime <- difftime(Sys.time(), t0, units = "mins")
    saveRDS(list(chains = chains5, runtime = chainTime), chainFile)
  }
  stopCluster(cl)
  saveRDS(list(model = model5, priors = list(meanlog = meanlog5, sdlog = sdlog5),
               crown = crown, nExtant = nExtant, nFossil = nFossil,
               sims = simStats, simTries = sum(vapply(sims, `[[`, 0, "tries")),
               simExploded = sum(vapply(sims, `[[`, 0, "exploded")),
               simTime = simTime, startTree = startTree, chains = chains5,
               nIter = nIter, burnin = burnin, runtime = chainTime),
          file.path(outDir, "rung5.rds"))
}
