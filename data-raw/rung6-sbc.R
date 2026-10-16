# Rung 6: simulation-based calibration (Talts et al. 2018) per (beta, gamma)
# cell, on fixed time trees of 20-40 tips with binary characters.
#
# Each replicate draws parameters from the prior, simulates a crown clade with
# the exact dCET simulator (both crown daughters with sampled extant
# descendants), keeps it if it has 20-40 tips (sampled ancestors excluded),
# simulates characters on the realized tempo-weighted branch lengths, and runs
# TempoMCMC() on the time tree and characters from an independent prior draw.
# Rank statistics: each parameter, the root bin (ties broken at random) and
# the path functional sum_e int_e e^x dt.  Rejection on tree size depends on
# the data only, so posterior ranks remain uniform; survival is instead imposed
# by resimulating at fixed parameters, matching the conditioned likelihood.
# Run from the package root: Rscript data-raw/rung6-sbc.R [nRep] [nIter]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
nRep <- if (length(args)) as.integer(args[[1]]) else 25L
nIter <- if (length(args) > 1) as.integer(args[[2]]) else 2400L
burnin <- 500L
nDraw <- 50L
thin <- nIter %/% nDraw
nChar <- 30L
crown <- 8
tipRange <- c(20, 40)
outDir <- "data-raw/rung6-results"
dir.create(file.path(outDir, "replicates"), recursive = TRUE, showWarnings = FALSE)

meanlog <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.05), theta = log(0.2),
             s = log(0.6), rate = log(0.1))
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.25, rate = 0.4)
priors <- TempoPriors(meanlog, sdlog)
PriorDraw <- function() {
  as.list(setNames(stats::rlnorm(6, meanlog, sdlog), names(meanlog)))
}
cells <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))

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

Replicate <- function(job) {
  cell <- cells[job$cell, ]
  model <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5, beta = cell$beta,
                     gamma = cell$gamma)
  set.seed(job$seed)
  tries <- 0
  repeat {
    tries <- tries + 1
    truth <- PriorDraw()
    set <- .Setup(truth, model)
    # The likelihood is conditioned on survival, so resimulate at the same
    # parameters until both crown daughters have sampled extant descendants;
    # redrawing the parameters instead would weight them by P(survival)
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
  chars <- SimulateCharacters(tree, nChar, 2)
  clockTrue <- sum(attr(ReconstructedTree(history, gamma = 1, rate = 1),
                        "characterLength"))
  tree$root.edge <- NULL
  tt <- TempoTree(tree, present = crown)
  t0 <- Sys.time()
  run <- TempoMCMC(tt, model, PriorDraw(), priors, nIter = nIter, burnin = burnin,
                   thin = thin, nSlice = 16,
                   moves = c(tree = 1, node = 3, tip = 3, branch = 3),
                   characters = chars)
  post <- run$samples
  rankOf <- function(draws, value) {
    sum(draws < value) + sample.int(sum(draws == value) + 1L, 1L) - 1L
  }
  ranks <- c(vapply(names(meanlog), function(p) rankOf(post[[p]], truth[[p]]), 0),
             rootBin = rankOf(post$rootBin, k), clock = rankOf(post$clock, clockTrue))
  out <- list(job = job, truth = c(unlist(truth), rootBin = k, clock = clockTrue),
              ranks = ranks, nTip = nTip, nAncestor = sum(attr(tree, "tipType") == "ancestor"),
              tries = tries, acceptance = run$acceptance,
              ess = vapply(post[c(names(meanlog), "clock")], function(x) {
                coda::effectiveSize(x)
              }, 0),
              posteriorMean = colMeans(post[c(names(meanlog), "clock")]),
              runtime = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, file.path(outDir, "replicates",
                         sprintf("cell%d-rep%03d.rds", job$cell, job$rep)))
  out
}

jobs <- expand.grid(rep = seq_len(nRep), cell = seq_len(nrow(cells)))
jobs$seed <- 6000 + seq_len(nrow(jobs))
# Interleave cells, so that an interrupted run leaves balanced cells
jobs <- jobs[order(jobs$rep, jobs$cell), ]
jobList <- split(jobs, seq_len(nrow(jobs)))
done <- file.exists(file.path(outDir, "replicates",
                              sprintf("cell%d-rep%03d.rds", jobs$cell, jobs$rep)))
cl <- makeCluster(16)
invisible(clusterCall(cl, function() {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
}))
clusterExport(cl, c("cells", "meanlog", "sdlog", "priors", "PriorDraw", "Join",
                    "crown", "tipRange", "nChar", "nIter", "burnin", "thin",
                    "outDir"))
t0 <- Sys.time()
invisible(clusterApplyLB(cl, jobList[!done], Replicate))
stopCluster(cl)
files <- list.files(file.path(outDir, "replicates"), full.names = TRUE)
results <- lapply(files, readRDS)
saveRDS(list(results = results, cells = cells, nDraw = nDraw, nIter = nIter,
             burnin = burnin, thin = thin, nChar = nChar, crown = crown,
             priors = list(meanlog = meanlog, sdlog = sdlog),
             runtime = difftime(Sys.time(), t0, units = "mins")),
        file.path(outDir, "rung6.rds"))
