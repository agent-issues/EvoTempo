# Rung 4: data augmentation vs the exact marginal likelihood (no characters).
#   4a(i)  DA posterior of (lambda, mu, psi, theta, s) on a fixed tree vs MCMC on
#          the marginal likelihood of rung 3;
#   4a(ii) at fixed parameters, path summaries from a chain of local path moves
#          vs importance-weighted FFBS stochastic maps;
#   4c     bridge Laplace transform (Lepage et al. 2006 identity);
#   4d     sum over paths (Proposition 3.3) by prior-path Monte Carlo.
# 4b (one character, product-space oracle) is in data-raw/rung5-validation.R.
# Run from the package root: Rscript data-raw/rung4-validation.R [nIter] [joint]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args)) as.integer(args[[1]]) else 5000L
# Second argument "joint": DA chains add joint parameter-path moves (#20)
useJoint <- length(args) > 1 && args[[2]] == "joint"
burnin <- 500L
nChain <- 8L
outDir <- "data-raw/rung4-results"
dir.create(outDir, showWarnings = FALSE)
# EVOTEMPO_GRID=default runs at the package default grid (K = 26, z in [-5, 5])
# and writes "-K26" result files (issue #17)
defaultGrid <- identical(Sys.getenv("EVOTEMPO_GRID"), "default")
suffix <- if (defaultGrid) "-K26" else ""

# K = 16 on z in [-3, 3] by default to keep the marginal chains affordable;
# both samplers use the same discretized model.
model <- if (defaultGrid) {
  DCETModel(rho = 0.5)
} else {
  DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5)
}
truth <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.1, s = 0.8)

# Data tree: crown clade, two daughters of the root bin simulated forwards
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
crown <- 6
set.seed(41)
setTrue <- .Setup(truth, model)
# At the default grid, reuse the K = 16 data tree
if (defaultGrid) tree <- readRDS(file.path(outDir, "rung4.rds"))$tree
while (!defaultGrid) {
  k <- .Sample1(setTrue$prior)
  h <- lapply(1:2, function(i) {
    SimulateDCET(crown, setTrue$grid$x, setTrue$Q, setTrue$lambda, setTrue$mu,
                 setTrue$psi, start = k, rho = model$rho)
  })
  if (any(vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0) == 0)) next
  tree <- ReconstructedTree(Join(h[[1]], h[[2]], setTrue$grid$x[k]))
  types <- attr(tree, "tipType")
  if (length(types) %in% 8:11 && sum(types == "fossil") >= 1 &&
      sum(types == "ancestor") >= 1) {
    break
  }
}
tree$root.edge <- NULL
tt <- TempoTree(tree, present = crown)
priors <- TempoPriors()

# 4a(i) -----------------------------------------------------------------------
Chain <- function(i, kind, nIter, burnin) {
  set.seed(1000 * i + (kind == "DA"))
  if (kind == "DA") {
    TempoMCMC(tt, model, truth, priors, nIter = nIter, burnin = burnin,
              nSlice = 20, moves = c(tree = 1, node = 1, tip = 1, branch = 1,
                                     joint = as.numeric(useJoint)))[
                c("samples", "acceptance")]
  } else {
    list(samples = MarginalMCMC(tt, model, truth, priors, nIter = nIter,
                                burnin = burnin))
  }
}
jobs <- expand.grid(i = seq_len(nChain), kind = c("DA", "marginal"),
                    stringsAsFactors = FALSE)
cl <- makeCluster(nrow(jobs))
invisible(clusterCall(cl, function() {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
}))
clusterExport(cl, c("tt", "model", "truth", "priors", "Chain", "nIter", "burnin", "jobs",
                    "suffix", "useJoint"))
t0 <- Sys.time()
chains <- parLapply(cl, seq_len(nrow(jobs)), function(j) {
  out <- Chain(jobs$i[j], jobs$kind[j], nIter, burnin)
  cat(format(Sys.time()), "chain", j, "done\n", file = paste0("rung4-progress", suffix, ".txt"),
      append = TRUE)
  out
})
runtime <- difftime(Sys.time(), t0, units = "mins")
saveRDS(list(jobs = jobs, chains = chains, runtime = runtime),
        file.path(outDir, paste0("rung4-chains", suffix, ".rds")))

# 4a(ii): local path moves at fixed parameters vs FFBS maps ----------------
mapJob <- parLapply(cl, 1:2, function(j) {
  set.seed(77 + j)
  if (j == 1) {
    TempoMCMC(tt, model, truth, list(), nIter = 4000, burnin = 200, nSlice = 20,
              moves = c(tree = 0, node = 2, tip = 2, branch = 2))[
                c("samples", "occupancy", "branchX", "acceptance")]
  } else {
    SampleTempoPaths(tt, truth, model, n = 4000, nSlice = 20)
  }
})
stopCluster(cl)

# 4c: bridge Laplace transform at a high-tempo frozen slice --------------
set <- .AddSlices(tt, .Precompute(tt, truth, model), 20)
ExpM <- function(A) {
  e <- eigen(A)
  Re(e$vectors %*% diag(exp(e$values)) %*% solve(e$vectors))
}
s <- 10
span <- 0.6
sl <- set$slices[[s]]
G <- .AddDiag(set$Q, 2 * set$lambda * sl$eBar - set$sigma)
clock <- exp(set$grid$x)
set.seed(5)
# Bridge end bins of the K = 16 design, mapped to the same z on other grids
BinAt <- function(k) round(1 + (k - 1) * (model$nBin - 1) / 15)
laplace <- do.call(rbind, lapply(list(BinAt(c(8, 8)), BinAt(c(7, 10)), BinAt(c(12, 10))), function(ends) {
  ell <- replicate(20000, {
    br <- .Bridge(set, s, span, ends[1], ends[2])
    sum(clock[br$k] * diff(c(0, br$t, span)))
  })
  do.call(rbind, lapply(c(0.2, 1, 5), function(u) {
    data.frame(from = ends[1], to = ends[2], u = u,
               empirical = mean(exp(-u * ell)),
               se = sd(exp(-u * ell)) / sqrt(length(ell)),
               exact = ExpM((G - u * diag(clock)) * span)[ends[1], ends[2]] /
                 ExpM(G * span)[ends[1], ends[2]])
  }))
}))
laplace$z <- (laplace$empirical - laplace$exact) / laplace$se

# 4d: Proposition 3.3 by prior paths, K = 4 ------------------------------------
set.seed(6)
K <- 4
Qr <- matrix(runif(K * K), K)
diag(Qr) <- 0
diag(Qr) <- -rowSums(Qr)
g <- -runif(K, 0, 1.5)
delta <- 0.7
nPath <- 50000
sumPaths <- matrix(0, K, K)
sumSq <- matrix(0, K, K)
for (k in seq_len(K)) {
  for (r in seq_len(nPath)) {
    state <- k
    now <- 0
    integral <- 0
    repeat {
      wait <- rexp(1, -Qr[state, state])
      if (now + wait >= delta) {
        integral <- integral + g[state] * (delta - now)
        break
      }
      integral <- integral + g[state] * wait
      now <- now + wait
      rates <- Qr[state, ]
      rates[state] <- 0
      state <- sample.int(K, 1, prob = rates)
    }
    w <- exp(integral)
    sumPaths[k, state] <- sumPaths[k, state] + w
    sumSq[k, state] <- sumSq[k, state] + w * w
  }
}
mcM <- sumPaths / nPath
mcSE <- sqrt(sumSq / nPath - mcM^2) / sqrt(nPath)
exactM <- ExpM((Qr + diag(g)) * delta)

saveRDS(list(tree = tree, tt = tt, truth = truth, model = model, priors = priors,
             jobs = jobs, chains = chains, runtime = runtime, nIter = nIter, joint = useJoint,
             burnin = burnin, mapChain = mapJob[[1]], maps = mapJob[[2]],
             laplace = laplace, prop33 = list(mc = mcM, se = mcSE, exact = exactM)),
        file.path(outDir, paste0("rung4", suffix, ".rds")))
cat("runtime", format(runtime), "\n")
print(laplace)
print(round((mcM - exactM) / mcSE, 2))
