# Exact-target oracles for TempoMCMC() paths that no other test runs (#51):
# psi compensation by `shift` (delta != 0), `shift` with `s` free, and the
# gamma = 0 default schedule.  Each returns batch-means z-scores of pooled
# chains against an exact posterior from TreeLikelihood() or MkLikelihood()
# on a grid.  Also sourced by data-raw/oracles-51.R, which sweeps seeds and
# mutations.
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
oracleTree <- TempoTree(structure(
  list(
    edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9), c(9, 4), c(9, 5)),
    edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
    tip.label = c("a", "b", "f", "c", "d"), Nnode = 4
  ),
  class = "phylo"
))
oraclePar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.4)
oracleChar <- c(0, 1, 1, 0, 1)
oracleChars <- matrix(c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
                        1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1), 5)

# Batch-means variance of a chain mean
BatchVar <- function(x, nBatch = 25) {
  m <- floor(length(x) / nBatch) * nBatch
  stats::var(colMeans(matrix(x[seq_len(m)], ncol = nBatch))) / nBatch
}

# z of pooled chain means against exact values, MCSE by batch means
PooledZ <- function(chains, exact) {
  n <- vapply(chains, length, 0)
  est <- sum(vapply(chains, sum, 0)) / sum(n)
  se <- sqrt(sum(vapply(chains, BatchVar, 0) * n^2)) / sum(n)
  (est - exact) / se
}

# Normalized posterior on an equally spaced log-scale grid (trapezoid rule)
GridPosterior <- function(logGrid, logPost) {
  w <- exp(logPost - max(logPost))
  h <- diff(logGrid)
  mass <- c(0, cumsum((w[-1] + w[-length(w)]) / 2 * h))
  list(x = logGrid, cdf = mass / mass[[length(mass)]],
       mean = sum((logGrid[-1] * w[-1] + logGrid[-length(w)] * w[-length(w)]) / 2 * h) /
         mass[[length(mass)]])
}

# z-scores of chains of a positive parameter against its exact posterior:
# the mean of log(p), and P(log p <= q) at the exact quartiles q
MarginalZ <- function(chains, post, probs = c(0.25, 0.5, 0.75)) {
  logChains <- lapply(chains, log)
  q <- stats::approx(post$cdf, post$x, probs, ties = "ordered")$y
  c(mean = PooledZ(logChains, post$mean),
    setNames(vapply(seq_along(probs), function(i) {
      PooledZ(lapply(logChains, function(x) as.numeric(x <= q[[i]])), probs[[i]])
    }, 0), paste0("q", 100 * probs)))
}

# Chains of different seeds share no random numbers
OracleChains <- function(seed, nChain, ...) {
  lapply(seq_len(nChain), function(i) {
    set.seed(1000 * seed + i)
    TempoMCMC(oracleTree, ...)
  })
}

# Only psi free, delta = 2: every accepted shift multiplies psi by exp(-2 d).
# Exact target: TreeLikelihood() x prior on a grid of psi.
OraclePsi <- function(seed = 1, nIter = 3000, nChain = 1, burnin = 200) {
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, delta = 2)
  priors <- TempoPriors(c(psi = log(0.1)), c(psi = 1))
  logGrid <- seq(log(0.1) - 4.5, log(0.1) + 4.5, length.out = 161)
  logPost <- vapply(logGrid, function(lp) {
    TreeLikelihood(oracleTree, replace(oraclePar, "psi", exp(lp)), model) +
      priors$psi(exp(lp))
  }, 0) + logGrid
  runs <- OracleChains(seed, nChain, model, oraclePar, priors, nIter = nIter,
                       burnin = burnin, nSlice = 4,
                       moves = c(node = 1, tip = 1, branch = 1, shift = 3))
  list(z = MarginalZ(lapply(runs, function(r) r$samples$psi),
                     GridPosterior(logGrid, logPost)),
       shift = mean(vapply(runs, function(r) r$acceptance[["shift"]], 0)))
}

# Marginals of a posterior on an equally spaced log-scale grid (rows: first
# parameter, columns: second), in the form of GridPosterior()
GridMarginals <- function(logGrid1, logGrid2, logPost) {
  top <- max(logPost)
  list(GridPosterior(logGrid1, log(rowSums(exp(logPost - top)))),
       GridPosterior(logGrid2, log(colSums(exp(logPost - top)))))
}

# `s` and `rate` free with one character, gamma = 2: every accepted shift
# multiplies `rate` by exp(-2 d), d = m h s at the current `s`, which keeps the
# character lengths and so the stored character likelihood.
# Exact target: TreeLikelihood(character = ) x prior on a grid of (s, rate).
OracleS <- function(seed = 1, nIter = 6000, nChain = 2, burnin = 200, nGrid = 33,
                    moves = c(node = 1, tip = 1, branch = 1, shift = 3)) {
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, gamma = 2)
  priors <- TempoPriors(c(s = log(0.6), rate = log(0.3)), c(s = 0.3, rate = 1))
  logS <- log(0.6) + seq(-4, 4, length.out = nGrid) * 0.3
  logRate <- log(0.3) + seq(-4, 4, length.out = nGrid)
  logPost <- outer(seq_len(nGrid), seq_len(nGrid), Vectorize(function(i, j) {
    par <- replace(oraclePar, c("s", "rate"), list(exp(logS[[i]]), exp(logRate[[j]])))
    TreeLikelihood(oracleTree, par, model, character = oracleChar) +
      priors$s(par$s) + priors$rate(par$rate) + logS[[i]] + logRate[[j]]
  }))
  post <- GridMarginals(logS, logRate, logPost)
  runs <- OracleChains(seed, nChain, model, replace(oraclePar, "s", 0.6), priors,
                       nIter = nIter, burnin = burnin, nSlice = 4, moves = moves,
                       characters = matrix(oracleChar))
  z <- c(MarginalZ(lapply(runs, function(r) r$samples$s), post[[1]]),
         MarginalZ(lapply(runs, function(r) r$samples$rate), post[[2]]))
  names(z) <- paste0(rep(c("s.", "rate."), each = 4), names(z))
  list(z = z, shift = mean(vapply(runs, function(r) r$acceptance[["shift"]], 0)))
}

# Largest difference, over `nShift` shift moves with `s` and `rate` free at
# s = 0.6, between the character likelihood that a shift keeps and its value
# recomputed from the shifted state, with six characters.  One character is
# too little data for OracleS() to see an error in d; this sees any.
ShiftIdentityError <- function(seed = 1, nShift = 40) {
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)
  par <- replace(oraclePar, "s", 0.6)
  priors <- TempoPriors(c(s = log(0.6), rate = log(0.3)), c(s = 0.3, rate = 0.5))
  CharLL <- function(tr, paths, p, set) {
    MkLikelihood(tr, oracleChars,
                 .CharacterLengths(tr, set$grid$x, paths, p$rate, model$gamma))
  }
  set.seed(seed)
  set <- .AddSlices(oracleTree, .Precompute(oracleTree, par, model), 4)
  msg <- .Messages(oracleTree, set)
  d <- .SampleTree(oracleTree, set, msg)
  st <- list(tree = oracleTree, par = par, set = set, msg = msg, paths = d$paths,
             bin = d$bin)
  st$charLL <- CharLL(oracleTree, st$paths, par, set)
  err <- 0
  nAccepted <- 0
  for (i in seq_len(nShift)) {
    res <- .ShiftMove(st, model, c("s", "rate"), priors, 4, CharLL)
    st <- res$state
    nAccepted <- nAccepted + res$accepted
    err <- max(err, abs(st$charLL - CharLL(st$tree, st$paths, st$par, st$set)))
  }
  c(error = err, accepted = nAccepted)
}

# The default schedule at gamma = 0 (`tree`, `node`, `tip`, `branch`, `joint`,
# `shift`), with `theta` and `rate` free and six characters: the characters
# depend on `rate` only, and `shift` compensates `theta`, not `rate`.
# Exact target: the posterior factorizes into theta | tree, from
# TreeLikelihood(), and rate | characters, from MkLikelihood() on the
# path-independent lengths.
OracleGamma0 <- function(seed = 1, nIter = 3000, nChain = 2, burnin = 200) {
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, gamma = 0)
  priors <- TempoPriors(c(theta = log(0.1), rate = log(0.3)), c(theta = 0.75, rate = 0.5))
  runs <- OracleChains(seed, nChain, model, oraclePar, priors, nIter = nIter,
                       burnin = burnin, nSlice = 4, characters = oracleChars)
  # Lengths per unit rate at gamma = 0: the durations of the branches with paths
  nodes <- .PathNodes(oracleTree)
  age <- oracleTree$age
  unitLengths <- replace(numeric(length(age)), nodes, age[oracleTree$parent[nodes]] - age[nodes])
  logTheta <- log(0.1) + seq(-4, 4, length.out = 161) * 0.75
  logRate <- log(0.3) + seq(-4, 4, length.out = 161) * 0.5
  postTheta <- vapply(logTheta, function(lt) {
    TreeLikelihood(oracleTree, replace(oraclePar, "theta", exp(lt)), model) +
      priors$theta(exp(lt))
  }, 0) + logTheta
  postRate <- vapply(logRate, function(lr) {
    MkLikelihood(oracleTree, oracleChars, exp(lr) * unitLengths) + priors$rate(exp(lr))
  }, 0) + logRate
  z <- c(MarginalZ(lapply(runs, function(r) r$samples$theta),
                   GridPosterior(logTheta, postTheta)),
         MarginalZ(lapply(runs, function(r) r$samples$rate),
                   GridPosterior(logRate, postRate)))
  names(z) <- paste0(rep(c("theta.", "rate."), each = 4), names(z))
  list(z = z, moves = runs[[1]]$moves,
       shift = mean(vapply(runs, function(r) r$acceptance[["shift"]], 0)))
}
