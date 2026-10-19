# Rung 7 (issue #7): power of (beta, gamma) model selection by marginal
# likelihoods, on fixed true time trees with characters.
#
# Evidence of cell (beta, gamma) for tree T and characters D:
#   log Z(beta, gamma) = log Z_T(beta) + log R(beta, gamma)
# - Z_T(beta) = int pi(Theta) p(T | Theta) dTheta, the tree-only evidence, with
#   tempo integrated out exactly by TreeLikelihood().  beta = 0: closed-form FBD
#   likelihood over (lambda, mu, psi), integrated by adaptive Gauss-Hermite
#   quadrature (theta, s do not enter).  beta = 1: bridge sampling (Meng & Wong
#   1996; Gronau et al. 2017) on MarginalMCMC() draws over (lambda, mu, psi,
#   theta, s).
# - R(beta, 0) = int pi(c) Mk(D | c t_e) dc: one-dimensional quadrature, the
#   same for both beta.
# - R(beta, 1) = E[Mk(D | paths, c)] under p(Theta, paths | T): stepping-stone
#   sampling (Xie et al. 2011) over power posteriors
#   p(Theta, paths | T) Mk(D | paths, c)^t, t = (k / K)^(1 / 0.3), by
#   TempoMCMC() on the fixed tree.  Thermodynamic integration on the same
#   draws is reported as a ladder diagnostic.
#
# TempoMCMC() has no power argument and does not report the character
# log-likelihood per sample, so this script (not the package) rebinds two
# functions inside the EvoTempo namespace for the duration of a run:
# - MkLikelihood -> power * MkLikelihood (every acceptance ratio of TempoMCMC()
#   reaches the characters through it);
# - .TreeStats -> a wrapper that records the per-node clock lengths
#   int_e e^x dt that TempoMCMC() passes to it at every stored sample, from
#   which the untempered Mk log-likelihood of each draw is recomputed.
# Both are checked on every chain: the recorded log-likelihood of the last
# draw must equal the chain's final charLL / power.  The engine itself is not
# modified; see issue #7 for the engine support that would replace this.
#
# Usage, from a checkout with EvoTempo installed (or EVOTEMPO_LOAD_ALL=1):
#   Rscript data-raw/rung7-power.R pilot <id>        (ids 1-7, see Pilot())
#   Rscript data-raw/rung7-power.R task <id> <nPer>  (production; design below)
# Environment:
#   EVOTEMPO_POWER_OUT  output folder (default data-raw/rung7-results)
#   EVOTEMPO_LOAD_ALL   if set, pkgload::load_all(".") instead of library()
if (nzchar(Sys.getenv("EVOTEMPO_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
ns <- asNamespace("EvoTempo")
for (fn in c(".Setup", ".Sample1")) {
  assign(fn, get(fn, ns))
}
args <- commandArgs(trailingOnly = TRUE)
outDir <- Sys.getenv("EVOTEMPO_POWER_OUT", "data-raw/rung7-results")
dir.create(file.path(outDir, "replicates"), recursive = TRUE, showWarnings = FALSE)

# ---- Tempering and per-draw recording (script-level; see header) -----------
.mkOriginal <- get("MkLikelihood", ns)
.statsOriginal <- get(".TreeStats", ns)
.temper <- new.env()
.temper$power <- 1
.temper$lengths <- list()
.TemperedMk <- function(tree, characters, lengths, nState = 2L,
                        engine = getOption("EvoTempo.engine", "cpp")) {
  if (.temper$power == 0) {
    return(0)
  }
  .temper$power * .mkOriginal(tree, characters, lengths, nState, engine)
}
.RecordingStats <- function(tree, lengths = NULL, inverse = NULL) {
  if (!is.null(lengths)) {
    .temper$lengths[[length(.temper$lengths) + 1L]] <- lengths
  }
  .statsOriginal(tree, lengths, inverse)
}
utils::assignInNamespace("MkLikelihood", .TemperedMk, ns = "EvoTempo")
utils::assignInNamespace(".TreeStats", .RecordingStats, ns = "EvoTempo")
.TreeStats <- .statsOriginal

# ---- Model and priors --------------------------------------------------------
grid0 <- list(nBin = 16, zRange = c(-3, 3))   # as rung 6 (#21)
rho <- 0.5
nSlice <- 16L
cells <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
CellLabel <- function(i) sprintf("(%d,%d)", cells$beta[i], cells$gamma[i])
psiLevels <- c(low = log(0.02), high = log(0.2))
sizeLevels <- list(small = list(crown = 8, extant = c(15, 30)),
                   large = list(crown = 11, extant = c(40, 80)))
nCharLevels <- c(30L, 100L, 300L)
nCharMax <- max(nCharLevels)
treePar <- c("lambda", "mu", "psi", "theta", "s")
Meanlog <- function(psiLevel) {
  c(lambda = log(0.5), mu = log(0.2), psi = psiLevels[[psiLevel]],
    theta = log(0.2), s = log(0.6), rate = log(0.1))
}
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.5, rate = 0.4)
Model <- function(beta, gamma) {
  DCETModel(nBin = grid0$nBin, zRange = grid0$zRange, rho = rho,
            beta = beta, gamma = gamma)
}
moves <- c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1)

# ---- Utilities ----------------------------------------------------------------
LogSumExp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}
LogMeanExp <- function(x) LogSumExp(x) - log(length(x))

# Geyer initial monotone sequence ESS (as rung 6)
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

Secs <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, secs = proc.time()[["elapsed"]] - t0)
}

# Branch durations above each node (tips, then internal nodes); 0 at the root
Durations <- function(tt) {
  d <- tt$age[tt$parent] - tt$age
  d[is.na(d)] <- 0
  d
}

# ---- Simulation (as rung 6, #21) ----------------------------------------------
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

# Survival imposed by resimulating at fixed parameters; rejection on the number
# of sampled extant tips (a data feature) redraws the parameters.  `tries`
# counts parameter draws, so 1 / mean(tries) estimates P(size in range | cell).
Simulate <- function(cell, meanlog, size) {
  model <- Model(cells$beta[cell], cells$gamma[cell])
  tries <- 0
  repeat {
    tries <- tries + 1
    truth <- as.list(setNames(stats::rlnorm(6, meanlog, sdlog), names(meanlog)))
    set <- .Setup(truth, model)
    repeat {
      k <- .Sample1(set$prior)
      h <- tryCatch(lapply(1:2, function(i) {
        SimulateDCET(size$crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                     start = k, rho = rho, maxEdges = 8000)
      }), error = function(e) NULL)
      if (is.null(h) ||
          all(vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0) > 0)) {
        break
      }
    }
    if (is.null(h)) next
    history <- Join(h[[1]], h[[2]], set$grid$x[k])
    tree <- ReconstructedTree(history, gamma = cells$gamma[cell], rate = truth$rate)
    nExtant <- sum(attr(tree, "tipType") == "extant")
    if (nExtant >= size$extant[1] && nExtant <= size$extant[2]) break
  }
  chars <- SimulateCharacters(tree, nCharMax, 2)
  tree$root.edge <- NULL
  tt <- TempoTree(tree, present = size$crown)
  list(truth = truth, rootBin = k, tt = tt, chars = chars, tries = tries,
       nExtant = sum(tt$tipType == "extant"), nFossil = sum(tt$tipType == "fossil"),
       nAncestor = sum(tt$tipType == "ancestor"))
}

# ---- Tree-only evidence, beta = 0: adaptive Gauss-Hermite on log scale -------
GaussHermite <- function(n) {
  # Golub-Welsch for the probabilists' weight exp(-x^2 / 2)
  b <- sqrt(seq_len(n - 1))
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1), 2:n)] <- b
  J[cbind(2:n, seq_len(n - 1))] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

# log int exp(f(u)) du over R^d, f in log space, by Laplace-centred GH
LogIntegrateGH <- function(f, u0, nNode = 20) {
  opt <- stats::optim(u0, function(u) -f(u), method = "BFGS", hessian = TRUE,
                      control = list(reltol = 1e-12, maxit = 500))
  H <- opt$hessian
  S <- solve(H)
  L <- t(chol(S))
  d <- length(u0)
  gh <- GaussHermite(nNode)
  idx <- as.matrix(expand.grid(rep(list(seq_len(nNode)), d)))
  z <- matrix(gh$x[idx], ncol = d)
  logw <- rowSums(matrix(log(gh$w[idx]), ncol = d))
  u <- sweep(z %*% t(L), 2, opt$par, `+`)
  # integrand / N(0, I) density of z, times the Jacobian |L|
  logPhi <- -0.5 * rowSums(z^2) - d / 2 * log(2 * pi)
  fu <- apply(u, 1, f)
  LogSumExp(logw + fu - logPhi) + sum(log(diag(L)))
}

LogZTree0 <- function(tt, meanlog, nNode = 20) {
  p <- c("lambda", "mu", "psi")
  f <- function(u) {
    v <- exp(u)
    FBDTreeLikelihood(tt, v[1], v[2], v[3], rho = rho) +
      sum(stats::dnorm(u, meanlog[p], sdlog[p], log = TRUE))
  }
  LogIntegrateGH(f, unname(meanlog[p]), nNode)
}

# ---- Bridge sampling ---------------------------------------------------------
# draws: matrix of posterior draws on the log scale; logPostDraws: unnormalised
# log posterior (log scale, Jacobian included) at the draws; LogPost: function.
# The first half of the draws fits a normal proposal; the second half bridges.
BridgeSample <- function(draws, logPostDraws, LogPost, nProp = NULL) {
  n <- nrow(draws)
  fitIdx <- seq_len(n %/% 2)
  post <- draws[-fitIdx, , drop = FALSE]
  l1q <- logPostDraws[-fitIdx]
  nProp <- if (is.null(nProp)) nrow(post) else nProp
  mu <- colMeans(draws[fitIdx, , drop = FALSE])
  R <- chol(stats::cov(draws[fitIdx, , drop = FALSE]))
  d <- ncol(draws)
  LogG <- function(u) {
    z <- backsolve(R, t(u) - mu, transpose = TRUE)
    -0.5 * colSums(z^2) - sum(log(diag(R))) - d / 2 * log(2 * pi)
  }
  prop <- sweep(matrix(stats::rnorm(nProp * d), nProp) %*% R, 2, mu, `+`)
  l2q <- apply(prop, 1, LogPost)
  l2q[is.na(l2q)] <- -Inf
  l1 <- l1q - LogG(post)
  l2 <- l2q - LogG(prop)
  n1 <- nrow(post)
  s1 <- n1 / (n1 + nProp)
  s2 <- nProp / (n1 + nProp)
  lstar <- stats::median(l1)
  e1 <- exp(l1 - lstar)
  e2 <- exp(l2 - lstar)
  r <- 1
  for (it in seq_len(2000)) {
    rOld <- r
    r <- mean(e2 / (s1 * e2 + s2 * r)) / mean(1 / (s1 * e1 + s2 * r))
    if (abs(log(r) - log(rOld)) < 1e-12) break
  }
  # Relative MSE (Fruhwirth-Schnatter 2004), autocorrelation of f2 by ESS
  f1 <- e2 / r / (s1 * e2 / r + s2)
  f2 <- 1 / (s1 * e1 / r + s2)
  re2 <- stats::var(f1) / nProp / mean(f1)^2 +
    stats::var(f2) / Ess(f2) / mean(f2)^2
  list(logZ = log(r) + lstar, se = sqrt(re2), iterations = it,
       nPost = n1, nProp = nProp, nFailed = sum(!is.finite(l2q)))
}

# Evidence of a fixed-tree target sampled by MarginalMCMC():
# free parameters `free`; `character` optional (product-space oracle).
LogZMarginal <- function(tt, model, meanlog, free, start, nIter, burnin,
                         character = NULL) {
  priors <- TempoPriors(meanlog[free], sdlog[free])
  run <- MarginalMCMC(tt, model, start, priors, nIter = nIter, burnin = burnin,
                      character = character)
  u <- log(as.matrix(run[free]))
  LogPrior <- function(u) sum(stats::dnorm(u, meanlog[free], sdlog[free], log = TRUE))
  # Log-normal prior density on the original scale times the Jacobian exp(u)
  # equals the normal density of u
  logPostDraws <- run$logLik + apply(u, 1, LogPrior)
  LogPost <- function(v) {
    par <- start
    par[free] <- as.list(exp(v))
    ll <- suppressWarnings(tryCatch(TreeLikelihood(tt, par, model, character),
                                    error = function(e) -Inf))
    # NaN (survival probability underflow at extreme proposals): zero density;
    # counted in nFailed
    if (!is.finite(ll)) -Inf else ll + LogPrior(v)
  }
  bs <- BridgeSample(u, logPostDraws, LogPost)
  c(bs, list(posteriorMean = colMeans(run[free]), ess = vapply(run[free], Ess, 0),
             last = as.list(run[nrow(run), free])))
}

# ---- Strict-clock character evidence: 1-D quadrature over log c ------------
LogRStrict <- function(tt, chars, meanlog) {
  d <- Durations(tt)
  m <- meanlog[["rate"]]
  s <- sdlog[["rate"]]
  f <- function(u) {
    vapply(u, function(ui) .mkOriginal(tt, chars, exp(ui) * d, 2L), 0) +
      stats::dnorm(u, m, s, log = TRUE)
  }
  coarse <- seq(m - 8 * s, m + 8 * s, length.out = 801)
  fc <- f(coarse)
  keep <- range(coarse[fc > max(fc) - 50])
  h <- diff(coarse[1:2])
  fine <- seq(max(keep[1] - h, m - 8 * s), min(keep[2] + h, m + 8 * s),
              length.out = 4001)
  ff <- f(fine)
  w <- rep(1, length(fine))
  w[c(1, length(fine))] <- 0.5
  LogSumExp(ff + log(w)) + log(diff(fine[1:2]))
}

# ---- Tempered chains and stepping-stone -------------------------------------
# One TempoMCMC() chain at power t; returns the untempered character
# log-likelihood of every stored draw, and the chain's diagnostics
TemperedChain <- function(tt, model, chars, start, priors, power, nIter, burnin,
                          thin, seed) {
  .temper$power <- power
  for (attempt in 0:2) {
    set.seed(seed + 1000003L * attempt)
    .temper$lengths <- list()
    t0 <- proc.time()[["elapsed"]]
    run <- tryCatch(
      TempoMCMC(tt, model, start, priors, nIter = nIter, burnin = burnin,
                thin = thin, nSlice = nSlice, moves = moves, characters = chars),
      error = function(e) e)
    secs <- proc.time()[["elapsed"]] - t0
    if (!inherits(run, "error")) break
    if (attempt == 2) stop(run)  # nocov
  }
  len <- .temper$lengths
  smp <- run$samples
  if (length(len) != nrow(smp)) {
    stop("Recorded ", length(len), " length vectors for ", nrow(smp), " draws")
  }
  d <- Durations(tt)
  logL <- vapply(seq_along(len), function(i) {
    ell <- if (model$gamma == 1) smp$rate[i] * len[[i]] else smp$rate[i] * d
    .mkOriginal(tt, chars, ell, 2L)
  }, 0)
  # Check the recording against the chain's own final state
  check <- if (power > 0) run$state$charLL / power - logL[length(logL)] else 0
  if (abs(check) > 1e-6 * max(1, abs(logL[length(logL)]))) {
    stop("Recorded log-likelihood disagrees with the chain state by ", check)
  }
  list(logL = logL, samples = smp, acceptance = run$acceptance,
       state = run$state, secs = secs, iterations = nIter + burnin,
       attempts = attempt + 1L, check = check)
}

Ladder <- function(K, a = 0.3) (seq(0, K) / K)^(1 / a)

# Stepping-stone and thermodynamic integration over a ladder of chains
SteppingStone <- function(ladder, logL) {
  K <- length(ladder) - 1L
  inc <- vapply(seq_len(K), function(k) {
    x <- (ladder[k + 1] - ladder[k]) * logL[[k]]
    w <- exp(x - max(x))
    c(LogMeanExp(x), stats::var(w) / (Ess(w) * mean(w)^2), Ess(w))
  }, numeric(3))
  inc[2, !is.finite(inc[2, ])] <- 0
  means <- vapply(logL, mean, 0)
  ti <- sum(diff(ladder) * (means[-1] + means[-length(means)]) / 2)
  list(logR = sum(inc[1, ]), se = sqrt(sum(inc[2, ])), increments = inc[1, ],
       incrementVar = inc[2, ], incrementEss = inc[3, ], logRTI = ti,
       meanLogL = means, essLogL = vapply(logL, Ess, 0))
}

# Full stepping-stone run for a gamma = 1 (or, for validation, gamma = 0) cell:
# chains at every rung, each started from the previous rung's final parameters
LogRLadder <- function(tt, model, chars, meanlog, start, K, nIter, burnin, thin,
                       seed) {
  priors <- TempoPriors(meanlog, sdlog)
  ladder <- Ladder(K)
  par <- start
  chains <- vector("list", K + 1L)
  for (k in seq_along(ladder)) {
    chains[[k]] <- TemperedChain(tt, model, chars, par, priors, ladder[k], nIter,
                                 burnin, thin, seed + k)
    par <- chains[[k]]$state$par
  }
  ss <- SteppingStone(ladder, lapply(chains, `[[`, "logL"))
  top <- chains[[K + 1L]]$samples
  c(ss, list(ladder = ladder, secs = vapply(chains, `[[`, 0, "secs"),
             attempts = vapply(chains, `[[`, 0L, "attempts"),
             posteriorMean = colMeans(top[names(meanlog)]),
             posteriorEss = vapply(top[names(meanlog)], Ess, 0),
             acceptance = do.call(rbind, lapply(chains, `[[`, "acceptance"))))
}

# ---- Pilot ------------------------------------------------------------------
# 1-4: timing on one tree each (size x psi), truth (1,1)
# 5-6: known-answer check with one character (product-space oracle) on a small
#      tree, truth (1,1) and (0,1): all four cells by two independent routes
# 7:   ladder resolution with many characters: gamma = 0 cell by stepping-stone
#      vs exact quadrature (nChar 100, 300), and a repeated gamma = 1 ladder
Pilot <- function(id) {
  file <- file.path(outDir, "pilot", sprintf("pilot%02d.rds", id))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  spec <- switch(id,
                 list(size = "small", psi = "low", cell = 4L),
                 list(size = "small", psi = "high", cell = 4L),
                 list(size = "large", psi = "low", cell = 4L),
                 list(size = "large", psi = "high", cell = 4L),
                 list(size = "small", psi = "high", cell = 4L),
                 list(size = "small", psi = "high", cell = 2L),
                 list(size = "small", psi = "high", cell = 3L))
  meanlog <- Meanlog(spec$psi)
  set.seed(79000L + id)
  simT <- Secs(Simulate(spec$cell, meanlog, sizeLevels[[spec$size]]))
  sim <- simT$value
  tt <- sim$tt
  out <- list(id = id, spec = spec, truth = sim$truth, rootBin = sim$rootBin,
              nExtant = sim$nExtant, nFossil = sim$nFossil,
              nAncestor = sim$nAncestor, tries = sim$tries, simSecs = simT$secs)
  start0 <- as.list(exp(meanlog))
  m1 <- Model(1, 1)
  # FBD closed form vs TreeLikelihood at beta = 0 (same conventions)
  out$fbdCheck <- FBDTreeLikelihood(tt, 0.5, 0.2, 0.05, rho = rho) -
    TreeLikelihood(tt, list(lambda = 0.5, mu = 0.2, psi = 0.05, theta = 0.2, s = 0.6),
                   Model(0, 0))
  z0 <- Secs(LogZTree0(tt, meanlog))
  out$logZT0 <- z0$value
  out$secsZT0 <- z0$secs
  if (id <= 4) {
    # Timing of each component
    zt1 <- Secs(LogZMarginal(tt, m1, meanlog, treePar, start0, nIter = 600,
                             burnin = 100))
    out$zT1 <- zt1$value
    out$secsZT1 <- zt1$secs
    out$secsTreeLik <- Secs(TreeLikelihood(tt, start0, m1))$secs
    out$logRStrict <- vapply(c(30L, nCharMax), function(n) {
      LogRStrict(tt, sim$chars[, seq_len(n), drop = FALSE], meanlog)
    }, 0)
    priors <- TempoPriors(meanlog, sdlog)
    grid <- expand.grid(power = c(0, 0.01, 0.2, 1), nChar = c(30L, nCharMax))
    out$chains <- lapply(seq_len(nrow(grid)), function(i) {
      ch <- TemperedChain(tt, m1, sim$chars[, seq_len(grid$nChar[i]), drop = FALSE],
                          start0, priors, grid$power[i], nIter = 1000,
                          burnin = 200, thin = 2, seed = 79100L + 10L * id + i)
      half <- length(ch$logL) %/% 2
      list(power = grid$power[i], nChar = grid$nChar[i], secs = ch$secs,
           secsPerIter = ch$secs / ch$iterations,
           essLogL = Ess(ch$logL), meanFirst = mean(ch$logL[seq_len(half)]),
           meanSecond = mean(ch$logL[-seq_len(half)]), sdLogL = stats::sd(ch$logL),
           essPar = vapply(ch$samples[names(meanlog)], Ess, 0),
           acceptance = ch$acceptance, logL = ch$logL)
    })
  } else if (id <= 6) {
    # Route 1 (oracle): MarginalMCMC on the exact one-character likelihood
    ch1 <- sim$chars[, 1]
    out$oracle <- lapply(1:4, function(cell) {
      b <- cells$beta[cell]
      g <- cells$gamma[cell]
      free <- c(if (b == 1) treePar else c("lambda", "mu", "psi"),
                if (g == 1 && b == 0) c("theta", "s"), "rate")
      r <- Secs(LogZMarginal(tt, Model(b, g), meanlog, free, start0, nIter = 4000,
                             burnin = 500, character = ch1))
      c(r$value, list(secs = r$secs, free = free))
    })
    # Route 2: tree evidence + character factor
    zt1 <- Secs(LogZMarginal(tt, m1, meanlog, treePar, start0, nIter = 4000,
                             burnin = 500))
    out$zT1 <- zt1$value
    out$secsZT1 <- zt1$secs
    chars1 <- sim$chars[, 1, drop = FALSE]
    out$logRStrict <- LogRStrict(tt, chars1, meanlog)
    out$ladders <- lapply(c(0, 1), function(b) {
      st <- start0
      if (b == 1) st[treePar] <- zt1$value$last
      r <- Secs(LogRLadder(tt, Model(b, 1), chars1, meanlog, st, K = 16,
                           nIter = 1000, burnin = 200, thin = 2,
                           seed = 79200L + 100L * id + 10L * b))
      c(r$value, list(secsTotal = r$secs))
    })
    route2 <- c(out$logZT0 + out$logRStrict, out$logZT0 + out$ladders[[1]]$logR,
                out$zT1$logZ + out$logRStrict, out$zT1$logZ + out$ladders[[2]]$logR)
    se2 <- c(0, out$ladders[[1]]$se, out$zT1$se, sqrt(out$zT1$se^2 + out$ladders[[2]]$se^2))
    out$compare <- data.frame(
      cell = vapply(1:4, CellLabel, ""),
      oracle = vapply(out$oracle, `[[`, 0, "logZ"),
      oracleSE = vapply(out$oracle, `[[`, 0, "se"),
      route2 = route2, route2SE = se2)
    out$compare$z <- (out$compare$route2 - out$compare$oracle) /
      sqrt(out$compare$oracleSE^2 + out$compare$route2SE^2)
  } else {
    # Ladder resolution: gamma = 0 by stepping-stone vs quadrature
    zt1 <- Secs(LogZMarginal(tt, m1, meanlog, treePar, start0, nIter = 2000,
                             burnin = 300))
    st <- start0
    st[treePar] <- zt1$value$last
    out$strict <- lapply(c(100L, nCharMax), function(n) {
      chars <- sim$chars[, seq_len(n), drop = FALSE]
      exact <- LogRStrict(tt, chars, meanlog)
      r <- Secs(LogRLadder(tt, Model(1, 0), chars, meanlog, st, K = 24,
                           nIter = 1000, burnin = 200, thin = 2,
                           seed = 79300L + n))
      c(r$value, list(exact = exact, nChar = n, secsTotal = r$secs,
                      z = (r$value$logR - exact) / r$value$se))
    })
    # Repeatability of a gamma = 1 ladder (two seeds), nChar 100
    out$repeat11 <- lapply(1:2, function(rep) {
      r <- Secs(LogRLadder(tt, Model(1, 1), sim$chars[, 1:100], meanlog, st,
                           K = 24, nIter = 1000, burnin = 200, thin = 2,
                           seed = 79400L + 1000L * rep))
      c(r$value, list(secsTotal = r$secs))
    })
  }
  out$sha <- Sys.getenv("EVOTEMPO_SHA")
  saveRDS(out, file)
  invisible(out)
}

# Every component at toy chain lengths (minutes)
Smoke <- function() {
  meanlog <- Meanlog("high")
  set.seed(79999L)
  sim <- Simulate(4L, meanlog, sizeLevels$small)
  tt <- sim$tt
  start0 <- as.list(exp(meanlog))
  chars <- sim$chars[, 1:30]
  cat("extant", sim$nExtant, "fossil", sim$nFossil, "SA", sim$nAncestor,
      "tries", sim$tries, "\n")
  cat("FBD - TreeLikelihood(beta = 0):",
      FBDTreeLikelihood(tt, 0.5, 0.2, 0.05, rho = rho) -
        TreeLikelihood(tt, list(lambda = 0.5, mu = 0.2, psi = 0.05, theta = 0.2,
                                s = 0.6), Model(0, 0)), "\n")
  cat("logZT0 GH 12 / 20:", LogZTree0(tt, meanlog, 12), LogZTree0(tt, meanlog, 20), "\n")
  z <- LogZMarginal(tt, Model(1, 1), meanlog, treePar, start0, nIter = 60, burnin = 20)
  cat("logZT1 (toy):", z$logZ, "se", z$se, "\n")
  cat("logRStrict:", LogRStrict(tt, chars, meanlog), "\n")
  priors <- TempoPriors(meanlog, sdlog)
  for (g in 0:1) for (p in c(0, 0.5, 1)) {
    ch <- TemperedChain(tt, Model(1, g), chars, start0, priors, p, nIter = 20,
                        burnin = 10, thin = 2, seed = 1L)
    cat(sprintf("gamma %d power %.1f: %.2f s/iter, mean logL %.2f, check %.2g\n",
                g, p, ch$secs / ch$iterations, mean(ch$logL), ch$check))
  }
  o <- LogZMarginal(tt, Model(1, 1), meanlog, c(treePar, "rate"), start0, nIter = 40,
                    burnin = 10, character = sim$chars[, 1])
  cat("oracle (toy):", o$logZ, "\n")
  l <- LogRLadder(tt, Model(1, 1), chars, meanlog, start0, K = 2, nIter = 20,
                  burnin = 10, thin = 2, seed = 5L)
  cat("ladder (toy): SS", l$logR, "TI", l$logRTI, "\n")
}

if (length(args) && args[[1]] == "smoke") {
  Smoke()
} else if (length(args) && args[[1]] == "pilot") {
  res <- Pilot(as.integer(args[[2]]))
  cat("pilot", res$id, "done; extant", res$nExtant, "fossil", res$nFossil, "\n")
} else if (length(args) && args[[1]] == "task") {
  stop("Production mode is added once the design is approved (issue #7)")
} else if (!interactive()) {
  stop("Usage: rung7-power.R smoke | pilot <id> | task <id> <nPerTask>")
}
