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
# - log Z(1,1) - log Z(0,1): bridge between the two posteriors (below), so
#   that only R(1, 1) needs a ladder.
# - R(beta, 1) = E[Mk(D | paths, c)] under p(Theta, paths | T): stepping-stone
#   sampling (Xie et al. 2011) over power posteriors
#   p(Theta, paths | T) Mk(D | paths, c)^t, t = (k / K)^(1 / 0.3), by
#   TempoMCMC() on the fixed tree.  Thermodynamic integration on the same
#   draws is reported as a ladder diagnostic.
#
# TempoMCMC() has no power argument and does not report the character
# log-likelihood per sample, so this script (not the package) rebinds two
# functions inside the EvoTempo namespace for the duration of a run:
# - .CharacterLikelihood -> power * its character log-likelihood (every
#   acceptance ratio of TempoMCMC() reaches the characters through the function
#   it returns; #26 replaced direct MkLikelihood() calls by it);
# - .TreeStats -> a wrapper that records the per-node clock lengths
#   int_e e^x dt that TempoMCMC() passes to it at every stored sample, from
#   which the untempered Mk log-likelihood of each draw is recomputed.
# Both are checked on every chain: the recorded log-likelihood of the last
# draw must equal the chain's final charLL / power.  For the beta contrast at
# gamma = 1 the wrapper also records each draw's (parameters, paths, root and
# node bins) from the calling TempoMCMC() frame, and the internal exact
# augmented density .LogAugmented() is evaluated at them under the other beta:
#   Z(1,1) / Z(0,1) by two-sample bridge sampling (Meng & Wong 1996) between
#   the two full posteriors, on the common space (Theta, paths): the character
#   likelihood is the same function of (paths, c) in both, so only the ratio
#   p(T, paths | Theta, beta = 1) / p(T, paths | Theta, beta = 0) enters.
# The recorded states are checked against the chain's own log posterior.  The engine itself is not
# modified; see issue #7 for the engine support that would replace this.
#
# Usage, from a checkout with EvoTempo installed (or EVOTEMPO_LOAD_ALL=1):
#   Rscript data-raw/rung7-power.R pilot <id>        (ids 1-7, see Pilot())
#   Rscript data-raw/rung7-power.R task <id> <nPer>  (trees (id-1)*nPer+1..id*nPer)
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
.temper$recordState <- FALSE
.temper$states <- list()
.clOriginal <- get(".CharacterLikelihood", ns)
.TemperedCharacters <- function(characters, nState, gamma) {
  Inner <- .clOriginal(characters, nState, gamma)
  Tempered <- function(tr, paths, p, set) {
    if (.temper$power == 0) {
      return(0)
    }
    .temper$power * Inner(tr, paths, p, set)
  }
  # Keeps the incremental cache's resolve hook (#26)
  attributes(Tempered) <- attributes(Inner)
  Tempered
}
.RecordingStats <- function(tree, lengths = NULL, inverse = NULL) {
  if (!is.null(lengths)) {
    .temper$lengths[[length(.temper$lengths) + 1L]] <- lengths
    if (.temper$recordState) {
      # The chain state of the calling TempoMCMC() frame at this stored draw
      st <- get("st", envir = parent.frame(), inherits = FALSE)
      .temper$states[[length(.temper$states) + 1L]] <- st[c("par", "paths", "bin")]
    }
  }
  .statsOriginal(tree, lengths, inverse)
}
utils::assignInNamespace(".CharacterLikelihood", .TemperedCharacters, ns = "EvoTempo")
utils::assignInNamespace(".TreeStats", .RecordingStats, ns = "EvoTempo")
.TreeStats <- .statsOriginal
.LogAugmented <- get(".LogAugmented", ns)
.Precompute <- get(".Precompute", ns)

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
# Log-normal priors truncated at +/- zTrunc sdlog and renormalised (outside mass
# 6.3e-5 per parameter).  TempoMCMC() skips a proposal with -Inf prior before
# evaluating the likelihood, so the truncation keeps its burn-in adaptation from
# proposing values (e.g. s > 10) at which the E equations fail (#37).  The same
# truncated prior is used to simulate and in every evidence calculation.
zTrunc <- 4
logMass <- log(stats::pnorm(zTrunc) - stats::pnorm(-zTrunc))
TruncPriors <- function(meanlog, sdlog) {
  out <- lapply(names(meanlog), function(p) {
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
  names(out) <- names(meanlog)
  out
}
# Density of u = log(x) (log-normal density times the Jacobian x)
LogPriorU <- function(u, m, s) {
  if (any(abs(u - m) > zTrunc * s)) -Inf else
    sum(stats::dnorm(u, m, s, log = TRUE)) - length(u) * logMass
}
PriorDraw <- function(meanlog) {
  z <- stats::rnorm(length(meanlog))
  while (any(out <- abs(z) > zTrunc)) {
    z[out] <- stats::rnorm(sum(out))
  }
  as.list(setNames(exp(meanlog + sdlog[names(meanlog)] * z), names(meanlog)))
}
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
    truth <- PriorDraw(meanlog)
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
  opt <- stats::optim(u0, function(u) min(-f(u), 1e300), method = "BFGS",
                      hessian = TRUE,
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
    lp <- LogPriorU(u, meanlog[p], sdlog[p])
    if (!is.finite(lp)) -Inf else FBDTreeLikelihood(tt, v[1], v[2], v[3], rho = rho) + lp
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
# The sampler is MarginalMCMC()'s algorithm (single-parameter multiplicative
# random walk, step sizes adapted every 20 burn-in iterations to 30%
# acceptance), written here on the log scale so that a solver failure at an
# extreme proposal counts as a rejection: MarginalMCMC() stops with an error
# there (as TempoMCMC() does, #37), and restarting would condition on not
# meeting one.
ParameterChain <- function(LogLik, u0, LogPrior, nIter, burnin) {
  d <- length(u0)
  u <- u0
  ll <- LogLik(u)
  lp <- LogPrior(u)
  step <- rep(1, d)
  hits <- numeric(d)
  out <- matrix(NA_real_, nIter, d, dimnames = list(NULL, names(u0)))
  outLL <- numeric(nIter)
  nFail <- 0
  for (i in seq_len(burnin + nIter)) {
    for (j in seq_len(d)) {
      v <- u
      v[j] <- u[j] + step[j] * (stats::runif(1) - 0.5)
      lpv <- LogPrior(v)
      if (!is.finite(lpv)) next
      llv <- LogLik(v)
      if (!is.finite(llv)) {
        nFail <- nFail + 1
        next
      }
      if (log(stats::runif(1)) < llv + lpv - ll - lp) {
        u <- v
        ll <- llv
        lp <- lpv
        hits[j] <- hits[j] + 1
      }
    }
    if (i <= burnin && i %% 20 == 0) {
      step <- step * exp(hits / 20 - 0.3)
      hits <- 0 * hits
    }
    if (i > burnin) {
      out[i - burnin, ] <- u
      outLL[i - burnin] <- ll
    }
  }
  list(u = out, logLik = outLL, nFail = nFail)
}

LogZMarginal <- function(tt, model, meanlog, free, start, nIter, burnin,
                         character = NULL) {
  # Log-normal prior density on the original scale times the Jacobian exp(u)
  # equals the normal density of u
  LogPrior <- function(u) LogPriorU(u, meanlog[free], sdlog[free])
  LogLik <- function(v) {
    par <- start
    par[free] <- as.list(exp(v))
    ll <- suppressWarnings(tryCatch(TreeLikelihood(tt, par, model, character),
                                    error = function(e) -Inf))
    # NaN (survival probability underflow, solver failure): zero density
    if (is.finite(ll)) ll else -Inf
  }
  u0 <- setNames(log(unlist(start[free])), free)
  run <- ParameterChain(LogLik, u0, LogPrior, nIter, burnin)
  u <- run$u
  logPostDraws <- run$logLik + apply(u, 1, LogPrior)
  LogPost <- function(v) {
    lp <- LogPrior(v)
    if (!is.finite(lp)) return(-Inf)
    ll <- LogLik(v)
    if (is.finite(ll)) ll + lp else -Inf
  }
  bs <- BridgeSample(u, logPostDraws, LogPost)
  x <- exp(u)
  c(bs, list(posteriorMean = colMeans(x), ess = apply(u, 2, Ess),
             last = as.list(x[nrow(x), ]), chainFailures = run$nFail))
}

# ---- Strict-clock character evidence: 1-D quadrature over log c ------------
LogRStrict <- function(tt, chars, meanlog) {
  d <- Durations(tt)
  m <- meanlog[["rate"]]
  s <- sdlog[["rate"]]
  f <- function(u) {
    vapply(u, function(ui) .mkOriginal(tt, chars, exp(ui) * d, 2L), 0) +
      stats::dnorm(u, m, s, log = TRUE) - logMass
  }
  lo <- m - zTrunc * s
  hi <- m + zTrunc * s
  coarse <- seq(lo, hi, length.out = 801)
  fc <- f(coarse)
  keep <- range(coarse[fc > max(fc) - 50])
  h <- diff(coarse[1:2])
  fine <- seq(max(keep[1] - h, lo), min(keep[2] + h, hi), length.out = 4001)
  ff <- f(fine)
  w <- rep(1, length(fine))
  w[c(1, length(fine))] <- 0.5
  LogSumExp(ff + log(w)) + log(diff(fine[1:2]))
}

# ---- Tempered chains and stepping-stone -------------------------------------
# One TempoMCMC() chain at power t; returns the untempered character
# log-likelihood of every stored draw, and the chain's diagnostics
TemperedChain <- function(tt, model, chars, start, priors, power, nIter, burnin,
                          thin, seed, recordState = FALSE) {
  .temper$power <- power
  .temper$recordState <- recordState
  on.exit(.temper$recordState <- FALSE)
  for (attempt in 0:2) {
    set.seed(seed + 1000003L * attempt)
    .temper$lengths <- list()
    .temper$states <- list()
    t0 <- proc.time()[["elapsed"]]
    run <- tryCatch(
      TempoMCMC(tt, model, start, priors, nIter = nIter, burnin = burnin,
                thin = thin, nSlice = nSlice, moves = moves, characters = chars),
      error = function(e) e)
    secs <- proc.time()[["elapsed"]] - t0
    if (!inherits(run, "error")) break
    message(sprintf("TempoMCMC error at power %.4g, attempt %d: %s; start %s",
                    power, attempt + 1L, conditionMessage(run),
                    paste(names(start), signif(unlist(start), 4), sep = "=",
                          collapse = " ")))
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
       attempts = attempt + 1L, check = check,
       states = if (recordState) .temper$states)
}

# Log augmented tree density log p(T, paths | Theta) of a recorded state under
# `model`; -Inf where the model's E equations cannot be solved
LogTreeDensity <- function(tt, model, state) {
  tryCatch(.LogAugmented(tt, .Precompute(tt, state$par, model), state$paths,
                         state$bin), error = function(e) -Inf)
}

# Two-sample bridge sampling for log(Z1 / Z0) from l = log(q1 / q0) at draws
# from posterior 0 (l0) and posterior 1 (l1); optimal bridge, iterated
BridgeTwo <- function(l0, l1) {
  n0 <- length(l0)
  n1 <- length(l1)
  s0 <- n0 / (n0 + n1)
  s1 <- n1 / (n0 + n1)
  shift <- stats::median(c(l0[is.finite(l0)], l1[is.finite(l1)]))
  e0 <- exp(l0 - shift)
  e1 <- exp(l1 - shift)
  r <- 1
  for (it in seq_len(5000)) {
    rOld <- r
    r <- mean(e0 / (s1 * e0 + s0 * r)) / mean(1 / (s1 * e1 + s0 * r))
    if (abs(log(r) - log(rOld)) < 1e-12) break
  }
  f0 <- e0 / (s1 * e0 + s0 * r)
  f1 <- 1 / (s1 * e1 + s0 * r)
  re2 <- stats::var(f0) / (Ess(f0) * mean(f0)^2) +
    stats::var(f1) / (Ess(f1) * mean(f1)^2)
  list(logBF = log(r) + shift, se = sqrt(re2), iterations = it,
       essF0 = Ess(f0), essF1 = Ess(f1),
       # importance-sampling estimates from one side only (diagnostics)
       logBFfrom0 = LogMeanExp(l0), logBFfrom1 = -LogMeanExp(-l1),
       sdL0 = stats::sd(l0), sdL1 = stats::sd(l1), meanL0 = mean(l0),
       meanL1 = mean(l1), nInfinite = sum(!is.finite(c(l0, l1))))
}

# log Z(beta = 1, gamma) - log Z(beta = 0, gamma) from recorded posterior chains
# (chain0 under beta = 0, chain1 under beta = 1, both at power 1)
BetaBridge <- function(tt, chain0, chain1, gamma, meanlog) {
  m0 <- Model(0, gamma)
  m1 <- Model(1, gamma)
  Side <- function(chain, own) {
    a0 <- vapply(chain$states, function(st) LogTreeDensity(tt, m0, st), 0)
    a1 <- vapply(chain$states, function(st) LogTreeDensity(tt, m1, st), 0)
    # Own-model density must equal logPost - log prior - character logL
    priors <- TruncPriors(meanlog, sdlog)
    lp <- vapply(chain$states, function(st) {
      sum(vapply(names(priors), function(q) priors[[q]](st$par[[q]]), 0))
    }, 0)
    resid <- (if (own == 0) a0 else a1) - (chain$samples$logPost - lp - chain$logL)
    list(l = a1 - a0, check = max(abs(resid)))
  }
  t0 <- proc.time()[["elapsed"]]
  s0 <- Side(chain0, 0)
  s1 <- Side(chain1, 1)
  c(BridgeTwo(s0$l, s1$l), list(check0 = s0$check, check1 = s1$check,
                                secs = proc.time()[["elapsed"]] - t0))
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
# The top rung (power 1) is the posterior; it can run longer (nIterTop,
# thinTop) and record its states for BetaBridge()
LogRLadder <- function(tt, model, chars, meanlog, start, K, nIter, burnin, thin,
                       seed, nIterTop = nIter, thinTop = thin, recordTop = FALSE) {
  priors <- TruncPriors(meanlog, sdlog)
  ladder <- Ladder(K)
  par <- start
  chains <- vector("list", K + 1L)
  for (k in seq_along(ladder)) {
    top <- k == length(ladder)
    chains[[k]] <- TemperedChain(tt, model, chars, par, priors, ladder[k],
                                 if (top) nIterTop else nIter, burnin,
                                 if (top) thinTop else thin, seed + k,
                                 recordState = top && recordTop)
    par <- chains[[k]]$state$par
  }
  ss <- SteppingStone(ladder, lapply(chains, `[[`, "logL"))
  top <- chains[[K + 1L]]$samples
  c(ss, list(ladder = ladder, top = chains[[K + 1L]],
             secs = vapply(chains, `[[`, 0, "secs"),
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
#      vs exact quadrature (nChar 100, 300)
# 8-10: beta contrast at gamma = 1 (nChar 30, 100, 300): bridge between the
#      two posteriors vs difference of two ladders; bridge repeatability
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
                 list(size = "small", psi = "high", cell = 3L),
                 list(size = "small", psi = "high", cell = 4L),
                 list(size = "small", psi = "high", cell = 4L),
                 list(size = "small", psi = "high", cell = 4L))
  meanlog <- Meanlog(spec$psi)
  set.seed(79000L + min(id, 8L))  # pilots 8-10 share one tree
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
    priors <- TruncPriors(meanlog, sdlog)
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
    # Route 1 (oracle): parameter chain on the exact one-character likelihood
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
    # Route 2: tree evidence + character factor; beta contrast by bridge too
    zt1 <- Secs(LogZMarginal(tt, m1, meanlog, treePar, start0, nIter = 4000,
                             burnin = 500))
    out$zT1 <- zt1$value
    out$secsZT1 <- zt1$secs
    chars1 <- sim$chars[, 1, drop = FALSE]
    out$logRStrict <- LogRStrict(tt, chars1, meanlog)
    ladders <- lapply(c(0, 1), function(b) {
      st <- start0
      if (b == 1) st[treePar] <- zt1$value$last
      r <- Secs(LogRLadder(tt, Model(b, 1), chars1, meanlog, st, K = 16,
                           nIter = 1000, burnin = 200, thin = 2,
                           seed = 79200L + 100L * id + 10L * b,
                           nIterTop = 4000, thinTop = 4, recordTop = TRUE))
      c(r$value, list(secsTotal = r$secs))
    })
    out$bridge <- BetaBridge(tt, ladders[[1]]$top, ladders[[2]]$top, 1, meanlog)
    out$ladders <- lapply(ladders, function(l) { l$top <- NULL; l })
    route2 <- c(out$logZT0 + out$logRStrict, out$logZT0 + out$ladders[[1]]$logR,
                out$zT1$logZ + out$logRStrict, out$zT1$logZ + out$ladders[[2]]$logR)
    se2 <- c(0, out$ladders[[1]]$se, out$zT1$se,
             sqrt(out$zT1$se^2 + out$ladders[[2]]$se^2))
    out$compare <- data.frame(
      cell = vapply(1:4, CellLabel, ""),
      oracle = vapply(out$oracle, `[[`, 0, "logZ"),
      oracleSE = vapply(out$oracle, `[[`, 0, "se"),
      route2 = route2, route2SE = se2)
    out$compare$z <- (out$compare$route2 - out$compare$oracle) /
      sqrt(out$compare$oracleSE^2 + out$compare$route2SE^2)
    out$betaContrast <- c(
      oracle = out$compare$oracle[4] - out$compare$oracle[2],
      oracleSE = sqrt(sum(out$compare$oracleSE[c(2, 4)]^2)),
      ladders = route2[4] - route2[2], laddersSE = sqrt(sum(se2[c(2, 4)]^2)),
      bridge = out$bridge$logBF, bridgeSE = out$bridge$se)
  } else if (id == 7) {
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
      r$value$top <- NULL
      c(r$value, list(exact = exact, nChar = n, secsTotal = r$secs,
                      z = (r$value$logR - exact) / r$value$se))
    })
  } else {
    # 8-10: beta contrast at gamma = 1 with many characters (30, 100, 300):
    # bridge vs difference of two ladders; repeatability of the bridge from a
    # second, independent pair of posterior chains
    n <- c(30L, 100L, nCharMax)[[id - 7L]]
    chars <- sim$chars[, seq_len(n), drop = FALSE]
    zt1 <- Secs(LogZMarginal(tt, m1, meanlog, treePar, start0, nIter = 3000,
                             burnin = 500))
    out$zT1 <- zt1$value
    out$logRStrict <- LogRStrict(tt, chars, meanlog)
    priors <- TruncPriors(meanlog, sdlog)
    st1 <- start0
    st1[treePar] <- zt1$value$last
    ladders <- lapply(c(0, 1), function(b) {
      r <- Secs(LogRLadder(tt, Model(b, 1), chars, meanlog, if (b) st1 else start0,
                           K = 16, nIter = 1000, burnin = 200, thin = 2,
                           seed = 79500L + 100L * id + 10L * b,
                           nIterTop = 8000, thinTop = 8, recordTop = TRUE))
      c(r$value, list(secsTotal = r$secs))
    })
    out$bridge <- BetaBridge(tt, ladders[[1]]$top, ladders[[2]]$top, 1, meanlog)
    tops2 <- lapply(c(0, 1), function(b) {
      TemperedChain(tt, Model(b, 1), chars,
                    ladders[[b + 1L]]$top$state$par, priors, 1, nIter = 8000,
                    burnin = 500, thin = 8, seed = 79600L + 100L * id + 10L * b,
                    recordState = TRUE)
    })
    out$bridge2 <- BetaBridge(tt, tops2[[1]], tops2[[2]], 1, meanlog)
    out$topEss <- lapply(ladders, function(l) {
      c(logL = Ess(l$top$logL), vapply(l$top$samples[names(meanlog)], Ess, 0))
    })
    out$topSecs <- vapply(tops2, `[[`, 0, "secs")
    out$ladders <- lapply(ladders, function(l) { l$top <- NULL; l })
    out$betaContrast <- c(
      ladders = out$zT1$logZ + out$ladders[[2]]$logR -
        (out$logZT0 + out$ladders[[1]]$logR),
      laddersSE = sqrt(out$zT1$se^2 + out$ladders[[1]]$se^2 + out$ladders[[2]]$se^2),
      bridge = out$bridge$logBF, bridgeSE = out$bridge$se,
      bridge2 = out$bridge2$logBF, bridge2SE = out$bridge2$se)
  }
  out$sha <- Sys.getenv("EVOTEMPO_SHA")
  saveRDS(out, file)
  invisible(out)
}

# ---- Production ---------------------------------------------------------------
# Tree g (1, 2, ...): true cell (g - 1) %% 4 + 1; fossil-density x size combination
# ((g - 1) %/% 4) %% 4 + 1 of `combos`; replicate block (g - 1) %/% 16 + 1;
# seed 70000 + g.  One tree carries nCharMax characters; the nChar levels use its
# first nChar, so the levels are nested and paired.  One task = one tree, all
# nChar levels (the tree evidence is shared); one output file per level, and
# levels already done are skipped.  In replicate block 1 the (0,1) cell is also
# estimated by its own stepping-stone ladder, to check the bridge at production
# settings.
combos <- expand.grid(psi = names(psiLevels), size = names(sizeLevels),
                      stringsAsFactors = FALSE)
EnvNum <- function(name, default) {
  v <- Sys.getenv(name)
  if (nzchar(v)) as.numeric(strsplit(v, ",")[[1]]) else default
}
prodLevels <- as.integer(EnvNum("EVOTEMPO_POWER_NCHAR", c(10, 30, 100)))
ladderK <- EnvNum("EVOTEMPO_POWER_K", 16)
ladderIter <- EnvNum("EVOTEMPO_POWER_ITER", 1000)
ladderBurn <- EnvNum("EVOTEMPO_POWER_BURN", 200)
ladderThin <- EnvNum("EVOTEMPO_POWER_THIN", 2)
treeIter <- EnvNum("EVOTEMPO_POWER_TREEITER", 3000)
topIter <- EnvNum("EVOTEMPO_POWER_TOPITER", 8000)
topThin <- EnvNum("EVOTEMPO_POWER_TOPTHIN", 8)
topBurn <- EnvNum("EVOTEMPO_POWER_TOPBURN", 500)

Task <- function(g) {
  files <- file.path(outDir, "replicates", sprintf("tree%04d-n%03d.rds", g, prodLevels))
  if (all(file.exists(files))) {
    return(invisible(NULL))
  }
  cell <- (g - 1L) %% 4L + 1L
  combo <- combos[((g - 1L) %/% 4L) %% nrow(combos) + 1L, ]
  block <- (g - 1L) %/% 16L + 1L
  meanlog <- Meanlog(combo$psi)
  t0 <- proc.time()[["elapsed"]]
  set.seed(70000L + g)
  sim <- Simulate(cell, meanlog, sizeLevels[[combo$size]])
  tt <- sim$tt
  start0 <- as.list(exp(meanlog))
  set.seed(71000L + g)
  zT0 <- Secs(LogZTree0(tt, meanlog))
  zT1 <- Secs(LogZMarginal(tt, Model(1, 1), meanlog, treePar, start0,
                           nIter = treeIter, burnin = 500))
  treeSecs <- proc.time()[["elapsed"]] - t0
  topSummary <- function(ch) {
    list(logL = ch$logL, essLogL = Ess(ch$logL),
         posteriorMean = colMeans(ch$samples[names(meanlog)]),
         posteriorEss = vapply(ch$samples[names(meanlog)], Ess, 0),
         acceptance = ch$acceptance, attempts = ch$attempts)
  }
  out <- NULL
  for (level in seq_along(prodLevels)) {
    if (file.exists(files[level])) next
    tl <- proc.time()[["elapsed"]]
    nChar <- prodLevels[[level]]
    chars <- sim$chars[, seq_len(nChar), drop = FALSE]
    rS <- Secs(LogRStrict(tt, chars, meanlog))
    # gamma contrast at beta = 1: stepping-stone ladder, recorded posterior top
    st1 <- start0
    st1[treePar] <- zT1$value$last
    seedBase <- 72000L + 10L * g + 1000000L * nChar
    lad <- Secs(LogRLadder(tt, Model(1, 1), chars, meanlog, st1, K = ladderK,
                           nIter = ladderIter, burnin = ladderBurn,
                           thin = ladderThin, seed = seedBase, nIterTop = topIter,
                           thinTop = topThin, recordTop = TRUE))
    # beta contrast at gamma = 1: posterior under (0,1), bridge to (1,1)
    top0 <- Secs(TemperedChain(tt, Model(0, 1), chars, lad$value$top$state$par,
                               TruncPriors(meanlog, sdlog), 1, nIter = topIter,
                               burnin = topBurn, thin = topThin,
                               seed = seedBase + 500000L, recordState = TRUE))
    br <- BetaBridge(tt, top0$value, lad$value$top, 1, meanlog)
    check01 <- if (block == 1L) {
      r <- Secs(LogRLadder(tt, Model(0, 1), chars, meanlog, start0, K = ladderK,
                           nIter = ladderIter, burnin = ladderBurn,
                           thin = ladderThin, seed = seedBase + 700000L))
      r$value$top <- topSummary(r$value$top)
      c(r$value, list(secs = r$secs, logZ = zT0$value + r$value$logR))
    }
    logZ11 <- zT1$value$logZ + lad$value$logR
    logZ <- c(zT0$value + rS$value, logZ11 - br$logBF,
              zT1$value$logZ + rS$value, logZ11)
    se11 <- sqrt(zT1$value$se^2 + lad$value$se^2)
    se <- c(0, sqrt(se11^2 + br$se^2), zT1$value$se, se11)
    logZTI <- c(NA, NA, NA, zT1$value$logZ + lad$value$logRTI)
    names(logZ) <- names(se) <- names(logZTI) <- vapply(1:4, CellLabel, "")
    ladder <- lad$value
    ladder$top <- topSummary(ladder$top)
    out <- list(g = g, block = block, cell = cell, beta = cells$beta[cell],
                gamma = cells$gamma[cell], psiLevel = combo$psi, size = combo$size,
                nChar = nChar, seed = 70000L + g, truth = sim$truth,
                rootBin = sim$rootBin, nExtant = sim$nExtant, nFossil = sim$nFossil,
                nAncestor = sim$nAncestor, tries = sim$tries,
                logZ = logZ, se = se, logZTI = logZTI,
                logZT0 = zT0$value, zT1 = zT1$value, logRStrict = rS$value,
                ladder = ladder, top01 = topSummary(top0$value), bridge = br,
                check01 = check01,
                secs = c(tree = if (level == 1L) treeSecs else 0,
                         strict = rS$secs, ladder11 = lad$secs, top01 = top0$secs,
                         bridge = br$secs,
                         check01 = if (is.null(check01)) 0 else check01$secs,
                         level = proc.time()[["elapsed"]] - tl),
                treeSecs = treeSecs,
                settings = list(K = ladderK, iter = ladderIter, burn = ladderBurn,
                                thin = ladderThin, topIter = topIter,
                                topThin = topThin, topBurn = topBurn,
                                treeIter = treeIter, grid = grid0, moves = moves,
                                nSlice = nSlice, sdlog = sdlog, meanlog = meanlog,
                                zTrunc = zTrunc),
                packageVersion = as.character(utils::packageVersion("EvoTempo")),
                sha = Sys.getenv("EVOTEMPO_SHA"))
    tmp <- paste0(files[level], ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, files[level])
    cat(sprintf("tree %d (cell %d, %s/%s) nChar %d: %.1f min; logZ %s\n", g, cell,
                combo$psi, combo$size, nChar, out$secs[["level"]] / 60,
                paste(sprintf("%.2f", logZ), collapse = " ")))
  }
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
  priors <- TruncPriors(meanlog, sdlog)
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
  tops <- lapply(0:1, function(b) {
    TemperedChain(tt, Model(b, 1), chars, start0, priors, 1, nIter = 30,
                  burnin = 10, thin = 2, seed = 7L + b, recordState = TRUE)
  })
  br <- BetaBridge(tt, tops[[1]], tops[[2]], 1, meanlog)
  cat("beta bridge (toy):", br$logBF, "se", br$se, "checks", br$check0, br$check1,
      "\n")
}

options(error = function() {
  calls <- sys.calls()
  message("Call stack:\n", paste(utils::tail(vapply(calls, function(x)
    paste(deparse(x, nlines = 1L), collapse = ""), ""), 15), collapse = "\n"))
  q(status = 1)
})
if (length(args) && args[[1]] == "smoke") {
  Smoke()
} else if (length(args) && args[[1]] == "dumptree") {
  # The simulated tree and characters of pilot <id>, for reproductions
  id <- as.integer(args[[2]])
  spec <- list(list(size = "small", psi = "low", cell = 4L),
               list(size = "small", psi = "high", cell = 4L),
               list(size = "large", psi = "low", cell = 4L),
               list(size = "large", psi = "high", cell = 4L),
               list(size = "small", psi = "high", cell = 4L),
               list(size = "small", psi = "high", cell = 2L),
               list(size = "small", psi = "high", cell = 3L),
               list(size = "small", psi = "high", cell = 4L))[[min(id, 8L)]]
  set.seed(79000L + min(id, 8L))
  sim <- Simulate(spec$cell, Meanlog(spec$psi), sizeLevels[[spec$size]])
  saveRDS(sim, args[[3]])
} else if (length(args) && args[[1]] == "pilot") {
  res <- Pilot(as.integer(args[[2]]))
  cat("pilot", res$id, "done; extant", res$nExtant, "fossil", res$nFossil, "\n")
} else if (length(args) && args[[1]] == "task") {
  n <- as.integer(args[[3]])
  for (g in (as.integer(args[[2]]) - 1L) * n + seq_len(n)) {
    Task(g)
  }
} else if (!interactive()) {
  stop("Usage: rung7-power.R smoke | pilot <id> | task <id> <nPerTask>")
}
