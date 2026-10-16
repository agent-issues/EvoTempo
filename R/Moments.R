#' Diversity moments
#'
#' Expected clade size \eqn{N_k(t) = E[n_t]}, and optionally its second raw
#' moment \eqn{E[n_t^2]}, for a clade founded by one lineage in bin \eqn{k}.
#'
#' Solves \eqn{\dot N = (Q + R) N} with \eqn{R = \mathrm{diag}(\lambda_k -
#' \mu_k)}, and \eqn{\dot m = (Q + R) m + 2 \Lambda N^2}, with
#' \eqn{\Lambda = \mathrm{diag}(\lambda_k)}; \eqn{N(0) = m(0) = 1}.
#' Divide by survival, `1 - ExtinctionProbability()`, to condition on
#' \eqn{n_t > 0}; \eqn{E[n^2]/E[n]} is the mean *experienced* clade size.
#'
#' @inheritParams ExtinctionProbability
#' @param second Logical; if `TRUE`, also return the second moment.
#' @returns `DiversityMoments()` returns a list with element `mean`, and
#' `meanSquare` if `second = TRUE`: matrices with one row per element of
#' `times` and one column per founding bin.
#' @examples
#' grid <- TempoGrid(41, c(-4, 4))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' DiversityMoments(100, Q, rates$lambda, rates$mu)$mean[, 21]
#' @family dCET numerics
#' @export
DiversityMoments <- function(times, Q, lambda, mu, second = FALSE,
                             rtol = 1e-10, atol = 1e-10) {
  nBin <- nrow(Q)
  lambda <- rep_len(lambda, nBin)
  growth <- .AddDiag(Q, lambda - rep_len(mu, nBin))
  if (!second) {
    return(list(mean = .Linear(times, growth, rep(1, nBin), rtol, atol)))
  }
  first <- seq_len(nBin)
  out <- .Integrate(
    rep(1, 2 * nBin), times,
    function(t, y) {
      n <- y[first]
      c(growth %*% n, growth %*% y[-first] + 2 * lambda * n * n)
    },
    function(t, y) {
      jac <- matrix(0, 2 * nBin, 2 * nBin)
      jac[first, first] <- growth
      jac[-first, -first] <- growth
      jac[cbind(first + nBin, first)] <- 4 * lambda * y[first]
      jac
    },
    rtol, atol
  )
  list(mean = out[, first, drop = FALSE], meanSquare = out[, -first, drop = FALSE])
}

#' Diversity of surviving clades through time
#'
#' For clades that survive to `horizon`, the expected number of species
#' \eqn{E(n_t \mid n_T > 0)} and of lineages with descendants at the horizon,
#' at each time \eqn{t} after the origin (Budd & Mann, eqs A.13-A.17).
#'
#' \eqn{E(n_t \mid n_T > 0) = (N(t) - C(t)) / S(T)}, where
#' \eqn{C(t) = E[n_t 1\{n_T = 0\}]} solves the tangent of the
#' generating-function equation along \eqn{G(s) = E(T - t + s)}, a separate
#' solve for each \eqn{t}.
#'
#' @inheritParams ExtinctionProbability
#' @param times Numeric vector giving times since the clade's origin, each no
#' greater than `horizon`.
#' @param horizon Numeric giving the time \eqn{T} at which survival is
#' required.
#' @returns `SurvivorDiversity()` returns a list of matrices `species` and
#' `lineages`, with one row per element of `times` and one column per founding
#' bin.
#' @examples
#' grid <- TempoGrid(41, c(-4, 4))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' SurvivorDiversity(c(50, 100), Q, rates$lambda, rates$mu, 100)$species[, 21]
#' @family dCET numerics
#' @export
SurvivorDiversity <- function(times, Q, lambda, mu, horizon,
                              rtol = 1e-10, atol = 1e-12) {
  if (any(times > horizon)) {
    stop("`times` must not exceed `horizon`")
  }
  nBin <- nrow(Q)
  lambda <- rep_len(lambda, nBin)
  mu <- rep_len(mu, nBin)
  sigma <- lambda + mu
  growth <- .AddDiag(Q, lambda - mu)
  ages <- sort(unique(c(horizon - times, horizon)))
  extinct <- ExtinctionProbability(ages, Q, lambda, mu, rtol = rtol, atol = atol)
  survival <- 1 - extinct[length(ages), ]
  meanN <- DiversityMoments(times, Q, lambda, mu, rtol = rtol, atol = atol)$mean
  first <- seq_len(nBin)
  rows <- lapply(seq_along(times), function(i) {
    t <- times[[i]]
    eStart <- extinct[match(horizon - t, ages), ]
    tangent <- .Integrate(
      c(eStart, eStart), t,
      function(s, y) {
        g <- y[first]
        c(mu - sigma * g + lambda * g * g + Q %*% g,
          Q %*% y[-first] + (2 * lambda * g - sigma) * y[-first])
      },
      function(s, y) {
        g <- y[first]
        jac <- matrix(0, 2 * nBin, 2 * nBin)
        jac[first, first] <- .AddDiag(Q, 2 * lambda * g - sigma)
        jac[-first, -first] <- .AddDiag(Q, 2 * lambda * g - sigma)
        jac[cbind(first + nBin, first)] <- 2 * lambda * y[-first]
        jac
      },
      rtol, atol
    )[1, -first]
    lineages <- .Linear(t, growth, 1 - eStart, rtol, atol)[1, ]
    c((meanN[i, ] - tangent) / survival, lineages / survival)
  })
  out <- do.call(rbind, rows)
  list(species = out[, first, drop = FALSE], lineages = out[, -first, drop = FALSE])
}

#' Tempo distribution of a clade
#'
#' Expected number of species in each tempo bin at time \eqn{t}, from the
#' discrete replicator-mutation equation \eqn{\dot n = (Q^\top + R) n}
#' (Budd & Mann, eq. 5).
#'
#' @inheritParams ExtinctionProbability
#' @param start Integer giving the founding bin, or numeric vector giving the
#' initial expected count in each bin.
#' @param normalize Logical; if `TRUE`, rescale each row to sum to one.
#' @returns `TempoDistribution()` returns a matrix with one row per element of
#' `times` and one column per bin.
#' @examples
#' grid <- TempoGrid(41, c(-4, 4))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' p <- TempoDistribution(100, Q, rates$lambda, rates$mu, start = 21)
#' sum(p * grid$x) # Mean log-tempo
#' @family dCET numerics
#' @export
TempoDistribution <- function(times, Q, lambda, mu, start, normalize = TRUE,
                              rtol = 1e-10, atol = 1e-14) {
  nBin <- nrow(Q)
  if (length(start) == 1L) {
    start <- replace(numeric(nBin), start, 1)
  }
  forward <- .AddDiag(t(Q), rep_len(lambda, nBin) - rep_len(mu, nBin))
  out <- .Linear(times, forward, start, rtol, atol)
  if (normalize) out / rowSums(out) else out
}

#' Equilibrium tempo distribution
#'
#' The normalized Perron eigenvector of \eqn{Q^\top + R}: the long-run tempo
#' distribution of contemporaneous species, whatever the founding tempo
#' (Budd & Mann, eq. A.25).
#'
#' The dominant eigenvalue is the long-run exponential growth rate of mean
#' diversity, and equals \eqn{\sum_k r_k p_k}.
#' For a valid tridiagonal generator, the eigenproblem is symmetrized via
#' detailed balance and solved in log space; with no net diversification the
#' detailed-balance law is returned directly.
#'
#' @inheritParams ExtinctionProbability
#' @returns `TempoEquilibrium()` returns a list with elements `p`
#' (equilibrium distribution) and `growth` (dominant eigenvalue).
#' @examples
#' grid <- TempoGrid(81, c(-8, 8))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' eq <- TempoEquilibrium(Q, rates$lambda, rates$mu)
#' eq$growth
#' sum(eq$p * grid$x)
#' @family dCET numerics
#' @export
TempoEquilibrium <- function(Q, lambda = 0, mu = 0) {
  nBin <- nrow(Q)
  r <- rep_len(lambda, nBin) - rep_len(mu, nBin)
  lower <- seq_len(nBin - 1L)
  up <- Q[cbind(lower, lower + 1L)]
  down <- Q[cbind(lower + 1L, lower)]
  if (all(Q[abs(row(Q) - col(Q)) > 1L] == 0) && all(up > 0) && all(down > 0)) {
    logPi <- c(0, cumsum(log(up) - log(down)))
    if (all(r == 0)) {
      p <- exp(logPi - max(logPi))
      return(list(p = p / sum(p), growth = 0))
    }
    off <- sqrt(up * down)
    sym <- diag(diag(Q) + r)
    sym[cbind(lower, lower + 1L)] <- off
    sym[cbind(lower + 1L, lower)] <- off
    eig <- eigen(sym, symmetric = TRUE)
    logP <- log(abs(eig[["vectors"]][, 1])) + logPi / 2
    p <- exp(logP - max(logP))
    growth <- eig[["values"]][[1]]
  } else {
    eig <- eigen(t(Q) + diag(r))
    dominant <- which.max(Re(eig[["values"]]))
    p <- abs(Re(eig[["vectors"]][, dominant]))
    growth <- Re(eig[["values"]][[dominant]])
  }
  list(p = p / sum(p), growth = growth)
}

#' Ancestral tempo
#'
#' Mean and standard deviation of the log-tempo of the ancestor, a time
#' `times` ago, of a randomly chosen modern species, when lineages at that time
#' have tempo distribution `prior` (Budd & Mann, eq. 7).
#'
#' The ancestor's law is proportional to `prior[k] * N_k(t)`, with \eqn{N}
#' from [DiversityMoments()].
#'
#' @inheritParams ExtinctionProbability
#' @param times Numeric vector giving look-back times.
#' @param x Numeric vector giving the log-tempo of each bin.
#' @param prior Numeric vector giving the tempo distribution at the earlier
#' time, e.g. `TempoEquilibrium()$p`.
#' @returns `AncestralTempo()` returns a matrix with columns `mean` and `sd`,
#' and one row per element of `times`.
#' @examples
#' grid <- TempoGrid(81, c(-8, 8))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' prior <- TempoEquilibrium(Q, rates$lambda, rates$mu)$p
#' AncestralTempo(100, grid$x, Q, rates$lambda, rates$mu, prior)
#' @family dCET numerics
#' @export
AncestralTempo <- function(times, x, Q, lambda, mu, prior,
                           rtol = 1e-10, atol = 1e-12) {
  n <- DiversityMoments(times, Q, lambda, mu, rtol = rtol, atol = atol)$mean
  w <- t(t(n) * prior)
  w <- w / rowSums(w)
  m <- as.vector(w %*% x)
  cbind(mean = m, sd = sqrt(pmax(as.vector(w %*% (x * x)) - m * m, 0)))
}

#' Branch moments and branch-duration density
#'
#' A branch runs from its birth to its first speciation or extinction, which
#' occur at total rate `rateOut` in each bin.
#' `BranchExpectation()` gives \eqn{E_k[\int_0^{T} f(X_t)\,dt]} for a branch
#' starting in bin \eqn{k}: `f = 1` gives the mean duration; `f = exp(x)`
#' gives the mean tempo-weighted length \eqn{w}, which is exactly
#' \eqn{1/(\lambda + \mu)} in every bin when
#' \eqn{\alpha = \beta = \gamma = 1} (Budd & Mann, eq. 12).
#' `BranchDurationDensity()` gives the density of the duration
#' (Budd & Mann, eq. 9).
#'
#' @inheritParams ExtinctionProbability
#' @param rateOut Numeric vector giving the per-bin rate at which a branch
#' ends, e.g. `lambda + mu`.
#' @param f Numeric vector giving the per-bin integrand.
#' @param times Numeric vector giving branch durations.
#' @returns `BranchExpectation()` returns a vector with one entry per
#' starting bin; `BranchDurationDensity()` returns a matrix with one row per
#' element of `times` and one column per starting bin.
#' @examples
#' grid <- TempoGrid(81, c(-8, 8))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' out <- rates$lambda + rates$mu
#' range(BranchExpectation(Q, out, exp(grid$x)))
#' BranchExpectation(Q, out)[c(21, 41, 61)]
#' @family dCET numerics
#' @export
BranchExpectation <- function(Q, rateOut, f = 1) {
  nBin <- nrow(Q)
  as.vector(solve(diag(rep_len(rateOut, nBin)) - Q, rep_len(f, nBin)))
}

#' @rdname BranchExpectation
#' @export
BranchDurationDensity <- function(times, Q, rateOut, rtol = 1e-10,
                                  atol = 1e-14) {
  rateOut <- rep_len(rateOut, nrow(Q))
  .Linear(times, .AddDiag(Q, -rateOut), rateOut, rtol, atol)
}
