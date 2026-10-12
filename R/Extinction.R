#' Extinction probabilities under dCET
#'
#' Solves the backward equation for \eqn{E_k(\tau)}, the probability that a
#' lineage in tempo bin \eqn{k} at age \eqn{\tau} leaves no sampled
#' descendant:
#' \deqn{\frac{dE_k}{d\tau} = \mu_k - (\lambda_k + \mu_k + \psi_k) E_k +
#' \lambda_k E_k^2 + \sum_j Q_{kj} E_j, \qquad E_k(0) = 1 - \rho.}
#'
#' Fossil samples count as sampled descendants; removal after sampling does
#' not enter \eqn{E}.
#' For the probability of no sampled *extant* descendant, set `psi = 0`.
#' With `rho = 1` and `psi = 0`, \eqn{E_k(t)} is the probability that a clade
#' founded in bin \eqn{k} is extinct after time \eqn{t} (Budd & Mann's
#' \eqn{G_x(t, 0)}).
#'
#' @param times Numeric vector giving ages (or clade durations) at which to
#' report \eqn{E}, in increasing order.
#' @param Q Numeric matrix giving the tempo generator, e.g. from
#' [TempoGenerator()].
#' @param lambda,mu Numeric vectors giving per-bin speciation and extinction
#' rates, e.g. from [BinRates()]; recycled to `nrow(Q)`.
#' @param psi Numeric vector giving per-bin fossil-sampling rates.
#' @param rho Numeric giving the probability of sampling an extant lineage.
#' @param rtol,atol Numerics giving the relative and absolute tolerances of
#' the stiff solver.
#' @returns `ExtinctionProbability()` returns a matrix with one row per
#' element of `times` and one column per tempo bin.
#' @examples
#' grid <- TempoGrid(41, c(-4, 4))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' E <- ExtinctionProbability(c(10, 100), Q, rates$lambda, rates$mu)
#' 1 - E[, 21] # Survival from x = 0
#' @family dCET numerics
#' @export
ExtinctionProbability <- function(times, Q, lambda, mu, psi = 0, rho = 1,
                                  rtol = 1e-10, atol = 1e-14) {
  if (length(rho) != 1L || !isTRUE(rho >= 0 && rho <= 1)) {
    stop("`rho` must lie in [0, 1]")
  }
  nBin <- nrow(Q)
  lambda <- rep_len(lambda, nBin)
  mu <- rep_len(mu, nBin)
  sigma <- lambda + mu + rep_len(psi, nBin)
  .Integrate(
    rep(1 - rho, nBin), times,
    function(t, y) mu - sigma * y + lambda * y * y + as.vector(Q %*% y),
    function(t, y) .AddDiag(Q, 2 * lambda * y - sigma),
    rtol, atol
  )
}

#' Extinction probability under the constant-rate fossilized birth-death
#' process
#'
#' Closed-form \eqn{p_0(t)} of Stadler (2010): the probability that a lineage
#' at age \eqn{t} leaves no sampled descendant.
#'
#' @param t Numeric vector giving ages.
#' @param lambda,mu,psi Non-negative numerics giving speciation, extinction
#' and fossil-sampling rates.
#' @inheritParams ExtinctionProbability
#' @returns `FBDExtinction()` returns a numeric vector of probabilities.
#' @references
#' Stadler T (2010). Sampling-through-time in birth-death trees. Journal of
#' Theoretical Biology 267:396--404. \doi{10.1016/j.jtbi.2010.09.010}
#' @examples
#' FBDExtinction(c(1, 10), lambda = 1, mu = 0.5, psi = 0.1, rho = 0.5)
#' @family dCET numerics
#' @export
FBDExtinction <- function(t, lambda, mu, psi = 0, rho = 1) {
  c1 <- sqrt((lambda - mu - psi)^2 + 4 * lambda * psi)
  if (c1 == 0) {
    return(1 - rho / (1 + rho * lambda * t))
  }
  c2 <- -(lambda - mu - 2 * lambda * rho - psi) / c1
  decay <- exp(-c1 * t) * (1 - c2)
  (lambda + mu + psi + c1 * (decay - (1 + c2)) / (decay + (1 + c2))) /
    (2 * lambda)
}
