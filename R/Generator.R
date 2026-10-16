#' Log-tempo grid
#'
#' Equally spaced log-tempo bins in standardized units \eqn{z = x / s}.
#'
#' Bins are equally spaced in \eqn{z}, so the Ornstein-Uhlenbeck block of the
#' generator does not depend on \eqn{s}, and a path stored as bin indices
#' stays valid when \eqn{s} changes.
#'
#' @param nBin Integer giving the number of bins, \eqn{K}.
#' @param zRange Numeric vector giving the lowest and highest bin centres, in
#' units of `s`.
#' @param s Positive numeric giving the stationary standard deviation of
#' log-tempo.
#' @returns `TempoGrid()` returns an object of class `TempoGrid`: a list with
#' elements `z` (bin centres in units of `s`), `x` (log-tempo, `s * z`), `h`
#' (spacing in `z`) and `s`.
#' @examples
#' grid <- TempoGrid(51, c(-5, 5))
#' range(grid$x)
#' @family dCET numerics
#' @export
TempoGrid <- function(nBin, zRange = c(-5, 5), s = 1) {
  nBin <- as.integer(nBin)
  if (length(nBin) != 1L || is.na(nBin) || nBin < 2L) {
    stop("`nBin` must be a single integer >= 2")
  }
  if (length(zRange) != 2L || !isTRUE(zRange[[2]] > zRange[[1]])) {
    stop("`zRange` must be an increasing pair")
  }
  if (length(s) != 1L || !isTRUE(s > 0)) {
    stop("`s` must be a single positive number")
  }
  z <- seq(zRange[[1]], zRange[[2]], length.out = nBin)
  structure(
    list(z = z, x = s * z, h = (zRange[[2]] - zRange[[1]]) / (nBin - 1L), s = s),
    class = "TempoGrid"
  )
}

#' Discretized tempo generator
#'
#' Constructs the generator \eqn{Q = \mathrm{diag}(e^{\alpha x_k}) A} of
#' discretized covariant tempo, where the tridiagonal \eqn{A} discretizes the
#' Ornstein-Uhlenbeck operator \eqn{-\theta z f' + \theta f''} with reflecting
#' ends.
#'
#' Central differences are second order, but have non-negative rates only if
#' \eqn{h \max_k |z_k| \le 2}; this is enforced.
#' Upwind differences are always valid, but first order, adding numerical
#' diffusion \eqn{\theta |z| h / 2}.
#' The square-root approximation is second order and always valid, and its
#' stationary law is exactly the Gaussian on the grid.
#' Central and upwind schemes satisfy \eqn{(A z)_k = -\theta z_k} exactly in
#' the interior.
#'
#' @param grid A [`TempoGrid`][TempoGrid()].
#' @param theta Positive numeric giving the mean-reversion rate per unit
#' tempo-time.
#' @param alpha Numeric giving the tempo self-rate exponent.
#' @param scheme Character string specifying the discretization:
#'  - `central`: central differences;
#'  - `upwind`: first-order upwind differences;
#'  - `sqra`: square-root approximation.
#' @returns `TempoGenerator()` returns a \eqn{K \times K} generator matrix:
#' non-negative off-diagonal entries; rows summing to zero.
#' @references
#' Budd GE, Mann RP (2025). Evolutionary tempo, supertaxa, and living fossils.
#' Systematic Biology 74:917--934. \doi{10.1093/sysbio/syaf020}
#' @examples
#' Q <- TempoGenerator(TempoGrid(21, c(-4, 4)), theta = 0.01)
#' range(rowSums(Q))
#' @family dCET numerics
#' @export
TempoGenerator <- function(grid, theta, alpha = 1,
                           scheme = c("central", "upwind", "sqra")) {
  scheme <- match.arg(scheme)
  if (!inherits(grid, "TempoGrid")) {
    stop("`grid` must be a TempoGrid")
  }
  if (length(theta) != 1L || !isTRUE(theta > 0)) {
    stop("`theta` must be a single positive number")
  }
  .OUBlock(grid[["z"]], grid[["h"]], theta, scheme) * exp(alpha * grid[["x"]])
}

.OUBlock <- function(z, h, theta, scheme) {
  nBin <- length(z)
  base <- theta / (h * h)
  rates <- switch(
    scheme,
    central = {
      if (h * max(abs(z)) > 2 * (1 + 1e-12)) {
        stop("Central scheme needs h * max|z| <= 2 (here ",
             signif(h * max(abs(z)), 4), "); increase `nBin`")
      }
      list(up = base - theta * z / (2 * h), down = base + theta * z / (2 * h))
    },
    upwind = list(up = base + theta * pmax(-z, 0) / h,
                  down = base + theta * pmax(z, 0) / h),
    sqra = {
      u <- z * z / 2
      list(up = base * exp(-(c(u[-1], 0) - u) / 2),
           down = base * exp(-(c(0, u[-nBin]) - u) / 2))
    }
  )
  .Tridiagonal(rates[["up"]], rates[["down"]])
}

# Up-rates to k + 1 and down-rates to k - 1; reflecting ends.
.Tridiagonal <- function(up, down) {
  nBin <- length(up)
  up[nBin] <- 0
  down[1] <- 0
  ret <- matrix(0, nBin, nBin)
  lower <- seq_len(nBin - 1L)
  ret[cbind(lower, lower + 1L)] <- up[lower]
  ret[cbind(lower + 1L, lower)] <- down[-1]
  diag(ret) <- -(up + down)
  ret
}

#' Budd & Mann's discretization
#'
#' Reproduces the method-of-lines operator of Budd & Mann's supplementary
#' code: ReacTran `tran.1D()` on `nCell` finite-volume cells, with centred
#' second differences for diffusion and a backward difference for the drift
#' at every \eqn{x}.
#'
#' The fixed drift difference is upwind for \eqn{x > 0} and downwind for
#' \eqn{x < 0}, giving effective diffusion \eqn{\theta(s^2 + x h / 2)}.
#' Down-rates are negative for \eqn{x < -s^2 / h}, so `Q` is not a valid
#' generator there.  Use it only to reproduce published values.
#'
#' @inheritParams TempoGenerator
#' @inheritParams TempoGrid
#' @param nCell Integer giving the number of grid cells.
#' @param xRange Numeric vector giving the limits of the log-tempo domain.
#' @returns `BuddMannGenerator()` returns a list with elements `x` (cell
#' centres), `h` (cell width) and `Q` (rate matrix with zero row sums).
#' @examples
#' bm <- BuddMannGenerator()
#' bm$h
#' @family dCET numerics
#' @export
BuddMannGenerator <- function(nCell = 101, xRange = c(-10, 10),
                              theta = 0.01, s = 1, alpha = 1) {
  h <- (xRange[[2]] - xRange[[1]]) / nCell
  x <- xRange[[1]] + h * (seq_len(nCell) - 0.5)
  diffusion <- theta * s * s / (h * h)
  q <- .Tridiagonal(rep(diffusion, nCell), diffusion + theta * x / h)
  list(x = x, h = h, Q = q * exp(alpha * x))
}

#' Per-bin rates
#'
#' @param grid A [`TempoGrid`][TempoGrid()], or numeric vector of log-tempos.
#' @param lambda,mu,psi Non-negative numerics giving speciation, extinction
#' and fossil-sampling rates at unit tempo.
#' @param beta,delta Numerics giving the coupling of diversification and of
#' fossil sampling to tempo.
#' @returns `BinRates()` returns a list of vectors `lambda`, `mu` and `psi`,
#' with \eqn{\lambda_k = \lambda e^{\beta x_k}},
#' \eqn{\mu_k = \mu e^{\beta x_k}} and \eqn{\psi_k = \psi e^{\delta x_k}}.
#' @examples
#' BinRates(TempoGrid(5, c(-2, 2)), lambda = 0.51, mu = 0.5)
#' @family dCET numerics
#' @export
BinRates <- function(grid, lambda, mu, psi = 0, beta = 1, delta = 0) {
  x <- if (inherits(grid, "TempoGrid")) grid[["x"]] else grid
  list(lambda = lambda * exp(beta * x), mu = mu * exp(beta * x),
       psi = psi * exp(delta * x))
}
