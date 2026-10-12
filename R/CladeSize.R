#' Clade-size distribution
#'
#' Probability \eqn{P_k(n_t = n)} that a clade founded by one lineage in bin
#' \eqn{k} has \eqn{n} species after time \eqn{t}, for
#' \eqn{n = 0, \dots, n_\mathrm{max}}.
#'
#' The probability generating function \eqn{G_k(t, z)} (Budd & Mann, eq. 2)
#' is solved at `nPoint` points on a circle of radius `radius` < 1, and
#' inverted by FFT (Cauchy's formula).
#' The radius damps aliasing: mass at \eqn{n + m \cdot} `nPoint` contributes
#' at most `radius ^ nPoint` relative error, while solver error is amplified
#' by `radius ^ -n`.  The defaults bound aliasing at \eqn{10^{-10}} with
#' `nPoint` \eqn{\ge 4 (n_\mathrm{max} + 1)}.
#'
#' @inheritParams ExtinctionProbability
#' @param time Numeric giving the clade age.
#' @param start Integer vector giving the founding bins to report.
#' @param nMax Integer giving the largest clade size to report.
#' @param nPoint Integer giving the number of contour points.
#' @param radius Numeric giving the contour radius.
#' @param cores Integer giving the number of parallel worker processes.
#' @returns `CladeSizeDistribution()` returns a matrix with rows
#' \eqn{n = 0, \dots, n_\mathrm{max}} and one column per element of `start`.
#' @examples
#' grid <- TempoGrid(21, c(-3, 3))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' p <- CladeSizeDistribution(10, Q, rates$lambda, rates$mu, start = 11,
#'                            nMax = 20)
#' sum(p) # Close to one: few clades exceed 20 species after 10 Myr
#' @family dCET numerics
#' @importFrom deSolve zvode
#' @importFrom parallel clusterCall makeCluster parLapply stopCluster
#' @importFrom stats fft
#' @export
CladeSizeDistribution <- function(time, Q, lambda, mu, start, nMax,
                                  nPoint = 2^ceiling(log2(4 * (nMax + 1))),
                                  radius = 1e-10^(1 / nPoint), cores = 1L,
                                  rtol = 1e-10, atol = 1e-14) {
  if (nPoint <= nMax) {
    stop("`nPoint` must exceed `nMax`")
  }
  nBin <- nrow(Q)
  lambda <- rep_len(lambda, nBin)
  mu <- rep_len(mu, nBin)
  sigma <- lambda + mu
  Solve <- function(z0) {
    deSolve::zvode(
      rep(z0, nBin), c(0, time),
      function(t, y, p) list(mu - sigma * y + lambda * y * y + as.vector(Q %*% y)),
      parms = NULL,
      jacfunc = function(t, y, p) Q + diag(2 * lambda * y - sigma),
      jactype = "fullusr", rtol = rtol, atol = atol, maxsteps = 1e6
    )[2, -1L][start]
  }
  z <- radius * exp(2i * pi * (seq_len(nPoint) - 1) / nPoint)
  values <- if (cores > 1L) {
    cl <- makeCluster(cores)
    on.exit(stopCluster(cl))
    SetPaths <- function(paths) .libPaths(paths)
    environment(SetPaths) <- globalenv()
    clusterCall(cl, SetPaths, .libPaths())
    parLapply(cl, z, Solve)
  } else {
    lapply(z, Solve)
  }
  values <- matrix(unlist(values), nPoint, length(start), byrow = TRUE)
  n <- 0:nMax
  apply(values, 2, function(g) Re(fft(g))[n + 1L]) / nPoint / radius^n
}
