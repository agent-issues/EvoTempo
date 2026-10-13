# Stiff integration from t = 0; rows of the result correspond to `times`.
#' @importFrom deSolve lsode
.Integrate <- function(y0, times, Func, Jac, rtol, atol) {
  if (any(times < 0) || is.unsorted(times)) {
    stop("`times` must be non-negative and sorted")
  }
  tOut <- unique(c(0, times))
  out <- if (length(tOut) == 1L) {
    matrix(y0, 1L)
  } else {
    sol <- lsode(y0, tOut, function(t, y, p) list(Func(t, y)), parms = NULL,
                 jacfunc = function(t, y, p) Jac(t, y), jactype = "fullusr",
                 rtol = rtol, atol = atol, maxsteps = 1e6)
    if (attr(sol, "istate")[[1]] < 0) {
      stop("ODE solver failed")  # nocov
    }
    unname(unclass(sol)[, -1L, drop = FALSE])
  }
  out[match(times, tOut), , drop = FALSE]
}

.AddDiag <- function(m, d) {
  diag(m) <- diag(m) + d
  m
}

.Linear <- function(times, m, y0, rtol, atol) {
  .Integrate(y0, times, function(t, y) as.vector(m %*% y),
             function(t, y) m, rtol, atol)
}
