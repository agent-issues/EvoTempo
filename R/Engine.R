#' Computation engine
#'
#' EvoTempo computes likelihoods, extinction probabilities and tempo paths
#' either in compiled C++ (`"cpp"`) or in pure R (`"R"`).
#' The R implementation is the reference against which the C++ engine is
#' tested; both give the same results to numerical tolerance, and the path
#' samplers draw the same random numbers in the same order.
#'
#' Functions with an `engine` argument use it for everything that they call.
#' Elsewhere, the engine is set by `options(EvoTempo.engine = )`.
#'
#' @param engine Character string selecting the engine:
#'  - `cpp`: compiled C++, with a stiff Radau IIA solver for the backward
#'    equations and MkPrime kernels for characters;
#'  - `R`: the pure-R reference, with `deSolve::lsode()`.
#' @returns `TempoEngine()` returns the engine that would be used, as a
#' character string.
#' @examples
#' TempoEngine()
#' @family dCET inference
#' @export
TempoEngine <- function(engine = getOption("EvoTempo.engine", "cpp")) {
  match.arg(engine, c("cpp", "R"))
}

.UseCpp <- function() {
  identical(getOption("EvoTempo.engine", "cpp"), "cpp")
}

# Signal an infeasible proposal returned by the C++ engine, else return x
.CheckInfeasible <- function(x) {
  if (!is.null(x[["infeasible"]])) {
    .Infeasible(x[["infeasible"]])
  }
  x
}
