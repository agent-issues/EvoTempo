#' Simulate a clade under discretized tempo
#'
#' Exact (Gillespie) forward simulation of the dCET process from a single
#' lineage, for time `time`.
#' Each lineage jumps between tempo bins at the rates in `Q`, speciates,
#' goes extinct, and is sampled as a fossil, all at its current bin's rates.
#' A sampled fossil lineage is removed with probability `removal`; extant
#' lineages are sampled with probability `rho`.
#'
#' Lineages evolve independently, so each is simulated in turn; the result is
#' exact for any `Q`.
#'
#' @inheritParams ExtinctionProbability
#' @param time Numeric giving the duration of the simulation (the origin age).
#' @param x Numeric vector giving the log-tempo of each bin.
#' @param Q Numeric matrix giving a tridiagonal tempo generator, e.g. from
#' [TempoGenerator()].
#' @param start Integer giving the founding lineage's bin.
#' @param removal Numeric giving the probability that a sampled fossil
#' lineage is removed.
#' @param maxEdges Integer giving the largest number of branches to simulate
#' before stopping with an error.
#' @returns `SimulateDCET()` returns an object of class `TempoHistory`: a
#' list with elements
#' - `edges`: data frame with one row per branch of the complete tree (from
#'   birth to speciation, extinction, removal or the present), giving
#'   `parent` branch, `start` and `end` times, end `type`, and whether an
#'   extant lineage is `sampled`;
#' - `paths`: list giving, for each branch, the tempo path as knots `t`
#'   (times) and `x` (log-tempo from each knot to the next);
#' - `fossils`: data frame of fossil samples (`edge`, `time`, `x`);
#' - `time`: the simulation duration.
#' @seealso [ReconstructedTree()], [SimulateCET()], [BranchClock()].
#' @examples
#' grid <- TempoGrid(41, c(-4, 4))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5, psi = 0.05)
#' set.seed(1)
#' history <- SimulateDCET(50, grid$x, Q, rates$lambda, rates$mu, rates$psi,
#'                         start = 21)
#' head(history$edges)
#' @family simulation
#' @importFrom stats rbinom rexp runif
#' @export
SimulateDCET <- function(time, x, Q, lambda, mu, psi = 0, start, rho = 1,
                         removal = 0, maxEdges = 1e5) {
  nBin <- nrow(Q)
  lower <- seq_len(nBin - 1L)
  up <- c(Q[cbind(lower, lower + 1L)], 0)
  down <- c(0, Q[cbind(lower + 1L, lower)])
  if (any(Q[abs(row(Q) - col(Q)) > 1L] != 0) || any(up < 0) || any(down < 0)) {
    stop("`Q` must be a tridiagonal generator with non-negative rates")
  }
  rates <- cbind(down, up, rep_len(lambda, nBin), rep_len(mu, nBin),
                 rep_len(psi, nBin))
  total <- rowSums(rates)
  cumProb <- t(apply(rates, 1, cumsum)) / total
  capacity <- 64L
  parent <- rep(NA_integer_, capacity)
  startTime <- endTime <- numeric(capacity)
  startBin <- endType <- integer(capacity)
  startBin[1] <- start
  paths <- vector("list", capacity)
  fossils <- list()
  nEdge <- 1L
  e <- 0L
  while (e < nEdge) {
    e <- e + 1L
    t <- startTime[[e]]
    k <- startBin[[e]]
    knotT <- t
    knotK <- k
    repeat {
      t <- t + if (total[[k]] > 0) rexp(1, total[[k]]) else Inf
      if (t >= time) {
        t <- time
        type <- 1L
        break
      }
      event <- which(runif(1) < cumProb[k, ])[1]
      if (event < 3L) {
        k <- k + 2L * event - 3L
        knotT <- c(knotT, t)
        knotK <- c(knotK, k)
        next
      }
      if (event == 5L) {
        fossils[[length(fossils) + 1L]] <- c(e, t, x[[k]])
        if (runif(1) < removal) {
          type <- 4L
          break
        }
        next
      }
      type <- event - 1L # 2: speciation; 3: extinction
      break
    }
    endTime[e] <- t
    endType[e] <- type
    paths[[e]] <- list(t = knotT, x = x[knotK])
    if (type == 2L) {
      .CheckEdges(nEdge, maxEdges)
      if (nEdge + 2L > capacity) {
        capacity <- 2L * capacity
        length(parent) <- length(startTime) <- length(endTime) <- capacity
        length(startBin) <- length(endType) <- length(paths) <- capacity
      }
      daughters <- nEdge + 1:2
      parent[daughters] <- e
      startTime[daughters] <- t
      startBin[daughters] <- k
      nEdge <- nEdge + 2L
    }
  }
  .TempoHistory(parent, startTime, endTime, endType, paths, nEdge, fossils, rho,
                time, "discrete")
}

.CheckEdges <- function(nEdge, maxEdges) {
  if (nEdge + 2L > maxEdges) {
    stop("More than `maxEdges` branches; reduce `time` or raise `maxEdges`")
  }
}

.TempoHistory <- function(parent, startTime, endTime, endType, paths, nEdge,
                          fossils, rho, time, kind) {
  used <- seq_len(nEdge)
  types <- c("present", "speciation", "extinction", "removed")[endType[used]]
  present <- types == "present"
  sampled <- rep(NA, length(types))
  sampled[present] <- rbinom(sum(present), 1, rho) == 1
  fossils <- if (length(fossils)) {
    f <- do.call(rbind, fossils)
    data.frame(edge = as.integer(f[, 1]), time = f[, 2], x = f[, 3])
  } else {
    data.frame(edge = integer(0), time = numeric(0), x = numeric(0))
  }
  structure(
    list(
      edges = data.frame(parent = parent[used], start = startTime[used],
                         end = endTime[used], type = types, sampled = sampled),
      paths = paths[used], fossils = fossils, time = time
    ),
    class = "TempoHistory", kind = kind
  )
}

#' Simulate a clade under continuous covariant tempo
#'
#' Forward simulation of Budd & Mann's continuous CET model
#' (\eqn{\alpha = \beta = 1}), exact in tempo-time \eqn{w = \int e^x dt}.
#'
#' In tempo-time, speciation and extinction occur at constant rates
#' \eqn{\lambda} and \eqn{\mu}, and log-tempo is an Ornstein-Uhlenbeck process
#' with exact transitions on steps of `dw`.
#' Real time, \eqn{\int e^{-x} dw}, is integrated by the trapezoid rule on
#' those steps; it places fossils (when `delta = 0`) and the present.
#' With `delta = 1`, fossil sampling has constant rate in tempo-time too.
#'
#' @inheritParams SimulateDCET
#' @inheritParams TempoGenerator
#' @param s Positive numeric giving the stationary standard deviation of
#' log-tempo.
#' @param lambda,mu,psi Non-negative numerics giving speciation, extinction
#' and fossil-sampling rates at unit tempo.
#' @param delta Integer specifying the coupling of fossil sampling:
#'  - `0`: constant rate `psi` in real time;
#'  - `1`: rate `psi * exp(x)` in real time.
#' @param x0 Numeric giving the founding lineage's log-tempo.
#' @param dw Positive numeric giving the tempo-time step of the path.
#' @returns `SimulateCET()` returns a `TempoHistory` (see [SimulateDCET()]),
#' whose paths are piecewise linear between knots, with an extra element `w`
#' (tempo-time at each knot).
#' @examples
#' set.seed(1)
#' history <- SimulateCET(50, theta = 0.01, s = 1, lambda = 0.51, mu = 0.5)
#' table(history$edges$type)
#' @family simulation
#' @importFrom stats filter rnorm rpois
#' @export
SimulateCET <- function(time, theta, s, lambda, mu, psi = 0, delta = 0,
                        x0 = 0, rho = 1, removal = 0, dw = 0.01,
                        maxEdges = 1e5) {
  if (!delta %in% c(0, 1)) {
    stop("`delta` must be 0 or 1")
  }
  wRate <- lambda + mu + delta * psi
  decay <- exp(-theta * dw)
  noise <- s * sqrt(1 - decay^2)
  capacity <- 64L
  parent <- rep(NA_integer_, capacity)
  startTime <- endTime <- startX <- numeric(capacity)
  endType <- integer(capacity)
  startX[1] <- x0
  paths <- vector("list", capacity)
  fossils <- list()
  nEdge <- 1L
  e <- 0L
  while (e < nEdge) {
    e <- e + 1L
    knotT <- startTime[[e]]
    knotX <- startX[[e]]
    knotW <- 0
    repeat {
      seg <- .OUSegment(rexp(1, wRate), knotX[[length(knotX)]], theta, s, dw,
                        decay, noise)
      segT <- knotT[[length(knotT)]] + c(0, cumsum(seg$dt))
      segW <- knotW[[length(knotW)]] + c(0, cumsum(seg$dw))
      cross <- which(segT >= time)[1]
      if (!is.na(cross)) {
        f <- (time - segT[cross - 1L]) / (segT[cross] - segT[cross - 1L])
        keep <- seq_len(cross - 1L)[-1]
        knotT <- c(knotT, segT[keep], time)
        knotX <- c(knotX, seg$x[keep],
                   seg$x[cross - 1L] + f * (seg$x[cross] - seg$x[cross - 1L]))
        knotW <- c(knotW, segW[keep],
                   segW[cross - 1L] + f * (segW[cross] - segW[cross - 1L]))
        type <- 1L
        break
      }
      knotT <- c(knotT, segT[-1])
      knotX <- c(knotX, seg$x[-1])
      knotW <- c(knotW, segW[-1])
      u <- runif(1) * wRate
      if (u < lambda) {
        type <- 2L
        break
      }
      if (u < lambda + mu) {
        type <- 3L
        break
      }
      fossils[[length(fossils) + 1L]] <- c(e, knotT[[length(knotT)]],
                                           knotX[[length(knotX)]])
      if (runif(1) < removal) {
        type <- 4L
        break
      }
    }
    if (delta == 0 && psi > 0) {
      span <- c(knotT[[1]], knotT[[length(knotT)]])
      fTimes <- sort(runif(rpois(1, psi * diff(span)), span[1], span[2]))
      for (ft in fTimes) {
        fx <- .Interpolate(knotT, knotX, ft)
        fossils[[length(fossils) + 1L]] <- c(e, ft, fx)
        if (runif(1) < removal) {
          keep <- knotT < ft
          knotW <- c(knotW[keep], .Interpolate(knotT, knotW, ft))
          knotX <- c(knotX[keep], fx)
          knotT <- c(knotT[keep], ft)
          type <- 4L
          break
        }
      }
    }
    endTime[e] <- knotT[[length(knotT)]]
    endType[e] <- type
    paths[[e]] <- list(t = knotT, x = knotX, w = knotW)
    if (type == 2L) {
      .CheckEdges(nEdge, maxEdges)
      if (nEdge + 2L > capacity) {
        capacity <- 2L * capacity
        length(parent) <- length(startTime) <- length(endTime) <- capacity
        length(startX) <- length(endType) <- length(paths) <- capacity
      }
      daughters <- nEdge + 1:2
      parent[daughters] <- e
      startTime[daughters] <- endTime[[e]]
      startX[daughters] <- knotX[[length(knotX)]]
      nEdge <- nEdge + 2L
    }
  }
  .TempoHistory(parent, startTime, endTime, endType, paths, nEdge, fossils, rho,
                time, "continuous")
}

# Exact OU over tempo-time `width` from `x0`, in steps of `dw`.
.OUSegment <- function(width, x0, theta, s, dw, decay, noise) {
  nStep <- ceiling(width / dw)
  steps <- c(rep(dw, nStep - 1L), width - dw * (nStep - 1L))
  xs <- if (nStep > 1L) {
    c(x0, as.vector(filter(noise * rnorm(nStep - 1L), decay, "recursive",
                           init = x0)))
  } else {
    x0
  }
  last <- xs[[nStep]]
  lastDecay <- exp(-theta * steps[[nStep]])
  xs <- c(xs, last * lastDecay + s * sqrt(1 - lastDecay^2) * rnorm(1))
  ex <- exp(-xs)
  list(x = xs, dw = steps, dt = steps * (ex[-1] + ex[-length(ex)]) / 2)
}

.Interpolate <- function(t, y, at) {
  i <- findInterval(at, t, rightmost.closed = TRUE)
  if (i >= length(t)) {
    return(y[[length(y)]])
  }
  y[[i]] + (at - t[[i]]) / (t[[i + 1L]] - t[[i]]) * (y[[i + 1L]] - y[[i]])
}

#' Character clock along a branch
#'
#' \eqn{c \int_a^b e^{\gamma x(t)} dt} along a branch of a simulated
#' history: the expected number of character changes for a character of rate
#' `rate` at unit tempo.
#' With `gamma = 1` and `rate = 1` this is the tempo-weighted length \eqn{w}.
#'
#' Exact for discretized histories; for continuous histories, exact when
#' `gamma` is 0 or 1, and otherwise integrated by the trapezoid rule.
#'
#' @param history A `TempoHistory`, from [SimulateDCET()] or [SimulateCET()].
#' @param edge Integer giving the branch.
#' @param from,to Numerics giving the time interval; `NULL` uses the whole
#' branch.
#' @param gamma Numeric giving the coupling of character rate to tempo.
#' @param rate Numeric giving the character rate at unit tempo, \eqn{c}.
#' @returns `BranchClock()` returns a numeric.
#' @examples
#' set.seed(1)
#' history <- SimulateCET(20, theta = 0.01, s = 1, lambda = 0.51, mu = 0.5)
#' BranchClock(history, 1)
#' @family simulation
#' @export
BranchClock <- function(history, edge, from = NULL, to = NULL, gamma = 1,
                        rate = 1) {
  path <- history$paths[[edge]]
  t <- path$t
  end <- history$edges$end[[edge]]
  from <- if (is.null(from)) t[[1]] else from
  to <- if (is.null(to)) end else to
  if (identical(attr(history, "kind"), "discrete")) {
    breaks <- c(t, end)
    lo <- pmax(breaks[-length(breaks)], from)
    hi <- pmin(breaks[-1], to)
    return(rate * sum(pmax(hi - lo, 0) * exp(gamma * path$x)))
  }
  cumulative <- if (gamma == 1) {
    path$w
  } else {
    ey <- exp(gamma * path$x)
    c(0, cumsum(diff(t) * (ey[-1] + ey[-length(ey)]) / 2))
  }
  rate * (.Interpolate(t, cumulative, to) - .Interpolate(t, cumulative, from))
}
