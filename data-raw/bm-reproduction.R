# Reproduction of Budd & Mann (2025) numerics: (A) their code, verbatim;
# (B) EvoTempo with BuddMannGenerator(); (C) converged dCET; (D) convergence.
# Run from the package root with EvoTempo, deSolve and ReacTran installed.
# Writes data-raw/bm-results/*.csv; figures via bm-figures.R.
library("EvoTempo")
library("deSolve")
library("ReacTran")
outDir <- "data-raw/bm-results"
dir.create(outDir, showWarnings = FALSE)
lambda <- 0.51
mu <- 0.5
theta <- 0.01
horizon <- 500

# (A) Budd & Mann's code: evaluate definitions only (up to the first `if`)
bmEnv <- new.env()
exprs <- parse("data-raw/bm/BuddMannTempoSI.R")
firstIf <- which(vapply(exprs, function(e) is.call(e) && identical(e[[1]], as.name("if")), TRUE))[1]
for (e in exprs[seq_len(firstIf - 1L)]) eval(e, bmEnv)
BuddMannRun <- function(env) {
  local({
    M <- evalM(r, theta, s)
    Mhighr <- evalM(2 * r, theta, s)
    Mlowr <- evalM(0, theta, s)
    D <- vode(rep(0, lx), seq(0, 500, 1), pdeG, parms = c(lambda, mu, theta, s),
              atol = 1e-15)[, 2:(lx + 1)]
    S <- 1 - D
    meanNfromx <- vode(rep(1, lx), seq(0, 500, 1), pdeM, parms = c(r, theta, s),
                       atol = 1e-15)[, 2:(lx + 1)]
    meanNN2fromx <- vode(rep(1, 2 * lx), seq(0, 500, 1), pdeMM2,
                         parms = c(lambda, mu, theta, s), atol = 1e-15)[, 2:(2 * lx + 1)]
    meanclade <- meanNN2fromx[501, x0idx] / (1 - D[501, x0idx])
    meanexclade <- (meanNN2fromx[501, x0idx + lx] / (1 - D[501, x0idx])) / meanclade
    # Fig 2B at selected times, code as published and with the reaction fixed
    pdeCfixed <- pdeC
    stopifnot(grepl("lambda * C2 * (C2 - 1)", paste(deparse(body(pdeC)), collapse = ""),
                    fixed = TRUE))
    body(pdeCfixed) <- parse(text = sub(
      "lambda * C2 * (C2 - 1)", "lambda * C2 * (C2 - 1 + mu/lambda)",
      paste(deparse(body(pdeC)), collapse = "\n"), fixed = TRUE))[[1]]
    Survivor <- function(i, f) {
      m0 <- c(D[501 - (i - 1), ], D[501 - (i - 1), ] - mu / lambda)
      corr <- vode(m0, c(500 - (i - 1), 500), f, parms = c(lambda, mu, theta, s),
                   atol = 1e-15)[2, 2:(lx + 1)]
      (meanNfromx[i, x0idx] - corr[x0idx]) / (1 - D[501, x0idx])
    }
    fig2bTimes <- c(100, 250, 400)
    fig2bCode <- vapply(fig2bTimes + 1, Survivor, 0, f = pdeC)
    fig2bFixed <- vapply(fig2bTimes + 1, Survivor, 0, f = pdeCfixed)
    lineages <- vapply(fig2bTimes + 1, function(i) {
      sum(M[i, , x0idx] * (1 - D[501 - (i - 1), ])) / S[501, x0idx]
    }, 0)
    Mean <- function(m, it, x0) {
      p <- m[it, , which.min(abs(Grid$x.mid - x0))]
      p <- p / sum(p)
      c(sum(p * x), sqrt(sum(p * x^2) - sum(p * x)^2))
    }
    Px <- M[501, , lx]
    Px <- Px / sum(Px)
    fig4a <- vapply(c(100, 200, 500) + 1, function(i) {
      w <- Px * meanNfromx[i, ]
      sum(w * x) / sum(w)
    }, 0)
    fig5 <- vode(ex$mid * (lambda + mu), seq(0, 20, 0.01), pdeM,
                 parms = c(-(lambda + mu), theta, s), atol = 1e-15)[, 2:(lx + 1)]
    ts <- seq(0, 20, 0.01)
    fig5mean <- vapply(c(-2, 0, 2), function(x0) {
      f <- fig5[, which.min(abs(x - x0))]
      sum(f * ts) / sum(f)
    }, 0)
    idx <- vapply(c(-2, 0, 2), function(x0) which.min(abs(x - x0)), 1L)
    list(
      S0 = S[501, x0idx], meanClade = meanclade, experienced = meanexclade,
      N100 = meanNfromx[101, idx], N500 = meanNfromx[501, idx],
      slope = diff(log(meanNfromx[c(401, 501), x0idx])) / 100,
      fig2bCode = fig2bCode, fig2bFixed = fig2bFixed, lineages = lineages,
      fig3aMinus = vapply(c(10, 50, 100, 200, 500) + 1, function(i) Mean(M, i, -2)[1], 0),
      fig3aPlus = vapply(c(10, 50, 100, 200, 500) + 1, function(i) Mean(M, i, 2)[1], 0),
      fig3aZero = Mean(M, 501, 0),
      fig3bLow = Mean(Mlowr, 501, 0)[1], fig3bHigh = Mean(Mhighr, 501, 0)[1],
      fig4a = fig4a, pxMean = sum(Px * x), fig5mean20 = fig5mean,
      fig5 = fig5[, idx]
    )
  }, envir = env)
}

# (B, C) The same quantities from EvoTempo, for any x-grid and rate matrix
EvoTempoRun <- function(x, Q, start = which.min(abs(x)), truncate5 = FALSE) {
  rates <- BinRates(x, lambda, mu)
  idx <- vapply(c(-2, 0, 2), function(x0) which.min(abs(x - x0)), 1L)
  E <- ExtinctionProbability(horizon, Q, rates$lambda, rates$mu)
  mom <- DiversityMoments(c(100, 400, 500), Q, rates$lambda, rates$mu, second = TRUE)
  S0 <- 1 - E[1, start]
  surv <- SurvivorDiversity(c(100, 250, 400), Q, rates$lambda, rates$mu, horizon)
  MeanFrom <- function(x0, times, r = theta) {
    rr <- BinRates(x, lambda = r, mu = 0)
    p <- TempoDistribution(times, Q, rr$lambda, 0, start = which.min(abs(x - x0)))
    as.vector(p %*% x)
  }
  pZero <- TempoDistribution(500, Q, theta * exp(x), 0, start)
  m0 <- sum(pZero * x)
  top <- TempoDistribution(500, Q, theta * exp(x), 0, start = length(x))
  fig4a <- AncestralTempo(c(100, 200, 500), x, Q, rates$lambda, rates$mu,
                          as.vector(top))[, "mean"]
  out <- rates$lambda + rates$mu
  list(
    S0 = S0, meanClade = mom$mean[3, start] / S0,
    experienced = mom$meanSquare[3, start] / mom$mean[3, start],
    N100 = mom$mean[1, idx], N500 = mom$mean[3, idx],
    slope = diff(log(mom$mean[2:3, start])) / 100,
    survivor = surv$species[, start], lineages = surv$lineages[, start],
    fig3aMinus = MeanFrom(-2, c(10, 50, 100, 200, 500)),
    fig3aPlus = MeanFrom(2, c(10, 50, 100, 200, 500)),
    fig3aZero = c(m0, sqrt(sum(pZero * x^2) - m0^2)),
    fig3bLow = MeanFrom(0, 500, 0), fig3bHigh = MeanFrom(0, 500, 2 * theta),
    fig4a = fig4a, pxMean = sum(top * x),
    eqGrowth = TempoEquilibrium(Q, rates$lambda, rates$mu)$growth,
    durations = BranchExpectation(Q, out)[idx],
    fig5 = BranchDurationDensity(seq(0, 20, 0.01), Q, out)[, idx],
    w = range(BranchExpectation(Q, out, exp(x)))
  )
}

bmCode <- BuddMannRun(bmEnv)
bmGen <- BuddMannGenerator()
bmOurs <- EvoTempoRun(bmGen$x, bmGen$Q)
fine <- TempoGrid(401, c(-10, 10))
converged <- EvoTempoRun(fine$x, TempoGenerator(fine, theta))
saveRDS(list(code = bmCode, ours = bmOurs, converged = converged),
        file.path(outDir, "bm-columns.rds"))
if (exists("skipConvergence") && skipConvergence) stop("Skipping convergence")

# (D) Convergence in h and in z-range
Quantities <- function(x, Q) {
  rates <- BinRates(x, lambda, mu)
  start <- which.min(abs(x))
  E <- ExtinctionProbability(horizon, Q, rates$lambda, rates$mu)[1, start]
  mom <- DiversityMoments(horizon, Q, rates$lambda, rates$mu, second = TRUE)
  top <- TempoDistribution(500, Q, theta * exp(x), 0, start = length(x))
  anc <- AncestralTempo(c(100, 500), x, Q, rates$lambda, rates$mu,
                        as.vector(top))[, "mean"]
  r0 <- TempoEquilibrium(Q)$p
  c(S0 = 1 - E, meanClade = mom$mean[1, start] / (1 - E),
    experienced = mom$meanSquare[1, start] / mom$mean[1, start],
    fig4a100 = anc[[1]], fig4a500 = anc[[2]], a26 = sum(r0 * x),
    growth = TempoEquilibrium(Q, rates$lambda, rates$mu)$growth)
}
rows <- list()
for (scheme in c("central", "upwind", "sqra")) {
  for (zMax in c(5, 10)) {
    for (h in c(0.4, 0.2, 0.1, 0.05)) {
      if (scheme == "central" && h * zMax > 2) next
      grid <- TempoGrid(round(2 * zMax / h) + 1, c(-zMax, zMax))
      q <- Quantities(grid$x, TempoGenerator(grid, theta, scheme = scheme))
      rows[[length(rows) + 1]] <- data.frame(scheme, zMax, h, K = length(grid$x), t(q))
    }
  }
}
for (nCell in c(51, 101, 201, 401, 801)) {
  g <- BuddMannGenerator(nCell)
  q <- Quantities(g$x, g$Q)
  rows[[length(rows) + 1]] <- data.frame(scheme = "budd-mann", zMax = 10, h = g$h,
                                         K = nCell, t(q))
}
for (zMax in c(3, 4, 5, 6, 8)) {
  grid <- TempoGrid(round(2 * zMax / 0.1) + 1, c(-zMax, zMax))
  q <- Quantities(grid$x, TempoGenerator(grid, theta))
  rows[[length(rows) + 1]] <- data.frame(scheme = "central-range", zMax, h = 0.1,
                                         K = length(grid$x), t(q))
}
convergence <- do.call(rbind, rows)
write.csv(convergence, file.path(outDir, "convergence.csv"), row.names = FALSE)
