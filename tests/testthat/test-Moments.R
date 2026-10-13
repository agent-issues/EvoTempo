coreGrid <- function(nBin = 161, zMax = 8) TempoGrid(nBin, c(-zMax, zMax))

test_that("Diversity moments match the constant-rate birth-death process", {
  Q <- TempoGenerator(TempoGrid(11, c(-2, 2)), 0.1)
  t <- c(0, 10, 100)
  lambda <- 0.6
  mu <- 0.4
  r <- lambda - mu
  mom <- DiversityMoments(t, Q, lambda, mu, second = TRUE)
  expect_equal(mom$mean, matrix(exp(r * t), 3, 11), tolerance = 1e-9)
  m2 <- 2 * lambda / r * exp(2 * r * t) + (1 - 2 * lambda / r) * exp(r * t)
  expect_equal(mom$meanSquare, matrix(m2, 3, 11), tolerance = 1e-7)
  expect_named(DiversityMoments(t, Q, lambda, mu), "mean")
})

test_that("Equilibrium growth is the Perron eigenvalue and the log-slope", {
  grid <- coreGrid()
  Q <- TempoGenerator(grid, 0.01)
  rates <- BinRates(grid, 0.51, 0.5)
  eq <- TempoEquilibrium(Q, rates$lambda, rates$mu)
  expect_equal(eq$growth, sum((rates$lambda - rates$mu) * eq$p))
  expect_equal(eq$growth, 0.014332, tolerance = 1e-4)
  n <- DiversityMoments(c(400, 500), Q, rates$lambda, rates$mu)$mean[, 81]
  expect_lt(abs(diff(log(n)) / 100 - eq$growth), 1e-6)
  # Equilibrium is stationary under the normalized replicator dynamics
  p <- TempoDistribution(50, Q, rates$lambda, rates$mu, start = eq$p)
  expect_equal(as.vector(p), eq$p, tolerance = 1e-7)
})

test_that("Equilibrium mean is 0 when r = theta (corrected A.28)", {
  grid <- coreGrid()
  for (scheme in c("central", "upwind")) {
    Q <- TempoGenerator(grid, 0.01, scheme = scheme)
    p <- TempoEquilibrium(Q, 0.51 * exp(grid$x), 0.5 * exp(grid$x))$p
    expect_equal(sum(p * grid$x), 0, tolerance = 1e-10)
  }
  # Corrected A.28 at r != theta: (r - theta)<x e^x> - r <x><e^x> = 0
  Q <- TempoGenerator(grid, 0.01)
  r <- 0.02
  p <- TempoEquilibrium(Q, r * exp(grid$x))$p
  ex <- exp(grid$x)
  expect_equal((r - 0.01) * sum(p * grid$x * ex) -
                 r * sum(p * grid$x) * sum(p * ex), 0, tolerance = 1e-10)
  expect_equal(sum(p * grid$x), 0.634, tolerance = 1e-3)
})

test_that("Equilibrium mean is -alpha s^2 when r = 0 (A.26)", {
  for (s in c(0.5, 1, 2)) {
    grid <- TempoGrid(201, c(-10, 10), s = s)
    Q <- TempoGenerator(grid, 0.01, scheme = "sqra")
    p <- TempoEquilibrium(Q)$p
    expect_equal(sum(p * grid$x), -s^2, tolerance = 1e-8)
  }
  grid <- TempoGrid(161, c(-8, 8))
  p <- TempoEquilibrium(TempoGenerator(grid, 0.01, alpha = 0.5,
                                       scheme = "sqra"))$p
  expect_equal(sum(p * grid$x), -0.5, tolerance = 1e-8)
})

test_that("TempoEquilibrium() handles non-generator rate matrices", {
  bm <- BuddMannGenerator()
  r <- 0.01 * exp(bm$x)
  eq <- TempoEquilibrium(bm$Q, r)
  expect_true(all(eq$p > 0))
  expect_equal(eq$growth, sum(r * eq$p), tolerance = 1e-8)
})

test_that("Mean log-tempo from x0 = 0 is exactly 0 when r = theta", {
  grid <- coreGrid()
  Q <- TempoGenerator(grid, 0.01)
  rates <- BinRates(grid, 0.51, 0.5)
  p <- TempoDistribution(c(10, 100, 200), Q, rates$lambda, rates$mu,
                         start = 81)
  expect_equal(as.vector(p %*% grid$x), numeric(3), tolerance = 1e-9)
  # Forward counts sum to the backward mean diversity (duality)
  n <- TempoDistribution(100, Q, rates$lambda, rates$mu, start = 61,
                         normalize = FALSE)
  expect_equal(sum(n), DiversityMoments(100, Q, rates$lambda,
                                        rates$mu)$mean[1, 61],
               tolerance = 1e-8)
})

test_that("Branch w is Exp(lambda + mu) in every bin (eq. 12)", {
  grid <- coreGrid()
  for (scheme in c("central", "upwind", "sqra")) {
    Q <- TempoGenerator(grid, 0.01, scheme = scheme)
    rates <- BinRates(grid, 0.51, 0.5)
    w <- BranchExpectation(Q, rates$lambda + rates$mu, exp(grid$x))
    expect_equal(w, rep(1 / 1.01, 161), tolerance = 1e-10)
  }
  Q <- TempoGenerator(grid, 0.01)
  durations <- BranchExpectation(Q, rates$lambda + rates$mu)
  expect_equal(durations[c(61, 81, 101)], c(7.244, 0.9998, 0.138),
               tolerance = 2e-3)
})

test_that("Branch-duration density integrates to one, with the right mean", {
  grid <- TempoGrid(41, c(-4, 4))
  Q <- TempoGenerator(grid, 0.01)
  rates <- BinRates(grid, 0.51, 0.5)
  out <- rates$lambda + rates$mu
  t <- seq(0, 60, by = 0.002)
  f <- BranchDurationDensity(t, Q, out)[, c(21, 26)]
  Trapz <- function(y) sum((y[-1] + y[-length(y)]) / 2) * 0.002
  expect_equal(apply(f, 2, Trapz), c(1, 1), tolerance = 1e-5)
  expect_equal(apply(f * t, 2, Trapz),
               BranchExpectation(Q, out)[c(21, 26)], tolerance = 1e-4)
})

test_that("SurvivorDiversity() matches the birth-death generating function", {
  Q <- TempoGenerator(TempoGrid(7, c(-2, 2)), 0.1)
  lambda <- 0.6
  mu <- 0.4
  horizon <- 20
  t <- c(0, 5, 12, 20)
  out <- SurvivorDiversity(t, Q, lambda, mu, horizon)
  Pgf <- function(time, z) {
    decay <- exp(-(lambda - mu) * time)
    (mu * (1 - z) - (mu - lambda * z) * decay) /
      (lambda * (1 - z) - (mu - lambda * z) * decay)
  }
  survival <- 1 - Pgf(horizon, 0)
  eps <- 1e-6
  expected <- vapply(t, function(time) {
    q <- Pgf(horizon - time, 0)
    deriv <- (Pgf(time, q + eps) - Pgf(time, q - eps)) / (2 * eps)
    (exp((lambda - mu) * time) - q * deriv) / survival
  }, 0)
  expect_equal(out$species[, 4], expected, tolerance = 1e-6)
  expect_equal(out$lineages[c(1, 4), 4],
               c(1, exp((lambda - mu) * horizon) / survival),
               tolerance = 1e-8)
  expect_error(SurvivorDiversity(30, Q, lambda, mu, horizon), "horizon")
})

test_that("AncestralTempo() weights the prior by expected descendants", {
  grid <- TempoGrid(41, c(-4, 4))
  Q <- TempoGenerator(grid, 0.01)
  rates <- BinRates(grid, 0.51, 0.5)
  oneHot <- replace(numeric(41), 25, 1)
  anc <- AncestralTempo(c(0, 50), grid$x, Q, rates$lambda, rates$mu, oneHot)
  expect_equal(anc[, "mean"], rep(grid$x[25], 2))
  expect_equal(anc[, "sd"], c(0, 0))
  prior <- TempoEquilibrium(Q, rates$lambda, rates$mu)$p
  anc <- AncestralTempo(c(0, 100), grid$x, Q, rates$lambda, rates$mu, prior)
  expect_equal(anc[1, "mean"], sum(prior * grid$x), ignore_attr = TRUE)
  expect_gt(anc[2, "mean"], anc[1, "mean"])
})
