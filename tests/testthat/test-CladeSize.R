test_that("CladeSizeDistribution() matches the birth-death geometric law", {
  Q <- TempoGenerator(TempoGrid(5, c(-1, 1)), 0.1)
  lambda <- 0.51
  mu <- 0.5
  time <- 5
  e <- exp((lambda - mu) * time)
  a <- mu * (e - 1) / (lambda * e - mu)
  b <- lambda * (e - 1) / (lambda * e - mu)
  expected <- c(a, (1 - a) * (1 - b) * b^(0:29))
  p <- CladeSizeDistribution(time, Q, lambda, mu, start = c(1, 3), nMax = 30)
  expect_equal(dim(p), c(31, 2))
  expect_equal(p[, 2], expected, tolerance = 1e-8)
  expect_error(CladeSizeDistribution(time, Q, lambda, mu, 1, 30, nPoint = 16),
               "nPoint")
})

test_that("CladeSizeDistribution() agrees with moments and runs in parallel", {
  grid <- TempoGrid(21, c(-3, 3))
  Q <- TempoGenerator(grid, 0.05)
  rates <- BinRates(grid, 0.51, 0.5)
  p <- CladeSizeDistribution(10, Q, rates$lambda, rates$mu, start = 11,
                             nMax = 200, cores = 2L)
  n <- 0:200
  extinct <- ExtinctionProbability(10, Q, rates$lambda, rates$mu)[1, 11]
  mom <- DiversityMoments(10, Q, rates$lambda, rates$mu, second = TRUE)
  expect_equal(p[1], extinct, tolerance = 1e-8)
  expect_equal(sum(n * p), mom$mean[1, 11], tolerance = 1e-4)
  expect_equal(sum(n * n * p), mom$meanSquare[1, 11], tolerance = 1e-3)
  expect_equal(sum(p), 1, tolerance = 1e-6) # Truncated at nMax
})
