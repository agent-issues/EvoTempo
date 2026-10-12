test_that("E collapses to the analytic FBD p0 when beta = delta = 0", {
  grid <- TempoGrid(21, c(-3, 3))
  Q <- TempoGenerator(grid, 0.05)
  ages <- c(0, 0.5, 5, 50)
  E <- ExtinctionProbability(ages, Q, lambda = 0.8, mu = 0.3, psi = 0.1,
                             rho = 0.4)
  expected <- FBDExtinction(ages, 0.8, 0.3, 0.1, 0.4)
  expect_equal(E, matrix(expected, length(ages), 21), tolerance = 1e-9)
  expect_equal(FBDExtinction(0, 0.8, 0.3, 0.1, 0.4), 0.6)
})

test_that("Critical FBD branch of p0 is correct", {
  Q <- TempoGenerator(TempoGrid(5, c(-1, 1)), 0.1)
  E <- ExtinctionProbability(c(1, 10), Q, 0.5, 0.5, rho = 0.7)
  expect_equal(E[, 1], FBDExtinction(c(1, 10), 0.5, 0.5, 0, 0.7),
               tolerance = 1e-9)
})

test_that("E tends to the tempo-invariant fixed point mu / lambda", {
  grid <- TempoGrid(21, c(-2, 2))
  Q <- TempoGenerator(grid, 0.1)
  rates <- BinRates(grid, 1, 0.5)
  E <- ExtinctionProbability(c(1, 400), Q, rates$lambda, rates$mu)
  expect_equal(E[2, ], rep(0.5, 21), tolerance = 1e-8)
  expect_true(all(E >= 0 & E <= 1))
  expect_true(sd(E[1, ]) > 0.01)
})

test_that("ExtinctionProbability() validates input", {
  Q <- TempoGenerator(TempoGrid(5, c(-1, 1)), 0.1)
  expect_error(ExtinctionProbability(1, Q, 1, 0.5, rho = 2), "rho")
  expect_error(ExtinctionProbability(c(2, 1), Q, 1, 0.5), "sorted")
  expect_equal(ExtinctionProbability(0, Q, 1, 0.5, rho = 0.3),
               matrix(0.7, 1, 5))
})
