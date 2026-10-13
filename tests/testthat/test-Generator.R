test_that("TempoGrid() validates input", {
  expect_error(TempoGrid(1), "nBin")
  expect_error(TempoGrid(NA), "nBin")
  expect_error(TempoGrid(5, c(1, -1)), "increasing")
  expect_error(TempoGrid(5, 1), "increasing")
  expect_error(TempoGrid(5, s = 0), "positive")
  grid <- TempoGrid(5, c(-2, 2), s = 2)
  expect_s3_class(grid, "TempoGrid")
  expect_equal(grid$x, 2 * grid$z)
  expect_equal(grid$h, 1)
})

test_that("TempoGenerator() returns valid generators", {
  grid <- TempoGrid(41, c(-4, 4))
  expect_error(TempoGenerator(list(), 0.01), "TempoGrid")
  expect_error(TempoGenerator(grid, -1), "theta")
  expect_error(TempoGenerator(TempoGrid(11, c(-8, 8)), 0.01), "h \\* max")
  for (scheme in c("central", "upwind", "sqra")) {
    Q <- TempoGenerator(grid, 0.01, scheme = scheme)
    expect_equal(rowSums(Q), numeric(41), tolerance = 1e-12)
    expect_true(all(Q[row(Q) != col(Q)] >= 0))
    expect_true(all(Q[abs(row(Q) - col(Q)) > 1] == 0))
  }
})

test_that("Central and upwind drift is exact in the interior", {
  grid <- TempoGrid(41, c(-4, 4), s = 1.5)
  interior <- 2:40
  for (scheme in c("central", "upwind")) {
    A <- TempoGenerator(grid, 0.02, alpha = 0, scheme = scheme)
    expect_equal(as.vector(A %*% grid$z)[interior], -0.02 * grid$z[interior],
                 tolerance = 1e-12)
  }
  # alpha scales rows by exp(alpha * x)
  A <- TempoGenerator(grid, 0.02, alpha = 0)
  Q <- TempoGenerator(grid, 0.02, alpha = 1)
  expect_equal(Q, A * exp(grid$x))
})

test_that("SqRA stationary law is the grid Gaussian", {
  grid <- TempoGrid(33, c(-4, 4))
  A <- TempoGenerator(grid, 0.01, alpha = 0, scheme = "sqra")
  p <- TempoEquilibrium(A)$p
  gauss <- exp(-grid$z^2 / 2)
  expect_equal(p, gauss / sum(gauss), tolerance = 1e-10)
})

test_that("BuddMannGenerator() reproduces ReacTran's tran.1D operator", {
  bm <- BuddMannGenerator()
  expect_equal(bm$h, 20 / 101)
  expect_equal(bm$x[51], 0, tolerance = 1e-12)
  expect_equal(rowSums(bm$Q), numeric(101), tolerance = 1e-9)
  down <- bm$Q[cbind(2:101, 1:100)]
  expect_equal(which(down < 0) + 1L, which(bm$x < -1 / bm$h)[-1])
  skip_if_not_installed("ReacTran")
  theta <- 0.01
  s <- 1.3
  bm <- BuddMannGenerator(theta = theta, s = s)
  grid <- ReacTran::setup.grid.1D(x.up = -10, x.down = 10, N = 101)
  set.seed(1)
  y <- runif(101)
  diffusion <- ReacTran::tran.1D(C = y, A = 1, D = 1, v = 0, dx = grid)$dC
  drift <- ReacTran::tran.1D(C = y, A = 1, D = 0, v = 1, dx = grid)$dC
  reference <- (diffusion * theta * s^2 + drift * theta * grid$x.mid) *
    exp(grid$x.mid)
  expect_equal(as.vector(bm$Q %*% y), reference, tolerance = 1e-10)
})

test_that("BinRates() scales rates by tempo", {
  rates <- BinRates(c(-1, 0, 1), 2, 1, psi = 0.5, beta = 1, delta = 0)
  expect_equal(rates$lambda, 2 * exp(c(-1, 0, 1)))
  expect_equal(rates$mu, exp(c(-1, 0, 1)))
  expect_equal(rates$psi, rep(0.5, 3))
  expect_equal(BinRates(TempoGrid(3, c(-1, 1)), 2, 1)$mu, exp(c(-1, 0, 1)))
})
