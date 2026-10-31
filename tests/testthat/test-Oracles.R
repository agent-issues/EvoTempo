# Exact-target oracles (helper-oracles.R; #51).  Each pools batch-means
# z-scores of the chain's mean log-parameter and its CDF at the exact
# quartiles.  data-raw/oracles-51.R sweeps seeds, and mutations of the shift
# move that each oracle must detect, at these chain lengths.
zMax <- 4

test_that("Shift moves compensate psi exactly (delta != 0)", {
  skip_on_cran()
  res <- OraclePsi(seed = 100, nIter = 3000)
  expect_gt(res$shift, 0.05)
  expect_lt(max(abs(res$z)), zMax)
})

test_that("Shift moves keep the character likelihood exactly with s free", {
  res <- ShiftIdentityError()
  expect_gt(res[["accepted"]], 5)
  expect_lt(res[["error"]], 1e-10)
})

test_that("Shift moves are exact with s free", {
  skip_on_cran()
  res <- OracleS(seed = 100)
  expect_gt(res$shift, 0.05)
  expect_lt(max(abs(res$z)), zMax)
})

test_that("The gamma = 0 default schedule is exact", {
  skip_on_cran()
  res <- OracleGamma0(seed = 100)
  # Whole-tree draws are kept after burn-in at gamma = 0
  expect_equal(res$moves, c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1, shift = 1))
  expect_gt(res$shift, 0.05)
  expect_lt(max(abs(res$z)), zMax)
})
