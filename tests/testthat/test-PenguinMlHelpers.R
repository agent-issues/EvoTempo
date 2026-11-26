# Helpers of the penguin common-footing summary (data-raw; review of #111).
# data-raw is not in the built package, so these run only from the source tree.
LoadMlHelpers <- function() {
  f <- test_path("..", "..", "data-raw", "empirical", "penguins", "cet", "ml-helpers.R")
  skip_if_not(file.exists(f))
  env <- new.env()
  sys.source(f, envir = env)
  env
}

test_that("CetPooledSe() is sd / sqrt(n) over replicates, floored at the estimator SE", {
  env <- LoadMlHelpers()
  # The four CET (0, 0) replicates of the #111 summary; the formula it replaced,
  # sqrt(se^2 + sd^2 / 2), gave 0.670 here (R111-1)
  reps <- c(-1753.2563, -1751.2199, -1751.6697, -1752.5823)
  expect_equal(env$CetPooledSe(0.179, reps), sd(reps) / 2)
  expect_lt(env$CetPooledSe(0.179, reps), 0.5)
  # The SE shrinks with more replicates
  expect_equal(env$CetPooledSe(0.1, c(reps, reps)), sd(c(reps, reps)) / sqrt(8))
  expect_equal(env$CetPooledSe(0.9, reps), 0.9)
  # One replicate: the estimator SE, not NA (R111-9)
  expect_equal(env$CetPooledSe(0.37, reps[1]), 0.37)
  expect_equal(env$CetPooledSe(0.37, c(reps[1], NA)), 0.37)
})

test_that("PoolProportion() weights mean and SE alike", {
  env <- LoadMlHelpers()
  p <- c(0.1, 0.2, 0.3, 0.2)
  expect_equal(env$PoolProportion(p, rep(10, 4)), c(p = mean(p), se = sd(p) / 2))
  w <- env$PoolProportion(c(0.1, 0.3), c(300, 100))
  expect_equal(w[["p"]], 0.15)
  expect_equal(w[["se"]], sqrt(2 * (0.75^2 * 0.05^2 + 0.25^2 * 0.15^2)))
  expect_true(is.na(env$PoolProportion(0.2, 5)[["se"]]))
})

test_that("CetLadderProblems() flags missing rungs, powers and chains", {
  env <- LoadMlHelpers()
  powers <- (0:3 / 4)^5
  Chains <- function(cell = "b0g0", reps = 1:2) {
    x <- expand.grid(power = c(powers, 1), rep = reps)
    data.frame(cell = cell, power = x$power, rep = x$rep, exists = TRUE, complete = TRUE)
  }
  ok <- Chains()
  expect_identical(env$CetLadderProblems(ok, powers), character())
  noTop <- ok[ok$power != 1, ]
  expect_match(env$CetLadderProblems(noTop, powers), "no t = 1 rung")
  gap <- ok[!(ok$rep == 2 & ok$power == powers[3]), ]
  expect_match(env$CetLadderProblems(gap, powers), "rep 2: 1 of 4 powers missing")
  dangling <- ok; dangling$exists[nrow(ok)] <- FALSE
  expect_match(env$CetLadderProblems(dangling, powers), "dangling")
  partial <- ok; partial$complete[2] <- FALSE
  expect_match(env$CetLadderProblems(partial, powers), "rep 1 power .* not complete")
})

test_that("EvenPowers() keeps t = 1 when a power is missing", {
  env <- LoadMlHelpers()
  powers <- (0:3 / 4)^5
  expect_identical(env$EvenPowers(c(powers, 1), powers), c(TRUE, FALSE, TRUE, FALSE, TRUE))
  expect_identical(env$EvenPowers(c(powers[-2], 1), powers), c(TRUE, TRUE, FALSE, TRUE))
})
