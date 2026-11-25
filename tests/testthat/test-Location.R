# The mean-rate tempo location (#116; dev/spec/tempo-location.md)

LocPhylo <- function(edge, len, labels) {
  structure(list(edge = edge, edge.length = len, tip.label = labels,
                 Nnode = length(labels) - 1L), class = "phylo")
}
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
locTree <- LocPhylo(rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3),
                          c(6, 9), c(9, 4), c(9, 5)),
                    c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
                    c("a", "b", "f", "c", "d"))
locStem <- locTree
locStem$root.edge <- 0.8
locStem$tip.label[3] <- "g" # an ordinary tip on a stemmed tree
locStem$edge.length[5] <- 0.4

# Budd & Mann parameters equivalent to `par` under the mean-rate location
# (spec L.3): every rate times exp(exponent * m), m = (alpha - 1/2) s^2
ToBuddMann <- function(par, model) {
  m <- (model$alpha - 0.5) * par$s^2
  expo <- c(theta = model$alpha, lambda = model$beta, mu = model$beta,
            psi = model$delta, rate = model$gamma)
  for (p in intersect(names(expo), names(par))) {
    par[[p]] <- par[[p]] * exp(expo[[p]] * m)
  }
  par
}

test_that("TempoGrid() places log-tempo by location", {
  grid <- TempoGrid(5, c(-2, 2), s = 2)
  expect_equal(grid$location, "mean")
  expect_equal(grid$x, 2 * grid$z - 2)
  expect_equal(TempoGrid(5, c(-2, 2), 2, "bm")$x, 2 * grid$z)
  expect_error(TempoGrid(5, location = "median"))
})

test_that("Default schemes and the centred positivity bound", {
  mean <- TempoGrid(26, c(-5, 5), s = 1)
  expect_equal(TempoGenerator(mean, 0.1),
               TempoGenerator(mean, 0.1, scheme = "sqra"))
  bm <- TempoGrid(26, c(-5, 5), s = 1, "bm")
  expect_equal(TempoGenerator(bm, 0.1),
               TempoGenerator(bm, 0.1, scheme = "central"))
  # Drift centre alpha s: h * max|z - s| = 0.4 * 6 > 2
  expect_error(TempoGenerator(mean, 0.1, scheme = "central"), "z - centre")
  expect_silent(TempoGenerator(mean, 0.1, alpha = 0, scheme = "central"))
  expect_equal(DCETModel()$scheme, "sqra")
  expect_equal(DCETModel(location = "bm")$scheme, "central")
  expect_equal(DCETModel(scheme = "upwind")$scheme, "upwind")
  expect_error(DCETModel(location = "zero"))
})

test_that("A grid without a location is Budd & Mann's", {
  bm <- TempoGrid(21, c(-4, 4), s = 1.3, "bm")
  legacy <- structure(bm[c("z", "x", "h", "s")], class = "TempoGrid")
  for (scheme in c("central", "sqra")) {
    expect_identical(TempoGenerator(legacy, 0.2, scheme = scheme),
                     TempoGenerator(bm, 0.2, scheme = scheme))
  }
})

test_that("The stationary law is the grid Gaussian, centred on the grid", {
  for (alpha in c(1, 0.5, 0)) {
    for (s in c(0.3, 1, 2)) {
      grid <- TempoGrid(41, c(-5, 5), s = s)
      p <- TempoEquilibrium(TempoGenerator(grid, 0.05, alpha = alpha))$p
      gauss <- dnorm(grid$z)
      expect_equal(p, gauss / sum(gauss), tolerance = 1e-12)
      expect_equal(sum(p * grid$x), -s^2 / 2, tolerance = 1e-12)
      # Mean tempo is 1 up to truncation of the tilted law N(s, 1) at z = 5
      expect_equal(sum(p * exp(grid$x)), 1,
                   tolerance = 2 * pnorm(5 - s, lower.tail = FALSE) + 1e-12)
    }
  }
})

test_that("Mean-rate model = Budd & Mann on the shifted window (spec L.3)", {
  cases <- expand.grid(alpha = c(1, 0), beta = c(1, 0), gamma = c(1, 0),
                       delta = c(0, 1), start = c("root", "origin"),
                       scheme = c("sqra", "upwind", "central"),
                       stringsAsFactors = FALSE)
  cases <- cases[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 47, 48), ]
  par <- list(lambda = 0.6, mu = 0.25, psi = 0.15, theta = 0.4, s = 0.9,
              rate = 0.3)
  character <- c(0, 1, 1, 0, 1)
  for (i in seq_len(nrow(cases))) {
    cs <- cases[i, ]
    tree <- TempoTree(if (cs$start == "origin") locStem else locTree)
    # Central needs a fine enough grid for the centred drift
    nBin <- if (cs$scheme == "central") 41 else 21
    shift <- cs$alpha * par$s
    Model <- function(location, zRange) {
      DCETModel(nBin = nBin, zRange = zRange, scheme = cs$scheme,
                alpha = cs$alpha, beta = cs$beta, gamma = cs$gamma,
                delta = cs$delta, rho = 0.6, start = cs$start,
                condition = "survival", location = location)
    }
    mean <- Model("mean", c(-3, 3))
    bm <- Model("bm", c(-3, 3) - shift)
    for (engine in c("R", "cpp")) {
      old <- options(EvoTempo.engine = engine)
      expect_equal(TreeLikelihood(tree, par, mean),
                   TreeLikelihood(tree, ToBuddMann(par, mean), bm),
                   tolerance = 1e-10, info = paste(i, engine))
      expect_equal(
        TreeLikelihood(tree, par, mean, character = character, nState = 2),
        TreeLikelihood(tree, ToBuddMann(par, mean), bm, character = character,
                       nState = 2),
        tolerance = 1e-10, info = paste(i, engine, "character")
      )
      options(old)
    }
  }
})

test_that("A model without a location is Budd & Mann's", {
  par <- list(lambda = 0.6, mu = 0.25, psi = 0.15, theta = 0.4, s = 0.9)
  bm <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.6, location = "bm")
  legacy <- bm
  legacy$location <- NULL
  expect_identical(.Setup(par, legacy), .Setup(par, bm))
  expect_identical(TreeLikelihood(TempoTree(locTree), par, legacy),
                   TreeLikelihood(TempoTree(locTree), par, bm))
})

test_that("SimulateCET() reverts to the location's mean", {
  m <- 0.7^2 / 2
  set.seed(4)
  bm <- .OUSegment(3, 0.5, 0.2, 0.7, 0.1, exp(-0.02), 0.7 * sqrt(1 - exp(-0.04)))
  set.seed(4)
  mean <- .OUSegment(3, 0.5 + m, 0.2, 0.7, 0.1, exp(-0.02),
                     0.7 * sqrt(1 - exp(-0.04)), ouMean = m)
  expect_equal(mean$x, bm$x + m)
  expect_equal(mean$dw, bm$dw)
  ex <- exp(-mean$x)
  expect_equal(mean$dt, mean$dw * (ex[-1] + ex[-length(ex)]) / 2)
  # Long-run lineage law in tempo-time: N(s^2/2, s^2) under "mean"
  set.seed(5)
  long <- .OUSegment(2e4, 0, 1, 0.7, 0.5, exp(-0.5), 0.7 * sqrt(1 - exp(-1)),
                     ouMean = m)
  expect_lt(abs(mean(long$x) - m), 0.03)
  expect_lt(abs(sd(long$x) - 0.7), 0.03)
  set.seed(6)
  history <- SimulateCET(5, theta = 0.5, s = 0.7, lambda = 0.3, mu = 0.1,
                         location = "bm")
  expect_s3_class(history, "TempoHistory")
})
