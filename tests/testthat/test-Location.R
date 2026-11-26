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
  rm(list = ls(.legacyWarned), envir = .legacyWarned)
  expect_warning(Q <- TempoGenerator(legacy, 0.2, scheme = "central"),
                 class = "legacyLocation")
  expect_identical(Q, TempoGenerator(bm, 0.2, scheme = "central"))
  expect_identical(TempoGenerator(legacy, 0.2, scheme = "sqra"),
                   TempoGenerator(bm, 0.2, scheme = "sqra"))
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
  # Rows 1-32 SqRA, 33-64 upwind, 65-96 central; within each, 1-16 root and
  # 17-32 origin start
  cases <- cases[c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 52, 59, 65, 70, 75, 80, 86, 91, 96), ]
  expect_setequal(paste(cases$scheme, cases$start),
                  paste(rep(c("sqra", "upwind", "central"), each = 2), c("root", "origin")))
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
  rm(list = ls(.legacyWarned), envir = .legacyWarned)
  expect_warning(legacySet <- .Setup(par, legacy), class = "legacyLocation")
  expect_identical(legacySet, .Setup(par, bm))
  # Once per session and kind of object
  expect_no_warning(.Setup(par, legacy))
  grid <- TempoGrid(11, c(-2, 2), 0.9, "bm")
  grid$location <- NULL
  expect_warning(legacyQ <- TempoGenerator(grid, 0.4), class = "legacyLocation")
  expect_identical(legacyQ, TempoGenerator(TempoGrid(11, c(-2, 2), 0.9, "bm"), 0.4))
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
  # SimulateCET() passes the location's OU mean to every segment
  ouMeans <- numeric(0)
  segment <- .OUSegment
  local_mocked_bindings(.OUSegment = function(width, x0, theta, s, dw, decay, noise,
                                              ouMean = 0) {
    ouMeans <<- c(ouMeans, ouMean)
    segment(width, x0, theta, s, dw, decay, noise, ouMean)
  })
  for (location in c("mean", "bm")) {
    ouMeans <- numeric(0)
    set.seed(6)
    history <- SimulateCET(5, theta = 0.5, s = 0.7, lambda = 0.3, mu = 0.1,
                           location = location)
    expect_s3_class(history, "TempoHistory")
    expect_gt(length(ouMeans), 0)
    expect_equal(unique(ouMeans), if (location == "mean") m else 0)
  }
})

test_that("location = \"bm\" reproduces the chains of the pre-#116 default", {
  # Column sums of a seeded chain with every move, free s and characters, as
  # produced by origin/main before #116 (e563e0c; identical() samples there and
  # here, in both engines). The tolerance allows for platform round-off only.
  tree <- TempoTree(locTree)
  chars <- matrix(c(0, 1, 1, 0, 1, 1, 1, 0, 0, 0), 5,
                  dimnames = list(locTree$tip.label, NULL))
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.6, gamma = 1, location = "bm")
  par <- list(lambda = 0.6, mu = 0.25, psi = 0.15, theta = 0.4, s = 0.9, rate = 0.3)
  priors <- c(TempoPriors(), rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  expected <- list(
    R = c(lambda = 8.7267744725229, s = 19.8004214405854, rate = 20.6907687532836,
          logPost = -635.747520580079, meanX = -9.46748009258043,
          treeLength = 294.830494052661, rootBin = 152),
    cpp = c(lambda = 8.7267744725229, s = 19.8004214405854, rate = 20.6907687532836,
            logPost = -635.74752057641, meanX = -9.46748009257895,
            treeLength = 294.830494052661, rootBin = 152))
  for (engine in names(expected)) {
    old <- options(EvoTempo.engine = engine)
    set.seed(11)
    run <- TempoMCMC(tree, model, par, priors, nIter = 30, nSlice = 6, characters = chars,
                     moves = c(tree = 1, node = 1, tip = 1, branch = 1, joint = 1, shift = 1,
                               nodeAge = 1, fossilAge = 1, slide = 1, saJump = 1))
    options(old)
    expect_equal(colSums(run$samples[names(expected[[engine]])]), expected[[engine]],
                 tolerance = 1e-9, info = engine)
  }
})

test_that("The central scheme's bound is an infeasible proposal, not an error", {
  # Under "mean" the bound h * max|z - alpha s| <= 2 depends on s: here s <= 3
  # ((a:1,b:1):1,c:2);
  tree <- TempoTree(LocPhylo(rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)), c(1, 1, 1, 2),
                             c("a", "b", "c")))
  par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1)
  cm <- DCETModel(nBin = 41, zRange = c(-5, 5), scheme = "central", rho = 1)
  expect_true(is.finite(TreeLikelihood(tree, modifyList(par, list(s = 2.9)), cm)))
  expect_error(TreeLikelihood(tree, modifyList(par, list(s = 3.5)), cm),
               class = "infeasibleProposal")
  priors <- sapply(names(par), function(n) function(x) 0, simplify = FALSE)
  set.seed(1)
  run <- MarginalMCMC(tree, cm, par, priors, nIter = 60)
  expect_true(all(run$s <= 3))
  expect_gt(max(run$s), 1.5)
})
