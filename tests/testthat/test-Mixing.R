# Moves and schedule from the mixing study (#19): shift move, TempoMoves(),
# fractional counts and the burn-in drop of whole-tree draws.
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
mixTree <- TempoTree(structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4), class = "phylo"))
mixPar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.4)
mixModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)
mixPriors <- c(TempoPriors(),
               rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
set.seed(31)
mixChars <- matrix(sample(0:1, 5 * 6, replace = TRUE), 5)

MixState <- function(model = mixModel, seed = 1) {
  set.seed(seed)
  set <- .AddSlices(mixTree, .Precompute(mixTree, mixPar, model), 4)
  msg <- .Messages(mixTree, set)
  d <- .SampleTree(mixTree, set, msg)
  st <- list(tree = mixTree, par = mixPar, set = set, msg = msg,
             paths = d[["paths"]], bin = d[["bin"]])
  st[["charLL"]] <- MixCharLL(model)(mixTree, st[["paths"]], mixPar, set)
  st
}

MixCharLL <- function(model = mixModel) {
  function(tr, paths, p, set) {
    MkLikelihood(tr, mixChars,
                 .CharacterLengths(tr, set[["grid"]][["x"]], paths, p[["rate"]],
                                   model[["gamma"]]))
  }
}

test_that("TempoMoves() scales local moves to the tree with characters", {
  base <- c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1, shift = 1)
  expect_equal(TempoMoves(mixTree), base)
  expect_equal(TempoMoves(mixTree, mixChars[, 0, drop = FALSE]), base)
  noGamma <- mixModel
  noGamma$gamma <- 0
  expect_equal(TempoMoves(mixTree, mixChars, noGamma), base)
  # 4 internal nodes, 4 tips with paths (f is a sampled ancestor), 7 branches
  expect_equal(TempoMoves(mixTree, mixChars, mixModel),
               c(node = 1, tip = 1, branch = 2, shift = 1))
  expect_equal(TempoMoves(mixTree, mixChars, local = 1),
               c(node = 4, tip = 4, branch = 7, shift = 1))
})

test_that("Fractional move counts are drawn independently of the state", {
  set.seed(2)
  expect_equal(.MoveCount(3), 3)
  expect_equal(stats::runif(1), {set.seed(2); stats::runif(1)})
  counts <- replicate(4000, .MoveCount(1.25))
  expect_setequal(unique(counts), c(1, 2))
  expect_equal(mean(counts), 1.25, tolerance = 0.05)
})

test_that("The shift map is an involution that keeps the character lengths", {
  st <- MixState()
  free <- c("rate", "lambda", "mu", "theta", "s")
  for (m in c(-1L, 1L)) {
    fwd <- .Shifted(st, m, mixModel, free)
    if (is.null(fwd)) next
    expect_equal(fwd$comp, c("rate", "lambda", "mu", "theta"))
    d <- m * st$set$grid$h * st$par$s
    expect_equal(fwd$par$rate, st$par$rate * exp(-d))
    expect_equal(fwd$par$s, st$par$s)
    expect_equal(fwd$par$psi, st$par$psi)
    expect_equal(fwd$logJac, -4 * d)
    expect_identical(fwd$bin, st$bin + m)
    nodes <- .PathNodes(mixTree)
    for (v in nodes) {
      expect_identical(fwd$paths[[v]]$t, st$paths[[v]]$t)
      expect_identical(fwd$paths[[v]]$k, st$paths[[v]]$k + m)
    }
    x <- st$set$grid$x
    expect_equal(.CharacterLengths(mixTree, x, fwd$paths, fwd$par$rate, 1),
                 .CharacterLengths(mixTree, x, st$paths, st$par$rate, 1))
    back <- .Shifted(c(fwd[c("par", "bin", "paths")], list(tree = mixTree, set = st$set)),
                     -m, mixModel, free)
    expect_equal(back$par, st$par)
    expect_identical(back$bin, st$bin)
    expect_identical(back$paths, st$paths)
    expect_equal(back$logJac, -fwd$logJac)
  }
  # Exponents select what is compensated
  model <- mixModel
  model[c("beta", "delta", "alpha")] <- list(0, 1, 0)
  expect_equal(.Shifted(st, -1L, model, c("rate", "psi", "theta"))$comp, c("rate", "psi"))
  # A bin at each end of the grid: every shift leaves it
  edge <- st
  v <- .PathNodes(mixTree)[[1]]
  edge$paths[[v]]$k[[1]] <- 1L
  edge$bin[[2]] <- 11L
  expect_null(.Shifted(edge, 1L, mixModel, free))
  expect_null(.Shifted(edge, -1L, mixModel, free))
  res <- .ShiftMove(edge, mixModel, free, mixPriors, 4, MixCharLL())
  expect_false(res$accepted)
  expect_identical(res$state, edge)
})

test_that("Shift moves are accepted on the exact target", {
  st <- MixState(seed = 3)
  CharLL <- MixCharLL()
  set.seed(4)
  acc <- 0
  for (i in 1:40) {
    res <- .ShiftMove(st, mixModel, c("rate", "lambda", "mu", "theta"), mixPriors,
                      4, CharLL)
    acc <- acc + res$accepted
    st <- res$state
    # The kept character likelihood is the recomputed one
    expect_equal(st$charLL, CharLL(st$tree, st$paths, st$par, st$set), tolerance = 1e-10)
    expect_true(.SameSet(st$set, st$tree, st$par))
  }
  expect_gt(acc, 0)
  # `rate` fixed: the character likelihood is recomputed
  for (i in 1:20) {
    res <- .ShiftMove(st, mixModel, "lambda", mixPriors, 4, CharLL)
    st <- res$state
    expect_equal(st$charLL, CharLL(st$tree, st$paths, st$par, st$set), tolerance = 1e-10)
  }
  # Only `rate` compensated: the frozen model is reused
  local_mocked_bindings(.Precompute = function(...) stop("not reused"))
  for (i in 1:5) {
    st <- .ShiftMove(st, mixModel, "rate", mixPriors, 4, CharLL)$state
  }
  expect_true(.SameSet(st$set, st$tree, st$par))
})

test_that("Shift moves reject infeasible solves and zero priors", {
  st <- MixState(seed = 5)
  set.seed(6)
  local_mocked_bindings(.Precompute = function(...) .Infeasible("ODE solver failed"))
  res <- .ShiftMove(st, mixModel, "lambda", mixPriors, 4, MixCharLL())
  expect_false(res$accepted)
  flat <- mixPriors
  flat$lambda <- function(x) if (x == 0.5) 0 else -Inf
  set.seed(1)
  expect_false(.ShiftMove(st, mixModel, "lambda", flat, 4, MixCharLL())$accepted)
})

test_that("TempoMCMC() runs the default schedule with characters", {
  set.seed(7)
  run <- TempoMCMC(mixTree, mixModel, mixPar, mixPriors[c("rate", "theta")],
                   nIter = 30, burnin = 4, nSlice = 4, characters = mixChars)
  expect_equal(run$moves, TempoMoves(mixTree, mixChars, mixModel))
  expect_gt(run$acceptance[["shift"]], 0)
  # Parameters fixed: shift moves bins only
  fixed <- TempoMCMC(mixTree, mixModel, mixPar, mixPriors[c("rate", "theta")],
                     nIter = 20, nSlice = 4, characters = mixChars,
                     moves = c(shift = 2, branch = 0.5), updateParameters = FALSE)
  expect_equal(unique(fixed$samples$rate), mixPar$rate)
  expect_equal(unique(fixed$samples$theta), mixPar$theta)
  expect_gt(fixed$acceptance[["shift"]], 0)
})

test_that("Whole-tree draws with low burn-in acceptance are dropped", {
  Run <- function(minAcceptance) {
    set.seed(8)
    TempoMCMC(mixTree, mixModel, mixPar, mixPriors["rate"], nIter = 3, burnin = 3,
              nSlice = 4, characters = mixChars,
              moves = c(tree = 1, joint = 1, branch = 1, shift = 1),
              minAcceptance = minAcceptance)
  }
  expect_equal(Run(0)$moves, c(tree = 1, branch = 1, joint = 1, shift = 1))
  expect_equal(Run(2)$moves, c(branch = 1, shift = 1))
})
