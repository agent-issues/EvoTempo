# Metropolis-coupled chains in TempoMCMC() (#110; dev/plans/mc3.md).
mc3Model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)
mc3Priors <- c(TempoPriors()[c("lambda", "theta")],
               rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
Mc3Run <- function(engine = "cpp", seed = 3, ...) {
  set.seed(seed)
  TempoMCMC(oracleTree, mc3Model, oraclePar, mc3Priors, nSlice = 4,
            characters = oracleChars, engine = engine, ...)
}

test_that("Coupled chains are validated", {
  for (bad in list(0, 1.5, NA, Inf, c(2, 2), "2")) {
    expect_error(Mc3Run(nIter = 1, nChain = bad), "positive whole numbers")
    expect_error(Mc3Run(nIter = 1, swapInterval = bad), "positive whole numbers")
  }
  for (bad in list(0, -1, NA, Inf, c(0.1, 0.2), "a")) {
    expect_error(Mc3Run(nIter = 1, deltaHeat = bad), "single positive number")
  }
  expect_error(
    TempoMCMC(oracleTree, mc3Model, oraclePar, mc3Priors, nIter = 1, nChain = 2),
    "needs `characters`"
  )
})

test_that("One chain draws the numbers of an uncoupled chain", {
  one <- Mc3Run(nIter = 6, burnin = 20, nChain = 1, deltaHeat = 3, swapInterval = 2)
  after <- stats::runif(1)
  ref <- Mc3Run(nIter = 6, burnin = 20)
  expect_identical(one, ref)
  expect_identical(after, stats::runif(1))
  expect_false(any(c("heats", "swaps") %in% names(one)))
})

test_that("Coupled chains match between engines", {
  moves <- c(node = 1.5, tip = 1, branch = 1, shift = 1, nodeAge = 1, slide = 1,
             fnpr = 1, saJump = 1)
  Run <- function(engine) {
    Mc3Run(engine, nIter = 12, burnin = 20, thin = 2, nChain = 3, deltaHeat = 0.5,
           moves = moves)
  }
  cpp <- Run("cpp")
  ref <- Run("R")
  expect_equal(cpp$samples, ref$samples, tolerance = 1e-8)
  expect_equal(cpp$occupancy, ref$occupancy, tolerance = 1e-8)
  expect_identical(cpp$acceptance, ref$acceptance)
  expect_identical(cpp$swaps, ref$swaps)
  expect_identical(cpp$heats, c(1, 1 / 1.5, 0.5))
  expect_named(cpp$swaps, c("hotter", "colder", "tries", "accepted", "rate"))
  expect_identical(cpp$swaps$colder, cpp$heats[1:2])
  expect_identical(cpp$swaps$hotter, cpp$heats[2:3])
  # One attempt per iteration, counted after burn-in
  expect_identical(sum(cpp$swaps$tries), 12)
  expect_gt(sum(cpp$swaps$accepted), 0)
  expect_identical(cpp$swaps$rate, cpp$swaps$accepted / pmax(cpp$swaps$tries, 1))
  # The cold state's stored character likelihood is unheated
  st <- cpp$state
  expect_equal(
    st$charLL,
    MkLikelihood(st$tree, oracleChars,
                 .CharacterLengths(st$tree, st$set$grid$x, st$paths, st$par$rate,
                                   mc3Model$gamma))
  )
  expect_identical(sum(Mc3Run(nIter = 12, nChain = 2, swapInterval = 5)$swaps$tries), 2)
})

test_that("Coupled chains with an origin and age moves match between engines", {
  phylo <- structure(
    list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9), c(9, 4), c(9, 5)),
         edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
         tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
    class = "phylo"
  )
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, start = "origin")
  Run <- function(engine) {
    set.seed(14)
    TempoMCMC(phylo, model, oraclePar, mc3Priors, nIter = 8, burnin = 20, nSlice = 4,
              characters = oracleChars, nChain = 2, deltaHeat = 1, engine = engine,
              moves = c(node = 1, tip = 1, branch = 1, shift = 1, nodeAge = 1, originAge = 1),
              agePriors = list(origin = function(a) stats::dunif(a, 4, 8, log = TRUE)))
  }
  cpp <- Run("cpp")
  ref <- Run("R")
  expect_equal(cpp$samples, ref$samples, tolerance = 1e-8)
  expect_identical(cpp$swaps, ref$swaps)
  expect_gt(length(unique(cpp$samples$originAge)), 1)
})

# Heat h at a level must act in every move exactly as a power h on the
# character likelihood does in an uncoupled chain.  A ladder (h1, h2) on the
# likelihood L is the same pair of targets as the ladder (1, h2 / h1) on
# L^h1, which a rebound MkLikelihood() supplies without any heat at level 1;
# a move that ignored its heat would be right in the second and wrong in the
# first.
test_that("Heat acts as a power on the character likelihood in every move", {
  phylo <- structure(
    list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9), c(9, 4), c(9, 5)),
         edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
         tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
    class = "phylo"
  )
  ranges <- cbind(c(0.2, 0, 0.5, 0.5, 0), c(2.5, 0, 3, 3.5, 0))
  shared <- c(tree = 1, node = 1, tip = 1, branch = 1, joint = 1, shift = 1, nodeAge = 1,
              fossilAge = 1, slide = 1, fnpr = 1, saJump = 1)
  cases <- list(
    list(model = mc3Model, moves = c(shared, rootAge = 1), priors = mc3Priors,
         agePriors = list(root = function(a) stats::dunif(a, 3.5, 6, log = TRUE))),
    list(model = DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, start = "origin"),
         moves = c(shared, originAge = 1), priors = mc3Priors,
         agePriors = list(origin = function(a) stats::dunif(a, 4, 8, log = TRUE))),
    # With `rate` fixed a shift changes the character likelihood
    list(model = mc3Model, moves = c(node = 1, tip = 1, shift = 3),
         priors = mc3Priors[c("lambda", "theta")], agePriors = list()),
    # Gamma rates among characters: the move on their shape
    list(model = mc3Model, moves = c(node = 1, tip = 1, slide = 1), priors = mc3Priors,
         agePriors = list(), characters = MkCharacters(oracleChars, nRateCat = 4),
         start = c(oraclePar, alpha = 0.8),
         charPriors = list(alpha = function(a) stats::dexp(a, 1, log = TRUE)))
  )
  Mk <- MkLikelihood
  for (case in cases) {
    Run <- function(ladder, power) {
      local_mocked_bindings(
        .HeatLadder = function(nChain, deltaHeat) ladder,
        MkLikelihood = function(...) power * Mk(...)
      )
      set.seed(21)
      TempoMCMC(TempoTree(phylo, ranges = ranges), case$model,
                if (is.null(case$start)) oraclePar else case$start, case$priors,
                nIter = 30, burnin = 20, nSlice = 4,
                characters = if (is.null(case$characters)) oracleChars else case$characters,
                moves = case$moves, agePriors = case$agePriors,
                charPriors = as.list(case$charPriors), minAcceptance = 0, nChain = 2)
    }
    heated <- Run(c(0.5, 0.2), 1)
    powered <- Run(c(1, 0.4), 0.5)
    expect_equal(heated$samples, powered$samples, tolerance = 1e-10)
    expect_identical(heated$acceptance, powered$acceptance)
    expect_identical(heated$swaps[c("tries", "accepted")], powered$swaps[c("tries", "accepted")])
    # Every move ran and was accepted at least once, and levels were exchanged
    expect_true(all(heated$acceptance[c(names(case$moves), names(case$charPriors))] > 0))
    expect_gt(sum(heated$swaps$accepted), 0)
    # The likelihood matters: the unheated chain differs
    expect_false(isTRUE(all.equal(heated$samples, Run(c(1, 0.4), 1)$samples)))
  }
})

test_that("Coupled chains do not depend on the cached character likelihood", {
  # The cache holds the last tree scored, whichever level it belongs to, and
  # diffs every call against it (dev/plans/mk-likelihood-cache.md); the levels
  # and the swaps interleave states, which must cost time only
  moves <- c(node = 1, tip = 1, branch = 1, shift = 1, nodeAge = 1, fossilAge = 1,
             slide = 1, fnpr = 2, saJump = 1)
  Run <- function(cache, engine = "cpp") {
    old <- options(EvoTempo.mkCache = cache)
    oldCheck <- Sys.getenv("EVOTEMPO_MK_CACHE_CHECK", NA)
    # With the cache on, every cached value is checked against a one-shot one
    Sys.setenv(EVOTEMPO_MK_CACHE_CHECK = if (cache) "1" else "0")
    on.exit({
      options(old)
      if (is.na(oldCheck)) {
        Sys.unsetenv("EVOTEMPO_MK_CACHE_CHECK")
      } else {
        Sys.setenv(EVOTEMPO_MK_CACHE_CHECK = oldCheck)
      }
    })
    Mc3Run(engine, nIter = 20, burnin = 20, nChain = 3, deltaHeat = 0.5, moves = moves,
           minAcceptance = 0)
  }
  .MkCache(clear = TRUE)
  cached <- Run(TRUE)
  expect_identical(.MkCache()$checkFailures, 0L)
  expect_gt(.MkCache()$recomputed, 0)
  plain <- Run(FALSE)
  expect_identical(cached$samples, plain$samples)
  expect_identical(cached$acceptance, plain$acceptance)
  expect_identical(cached$swaps, plain$swaps)
  expect_equal(cached$samples, Run(FALSE, "R")$samples, tolerance = 1e-8)
  # Tree moves and swaps were accepted, so states of different trees interleave
  expect_true(all(cached$acceptance[c("slide", "fnpr", "nodeAge")] > 0))
  expect_gt(sum(cached$swaps$accepted), 0)
})

test_that("Heated levels drop moves silently", {
  wholeTree <- c(tree = 1, joint = 1, node = 1, tip = 1, shift = 1)
  warned <- character(0)
  run <- withCallingHandlers(
    Mc3Run(nIter = 2, burnin = 100, nChain = 2, moves = wholeTree, minAcceptance = 2),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # One warning per dropped move, from the cold level only
  expect_length(warned, 2)
  expect_match(warned, "dropped after burn-in")
  expect_false(any(c("tree", "joint") %in% names(run$moves)))
})

test_that("Diagnostics and the move schedule are the cold level's", {
  # Whole-tree draws ignore the characters: accepted about 0.2-0.45 of the
  # time at heat 1 here, and nearly always at heat 0.001
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, gamma = 2)
  chars <- cbind(oracleChars, oracleChars, oracleChars)
  set.seed(1)
  expect_warning(
    run <- TempoMCMC(oracleTree, model, oraclePar, mc3Priors[c("lambda", "theta")],
                     nIter = 20, burnin = 100, nSlice = 4, characters = chars,
                     moves = c(tree = 1, node = 1, tip = 1), minAcceptance = 0.7,
                     nChain = 2, deltaHeat = 999),
    "`tree` dropped after burn-in"
  )
  expect_lt(run$burninAcceptance[["tree"]], 0.7)
  expect_named(run$moves, c("node", "tip"))
  expect_identical(run$acceptance[["tree"]], 0)
  # The hot level keeps `tree`: as a single chain at that power it is not dropped
  Mk <- MkLikelihood
  local_mocked_bindings(MkLikelihood = function(...) Mk(...) / 1000)
  set.seed(1)
  hot <- TempoMCMC(oracleTree, model, oraclePar, mc3Priors[c("lambda", "theta")],
                   nIter = 20, burnin = 100, nSlice = 4, characters = chars,
                   moves = c(tree = 1, node = 1, tip = 1), minAcceptance = 0.7)
  expect_gt(hot$burninAcceptance[["tree"]], 0.9)
  expect_named(hot$moves, c("tree", "node", "tip"))
})

test_that(".SwapHeats() exchanges levels, not states", {
  State <- function(logLik, heat) {
    list(tree = logLik, paths = NULL, par = NULL, set = NULL, charLL = heat * logLik)
  }
  CharLL <- function(tr, paths, p, set, h) h * tr
  Ladder <- function(logLik) {
    list(heats = c(1, 0.5), chain = list(State(logLik[[1]], 1), State(logLik[[2]], 0.5)),
         slot = 1:2, swapTry = 0, swapAcc = 0)
  }
  # The hotter state is the likelier: always accepted
  set.seed(1)
  up <- .SwapHeats(Ladder(c(-10, -4)), CharLL)
  expect_identical(up$slot, 2:1)
  expect_identical(c(up$swapTry, up$swapAcc), c(1, 1))
  expect_identical(vapply(up$chain, `[[`, 0, "charLL"), c(-5, -4))
  expect_identical(vapply(up$chain, `[[`, 0, "tree"), c(-10, -4))
  # Probability exp(0.5 * (-400 + 10)): rejected
  down <- .SwapHeats(Ladder(c(-10, -400)), CharLL)
  expect_identical(down$slot, 1:2)
  expect_identical(c(down$swapTry, down$swapAcc), c(1, 0))
  expect_identical(vapply(down$chain, `[[`, 0, "charLL"), c(-10, -200))
  # Acceptance frequency exp((h_a - h_b) (l_b - l_a)) = exp(-0.5)
  set.seed(2)
  ladder <- Ladder(c(-10, -11))
  accepted <- mean(vapply(1:4000, function(i) .SwapHeats(ladder, CharLL)$swapAcc, 0))
  expect_lt(abs(accepted - exp(-0.5)), 4 * sqrt(exp(-0.5) * (1 - exp(-0.5)) / 4000))
  # Two impossible states do not swap; an impossible colder state always does
  none <- .SwapHeats(Ladder(c(-Inf, -Inf)), CharLL)
  expect_identical(none$swapAcc, 0)
  expect_identical(.SwapHeats(Ladder(c(-Inf, -4)), CharLL)$slot, 2:1)
  expect_identical(.SwapHeats(Ladder(c(-4, -Inf)), CharLL)$slot, 1:2)
  # Both uniforms are drawn whatever the outcome
  set.seed(3)
  .SwapHeats(Ladder(c(-Inf, -Inf)), CharLL)
  after <- stats::runif(1)
  set.seed(3)
  stats::runif(2)
  expect_identical(after, stats::runif(1))
})

test_that("The cold level of coupled chains is exact", {
  skip_on_cran()
  res <- OracleMc3(seed = 100, nIter = if (identical(Sys.getenv("R_COVR"), "true")) 600 else 4000)
  expect_gt(min(res$swaps), 0.05)
  expect_lt(max(abs(res$z)), 4)
})
