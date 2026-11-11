# MkCharacters(): states per character, ordered characters, ambiguity,
# Mkv per missing-data mask and discrete Gamma rates (#47, #97).
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
mkPhylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
  class = "phylo")
mkTree <- TempoTree(mkPhylo)
mkPar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.4)
mkModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)
set.seed(11)
mkLengths <- runif(9, 0.05, 0.6)
mkLengths[3] <- 0  # the sampled ancestor
# Five characters: binary, ordered three-state with a set, unordered
# four-state, ordered four-state with missing cells, binary with {01}
mkCells <- matrix(c("0", "1", "1", "0", "?",
                    "0", "{01}", "2", "1", "2",
                    "3", "0", "1", "1", "2",
                    "0", "3", "-", "2", "1",
                    "1", "{01}", "0", "0", "1"), 5)
mkOrdered <- c(FALSE, TRUE, FALSE, TRUE, FALSE)

# Matrix exponential by scaling and squaring of a Taylor series
MatExp <- function(A) {
  j <- max(0, ceiling(log2(max(abs(A)) * 10)))
  A <- A / 2^j
  out <- term <- diag(nrow(A))
  for (k in 1:25) {
    term <- term %*% A / k
    out <- out + term
  }
  for (i in seq_len(j)) out <- out %*% out
  out
}

test_that("GammaRates() matches fnDiscretizeGamma", {
  # Yang (1994), Table 1: alpha = 0.5, four categories
  expect_equal(GammaRates(0.5), c(0.0334, 0.2519, 0.8203, 2.8944), tolerance = 2e-3)
  for (a in c(0.1, 1.1, 20)) {
    expect_equal(mean(GammaRates(a, 6)), 1, tolerance = 1e-12)
    expect_true(all(diff(GammaRates(a, 6)) > 0))
  }
  expect_equal(GammaRates(3, 1), 1)
  expect_error(GammaRates(0), "positive")
  expect_error(GammaRates(c(1, 2)), "positive")
})

test_that("MkCharacters() reads states, sets and unknowns", {
  num <- matrix(c(0, 1, NA, 2, 1, 0, 1, 1, 0, 0), 5)
  mk <- MkCharacters(num)
  expect_s3_class(mk, "MkCharacters")
  expect_equal(mk$nState, c(3L, 2L))
  expect_equal(mk$nChar, 2L)
  expect_equal(.NCharacters(mk), 2L)
  expect_equal(.NCharacters(num), 2L)
  chr <- MkCharacters(matrix(c("0", "1", "?", "2", "1", "0", "1", "1", "0", "0"), 5))
  expect_equal(chr$partitions, mk$partitions)
  # A full set is unknown; a partial set sums over its states
  full <- MkCharacters(matrix(c("0", "1", "{01}", "0", "1"), 5))
  gap <- MkCharacters(matrix(c("0", "1", "?", "0", "1"), 5))
  expect_equal(full$partitions, gap$partitions)
  part <- MkCharacters(matrix(c("0", "1", "(01)", "0", "2"), 5))
  expect_equal(array(part$partitions[[1]]$tip, c(3, 1, 5))[, 1, 3], c(1, 1, 0))
  # A two-state ordered character is unordered
  expect_false(MkCharacters(num, ordered = TRUE)$ordered[[2]])
  expect_true(MkCharacters(num, ordered = TRUE)$ordered[[1]])
  expect_error(MkCharacters(matrix(c("0", "x", "1"), 3)), "Cannot read")
  expect_error(MkCharacters(matrix(c(0, -1, 1), 3)), "whole numbers")
  expect_error(MkCharacters(matrix(0:1, 1)), "two rows")
  expect_error(MkCharacters(num, nState = 2), "exceed")
  expect_error(MkCharacters(num, nState = c(3, 2, 2)), "one per character")
  expect_error(MkCharacters(num, nState = 1.5), "whole number")
  expect_error(MkCharacters(num, ordered = NA), "TRUE or FALSE")
  expect_error(MkCharacters(num, nRateCat = 0), "nRateCat")
  expect_error(MkCharacters(num, nRateCat = 2.5), "nRateCat")
  expect_error(MkCharacters(num, coding = "informative"))
  # Not variable: every scored cell holds the same state (or set)
  expect_error(MkCharacters(matrix(c(1, 1, NA, 1, 1, 0, 1, 1, 0, 0), 5),
                            coding = "variable"), "Character\\(s\\) 1 cannot")
  expect_error(MkCharacters(matrix(c("{01}", "{01}", "?", "{01}", "{01}"), 5, 1)[, c(1, 1)],
                            nState = 3, coding = "variable"), "1, 2 cannot")
  # Consistent with a constant character only through a set: kept, as RevBayes
  expect_s3_class(MkCharacters(matrix(c("{01}", "0", "?", "0", "0"), 5), nState = 3,
                               coding = "variable"), "MkCharacters")
})

test_that("Plain Mk on MkCharacters equals MkLikelihood() on a matrix", {
  set.seed(3)
  chars <- matrix(sample(0:2, 5 * 6, replace = TRUE), 5)
  chars[2, 3] <- NA
  mk <- MkCharacters(chars, nState = 3)
  for (engine in c("cpp", "R")) {
    expect_equal(MkLikelihood(mkTree, mk, mkLengths, engine = engine),
                 MkLikelihood(mkTree, chars, mkLengths, 3, engine = "R"),
                 tolerance = 1e-12)
  }
  # nRateCat = 1 ignores alpha
  expect_equal(MkLikelihood(mkTree, mk, mkLengths, alpha = 3),
               MkLikelihood(mkTree, mk, mkLengths))
})

test_that("Ordered characters and sets follow the rate matrix", {
  cherry <- structure(list(edge = rbind(c(3, 1), c(3, 2)), edge.length = c(1, 1),
                           tip.label = c("a", "b"), Nnode = 1), class = "phylo")
  lengths <- c(0.3, 0.5, 0)
  for (K in 3:5) {
    Q <- .MkRateMatrix(K, TRUE)
    expect_equal(-sum(diag(Q)) / K, 1)
    expect_true(all(Q[abs(row(Q) - col(Q)) > 1] == 0))
    P <- MatExp(Q * 0.8)
    states <- matrix(c(0, K - 1, 1, 1, 0, 2), 2)
    mk <- MkCharacters(states, nState = K, ordered = TRUE)
    expected <- sum(log(P[cbind(states[1, ] + 1, states[2, ] + 1)] / K))
    for (engine in c("cpp", "R")) {
      expect_equal(MkLikelihood(cherry, mk, lengths, engine = engine), expected,
                   tolerance = 1e-10)
    }
    # A set sums over its states
    set <- MkCharacters(matrix(c("0", "{12}"), 2), nState = K, ordered = TRUE)
    expect_equal(MkLikelihood(cherry, set, lengths),
                 log(sum(P[1, 2:3]) / K), tolerance = 1e-10)
  }
})

test_that("MkCharacters likelihoods agree between engines", {
  for (coding in c("all", "variable")) {
    for (nCat in c(1, 4)) {
      mk <- MkCharacters(mkCells, ordered = mkOrdered, coding = coding,
                         nRateCat = nCat)
      for (alpha in c(0.2, 1.1, 8)) {
        cpp <- MkLikelihood(mkTree, mk, mkLengths, alpha = alpha, engine = "cpp")
        ref <- MkLikelihood(mkTree, mk, mkLengths, alpha = alpha, engine = "R")
        expect_true(is.finite(cpp))
        expect_equal(cpp, ref, tolerance = 1e-10)
      }
    }
  }
  mk <- MkCharacters(mkCells, nRateCat = 4)
  expect_error(MkLikelihood(mkTree, mk, mkLengths), "`alpha` is needed")
  expect_error(MkLikelihood(TempoTree(structure(
    list(edge = rbind(c(3, 1), c(3, 2)), edge.length = c(1, 1),
         tip.label = c("a", "b"), Nnode = 1), class = "phylo")),
    mk, c(1, 1, 0), alpha = 1), "one row per tip")
})

test_that("The Mkv correction conditions on variability per mask", {
  # Every variable pattern of a mask: conditional probabilities sum to 1
  Patterns <- function(K, nScored) {
    grid <- as.matrix(expand.grid(rep(list(seq_len(K) - 1L), nScored)))
    grid[apply(grid, 1, function(x) length(unique(x)) > 1), , drop = FALSE]
  }
  cases <- list(list(K = 2, ordered = FALSE, missing = integer(0), nCat = 4),
                list(K = 2, ordered = FALSE, missing = 3L, nCat = 1),
                list(K = 3, ordered = TRUE, missing = 5L, nCat = 4),
                list(K = 3, ordered = FALSE, missing = integer(0), nCat = 3))
  for (case in cases) {
    scored <- setdiff(1:5, case$missing)
    pats <- Patterns(case$K, length(scored))
    total <- 0
    for (i in seq_len(nrow(pats))) {
      x <- rep(NA_integer_, 5)
      x[scored] <- pats[i, ]
      mk <- MkCharacters(matrix(x, 5), nState = case$K, ordered = case$ordered,
                         coding = "variable", nRateCat = case$nCat)
      total <- total + exp(MkLikelihood(mkTree, mk, mkLengths, alpha = 0.6,
                                        engine = "R"))
    }
    expect_equal(total, 1, tolerance = 1e-12)
  }
  # RevBayes's mixture correction: log mean_c L_c - log(1 - mean_c P_c)
  x <- matrix(c(0, 1, NA, 1, 2), 5)
  rates <- GammaRates(0.6, 4)
  L <- vapply(rates, function(r) {
    exp(MkLikelihood(mkTree, MkCharacters(x, nState = 3), mkLengths * r))
  }, 0)
  P <- vapply(rates, function(r) {
    3 * exp(MkLikelihood(mkTree, MkCharacters(matrix(c(0, 0, NA, 0, 0), 5), nState = 3),
                         mkLengths * r))
  }, 0)
  mk <- MkCharacters(x, nState = 3, coding = "variable", nRateCat = 4)
  for (engine in c("cpp", "R")) {
    expect_equal(MkLikelihood(mkTree, mk, mkLengths, alpha = 0.6, engine = engine),
                 log(mean(L)) - log1p(-mean(P)), tolerance = 1e-12)
  }
  # Characters with one mask share one correction
  two <- MkCharacters(cbind(x, c(1, 0, NA, 1, 1)), nState = 3, coding = "variable")
  expect_equal(two$partitions[[1]]$count, 2L)
  expect_equal(two$partitions[[1]]$nCol, 3L)
})

test_that("TempoMCMC() samples the Gamma shape, identically in both engines", {
  mk <- MkCharacters(mkCells, ordered = mkOrdered, coding = "variable", nRateCat = 4)
  priors <- c(TempoPriors()["theta"],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  aPrior <- list(alpha = function(a) stats::dexp(a, 1, log = TRUE))
  start <- c(mkPar, alpha = 0.8)
  Run <- function(engine, charPriors = aPrior, ...) {
    set.seed(5)
    TempoMCMC(mkTree, mkModel, start, priors, nIter = 6, burnin = 20, nSlice = 4,
              characters = mk, charPriors = charPriors,
              moves = c(node = 1, tip = 1, branch = 1, shift = 1, nodeAge = 1),
              minAcceptance = 0, engine = engine, ...)
  }
  cpp <- Run("cpp")
  ref <- Run("R")
  expect_equal(cpp$samples, ref$samples, tolerance = 1e-7)
  expect_equal(cpp$acceptance, ref$acceptance)
  expect_equal(names(cpp$samples)[1:4], c("theta", "rate", "alpha", "logPost"))
  expect_true("alpha" %in% names(cpp$acceptance))
  expect_gt(length(unique(cpp$samples$alpha)), 1)
  # logPost includes the character likelihood and the prior of alpha
  st <- cpp$state
  lengths <- .CharacterLengths(st$tree, st$set$grid$x, st$paths, st$par$rate,
                               mkModel$gamma)
  expect_equal(
    cpp$samples$logPost[[6]],
    .LogAugmented(st$tree, st$set, st$paths, st$bin) +
      MkLikelihood(st$tree, mk, lengths, alpha = st$par$alpha) +
      .LogPrior(st$par, priors) + aPrior$alpha(st$par$alpha),
    tolerance = 1e-10)
  expect_equal(st$charLL, MkLikelihood(st$tree, mk, lengths, alpha = st$par$alpha))
  # Fixed alpha: no move, the start value throughout
  fixed <- Run("cpp", charPriors = list())
  expect_equal(fixed$samples$alpha, rep(0.8, 6))
  expect_false("alpha" %in% names(fixed$acceptance))
  held <- Run("cpp", updateParameters = FALSE)
  expect_equal(held$samples$alpha, rep(0.8, 6))
  # Without rate categories, MkCharacters gives the matrix chain
  chars <- matrix(as.integer(mkCells[, c(1, 3)] == "1"), 5)
  chars[is.na(suppressWarnings(as.integer(mkCells[, c(1, 3)])))] <- NA
  RunPlain <- function(characters) {
    set.seed(8)
    TempoMCMC(mkTree, mkModel, mkPar, priors, nIter = 5, nSlice = 4,
              characters = characters,
              moves = c(node = 1, tip = 1, branch = 1, shift = 1), minAcceptance = 0)
  }
  expect_equal(RunPlain(MkCharacters(chars, nState = 2))$samples,
               RunPlain(chars)$samples, tolerance = 1e-10)
  expect_equal(TempoMoves(mkTree, mk, mkModel), TempoMoves(mkTree, chars, mkModel))
})

test_that("TempoMCMC() checks the Gamma shape and its prior", {
  mk <- MkCharacters(mkCells, nRateCat = 4)
  priors <- TempoPriors()["theta"]
  Run <- function(start = c(mkPar, alpha = 1), charPriors = list(),
                  characters = mk, pri = priors) {
    TempoMCMC(mkTree, mkModel, start, pri, nIter = 1, nSlice = 4,
              characters = characters, charPriors = charPriors,
              moves = c(node = 1), minAcceptance = 0)
  }
  expect_error(Run(start = mkPar), "start\\$alpha")
  expect_error(Run(start = c(mkPar, alpha = -1)), "start\\$alpha")
  expect_error(Run(charPriors = list(beta = function(a) 0)), "charPriors")
  expect_error(Run(charPriors = list(function(a) 0)), "charPriors")
  expect_error(Run(pri = c(priors, alpha = function(a) 0)), "goes in `charPriors`")
  expect_error(Run(charPriors = list(alpha = function(a) -Inf)), "outside the support")
  expect_error(Run(characters = MkCharacters(mkCells),
                   charPriors = list(alpha = function(a) 0)), "nRateCat")
  expect_equal(nrow(Run()$samples), 1)
  # A proposal outside the prior's support is rejected
  set <- .AddSlices(mkTree, .Precompute(mkTree, mkPar, mkModel), 4)
  st <- list(tree = mkTree, par = c(mkPar, alpha = 1), set = set, charLL = 0)
  res <- .AlphaMove(st, 1, function(a) if (a == 1) 0 else -Inf,
                    function(...) stop("not called"))
  expect_false(res$accepted)
  expect_identical(res$state, st)
})
