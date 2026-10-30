# Incremental character pruning (#26, item 3): the persistent MkPrime pruner
# behind .CharacterLikelihood() against full pruning.  Fixture as
# test-Engine.R: ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils;
# f sampled ancestor; fossil ages uncertain, so that every tree move applies.
cachePhylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
  class = "phylo")
cacheTree <- TempoTree(
  cachePhylo,
  ranges = cbind(c(0.2, 0, 0.5, 0.5, 0), c(2.5, 0, 3, 3.5, 0))
)
cachePar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.4)
cacheModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)

# Balanced tree with 2^depth extant tips and unit branch lengths
BalancedPhylo <- function(depth) {
  nTip <- 2L^depth
  edge <- matrix(0L, 0, 2)
  nextNode <- nTip + 1L
  nextTip <- 1L
  Grow <- function(node, d) {
    for (k in 1:2) {
      if (d == 1) {
        child <- nextTip
        nextTip <<- nextTip + 1L
      } else {
        nextNode <<- nextNode + 1L
        child <- nextNode
      }
      edge <<- rbind(edge, c(node, child))
      if (d > 1) Grow(child, d - 1)
    }
  }
  Grow(nTip + 1L, depth)
  structure(list(edge = edge, edge.length = rep(1, nrow(edge)),
                 tip.label = paste0("t", seq_len(nTip)), Nnode = nTip - 1L),
            class = "phylo")
}

WithCpp <- function(expr, incremental = TRUE) {
  old <- options(EvoTempo.engine = "cpp", EvoTempo.incremental = incremental)
  on.exit(options(old))
  expr
}

test_that("Incremental pruning matches full pruning through accepted and rejected moves", {
  WithCpp({
    set.seed(2601)
    tree <- cacheTree
    chars <- matrix(sample(c(0:2, NA), 5 * 7, replace = TRUE), 5)
    states <- matrix(as.integer(chars), 5)
    CL <- .CharacterLikelihood(chars, 3L, 1)
    expect_false(is.null(attr(CL, "resolve")))
    nCall <- 0
    Full <- function(tr, paths, p, set) {
      .MkLikelihoodCpp(tr, states,
                       .CharacterLengths(tr, set$grid$x, paths, p$rate, 1), 3L)
    }
    # Every evaluation, of proposals and of current states alike, is checked
    Checked <- structure(function(tr, paths, p, set) {
      value <- CL(tr, paths, p, set)
      len <- .CharacterLengths(tr, set$grid$x, paths, p$rate, 1)
      expect_identical(value, .MkLikelihoodCpp(tr, states, len, 3L))
      expect_equal(value, MkLikelihood(tr, chars, len, 3L, engine = "R"),
                   tolerance = 1e-10)
      nCall <<- nCall + 1
      value
    }, resolve = attr(CL, "resolve"))
    set <- .AddSlices(tree, .Precompute(tree, cachePar, cacheModel), 4)
    msg <- .Messages(tree, set)
    draw <- .SampleTree(tree, set, msg)
    st <- list(tree = tree, par = cachePar, set = set, msg = msg,
               paths = draw[["paths"]], bin = draw[["bin"]])
    st$charLL <- Checked(tree, st$paths, cachePar, set)
    .ResolveCharacters(Checked, TRUE)
    treeMoves <- c("nodeAge", "fossilAge", "slide", "saJump")
    moves <- c("node", "tip", "branch", "tree", treeMoves)
    accepted <- setNames(numeric(length(moves)), moves)
    rejected <- accepted
    for (i in seq_len(120)) {
      move <- sample(moves, 1)
      res <- if (move %in% treeMoves) {
        .TreeMove(move, st, Checked)
      } else {
        .PathMove(move, st, Checked)
      }
      st <- res$state
      if (res$accepted) {
        accepted[[move]] <- accepted[[move]] + 1
      } else {
        rejected[[move]] <- rejected[[move]] + 1
      }
      expect_identical(st$charLL, Full(st$tree, st$paths, st$par, st$set))
      if (i %% 10 == 0) {
        # A proposal left unresolved is diffed against at the next call
        other <- .SampleTree(st$tree, st$set, st$msg)
        Checked(st$tree, other$paths, st$par, st$set)
        expect_identical(Checked(st$tree, st$paths, st$par, st$set), st$charLL)
        .ResolveCharacters(Checked, TRUE)
      }
    }
    expect_gt(nCall, 120)
    expect_true(all(accepted[c("node", "tip", "branch", "slide")] > 0))
    expect_gt(sum(rejected), 0)
    expect_gt(sum(accepted[treeMoves]), 0)
    # A shift with `rate` fixed changes every length; a rejection restores them
    for (i in seq_len(6)) {
      st <- .ShiftMove(st, cacheModel, "lambda", TempoPriors()["lambda"], 4,
                       Checked)$state
      expect_identical(st$charLL, Full(st$tree, st$paths, st$par, st$set))
    }
  })
})

test_that("The character cache recomputes only the changed clades", {
  WithCpp({
    tree <- TempoTree(BalancedPhylo(4))
    nAll <- tree$nTip + tree$nNode
    set.seed(2602)
    chars <- matrix(sample(0:1, tree$nTip * 9, replace = TRUE), tree$nTip)
    states <- matrix(as.integer(chars), tree$nTip)
    cache <- .CharacterCacheCpp(states, 2L)
    len <- stats::runif(nAll, 0.01, 0.3)
    Eval <- function(tr, l) {
      value <- .CachedMkLikelihoodCpp(cache, tr$parent, l)
      expect_identical(value, .MkLikelihoodCpp(tr, states, l, 2L))
      value
    }
    base <- Eval(tree, len)
    expect_equal(.CharacterCacheRecomputedCpp(cache), tree$nNode)
    .ResolveCharacterCacheCpp(cache, TRUE)
    # A tip branch: the four nodes from its parent to the root
    tip <- 3L
    depth <- 0L
    v <- tip
    while (!is.na(tree$parent[[v]])) {
      depth <- depth + 1L
      v <- tree$parent[[v]]
    }
    expect_equal(depth, 4L)
    moved <- replace(len, tip, len[[tip]] * 2)
    Eval(tree, moved)
    expect_equal(.CharacterCacheRecomputedCpp(cache), depth)
    # Rejected: the journal restores the old CLs, so nothing is recomputed
    .ResolveCharacterCacheCpp(cache, FALSE)
    expect_identical(Eval(tree, len), base)
    expect_equal(.CharacterCacheRecomputedCpp(cache), 0L)
    .ResolveCharacterCacheCpp(cache, TRUE)
    # Accepted: kept
    Eval(tree, moved)
    .ResolveCharacterCacheCpp(cache, TRUE)
    Eval(tree, moved)
    expect_equal(.CharacterCacheRecomputedCpp(cache), 0L)
    # Unresolved proposals, then a change on the other side of the root
    Eval(tree, len)
    other <- replace(len, 16L, 0.5)
    Eval(tree, other)
    .ResolveCharacterCacheCpp(cache, FALSE)
    Eval(tree, other)
    .ResolveCharacterCacheCpp(cache, TRUE)
    # A regraft: tips 1 and j swap cherries; their parents and ancestors only
    Ancestors <- function(v) {
      out <- v
      while (!is.na(tree$parent[[v]])) {
        v <- tree$parent[[v]]
        out <- c(out, v)
      }
      out
    }
    p1 <- tree$parent[[1]]
    j <- which(tree$parent[seq_len(tree$nTip)] != p1)[[1]]
    pj <- tree$parent[[j]]
    swapped <- tree
    swapped$parent[c(1, j)] <- c(pj, p1)
    swapped <- .Rebuild(swapped)
    Eval(swapped, other)
    expect_equal(.CharacterCacheRecomputedCpp(cache),
                 length(union(Ancestors(p1), Ancestors(pj))))
    .ResolveCharacterCacheCpp(cache, FALSE)
    Eval(tree, other)
    expect_equal(.CharacterCacheRecomputedCpp(cache), 0L)
    Eval(swapped, other)
    .ResolveCharacterCacheCpp(cache, TRUE)
    expect_identical(Eval(swapped, other),
                     .MkLikelihoodCpp(swapped, states, other, 2L))
    # Resolving twice, or with nothing pending, is harmless
    .ResolveCharacterCacheCpp(cache, TRUE)
    .ResolveCharacterCacheCpp(cache, FALSE)
    expect_identical(Eval(swapped, other),
                     .MkLikelihoodCpp(swapped, states, other, 2L))
  })
})

test_that(".CharacterLikelihood() falls back to full pruning where it must", {
  set <- .AddSlices(cacheTree, .Precompute(cacheTree, cachePar, cacheModel), 4)
  set.seed(2603)
  draw <- .SampleTree(cacheTree, set, .Messages(cacheTree, set))
  chars <- matrix(c(0, 1, 1, 0, NA), 5)
  ref <- MkLikelihood(cacheTree, chars,
                      .CharacterLengths(cacheTree, set$grid$x, draw$paths, 0.4, 1),
                      engine = "R")
  for (engine in c("R", "cpp")) {
    for (incremental in c(TRUE, FALSE)) {
      old <- options(EvoTempo.engine = engine, EvoTempo.incremental = incremental)
      CL <- .CharacterLikelihood(chars, 2L, 1)
      expect_equal(is.null(attr(CL, "resolve")), engine == "R" || !incremental)
      expect_equal(CL(cacheTree, draw$paths, cachePar, set), ref, tolerance = 1e-10)
      expect_equal(.CharacterLikelihood(chars[, 0, drop = FALSE], 2L, 1)(
        cacheTree, draw$paths, cachePar, set), 0)
      options(old)
    }
  }
  expect_equal(.CharacterLikelihood(NULL, 2L, 1)(cacheTree, draw$paths, cachePar, set), 0)
  expect_null(.ResolveCharacters(function(...) 0, TRUE))
  WithCpp({
    CL <- .CharacterLikelihood(matrix(0, 3, 1), 2L, 1)
    expect_error(CL(cacheTree, draw$paths, cachePar, set), "one row per tip")
    expect_error(.CharacterLikelihood(matrix(2, 5, 1), 2L, 1), "n_states")
    # An external pointer does not survive serialization
    stale <- unserialize(serialize(attr(.CharacterLikelihood(chars, 2L, 1), "cache"),
                                   NULL))
    expect_error(.CachedMkLikelihoodCpp(stale, cacheTree$parent, numeric(9)),
                 "no longer valid")
  })
})

test_that("Seeded chains are identical with and without incremental pruning", {
  set.seed(2604)
  chars <- matrix(sample(0:1, 5 * 6, replace = TRUE), 5)
  priors <- c(TempoPriors()[c("lambda", "theta")],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  Run <- function(incremental) {
    WithCpp({
      set.seed(2605)
      TempoMCMC(cacheTree, cacheModel, cachePar, priors, nIter = 10, burnin = 4,
                nSlice = 4, characters = chars,
                moves = c(TempoMoves(cacheTree, chars, local = 0.6), joint = 0.5,
                          nodeAge = 1, fossilAge = 1, slide = 1, saJump = 1),
                engine = "cpp")
    }, incremental)
  }
  inc <- Run(TRUE)
  full <- Run(FALSE)
  expect_identical(inc$samples, full$samples)
  expect_identical(inc$acceptance, full$acceptance)
  expect_identical(inc$state$charLL, full$state$charLL)
  expect_gt(sum(inc$acceptance[c("node", "tip", "branch", "slide", "nodeAge")] > 0), 2)
})
