# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
movePhylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4), class = "phylo")
moveRanges <- cbind(c(0.2, 0, 0.5, 0.5, 0), c(2.5, 0, 3, 3.5, 0))
moveTree <- TempoTree(movePhylo, ranges = moveRanges)
movePar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.3)
moveModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)

test_that("MkLikelihood() matches a two-tip closed form", {
  cherry <- structure(list(edge = rbind(c(3, 1), c(3, 2)), edge.length = c(1, 1),
                           tip.label = c("a", "b"), Nnode = 1), class = "phylo")
  ell <- c(0.3, 0.5, 0)
  for (k in 2:4) {
    stay <- function(l) 1 / k + (k - 1) / k * exp(-k * l / (k - 1))
    same <- stay(0.8)
    expect_equal(MkLikelihood(cherry, matrix(c(0, 0, 1, 0), 2), ell, k),
                 log(same / k) + log((1 - same) / (k - 1) / k))
  }
  expect_equal(MkLikelihood(cherry, matrix(c(NA, 1), 2), ell), log(1 / 2))
  expect_error(MkLikelihood(cherry, matrix(0, 3, 1), ell), "one row per tip")
})

test_that("Product-space likelihood factorizes when tempo is irrelevant", {
  character <- c(0, 1, 1, NA, 0)
  m0 <- DCETModel(nBin = 11, zRange = c(-2, 2), beta = 0, gamma = 0, rho = 0.7)
  nodes <- .PathNodes(moveTree)
  len <- numeric(9)
  len[nodes] <- 0.3 * (moveTree$age[moveTree$parent[nodes]] - moveTree$age[nodes])
  expect_equal(TreeLikelihood(moveTree, movePar, m0, character = character),
               FBDTreeLikelihood(moveTree, 0.5, 0.2, 0.1, 0.7) +
                 MkLikelihood(moveTree, matrix(character), len),
               tolerance = 1e-8)
  expect_error(TreeLikelihood(moveTree, movePar, m0, character = 1:2), "one state")
  # With tempo, the product-space likelihood equals the expectation over paths
  m1 <- moveModel
  set.seed(1)
  maps <- lapply(1:300, function(i) {
    set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, m1), 6)
    msg <- .Messages(moveTree, set)
    d <- .SampleTree(moveTree, set, msg)
    lw <- .LogWeight(moveTree, set, msg, d$paths, .PathNodes(moveTree))
    c(lw, MkLikelihood(moveTree, matrix(character),
                       .CharacterLengths(moveTree, set$grid$x, d$paths, 0.3, 1)),
      .FrozenLikelihood(moveTree, set, msg))
  })
  m <- do.call(rbind, maps)
  w <- exp(m[, 1] + m[, 2])
  estimate <- m[1, 3] + log(mean(w))
  expect_equal(estimate, TreeLikelihood(moveTree, movePar, m1, character = character),
               tolerance = 4 * sd(w) / sqrt(300) / mean(w))
})

test_that("Rebuild and helpers keep the tree consistent", {
  rebuilt <- .Rebuild(moveTree)
  expect_equal(rebuilt$kids, moveTree$kids)
  expect_equal(rebuilt$postorder, moveTree$postorder)
  expect_equal(sort(.Descendants(moveTree, 7)), c(1, 2, 3, 7, 8))
  expect_equal(sort(.SlideNodes(moveTree)), c(1, 4, 5, 8))
  expect_equal(sort(.Crossing(moveTree, 2.5, 9, 4)), c(1, 5, 8))
  expect_equal(.LogSumExp(c(log(2), log(3))), log(5))
  set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, moveModel), 6)
  pc <- .Pieces(moveTree, set, 7)
  expect_equal(.LogTransferRow(pc, 3)[5], .LogTransferCol(pc, 5)[3], tolerance = 1e-10)
  expect_error(TempoTree(movePhylo, ranges = cbind(c(2, 0, 0, 0, 0), 3)), "within its range")
})

test_that("Each tree move proposes valid trees", {
  set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, moveModel), 6)
  msg <- .Messages(moveTree, set)
  set.seed(2)
  d <- .SampleTree(moveTree, set, msg)
  st <- list(tree = moveTree, set = set, paths = d$paths, bin = d$bin)
  seen <- character(0)
  for (i in 1:60) {
    for (move in c("nodeAge", "fossilAge", "slide", "saJump")) {
      prop <- .ProposeTreeMove(move, st$tree, set, st)
      if (is.null(prop)) next
      seen <- c(seen, move)
      tr <- prop$tree
      expect_true(is.finite(prop$logHastings))
      expect_true(all(tr$age[tr$parent[.PathNodes(tr)]] >= tr$age[.PathNodes(tr)]))
      expect_equal(tr$age[[tr$root]], moveTree$age[[moveTree$root]])
      for (v in .PathNodes(tr)) {
        path <- prop$paths[[v]]
        expect_equal(path$k[[1]], prop$bin[[tr$parent[[v]]]])
        expect_equal(path$k[[length(path$k)]], prop$bin[[v]])
      }
      tips <- tr$age[seq_len(tr$nTip)]
      expect_true(all(tips >= tr$range[, 1] - 1e-12 & tips <= tr$range[, 2] + 1e-12))
      if (stats::runif(1) < 0.5) {
        st$tree <- tr
        st$paths <- prop$paths
        st$bin <- prop$bin
      }
    }
  }
  expect_setequal(unique(seen), c("nodeAge", "fossilAge", "slide", "saJump"))
  expect_true(any(st$tree$tipType == "ancestor") || any(seen == "saJump"))
})

test_that("Tree moves with no candidates return NULL", {
  cherry <- TempoTree(structure(list(edge = rbind(c(3, 1), c(3, 2)),
                                     edge.length = c(1, 1), tip.label = c("a", "b"),
                                     Nnode = 1), class = "phylo"))
  set <- .AddSlices(cherry, .Precompute(cherry, movePar, moveModel), 4)
  d <- .SampleTree(cherry, set, .Messages(cherry, set))
  st <- list(tree = cherry, set = set, paths = d$paths, bin = d$bin)
  for (move in c("nodeAge", "fossilAge", "slide")) {
    expect_null(.ProposeTreeMove(move, cherry, set, st))
  }
  # Fossil tip whose parent is the root
  rootFossil <- TempoTree(structure(list(edge = rbind(c(4, 1), c(4, 5), c(5, 2), c(5, 3)),
                                         edge.length = c(1, 1, 2, 2),
                                         tip.label = c("f", "a", "b"), Nnode = 2),
                                    class = "phylo"))
  set <- .AddSlices(rootFossil, .Precompute(rootFossil, movePar, moveModel), 4)
  d <- .SampleTree(rootFossil, set, .Messages(rootFossil, set))
  st <- list(tree = rootFossil, set = set, paths = d$paths, bin = d$bin)
  expect_null(.ProposeTreeMove("saJump", rootFossil, set, st))
})

test_that("saJump leaves a tree without fossils unchanged", {
  extantOnly <- structure(
    list(
      edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
      edge.length = c(1, 2, 2, 3),
      tip.label = c("a", "b", "c"),
      Nnode = 2
    ),
    class = "phylo"
  )
  set.seed(1)
  run <- TempoMCMC(extantOnly, DCETModel(nBin = 11, zRange = c(-2, 2)),
                   list(lambda = 0.5, mu = 0.2, theta = 0.1, s = 1),
                   TempoPriors()[c("lambda", "theta")], nIter = 1, nSlice = 4,
                   moves = c(saJump = 1))
  st <- run$state
  expect_null(.ProposeTreeMove("saJump", st$tree, st$set, st))
  expect_identical(run$state$tree, TempoTree(extantOnly))
  expect_equal(run$acceptance[["saJump"]], 0)
})

test_that("TempoMCMC() runs every move with characters", {
  set.seed(3)
  chars <- matrix(sample(0:1, 5 * 6, TRUE), 5)
  run <- TempoMCMC(moveTree, moveModel, movePar,
                   TempoPriors()[c("lambda", "theta")], nIter = 6, nSlice = 4,
                   moves = c(tree = 1, node = 1, tip = 1, branch = 1, nodeAge = 1,
                             fossilAge = 1, slide = 1, saJump = 1, joint = 1),
                   characters = chars)
  expect_true(all(run$acceptance > 0))
  expect_true(all(is.finite(run$samples$logPost)))
  expect_true(all(run$samples$clock > 0))
  expect_error(TempoMCMC(moveTree, moveModel, movePar, list(), nIter = 1,
                         moves = c(jump = 1)), "Unknown move")
  stats <- .TreeStats(moveTree)
  expect_equal(stats[["nAncestor"]], 1)
  expect_equal(stats[["cherries"]], 1)
  expect_true(is.na(stats[["clock"]]))
  extantOnly <- TempoTree(structure(list(edge = rbind(c(3, 1), c(3, 2)),
                                         edge.length = c(1, 1), tip.label = c("a", "b"),
                                         Nnode = 1), class = "phylo"))
  expect_true(is.na(.TreeStats(extantOnly)[["oldestFossil"]]))
  expect_true(is.na(.TreeStats(extantOnly)[["meanNodeAge"]]))
})

test_that("MarginalMCMC() accepts one character", {
  set.seed(4)
  out <- MarginalMCMC(moveTree, moveModel, movePar,
                      c(TempoPriors()["lambda"],
                        rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE)),
                      nIter = 2, character = c(0, 1, 1, NA, 0))
  expect_equal(dim(out), c(2, 3))
})

test_that("Tree moves with infeasible bridges are rejected", {
  tiny <- TempoTree(structure(list(edge = rbind(c(5, 6), c(6, 7), c(7, 1), c(7, 2),
                                                c(6, 3), c(5, 4)),
                                   edge.length = c(1, 1, 1, 1, 2, 3) * 1e-7,
                                   tip.label = c("a", "b", "c", "d"), Nnode = 3),
                              class = "phylo"))
  set <- .AddSlices(tiny, .Precompute(tiny, movePar, moveModel), 2)
  bin <- c(11L, 11L, 1L, 1L, 1L, 11L, 1L)
  st <- list(tree = tiny, par = movePar, set = set, msg = .Messages(tiny, set),
             paths = list(), bin = bin, charLL = 0)
  res <- .TreeMove("nodeAge", st, function(...) 0)
  expect_false(res$accepted)
})

test_that("Joint parameter-path moves run and adapt", {
  set.seed(4)
  chars <- matrix(sample(0:1, 5 * 3, TRUE), 5)
  priors <- c(TempoPriors()["theta"],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  run <- TempoMCMC(moveTree, moveModel, c(movePar, rate = 0.3), priors,
                   nIter = 4, burnin = 20, nSlice = 4,
                   moves = c(tree = 1, joint = 1), characters = chars)
  expect_equal(nrow(run$samples), 4)
  expect_gt(run$acceptance[["joint"]], 0)
  # Proposals outside the prior support are skipped
  bounded <- list(theta = function(x) if (x > 1e6) 0 else -Inf)
  fixed <- TempoMCMC(moveTree, moveModel, movePar, bounded, nIter = 2,
                     nSlice = 4, moves = c(joint = 1))
  expect_equal(fixed$acceptance[["joint"]], 0)
  # No joint moves when parameters are held fixed
  held <- TempoMCMC(moveTree, moveModel, movePar, TempoPriors()["theta"], nIter = 2,
                    nSlice = 4, moves = c(joint = 1), updateParameters = FALSE)
  expect_equal(held$acceptance[["joint"]], 0)
  # Infeasible FFBS draws are rejected
  set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, moveModel), 4)
  msg <- .Messages(moveTree, set)
  d <- .SampleTree(moveTree, set, msg)
  st <- list(tree = moveTree, par = movePar, set = set, msg = msg, paths = d$paths,
             bin = d$bin, charLL = 0)
  res <- .JointSweep(st, TempoPriors()["theta"], c(theta = 1), moveModel, 4,
                     function(...) 0, Draw = function(...) .Infeasible())
  expect_false(res$accepted[["theta"]])
})
