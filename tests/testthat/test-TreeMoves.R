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

# ((((a:1,b:1):1,c:2):2,(d:3,e:3):1):1,(g:1.5,h:1.5):3.5); tips a-e, g, h are 1-7
sevenPhylo <- structure(
  list(edge = rbind(c(8, 9), c(9, 10), c(10, 11), c(11, 1), c(11, 2), c(10, 3),
                    c(9, 12), c(12, 4), c(12, 5), c(8, 13), c(13, 6), c(13, 7)),
       edge.length = c(1, 2, 1, 1, 1, 2, 1, 3, 3, 3.5, 1.5, 1.5),
       tip.label = c("a", "b", "c", "d", "e", "g", "h"), Nnode = 6), class = "phylo")

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
    for (move in c("nodeAge", "fossilAge", "slide", "fnpr", "saJump")) {
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
  expect_setequal(unique(seen), c("nodeAge", "fossilAge", "slide", "fnpr", "saJump"))
  expect_true(any(st$tree$tipType == "ancestor") || any(seen == "saJump"))
})

test_that("Tree moves with no candidates return NULL", {
  cherry <- TempoTree(structure(list(edge = rbind(c(3, 1), c(3, 2)),
                                     edge.length = c(1, 1), tip.label = c("a", "b"),
                                     Nnode = 1), class = "phylo"))
  set <- .AddSlices(cherry, .Precompute(cherry, movePar, moveModel), 4)
  d <- .SampleTree(cherry, set, .Messages(cherry, set))
  st <- list(tree = cherry, set = set, paths = d$paths, bin = d$bin)
  for (move in c("nodeAge", "fossilAge", "slide", "fnpr")) {
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
                             fossilAge = 1, slide = 1, fnpr = 1, saJump = 1, joint = 1),
                   characters = chars)
  expect_true(all(run$acceptance[c(names(run$moves), "parameters")] > 0))
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
  expect_gt(run$burninAcceptance[["joint"]], 0)
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

test_that("fnpr targets are symmetric and local", {
  # All extant
  tr <- TempoTree(sevenPhylo)
  nTip <- tr$nTip
  ab <- tr$parent[[1]]
  abc <- tr$parent[[ab]]
  gh <- tr$parent[[6]]
  # Prune c (parent abc, sibling ab): lineages crossing age(abc) = 2 in the
  # pruned tree are ab, d, e and gh; ab's MRCA with d and e is 4, with gh 5
  clade <- .Descendants(tr, 3)
  expect_setequal(.FnprTargets(tr, abc, ab, clade, Inf), c(4, 5, gh))
  expect_setequal(.FnprTargets(tr, abc, ab, clade, 4), c(4, 5))
  expect_length(.FnprTargets(tr, abc, ab, clade, 3.9), 0)
  pruned <- .PrunedParent(tr, abc)
  expect_equal(pruned[[ab]], tr$parent[[abc]])
  expect_equal(.MRCAAge(tr, pruned, ab, c(4, gh)), c(4, 5))
  # Regraft onto d: the reverse target set, around d, has the same size
  newTree <- tr
  newTree$parent[ab] <- tr$parent[[abc]]
  newTree$parent[abc] <- tr$parent[[4]]
  newTree$parent[4] <- abc
  newTree <- .Rebuild(newTree)
  for (h in c(4, 4.5, 5, Inf)) {
    expect_equal(length(.FnprTargets(newTree, abc, 4, clade, h)),
                 length(.FnprTargets(tr, abc, ab, clade, h)))
  }
  # No other lineage alive at age(p): a and d are fossils older than (b, c)
  # ((a:0.4,(b:0.5,c:0.5):0.5):1,d:0.2);
  lone <- TempoTree(structure(
    list(edge = rbind(c(5, 6), c(6, 1), c(6, 7), c(7, 2), c(7, 3), c(5, 4)),
         edge.length = c(1, 0.4, 0.5, 0.5, 0.5, 0.2),
         tip.label = c("a", "b", "c", "d"), Nnode = 3), class = "phylo"))
  bc <- lone$parent[[2]]
  expect_length(.FnprTargets(lone, bc, 3, .Descendants(lone, 2), Inf), 0)
  # With no target there is no proposal
  set <- .AddSlices(tr, .Precompute(tr, movePar, moveModel), 4)
  d <- .SampleTree(tr, set, .Messages(tr, set))
  st <- list(tree = tr, set = set, paths = d$paths, bin = d$bin)
  local_mocked_bindings(.FnprTargets = function(...) integer(0))
  expect_null(.ProposeTreeMove("fnpr", tr, set, st))
})

test_that("slide with no lineage at the regraft height proposes nothing", {
  set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, moveModel), 4)
  d <- .SampleTree(moveTree, set, .Messages(moveTree, set))
  st <- list(tree = moveTree, set = set, paths = d$paths, bin = d$bin)
  local_mocked_bindings(.Crossing = function(...) integer(0))
  expect_null(.ProposeTreeMove("slide", moveTree, set, st))
})

# Fixtures for fnpr: moveTree with a crown start, and with an origin 0.6 above
FnprStates <- function(seed) {
  originModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, start = "origin")
  stemmed <- structure(c(unclass(movePhylo), root.edge = 0.6), class = "phylo")
  fixtures <- list(
    list(tree = moveTree, model = moveModel),
    list(tree = .ForModel(TempoTree(stemmed, ranges = moveRanges), originModel),
         model = originModel)
  )
  lapply(fixtures, function(f) {
    set <- .AddSlices(f$tree, .Precompute(f$tree, movePar, f$model), 6)
    set.seed(seed)
    d <- .SampleTree(f$tree, set, .Messages(f$tree, set))
    list(tree = f$tree, model = f$model, set = set, paths = d$paths, bin = d$bin)
  })
}

# Clades of a tree, as sorted tip sets
Clades <- function(tr) {
  sort(vapply(tr$postorder, function(u) {
    tips <- .Descendants(tr, u)
    paste(sort(tips[tips <= tr$nTip]), collapse = ",")
  }, ""))
}

test_that("fnpr Hastings ratio matches the proposal densities in both directions", {
  # Recomputes log q(reverse) - log q(forward) from the definition of the move
  # (dev/plans/local-fnpr.md): choice of v, uniform target among the targets
  # at the drawn threshold, bin of p from its frozen conditional, and frozen
  # bridges on the four branches whose top changes.  A wrong or missing term
  # in .ProposeTreeMove()'s ratio fails this test
  nChecked <- 0
  for (st in FnprStates(7)) {
    set <- st$set
    for (i in 1:40) {
      prop <- .ProposeTreeMove("fnpr", st$tree, set, st)
      if (is.null(prop)) next
      old <- st$tree
      new <- prop$tree
      v <- prop$fnpr$v
      y <- prop$fnpr$y
      h <- prop$fnpr$threshold
      p <- old$parent[[v]]
      pp <- old$parent[[p]]
      sib <- setdiff(old$kids[p - old$nTip, ], v)
      x <- old$parent[[y]]
      clade <- .Descendants(old, v)
      expect_equal(new$parent[c(p, y, sib)], c(x, p, pp))
      expect_equal(new$age, old$age)
      fwdT <- .FnprTargets(old, p, sib, clade, h)
      revT <- .FnprTargets(new, p, y, clade, h)
      expect_true(y %in% fwdT)
      expect_true(sib %in% revT)
      BinQ <- function(tr, kTop, kids, bins, k) {
        lw <- .NodeLogWeights(tr, set, p, set$lambda, kTop, kids, bins[kids])
        lw[[k]] - .LogSumExp(lw)
      }
      Bridges <- function(tr, paths) {
        sum(vapply(c(p, y, v, sib), function(u) .LogBridge(tr, set, u, paths[[u]]), 0))
      }
      logFwd <- -log(length(.SlideNodes(old))) - log(length(fwdT)) +
        BinQ(new, st$bin[[x]], c(y, v), prop$bin, prop$bin[[p]]) +
        Bridges(new, prop$paths)
      logRev <- -log(length(.SlideNodes(new))) - log(length(revT)) +
        BinQ(old, st$bin[[pp]], c(sib, v), st$bin, st$bin[[p]]) +
        Bridges(old, st$paths)
      expect_equal(prop$logHastings, logRev - logFwd, tolerance = 1e-8)
      expect_false(identical(Clades(new), Clades(old)))
      nChecked <- nChecked + 1
      # Move on half the time, so later checks start from other states
      if (i %% 2 == 0) {
        st$tree <- new
        st$paths <- prop$paths
        st$bin <- prop$bin
      }
    }
  }
  expect_gt(nChecked, 20)
})

test_that("fnpr runs from an origin start", {
  st <- FnprStates(8)[[2]]
  set.seed(9)
  run <- TempoMCMC(st$tree, st$model, movePar, TempoPriors()["lambda"], nIter = 20,
                   nSlice = 4, moves = c(node = 1, fnpr = 2), minAcceptance = 0)
  expect_gt(run$acceptance[["fnpr"]], 0)
  expect_equal(run$state$tree$age, st$tree$age)
  expect_equal(.StartAge(run$state$tree), .StartAge(st$tree))
})
