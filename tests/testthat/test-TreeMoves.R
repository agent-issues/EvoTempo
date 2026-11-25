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
    for (move in c("nodeAge", "fossilAge", "slide", "fnpr", "wfnpr", "saJump")) {
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
  expect_setequal(unique(seen), c("nodeAge", "fossilAge", "slide", "fnpr", "wfnpr", "saJump"))
  expect_true(any(st$tree$tipType == "ancestor") || any(seen == "saJump"))
})

test_that("Tree moves with no candidates return NULL", {
  cherry <- TempoTree(structure(list(edge = rbind(c(3, 1), c(3, 2)),
                                     edge.length = c(1, 1), tip.label = c("a", "b"),
                                     Nnode = 1), class = "phylo"))
  set <- .AddSlices(cherry, .Precompute(cherry, movePar, moveModel), 4)
  d <- .SampleTree(cherry, set, .Messages(cherry, set))
  st <- list(tree = cherry, set = set, paths = d$paths, bin = d$bin)
  for (move in c("nodeAge", "fossilAge", "slide", "fnpr", "wfnpr")) {
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
  set.seed(4)
  chars <- matrix(sample(0:1, 5 * 6, TRUE), 5)
  run <- TempoMCMC(moveTree, moveModel, movePar,
                   TempoPriors()[c("lambda", "theta")], nIter = 6, nSlice = 4,
                   moves = c(tree = 1, node = 1, tip = 1, branch = 1, nodeAge = 1,
                             fossilAge = 1, slide = 1, fnpr = 1, wfnpr = 1, saJump = 1, joint = 1),
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
  # Tips a-e, g, h are 1-7; all extant
  phy <- ape::read.tree(text = "((((a:1,b:1):1,c:2):2,(d:3,e:3):1):1,(g:1.5,h:1.5):3.5);")
  tr <- TempoTree(phy)
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
  lone <- TempoTree(ape::read.tree(text = "((a:0.4,(b:0.5,c:0.5):0.5):1,d:0.2);"))
  bc <- lone$parent[[2]]
  expect_length(.FnprTargets(lone, bc, 3, .Descendants(lone, 2), Inf), 0)
  # A threshold just above the parent leaves no target: no proposal
  set <- .AddSlices(tr, .Precompute(tr, movePar, moveModel), 4)
  d <- .SampleTree(tr, set, .Messages(tr, set))
  st <- list(tree = tr, set = set, paths = d$paths, bin = d$bin)
  local_mocked_bindings(.fnprSpan = 1e9)
  set.seed(1)
  nulls <- vapply(1:20, function(i) is.null(.ProposeTreeMove("fnpr", tr, set, st)), TRUE)
  expect_true(any(nulls))
})

test_that("fnpr keeps node ages and changes the topology", {
  set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, moveModel), 6)
  set.seed(5)
  d <- .SampleTree(moveTree, set, .Messages(moveTree, set))
  st <- list(tree = moveTree, set = set, paths = d$paths, bin = d$bin)
  changed <- 0
  for (i in 1:30) {
    prop <- .ProposeTreeMove("fnpr", moveTree, set, st)
    if (is.null(prop)) next
    expect_equal(sort(prop$tree$age), sort(moveTree$age))
    changed <- changed + !identical(prop$tree$parent, moveTree$parent)
  }
  expect_gt(changed, 0)
})

# State of `moveTree` with characters that depend on tempo, for wfnpr
WfnprState <- function(seed) {
  set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, moveModel), 6)
  set.seed(seed)
  msg <- .Messages(moveTree, set)
  d <- .SampleTree(moveTree, set, msg)
  chars <- matrix(sample(0:1, 5 * 4, TRUE), 5)
  CharLL <- function(tr, paths, p, s) {
    MkLikelihood(tr, chars, .CharacterLengths(tr, s$grid$x, paths, p$rate, 1))
  }
  list(tree = moveTree, par = movePar, set = set, msg = msg, paths = d$paths,
       bin = d$bin, CharLL = CharLL, charLL = CharLL(moveTree, d$paths, movePar, set))
}

test_that("wfnpr weights satisfy detailed balance between every pair of members", {
  st <- WfnprState(6)
  set <- st$set
  # Exact target of a member, and the fnpr proposal density of its four
  # redrawn branches when it is reached from the member attached at `from`
  Member <- function(m, v, from) {
    tr <- m$tree
    p <- tr$parent[[v]]
    lw <- .NodeLogWeights(tr, set, p, set$lambda, m$bin[[tr$parent[[p]]]],
                          c(m$y, v), m$bin[c(m$y, v)])
    c(target = .LogAugmented(tr, set, m$paths, m$bin) +
        st$CharLL(tr, m$paths, movePar, set),
      q = .LogBridges(tr, set, c(p, m$y, v, from), m$paths) + lw[[m$bin[[p]]]] -
        .LogSumExp(lw))
  }
  nMember <- integer(0)
  # v: fossil tip a (1), tips c and d (4, 5), sampled-ancestor node (8)
  for (v in .SlideNodes(moveTree)) {
    p <- moveTree$parent[[v]]
    sib <- setdiff(moveTree$kids[p - moveTree$nTip, ], v)
    targets <- sort(.FnprTargets(moveTree, p, sib, .Descendants(moveTree, v), Inf))
    members <- c(
      list(list(y = sib, tree = moveTree, bin = st$bin, paths = st$paths, logW = 0)),
      .WfnprEnsemble(moveTree, set, st, v, targets)
    )
    nMember <- c(nMember, length(members))
    for (a in seq_along(members)) {
      for (b in seq_along(members)[-a]) {
        # The shared path on the old sibling's branch puts every member one
        # fnpr move away from every other
        ma <- Member(members[[a]], v, members[[b]]$y)
        mb <- Member(members[[b]], v, members[[a]]$y)
        expect_equal(members[[b]]$logW - members[[a]]$logW,
                     mb[["target"]] - ma[["target"]] + ma[["q"]] - mb[["q"]],
                     tolerance = 1e-8)
      }
    }
  }
  expect_gt(max(nMember), 2)
})

test_that("wfnpr bounds its candidates and rejects when none has weight", {
  tr <- TempoTree(ape::read.tree(
    text = "((((a:1,b:1):1,c:2):2,(d:3,e:3):1):1,(g:1.5,h:1.5):3.5);"
  ))
  set <- .AddSlices(tr, .Precompute(tr, movePar, moveModel), 4)
  set.seed(7)
  d <- .SampleTree(tr, set, .Messages(tr, set))
  st <- list(tree = tr, set = set, paths = d$paths, bin = d$bin)
  nTarget <- integer(0)
  Ensemble <- .WfnprEnsemble
  local_mocked_bindings(
    .wfnprMax = 3L,
    .WfnprEnsemble = function(tree, set, st, v, targets) {
      nTarget <<- c(nTarget, length(targets))
      expect_false(is.unsorted(targets))
      Ensemble(tree, set, st, v, targets)
    }
  )
  props <- lapply(1:40, function(i) .ProposeTreeMove("wfnpr", tr, set, st))
  expect_equal(max(nTarget), 2)
  expect_true(any(vapply(props, is.null, TRUE)))
  kept <- Filter(Negate(is.null), props)
  expect_true(all(vapply(kept, function(prop) {
    is.finite(prop$logRatio) && identical(sort(prop$tree$age), sort(tr$age)) &&
      !identical(prop$tree$parent, tr$parent)
  }, TRUE)))
  # Candidates of zero target density are never chosen
  local_mocked_bindings(.CandidateLogLik = function(...) -Inf)
  full <- WfnprState(8)
  set.seed(9)
  for (i in 1:10) {
    expect_null(.ProposeTreeMove("wfnpr", moveTree, full$set, full))
    expect_false(.TreeMove("wfnpr", full, full$CharLL)$accepted)
  }
})

test_that("wfnpr samples topologies with their exact probabilities", {
  # Ages and parameters fixed, characters independent of tempo: the topology
  # has probability proportional to TreeLikelihood() times the Mk likelihood
  model <- DCETModel(nBin = 5, zRange = c(-2, 2), gamma = 0, rho = 0.7)
  set.seed(11)
  chars <- matrix(sample(0:1, 5 * 2, TRUE), 5)
  Mk <- function(tr) {
    nodes <- .PathNodes(tr)
    len <- numeric(9)
    len[nodes] <- movePar$rate * (tr$age[tr$parent[nodes]] - tr$age[nodes])
    MkLikelihood(tr, chars, len)
  }
  Key <- function(tr) paste(tr$parent, collapse = " ")
  # Every topology that fixed-height regrafts reach from `moveTree`
  trees <- list(moveTree)
  names(trees) <- Key(moveTree)
  todo <- trees
  while (length(todo)) {
    tr <- todo[[1]]
    todo <- todo[-1]
    for (v in .SlideNodes(tr)) {
      p <- tr$parent[[v]]
      sib <- setdiff(tr$kids[p - tr$nTip, ], v)
      for (y in .FnprTargets(tr, p, sib, .Descendants(tr, v), Inf)) {
        new <- tr
        new$parent[sib] <- tr$parent[[p]]
        new$parent[p] <- new$parent[[y]]
        new$parent[y] <- p
        new <- .Rebuild(new)
        if (is.null(trees[[Key(new)]])) {
          trees[[Key(new)]] <- todo[[Key(new)]] <- new
        }
      }
    }
  }
  expect_gt(length(trees), 4)
  logP <- vapply(trees, function(tr) TreeLikelihood(tr, movePar, model) + Mk(tr), 0)
  exact <- exp(logP - .LogSumExp(logP))
  CharLL <- function(tr, paths, p, s) Mk(tr)
  set <- .AddSlices(moveTree, .Precompute(moveTree, movePar, model), 4)
  msg <- .Messages(moveTree, set)
  set.seed(12)
  d <- .SampleTree(moveTree, set, msg)
  st <- list(tree = moveTree, par = movePar, set = set, msg = msg, paths = d$paths,
             bin = d$bin, charLL = Mk(moveTree))
  nIter <- 1500
  visits <- character(nIter)
  accepted <- 0
  for (i in seq_len(nIter)) {
    st <- .PathMove("tree", st, CharLL)$state
    res <- .TreeMove("wfnpr", st, CharLL)
    st <- res$state
    accepted <- accepted + res$accepted
    visits[[i]] <- Key(st$tree)
  }
  expect_true(all(visits %in% names(trees)))
  expect_gt(accepted / nIter, 0.2)
  observed <- table(factor(visits, names(trees))) / nIter
  # Four standard errors, with the variance inflated threefold for autocorrelation
  possible <- exact > 0
  expect_true(any(possible) && all(observed[!possible] == 0))
  z <- (observed - exact) / sqrt(3 * exact * (1 - exact) / nIter)
  expect_lt(max(abs(z[possible])), 4)
})

test_that("wfnpr leaves the prior over trees invariant", {
  # No characters, fixed parameters: tree statistics with and without the move
  Run <- function(seed, moves) {
    set.seed(seed)
    TempoMCMC(moveTree, DCETModel(nBin = 5, zRange = c(-2, 2), rho = 0.7), movePar,
              list(), nIter = 1200, burnin = 100, thin = 2, nSlice = 4,
              moves = c(tree = 1, nodeAge = 1, fossilAge = 1, saJump = 1, moves),
              minAcceptance = 0)
  }
  with <- Run(13, c(wfnpr = 2))
  without <- Run(14, c(slide = 1, fnpr = 1))
  expect_gt(with$acceptance[["wfnpr"]], 0.1)
  for (stat in c("nAncestor", "cherries", "treeLength", "youngestFossil")) {
    # Batch means: 20 batches of 30 draws
    Batch <- function(run) colMeans(matrix(run$samples[[stat]], 30))
    a <- Batch(with)
    b <- Batch(without)
    z <- (mean(a) - mean(b)) / sqrt(var(a) / 20 + var(b) / 20)
    expect_lt(abs(z), 3.5)
  }
})
