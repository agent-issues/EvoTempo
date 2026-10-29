# Origin start (stem path) and moves of the start age (#18).
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7) with a stem of 0.6 above the
# root; a, c fossils; f sampled ancestor
orPhylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
  class = "phylo")
orRanged <- TempoTree(orPhylo,
                      ranges = cbind(c(0.2, 0, 0.5, 0.5, 0), c(2.5, 0, 3, 3.5, 0)))
orPar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.4)
orModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, start = "origin")
crModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)
orTree <- .ForModel(TempoTree(orPhylo), orModel)
set.seed(41)
orChars <- matrix(sample(0:1, 5 * 4, replace = TRUE), 5)
orPriors <- c(TempoPriors()[c("lambda", "theta")],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
nMCo <- if (identical(Sys.getenv("R_COVR"), "true")) 150 else 3000

WithEngineO <- function(engine, expr) {
  old <- options(EvoTempo.engine = engine)
  on.exit(options(old))
  expr
}

test_that(".ForModel() marks the origin and checks the stem", {
  expect_true(.Origin(orTree))
  expect_false(.Origin(.ForModel(TempoTree(orPhylo), crModel)))
  expect_equal(.StartAge(orTree), 4.6)
  expect_equal(.PathNodes(orTree), c(1, 2, 4, 5, 6, 7, 8, 9))
  expect_equal(.TopAge(orTree, c(6, 7)), c(4.6, 4))
  noStem <- orPhylo
  noStem$root.edge <- NULL
  expect_error(.ForModel(TempoTree(noStem), orModel), "positive root edge")
  # A sampled ancestor at the root
  saRoot <- TempoTree(structure(list(edge = rbind(c(4, 1), c(4, 5), c(5, 2), c(5, 3)),
                                     edge.length = c(0, 1, 2, 2),
                                     tip.label = c("s", "b", "c"), Nnode = 2,
                                     root.edge = 1), class = "phylo"))
  expect_error(.ForModel(saRoot, orModel), "sampled ancestor at the root")
  expect_error(TempoMCMC(saRoot, orModel, orPar, list(), nIter = 1), "sampled ancestor")
})

test_that("The stem path recovers the origin likelihood", {
  set <- .AddSlices(orTree, .Precompute(orTree, orPar, orModel), 8)
  expect_equal(max(set$bounds), 4.6)
  expect_equal(set$survive,
               1 - ExtinctionProbability(4.6, set$Q, set$lambda, set$mu, 0, 0.7)[1, ],
               tolerance = 1e-8)
  msg <- .Messages(orTree, set)
  expect_length(msg$top, 11)
  exact <- TreeLikelihood(orTree, orPar, orModel)
  expect_lt(abs(exact - .FrozenLikelihood(orTree, set, msg)), 0.01)
  set.seed(2)
  maps <- SampleTempoPaths(orPhylo, orPar, orModel, n = 60, nSlice = 8)
  expect_equal(maps$logLik, exact, tolerance = 5 * sd(exp(maps$logWeight)) / sqrt(60))
  expect_equal(ncol(maps$branchX), 8)
  draw <- .SampleTree(orTree, set, msg)
  expect_length(draw$bin, 10)
  stem <- draw$paths[[6]]
  expect_equal(stem$k[[1]], draw$bin[[10]])
  expect_equal(stem$k[[length(stem$k)]], draw$bin[[6]])
  expect_true(all(diff(c(4.6, stem$t, 4)) <= 0))
  expect_true(is.finite(.LogAugmented(orTree, set, draw$paths, draw$bin)))
  # The stem's tempo-weighted length is that of the root, which pruning ignores
  len <- .CharacterLengths(orTree, set$grid$x, draw$paths, 1, 0)
  expect_equal(len[[6]], 0.6)
  # Joint conditioning
  joint <- orModel
  joint$joint <- TRUE
  setJ <- .AddSlices(orTree, .Precompute(orTree, orPar, joint), 8)
  msgJ <- .Messages(orTree, setJ)
  expect_lt(abs(TreeLikelihood(orTree, orPar, joint) -
                  .FrozenLikelihood(orTree, setJ, msgJ, TRUE)), 0.01)
  # Extant conditioning with an origin needs one extant tip anywhere
  fossilOnly <- orTree
  fossilOnly$tipType[c(2, 5)] <- "fossil"
  expect_true(.Compatible(fossilOnly, "survival", TRUE))
  expect_false(.Compatible(fossilOnly, "extant", TRUE))
  expect_equal(.LogAugmented(fossilOnly, set, draw$paths, draw$bin), -Inf)
  expect_equal(WithEngineO("cpp", .LogAugmented(fossilOnly, set, draw$paths, draw$bin)),
               -Inf)
})

test_that("Path moves with an origin leave the exact conditional law invariant", {
  set.seed(3)
  run <- TempoMCMC(orPhylo, orModel, orPar, list(), nIter = nMCo / 75, nSlice = 6,
                   moves = c(node = 3, tip = 2, branch = 2))
  expect_true(all(run$acceptance[c("node", "tip", "branch")] > 0.9))
  st <- run$state
  for (v in .PathNodes(orTree)) {
    path <- st$paths[[v]]
    expect_equal(path$k[[1]], .TopBin(orTree, st$bin, v))
    expect_equal(path$k[[length(path$k)]], st$bin[[v]])
  }
  expect_equal(unique(run$samples$originAge), 4.6)
  maps <- SampleTempoPaths(orPhylo, orPar, orModel, n = nMCo / 30, nSlice = 6)
  w <- exp(maps$logWeight - max(maps$logWeight))
  expect_lt(abs(mean(run$samples$meanX) - sum(w * maps$meanX) / sum(w)), 1)
})

test_that("Origin kernels agree between engines", {
  set <- .AddSlices(orTree, .Precompute(orTree, orPar, orModel), 6)
  msgR <- WithEngineO("R", .Messages(orTree, set))
  msgC <- WithEngineO("cpp", .Messages(orTree, set))
  expect_equal(msgC, msgR, tolerance = 1e-10)
  expect_equal(WithEngineO("cpp", .Pieces(orTree, set, 6)),
               WithEngineO("R", .Pieces(orTree, set, 6)), tolerance = 1e-10)
  for (seed in 1:3) {
    set.seed(seed)
    dR <- WithEngineO("R", .SampleTree(orTree, set, msgR))
    set.seed(seed)
    dC <- WithEngineO("cpp", .SampleTree(orTree, set, msgR))
    expect_identical(dC$bin, as.integer(dR$bin))
    expect_equal(dC$paths, dR$paths, tolerance = 1e-9)
    expect_equal(WithEngineO("cpp", .LogAugmented(orTree, set, dR$paths, dR$bin)),
                 WithEngineO("R", .LogAugmented(orTree, set, dR$paths, dR$bin)),
                 tolerance = 1e-10)
    nodes <- .PathNodes(orTree)
    expect_equal(WithEngineO("cpp", .LogWeight(orTree, set, msgR, dR$paths, nodes)),
                 WithEngineO("R", .LogWeight(orTree, set, msgR, dR$paths, nodes)),
                 tolerance = 1e-10)
    expect_equal(WithEngineO("cpp", .CharacterLengths(orTree, set$grid$x, dR$paths, 0.4, 1)),
                 WithEngineO("R", .CharacterLengths(orTree, set$grid$x, dR$paths, 0.4, 1)),
                 tolerance = 1e-12)
    expect_equal(WithEngineO("cpp", .LogBridges(orTree, set, c(6, 7), dR$paths)),
                 WithEngineO("R", .LogBridges(orTree, set, c(6, 7), dR$paths)),
                 tolerance = 1e-10)
    expect_equal(WithEngineO("cpp", .LogPathExact(orTree, set, 6, dR$paths[[6]])),
                 WithEngineO("R", .LogPathExact(orTree, set, 6, dR$paths[[6]])),
                 tolerance = 1e-10)
    set.seed(seed)
    bR <- WithEngineO("R", .DrawBridges(orTree, set, c(6, 9), dR$bin, dR$paths))
    set.seed(seed)
    bC <- WithEngineO("cpp", .DrawBridges(orTree, set, c(6, 9), dR$bin, dR$paths))
    expect_equal(bC, bR, tolerance = 1e-9)
    set.seed(seed)
    cR <- WithEngineO("R", .SampleClade(orTree, set, msgR, dR$paths, dR$bin, 7, dR$bin[[6]]))
    set.seed(seed)
    cC <- WithEngineO("cpp", .SampleClade(orTree, set, msgR, dR$paths, dR$bin, 7, dR$bin[[6]]))
    expect_identical(cC$bin, as.integer(cR$bin))
    expect_length(cC$bin, 10)
  }
  ranged <- .ForModel(orRanged, orModel)
  expect_identical(WithEngineO("cpp", .Rebuild(ranged)), WithEngineO("R", .Rebuild(ranged)))
})

test_that("Chains with an origin and age moves match between engines", {
  Run <- function(engine, model, tree, moves, agePriors) {
    set.seed(14)
    TempoMCMC(tree, model, orPar, orPriors, nIter = 10, burnin = 4, nSlice = 4,
              characters = orChars, moves = moves, agePriors = agePriors,
              engine = engine)
  }
  originPrior <- list(origin = function(a) stats::dunif(a, 4, 8, log = TRUE))
  rootPrior <- list(root = function(a) stats::dunif(a, 3.5, 6, log = TRUE))
  cases <- list(
    list(orModel, orRanged, c(tree = 1, node = 2, tip = 1, branch = 1, joint = 1,
                              shift = 1, nodeAge = 2, fossilAge = 1, slide = 1,
                              saJump = 1, originAge = 2), originPrior),
    list(crModel, orRanged, c(tree = 1, node = 1, branch = 1, shift = 1, nodeAge = 1,
                              slide = 1, rootAge = 2), rootPrior)
  )
  for (case in cases) {
    cpp <- Run("cpp", case[[1]], case[[2]], case[[3]], case[[4]])
    ref <- Run("R", case[[1]], case[[2]], case[[3]], case[[4]])
    expect_equal(cpp$samples, ref$samples, tolerance = 1e-7)
    expect_equal(cpp$acceptance, ref$acceptance)
    age <- names(case[[3]])[grep("Age$", names(case[[3]]))]
    expect_gt(cpp$acceptance[[age[length(age)]]], 0)
  }
  expect_gt(length(unique(cpp$samples$rootAge)), 1)
})

test_that("Age moves check their set-up and reject what they cannot do", {
  p <- list(root = function(a) 0)
  expect_error(TempoMCMC(orPhylo, orModel, orPar, list(), nIter = 1,
                         moves = c(rootAge = 1), agePriors = p), "nodeAge")
  expect_error(TempoMCMC(orPhylo, crModel, orPar, list(), nIter = 1,
                         moves = c(originAge = 1), agePriors = list(origin = p$root)),
               "origin")
  expect_error(TempoMCMC(orPhylo, crModel, orPar, list(), nIter = 1,
                         moves = c(rootAge = 1)), "agePriors\\$root")
  tree <- .ForModel(TempoTree(orPhylo), crModel)
  set <- .AddSlices(tree, .Precompute(tree, orPar, crModel), 4)
  msg <- .Messages(tree, set)
  set.seed(5)
  d <- .SampleTree(tree, set, msg)
  st <- list(tree = tree, par = orPar, set = set, msg = msg, paths = d$paths,
             bin = d$bin, charLL = 0)
  NoChar <- function(...) 0
  # Zero prior density at the proposal
  res <- .AgeMove("rootAge", st, crModel, 4, NoChar, function(a) if (a == 4) 0 else -Inf, 1)
  expect_false(res$accepted)
  expect_identical(res$state, st)
  # Infeasible frozen model at the new age
  local_mocked_bindings(.Precompute = function(...) .Infeasible("ODE solver failed"))
  expect_false(.AgeMove("rootAge", st, crModel, 4, NoChar, function(a) 0, 1)$accepted)
})

test_that("Root-age moves keep the root above its children and rebuild the model", {
  tree <- .ForModel(TempoTree(orPhylo), crModel)
  set.seed(6)
  run <- TempoMCMC(tree, crModel, orPar, list(), nIter = 40, burnin = 20, nSlice = 4,
                   moves = c(node = 1, branch = 2, rootAge = 2),
                   agePriors = list(root = function(a) stats::dexp(a, 0.25, log = TRUE)))
  expect_true(all(run$samples$rootAge > 3.3))
  expect_gt(stats::sd(run$samples$rootAge), 0)
  st <- run$state
  expect_equal(max(st$set$bounds), st$tree$age[[st$tree$root]])
  # Stem length is fixed while the root moves under an origin: origin age fixed
  run <- TempoMCMC(orRanged, orModel, orPar, list(), nIter = 30, nSlice = 4,
                   moves = c(branch = 1, nodeAge = 4))
  expect_equal(unique(run$samples$originAge), 4.6)
  expect_gt(length(unique(run$samples$rootAge)), 1)
  expect_true(all(run$samples$rootAge < 4.6))
})
