# The C++ engine against the R reference.  Fixtures as test-Paths.R:
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
engineTree <- TempoTree(structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4, root.edge = 0.6),
  class = "phylo"))
enginePar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1,
                  rate = 0.4)
engineModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)

# Run `expr` with the given engine
WithEngine <- function(engine, expr) {
  old <- options(EvoTempo.engine = engine)
  on.exit(options(old))
  expr
}

test_that("TempoEngine() validates its choice", {
  expect_equal(TempoEngine("R"), "R")
  expect_equal(WithEngine("R", TempoEngine()), "R")
  expect_equal(WithEngine(NULL, TempoEngine()), "cpp")
  expect_error(TempoEngine("fortran"))
})

test_that("Extinction probabilities agree between engines", {
  for (grid in list(c(16, 3), c(26, 5))) {
    g <- TempoGrid(grid[1], c(-grid[2], grid[2]), 1.3)
    Q <- TempoGenerator(g, 0.8)
    r <- BinRates(g, 0.5, 0.2, 0.1)
    times <- c(0, 0.1, 0.1, 2, 7.5)
    cpp <- ExtinctionProbability(times, Q, r$lambda, r$mu, r$psi, 0.6, engine = "cpp")
    ref <- ExtinctionProbability(times, Q, r$lambda, r$mu, r$psi, 0.6, engine = "R")
    expect_equal(cpp, ref, tolerance = 1e-9)
    expect_equal(cpp[1, ], rep(0.4, grid[1]))
  }
  expect_error(ExtinctionProbability(c(2, 1), Q, 1, 0.5), "sorted")
})

test_that("Tree likelihoods agree between engines", {
  for (cond in c("extant", "survival", "none")) {
    for (start in c("root", "origin")) {
      for (joint in c(FALSE, TRUE)) {
        model <- engineModel
        model[c("condition", "start", "joint", "removal")] <-
          list(cond, start, joint, if (joint) 0.3 else 0)
        cpp <- TreeLikelihood(engineTree, enginePar, model, engine = "cpp")
        ref <- TreeLikelihood(engineTree, enginePar, model, engine = "R")
        expect_equal(cpp, ref, tolerance = 1e-9)
      }
    }
  }
  # Solver counters (profiling) record the steps taken, then reset
  invisible(.SolverCountsCpp())
  TreeLikelihood(engineTree, enginePar, engineModel, engine = "cpp")
  counts <- .SolverCountsCpp()
  expect_gt(counts[["steps"]], 0)
  expect_gte(counts[["newton"]], counts[["steps"]])
  expect_equal(.SolverCountsCpp()[["steps"]], 0)
  # Stiff generator: high tempo at the edges of a wide grid
  stiff <- DCETModel(nBin = 26, zRange = c(-5, 5), rho = 0.7)
  par <- replace(enginePar, c("theta", "s"), list(1, 1.5))
  expect_equal(TreeLikelihood(engineTree, par, stiff, engine = "cpp"),
               TreeLikelihood(engineTree, par, stiff, engine = "R"),
               tolerance = 1e-9)
  # Product space with one character, including an unknown state
  chr <- c(0, 1, NA, 1, 0)
  expect_equal(
    TreeLikelihood(engineTree, enginePar, engineModel, character = chr,
                   engine = "cpp"),
    TreeLikelihood(engineTree, enginePar, engineModel, character = chr,
                   engine = "R"),
    tolerance = 1e-9)
})

test_that("Pruning with a dense generator agrees between engines", {
  set.seed(1)
  K <- 7
  Q <- matrix(runif(K * K), K)
  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)
  prior <- runif(K)
  prior <- prior / sum(prior)
  lam <- runif(K, 0.3, 0.8)
  Prune <- function(engine) {
    WithEngine(engine, .PruneTree(engineTree, Q, lam, 0.2, 0.1, 0.6, 0.4, prior,
                                  "survival", FALSE, "origin"))
  }
  expect_equal(Prune("cpp"), Prune("R"), tolerance = 1e-9)
  # A zero likelihood (extant tip, rho = 0)
  expect_equal(WithEngine("cpp", .PruneTree(engineTree, Q, lam, 0.2, 0.1, 0, 0,
                                            prior, "none", FALSE, "root")),
               -Inf)
})

test_that("E grids and frozen messages agree between engines", {
  setR <- WithEngine("R", .AddSlices(engineTree, .Precompute(engineTree, enginePar, engineModel), 6))
  setC <- WithEngine("cpp", .AddSlices(engineTree, .Precompute(engineTree, enginePar, engineModel), 6))
  # lsode's error at these tolerances exceeds that of the C++ solver; dE
  # amplifies it by the generator's rates
  expect_equal(setC$eg[c("ages", "E", "G")], setR$eg[c("ages", "E", "G")],
               tolerance = 1e-8)
  expect_equal(setC$eg$dE, setR$eg$dE, tolerance = 1e-6)
  expect_equal(setC$survive, setR$survive, tolerance = 1e-9)
  msgR <- WithEngine("R", .Messages(engineTree, setR))
  msgC <- WithEngine("cpp", .Messages(engineTree, setR))
  expect_equal(msgC, msgR, tolerance = 1e-10)
  expect_equal(WithEngine("cpp", .Pieces(engineTree, setR, 7)),
               WithEngine("R", .Pieces(engineTree, setR, 7)), tolerance = 1e-10)
  expect_equal(WithEngine("cpp", .SliceExp(setR, 3, 0.4)),
               WithEngine("R", .SliceExp(setR, 3, 0.4)), tolerance = 1e-10)
  # Slice eigendecompositions from the same E grid; eigenvectors are unique up
  # to sign, so compare the matrices they rebuild
  slicesC <- WithEngine("cpp", .AddSlices(engineTree, setR, 6))
  expect_equal(slicesC[c("halfPi", "bounds")], setR[c("halfPi", "bounds")])
  Rebuild <- function(sl) sl$vectors %*% (sl$values * t(sl$vectors))
  for (s in seq_along(setR$slices)) {
    expect_equal(slicesC$slices[[s]][c("eBar", "values")],
                 setR$slices[[s]][c("eBar", "values")], tolerance = 1e-12)
    expect_equal(Rebuild(slicesC$slices[[s]]), Rebuild(setR$slices[[s]]),
                 tolerance = 1e-12)
    expect_equal(abs(slicesC$slices[[s]]$vectors), abs(setR$slices[[s]]$vectors),
                 tolerance = 1e-10)
  }
  # A grid at age zero needs no integration
  Q <- setR$Q
  zero <- .EGridCpp(c(0, 0), Q, setR$lambda, setR$mu, setR$psi, 0.3, 1e-10, 1e-14)
  expect_equal(zero$E, matrix(0.7, 2, nrow(Q)))
  expect_equal(zero$G, matrix(0, 2, nrow(Q)))
})

test_that("Parameter moves reuse the frozen model when only `rate` changes", {
  set.seed(4)
  chars <- matrix(sample(0:1, 5 * 3, replace = TRUE), 5)
  priors <- c(TempoPriors()["theta"],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  Run <- function() {
    set.seed(8)
    TempoMCMC(engineTree, engineModel, enginePar, priors, nIter = 8, nSlice = 4,
              characters = chars, moves = c(tree = 1, joint = 1, nodeAge = 1))
  }
  nCall <- 0
  Precompute <- .Precompute
  local_mocked_bindings(.Precompute = function(...) {
    nCall <<- nCall + 1
    Precompute(...)
  })
  reuse <- Run()
  nReuse <- nCall
  # Recomputing every set from scratch leaves the chain unchanged
  nCall <- 0
  local_mocked_bindings(.SameSet = function(...) FALSE)
  fresh <- Run()
  expect_lt(nReuse, nCall)
  expect_identical(reuse$samples, fresh$samples)
  expect_identical(reuse$acceptance, fresh$acceptance)
})

test_that("Path samplers draw the same paths in both engines", {
  set <- .AddSlices(engineTree, .Precompute(engineTree, enginePar, engineModel), 6)
  msg <- .Messages(engineTree, set)
  Draw <- function(engine, seed, f) {
    set.seed(seed)
    WithEngine(engine, f())
  }
  for (seed in 1:5) {
    drawR <- Draw("R", seed, function() .SampleTree(engineTree, set, msg))
    drawC <- Draw("cpp", seed, function() .SampleTree(engineTree, set, msg))
    expect_identical(drawC$bin, drawR$bin)
    expect_equal(drawC$paths, drawR$paths, tolerance = 1e-9)
    # Densities of the drawn configuration
    nodes <- .PathNodes(engineTree)
    expect_equal(WithEngine("cpp", .LogAugmented(engineTree, set, drawR$paths, drawR$bin)),
                 WithEngine("R", .LogAugmented(engineTree, set, drawR$paths, drawR$bin)),
                 tolerance = 1e-10)
    expect_equal(WithEngine("cpp", .LogWeight(engineTree, set, msg, drawR$paths, nodes)),
                 WithEngine("R", .LogWeight(engineTree, set, msg, drawR$paths, nodes)),
                 tolerance = 1e-10)
    expect_equal(
      WithEngine("cpp", .PathWeight(set, msg$pieces[[4]], drawR$paths[[4]])),
      WithEngine("R", .PathWeight(set, msg$pieces[[4]], drawR$paths[[4]])),
      tolerance = 1e-10)
  }
  # Fixed end bin, and a clade below a branch
  for (engine in c("R", "cpp")) {
    br <- Draw(engine, 9, function() .SampleBranch(set, msg$pieces[[5]], NULL, 4L, 8L))
    expect_equal(br$k[c(1, length(br$k))], c(4, 8))
  }
  expect_equal(Draw("cpp", 9, function() .SampleBranch(set, msg$pieces[[5]], NULL, 4L, 8L)),
               Draw("R", 9, function() .SampleBranch(set, msg$pieces[[5]], NULL, 4L, 8L)),
               tolerance = 1e-9)
  base <- Draw("R", 3, function() .SampleTree(engineTree, set, msg))
  cladeR <- Draw("R", 4, function() {
    .SampleClade(engineTree, set, msg, base$paths, base$bin, 7, base$bin[[6]])
  })
  cladeC <- Draw("cpp", 4, function() {
    .SampleClade(engineTree, set, msg, base$paths, base$bin, 7, base$bin[[6]])
  })
  expect_identical(cladeC$bin, as.integer(cladeR$bin))
  expect_equal(cladeC$paths, cladeR$paths, tolerance = 1e-9)
  expect_error(WithEngine("cpp", .SampleClade(engineTree, set, msg, base$paths,
                                              base$bin, 6, 1L)), "No pieces")
  # Bridges
  for (ends in list(c(6, 6), c(5, 8))) {
    bR <- Draw("R", 11, function() .Bridge(set, 2, 0.8, ends[1], ends[2]))
    bC <- Draw("cpp", 11, function() .Bridge(set, 2, 0.8, ends[1], ends[2]))
    expect_equal(bC, bR, tolerance = 1e-9)
  }
})

test_that("Joint conditioning and infeasible proposals are handled in C++", {
  joint <- engineModel
  joint$joint <- TRUE
  set <- .AddSlices(engineTree, .Precompute(engineTree, enginePar, joint), 4)
  msg <- .Messages(engineTree, set)
  set.seed(2)
  dR <- WithEngine("R", .SampleTree(engineTree, set, msg))
  set.seed(2)
  dC <- WithEngine("cpp", .SampleTree(engineTree, set, msg))
  expect_identical(dC$bin, dR$bin)
  expect_equal(WithEngine("cpp", .LogAugmented(engineTree, set, dR$paths, dR$bin)),
               WithEngine("R", .LogAugmented(engineTree, set, dR$paths, dR$bin)),
               tolerance = 1e-10)
  expect_error(WithEngine("cpp", .Bridge(set, 1, 1e-9, 1, 11)),
               class = "infeasibleProposal")
  tiny <- list(edges = c(1e-9, 0), slice = 1L, M = list(.SliceExp(set, 1, 1e-9)))
  expect_error(WithEngine("cpp", .SampleBranch(set, tiny, NULL, 1L, 11L)),
               class = "infeasibleProposal")
  bad <- msg
  bad$D[] <- 0
  expect_error(WithEngine("cpp", .SampleTree(engineTree, set, bad)), "No positive",
               class = "infeasibleProposal")
  noBeta <- msg
  noBeta$beta[[7]] <- lapply(noBeta$beta[[7]], function(b) 0 * b)
  expect_error(WithEngine("cpp", .SampleClade(engineTree, set, noBeta, dR$paths,
                                              dR$bin, 7, dR$bin[[6]])),
               class = "infeasibleProposal")
  # Outside the crown conditioning
  fossilSide <- TempoTree(structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
                                         edge.length = c(1, 2, 2, 1),
                                         tip.label = c("a", "b", "c"), Nnode = 2),
                                    class = "phylo"),
                          tipType = c("extant", "extant", "fossil"))
  fs <- .AddSlices(fossilSide, .Precompute(fossilSide, enginePar, engineModel), 4)
  d <- .SampleTree(fossilSide, fs, .Messages(fossilSide, fs))
  expect_equal(WithEngine("cpp", .LogAugmented(fossilSide, fs, d$paths, d$bin)), -Inf)
})

test_that("Character likelihoods agree between engines", {
  set.seed(3)
  chars <- matrix(sample(0:2, 5 * 6, replace = TRUE), 5)
  chars[2, 3] <- NA
  lengths <- runif(9, 0, 0.6)
  lengths[3] <- 0
  expect_equal(MkLikelihood(engineTree, chars, lengths, 3, engine = "cpp"),
               MkLikelihood(engineTree, chars, lengths, 3, engine = "R"),
               tolerance = 1e-12)
  expect_equal(MkLikelihood(engineTree, chars[, 0, drop = FALSE], lengths, 3,
                            engine = "cpp"), 0)
  expect_error(MkLikelihood(engineTree, chars, lengths, 2, engine = "cpp"),
               "n_states")
  set <- .AddSlices(engineTree, .Precompute(engineTree, enginePar, engineModel), 4)
  d <- .SampleTree(engineTree, set, .Messages(engineTree, set))
  x <- set$grid$x
  expect_equal(WithEngine("cpp", .CharacterLengths(engineTree, x, d$paths, 0.4, 0.7)),
               WithEngine("R", .CharacterLengths(engineTree, x, d$paths, 0.4, 0.7)),
               tolerance = 1e-12)
})

test_that("Chains match between engines", {
  set.seed(4)
  chars <- matrix(sample(0:1, 5 * 3, replace = TRUE), 5)
  priors <- c(TempoPriors()[c("lambda", "theta")],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  model <- engineModel
  Run <- function(engine) {
    set.seed(5)
    TempoMCMC(engineTree, model, enginePar, priors, nIter = 6, burnin = 0,
              nSlice = 4, characters = chars,
              moves = c(tree = 1, node = 1, tip = 1, branch = 1, joint = 1,
                        nodeAge = 1, slide = 1), engine = engine)
  }
  cpp <- Run("cpp")
  ref <- Run("R")
  expect_equal(cpp$samples, ref$samples, tolerance = 1e-7)
  expect_equal(cpp$acceptance, ref$acceptance)
  set.seed(6)
  mC <- MarginalMCMC(engineTree, model, enginePar, priors["lambda"], nIter = 3,
                     engine = "cpp")
  set.seed(6)
  mR <- MarginalMCMC(engineTree, model, enginePar, priors["lambda"], nIter = 3,
                     engine = "R")
  expect_equal(mC, mR, tolerance = 1e-8)
  set.seed(7)
  pC <- SampleTempoPaths(engineTree, enginePar, model, n = 4, nSlice = 4, engine = "cpp")
  set.seed(7)
  pR <- SampleTempoPaths(engineTree, enginePar, model, n = 4, nSlice = 4, engine = "R")
  expect_equal(pC[c("occupancy", "branchX", "meanX")], pR[c("occupancy", "branchX", "meanX")],
               tolerance = 1e-8)
  # Weights near zero: compare on the absolute scale
  expect_lt(max(abs(pC$logWeight - pR$logWeight)), 1e-8)
  expect_equal(pC$logLik, pR$logLik, tolerance = 1e-9)
})
