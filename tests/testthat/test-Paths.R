# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
pathTree <- TempoTree(structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                    c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4), class = "phylo"))
pathPar <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1)
pathModel <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)

# Smaller Monte Carlo samples under coverage instrumentation
nMC <- if (identical(Sys.getenv("R_COVR"), "true")) 150 else 3000

ExpM <- function(A) {
  e <- eigen(A)
  Re(e$vectors %*% diag(exp(e$values)) %*% solve(e$vectors))
}

test_that("E and its integral are interpolated accurately", {
  set <- .Precompute(pathTree, pathPar, pathModel)
  at <- c(0.05, 0.77, 1.9, 3.95)
  exact <- ExtinctionProbability(at, set$Q, set$lambda, set$mu, set$psi, 0.7)
  approx <- t(vapply(at, function(a) .EVector(set$eg, a), numeric(11)))
  expect_equal(approx, exact, tolerance = 1e-8)
  # G is the integral of E
  k <- 4
  G <- .Hermite(set$eg, 2.5, k) - .Hermite(set$eg, 0.3, k)
  E <- function(a) vapply(a, function(x) .Hermite(set$eg, x, k, "E"), 0)
  expect_equal(G, integrate(E, 0.3, 2.5, rel.tol = 1e-10)$value, tolerance = 1e-8)
})

test_that("Survival terms follow the conditioning", {
  extant <- .Precompute(pathTree, pathPar, pathModel)
  noPsi <- pathPar
  noPsi$psi <- 0
  ultra <- TempoTree(structure(list(edge = rbind(c(3, 1), c(3, 2)),
                                    edge.length = c(2, 2), tip.label = c("a", "b"),
                                    Nnode = 1), class = "phylo"))
  a <- .Precompute(ultra, noPsi, pathModel)
  expect_equal(a$survive, (1 - ExtinctionProbability(2, a$Q, a$lambda, a$mu, 0, 0.7)[1, ])^2,
               tolerance = 1e-8)
  surv <- pathModel
  surv$condition <- "survival"
  b <- .Precompute(pathTree, pathPar, surv)
  expect_true(all(b$survive > extant$survive))
  surv$condition <- "none"
  expect_equal(.Precompute(pathTree, pathPar, surv)$survive, rep(1, 11))
})

test_that("Bridges have the endpoint-conditioned Laplace transform (rung 4c)", {
  set <- .AddSlices(pathTree, .Precompute(pathTree, pathPar, pathModel), 4)
  s <- 2
  span <- 0.8
  sl <- set$slices[[s]]
  G <- .AddDiag(set$Q, 2 * set$lambda * sl$eBar - set$sigma)
  clock <- exp(set$grid$x)
  set.seed(1)
  for (ends in list(c(6, 6), c(5, 8))) {
    bridges <- replicate(nMC, .Bridge(set, s, span, ends[1], ends[2]),
                         simplify = FALSE)
    expect_true(all(vapply(bridges, function(br) {
      identical(br$k[c(1, length(br$k))], as.numeric(ends)) ||
        identical(br$k[c(1, length(br$k))], as.integer(ends))
    }, TRUE)))
    expect_true(all(vapply(bridges, function(br) all(abs(diff(br$k)) == 1), TRUE)))
    ell <- vapply(bridges, function(br) sum(clock[br$k] * diff(c(0, br$t, span))), 0)
    for (u in c(0.5, 2)) {
      expected <- ExpM((G - u * diag(clock)) * span)[ends[1], ends[2]] /
        ExpM(G * span)[ends[1], ends[2]]
      expect_equal(mean(exp(-u * ell)), expected, tolerance = 4 * sd(exp(-u * ell)) / sqrt(nMC) / expected)
    }
  }
  # Frozen transfer matches the matrix exponential
  expect_equal(.SliceExp(set, s, span), ExpM(G * span), tolerance = 1e-10)
})

test_that("Exact / frozen importance weights recover the marginal likelihood", {
  set <- .AddSlices(pathTree, .Precompute(pathTree, pathPar, pathModel), 8)
  msg <- .Messages(pathTree, set)
  exact <- TreeLikelihood(pathTree, pathPar, pathModel)
  frozen <- .FrozenLikelihood(pathTree, set, msg)
  expect_lt(abs(exact - frozen), 0.01)
  set.seed(2)
  maps <- SampleTempoPaths(pathTree, pathPar, pathModel, n = 60, nSlice = 8)
  expect_equal(maps$logLik, exact, tolerance = 5 * sd(exp(maps$logWeight)) / sqrt(60))
  expect_equal(dim(maps$occupancy), c(60, 11))
  expect_equal(rowSums(maps$occupancy), rep(sum(pathTree$age[pathTree$parent[.PathNodes(pathTree)]] -
                                                  pathTree$age[.PathNodes(pathTree)]), 60))
  # A sampled configuration has the stated bins at its path ends
  draw <- .SampleTree(pathTree, set, msg)
  for (v in .PathNodes(pathTree)) {
    path <- draw$paths[[v]]
    expect_equal(path$k[[1]], draw$bin[[pathTree$parent[[v]]]])
    expect_equal(path$k[[length(path$k)]], draw$bin[[v]])
    expect_true(all(diff(c(pathTree$age[[pathTree$parent[[v]]]], path$t,
                           pathTree$age[[v]])) <= 0))
  }
  expect_equal(draw$bin[3], draw$bin[8]) # Sampled ancestor shares its node's bin
  expect_true(is.finite(.LogAugmented(pathTree, set, draw$paths, draw$bin)))
  # Joint conditioning
  joint <- pathModel
  joint$joint <- TRUE
  setJ <- .AddSlices(pathTree, .Precompute(pathTree, pathPar, joint), 8)
  msgJ <- .Messages(pathTree, setJ)
  expect_lt(abs(TreeLikelihood(pathTree, pathPar, joint) -
                  .FrozenLikelihood(pathTree, setJ, msgJ, TRUE)), 0.01)
  drawJ <- .SampleTree(pathTree, setJ, msgJ)
  expect_true(is.finite(.LogAugmented(pathTree, setJ, drawJ$paths, drawJ$bin)))
})

test_that("Path moves leave the exact conditional law invariant", {
  set.seed(3)
  run <- TempoMCMC(pathTree, pathModel, pathPar, list(), nIter = nMC / 75, nSlice = 6,
                   moves = c(tree = 0, node = 2, tip = 2, branch = 2))
  expect_true(all(run$acceptance[c("node", "tip", "branch")] > 0.9))
  st <- run$state
  for (v in .PathNodes(pathTree)) {
    path <- st$paths[[v]]
    expect_equal(path$k[[1]], st$bin[[pathTree$parent[[v]]]])
    expect_equal(path$k[[length(path$k)]], st$bin[[v]])
  }
  expect_equal(st$bin[3], st$bin[8])
  # Mean log-tempo agrees with importance-weighted FFBS draws
  maps <- SampleTempoPaths(pathTree, pathPar, pathModel, n = nMC / 30, nSlice = 6)
  w <- exp(maps$logWeight - max(maps$logWeight))
  expect_lt(abs(mean(run$samples$meanX) - sum(w * maps$meanX) / sum(w)), 1)
})

test_that("TempoMCMC() updates parameters and records summaries", {
  set.seed(4)
  priors <- TempoPriors()[c("lambda", "theta")]
  priors$s <- function(x) stats::dunif(x, 0.5, 1.5, log = TRUE)
  run <- TempoMCMC(pathTree, pathModel, pathPar, priors, nIter = 6, burnin = 20,
                   thin = 2, nSlice = 4, moves = c(tree = 1))
  expect_equal(nrow(run$samples), 3)
  expect_named(run$samples[1:6], c("lambda", "theta", "s", "logPost", "rootBin", "meanX"))
  expect_equal(ncol(run$branchX), length(.PathNodes(pathTree)))
  expect_gt(run$acceptance[["tree"]], 0.9)
  fixed <- TempoMCMC(pathTree, pathModel, pathPar, priors, nIter = 2, nSlice = 4,
                     updateParameters = FALSE)
  expect_equal(fixed$samples$lambda, rep(0.5, 2))
  origin <- pathModel
  origin$start <- "origin"
  expect_error(TempoMCMC(pathTree, origin, pathPar, priors, nIter = 1), "root")
})

test_that("MarginalMCMC() samples parameters", {
  set.seed(5)
  priors <- TempoPriors()["lambda"]
  out <- MarginalMCMC(pathTree, pathModel, pathPar, priors, nIter = 4,
                      burnin = 20, thin = 2)
  expect_equal(dim(out), c(2, 2))
  expect_true(all(is.finite(out$logLik)))
  expect_error(TempoPriors(c(a = 1), c(b = 1)), "same parameters")
})

test_that("Path moves reject when the character likelihood forbids them", {
  set <- .AddSlices(pathTree, .Precompute(pathTree, pathPar, pathModel), 4)
  msg <- .Messages(pathTree, set)
  set.seed(6)
  d <- .SampleTree(pathTree, set, msg)
  st <- list(tree = pathTree, par = pathPar, set = set, msg = msg,
             paths = d$paths, bin = d$bin, charLL = 0)
  res <- .PathMove("branch", st, function(...) -Inf)
  expect_false(res$accepted)
  expect_identical(res$state$paths, st$paths)
})

test_that("Infeasible bridges are signalled and their moves rejected", {
  set <- .AddSlices(pathTree, .Precompute(pathTree, pathPar, pathModel), 4)
  expect_error(.Bridge(set, 1, 1e-9, 1, 11), class = "infeasibleProposal")
  # A tree 1e-7 Myr deep whose node bins alternate between the extreme bins
  tiny <- TempoTree(structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
                                   edge.length = c(1, 1, 1, 2) * 1e-7,
                                   tip.label = c("a", "b", "c"), Nnode = 2),
                              class = "phylo"))
  set <- .AddSlices(tiny, .Precompute(tiny, pathPar, pathModel), 2)
  bin <- c(11L, 11L, 1L, 11L, 1L)
  st <- list(tree = tiny, par = pathPar, set = set, msg = .Messages(tiny, set),
             paths = list(), bin = bin, charLL = 0)
  res <- .PathMove("branch", st, function(...) 0)
  expect_false(res$accepted)
})

test_that(".Sample1() draws in proportion and rejects empty weights", {
  set.seed(9)
  draws <- vapply(1:4000, function(i) .Sample1(c(0.2, 0, 0.5, 0.3)), 0L)
  expect_equal(as.vector(table(factor(draws, 1:4))) / 4000, c(0.2, 0, 0.5, 0.3),
               tolerance = 0.1)
  expect_equal(.Sample1(5), 1L)
  expect_error(.Sample1(c(0, -1e-18)), "No positive")
  expect_error(.Sample1(c(NaN, 1)), "No positive")
})
