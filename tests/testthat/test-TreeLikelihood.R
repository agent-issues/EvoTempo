Phylo <- function(edge, len, labels = paste0("t", seq_len(nrow(edge) / 2 + 1)),
                  nNode = length(labels) - 1L) {
  structure(list(edge = edge, edge.length = len, tip.label = labels,
                 Nnode = nNode), class = "phylo")
}
# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
saTree <- Phylo(rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9),
                      c(9, 4), c(9, 5)),
                c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3), c("a", "b", "f", "c", "d"))

test_that("TempoTree() reads ages and tip types", {
  tt <- TempoTree(saTree)
  expect_s3_class(tt, "TempoTree")
  expect_identical(TempoTree(tt), tt)
  expect_equal(tt$tipType, c("fossil", "extant", "ancestor", "fossil", "extant"))
  expect_equal(tt$age, c(1, 0, 1.5, 2.1, 0, 4, 3, 1.5, 3.3))
  expect_equal(tt$sampled, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(tt$root, 6L)
  expect_true(all(match(tt$parent[tt$postorder][-length(tt$postorder)],
                        tt$postorder) > seq_len(length(tt$postorder) - 1)))
  expect_equal(tt$stem, 0)
  rooted <- saTree
  rooted$root.edge <- 0.5
  expect_equal(TempoTree(rooted)$stem, 0.5)

  # Ages and types from attributes, as set by ReconstructedTree()
  attributed <- structure(saTree, nodeTime = 4 - tt$age,
                          tipType = c("extinct", "extant", "ancestor", "fossil",
                                      "extant"))
  expect_equal(TempoTree(attributed, present = 4)$age, tt$age)
  expect_equal(TempoTree(attributed)$tipType[1], "fossil")
  expect_equal(TempoTree(saTree, ages = tt$age + 1e-10)$age, tt$age, tolerance = 1e-9)
})

test_that("TempoTree() validates input", {
  expect_error(TempoTree(list()), "phylo")
  expect_error(TempoTree(Phylo(rbind(c(4, 1), c(4, 2), c(4, 3)), c(1, 1, 1),
                               c("a", "b", "c"), nNode = 1)), "binary")
  expect_error(TempoTree(saTree, ages = rev(TempoTree(saTree)$age)), "decrease")
  expect_error(TempoTree(saTree, tipType = rep("tip", 5)), "tipType")
  expect_error(TempoTree(saTree, tipType = rep("extant", 5)), "age zero")
  twoSA <- Phylo(rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)), c(1, 0, 0, 2))
  expect_error(TempoTree(twoSA, tipType = c("ancestor", "ancestor", "extant")),
               "two sampled")
})

test_that("DCETModel() validates and stores settings", {
  m <- DCETModel(nBin = 11, rootPrior = "unif", condition = "none")
  expect_s3_class(m, "DCETModel")
  expect_equal(m$rootPrior, "uniform")
  expect_equal(DCETModel(nBin = 3, rootPrior = c(1, 2, 1))$rootPrior, c(1, 2, 1))
  expect_error(DCETModel(nBin = 3, rootPrior = c(1, 2)), "one non-negative")
  expect_error(DCETModel(rootPrior = "flat"))
})

test_that("Start-bin priors", {
  par <- list(lambda = 0.6, mu = 0.2, theta = 0.3, s = 1)
  for (prior in list("equilibrium", "stationary", "uniform", 1:11)) {
    set <- .Setup(par, DCETModel(nBin = 11, zRange = c(-2, 2), rootPrior = prior))
    expect_equal(sum(set$prior), 1)
  }
  expect_equal(.Setup(par, DCETModel(nBin = 11, zRange = c(-2, 2),
                                     rootPrior = 1:11))$prior, (1:11) / 66)
  eq <- .Setup(par, DCETModel(nBin = 11, zRange = c(-2, 2)))
  expect_equal(eq$prior, TempoEquilibrium(eq$Q, eq$lambda, eq$mu)$p)
})

test_that("beta = 0 collapses to the analytic FBD likelihood for any Q", {
  rooted <- saTree
  rooted$root.edge <- 0.8
  tt <- TempoTree(rooted)
  set.seed(1)
  for (condition in c("none", "extant", "survival")) {
    for (removal in c(0, 0.4)) {
      for (start in c("root", "origin")) {
        nBin <- 6
        Q <- matrix(runif(nBin * nBin), nBin)
        diag(Q) <- 0
        diag(Q) <- -rowSums(Q)
        prior <- runif(nBin)
        expect_equal(
          .PruneTree(tt, Q, 0.5, 0.2, 0.1, 0.6, removal, prior / sum(prior),
                     condition, FALSE, start),
          FBDTreeLikelihood(tt, 0.5, 0.2, 0.1, 0.6, removal, condition, start),
          tolerance = 1e-8
        )
      }
    }
  }
  # Through the public interface, on a tempo generator
  par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.4, s = 1.3)
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), beta = 0, rho = 0.6)
  expect_equal(TreeLikelihood(saTree, par, model),
               FBDTreeLikelihood(saTree, 0.5, 0.2, 0.1, 0.6), tolerance = 1e-8)
  # Extant-only, no fossil sampling
  extant <- Phylo(rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)), c(1, 2, 2, 3))
  par$psi <- NULL
  expect_equal(TreeLikelihood(extant, par, model),
               FBDTreeLikelihood(extant, 0.5, 0.2, 0, 0.6), tolerance = 1e-8)
})

test_that("Joint and marginal conditioning coincide when tempo is irrelevant", {
  par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.4, s = 1)
  joint <- DCETModel(nBin = 11, zRange = c(-2, 2), beta = 0, joint = TRUE)
  marginal <- DCETModel(nBin = 11, zRange = c(-2, 2), beta = 0)
  expect_equal(TreeLikelihood(saTree, par, joint),
               TreeLikelihood(saTree, par, marginal), tolerance = 1e-8)
  # but not otherwise
  joint$beta <- marginal$beta <- 1
  expect_gt(abs(TreeLikelihood(saTree, par, joint) -
                  TreeLikelihood(saTree, par, marginal)), 1e-3)
})

test_that("Tempo-dependent likelihood is finite and responds to beta", {
  par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.4, s = 1)
  rooted <- saTree
  rooted$root.edge <- 0.8
  m1 <- DCETModel(nBin = 11, zRange = c(-2, 2), start = "origin")
  ll1 <- TreeLikelihood(rooted, par, m1)
  expect_true(is.finite(ll1))
  m0 <- m1
  m0$beta <- 0
  expect_gt(abs(ll1 - TreeLikelihood(rooted, par, m0)), 1e-3)
  # Survival conditioning with removal
  m1$removal <- 0.3
  m1$condition <- "survival"
  expect_true(is.finite(TreeLikelihood(rooted, par, m1)))
})

test_that("TreeLikelihood() rejects impossible settings", {
  par <- list(lambda = 0.5, mu = 0.2, theta = 0.4, s = 1)
  model <- DCETModel(nBin = 11, zRange = c(-2, 2))
  expect_error(TreeLikelihood(saTree, par, model), "psi > 0")
  par$psi <- 0.1
  expect_error(TreeLikelihood(saTree, par, DCETModel(nBin = 11, zRange = c(-2, 2),
                                                    start = "origin")),
               "root edge")
  rootSA <- Phylo(rbind(c(4, 1), c(4, 5), c(5, 2), c(5, 3)), c(0, 1, 2, 2))
  expect_error(TreeLikelihood(rootSA, par, model), "root cannot")
  rootSA$root.edge <- 1
  model$start <- "origin"
  expect_true(is.finite(TreeLikelihood(rootSA, par, model)))
  model$rho <- 0
  model$condition <- "none"
  expect_equal(TreeLikelihood(rootSA, par, model), -Inf)
})

test_that("FBDTreeLikelihood() matches quadrature of the D equation", {
  # Two extant tips: D = rho * exp(int g); root carries no lambda
  lambda <- 0.7
  mu <- 0.3
  rho <- 0.5
  g <- function(t) -(lambda + mu) + 2 * lambda * FBDExtinction(t, lambda, mu, 0, rho)
  cherry <- Phylo(rbind(c(3, 1), c(3, 2)), c(2, 2))
  expect_equal(FBDTreeLikelihood(cherry, lambda, mu, rho = rho, condition = "none"),
               2 * (log(rho) + integrate(g, 0, 2, rel.tol = 1e-12)$value),
               tolerance = 1e-9)
})

test_that("Trees incompatible with crown conditioning have zero likelihood", {
  # Root daughters: (a, b) with an extant tip, and a lone fossil c
  fossilSide <- Phylo(rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)), c(1, 2, 2, 1),
                      c("a", "b", "c"))
  tt <- TempoTree(fossilSide, tipType = c("extant", "extant", "fossil"))
  par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.4, s = 1)
  model <- DCETModel(nBin = 11, zRange = c(-2, 2))
  expect_equal(TreeLikelihood(tt, par, model), -Inf)
  expect_equal(FBDTreeLikelihood(tt, 0.5, 0.2, 0.1), -Inf)
  model$condition <- "survival"
  expect_true(is.finite(TreeLikelihood(tt, par, model)))
  expect_true(.Compatible(tt, "survival"))
  # With an origin, at least one extant tip is needed
  fossilsOnly <- TempoTree(fossilSide, ages = c(1, 1, 2, 3, 2),
                           tipType = rep("fossil", 3))
  expect_false(.Compatible(fossilsOnly, "extant", origin = TRUE))
  expect_true(.Compatible(tt, "extant", origin = TRUE))
})
