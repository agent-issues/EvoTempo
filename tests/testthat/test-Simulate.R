grid <- TempoGrid(21, c(-3, 3))
Q <- TempoGenerator(grid, theta = 0.05)
rates <- BinRates(grid, lambda = 0.51, mu = 0.5, psi = 0.05)

Surviving <- function(Simulate, minExtant = 4,
                      Accept = function(h) sum(h$edges$type == "present") >= minExtant) {
  for (i in seq_len(500)) {
    history <- Simulate()
    if (Accept(history)) return(history)
  }
  stop("No accepted history")
}
Removes <- function(h) "removed" %in% h$edges$type && nrow(h$fossils) > 1

test_that("SimulateDCET() validates input", {
  expect_error(SimulateDCET(1, BuddMannGenerator()$x, BuddMannGenerator()$Q,
                            0.51, 0.5, start = 51), "generator")
  full <- matrix(1, 3, 3)
  diag(full) <- -2
  expect_error(SimulateDCET(1, 1:3, full, 1, 0, start = 1), "tridiagonal")
  set.seed(1)
  expect_error(SimulateDCET(200, grid$x, Q, 2, 0, start = 11, maxEdges = 10),
               "maxEdges")
})

test_that("SimulateDCET() histories are consistent", {
  set.seed(2)
  history <- Surviving(function() {
    SimulateDCET(60, grid$x, Q, rates$lambda, rates$mu, rates$psi, start = 11,
                 rho = 0.5)
  })
  expect_s3_class(history, "TempoHistory")
  edges <- history$edges
  expect_true(all(edges$end >= edges$start))
  daughters <- which(!is.na(edges$parent))
  expect_equal(edges$start[daughters], edges$end[edges$parent[daughters]])
  expect_true(all(edges$type[unique(edges$parent[daughters])] == "speciation"))
  expect_true(all(edges$end[edges$type == "present"] == 60))
  expect_true(all(is.na(edges$sampled) == (edges$type != "present")))
  inEdge <- history$fossils$time >= edges$start[history$fossils$edge] &
    history$fossils$time <= edges$end[history$fossils$edge]
  expect_true(all(inEdge))
  expect_true(all(history$fossils$x %in% grid$x))
  # Paths start in the parent's final bin
  lastX <- vapply(history$paths, function(p) p$x[[length(p$x)]], 0)
  firstX <- vapply(history$paths, function(p) p$x[[1]], 0)
  expect_equal(firstX[daughters], lastX[edges$parent[daughters]])
})

test_that("A lineage with no events persists to the present", {
  still <- matrix(0, 3, 3)
  history <- SimulateDCET(10, c(-1, 0, 1), still, 0, 0, start = 2)
  expect_equal(history$edges$type, "present")
  expect_equal(BranchClock(history, 1, gamma = 0), 10)
  expect_equal(BranchClock(history, 1, from = 2, to = 5, rate = 3), 9)
  expect_null(ReconstructedTree(history))
})

test_that("BranchClock() integrates exp(gamma x) along discrete paths", {
  set.seed(3)
  history <- SimulateDCET(20, grid$x, Q * 50, 0, 0.05, start = 11)
  path <- history$paths[[1]]
  dwell <- diff(c(path$t, history$edges$end[1]))
  expect_equal(BranchClock(history, 1), sum(dwell * exp(path$x)))
  expect_equal(BranchClock(history, 1, gamma = 0), sum(dwell))
  expect_gt(length(path$t), 2)
})

test_that("Removal ends sampled lineages", {
  set.seed(4)
  history <- Surviving(function() {
    SimulateDCET(60, grid$x, Q, rates$lambda, rates$mu, rates$psi * 4,
                 start = 11, removal = 1)
  }, Accept = Removes)
  removed <- history$edges$type == "removed"
  expect_true(any(removed))
  expect_true(all(history$edges$end[removed] %in% history$fossils$time))
  tree <- ReconstructedTree(history)
  expect_false("ancestor" %in% attr(tree, "tipType"))
  expect_true("fossil" %in% attr(tree, "tipType"))
})

test_that("SimulateCET() is exact in tempo-time", {
  expect_error(SimulateCET(1, 0.01, 1, 1, 1, delta = 2), "delta")
  set.seed(5)
  flat <- SimulateCET(30, theta = 0.01, s = 1e-9, lambda = 0, mu = 0.1)
  expect_equal(flat$edges$type, "extinction")
  expect_equal(BranchClock(flat, 1), flat$edges$end - flat$edges$start,
               tolerance = 1e-6)
  expect_equal(BranchClock(flat, 1, gamma = 0.5), BranchClock(flat, 1),
               tolerance = 1e-6)
  set.seed(6)
  history <- Surviving(function() {
    SimulateCET(40, 0.01, 1, 0.51, 0.5, psi = 0.05, rho = 0.7)
  })
  edges <- history$edges
  expect_true(all(edges$end[edges$type == "present"] == 40))
  w <- vapply(seq_len(nrow(edges)), function(e) BranchClock(history, e), 0)
  pathW <- vapply(history$paths, function(p) p$w[[length(p$w)]], 0)
  expect_equal(w, pathW)
  daughters <- which(!is.na(edges$parent))
  expect_equal(edges$start[daughters], edges$end[edges$parent[daughters]])
  expect_equal(BranchClock(flat, 1, to = flat$edges$end + 5), BranchClock(flat, 1))
})

test_that("SimulateCET() samples fossils in real or tempo time", {
  set.seed(7)
  removedReal <- Surviving(function() {
    SimulateCET(40, 0.01, 1, 0.51, 0.5, psi = 0.3, removal = 0.5)
  }, Accept = Removes)
  expect_true("removed" %in% removedReal$edges$type)
  set.seed(8)
  tempoTime <- Surviving(function() {
    SimulateCET(40, 0.01, 1, 0.51, 0.5, psi = 0.3, delta = 1, removal = 0.5)
  }, Accept = function(h) Removes(h) && "present" %in% h$edges$type)
  expect_true(all(c("removed", "present") %in% tempoTime$edges$type))
  expect_gt(nrow(tempoTime$fossils), 0)
  set.seed(9)
  expect_error(SimulateCET(400, 0.01, 1, 2, 0, maxEdges = 10), "maxEdges")
})

test_that("ReconstructedTree() returns valid sampled and complete trees", {
  set.seed(10)
  history <- Surviving(function() {
    SimulateDCET(60, grid$x, Q, rates$lambda, rates$mu, rates$psi, start = 11,
                 rho = 0.5)
  })
  tree <- ReconstructedTree(history)
  expect_s3_class(tree, "phylo")
  nTip <- length(tree$tip.label)
  expect_equal(nrow(tree$edge), nTip + tree$Nnode - 1)
  expect_equal(tabulate(tree$edge[, 1])[nTip + seq_len(tree$Nnode)],
               rep(2L, tree$Nnode))
  nodeTime <- attr(tree, "nodeTime")
  expect_equal(nodeTime[tree$edge[, 2]] - nodeTime[tree$edge[, 1]],
               tree$edge.length)
  type <- attr(tree, "tipType")
  expect_equal(sum(type == "extant"), sum(history$edges$sampled, na.rm = TRUE))
  expect_equal(sum(type %in% c("fossil", "ancestor")), nrow(history$fossils))
  expect_true(all(tree$edge.length[match(which(type == "ancestor"),
                                         tree$edge[, 2])] == 0))
  expect_equal(nodeTime[nTip + 1] - tree$root.edge, 0)
  expect_true(all(attr(tree, "characterLength") >= 0))
  complete <- ReconstructedTree(history, complete = TRUE)
  expect_equal(length(complete$tip.label),
               sum(history$edges$type %in% c("present", "extinction")))
  expect_setequal(unique(attr(complete, "tipType")), c("extant", "extinct"))
  continuous <- Surviving(function() {
    SimulateCET(40, 0.01, 1, 0.51, 0.5, psi = 0.05, rho = 0.5)
  })
  ctree <- ReconstructedTree(continuous, gamma = 0, rate = 2)
  expect_equal(attr(ctree, "characterLength"), 2 * ctree$edge.length,
               tolerance = 1e-9)
  unsampled <- SimulateDCET(5, grid$x, Q, rates$lambda, rates$mu, start = 11,
                            rho = 0)
  expect_null(ReconstructedTree(unsampled))
})

test_that("SimulateCharacters() follows branch lengths", {
  set.seed(11)
  tree <- ReconstructedTree(Surviving(function() {
    SimulateDCET(60, grid$x, Q, rates$lambda, rates$mu, start = 11)
  }))
  nTip <- length(tree$tip.label)
  static <- SimulateCharacters(tree, 50, 3, lengths = rep(0, nrow(tree$edge)))
  expect_equal(dim(static), c(nTip, 50))
  expect_equal(rownames(static), tree$tip.label)
  expect_true(all(apply(static, 2, function(x) length(unique(x)) == 1)))
  saturated <- SimulateCharacters(tree, 2000, 4, lengths = rep(50, nrow(tree$edge)))
  expect_equal(as.vector(table(saturated)) / length(saturated), rep(0.25, 4),
               tolerance = 0.05)
  expect_true(all(SimulateCharacters(tree, 3) %in% 0:1))
})
