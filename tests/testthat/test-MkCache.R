# Cached character likelihood (dev/plans/mk-likelihood-cache.md).  Fixture as
# test-Engine.R: a, c fossils with uncertain ages; f a sampled ancestor
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

# Run `expr` with the cache and its exact check on or off
WithCache <- function(on, expr, check = on) {
  old <- options(EvoTempo.mkCache = on)
  oldCheck <- Sys.getenv("EVOTEMPO_MK_CACHE_CHECK", NA)
  Sys.setenv(EVOTEMPO_MK_CACHE_CHECK = if (check) "1" else "0")
  on.exit({
    options(old)
    if (is.na(oldCheck)) {
      Sys.unsetenv("EVOTEMPO_MK_CACHE_CHECK")
    } else {
      Sys.setenv(EVOTEMPO_MK_CACHE_CHECK = oldCheck)
    }
  })
  expr
}

test_that("Cached likelihoods equal one-shot evaluations, recomputing only what changed", {
  set.seed(1)
  chars <- matrix(sample(c(0:2, NA), 5 * 7, replace = TRUE), 5)
  tr <- cacheTree
  n <- tr$nTip + tr$nNode
  .MkCache(clear = TRUE)
  Score <- function(len) {
    before <- .MkCache()$recomputed
    value <- WithCache(TRUE, MkLikelihood(tr, chars, len, 3L))
    expect_identical(value, WithCache(FALSE, MkLikelihood(tr, chars, len, 3L)))
    expect_equal(value, MkLikelihood(tr, chars, len, 3L, engine = "R"), tolerance = 1e-12)
    .MkCache()$recomputed - before
  }
  a <- runif(n, 0, 0.5)
  b <- replace(a, 4, 0.9)  # tip c, below node 9 below the root
  c2 <- replace(a, 2, 0.1)  # tip b, three nodes up
  expect_equal(Score(a), tr$nNode)
  expect_equal(.MkCache()$entries, 1)
  expect_equal(Score(b), 2)  # proposal
  expect_equal(Score(a), 0)  # rejected: restored, not recomputed
  expect_equal(Score(b), 2)  # proposal
  expect_equal(Score(c2), 3)  # nearer to `a`: reject, then update
  expect_equal(Score(c2), 0)  # unchanged
  d <- replace(c2, 5, 0.2)
  expect_equal(Score(d), 2)  # nearer to `c2`: accept, then update
  expect_equal(Score(a), 3)  # nearer to `c2`: reject, then update
  expect_equal(Score(a * 2), tr$nNode)  # `rate` changes every length
  expect_equal(.MkCache()$entries, 1)

  sites <- .MkSiteLogLik(tr, chars, a, 3L)
  oneShot <- vapply(seq_len(ncol(chars)), function(j) {
    MkLikelihood(tr, chars[, j, drop = FALSE], a, 3L, engine = "R")
  }, 0)
  expect_equal(sites, oneShot, tolerance = 1e-12)
  expect_identical(WithCache(FALSE, .MkSiteLogLik(tr, chars, a, 3L)), sites)
  expect_identical(.MkSiteLogLik(tr, chars * 1, a, 3L), sites)
  expect_identical(.MkSiteLogLik(tr, chars[, 0], a, 3L), numeric(0))
  expect_identical(MkLikelihood(tr, chars[, 0], a, 3L), 0)
})

test_that("Cache entries are keyed by data and number of states, and bounded", {
  .MkCache(clear = TRUE)
  tr <- cacheTree
  len <- rep(0.3, tr$nTip + tr$nNode)
  chars <- matrix(c(0L, 1L, 0L, 1L, 1L), 5)
  MkLikelihood(tr, chars, len, 2L)
  MkLikelihood(tr, chars, len, 3L)
  MkLikelihood(tr, 1L - chars, len, 2L)
  expect_equal(.MkCache()$entries, 3)
  for (i in 1:40) {
    MkLikelihood(tr, cbind(chars, i %% 2L, i %/% 2L %% 2L, i %/% 4L %% 2L,
                           i %/% 8L %% 2L, i %/% 16L %% 2L, i %/% 32L %% 2L), len, 2L)
  }
  expect_equal(.MkCache()$entries, 32)
  expect_equal(MkLikelihood(tr, chars, len, 2L),
               MkLikelihood(tr, chars, len, 2L, engine = "R"), tolerance = 1e-12)
  expect_error(MkLikelihood(tr, chars + 1L, len, 2L), "exceeds")
  expect_error(MkLikelihood(tr, chars, len, 1L))
  n <- .MkCache(clear = TRUE)$entries
  # A tree the pruner rejects leaves no entry behind
  bad <- replace(tr, "parent", list(replace(tr$parent, 1, 2L)))
  expect_error(MkLikelihood(bad, chars * 0L, len, 2L), "tip has children")
  expect_equal(.MkCache()$entries, n)
  expect_equal(.MkCache(clear = TRUE)$entries, 0)
})

test_that("Cached chains match uncached chains, exactly, in every (beta, gamma) cell", {
  set.seed(2)
  chars <- matrix(sample(c(0:1, NA), 5 * 4, replace = TRUE, prob = c(4, 4, 1)), 5)
  priors <- c(TempoPriors()[c("lambda", "theta")],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  for (beta in 0:1) for (gamma in 0:1) {
    model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, beta = beta,
                       gamma = gamma)
    moves <- c(node = 1, tip = 1, branch = 1, shift = gamma, tree = 1 - gamma,
               joint = 1 - gamma, nodeAge = 2, fossilAge = 2, slide = 2, saJump = 2)
    Run <- function(on) {
      .MkCache(clear = TRUE)
      set.seed(3)
      WithCache(on, TempoMCMC(cacheTree, model, cachePar, priors, nIter = 8,
                              burnin = 0, nSlice = 4, characters = chars,
                              moves = moves, minAcceptance = 0))
    }
    cached <- Run(TRUE)
    expect_gt(.MkCache()$entries, 0)
    expect_identical(cached$samples, Run(FALSE)$samples)
    expect_gt(sum(cached$acceptance[c("nodeAge", "slide", "saJump")]), 0)
  }
})

test_that("A rebound MkLikelihood() that splits the characters stays exact", {
  set.seed(4)
  chars <- matrix(sample(0:2, 5 * 6, replace = TRUE), 5)
  original <- MkLikelihood
  Split <- function(tree, characters, lengths, nState = 2L, engine = "cpp") {
    original(tree, characters[, 1:3], lengths, 2L, engine) +
      original(tree, characters[, 4:6], lengths, 3L, engine) +
      sum(.MkSiteLogLik(tree, characters[, 1:2] * 0L, lengths, 3L))
  }
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, gamma = 1)
  priors <- c(TempoPriors()["lambda"],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  chars[, 1:3] <- chars[, 1:3] %% 2L
  Run <- function(on) {
    .MkCache(clear = TRUE)
    set.seed(5)
    local_mocked_bindings(MkLikelihood = Split)
    WithCache(on, TempoMCMC(cacheTree, model, cachePar, priors, nIter = 8,
                            burnin = 0, nSlice = 4, characters = chars,
                            moves = c(branch = 1, nodeAge = 1, slide = 2, saJump = 2)))
  }
  cached <- Run(TRUE)
  expect_equal(.MkCache()$entries, 3)
  expect_identical(cached$samples, Run(FALSE)$samples)
})

test_that("Cached chains match uncached chains with a free origin", {
  set.seed(6)
  chars <- matrix(sample(0:1, 5 * 4, replace = TRUE), 5)
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, gamma = 1,
                     start = "origin")
  priors <- c(rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  Run <- function(on) {
    set.seed(7)
    WithCache(on, TempoMCMC(cacheTree, model, cachePar, priors, nIter = 8,
                            burnin = 0, nSlice = 4, characters = chars,
                            moves = c(branch = 1, nodeAge = 2, originAge = 2),
                            agePriors = list(origin = function(a) -a)))
  }
  cached <- Run(TRUE)
  expect_gt(cached$acceptance[["originAge"]], 0)
  expect_identical(cached$samples, Run(FALSE)$samples)
})

test_that("Malformed input is an error, and the self-check fires on a stale result", {
  tr <- cacheTree
  len <- rep(0.3, tr$nTip + tr$nNode)
  chars <- matrix(c(0L, 1L, 0L, 1L, 1L), 5)
  expect_error(.MkSiteLogLik(tr, matrix(c(0, 1, 0, 1), 2), len, 2L), "one row per tip")
  expect_error(.MkCpp(list(nTip = 0L, parent = integer(0)), matrix(0L, 0, 1),
                      numeric(0), 2L, FALSE, TRUE), "two tips")
  expect_error(MkLikelihood(tr, chars, len[-1], 2L), "one value per node")
  zeroRoot <- replace(tr, "parent", list(replace(tr$parent, tr$root, 0L)))
  expect_identical(MkLikelihood(zeroRoot, chars, len, 2L), MkLikelihood(tr, chars, len, 2L))

  # NaN lengths give NaN, not a false alarm
  nanTree <- list(nTip = 3L, parent = c(4L, 4L, 5L, 5L, NA))
  WithCache(TRUE, expect_identical(
    .MkCpp(nanTree, matrix(c(0L, 1L, 0L), 3), c(NaN, 0.2, 0.3, 0.1, 0), 2L, FALSE, TRUE),
    NaN
  ))

  .MkCache(clear = TRUE)
  first <- WithCache(TRUE, MkLikelihood(tr, chars, len, 2L))
  .MkCacheCpp(FALSE, TRUE)
  expect_equal(WithCache(TRUE, MkLikelihood(tr, chars, len, 2L), check = FALSE),
                   first + log(0.5) * ncol(chars))
  .MkCacheCpp(FALSE, TRUE)
  expect_error(WithCache(TRUE, MkLikelihood(tr, chars, len, 2L)), "differs")
  expect_equal(.MkCache()$checkFailures, 1)
  expect_equal(.MkCache()$entries, 0)  # the failed entry is dropped
  expect_identical(WithCache(TRUE, MkLikelihood(tr, chars, len, 2L)), first)
})

test_that("Random call sequences give the one-shot result, bitwise", {
  RandomTree <- function(n) {
    parent <- rep(NA_integer_, 2 * n - 1)
    roots <- seq_len(n)
    for (lab in sample(n + seq_len(n - 1))) {
      pick <- sample(length(roots), 2)
      parent[roots[pick]] <- lab
      roots <- c(roots[-pick], lab)
    }
    parent
  }
  Subtree <- function(parent, v) {
    out <- v
    repeat {
      kids <- setdiff(which(parent %in% out), out)
      if (!length(kids)) return(out)
      out <- c(out, kids)
    }
  }
  # Prune v with its parent p and regraft p above a node t outside the subtree;
  # moving the root's child changes which node is the root
  Spr <- function(parent) {
    cand <- which(!is.na(parent))
    v <- cand[sample(length(cand), 1)]
    p <- parent[v]
    parent[setdiff(which(parent == p), v)] <- parent[p]
    parent[p] <- NA
    ok <- setdiff(seq_along(parent), c(Subtree(parent, v), p))
    t <- ok[sample(length(ok), 1)]
    parent[p] <- parent[t]
    parent[t] <- p
    parent[v] <- p
    parent
  }
  Relabel <- function(st, n) {
    perm <- c(seq_len(n), n + sample(n - 1))
    parent <- st$parent
    parent[perm] <- perm[st$parent]
    len <- st$len
    len[perm] <- st$len
    list(parent = parent, len = len)
  }
  Propose <- function(st, n) {
    N <- length(st$parent)
    switch(sample(c("len1", "all", "spr", "relabel", "zero"), 1, prob = c(4, 1, 4, 1, 1)),
           len1 = st$len[sample(N, 1)] <- stats::rexp(1, 3),
           all = st$len <- st$len * exp(stats::runif(1, -1, 1)),
           spr = st$parent <- Spr(st$parent),
           relabel = st <- Relabel(st, n),
           zero = st$len[sample(N, 1)] <- 0)
    st
  }
  nRoot <- 0
  Eval <- function(st, chars, n, k) {
    tree <- list(nTip = n, parent = as.integer(st$parent))
    expect_identical(.MkCpp(tree, chars, st$len, k, TRUE, TRUE),
                     .MkCpp(tree, chars, st$len, k, TRUE, FALSE))
  }
  set.seed(9)
  .MkCache(clear = TRUE)
  WithCache(TRUE, for (n in c(4, 12)) for (k in 2:3) {
    chars <- matrix(sample(c(0:(k - 1), NA), n * 5, replace = TRUE), n)
    storage.mode(chars) <- "integer"
    # Two chains interleaved on one cache entry
    chains <- lapply(1:2, function(i) {
      list(parent = RandomTree(n), len = stats::rexp(2 * n - 1, 3))
    })
    history <- list()
    for (i in 1:60) {
      ch <- if (i %% 20 < 14) 1 else sample(2, 1)
      cur <- chains[[ch]]
      prop <- Propose(cur, n)
      nRoot <- nRoot + (which(is.na(prop$parent)) != which(is.na(cur$parent)))
      Eval(prop, chars, n, k)
      u <- stats::runif(1)
      if (u < 0.35) {
        chains[[ch]] <- prop
        history[[length(history) + 1]] <- cur
      } else if (u < 0.5) {
        Eval(cur, chars, n, k)
      } else if (u < 0.6 && length(history)) {
        Eval(history[[sample(length(history), 1)]], chars, n, k)
      }
    }
  })
  expect_gt(nRoot, 0)
  expect_equal(.MkCache()$checkFailures, 0)
})

test_that("Cached chains with fnpr match uncached chains and the R engine", {
  set.seed(12)
  chars <- matrix(sample(c(0:2, NA), 5 * 6, replace = TRUE, prob = c(3, 3, 3, 1)), 5)
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7, gamma = 1)
  priors <- c(TempoPriors()["lambda"],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  # Every iteration is stored; with the cache on, each evaluation is also
  # checked bitwise against a fresh one
  Run <- function(on, engine = "cpp") {
    .MkCache(clear = TRUE)
    set.seed(13)
    WithCache(on, TempoMCMC(cacheTree, model, cachePar, priors, nIter = 40,
                            burnin = 0, thin = 1, nSlice = 4, characters = chars,
                            nState = 3L, moves = c(branch = 1, fnpr = 3, nodeAge = 1),
                            minAcceptance = 0, engine = engine))
  }
  cached <- Run(TRUE)
  expect_equal(.MkCache()$checkFailures, 0)
  expect_gt(.MkCache()$entries, 0)
  # fnpr is both accepted (the cache keeps the regrafted tree) and rejected
  # (it returns to the tree before)
  expect_gt(cached$acceptance[["fnpr"]], 0)
  expect_lt(cached$acceptance[["fnpr"]], 1)
  expect_equal(nrow(cached$samples), 40)
  expect_identical(cached$samples, Run(FALSE)$samples)
  expect_equal(cached$samples, Run(FALSE, "R")$samples, tolerance = 1e-8)
})
