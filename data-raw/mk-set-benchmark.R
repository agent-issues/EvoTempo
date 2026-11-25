# Fixed cost of one character-likelihood evaluation on the penguin data (b1g1),
# and its share of a production-schedule iteration (#110).
#
#   Rscript data-raw/mk-set-benchmark.R <mode> <n> [variant] [common.R]
#
# Run from the package root with EvoTempo installed (and on .libPaths()).
# mode:
#   micro   ms per PartitionedMk() call over `n` calls, with no edge, one tip
#           edge or every edge changed since the last call (each call is then
#           an accepted proposal); beside it the same evaluation by twelve
#           direct .MkCpp() calls behind a thin R wrapper (`thinTwelve`), by the
#           bare C++ calls (`direct12`, `directSet`), and the cost of making the
#           lengths (`lengthsOnly`)
#   iter    ms per iteration of the production schedule over `n` iterations,
#           and a checksum of the samples, which must agree between variants
#   phases  ms per evaluation inside MkPrime::Pruner, by phase, over `n`
#           proposals (data-raw/mk-set-phases.cpp; needs a compiler)
# variant "twelve" rebinds the likelihood to `thinTwelve`: the thinner R path
# alone, without the single C++ entry point.  `thinTwelve` is thinner than the
# real fast path, which keeps TempoTree(), the attribute lookups, MkSetFor() and
# the power, so the difference between the two understates what the single
# entry point saves.  `common.R` is the pipeline file to
# source: give that of an earlier commit, with a build of it, for a baseline.
a <- commandArgs(trailingOnly = TRUE)
mode <- a[[1]]
n <- as.integer(a[[2]])
variant <- if (is.na(a[3])) "default" else a[[3]]
suppressMessages(library("EvoTempo"))
source(if (is.na(a[4])) "data-raw/empirical/penguins/cet/common.R" else a[[4]])

dat <- LoadPenguins()
tree <- dat$tree
chars <- dat$characters
parts <- attr(chars, "parts")
pat <- attr(chars, "mkv")
kU <- unique(pat$k)
ofK <- lapply(kU, function(k) which(pat$k == k))
consts <- unname(pat$const[as.character(kU)])
Mk <- get(".MkCpp", .ns)
Twelve <- function(tree, characters, lengths, ...) {
  ll <- 0
  for (p in parts) {
    ll <- ll + Mk(tree, p$x, lengths, p$k, FALSE, TRUE)
  }
  site <- numeric(length(pat$k))
  for (j in seq_along(kU)) {
    site[ofK[[j]]] <- Mk(tree, consts[[j]], lengths, kU[[j]], TRUE, TRUE)
  }
  ll + -sum(pat$count * log1p(-pat$k * exp(site)))
}
if (variant == "twelve") {
  PartitionedMk <- Twelve
  utils::assignInNamespace("MkLikelihood", Twelve, ns = "EvoTempo")
}
set.seed(1)
len <- runif(tree$nTip + tree$nNode, 0.01, 0.3)

if (mode == "iter") {
  model <- CetModel(1, 1)
  par <- list(lambda = 0.2, mu = 0.15, psi = 0.05, theta = 0.05, s = 0.2, rate = 0.017)
  Run <- function(nIter) {
    set.seed(1)
    .cet$records <- list()
    t0 <- proc.time()[[3]]
    run <- TempoMCMC(tree, model, par, CetPriors(), nIter = nIter, burnin = 0, thin = 1,
                     nSlice = 20, moves = CetMoves(tree, chars, model), characters = chars,
                     agePriors = CetAgePriors(), minAcceptance = 0)
    list(secs = proc.time()[[3]] - t0,
         checksum = format(sum(vapply(run$samples, sum, 0)), digits = 17))
  }
  Run(2)
  res <- Run(n)
  cat(variant, round(1000 * res$secs / n, 1), "ms per iteration; checksum",
      res$checksum, "\n")
} else if (mode == "micro") {
  Time <- function(Call) {
    Call(1L)
    t0 <- proc.time()[[3]]
    for (i in seq_len(n)) Call(i)
    1000 * (proc.time()[[3]] - t0) / n
  }
  # `oneTipEdge` changes one edge of the lengths of the previous call
  running <- len
  Lengths <- list(
    unchanged = function(i) len,
    oneTipEdge = function(i) {
      running[1 + i %% 50] <<- len[1 + i %% 50] * (1 + 1e-3 * (1 + i %% 7))
      running
    },
    allEdges = function(i) len * (1 + 1e-6 * i)
  )
  Direct12 <- function(l) {
    for (p in parts) Mk(tree, p$x, l, p$k, FALSE, TRUE)
    for (j in seq_along(kU)) Mk(tree, consts[[j]], l, kU[[j]], TRUE, TRUE)
  }
  res <- t(vapply(Lengths, function(L) {
    row <- c(lengthsOnly = Time(function(i) L(i)),
             PartitionedMk = Time(function(i) PartitionedMk(tree, chars, L(i))),
             thinTwelve = Time(function(i) Twelve(tree, chars, L(i))),
             direct12 = Time(function(i) Direct12(L(i))), directSet = NA)
    if (exists(".MkSetLogLikCpp", .ns)) {
      s <- MkSetFor(parts, pat, tree$nTip)$set
      SetCpp <- get(".MkSetLogLikCpp", .ns)
      row[["directSet"]] <- Time(function(i) SetCpp(s$ptr, tree$parent, L(i), s$matrices, s$nState))
    }
    row
  }, numeric(5)))
  print(round(res, 4))
} else if (mode == "phases") {
  Rcpp::sourceCpp("data-raw/mk-set-phases.cpp")
  mats <- c(lapply(parts, `[[`, "x"), consts)
  k <- as.integer(c(vapply(parts, `[[`, 0, "k"), kU))
  res <- rbind(
    oneTipEdge = Phases(tree$parent, mats, k, len, n, 0, TRUE),
    oneTipEdgeAccepted = Phases(tree$parent, mats, k, len, n, 0, FALSE),
    allEdges = Phases(tree$parent, mats, k, len, n, 1, TRUE),
    oneTipEdgeBySetEdgeLength = Phases(tree$parent, mats, k, len, n, 2, TRUE)
  )
  print(round(cbind(res, total = rowSums(res[, 1:4])), 4))
}
