# Validation of the likelihood-weighted regraft `wfnpr` (#110;
# dev/spec/weighted-regraft.md) on a five-tip tree with a sampled ancestor, at
# fixed node ages and parameters.
#
# 1. gamma = 0 (characters independent of tempo), with and without characters:
#    the topology has probability proportional to TreeLikelihood() times the Mk
#    likelihood.  Every topology that fixed-height regrafts reach is enumerated,
#    and the frequencies in independent chains are compared with the exact
#    probabilities (z from the between-chain standard error).  `fnpr` chains
#    run alongside as a control.
# 2. gamma = 1 (no exact oracle): `wfnpr` chains against `fnpr` chains.
#
# Run from the package root with the package installed:
#   Rscript data-raw/wfnpr-validation.R [nChain] [nIter]
# Results: data-raw/wfnpr-results/wfnpr-validation.csv
suppressMessages(library("EvoTempo"))
args <- as.numeric(commandArgs(TRUE))
nChain <- if (length(args) > 0) args[[1]] else 8
nIter <- if (length(args) > 1) args[[2]] else 20000
ns <- asNamespace("EvoTempo")
for (f in c(".PathNodes", ".SlideNodes", ".FnprTargets", ".Descendants", ".Rebuild",
            ".AddSlices", ".Precompute", ".Messages", ".SampleTree", ".PathMove",
            ".TreeMove", ".LogSumExp", ".CharacterLengths")) {
  assign(f, get(f, ns))
}

# ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
phylo <- structure(
  list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9), c(9, 4), c(9, 5)),
       edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
       tip.label = c("a", "b", "f", "c", "d"), Nnode = 4), class = "phylo")
start <- TempoTree(phylo)
par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.3)
set.seed(1)
chars <- matrix(sample(0:1, 5 * 4, TRUE), 5)
Key <- function(tr) paste(tr$parent, collapse = " ")

trees <- list(start)
names(trees) <- Key(start)
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

Chain <- function(move, model, CharLL, seed) {
  set.seed(seed)
  set <- .AddSlices(start, .Precompute(start, par, model), 6)
  msg <- .Messages(start, set)
  d <- .SampleTree(start, set, msg)
  st <- list(tree = start, par = par, set = set, msg = msg, paths = d$paths,
             bin = d$bin)
  st$charLL <- CharLL(start, st$paths, par, set)
  visits <- character(nIter)
  accepted <- 0
  for (i in seq_len(nIter)) {
    # gamma = 0: an exact draw of every path; gamma = 1: local path moves
    for (pathMove in if (model$gamma == 0) "tree" else c("node", "tip", "branch")) {
      st <- .PathMove(pathMove, st, CharLL)$state
    }
    res <- .TreeMove(move, st, CharLL)
    st <- res$state
    accepted <- accepted + res$accepted
    visits[[i]] <- Key(st$tree)
  }
  c(table(factor(visits, names(trees))) / nIter, acceptance = accepted / nIter)
}

cases <- list(
  prior = list(gamma = 0, chars = FALSE),
  characters = list(gamma = 0, chars = TRUE),
  tempoCharacters = list(gamma = 1, chars = TRUE)
)
out <- do.call(rbind, lapply(names(cases), function(case) {
  model <- DCETModel(nBin = 7, zRange = c(-2, 2), gamma = cases[[case]]$gamma, rho = 0.7)
  CharLL <- function(tr, paths, p, s) {
    if (!cases[[case]]$chars) {
      return(0)
    }
    if (model$gamma == 0) {
      nodes <- .PathNodes(tr)
      len <- numeric(length(tr$parent))
      len[nodes] <- p$rate * (tr$age[tr$parent[nodes]] - tr$age[nodes])
      return(MkLikelihood(tr, chars, len))
    }
    MkLikelihood(tr, chars, .CharacterLengths(tr, s$grid$x, paths, p$rate, model$gamma))
  }
  runs <- lapply(c(wfnpr = "wfnpr", fnpr = "fnpr"), function(move) {
    vapply(seq_len(nChain), function(i) {
      Chain(move, model, CharLL, 1000 * match(case, names(cases)) + 10 * i + (move == "fnpr"))
    }, numeric(length(trees) + 1L))
  })
  Mean <- function(move) rowMeans(runs[[move]])
  Se <- function(move) apply(runs[[move]], 1, sd) / sqrt(nChain)
  exact <- if (model$gamma == 0) {
    logP <- vapply(trees, function(tr) {
      TreeLikelihood(tr, par, model) + CharLL(tr, NULL, par, NULL)
    }, 0)
    c(exp(logP - .LogSumExp(logP)), NA)
  } else {
    rep(NA_real_, length(trees) + 1L)
  }
  data.frame(
    case = case, topology = c(sprintf("t%02d", seq_along(trees)), "acceptance"),
    exact = exact, wfnpr = Mean("wfnpr"), wfnprSe = Se("wfnpr"), fnpr = Mean("fnpr"),
    fnprSe = Se("fnpr"),
    zExact = (Mean("wfnpr") - exact) / Se("wfnpr"),
    zFnpr = (Mean("wfnpr") - Mean("fnpr")) / sqrt(Se("wfnpr")^2 + Se("fnpr")^2),
    row.names = NULL
  )
}))
out[-(1:2)] <- signif(out[-(1:2)], 4)
write.csv(out, "data-raw/wfnpr-results/wfnpr-validation.csv", row.names = FALSE)
print(out)
