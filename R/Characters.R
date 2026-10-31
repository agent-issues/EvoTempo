#' Mk likelihood of discrete characters
#'
#' Log-likelihood of discrete characters under the Mk (Jukes-Cantor) model
#' with `nState` states and a uniform root, by pruning on given branch lengths
#' in expected changes per character.
#' The C++ engine uses the MkPrime pruning kernels; used by [TempoMCMC()].
#'
#' Sampled ancestors are tips on zero-length branches, so their states are
#' observed at the node that bears them.
#'
#' @param tree A [`TempoTree`][TempoTree()], or a tree that it accepts.
#' @param characters Integer matrix with one row per tip (in the order of the
#' tree's tips) and one column per character, giving states `0` to
#' `nState - 1`, or `NA` where unknown.
#' @param lengths Numeric vector giving the length of the branch above each
#' node (tips, then internal nodes), in expected changes; ignored at the root.
#' @param nState Integer giving the number of states.
#' @inheritParams TempoEngine
#' @returns `MkLikelihood()` returns the natural logarithm of the likelihood.
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' chars <- matrix(c(0, 0, 1, 1, NA, 1), 3)
#' MkLikelihood(tree, chars, c(0.2, 0.2, 0.3, 0, 0.1))
#' @family dCET inference
#' @export
MkLikelihood <- function(tree, characters, lengths, nState = 2L,
                         engine = getOption("EvoTempo.engine", "cpp")) {
  tree <- TempoTree(tree)
  nTip <- tree[["nTip"]]
  characters <- as.matrix(characters)
  if (nrow(characters) != nTip) {
    stop("`characters` needs one row per tip")
  }
  if (TempoEngine(engine) == "cpp") {
    states <- matrix(as.integer(characters), nTip)
    return(.MkLikelihoodCpp(tree, states, as.double(lengths), as.integer(nState)))
  }
  nChar <- ncol(characters)
  partial <- vector("list", nTip + tree[["nNode"]])
  for (i in seq_len(nTip)) {
    p <- matrix(0, nState, nChar)
    obs <- characters[i, ]
    p[, is.na(obs)] <- 1
    known <- which(!is.na(obs))
    p[cbind(obs[known] + 1L, known)] <- 1
    partial[[i]] <- p
  }
  logScale <- 0
  for (v in tree[["postorder"]]) {
    value <- 1
    for (kid in tree[["kids"]][v - nTip, ]) {
      ell <- lengths[[kid]]
      stay <- 1 / nState + (nState - 1) / nState * exp(-nState * ell / (nState - 1))
      move <- (1 - stay) / (nState - 1)
      L <- partial[[kid]]
      value <- value * (rep(move * colSums(L), each = nState) + (stay - move) * L)
    }
    scale <- value[cbind(max.col(t(value), "first"), seq_len(nChar))]
    partial[[v]] <- value / rep(scale, each = nState)
    logScale <- logScale + sum(log(scale))
  }
  sum(log(colSums(partial[[tree[["root"]]]]) / nState)) + logScale
}

# Tempo-weighted branch lengths c * sum_k exp(gamma x_k) dwell_k, per node
# (with an origin, the root's is the stem's, which MkLikelihood() ignores)
.CharacterLengths <- function(tree, x, paths, rate, gamma) {
  if (.UseCpp()) {
    return(.CharacterLengthsCpp(tree, x, paths, rate, gamma))
  }
  age <- tree[["age"]]
  out <- numeric(tree[["nTip"]] + tree[["nNode"]])
  for (v in .PathNodes(tree)) {
    int <- .Intervals(paths[[v]], .TopAge(tree, v), age[[v]])
    out[v] <- rate * sum(exp(gamma * x[int[["bin"]]]) * (int[["hi"]] - int[["lo"]]))
  }
  out
}

# Character log-likelihood of a chain state, for TempoMCMC(): returns
# function(tree, paths, par, set).  With the C++ engine the function holds a
# persistent MkPrime pruner (src/characters.cpp; #26), so that a call
# recomputes only the branch lengths whose path or end ages changed, and the
# nodes below which a branch length or the topology changed, since the
# previous call; it gives the same value, bit for bit, as full pruning.
# Moves report each Metropolis decision with .ResolveCharacters(), which lets
# a rejected proposal be undone from the pruner's journal; correctness does
# not depend on it.  The R engine, and
# option EvoTempo.incremental = FALSE, prune afresh at every call.
.CharacterLikelihood <- function(characters, nState, gamma) {
  if (is.null(characters)) {
    return(function(tr, paths, p, set) 0)
  }
  Lengths <- function(tr, paths, p, set) {
    .CharacterLengths(tr, set[["grid"]][["x"]], paths, p[["rate"]], gamma)
  }
  if (!.UseCpp() || !ncol(characters) ||
      !isTRUE(getOption("EvoTempo.incremental", TRUE))) {
    return(function(tr, paths, p, set) {
      MkLikelihood(tr, characters, Lengths(tr, paths, p, set), nState)
    })
  }
  nTip <- nrow(characters)
  cache <- .CharacterCacheCpp(matrix(as.integer(characters), nTip), as.integer(nState))
  # The cache also memoizes each branch's length (.CharacterLengths())
  structure(function(tr, paths, p, set) {
    if (tr[["nTip"]] != nTip) {
      stop("`characters` needs one row per tip")
    }
    .CachedCharacterLikelihoodCpp(cache, tr, set[["grid"]][["x"]], paths, p[["rate"]],
                                  gamma)
  }, resolve = function(accepted) .ResolveCharacterCacheCpp(cache, accepted),
  cache = cache)
}

# Report the Metropolis decision on the last proposal to a character
# likelihood from .CharacterLikelihood(); a no-op for any other function
.ResolveCharacters <- function(CharLL, accepted) {
  Resolve <- attr(CharLL, "resolve")
  if (!is.null(Resolve)) {
    Resolve(isTRUE(accepted))
  }
  invisible()
}
