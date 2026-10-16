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
.CharacterLengths <- function(tree, x, paths, rate, gamma) {
  if (.UseCpp()) {
    return(.CharacterLengthsCpp(tree, x, paths, rate, gamma))
  }
  age <- tree[["age"]]
  out <- numeric(tree[["nTip"]] + tree[["nNode"]])
  for (v in .PathNodes(tree)) {
    int <- .Intervals(paths[[v]], age[[tree[["parent"]][[v]]]], age[[v]])
    out[v] <- rate * sum(exp(gamma * x[int[["bin"]]]) * (int[["hi"]] - int[["lo"]]))
  }
  out
}
