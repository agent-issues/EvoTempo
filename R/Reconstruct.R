#' Trees from a simulated history
#'
#' Extracts the sampled (reconstructed) tree from a simulated history, or the
#' complete tree.
#'
#' In the sampled tree, tips are sampled extant lineages and fossils with no
#' sampled descendant; a fossil with sampled descendants is a sampled
#' ancestor, represented as a tip on a zero-length branch.
#' Branch lengths are in time; attribute `characterLength` gives each branch's
#' [BranchClock()], summed over the history's branches that it spans.
#'
#' @inheritParams BranchClock
#' @param complete Logical; if `TRUE`, return the complete tree, with a tip
#' for every extinct and surviving lineage and no fossils.
#' @returns `ReconstructedTree()` returns a tree of class `phylo`, with a
#' `root.edge` (the stem, from the origin), and attributes `tipType`
#' (`"extant"`, `"extinct"`, `"fossil"` or `"ancestor"`), `nodeTime`
#' (time since the origin, for tips then internal nodes) and
#' `characterLength`; or `NULL` if fewer than two tips are sampled.
#' @examples
#' grid <- TempoGrid(41, c(-4, 4))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5, psi = 0.05)
#' set.seed(3)
#' repeat { # Until at least two lineages are sampled
#'   history <- SimulateDCET(60, grid$x, Q, rates$lambda, rates$mu, rates$psi,
#'                           start = 21)
#'   tree <- ReconstructedTree(history)
#'   if (!is.null(tree)) break
#' }
#' table(attr(tree, "tipType"))
#' @family simulation
#' @export
ReconstructedTree <- function(history, gamma = 1, rate = 1, complete = FALSE) {
  edges <- history$edges
  nEdge <- nrow(edges)
  children <- split(seq_len(nEdge)[-1], edges$parent[-1])
  Kids <- function(e) children[[as.character(e)]]
  fossilTimes <- if (complete) {
    rep(list(numeric(0)), nEdge)
  } else {
    lapply(seq_len(nEdge), function(e) sort(history$fossils$time[history$fossils$edge == e]))
  }
  hasSample <- logical(nEdge)
  for (e in rev(seq_len(nEdge))) {
    hasSample[e] <- length(fossilTimes[[e]]) > 0 ||
      (if (complete) edges$type[e] != "removed" else isTRUE(edges$sampled[e])) ||
      any(hasSample[Kids(e)])
  }
  if (!hasSample[1]) {
    return(NULL)
  }
  Tip <- function(kind, at, from, e) {
    list(kind = kind, time = at, len = at - from,
         clock = BranchClock(history, e, from, at, gamma, rate))
  }
  Build <- function(e, from) {
    later <- fossilTimes[[e]][fossilTimes[[e]] > from]
    descends <- any(hasSample[Kids(e)]) ||
      (if (complete) edges$type[e] != "removed" else isTRUE(edges$sampled[e]))
    if (length(later)) {
      f <- later[[1]]
      if (length(later) > 1L || descends) {
        node <- Tip("node", f, from, e)
        node$kids <- list(Tip("ancestor", f, f, e), Build(e, f))
        return(node)
      }
      return(Tip("fossil", f, from, e))
    }
    end <- edges$end[[e]]
    switch(
      edges$type[[e]],
      present = Tip("extant", end, from, e),
      extinction = Tip("extinct", end, from, e),
      speciation = {
        kids <- Kids(e)[hasSample[Kids(e)]]
        if (length(kids) == 2L) {
          node <- Tip("node", end, from, e)
          node$kids <- lapply(kids, Build, end)
          node
        } else {
          below <- Build(kids, end)
          below$len <- below$len + end - from
          below$clock <- below$clock + BranchClock(history, e, from, end, gamma, rate)
          below
        }
      }
    )
  }
  .AsPhylo(Build(1L, 0))
}

.AsPhylo <- function(top) {
  acc <- new.env()
  acc$tips <- list()
  acc$nInternal <- 0L
  Collect <- function(node) {
    if (is.null(node$kids)) {
      acc$tips[[length(acc$tips) + 1L]] <- node
    } else {
      acc$nInternal <- acc$nInternal + 1L
      lapply(node$kids, Collect)
    }
    NULL
  }
  Collect(top)
  nTip <- length(acc$tips)
  if (nTip < 2L) {
    return(NULL)
  }
  nEdge <- nTip + acc$nInternal - 1L
  acc$edge <- matrix(0L, nEdge, 2)
  acc$len <- acc$clock <- numeric(nEdge)
  acc$tipType <- character(nTip)
  acc$nodeTime <- numeric(nTip + acc$nInternal)
  acc$tip <- 0L
  acc$node <- nTip
  acc$row <- 0L
  Number <- function(node) {
    if (is.null(node$kids)) {
      acc$tip <- acc$tip + 1L
      acc$tipType[acc$tip] <- node$kind
      acc$nodeTime[acc$tip] <- node$time
      return(acc$tip)
    }
    acc$node <- acc$node + 1L
    me <- acc$node
    acc$nodeTime[me] <- node$time
    for (kid in node$kids) {
      acc$row <- acc$row + 1L
      row <- acc$row
      acc$edge[row, ] <- c(me, Number(kid))
      acc$len[row] <- kid$len
      acc$clock[row] <- kid$clock
    }
    me
  }
  Number(top)
  structure(
    list(edge = acc$edge, edge.length = acc$len,
         tip.label = paste0(substr(acc$tipType, 1, 1), seq_len(nTip)),
         Nnode = acc$nInternal, root.edge = top$len),
    class = "phylo", order = "cladewise", tipType = acc$tipType,
    nodeTime = acc$nodeTime, characterLength = acc$clock
  )
}

#' Simulate characters on a tree
#'
#' Simulates discrete characters under the Mk (Jukes-Cantor) model with
#' `nState` states, on given branch lengths in expected changes per character.
#'
#' @param tree A tree of class `phylo`.
#' @param nChar Integer giving the number of characters.
#' @param nState Integer giving the number of states.
#' @param lengths Numeric vector giving the length of each edge of `tree`, in
#' expected changes.
#' @returns `SimulateCharacters()` returns an integer matrix, with one row per
#' tip and one column per character, of states `0` to `nState - 1`.
#' @examples
#' grid <- TempoGrid(41, c(-4, 4))
#' Q <- TempoGenerator(grid, theta = 0.01)
#' rates <- BinRates(grid, lambda = 0.51, mu = 0.5)
#' set.seed(2)
#' repeat { # Until the clade survives
#'   tree <- ReconstructedTree(SimulateDCET(40, grid$x, Q, rates$lambda,
#'                                          rates$mu, start = 21))
#'   if (!is.null(tree)) break
#' }
#' SimulateCharacters(tree, 5)
#' @family simulation
#' @importFrom stats runif
#' @export
SimulateCharacters <- function(tree, nChar, nState = 2L,
                               lengths = attr(tree, "characterLength")) {
  nTip <- length(tree$tip.label)
  states <- matrix(0L, nTip + tree$Nnode, nChar)
  root <- nTip + 1L
  states[root, ] <- sample.int(nState, nChar, replace = TRUE) - 1L
  stay <- 1 / nState + (nState - 1) / nState *
    exp(-nState * lengths / (nState - 1))
  for (i in order(.Depth(tree))) {
    from <- states[tree$edge[i, 1], ]
    change <- runif(nChar) > stay[[i]]
    shift <- sample.int(nState - 1L, nChar, replace = TRUE)
    states[tree$edge[i, 2], ] <- ifelse(change, (from + shift) %% nState, from)
  }
  ret <- states[seq_len(nTip), , drop = FALSE]
  rownames(ret) <- tree$tip.label
  ret
}

# Number of edges between each edge's parent and the root.
.Depth <- function(tree) {
  parentEdge <- match(tree$edge[, 1], tree$edge[, 2])
  depth <- integer(nrow(tree$edge))
  for (i in seq_along(depth)) {
    j <- parentEdge[[i]]
    while (!is.na(j)) {
      depth[i] <- depth[i] + 1L
      j <- parentEdge[[j]]
    }
  }
  depth
}
