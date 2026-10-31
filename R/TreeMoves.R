# Time-tree moves for TempoMCMC().  Each move changes the tree, redraws the
# bins of moved nodes from their frozen conditional law and the paths on
# affected branches from frozen bridges, and is accepted on the exact target
# with explicit proposal densities in both directions.  These moves keep the
# start age (the root's, or with an origin the origin's) fixed; .AgeMove() in
# R/MCMC.R moves it.

# Recompute derived fields after a change to `parent`, `age` or `tipType`
.Rebuild <- function(tree) {
  if (.UseCpp()) {
    return(.RebuildCpp(tree))
  }
  nTip <- tree[["nTip"]]
  nAll <- nTip + tree[["nNode"]]
  parent <- tree[["parent"]]
  kids <- matrix(0L, tree[["nNode"]], 2)
  for (i in seq_len(tree[["nNode"]])) {
    kids[i, ] <- which(parent == nTip + i)
  }
  tree[["kids"]] <- kids
  isAncestor <- c(tree[["tipType"]] == "ancestor", logical(tree[["nNode"]]))
  tree[["sampled"]] <- isAncestor[kids[, 1]] | isAncestor[kids[, 2]]
  sa <- which(isAncestor)
  tree[["age"]][sa] <- tree[["age"]][parent[sa]]
  order <- .Preorder(parent, tree[["root"]])
  tree[["postorder"]] <- rev(order[order > nTip])
  tree
}

# Log of the frozen transfer row from bin k (vector over end bins)
.LogTransferRow <- function(pc, k) {
  if (.UseCpp()) {
    return(.LogTransferRowCpp(pc, k))
  }
  r <- replace(numeric(nrow(pc[["M"]][[1]])), k, 1)
  logScale <- 0
  for (M in pc[["M"]]) {
    r <- as.vector(r %*% M)
    top <- max(r)
    r <- r / top
    logScale <- logScale + log(top)
  }
  log(pmax(r, 0)) + logScale
}

.LogTransferCol <- function(pc, j) {
  if (.UseCpp()) {
    return(.LogTransferColCpp(pc, j))
  }
  r <- replace(numeric(nrow(pc[["M"]][[1]])), j, 1)
  logScale <- 0
  for (M in rev(pc[["M"]])) {
    r <- as.vector(M %*% r)
    top <- max(r)
    r <- r / top
    logScale <- logScale + log(top)
  }
  log(pmax(r, 0)) + logScale
}

# Exact log density of one branch's path (CTMC part and exp(int g))
.LogPathExact <- function(tree, set, v, path) {
  if (.UseCpp()) {
    return(.LogPathExactCpp(tree, set, v, path))
  }
  int <- .Intervals(path, .TopAge(tree, v), tree[["age"]][[v]])
  k <- int[["bin"]]
  Q <- set[["Q"]]
  kk <- path[["k"]]
  sum((diag(Q)[k] - set[["sigma"]][k]) * (int[["hi"]] - int[["lo"]]) +
        2 * set[["lambda"]][k] * (.Hermite(set[["eg"]], int[["hi"]], k) -
                                    .Hermite(set[["eg"]], int[["lo"]], k))) +
    (if (length(kk) > 1L) sum(log(Q[cbind(kk[-length(kk)], kk[-1])])) else 0)
}

# Log density of a path under the frozen bridge between its end bins
.LogBridge <- function(tree, set, v, path, pc = .Pieces(tree, set, v)) {
  k <- path[["k"]]
  .LogPathExact(tree, set, v, path) - .PathWeight(set, pc, path) -
    .LogTransferRow(pc, k[[1]])[[k[[length(k)]]]]
}

# Draw a bridge on branch v of `tree`; returns the path and its log density
.DrawBridge <- function(tree, set, v, kTop, kEnd) {
  pc <- .Pieces(tree, set, v)
  path <- .SampleBranch(set, pc, NULL, kTop, kEnd)
  list(path = path, logQ = .LogBridge(tree, set, v, path, pc))
}

# Redraw the paths on branches `nodes` between the fixed bins `bin`; returns
# the updated paths and their log proposal density
.DrawBridges <- function(tree, set, nodes, bin, paths) {
  if (.UseCpp()) {
    return(.CheckInfeasible(.DrawBridgesCpp(tree, set, nodes, bin, paths)))
  }
  logQ <- 0
  for (u in nodes) {
    d <- .DrawBridge(tree, set, u, .TopBin(tree, bin, u), bin[[u]])
    paths[[u]] <- d[["path"]]
    logQ <- logQ + d[["logQ"]]
  }
  list(paths = paths, logQ = logQ)
}

# Log bridge density of the paths on branches `nodes`
.LogBridges <- function(tree, set, nodes, paths) {
  if (.UseCpp()) {
    return(.LogBridgesCpp(tree, set, nodes, paths))
  }
  sum(vapply(nodes, function(u) .LogBridge(tree, set, u, paths[[u]]), 0))
}

# Fossil tip: log weight of each tip bin given the top bin, and draws
.TipLogWeights <- function(tree, set, v, kTop, pc = .Pieces(tree, set, v)) {
  if (.UseCpp()) {
    return(.TipLogWeightsCpp(tree, set, v, kTop, pc))
  }
  .LogTransferRow(pc, kTop) + log(.TipValues(tree, set)[v, ])
}

.DrawTip <- function(tree, set, v, kTop) {
  if (.UseCpp()) {
    return(.CheckInfeasible(.DrawTipCpp(tree, set, v, kTop)))
  }
  pc <- .Pieces(tree, set, v)
  lw <- .TipLogWeights(tree, set, v, kTop, pc)
  kTip <- .Sample1(exp(lw - max(lw)))
  path <- .SampleBranch(set, pc, NULL, kTop, kTip)
  list(path = path, bin = kTip,
       logQ = .LogBridge(tree, set, v, path, pc) + lw[[kTip]] - .LogSumExp(lw))
}

.LogTipQ <- function(tree, set, v, path) {
  if (.UseCpp()) {
    return(.LogTipQCpp(tree, set, v, path))
  }
  pc <- .Pieces(tree, set, v)
  k <- path[["k"]]
  lw <- .TipLogWeights(tree, set, v, k[[1]], pc)
  .LogBridge(tree, set, v, path, pc) + lw[[k[[length(k)]]]] - .LogSumExp(lw)
}

.LogSumExp <- function(x) {
  top <- max(x)
  top + log(sum(exp(x - top)))
}

# Log weights of a node's bin given its neighbours: `factor` per bin, the
# branch above (from bin kTop), and branches below to fixed bins `ends`
# (NA: a fossil tip whose bin is marginalized)
.NodeLogWeights <- function(tree, set, v, factor, kTop, kids, ends) {
  if (.UseCpp()) {
    return(.NodeLogWeightsCpp(tree, set, v, factor, kTop, kids, ends))
  }
  lw <- log(factor)
  lw <- lw + .LogTransferRow(.Pieces(tree, set, v), kTop)
  for (i in seq_along(kids)) {
    pc <- .Pieces(tree, set, kids[[i]])
    if (is.na(ends[[i]])) {
      tipD <- log(.TipValues(tree, set)[kids[[i]], ])
      nBin <- length(tipD)
      lw <- lw + vapply(seq_len(nBin), function(l) {
        .LogSumExp(.LogTransferRow(pc, l) + tipD)
      }, 0)
    } else {
      lw <- lw + .LogTransferCol(pc, ends[[i]])
    }
  }
  lw
}

.DrawFromLog <- function(lw) {
  k <- .Sample1(exp(lw - max(lw)))
  list(k = k, logQ = lw[[k]] - .LogSumExp(lw))
}

# Nodes eligible for the subtree slide: parent is a non-root bifurcation
.SlideNodes <- function(tree) {
  nTip <- tree[["nTip"]]
  p <- tree[["parent"]]
  v <- setdiff(.PathNodes(tree), tree[["root"]])
  pv <- p[v]
  v[pv != tree[["root"]] & !tree[["sampled"]][pv - nTip]]
}

# Branches of `tree` (excluding `exclude`) that cross `at`, when node `p` is
# pruned (its other child hangs from `p`'s parent)
.Crossing <- function(tree, at, p, exclude) {
  parent <- tree[["parent"]]
  nodes <- setdiff(.PathNodes(tree), c(exclude, p, tree[["root"]]))
  top <- tree[["age"]][parent[nodes]]
  hangs <- parent[nodes] == p
  top[hangs] <- tree[["age"]][[parent[[p]]]]
  nodes[top > at & tree[["age"]][nodes] < at]
}

.Descendants <- function(tree, v) {
  out <- v
  frontier <- v
  while (length(frontier)) {
    frontier <- which(tree[["parent"]] %in% frontier)
    out <- c(out, frontier)
  }
  out
}

# Propose a tree move; returns list(tree, bin, paths, logHastings) with
# logHastings = log q(reverse) - log q(forward), or NULL if not applicable
.ProposeTreeMove <- function(move, tree, set, st) {
  bin <- st[["bin"]]
  paths <- st[["paths"]]
  parent <- tree[["parent"]]
  age <- tree[["age"]]
  nTip <- tree[["nTip"]]
  Redraw <- function(newTree, nodes, newBin) {
    .DrawBridges(newTree, set, nodes, newBin, paths)
  }
  OldQ <- function(nodes) {
    .LogBridges(tree, set, nodes, paths)
  }
  switch(
    move,
    nodeAge = {
      # With an origin the root is a bifurcation below the stem, and its age
      # moves too (the origin age stays fixed)
      cand <- tree[["postorder"]][!tree[["sampled"]][tree[["postorder"]] - nTip] &
                                     (tree[["postorder"]] != tree[["root"]] |
                                        .Origin(tree))]
      if (!length(cand)) return(NULL)
      v <- cand[[.Sample1(rep(1, length(cand)))]]
      kids <- tree[["kids"]][v - nTip, ]
      newTree <- tree
      top <- .TopAge(tree, v)
      newTree[["age"]][v] <- stats::runif(1, max(age[kids]), top)
      if (v == tree[["root"]]) {
        newTree[["stem"]] <- top - newTree[["age"]][[v]]
      }
      newTree <- .Rebuild(newTree)
      nodes <- c(v, kids)
      fwd <- Redraw(newTree, nodes, bin)
      list(tree = newTree, bin = bin, paths = fwd[["paths"]],
           logHastings = OldQ(nodes) - fwd[["logQ"]])
    },
    fossilAge = {
      range <- tree[["range"]]
      cand <- which(range[, 2] > range[, 1])
      if (!length(cand)) return(NULL)
      f <- cand[[.Sample1(rep(1, length(cand)))]]
      newTree <- tree
      if (tree[["tipType"]][[f]] == "fossil") {
        newTree[["age"]][f] <- stats::runif(1, range[f, 1],
                                            min(range[f, 2], age[[parent[[f]]]]))
        nodes <- f
      } else {
        host <- parent[[f]]
        cont <- setdiff(.PathKids(tree, host), f)
        lo <- max(range[f, 1], age[[cont]])
        hi <- min(range[f, 2], age[[parent[[host]]]])
        newTree[["age"]][host] <- stats::runif(1, lo, hi)
        nodes <- c(host, cont)
      }
      newTree <- .Rebuild(newTree)
      fwd <- Redraw(newTree, nodes, bin)
      list(tree = newTree, bin = bin, paths = fwd[["paths"]],
           logHastings = OldQ(nodes) - fwd[["logQ"]])
    },
    slide = {
      cand <- .SlideNodes(tree)
      if (!length(cand)) return(NULL)
      v <- cand[[.Sample1(rep(1, length(cand)))]]
      p <- parent[[v]]
      pp <- parent[[p]]
      sib <- setdiff(tree[["kids"]][p - nTip, ], v)
      clade <- .Descendants(tree, v)
      rootAge <- age[[tree[["root"]]]]
      at <- stats::runif(1, age[[v]], rootAge)
      fwdTargets <- .Crossing(tree, at, p, clade)
      if (!length(fwdTargets)) return(NULL)
      y <- fwdTargets[[.Sample1(rep(1, length(fwdTargets)))]]
      revTargets <- .Crossing(tree, age[[p]], p, clade)
      newTree <- tree
      newTree[["parent"]][sib] <- pp
      x <- newTree[["parent"]][[y]]
      newTree[["parent"]][p] <- x
      newTree[["parent"]][y] <- p
      newTree[["age"]][p] <- at
      newTree <- .Rebuild(newTree)
      # Bin of p: forward in the new geometry, reverse in the old
      lwNew <- .NodeLogWeights(newTree, set, p, set[["lambda"]], bin[[x]],
                               c(y, v), bin[c(y, v)])
      lwOld <- .NodeLogWeights(tree, set, p, set[["lambda"]], bin[[pp]],
                               c(sib, v), bin[c(sib, v)])
      draw <- .DrawFromLog(lwNew)
      newBin <- bin
      newBin[p] <- draw[["k"]]
      nodesNew <- unique(c(p, y, v, sib))
      fwd <- Redraw(newTree, nodesNew, newBin)
      list(tree = newTree, bin = newBin, paths = fwd[["paths"]],
           logHastings = OldQ(unique(c(p, sib, v, y))) + lwOld[[bin[[p]]]] - .LogSumExp(lwOld) -
             fwd[["logQ"]] - draw[["logQ"]] +
             log(length(fwdTargets)) - log(length(revTargets)) +
             log(length(cand)) - log(length(.SlideNodes(newTree))))
    },
    saJump = {
      fossils <- which(tree[["tipType"]] != "extant")
      if (!length(fossils)) return(NULL)
      f <- fossils[[.Sample1(rep(1, length(fossils)))]]
      p <- parent[[f]]
      if (p == tree[["root"]]) return(NULL)
      pp <- parent[[p]]
      newTree <- tree
      if (tree[["tipType"]][[f]] == "fossil") {
        # Tip to sampled ancestor: p moves down to the fossil's age
        sib <- setdiff(tree[["kids"]][p - nTip, ], f)
        if (tree[["sampled"]][p - nTip] || !(age[[f]] > age[[sib]])) return(NULL)
        newTree[["tipType"]][f] <- "ancestor"
        newTree[["age"]][p] <- age[[f]]
        newTree <- .Rebuild(newTree)
        lwNew <- .NodeLogWeights(newTree, set, p, set[["psi"]], bin[[pp]], sib, bin[[sib]])
        lwOld <- .NodeLogWeights(tree, set, p, set[["lambda"]], bin[[pp]],
                                 c(sib, f), c(bin[[sib]], NA))
        draw <- .DrawFromLog(lwNew)
        newBin <- bin
        newBin[c(p, f)] <- draw[["k"]]
        fwd <- Redraw(newTree, c(p, sib), newBin)
        newPaths <- fwd[["paths"]]
        newPaths[f] <- list(NULL)
        list(tree = newTree, bin = newBin, paths = newPaths,
             logHastings = OldQ(c(p, sib)) + .LogTipQ(tree, set, f, paths[[f]]) +
               lwOld[[bin[[p]]]] - .LogSumExp(lwOld) - log(age[[pp]] - age[[f]]) -
               fwd[["logQ"]] - draw[["logQ"]])
      } else {
        # Sampled ancestor to tip: p rises uniformly towards its parent
        sib <- setdiff(tree[["kids"]][p - nTip, ], f)
        newTree[["tipType"]][f] <- "fossil"
        newTree[["age"]][p] <- stats::runif(1, age[[f]], age[[pp]])
        newTree[["age"]][f] <- age[[f]]
        newTree <- .Rebuild(newTree)
        newTree[["age"]][f] <- age[[f]]
        lwNew <- .NodeLogWeights(newTree, set, p, set[["lambda"]], bin[[pp]],
                                 c(sib, f), c(bin[[sib]], NA))
        lwOld <- .NodeLogWeights(tree, set, p, set[["psi"]], bin[[pp]], sib, bin[[sib]])
        draw <- .DrawFromLog(lwNew)
        newBin <- bin
        newBin[p] <- draw[["k"]]
        fwd <- Redraw(newTree, c(p, sib), newBin)
        tip <- .DrawTip(newTree, set, f, newBin[[p]])
        newBin[f] <- tip[["bin"]]
        newPaths <- fwd[["paths"]]
        newPaths[[f]] <- tip[["path"]]
        list(tree = newTree, bin = newBin, paths = newPaths,
             logHastings = OldQ(c(p, sib)) + lwOld[[bin[[p]]]] - .LogSumExp(lwOld) +
               log(age[[pp]] - age[[f]]) - fwd[["logQ"]] - draw[["logQ"]] -
               tip[["logQ"]])
      }
    }
  )
}

# Summaries of a tree and its paths, comparable with simulated trees;
# `lengths` and `inverse` cover the branches below the root
.TreeStats <- function(tree, lengths = NULL, inverse = NULL) {
  nTip <- tree[["nTip"]]
  age <- tree[["age"]]
  # Branches below the root (with an origin, the stem is excluded)
  nodes <- setdiff(.PathNodes(tree), tree[["root"]])
  internal <- tree[["postorder"]]
  bif <- internal[!tree[["sampled"]][internal - nTip] & internal != tree[["root"]]]
  kids <- tree[["kids"]]
  cherries <- sum(kids[, 1] <= nTip & kids[, 2] <= nTip & !tree[["sampled"]])
  fossilAge <- age[seq_len(nTip)][tree[["tipType"]] != "extant"]
  c(nAncestor = sum(tree[["tipType"]] == "ancestor"),
    treeLength = sum(age[tree[["parent"]][nodes]] - age[nodes]),
    meanNodeAge = if (length(bif)) mean(age[bif]) else NA,
    cherries = cherries,
    oldestFossil = if (length(fossilAge)) max(fossilAge) else NA,
    youngestFossil = if (length(fossilAge)) min(fossilAge) else NA,
    clock = if (is.null(lengths)) NA else sum(lengths),
    inverseClock = if (is.null(inverse)) NA else sum(inverse),
    rootAge = age[[tree[["root"]]]])
}
