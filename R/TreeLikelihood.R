#' Time tree for dCET inference
#'
#' Converts a rooted, binary time tree into the representation used by the
#' dCET likelihood and samplers: an age (time before present) for every node
#' and a type for every tip.
#'
#' Sampled ancestors are tips on zero-length branches; the node that bears
#' one is a sampled node rather than a bifurcation.
#' Ages are taken, in order of preference, from `ages`; from attributes
#' `nodeTime` and `present` (as set by [ReconstructedTree()], with `present`
#' the simulation duration); or from branch lengths, with the deepest tip at
#' age zero.
#'
#' @param tree A rooted binary tree of class `phylo`, with branch lengths.
#' @param ages Optional numeric vector giving the age of each node (tips, then
#' internal nodes, in the numbering of `tree`).
#' @param tipType Optional character vector giving the type of each tip:
#' `"extant"`, `"fossil"` or `"ancestor"`.
#' Defaults to attribute `tipType` of `tree` (with `"extinct"` read as
#' `"fossil"`), or else to `"ancestor"` for tips on zero-length branches,
#' `"extant"` for tips of age zero, and `"fossil"` otherwise.
#' @param present Numeric giving the time of the present on the scale of
#' attribute `nodeTime`.
#' @param tol Numeric giving the tolerance for treating an age or a branch
#' length as zero.
#' @returns `TempoTree()` returns an object of class `TempoTree`: a list with
#' elements `nTip`, `nNode`, `parent` (parent of each node; `NA` at the root),
#' `kids` (matrix of the two children of each internal node), `age`,
#' `tipType`, `sampled` (whether each internal node bears a sampled ancestor),
#' `root` (index of the root node), `postorder` (internal nodes, children
#' before parents) and `stem` (length of the root edge, or zero).
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' TempoTree(tree)$age
#' @family dCET inference
#' @export
TempoTree <- function(tree, ages = NULL, tipType = NULL, present = NULL,
                      tol = 1e-8) {
  if (inherits(tree, "TempoTree")) {
    return(tree)
  }
  if (!inherits(tree, "phylo") || is.null(tree[["edge.length"]])) {
    stop("`tree` must be a phylo object with branch lengths")
  }
  nTip <- length(tree[["tip.label"]])
  nNode <- tree[["Nnode"]]
  edge <- tree[["edge"]]
  nAll <- nTip + nNode
  parent <- rep(NA_integer_, nAll)
  parent[edge[, 2]] <- edge[, 1]
  root <- setdiff(seq_len(nAll), edge[, 2])
  counts <- tabulate(edge[, 1], nAll)
  if (length(root) != 1L || nNode != nTip - 1L ||
      any(counts[nTip + seq_len(nNode)] != 2L)) {
    stop("`tree` must be rooted and binary")
  }
  len <- numeric(nAll)
  len[edge[, 2]] <- tree[["edge.length"]]
  order <- .Preorder(parent, root)
  if (is.null(ages)) {
    nodeTime <- attr(tree, "nodeTime")
    if (!is.null(nodeTime) && !is.null(present)) {
      ages <- present - nodeTime
    } else {
      depth <- numeric(nAll)
      for (v in order[-1]) {
        depth[v] <- depth[parent[v]] + len[v]
      }
      ages <- max(depth[seq_len(nTip)]) - depth
    }
  }
  ages[abs(ages) < tol] <- 0
  if (length(ages) != nAll || any(ages < 0) ||
      any(ages[edge[, 1]] - ages[edge[, 2]] < -tol)) {
    stop("`ages` must be non-negative and decrease from root to tips")
  }
  if (is.null(tipType)) {
    tipType <- attr(tree, "tipType")
    if (is.null(tipType)) {
      tipType <- ifelse(len[seq_len(nTip)] < tol, "ancestor",
                        ifelse(ages[seq_len(nTip)] == 0, "extant", "fossil"))
    }
    tipType[tipType == "extinct"] <- "fossil"
  }
  if (length(tipType) != nTip ||
      !all(tipType %in% c("extant", "fossil", "ancestor"))) {
    stop("`tipType` must give \"extant\", \"fossil\" or \"ancestor\" per tip")
  }
  if (any(ages[seq_len(nTip)][tipType == "extant"] != 0)) {
    stop("Extant tips must have age zero")
  }
  kids <- matrix(0L, nNode, 2)
  for (i in seq_len(nNode)) {
    kids[i, ] <- edge[edge[, 1] == nTip + i, 2]
  }
  isAncestor <- c(tipType == "ancestor", logical(nNode))
  sampled <- isAncestor[kids[, 1]] | isAncestor[kids[, 2]]
  if (any(isAncestor[kids[, 1]] & isAncestor[kids[, 2]])) {
    stop("A node cannot bear two sampled ancestors")
  }
  # Sampled-ancestor tips sit at their node's age
  saTip <- which(isAncestor)
  ages[saTip] <- ages[parent[saTip]]
  structure(
    list(nTip = nTip, nNode = nNode, parent = parent, kids = kids, age = ages,
         tipType = tipType, sampled = sampled, root = root,
         postorder = rev(order[order > nTip]),
         stem = if (is.null(tree[["root.edge"]])) 0 else tree[["root.edge"]],
         tipLabel = tree[["tip.label"]]),
    class = "TempoTree"
  )
}

.Preorder <- function(parent, root) {
  children <- split(seq_along(parent), factor(parent, levels = seq_along(parent)))
  out <- integer(length(parent))
  stack <- root
  i <- 0L
  while (length(stack)) {
    v <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    i <- i + 1L
    out[i] <- v
    stack <- c(stack, rev(children[[v]]))
  }
  out
}

#' dCET model specification
#'
#' Fixed settings of a dCET analysis: the tempo grid and discretization, the
#' couplings, the sampling scheme, the start-bin prior and the conditioning.
#'
#' The likelihood is the density of the oriented reconstructed tree with
#' labelled tips (spec section 3.4): a factor \eqn{\lambda_k} at each
#' bifurcation, \eqn{\rho} at each extant tip,
#' \eqn{\psi_k[r + (1 - r)E_k(\tau)]} at each fossil tip, and
#' \eqn{\psi_k(1 - r)} at each sampled ancestor, with \eqn{r} the removal
#' probability.
#' With `start = "root"` the root bifurcation carries no \eqn{\lambda} factor
#' and both daughters start in the root bin; with `start = "origin"` a stem
#' branch of length `root.edge` precedes the root, which carries a
#' \eqn{\lambda} factor.
#'
#' @param nBin Integer giving the number of tempo bins.
#' @param zRange Numeric vector giving the range of bin centres in units of
#' `s`.
#' @inheritParams TempoGenerator
#' @param beta,gamma,delta Numerics giving the coupling to tempo of
#' diversification, of character change and of fossil sampling.
#' @param rho Numeric giving the probability of sampling an extant lineage.
#' @param removal Numeric giving the probability that a sampled fossil lineage
#' is removed.
#' @param rootPrior Character string or numeric vector specifying the prior
#' on the start bin:
#'  - `equilibrium`: Perron vector of \eqn{Q^\top + R}, the tempo distribution
#'    of contemporaneous species at equilibrium (spec section 5.4, option 1);
#'  - `stationary`: lineage-stationary law of \eqn{Q};
#'  - `uniform`: equal weights;
#'  - a numeric vector of weights, one per bin.
#' @param condition Character string specifying the survival conditioning:
#'  - `extant`: each root daughter (or the stem) has a sampled extant
#'    descendant (trees without one have likelihood zero);
#'  - `survival`: each has any sampled descendant;
#'  - `none`: no conditioning.
#' @param joint Logical; if `TRUE`, condition on survival given the start bin
#' (joint conditioning); otherwise divide by the prior-averaged survival
#' probability (marginal conditioning).
#' @param start Character string: `root` or `origin` (see Details).
#' @returns `DCETModel()` returns an object of class `DCETModel`.
#' @examples
#' DCETModel(beta = 0)
#' @family dCET inference
#' @export
DCETModel <- function(nBin = 26, zRange = c(-5, 5),
                      scheme = c("central", "upwind", "sqra"), alpha = 1,
                      beta = 1, gamma = 1, delta = 0, rho = 1, removal = 0,
                      rootPrior = "equilibrium",
                      condition = c("extant", "survival", "none"),
                      joint = FALSE, start = c("root", "origin")) {
  if (is.character(rootPrior)) {
    rootPrior <- match.arg(rootPrior, c("equilibrium", "stationary", "uniform"))
  } else if (length(rootPrior) != nBin || any(rootPrior < 0)) {
    stop("A numeric `rootPrior` needs one non-negative weight per bin")
  }
  structure(
    list(nBin = as.integer(nBin), zRange = zRange, scheme = match.arg(scheme),
         alpha = alpha, beta = beta, gamma = gamma, delta = delta, rho = rho,
         removal = removal, rootPrior = rootPrior,
         condition = match.arg(condition), joint = joint,
         start = match.arg(start)),
    class = "DCETModel"
  )
}

# Grid, generator, per-bin rates and start prior at parameters `par`
.Setup <- function(par, model) {
  grid <- TempoGrid(model[["nBin"]], model[["zRange"]], par[["s"]])
  Q <- TempoGenerator(grid, par[["theta"]], model[["alpha"]], model[["scheme"]])
  rates <- BinRates(grid, par[["lambda"]], par[["mu"]],
                    if (is.null(par[["psi"]])) 0 else par[["psi"]],
                    model[["beta"]], model[["delta"]])
  prior <- .RootPrior(model[["rootPrior"]], Q, rates, model[["removal"]])
  c(list(grid = grid, Q = Q, prior = prior), rates)
}

.RootPrior <- function(rootPrior, Q, rates, removal = 0) {
  if (is.numeric(rootPrior)) {
    return(rootPrior / sum(rootPrior))
  }
  switch(rootPrior,
         equilibrium = TempoEquilibrium(
           Q, rates[["lambda"]], rates[["mu"]] + removal * rates[["psi"]])$p,
         stationary = TempoEquilibrium(Q)$p,
         uniform = rep(1 / nrow(Q), nrow(Q)))
}

#' Marginal dCET tree likelihood
#'
#' Exact tempo-marginal likelihood of a fixed time tree under dCET, with
#' extant sampling, fossil tips, sampled ancestors and fossil sampling.
#' The backward equations for \eqn{E} and \eqn{D} (spec equations 2.1 and 3.1)
#' are integrated jointly along each branch with a stiff solver, and combined
#' at nodes by the SSE pruning recursion (spec section 3.4).
#'
#' @param tree A [`TempoTree`][TempoTree()], or a tree that it accepts.
#' @param par Named list or numeric vector of parameters: `lambda`, `mu`,
#' `theta`, `s`, and optionally `psi` (default 0).
#' @param model A [`DCETModel`][DCETModel()].
#' @param rtol,atol Numerics giving the tolerances of the stiff solver.
#' @returns `TreeLikelihood()` returns the natural logarithm of the
#' likelihood.
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' TreeLikelihood(tree, list(lambda = 0.5, mu = 0.2, theta = 0.1, s = 1),
#'                DCETModel(nBin = 11, zRange = c(-2, 2)))
#' @family dCET inference
#' @export
TreeLikelihood <- function(tree, par, model = DCETModel(), rtol = 1e-10,
                           atol = 1e-14) {
  tree <- TempoTree(tree)
  par <- as.list(par)
  set <- .Setup(par, model)
  .PruneTree(tree, set[["Q"]], set[["lambda"]], set[["mu"]], set[["psi"]],
             model[["rho"]], model[["removal"]], set[["prior"]],
             model[["condition"]], model[["joint"]], model[["start"]],
             rtol, atol)
}

# Whether the tree is possible under the conditioning: with extant (crown)
# conditioning each root daughter, or with an origin the tree, needs a
# sampled extant descendant
.Compatible <- function(tree, condition, origin = FALSE) {
  if (condition != "extant") {
    return(TRUE)
  }
  extant <- which(tree[["tipType"]] == "extant")
  if (origin) {
    return(length(extant) > 0)
  }
  root <- tree[["root"]]
  # Child of the root above each extant tip
  top <- vapply(extant, function(v) {
    while (tree[["parent"]][[v]] != root) {
      v <- tree[["parent"]][[v]]
    }
    v
  }, 0)
  all(tree[["kids"]][root - tree[["nTip"]], ] %in% top)
}

# Exact pruning for arbitrary generator `Q` and per-bin rates
.PruneTree <- function(tree, Q, lambda, mu, psi, rho, removal, prior,
                       condition, joint, start, rtol = 1e-10, atol = 1e-14) {
  nBin <- nrow(Q)
  lambda <- rep_len(lambda, nBin)
  mu <- rep_len(mu, nBin)
  psi <- rep_len(psi, nBin)
  nTip <- tree[["nTip"]]
  age <- tree[["age"]]
  origin <- start == "origin"
  if (origin && !(tree[["stem"]] > 0)) {
    stop("`start = \"origin\"` needs a tree with a positive root edge")
  }
  if (tree[["sampled"]][tree[["root"]] - nTip] && !origin) {
    stop("The root cannot bear a sampled ancestor unless `start = \"origin\"`")
  }
  if (any(tree[["tipType"]] != "extant") && !all(psi > 0)) {
    stop("Fossil tips and sampled ancestors need `psi > 0`")
  }
  if (!.Compatible(tree, condition, origin)) {
    return(-Inf)
  }
  topAge <- age[tree[["root"]]] + if (origin) tree[["stem"]] else 0
  sorted <- sort(unique(c(age, topAge)))
  eAt <- ExtinctionProbability(sorted, Q, lambda, mu, psi, rho, rtol, atol)
  E <- function(a) eAt[match(a, sorted), ]
  D <- matrix(0, tree[["nTip"]] + tree[["nNode"]], nBin)
  logScale <- 0
  tipType <- tree[["tipType"]]
  D[which(tipType == "extant"), ] <- rho
  for (f in which(tipType == "fossil")) {
    D[f, ] <- psi * (removal + (1 - removal) * E(age[[f]]))
  }
  sigma <- lambda + mu + psi
  Branch <- function(v, toAge) {
    .BranchODE(D[v, ], E(age[[v]]), toAge - age[[v]], Q, lambda, mu, sigma,
               rtol, atol)
  }
  for (v in tree[["postorder"]]) {
    kids <- tree[["kids"]][v - nTip, ]
    value <- if (tree[["sampled"]][v - nTip]) {
      isSA <- kids <= nTip & tipType[pmin(kids, nTip)] == "ancestor"
      psi * (1 - removal) * Branch(kids[!isSA], age[[v]])
    } else {
      Branch(kids[[1]], age[[v]]) * Branch(kids[[2]], age[[v]]) *
        (if (v == tree[["root"]] && !origin) 1 else lambda)
    }
    scale <- max(value)
    if (!(scale > 0)) {
      return(-Inf)
    }
    D[v, ] <- value / scale
    logScale <- logScale + log(scale)
  }
  top <- D[tree[["root"]], ]
  if (origin) {
    top <- Branch(tree[["root"]], topAge)
  }
  survive <- if (condition == "none") {
    rep(1, nBin)
  } else {
    pr <- if (condition == "extant") {
      1 - ExtinctionProbability(topAge, Q, lambda, mu + removal * psi, 0, rho,
                                rtol, atol)[1, ]
    } else {
      1 - E(topAge)
    }
    if (origin) pr else pr * pr
  }
  if (joint) {
    log(sum(prior * top / survive)) + logScale
  } else {
    log(sum(prior * top)) + logScale - log(sum(prior * survive))
  }
}

# Integrate (E, D) jointly over age span `span` from D0, E0
.BranchODE <- function(D0, E0, span, Q, lambda, mu, sigma, rtol, atol) {
  if (span <= 0 || !(max(D0) > 0)) {
    return(D0)
  }
  nBin <- length(D0)
  eIdx <- seq_len(nBin)
  dIdx <- nBin + eIdx
  scale <- max(D0)
  out <- .Integrate(
    c(E0, D0 / scale), span,
    function(t, y) {
      e <- y[eIdx]
      d <- y[dIdx]
      g <- 2 * lambda * e - sigma
      c(mu - sigma * e + lambda * e * e + as.vector(Q %*% e),
        g * d + as.vector(Q %*% d))
    },
    function(t, y) {
      e <- y[eIdx]
      jac <- matrix(0, 2 * nBin, 2 * nBin)
      g <- 2 * lambda * e - sigma
      jac[eIdx, eIdx] <- .AddDiag(Q, g)
      jac[dIdx, dIdx] <- .AddDiag(Q, g)
      jac[cbind(dIdx, eIdx)] <- 2 * lambda * y[dIdx]
      jac
    },
    rtol, atol
  )
  pmax(out[1, dIdx], 0) * scale
}

#' Constant-rate fossilized birth-death tree likelihood
#'
#' Closed-form likelihood of a time tree under the constant-rate fossilized
#' birth-death process (Stadler 2010; Gavryushkina et al. 2014), in the same
#' labelling and conditioning conventions as [TreeLikelihood()].
#' The dCET likelihood with `beta = 0` and `delta = 0` equals this for any
#' tempo generator (spec section 7.2).
#'
#' Each branch from age \eqn{\tau_c} to \eqn{\tau_p} contributes
#' \eqn{q(\tau_c) / q(\tau_p)}, with
#' \eqn{q(t) = 2(1 - c_2^2) + e^{-c_1 t}(1 - c_2)^2 + e^{c_1 t}(1 + c_2)^2}.
#'
#' @inheritParams TreeLikelihood
#' @inheritParams FBDExtinction
#' @inheritParams DCETModel
#' @returns `FBDTreeLikelihood()` returns the natural logarithm of the
#' likelihood.
#' @references
#' Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014). Bayesian inference
#' of sampled ancestor trees for epidemiology and fossil calibration. PLoS
#' Computational Biology 10:e1003919. \doi{10.1371/journal.pcbi.1003919}
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' FBDTreeLikelihood(tree, lambda = 0.5, mu = 0.2)
#' @family dCET inference
#' @export
FBDTreeLikelihood <- function(tree, lambda, mu, psi = 0, rho = 1, removal = 0,
                              condition = c("extant", "survival", "none"),
                              start = c("root", "origin")) {
  tree <- TempoTree(tree)
  condition <- match.arg(condition)
  origin <- match.arg(start) == "origin"
  if (!.Compatible(tree, condition, origin)) {
    return(-Inf)
  }
  nTip <- tree[["nTip"]]
  age <- tree[["age"]]
  tipType <- tree[["tipType"]]
  LogQ <- function(t) {
    c1 <- sqrt((lambda - mu - psi)^2 + 4 * lambda * psi)
    c2 <- -(lambda - mu - 2 * lambda * rho - psi) / c1
    log(2 * (1 - c2^2) + exp(-c1 * t) * (1 - c2)^2 + exp(c1 * t) * (1 + c2)^2)
  }
  child <- which(!is.na(tree[["parent"]]) &
                   !(seq_along(age) <= nTip & c(tipType, character(tree[["nNode"]])) == "ancestor"))
  root <- tree[["root"]]
  topAge <- age[[root]] + if (origin) tree[["stem"]] else 0
  branches <- sum(LogQ(age[child]) - LogQ(age[tree[["parent"]][child]])) +
    if (origin) LogQ(age[[root]]) - LogQ(topAge) else 0
  fossil <- age[seq_len(nTip)][tipType == "fossil"]
  nBif <- sum(!tree[["sampled"]]) - if (origin) 0 else 1
  survive <- switch(condition,
                    none = 1,
                    extant = 1 - FBDExtinction(topAge, lambda, mu + removal * psi, 0, rho),
                    survival = 1 - FBDExtinction(topAge, lambda, mu, psi, rho))
  branches + nBif * log(lambda) + sum(tipType == "extant") * log(rho) +
    sum(log(psi * (removal + (1 - removal) *
                     FBDExtinction(fossil, lambda, mu, psi, rho)))) +
    (if (any(tipType == "ancestor")) {
      sum(tipType == "ancestor") * log(psi * (1 - removal))
    } else {
      0
    }) -
    (if (condition == "none") 0 else if (origin) log(survive) else 2 * log(survive))
}
