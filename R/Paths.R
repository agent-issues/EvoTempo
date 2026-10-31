# Tempo paths on a fixed time tree: exact augmented density, and proposals
# from a model whose E is frozen on time slices (spec section 8).
#
# A path on the branch above node v is list(t, k): jump ages `t` (decreasing,
# strictly inside the branch) and bins `k` (length(t) + 1), with k[1] the bin
# of the parent node and k[length(k)] that of v.

# E and G = integral of E, on an age grid including `extra`, with derivatives
# for cubic Hermite interpolation
.EGrid <- function(Q, lambda, mu, psi, rho, maxAge, extra = numeric(0),
                   rtol = 1e-10, atol = 1e-14) {
  nBin <- nrow(Q)
  ages <- sort(unique(c(seq(0, maxAge, length.out = 401),
                        maxAge * 10^seq(-6, 0, length.out = 61), extra)))
  if (.UseCpp()) {
    eg <- .SolverInfeasible(
      .EGridCpp(ages, Q, lambda, mu, psi, rho, rtol, atol)
    )
    return(c(list(ages = ages), eg))
  }
  sigma <- lambda + mu + psi
  Rhs <- function(e) mu - sigma * e + lambda * e * e + as.vector(Q %*% e)
  eIdx <- seq_len(nBin)
  out <- .SolverInfeasible(.Integrate(
    c(rep(1 - rho, nBin), numeric(nBin)), ages,
    function(t, y) c(Rhs(y[eIdx]), y[eIdx]),
    function(t, y) {
      jac <- matrix(0, 2 * nBin, 2 * nBin)
      jac[eIdx, eIdx] <- .AddDiag(Q, 2 * lambda * y[eIdx] - sigma)
      jac[cbind(nBin + eIdx, eIdx)] <- 1
      jac
    },
    rtol, atol
  ))
  E <- out[, eIdx, drop = FALSE]
  list(ages = ages, E = E, G = out[, nBin + eIdx, drop = FALSE],
       dE = t(apply(E, 1, Rhs)))
}

# Interpolated E or G at ages `at` in bins `bin`
.Hermite <- function(eg, at, bin, value = "G") {
  ages <- eg[["ages"]]
  i <- findInterval(at, ages, all.inside = TRUE)
  h <- ages[i + 1L] - ages[i]
  u <- (at - ages[i]) / h
  if (value == "G") {
    y <- eg[["G"]]
    dy <- eg[["E"]]
  } else {
    y <- eg[["E"]]
    dy <- eg[["dE"]]
  }
  i0 <- cbind(i, bin)
  i1 <- cbind(i + 1L, bin)
  u2 <- u * u
  u3 <- u2 * u
  (2 * u3 - 3 * u2 + 1) * y[i0] + (u3 - 2 * u2 + u) * h * dy[i0] +
    (3 * u2 - 2 * u3) * y[i1] + (u3 - u2) * h * dy[i1]
}

.EVector <- function(eg, at) {
  nBin <- ncol(eg[["E"]])
  .Hermite(eg, rep(at, nBin), seq_len(nBin), "E")
}

# Whether the paths start at an origin above the root, with a path on the
# stem (`start = "origin"`; set by .ForModel())
.Origin <- function(tree) {
  isTRUE(tree[["origin"]])
}

# The tree as the samplers use it under `model`: with `start = "origin"`, the
# root edge (`stem`) carries a path from the origin, whose bin is stored as
# element nTip + nNode + 1 of the bin vector
.ForModel <- function(tree, model) {
  origin <- model[["start"]] == "origin"
  if (origin) {
    if (!(tree[["stem"]] > 0)) {
      stop("`start = \"origin\"` needs a tree with a positive root edge")
    }
    if (tree[["sampled"]][tree[["root"]] - tree[["nTip"]]]) {
      stop("The samplers do not support a sampled ancestor at the root")
    }
  }
  tree[["origin"]] <- if (origin) TRUE
  tree
}

# Age at which the process starts: the origin, or the root
.StartAge <- function(tree) {
  tree[["age"]][[tree[["root"]]]] + if (.Origin(tree)) tree[["stem"]] else 0
}

# Age of the top of the branch above each node in `v`: its parent's, or for
# the root (with an origin) the origin's
.TopAge <- function(tree, v) {
  root <- tree[["root"]]
  top <- tree[["age"]][tree[["parent"]][v]]
  top[v == root] <- tree[["age"]][[root]] + tree[["stem"]]
  top
}

# Bin at the top of the branch above node v: its parent's, or for the root
# the origin bin
.TopBin <- function(tree, bin, v) {
  if (v == tree[["root"]]) {
    bin[[tree[["nTip"]] + tree[["nNode"]] + 1L]]
  } else {
    bin[[tree[["parent"]][[v]]]]
  }
}

# Model quantities at parameters `par`, for the augmented density
.Precompute <- function(tree, par, model) {
  set <- .Setup(par, model)
  Q <- set[["Q"]]
  nBin <- nrow(Q)
  rho <- model[["rho"]]
  removal <- model[["removal"]]
  rootAge <- .StartAge(tree)
  # Conditioning on both root daughters, or on the single stem lineage
  power <- if (.Origin(tree)) 1 else 2
  eg <- .EGrid(Q, set[["lambda"]], set[["mu"]], set[["psi"]], rho, rootAge,
               tree[["age"]])
  set[["eg"]] <- eg
  set[["joint"]] <- model[["joint"]]
  set[["condition"]] <- model[["condition"]]
  set[["removal"]] <- removal
  set[["rho"]] <- rho
  set[["sigma"]] <- set[["lambda"]] + set[["mu"]] + set[["psi"]]
  set[["survive"]] <- switch(
    model[["condition"]],
    none = rep(1, nBin),
    extant = if (all(set[["psi"]] == 0)) {
      (1 - .EVector(eg, rootAge))^power
    } else {
      eRoot <- .SolverInfeasible(ExtinctionProbability(
        rootAge, Q, set[["lambda"]], set[["mu"]] + removal * set[["psi"]], 0,
        rho
      ))
      (1 - eRoot[1, ])^power
    },
    survival = (1 - .EVector(eg, rootAge))^power
  )
  set[["key"]] <- .SetKey(tree, par)
  set
}

# What a set from .Precompute() depends on besides the model: the parameters
# that enter E and the generator (not `rate`), the node ages, which are
# points of the E grid, and the start age.  Equal keys mean identical sets.
.SetKey <- function(tree, par) {
  list(par = vapply(c("lambda", "mu", "psi", "theta", "s"), function(p) {
    if (is.null(par[[p]])) 0 else as.numeric(par[[p]])
  }, 0), age = tree[["age"]], start = .StartAge(tree))
}

# Whether `set` is the set of `tree` at parameters `par`
.SameSet <- function(set, tree, par) {
  identical(set[["key"]], .SetKey(tree, par))
}

# Frozen-E slices: symmetrized generator Q + diag(g) per slice
.AddSlices <- function(tree, set, nSlice) {
  Q <- set[["Q"]]
  nBin <- nrow(Q)
  lower <- seq_len(nBin - 1L)
  up <- Q[cbind(lower, lower + 1L)]
  down <- Q[cbind(lower + 1L, lower)]
  logPi <- c(0, cumsum(log(up) - log(down)))
  set[["halfPi"]] <- exp((logPi - max(logPi)) / 2)
  off <- sqrt(up * down)
  bounds <- seq(0, .StartAge(tree), length.out = nSlice + 1L)
  set[["bounds"]] <- bounds
  if (.UseCpp()) {
    set[["slices"]] <- .SlicesCpp(set, (bounds[-1] + bounds[-length(bounds)]) / 2,
                                  off)
    return(set)
  }
  set[["slices"]] <- lapply(seq_len(nSlice), function(s) {
    eBar <- .EVector(set[["eg"]], (bounds[s] + bounds[s + 1L]) / 2)
    sym <- diag(diag(Q) + 2 * set[["lambda"]] * eBar - set[["sigma"]])
    sym[cbind(lower, lower + 1L)] <- off
    sym[cbind(lower + 1L, lower)] <- off
    eig <- eigen(sym, symmetric = TRUE)
    list(eBar = eBar, values = eig[["values"]], vectors = eig[["vectors"]])
  })
  set
}

# Frozen transfer matrix of slice `s` over `span`
.SliceExp <- function(set, s, span) {
  if (.UseCpp()) {
    return(.SliceExpCpp(set, s, span))
  }
  sl <- set[["slices"]][[s]]
  hp <- set[["halfPi"]]
  v <- sl[["vectors"]]
  (v %*% (exp(sl[["values"]] * span) * t(v))) * outer(1 / hp, hp)
}

.TipValues <- function(tree, set) {
  rho <- set[["rho"]]
  removal <- set[["removal"]]
  nBin <- nrow(set[["Q"]])
  out <- matrix(0, tree[["nTip"]], nBin)
  for (i in seq_len(tree[["nTip"]])) {
    out[i, ] <- switch(
      tree[["tipType"]][[i]],
      extant = rho,
      fossil = set[["psi"]] * (removal + (1 - removal) *
                                 .EVector(set[["eg"]], tree[["age"]][[i]])),
      ancestor = 0
    )
  }
  out
}

# Branches carrying paths: every non-root node except sampled-ancestor tips,
# and the root when the stem starts at an origin
.PathNodes <- function(tree) {
  nAll <- tree[["nTip"]] + tree[["nNode"]]
  v <- seq_len(nAll)
  v[(!is.na(tree[["parent"]]) | (.Origin(tree) & v == tree[["root"]])) &
      !(v <= tree[["nTip"]] & c(tree[["tipType"]] == "ancestor",
                                logical(tree[["nNode"]])))]
}

# Pieces of the branch above v: boundaries (old to young) and slice indices
.Pieces <- function(tree, set, v) {
  if (.UseCpp()) {
    return(.PiecesCpp(tree, set, v))
  }
  top <- .TopAge(tree, v)
  bottom <- tree[["age"]][[v]]
  b <- set[["bounds"]]
  inside <- b[b > bottom + 1e-12 & b < top - 1e-12]
  edges <- c(top, rev(inside), bottom)
  mids <- (edges[-1] + edges[-length(edges)]) / 2
  slice <- pmin(findInterval(mids, b), length(set[["slices"]]))
  list(edges = edges, slice = slice,
       M = lapply(seq_along(slice), function(j) {
         .SliceExp(set, slice[[j]], edges[[j]] - edges[[j + 1L]])
       }))
}

# Postorder frozen messages: D at each node, and per branch the backward
# message at each piece boundary; log of the frozen likelihood.  With an
# origin, `top` is the message at the origin.
.Messages <- function(tree, set) {
  if (.UseCpp()) {
    return(.MessagesCpp(tree, set))
  }
  nTip <- tree[["nTip"]]
  nAll <- nTip + tree[["nNode"]]
  nBin <- nrow(set[["Q"]])
  D <- matrix(0, nAll, nBin)
  D[seq_len(nTip), ] <- .TipValues(tree, set)
  pieces <- vector("list", nAll)
  beta <- vector("list", nAll)
  logScale <- 0
  removal <- set[["removal"]]
  origin <- .Origin(tree)
  Up <- function(v) {
    pc <- .Pieces(tree, set, v)
    m <- length(pc[["M"]])
    b <- vector("list", m + 1L)
    b[[m + 1L]] <- D[v, ]
    for (j in rev(seq_len(m))) {
      x <- as.vector(pc[["M"]][[j]] %*% b[[j + 1L]])
      b[[j]] <- x / max(x)
      logScale <<- logScale + log(max(x))
    }
    pieces[[v]] <<- pc
    beta[[v]] <<- b
    b[[1]]
  }
  for (v in tree[["postorder"]]) {
    kids <- tree[["kids"]][v - nTip, ]
    if (tree[["sampled"]][v - nTip]) {
      cont <- kids[!(kids <= nTip & tree[["tipType"]][pmin(kids, nTip)] == "ancestor")]
      value <- set[["psi"]] * (1 - removal) * Up(cont)
    } else {
      value <- Up(kids[[1]]) * Up(kids[[2]]) *
        (if (v == tree[["root"]] && !origin) 1 else set[["lambda"]])
    }
    D[v, ] <- value / max(value)
    logScale <- logScale + log(max(value))
  }
  top <- if (origin) Up(tree[["root"]])
  out <- list(D = D, pieces = pieces, beta = beta, logScale = logScale)
  if (origin) {
    out[["top"]] <- top
  }
  out
}

# Frozen message at the start: at the origin, or at the root
.StartMessage <- function(tree, msg) {
  if (.Origin(tree)) msg[["top"]] else msg[["D"]][tree[["root"]], ]
}

# Frozen likelihood (for testing the importance identity)
.FrozenLikelihood <- function(tree, set, msg, joint = FALSE) {
  top <- .StartMessage(tree, msg)
  prior <- set[["prior"]]
  if (joint) {
    log(sum(prior * top / set[["survive"]])) + msg[["logScale"]]
  } else {
    log(sum(prior * top)) + msg[["logScale"]] - log(sum(prior * set[["survive"]]))
  }
}

# Signal a proposal that cannot be drawn numerically; callers reject it
.Infeasible <- function(message = "Infeasible bridge") {
  stop(structure(class = c("infeasibleProposal", "error", "condition"),
                 list(message = message, call = NULL)))
}

# A failed stiff solve as an infeasible proposal; other errors pass unchanged
.SolverInfeasible <- function(expr) {
  tryCatch(expr, error = function(e) {
    if (identical(conditionMessage(e), "ODE solver failed")) {
      .Infeasible(conditionMessage(e))
    }
    stop(e)
  })
}

# One draw from weights `w` (negative round-off treated as zero)
.Sample1 <- function(w) {
  cw <- cumsum(pmax(w, 0))
  total <- cw[[length(cw)]]
  if (!is.finite(total) || !(total > 0)) {
    # Weights underflow when fixed end bins cannot be connected numerically
    .Infeasible("No positive finite weight to sample from")
  }
  findInterval(stats::runif(1) * total, cw) + 1L
}

# Endpoint-conditioned path under the frozen generator of a slice, over
# `span`, from bin i to bin j, by direct sampling (Hobolth & Stone 2009) on
# the symmetrized eigendecomposition.  Jump times are measured forward from
# the old end.
.Bridge <- function(set, s, span, i, j) {
  if (.UseCpp()) {
    return(.CheckInfeasible(.BridgeCpp(set, s, span, i, j)))
  }
  sl <- set[["slices"]][[s]]
  lam <- sl[["values"]]
  V <- sl[["vectors"]]
  hp <- set[["halfPi"]]
  Q <- set[["Q"]]
  diagG <- diag(Q) + (2 * set[["lambda"]] * sl[["eBar"]] - set[["sigma"]])
  times <- numeric(0)
  bins <- i
  a <- i
  now <- 0
  repeat {
    remain <- span - now
    gaa <- diagG[[a]]
    rateDiff <- gaa - lam
    # Integral over [0, tau] of exp(gaa u + lam (remain - u)), per eigenvalue
    J <- function(tau) {
      ifelse(rateDiff > 0,
             exp(gaa * tau + lam * (remain - tau)) * -expm1(-rateDiff * tau) / rateDiff,
             ifelse(rateDiff < 0,
                    exp(lam * remain) * -expm1(rateDiff * tau) / -rateDiff,
                    tau * exp(lam * remain)))
    }
    nbr <- which(Q[a, ] > 0)
    coef <- lapply(nbr, function(c) Q[a, c] * hp[[j]] / hp[[c]] * V[c, ] * V[j, ])
    jAll <- J(remain)
    w <- c(if (a == j) exp(gaa * remain) else 0,
           vapply(coef, function(cf) sum(cf * jAll), 0))
    if (!(sum(pmax(w, 0)) > 0)) {
      # Endpoints too far apart for the remaining time: the bridge density
      # underflows and cannot be sampled
      .Infeasible()
    }
    pick <- .Sample1(w)
    if (pick == 1L) {
      break
    }
    cf <- coef[[pick - 1L]]
    target <- stats::runif(1) * w[[pick]]
    tau <- stats::uniroot(function(x) sum(cf * J(x)) - target, c(0, remain),
                          tol = 1e-12 * span)[["root"]]
    now <- now + tau
    a <- nbr[[pick - 1L]]
    times <- c(times, now)
    bins <- c(bins, a)
  }
  list(t = times, k = bins)
}

# Sample the path on branch v given the bin at its top; the bottom bin is
# drawn from the backward messages (`end = NULL`) or fixed.
.SampleBranch <- function(set, pc, beta, kTop, end = NULL) {
  if (.UseCpp()) {
    return(.CheckInfeasible(.SampleBranchCpp(set, pc, beta, kTop, end)))
  }
  m <- length(pc[["M"]])
  if (!is.null(end)) {
    # Backward messages from a fixed end bin
    nBin <- nrow(set[["Q"]])
    beta <- vector("list", m + 1L)
    beta[[m + 1L]] <- replace(numeric(nBin), end, 1)
    for (j in rev(seq_len(m))) {
      x <- as.vector(pc[["M"]][[j]] %*% beta[[j + 1L]])
      beta[[j]] <- x / max(x)
    }
  }
  b <- integer(m + 1L)
  b[1] <- kTop
  for (j in seq_len(m)) {
    b[j + 1L] <- .Sample1(pc[["M"]][[j]][b[j], ] * beta[[j + 1L]])
  }
  jumpT <- numeric(0)
  bins <- kTop
  edges <- pc[["edges"]]
  for (j in seq_len(m)) {
    span <- edges[[j]] - edges[[j + 1L]]
    br <- .Bridge(set, pc[["slice"]][[j]], span, b[j], b[j + 1L])
    jumpT <- c(jumpT, edges[[j]] - br[["t"]])
    bins <- c(bins, br[["k"]][-1])
  }
  list(t = jumpT, k = bins)
}

# Sample bins and paths in the clade below the top of branch v
.SampleClade <- function(tree, set, msg, paths, bin, v, kTop) {
  if (.UseCpp()) {
    return(.CheckInfeasible(.SampleCladeCpp(tree, set, msg, paths, bin, v, kTop)))
  }
  nTip <- tree[["nTip"]]
  path <- .SampleBranch(set, msg[["pieces"]][[v]], msg[["beta"]][[v]], kTop)
  paths[[v]] <- path
  bin[v] <- path[["k"]][[length(path[["k"]])]]
  if (v > nTip) {
    for (kid in tree[["kids"]][v - nTip, ]) {
      if (kid <= nTip && tree[["tipType"]][[kid]] == "ancestor") {
        bin[kid] <- bin[v]
      } else {
        out <- .SampleClade(tree, set, msg, paths, bin, kid, bin[v])
        paths <- out[["paths"]]
        bin <- out[["bin"]]
      }
    }
  }
  list(paths = paths, bin = bin)
}

# Draw every bin and path from the frozen tree-conditional law
.SampleTree <- function(tree, set, msg) {
  if (.UseCpp()) {
    return(.CheckInfeasible(.SampleTreeCpp(tree, set, msg)))
  }
  nAll <- tree[["nTip"]] + tree[["nNode"]]
  root <- tree[["root"]]
  origin <- .Origin(tree)
  paths <- vector("list", nAll)
  bin <- integer(nAll + origin)
  w <- set[["prior"]] * .StartMessage(tree, msg)
  if (set[["joint"]]) {
    w <- w / set[["survive"]]
  }
  if (origin) {
    bin[nAll + 1L] <- .Sample1(w)
    return(.SampleClade(tree, set, msg, paths, bin, root, bin[[nAll + 1L]]))
  }
  bin[root] <- .Sample1(w)
  out <- list(paths = paths, bin = bin)
  for (kid in tree[["kids"]][root - tree[["nTip"]], ]) {
    out <- .SampleClade(tree, set, msg, out[["paths"]], out[["bin"]], kid,
                        out[["bin"]][root])
  }
  out
}

# Intervals of a path: bins with old (hi) and young (lo) ends
.Intervals <- function(path, top, bottom) {
  list(bin = path[["k"]], hi = c(top, path[["t"]]), lo = c(path[["t"]], bottom))
}

# log(exact / frozen) weight of one path, given the pieces of its branch
.PathWeight <- function(set, pc, path) {
  if (.UseCpp()) {
    return(.PathWeightCpp(set, pc, path))
  }
  b <- set[["bounds"]]
  edges <- pc[["edges"]]
  cuts <- sort(unique(c(edges, path[["t"]])), decreasing = TRUE)
  hi <- cuts[-length(cuts)]
  lo <- cuts[-1]
  mid <- (hi + lo) / 2
  k <- path[["k"]][findInterval(-mid, -c(edges[[1]], path[["t"]]))]
  s <- pmin(findInterval(mid, b), length(set[["slices"]]))
  eBar <- vapply(seq_along(s), function(i) set[["slices"]][[s[i]]][["eBar"]][[k[i]]], 0)
  sum(2 * set[["lambda"]][k] * (.Hermite(set[["eg"]], hi, k) -
                                  .Hermite(set[["eg"]], lo, k) - eBar * (hi - lo)))
}

# log(exact / frozen) weight of the paths on branches `nodes`
.LogWeight <- function(tree, set, msg, paths, nodes) {
  if (.UseCpp()) {
    return(.LogWeightCpp(tree, set, msg, paths, nodes))
  }
  sum(vapply(nodes, function(v) .PathWeight(set, msg[["pieces"]][[v]], paths[[v]]), 0))
}

# Exact log augmented density (spec equation 3.4), start at the root or, with
# an origin, at the top of the stem
.LogAugmented <- function(tree, set, paths, bin) {
  if (.UseCpp()) {
    return(.LogAugmentedCpp(tree, set, paths, bin))
  }
  origin <- .Origin(tree)
  if (!.Compatible(tree, set[["condition"]], origin)) {
    return(-Inf)
  }
  nTip <- tree[["nTip"]]
  age <- tree[["age"]]
  nodes <- .PathNodes(tree)
  top <- .TopAge(tree, nodes)
  ints <- lapply(seq_along(nodes), function(i) {
    .Intervals(paths[[nodes[[i]]]], top[[i]], age[[nodes[[i]]]])
  })
  k <- unlist(lapply(ints, `[[`, "bin"))
  hi <- unlist(lapply(ints, `[[`, "hi"))
  lo <- unlist(lapply(ints, `[[`, "lo"))
  Q <- set[["Q"]]
  jumps <- do.call(rbind, lapply(nodes, function(v) {
    kk <- paths[[v]][["k"]]
    if (length(kk) > 1L) cbind(kk[-length(kk)], kk[-1]) else NULL
  }))
  lambda <- set[["lambda"]]
  pathPart <- sum((diag(Q)[k] - set[["sigma"]][k]) * (hi - lo) +
                    2 * lambda[k] * (.Hermite(set[["eg"]], hi, k) -
                                       .Hermite(set[["eg"]], lo, k))) +
    (if (is.null(jumps)) 0 else sum(log(Q[jumps])))
  internal <- tree[["postorder"]]
  bif <- internal[!tree[["sampled"]][internal - nTip] &
                    (internal != tree[["root"]] | origin)]
  sa <- internal[tree[["sampled"]][internal - nTip]]
  tips <- seq_len(nTip)
  tipPart <- sum(log(.TipValues(tree, set)[cbind(tips, bin[tips])][tree[["tipType"]] != "ancestor"]))
  start <- bin[[if (origin) nTip + tree[["nNode"]] + 1L else tree[["root"]]]]
  survive <- set[["survive"]]
  prior <- set[["prior"]]
  condition <- if (set[["joint"]]) {
    log(survive[[start]])
  } else {
    log(sum(prior * survive))
  }
  log(prior[[start]]) + pathPart + sum(log(lambda[bin[bif]])) +
    sum(log(set[["psi"]][bin[sa]] * (1 - set[["removal"]]))) + tipPart - condition
}
