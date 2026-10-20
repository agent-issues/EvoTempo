# Mixing study for TempoMCMC() with many characters (issue #19): shared code.
#
# A stand-alone driver that reproduces the move schedule of TempoMCMC() on a
# fixed time tree (path moves `tree`, `node`, `tip`, `branch`, `joint`, then
# the random-walk parameter sweep, adapting step sizes every 20 burn-in
# iterations), calling the engine's internals through `:::`, so that candidate
# moves can be added WITHOUT touching the engine.  With no extra moves the
# driver consumes random numbers exactly as TempoMCMC() does and returns the
# same chain (`CheckDriver()`).
#
# Candidate moves (specification in dev/mixing/README.md):
# - ShiftMove(): every bin moves by m = +/-1 and the multipliers of the
#   tempo-scaled rates are divided by exp(m h s * exponent), so that the
#   character lengths are exactly unchanged;
# - AnchoredJoint(): a parameter is proposed together with fresh frozen
#   bridges of every path between "anchor" bins (node bins, and the path's
#   bins at a subset of slice boundaries), accepted on the exact target with
#   the anchored bridge densities as proposal densities.
#
# Sourced by the scripts in dev/mixing/ from the package root.

if (nzchar(Sys.getenv("EVOTEMPO_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = FALSE))
} else {
  suppressMessages(library(EvoTempo))
}
ns <- asNamespace("EvoTempo")
for (fn in c(".Setup", ".Sample1", ".Precompute", ".AddSlices", ".Messages",
             ".SampleTree", ".SampleBranch", ".LogAugmented", ".LogWeight",
             ".PathWeight", ".FrozenLikelihood", ".JointSweep", ".PathMove",
             ".ParameterSweep", ".Adapt", ".SameSet", ".PathNodes",
             ".CharacterLengths", ".PathSummary", ".LogPathExact",
             ".LogTransferRow", ".LogBridges", ".Intervals", ".LogPrior",
             ".DrawBridges")) {
  assign(fn, get(fn, ns))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Model and priors (as rung 7, #7) ------------------------------------------
grid0 <- list(nBin = 16, zRange = c(-3, 3))
rho <- 0.5
nSlice <- 16L
psiLevels <- c(low = log(0.02), high = log(0.2))
sizeLevels <- list(small = list(crown = 8, extant = c(15, 30)),
                   large = list(crown = 11, extant = c(40, 80)))
Meanlog <- function(psiLevel = "high") {
  c(lambda = log(0.5), mu = log(0.2), psi = psiLevels[[psiLevel]],
    theta = log(0.2), s = log(0.6), rate = log(0.1))
}
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.5, rate = 0.4)
zTrunc <- 4
logMass <- log(stats::pnorm(zTrunc) - stats::pnorm(-zTrunc))
TruncPriors <- function(meanlog, sdlog) {
  out <- lapply(names(meanlog), function(p) {
    m <- meanlog[[p]]
    s <- sdlog[[p]]
    force(m)
    force(s)
    function(x) {
      if (!(x > 0) || abs(log(x) - m) > zTrunc * s) {
        return(-Inf)
      }
      stats::dlnorm(x, m, s, log = TRUE) - logMass
    }
  })
  names(out) <- names(meanlog)
  out
}
PriorDraw <- function(meanlog) {
  z <- stats::rnorm(length(meanlog))
  while (any(out <- abs(z) > zTrunc)) {
    z[out] <- stats::rnorm(sum(out))
  }
  as.list(setNames(exp(meanlog + sdlog[names(meanlog)] * z), names(meanlog)))
}
Model <- function(beta = 1, gamma = 1) {
  DCETModel(nBin = grid0$nBin, zRange = grid0$zRange, rho = rho,
            beta = beta, gamma = gamma)
}
baseMoves <- c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1)

# ---- Simulation (as rung 7) ------------------------------------------------------
Join <- function(h1, h2, x0) {
  n1 <- nrow(h1$edges)
  e1 <- h1$edges
  e2 <- h2$edges
  e1$parent <- ifelse(is.na(e1$parent), 1L, e1$parent + 1L)
  e2$parent <- ifelse(is.na(e2$parent), 1L, e2$parent + 1L + n1)
  f1 <- h1$fossils
  f2 <- h2$fossils
  f1$edge <- f1$edge + 1L
  f2$edge <- f2$edge + 1L + n1
  structure(list(edges = rbind(data.frame(parent = NA_integer_, start = 0, end = 0,
                                          type = "speciation", sampled = NA), e1, e2),
                 paths = c(list(list(t = 0, x = x0)), h1$paths, h2$paths),
                 fossils = rbind(f1, f2), time = h1$time),
            class = "TempoHistory", kind = "discrete")
}

# Survival imposed by resimulating at fixed parameters (AGENTS.md); rejection
# on the number of extant tips redraws the parameters.  nCharMax characters
# are simulated on the true tempo-weighted lengths; analyses use the first n,
# so the character sets are nested.
Simulate <- function(model, meanlog, size, nCharMax = 300L) {
  repeat {
    truth <- PriorDraw(meanlog)
    set <- .Setup(truth, model)
    repeat {
      k <- .Sample1(set$prior)
      h <- tryCatch(lapply(1:2, function(i) {
        SimulateDCET(size$crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                     start = k, rho = rho, maxEdges = 8000)
      }), error = function(e) NULL)
      if (is.null(h) ||
          all(vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0) > 0)) {
        break
      }
    }
    if (is.null(h)) next
    history <- Join(h[[1]], h[[2]], set$grid$x[k])
    tree <- ReconstructedTree(history, gamma = model$gamma, rate = truth$rate)
    nExtant <- sum(attr(tree, "tipType") == "extant")
    if (nExtant >= size$extant[1] && nExtant <= size$extant[2]) break
  }
  chars <- SimulateCharacters(tree, nCharMax, 2)
  tree$root.edge <- NULL
  tt <- TempoTree(tree, present = size$crown)
  list(truth = truth, rootBin = k, tt = tt, chars = chars,
       nExtant = sum(tt$tipType == "extant"), nFossil = sum(tt$tipType == "fossil"),
       nAncestor = sum(tt$tipType == "ancestor"))
}

# ---- Utilities ------------------------------------------------------------------
# Geyer initial monotone sequence ESS (as rungs 6-7)
Ess <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4 || stats::var(x) == 0) {
    return(NA_real_)
  }
  ac <- stats::acf(x, lag.max = n - 1, plot = FALSE, demean = TRUE)$acf[, 1, 1]
  nPair <- length(ac) %/% 2
  pairSum <- ac[2 * seq_len(nPair) - 1] + ac[2 * seq_len(nPair)]
  stop1 <- which(pairSum <= 0)[1]
  if (!is.na(stop1)) {
    pairSum <- pairSum[seq_len(stop1 - 1)]
  }
  tau <- if (length(pairSum)) -1 + 2 * sum(cummin(pairSum)) else 1
  min(n, n / max(tau, 1e-8))
}

CharLLFun <- function(characters, model, nState = 2L) {
  force(characters)
  function(tr, paths, p, set) {
    if (is.null(characters)) {
      return(0)
    }
    MkLikelihood(tr, characters,
                 .CharacterLengths(tr, set[["grid"]][["x"]], paths, p[["rate"]],
                                   model[["gamma"]]), nState)
  }
}

# Unnormalised log target of a state (augmented density x characters x prior)
LogTarget <- function(st, priors) {
  .LogAugmented(st$tree, st$set, st$paths, st$bin) + st$charLL +
    .LogPrior(st$par, priors)
}

NJumps <- function(st) {
  sum(vapply(.PathNodes(st$tree), function(v) length(st$paths[[v]]$t), 0))
}

# Per-branch clock int_e e^x dt (gamma = 1, rate = 1), path nodes only
BranchClock <- function(st, nodes = .PathNodes(st$tree)) {
  .CharacterLengths(st$tree, st$set$grid$x, st$paths, 1, 1)[nodes]
}

# ---- Candidate move 1: global bin shift -----------------------------------------
# Every bin (nodes and path pieces) moves by m in {-1, +1}, so x moves by
# d = m h s everywhere.  Each compensated parameter p is multiplied by
# exp(-e_p d), e_p its tempo exponent: rate (gamma), lambda and mu (beta), psi
# (delta), theta (alpha).  With rate compensated the character lengths are
# unchanged exactly.  T_m is a bijection with inverse T_-m; m is drawn with
# probability 1/2 each way; Jacobian in natural scale prod_p p' / p.  A shift
# that leaves the grid is rejected (its reverse would leave it too).
ShiftMove <- function(st, model, priors, comp, CharLL) {
  K <- model$nBin
  m <- if (stats::runif(1) < 0.5) -1L else 1L
  nodes <- .PathNodes(st$tree)
  allBins <- c(st$bin[st$bin > 0], unlist(lapply(st$paths[nodes], `[[`, "k")))
  if (any(allBins + m < 1L | allBins + m > K)) {
    return(list(state = st, accepted = FALSE))
  }
  d <- m * st$set$grid$h * st$par$s
  expo <- c(rate = model$gamma, lambda = model$beta, mu = model$beta,
            psi = model$delta, theta = model$alpha)
  comp <- intersect(comp, names(priors))
  par <- st$par
  logJac <- 0
  lp <- 0
  for (p in comp) {
    par[[p]] <- st$par[[p]] * exp(-expo[[p]] * d)
    logJac <- logJac - expo[[p]] * d
    lp <- lp + priors[[p]](par[[p]]) - priors[[p]](st$par[[p]])
  }
  if (!is.finite(lp)) {
    return(list(state = st, accepted = FALSE))
  }
  bin <- st$bin
  bin[bin > 0] <- bin[bin > 0] + m
  paths <- st$paths
  for (v in nodes) {
    paths[[v]]$k <- paths[[v]]$k + m
  }
  set <- tryCatch(.Precompute(st$tree, par, model), error = function(e) NULL)
  if (is.null(set)) {
    return(list(state = st, accepted = FALSE))
  }
  charLL <- CharLL(st$tree, paths, par, set)
  logRatio <- .LogAugmented(st$tree, set, paths, bin) + charLL -
    .LogAugmented(st$tree, st$set, st$paths, st$bin) - st$charLL + lp + logJac
  if (is.finite(logRatio) && log(stats::runif(1)) < logRatio) {
    st$par <- par
    st$set <- .AddSlices(st$tree, set, length(st$set$slices))
    st$msg <- .Messages(st$tree, st$set)
    st$bin <- bin
    st$paths <- paths
    st$charLL <- charLL
    return(list(state = st, accepted = TRUE))
  }
  list(state = st, accepted = FALSE)
}

# ---- Candidate move 2: anchored joint parameter-path move ---------------------
# Anchors on branch v: its end bins, and the path's bin at each interior piece
# boundary (a slice bound b_i, i = 0..nSlice) with i %% every == phase.
# every = Inf keeps the node bins only.  The anchor times depend only on the
# tree and the (state-independent) phase, and the proposal leaves the anchor
# bins unchanged, so q(P' | anchors(P)) and q(P | anchors(P')) share them.

# Bin of a path at age `a` strictly inside its branch (t: decreasing jump ages)
BinAt <- function(path, a) {
  path$k[[sum(path$t > a) + 1L]]
}

# Anchored segments of the pieces of one branch: list of piece-index vectors
Segments <- function(pc, set, every, phase) {
  n <- length(pc$M)
  if (n == 1L || !is.finite(every)) {
    return(list(seq_len(n)))
  }
  b <- set$bounds
  idx <- vapply(pc$edges[2:n], function(e) which.min(abs(b - e)), 0L) - 1L
  cut <- which(idx %% every == phase)   # piece j ends on an anchor
  unname(split(seq_len(n), cumsum(c(0L, as.integer(seq_len(n - 1L) %in% cut)))))
}

SubPieces <- function(pc, js) {
  list(edges = pc$edges[c(js, max(js) + 1L)], slice = pc$slice[js], M = pc$M[js])
}

SegmentEnds <- function(pc, path, js) {
  lastJ <- js[[length(js)]]
  c(if (js[[1]] == 1L) path$k[[1]] else BinAt(path, pc$edges[[js[[1]]]]),
    if (lastJ == length(pc$M)) path$k[[length(path$k)]] else
      BinAt(path, pc$edges[[lastJ + 1L]]))
}

# Log density of all paths under anchored frozen bridges:
# frozen path density / prod over segments of the frozen transfer between anchors
LogAnchored <- function(tree, set, msg, paths, every, phase) {
  out <- 0
  for (v in .PathNodes(tree)) {
    pc <- msg$pieces[[v]]
    path <- paths[[v]]
    out <- out + .LogPathExact(tree, set, v, path) - .PathWeight(set, pc, path)
    for (js in Segments(pc, set, every, phase)) {
      ab <- SegmentEnds(pc, path, js)
      out <- out - .LogTransferRow(SubPieces(pc, js), ab[[1]])[[ab[[2]]]]
    }
  }
  out
}

# Draw anchored frozen bridges under set/msg, anchors read from `paths`
DrawAnchored <- function(tree, set, msg, paths, every, phase) {
  new <- paths
  for (v in .PathNodes(tree)) {
    pc <- msg$pieces[[v]]
    path <- paths[[v]]
    tAll <- numeric(0)
    kAll <- path$k[[1]]
    for (js in Segments(pc, set, every, phase)) {
      ab <- SegmentEnds(pc, path, js)
      seg <- .SampleBranch(set, SubPieces(pc, js), NULL, ab[[1]], ab[[2]])
      tAll <- c(tAll, seg$t)
      kAll <- c(kAll, seg$k[-1])
    }
    new[[v]] <- list(t = tAll, k = kAll)
  }
  new
}

# One anchored joint update of each parameter in `which`.  Target ratio x
# q(P | anchors, Theta) / q(P' | anchors, Theta') x prior ratio x factor.
AnchoredJoint <- function(st, model, priors, step, which, every, CharLL) {
  tree <- st$tree
  nSl <- length(st$set$slices)
  accepted <- setNames(logical(length(which)), which)
  for (p in which) {
    phase <- if (is.finite(every)) sample.int(every, 1L) - 1L else 0L
    proposal <- st$par
    factor <- exp(step[[p]] * (stats::runif(1) - 0.5))
    proposal[[p]] <- proposal[[p]] * factor
    lp <- priors[[p]](proposal[[p]])
    if (!is.finite(lp)) next
    if (.SameSet(st$set, tree, proposal)) {
      set <- st$set
      msg <- st$msg
    } else {
      set <- tryCatch(.AddSlices(tree, .Precompute(tree, proposal, model), nSl),
                      error = function(e) NULL)
      if (is.null(set)) next
      msg <- .Messages(tree, set)
    }
    paths <- tryCatch(DrawAnchored(tree, set, msg, st$paths, every, phase),
                      infeasibleProposal = function(e) NULL)
    if (is.null(paths)) next
    logQnew <- LogAnchored(tree, set, msg, paths, every, phase)
    logQold <- LogAnchored(tree, st$set, st$msg, st$paths, every, phase)
    charLL <- CharLL(tree, paths, proposal, set)
    logRatio <- .LogAugmented(tree, set, paths, st$bin) + charLL + lp - logQnew -
      (.LogAugmented(tree, st$set, st$paths, st$bin) + st$charLL +
         priors[[p]](st$par[[p]]) - logQold) + log(factor)
    if (is.finite(logRatio) && log(stats::runif(1)) < logRatio) {
      st$par <- proposal
      st$set <- set
      st$msg <- msg
      st$paths <- paths
      st$charLL <- charLL
      accepted[[p]] <- TRUE
    }
  }
  list(state = st, accepted = accepted)
}

# ---- Driver ------------------------------------------------------------------------
# `extra`: list of extra moves, run after the TempoMCMC() path/joint schedule
# and before the parameter sweep.  Each is list(name, n, Fun, which, adapt,
# target, parameter), Fun(st, ctx) -> list(state, accepted), with ctx holding
# model, priors, CharLL and the move's own adaptive step vector.
# Stored per draw: parameters, logPost, charLL, rootBin, meanX, clock, nJumps;
# per branch: log clock and time-averaged log-tempo (branchX).
Driver <- function(tt, model, start, priors, nIter, burnin = 0, thin = 1,
                   nSlice = 16L, moves = baseMoves, characters = NULL,
                   updateParameters = TRUE, extra = list(), sweep = TRUE,
                   initial = NULL) {
  old <- options(EvoTempo.engine = TempoEngine())
  on.exit(options(old))
  tree <- TempoTree(tt)
  moveNames <- c("tree", "node", "tip", "branch", "joint")
  moves <- replace(setNames(numeric(length(moveNames)), moveNames), names(moves), moves)
  CharLL <- CharLLFun(characters, model)
  if (is.null(initial)) {
    par <- as.list(start)
    set <- .AddSlices(tree, .Precompute(tree, par, model), nSlice)
    msg <- .Messages(tree, set)
    draw <- .SampleTree(tree, set, msg)
    st <- list(tree = tree, par = par, set = set, msg = msg,
               paths = draw[["paths"]], bin = draw[["bin"]])
  } else {
    st <- initial[c("tree", "par", "set", "msg", "paths", "bin")]
  }
  st[["charLL"]] <- CharLL(tree, st[["paths"]], st[["par"]], st[["set"]])
  Target <- function(p, current) {
    s <- if (.SameSet(current[["set"]], st[["tree"]], p)) {
      current[["set"]]
    } else {
      .Precompute(st[["tree"]], p, model)
    }
    list(value = .LogAugmented(st[["tree"]], s, st[["paths"]], st[["bin"]]) +
           CharLL(st[["tree"]], st[["paths"]], p, s), set = s)
  }
  step <- setNames(rep(1, length(priors)), names(priors))
  hits <- numeric(length(priors))
  jointStep <- step
  jointHits <- hits
  extraNames <- vapply(extra, `[[`, "", "name")
  extraStep <- lapply(extra, function(e) e$step0 %||% step)
  extraHits <- lapply(extra, function(e) setNames(numeric(length(priors)), names(priors)))
  pathNodes <- .PathNodes(tree)
  nKeep <- nIter %/% thin
  cols <- c(names(priors), "logPost", "charLL", "rootBin", "meanX", "clock", "nJumps")
  samples <- matrix(NA_real_, nKeep, length(cols), dimnames = list(NULL, cols))
  logClock <- branchX <- matrix(NA_real_, nKeep, length(pathNodes))
  allNames <- c(moveNames, extraNames, "parameters")
  moveTry <- moveAcc <- moveSec <- setNames(numeric(length(allNames)), allNames)
  Tick <- function() proc.time()[["elapsed"]]
  for (i in seq_len(burnin + nIter)) {
    for (move in moveNames[moves > 0]) {
      for (rep in seq_len(moves[[move]])) {
        t0 <- Tick()
        if (move == "joint") {
          if (!updateParameters || !length(priors)) next
          res <- .JointSweep(st, priors, jointStep, model, nSlice, CharLL)
          jointHits <- jointHits + res[["accepted"]]
          moveTry[[move]] <- moveTry[[move]] + length(priors)
          moveAcc[[move]] <- moveAcc[[move]] + sum(res[["accepted"]])
          st <- res[["state"]]
        } else {
          res <- .PathMove(move, st, CharLL)
          st <- res[["state"]]
          moveTry[[move]] <- moveTry[[move]] + 1
          moveAcc[[move]] <- moveAcc[[move]] + res[["accepted"]]
        }
        moveSec[[move]] <- moveSec[[move]] + Tick() - t0
      }
    }
    for (e in seq_along(extra)) {
      ex <- extra[[e]]
      if (!updateParameters && isTRUE(ex$parameter)) next
      for (rep in seq_len(ex$n)) {
        t0 <- Tick()
        res <- ex$Fun(st, list(model = model, priors = priors, CharLL = CharLL,
                               step = extraStep[[e]]))
        st <- res$state
        acc <- res$accepted
        if (!is.null(names(acc))) {
          extraHits[[e]][names(acc)] <- extraHits[[e]][names(acc)] + acc
        }
        moveTry[[ex$name]] <- moveTry[[ex$name]] + length(acc)
        moveAcc[[ex$name]] <- moveAcc[[ex$name]] + sum(acc)
        moveSec[[ex$name]] <- moveSec[[ex$name]] + Tick() - t0
      }
    }
    if (sweep && updateParameters && length(priors)) {
      t0 <- Tick()
      current <- list(value = .LogAugmented(st[["tree"]], st[["set"]],
                                            st[["paths"]], st[["bin"]]) +
                        st[["charLL"]], set = st[["set"]])
      sw <- .ParameterSweep(st[["par"]], priors, step, current, Target)
      moveTry[["parameters"]] <- moveTry[["parameters"]] + length(priors)
      moveAcc[["parameters"]] <- moveAcc[["parameters"]] + sum(sw[["accepted"]])
      hits <- hits + sw[["accepted"]]
      if (any(sw[["accepted"]])) {
        st[["par"]] <- sw[["par"]]
        if (!.SameSet(st[["set"]], st[["tree"]], st[["par"]])) {
          st[["set"]] <- .AddSlices(st[["tree"]], sw[["current"]][["set"]], nSlice)
          st[["msg"]] <- .Messages(st[["tree"]], st[["set"]])
        }
        st[["charLL"]] <- CharLL(st[["tree"]], st[["paths"]], st[["par"]], st[["set"]])
      }
      moveSec[["parameters"]] <- moveSec[["parameters"]] + Tick() - t0
    }
    if (updateParameters && length(priors) && i <= burnin && i %% 20 == 0) {
      if (sweep) {
        step <- .Adapt(step, hits / 20)
        hits <- 0 * hits
      }
      if (moves[["joint"]] > 0) {
        jointStep <- .Adapt(jointStep, jointHits / (20 * moves[["joint"]]))
        jointHits <- 0 * jointHits
      }
      for (e in seq_along(extra)) {
        w <- extra[[e]]$which
        if (isTRUE(extra[[e]]$adapt) && length(w)) {
          extraStep[[e]][w] <- .Adapt(extraStep[[e]][w],
                                      extraHits[[e]][w] / (20 * extra[[e]]$n),
                                      target = extra[[e]]$target %||% 0.3)
          extraHits[[e]][w] <- 0
        }
      }
    }
    j <- i - burnin
    if (j > 0 && j %% thin == 0) {
      summ <- .PathSummary(st$tree, st, pathNodes)
      clk <- BranchClock(st, pathNodes)
      samples[j / thin, ] <- c(unlist(st$par[names(priors)]),
                               LogTarget(st, priors), st$charLL,
                               st$bin[[st$tree$root]], summ$meanX, sum(clk), NJumps(st))
      logClock[j / thin, ] <- log(clk)
      branchX[j / thin, ] <- summ$branchX
    }
  }
  list(samples = as.data.frame(samples), logClock = logClock, branchX = branchX,
       acceptance = moveAcc / pmax(moveTry, 1), tries = moveTry, seconds = moveSec,
       state = st, step = step, jointStep = jointStep, extraStep = extraStep)
}

# Extra-move constructors for Driver()
ShiftExtra <- function(comp, n = 1) {
  list(name = paste0("shift[", paste(comp, collapse = ","), "]"), n = n,
       parameter = TRUE, which = character(0),
       Fun = function(st, ctx) ShiftMove(st, ctx$model, ctx$priors, comp, ctx$CharLL))
}
AnchoredExtra <- function(which, every, n = 1, target = 0.3) {
  list(name = paste0("anchored", every, "[", paste(which, collapse = ","), "]"), n = n,
       parameter = TRUE, which = which, adapt = TRUE, target = target,
       Fun = function(st, ctx) {
         AnchoredJoint(st, ctx$model, ctx$priors, ctx$step, which, every, ctx$CharLL)
       })
}

# Driver() reproduces TempoMCMC() exactly when no extra moves are given
CheckDriver <- function(tt, model, start, priors, chars, nIter = 30, burnin = 20) {
  set.seed(1)
  a <- TempoMCMC(tt, model, start, priors, nIter = nIter, burnin = burnin,
                 nSlice = nSlice, moves = baseMoves, characters = chars)
  set.seed(1)
  b <- Driver(tt, model, start, priors, nIter = nIter, burnin = burnin,
              nSlice = nSlice, moves = baseMoves, characters = chars)
  isTRUE(all.equal(as.matrix(a$samples[names(priors)]),
                   as.matrix(b$samples[names(priors)]), tolerance = 0)) &&
    isTRUE(all.equal(a$samples$logPost, b$samples$logPost))
}

# ESS summary of a Driver() run
Summarise <- function(run) {
  ess <- vapply(run$samples, Ess, 0)
  bc <- apply(run$logClock, 2, Ess)
  bx <- apply(run$branchX, 2, Ess)
  list(ess = ess, essBranchClock = bc, essBranchX = bx,
       seconds = sum(run$seconds), acceptance = run$acceptance,
       moveSeconds = run$seconds, n = nrow(run$samples))
}
