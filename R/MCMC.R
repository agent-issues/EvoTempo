#' Default parameter priors
#'
#' Log-normal priors on the positive dCET parameters, as a named list of
#' functions returning log densities.
#'
#' @param meanlog,sdlog Named numeric vectors giving the parameters of each
#' log-normal prior; names select the parameters.
#' @returns `TempoPriors()` returns a named list of functions, each giving the
#' log prior density of one parameter.
#' @examples
#' priors <- TempoPriors()
#' priors$theta(0.1)
#' @family dCET inference
#' @export
TempoPriors <- function(meanlog = c(lambda = log(0.5), mu = log(0.2),
                                    psi = log(0.1), theta = log(0.1),
                                    s = log(0.8)),
                        sdlog = c(lambda = 0.5, mu = 0.5, psi = 0.5,
                                  theta = 0.75, s = 0.3)) {
  if (!setequal(names(meanlog), names(sdlog))) {
    stop("`meanlog` and `sdlog` must name the same parameters")
  }
  sdlog <- sdlog[names(meanlog)]
  out <- lapply(names(meanlog), function(p) {
    m <- meanlog[[p]]
    s <- sdlog[[p]]
    force(m)
    force(s)
    function(x) stats::dlnorm(x, m, s, log = TRUE)
  })
  names(out) <- names(meanlog)
  out
}

.LogPrior <- function(par, priors) {
  sum(vapply(names(priors), function(p) priors[[p]](par[[p]]), 0))
}

# Shared random-walk parameter update: multiplicative proposals
.ParameterSweep <- function(par, priors, step, current, Target) {
  accepted <- setNames(logical(length(priors)), names(priors))
  for (p in names(priors)) {
    proposal <- par
    factor <- exp(step[[p]] * (stats::runif(1) - 0.5))
    proposal[[p]] <- par[[p]] * factor
    lp <- priors[[p]](proposal[[p]])
    if (!is.finite(lp)) {
      next
    }
    candidate <- Target(proposal)
    logRatio <- candidate[["value"]] + lp - current[["value"]] -
      priors[[p]](par[[p]]) + log(factor)
    if (is.finite(logRatio) && log(stats::runif(1)) < logRatio) {
      par <- proposal
      current <- candidate
      accepted[[p]] <- TRUE
    }
  }
  list(par = par, current = current, accepted = accepted)
}

.Adapt <- function(step, rate, target = 0.3) {
  step * exp(rate - target)
}

#' MCMC on the marginal dCET tree likelihood
#'
#' Random-walk Metropolis on the parameters of a fixed time tree, using the
#' exact tempo-marginal likelihood of [TreeLikelihood()].
#' The reference against which [TempoMCMC()] is validated (rung 4).
#'
#' @inheritParams TreeLikelihood
#' @param start Named list giving starting values of all parameters; those
#' without a prior are held fixed.
#' @param priors Named list of functions giving log prior densities, as from
#' [TempoPriors()]; names select the free parameters.
#' @param nIter Integer giving the number of iterations after burn-in.
#' @param burnin Integer giving the number of burn-in iterations, during which
#' step sizes adapt.
#' @param thin Integer giving the interval between stored samples.
#' @inheritParams TreeLikelihood
#' @returns `MarginalMCMC()` returns a data frame of samples: one column per
#' free parameter, and `logLik`.
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' set.seed(1)
#' MarginalMCMC(tree, DCETModel(nBin = 11, zRange = c(-2, 2)),
#'              list(lambda = 0.5, mu = 0.2, theta = 0.1, s = 1),
#'              TempoPriors()[c("lambda", "mu")], nIter = 20, burnin = 10)
#' @family dCET inference
#' @export
MarginalMCMC <- function(tree, model, start, priors, nIter, burnin = 0,
                         thin = 1, character = NULL, nState = 2L) {
  tree <- TempoTree(tree)
  Target <- function(par) {
    list(value = TreeLikelihood(tree, par, model, character, nState))
  }
  .RunParameterChain(start, priors, nIter, burnin, thin, Target)
}

.RunParameterChain <- function(par, priors, nIter, burnin, thin, Target) {
  par <- as.list(par)
  current <- Target(par)
  step <- setNames(rep(1, length(priors)), names(priors))
  nKeep <- nIter %/% thin
  out <- matrix(NA_real_, nKeep, length(priors) + 1L,
                dimnames = list(NULL, c(names(priors), "logLik")))
  hits <- numeric(length(priors))
  for (i in seq_len(burnin + nIter)) {
    sweep <- .ParameterSweep(par, priors, step, current, Target)
    par <- sweep[["par"]]
    current <- sweep[["current"]]
    hits <- hits + sweep[["accepted"]]
    if (i <= burnin && i %% 20 == 0) {
      step <- .Adapt(step, hits / 20)
      hits <- 0 * hits
    }
    j <- i - burnin
    if (j > 0 && j %% thin == 0) {
      out[j / thin, ] <- c(unlist(par[names(priors)]), current[["value"]])
    }
  }
  as.data.frame(out)
}

#' Data-augmented dCET MCMC
#'
#' Samples parameters jointly with the tempo path on every branch of a time
#' tree, under the exact augmented density (spec section 3.4), optionally with
#' discrete characters scored on tempo-weighted branch lengths and with moves
#' on node ages, fossil ages and topology.
#'
#' Path moves propose from a frozen-extinction model, in which
#' \eqn{E(\tau)} is held at its value at the midpoint of each of `nSlice`
#' time slices: whole-tree forward-filtering backward-sampling (`tree`);
#' `node`, which redraws a node's bin and its incident paths from their frozen
#' conditional law; `tip`, which redraws a terminal branch and its tip bin;
#' and `branch`, which redraws one path between fixed end bins.
#' Within each slice, paths are drawn exactly from the endpoint-conditioned
#' frozen law (Hobolth & Stone 2009, direct sampling).
#' Every path move is accepted by Metropolis-Hastings on the ratio of exact to
#' frozen densities, times the character likelihood ratio, so the chain
#' targets the exact model.
#'
#' Tree moves (root age fixed): `nodeAge` redraws the age of a bifurcation
#' uniformly between its children and parent; `fossilAge` redraws a fossil's
#' age within its stratigraphic range (see [TempoTree()]); `slide` prunes a
#' subtree and regrafts it at a uniform height onto a uniformly chosen branch
#' crossing that height; `saJump` converts a fossil tip to a sampled ancestor
#' or back (Gavryushkina et al. 2014).
#' Each redraws affected node bins and paths from the frozen model, with
#' proposal densities in the Hastings ratio.
#' The target is the density of the unoriented labelled tree, i.e. the
#' oriented density of [TreeLikelihood()] times \eqn{2^{|V|}}, \eqn{|V|} the
#' number of bifurcations.
#'
#' Characters follow the Mk model on branch lengths
#' \eqn{\ell_e = c \int_e e^{\gamma x} dt}, with \eqn{c} the parameter `rate`
#' and \eqn{\gamma} from `model`; see [MkLikelihood()].
#'
#' @inheritParams MarginalMCMC
#' @inheritParams MkLikelihood
#' @param nSlice Integer giving the number of time slices of the proposal.
#' @param moves Named numeric vector giving the number of each move per
#' iteration (see Details); unnamed moves are not used.
#' @param characters Optional matrix of character states, as for
#' [MkLikelihood()].
#' @param updateParameters Logical; if `FALSE`, hold parameters fixed.
#' @returns `TempoMCMC()` returns a list with elements
#' - `samples`: data frame with one column per free parameter, `logPost`
#'   (log augmented density, plus character log-likelihood and log prior),
#'   `rootBin`, `meanX` (lineage-time-averaged log-tempo over the tree) and
#'   tree statistics: `nAncestor`, `treeLength`, `meanNodeAge` (of non-root
#'   bifurcations), `cherries`, `oldestFossil`, `youngestFossil`, `clock`
#'   (\eqn{\sum_e \int_e e^x dt}) and `inverseClock`
#'   (\eqn{\sum_e \int_e e^{-x} dt});
#' - `occupancy`: matrix of time spent in each bin (columns), per sample;
#' - `branchX`: matrix of the time-averaged log-tempo on each branch (for
#'   chains without tree moves);
#' - `acceptance`: acceptance rate of each move;
#' - `state`: the final state.
#' @references
#' Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014). Bayesian inference
#' of sampled ancestor trees for epidemiology and fossil calibration. PLoS
#' Computational Biology 10:e1003919. \doi{10.1371/journal.pcbi.1003919}
#'
#' Hobolth A, Stone EA (2009). Simulation from endpoint-conditioned,
#' continuous-time Markov chains on a finite state space, with applications to
#' molecular evolution. Annals of Applied Statistics 3:1204--1231.
#' \doi{10.1214/09-AOAS247}
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' set.seed(1)
#' run <- TempoMCMC(tree, DCETModel(nBin = 11, zRange = c(-2, 2)),
#'                  list(lambda = 0.5, mu = 0.2, theta = 0.1, s = 1),
#'                  TempoPriors()[c("lambda", "theta")], nIter = 5, nSlice = 4)
#' run$samples
#' @family dCET inference
#' @export
TempoMCMC <- function(tree, model, start, priors, nIter, burnin = 0, thin = 1,
                      nSlice = 20, moves = c(tree = 1), characters = NULL,
                      nState = 2L, updateParameters = TRUE) {
  tree <- TempoTree(tree)
  if (model[["start"]] != "root") {
    stop("TempoMCMC() supports `start = \"root\"` only")
  }
  moveNames <- c("tree", "node", "tip", "branch", "nodeAge", "fossilAge",
                 "slide", "saJump")
  if (!all(names(moves) %in% moveNames)) {
    stop("Unknown move: ", paste(setdiff(names(moves), moveNames), collapse = ", "))
  }
  moves <- replace(setNames(numeric(length(moveNames)), moveNames), names(moves),
                   moves)
  par <- as.list(start)
  CharLL <- function(tr, paths, p, set) {
    if (is.null(characters)) {
      return(0)
    }
    MkLikelihood(tr, characters,
                 .CharacterLengths(tr, set[["grid"]][["x"]], paths, p[["rate"]],
                                   model[["gamma"]]), nState)
  }
  set <- .AddSlices(tree, .Precompute(tree, par, model), nSlice)
  msg <- .Messages(tree, set)
  draw <- .SampleTree(tree, set, msg)
  st <- list(tree = tree, par = par, set = set, msg = msg,
             paths = draw[["paths"]], bin = draw[["bin"]])
  st[["charLL"]] <- CharLL(tree, st[["paths"]], par, set)
  Target <- function(p) {
    s <- .Precompute(st[["tree"]], p, model)
    list(value = .LogAugmented(st[["tree"]], s, st[["paths"]], st[["bin"]]) +
           CharLL(st[["tree"]], st[["paths"]], p, s), set = s)
  }
  step <- setNames(rep(1, length(priors)), names(priors))
  hits <- numeric(length(priors))
  pathNodes <- .PathNodes(tree)
  nKeep <- nIter %/% thin
  statNames <- names(.TreeStats(tree))
  samples <- matrix(NA_real_, nKeep, length(priors) + 3L + length(statNames),
                    dimnames = list(NULL, c(names(priors), "logPost", "rootBin",
                                            "meanX", statNames)))
  occupancy <- matrix(NA_real_, nKeep, model[["nBin"]])
  branchX <- matrix(NA_real_, nKeep, length(pathNodes))
  moveTry <- moveAcc <- setNames(numeric(length(moveNames) + 1L),
                                 c(moveNames, "parameters"))
  treeMoves <- c("nodeAge", "fossilAge", "slide", "saJump")
  for (i in seq_len(burnin + nIter)) {
    for (move in moveNames[moves > 0]) {
      for (rep in seq_len(moves[[move]])) {
        res <- if (move %in% treeMoves) {
          .TreeMove(move, st, CharLL)
        } else {
          .PathMove(move, st, CharLL)
        }
        st <- res[["state"]]
        moveTry[[move]] <- moveTry[[move]] + 1
        moveAcc[[move]] <- moveAcc[[move]] + res[["accepted"]]
      }
    }
    if (updateParameters && length(priors)) {
      current <- list(value = .LogAugmented(st[["tree"]], st[["set"]],
                                            st[["paths"]], st[["bin"]]) +
                        st[["charLL"]], set = st[["set"]])
      sweep <- .ParameterSweep(st[["par"]], priors, step, current, Target)
      moveTry[["parameters"]] <- moveTry[["parameters"]] + length(priors)
      moveAcc[["parameters"]] <- moveAcc[["parameters"]] + sum(sweep[["accepted"]])
      hits <- hits + sweep[["accepted"]]
      if (any(sweep[["accepted"]])) {
        st[["par"]] <- sweep[["par"]]
        st[["set"]] <- .AddSlices(st[["tree"]], sweep[["current"]][["set"]], nSlice)
        st[["msg"]] <- .Messages(st[["tree"]], st[["set"]])
        st[["charLL"]] <- CharLL(st[["tree"]], st[["paths"]], st[["par"]], st[["set"]])
      }
      if (i <= burnin && i %% 20 == 0) {
        step <- .Adapt(step, hits / 20)
        hits <- 0 * hits
      }
    }
    j <- i - burnin
    if (j > 0 && j %% thin == 0) {
      tr <- st[["tree"]]
      summ <- .PathSummary(tr, st, .PathNodes(tr))
      x <- st[["set"]][["grid"]][["x"]]
      samples[j / thin, ] <- c(
        unlist(st[["par"]][names(priors)]),
        .LogAugmented(tr, st[["set"]], st[["paths"]], st[["bin"]]) +
          st[["charLL"]] + .LogPrior(st[["par"]], priors),
        st[["bin"]][[tr[["root"]]]], summ[["meanX"]],
        .TreeStats(tr, .CharacterLengths(tr, x, st[["paths"]], 1, 1),
                   .CharacterLengths(tr, x, st[["paths"]], 1, -1)))
      occupancy[j / thin, ] <- summ[["occupancy"]]
      if (length(summ[["branchX"]]) == ncol(branchX)) {
        branchX[j / thin, ] <- summ[["branchX"]]
      }
    }
  }
  list(samples = as.data.frame(samples), occupancy = occupancy,
       branchX = branchX, acceptance = moveAcc / pmax(moveTry, 1),
       state = st)
}

# One tree move, accepted on the exact target
.TreeMove <- function(move, st, CharLL) {
  tree <- st[["tree"]]
  set <- st[["set"]]
  prop <- tryCatch(.ProposeTreeMove(move, tree, set, st),
                   infeasibleProposal = function(e) NULL)
  if (is.null(prop)) {
    return(list(state = st, accepted = FALSE))
  }
  newTree <- prop[["tree"]]
  newChar <- CharLL(newTree, prop[["paths"]], st[["par"]], set)
  nBif <- function(tr) sum(!tr[["sampled"]])
  logRatio <- .LogAugmented(newTree, set, prop[["paths"]], prop[["bin"]]) -
    .LogAugmented(tree, set, st[["paths"]], st[["bin"]]) +
    (nBif(newTree) - nBif(tree)) * log(2) + newChar - st[["charLL"]] +
    prop[["logHastings"]]
  if (is.finite(logRatio) && log(stats::runif(1)) < logRatio) {
    st[["tree"]] <- newTree
    st[["paths"]] <- prop[["paths"]]
    st[["bin"]] <- prop[["bin"]]
    st[["charLL"]] <- newChar
    st[["msg"]] <- .Messages(newTree, set)
    return(list(state = st, accepted = TRUE))
  }
  list(state = st, accepted = FALSE)
}

# Dwell per bin, and time-averaged log-tempo per branch and overall
.PathSummary <- function(tree, st, pathNodes = .PathNodes(tree)) {
  age <- tree[["age"]]
  x <- st[["set"]][["grid"]][["x"]]
  nBin <- length(x)
  occupancy <- numeric(nBin)
  branchX <- numeric(length(pathNodes))
  for (i in seq_along(pathNodes)) {
    v <- pathNodes[[i]]
    int <- .Intervals(st[["paths"]][[v]], age[[tree[["parent"]][[v]]]], age[[v]])
    dwell <- int[["hi"]] - int[["lo"]]
    grouped <- rowsum(dwell, int[["bin"]])
    idx <- as.integer(rownames(grouped))
    occupancy[idx] <- occupancy[idx] + grouped[, 1]
    branchX[i] <- sum(dwell * x[int[["bin"]]]) / sum(dwell)
  }
  list(occupancy = occupancy, branchX = branchX,
       meanX = sum(occupancy * x) / sum(occupancy))
}

# Transfer of a branch: row from a fixed top bin, or column to a fixed end bin
.TransferRow <- function(pc, top) {
  r <- replace(numeric(nrow(pc[["M"]][[1]])), top, 1)
  for (M in pc[["M"]]) {
    r <- as.vector(r %*% M)
    r <- r / max(r)
  }
  r
}

.TransferCol <- function(pc, end) {
  r <- replace(numeric(nrow(pc[["M"]][[1]])), end, 1)
  for (M in rev(pc[["M"]])) {
    r <- as.vector(M %*% r)
    r <- r / max(r)
  }
  r
}

.PathKids <- function(tree, v) {
  nTip <- tree[["nTip"]]
  kids <- tree[["kids"]][v - nTip, ]
  kids[!(kids <= nTip & tree[["tipType"]][pmin(kids, nTip)] == "ancestor")]
}

# One path move; returns the new state and whether it was accepted
.PathMove <- function(move, st, CharLL) {
  tree <- st[["tree"]]
  set <- st[["set"]]
  msg <- st[["msg"]]
  nTip <- tree[["nTip"]]
  parent <- tree[["parent"]]
  bin <- st[["bin"]]
  paths <- st[["paths"]]
  nodes <- .PathNodes(tree)
  proposal <- tryCatch({
    switch(
      move,
      tree = {
        draw <- .SampleTree(tree, set, msg)
        changed <- nodes
        newPaths <- draw[["paths"]]
        newBin <- draw[["bin"]]
      },
      tip = {
        tips <- nodes[nodes <= nTip]
        v <- tips[[.Sample1(rep(1, length(tips)))]]
        draw <- .SampleClade(tree, set, msg, paths, bin, v, bin[[parent[[v]]]])
        changed <- v
        newPaths <- draw[["paths"]]
        newBin <- draw[["bin"]]
      },
      branch = {
        v <- nodes[[.Sample1(rep(1, length(nodes)))]]
        changed <- v
        newPaths <- paths
        newPaths[[v]] <- .SampleBranch(set, msg[["pieces"]][[v]], NULL,
                                       bin[[parent[[v]]]], bin[[v]])
        newBin <- bin
      },
      node = {
        internal <- tree[["postorder"]]
        v <- internal[[.Sample1(rep(1, length(internal)))]]
        kids <- .PathKids(tree, v)
        w <- if (v == tree[["root"]]) {
          set[["prior"]] / (if (set[["joint"]]) set[["survive"]] else 1)
        } else if (tree[["sampled"]][v - nTip]) {
          set[["psi"]] * .TransferRow(msg[["pieces"]][[v]], bin[[parent[[v]]]])
        } else {
          set[["lambda"]] * .TransferRow(msg[["pieces"]][[v]], bin[[parent[[v]]]])
        }
        for (kid in kids) {
          w <- w * .TransferCol(msg[["pieces"]][[kid]], bin[[kid]])
        }
        newBin <- bin
        newBin[v] <- .Sample1(w)
        sa <- tree[["kids"]][v - nTip, ]
        newBin[sa[sa <= nTip & tree[["tipType"]][pmin(sa, nTip)] == "ancestor"]] <- newBin[v]
        changed <- c(if (v != tree[["root"]]) v, kids)
        newPaths <- paths
        for (u in changed) {
          newPaths[[u]] <- .SampleBranch(set, msg[["pieces"]][[u]], NULL,
                                         newBin[[parent[[u]]]], newBin[[u]])
        }
      }
    )
    TRUE
  }, infeasibleProposal = function(e) NULL)
  if (is.null(proposal)) {
    return(list(state = st, accepted = FALSE))
  }
  newChar <- CharLL(tree, newPaths, st[["par"]], set)
  logRatio <- .LogWeight(tree, set, msg, newPaths, changed) -
    .LogWeight(tree, set, msg, paths, changed) + newChar - st[["charLL"]]
  if (log(stats::runif(1)) < logRatio) {
    st[["paths"]] <- newPaths
    st[["bin"]] <- newBin
    st[["charLL"]] <- newChar
    return(list(state = st, accepted = TRUE))
  }
  list(state = st, accepted = FALSE)
}

#' Stochastic mapping of tempo paths
#'
#' Draws tempo paths on every branch of a fixed time tree from the frozen
#' tree-conditional law used as a proposal by [TempoMCMC()], with importance
#' weights that correct them to the exact law.
#'
#' @inheritParams TempoMCMC
#' @param par Named list giving all parameters.
#' @param n Integer giving the number of draws.
#' @returns `SampleTempoPaths()` returns a list with elements `occupancy`,
#' `branchX` and `meanX` (as in [TempoMCMC()], one row or element per draw),
#' `logWeight` (log importance weight of each draw), and `logLik`, the
#' importance-sampling estimate of the log marginal likelihood.
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' set.seed(1)
#' maps <- SampleTempoPaths(tree, list(lambda = 0.5, mu = 0.2, theta = 0.1, s = 1),
#'                          DCETModel(nBin = 11, zRange = c(-2, 2)), n = 5,
#'                          nSlice = 4)
#' maps$logLik
#' @family dCET inference
#' @export
SampleTempoPaths <- function(tree, par, model, n, nSlice = 20) {
  tree <- TempoTree(tree)
  set <- .AddSlices(tree, .Precompute(tree, as.list(par), model), nSlice)
  msg <- .Messages(tree, set)
  pathNodes <- .PathNodes(tree)
  draws <- lapply(seq_len(n), function(i) {
    d <- .SampleTree(tree, set, msg)
    st <- list(set = set, paths = d[["paths"]], bin = d[["bin"]])
    c(.PathSummary(tree, st, pathNodes),
      logWeight = .LogWeight(tree, set, msg, d[["paths"]], pathNodes))
  })
  logWeight <- vapply(draws, `[[`, 0, "logWeight")
  top <- max(logWeight)
  list(occupancy = do.call(rbind, lapply(draws, `[[`, "occupancy")),
       branchX = do.call(rbind, lapply(draws, `[[`, "branchX")),
       meanX = vapply(draws, `[[`, 0, "meanX"),
       logWeight = logWeight,
       logLik = .FrozenLikelihood(tree, set, msg, model[["joint"]]) + top +
         log(mean(exp(logWeight - top))))
}
