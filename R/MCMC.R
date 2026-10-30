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
.ParameterSweep <- function(par, priors, step, current, Target,
                            Resolve = function(accepted) NULL) {
  accepted <- setNames(logical(length(priors)), names(priors))
  for (p in names(priors)) {
    proposal <- par
    factor <- exp(step[[p]] * (stats::runif(1) - 0.5))
    proposal[[p]] <- par[[p]] * factor
    lp <- priors[[p]](proposal[[p]])
    if (!is.finite(lp)) {
      next
    }
    candidate <- tryCatch(
      Target(proposal, current),
      infeasibleProposal = function(e) list(value = -Inf)
    )
    logRatio <- candidate[["value"]] + lp - current[["value"]] -
      priors[[p]](par[[p]]) + log(factor)
    if (is.finite(logRatio) && log(stats::runif(1)) < logRatio) {
      par <- proposal
      current <- candidate
      accepted[[p]] <- TRUE
    }
    Resolve(accepted[[p]])
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
#' @inheritParams TempoEngine
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
                         thin = 1, character = NULL, nState = 2L,
                         engine = getOption("EvoTempo.engine", "cpp")) {
  old <- options(EvoTempo.engine = TempoEngine(engine))
  on.exit(options(old))
  tree <- TempoTree(tree)
  Target <- function(par, current) {
    list(value = .SolverInfeasible(
      TreeLikelihood(tree, par, model, character, nState)
    ))
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
#' `joint` proposes each free parameter in turn together with a fresh FFBS
#' draw of every path under the proposed value; its acceptance ratio uses the
#' frozen likelihood and the exact-to-frozen weight, so it mixes like MCMC on
#' the marginal likelihood (at the cost of one FFBS draw per parameter).
#' Without it, parameters such as `theta` that the paths determine closely
#' mix slowly, at least without characters.
#'
#' `shift` moves every bin (nodes, tips and every piece of every path) by
#' \eqn{m = \pm 1} with probability 1/2 each, keeping the jump times, so
#' that log-tempo moves by \eqn{d = m h s} everywhere (\eqn{h} the bin
#' spacing in units of `s`).
#' Each free parameter that tempo scales is multiplied by
#' \eqn{e^{-e_p d}}, \eqn{e_p} its exponent: `rate` by \eqn{e^{-\gamma d}},
#' `lambda` and `mu` by \eqn{e^{-\beta d}}, `psi` by \eqn{e^{-\delta d}}
#' and `theta` by \eqn{e^{-\alpha d}}.
#' The map is a bijection whose inverse is the shift by \eqn{-m}; the
#' Hastings ratio is its Jacobian, \eqn{\prod_p p'/p}, and a shift that
#' would leave the grid is rejected.
#' With `rate` free, the tempo-weighted lengths, and so the character
#' likelihood, are unchanged exactly, so the move travels along the ridge
#' between `rate` and the overall tempo level that the characters leave
#' (spec section 9.1; issue #19).
#'
#' With characters, the whole-tree draws of `tree` and `joint` ignore the
#' characters and are almost never accepted, so the default schedule
#' ([TempoMoves()]) omits them and instead runs local path moves in
#' proportion to the size of the tree.
#' A fractional count \eqn{n + f} runs \eqn{n} moves and, with probability
#' \eqn{f}, one more.
#' At the end of burn-in, `tree` and `joint` are dropped from the sampling
#' phase if their acceptance during burn-in was below `minAcceptance`; move
#' frequencies change only during burn-in, so the sampling phase is a fixed
#' mixture of kernels.
#'
#' Tree moves: `nodeAge` redraws the age of a bifurcation
#' uniformly between its children and parent (with `start = "origin"`, the
#' root's too, between its children and the origin); `fossilAge` redraws a fossil's
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
#' The start age is fixed unless it has a prior in `agePriors`.
#' With `start = "root"`, `rootAge` scales the distance between the root and
#' its older child by \eqn{e^{u}}, \eqn{u \sim U(-\epsilon/2, \epsilon/2)},
#' rebuilds the frozen model at the new root age (its E grid and slices end
#' there) and redraws the paths on the two root branches between their
#' fixed end bins.
#' The survival conditioning is evaluated at the current root age, so the
#' target is the root-age prior times the likelihood conditioned on survival
#' from that age (spec section 6).
#' With `start = "origin"`, the stem carries a path from the origin bin
#' (drawn from the start-bin prior) to the root, the root carries a
#' \eqn{\lambda} factor, survival is that of the single stem lineage, and
#' `originAge` scales the stem length in the same way, redrawing the stem's
#' path.
#' Step sizes \eqn{\epsilon} adapt during burn-in.
#'
#' Characters follow the Mk model on branch lengths
#' \eqn{\ell_e = c \int_e e^{\gamma x} dt}, with \eqn{c} the parameter `rate`
#' and \eqn{\gamma} from `model`; see [MkLikelihood()].
#'
#' @inheritParams MarginalMCMC
#' @inheritParams MkLikelihood
#' @param nSlice Integer giving the number of time slices of the proposal.
#' @param moves Named numeric vector giving the number of each move per
#' iteration (see Details); unnamed moves are not used.  Defaults to
#' [TempoMoves()] for `tree` and `characters`.
#' @param characters Optional matrix of character states, as for
#' [MkLikelihood()].
#' @param updateParameters Logical; if `FALSE`, hold parameters fixed (`shift`
#' then moves the bins only).
#' @param minAcceptance Numeric: `tree` and `joint` moves whose acceptance
#' during burn-in falls below this are not used after burn-in.
#' @param agePriors Named list of functions giving log prior densities of
#' the root age (`root`; with `start = "root"`) and of the origin age
#' (`origin`; with `start = "origin"`), needed by the `rootAge` and
#' `originAge` moves.
#' @returns `TempoMCMC()` returns a list with elements
#' - `samples`: data frame with one column per free parameter, `logPost`
#'   (log augmented density, plus character log-likelihood and log prior),
#'   `rootBin`, `meanX` (lineage-time-averaged log-tempo over the tree) and
#'   tree statistics: `nAncestor`, `treeLength`, `meanNodeAge` (of non-root
#'   bifurcations), `cherries`, `oldestFossil`, `youngestFossil`, `clock`
#'   (\eqn{\sum_e \int_e e^x dt}) and `inverseClock`
#'   (\eqn{\sum_e \int_e e^{-x} dt}), over the branches below the root,
#'   `rootAge`, and with `start = "origin"`, `originAge`;
#' - `occupancy`: matrix of time spent in each bin (columns), per sample;
#' - `branchX`: matrix of the time-averaged log-tempo on each branch (for
#'   chains without tree moves);
#' - `acceptance`: acceptance rate of each move;
#' - `moves`: the move counts used after burn-in;
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
#' @importFrom stats setNames
#' @export
TempoMCMC <- function(tree, model, start, priors, nIter, burnin = 0, thin = 1,
                      nSlice = 20, moves = NULL, characters = NULL,
                      nState = 2L, updateParameters = TRUE, minAcceptance = 0.01,
                      agePriors = list(),
                      engine = getOption("EvoTempo.engine", "cpp")) {
  old <- options(EvoTempo.engine = TempoEngine(engine))
  on.exit(options(old))
  tree <- .ForModel(TempoTree(tree), model)
  origin <- .Origin(tree)
  moveNames <- c("tree", "node", "tip", "branch", "joint", "shift", "nodeAge",
                 "fossilAge", "slide", "saJump", "rootAge", "originAge")
  if (is.null(moves)) {
    moves <- TempoMoves(tree, characters, model)
  }
  if (!all(names(moves) %in% moveNames)) {
    stop("Unknown move: ", paste(setdiff(names(moves), moveNames), collapse = ", "))
  }
  moves <- replace(setNames(numeric(length(moveNames)), moveNames), names(moves),
                   moves)
  ageMoves <- c(rootAge = "root", originAge = "origin")
  for (move in names(ageMoves)[moves[names(ageMoves)] > 0]) {
    if ((move == "originAge") != origin) {
      stop("`", move, "` needs `start = \"", if (origin) "root" else "origin", "\"`",
           if (origin) "; with an origin, `nodeAge` moves the root" else "")
    }
    if (!is.function(agePriors[[ageMoves[[move]]]])) {
      stop("`", move, "` needs a prior `agePriors$", ageMoves[[move]], "`")
    }
  }
  ageStep <- c(rootAge = 1, originAge = 1)
  ageHits <- 0 * ageStep
  par <- as.list(start)
  if (!is.null(characters)) {
    characters <- as.matrix(characters)
  }
  CharLL <- .CharacterLikelihood(characters, nState, model[["gamma"]])
  set <- .AddSlices(tree, .Precompute(tree, par, model), nSlice)
  msg <- .Messages(tree, set)
  draw <- .SampleTree(tree, set, msg)
  st <- list(tree = tree, par = par, set = set, msg = msg,
             paths = draw[["paths"]], bin = draw[["bin"]])
  st[["charLL"]] <- CharLL(tree, st[["paths"]], par, set)
  .ResolveCharacters(CharLL, TRUE)
  # Reuses the current set when only parameters outside E change (`rate`)
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
  pathNodes <- .PathNodes(tree)
  nKeep <- nIter %/% thin
  statNames <- c(names(.TreeStats(tree)), if (origin) "originAge")
  samples <- matrix(NA_real_, nKeep, length(priors) + 3L + length(statNames),
                    dimnames = list(NULL, c(names(priors), "logPost", "rootBin",
                                            "meanX", statNames)))
  occupancy <- matrix(NA_real_, nKeep, model[["nBin"]])
  branchX <- matrix(NA_real_, nKeep, length(pathNodes))
  moveTry <- moveAcc <- setNames(numeric(length(moveNames) + 1L),
                                 c(moveNames, "parameters"))
  treeMoves <- c("nodeAge", "fossilAge", "slide", "saJump")
  # Character lengths below the root, for the clock statistics
  Below <- function(tr, gamma) {
    replace(.CharacterLengths(tr, st[["set"]][["grid"]][["x"]], st[["paths"]], 1,
                              gamma), tr[["root"]], 0)
  }
  free <- if (updateParameters) names(priors) else character(0)
  for (i in seq_len(burnin + nIter)) {
    for (move in moveNames[moves > 0]) {
      for (rep in seq_len(.MoveCount(moves[[move]]))) {
        if (move == "shift") {
          res <- .ShiftMove(st, model, free, priors, nSlice, CharLL)
          st <- res[["state"]]
          moveTry[[move]] <- moveTry[[move]] + 1
          moveAcc[[move]] <- moveAcc[[move]] + res[["accepted"]]
          next
        }
        if (move %in% names(ageMoves)) {
          res <- .AgeMove(move, st, model, nSlice, CharLL,
                          agePriors[[ageMoves[[move]]]], ageStep[[move]])
          st <- res[["state"]]
          moveTry[[move]] <- moveTry[[move]] + 1
          moveAcc[[move]] <- moveAcc[[move]] + res[["accepted"]]
          ageHits[[move]] <- ageHits[[move]] + res[["accepted"]]
          next
        }
        if (move == "joint") {
          if (!updateParameters || !length(priors)) next
          res <- .JointSweep(st, priors, jointStep, model, nSlice, CharLL)
          jointHits <- jointHits + res[["accepted"]]
          moveTry[[move]] <- moveTry[[move]] + length(priors)
          moveAcc[[move]] <- moveAcc[[move]] + sum(res[["accepted"]])
          st <- res[["state"]]
          next
        }
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
      sweep <- .ParameterSweep(st[["par"]], priors, step, current, Target,
                               function(a) .ResolveCharacters(CharLL, a))
      moveTry[["parameters"]] <- moveTry[["parameters"]] + length(priors)
      moveAcc[["parameters"]] <- moveAcc[["parameters"]] + sum(sweep[["accepted"]])
      hits <- hits + sweep[["accepted"]]
      if (any(sweep[["accepted"]])) {
        st[["par"]] <- sweep[["par"]]
        if (!.SameSet(st[["set"]], st[["tree"]], st[["par"]])) {
          st[["set"]] <- .AddSlices(st[["tree"]], sweep[["current"]][["set"]], nSlice)
          st[["msg"]] <- .Messages(st[["tree"]], st[["set"]])
        }
        st[["charLL"]] <- CharLL(st[["tree"]], st[["paths"]], st[["par"]], st[["set"]])
        .ResolveCharacters(CharLL, TRUE)
      }
      if (i <= burnin && i %% 20 == 0) {
        step <- .Adapt(step, hits / 20)
        hits <- 0 * hits
        if (moves[["joint"]] > 0) {
          jointStep <- .Adapt(jointStep, jointHits / (20 * moves[["joint"]]))
          jointHits <- 0 * jointHits
        }
      }
    }
    if (i <= burnin && i %% 20 == 0) {
      for (move in names(ageMoves)[moves[names(ageMoves)] > 0]) {
        ageStep[[move]] <- .Adapt(ageStep[[move]], ageHits[[move]] / (20 * moves[[move]]))
      }
      ageHits <- 0 * ageHits
    }
    if (i == burnin) {
      # Whole-tree draws that the characters reject are dropped after burn-in
      for (move in c("tree", "joint")) {
        if (moves[[move]] > 0 && moveTry[[move]] > 0 &&
            moveAcc[[move]] / moveTry[[move]] < minAcceptance) {
          moves[[move]] <- 0
        }
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
        .TreeStats(tr, Below(tr, 1), Below(tr, -1)),
        if (origin) .StartAge(tr))
      occupancy[j / thin, ] <- summ[["occupancy"]]
      if (length(summ[["branchX"]]) == ncol(branchX)) {
        branchX[j / thin, ] <- summ[["branchX"]]
      }
    }
  }
  list(samples = as.data.frame(samples), occupancy = occupancy,
       branchX = branchX, acceptance = moveAcc / pmax(moveTry, 1),
       moves = moves[moves > 0], state = st)
}

#' Default move schedule of TempoMCMC()
#'
#' Move counts per iteration for [TempoMCMC()], following the mixing study of
#' issue #19 (`dev/mixing/README.md`).
#'
#' Without characters, or with \eqn{\gamma = 0} (when the characters do not
#' depend on the paths), the schedule is one whole-tree draw (`tree`), two
#' each of `node`, `tip` and `branch`, one `joint` sweep and one `shift`.
#' With characters that depend on tempo, whole-tree draws from the frozen
#' model are almost never accepted, so `tree` and `joint` are omitted, and
#' `node`, `tip` and `branch` moves are run for a fraction `local` of the
#' internal nodes, the tips with paths and the branches with paths
#' respectively (rounded up), with one `shift`.
#' Tree moves (`nodeAge`, `fossilAge`, `slide`, `saJump`) are not included;
#' add them to the returned vector.
#'
#' @param tree A [`TempoTree`][TempoTree()], or a tree that it accepts.
#' @param characters Optional matrix of character states, as for
#' [TempoMCMC()].
#' @param model Optional [`DCETModel`][DCETModel()]; only its `gamma` is used.
#' @param local Numeric giving the number of each local path move per
#' iteration, as a fraction of the elements it updates.
#' @returns `TempoMoves()` returns a named numeric vector of move counts, for
#' the `moves` argument of [TempoMCMC()].
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' TempoMoves(tree)
#' TempoMoves(tree, characters = matrix(c(0, 1, 1), 3))
#' c(TempoMoves(tree), nodeAge = 2, slide = 1)
#' @family dCET inference
#' @export
TempoMoves <- function(tree, characters = NULL, model = NULL, local = 0.25) {
  tree <- TempoTree(tree)
  scored <- !is.null(characters) && NCOL(characters) > 0 &&
    (is.null(model) || model[["gamma"]] != 0)
  if (!scored) {
    return(c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1, shift = 1))
  }
  nodes <- .PathNodes(tree)
  c(node = ceiling(local * length(tree[["postorder"]])),
    tip = ceiling(local * sum(nodes <= tree[["nTip"]])),
    branch = ceiling(local * length(nodes)),
    shift = 1)
}

# Number of repeats of a move this iteration: a fractional count n + f gives n,
# or n + 1 with probability f, independently of the state
.MoveCount <- function(n) {
  if (n == floor(n)) {
    return(n)
  }
  floor(n) + (stats::runif(1) < n - floor(n))
}

# Shift move (#19; dev/mixing/README.md): every bin moves by m = +/-1, so the
# log-tempo moves by d = m h s everywhere, and each free parameter p that
# tempo scales is multiplied by exp(-e_p d).  The map is a bijection with
# inverse the shift by -m, its Jacobian in natural scale is prod p' / p, and
# with `rate` free (or gamma = 0) the character lengths, and so the character
# likelihood, are unchanged exactly, so the stored value is kept.  A shift
# that would leave the grid is rejected; so would its reverse.
.ShiftMove <- function(st, model, free, priors, nSlice, CharLL) {
  reject <- list(state = st, accepted = FALSE)
  m <- if (stats::runif(1) < 0.5) -1L else 1L
  new <- .Shifted(st, m, model, free)
  if (is.null(new)) {
    return(reject)
  }
  tree <- st[["tree"]]
  logPrior <- sum(vapply(new[["comp"]], function(p) {
    priors[[p]](new[["par"]][[p]]) - priors[[p]](st[["par"]][[p]])
  }, 0))
  if (!is.finite(logPrior)) {
    return(reject)
  }
  same <- .SameSet(st[["set"]], tree, new[["par"]])
  set <- if (same) {
    st[["set"]]
  } else {
    tryCatch(.Precompute(tree, new[["par"]], model),
             infeasibleProposal = function(e) NULL)
  }
  if (is.null(set)) {
    return(reject)
  }
  invariant <- model[["gamma"]] == 0 || "rate" %in% new[["comp"]]
  charLL <- if (invariant) {
    st[["charLL"]]
  } else {
    CharLL(tree, new[["paths"]], new[["par"]], set)
  }
  logRatio <- .LogAugmented(tree, set, new[["paths"]], new[["bin"]]) + charLL -
    .LogAugmented(tree, st[["set"]], st[["paths"]], st[["bin"]]) - st[["charLL"]] +
    logPrior + new[["logJac"]]
  accept <- is.finite(logRatio) && log(stats::runif(1)) < logRatio
  .ResolveCharacters(CharLL, accept)
  if (accept) {
    if (!same) {
      st[["set"]] <- .AddSlices(tree, set, nSlice)
      st[["msg"]] <- .Messages(tree, st[["set"]])
    }
    st[["par"]] <- new[["par"]]
    st[["bin"]] <- new[["bin"]]
    st[["paths"]] <- new[["paths"]]
    st[["charLL"]] <- charLL
    return(list(state = st, accepted = TRUE))
  }
  reject
}

# The state shifted by m bins, with the free parameters in `free` that tempo
# scales compensated; NULL if a bin would leave the grid.  Returns the new
# parameters, bins and paths, the compensated parameters `comp` and the log
# Jacobian (in natural scale) of the map.
.Shifted <- function(st, m, model, free) {
  nodes <- .PathNodes(st[["tree"]])
  bin <- st[["bin"]]
  paths <- st[["paths"]]
  used <- bin > 0
  allBins <- c(bin[used], unlist(lapply(paths[nodes], `[[`, "k")))
  if (any(allBins + m < 1L | allBins + m > model[["nBin"]])) {
    return(NULL)
  }
  d <- m * st[["set"]][["grid"]][["h"]] * st[["par"]][["s"]]
  expo <- c(rate = model[["gamma"]], lambda = model[["beta"]], mu = model[["beta"]],
            psi = model[["delta"]], theta = model[["alpha"]])
  comp <- intersect(names(expo)[expo != 0], free)
  par <- st[["par"]]
  for (p in comp) {
    par[[p]] <- par[[p]] * exp(-expo[[p]] * d)
  }
  bin[used] <- bin[used] + m
  for (v in nodes) {
    paths[[v]][["k"]] <- paths[[v]][["k"]] + m
  }
  list(par = par, bin = bin, paths = paths, comp = comp,
       logJac = -sum(expo[comp]) * d)
}

# Joint parameter-path moves: for each free parameter, propose a new value and
# redraw every path by FFBS under it.  Because exact / proposal = frozen
# likelihood x exact-to-frozen weight, the acceptance ratio needs no proposal
# density for the paths; the move mixes like MCMC on the marginal likelihood.
.JointSweep <- function(st, priors, step, model, nSlice, CharLL,
                        Draw = .SampleTree) {
  tree <- st[["tree"]]
  nodes <- .PathNodes(tree)
  accepted <- setNames(logical(length(priors)), names(priors))
  Value <- function(set, msg, paths, charLL) {
    .FrozenLikelihood(tree, set, msg, set[["joint"]]) +
      .LogWeight(tree, set, msg, paths, nodes) + charLL
  }
  current <- Value(st[["set"]], st[["msg"]], st[["paths"]], st[["charLL"]])
  for (p in names(priors)) {
    proposal <- st[["par"]]
    factor <- exp(step[[p]] * (stats::runif(1) - 0.5))
    proposal[[p]] <- proposal[[p]] * factor
    lp <- priors[[p]](proposal[[p]])
    if (!is.finite(lp)) {
      next
    }
    draw <- tryCatch({
      if (.SameSet(st[["set"]], tree, proposal)) {
        # `rate` does not enter the frozen model: only the paths are redrawn
        set <- st[["set"]]
        msg <- st[["msg"]]
      } else {
        set <- .AddSlices(tree, .Precompute(tree, proposal, model), nSlice)
        msg <- .Messages(tree, set)
      }
      Draw(tree, set, msg)
    }, infeasibleProposal = function(e) NULL)
    if (is.null(draw)) {
      next
    }
    charLL <- CharLL(tree, draw[["paths"]], proposal, set)
    candidate <- Value(set, msg, draw[["paths"]], charLL)
    logRatio <- candidate + lp - current - priors[[p]](st[["par"]][[p]]) + log(factor)
    accept <- is.finite(logRatio) && log(stats::runif(1)) < logRatio
    .ResolveCharacters(CharLL, accept)
    if (accept) {
      st[["par"]] <- proposal
      st[["set"]] <- set
      st[["msg"]] <- msg
      st[["paths"]] <- draw[["paths"]]
      st[["bin"]] <- draw[["bin"]]
      st[["charLL"]] <- charLL
      current <- candidate
      accepted[[p]] <- TRUE
    }
  }
  list(state = st, accepted = accepted)
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
  accept <- is.finite(logRatio) && log(stats::runif(1)) < logRatio
  .ResolveCharacters(CharLL, accept)
  if (accept) {
    st[["tree"]] <- newTree
    st[["paths"]] <- prop[["paths"]]
    st[["bin"]] <- prop[["bin"]]
    st[["charLL"]] <- newChar
    st[["msg"]] <- .Messages(newTree, set)
    return(list(state = st, accepted = TRUE))
  }
  list(state = st, accepted = FALSE)
}

# Move of the start age (`rootAge` with start = "root", `originAge` with an
# origin): the distance from the start to its older child is scaled by
# exp(step (U - 1/2)); the frozen model is rebuilt at the new tree (its E grid
# and slices end at the start age); the paths on the branches from the start
# (the root's two daughters, or the stem) are redrawn as frozen bridges
# between their fixed end bins.  Accepted on the exact target, which includes
# the survival conditioning at the new start age, with the age prior, the
# Jacobian of the scaling and the bridge densities under the frozen model of
# each state in the Hastings ratio.
.AgeMove <- function(move, st, model, nSlice, CharLL, prior, step) {
  reject <- list(state = st, accepted = FALSE)
  tree <- st[["tree"]]
  root <- tree[["root"]]
  par <- st[["par"]]
  oldAge <- .StartAge(tree)
  newTree <- tree
  factor <- exp(step * (stats::runif(1) - 0.5))
  if (move == "rootAge") {
    nodes <- tree[["kids"]][root - tree[["nTip"]], ]
    bottom <- max(tree[["age"]][nodes])
    newTree[["age"]][root] <- bottom + (oldAge - bottom) * factor
  } else {
    nodes <- root
    newTree[["stem"]] <- tree[["stem"]] * factor
  }
  newAge <- .StartAge(newTree)
  logPrior <- prior(newAge) - prior(oldAge)
  if (!is.finite(logPrior)) {
    return(reject)
  }
  prop <- tryCatch({
    # The reverse move would draw from a frozen model built afresh at the
    # current tree
    oldSet <- if (.SameSet(st[["set"]], tree, par)) {
      st[["set"]]
    } else {
      .AddSlices(tree, .Precompute(tree, par, model), nSlice)
    }
    set <- .AddSlices(newTree, .Precompute(newTree, par, model), nSlice)
    list(oldSet = oldSet, set = set,
         fwd = .DrawBridges(newTree, set, nodes, st[["bin"]], st[["paths"]]))
  }, infeasibleProposal = function(e) NULL)
  if (is.null(prop)) {
    return(reject)
  }
  set <- prop[["set"]]
  newPaths <- prop[["fwd"]][["paths"]]
  newChar <- CharLL(newTree, newPaths, par, set)
  logRatio <- .LogAugmented(newTree, set, newPaths, st[["bin"]]) -
    .LogAugmented(tree, prop[["oldSet"]], st[["paths"]], st[["bin"]]) +
    newChar - st[["charLL"]] + logPrior + log(factor) +
    .LogBridges(tree, prop[["oldSet"]], nodes, st[["paths"]]) - prop[["fwd"]][["logQ"]]
  accept <- is.finite(logRatio) && log(stats::runif(1)) < logRatio
  .ResolveCharacters(CharLL, accept)
  if (accept) {
    st[["tree"]] <- newTree
    st[["set"]] <- set
    st[["msg"]] <- .Messages(newTree, set)
    st[["paths"]] <- newPaths
    st[["charLL"]] <- newChar
    return(list(state = st, accepted = TRUE))
  }
  reject
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
    int <- .Intervals(st[["paths"]][[v]], .TopAge(tree, v), age[[v]])
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
                                       .TopBin(tree, bin, v), bin[[v]])
        newBin <- bin
      },
      node = {
        # Internal nodes, and with an origin the origin bin
        internal <- tree[["postorder"]]
        origin <- .Origin(tree)
        root <- tree[["root"]]
        pick <- .Sample1(rep(1, length(internal) + origin))
        newBin <- bin
        newPaths <- paths
        startW <- set[["prior"]] / (if (set[["joint"]]) set[["survive"]] else 1)
        if (pick > length(internal)) {
          top <- length(bin)
          newBin[top] <- .Sample1(startW * .TransferCol(msg[["pieces"]][[root]],
                                                        bin[[root]]))
          changed <- root
          newPaths[[root]] <- .SampleBranch(set, msg[["pieces"]][[root]], NULL,
                                            newBin[[top]], bin[[root]])
        } else {
          v <- internal[[pick]]
          kids <- .PathKids(tree, v)
          w <- if (v == root && !origin) {
            startW
          } else if (tree[["sampled"]][v - nTip]) {
            set[["psi"]] * .TransferRow(msg[["pieces"]][[v]], .TopBin(tree, bin, v))
          } else {
            set[["lambda"]] * .TransferRow(msg[["pieces"]][[v]], .TopBin(tree, bin, v))
          }
          for (kid in kids) {
            w <- w * .TransferCol(msg[["pieces"]][[kid]], bin[[kid]])
          }
          newBin[v] <- .Sample1(w)
          sa <- tree[["kids"]][v - nTip, ]
          newBin[sa[sa <= nTip & tree[["tipType"]][pmin(sa, nTip)] == "ancestor"]] <- newBin[v]
          changed <- c(if (v != root || origin) v, kids)
          for (u in changed) {
            newPaths[[u]] <- .SampleBranch(set, msg[["pieces"]][[u]], NULL,
                                           .TopBin(tree, newBin, u), newBin[[u]])
          }
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
  accept <- log(stats::runif(1)) < logRatio
  .ResolveCharacters(CharLL, accept)
  if (accept) {
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
SampleTempoPaths <- function(tree, par, model, n, nSlice = 20,
                             engine = getOption("EvoTempo.engine", "cpp")) {
  old <- options(EvoTempo.engine = TempoEngine(engine))
  on.exit(options(old))
  tree <- .ForModel(TempoTree(tree), model)
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
