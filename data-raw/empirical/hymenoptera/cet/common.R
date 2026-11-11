# Hymenoptera CET pipeline: definitions shared by run-chain.R, check-mkv.R,
# bayes-factors.R and check-prior.R.  Source from the package root after
# loading EvoTempo.  Follows data-raw/empirical/penguins/cet/common.R; the
# differences are marked "Hymenoptera:".
#
# The engine (EvoTempo main @ 48e6bf7) lacks several things this analysis
# needs.  As in the penguin pipeline, this file rebinds functions inside the
# EvoTempo namespace for the duration of a run, without changing the package:
#
# 1. MkLikelihood() scores unordered (Jukes-Cantor) characters with one
#    number of states, one state or NA per cell, and no ascertainment
#    correction.  Hymenoptera: the matrix has 43 ordered characters (K = 3-5)
#    and 22 ambiguous cells, so the rebound MkLikelihood() (HymMk()) scores
#    each (ordering, K) partition with its own rate matrix through the C++
#    kernel mk-columns.cpp (any 0/1 tip partials, so ambiguity sets are summed
#    over, as in RevBayes), and applies the Mkv correction per missing-data
#    mask, as RevBayes dnPhyloCTMC(coding = "variable") (#65, #80).  It
#    raises the character likelihood, correction included, to a power
#    (`.cet$power`) for stepping-stone runs.
# 2. Hymenoptera: the Holometabola constraint and calibration of the
#    comparator (a clade constraint and an offset-exponential density on the
#    clade's age) are tree-only terms with no place in TempoMCMC().  HymMk()
#    adds them, unpowered, to its return value (TreeTerm()).  This is exact:
#    TempoMCMC() rescores the characters on every proposed tree (tree and
#    start-age moves) and reuses the stored value only when the tree is
#    unchanged (path, shift and parameter moves).  The recorded charLL
#    excludes the term, and RecordingStats() checks the decomposition at every
#    sample.
# 3. Hymenoptera: the comparator's prior on the FBD rates is on net
#    diversification d = lambda - mu ~ Exp(100), turnover mu / lambda ~ U(0, 1)
#    and sampling proportion psi / (mu + psi) ~ U(0, 1), which is not a
#    product of priors on lambda, mu, psi.  TempoMCMC() evaluates priors one
#    parameter at a time, as priors[[p]](value), at four call sites:
#    .ParameterSweep() and .JointSweep() (one parameter changes, the others
#    are at their current values), .ShiftMove() (lambda and mu are rescaled
#    together when beta != 0) and .LogPrior() (all parameters, for logPost).
#    FbdConditional() gives lambda, mu and psi priors that evaluate the joint
#    density at the caller's current parameters with one entry replaced; they
#    identify their caller and stop() at any other call site.  .ShiftMove()
#    is wrapped (ShiftJoint()) so that the rescaled pair is scored jointly,
#    and .LogPrior() is rebound (LogPriorJoint()).  check-prior.R checks the
#    result against an independent sampler.
# 4. TempoMCMC() returns no per-sample trees or states and cannot checkpoint
#    or resume; the rebound .TreeStats() (RecordingStats(), as in the penguin
#    pipeline) records them and continues a captured chain exactly.
# See README.md and the engine issues #47 and #48.
suppressMessages(library(ape))
.ns <- asNamespace("EvoTempo")
.cet <- new.env()
.cet$power <- 1
.cet$mkv <- Sys.getenv("CET_MKV", "pattern")
stopifnot("CET_MKV must be \"pattern\" or \"none\"" = .cet$mkv %in% c("pattern", "none"))
MkvSuffix <- function(mkv = .cet$mkv) if (mkv == "pattern") "" else paste0("_mkv", mkv)
.cet$records <- list()
.cet$resume <- NULL
.cet$captured <- NULL
.cet$injected <- FALSE
if (!exists(".mkOriginal")) {
  .mkOriginal <- get("MkLikelihood", .ns)
  .statsOriginal <- get(".TreeStats", .ns)
  .logPriorOriginal <- get(".LogPrior", .ns)
  .shiftOriginal <- get(".ShiftMove", .ns)
  .sweepOriginal <- get(".ParameterSweep", .ns)
  .jointOriginal <- get(".JointSweep", .ns)
}
cetDir <- "data-raw/empirical/hymenoptera/cet"
Rcpp::sourceCpp(file.path(cetDir, "mk-columns.cpp"),
                cacheDir = Sys.getenv("CET_RCPP_CACHE", file.path(tempdir(), "cet-rcpp")))

# ---- Characters -----------------------------------------------------------------
# Rate matrix of a partition, rescaled to mean rate 1 under its (uniform)
# stationary distribution, and its eigen decomposition (symmetric, so
# P(t) = V diag(exp(ev t)) V').  Unordered: Jukes-Cantor (RevBayes fnJC(K)).
# Ordered: changes between adjacent states only, at equal rates
# (fnOrderedRateMatrix(maxState = K), RevBayes 1.4 / dev); mean rate
# 2 (K - 1) r / K = 1, so r = K / (2 (K - 1)).
RateMatrix <- function(K, ordered) {
  if (ordered) {
    Q <- matrix(0, K, K)
    r <- K / (2 * (K - 1))
    Q[cbind(1:(K - 1), 2:K)] <- r
    Q[cbind(2:K, 1:(K - 1))] <- r
  } else {
    Q <- matrix(1 / (K - 1), K, K)
  }
  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)
  stopifnot(isSymmetric(Q), abs(-mean(diag(Q)) - 1) < 1e-12)
  e <- eigen(Q, symmetric = TRUE)
  list(Q = Q, V = e$vectors, ev = e$values)
}

# Tip partials of one cell: a state, a set "{01}", or "" (unknown)
CellPartial <- function(cell, K) {
  p <- numeric(K)
  if (!nzchar(cell)) return(p + 1)
  p[as.integer(strsplit(gsub("[^0-9]", "", cell), "")[[1]]) + 1L] <- 1
  stopifnot(sum(p) >= 1)
  p
}

# One entry per (ordering, K) partition: the kernel inputs for its data
# columns followed by its constant columns.  Mkv (#65): a missing-data mask
# is the set of tips with a scored cell ('?' and '-' unknown; an ambiguity set
# counts as scored, as for RevBayes, which treats only '?' and '-' as gaps);
# for each distinct mask within a partition, the constant columns hold state
# j at the scored tips and unknown elsewhere.  Under Jukes-Cantor with a
# uniform root P(all scored tips in state j) is the same for every j, so
# unordered partitions need one constant column per mask (state 0) and use
# K x its probability; ordered partitions have one column per state.
BuildPartitions <- function(cells, partitions) {
  nTip <- nrow(cells)
  labels <- unique(partitions$partition[order(partitions$ordered, partitions$nState)])
  lapply(labels, function(lab) {
    cols <- which(partitions$partition == lab)
    K <- partitions$nState[cols[1]]
    ordered <- partitions$ordered[cols[1]]
    x <- cells[, cols, drop = FALSE]
    # A cell is scored unless unknown or an ambiguity set covering all K
    # states (in a binary character, {01}): RevBayes treats both as gaps for
    # the Mkv masks (check-mkv.R: counting {01} as scored made u2 differ by
    # up to 6.8e-3 while every K >= 3 partition matched to 1e-10)
    full <- matrix(nchar(gsub("[^0-9]", "", x)) == K & grepl("^\\{", x), nrow(x))
    observed <- matrix(nzchar(x), nrow(x)) & !full
    key <- apply(observed, 2, function(o) paste(which(o), collapse = ","))
    first <- !duplicated(key)
    maskObs <- observed[, first, drop = FALSE]
    count <- tabulate(match(key, key[first]), sum(first))
    stopifnot(all(colSums(maskObs) >= 2), sum(count) == length(cols))
    nConstState <- if (ordered) K else 1L
    nData <- length(cols)
    nMask <- ncol(maskObs)
    nCol <- nData + nMask * nConstState
    tip <- array(0, c(K, nCol, nTip))
    for (i in seq_len(nTip)) {
      for (j in seq_len(nData)) tip[, j, i] <- CellPartial(x[i, j], K)
      for (m in seq_len(nMask)) for (s in seq_len(nConstState)) {
        tip[, nData + (m - 1L) * nConstState + s, i] <-
          if (maskObs[i, m]) as.numeric(seq_len(K) == s) else rep(1, K)
      }
    }
    c(list(label = lab, K = K, ordered = ordered, columns = cols, nData = nData,
           nMask = nMask, nConstState = nConstState, count = count, nCol = nCol,
           tip = as.vector(tip), maskObserved = maskObs),
      RateMatrix(K, ordered))
  })
}

# Pure-R reference of MkColumnLogLik() (check-mkv.R): plain pruning with
# per-node rescaling
ColumnLogLikR <- function(tree, g, lengths) {
  nTip <- tree[["nTip"]]
  K <- g$K
  tip <- array(g$tip, c(K, g$nCol, nTip))
  part <- vector("list", nTip + tree[["nNode"]])
  for (i in seq_len(nTip)) part[[i]] <- matrix(tip[, , i], K)
  logScale <- numeric(g$nCol)
  for (v in tree[["postorder"]]) {
    a <- matrix(1, K, g$nCol)
    for (kid in tree[["kids"]][v - nTip, ]) {
      P <- g$V %*% diag(exp(g$ev * lengths[[kid]]), K) %*% t(g$V)
      a <- a * (P %*% part[[kid]])
    }
    mx <- apply(a, 2, max)
    part[[v]] <- sweep(a, 2, mx, "/")
    logScale <- logScale + log(mx)
  }
  log(colSums(part[[tree[["root"]]]]) / K) + logScale
}

# Data log-likelihood and Mkv correction of one partition
PartitionLogLik <- function(tree, g, lengths, impl = "cpp") {
  ll <- if (impl == "cpp") {
    MkColumnLogLik(as.integer(tree[["postorder"]]), tree[["kids"]], tree[["nTip"]],
                   tree[["root"]], as.double(lengths), g$V, g$ev, g$tip, g$nCol)
  } else {
    ColumnLogLikR(tree, g, lengths)
  }
  data <- sum(ll[seq_len(g$nData)])
  pConst <- exp(ll[-seq_len(g$nData)])
  pMask <- if (g$ordered) colSums(matrix(pConst, g$nConstState)) else g$K * pConst
  c(data = data, correction = -sum(g$count * log1p(-pMask)))
}

# Character log-likelihood, Mkv-corrected, raised to `power`, plus the
# unpowered tree-only term (constraint and calibration) if `treeTerm`
HymMk <- function(tree, characters, lengths, nState = 2L,
                  engine = getOption("EvoTempo.engine", "cpp"),
                  power = .cet$power, mkv = .cet$mkv, treeTerm = TRUE, impl = "cpp") {
  parts <- attr(characters, "hym")
  if (is.null(parts)) {
    return(.mkOriginal(tree, characters, lengths, nState, engine))
  }
  tree <- TempoTree(tree)
  tt <- if (treeTerm) TreeTerm(tree) else 0
  if (power == 0 || !is.finite(tt)) {
    return(tt)
  }
  ll <- 0
  for (g in parts) {
    v <- PartitionLogLik(tree, g, lengths, impl)
    ll <- ll + v[["data"]] + switch(mkv, pattern = v[["correction"]], none = 0,
                                    stop("`mkv` must be \"pattern\" or \"none\""))
  }
  power * ll + tt
}

# ---- Tree-only terms: Holometabola constraint and calibration -------------------
# Age of the node whose descendant tips are exactly `members` (NA if they do
# not form a clade), or with `mrca = TRUE` the age of their most recent common
# ancestor
CladeAge <- function(tree, members, mrca = FALSE) {
  nTip <- tree[["nTip"]]
  isM <- tree[["tipLabel"]] %in% members
  nM <- c(as.integer(isM), integer(tree[["nNode"]]))
  nT <- c(rep(1L, nTip), integer(tree[["nNode"]]))
  kids <- tree[["kids"]]
  for (v in tree[["postorder"]]) {
    k <- kids[v - nTip, ]
    nM[v] <- nM[k[1]] + nM[k[2]]
    nT[v] <- nT[k[1]] + nT[k[2]]
  }
  n <- sum(isM)
  if (n < 2L) return(NA_real_)
  hit <- if (mrca) which(nM == n)[which.min(nT[nM == n])] else which(nM == n & nT == n)
  if (length(hit)) tree[["age"]][hit[1]] else NA_real_
}
# Holometabola (with all fossils; 111 of 113 tips) monophyletic, as the
# comparator's dnConstrainedTopology(constraints = v(clade_holo)), and its age
# t ~ offset exponential, minimum 302 Ma, mean 396 Ma (the comparator's
# clamped-exponential calibration, MrBayes `calibrate`): exp(-(t - 302)/94)/94.
# Unnormalised, as in RevBayes (no renormalisation of the tree prior)
HoloCalibration <- function(age) if (age > 302) -log(94) - (age - 302) / 94 else -Inf
TreeTerm <- function(tree) {
  a <- CladeAge(tree, .cet$holo)
  if (is.na(a)) -Inf else HoloCalibration(a)
}

# ---- FBD-rate prior of the comparator, as conditional priors --------------------
# log density of (lambda, mu, psi) under d = lambda - mu ~ Exp(100),
# r = mu / lambda ~ U(0, 1), s = psi / (mu + psi) ~ U(0, 1):
#   p(lambda, mu, psi) = 100 exp(-100 d) |d(d, r, s) / d(lambda, mu, psi)|,
#   |J| = (lambda - mu) mu / (lambda^2 (mu + psi)^2).
FbdLogPrior <- function(par) {
  l <- par[["lambda"]]
  m <- par[["mu"]]
  p <- par[["psi"]]
  if (!(is.finite(l) && is.finite(m) && is.finite(p) && m > 0 && l > m && p > 0)) {
    return(-Inf)
  }
  d <- l - m
  log(100) - 100 * d + log(d) + log(m) - 2 * log(l) - 2 * log(m + p)
}
fbdNames <- c("lambda", "mu", "psi")
# Prior for parameter `p` in TempoMCMC()'s `priors`: the joint density at the
# caller's current parameters with `p` set to x.  The current parameters are
# `par` in .ParameterSweep() and `st$par` in .JointSweep(); in both, every
# call changes `p` only.  Any other caller is an error (.ShiftMove() and
# .LogPrior() never call these functions: see ShiftJoint(), LogPriorJoint())
FbdConditional <- function(p) {
  force(p)
  f <- function(x) {
    caller <- sys.function(sys.parent())
    pf <- parent.frame()
    cur <- if (identical(caller, .sweepOriginal)) {
      get("par", envir = pf, inherits = FALSE)
    } else if (identical(caller, .jointOriginal)) {
      get("st", envir = pf, inherits = FALSE)[["par"]]
    } else {
      stop("FBD-rate prior on `", p, "` called from an unexpected function; ",
           "common.R's joint prior supports only .ParameterSweep() and .JointSweep()")
    }
    cur[[p]] <- x
    FbdLogPrior(cur)
  }
  attr(f, "cetJoint") <- p
  f
}
IsJoint <- function(priors) {
  any(vapply(priors, function(f) !is.null(attr(f, "cetJoint")), TRUE))
}
# .LogPrior() for logPost: the joint FBD-rate density plus the other priors
LogPriorJoint <- function(par, priors) {
  if (!IsJoint(priors)) return(.logPriorOriginal(par, priors))
  stopifnot(all(fbdNames %in% names(priors)))
  other <- setdiff(names(priors), fbdNames)
  FbdLogPrior(par) + sum(vapply(other, function(p) priors[[p]](par[[p]]), 0))
}
# .ShiftMove() rescales the free parameters with a non-zero tempo exponent,
# p -> p exp(-e_p d) (lambda and mu share e = beta), and scores
# sum over those p of priors[[p]](new) - priors[[p]](old).  With the joint
# prior, the first rescaled FBD rate (the carrier) gets a prior that returns
# the joint density at the old state for the old value and at the rescaled
# state for any other value, and the other rescaled FBD rates get 0, so the
# sum is FbdLogPrior(new) - FbdLogPrior(old) exactly
ShiftJoint <- function(st, model, free, priors, nSlice, CharLL) {
  if (IsJoint(priors)) {
    expo <- c(rate = model[["gamma"]], lambda = model[["beta"]], mu = model[["beta"]],
              psi = model[["delta"]], theta = model[["alpha"]])
    comp <- intersect(names(expo)[expo != 0], free)
    jc <- intersect(comp, fbdNames)
    if (length(jc)) {
      old <- st[["par"]]
      carrier <- jc[1]
      priors[[carrier]] <- function(x) {
        if (identical(x, old[[carrier]])) return(FbdLogPrior(old))
        new <- old
        for (q in jc) new[[q]] <- old[[q]] * (x / old[[carrier]])^(expo[[q]] / expo[[carrier]])
        FbdLogPrior(new)
      }
      for (q in setdiff(jc, carrier)) priors[[q]] <- function(x) 0
    }
  }
  res <- .shiftOriginal(st, model, free, priors, nSlice, CharLL)
  if (IsJoint(priors) && all(fbdNames %in% free) && .cet$blockRounds > 0) {
    tm <- parent.frame()
    burn <- get("i", envir = tm, inherits = FALSE) <= get("burnin", envir = tm, inherits = FALSE)
    for (k in seq_len(.cet$blockRounds)) res$state <- FbdBlockMoves(res$state, model, nSlice, burn)
  }
  res
}

# FBD block moves (README, "Moves").  With turnover near 1 (the comparator's
# posterior: 0.9999), lambda and mu lie on a ridge lambda - mu = d << lambda,
# along which the engine's one-parameter moves cannot travel (check-prior.R
# without these moves: beta = 0 chains stuck, mean log lambda -1.1 against
# 2.5).  Three random-walk Metropolis moves, as the comparator's
# mvScale(diversification), mvSlide(turnover), mvSlide(sampling_prop), act on
# u = (log d, logit turnover, logit sampling proportion), with
#   lambda = d (1 + e^u2), mu = d e^u2, psi = mu e^u3,
# so each keeps the other two coordinates fixed.  The target in u is the
# paths-conditional density TempoMCMC()'s parameter sweep uses (augmented
# density at the current paths and bins; characters do not depend on these
# rates) times the prior in u: Exp(100) on d and flat turnover and sampling
# proportion, with the log / logit Jacobians.  This prior is written
# independently of FbdLogPrior(); the two agree only if FbdLogPrior()'s
# Jacobian is right (check-prior.R).  They run in the `shift` slot (once per
# iteration in every cell), after the shift move, `CET_FBD_BLOCK` rounds
# (default 1; 0 disables them); step sizes adapt during burn-in as the
# engine's (.Adapt(), target acceptance 0.3, every 20 calls) and are part of
# the checkpoint.
.cet$blockRounds <- as.integer(Sys.getenv("CET_FBD_BLOCK", "1"))
blockNames <- c("d", "turnover", "sampling")
ResetBlock <- function() {
  .cet$blockStep <- setNames(c(1, 1, 1), blockNames)
  .cet$blockHits <- setNames(numeric(3), blockNames)
  .cet$blockCalls <- 0L
  .cet$blockTry <- .cet$blockAcc <- setNames(numeric(3), blockNames)
}
ResetBlock()
FbdToU <- function(par) {
  d <- par[["lambda"]] - par[["mu"]]
  c(log(d), log(par[["mu"]]) - log(d), log(par[["psi"]]) - log(par[["mu"]]))
}
FbdFromU <- function(u) {
  d <- exp(u[1])
  mu <- d * exp(u[2])
  list(lambda = d + mu, mu = mu, psi = mu * exp(u[3]))
}
FbdUPrior <- function(u) {
  log(100) - 100 * exp(u[1]) + u[1] +
    stats::plogis(u[2], log.p = TRUE) + stats::plogis(-u[2], log.p = TRUE) +
    stats::plogis(u[3], log.p = TRUE) + stats::plogis(-u[3], log.p = TRUE)
}
FbdBlockMoves <- function(st, model, nSlice, burn) {
  tree <- st[["tree"]]
  cur <- get(".LogAugmented", .ns)(tree, st[["set"]], st[["paths"]], st[["bin"]])
  for (j in seq_along(blockNames)) {
    u <- FbdToU(st[["par"]])
    v <- u
    v[j] <- u[j] + .cet$blockStep[[j]] * (stats::runif(1) - 0.5)
    prop <- st[["par"]]
    prop[fbdNames] <- FbdFromU(v)
    .cet$blockTry[[j]] <- .cet$blockTry[[j]] + 1
    set <- tryCatch(get(".Precompute", .ns)(tree, prop, model),
                    infeasibleProposal = function(e) NULL)
    if (is.null(set)) next
    new <- tryCatch(get(".LogAugmented", .ns)(tree, set, st[["paths"]], st[["bin"]]),
                    infeasibleProposal = function(e) -Inf)
    logRatio <- new - cur + FbdUPrior(v) - FbdUPrior(u)
    if (is.finite(logRatio) && log(stats::runif(1)) < logRatio) {
      st[["par"]] <- prop
      st[["set"]] <- get(".AddSlices", .ns)(tree, set, nSlice)
      st[["msg"]] <- get(".Messages", .ns)(tree, st[["set"]])
      cur <- new
      .cet$blockAcc[[j]] <- .cet$blockAcc[[j]] + 1
      if (burn) .cet$blockHits[[j]] <- .cet$blockHits[[j]] + 1
    }
  }
  if (burn) {
    .cet$blockCalls <- .cet$blockCalls + 1L
    if (.cet$blockCalls %% 20L == 0L) {
      .cet$blockStep <- get(".Adapt", .ns)(.cet$blockStep, .cet$blockHits / 20)
      .cet$blockHits[] <- 0
    }
  }
  st
}

# ---- Recording and exact resume (as the penguin pipeline) -----------------------
RecordingStats <- function(tree, lengths = NULL, inverse = NULL) {
  pf <- parent.frame()
  if (is.null(lengths)) {
    r <- .cet$resume
    if (!is.null(r)) {
      for (v in c("st", "step", "jointStep", "ageStep", "moves")) {
        assign(v, r[[v]], envir = pf)
      }
      assign(".Random.seed", r[["seed"]], envir = globalenv())
      .cet$blockStep <- r[["blockStep"]]
      .cet$blockCalls <- r[["blockCalls"]]
      .cet$resume <- NULL
      .cet$injected <- TRUE
    }
  } else {
    st <- get("st", envir = pf, inherits = FALSE)
    model <- get("model", envir = pf, inherits = FALSE)
    chars <- get("characters", envir = pf, inherits = FALSE)
    tr <- st[["tree"]]
    if (is.null(chars)) {
      # No characters (check-prior.R): TempoMCMC() scores none, so no tree term
      charLL <- 0
      treeTerm <- 0
    } else {
      len <- get(".CharacterLengths", .ns)(tr, st[["set"]][["grid"]][["x"]], st[["paths"]],
                                            st[["par"]][["rate"]], model[["gamma"]])
      charLL <- HymMk(tr, chars, len, power = 1, treeTerm = FALSE)
      treeTerm <- TreeTerm(tr)
    }
    # The chain's stored value must decompose into the two (see 2. above)
    if (abs(st[["charLL"]] - (.cet$power * charLL + treeTerm)) > 1e-6 * (1 + abs(charLL))) {
      stop(sprintf("stored charLL %.10g != power x charLL + treeTerm = %.10g",
                   st[["charLL"]], .cet$power * charLL + treeTerm))
    }
    .cet$records[[length(.cet$records) + 1L]] <- list(
      newick = TempoNewick(tr), charLL = charLL, treeTerm = treeTerm,
      ages = c(holoAge = CladeAge(tr, .cet$holo),
               hymCrownAge = CladeAge(tr, .cet$hymCrown, mrca = TRUE),
               hymTotalAge = CladeAge(tr, .cet$hymAll, mrca = TRUE)),
      state = st[c("tree", "par", "paths", "bin")])
    .cet$captured <- list(st = st, step = get("step", pf), jointStep = get("jointStep", pf),
                          ageStep = get("ageStep", pf), moves = get("moves", pf),
                          seed = get(".Random.seed", globalenv()),
                          blockStep = .cet$blockStep, blockCalls = .cet$blockCalls)
  }
  .statsOriginal(tree, lengths, inverse)
}

utils::assignInNamespace("MkLikelihood", HymMk, ns = "EvoTempo")
utils::assignInNamespace(".TreeStats", RecordingStats, ns = "EvoTempo")
utils::assignInNamespace(".LogPrior", LogPriorJoint, ns = "EvoTempo")
utils::assignInNamespace(".ShiftMove", ShiftJoint, ns = "EvoTempo")

# Newick of a TempoTree, with branch lengths from ages; sampled ancestors are
# zero-length tips; a stem, if any, is the root edge
TempoNewick <- function(tt, digits = 6) {
  nTip <- tt[["nTip"]]
  age <- tt[["age"]]
  lab <- tt[["tipLabel"]]
  Sub <- function(v) {
    len <- if (is.na(tt[["parent"]][v])) NA else age[tt[["parent"]][v]] - age[v]
    s <- if (v <= nTip) {
      lab[v]
    } else {
      k <- tt[["kids"]][v - nTip, ]
      paste0("(", Sub(k[1]), ",", Sub(k[2]), ")")
    }
    if (is.na(len)) s else paste0(s, ":", formatC(max(len, 0), digits = digits, format = "f"))
  }
  stem <- tt[["stem"]]
  paste0(Sub(tt[["root"]]),
         if (!is.null(stem) && stem > 0) paste0(":", formatC(stem, digits = digits, format = "f")),
         ";")
}

# ---- Data ----------------------------------------------------------------------
hymDir <- "data-raw/empirical/hymenoptera"
local({
  cons <- read.csv(file.path(hymDir, "hymenoptera_constraints.csv"))
  tips <- read.csv(file.path(hymDir, "hymenoptera_tips.csv"))
  Set <- function(name) strsplit(cons$taxa[cons$constraint == name], " ")[[1]]
  .cet$holo <- Set("holometabola_withFossils")
  .cet$hymCrown <- Set("hymenoptera")
  .cet$hymAll <- tips$taxon[tips$group == "Hymenoptera"]
  stopifnot(length(.cet$holo) == 111L, length(.cet$hymCrown) == 60L, length(.cet$hymAll) == 105L)
})

LoadHymenoptera <- function(startTree = file.path(cetDir, "data", "start-tree.nwk")) {
  tips <- read.csv(file.path(hymDir, "hymenoptera_tips.csv"))
  tr <- read.tree(startTree)
  nTip <- Ntip(tr)
  info <- tips[match(tr$tip.label, tips$taxon), ]
  stopifnot(!anyNA(info$taxon))
  depth <- node.depth.edgelength(tr)
  extant <- info$status == "extant"
  age <- max(depth[which(extant)]) - depth
  age[which(extant)] <- 0
  len <- numeric(nTip + tr$Nnode)
  len[tr$edge[, 2]] <- tr$edge.length
  tipType <- ifelse(extant, "extant", ifelse(len[seq_len(nTip)] < 1e-9, "ancestor", "fossil"))
  ranges <- cbind(info$min_age, info$max_age)
  ranges[extant, ] <- 0
  fos <- which(tipType == "fossil")
  age[fos] <- pmin(pmax(age[fos], ranges[fos, 1]), ranges[fos, 2])
  tt <- TempoTree(tr, ages = age, tipType = tipType, ranges = ranges)
  chars <- read.csv(file.path(cetDir, "data", "characters.csv"), check.names = FALSE,
                    colClasses = "character", na.strings = character(0))
  parts <- read.csv(file.path(cetDir, "data", "partitions.csv"))
  cells <- as.matrix(chars[match(tr$tip.label, chars$taxon), -1])
  stopifnot(!anyNA(cells), ncol(cells) == nrow(parts))
  rownames(cells) <- tr$tip.label
  # TempoMCMC() needs a matrix with one row per tip; the rebound MkLikelihood()
  # scores the partitions in attribute "hym".  Unambiguous states, NA otherwise
  x <- matrix(suppressWarnings(as.integer(cells)), nTip, dimnames = dimnames(cells))
  attr(x, "hym") <- BuildPartitions(cells, parts)
  list(tree = tt, characters = x, cells = cells, partitions = parts, tips = info)
}

# ---- Model, priors, moves -------------------------------------------------------
cells <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
CellLabel <- function(beta, gamma) sprintf("b%sg%s", format(beta), format(gamma))

# Hymenoptera: crown start (the comparator conditions on the root age) with
# conditioning on a sampled extant descendant of each root daughter (RevBayes
# condition = "survival" with rho > 0); rho = 0.0005 as the comparator
CetModel <- function(beta, gamma,
                     nBin = as.integer(Sys.getenv("CET_NBIN", "26")),
                     zMax = as.numeric(Sys.getenv("CET_ZMAX", "5"))) {
  DCETModel(nBin = nBin, zRange = c(-zMax, zMax), beta = beta, gamma = gamma,
            rho = 5e-4, condition = "extant", start = "root")
}

LnTrunc <- function(meanlog, sdlog, zTrunc = 4) {
  logMass <- log(2 * stats::pnorm(zTrunc) - 1)
  force(meanlog)
  force(sdlog)
  function(x) {
    if (!(x > 0) || abs(log(x) - meanlog) > zTrunc * sdlog) return(-Inf)
    stats::dlnorm(x, meanlog, sdlog, log = TRUE) - logMass
  }
}
# Exponential prior truncated to (0, hi], renormalised
ExpTrunc <- function(rate, hi) {
  logMass <- stats::pexp(hi, rate, log.p = TRUE)
  function(x) {
    if (!(x > 0 && x <= hi)) return(-Inf)
    stats::dexp(x, rate, log = TRUE) - logMass
  }
}
# Rates per Myr.  lambda, mu, psi: the comparator's prior on (d, turnover,
# sampling proportion), through FbdConditional().  theta: reversion time
# 1/theta with median 120 Myr (one third of the comparator's posterior root
# age, the ratio of the penguin prior's 20 Myr to the penguin root of ~60 Ma),
# 95% 17-850 Myr.  s: the comparator's UCLN log-sd prior Exp(mean 0.587405),
# truncated at 3 (numerical bound of the penguin pipeline; removes 0.61% of
# the mass).  rate: the comparator's clock prior LogNormal(-7.1, 0.5),
# truncated at +-4 sdlog.
CetPriors <- function() {
  list(lambda = FbdConditional("lambda"), mu = FbdConditional("mu"),
       psi = FbdConditional("psi"), theta = LnTrunc(log(1 / 120), 1),
       s = ExpTrunc(1 / 0.587405, 3), rate = LnTrunc(-7.1, 0.5))
}
# Log prior of all parameters (for outputs outside TempoMCMC())
CetLogPrior <- function(par, priors = CetPriors()) LogPriorJoint(par, priors)
# Root age ~ offset exponential, minimum 315 Ma, mean 396 Ma (the comparator's
# root_time ~ dnExponential(1/81, offset = 315))
CetAgePriors <- function() {
  list(root = function(a) if (a > 315) -log(81) - (a - 315) / 81 else -Inf)
}
# Near the comparator's posterior medians (parameters_morph.csv: lambda 42.48,
# mu 42.48, d 0.0055, psi 2.18e-6, clock_rate 0.0019); theta at its prior
# median, s 0.3.  Replicates > 1 are jittered by exp(N(0, 0.2^2)) on d,
# lambda, psi, theta, s and rate, keeping mu = lambda - d
CetStart <- function(rep) {
  d <- 0.0055
  start <- list(lambda = 42.48, psi = 2.18e-6, theta = 1 / 120, s = 0.3, rate = 0.0019)
  if (rep > 1) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
    set.seed(9000 + rep)
    start <- lapply(start, function(v) v * exp(stats::rnorm(1, 0, 0.2)))
    d <- d * exp(stats::rnorm(1, 0, 0.2))
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  }
  c(start["lambda"], list(mu = start$lambda - d), start[c("psi", "theta", "s", "rate")])
}
CetMoves <- function(tt, characters, model,
                     scale = as.numeric(Sys.getenv("CET_TREE_MOVES", "1"))) {
  nBif <- sum(!tt$sampled)
  nFossil <- sum(tt$tipType != "extant")
  c(TempoMoves(tt, characters, model),
    nodeAge = ceiling(scale * 0.25 * nBif), fossilAge = ceiling(scale * 0.25 * nFossil),
    slide = ceiling(scale * 0.15 * nBif), saJump = ceiling(scale * 0.1 * nFossil),
    rootAge = 1)
}

# ---- Utilities -------------------------------------------------------------------
LogSumExp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}
LogMeanExp <- function(x) LogSumExp(x) - log(length(x))
# Geyer initial monotone sequence ESS (for the bridge-sampling errors, as the
# penguin pipeline; convergence summaries use MkPrime's rank-normalised R-hat
# and ESS)
Ess <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4 || stats::var(x) == 0) return(NA_real_)
  ac <- stats::acf(x, lag.max = n - 1, plot = FALSE, demean = TRUE)$acf[, 1, 1]
  nPair <- length(ac) %/% 2
  pairSum <- ac[2 * seq_len(nPair) - 1] + ac[2 * seq_len(nPair)]
  stop1 <- which(pairSum <= 0)[1]
  if (!is.na(stop1)) pairSum <- pairSum[seq_len(stop1 - 1)]
  tau <- if (length(pairSum)) -1 + 2 * sum(cummin(pairSum)) else 1
  min(n, n / max(tau, 1e-8))
}
