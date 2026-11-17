# Hymenoptera CET pipeline: definitions shared by run-chain.R, check-mkv.R,
# bayes-factors.R, check-prior.R and the check scripts.  Source from the
# package root after loading EvoTempo.  Follows
# data-raw/empirical/penguins/cet/common.R; the differences are marked
# "Hymenoptera:".
#
# The characters use the package's model (#97, #98): MkCharacters() with the
# number of states and the ordering of each character, ambiguity sets, the
# Mkv correction per missing-data mask (coding = "variable", as RevBayes
# dnPhyloCTMC(coding = "variable")) and Gamma(4) rates among characters, whose
# shape alpha is a parameter of TempoMCMC() (start$alpha) with the
# comparator's prior Exp(1) (CetCharPriors()).  This is the comparator's
# character model (../revbayes/hymenoptera_fbd.Rev).
#
# What the engine still lacks, this file adds by rebinding functions inside
# the EvoTempo namespace for the duration of a run, without changing the
# package:
#
# 1. Power and tree term.  The rebound MkLikelihood() (CetMk()) calls the
#    package's MkLikelihood() and raises the character likelihood, Mkv
#    correction included, to a power (`.cet$power`) for stepping-stone runs.
# 2. Hymenoptera: the Holometabola constraint and calibration of the
#    comparator (a clade constraint and an offset-exponential density on the
#    clade's age) are tree-only terms with no place in TempoMCMC().  CetMk()
#    adds them, unpowered, to its return value (TreeTerm()).  This is exact:
#    TempoMCMC() rescores the characters on every proposed tree (tree and
#    start-age moves) and reuses the stored value only when the tree is
#    unchanged (path, shift, parameter and alpha moves; in the alpha move the
#    term cancels).  The recorded charLL excludes the term, and
#    RecordingStats() checks the decomposition at every sample.
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
#    pipeline) records them and continues a captured chain exactly (step
#    sizes of the parameter, joint, start-age and alpha moves and of the FBD
#    block moves included).
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
stopifnot("EvoTempo lacks MkCharacters() (#97): install the gamma-acrv-97 branch or later" =
            exists("MkCharacters", .ns))
cetDir <- "data-raw/empirical/hymenoptera/cet"

# ---- Characters -----------------------------------------------------------------
# Number of discrete Gamma rate categories among characters, as the comparator
# (fnDiscretizeGamma(alpha_morph, alpha_morph, 4)).  CET_NRATECAT = 1 (no rate
# variation; alpha is then not a parameter) is a diagnostic variant
# (README, "Root ages"); run-chain.R records it in the checkpoint.
.cet$nRateCat <- as.integer(Sys.getenv("CET_NRATECAT", "4"))
stopifnot("CET_NRATECAT must be a whole number >= 1" = isTRUE(.cet$nRateCat >= 1L))
# Holometabola calibration on (default) or off (CET_HOLO_CAL = false: the
# constraint only), a diagnostic variant for prior-only comparisons with
# RevBayes, whose ignoreAllData() also ignores the clamped calibration node
.cet$holoCal <- Sys.getenv("CET_HOLO_CAL", "true") == "true"
# Coding of the Mkv variant: "pattern" is RevBayes coding = "variable" (one
# correction per missing-data mask, #65, #80); "none" is coding = "all"
MkvCoding <- function(mkv = .cet$mkv) switch(mkv, pattern = "variable", none = "all",
                                              stop("`mkv` must be \"pattern\" or \"none\""))
# Characters of LoadHymenoptera(): `cells` (taxon x character, row names the
# tip labels; a state, a set such as "{01}", or "" for unknown) and
# `partitions` (data/partitions.csv: nState, ordered per character)
CetCharacters <- function(cells, partitions, mkv = .cet$mkv, nRateCat = .cet$nRateCat) {
  stopifnot(!is.null(rownames(cells)), ncol(cells) == nrow(partitions))
  MkCharacters(cells, nState = partitions$nState, ordered = partitions$ordered,
               coding = MkvCoding(mkv), nRateCat = nRateCat)
}
# The rebound MkLikelihood(): power x the package's character log-likelihood
# (Mkv correction and Gamma mixture included) + the unpowered tree term.
# Input other than MkCharacters() goes to the package's function unchanged.
CetMk <- function(tree, characters, lengths, nState = 2L,
                  engine = getOption("EvoTempo.engine", "cpp"), alpha = NULL,
                  power = .cet$power, treeTerm = TRUE) {
  if (!inherits(characters, "MkCharacters")) {
    return(.mkOriginal(tree, characters, lengths, nState, engine, alpha))
  }
  tree <- TempoTree(tree)
  tt <- if (treeTerm) TreeTerm(tree) else 0
  if (power == 0 || !is.finite(tt)) {
    return(tt)
  }
  power * .mkOriginal(tree, characters, lengths, nState, engine, alpha) + tt
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
TreeTerm <- function(tree, calibration = .cet$holoCal) {
  a <- CladeAge(tree, .cet$holo)
  if (is.na(a)) -Inf else if (calibration) HoloCalibration(a) else 0
}

# ---- FBD-rate prior of the comparator, as conditional priors --------------------
# log density of (lambda, mu, psi) under d = lambda - mu ~ Exp(100),
# r = mu / lambda ~ U(0, 1), s = psi / (mu + psi) ~ U(0, 1):
#   p(lambda, mu, psi) = 100 exp(-100 d) |d(d, r, s) / d(lambda, mu, psi)|,
#   |J| = (lambda - mu) mu / (lambda^2 (mu + psi)^2).
FbdLogPrior <- function(par, beta = 0) {
  l <- par[["lambda"]]
  m <- par[["mu"]]
  p <- par[["psi"]]
  if (!(is.finite(l) && is.finite(m) && is.finite(p) && m > 0 && l > m && p > 0)) {
    return(-Inf)
  }
  # Lineage-averaged rates (CET_FBD_SCALE = "lineage"; README, "Cell (1, 0)"):
  # the density of (lambda, mu, psi) is that of (k lambda, k mu, psi) under the
  # comparator's prior, times k^2 (the Jacobian), with k = LineageMean(s, beta)
  k <- if (beta != 0 && .cet$fbdScale == "lineage") LineageMean(par[["s"]], beta) else 1
  l <- k * l
  m <- k * m
  d <- l - m
  log(100) - 100 * d + log(d) + log(m) - 2 * log(l) - 2 * log(m + p) + 2 * log(k)
}
# ---- Lineage-averaged rates and the bound on s (diagnosis of cell (1, 0)) --------
# CET_FBD_SCALE = "none" (default): the comparator's prior is on lambda, mu, the
# rates at log-tempo x = 0.  "lineage": it is on the lineage-averaged rates
# k lambda, k mu, k = E_pi[e^(beta x)], pi the stationary law of the tempo
# generator Q = diag(e^(alpha x)) A (mean x = -s^2, k about e^(-s^2 / 2) for
# alpha = beta = 1; diag-ridge.R); s then enters the joint prior.  k = 1 when
# beta = 0, so cells with beta = 0 are unchanged.  CET_S_MAX (default 3) is the
# upper bound of the prior on s.
.cet$fbdScale <- Sys.getenv("CET_FBD_SCALE", "none")
stopifnot("CET_FBD_SCALE must be \"none\" or \"lineage\"" = .cet$fbdScale %in% c("none", "lineage"))
.cet$sMax <- as.numeric(Sys.getenv("CET_S_MAX", "3"))
stopifnot("CET_S_MAX must be a positive number" = isTRUE(.cet$sMax > 0))
LineageMean <- function(s, beta, nBin = as.integer(Sys.getenv("CET_NBIN", "26")),
                        zMax = as.numeric(Sys.getenv("CET_ZMAX", "5"))) {
  if (!is.finite(s) || s <= 0) return(NA_real_)
  grid <- TempoGrid(nBin, c(-zMax, zMax), s)
  p <- TempoEquilibrium(TempoGenerator(grid, 1, 1, "central"))$p
  sum(p * exp(beta * grid$x))
}
fbdNames <- c("lambda", "mu", "psi")
# Prior for parameter `p` in TempoMCMC()'s `priors`: the joint density at the
# caller's current parameters with `p` set to x.  The current parameters are
# `par` in .ParameterSweep() and `st$par` in .JointSweep(); in both, every
# call changes `p` only.  Any other caller is an error (.ShiftMove() and
# .LogPrior() never call these functions: see ShiftJoint(), LogPriorJoint())
FbdConditional <- function(p, JointPrior = FbdLogPrior) {
  force(p)
  force(JointPrior)
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
    JointPrior(cur)
  }
  attr(f, "cetJoint") <- p
  attr(f, "cetJointPrior") <- JointPrior
  f
}
IsJoint <- function(priors) {
  any(vapply(priors, function(f) !is.null(attr(f, "cetJoint")), TRUE))
}
JointNames <- function(priors) {
  names(priors)[vapply(priors, function(f) !is.null(attr(f, "cetJoint")), TRUE)]
}
JointPriorOf <- function(priors) attr(priors[[JointNames(priors)[1]]], "cetJointPrior")
# .LogPrior() for logPost: the joint density (FBD rates, and s with
# CET_FBD_SCALE = "lineage") plus the other priors
LogPriorJoint <- function(par, priors) {
  if (!IsJoint(priors)) return(.logPriorOriginal(par, priors))
  stopifnot(all(fbdNames %in% names(priors)))
  joint <- JointNames(priors)
  other <- setdiff(names(priors), joint)
  JointPriorOf(priors)(par) + sum(vapply(other, function(p) priors[[p]](par[[p]]), 0))
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
    jc <- intersect(comp, JointNames(priors))
    if (length(jc)) {
      old <- st[["par"]]
      carrier <- jc[1]
      JP <- JointPriorOf(priors)
      priors[[carrier]] <- function(x) {
        if (identical(x, old[[carrier]])) return(JP(old))
        new <- old
        for (q in jc) new[[q]] <- old[[q]] * (x / old[[carrier]])^(expo[[q]] / expo[[carrier]])
        JP(new)
      }
      for (q in setdiff(jc, carrier)) priors[[q]] <- function(x) 0
    }
  }
  res <- .shiftOriginal(st, model, free, priors, nSlice, CharLL)
  if (IsJoint(priors) && all(fbdNames %in% free) && .cet$blockRounds > 0) {
    tm <- parent.frame()
    burn <- get("i", envir = tm, inherits = FALSE) <= get("burnin", envir = tm, inherits = FALSE)
    # The comparator's prior on (lambda, mu, psi) has its own form in u
    # (FbdUPrior()); a joint prior that includes s (CET_FBD_SCALE = "lineage")
    # goes through the general target
    JP <- if ("s" %in% JointNames(priors)) JointPriorOf(priors)
    for (k in seq_len(.cet$blockRounds)) {
      res$state <- FbdBlockMoves(res$state, model, nSlice, burn, JP)
    }
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
# With a joint prior JP other than the comparator's on (lambda, mu, psi)
# (CET_FBD_SCALE = "lineage"), the target in u is JP + log|d(lambda, mu, psi) /
# du|, and log|d(lambda, mu, psi) / du| = FbdUPrior(u) - FbdLogPrior(par(u))
# (beta = 0: the comparator's density on (lambda, mu, psi) is its density in u
# divided by that Jacobian); with the comparator's prior this is FbdUPrior(u)
FbdBlockMoves <- function(st, model, nSlice, burn, JP = NULL) {
  UTarget <- function(u, par) {
    if (is.null(JP)) return(FbdUPrior(u))
    JP(par) + FbdUPrior(u) - FbdLogPrior(par)
  }
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
    logRatio <- new - cur + UTarget(v, prop) - UTarget(u, st[["par"]])
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
# Loop variables of TempoMCMC() that a resumed chain must continue with: step
# sizes of the parameter sweep, joint sweep, start-age moves and (with Gamma
# rates) the alpha move, and the move counts after burn-in; plus the chain
# state and the random-number state.  The FBD block-move step sizes are in .cet.
resumeVars <- c("st", "step", "jointStep", "ageStep", "alphaStep", "moves")
RecordingStats <- function(tree, lengths = NULL, inverse = NULL) {
  pf <- parent.frame()
  if (is.null(lengths)) {
    r <- .cet$resume
    if (!is.null(r)) {
      # TempoMCMC() calls .TreeStats(tree) once, after setting up these
      # variables and before its first iteration
      stopifnot(all(resumeVars %in% names(r)), all(resumeVars %in% ls(pf)))
      for (v in resumeVars) {
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
      charLL <- CetMk(tr, chars, len, alpha = st[["par"]][["alpha"]], power = 1,
                      treeTerm = FALSE)
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
    .cet$captured <- c(mget(resumeVars, envir = pf, inherits = FALSE),
                       list(seed = get(".Random.seed", globalenv()),
                            blockStep = .cet$blockStep, blockCalls = .cet$blockCalls))
  }
  .statsOriginal(tree, lengths, inverse)
}

utils::assignInNamespace("MkLikelihood", CetMk, ns = "EvoTempo")
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

# Reference start tree (prepare.R); production chains start from their own
# draws (draw-start-trees.R; StartTreeFor())
referenceTree <- file.path(cetDir, "data", "start-tree.nwk")
startTreeDir <- file.path(cetDir, "data", "start-trees")
LoadHymenoptera <- function(startTree = referenceTree, mkv = .cet$mkv,
                            nRateCat = .cet$nRateCat) {
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
  # One row per tip, in the tree's tip order, with the tip labels as row names
  # (MkCharacters() keeps them and they must equal the tree's tip labels)
  cells <- as.matrix(chars[match(tr$tip.label, chars$taxon), -1])
  stopifnot(!anyNA(cells), ncol(cells) == nrow(parts))
  rownames(cells) <- tr$tip.label
  x <- CetCharacters(cells, parts, mkv, nRateCat)
  stopifnot(identical(x$tipLabel, tt$tipLabel))
  list(tree = tt, characters = x, cells = cells, partitions = parts, tips = info)
}

# ---- Model, priors, moves -------------------------------------------------------
# Chain labels, the expected chain set, segment caps and eval freshness: as the
# penguin pipeline (#75, #101, #103, #104, #107 review)
cells <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
# A number as it appears in chain labels and cross.csv column names: plain
# decimal notation (never scientific, which ParseLabel() cannot read), with
# the fewest significant digits that read back as exactly `x` (#75)
FormatValue <- function(x) {
  stopifnot(length(x) == 1L, is.finite(x), x >= 0)
  for (d in 1:17) {
    s <- format(x, scientific = FALSE, digits = d, trim = TRUE)
    if (as.numeric(s) == x) return(s)
  }
  stop("cannot write ", x, " in plain decimal notation")
}
CellLabel <- function(beta, gamma) sprintf("b%sg%s", FormatValue(beta), FormatValue(gamma))

# Chain labels (run-chain.R): b<beta>g<gamma>[_mkv<variant>][_t<power>]_run<rep>.
# ParseLabel() stops on a label it cannot read (#75); before, it returned NULL,
# and Chains() dropped the directory without a message.
ParseLabel <- function(label) {
  m <- regmatches(label, regexec(
    "^b([0-9.]+)g([0-9.]+)(_mkv(none))?(_t([0-9.]+))?_run([0-9]+)$", label))[[1]]
  out <- if (length(m)) {
    data.frame(label = label, beta = as.numeric(m[2]), gamma = as.numeric(m[3]),
               mkv = if (nzchar(m[5])) m[5] else "pattern",
               power = if (nzchar(m[7])) as.numeric(m[7]) else 1, rep = as.integer(m[8]))
  }
  if (is.null(out) || anyNA(out[, c("beta", "gamma", "power", "rep")])) {
    stop(sprintf(paste("cannot parse chain label '%s': expected",
                       "b<beta>g<gamma>[_mkv<variant>][_t<power>]_run<rep>",
                       "with numbers in plain decimal notation"), label))
  }
  out
}
# Chain directories under `root`: every subdirectory holding run-chain.R
# output (seg*_samples.csv, checkpoint.rds or progress.txt) must have a
# readable label; other subdirectories (e.g. summary/) are ignored.  Only
# chains with at least one complete segment are returned.
Chains <- function(root) {
  dirs <- list.dirs(root, recursive = FALSE)
  isChain <- vapply(dirs, function(d) {
    length(Sys.glob(file.path(d, "seg*_samples.csv"))) > 0L ||
      file.exists(file.path(d, "checkpoint.rds")) || file.exists(file.path(d, "progress.txt"))
  }, TRUE)
  out <- do.call(rbind, lapply(dirs[isChain], function(d) {
    p <- ParseLabel(basename(d))
    if (!length(Sys.glob(file.path(d, "seg*_samples.csv")))) return(NULL)
    cbind(p, dir = d)
  }))
  if (is.null(out)) {
    out <- data.frame(label = character(), beta = numeric(), gamma = numeric(),
                      mkv = character(), power = numeric(), rep = integer(), dir = character())
  }
  out
}

# ---- Expected chain set (#101) ----------------------------------------------------
# `bayes-factors.R combine` must use every chain that was launched, not just the
# chains it happens to find with a cross.csv: eval writes cross.csv only after
# its checks pass, so a chain whose eval failed (or was never run) used to drop
# out of the Bayes factors without a message (found in the trilobite pipeline,
# #91 delta review, N2; this pipeline had the same gap).
#
# The expected set comes from the launch record <run root>/launch.csv, one row
# per launched group of chains (cells x replicates x rungs):
#   cell   the (beta, gamma) label as run-chain.R writes it: b0g0, b1g1, b0.5g0
#   reps   number of replicates; the group is runs 1..reps
#   power  character power (optional; default 1)
#   mkv    Mkv variant (optional; default "pattern"; or "none")
# For example, the production launch (run-array.sh --array=1-16) plus one rung
# of 4 replicates (CET_CELL=b0.5g0, --array=1-4):
#   cell,reps
#   b0g0,4
#   b0g1,4
#   b1g0,4
#   b1g1,4
#   b0.5g0,4
# A rung launched with its own replicate count, a power ladder (power = t) or
# another Mkv variant (mkv = none) gets its own row.
# Without launch.csv the expected set is the production launch alone: the four
# cells x CET_BF_REPS (default 4) replicates at power 1.  `reps` must be a
# whole number >= 1 (2.7 used to be truncated to 2), and a record that selects
# no chain of the CET_MKV variant stops (#103 review, R103-8, R103-10).
ExpectedChains <- function(root, mkv = .cet$mkv) {
  f <- file.path(root, "launch.csv")
  if (file.exists(f)) {
    L <- utils::read.csv(f, colClasses = "character", strip.white = TRUE)
    if (!all(c("cell", "reps") %in% names(L)) || !nrow(L)) {
      stop(f, " needs columns cell and reps (optional: power, mkv) and at least one row")
    }
    source <- f
  } else {
    L <- data.frame(cell = mapply(CellLabel, cells$beta, cells$gamma),
                    reps = Sys.getenv("CET_BF_REPS", "4"))
    source <- sprintf("no launch.csv: production default, 4 cells x %s replicates", L$reps[1])
  }
  if (is.null(L$power)) L$power <- "1"
  if (is.null(L$mkv)) L$mkv <- "pattern"
  L$power[is.na(L$power) | !nzchar(L$power)] <- "1"
  L$mkv[is.na(L$mkv) | !nzchar(L$mkv)] <- "pattern"
  out <- do.call(rbind, lapply(seq_len(nrow(L)), function(i) {
    reps <- if (grepl("^[0-9]+$", trimws(L$reps[i]))) as.integer(L$reps[i]) else NA_integer_
    power <- suppressWarnings(as.numeric(L$power[i]))
    if (is.na(reps) || reps < 1L || is.na(power)) {
      stop(sprintf("%s row %d: reps must be a whole number >= 1 (got '%s') and power a number",
                   source, i, L$reps[i]))
    }
    cell <- ParseLabel(paste0(L$cell[i], "_run1"))
    label <- paste0(CellLabel(cell$beta, cell$gamma), MkvSuffix(L$mkv[i]),
                    if (power != 1) paste0("_t", FormatValue(power)), "_run", seq_len(reps))
    do.call(rbind, lapply(label, ParseLabel))
  }))
  out <- out[out$mkv == mkv, ]
  if (!nrow(out)) stop(sprintf("%s: no expected chain of Mkv variant %s (CET_MKV)", source, mkv))
  if (anyDuplicated(ChainKey(out))) stop(source, " lists a chain twice")
  attr(out, "source") <- source
  out
}
ChainKey <- function(x) {
  if (!nrow(x)) return(character())
  paste(vapply(x$beta, FormatValue, ""), vapply(x$gamma, FormatValue, ""), x$mkv,
        vapply(x$power, FormatValue, ""), x$rep)
}
# Every chain directory under `root` (run-chain.R output: seg*_samples.csv,
# checkpoint.rds or progress.txt), with or without a complete segment, so that a
# launched chain missing from the record is flagged even before it has written
# a segment (R103-9)
ChainDirs <- function(root) {
  dirs <- list.dirs(root, recursive = FALSE)
  isChain <- vapply(dirs, function(d) {
    length(Sys.glob(file.path(d, "seg*_samples.csv"))) > 0L ||
      file.exists(file.path(d, "checkpoint.rds")) || file.exists(file.path(d, "progress.txt"))
  }, TRUE)
  out <- do.call(rbind, lapply(dirs[isChain], function(d) cbind(ParseLabel(basename(d)), dir = d)))
  if (is.null(out)) {
    out <- data.frame(label = character(), beta = numeric(), gamma = numeric(),
                      mkv = character(), power = numeric(), rep = integer(), dir = character())
  }
  out
}
# Segment cap (#107 review, R107-9): so that a summary can be reproduced while
# its chains are extended, CET_MAX_SEG = m > 0 makes eval, combine and the
# fingerprint read only segments 1..m of each chain, and CET_MAX_SEG_RUNG
# (default CET_MAX_SEG) does the same for chains other than the four
# production cells at power 1 (ladder rungs, power ladders).  For example,
# the penguin results of #107 read CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20.
# 0 (default) reads every segment.
SegCap <- function(label) {
  m <- suppressWarnings(as.integer(Sys.getenv("CET_MAX_SEG", "0")))
  if (is.na(m)) stop("CET_MAX_SEG must be a whole number")
  p <- tryCatch(ParseLabel(label), error = function(e) NULL)
  production <- !is.null(p) && p$power == 1 && p$beta %in% c(0, 1) && p$gamma %in% c(0, 1)
  if (!production && nzchar(Sys.getenv("CET_MAX_SEG_RUNG"))) {
    m <- suppressWarnings(as.integer(Sys.getenv("CET_MAX_SEG_RUNG")))
    if (is.na(m)) stop("CET_MAX_SEG_RUNG must be a whole number")
  }
  m
}
SegFiles <- function(dir, what = "samples.csv") {
  f <- sort(Sys.glob(file.path(dir, paste0("seg*_", what))))
  m <- SegCap(basename(dir))
  if (m > 0L) f <- utils::head(f, m)
  f
}
# A chain is complete when the last line of its progress.txt says so
# (run-chain.R writes "<time> complete" after the last segment, and a resumed
# or extended chain writes later lines after any earlier "complete")
ChainComplete <- function(dir) {
  # Under a segment cap, complete means the capped set is full
  m <- SegCap(basename(dir))
  if (m > 0L) return(length(SegFiles(dir)) >= m)
  f <- file.path(dir, "progress.txt")
  if (!file.exists(f)) return(FALSE)
  x <- trimws(readLines(f, warn = FALSE))
  x <- x[nzchar(x)]
  length(x) > 0L && grepl(" complete$", x[length(x)])
}
# What eval read (R103-2): the number of sample rows and an md5 over the
# chain's seg*_samples.csv files (not the states files, which may be deleted
# once eval has run).  eval records both in cross.csv (evalSamples,
# evalSegMd5), and combine refuses a cross.csv whose chain has changed since.
SegmentFingerprint <- function(dir) {
  f <- SegFiles(dir)
  n <- sum(vapply(f, function(x) nrow(utils::read.csv(x)), 0L))
  tmp <- tempfile()
  on.exit(unlink(tmp))
  writeLines(paste(basename(f), unname(tools::md5sum(f))), tmp)
  # "md5:" keeps read.csv from reading an all-digit hash as a number
  list(samples = n, md5 = paste0("md5:", unname(tools::md5sum(tmp))))
}
# Problems with the cross.csv files combine is about to use (`cross`, named by
# chain label; `ch`, their rows of Chains()): stale (eval read other sample
# files than the chain has now, or recorded nothing; R103-2), and, when `tol`
# is given, eval's checks: `check` and `checkChar` must be within the stricter
# of `tol` and the tolerance eval recorded (evalTol) (review of #91, A2;
# R103-2).  Returns one line per bad chain.
CrossProblems <- function(ch, cross, tol = NULL) {
  unlist(lapply(seq_along(cross), function(i) {
    x <- cross[[i]]
    lab <- names(cross)[i]
    if (!all(c("evalSamples", "evalSegMd5") %in% names(x))) {
      return(sprintf("%s: cross.csv does not record what eval read (written before #103); rerun eval", lab))
    }
    fp <- SegmentFingerprint(ch$dir[ch$label == lab])
    if (x$evalSamples[1] != fp$samples || x$evalSegMd5[1] != fp$md5) {
      return(sprintf("%s: cross.csv is stale: eval read %d samples (md5 %s), the chain now has %d (md5 %s); rerun eval",
                     lab, x$evalSamples[1], x$evalSegMd5[1], fp$samples, fp$md5))
    }
    if (is.null(tol)) return(NULL)
    cols <- intersect(c("check", "checkChar"), names(x))
    if (length(cols) < 2L) return(sprintf("%s: cross.csv has no checkChar column; rerun eval", lab))
    use <- min(tol, if (!is.null(x$evalTol)) x$evalTol[1] else tol)
    v <- as.matrix(x[, cols])
    if (!all(is.finite(v)) || !all(abs(v) <= use)) {
      return(sprintf("%s: cross.csv fails eval's checks (max |check|, |checkChar| %.3g > tolerance %g)",
                     lab, max(abs(v)), use))
    }
    NULL
  }))
}
# Compare the chains under `root` (`found`: Chains(), chains with a complete
# segment) with the expected set.  For the CET_MKV variant, combine stops if:
# - an expected chain has no complete segment;
# - an expected chain at power 1 (the chains the bridges use) has no cross.csv;
# - an expected chain is not complete (last progress.txt line; R103-2);
# - a chain directory, with or without a segment, is not in the record (R103-9);
# - two directories hold the same chain (e.g. b0g0_run1 and b0g0_run01; R103-7).
# The opt-outs are explicit and per chain: CET_BF_ALLOW_MISSING (labels such as
# b0g0_run3) lets combine leave a missing chain out, everywhere: bridges and
# stepping-stone; CET_BF_ALLOW_INCOMPLETE lets it use a chain that is not
# complete.  Each named chain must be expected.  Returns the labels left out,
# with attribute "incomplete" (the incomplete chains used on request).
CheckChainSet <- function(root, found, mkv = .cet$mkv,
                          allow = strsplit(Sys.getenv("CET_BF_ALLOW_MISSING"), ",")[[1]],
                          allowIncomplete = strsplit(Sys.getenv("CET_BF_ALLOW_INCOMPLETE"), ",")[[1]]) {
  allow <- trimws(allow[nzchar(trimws(allow))])
  allowIncomplete <- trimws(allowIncomplete[nzchar(trimws(allowIncomplete))])
  expected <- ExpectedChains(root, mkv)
  found <- found[found$mkv == mkv, ]
  dirs <- ChainDirs(root)
  dirs <- dirs[dirs$mkv == mkv, ]
  dKey <- ChainKey(dirs)
  fKey <- ChainKey(found)
  eKey <- ChainKey(expected)
  hasCross <- file.exists(file.path(found$dir, "cross.csv"))
  absent <- expected$label[!eKey %in% fKey]
  noCross <- found$label[fKey %in% eKey & found$power == 1 & !hasCross]
  unknown <- dirs$label[!dKey %in% eKey]
  unknown <- paste0(unknown, ifelse(dirs$label[!dKey %in% eKey] %in% found$label, "",
                                    " (no complete segment yet)"))
  dup <- dirs$label[dKey %in% dKey[duplicated(dKey)]]
  missing <- c(absent, noCross)
  leftOut <- intersect(missing, allow)
  # Completeness matters for every chain that will be used (not for one left out)
  incomplete <- setdiff(found$label[fKey %in% eKey & !vapply(found$dir, ChainComplete, TRUE)], leftOut)
  badAllow <- setdiff(c(allow, allowIncomplete), expected$label)
  cat(sprintf(paste("Expected chains (mkv %s): %d (%s); found %d; missing %d; not complete %d;",
                    "not in the launch record %d\n"),
              mkv, nrow(expected), attr(expected, "source"), nrow(found), length(missing),
              length(incomplete), length(unknown)))
  List <- function(what, x) if (length(x)) paste0(what, ": ", paste(x, collapse = ", "))
  problems <- c(
    List("no complete segment (chain directory missing or empty)", setdiff(absent, allow)),
    List("no cross.csv (eval failed or was not run)", setdiff(noCross, allow)),
    List("not complete (last line of progress.txt is not 'complete'; CET_BF_ALLOW_INCOMPLETE to use anyway)",
         setdiff(incomplete, allowIncomplete)),
    List("found but not in the launch record (add to <run root>/launch.csv)", unknown),
    List("two directories hold the same chain", dup),
    List("CET_BF_ALLOW_MISSING / CET_BF_ALLOW_INCOMPLETE name chains that are not expected", badAllow))
  if (length(problems)) {
    stop("combine: the chain set does not match the launch record; no Bayes factor is reported.\n  ",
         paste(problems, collapse = "\n  "),
         "\n  To leave a chain out on purpose, name it in CET_BF_ALLOW_MISSING.", call. = FALSE)
  }
  if (length(leftOut)) {
    cat(sprintf("Left out on request (CET_BF_ALLOW_MISSING): %s (from bridges and stepping-stone)\n",
                paste(leftOut, collapse = ", ")))
  }
  usedIncomplete <- intersect(incomplete, allowIncomplete)
  if (length(usedIncomplete)) {
    cat(sprintf("Used although not complete (CET_BF_ALLOW_INCOMPLETE): %s\n",
                paste(usedIncomplete, collapse = ", ")))
  }
  structure(leftOut, incomplete = usedIncomplete)
}
# chains-used.csv (R103-12): every chain of the CET_MKV variant combine found
# or expected, with its role in the outputs beside it
ChainsUsed <- function(chAll, leftOut, cross, mkv = .cet$mkv) {
  x <- chAll[chAll$mkv == mkv, c("label", "beta", "gamma", "mkv", "power", "rep")]
  if (length(leftOut)) {
    x <- rbind(x, do.call(rbind, lapply(leftOut, ParseLabel))[, names(x)])
    x <- x[!duplicated(x$label), ]
  }
  x$role <- ifelse(x$label %in% leftOut, "left out (CET_BF_ALLOW_MISSING)",
                   ifelse(x$power == 1, "bridge", "stepping-stone"))
  x$incompleteAllowed <- x$label %in% attr(leftOut, "incomplete")
  x$evalSamples <- vapply(x$label, function(l) if (is.null(cross[[l]])) NA_real_ else
    as.numeric(cross[[l]]$evalSamples[1]), 0)
  x$evalSegMd5 <- vapply(x$label, function(l) if (is.null(cross[[l]]$evalSegMd5)) NA_character_ else
    as.character(cross[[l]]$evalSegMd5[1]), "")
  x
}

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
CetPriors <- function(beta = 0) {
  sPrior <- ExpTrunc(1 / 0.587405, .cet$sMax)
  scaled <- beta != 0 && .cet$fbdScale == "lineage"
  force(beta)
  JP <- if (scaled) {
    function(par) FbdLogPrior(par, beta) + sPrior(par[["s"]])
  } else {
    FbdLogPrior
  }
  pri <- list(lambda = NULL, mu = NULL, psi = NULL, theta = LnTrunc(log(1 / 120), 1),
              s = sPrior, rate = LnTrunc(-7.1, 0.5))
  for (p in c(fbdNames, if (scaled) "s")) pri[[p]] <- FbdConditional(p, JP)
  pri
}
# Hymenoptera: shape of the Gamma rates among characters, alpha ~ Exp(1), as
# the comparator's alpha_morph ~ dnExponential(1); TempoMCMC(charPriors = )
CetCharPriors <- function() {
  list(alpha = function(a) if (is.finite(a) && a > 0) -a else -Inf)
}
# Log prior of all parameters (for outputs outside TempoMCMC()): the
# parameters of `priors` (joint FBD-rate prior included) and alpha, if present
CetLogPrior <- function(par, priors = CetPriors(), charPriors = CetCharPriors()) {
  LogPriorJoint(par, priors) + if (is.null(par[["alpha"]])) 0 else charPriors$alpha(par[["alpha"]])
}
# Root age ~ offset exponential, minimum 315 Ma, mean 396 Ma (the comparator's
# root_time ~ dnExponential(1/81, offset = 315))
CetAgePriors <- function() {
  list(root = function(a) if (a > 315) -log(81) - (a - 315) / 81 else -Inf)
}

# ---- Starting trees and parameters (#89 review, R89-5; as the penguin pipeline) --
# Production chain k = 4 (rep - 1) + cell (cells 1-4) starts from its own
# draw from the comparator's posterior, data/start-trees/start-tree-<kk>.nwk
# (draw-start-trees.R: 16 post-burn-in samples drawn at random without
# replacement, seed 20261111), when data/start-trees/draws.csv lists k.  Other
# chains (ladder rungs such as b0.5g1, power chains of a cell label,
# replicates beyond the draws) start from the reference tree
# data/start-tree.nwk (prepare.R: the highest posterior sample).
# CET_START_TREE overrides both.  Returns the tree file and the comparator
# sample it came from (a one-row data frame; NULL for an override).
ChainIndex <- function(cell, rep) if (cell > 0L) 4L * (rep - 1L) + cell else NA_integer_
StartTreeFor <- function(chainK, override = Sys.getenv("CET_START_TREE", "")) {
  if (nzchar(override)) return(list(file = override, sample = NULL))
  drawsFile <- file.path(startTreeDir, "draws.csv")
  if (!is.na(chainK) && file.exists(drawsFile)) {
    draws <- utils::read.csv(drawsFile)
    if (chainK %in% draws$k) {
      own <- file.path(startTreeDir, sprintf("start-tree-%02d.nwk", chainK))
      if (!file.exists(own)) stop("draws.csv lists chain ", chainK, " but ", own, " is missing")
      return(list(file = own, sample = draws[draws$k == chainK, ]))
    }
  }
  refFile <- file.path(cetDir, "data", "start-tree-sample.csv")
  list(file = referenceTree, sample = if (file.exists(refFile)) utils::read.csv(refFile))
}
# Log posterior of a start (as TempoMCMC() builds its first state: tempo-path
# draw from .Precompute(), .AddSlices(), .Messages(), .SampleTree()); -Inf with
# attribute "error" if it cannot be computed
CetStartLogPost <- function(par, tree, model, characters, priors = CetPriors(),
                            nSlice = 20) {
  tryCatch({
    lp <- CetLogPrior(par, priors)
    if (!is.finite(lp)) return(-Inf)
    tt <- get(".NoRootAncestor", .ns)(get(".ForModel", .ns)(TempoTree(tree), model))
    la <- CetAgePriors()$root(get(".StartAge", .ns)(tt))
    set <- get(".AddSlices", .ns)(tt, get(".Precompute", .ns)(tt, par, model), nSlice)
    msg <- get(".Messages", .ns)(tt, set)
    draw <- get(".SampleTree", .ns)(tt, set, msg)
    aug <- get(".LogAugmented", .ns)(tt, set, draw[["paths"]], draw[["bin"]])
    len <- get(".CharacterLengths", .ns)(tt, set[["grid"]][["x"]], draw[["paths"]],
                                         par[["rate"]], model[["gamma"]])
    chr <- CetMk(tt, characters, len, alpha = par[["alpha"]])
    lp + la + aug + chr
  }, error = function(e) structure(-Inf, error = conditionMessage(e)))
}
# Starting parameters.  Hymenoptera: lambda, mu, psi, rate (clock_rate) and
# alpha (alpha_morph) are those of the comparator sample that gave the start
# tree, so each production chain starts at its own joint draw from the
# comparator's posterior; the reference tree takes those of its sample
# (data/start-tree-sample.csv).  theta and s, which the comparator lacks, are
# drawn from their priors, theta ~ LogNormal(log(1/120), 1) within +-4 sdlog
# and s ~ Exp(mean 0.587405) below 3, from a stream of their own,
# set.seed(9000 + index), which leaves the chain's random-number state
# untouched; they are redrawn while the start log posterior is not finite
# (with `tree` given; CetStartLogPost()).  Returns the parameters, with
# attributes "attempts" and "logPost".
CetStart <- function(index, sample, tree = NULL, model = NULL, characters = NULL,
                     priors = CetPriors(), nSlice = 20, maxAttempts = 100L) {
  stopifnot(is.data.frame(sample), nrow(sample) == 1L)
  fixed <- list(lambda = sample$lambda, mu = sample$mu, psi = sample$psi,
                rate = sample$clock_rate, alpha = sample$alpha_morph)
  stopifnot(all(vapply(fixed, function(v) is.finite(v) && v > 0, TRUE)))
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(9000 + index)
  lastError <- NULL
  for (attempt in seq_len(maxAttempts)) {
    repeat {
      theta <- exp(stats::rnorm(1, log(1 / 120), 1))
      if (abs(log(theta) - log(1 / 120)) <= 4) break
    }
    repeat {
      s <- stats::rexp(1, 1 / 0.587405)
      if (s <= .cet$sMax) break
    }
    start <- c(fixed[c("lambda", "mu", "psi")], list(theta = theta, s = s),
               fixed[c("rate", "alpha")])
    lp <- if (is.null(tree)) NA_real_ else CetStartLogPost(start, tree, model, characters,
                                                           priors, nSlice)
    if (is.null(tree) || is.finite(lp)) {
      return(structure(start, attempts = attempt, logPost = as.numeric(lp)))
    }
    lastError <- attr(lp, "error")
  }
  stop("CetStart(): no start with a finite log posterior in ", maxAttempts, " draws",
       if (!is.null(lastError)) paste0("; last error: ", lastError))
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
