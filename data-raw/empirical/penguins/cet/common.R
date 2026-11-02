# Penguin CET pipeline: definitions shared by run-chain.R, bayes-factors.R and
# summarise.R.  Source from the package root after loading EvoTempo.
#
# The engine (EvoTempo on origin-root-age-18, PR #43) lacks four things this
# empirical analysis needs; as in rung 7 (data-raw/rung7-power.R), this file
# rebinds two functions inside the EvoTempo namespace for the duration of a
# run, without changing the package:
#
# 1. MkLikelihood() takes one number of states for all characters, and no
#    ascertainment correction.  The penguin matrix has partitions with K = 2-7
#    states and only variable characters, so the rebound MkLikelihood() scores
#    each K-partition separately and applies the Mkv correction (Lewis 2001):
#    log L - n_K log(1 - K P(all tips in state 0)), P from a constant
#    dummy column, the BEAST2 / RevBayes coding = "variable" convention.  It
#    also raises the character likelihood to a power (`.cet$power`) for
#    character-power stepping-stone runs.  It applies only to character
#    matrices that carry the attribute "parts" (from LoadPenguins()).
# 2. TempoMCMC() returns no per-sample trees or states.  The rebound
#    .TreeStats(), which TempoMCMC() calls at every stored sample, reads the
#    chain state `st` from its caller's frame and records the tree (Newick),
#    the state (tree, parameters, paths, bins) and the character
#    log-likelihood.
# 3. TempoMCMC() cannot checkpoint or resume.  The same hook captures, at
#    every stored sample, everything that the sampling phase carries from one
#    iteration to the next (st, step sizes, move counts, the RNG state); on
#    the first call of a run (TempoMCMC() calls .TreeStats() once before its
#    loop) it can inject a captured state, so a chain runs as a series of
#    segments that continue each other exactly (run-chain.R; tested in
#    smoke.sh).
# 4. Priors are independent per parameter, so the FBD prior of the source
#    analysis (on net diversification, turnover and sampling proportion) cannot
#    be expressed; independent priors on lambda, mu, psi are used instead, and
#    the RevBayes morphology-only comparator uses the same ones.
# See the engine issues linked from README.md.
suppressMessages(library(ape))
.ns <- asNamespace("EvoTempo")
.cet <- new.env()
.cet$power <- 1
.cet$mkv <- TRUE
.cet$records <- list()
.cet$resume <- NULL
.cet$captured <- NULL
.cet$injected <- FALSE
if (!exists(".mkOriginal")) {
  .mkOriginal <- get("MkLikelihood", .ns)
  .statsOriginal <- get(".TreeStats", .ns)
}

PartitionedMk <- function(tree, characters, lengths, nState = 2L,
                          engine = getOption("EvoTempo.engine", "cpp"),
                          power = .cet$power) {
  parts <- attr(characters, "parts")
  if (is.null(parts)) {
    return(.mkOriginal(tree, characters, lengths, nState, engine))
  }
  if (power == 0) {
    return(0)
  }
  ll <- 0
  for (p in parts) {
    ll <- ll + .mkOriginal(tree, p$x, lengths, p$k, engine)
    if (.cet$mkv) {
      lc <- .mkOriginal(tree, p$const, lengths, p$k, engine)
      ll <- ll - ncol(p$x) * log1p(-p$k * exp(lc))
    }
  }
  power * ll
}

RecordingStats <- function(tree, lengths = NULL, inverse = NULL) {
  pf <- parent.frame()
  if (is.null(lengths)) {
    r <- .cet$resume
    if (!is.null(r)) {
      # First call of a resumed segment: continue the captured chain
      for (v in c("st", "step", "jointStep", "ageStep", "moves")) {
        assign(v, r[[v]], envir = pf)
      }
      assign(".Random.seed", r[["seed"]], envir = globalenv())
      .cet$resume <- NULL
      .cet$injected <- TRUE
    }
  } else {
    st <- get("st", envir = pf, inherits = FALSE)
    model <- get("model", envir = pf, inherits = FALSE)
    charLL <- if (.cet$power == 1) {
      st[["charLL"]]
    } else {
      chars <- get("characters", envir = pf, inherits = FALSE)
      if (is.null(chars)) 0 else PartitionedMk(
        st[["tree"]], chars,
        get(".CharacterLengths", .ns)(st[["tree"]], st[["set"]][["grid"]][["x"]],
                                      st[["paths"]], st[["par"]][["rate"]],
                                      model[["gamma"]]), power = 1)
    }
    .cet$records[[length(.cet$records) + 1L]] <- list(
      newick = TempoNewick(st[["tree"]]), charLL = charLL,
      state = st[c("tree", "par", "paths", "bin")])
    .cet$captured <- list(st = st, step = get("step", pf), jointStep = get("jointStep", pf),
                          ageStep = get("ageStep", pf), moves = get("moves", pf),
                          seed = get(".Random.seed", globalenv()))
  }
  .statsOriginal(tree, lengths, inverse)
}

utils::assignInNamespace("MkLikelihood", PartitionedMk, ns = "EvoTempo")
utils::assignInNamespace(".TreeStats", RecordingStats, ns = "EvoTempo")

# Newick of a TempoTree, with branch lengths from ages; sampled ancestors are
# zero-length tips; the stem is the root edge
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
  paste0(Sub(tt[["root"]]), ":", formatC(tt[["stem"]], digits = digits, format = "f"), ";")
}

# ---- Data ----------------------------------------------------------------------
cetDir <- "data-raw/empirical/penguins/cet"

LoadPenguins <- function(startTree = file.path(cetDir, "data", "start-tree.nwk")) {
  tips <- read.csv("data-raw/empirical/penguins/penguins_tips.csv")
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
  # Rounding in the Newick: keep fossil tips inside their ranges
  fos <- which(tipType == "fossil")
  age[fos] <- pmin(pmax(age[fos], ranges[fos, 1]), ranges[fos, 2])
  tt <- TempoTree(tr, ages = age, tipType = tipType, ranges = ranges)
  chars <- read.csv(file.path(cetDir, "data", "characters.csv"), check.names = FALSE)
  parts <- read.csv(file.path(cetDir, "data", "partitions.csv"))
  x <- as.matrix(chars[match(tr$tip.label, chars$taxon), -1])
  storage.mode(x) <- "integer"
  rownames(x) <- tr$tip.label
  attr(x, "parts") <- lapply(sort(unique(parts$nState)), function(k) {
    cols <- which(parts$nState == k)
    list(k = k, x = x[, cols, drop = FALSE], const = matrix(0L, nTip, 1))
  })
  list(tree = tt, characters = x, tips = info)
}

# ---- Model, priors, moves -------------------------------------------------------
cells <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
CellLabel <- function(beta, gamma) sprintf("b%sg%s", format(beta), format(gamma))

CetModel <- function(beta, gamma,
                     nBin = as.integer(Sys.getenv("CET_NBIN", "26")),
                     zMax = as.numeric(Sys.getenv("CET_ZMAX", "5"))) {
  DCETModel(nBin = nBin, zRange = c(-zMax, zMax), beta = beta, gamma = gamma,
            rho = 1, condition = "extant", start = "origin")
}

# Log-normal priors truncated at +-4 sdlog (as rung 7; keeps the stiff E solver
# away from extreme rates), renormalised
LnTrunc <- function(meanlog, sdlog, zTrunc = 4) {
  logMass <- log(2 * stats::pnorm(zTrunc) - 1)
  force(meanlog)
  force(sdlog)
  function(x) {
    if (!(x > 0) || abs(log(x) - meanlog) > zTrunc * sdlog) return(-Inf)
    stats::dlnorm(x, meanlog, sdlog, log = TRUE) - logMass
  }
}
# Rates per Myr.  lambda, mu, psi: broad, the same in the RevBayes comparator
# (revbayes/penguins_morph_fbd.Rev).  theta: reversion time 1/theta with median
# 20 Myr, 95% 2.8-140 Myr.  s: the UCLN log-sd prior of the source analysis
# and the RevBayes baselines (Gamma, shape 0.5396, scale 0.3819).  rate: the
# morphological clock prior of the source analysis, LogNormal(-5.5, 2).
CetPriors <- function() {
  list(lambda = LnTrunc(log(0.1), 1.5), mu = LnTrunc(log(0.1), 1.5),
       psi = LnTrunc(log(0.05), 1.5), theta = LnTrunc(log(0.05), 1),
       s = GammaTrunc(0.5396, 0.3819, 3), rate = LnTrunc(-5.5, 2))
}
# Gamma prior truncated to (0, hi], renormalised.  The bound matters for s: in
# the (0, 0) cell, where s enters only the path density, the burn-in step size
# for s grew to 12.8 and a proposal of s = 86.8 reached the generator, whose
# self-rate theta e^x (alpha = 1) overflows at x = 5 s = 434, and the eigen()
# failure in TempoEquilibrium() escapes TempoMCMC() as an error (12-minute
# pilot, 2026-10-31).  Above 3 lies 8.7e-5 of the untruncated prior mass.
GammaTrunc <- function(shape, scale, hi) {
  logMass <- stats::pgamma(hi, shape, scale = scale, log.p = TRUE)
  function(x) {
    if (!(x > 0 && x <= hi)) return(-Inf)
    stats::dgamma(x, shape, scale = scale, log = TRUE) - logMass
  }
}
# Origin ~ U(0, 160) Ma, as the source analysis and the RevBayes baselines
CetAgePriors <- function() {
  list(origin = function(a) if (a > 0 && a < 160) -log(160) else -Inf)
}
# Near the RevBayes total-evidence round-1 posterior means (lambda 0.20,
# mu 0.16; psi from the sampling proportion; clock_morph 0.017), jittered per
# replicate by exp(N(0, 0.2^2))
CetStart <- function(rep) {
  start <- list(lambda = 0.2, mu = 0.15, psi = 0.05, theta = 0.05, s = 0.2, rate = 0.017)
  if (rep > 1) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
    set.seed(9000 + rep)
    start <- lapply(start, function(v) v * exp(stats::rnorm(1, 0, 0.2)))
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  }
  start
}
# Default schedule (TempoMoves(): local path moves for 0.25 of the tree and a
# shift with gamma != 0; whole-tree draws and joint sweeps with gamma = 0)
# plus tree moves in proportion to the tree, scaled by CET_TREE_MOVES
CetMoves <- function(tt, characters, model,
                     scale = as.numeric(Sys.getenv("CET_TREE_MOVES", "1"))) {
  nBif <- sum(!tt$sampled)
  nFossil <- sum(tt$tipType != "extant")
  c(TempoMoves(tt, characters, model),
    nodeAge = ceiling(scale * 0.25 * nBif), fossilAge = ceiling(scale * 0.25 * nFossil),
    slide = ceiling(scale * 0.15 * nBif), saJump = ceiling(scale * 0.1 * nFossil),
    originAge = 1)
}

# ---- Utilities -------------------------------------------------------------------
LogSumExp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}
LogMeanExp <- function(x) LogSumExp(x) - log(length(x))
# Geyer initial monotone sequence ESS (as rungs 6-7)
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
# Gelman-Rubin potential scale reduction over a list of chains
Psrf <- function(chains) {
  chains <- lapply(chains, function(x) x[is.finite(x)])
  n <- min(lengths(chains))
  if (length(chains) < 2 || n < 4) return(NA_real_)
  m <- vapply(chains, function(x) mean(utils::tail(x, n)), 0)
  w <- mean(vapply(chains, function(x) stats::var(utils::tail(x, n)), 0))
  b <- n * stats::var(m)
  if (w == 0) return(NA_real_)
  sqrt(((n - 1) / n * w + b / n) / w)
}
