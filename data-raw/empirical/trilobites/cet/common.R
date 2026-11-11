# Trilobite CET pipeline: definitions shared by run-chain.R, bayes-factors.R,
# check-*.R and summarise.R.  Source from the package root after loading
# EvoTempo.  Adapted from data-raw/empirical/penguins/cet/common.R; see
# README.md for every difference from the penguin pipeline.
#
# Time is in Myr before the end of the Cambrian sampling window, the
# Cambrian-Ordovician boundary at 485.4 Ma (as the RevBayes comparator,
# ../revbayes): age = Ma - 485.4, so all 107 tips are fossils with ages
# 0.6-33.6 and nothing is sampled at age 0 (rho = 0).  EvoTempo has no offset
# argument: age 0 is the end of the process, and the shift is the whole
# "offset".  The tree is conditioned on the stem lineage leaving at least one
# sample (condition = "survival", start = "origin"), RevBayes's
# condition = "sampling" in origin mode (README, "Conditioning").
#
# As in the penguin pipeline, this file rebinds functions inside the EvoTempo
# namespace for the duration of a run, without changing the package:
#
# 1. MkLikelihood() takes one number of states, symmetric Mk, one state or NA
#    per cell, and no ascertainment correction.  The rebound PartitionedMk()
#    scores the RevBayes comparator's character model minus its Gamma rate
#    variation: unordered characters by number of states K (C++ engine for
#    characters without ambiguous cells), characters with polymorphic cells
#    with those cells ambiguous (R pruning, SymLogLik()), characters 45 and 47
#    under the ordered (adjacent-state) model rescaled to mean rate 1
#    (OrdLogLik()), and the Mkv correction once per missing-data pattern and
#    partition, as RevBayes dnPhyloCTMC(coding = "variable") (#65).  Power
#    `.cet$power` multiplies the corrected log-likelihood.  It applies only to
#    character matrices carrying the attribute "cet" (LoadTrilobites()).
# 2. .TreeStats(): per-sample trees and states, and exact resume (as the
#    penguin pipeline; RecordingStats()).
# 3. .ProposeTreeMove(): the RevBayes comparator (and the source analyses)
#    constrain the ingroup (all but the three fallotaspidoids) to be
#    monophyletic; EvoTempo's tree moves know no constraints.  The rebound
#    proposal returns NULL (a rejection) for a proposed tree that breaks the
#    constraint, which makes the chain target the posterior restricted to
#    constrained trees (the proposal kernel is unchanged, so the Hastings
#    ratio of every allowed move is unchanged).
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
.cet$ingroup <- NULL
.cet$constraintRejections <- list()
if (!exists(".mkOriginal")) {
  .mkOriginal <- get("MkLikelihood", .ns)
  .statsOriginal <- get(".TreeStats", .ns)
  .proposeOriginal <- get(".ProposeTreeMove", .ns)
}
triDir <- "data-raw/empirical/trilobites"
cetDir <- file.path(triDir, "cet")
windowEnd <- 485.4
# CetPrune(): C++ pruning of the columns MkLikelihood() cannot score (prune.cpp)
if (!exists("CetPrune")) Rcpp::sourceCpp(file.path(cetDir, "prune.cpp"))

# ---- Character likelihood -----------------------------------------------------
PartitionedMk <- function(tree, characters, lengths, nState = 2L,
                          engine = getOption("EvoTempo.engine", "cpp"),
                          power = .cet$power, mkv = .cet$mkv) {
  cm <- attr(characters, "cet")
  if (is.null(cm)) {
    return(.mkOriginal(tree, characters, lengths, nState, engine))
  }
  if (power == 0) {
    return(0)
  }
  tree <- TempoTree(tree)
  ll <- 0
  for (p in cm$cpp) {
    ll <- ll + .mkOriginal(tree, p$x, lengths, p$k, engine)
  }
  if (TempoEngine(engine) == "cpp") {
    # C++ route: every column of cm$pack in one pass
    pk <- cm$pack
    kids <- tree[["kids"]]
    storage.mode(kids) <- "integer"
    col <- CetPrune(as.integer(tree[["postorder"]]), kids, tree[["nTip"]], tree[["root"]],
                    as.double(lengths), pk$tip, pk$k, pk$ordId, pk$ordV, pk$ordLam)
    ll <- ll + sum(col[pk$role == "data"])
    if (mkv == "pattern") {
      ll <- ll - sum(cm$mkv$count * log1p(-cm$mkv$k * exp(col[pk$role == "mkv"]))) -
        sum(vapply(split(col[pk$role == "ordConst"], pk$ordChar[pk$role == "ordConst"]),
                   function(x) log1p(-sum(exp(x))), 0))
    } else if (mkv != "none") {
      stop("`.cet$mkv` must be \"pattern\" or \"none\"")
    }
    return(power * ll)
  }
  # R route (the reference for check-mkv.R)
  if (!is.null(cm$sym)) {
    ll <- ll + sum(SymLogLik(tree, cm$sym, lengths))
  }
  ord <- lapply(cm$ord, function(o) OrdLogLik(tree, o, lengths))
  ll <- ll + sum(vapply(ord, `[[`, 0, 1))
  ll <- ll + switch(
    mkv,
    pattern = {
      pat <- cm$mkv
      stopifnot(!is.null(pat))
      (if (length(pat$k)) {
        -sum(pat$count * log1p(-pat$k * exp(ConstantPatternLogLik(tree, pat, lengths))))
      } else 0) -
        sum(vapply(ord, function(x) log1p(-sum(exp(x[-1]))), 0))
    },
    none = 0,
    stop("`.cet$mkv` must be \"pattern\" or \"none\""))
  power * ll
}

# Pruning of characters with ambiguous cells under symmetric Mk (uniform root),
# all K at once: partial likelihoods are kMax x nChar matrices, with the rows
# above each character's K held at zero (sym$valid).  Returns the
# log-likelihood of each column.  Same recursion and per-node rescaling as the
# R engine of MkLikelihood().
SymLogLik <- function(tree, sym, lengths) {
  nTip <- tree[["nTip"]]
  kids <- tree[["kids"]]
  k <- sym$k
  kMax <- nrow(sym$valid)
  nChar <- length(k)
  partial <- vector("list", nTip + tree[["nNode"]])
  partial[seq_len(nTip)] <- sym$tip
  logScale <- numeric(nChar)
  kRow <- rep(k, each = kMax)
  for (v in tree[["postorder"]]) {
    value <- sym$valid
    for (kid in kids[v - nTip, ]) {
      ell <- lengths[[kid]]
      stay <- 1 / kRow + (kRow - 1) / kRow * exp(-kRow * ell / (kRow - 1))
      move <- (1 - stay) / (kRow - 1)
      L <- partial[[kid]]
      value <- value * (move * rep(colSums(L), each = kMax) + (stay - move) * L)
    }
    scale <- value[cbind(max.col(t(value), "first"), seq_len(nChar))]
    partial[[v]] <- value / rep(scale, each = kMax)
    logScale <- logScale + log(scale)
  }
  log(colSums(partial[[tree[["root"]]]]) / k) + logScale
}

# Ordered character (adjacent-state changes at equal rates, rescaled to mean
# rate 1, uniform root): columns of o$tip are the data column and the K
# constant columns of its Mkv correction (observed tips in state j, others
# unknown).  Returns the log-likelihood of each column.
OrdLogLik <- function(tree, o, lengths) {
  nTip <- tree[["nTip"]]
  kids <- tree[["kids"]]
  nCol <- ncol(o$tip[[1]])
  partial <- vector("list", nTip + tree[["nNode"]])
  partial[seq_len(nTip)] <- o$tip
  logScale <- numeric(nCol)
  for (v in tree[["postorder"]]) {
    value <- 1
    for (kid in kids[v - nTip, ]) {
      P <- o$vectors %*% (exp(o$values * lengths[[kid]]) * o$inverse)
      value <- value * (P %*% partial[[kid]])
    }
    scale <- apply(value, 2, max)
    partial[[v]] <- value / rep(scale, each = o$k)
    logScale <- logScale + log(scale)
  }
  log(colSums(partial[[tree[["root"]]]]) / o$k) + logScale
}

OrderedQ <- function(k) {
  r <- k / (2 * (k - 1))
  Q <- matrix(0, k, k)
  Q[cbind(seq_len(k - 1), 2:k)] <- r
  Q[cbind(2:k, seq_len(k - 1))] <- r
  diag(Q) <- -rowSums(Q)
  Q
}

# Distinct missing-data patterns within each unordered partition (K): for
# every pattern, K, the tips observed (logical, one row per tip; a partly
# ambiguous cell counts as observed; ?, - and a cell ambiguous over all K
# states as unobserved, as RevBayes with treatUnknownAsGap) and the number of
# characters that share it
MkvPatterns <- function(observed, k) {
  key <- paste(k, apply(observed, 2, function(o) paste(which(o), collapse = ",")))
  first <- !duplicated(key)
  pat <- list(k = k[first], observed = observed[, first, drop = FALSE],
              count = tabulate(match(key, key[first]), sum(first)))
  stopifnot(all(colSums(pat$observed) >= 2), sum(pat$count) == ncol(observed))
  pat
}

# log P(every observed tip in state 0) under symmetric Mk, per pattern (as the
# penguin pipeline, which checks it against MkLikelihood() column by column)
ConstantPatternLogLik <- function(tree, pat, lengths) {
  nTip <- tree[["nTip"]]
  k <- pat$k
  kU <- unique(k)
  ki <- match(k, kU)
  L0 <- matrix(1, length(k), nTip + tree[["nNode"]])
  L1 <- L0
  L1[, seq_len(nTip)] <- t(1 - pat$observed)
  kids <- tree[["kids"]]
  logScale <- 0
  for (v in tree[["postorder"]]) {
    a0 <- 1
    a1 <- 1
    for (kid in kids[v - nTip, ]) {
      stay <- 1 / kU + (kU - 1) / kU * exp(-kU * lengths[[kid]] / (kU - 1))
      off <- 1 - stay
      same <- stay + (kU - 2) / (kU - 1) * off
      l0 <- L0[, kid]
      l1 <- L1[, kid]
      a0 <- a0 * (stay[ki] * l0 + off[ki] * l1)
      a1 <- a1 * ((off / (kU - 1))[ki] * l0 + same[ki] * l1)
    }
    scale <- pmax(a0, a1)
    L0[, v] <- a0 / scale
    L1[, v] <- a1 / scale
    logScale <- logScale + log(scale)
  }
  root <- tree[["root"]]
  log((L0[, root] + (k - 1) * L1[, root]) / k) + logScale
}

# Character model from the coded cells (tips in tree order) and partitions.csv
BuildCharacters <- function(cells, parts) {
  nTip <- nrow(cells)
  stopifnot(ncol(cells) == nrow(parts))
  States <- function(cell) if (cell == "") integer(0) else as.integer(strsplit(cell, "/")[[1]])
  ambiguous <- parts$nAmbiguous > 0
  single <- cells != "" & !grepl("/", cells)
  first <- matrix(NA_integer_, nTip, ncol(cells), dimnames = dimnames(cells))
  first[single] <- as.integer(cells[single])
  # A cell is observed for the Mkv masks unless it is ? or -, or ambiguous
  # over all K states, which RevBayes also treats as missing (check-mkv.R:
  # character 12 has one {01} cell; counting it as observed moves partition
  # u2 by 1.4e-3 nats at rate 0.005)
  nStates <- matrix(lengths(strsplit(cells, "/")), nTip)
  observed <- cells != "" & nStates < rep(parts$k, each = nTip)
  # Every character must vary among its cells (the Mkv conditional)
  varies <- vapply(seq_len(ncol(cells)), function(j) {
    length(unique(unlist(lapply(cells[, j], States)))) > 1
  }, TRUE)
  if (!all(varies)) stop("characters constant among their scored tips: ", which(!varies))
  TipPartial <- function(cell, k) {
    p <- numeric(k)
    s <- States(cell)
    if (length(s)) p[s + 1L] <- 1 else p[] <- 1
    p
  }
  cppCols <- which(!parts$ordered & !ambiguous)
  cpp <- lapply(sort(unique(parts$k[cppCols])), function(k) {
    cols <- cppCols[parts$k[cppCols] == k]
    list(k = k, x = first[, cols, drop = FALSE], chars = parts$char[cols])
  })
  symCols <- which(!parts$ordered & ambiguous)
  sym <- if (length(symCols)) {
    k <- parts$k[symCols]
    kMax <- max(k)
    valid <- outer(seq_len(kMax), k, "<=") * 1
    list(k = k, valid = valid, chars = parts$char[symCols],
         tip = lapply(seq_len(nTip), function(i) {
           vapply(seq_along(symCols), function(j) {
             c(TipPartial(cells[i, symCols[j]], k[j]), numeric(kMax - k[j]))
           }, numeric(kMax))
         }))
  }
  ord <- lapply(which(parts$ordered), function(j) {
    k <- parts$k[j]
    e <- eigen(OrderedQ(k), symmetric = TRUE)
    list(k = k, char = parts$char[j], vectors = e$vectors, values = e$values,
         inverse = t(e$vectors),
         tip = lapply(seq_len(nTip), function(i) {
           cbind(TipPartial(cells[i, j], k),
                 if (observed[i, j]) diag(k) else matrix(1, k, k))
         }))
  })
  unordered <- which(!parts$ordered)
  mkv <- if (length(unordered)) {
    MkvPatterns(observed[, unordered, drop = FALSE], parts$k[unordered])
  } else {
    list(k = integer(0), observed = matrix(FALSE, nTip, 0), count = integer(0))
  }
  # Columns for CetPrune(): the ambiguous characters, the unordered Mkv
  # patterns (observed tips in state 0, others unknown), and per ordered
  # character its data column and its K constant columns
  colK <- integer(0)
  colOrd <- integer(0)
  role <- character(0)
  ordChar <- integer(0)
  tipCols <- vector("list", nTip)
  Add <- function(partials, k, ordId, r, oc = NA_integer_) {
    colK <<- c(colK, k)
    colOrd <<- c(colOrd, ordId)
    role <<- c(role, r)
    ordChar <<- c(ordChar, oc)
    for (i in seq_len(nTip)) tipCols[[i]] <<- c(tipCols[[i]], partials[[i]])
  }
  for (j in seq_along(symCols)) {
    Add(lapply(seq_len(nTip), function(i) TipPartial(cells[i, symCols[j]], parts$k[symCols[j]])),
        parts$k[symCols[j]], -1L, "data")
  }
  for (j in seq_along(mkv$k)) {
    Add(lapply(seq_len(nTip), function(i) {
      if (mkv$observed[i, j]) c(1, numeric(mkv$k[j] - 1L)) else rep(1, mkv$k[j])
    }), mkv$k[j], -1L, "mkv")
  }
  for (o in seq_along(ord)) {
    oo <- ord[[o]]
    for (cc in seq_len(oo$k + 1L)) {
      Add(lapply(seq_len(nTip), function(i) oo$tip[[i]][, cc]), oo$k, o - 1L,
          if (cc == 1L) "data" else "ordConst", oo$char)
    }
  }
  pack <- list(k = as.integer(colK), ordId = as.integer(colOrd), role = role, ordChar = ordChar,
               tip = as.double(unlist(tipCols)),
               ordV = lapply(ord, `[[`, "vectors"), ordLam = lapply(ord, `[[`, "values"))
  list(cpp = cpp, sym = sym, ord = ord, mkv = mkv, pack = pack)
}

# Every character by the R pruners (SymLogLik() for all unordered
# characters, OrdLogLik(), ConstantPatternLogLik()), for check-mkv.R: a route
# to PartitionedMk() that shares no code with MkLikelihood() or CetPrune()
AllByR <- function(tree, cells, parts, lengths, mkv = "pattern") {
  p2 <- parts
  p2$nAmbiguous[!p2$ordered] <- pmax(p2$nAmbiguous[!p2$ordered], 1L)
  cm <- BuildCharacters(cells, p2)
  stopifnot(!length(cm$cpp))
  x <- matrix(0L, nrow(cells), ncol(cells))
  attr(x, "cet") <- cm
  PartitionedMk(tree, x, lengths, engine = "R", power = 1, mkv = mkv)
}

utils::assignInNamespace("MkLikelihood", PartitionedMk, ns = "EvoTempo")

# ---- Ingroup constraint ----------------------------------------------------------
# TRUE if some node has exactly the ingroup tips below it (sampled ancestors
# count as tips)
IngroupMonophyletic <- function(tt, ingroup) {
  nTip <- tt[["nTip"]]
  nIn <- c(as.integer(ingroup), integer(tt[["nNode"]]))
  nAll <- c(rep(1L, nTip), integer(tt[["nNode"]]))
  for (v in tt[["postorder"]]) {
    k <- tt[["kids"]][v - nTip, ]
    nIn[v] <- nIn[k[1]] + nIn[k[2]]
    nAll[v] <- nAll[k[1]] + nAll[k[2]]
  }
  target <- sum(ingroup)
  any(nIn == target & nAll == target)
}
ConstrainedProposal <- function(move, tree, set, st) {
  prop <- .proposeOriginal(move, tree, set, st)
  if (!is.null(prop) && !is.null(.cet$ingroup) &&
      !IngroupMonophyletic(prop[["tree"]], .cet$ingroup)) {
    .cet$constraintRejections[[move]] <- (.cet$constraintRejections[[move]] %||% 0) + 1
    return(NULL)
  }
  prop
}
utils::assignInNamespace(".ProposeTreeMove", ConstrainedProposal, ns = "EvoTempo")

# ---- Per-sample records and exact resume (as the penguin pipeline) ---------------
RecordingStats <- function(tree, lengths = NULL, inverse = NULL) {
  pf <- parent.frame()
  if (is.null(lengths)) {
    r <- .cet$resume
    if (!is.null(r)) {
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

# ---- Data -------------------------------------------------------------------------
StartTreeFile <- function(rep) {
  file.path(cetDir, "data", sprintf("start-tree-%d.nwk", (rep - 1L) %% 6L + 1L))
}

LoadTrilobites <- function(startTree = StartTreeFile(1L)) {
  tips <- read.csv(file.path(triDir, "trilobites_tips.csv"))
  tr <- read.tree(startTree)
  nTip <- Ntip(tr)
  info <- tips[match(tr$tip.label, tips$taxon), ]
  stopifnot(!anyNA(info$taxon), all(info$min_age == info$max_age))
  tipAge <- info$min_age - windowEnd
  depth <- node.depth.edgelength(tr)
  desc <- prop.part(tr)
  age <- c(tipAge, vapply(seq_len(tr$Nnode), function(i) {
    d <- desc[[i]]
    mean(tipAge[d] + depth[d] - depth[nTip + i])
  }, 0))
  len <- numeric(nTip + tr$Nnode)
  len[tr$edge[, 2]] <- tr$edge.length
  parent <- integer(nTip + tr$Nnode)
  parent[tr$edge[, 2]] <- tr$edge[, 1]
  tipType <- ifelse(len[seq_len(nTip)] < 1e-9, "ancestor", "fossil")
  sa <- which(tipType == "ancestor")
  age[parent[sa]] <- tipAge[sa]
  # Point ages: each range is the age itself, so no move changes a tip age
  tt <- TempoTree(tr, ages = age, tipType = tipType, ranges = cbind(tipAge, tipAge))
  chars <- read.csv(file.path(cetDir, "data", "characters.csv"), check.names = FALSE,
                    colClasses = "character")
  parts <- read.csv(file.path(cetDir, "data", "partitions.csv"))
  cells <- as.matrix(chars[match(tr$tip.label, chars$taxon), -1])
  cells[is.na(cells)] <- ""
  rownames(cells) <- tr$tip.label
  stopifnot(identical(colnames(cells), paste0("c", parts$char)))
  cm <- BuildCharacters(cells, parts)
  x <- matrix(NA_integer_, nTip, ncol(cells), dimnames = dimnames(cells))
  single <- cells != "" & !grepl("/", cells)
  x[single] <- as.integer(cells[single])
  attr(x, "cet") <- cm
  cons <- read.csv(file.path(triDir, "trilobites_constraints.csv"))
  ingroupTaxa <- strsplit(cons$taxa[cons$constraint == "Ingroup_node"], " ")[[1]]
  ingroup <- tr$tip.label %in% ingroupTaxa
  stopifnot(sum(ingroup) == 104, IngroupMonophyletic(tt, ingroup))
  .cet$ingroup <- ingroup
  list(tree = tt, characters = x, cells = cells, parts = parts, tips = info, ingroup = ingroup)
}

# ---- Model, priors, moves -------------------------------------------------------------
cellTable <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
CellLabel <- function(beta, gamma) sprintf("b%sg%s", format(beta), format(gamma))
GridSettings <- function() {
  list(nBin = as.integer(Sys.getenv("CET_NBIN", "26")),
       zMax = as.numeric(Sys.getenv("CET_ZMAX", "5")))
}
CetModel <- function(beta, gamma, nBin = GridSettings()$nBin, zMax = GridSettings()$zMax) {
  DCETModel(nBin = nBin, zRange = c(-zMax, zMax), beta = beta, gamma = gamma,
            rho = 0, condition = "survival", start = "origin")
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

# Marginal priors on lambda, mu and psi induced by the RevBayes comparator's
# prior on (d, nu, s): net diversification d ~ LogNormal(-3.5, 1.5), turnover
# nu ~ U(0, 1), sampling proportion s ~ U(0, 1), with lambda = d / (1 - nu),
# mu = nu lambda, psi = s mu / (1 - s).  On the log scale
#   log lambda = log d + E,        E = -log(1 - nu) ~ Exp(1)
#   log mu     = log d + L1,       L1 = logit(nu) ~ Logistic(0, 1)
#   log psi    = log d + L1 + L2,  L2 = logit(s) ~ Logistic(0, 1),
# computed by discrete convolution of cell probabilities (step h), truncated
# to the central 1 - 2 pnorm(-4) of their mass (as the penguin pipeline's
# +-4 sdlog), renormalised, and interpolated by a spline of the log density.
# The joint prior is NOT a product of these marginals (it has mu < lambda and
# a dependence between all three); TempoMCMC() takes independent priors only
# (README, "Priors").  check-priors.R checks the marginals by simulation.
FbdMarginalPriors <- function(meanlog = -3.5, sdlog = 1.5, h = 0.002,
                              tail = 2 * stats::pnorm(-4)) {
  Cells <- function(lo, hi, cdf) {
    x <- seq(lo, hi, by = h)
    list(x0 = lo, p = diff(cdf(c(x - h / 2, hi + h / 2))))
  }
  N <- Cells(meanlog - 12 * sdlog, meanlog + 12 * sdlog,
             function(q) stats::pnorm(q, meanlog, sdlog))
  E <- Cells(0, 40, function(q) stats::pexp(pmax(q, 0)))
  L <- Cells(-40, 40, stats::plogis)
  Conv <- function(a, b) {
    p <- stats::convolve(a$p, rev(b$p), type = "open")
    list(x0 = a$x0 + b$x0, p = pmax(p, 0))
  }
  Prior <- function(d) {
    x <- d$x0 + h * (seq_along(d$p) - 1)
    cum <- cumsum(d$p) / sum(d$p)
    lo <- x[which(cum >= tail / 2)[1]]
    hi <- x[which(cum >= 1 - tail / 2)[1]]
    sel <- x >= lo - 0.1 & x <= hi + 0.1 & d$p > 0
    thin <- which(sel)[seq(1, sum(sel), by = 5)]
    mass <- sum(d$p[x >= lo & x <= hi]) / sum(d$p)
    spl <- stats::splinefun(x[thin], log(d$p[thin] / (h * sum(d$p))), method = "natural")
    force(lo)
    force(hi)
    f <- function(v) {
      if (!(v > 0)) return(-Inf)
      y <- log(v)
      if (y < lo || y > hi) return(-Inf)
      spl(y) - log(mass) - y
    }
    attr(f, "logRange") <- c(lo, hi)
    f
  }
  NL <- Conv(N, L)
  list(lambda = Prior(Conv(N, E)), mu = Prior(NL), psi = Prior(Conv(NL, L)))
}

# The comparator's joint prior on (lambda, mu, psi), untruncated: (d, nu, s)
# as above, with Jacobian |d(d, nu, s) / d(lambda, mu, psi)| =
# (lambda - mu) mu / (lambda^2 (mu + psi)^2); zero unless lambda > mu.
JointFbdLogPrior <- function(lambda, mu, psi, meanlog = -3.5, sdlog = 1.5) {
  if (!(lambda > mu && mu > 0 && psi > 0)) return(-Inf)
  d <- lambda - mu
  stats::dlnorm(d, meanlog, sdlog, log = TRUE) + log(d) + log(mu) - 2 * log(lambda) -
    2 * log(mu + psi)
}
# Log importance weight from the CET target (independent marginal priors) to
# the target with the comparator's joint prior: reweighting the samples by
# exp(JointPriorLogWeight()) gives posterior summaries and Bayes factors under
# the joint prior (bayes-factors.R, summarise.R), as far as the weights'
# effective sample size allows
JointPriorLogWeight <- function(par, priors = CetPriors()) {
  JointFbdLogPrior(par$lambda, par$mu, par$psi) - priors$lambda(par$lambda) -
    priors$mu(par$mu) - priors$psi(par$psi)
}

# Log-SD sigma of the UCLN branch rates when the coefficient of variation
# CV = sqrt(exp(sigma^2) - 1) ~ Exponential(rate 3), as the RevBayes
# comparator (ucld_cv ~ dnExponential(3), ucld_sd = sqrt(ln(1 + cv^2))):
# p(sigma) = 3 exp(-3 CV) dCV/dsigma, dCV/dsigma = sigma exp(sigma^2) / CV.
# Bounded at 3 for numerical safety (as the penguin pipeline: the self-rate
# theta e^x overflows for large s); the mass above 3 (CV > 90) is exp(-270).
UclnSdPrior <- function(rate = 3, hi = 3) {
  function(x) {
    if (!(x > 0 && x <= hi)) return(-Inf)
    cv <- sqrt(expm1(x * x))
    log(rate) - rate * cv + log(x) + x * x - log(cv)
  }
}

# Rates per Myr.  lambda, mu, psi: marginals of the comparator's (d, nu, s)
# prior.  theta: as the penguin pipeline (reversion time median 20 Myr; no
# counterpart in RevBayes).  s: the comparator's UCLN log-SD prior.  rate:
# U(0, 100), the comparator's clock_morph prior (its mean branch rate).
CetPriors <- function() {
  if (is.null(.cet$fbdPriors)) .cet$fbdPriors <- FbdMarginalPriors()
  c(.cet$fbdPriors,
    list(theta = LnTrunc(log(0.05), 1), s = UclnSdPrior(),
         rate = function(x) if (x > 0 && x < 100) -log(100) else -Inf))
}
# Origin ~ U(0.6, 1000.6) Myr before 485.4 Ma (U(0, 1000) above the youngest
# tip), as the comparator's origin_time
CetAgePriors <- function() {
  list(origin = function(a) if (a > 0.6 && a < 1000.6) -log(1000) else -Inf)
}
# Starting values: the RevBayes sample whose tree starts the replicate
# (prepare.R, start-values.csv): lambda, mu, psi, rate = clock_morph,
# s = ucld_sd; theta 0.05 (no counterpart)
CetStart <- function(rep) {
  sv <- read.csv(file.path(cetDir, "data", "start-values.csv"))
  r <- sv[sv$rep == (rep - 1L) %% 6L + 1L, ]
  list(lambda = r$lambda, mu = r$mu, psi = r$psi, theta = 0.05,
       s = min(r$ucld_sd, 2.5), rate = r$clock_morph)
}
# TempoMoves() plus tree moves in proportion to the tree; no fossilAge (point
# ages), CET_TREE_MOVES scales the tree moves
CetMoves <- function(tt, characters, model,
                     scale = as.numeric(Sys.getenv("CET_TREE_MOVES", "1"))) {
  nBif <- sum(!tt$sampled)
  c(TempoMoves(tt, characters, model),
    nodeAge = ceiling(scale * 0.25 * nBif), slide = ceiling(scale * 0.15 * nBif),
    saJump = ceiling(scale * 0.1 * tt$nTip), originAge = 1)
}

# ---- Utilities ---------------------------------------------------------------------
LogSumExp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}
LogMeanExp <- function(x) LogSumExp(x) - log(length(x))
# Convergence through MkPrime's exported ConvergenceDiagnostics(), which needs
# an MkPosterior: rank-normalised split-chain R-hat (Vehtari et al. 2021) over
# chains and ESS of the pooled samples.  `chains` is a list of numeric
# matrices with the same column names (one per chain).
MkpPosterior <- function(chains) {
  chains <- lapply(chains, as.matrix)
  structure(list(samples = do.call(rbind, chains),
                 per_run = lapply(chains, function(x) list(samples = x)),
                 nRuns = length(chains), burnin = 0L), class = "MkPosterior")
}
RhatEss <- function(chains) {
  d <- suppressWarnings(MkPrime::ConvergenceDiagnostics(MkpPosterior(chains)))
  q <- names(d$ess)
  data.frame(quantity = q, ess = unname(d$ess),
             rhat = if (is.null(d$rhat)) NA_real_ else unname(d$rhat[q]))
}
# ESS of one sequence (MkPrime, as above)
Ess <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || stats::var(x) == 0) return(NA_real_)
  RhatEss(list(matrix(x, ncol = 1, dimnames = list(NULL, "x"))))$ess[[1]]
}
