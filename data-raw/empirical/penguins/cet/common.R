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
#    each K-partition separately and applies the Mkv correction (Lewis 2001)
#    per missing-data pattern, as RevBayes dnPhyloCTMC(coding = "variable")
#    does (#65): for each partition and each distinct set of observed tips
#    among its characters, log L - n_pattern log(1 - K P(all observed tips in
#    state 0)), with P from a constant dummy column in which the unobserved
#    tips are unknown (ConstantPatternLogLik()).  Under symmetric Mk with a
#    uniform root, P(all observed tips in state j) is the same for every j.
#    `.cet$mkv = "allTips"` (environment CET_MKV; part of the chain label
#    and checkpoint, see run-chain.R) restores the earlier
#    variant (one constant column with every tip observed, which
#    under-corrects when tips are missing); "none" drops the correction.  The rebound function also raises the
#    character likelihood, correction included, to a power (`.cet$power`)
#    for character-power stepping-stone runs.  It applies only to character
#    matrices that carry the attributes "parts" and "mkv" (from
#    LoadPenguins()).
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
.cet$mkv <- Sys.getenv("CET_MKV", "pattern")
stopifnot("CET_MKV must be \"pattern\", \"allTips\" or \"none\"" =
            .cet$mkv %in% c("pattern", "allTips", "none"))
# Label suffix of a non-default Mkv variant (run-chain.R, bayes-factors.R,
# check-smoke.R): "" for the default, "_mkv<variant>" otherwise (#74)
MkvSuffix <- function(mkv = .cet$mkv) if (mkv == "pattern") "" else paste0("_mkv", mkv)
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
                          power = .cet$power, mkv = .cet$mkv) {
  parts <- attr(characters, "parts")
  if (is.null(parts)) {
    return(.mkOriginal(tree, characters, lengths, nState, engine))
  }
  if (power == 0) {
    return(0)
  }
  tree <- TempoTree(tree)
  ll <- 0
  for (p in parts) {
    ll <- ll + .mkOriginal(tree, p$x, lengths, p$k, engine)
  }
  ll <- ll + switch(
    mkv,
    pattern = {
      pat <- attr(characters, "mkv")
      # Without the attribute, ConstantPatternLogLik() would return numeric(0)
      # and the correction would silently become 0 (#74)
      stopifnot("characters with \"parts\" need the \"mkv\" attribute (MkvPatterns())" =
                  !is.null(pat))
      -sum(pat$count * log1p(-pat$k * exp(ConstantPatternLogLik(tree, pat, lengths))))
    },
    allTips = -sum(vapply(parts, function(p) {
      ncol(p$x) * log1p(-p$k * exp(.mkOriginal(tree, p$const, lengths, p$k, engine)))
    }, 0)),
    none = 0,
    stop("`.cet$mkv` must be \"pattern\", \"allTips\" or \"none\""))
  power * ll
}

# Distinct missing-data patterns within each K-partition: for every pattern,
# K, the tips observed (logical, one row per tip) and the number of
# characters that share it
MkvPatterns <- function(x, nState) {
  observed <- !is.na(x)
  # The Mkv conditional assumes that every character varies among its observed
  # tips; a character constant there (such as character 81, which prepare.R
  # drops) must be removed first (#74)
  varies <- apply(x, 2, function(col) length(unique(col[!is.na(col)])) > 1L)
  if (!all(varies)) {
    stop("characters constant among their observed tips: ",
         paste(utils::head(which(!varies), 10), collapse = ", "))
  }
  key <- paste(nState, apply(observed, 2, function(o) paste(which(o), collapse = ",")))
  first <- !duplicated(key)
  pat <- list(k = nState[first], observed = observed[, first, drop = FALSE],
              count = tabulate(match(key, key[first]), sum(first)))
  # A pattern with fewer than two observed tips is constant whatever the
  # states, so its correction would be log(0)
  stopifnot(all(colSums(pat$observed) >= 2), sum(pat$count) == ncol(x))
  pat
}

# log P(every observed tip in state 0) under Mk with K states and a uniform
# root, for each pattern of `pat` (MkvPatterns()): pruning of one constant
# column per pattern, with the tips not observed in it as unknown.  Under
# symmetric Mk every partial likelihood of such a column has the form
# (L0, L1, ..., L1), so two numbers per node and pattern suffice, and all the
# patterns of all partitions are pruned together.  Same recursion and
# per-node rescaling as the R engine of MkLikelihood(); check-mkv.R checks it
# against MkLikelihood() column by column.
ConstantPatternLogLik <- function(tree, pat, lengths) {
  nTip <- tree[["nTip"]]
  k <- pat$k
  kU <- unique(k)
  ki <- match(k, kU)
  # One column per node (contiguous), one row per pattern
  L0 <- matrix(1, length(k), nTip + tree[["nNode"]])
  L1 <- L0
  L1[, seq_len(nTip)] <- t(1 - pat$observed)
  kids <- tree[["kids"]]
  logScale <- 0
  for (v in tree[["postorder"]]) {
    a0 <- 1
    a1 <- 1
    for (kid in kids[v - nTip, ]) {
      # Per distinct K: P(stay), and (K - 1) P(move) = 1 - P(stay)
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
  attr(x, "mkv") <- MkvPatterns(x, parts$nState)
  list(tree = tt, characters = x, tips = info)
}

# ---- Model, priors, moves -------------------------------------------------------
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
    "^b([0-9.]+)g([0-9.]+)(_mkv(allTips|none))?(_t([0-9.]+))?_run([0-9]+)$", label))[[1]]
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
#   mkv    Mkv variant (optional; default "pattern")
# For example, the production launch (run-array.sh --array=1-16) plus one rung
# of 4 replicates (CET_CELL=b0.5g0, --array=1-4):
#   cell,reps
#   b0g0,4
#   b0g1,4
#   b1g0,4
#   b1g1,4
#   b0.5g0,4
# A rung launched with its own replicate count, a power ladder (power = t) or
# another Mkv variant (mkv = allTips or none) gets its own row.
# Without launch.csv the expected set is the production launch alone: the four
# cells x CET_BF_REPS (default 4) replicates at power 1.
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
    reps <- suppressWarnings(as.integer(L$reps[i]))
    power <- suppressWarnings(as.numeric(L$power[i]))
    if (is.na(reps) || reps < 1L || is.na(power)) {
      stop(sprintf("%s row %d: reps must be a positive integer and power a number", source, i))
    }
    cell <- ParseLabel(paste0(L$cell[i], "_run1"))
    label <- paste0(CellLabel(cell$beta, cell$gamma), MkvSuffix(L$mkv[i]),
                    if (power != 1) paste0("_t", FormatValue(power)), "_run", seq_len(reps))
    do.call(rbind, lapply(label, ParseLabel))
  }))
  out <- out[out$mkv == mkv, ]
  if (anyDuplicated(ChainKey(out))) stop(source, " lists a chain twice")
  attr(out, "source") <- source
  out
}
ChainKey <- function(x) {
  if (!nrow(x)) return(character())
  paste(vapply(x$beta, FormatValue, ""), vapply(x$gamma, FormatValue, ""), x$mkv,
        vapply(x$power, FormatValue, ""), x$rep)
}
# Compare the chains found under `root` (Chains(), one Mkv variant) with the
# expected set.  Stops if any expected chain has no complete segment, if any
# expected chain at power 1 (the chains the bridges use) has no cross.csv, or
# if a chain is found that the launch record does not list.  The only opt-out
# is explicit and per chain: CET_BF_ALLOW_MISSING, a comma-separated list of
# chain labels (e.g. b0g0_run3) that combine may leave out; each must be in the
# expected set.  Returns the labels left out under that opt-out.
CheckChainSet <- function(root, found, mkv = .cet$mkv,
                          allow = strsplit(Sys.getenv("CET_BF_ALLOW_MISSING"), ",")[[1]]) {
  allow <- trimws(allow[nzchar(trimws(allow))])
  expected <- ExpectedChains(root, mkv)
  found <- found[found$mkv == mkv, ]
  fKey <- ChainKey(found)
  eKey <- ChainKey(expected)
  hasCross <- file.exists(file.path(found$dir, "cross.csv"))
  absent <- expected$label[!eKey %in% fKey]
  noCross <- found$label[fKey %in% eKey & found$power == 1 & !hasCross]
  unknown <- found$label[!fKey %in% eKey]
  missing <- c(absent, noCross)
  badAllow <- setdiff(allow, expected$label)
  cat(sprintf("Expected chains (mkv %s): %d (%s); found %d; missing %d; not in the launch record %d\n",
              mkv, nrow(expected), attr(expected, "source"), nrow(found), length(missing),
              length(unknown)))
  List <- function(what, x) if (length(x)) paste0(what, ": ", paste(x, collapse = ", "))
  problems <- c(
    List("no complete segment (chain directory missing or empty)", setdiff(absent, allow)),
    List("no cross.csv (eval failed or was not run)", setdiff(noCross, allow)),
    List("found but not in the launch record (add to <run root>/launch.csv)", unknown),
    List("CET_BF_ALLOW_MISSING names chains that are not expected", badAllow))
  if (length(problems)) {
    stop("combine: the chain set does not match the launch record; no Bayes factor is reported.\n  ",
         paste(problems, collapse = "\n  "),
         "\n  To leave a chain out on purpose, name it in CET_BF_ALLOW_MISSING.", call. = FALSE)
  }
  leftOut <- intersect(missing, allow)
  if (length(leftOut)) {
    cat(sprintf("Left out on request (CET_BF_ALLOW_MISSING): %s\n", paste(leftOut, collapse = ", ")))
  }
  leftOut
}

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
# Starting parameters: the round-1 total-evidence posterior means (lambda
# 0.20, mu 0.16; psi from the sampling proportion; clock_morph 0.017), each
# multiplied by exp(N(0, 1)) for every chain, replicate 1 included (#89 review,
# R89-5: dispersed starts; before, exp(N(0, 0.2^2)) and none for replicate 1).
# The draw uses its own stream, set.seed(9000 + index), and leaves the chain's
# random-number state untouched.  A draw is rejected, and the next one from the
# same stream taken, while its log posterior is not finite: the log prior, the
# origin age prior, the log augmented density of a tempo-path draw from the
# start state (as TempoMCMC() builds it: .Precompute(), .AddSlices(),
# .Messages(), .SampleTree()) and the character log-likelihood.  Without
# `tree`, no check is made (as before).  Returns the parameters, with
# attributes "attempts" and "logPost".
CetStartLogPost <- function(par, tree, model, characters, priors, nSlice = 20) {
  tryCatch({
    lp <- sum(vapply(names(priors), function(p) priors[[p]](par[[p]]), 0))
    if (!is.finite(lp)) return(-Inf)
    tt <- get(".NoRootAncestor", .ns)(get(".ForModel", .ns)(TempoTree(tree), model))
    la <- CetAgePriors()$origin(get(".StartAge", .ns)(tt))
    set <- get(".AddSlices", .ns)(tt, get(".Precompute", .ns)(tt, par, model), nSlice)
    msg <- get(".Messages", .ns)(tt, set)
    draw <- get(".SampleTree", .ns)(tt, set, msg)
    aug <- get(".LogAugmented", .ns)(tt, set, draw[["paths"]], draw[["bin"]])
    chr <- get("MkLikelihood", .ns)(tt, characters,
                                    get(".CharacterLengths", .ns)(tt, set[["grid"]][["x"]],
                                                                 draw[["paths"]], par[["rate"]],
                                                                 model[["gamma"]]), 2L)
    lp + la + aug + chr
  }, error = function(e) structure(-Inf, error = conditionMessage(e)))
}
CetStart <- function(index, tree = NULL, model = NULL, characters = NULL, priors = CetPriors(),
                     nSlice = 20, sdlog = 1, maxAttempts = 100L) {
  centre <- list(lambda = 0.2, mu = 0.15, psi = 0.05, theta = 0.05, s = 0.2, rate = 0.017)
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(9000 + index)
  lastError <- NULL
  for (attempt in seq_len(maxAttempts)) {
    start <- lapply(centre, function(v) v * exp(stats::rnorm(1, 0, sdlog)))
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
