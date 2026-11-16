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
# 1. MkLikelihood(): the characters are the package's MkCharacters() (#97,
#    #98; LoadCharacters()): the comparator's model, Mk (fnJC(K)) by number of
#    states K, characters 45 and 47 ordered (fnOrderedRateMatrix(maxState =
#    K)), polymorphic cells as ambiguous sets, coding = "variable" (one Mkv
#    correction per missing-data mask and partition), and Gamma rates among
#    characters with 4 categories and a free shape alpha ~ Exponential(rate
#    2) (CetCharPriors()).  The rebinding PoweredMk() only multiplies the
#    package's log-likelihood of MkCharacters objects by the power
#    `.cet$power` (stepping-stone chains); at power 1 it returns the
#    package's value unchanged.
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
# Gamma rate categories among characters, as the comparator
# (fnDiscretizeGamma(alpha, alpha, 4))
cetRateCats <- 4L
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

# ---- Character likelihood -----------------------------------------------------
# The package's MkLikelihood() raised to the power `.cet$power` for
# MkCharacters objects; any other call (no characters, a plain matrix) passes
# through unchanged
PoweredMk <- function(tree, characters, lengths, nState = 2L,
                      engine = getOption("EvoTempo.engine", "cpp"), alpha = NULL,
                      power = .cet$power) {
  if (!inherits(characters, "MkCharacters")) {
    return(.mkOriginal(tree, characters, lengths, nState, engine, alpha))
  }
  if (power == 0) {
    return(0)
  }
  ll <- .mkOriginal(tree, characters, lengths, nState, engine, alpha)
  if (power == 1) ll else power * ll
}
utils::assignInNamespace("MkLikelihood", PoweredMk, ns = "EvoTempo")

# Prior on the Gamma shape alpha: the comparator's alpha_morph ~
# dnExponential(2) (mean 0.5; BEAST 1.8)
CetCharPriors <- function() {
  list(alpha = function(a) if (a > 0) stats::dexp(a, rate = 2, log = TRUE) else -Inf)
}

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
      # alphaStep: the Gamma shape's adapted step (#97), fixed after burn-in
      for (v in c("st", "step", "jointStep", "ageStep", "moves", "alphaStep")) {
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
      if (is.null(chars)) 0 else PoweredMk(
        st[["tree"]], chars,
        get(".CharacterLengths", .ns)(st[["tree"]], st[["set"]][["grid"]][["x"]],
                                      st[["paths"]], st[["par"]][["rate"]],
                                      model[["gamma"]]), alpha = st[["par"]][["alpha"]], power = 1)
    }
    .cet$records[[length(.cet$records) + 1L]] <- list(
      newick = TempoNewick(st[["tree"]]), charLL = charLL,
      state = st[c("tree", "par", "paths", "bin")])
    .cet$captured <- list(st = st, step = get("step", pf), jointStep = get("jointStep", pf),
                          ageStep = get("ageStep", pf), moves = get("moves", pf),
                          alphaStep = get("alphaStep", pf),
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
# Production chain k = 4 (rep - 1) + cell (cells 1-4 = b0g0, b0g1, b1g0, b1g1;
# run-array.sh task k) starts from draw k of prepare.R; a ladder rung or power
# chain (cell 0) uses draw <rep>.  Indices above the number of draws wrap.
StartIndex <- function(cell, rep, nDraw = nrow(read.csv(file.path(cetDir, "data", "start-values.csv")))) {
  k <- if (cell >= 1L) 4L * (rep - 1L) + cell else rep
  (k - 1L) %% nDraw + 1L
}
StartTreeFile <- function(k) file.path(cetDir, "data", sprintf("start-tree-%02d.nwk", k))

# Characters as MkCharacters() (rows = tips in the order `tipLabel`, which
# must be the tip order of the tree they are scored on: bayes-factors.R eval
# uses each chain's own; review of #91, F1): K and ordering per character
# from partitions.csv (the comparator's K, from the matrix symbols), cells as
# written by prepare.R ("0", "10", "0/1"; "" unknown), coding "variable" (or
# "all" with CET_MKV=none), `nRateCat` Gamma categories.  The row names are
# the tip labels, which MkLikelihood() checks against the tree.
LoadCharacters <- function(tipLabel, mkv = .cet$mkv, nRateCat = cetRateCats) {
  chars <- read.csv(file.path(cetDir, "data", "characters.csv"), check.names = FALSE,
                    colClasses = "character")
  parts <- read.csv(file.path(cetDir, "data", "partitions.csv"))
  stopifnot(setequal(chars$taxon, tipLabel), !anyDuplicated(tipLabel))
  cells <- as.matrix(chars[match(tipLabel, chars$taxon), -1])
  cells[is.na(cells)] <- ""
  rownames(cells) <- tipLabel
  stopifnot(identical(colnames(cells), paste0("c", parts$char)),
            identical(rownames(cells), tipLabel))
  coding <- switch(mkv, pattern = "variable", none = "all",
                   stop("CET_MKV must be \"pattern\" or \"none\""))
  mk <- MkCharacters(cells, nState = parts$k, ordered = parts$ordered, coding = coding,
                     nRateCat = nRateCat)
  stopifnot(identical(mk$tipLabel, tipLabel), identical(as.integer(mk$nState), as.integer(parts$k)),
            identical(mk$ordered, parts$ordered & parts$k > 2L))
  list(characters = mk, cells = cells, parts = parts)
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
  ch <- LoadCharacters(tt$tipLabel)
  cons <- read.csv(file.path(triDir, "trilobites_constraints.csv"))
  ingroupTaxa <- strsplit(cons$taxa[cons$constraint == "Ingroup_node"], " ")[[1]]
  ingroup <- tt$tipLabel %in% ingroupTaxa
  stopifnot(sum(ingroup) == 104, IngroupMonophyletic(tt, ingroup))
  .cet$ingroup <- ingroup
  list(tree = tt, characters = ch$characters, cells = ch$cells, parts = ch$parts, tips = info,
       ingroup = ingroup)
}

# ---- Model, priors, moves -------------------------------------------------------------
cellTable <- data.frame(beta = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
# A number as it appears in chain labels and cross.csv column names: plain
# decimal notation (never scientific, which ParseLabel() cannot read), with
# the fewest significant digits that read back as exactly `x`.  One number at
# a time: format() on a vector pads to a common width, so c(0, 0.5, 1) became
# "0.0", "0.5", "1.0" and the cross.csv columns of a ladder rung no longer
# matched the names combine looked up (#100; as the penguin pipeline, #75)
FormatValue <- function(x) {
  stopifnot(length(x) == 1L, is.finite(x), x >= 0)
  for (d in 1:17) {
    s <- format(x, scientific = FALSE, digits = d, trim = TRUE)
    if (as.numeric(s) == x) return(s)
  }
  stop("cannot write ", x, " in plain decimal notation")
}
CellLabel <- function(beta, gamma) sprintf("b%sg%s", FormatValue(beta), FormatValue(gamma))

# Chain labels (run-chain.R): b<beta>g<gamma>[_mkvnone][_t<power>]_run<rep>.
# ParseLabel() stops on a label it cannot read; Chains() used to drop such a
# directory without a message (#100; penguin #75).
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
# out of the Bayes factors without a message (#91 delta review, N2: removing one
# cross.csv moved the b0g0 -> b0g1 log BF from 7.13 to 7.54, exit 0).
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
    L <- data.frame(cell = mapply(CellLabel, cellTable$beta, cellTable$gamma),
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
# A chain is complete when the last line of its progress.txt says so
# (run-chain.R writes "<time> complete" after the last segment, and a resumed
# or extended chain writes later lines after any earlier "complete")
ChainComplete <- function(dir) {
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
  f <- sort(Sys.glob(file.path(dir, "seg*_samples.csv")))
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
# Starting values: the RevBayes sample whose tree starts chain k
# (prepare.R, start-values.csv): lambda, mu, psi, rate = clock_morph,
# s = ucld_sd, alpha = alpha_morph (the Gamma shape); theta 0.05 (no
# counterpart)
CetStart <- function(k) {
  sv <- read.csv(file.path(cetDir, "data", "start-values.csv"))
  r <- sv[sv$k == k, ]
  stopifnot(nrow(r) == 1)
  list(lambda = r$lambda, mu = r$mu, psi = r$psi, theta = 0.05,
       s = min(r$ucld_sd, 2.5), rate = r$clock_morph, alpha = r$alpha_morph)
}
# Move counts per iteration, the same for every chain whatever its start tree
# (review of #91, F7): the counts of TempoMoves() and of the tree moves for a
# 107-tip tree WITHOUT sampled ancestors (106 internal nodes, 106
# bifurcations, 212 branches with paths): with gamma != 0, node 27, tip 27,
# branch 53, shift 1; with gamma = 0, TempoMoves()'s fixed schedule (tree 1,
# node 2, tip 2, branch 2, joint 1, shift 1); plus nodeAge 27 (0.25 x
# bifurcations), slide 16 (0.15 x), saJump 11 (0.1 x tips) and originAge 1;
# no fossilAge (point ages).  CET_TREE_MOVES scales the tree moves.
CetMoves <- function(tt, characters, model,
                     scale = as.numeric(Sys.getenv("CET_TREE_MOVES", "1"))) {
  nTip <- tt$nTip
  nBif <- nTip - 1L
  path <- c(if (model[["gamma"]] != 0 && !is.null(characters)) {
    c(node = ceiling(0.25 * (nTip - 1L)), tip = ceiling(0.25 * nTip),
      branch = ceiling(0.25 * (2L * nTip - 2L)), shift = 1)
  } else {
    TempoMoves(tt, NULL, model)
  })
  c(path, nodeAge = ceiling(scale * 0.25 * nBif), slide = ceiling(scale * 0.15 * nBif),
    saJump = ceiling(scale * 0.1 * nTip), originAge = 1)
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
# Effective sample size of importance weights given on the log scale, one
# vector per chain: MkPrime's bulk ESS (rank-normalised split chains, which
# accounts for autocorrelation) of the normalised weights, tail-equalised over
# chains (review of #91, F3); and the iid Kish ESS, sum(w)^2 / sum(w^2)
WeightEss <- function(lwByChain) {
  lw <- unlist(lwByChain)
  fin <- is.finite(lw)
  if (!any(fin)) return(c(mkprime = 0, kish = 0))
  m <- max(lw[fin])
  w <- lapply(lwByChain, function(x) ifelse(is.finite(x), exp(x - m), 0))
  n <- min(lengths(w))
  mat <- vapply(w, function(x) utils::tail(x, n), numeric(n))
  if (n == 1L) mat <- matrix(mat, nrow = 1)
  wa <- unlist(w)
  c(mkprime = if (n >= 4 && stats::var(as.vector(mat)) > 0) MkPrime:::.EssBulk(mat) else NA_real_,
    kish = sum(wa)^2 / sum(wa^2))
}
# ESS of one sequence (MkPrime, as above)
Ess <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || stats::var(x) == 0) return(NA_real_)
  RhatEss(list(matrix(x, ncol = 1, dimnames = list(NULL, "x"))))$ess[[1]]
}
