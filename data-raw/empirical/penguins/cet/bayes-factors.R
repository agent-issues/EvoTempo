# Penguin CET pipeline: Bayes factors between the four (beta, gamma) cells.
#
#   Rscript data-raw/empirical/penguins/cet/bayes-factors.R eval <chain dir>
#   Rscript data-raw/empirical/penguins/cet/bayes-factors.R combine <run root>
#
# Procedure (README.md, "Bayes factors"): rung 7 (data-raw/rung7-power.R)
# compared cells on a fixed true tree, combining tree-only evidence by
# quadrature, a character-power stepping-stone ladder, and a two-sample bridge
# between the (0,1) and (1,1) posteriors.  With tree uncertainty the first two
# have no counterpart (the tree-only evidence would integrate over every
# labelled time tree), but the bridge does: all four cells share one state
# space x = (time tree, tip ages, tempo paths and bins, lambda, mu, psi, theta,
# s, rate) and one prior on it, so for cells A and B
#   log q_B(x) - log q_A(x) = [log p(T, paths | Theta, beta_B)
#                               - log p(T, paths | Theta, beta_A)]
#                            + [log Mkv(D | T, paths, rate, gamma_B)
#                               - log Mkv(D | T, paths, rate, gamma_A)],
# both terms computable at any recorded state, and two-sample bridge sampling
# (Meng & Wong 1996; BridgeTwo() as rung 7) on draws from both posteriors
# estimates log Z_B / Z_A, the Bayes factor for all the data (characters, tip
# ages, sampled ancestors) with the tree integrated out by the MCMC.
#
# - Edges: beta at gamma = 0 and 1; gamma at beta = 0 and 1.  The square
#   closes: e(00->10) + e(10->11) - e(00->01) - e(01->11) = 0 up to error, a
#   check that no single edge can give.
# - Error: the iterative-bridge relative-MSE SE (Fruhwirth-Schnatter 2004),
#   which rung 7 found ~1.6x too small, and the spread of the edge estimated
#   from each replicate pair separately (rep k of A with rep k of B).
# - Overlap: if an edge's bridge weights have ESS < 100, or the two one-sided
#   importance estimates differ by > 2 nats, or the square fails to close by
#   > 3 SE, run ladder rungs at intermediate model values along that edge
#   (DCETModel() takes numeric beta and gamma; run-chain.R cell label such as
#   b1g0.5).  combine then bridges consecutive rungs and sums them (a
#   model-index stepping-stone path).  Rungs found under <run root> are used
#   automatically.
# - Character-power chains (CET_POWER = t, labels <cell>_t<t>), if present,
#   give log p(characters | tip ages, cell) by stepping-stone, the quantity
#   that RevBayes's powerPosterior() estimates for the morphology-only UCLN
#   comparator; that is the only cross-framework comparison (the RevBayes and
#   EvoTempo tree densities differ in conventions, so full-data evidences are
#   not comparable across programs).
#
# `eval` writes <chain dir>/cross.csv: for every recorded state (after
# dropping CET_BF_BURN, default 0.25, of each chain) the log augmented density
# at each beta and the character log-likelihood at each gamma that occur among
# the chains under the same root (CET_BF_BETAS / CET_BF_GAMMAS override), and a
# check that the chain's own values reproduce its logPost.
#
# `combine` checks the chains it finds against the expected set: the launch
# record <run root>/launch.csv (cell, reps, optional power and mkv; one row per
# launched group, e.g. each ladder rung), or without one the production launch
# of 4 cells x CET_BF_REPS (default 4) replicates.  It stops if an expected
# chain has no complete segment, if an expected power-1 chain has no cross.csv,
# or if it finds a chain the record does not list; CET_BF_ALLOW_MISSING (chain
# labels) is the only, explicit, opt-out (#101).
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/penguins/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
burn <- as.numeric(Sys.getenv("CET_BF_BURN", "0.25"))

# ParseLabel() and Chains() are in common.R (shared with run-chain.R, #75)
LoadSamples <- function(dir) {
  smp <- do.call(rbind, lapply(sort(Sys.glob(file.path(dir, "seg*_samples.csv"))), utils::read.csv))
  keep <- seq_len(nrow(smp)) > floor(burn * nrow(smp))
  list(samples = smp, keep = keep)
}
LoadStates <- function(dir) {
  unlist(lapply(sort(Sys.glob(file.path(dir, "seg*_states.rds"))), readRDS), recursive = FALSE)
}
LogAugAt <- function(state, beta, gamma) {
  model <- CetModel(beta, gamma)
  tryCatch(get(".LogAugmented", .ns)(state$tree, get(".Precompute", .ns)(state$tree, state$par, model),
                                     state$paths, state$bin),
           error = function(e) -Inf)
}
CharLLAt <- function(state, gamma, characters) {
  model <- CetModel(0, gamma)
  x <- TempoGrid(model$nBin, model$zRange, state$par$s)$x
  PartitionedMk(state$tree, characters,
                get(".CharacterLengths", .ns)(state$tree, x, state$paths, state$par$rate, gamma),
                power = 1)
}

if (args[[1]] == "eval") {
  dir <- args[[2]]
  me <- ParseLabel(basename(dir))
  if (is.null(me)) stop("not a chain directory of run-chain.R: ", dir)
  # Evaluate under the chain's own Mkv variant (its label, and its checkpoint
  # if present), not this job's CET_MKV (#74)
  if (nzchar(Sys.getenv("CET_MKV")) && Sys.getenv("CET_MKV") != me$mkv) {
    stop(sprintf("CET_MKV is %s, but chain %s has mkv %s", Sys.getenv("CET_MKV"),
                 basename(dir), me$mkv))
  }
  ckpFile <- file.path(dir, "checkpoint.rds")
  ckp <- if (file.exists(ckpFile)) readRDS(ckpFile)
  if (!is.null(ckp$mkv) && ckp$mkv != me$mkv) {
    stop(sprintf("%s records mkv %s, but its label says %s", ckpFile, ckp$mkv, me$mkv))
  }
  .cet$mkv <- me$mkv
  all <- Chains(dirname(dir))
  all <- all[all$mkv == me$mkv, ]
  betas <- if (nzchar(Sys.getenv("CET_BF_BETAS"))) {
    as.numeric(strsplit(Sys.getenv("CET_BF_BETAS"), ",")[[1]])
  } else sort(unique(all$beta))
  gammas <- if (nzchar(Sys.getenv("CET_BF_GAMMAS"))) {
    as.numeric(strsplit(Sys.getenv("CET_BF_GAMMAS"), ",")[[1]])
  } else sort(unique(all$gamma))
  # The chain's own start tree (run-chain.R records it): the character rows
  # follow its tip order, and chains start from different trees (#89 review)
  startTree <- if (is.null(ckp$startTree)) file.path(cetDir, "data", "start-tree.nwk") else ckp$startTree
  if (!is.null(ckp$startTreeMd5) && unname(tools::md5sum(startTree)) != ckp$startTreeMd5) {
    stop(sprintf("%s has changed since chain %s started from it", startTree, basename(dir)))
  }
  dat <- LoadPenguins(startTree)
  priors <- CetPriors()
  s <- LoadSamples(dir)
  states <- LoadStates(dir)
  stopifnot(length(states) == nrow(s$samples))
  sameTips <- vapply(states, function(st) identical(st$tree$tipLabel, rownames(dat$characters)), TRUE)
  if (!all(sameTips)) {
    stop(sprintf("%s: %d of %d states have a tip order other than the character rows of %s",
                 basename(dir), sum(!sameTips), length(states), startTree))
  }
  idx <- which(s$keep)
  t0 <- proc.time()[["elapsed"]]
  out <- do.call(rbind, lapply(idx, function(i) {
    st <- states[[i]]
    aug <- vapply(betas, function(b) LogAugAt(st, b, me$gamma), 0)
    chr <- vapply(gammas, function(g) CharLLAt(st, g, dat$characters), 0)
    lp <- sum(vapply(names(priors), function(p) priors[[p]](st$par[[p]]), 0))
    # The chain's logPost includes the origin age prior (TempoMCMC() adds the
    # prior of a free start age since #56; CetMoves() sets originAge = 1).  It
    # is the same for every model point, so it cancels in the Bayes factors
    # and is not part of logPrior; the check needs it (#67).  logPrior means
    # the same in the chain's seg*_samples.csv, which has the origin age prior
    # as a separate column, logOriginPrior (#74).  Under U(0, 160) that prior
    # is a constant on its support, so this check cannot tell at which age
    # TempoMCMC() evaluated it.
    la <- CetAgePriors()$origin(get(".StartAge", .ns)(st$tree))
    own <- LogAugAt(st, me$beta, me$gamma) + me$power * CharLLAt(st, me$gamma, dat$characters) +
      lp + la
    c(sample = i, setNames(aug, paste0("aug_b", vapply(betas, FormatValue, ""))),
      setNames(chr, paste0("chr_g", vapply(gammas, FormatValue, ""))), logPrior = lp,
      check = own - s$samples$logPost[i])
  }))
  utils::write.csv(out, file.path(dir, "cross.csv"), row.names = FALSE)
  cat(sprintf("%s: %d states, %.1f s; max |own - logPost| = %.3g\n", basename(dir), nrow(out),
              proc.time()[["elapsed"]] - t0, max(abs(out[, "check"]))))
  quit(save = "no")
}

# ---- combine ---------------------------------------------------------------------
# Two-sample bridge for log(Z1 / Z0) from l = log(q1 / q0) at draws from
# posterior 0 (l0) and posterior 1 (l1); as rung 7
BridgeTwo <- function(l0, l1) {
  n0 <- length(l0)
  n1 <- length(l1)
  s0 <- n0 / (n0 + n1)
  s1 <- n1 / (n0 + n1)
  shift <- stats::median(c(l0[is.finite(l0)], l1[is.finite(l1)]))
  e0 <- exp(l0 - shift)
  e1 <- exp(l1 - shift)
  r <- 1
  for (it in seq_len(5000)) {
    rOld <- r
    r <- mean(e0 / (s1 * e0 + s0 * r)) / mean(1 / (s1 * e1 + s0 * r))
    if (abs(log(r) - log(rOld)) < 1e-12) break
  }
  f0 <- e0 / (s1 * e0 + s0 * r)
  f1 <- 1 / (s1 * e1 + s0 * r)
  re2 <- stats::var(f0) / (Ess(f0) * mean(f0)^2) + stats::var(f1) / (Ess(f1) * mean(f1)^2)
  list(logBF = log(r) + shift, se = sqrt(re2), essF0 = Ess(f0), essF1 = Ess(f1),
       logBFfrom0 = LogMeanExp(l0), logBFfrom1 = -LogMeanExp(-l1),
       nInfinite = sum(!is.finite(c(l0, l1))))
}
SteppingStone <- function(ladder, logL) {
  K <- length(ladder) - 1L
  inc <- vapply(seq_len(K), function(k) {
    x <- (ladder[k + 1] - ladder[k]) * logL[[k]]
    w <- exp(x - max(x))
    c(LogMeanExp(x), stats::var(w) / (Ess(w) * mean(w)^2))
  }, numeric(2))
  inc[2, !is.finite(inc[2, ])] <- 0
  list(logR = sum(inc[1, ]), se = sqrt(sum(inc[2, ])))
}

root <- args[[2]]
stopifnot("combine: <run root> does not exist" = dir.exists(root))
# Every cell with character-power chains of this variant must have a complete
# ladder, t = 0 to t = 1, before anything is combined.  The t = 1 rung is the
# cell's ordinary chains: without them under this root, the stepping-stone
# would estimate log Z(t_max) / Z(0) and report it as the marginal (#88)
LadderProblems <- function(ch) {
  pw <- unique(ch[ch$power != 1, c("beta", "gamma")])
  unlist(lapply(seq_len(nrow(pw)), function(i) {
    sel <- ch$beta == pw$beta[i] & ch$gamma == pw$gamma[i]
    ladder <- sort(unique(ch$power[sel]))
    missing <- c(if (ladder[1] != 0) "t = 0", if (ladder[length(ladder)] != 1) "t = 1")
    few <- length(ladder) < 3L
    if (!length(missing) && !few) return(NULL)
    sprintf("cell %s (mkv %s): power ladder %s %s", CellLabel(pw$beta[i], pw$gamma[i]), .cet$mkv,
            paste(vapply(ladder, FormatValue, ""), collapse = ", "),
            paste(c(if (length(missing)) paste("lacks", paste(missing, collapse = " and ")),
                    if (few) "has fewer than 3 rungs"), collapse = " and "))
  }))
}
chAll <- Chains(root)
problems <- LadderProblems(chAll[chAll$mkv == .cet$mkv, ])
if (length(problems)) {
  stop("combine: incomplete character-power ladder(s); no marginal is reported.\n  ",
       paste(problems, collapse = "\n  "),
       "\n  The t = 1 rung is the cell's chains at power 1 (production chains) under the same root.",
       call. = FALSE)
}
# Every expected chain (launch record <run root>/launch.csv, or the production
# default of 4 cells x 4 replicates; common.R, ExpectedChains()) must be present,
# and every expected power-1 chain must have a cross.csv, unless named in
# CET_BF_ALLOW_MISSING; a chain not in the record also stops combine.  Before,
# a chain without cross.csv was dropped from the bridges without a message (#101)
leftOut <- CheckChainSet(root, chAll)
# Bridge only chains of one Mkv variant: CET_MKV, default "pattern" (#74)
ch <- chAll[chAll$mkv == .cet$mkv & !chAll$label %in% leftOut &
              file.exists(file.path(chAll$dir, "cross.csv")), ]
cross <- lapply(ch$dir, function(d) utils::read.csv(file.path(d, "cross.csv"), check.names = FALSE))
names(cross) <- ch$label
LogQ <- function(x, beta, gamma) {
  x[[paste0("aug_b", FormatValue(beta))]] + x[[paste0("chr_g", FormatValue(gamma))]] + x$logPrior
}
# Bridge between model points A and B (beta, gamma), pooling replicates (or
# replicate `rep` only)
Bridge <- function(a, b, rep = NULL) {
  Draws <- function(p) {
    sel <- ch$power == 1 & ch$beta == p[1] & ch$gamma == p[2] & (is.null(rep) | ch$rep %in% rep)
    do.call(rbind, cross[ch$label[sel]])
  }
  xa <- Draws(a)
  xb <- Draws(b)
  if (is.null(xa) || is.null(xb)) return(NULL)
  BridgeTwo(LogQ(xa, b[1], b[2]) - LogQ(xa, a[1], a[2]),
            LogQ(xb, b[1], b[2]) - LogQ(xb, a[1], a[2]))
}
# An edge from a to b along beta (or gamma), through any ladder rungs present
Edge <- function(a, b) {
  along <- if (a[1] != b[1]) 1 else 2
  pts <- unique(ch[ch$power == 1, c("beta", "gamma")])
  other <- if (along == 1) pts$gamma == a[2] else pts$beta == a[1]
  v <- pts[other, along]
  v <- sort(unique(v[v >= min(a[along], b[along]) & v <= max(a[along], b[along])]))
  if (a[along] > b[along]) v <- rev(v)
  # The path must run from a to b: both end points need chains (#88)
  if (length(v) < 2L || v[1] != a[along] || v[length(v)] != b[along]) {
    cat(sprintf("Edge %s -> %s: no chains with cross.csv at %s; edge left out\n",
                CellLabel(a[1], a[2]), CellLabel(b[1], b[2]),
                paste(c(if (!length(v) || v[1] != a[along]) CellLabel(a[1], a[2]),
                        if (!length(v) || v[length(v)] != b[along]) CellLabel(b[1], b[2])),
                      collapse = " and ")))
    return(NULL)
  }
  P <- function(val) if (along == 1) c(val, a[2]) else c(a[1], val)
  steps <- lapply(seq_len(length(v) - 1L), function(k) Bridge(P(v[k]), P(v[k + 1])))
  if (any(vapply(steps, is.null, TRUE)) || !length(steps)) return(NULL)
  reps <- sort(unique(ch$rep))
  perRep <- vapply(reps, function(r) {
    s <- lapply(seq_len(length(v) - 1L), function(k) Bridge(P(v[k]), P(v[k + 1]), r))
    if (any(vapply(s, is.null, TRUE))) NA else sum(vapply(s, `[[`, 0, "logBF"))
  }, 0)
  data.frame(from = CellLabel(a[1], a[2]), to = CellLabel(b[1], b[2]), points = length(v),
             logBF = sum(vapply(steps, `[[`, 0, "logBF")),
             se = sqrt(sum(vapply(steps, `[[`, 0, "se")^2)),
             seReplicates = stats::sd(perRep, na.rm = TRUE) / sqrt(sum(is.finite(perRep))),
             minEss = min(vapply(steps, function(s) min(s$essF0, s$essF1), 0)),
             maxOneSidedGap = max(vapply(steps, function(s) abs(s$logBFfrom0 - s$logBFfrom1), 0)),
             nInfinite = sum(vapply(steps, `[[`, 0, "nInfinite")),
             perReplicate = paste(format(perRep, digits = 5), collapse = " "))
}
edges <- rbind(Edge(c(0, 0), c(1, 0)), Edge(c(0, 1), c(1, 1)),
               Edge(c(0, 0), c(0, 1)), Edge(c(1, 0), c(1, 1)))
if (is.null(edges)) {
  cat("No bridge edge has chains with cross.csv at both ends\n")
} else {
  print(edges, row.names = FALSE)
  utils::write.csv(edges, file.path(root, "bayes-factors.csv"), row.names = FALSE)
}
if (!is.null(edges) && nrow(edges) == 4) {
  e <- setNames(edges$logBF, paste(edges$from, edges$to))
  se <- edges$se
  closure <- e[["b0g0 b1g0"]] + e[["b1g0 b1g1"]] - e[["b0g0 b0g1"]] - e[["b0g1 b1g1"]]
  cat(sprintf("Square closure: %.3f nats (SE %.3f)\n", closure, sqrt(sum(se^2))))
  # log Z relative to (0,0); (1,1) as the mean of its two paths
  logZ <- c(b0g0 = 0, b1g0 = e[["b0g0 b1g0"]], b0g1 = e[["b0g0 b0g1"]],
            b1g1 = mean(c(e[["b0g0 b1g0"]] + e[["b1g0 b1g1"]], e[["b0g0 b0g1"]] + e[["b0g1 b1g1"]])))
  prob <- exp(logZ - LogSumExp(logZ))
  res <- data.frame(cell = names(logZ), logZ = logZ, posteriorProbability = prob)
  print(res, row.names = FALSE)
  utils::write.csv(res, file.path(root, "model-probabilities.csv"), row.names = FALSE)
  writeLines(sprintf("closure %.4f se %.4f", closure, sqrt(sum(se^2))), file.path(root, "closure.txt"))
}
# Character-power stepping-stone per cell, if power chains exist.  Only chains
# of the Mkv variant CET_MKV (default "pattern"), as for the bridges above:
# pooling the charLL of chains of two variants mixes two targets (#75).
# Ladders were checked to run from t = 0 to t = 1 at the start of combine.
ch <- chAll[chAll$mkv == .cet$mkv & !chAll$label %in% leftOut, ]
pw <- unique(ch[ch$power != 1, c("beta", "gamma")])
for (i in seq_len(nrow(pw))) {
  sel <- ch$beta == pw$beta[i] & ch$gamma == pw$gamma[i]
  ladder <- sort(unique(ch$power[sel]))
  stopifnot(ladder[1] == 0, ladder[length(ladder)] == 1, length(ladder) >= 3L)
  logL <- lapply(ladder, function(t) {
    unlist(lapply(ch$dir[sel & ch$power == t], function(d) {
      s <- LoadSamples(d)
      s$samples$charLL[s$keep]
    }))
  })
  ss <- SteppingStone(ladder, logL)
  cat(sprintf("%s (mkv %s): log p(characters | tip ages) = %.3f (SE %.3f), %d rungs, %d chains\n",
              CellLabel(pw$beta[i], pw$gamma[i]), .cet$mkv, ss$logR, ss$se, length(ladder),
              sum(sel)))
}
