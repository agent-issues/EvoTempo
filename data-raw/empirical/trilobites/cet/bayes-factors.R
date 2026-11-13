# Trilobite CET pipeline: Bayes factors between the four (beta, gamma) cells.
#
#   Rscript data-raw/empirical/trilobites/cet/bayes-factors.R eval <chain dir>
#   Rscript data-raw/empirical/trilobites/cet/bayes-factors.R combine <run root>
#
# The procedure of the penguin pipeline (../../penguins/cet/bayes-factors.R and
# its README, "Bayes factors"): all four cells share one state space (time
# tree with sampled ancestors, tempo paths and bins, lambda, mu, psi, theta, s,
# rate, origin) and one prior on it, so for cells A and B
#   log q_B(x) - log q_A(x) = [log p(T, paths | Theta, beta_B) - log p(T, paths | Theta, beta_A)]
#                            + [log Mkv(D | T, paths, rate, gamma_B) - log Mkv(D | T, paths, rate, gamma_A)]
# at every recorded state (`eval`), and two-sample bridge sampling (Meng & Wong
# 1996) on the draws of both posteriors estimates log Z_B / Z_A.  Four edges
# (beta at gamma = 0 and 1, gamma at beta = 0 and 1) must close the square;
# ladder rungs (labels such as b1g0.5) are bridged in sequence when present.
# Errors: the bridge relative-MSE SE (Fruhwirth-Schnatter 2004; rung 7 found it
# ~1.6x too small) and the spread over replicate pairs.  ESS is MkPrime's
# (common.R, Ess()).
#
# Differences from the penguin script:
# - `eval` refuses a chain whose checkpoint is missing or records no Mkv
#   variant, or records a grid other than this job's CET_NBIN / CET_ZMAX, and
#   stops when a state's recomputed log posterior, or its recomputed character
#   log-likelihood (which the power-0 check cannot see), differs from the
#   chain's by more than CET_BF_TOL (default 1e-6) (#76; review of #91, A4);
# - `eval` scores the characters in the chain's own tip order
#   (states[[1]]$tree$tipLabel; every chain has its own start tree), and
#   writes cross.csv only after the check passes (temporary file, then
#   rename); `combine` refuses any cross.csv whose recorded checks exceed
#   CET_BF_TOL (review of #91, F1, A1, A2);
# - cross.csv has logJointW, the log importance weight from the CET priors to
#   the comparator's joint (d, nu, s) prior on (lambda, mu, psi)
#   (JointPriorLogWeight()); `combine` reports each edge reweighted to that
#   prior, log BF' = log BF + log E_B[w] - log E_A[w] (self-normalised), with
#   the weights' ESS per cell: MkPrime's bulk ESS (autocorrelation-aware) and
#   the iid Kish ESS (README, "Priors"; review of #91, F3);
# - column names and chain labels write numbers with FormatValue() (common.R):
#   format() on the vector c(0, 0.5, 1) gave aug_b0.0 / aug_b0.5 / aug_b1.0,
#   which combine could not find, so any ladder rung crashed it (#100);
#   ParseLabel() and Chains() are in common.R and stop on a label they cannot
#   read (as the penguin pipeline, #75);
# - combine checks the chains it finds against the expected set (launch
#   record <run root>/launch.csv, or the production default of 4 cells x 4
#   replicates; common.R, ExpectedChains()) and stops if any expected chain has
#   no complete segment or no cross.csv, or if a chain is not in the record;
#   CET_BF_ALLOW_MISSING (chain labels) is the only, explicit, opt-out (#101).
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
burn <- as.numeric(Sys.getenv("CET_BF_BURN", "0.25"))
tol <- as.numeric(Sys.getenv("CET_BF_TOL", "1e-6"))

# ParseLabel(), Chains(), FormatValue() and the expected chain set are in common.R
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
  ckpFile <- file.path(dir, "checkpoint.rds")
  if (!file.exists(ckpFile)) stop("no checkpoint in ", dir, ": cannot verify the chain's target")
  ckp <- readRDS(ckpFile)
  if (is.null(ckp$mkv) || ckp$mkv != me$mkv) {
    stop(sprintf("%s records mkv %s, but its label says %s", ckpFile, format(ckp$mkv), me$mkv))
  }
  grid <- GridSettings()
  if (!identical(ckp$nBin, grid$nBin) || !identical(ckp$zMax, grid$zMax)) {
    stop(sprintf("%s has grid %s bins, z +-%s; this job has %d bins, z +-%g (set CET_NBIN, CET_ZMAX)",
                 ckpFile, format(ckp$nBin), format(ckp$zMax), grid$nBin, grid$zMax))
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
  priors <- CetPriors()
  s <- LoadSamples(dir)
  states <- LoadStates(dir)
  stopifnot(length(states) == nrow(s$samples))
  # Characters in the chain's own tip order (tree moves keep tip indices)
  tipLabel <- states[[1]]$tree$tipLabel
  stopifnot(all(vapply(states, function(st) identical(st$tree$tipLabel, tipLabel), TRUE)))
  dat <- LoadCharacters(tipLabel)
  stopifnot(identical(rownames(dat$characters), tipLabel))
  idx <- which(s$keep)
  t0 <- proc.time()[["elapsed"]]
  out <- do.call(rbind, lapply(idx, function(i) {
    st <- states[[i]]
    aug <- vapply(betas, function(b) LogAugAt(st, b, me$gamma), 0)
    chr <- vapply(gammas, function(g) CharLLAt(st, g, dat$characters), 0)
    lp <- sum(vapply(names(priors), function(p) priors[[p]](st$par[[p]]), 0))
    # The chain's logPost includes the origin age prior (constant on its
    # support, so it cancels in the Bayes factors and is not in logPrior)
    la <- CetAgePriors()$origin(get(".StartAge", .ns)(st$tree))
    ownChar <- CharLLAt(st, me$gamma, dat$characters)
    own <- LogAugAt(st, me$beta, me$gamma) + me$power * ownChar + lp + la
    c(sample = i, setNames(aug, paste0("aug_b", vapply(betas, FormatValue, ""))),
      setNames(chr, paste0("chr_g", vapply(gammas, FormatValue, ""))), logPrior = lp,
      logJointW = JointPriorLogWeight(st$par, priors), check = own - s$samples$logPost[i],
      checkChar = ownChar - s$samples$charLL[i])
  }))
  worst <- max(abs(out[, c("check", "checkChar")]))
  cat(sprintf("%s: %d states, %.1f s; max |own - logPost| = %.3g, max |own - charLL| = %.3g\n",
              basename(dir), nrow(out), proc.time()[["elapsed"]] - t0, max(abs(out[, "check"])),
              max(abs(out[, "checkChar"]))))
  crossFile <- file.path(dir, "cross.csv")
  if (!(worst <= tol)) {
    unlink(crossFile)
    stop(sprintf("recomputed log posterior or charLL differs from the chain's by %.3g > %g; no cross.csv written",
                 worst, tol))
  }
  tmp <- tempfile("cross-", tmpdir = dir, fileext = ".csv")
  utils::write.csv(out, tmp, row.names = FALSE)
  stopifnot(file.rename(tmp, crossFile))
  quit(save = "no")
}

# ---- combine ---------------------------------------------------------------------
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
chAll <- Chains(root)
# Every expected chain must be present with a cross.csv, unless named in
# CET_BF_ALLOW_MISSING (#101); only then are chains without one dropped
leftOut <- CheckChainSet(root, chAll)
ch <- chAll[chAll$mkv == .cet$mkv & !chAll$label %in% leftOut &
              file.exists(file.path(chAll$dir, "cross.csv")), ]
cross <- lapply(ch$dir, function(d) utils::read.csv(file.path(d, "cross.csv"), check.names = FALSE))
names(cross) <- ch$label
# Refuse any cross.csv that did not pass eval's checks (review of #91, A2)
bad <- vapply(cross, function(x) {
  cols <- intersect(c("check", "checkChar"), names(x))
  length(cols) < 2 || !all(abs(as.matrix(x[, cols])) <= tol)
}, TRUE)
if (any(bad)) {
  stop("cross.csv fails eval's checks (or predates checkChar) for: ",
       paste(ch$label[bad], collapse = ", "), "; rerun eval on those chains")
}
LogQ <- function(x, beta, gamma) {
  a <- paste0("aug_b", FormatValue(beta))
  g <- paste0("chr_g", FormatValue(gamma))
  if (!all(c(a, g) %in% names(x))) {
    stop(sprintf("cross.csv has no column %s; rerun eval on every chain after adding a rung",
                 paste(setdiff(c(a, g), names(x)), collapse = " or ")), call. = FALSE)
  }
  x[[a]] + x[[g]] + x$logPrior
}
Draws <- function(p, rep = NULL) {
  sel <- ch$power == 1 & ch$beta == p[1] & ch$gamma == p[2] & (is.null(rep) | ch$rep %in% rep)
  do.call(rbind, cross[ch$label[sel]])
}
Bridge <- function(a, b, rep = NULL) {
  xa <- Draws(a, rep)
  xb <- Draws(b, rep)
  if (is.null(xa) || is.null(xb)) return(NULL)
  BridgeTwo(LogQ(xa, b[1], b[2]) - LogQ(xa, a[1], a[2]),
            LogQ(xb, b[1], b[2]) - LogQ(xb, a[1], a[2]))
}
Edge <- function(a, b) {
  along <- if (a[1] != b[1]) 1 else 2
  pts <- unique(ch[ch$power == 1, c("beta", "gamma")])
  other <- if (along == 1) pts$gamma == a[2] else pts$beta == a[1]
  v <- pts[other, along]
  v <- sort(unique(v[v >= min(a[along], b[along]) & v <= max(a[along], b[along])]))
  if (a[along] > b[along]) v <- rev(v)
  # The path must run from a to b: both end points need chains (penguin #88);
  # possible here only if CET_BF_ALLOW_MISSING left out a whole cell
  if (length(v) < 2L || v[1] != a[along] || v[length(v)] != b[along]) {
    cat(sprintf("Edge %s -> %s: no chains with cross.csv at an end point; edge left out\n",
                CellLabel(a[1], a[2]), CellLabel(b[1], b[2])))
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
  wa <- Draws(a)$logJointW
  wb <- Draws(b)$logJointW
  WByChain <- function(p) {
    sel <- ch$power == 1 & ch$beta == p[1] & ch$gamma == p[2]
    lapply(cross[ch$label[sel]], `[[`, "logJointW")
  }
  ea <- WeightEss(WByChain(a))
  eb <- WeightEss(WByChain(b))
  logBF <- sum(vapply(steps, `[[`, 0, "logBF"))
  data.frame(from = CellLabel(a[1], a[2]), to = CellLabel(b[1], b[2]), points = length(v),
             logBF = logBF,
             se = sqrt(sum(vapply(steps, `[[`, 0, "se")^2)),
             seReplicates = stats::sd(perRep, na.rm = TRUE) / sqrt(sum(is.finite(perRep))),
             minEss = min(vapply(steps, function(s) min(s$essF0, s$essF1), 0)),
             maxOneSidedGap = max(vapply(steps, function(s) abs(s$logBFfrom0 - s$logBFfrom1), 0)),
             nInfinite = sum(vapply(steps, `[[`, 0, "nInfinite")),
             logBFjointPrior = logBF + LogMeanExp(wb) - LogMeanExp(wa),
             weightEssFrom = ea[["mkprime"]], weightEssTo = eb[["mkprime"]],
             weightKishFrom = ea[["kish"]], weightKishTo = eb[["kish"]],
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
  logZ <- c(b0g0 = 0, b1g0 = e[["b0g0 b1g0"]], b0g1 = e[["b0g0 b0g1"]],
            b1g1 = mean(c(e[["b0g0 b1g0"]] + e[["b1g0 b1g1"]], e[["b0g0 b0g1"]] + e[["b0g1 b1g1"]])))
  prob <- exp(logZ - LogSumExp(logZ))
  res <- data.frame(cell = names(logZ), logZ = logZ, posteriorProbability = prob)
  print(res, row.names = FALSE)
  utils::write.csv(res, file.path(root, "model-probabilities.csv"), row.names = FALSE)
  writeLines(sprintf("closure %.4f se %.4f", closure, sqrt(sum(se^2))), file.path(root, "closure.txt"))
}
ch <- chAll[chAll$mkv == .cet$mkv & !chAll$label %in% leftOut, ]
pw <- unique(ch[, c("beta", "gamma")])
for (i in seq_len(nrow(pw))) {
  sel <- ch$beta == pw$beta[i] & ch$gamma == pw$gamma[i]
  ladder <- sort(unique(ch$power[sel]))
  if (length(ladder) < 3 || ladder[1] != 0) next
  logL <- lapply(ladder, function(t) {
    unlist(lapply(ch$dir[sel & ch$power == t], function(d) {
      s <- LoadSamples(d)
      s$samples$charLL[s$keep]
    }))
  })
  ss <- SteppingStone(ladder, logL)
  cat(sprintf("%s: log p(characters | tip ages) = %.3f (SE %.3f), %d rungs\n",
              CellLabel(pw$beta[i], pw$gamma[i]), ss$logR, ss$se, length(ladder)))
}
