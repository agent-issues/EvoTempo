# Trilobite CET pipeline: free-parameter simulation-based calibration (SBC-lite)
# of the all-fossil set-up (review of #91, F4).
#
#   Rscript data-raw/empirical/trilobites/cet/check-sbc.R [nRep] [nIter]
#
# check-window.R holds the parameters fixed, so the parameter sweep, the joint
# sweep and the shift move were not tested against simulation at rho = 0 with
# beta = 1.  This script does, on small trees, at modest scale.
#
# Design: K = 16, z in [-3, 3], beta = 1, delta = 0, rho = 0,
# start = "origin", condition = "survival" (marginal), no characters, point
# fossil ages (as the trilobites: no fossilAge move).  Priors (narrow, to keep
# simulated trees small): lambda ~ LogNormal(log 0.5, 0.25), mu ~
# LogNormal(log 0.4, 0.25), psi ~ LogNormal(log 0.2, 0.3), theta ~
# LogNormal(log 0.2, 0.4), s ~ LogNormal(log 0.6, 0.25); origin age ~ U(3, 8).
# Per replicate:
#   1. draw the parameters and the origin age from the priors;
#   2. draw a start bin and simulate forward to age 0 with rho = 0, repeating
#      at the same parameters and origin age (new start bin each time) until
#      the stem leaves at least one fossil (the survival conditioning);
#   3. keep the tree if it has 3-12 fossil samples, at least two tips and no
#      sampled ancestor at the root; otherwise go back to 1.  This selection
#      depends on the data only, so the posterior given the kept data is
#      unchanged and SBC ranks stay uniform;
#   4. run TempoMCMC() on that tree with free parameters and origin age, from
#      an independent prior draw, moves: tree 1, node 2, tip 2, branch 2,
#      joint 1, shift 1, nodeAge 2, slide 2, saJump 2, originAge 2; nIter
#      iterations after 2000 burn-in, thinned to 99 draws;
#   5. rank of the true value among the 99 draws, for lambda, mu, psi, theta,
#      s and the origin age.
# Test, fixed before any run: for each of the six quantities, chi-square test
# of uniform ranks over 10 bins (ranks 0-99 in bins of 10); Holm-adjusted
# p >= 0.01 for all six passes.  The MkPrime ESS of each replicate's 99 draws
# is reported (thinning adequacy: dependence among draws makes ranks
# non-uniform, so a failure with low ESS points to mixing, not bias).
# Environment: CET_SBC_OUT (default data-raw/empirical/trilobites/cet/sbc-check),
#   EVOTEMPO_CORES (default 16), CET_SBC_SEED (seed offset, default 0).
suppressMessages(library(EvoTempo))
library(parallel)
ns <- asNamespace("EvoTempo")
args <- commandArgs(trailingOnly = TRUE)
nRep <- if (length(args) > 0) as.integer(args[[1]]) else 200L
nIter <- if (length(args) > 1) as.integer(args[[2]]) else 19800L
nCore <- as.integer(Sys.getenv("EVOTEMPO_CORES", "16"))
outDir <- Sys.getenv("CET_SBC_OUT", "data-raw/empirical/trilobites/cet/sbc-check")
seedOffset <- as.integer(Sys.getenv("CET_SBC_SEED", "0"))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
nDraw <- 99L
thin <- nIter %/% nDraw
stopifnot(thin * nDraw == nIter)
burnin <- 2000L
meanlog <- c(lambda = log(0.5), mu = log(0.4), psi = log(0.2), theta = log(0.2), s = log(0.6))
sdlog <- c(lambda = 0.25, mu = 0.25, psi = 0.3, theta = 0.4, s = 0.25)
originLo <- 3
originHi <- 8
nRange <- c(3L, 12L)
model <- DCETModel(nBin = 16, zRange = c(-3, 3), beta = 1, rho = 0, condition = "survival",
                   start = "origin")
moves <- c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1, shift = 1, nodeAge = 2,
           slide = 2, saJump = 2, originAge = 2)

OneReplicate <- function(i) {
  set.seed(91000 + i + seedOffset)
  Setup <- get(".Setup", ns)
  Sample1 <- get(".Sample1", ns)
  priors <- TempoPriors(meanlog, sdlog)
  agePriors <- list(origin = function(a) stats::dunif(a, originLo, originHi, log = TRUE))
  tries <- 0L
  repeat {
    tries <- tries + 1L
    par <- as.list(stats::setNames(stats::rlnorm(5, meanlog, sdlog), names(meanlog)))
    T0 <- stats::runif(1, originLo, originHi)
    set <- Setup(par, model)
    repeat {
      k <- Sample1(set$prior)
      h <- tryCatch(SimulateDCET(T0, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                                 start = k, rho = 0, maxEdges = 5000),
                    error = function(e) NULL)
      if (is.null(h) || nrow(h$fossils) > 0) break
    }
    if (is.null(h)) next
    n <- nrow(h$fossils)
    if (n < nRange[1] || n > nRange[2]) next
    tree <- ReconstructedTree(h, gamma = 1)
    if (is.null(tree)) next
    tt <- TempoTree(tree, present = T0)
    if (tt$sampled[tt$root - tt$nTip] || any(tt$tipType == "extant")) next
    break
  }
  start <- as.list(stats::setNames(stats::rlnorm(5, meanlog, sdlog), names(meanlog)))
  t0 <- proc.time()[["elapsed"]]
  run <- tryCatch(TempoMCMC(tt, model, start, priors, nIter = nIter, burnin = burnin, thin = thin,
                            nSlice = 16, moves = moves, agePriors = agePriors),
                  error = function(e) e)
  if (inherits(run, "error")) {
    return(data.frame(rep = i, error = conditionMessage(run)))
  }
  smp <- run$samples
  truth <- c(unlist(par), originAge = T0)
  q <- c(names(meanlog), "originAge")
  ranks <- vapply(q, function(p) sum(smp[[p]] < truth[[p]]), 0)
  ess <- vapply(q, function(p) {
    x <- smp[[p]]
    if (stats::var(x) == 0) return(NA_real_)
    MkPrime::ConvergenceDiagnostics(structure(
      list(samples = matrix(x, ncol = 1, dimnames = list(NULL, "x")), nRuns = 1L, burnin = 0L),
      class = "MkPosterior"))$ess[[1]]
  }, 0)
  data.frame(rep = i, error = NA_character_, nSamples = tt$nTip, nAncestor = sum(tt$tipType == "ancestor"),
             tries = tries, secs = proc.time()[["elapsed"]] - t0,
             t(stats::setNames(truth, paste0("true_", q))), t(stats::setNames(ranks, paste0("rank_", q))),
             t(stats::setNames(ess, paste0("ess_", q))))
}

cl <- makeCluster(nCore)
invisible(clusterCall(cl, function() {
  suppressMessages(library(EvoTempo))
  NULL
}))
clusterExport(cl, c("ns", "meanlog", "sdlog", "originLo", "originHi", "nRange", "model", "moves",
                    "nIter", "burnin", "thin", "seedOffset", "OneReplicate"))
t0 <- Sys.time()
res <- parLapply(cl, seq_len(nRep), OneReplicate)
stopCluster(cl)
cols <- unique(unlist(lapply(res, names)))
res <- do.call(rbind, lapply(res, function(r) {
  r[setdiff(cols, names(r))] <- NA
  r[cols]
}))
write.csv(res, file.path(outDir, "sbc-ranks.csv"), row.names = FALSE)
ok <- is.na(res$error)
cat(sprintf("%d replicates, %d failed; %s\n", nRep, sum(!ok), format(difftime(Sys.time(), t0, units = "mins"))))
q <- c(names(meanlog), "originAge")
summary <- do.call(rbind, lapply(q, function(p) {
  r <- res[[paste0("rank_", p)]][ok]
  counts <- tabulate(r %/% 10 + 1, 10)
  test <- stats::chisq.test(counts, p = rep(0.1, 10))
  data.frame(quantity = p, n = length(r), chisq = unname(test$statistic), p = test$p.value,
             meanRankOver99 = mean(r) / 99, medianEss = stats::median(res[[paste0("ess_", p)]][ok], na.rm = TRUE),
             bins = paste(counts, collapse = " "))
}))
summary$pHolm <- stats::p.adjust(summary$p, "holm")
summary$pass <- summary$pHolm >= 0.01
print(summary, row.names = FALSE, digits = 3)
cat(sprintf("trees: %.1f samples on average (%d-%d), %.2f sampled ancestors; %.1f s per chain\n",
            mean(res$nSamples[ok]), min(res$nSamples[ok]), max(res$nSamples[ok]),
            mean(res$nAncestor[ok]), mean(res$secs[ok])))
write.csv(summary, file.path(outDir, "sbc-summary.csv"), row.names = FALSE)
cat(sprintf("check-sbc.R: %s (min Holm p %.3g)\n", if (all(summary$pass)) "PASS" else "FAIL",
            min(summary$pHolm)))
