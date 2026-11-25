# Penguin CET: common-footing marginal likelihoods (README, "Common-footing
# marginal likelihoods").
#
#   Rscript data-raw/empirical/penguins/cet/ml-summary.R <CET ladder root> \
#     <RevBayes stepping-stone root> <RevBayes posterior root> <bayes-factors.csv> <out dir>
#
# <CET ladder root>: ml-root.sh output (ladder chains <cell>_t<t>_run<rep> from
#   ml-ladder.sh; links to the cell's production chains, t = 1).
# <RevBayes stepping-stone root>: <clock>_run<rep>/ from run-hamilton.sh with
#   ML=true ML_MODE=stone CATS=65 (output/*_pp_stone_<k>.out; stone 1 has
#   power 1, stone CATS power 0; monitors of stones 1 and CATS).
# <RevBayes posterior root>: the comparator's MCMC chains (<clock>_run<rep>/).
# <bayes-factors.csv>: `bayes-factors.R combine` output for the production
#   cells (edges b0g0 -> b0g1 and b1g0 -> b1g1).
#
# What each side estimates.  Both anneal only the character likelihood: CET
# raises PartitionedMk() to the power t (common.R), RevBayes's powerPosterior()
# only the clamped phy_morph nodes.  The tree (FBD with tip-age intervals),
# clock and FBD parameters stay at their prior, so each estimates
#   log p(characters | tip-age intervals, model) = log E_prior[ L(characters) ],
# the expectation under that program's normalised prior.  Normalising constants
# of the tree density (-log n!, the 1 / width of tip-age intervals) do not
# enter; the 2^|V| orientation factor varies with the number of sampled
# ancestors, so it does not cancel, but both programs carry it identically
# (R/MCMC.R; RevBayes lnProbTreeShape).  The programs' priors differ in one known
# respect: EvoTempo does not allow a sampled ancestor at the root, RevBayes
# does.  If the priors agree on every other tree, CET's prior is RevBayes's
# conditioned on A = {no sampled ancestor at the root}, and
#   log Z_CET = log Z_RB + log P_RB(A | characters) - log P_RB(A),
# with P_RB(A) from the power-0 stones and P_RB(A | characters) from the
# comparator's posterior.  The anchor checks this for CET (0, 0) against
# RevBayes strict; the remaining difference, CET's truncation of its priors
# (+-4 sdlog on lambda, mu, psi and rate, 6e-5 of prior mass each; s on (0, 3];
# theta has no RevBayes counterpart), is below 1e-3 nats.
#
# Estimators, the same for both programs: stepping-stone (Xie et al. 2011;
# SteppingStone() of bayes-factors.R, with its per-stone variance) per
# replicate and pooled over replicates (SE of the pooled CET estimate:
# max(estimator SE, sd / sqrt(n)) over replicates, CetPooledSe()), and path sampling (trapezoid rule over
# the mean character log-likelihood at each power).  Burn-in: the first 25%
# of each CET chain's samples (bayes-factors.R); RevBayes drops 25% of each
# stone itself (runOneStone(burninFraction = 0.25)).
#
# Outputs: ml-stones.csv (per program, cell or clock, replicate and power: mean,
# sd and ESS of the character log-likelihood), ml-estimates.csv, root-sa.csv,
# prior-check.csv, anchor.txt, common-footing.csv.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 5L)
cetRoot <- args[[1]]; rbmlRoot <- args[[2]]; rbRoot <- args[[3]]; bfFile <- args[[4]]
outDir <- args[[5]]
burn <- 0.25
source("data-raw/empirical/penguins/cet/ml-helpers.R")
cats <- as.integer(Sys.getenv("ML_CATS", "65"))
alpha <- as.numeric(Sys.getenv("ML_ALPHA", "0.2"))
partial <- nzchar(Sys.getenv("ML_ALLOW_PARTIAL"))
suppressMessages(library(ape))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
tips <- utils::read.csv("data-raw/empirical/penguins/penguins_tips.csv")
extant <- tips$taxon[tips$status == "extant"]
fossils <- tips$taxon[tips$status == "fossil"]

# ---- Estimators (as bayes-factors.R) ---------------------------------------------
LogSumExp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
LogMeanExp <- function(x) LogSumExp(x) - log(length(x))
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
# logL: list over powers t_1 < ... < t_K (t_K = 1 not used); each element a
# list of per-chain sample vectors (ESS per chain, then summed)
SteppingStone <- function(ladder, logL) {
  K <- length(ladder) - 1L
  inc <- vapply(seq_len(K), function(k) {
    x <- lapply(logL[[k]], function(v) (ladder[k + 1] - ladder[k]) * v)
    all <- unlist(x)
    m <- max(all)
    w <- lapply(x, function(v) exp(v - m))
    ess <- sum(vapply(w, Ess, 0), na.rm = TRUE)
    c(LogMeanExp(all), stats::var(unlist(w)) / (max(ess, 1) * mean(unlist(w))^2))
  }, numeric(2))
  inc[2, !is.finite(inc[2, ])] <- 0
  list(logZ = sum(inc[1, ]), se = sqrt(sum(inc[2, ])))
}
PathSampling <- function(ladder, meanL) sum(diff(ladder) * (head(meanL, -1) + tail(meanL, -1)) / 2)

# ---- Trees ------------------------------------------------------------------------
TreeFacts <- function(tr) {
  nTip <- length(tr$tip.label)
  d <- node.depth.edgelength(tr)
  age <- max(d[match(extant, tr$tip.label)]) - d
  tipEdge <- match(seq_len(nTip), tr$edge[, 2])
  sa <- tr$edge.length[tipEdge] < 1e-8 & tr$tip.label %in% fossils
  c(rootSA = any(sa & tr$edge[tipEdge, 1] == nTip + 1L), nSA = sum(sa),
    crownAge = age[getMRCA(tr, extant)], rootAgeTree = age[nTip + 1L])
}
Facts <- function(nwk) {
  if (!length(nwk)) return(NULL)
  trs <- read.tree(text = nwk)
  if (inherits(trs, "phylo")) trs <- list(trs)
  as.data.frame(do.call(rbind, lapply(trs, TreeFacts)))
}
Burn <- function(n) seq_len(n) > floor(burn * n)

# ---- CET ladders ------------------------------------------------------------------
# list.files, not list.dirs: a dangling link to a production chain must be seen
dirs <- list.files(cetRoot, full.names = TRUE)
lab <- basename(dirs)
m <- regmatches(lab, regexec("^(b[0-9.]+g[0-9.]+)(_t([0-9.]+))?_run([0-9]+)$", lab))
ok <- lengths(m) > 0
cet <- data.frame(dir = dirs[ok], cell = vapply(m[ok], `[`, "", 2),
                  power = vapply(m[ok], function(x) if (nzchar(x[4])) as.numeric(x[4]) else 1, 0),
                  rep = as.integer(vapply(m[ok], `[`, "", 5)))
cet$exists <- dir.exists(cet$dir)
cet$complete <- vapply(file.path(cet$dir, "progress.txt"), function(f) {
  if (!file.exists(f)) return(FALSE)
  x <- readLines(f, warn = FALSE)
  length(x) > 0 && grepl("complete[[:space:]]*$", x[length(x)])
}, TRUE)
LoadCet <- function(dir) {
  segs <- sort(Sys.glob(file.path(dir, "seg*_samples.csv")))
  if (!length(segs)) return(NULL)
  smp <- do.call(rbind, lapply(segs, utils::read.csv))
  smp[Burn(nrow(smp)), ]
}
expectPow <- (seq(0, cats - 2) / (cats - 1))^(1 / alpha)
stones <- list(); est <- list(); cetPrior <- NULL
ladderCells <- unique(cet$cell[cet$power != 1])
problems <- CetLadderProblems(cet[cet$cell %in% ladderCells, ], expectPow)
if (length(problems)) {
  msg <- paste0("CET ladder(s) incomplete:\n  ", paste(problems, collapse = "\n  "))
  if (partial) warning(msg, call. = FALSE) else stop(msg, call. = FALSE)
}
cet <- cet[cet$exists, ]
for (cell in ladderCells) {
  sel <- cet[cet$cell == cell, ]
  ladderReps <- sort(unique(sel$rep[sel$power != 1]))
  data <- lapply(seq_len(nrow(sel)), function(i) LoadCet(sel$dir[i]))
  have <- !vapply(data, is.null, TRUE)
  sel <- sel[have, ]; data <- data[have]
  for (i in seq_len(nrow(sel))) {
    L <- data[[i]]$charLL
    stones[[length(stones) + 1]] <- data.frame(program = "CET", model = cell, rep = sel$rep[i],
      power = sel$power[i], n = length(L), meanL = mean(L), sdL = stats::sd(L), essL = Ess(L))
  }
  pw <- sort(unique(sel$power))
  # Every replicate's ladder ends in all of the cell's production chains at t = 1.
  # The stepping-stone never uses the t = 1 draws, but path sampling does, so the
  # per-replicate path-sampling estimates are correlated (review of #111, N-1).
  repSets <- c(as.list(ladderReps), if (length(ladderReps) > 1) list(ladderReps))
  for (r in repSets) {
    lad <- lapply(pw, function(t) lapply(which(sel$power == t & (sel$power == 1 | sel$rep %in% r)),
                                        function(i) data[[i]]$charLL))
    keep <- lengths(lad) > 0
    ss <- SteppingStone(pw[keep], lad[keep])
    meanL <- vapply(lad[keep], function(x) mean(unlist(x)), 0)
    est[[length(est) + 1]] <- data.frame(program = "CET", model = cell,
      rep = if (length(r) > 1 || length(ladderReps) == 1) "pooled" else as.character(r), stones = sum(keep) - 1L,
      logZ_ss = ss$logZ, se_ss = ss$se, logZ_ps = PathSampling(pw[keep], meanL),
      logZ_ss_even = with(list(ev = EvenPowers(pw[keep], expectPow)),
                          SteppingStone(pw[keep][ev], lad[keep][ev])$logZ))
  }
  # Prior samples (t = 0) of this cell, for the prior check
  for (i in which(sel$power == 0)) {
    s <- data[[i]]
    nwk <- unlist(lapply(sort(Sys.glob(file.path(sel$dir[i], "seg*_trees.nwk"))), readLines))
    nwk <- nwk[Burn(length(nwk))]
    f <- Facts(nwk)
    cetPrior <- rbind(cetPrior, data.frame(source = paste0("CET ", cell, " t = 0"), rep = sel$rep[i],
      lambda = s$lambda, mu = s$mu, psi = s$psi, originAge = s$originAge, rootAge = s$rootAge,
      crownAge = f$crownAge, nSA = f$nSA, rootSA = f$rootSA))
  }
}

# ---- RevBayes stones --------------------------------------------------------------
ReadTab <- function(f) {
  x <- readLines(f)
  hdr <- strsplit(x[1], "\t")[[1]]
  body <- x[-1]
  body <- body[lengths(strsplit(body, "\t")) >= length(hdr) - 1L]
  utils::read.delim(text = c(x[1], body), quote = "", comment.char = "", check.names = FALSE)
}
rbPow <- ((cats - seq_len(cats)) / (cats - 1))^(1 / alpha)
rbPrior <- NULL; rbStone1 <- NULL
for (d in sort(Sys.glob(file.path(rbmlRoot, "*_run*")))) {
  clock <- sub("_run.*", "", basename(d)); r <- as.integer(sub(".*_run", "", basename(d)))
  stem <- file.path(d, "output", sprintf("penguins_morph_%s_run%d", clock, r))
  files <- sprintf("%s_pp_stone_%d.out", stem, seq_len(cats))
  have <- file.exists(files)
  if (!all(have) && !partial) stop(sprintf("RevBayes %s run %d: %d stones missing", clock, r, sum(!have)))
  if (!any(have)) next
  L <- lapply(files[have], function(f) ReadTab(f)$likelihood)
  pw <- rbPow[have]
  filePow <- vapply(files[have], function(f) ReadTab(f)$power[1], 0)
  stopifnot("stone powers differ from the schedule" = all(abs(filePow - pw) <= 1e-5 * pmax(pw, 1e-300) + 1e-12))
  for (k in seq_along(pw)) {
    stones[[length(stones) + 1]] <- data.frame(program = "RevBayes", model = clock, rep = r, power = pw[k],
      n = length(L[[k]]), meanL = mean(L[[k]]), sdL = stats::sd(L[[k]]), essL = Ess(L[[k]]))
  }
  o <- order(pw); pw <- pw[o]; L <- L[o]
  ss <- SteppingStone(pw, lapply(L, list))
  own <- if (file.exists(paste0(stem, "_ml.txt"))) scan(paste0(stem, "_ml.txt"), what = "", quiet = TRUE) else character()
  est[[length(est) + 1]] <- data.frame(program = "RevBayes", model = clock, rep = as.character(r),
    stones = length(pw) - 1L, logZ_ss = ss$logZ, se_ss = ss$se, logZ_ps = PathSampling(pw, vapply(L, mean, 0)),
    logZ_ss_even = SteppingStone(pw[c(TRUE, FALSE)], lapply(L, list)[c(TRUE, FALSE)])$logZ,
    rb_ss = if (length(own) >= 2) as.numeric(own[2]) else NA_real_,
    rb_ps = if (length(own) >= 4) as.numeric(own[4]) else NA_real_)
  # Prior (power 0) and power-1 stone monitors
  for (k in c(cats, 1L)) {
    lg <- sprintf("%s_stone_%d.log", stem, k); tf <- sprintf("%s_stone_%d.trees", stem, k)
    if (!file.exists(lg) || !file.exists(tf)) next
    p <- ReadTab(lg); t <- ReadTab(tf)
    t <- t[grepl(";\\s*$", t$fbd_tree), ]
    it <- intersect(p$Iteration, t$Iteration)
    p <- p[match(it, p$Iteration), ]; t <- t[match(it, t$Iteration), ]
    keep <- Burn(nrow(p)); p <- p[keep, ]; f <- Facts(t$fbd_tree[keep])
    row <- data.frame(source = sprintf("RevBayes %s power %g", clock, rbPow[k]), rep = r,
      lambda = p$lambda, mu = p$mu, psi = p$psi, originAge = p$origin_time, rootAge = p$root_age,
      crownAge = f$crownAge, nSA = f$nSA, rootSA = f$rootSA)
    if (k == cats) rbPrior <- rbind(rbPrior, row) else rbStone1 <- rbind(rbStone1, row)
  }
}
stones <- do.call(rbind, stones)
utils::write.csv(stones, file.path(outDir, "ml-stones.csv"), row.names = FALSE)

# Pooled RevBayes estimate: mean over replicates, SE from their spread
est <- do.call(rbind, lapply(est, function(e) { for (k in c("rb_ss", "rb_ps")) if (is.null(e[[k]])) e[[k]] <- NA_real_; e }))
for (clock in unique(est$model[est$program == "RevBayes"])) {
  e <- est[est$program == "RevBayes" & est$model == clock, ]
  est <- rbind(est, data.frame(program = "RevBayes", model = clock, rep = "mean", stones = e$stones[1],
    logZ_ss = mean(e$logZ_ss), se_ss = stats::sd(e$logZ_ss) / sqrt(nrow(e)), logZ_ps = mean(e$logZ_ps),
    logZ_ss_even = mean(e$logZ_ss_even), rb_ss = mean(e$rb_ss), rb_ps = mean(e$rb_ps)))
}
est$sdReplicates <- NA_real_
for (k in unique(paste(est$program, est$model))) {
  i <- paste(est$program, est$model) == k & !est$rep %in% c("pooled", "mean")
  est$sdReplicates[paste(est$program, est$model) == k & est$rep %in% c("pooled", "mean")] <- stats::sd(est$logZ_ss[i])
}
utils::write.csv(est, file.path(outDir, "ml-estimates.csv"), row.names = FALSE)

# ---- Sampled ancestor at the root ---------------------------------------------------
LoadRbTrees <- function(dir) {
  stem <- file.path(dir, "output", paste0("penguins_morph_", basename(dir)))
  t <- ReadTab(paste0(stem, ".trees"))
  t <- t[grepl(";\\s*$", t$fbd_tree), ]
  t$fbd_tree[Burn(nrow(t))]
}
rsa <- list()
for (clock in c("strict", "ucln")) {
  for (d in sort(Sys.glob(file.path(rbRoot, paste0(clock, "_run*"))))) {
    f <- Facts(LoadRbTrees(d))
    rsa[[length(rsa) + 1]] <- data.frame(source = paste("RevBayes", clock, "posterior (comparator MCMC)"),
                                         rep = basename(d), n = nrow(f), pRootSA = mean(f$rootSA))
  }
}
if (!is.null(rbPrior)) {
  for (k in unique(paste(rbPrior$source, rbPrior$rep))) {
    x <- rbPrior[paste(rbPrior$source, rbPrior$rep) == k, ]
    rsa[[length(rsa) + 1]] <- data.frame(source = x$source[1], rep = as.character(x$rep[1]),
                                         n = nrow(x), pRootSA = mean(x$rootSA))
  }
}
if (!is.null(rbStone1)) {
  for (k in unique(paste(rbStone1$source, rbStone1$rep))) {
    x <- rbStone1[paste(rbStone1$source, rbStone1$rep) == k, ]
    rsa[[length(rsa) + 1]] <- data.frame(source = x$source[1], rep = as.character(x$rep[1]),
                                         n = nrow(x), pRootSA = mean(x$rootSA))
  }
}
rsa <- do.call(rbind, rsa)
utils::write.csv(rsa, file.path(outDir, "root-sa.csv"), row.names = FALSE)
Pool <- function(pattern) {
  x <- rsa[grepl(pattern, rsa$source), ]
  PoolProportion(x$pRootSA, x$n)
}
pi0 <- Pool("power 0$")
post <- list(strict = Pool("strict posterior"), ucln = Pool("ucln posterior"))
Corr <- function(clock) log1p(-post[[clock]][["p"]]) - log1p(-pi0[["p"]])
CorrSe <- function(clock) sqrt((post[[clock]][["se"]] / (1 - post[[clock]][["p"]]))^2 +
                               (pi0[["se"]] / (1 - pi0[["p"]]))^2)

# ---- Prior check: CET t = 0 against RevBayes power 0 -------------------------------
if (!is.null(cetPrior) && !is.null(rbPrior)) {
  qs <- c("lambda", "mu", "psi", "originAge", "rootAge", "crownAge", "nSA")
  sets <- list()
  for (s in unique(cetPrior$source)) sets[[s]] <- cetPrior[cetPrior$source == s, ]
  sets[["RevBayes power 0, all"]] <- rbPrior
  sets[["RevBayes power 0, no root SA"]] <- rbPrior[!rbPrior$rootSA, ]
  ref <- sets[["RevBayes power 0, no root SA"]]
  pc <- do.call(rbind, lapply(names(sets), function(s) do.call(rbind, lapply(qs, function(q) {
    x <- sets[[s]][[q]]
    data.frame(sample = s, quantity = q, n = length(x), mean = mean(x),
               q05 = stats::quantile(x, 0.05, names = FALSE), median = stats::median(x),
               q95 = stats::quantile(x, 0.95, names = FALSE),
               ksD_vs_RBnoRootSA = suppressWarnings(stats::ks.test(x, ref[[q]])$statistic[[1]]))
  }))))
  utils::write.csv(pc, file.path(outDir, "prior-check.csv"), row.names = FALSE)
}

# ---- Anchor and common footing -------------------------------------------------------
Get <- function(program, model) {
  e <- est[est$program == program & est$model == model & est$rep %in% c("pooled", "mean"), ]
  if (!nrow(e)) return(c(NA, NA))
  perRep <- est$logZ_ss[est$program == program & est$model == model & !est$rep %in% c("pooled", "mean")]
  se <- if (program == "RevBayes") e$se_ss else CetPooledSe(e$se_ss, perRep)
  c(e$logZ_ss, se)
}
bf <- if (file.exists(bfFile)) utils::read.csv(bfFile) else NULL
Edge <- function(from, to) {
  if (is.null(bf)) return(c(NA, NA))
  e <- bf[bf$from == from & bf$to == to, ]
  if (!nrow(e)) return(c(NA, NA))
  c(e$logBF, max(e$se, e$seReplicates))
}
c00 <- Get("CET", "b0g0"); c10 <- Get("CET", "b1g0")
rs <- Get("RevBayes", "strict"); ru <- Get("RevBayes", "ucln")
e01 <- Edge("b0g0", "b0g1"); e11 <- Edge("b1g0", "b1g1")
pred <- Corr("strict")
obs <- c00[1] - rs[1]
anchor <- c(
  sprintf("P(root SA) under the RevBayes prior (power-0 stones, %d samples): %.4f (SE %.4f)",
          sum(rsa$n[grepl("power 0$", rsa$source)]), pi0[["p"]], pi0[["se"]]),
  sprintf("P(root SA | characters), RevBayes strict / UCLN posterior: %.4f / %.4f",
          post$strict[["p"]], post$ucln[["p"]]),
  sprintf("log Z, CET (0, 0):      %.3f (SE %.3f)", c00[1], c00[2]),
  sprintf("log Z, RevBayes strict: %.3f (SE %.3f)", rs[1], rs[2]),
  sprintf("observed  log Z_CET(0,0) - log Z_RB,strict = %.3f (SE %.3f)", obs, sqrt(c00[2]^2 + rs[2]^2)),
  sprintf("predicted log(1 - P(root SA | ch)) - log(1 - P(root SA)) = %.3f (SE %.3f)", pred, CorrSe("strict")),
  sprintf("gap (observed - predicted) = %.3f (SE %.3f)", obs - pred,
          sqrt(c00[2]^2 + rs[2]^2 + CorrSe("strict")^2)))
writeLines(anchor, file.path(outDir, "anchor.txt"))
cat(anchor, sep = "\n")
cf <- data.frame(
  model = c("CET (0, 0)", "CET (0, 1)", "CET (1, 0)", "CET (1, 1)",
            "RevBayes strict, on CET's support", "RevBayes UCLN, on CET's support",
            "RevBayes strict", "RevBayes UCLN"),
  how = c("ladder", "ladder (0, 0) + bridge (0, 0) -> (0, 1)", "ladder", "ladder (1, 0) + bridge (1, 0) -> (1, 1)",
          "stepping-stone + root-SA term", "stepping-stone + root-SA term", "stepping-stone", "stepping-stone"),
  logZ = c(c00[1], c00[1] + e01[1], c10[1], c10[1] + e11[1], rs[1] + Corr("strict"), ru[1] + Corr("ucln"), rs[1], ru[1]),
  se = c(c00[2], sqrt(c00[2]^2 + e01[2]^2), c10[2], sqrt(c10[2]^2 + e11[2]^2),
         sqrt(rs[2]^2 + CorrSe("strict")^2), sqrt(ru[2]^2 + CorrSe("ucln")^2), rs[2], ru[2]))
ref <- cf$logZ[cf$model == "RevBayes UCLN, on CET's support"]
cf$vsUclnOnCetSupport <- cf$logZ - ref
utils::write.csv(cf, file.path(outDir, "common-footing.csv"), row.names = FALSE)
print(cf, digits = 5, row.names = FALSE)
