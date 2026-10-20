# Rung 7 (issue #7): summaries of the power study run by data-raw/rung7-power.R.
#   Rscript data-raw/rung7-power-summary.R <outDir>
# reads <outDir>/replicates/tree*-n*.rds and writes the CSV files of
# data-raw/rung7-results/ (see its README).  Base R only.
args <- commandArgs(trailingOnly = TRUE)
outDir <- if (length(args)) args[[1]] else "data-raw/rung7-results"
resDir <- if (length(args) > 1) args[[2]] else "data-raw/rung7-results"
dir.create(resDir, showWarnings = FALSE, recursive = TRUE)
files <- list.files(file.path(outDir, "replicates"), "^tree.*\\.rds$", full.names = TRUE)
reps <- lapply(files, readRDS)
cellLab <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)")
Softmax <- function(x) {
  e <- exp(x - max(x))
  e / sum(e)
}
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

runs <- do.call(rbind, lapply(reps, function(r) {
  data.frame(
    tree = r$g, block = r$block, cell = r$cell, truth = cellLab[r$cell],
    psiLevel = r$psiLevel, size = r$size, nChar = r$nChar,
    nExtant = r$nExtant, nFossil = r$nFossil, nAncestor = r$nAncestor,
    tries = r$tries,
    lambda = r$truth$lambda, mu = r$truth$mu, psi = r$truth$psi,
    theta = r$truth$theta, s = r$truth$s, rate = r$truth$rate,
    logZ00 = r$logZ[[1]], logZ01 = r$logZ[[2]], logZ10 = r$logZ[[3]],
    logZ11 = r$logZ[[4]], se01 = r$se[[2]], se10 = r$se[[3]], se11 = r$se[[4]],
    logZ11TI = r$logZTI[[4]],
    logZT0 = r$logZT0, logZT1 = r$zT1$logZ, logZT1se = r$zT1$se,
    logRStrict = r$logRStrict, logR11 = r$ladder$logR, logR11se = r$ladder$se,
    bridgeLogBF = r$bridge$logBF, bridgeSE = r$bridge$se,
    bridgeEss0 = r$bridge$essF0, bridgeEss1 = r$bridge$essF1,
    bridgeFrom0 = r$bridge$logBFfrom0, bridgeFrom1 = r$bridge$logBFfrom1,
    bridgeCheck = max(r$bridge$check0, r$bridge$check1),
    essTop11 = r$ladder$top$essLogL, essTop01 = r$top01$essLogL,
    minEssTop = min(r$ladder$top$posteriorEss, r$top01$posteriorEss),
    attempts = sum(r$ladder$attempts) + r$top01$attempts,
    logZ01ladder = if (is.null(r$check01)) NA else r$check01$logZ,
    logZ01ladderSE = if (is.null(r$check01)) NA else r$check01$se,
    hours = (r$secs[["level"]] + r$secs[["tree"]]) / 3600,
    sha = r$sha, stringsAsFactors = FALSE)
}))
runs <- runs[order(runs$tree, runs$nChar), ]

# P(size in range | cell, psi, size) from the parameter draws per tree
# (1 / mean tries; one value per tree, so deduplicated over nChar levels)
trees <- runs[!duplicated(runs$tree), ]
acc <- aggregate(tries ~ cell + psiLevel + size, trees,
                 function(x) c(n = length(x), sum = sum(x)))
acc <- data.frame(acc[c("cell", "psiLevel", "size")],
                  nTree = acc$tries[, "n"], nDraw = acc$tries[, "sum"])
acc$pAccept <- acc$nTree / acc$nDraw

logZ <- as.matrix(runs[c("logZ00", "logZ01", "logZ10", "logZ11")])
post <- t(apply(logZ, 1, Softmax))
# Correction for rejection on tree size within a fixed true cell (README):
# posterior in the simulated ensemble is proportional to p(D | m) / P(A | m)
logA <- t(vapply(seq_len(nrow(runs)), function(i) {
  vapply(1:4, function(m) {
    a <- acc$pAccept[acc$cell == m & acc$psiLevel == runs$psiLevel[i] &
                       acc$size == runs$size[i]]
    if (length(a)) log(a) else 0
  }, 0)
}, numeric(4)))
postC <- t(apply(logZ - logA, 1, Softmax))
colnames(post) <- colnames(postC) <- cellLab
runs$selected <- cellLab[max.col(post, "first")]
runs$correct <- runs$selected == runs$truth
runs$pTrue <- post[cbind(seq_len(nrow(runs)), runs$cell)]
runs$pTrueCorrected <- postC[cbind(seq_len(nrow(runs)), runs$cell)]
for (m in 1:4) {
  runs[[paste0("p", gsub("[(,)]", "", cellLab[m]))]] <- post[, m]
  runs[[paste0("pc", gsub("[(,)]", "", cellLab[m]))]] <- postC[, m]
}
runs$logBFbeta1 <- runs$logZ11 - runs$logZ01
runs$logBFbeta0 <- runs$logZ10 - runs$logZ00
runs$logBFgamma0 <- runs$logZ01 - runs$logZ00
runs$logBFgamma1 <- runs$logZ11 - runs$logZ10
runs$seBeta1 <- runs$bridgeSE
runs$seGamma1 <- runs$logR11se
runs$seGamma0 <- sqrt(runs$se11^2 + runs$bridgeSE^2)
write.csv(runs, file.path(resDir, "rung7-runs.csv"), row.names = FALSE)
write.csv(acc, file.path(resDir, "rung7-acceptance.csv"), row.names = FALSE)

# Recovery of the true cell, per factor combination and true cell
Wilson <- function(k, n, z = 1.96) {
  p <- k / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(lo = centre - half, hi = centre + half)
}
key <- c("nChar", "psiLevel", "size", "truth")
groups <- split(runs, runs[key], drop = TRUE)
recovery <- do.call(rbind, lapply(groups, function(d) {
  ci <- Wilson(sum(d$correct), nrow(d))
  data.frame(d[1, key], n = nrow(d), pCorrect = mean(d$correct),
             lo = ci[["lo"]], hi = ci[["hi"]], meanPTrue = mean(d$pTrue),
             meanPTrueCorrected = mean(d$pTrueCorrected))
}))
recovery <- recovery[do.call(order, recovery[key]), ]
write.csv(recovery, file.path(resDir, "rung7-recovery.csv"), row.names = FALSE)

# Pairwise contrasts: decision by the sign of log BF, and "positive" evidence
# (|log BF| > 1, Kass & Raftery 1995: 2 ln BF > 2) for the true model
Contrast <- function(name, bf, se, cellA, cellB) {
  # cellA: the cell favoured by log BF > 0; cellB the alternative
  d <- runs[runs$cell %in% c(cellA, cellB), ]
  d$bf <- d[[bf]]
  d$se <- d[[se]]
  d$trueA <- d$cell == cellA
  sp <- split(d, d[c("nChar", "psiLevel", "size", "trueA")], drop = TRUE)
  do.call(rbind, lapply(sp, function(x) {
    right <- if (x$trueA[1]) x$bf > 0 else x$bf < 0
    strong <- if (x$trueA[1]) x$bf > 1 else x$bf < -1
    ci <- Wilson(sum(right), nrow(x))
    data.frame(contrast = name, nChar = x$nChar[1], psiLevel = x$psiLevel[1],
               size = x$size[1], truth = cellLab[if (x$trueA[1]) cellA else cellB],
               n = nrow(x), pRight = mean(right), lo = ci[["lo"]], hi = ci[["hi"]],
               pPositive = mean(strong), medianLogBF = stats::median(x$bf),
               fracUnresolved = mean(abs(x$bf) < 2 * x$se))
  }))
}
contrasts <- rbind(
  Contrast("beta | gamma = 1: (1,1) vs (0,1)", "logBFbeta1", "seBeta1", 4, 2),
  Contrast("beta | gamma = 0: (1,0) vs (0,0)", "logBFbeta0", "se10", 3, 1),
  Contrast("gamma | beta = 1: (1,1) vs (1,0)", "logBFgamma1", "seGamma1", 4, 3),
  Contrast("gamma | beta = 0: (0,1) vs (0,0)", "logBFgamma0", "seGamma0", 2, 1))
write.csv(contrasts, file.path(resDir, "rung7-contrasts.csv"), row.names = FALSE)

# Calibration: pooled (dataset, cell) pairs binned by posterior probability
Calibration <- function(P, label) {
  p <- as.vector(P)
  hit <- as.vector(outer(runs$cell, 1:4, `==`))
  bins <- cut(p, c(0, 0.01, seq(0.1, 0.9, 0.1), 0.99, 1), include.lowest = TRUE)
  do.call(rbind, lapply(split(data.frame(p, hit), bins, drop = TRUE), function(x) {
    ci <- Wilson(sum(x$hit), nrow(x))
    data.frame(probabilities = label, n = nrow(x), meanP = mean(x$p),
               observed = mean(x$hit), lo = ci[["lo"]], hi = ci[["hi"]])
  }))
}
calibration <- rbind(Calibration(post, "uniform prior"),
                     Calibration(postC, "size-corrected"))
calibration$bin <- rownames(calibration)
write.csv(calibration, file.path(resDir, "rung7-calibration.csv"), row.names = FALSE)

# Data-averaged posterior (Schad et al. 2022): E[P(m | D)] = P(m) = 1/4, with
# true cells balanced; reported per nChar level, with Monte Carlo SE
averaged <- do.call(rbind, lapply(split(seq_len(nrow(runs)), runs$nChar), function(i) {
  data.frame(nChar = runs$nChar[i[1]], cell = cellLab, n = length(i),
             meanP = colMeans(post[i, , drop = FALSE]),
             meanPCorrected = colMeans(postC[i, , drop = FALSE]),
             se = apply(postC[i, , drop = FALSE], 2, stats::sd) / sqrt(length(i)),
             fracTrue = vapply(1:4, function(m) mean(runs$cell[i] == m), 0))
}))
write.csv(averaged, file.path(resDir, "rung7-averaged.csv"), row.names = FALSE)

# Bridge vs a separate (0,1) ladder (replicate block 1)
chk <- runs[!is.na(runs$logZ01ladder), ]
chk$z <- (chk$logZ01 - chk$logZ01ladder) / sqrt(chk$se01^2 + chk$logZ01ladderSE^2)
write.csv(chk[c("tree", "nChar", "psiLevel", "size", "truth", "logZ01", "se01",
                "logZ01ladder", "logZ01ladderSE", "z")],
          file.path(resDir, "rung7-bridge-check.csv"), row.names = FALSE)

cat(sprintf("%d runs from %d trees; %.1f core-hours of replicate run time\n",
            nrow(runs), nrow(trees), sum(runs$hours)))
print(recovery[c(key, "n", "pCorrect")], row.names = FALSE)
print(contrasts[grepl("^beta \\| gamma = 1", contrasts$contrast),
                c("nChar", "psiLevel", "size", "truth", "n", "pRight", "pPositive")],
      row.names = FALSE)
