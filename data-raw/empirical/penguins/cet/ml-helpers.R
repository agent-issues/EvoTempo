# Helpers of ml-summary.R that are tested from the source tree
# (tests/testthat/test-PenguinMlHelpers.R; review of #111).

# Standard error of a CET ladder's pooled log Z (R111-1).  The replicate is the
# unit: the spread of the per-replicate stepping-stone estimates already
# contains their within-chain noise, so the SE is sd / sqrt(n), floored at the
# pooled estimator's own SE (as Edge() and bayes-factors.R take
# max(se, seReplicates)).  With one replicate the spread is unavailable and the
# estimator SE is returned (R111-9).
CetPooledSe <- function(sePooled, perRep) {
  perRep <- perRep[is.finite(perRep)]
  n <- length(perRep)
  if (n < 2L) return(sePooled)
  max(sePooled, stats::sd(perRep) / sqrt(n))
}

# A proportion pooled over runs, weighted by their sample counts n, with the SE
# of that weighted mean from the between-run spread (R111-8).  With equal
# weights the SE is sd(p) / sqrt(m).
PoolProportion <- function(p, n) {
  w <- n / sum(n)
  m <- length(p)
  pBar <- sum(w * p)
  se <- if (m < 2L) NA_real_ else sqrt(m / (m - 1) * sum(w^2 * (p - pBar)^2))
  c(p = pBar, se = se)
}

# Problems with a cell's CET ladder (R111-3).  `chains`: one row per chain
# directory with columns cell, power, rep, exists (FALSE for a dangling link)
# and complete (progress.txt ends with "complete").  `powers`: the expected
# powers t < 1.  Returns a character vector, empty if the ladder is whole:
# every replicate has every power, the t = 1 rung (production chains) is
# present, and every chain exists and is complete.
CetLadderProblems <- function(chains, powers) {
  out <- character()
  bad <- chains[!chains$exists, ]
  if (nrow(bad)) out <- c(out, sprintf("%s rep %d power %g: dangling link or missing directory",
                                       bad$cell, bad$rep, bad$power))
  inc <- chains[chains$exists & !chains$complete, ]
  if (nrow(inc)) out <- c(out, sprintf("%s rep %d power %g: chain not complete",
                                       inc$cell, inc$rep, inc$power))
  for (cell in unique(chains$cell)) {
    sel <- chains[chains$cell == cell & chains$exists, ]
    if (!any(sel$power == 1)) {
      out <- c(out, sprintf("%s: no t = 1 rung (production chains; ml-root.sh)", cell))
    }
    for (r in sort(unique(sel$rep[sel$power != 1]))) {
      have <- sel$power[sel$rep == r & sel$power != 1]
      missing <- setdiff(signif(powers, 12), signif(have, 12))
      if (length(missing)) out <- c(out, sprintf("%s rep %d: %d of %d powers missing", cell, r,
                                                 length(missing), length(powers)))
    }
  }
  out
}

# Every second power of the schedule plus t = 1, by position in the schedule
# rather than in the powers present, so a missing power cannot drop t = 1
# from the even sub-ladder (R111-3).
EvenPowers <- function(pw, powers) {
  j <- match(signif(pw, 12), signif(powers, 12))
  pw == 1 | (!is.na(j) & j %% 2L == 1L)
}
