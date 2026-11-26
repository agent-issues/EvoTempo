# Penguin CET pipeline: sensitivity of the gamma Bayes factor to the prior on s
# (README, "Model and priors", open decision; #107 review, R107-16).
#
#   Rscript data-raw/empirical/penguins/cet/s-prior-reweight.R <results dir> [out.csv]
#
# Reads <results dir>/samples.csv.gz (stage1-results.R: post-burn-in draws per
# chain of the four production cells).  Changing the prior on s from pi to pi'
# multiplies a cell's evidence by E_post[pi'(s) / pi(s)] (importance
# reweighting of that cell's posterior draws), so the gamma Bayes factor at
# fixed beta moves by
#   log E_post(beta, 1)[w] - log E_post(beta, 0)[w],   w = pi'(s) / pi(s).
# pi: Gamma(0.5396, scale 0.3819) on (0, 3] (CetPriors()).  pi': log-normal
# alternatives, also on (0, 3]: matched to pi's median and 97.5% quantile, and
# LN(log 0.2, 1), LN(log 0.2, 0.5), LN(log 0.5, 1).  Reports log E[w] per cell
# (and its range over chains) and the Kish ESS of the weights.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
s <- utils::read.csv(file.path(args[[1]], "samples.csv.gz"))
s <- s[s$analysis %in% c("b0g0", "b0g1", "b1g0", "b1g1"), c("analysis", "chain", "s")]
shape <- 0.5396; scale <- 0.3819; up <- 3
LogPi <- function(x) stats::dgamma(x, shape, scale = scale, log = TRUE) -
  stats::pgamma(up, shape, scale = scale, log.p = TRUE)
med <- stats::qgamma(0.5, shape, scale = scale)
q975 <- stats::qgamma(0.975, shape, scale = scale)
alts <- list(matched = c(log(med), (log(q975) - log(med)) / stats::qnorm(0.975)),
             ln_0.2_1 = c(log(0.2), 1), ln_0.2_0.5 = c(log(0.2), 0.5), ln_0.5_1 = c(log(0.5), 1))
cat(sprintf("Gamma prior: median %.4f, 97.5%% quantile %.4f, P(s < 0.01) %.3f\n", med, q975,
            stats::pgamma(0.01, shape, scale = scale)))
LogMeanExp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
out <- NULL
for (a in names(alts)) {
  p <- alts[[a]]
  lw <- stats::dlnorm(s$s, p[1], p[2], log = TRUE) - stats::plnorm(up, p[1], p[2], log.p = TRUE) - LogPi(s$s)
  for (cell in sort(unique(s$analysis))) {
    i <- s$analysis == cell
    w <- exp(lw[i] - max(lw[i]))
    perChain <- tapply(lw[i], s$chain[i], LogMeanExp)
    out <- rbind(out, data.frame(prior = a, meanlog = p[1], sdlog = p[2], cell = cell,
                                 logEw = LogMeanExp(lw[i]), kish = sum(w)^2 / sum(w^2), n = sum(i),
                                 chainMin = min(perChain), chainMax = max(perChain)))
  }
}
print(out, digits = 4, row.names = FALSE)
shift <- do.call(rbind, lapply(names(alts), function(a) {
  e <- stats::setNames(out$logEw[out$prior == a], out$cell[out$prior == a])
  data.frame(prior = a, gammaShiftBeta0 = e[["b0g1"]] - e[["b0g0"]], gammaShiftBeta1 = e[["b1g1"]] - e[["b1g0"]],
             betaShiftGamma1 = e[["b1g1"]] - e[["b0g1"]], betaShiftGamma0 = e[["b1g0"]] - e[["b0g0"]])
}))
cat("\nShift in log BF under each alternative prior on s:\n")
print(shift, digits = 3, row.names = FALSE)
if (length(args) >= 2L) {
  utils::write.csv(merge(out, shift, by = "prior"), args[[2]], row.names = FALSE)
}
