# Penguin CET pipeline, stage 1 figures, from the outputs of stage1-results.R.
#
#   Rscript data-raw/empirical/penguins/cet/stage1-figures.R <results dir> [<fig dir>]
#
# Reads samples.csv.gz (post-burn-in draws per chain) and writes PNGs:
#   fig-posteriors.png   pooled posterior densities of the main quantities,
#                        four CET cells and the two RevBayes comparators
#   fig-b0g0-strict.png  cell (0, 0) against the strict-clock comparator, the
#                        RevBayes counterpart of a constant-rate, constant-clock
#                        model, per chain (the lambda < mu check)
#   fig-traces.png       per-chain traces (post burn-in) of log10 theta, s,
#                        crown age and tree length in each CET cell (the
#                        quantities flagged at 200 000 iterations; charLL,
#                        flagged at 400 000, is in convergence.csv)
# Base graphics only.
args <- commandArgs(trailingOnly = TRUE)
resDir <- args[[1]]
figDir <- if (length(args) >= 2L) args[[2]] else resDir
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
s <- utils::read.csv(file.path(resDir, "samples.csv.gz"))
analyses <- c("b0g0", "b0g1", "b1g0", "b1g1", "revbayes_strict", "revbayes_ucln")
analyses <- intersect(analyses, unique(s$analysis))
col <- c(b0g0 = "#1b9e77", b0g1 = "#d95f02", b1g0 = "#7570b3", b1g1 = "#e7298a",
         revbayes_strict = "black", revbayes_ucln = "grey50")
lty <- c(b0g0 = 1, b0g1 = 1, b1g0 = 1, b1g1 = 1, revbayes_strict = 2, revbayes_ucln = 2)
lab <- c(crownAge = "crown age (Ma)", rootAge = "root age (Ma)", originAge = "origin age (Ma)",
         nSA = "sampled ancestors", lambda = "lambda (per Myr)", mu = "mu (per Myr)",
         diversification = "lambda - mu", turnover = "mu / lambda", psi = "psi (per Myr)",
         rate = "clock rate (per Myr)", theta = "theta", s = "s", charLL = "character log-likelihood")
Dens <- function(p, list = analyses, logx = FALSE) {
  ds <- lapply(list, function(a) {
    x <- s[s$analysis == a, p]
    x <- x[is.finite(x)]
    if (!length(x)) return(NULL)
    if (logx) x <- log10(x)
    if (p == "nSA") {
      t <- table(factor(x, levels = 0:55)) / length(x)
      return(list(x = 0:55, y = as.numeric(t)))
    }
    stats::density(x, n = 512)
  })
  names(ds) <- list
  ds <- ds[!vapply(ds, is.null, TRUE)]
  xr <- range(unlist(lapply(ds, `[[`, "x")))
  if (p == "nSA") xr <- c(0, 30)
  yr <- c(0, max(unlist(lapply(ds, `[[`, "y"))))
  plot(NA, xlim = xr, ylim = yr, xlab = if (logx) paste("log10", lab[[p]]) else lab[[p]],
       ylab = if (p == "nSA") "probability" else "density", main = "")
  for (a in names(ds)) {
    if (p == "nSA") {
      graphics::lines(ds[[a]]$x, ds[[a]]$y, col = col[[a]], lty = lty[[a]], type = "s")
    } else {
      graphics::lines(ds[[a]], col = col[[a]], lty = lty[[a]], lwd = 1.5)
    }
  }
}
Legend <- function(list = analyses) {
  graphics::legend("topright", legend = sub("revbayes_", "RevBayes ", list), col = col[list],
                   lty = lty[list], lwd = 1.5, bty = "n", cex = 0.8)
}

grDevices::png(file.path(figDir, "fig-posteriors.png"), width = 2400, height = 1800, res = 200)
graphics::par(mfrow = c(3, 4), mar = c(4, 4, 1, 1))
for (p in c("crownAge", "rootAge", "originAge", "nSA", "lambda", "mu", "turnover", "psi", "rate",
            "theta", "s", "charLL")) {
  Dens(p, logx = p %in% c("rate", "theta"))
  if (p == "crownAge") Legend()
}
grDevices::dev.off()

grDevices::png(file.path(figDir, "fig-b0g0-strict.png"), width = 2400, height = 1200, res = 200)
graphics::par(mfrow = c(2, 4), mar = c(4, 4, 1, 1))
for (p in c("lambda", "mu", "diversification", "turnover", "psi", "rate", "crownAge", "rootAge")) {
  xr <- range(stats::quantile(s[s$analysis %in% c("b0g0", "revbayes_strict"), p], c(0.001, 0.999)))
  plot(NA, xlim = xr, ylim = c(0, 1), xlab = lab[[p]], ylab = "cumulative probability")
  for (a in c("b0g0", "revbayes_strict")) {
    for (k in unique(s$chain[s$analysis == a])) {
      x <- s[s$analysis == a & s$chain == k, p]
      graphics::lines(stats::ecdf(x), col = grDevices::adjustcolor(col[[a]], 0.6), lty = lty[[a]],
                      do.points = FALSE, verticals = TRUE)
    }
  }
  if (p == "diversification") graphics::abline(v = 0, col = "grey70")
  if (p == "lambda") {
    graphics::legend("bottomright", legend = c("CET (0, 0), per chain", "RevBayes strict, per chain"),
                     col = col[c("b0g0", "revbayes_strict")], lty = c(1, 2), bty = "n", cex = 0.8)
  }
}
grDevices::dev.off()

cells <- intersect(c("b0g0", "b0g1", "b1g0", "b1g1"), analyses)
grDevices::png(file.path(figDir, "fig-traces.png"), width = 2400, height = 600 * length(cells), res = 200)
graphics::par(mfrow = c(length(cells), 4), mar = c(3, 4, 2, 1))
chainCol <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
for (a in cells) {
  for (p in c("theta", "s", "crownAge", "treeLength")) {
    x <- s[s$analysis == a, ]
    ks <- unique(x$chain)
    yr <- range(if (p == "theta") log10(x[[p]]) else x[[p]], na.rm = TRUE)
    plot(NA, xlim = range(x$draw), ylim = yr, xlab = "", ylab = if (p == "theta") "log10 theta" else p,
         main = paste(a, p), cex.main = 0.9)
    for (i in seq_along(ks)) {
      y <- x[x$chain == ks[i], ]
      graphics::lines(y$draw, if (p == "theta") log10(y[[p]]) else y[[p]],
                      col = grDevices::adjustcolor(chainCol[i], 0.6), lwd = 0.5)
    }
  }
}
grDevices::dev.off()
cat("figures written to", figDir, "\n")
