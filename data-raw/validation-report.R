# Builds data-raw/validation-report/: rung 1-2 report (markdown) and figures,
# from data-raw/bm-results/ and data-raw/rung2-results/.
library("EvoTempo")
repDir <- "data-raw/validation-report"
dir.create(repDir, showWarnings = FALSE)
bm <- readRDS("data-raw/bm-results/bm-columns.rds")
cv <- read.csv("data-raw/bm-results/convergence.csv")
r2 <- readRDS("data-raw/rung2-results/rung2.rds")
pnFiles <- c(code = "data-raw/bm-results/pn-budd-mann.rds",
             converged = "data-raw/bm-results/pn-central-101.rds")
pn <- lapply(pnFiles, function(f) if (file.exists(f)) readRDS(f)[, 1] else NULL)
lambda <- 0.51
mu <- 0.5
theta <- 0.01
Fmt <- function(v, d = 4) paste(formatC(v, digits = d, format = "fg", flag = "#"), collapse = " / ")
Pct <- function(a, b) sprintf("%+.1f%%", 100 * (a / b - 1))
Png <- function(name, width = 17, height = 8.5) {
  png(file.path(repDir, name), width = width, height = height, units = "cm", res = 200,
      pointsize = 9)
  par(mar = c(4.4, 4.8, 1.8, 0.8), mgp = c(3.3, 0.7, 0), las = 1, cex.main = 0.95)
}

# Clade-size summaries from P_n (start x = 0, T = 500)
PnSummary <- function(p) {
  if (is.null(p)) return(rep(NA, 5))
  s <- 1 - p[1]
  q <- p[-1] / s # P(n | n > 0), n = 1..nMax
  n <- seq_along(q)
  cdf <- cumsum(q)
  tail <- 1 - cdf[length(cdf)]
  # Sister groups: P(smaller share <= 0.05) = 2 sum_n q_n P(other >= 19 n)
  Above <- function(k) ifelse(k <= length(q), 1 - c(0, cdf)[k], NA)
  small <- n[19 * n <= length(q)]
  share <- 2 * sum(q[small] * Above(19 * small))
  rest <- 1 - cdf[max(small)] # mass with smaller sister > max(small)
  c(lt1000 = cdf[999], gt50000 = tail, share = share,
    shareUpper = share + 2 * rest * tail, mass = s)
}
pnSum <- lapply(pn, PnSummary)

# --- Figure 1: convergence of Fig 1A quantities ---
central <- cv[cv$scheme == "central" & cv$zMax == 10, ]
ref <- with(central[order(central$h), ], c(
  meanClade = meanClade[1] + (meanClade[1] - meanClade[2]) / 3,
  experienced = experienced[1] + (experienced[1] - experienced[2]) / 3))
Png("fig-convergence.png")
par(mfrow = c(1, 2))
schemes <- list(central = c("central", 10, "black", 16), upwind = c("upwind", 10, "#D55E00", 17),
                sqra = c("sqra", 10, "#0072B2", 15), "Budd & Mann" = c("budd-mann", 10, "#009E73", 18))
for (q in c("meanClade", "experienced")) {
  plot(NA, xlim = c(0.02, 0.45), ylim = c(1e-5, 1), log = "xy",
       xlab = "Grid spacing h (units of s)", ylab = "Relative error",
       main = c(meanClade = "Mean surviving clade size", experienced = "Experienced clade size")[q])
  for (nm in names(schemes)) {
    sc <- schemes[[nm]]
    d <- cv[cv$scheme == sc[1] & cv$zMax == as.numeric(sc[2]), ]
    d <- d[order(d$h), ]
    lines(d$h, abs(d[[q]] / ref[q] - 1), type = "b", col = sc[3], pch = as.numeric(sc[4]))
  }
  lines(c(0.05, 0.4), 0.02 * c(1, 8)^1, lty = 3, col = "grey40")
  lines(c(0.05, 0.4), 0.0006 * c(1, 64), lty = 2, col = "grey40")
  if (q == "meanClade") {
    legend("bottomright", c(names(schemes), "order 1", "order 2"), bty = "n",
           col = c(vapply(schemes, `[`, "", 3), "grey40", "grey40"),
           pch = c(vapply(schemes, function(s) as.numeric(s[4]), 0), NA, NA),
           lty = c(1, 1, 1, 1, 3, 2))
  }
}
dev.off()

# --- Figure 2: Fig 1A clade-size distribution ---
if (!is.null(pn$code) && !is.null(pn$converged)) {
  Png("fig-clade-size.png")
  par(mfrow = c(1, 2))
  n <- 1:50000
  plot(n, pn$converged[-1] / (1 - pn$converged[1]), type = "l", log = "y",
       xlab = "Species, n", ylab = "P(n | n > 0)", ylim = c(1e-7, 1e-3), xlim = c(0, 50000),
       main = "Clade size after 500 Myr, x0 = 0")
  lines(n, pn$code[-1] / (1 - pn$code[1]), col = "#009E73")
  lines(n, dgeom(n - 1, 1 / bm$converged$meanClade), lty = 2)
  abline(v = c(bm$converged$meanClade, bm$code$meanClade), col = c("black", "#009E73"), lty = 3)
  legend("topright", c("dCET central K = 101, z in +-5", "Budd & Mann scheme", "geometric, same mean"),
         col = c("black", "#009E73", "black"), lty = c(1, 1, 2), bty = "n")
  plot(n, cumsum(pn$converged[-1]) / (1 - pn$converged[1]), type = "l", log = "x",
       xlab = "Species, n", ylab = "P(N <= n | N > 0)", ylim = c(0, 1),
       main = "Cumulative")
  lines(n, cumsum(pn$code[-1]) / (1 - pn$code[1]), col = "#009E73")
  abline(v = c(1000, 50000), h = c(1 / 3, 3 / 4), lty = 3, col = "grey50")
  dev.off()
}

# --- Figure 3: B&M Figs 2A, 3A, 4A curves, B&M operator vs converged ---
bmGen <- BuddMannGenerator()
fine <- TempoGrid(401, c(-10, 10))
fineQ <- TempoGenerator(fine, theta)
t <- seq(0, 500, by = 5)
Curves <- function(x, Q) {
  r <- BinRates(x, lambda, mu)
  idx <- vapply(c(-2, 0, 2), function(x0) which.min(abs(x - x0)), 1L)
  n <- DiversityMoments(t, Q, r$lambda, r$mu)$mean[, idx]
  m <- vapply(idx, function(k) {
    as.vector(TempoDistribution(t, Q, theta * exp(x), 0, start = k) %*% x)
  }, numeric(length(t)))
  top <- TempoDistribution(500, Q, theta * exp(x), 0, start = length(x))
  anc <- AncestralTempo(t, x, Q, r$lambda, r$mu, as.vector(top))
  list(n = n, m = m, anc = anc)
}
cBM <- Curves(bmGen$x, bmGen$Q)
cConv <- Curves(fine$x, fineQ)
Png("fig-bm-curves.png", height = 6.5)
par(mfrow = c(1, 3))
cols <- c("#0072B2", "black", "#D55E00")
matplot(t, cConv$n, type = "l", lty = 1, col = cols, log = "y", xlab = "Time since origin (Myr)",
        ylab = "Expected species", main = "Fig 2A: mean diversity")
matlines(t, cBM$n, lty = 2, col = cols)
legend("topleft", c("x0 = -2", "x0 = 0", "x0 = 2", "converged", "B&M scheme"),
       col = c(cols, "grey30", "grey30"), lty = c(1, 1, 1, 1, 2), bty = "n")
matplot(t, cConv$m, type = "l", lty = 1, col = cols, xlab = "Time since origin (Myr)",
        ylab = "Mean log-tempo", main = "Fig 3A: r = theta", ylim = c(-2.2, 2.2))
matlines(t, cBM$m, lty = 2, col = cols)
plot(t, cConv$anc[, "mean"], type = "l", xlab = "Look-back time (Myr)",
     ylab = "Mean ancestral log-tempo", main = "Fig 4A: ancestors of modern species",
     ylim = c(0, 0.7))
lines(t, cBM$anc[, "mean"], lty = 2)
abline(h = 0.6, col = "grey60", lty = 3)
text(480, 0.62, "paper: c. 0.6", adj = 1, col = "grey40")
dev.off()

# --- Figure 4: branch w law (rung 2) ---
Png("fig-branch-w.png", height = 7)
par(mfrow = c(1, 3))
qq <- ppoints(500)
for (i in 1:3) {
  w <- r2$branches[[i]][, "w"]
  plot(qexp(qq, lambda + mu), quantile(w, qq), pch = 16, cex = 0.4,
       xlab = "Exp(lambda + mu) quantile", ylab = "Simulated branch w",
       main = sprintf("x0 = %g: KS p = %.2f; control p = %.1g", c(-2, 0, 2)[i],
                      r2$branchTable$ksP[i], r2$branchTable$ksPcontrol[i]))
  points(qexp(qq, lambda + mu), quantile(r2$control[[i]][, "w"], qq), pch = 1, cex = 0.4,
         col = "#D55E00")
  points(qexp(qq, lambda + mu), quantile(r2$continuousBranches[[i]][, "w"], qq), pch = 3,
         cex = 0.4, col = "#0072B2")
  abline(0, 1, col = "grey50")
  if (i == 1) legend("topleft", c("dCET (beta = gamma = 1)", "continuous CET",
                                  "control: beta = 0, gamma = 1"),
                     pch = c(16, 3, 1), col = c("black", "#0072B2", "#D55E00"), bty = "n")
}
dev.off()

# --- Figure 5: tempo distribution among extant species at T = 100 ---
Png("fig-tempo-distribution.png", width = 9, height = 7)
th <- r2$tempoHist
plot(th$x, th$simulated / sum(th$simulated), type = "h", lwd = 3, col = "grey60",
     xlab = "Log-tempo x", ylab = "Proportion of extant species", xlim = c(-3, 4),
     ylim = c(0, 1.3 * max(th$expected)),
     main = "Extant tempo, T = 100 Myr, x0 = 0")
lines(th$x, th$expected, type = "b", pch = 16, cex = 0.5)
legend("topleft", c("simulated (pooled)", "replicator-mutation ODE"), lwd = c(3, 1),
       col = c("grey60", "black"), bty = "n")
dev.off()

# --- Tables ---
cd <- bm$code
our <- bm$ours
# B&M's pdeMM2 adds dN/dt to the second-moment equation, so their C2 is
# E[n^2] + t dN/dt; adding t (Q + R) N / N to EvoTempo's value reproduces it.
bmRates <- BinRates(bmGen$x, lambda, mu)
bmN <- DiversityMoments(500, bmGen$Q, bmRates$lambda, bmRates$mu)$mean[1, ]
bmNdot <- as.vector((bmGen$Q + diag(bmRates$lambda - bmRates$mu)) %*% bmN)
experiencedAsCode <- our$experienced + 500 * bmNdot[51] / bmN[51]
co <- bm$converged
rung1 <- rbind(
  c("S0(500), survival from x0 = 0", "-", Fmt(cd$S0), Fmt(our$S0), Fmt(co$S0), Pct(cd$S0, co$S0)),
  c("Fig 1A mean clade size E(n | n > 0)", "c. 60 000", Fmt(cd$meanClade, 6), Fmt(our$meanClade, 6),
    Fmt(co$meanClade, 6), Pct(cd$meanClade, co$meanClade)),
  c("Fig 1A experienced size E[n^2]/E[n]", "c. 400 000", Fmt(cd$experienced, 7),
    paste0(Fmt(our$experienced, 7), " (+ t dN/dt / N, as in pdeMM2: ", Fmt(experiencedAsCode, 7), ")"),
    Fmt(co$experienced, 6), Pct(cd$experienced, co$experienced)),
  c("ratio experienced / mean", "c. 8", Fmt(cd$experienced / cd$meanClade, 3),
    Fmt(our$experienced / our$meanClade, 3), Fmt(co$experienced / co$meanClade, 3), ""),
  c("Fig 1A P(n < 1000 | n > 0)", "c. 1/3", "not run (evalQ)", Fmt(pnSum$code["lt1000"], 3),
    Fmt(pnSum$converged["lt1000"], 3), ""),
  c("Fig 1A P(n > 50 000 | n > 0)", "c. 1/4", "not run (evalQ)", Fmt(pnSum$code["gt50000"], 3),
    Fmt(pnSum$converged["gt50000"], 3), ""),
  c("Fig 1B P(smaller sister share <= 0.05)", "c. 0.5", "not run (evalQ)",
    sprintf("%s - %s", Fmt(pnSum$code["share"], 3), Fmt(pnSum$code["shareUpper"], 3)),
    sprintf("%s - %s", Fmt(pnSum$converged["share"], 3), Fmt(pnSum$converged["shareUpper"], 3)), ""),
  c("Fig 2A N_x(100), x0 = -2 / 0 / 2", "-", Fmt(cd$N100), Fmt(our$N100), Fmt(co$N100), ""),
  c("Fig 2A N_x(500), x0 = -2 / 0 / 2", "-", Fmt(cd$N500), Fmt(our$N500), Fmt(co$N500), ""),
  c("Fig 2A log-slope 400-500 Myr", "\"fixed rate\"", Fmt(cd$slope), Fmt(our$slope), Fmt(co$slope),
    Pct(cd$slope, co$slope)),
  c("Perron growth rate", "-", "-", Fmt(our$eqGrowth), Fmt(co$eqGrowth), ""),
  c("Fig 2B E(n_t | n_500 > 0), t = 100 / 250 / 400", "-",
    paste0(Fmt(cd$fig2bCode, 6), " (as published); ", Fmt(cd$fig2bFixed, 6), " (pdeC fixed)"),
    Fmt(our$survivor, 6), Fmt(co$survivor, 5), ""),
  c("Fig 2B lineages, t = 100 / 250 / 400", "-", Fmt(cd$lineages), Fmt(our$lineages), Fmt(co$lineages), ""),
  c("Fig 3A mean, x0 = -2, t = 10/50/100/200/500", "-", Fmt(cd$fig3aMinus, 3), Fmt(our$fig3aMinus, 3),
    Fmt(co$fig3aMinus, 3), ""),
  c("Fig 3A mean, x0 = +2, t = 10/50/100/200/500", "-", Fmt(cd$fig3aPlus, 3), Fmt(our$fig3aPlus, 3),
    Fmt(co$fig3aPlus, 3), ""),
  c("Fig 3A x0 = 0, t = 500: mean, sd", "mean 0", Fmt(c(round(cd$fig3aZero[1], 10), cd$fig3aZero[2]), 3),
    Fmt(c(round(our$fig3aZero[1], 10), our$fig3aZero[2]), 3), Fmt(c(round(co$fig3aZero[1], 10), co$fig3aZero[2]), 3), ""),
  c("Fig 3B mean at t = 500, x0 = 0: r = 0 / 2 theta", "-1 (A.26 limit)", Fmt(c(cd$fig3bLow, cd$fig3bHigh), 3),
    Fmt(c(our$fig3bLow, our$fig3bHigh), 3), Fmt(c(co$fig3bLow, co$fig3bHigh), 3), ""),
  c("Fig 4A mean ancestral x, look-back 100/200/500", "plateau c. 0.6", Fmt(cd$fig4a, 3), Fmt(our$fig4a, 3),
    Fmt(co$fig4a, 3), ""),
  c("Fig 5A mean branch duration, x0 = -2/0/2", "-", paste0(Fmt(cd$fig5mean20, 3), " (truncated at 20 Myr, renormalized)"),
    Fmt(our$durations, 4), Fmt(co$durations, 4), ""),
  c("Fig 5B mean branch w (every bin)", "1/(lambda+mu) = 0.990099", "analytic in code",
    Fmt(our$w, 7), Fmt(co$w, 7), "")
)
colnames(rung1) <- c("Quantity", "Paper text", "B&M code (verbatim re-run)",
                     "EvoTempo, B&M operator", "EvoTempo converged (central K = 401, z in +-10)",
                     "code vs converged")
MdTable <- function(m) {
  c(paste("|", paste(colnames(m), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(m)), collapse = " | "), "|"),
    apply(m, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
# Observed order from successive differences: |q_i - q_{i+1}| ~ h_i^p
Order <- function(d, q) {
  n <- nrow(d)
  if (n < 3L) return(rep(NA, n))
  e <- abs(diff(d[[q]]))
  c(NA, NA, log(e[-length(e)] / e[-1]) / log(d$h[seq_len(n - 2L)] / d$h[seq_len(n - 2L) + 1L]))
}
convRows <- do.call(rbind, lapply(split(cv, paste(cv$scheme, cv$zMax)), function(d) {
  d <- d[order(-d$h), ]
  data.frame(d[, c("scheme", "zMax", "h", "K", "S0", "meanClade", "experienced",
                   "fig4a100", "a26")],
             order = Order(d, "meanClade"))
}))
convTable <- as.matrix(transform(
  convRows, h = signif(h, 3), S0 = signif(S0, 5), meanClade = round(meanClade),
  experienced = round(experienced), fig4a100 = round(fig4a100, 4), a26 = round(a26, 4),
  order = round(order, 2)))
colnames(convTable) <- c("scheme", "z_max", "h", "K", "S0(500)", "mean clade", "experienced",
                         "Fig 4A @100", "A.26 mean (r = 0)", "observed order (mean clade)")

bt <- r2$branchTable
ct <- r2$cladeTable
colnames(ct) <- c("x0", "nClade", paste0("dcet.", colnames(ct)[3:9]),
                  paste0("cet.", colnames(ct)[10:16]))
nu <- r2$numerics
cn <- r2$convergedNumerics
Z <- function(sim, se, ref) sprintf("%.2f", (sim - ref) / se)
branchRows <- cbind(
  "x0" = bt$x0,
  "KS p, dCET" = signif(bt$ksP, 2), "KS p, continuous CET" = signif(bt$ksPcontinuous, 2),
  "KS p, control (beta = 0)" = signif(bt$ksPcontrol, 2),
  "mean w" = round(bt$meanW, 4),
  "mean duration: dCET sim" = sprintf("%.3f +- %.3f", bt$meanDuration, bt$seDuration),
  "numerics (K = 101)" = round(bt$expectedDuration, 4),
  "continuous sim" = sprintf("%.3f +- %.3f", bt$meanDurationContinuous, bt$seDurationContinuous),
  "converged (K = 401)" = round(cn$duration, 4)
)
cladeRows <- cbind(
  "x0" = ct[, "x0"], "clades" = ct[, "nClade"],
  "survival: dCET sim" = sprintf("%.4f +- %.4f", ct[, "dcet.survival"], ct[, "dcet.seSurvival"]),
  "S(100) numerics" = round(nu$survival, 4),
  "z" = Z(ct[, "dcet.survival"], ct[, "dcet.seSurvival"], nu$survival),
  "mean n: dCET sim" = sprintf("%.3f +- %.3f", ct[, "dcet.meanN"], ct[, "dcet.seMeanN"]),
  "N(100) numerics" = round(nu$meanN, 3),
  "z " = Z(ct[, "dcet.meanN"], ct[, "dcet.seMeanN"], nu$meanN),
  "extant mean x: sim" = sprintf("%.3f +- %.3f", ct[, "dcet.meanX"], ct[, "dcet.seMeanX"]),
  "numerics" = round(nu$meanX, 3),
  "z  " = Z(ct[, "dcet.meanX"], ct[, "dcet.seMeanX"], nu$meanX)
)
cetRows <- cbind(
  "x0" = ct[, "x0"],
  "survival: CET sim" = sprintf("%.4f +- %.4f", ct[, "cet.survival"], ct[, "cet.seSurvival"]),
  "dCET h = 0.4 / 0.2 / 0.1 / converged" = apply(cbind(r2$coarse[[1]]$survival, r2$coarse[[2]]$survival,
                                                    nu$survival, cn$survival), 1, Fmt, d = 4),
  "mean n: CET sim" = sprintf("%.3f +- %.3f", ct[, "cet.meanN"], ct[, "cet.seMeanN"]),
  "dCET h = 0.4 / 0.2 / 0.1 / converged " = apply(cbind(r2$coarse[[1]]$meanN, r2$coarse[[2]]$meanN,
                                                     nu$meanN, cn$meanN), 1, Fmt, d = 4),
  "z vs converged (survival, n)" = paste(Z(ct[, "cet.survival"], ct[, "cet.seSurvival"], cn$survival),
                                         Z(ct[, "cet.meanN"], ct[, "cet.seMeanN"], cn$meanN), sep = ", ")
)
dc <- r2$durationCheck
durationRows <- cbind(Simulator = dc$simulator, Branches = dc$n,
                      "Mean duration" = sprintf("%.5f +- %.5f", dc$mean, dc$se),
                      Expected = round(dc$expected, 5), z = Z(dc$mean, dc$se, dc$expected))
st <- r2$samplingTable
samplingRows <- cbind(Quantity = st$quantity,
                      Simulated = sprintf("%.4f +- %.4f", st$simulated, st$se),
                      Expected = round(st$expected, 4),
                      z = Z(st$simulated, st$se, st$expected))
saveRDS(list(rung1 = rung1, convergence = convTable, branch = branchRows, clade = cladeRows,
             cet = cetRows, sampling = samplingRows, durations = durationRows, pnSummary = pnSum),
        file.path(repDir, "tables.rds"))
writeLines(c(MdTable(rung1), "", MdTable(convTable), "", MdTable(branchRows), "",
             MdTable(durationRows), "", MdTable(cladeRows), "", MdTable(cetRows), "",
             MdTable(samplingRows)),
           file.path(repDir, "tables.md"))
