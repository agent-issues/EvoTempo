# #72 / verify43c: pre-registered analysis (prereg-verify43c.md). Written before any chain data.
args <- commandArgs(trailingOnly = TRUE)
pc <- read.csv(args[[1]]); pr <- read.csv(args[[2]]); fe <- read.csv(args[[3]])
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95); Fcol <- paste0("F", c("05", "25", "50", "75", "95"))
planned <- c(A = 40, B = 24, C = 16)
arms <- list()
for (arm in names(planned)) {
  d <- pc[pc$arm == arm, ]; f <- fe[fe$arm == arm, ]
  N <- nrow(d); nRun <- length(unique(d$offset))
  trapScript <- unique(signif(d$exactMean, 12))
  stopifnot(length(trapScript) == 1)
  check <- f$trapMean - d$exactMean[1]           # transcription check, should be ~0
  useSimpson <- abs(f$simpsonMean - d$exactMean[1]) > 0.0002
  ref <- if (useSimpson) f$simpsonMean else d$exactMean[1]
  b <- mean(d$mean) - ref; SE <- sd(d$mean) / sqrt(N); z <- b / SE
  sdx <- f$trapSd; r <- b / sdx; SEr <- SE / sdx
  zp <- vapply(seq_along(probs), function(k) {
    F <- d[[Fcol[k]]]; (mean(F) - probs[k]) / (sd(F) / sqrt(N)) }, 0)
  zr <- tapply(d$mean, d$offset, function(m) (mean(m) - ref) / (sd(m) / sqrt(length(m))))
  ks <- suppressWarnings(ks.test(as.numeric(zr), "pt", df = 31))
  arms[[arm]] <- data.frame(arm = arm, case = d$case[1], runsPlanned = planned[[arm]], runsDone = nRun,
    chains = N, underpowered = nRun < 0.75 * planned[[arm]], exactTrap = d$exactMean[1],
    exactSimpson = f$simpsonMean, simpsonMinusTrap = f$simpsonMean - d$exactMean[1],
    transcriptionCheck = check, refUsed = if (useSimpson) "simpson" else "trapezoid",
    exactSd = sdx, mcmcMean = mean(d$mean), b = b, SE = SE, z = z,
    ci99lo = b - 2.576 * SE, ci99hi = b + 2.576 * SE, r = r, SEr = SEr,
    z05 = zp[1], z25 = zp[2], z50 = zp[3], z75 = zp[4], z95 = zp[5],
    secondaryFlag = any(abs(zp) >= 3.5), ksD = unname(ks$statistic), ksP = ks$p.value,
    meanRunZ = mean(zr), ownFinding = abs(z) >= 3)
}
A <- arms$A; B <- arms$B; C <- arms$C
res <- do.call(rbind, arms)
# Rule 1
rule1 <- if (abs(A$z) >= 3) "REAL" else if (A$ci99lo > -0.0010 || A$ci99hi < -0.0010)
  "NO BIAS OF THE REPORTED SIZE" else "INCONCLUSIVE"
zDiff <- (A$r - B$r) / sqrt(A$SEr^2 + B$SEr^2)
bound <- if (abs(B$z) < 2 && abs(zDiff) >= 2.5) "localises to the bound region" else
  if (abs(B$z) >= 3 && sign(B$z) == sign(A$z)) "not the bound (target or move in general)" else
  "localisation inconclusive"
start <- if (abs(C$z) >= 3 && sign(C$z) == sign(A$z)) "also at the origin start (not specific to the rootAge move)" else
  if (abs(C$z) < 2) "root-start specific" else "inconclusive"
verdict <- data.frame(
  rule = c("1 primary (arm A)", "2a bound (A vs B)", "2b start type (C)", "3 secondary flags",
           "4 arm B own", "4 arm C own", "5 completion"),
  applies = c(TRUE, rule1 == "REAL", rule1 == "REAL", TRUE, TRUE, TRUE, TRUE),
  outcome = c(rule1, bound, start,
              paste(res$arm[res$secondaryFlag], collapse = ",") |> (\(x) if (x == "") "none" else x)(),
              if (B$ownFinding) "NEW FINDING |z|>=3" else "no finding",
              if (C$ownFinding) "NEW FINDING |z|>=3" else "no finding",
              paste(sprintf("%s %d/%d%s", res$arm, res$runsDone, res$runsPlanned,
                            ifelse(res$underpowered, " UNDERPOWERED", "")), collapse = "; ")),
  stat = c(sprintf("z_A %.3f; b_A %.6f +- %.6f; 99%% CI (%.6f, %.6f)", A$z, A$b, A$SE, A$ci99lo, A$ci99hi),
           sprintf("z_B %.3f; zDiff (r_A - r_B) %.3f", B$z, zDiff),
           sprintf("z_C %.3f", C$z), "", sprintf("z_B %.3f", B$z), sprintf("z_C %.3f", C$z), ""))
write.csv(res, "bias72-arms.csv", row.names = FALSE)
write.csv(verdict, "bias72-verdict.csv", row.names = FALSE)
print(t(res), quote = FALSE); print(verdict, right = FALSE)
