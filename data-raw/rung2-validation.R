# Validation rung 2: forward simulators vs dCET numerics (issue #2).
# Run from the package root with EvoTempo installed. Each section can run as a
# separate process: EVOTEMPO_SECTION = branches, dcet1..3, cet1..3, sampling;
# then EVOTEMPO_SECTION = combine assembles data-raw/rung2-results/rung2.rds.
library("EvoTempo")
outDir <- "data-raw/rung2-results"
dir.create(outDir, showWarnings = FALSE)
section <- Sys.getenv("EVOTEMPO_SECTION", "all")
Run <- function(name) section %in% c("all", name)
Progress <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
lambda <- 0.51
mu <- 0.5
theta <- 0.01
grid <- TempoGrid(101, c(-5, 5))
Q <- TempoGenerator(grid, theta)
rates <- BinRates(grid, lambda, mu)
starts <- c(-2, 0, 2)
startBin <- vapply(starts, function(x0) which.min(abs(grid$x - x0)), 1L)
out <- rates$lambda + rates$mu
horizon <- 100
nClade <- c(10000, 20000, 6000)
psi <- 0.05
rho <- 0.5

# 1-2. Complete-tree branches: w ~ Exp(lambda + mu); durations.
# A branch ends at its first speciation or extinction, so simulating with
# speciation folded into extinction gives the same branch law.
nBranch <- 4000
Branches <- function(n, k, beta = 1) {
  end <- (lambda + mu) * exp(beta * grid$x)
  t(vapply(seq_len(n), function(i) {
    h <- SimulateDCET(1e9, grid$x, Q, 0, end, start = k)
    c(w = BranchClock(h, 1), duration = h$edges$end[1])
  }, numeric(2)))
}
ContinuousBranches <- function(n, x0) {
  t(vapply(seq_len(n), function(i) {
    h <- SimulateCET(1e9, theta, 1, 0, lambda + mu, x0 = x0)
    c(w = BranchClock(h, 1), duration = h$edges$end[1])
  }, numeric(2)))
}
if (Run("branches")) {
  set.seed(1)
  branches <- lapply(startBin, Branches, n = nBranch)
  Progress("dCET branches")
  control <- lapply(startBin, Branches, n = nBranch, beta = 0)
  Progress("control branches")
  continuousBranches <- lapply(starts, ContinuousBranches, n = nBranch)
  Progress("continuous branches")
  saveRDS(list(branches = branches, control = control,
               continuousBranches = continuousBranches),
          file.path(outDir, "branches.rds"))
}

# Follow-up: x0 = 2 durations with 10x the branches (dCET) and 7.5x (CET)
if (Run("durations")) {
  set.seed(2)
  end <- (lambda + mu) * exp(grid$x)
  dcetD <- vapply(seq_len(40000), function(i) {
    SimulateDCET(1e9, grid$x, Q, 0, end, start = startBin[3])$edges$end[1]
  }, 0)
  cetD <- vapply(seq_len(30000), function(i) {
    SimulateCET(1e9, theta, 1, 0, lambda + mu, x0 = 2)$edges$end[1]
  }, 0)
  saveRDS(list(dcet = dcetD, cet = cetD), file.path(outDir, "durations.rds"))
}

# 3-4. Clades to T = 100: survival, mean diversity, tempo distribution
Clades <- function(n, Simulate) {
  lapply(seq_len(n), function(i) {
    if (i %% 1000 == 0) Progress(i)
    h <- Simulate()
    present <- h$edges$type == "present"
    endX <- vapply(h$paths[present], function(p) p$x[[length(p$x)]], 0)
    list(n = sum(present), x = endX)
  })
}
for (i in 1:3) {
  if (Run(paste0("dcet", i))) {
    set.seed(10 + i)
    clades <- Clades(nClade[i], function() {
      SimulateDCET(horizon, grid$x, Q, rates$lambda, rates$mu, start = startBin[i],
                   maxEdges = 1e7)
    })
    saveRDS(clades, file.path(outDir, paste0("dcet", i, ".rds")))
  }
  if (Run(paste0("cet", i))) {
    set.seed(20 + i)
    clades <- Clades(nClade[i], function() {
      SimulateCET(horizon, theta, 1, lambda, mu, x0 = starts[i], maxEdges = 1e7)
    })
    saveRDS(clades, file.path(outDir, paste0("cet", i, ".rds")))
  }
}

# 6. Sampling layer: fossils, rho, and the E-ODE
if (Run("sampling")) {
  set.seed(30)
  sampling <- do.call(rbind, lapply(seq_len(20000), function(i) {
    if (i %% 1000 == 0) Progress(i)
    h <- SimulateDCET(horizon, grid$x, Q, rates$lambda, rates$mu, psi,
                      start = startBin[2], rho = rho, maxEdges = 1e7)
    c(fossils = nrow(h$fossils), lineageTime = sum(h$edges$end - h$edges$start),
      anySample = nrow(h$fossils) > 0 || any(h$edges$sampled, na.rm = TRUE),
      anyExtant = any(h$edges$sampled, na.rm = TRUE))
  }))
  saveRDS(sampling, file.path(outDir, "sampling.rds"))
}

if (Run("combine")) {
fine <- TempoGrid(401, c(-10, 10))
fineQ <- TempoGenerator(fine, theta)
fineRates <- BinRates(fine, lambda, mu)
fineBin <- vapply(starts, function(x0) which.min(abs(fine$x - x0)), 1L)
Numerics <- function(x, Qm, bins) {
  r <- BinRates(x, lambda, mu)
  E <- ExtinctionProbability(horizon, Qm, r$lambda, r$mu)[1, bins]
  mom <- DiversityMoments(horizon, Qm, r$lambda, r$mu, second = TRUE)
  p <- vapply(bins, function(k) {
    d <- TempoDistribution(horizon, Qm, r$lambda, r$mu, start = k)
    c(sum(d * x), sqrt(sum(d * x^2) - sum(d * x)^2))
  }, numeric(2))
  durations <- BranchExpectation(Qm, r$lambda + r$mu)[bins]
  data.frame(survival = 1 - E, meanN = mom$mean[1, bins],
             sdN = sqrt(mom$meanSquare[1, bins] - mom$mean[1, bins]^2),
             meanX = p[1, ], sdX = p[2, ], duration = durations)
}
numerics <- Numerics(grid$x, Q, startBin)
convergedNumerics <- Numerics(fine$x, fineQ, fineBin)
coarse <- lapply(c(0.4, 0.2), function(h) {
  g <- TempoGrid(round(10 / h) + 1, c(-5, 5))
  Numerics(g$x, TempoGenerator(g, theta),
           vapply(starts, function(x0) which.min(abs(g$x - x0)), 1L))
})
Summarize <- function(clades) {
  n <- vapply(clades, `[[`, 0, "n")
  pooled <- unlist(lapply(clades, `[[`, "x"))
  boot <- replicate(200, {
    pick <- sample.int(length(clades), replace = TRUE)
    xs <- unlist(lapply(clades[pick], `[[`, "x"))
    mean(xs)
  })
  c(survival = mean(n > 0), seSurvival = sqrt(mean(n > 0) * (1 - mean(n > 0)) / length(n)),
    meanN = mean(n), seMeanN = sd(n) / sqrt(length(n)),
    meanX = mean(pooled), seMeanX = sd(boot), sdX = sd(pooled))
}
b <- readRDS(file.path(outDir, "branches.rds"))
branches <- b$branches
control <- b$control
continuousBranches <- b$continuousBranches
KS <- function(w) suppressWarnings(ks.test(w, "pexp", lambda + mu)$p.value)
branchTable <- data.frame(
  x0 = starts,
  ksP = vapply(branches, function(b) KS(b[, "w"]), 0),
  ksPpooled = KS(unlist(lapply(branches, function(b) b[, "w"]))),
  ksPcontrol = vapply(control, function(b) KS(b[, "w"]), 0),
  ksPcontinuous = vapply(continuousBranches, function(b) KS(b[, "w"]), 0),
  meanW = vapply(branches, function(b) mean(b[, "w"]), 0),
  meanDuration = vapply(branches, function(b) mean(b[, "duration"]), 0),
  seDuration = vapply(branches, function(b) sd(b[, "duration"]) / nBranch^0.5, 0),
  expectedDuration = BranchExpectation(Q, out)[startBin],
  meanDurationContinuous = vapply(continuousBranches, function(b) mean(b[, "duration"]), 0),
  seDurationContinuous = vapply(continuousBranches,
                                function(b) sd(b[, "duration"]) / nBranch^0.5, 0)
)
d2 <- readRDS(file.path(outDir, "durations.rds"))
durationCheck <- data.frame(
  simulator = c("dCET", "continuous CET"), n = c(length(d2$dcet), length(d2$cet)),
  mean = c(mean(d2$dcet), mean(d2$cet)),
  se = c(sd(d2$dcet) / sqrt(length(d2$dcet)), sd(d2$cet) / sqrt(length(d2$cet))),
  expected = c(BranchExpectation(Q, out)[startBin[3]],
               BranchExpectation(fineQ, fineRates$lambda + fineRates$mu)[fineBin[3]])
)
dcetClades <- lapply(1:3, function(i) readRDS(file.path(outDir, paste0("dcet", i, ".rds"))))
cetClades <- lapply(1:3, function(i) readRDS(file.path(outDir, paste0("cet", i, ".rds"))))
cladeTable <- cbind(
  x0 = starts, nClade = nClade,
  dcet = do.call(rbind, lapply(dcetClades, Summarize)),
  cet = do.call(rbind, lapply(cetClades, Summarize))
)
# Pooled tempo histogram from x0 = 0 vs numerics
tempoHist <- list(
  x = grid$x,
  simulated = tabulate(match(unlist(lapply(dcetClades[[2]], `[[`, "x")), grid$x),
                       length(grid$x)),
  expected = as.vector(TempoDistribution(horizon, Q, rates$lambda, rates$mu,
                                         start = startBin[2]))
)

sampling <- readRDS(file.path(outDir, "sampling.rds"))
eAny <- ExtinctionProbability(horizon, Q, rates$lambda, rates$mu, psi, rho)[1, startBin[2]]
eExtant <- ExtinctionProbability(horizon, Q, rates$lambda, rates$mu, 0, rho)[1, startBin[2]]
samplingTable <- data.frame(
  quantity = c("fossils / (psi x lineage time)", "P(any sample)", "P(sampled extant)"),
  simulated = c(sum(sampling[, "fossils"]) / (psi * sum(sampling[, "lineageTime"])),
                mean(sampling[, "anySample"]), mean(sampling[, "anyExtant"])),
  se = c(1 / sqrt(sum(sampling[, "fossils"])),
         sqrt(mean(sampling[, "anySample"]) * (1 - mean(sampling[, "anySample"])) / nrow(sampling)),
         sqrt(mean(sampling[, "anyExtant"]) * (1 - mean(sampling[, "anyExtant"])) / nrow(sampling))),
  expected = c(1, 1 - eAny, 1 - eExtant)
)
saveRDS(list(branchTable = branchTable, branches = branches, control = control,
             continuousBranches = continuousBranches, cladeTable = cladeTable,
             numerics = numerics, convergedNumerics = convergedNumerics,
             coarse = coarse, tempoHist = tempoHist, samplingTable = samplingTable,
             durationCheck = durationCheck),
        file.path(outDir, "rung2.rds"))
}
