# Validation 3 of #97: SBC-lite (Talts et al. 2018) of TempoMCMC() with
# Gamma rates among characters and the Mkv correction, in cells
# (beta, gamma) = (1, 0) and (1, 1).
#
#   Rscript data-raw/gamma-acrv-sbc.R run <task> <nTask> <outDir> [nRep] [nIter]
#   Rscript data-raw/gamma-acrv-sbc.R summary <outDir>
#
# As rung 6 (data-raw/rung6-sbc.R): each replicate draws the parameters from
# the prior, simulates a crown clade of 8 time units with the exact dCET
# simulator, resimulating at fixed parameters until both crown daughters have
# sampled extant descendants, and keeps trees with 15-40 tips.  Then, new
# here: alpha ~ LogNormal(0, 0.5) and 45 characters (30 binary, 15
# three-state, unordered) on the realized tempo-weighted lengths, each with a
# Gamma category drawn uniformly from GammaRates(alpha, 4), redrawn
# (category and states) until the character is variable, which is the
# generative model of MkCharacters(coding = "variable", nRateCat = 4).
# TempoMCMC() runs on the fixed true tree with every parameter and alpha
# free, from an independent prior draw.  Ranks of 50 thinned draws: each
# parameter, alpha, the root bin and the path functional sum_e int_e e^x dt.
# Replicates g = task, task + nTask, ...; cell = (g - 1) %% 2 + 1;
# seed 97000 + g.  `summary` writes ranks.csv and a chi-square test (10 bins)
# per statistic and cell.
# Options (environment): SBC_CELLS, cells to run (default "1,2");
# SBC_NCAT, Gamma categories (default 4; 1 = no rate variation: characters
# simulated at rate 1, nRateCat = 1, alpha not estimated, a matched control
# with the same trees, since the prior draw is the same).  The chain keeps
# every iteration; ranks use draws thin, 2 thin, ..., nIter (thin = nIter /
# 50), and ESS (Geyer) is reported on the whole chain (`ess`) and on the 50
# draws (`essDraws`).
suppressMessages(library(EvoTempo))
ns <- asNamespace("EvoTempo")
args <- commandArgs(TRUE)
cells <- data.frame(beta = c(1, 1), gamma = c(0, 1))
nDraw <- 50L
statNames <- c("lambda", "mu", "psi", "theta", "s", "rate", "alpha", "rootBin", "clock")
nCat <- as.integer(Sys.getenv("SBC_NCAT", "4"))
cellsRun <- as.integer(strsplit(Sys.getenv("SBC_CELLS", "1,2"), ",")[[1]])

if (args[[1]] == "summary") {
  outDir <- args[[2]]
  files <- list.files(file.path(outDir, "replicates"), "\\.rds$", full.names = TRUE)
  reps <- lapply(files, readRDS)
  ranks <- do.call(rbind, lapply(reps, function(r) {
    data.frame(g = r$g, cell = r$cell, beta = cells$beta[r$cell], gamma = cells$gamma[r$cell],
               nTip = r$nTip, t(r$ranks), nIter = if (is.null(r$nIter)) NA else r$nIter,
               essTheta = r$ess[["theta"]], essS = r$ess[["s"]],
               essAlpha = if ("alpha" %in% names(r$ess)) r$ess[["alpha"]] else NA,
               essMin = min(r$ess, na.rm = TRUE), minutes = r$runtime)
  }))
  ranks <- ranks[order(ranks$g), ]
  write.csv(ranks, file.path(outDir, "ranks.csv"), row.names = FALSE)
  tests <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ci) {
    sub <- ranks[ranks$cell == ci, ]
    do.call(rbind, lapply(statNames[colSums(!is.na(sub[statNames])) > 0], function(p) {
      bins <- cut(sub[[p]], seq(-0.5, nDraw + 0.5, length.out = 11))
      counts <- table(bins)
      chi <- suppressWarnings(stats::chisq.test(counts))
      data.frame(beta = cells$beta[ci], gamma = cells$gamma[ci], statistic = p,
                 nRep = nrow(sub), chisq = unname(chi$statistic), p = chi$p.value,
                 meanRank = mean(sub[[p]]) / nDraw, minBin = min(counts),
                 maxBin = max(counts))
    }))
  }))
  tests$pHolm <- stats::p.adjust(tests$p, "holm")
  print(tests, digits = 3, row.names = FALSE)
  write.csv(tests, file.path(outDir, "sbc-tests.csv"), row.names = FALSE)
  cat(sprintf("%d replicates; min Holm p %.3f; median ESS (whole chain) theta %.0f, s %.0f, alpha %.0f; %.1f core-h\n",
              nrow(ranks), min(tests$pHolm), stats::median(ranks$essTheta),
              stats::median(ranks$essS), stats::median(ranks$essAlpha, na.rm = TRUE),
              sum(ranks$minutes) / 60))
  for (ci in unique(ranks$cell)) {
    sub <- ranks[ranks$cell == ci, ]
    cat(sprintf("cell %d: ESS theta quantiles (5, 25, 50, 75, 95%%): %s\n", ci,
                paste(round(stats::quantile(sub$essTheta, c(.05, .25, .5, .75, .95))), collapse = " ")))
  }
  quit(save = "no")
}

task <- as.integer(args[[2]])
nTask <- as.integer(args[[3]])
outDir <- args[[4]]
nRep <- if (length(args) > 4) as.integer(args[[5]]) else 100L
nIter <- if (length(args) > 5) as.integer(args[[6]]) else 3000L
burnin <- 600L
thin <- nIter %/% nDraw
crown <- 8
tipRange <- c(15, 40)
dir.create(file.path(outDir, "replicates"), recursive = TRUE, showWarnings = FALSE)

meanlog <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.05), theta = log(0.2),
             s = log(0.6), rate = log(0.1))
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.25, rate = 0.4)
priors <- TempoPriors(meanlog, sdlog)
alphaPrior <- function(a) stats::dlnorm(a, 0, 0.5, log = TRUE)
PriorDraw <- function() {
  c(as.list(setNames(stats::rlnorm(6, meanlog, sdlog), names(meanlog))),
    alpha = stats::rlnorm(1, 0, 0.5))
}

Join <- function(h1, h2, x0) {
  n1 <- nrow(h1$edges)
  e1 <- h1$edges
  e2 <- h2$edges
  e1$parent <- ifelse(is.na(e1$parent), 1L, e1$parent + 1L)
  e2$parent <- ifelse(is.na(e2$parent), 1L, e2$parent + 1L + n1)
  f1 <- h1$fossils
  f2 <- h2$fossils
  f1$edge <- f1$edge + 1L
  f2$edge <- f2$edge + 1L + n1
  structure(list(edges = rbind(data.frame(parent = NA_integer_, start = 0, end = 0,
                                          type = "speciation", sampled = NA), e1, e2),
                 paths = c(list(list(t = 0, x = x0)), h1$paths, h2$paths),
                 fossils = rbind(f1, f2), time = h1$time),
            class = "TempoHistory", kind = "discrete")
}

# Geyer initial monotone sequence ESS (coda is not on Hamilton)
Ess <- function(x) {
  n <- length(x)
  if (stats::var(x) == 0) return(NA_real_)
  ac <- stats::acf(x, lag.max = min(n - 1, 4000), plot = FALSE)$acf[, 1, 1]
  m <- floor((length(ac) - 1) / 2)
  gam <- ac[2 * seq_len(m) - 1] + ac[2 * seq_len(m)]
  k <- which(gam <= 0)[1]
  gam <- gam[seq_len(if (is.na(k)) m else k - 1)]
  gam <- cummin(gam)
  n / (-1 + 2 * sum(gam))
}

# Variable characters with a Gamma category each, on per-edge lengths
SimulateGammaCharacters <- function(tree, nChar, nState, alpha) {
  rates <- if (nCat > 1) GammaRates(alpha, nCat) else rep(1, 4)
  out <- vapply(seq_len(nChar), function(j) {
    repeat {
      r <- rates[sample.int(4, 1)]
      x <- SimulateCharacters(tree, 1, nState, lengths = r * attr(tree, "characterLength"))
      if (length(unique(x[, 1])) > 1) return(x[, 1])
    }
  }, integer(length(tree$tip.label)))
  out
}

Replicate <- function(g) {
  ci <- (g - 1) %% nrow(cells) + 1
  cell <- cells[ci, ]
  model <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5, beta = cell$beta,
                     gamma = cell$gamma)
  set.seed(97000 + g)
  repeat {
    truth <- PriorDraw()
    set <- ns$.Setup(truth, model)
    repeat {
      k <- ns$.Sample1(set$prior)
      h <- tryCatch(lapply(1:2, function(i) {
        SimulateDCET(crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                     start = k, rho = model$rho, maxEdges = 5000)
      }), error = function(e) NULL)
      if (is.null(h) ||
          all(vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0) > 0)) {
        break
      }
    }
    if (is.null(h)) next
    history <- Join(h[[1]], h[[2]], set$grid$x[k])
    tree <- ReconstructedTree(history, gamma = cell$gamma, rate = truth$rate)
    nTip <- sum(attr(tree, "tipType") != "ancestor")
    if (nTip >= tipRange[1] && nTip <= tipRange[2]) break
  }
  chars <- cbind(SimulateGammaCharacters(tree, 30, 2L, truth$alpha),
                 SimulateGammaCharacters(tree, 15, 3L, truth$alpha))
  nState <- rep(c(2L, 3L), c(30, 15))
  mk <- MkCharacters(chars, nState = nState, coding = "variable", nRateCat = nCat)
  clockTrue <- sum(attr(ReconstructedTree(history, gamma = 1, rate = 1),
                        "characterLength"))
  tree$root.edge <- NULL
  tt <- TempoTree(tree, present = crown)
  t0 <- Sys.time()
  run <- TempoMCMC(tt, model, PriorDraw(), priors, nIter = nIter, burnin = burnin,
                   thin = 1, nSlice = 16,
                   moves = c(tree = 1, node = 3, tip = 3, branch = 3, shift = 1),
                   characters = mk,
                   charPriors = if (nCat > 1) list(alpha = alphaPrior) else list(),
                   minAcceptance = 0)
  full <- run$samples
  post <- full[seq(thin, nIter, by = thin), ]
  essNames <- c(names(meanlog), if (nCat > 1) "alpha", "clock")
  rankOf <- function(draws, value) {
    sum(draws < value) + sample.int(sum(draws == value) + 1L, 1L) - 1L
  }
  truthAll <- c(unlist(truth), rootBin = k, clock = clockTrue)
  ranks <- vapply(statNames, function(p) {
    if (p == "alpha" && nCat == 1) NA_real_ else rankOf(post[[p]], truthAll[[p]])
  }, 0)
  out <- list(g = g, cell = ci, truth = truthAll, ranks = ranks, nTip = nTip,
              nAncestor = sum(attr(tree, "tipType") == "ancestor"),
              nIter = nrow(full), nCat = nCat, acceptance = run$acceptance,
              ess = vapply(full[essNames], Ess, 0),
              essDraws = vapply(post[essNames], Ess, 0),
              posteriorMean = colMeans(post[essNames]),
              runtime = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, file.path(outDir, "replicates", sprintf("g%03d.rds", g)))
  out
}

gs <- seq(task, 2L * nRep, by = nTask)
gs <- gs[((gs - 1L) %% nrow(cells) + 1L) %in% cellsRun]
gs <- gs[!file.exists(file.path(outDir, "replicates", sprintf("g%03d.rds", gs)))]
invisible(parallel::mclapply(gs, function(g) {
  tryCatch(Replicate(g), error = function(e) {
    writeLines(conditionMessage(e), file.path(outDir, sprintf("error-g%03d.txt", g)))
  })
}, mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), mc.preschedule = FALSE))
cat("task", task, "done:", length(gs), "replicates\n")
