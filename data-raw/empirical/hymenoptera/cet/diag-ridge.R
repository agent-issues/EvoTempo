# Hymenoptera CET pipeline: diagnosis of the runaway of cell (1, 0) in the
# pilot (README, "Smoke run and cost projection"; PR #115).
#
#   Rscript data-raw/empirical/hymenoptera/cet/diag-ridge.R <pilot dir> <penguin runs dir> <out dir>
#
# Run from the package root with EvoTempo installed.  Writes to <out dir>:
#   ridge-location.csv  where the tempo process puts log-tempo x, as a function
#                       of s: the stationary law of Q = diag(e^(alpha x)) A
#                       (A the OU operator in z = x / s, mean 0) and the
#                       equilibrium root-bin prior (DCETModel(rootPrior =
#                       "equilibrium")) at the pilot's rates: mean x, and
#                       E[e^(beta x)] under each, against -s^2 and e^(-s^2/2)
#                       (the continuous-space values for alpha = beta = 1)
#   ridge-pilot.csv     the (1, 0) pilot chain by segment: lambda, mu, theta,
#                       s, mean path x, root bin, <e^x> = clock / treeLength
#                       (tree-length-weighted mean tempo factor, from
#                       .TreeStats()), lambda <e^x>, theta <e^x>, and the
#                       stationary-law effective rate lambda E_pi[e^x]
#   ridge-profile.csv   marginal tree likelihood (TreeLikelihood(), paths
#                       integrated out; atol 1e-60, #90) on one fixed tree
#                       (the chain's last recorded tree) along the straight
#                       line in log-parameter space from the segment-1 state
#                       to the last state (lambda, mu, psi, theta, s), with
#                       the log prior of the pipeline; and at default atol,
#                       for #90
#   ridge-penguin.csv   the penguin pipeline's (1, 0) production chains by
#                       window of segments: lambda, s, theta, mean x
# Each table is written as it is computed (a stalled evaluation keeps the rest).
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/hymenoptera/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
pilot <- args[[1]]
penguin <- args[[2]]
outDir <- args[[3]]
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
Out <- function(x, f) utils::write.csv(x, file.path(outDir, f), row.names = FALSE)
Samples <- function(d) do.call(rbind, lapply(sort(Sys.glob(file.path(d, "seg*_samples.csv"))), utils::read.csv))

# ---- 1. Location of the tempo process ---------------------------------------------
model <- CetModel(1, 0)
loc <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3), function(s) {
  grid <- TempoGrid(model$nBin, model$zRange, s)
  Q <- TempoGenerator(grid, 0.01, model$alpha, model$scheme)
  st <- TempoEquilibrium(Q)$p
  rates <- BinRates(grid, 40, 40 - 0.0055, 2e-6, 1, 0)
  eq <- get(".RootPrior", .ns)("equilibrium", Q, rates)
  data.frame(s = s, meanXStationary = sum(st * grid$x), minusS2 = -s^2,
             eXStationary = sum(st * exp(grid$x)), expMinusS2Half = exp(-s^2 / 2),
             meanXEquilibrium = sum(eq * grid$x), eXEquilibrium = sum(eq * exp(grid$x)),
             lowestBinMass = st[1])
}))
print(loc, digits = 4, row.names = FALSE)
Out(loc, "ridge-location.csv")
EStat <- function(s) {
  grid <- TempoGrid(model$nBin, model$zRange, s)
  st <- TempoEquilibrium(TempoGenerator(grid, 0.01, model$alpha, model$scheme))$p
  sum(st * exp(grid$x))
}

# ---- 2. The pilot chain by segment --------------------------------------------------
d10 <- file.path(pilot, "b1g0_run1")
x <- Samples(d10)
seg <- do.call(rbind, lapply(split(x, x$segment), function(y) {
  ex <- y$clock / y$treeLength
  data.frame(segment = y$segment[1], lambda = stats::median(y$lambda), mu = stats::median(y$mu),
             theta = stats::median(y$theta), s = stats::median(y$s), meanX = mean(y$meanX),
             rootBin = stats::median(y$rootBin), eX = stats::median(ex),
             lambdaEX = stats::median(y$lambda * ex), thetaEX = stats::median(y$theta * ex),
             lambdaStationary = stats::median(y$lambda) * EStat(stats::median(y$s)),
             rootAge = stats::median(y$rootAge))
}))
print(seg, digits = 4, row.names = FALSE)
Out(seg, "ridge-pilot.csv")

# ---- 4. Penguin (1, 0) production chains ---------------------------------------------
pen <- do.call(rbind, lapply(sort(Sys.glob(file.path(penguin, "b1g0_run*"))), function(d) {
  y <- Samples(d)
  if (is.null(y) || !nrow(y)) return(NULL)
  w <- cut(y$iteration, breaks = unique(c(0, seq(1e5, max(y$iteration) + 1e5, by = 1e5))))
  do.call(rbind, lapply(split(y, w, drop = TRUE), function(z) {
    data.frame(chain = basename(d), toIteration = max(z$iteration), n = nrow(z),
               lambdaMedian = stats::median(z$lambda), lambdaMax = max(z$lambda),
               s = stats::median(z$s), theta = stats::median(z$theta), meanX = mean(z$meanX),
               lambdaEX = stats::median(z$lambda * z$clock / z$treeLength))
  }))
}))
print(pen, digits = 4, row.names = FALSE)
Out(pen, "ridge-penguin.csv")

# ---- 3. Profile of the marginal tree likelihood on a fixed tree ----------------------
states <- unlist(lapply(sort(Sys.glob(file.path(d10, "seg*_states.rds"))), readRDS), recursive = FALSE)
last <- states[[length(states)]]
first <- states[[1]]
tree <- last$tree
pars <- c("lambda", "mu", "psi", "theta", "s")
p0 <- unlist(first$par[pars])
p1 <- unlist(last$par[pars])
priors <- CetPriors()
rows <- list()
for (t in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
  p <- as.list(exp((1 - t) * log(p0) + t * log(p1)))
  # keep lambda - mu on the line in log d, so that d stays positive
  d <- exp((1 - t) * log(p0[["lambda"]] - p0[["mu"]]) + t * log(p1[["lambda"]] - p1[["mu"]]))
  p$mu <- p$lambda - d
  t0 <- proc.time()[["elapsed"]]
  ll <- tryCatch(TreeLikelihood(tree, p, model, atol = 1e-60), error = function(e) NA_real_)
  secs <- proc.time()[["elapsed"]] - t0
  llDefault <- tryCatch(TreeLikelihood(tree, p, model), error = function(e) NA_real_)
  lp <- LogPriorJoint(c(p, list(rate = last$par$rate)), priors) - priors$rate(last$par$rate)
  rows[[length(rows) + 1L]] <- data.frame(t = t, lambda = p$lambda, mu = p$mu, psi = p$psi,
    theta = p$theta, s = p$s, lambdaStationary = p$lambda * EStat(p$s), logLik = ll,
    logPrior = lp, logLikPlusPrior = ll + lp, logLikDefaultAtol = llDefault, secs = secs)
  Out(do.call(rbind, rows), "ridge-profile.csv")
  print(utils::tail(do.call(rbind, rows), 1), digits = 8, row.names = FALSE)
}
