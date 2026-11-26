# Tempo location (#116), check V0 of dev/spec/tempo-location.md: the
# Hymenoptera (1, 0) ridge profile under the mean-rate model.
#
#   Rscript data-raw/tempo-location-v0.R <pilot chain dir> <t> <out dir>
#
# <pilot chain dir> holds the (1, 0) pilot's seg*_states.rds
# (hymenoptera-cet-pipeline, Slurm pilot; on Hamilton
# /nobackup/pjjg18/evotempo-hym-cet/pilot/b1g0_run1).  <t> in [0, 1] is the
# point on the straight line in log-parameter space from the first to the
# last recorded state, as in that branch's diag-ridge.R ("ridge-profile.csv",
# whose old-model timings this check is compared with).  The parameters there
# are Budd & Mann's (location "bm"); they are mapped to the mean-rate model by
# spec (L.3): lambda, mu and theta times exp(-s^2 / 2) (alpha = beta = 1,
# delta = 0).  Writes <out dir>/v0-<t>.csv with, for the pilot's model
# (K = 26, z in [-5, 5], rho = 5e-4, extant conditioning, root start):
#   logLik, secs            mean-rate model, default atol
#   logLik60, secs60        the same with atol = 1e-60 (as ridge-profile.csv)
#   logLikBm                Budd & Mann on the window shifted by -s (must equal
#                           logLik: spec L.3)
#   logLikRef, secsRef      reference grid: K = 51, z in [-7, 7]
#   topPrior, top3Prior     mass of the top bin / top three bins of the
#                           equilibrium root prior p*
#   topTilt, top3Tilt       the same for the e^x-tilted lineage law
#   lowPrior                mass of the bottom three bins of p*
suppressMessages(library(EvoTempo))
args <- commandArgs(trailingOnly = TRUE)
chain <- args[[1]]
t <- as.numeric(args[[2]])
outDir <- args[[3]]
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

states <- unlist(lapply(sort(Sys.glob(file.path(chain, "seg*_states.rds"))), readRDS),
                 recursive = FALSE)
first <- states[[1]]
last <- states[[length(states)]]
tree <- last$tree
pars <- c("lambda", "mu", "psi", "theta", "s")
p0 <- unlist(first$par[pars])
p1 <- unlist(last$par[pars])
p <- as.list(exp((1 - t) * log(p0) + t * log(p1)))
d <- exp((1 - t) * log(p0[["lambda"]] - p0[["mu"]]) + t * log(p1[["lambda"]] - p1[["mu"]]))
p$mu <- p$lambda - d
s <- p$s
pMean <- p
for (q in c("lambda", "mu", "theta")) pMean[[q]] <- p[[q]] * exp(-s^2 / 2)

Model <- function(nBin, zRange, location, scheme = NULL) {
  DCETModel(nBin = nBin, zRange = zRange, beta = 1, gamma = 0, rho = 5e-4,
            condition = "extant", start = "root", location = location,
            scheme = scheme)
}
Timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- tryCatch(expr, error = function(e) NA_real_)
  c(value, proc.time()[["elapsed"]] - t0)
}
model <- Model(26, c(-5, 5), "mean")
a <- Timed(TreeLikelihood(tree, pMean, model))
a60 <- Timed(TreeLikelihood(tree, pMean, model, atol = 1e-60))
bm <- TreeLikelihood(tree, p, Model(26, c(-5, 5) - s, "bm", "sqra"))
ref <- Timed(TreeLikelihood(tree, pMean, Model(51, c(-7, 7), "mean")))
set <- EvoTempo:::.Setup(pMean, model)
prior <- set$prior
tilt <- TempoEquilibrium(set$Q)$p * exp(set$grid$x)
tilt <- tilt / sum(tilt)
K <- length(prior)
out <- data.frame(
  t = t, s = s, lambdaBm = p$lambda, lambda = pMean$lambda, mu = pMean$mu,
  theta = pMean$theta, psi = p$psi,
  topRateBm = p$lambda * exp(5 * s), topRate = pMean$lambda * exp(5 * s - s^2 / 2),
  logLik = a[[1]], secs = a[[2]], logLik60 = a60[[1]], secs60 = a60[[2]],
  logLikBm = bm, logLikRef = ref[[1]], secsRef = ref[[2]],
  topPrior = prior[[K]], top3Prior = sum(prior[K - 0:2]), lowPrior = sum(prior[1:3]),
  topTilt = tilt[[K]], top3Tilt = sum(tilt[K - 0:2]),
  sha = Sys.getenv("EVOTEMPO_SHA")
)
print(out, digits = 8, row.names = FALSE)
utils::write.csv(out, file.path(outDir, sprintf("v0-%.1f.csv", t)), row.names = FALSE)
