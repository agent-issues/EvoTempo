# Shift move (#19) at fixed parameters: the root-bin posterior given the
# rung 4b tree and its character, from TempoMCMC() chains whose path moves
# include `shift` (bins only: no parameter is free, so nothing is compensated)
# and from reference chains without it, against the exact posterior.
# The exact root-bin posterior comes from the product-space likelihood with a
# point-mass start prior e_k: TreeLikelihood() then returns
# log p(T, c | k) - log S_k, so p(k | T, c) is proportional to
# pi_k S_k exp(L_k) under the model's own prior pi and marginal conditioning.
# Run from the package root:
#   Rscript data-raw/rung4b-shift-fixed.R [nIter] [nChain]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args)) as.integer(args[[1]]) else 20000L
nChain <- if (length(args) > 1) as.integer(args[[2]]) else 16L
burnin <- 500L
nSlice <- 20L
r4b <- readRDS("data-raw/rung5-results/rung4b.rds")
tt <- r4b$tt; model <- r4b$model; truth <- r4b$truth; character <- r4b$character
K <- model$nBin

# Exact posterior of the root bin
set <- .Precompute(tt, truth, model)
logL <- vapply(seq_len(K), function(k) {
  m <- model
  m$rootPrior <- replace(numeric(K), k, 1)
  TreeLikelihood(tt, truth, m, character = character)
}, 0)
lw <- log(set$prior) + log(set$survive) + logL
exact <- exp(lw - max(lw))
exact <- exact / sum(exact)

configs <- list(
  shift = c(node = 2, tip = 2, branch = 2, shift = 2),
  reference = c(tree = 1, node = 2, tip = 2, branch = 2)
)
jobs <- expand.grid(chain = seq_len(nChain), config = names(configs),
                    stringsAsFactors = FALSE)
cl <- makeCluster(min(16L, nrow(jobs)))
invisible(clusterCall(cl, function() suppressMessages(pkgload::load_all(".", quiet = TRUE))))
clusterExport(cl, c("tt", "model", "truth", "character", "configs", "jobs", "nIter",
                    "burnin", "nSlice"))
t0 <- Sys.time()
runs <- clusterApplyLB(cl, seq_len(nrow(jobs)), function(j) {
  set.seed(19400 + j)
  run <- TempoMCMC(tt, model, truth, list(), nIter = nIter, burnin = burnin,
                   nSlice = nSlice, moves = configs[[jobs$config[j]]],
                   characters = matrix(character), minAcceptance = 0)
  list(rootBin = run$samples$rootBin, meanX = run$samples$meanX,
       acceptance = run$acceptance)
})
stopCluster(cl)
runtime <- difftime(Sys.time(), t0, units = "mins")

# Batch-means variance of a chain mean
BatchVar <- function(x, nBatch = 25) {
  m <- floor(length(x) / nBatch) * nBatch
  stats::var(colMeans(matrix(x[seq_len(m)], ncol = nBatch))) / nBatch
}
Pooled <- function(x) {
  n <- vapply(x, length, 0)
  c(mean = sum(vapply(x, sum, 0)) / sum(n),
    se = sqrt(sum(vapply(x, BatchVar, 0) * n^2)) / sum(n))
}
rows <- list()
for (cfg in names(configs)) {
  ch <- runs[jobs$config == cfg]
  bins <- lapply(ch, `[[`, "rootBin")
  for (k in which(exact > 1e-3)) {
    pm <- Pooled(lapply(bins, function(b) as.numeric(b == k)))
    rows[[length(rows) + 1L]] <- data.frame(config = cfg, quantity = paste0("P(rootBin=", k, ")"),
                                            exact = exact[[k]], mcmc = pm[["mean"]], mcse = pm[["se"]])
  }
  pm <- Pooled(bins)
  rows[[length(rows) + 1L]] <- data.frame(config = cfg, quantity = "E(rootBin)",
                                          exact = sum(seq_len(K) * exact), mcmc = pm[["mean"]],
                                          mcse = pm[["se"]])
}
out <- do.call(rbind, rows)
out$z <- (out$mcmc - out$exact) / out$mcse
# meanX between the two samplers (no exact value)
mx <- lapply(names(configs), function(cfg) Pooled(lapply(runs[jobs$config == cfg], `[[`, "meanX")))
meanX <- data.frame(config = names(configs), mean = vapply(mx, `[[`, 0, "mean"),
                    mcse = vapply(mx, `[[`, 0, "se"))
meanXz <- diff(meanX$mean) / sqrt(sum(meanX$mcse^2))
acceptance <- do.call(rbind, lapply(names(configs), function(cfg) {
  a <- colMeans(do.call(rbind, lapply(runs[jobs$config == cfg], `[[`, "acceptance")))
  data.frame(config = cfg, t(a[a > 0]), check.names = FALSE)
}))
dir.create("data-raw/mixing-results", showWarnings = FALSE)
saveRDS(list(exact = exact, out = out, meanX = meanX, meanXz = meanXz, runs = runs,
             jobs = jobs, nIter = nIter, burnin = burnin, runtime = runtime),
        "data-raw/mixing-results/rung4b-shift-fixed.rds")
write.csv(out, "data-raw/mixing-results/rung4b-shift-fixed.csv", row.names = FALSE)
print(out, digits = 3)
print(meanX, digits = 4)
cat("meanX z (reference - shift):", round(meanXz, 2), "\n")
print(acceptance, digits = 3)
cat("runtime", format(runtime), "\n")
