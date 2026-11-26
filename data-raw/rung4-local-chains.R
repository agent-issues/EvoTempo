# Rung 4a(ii): eight independent chains of local path moves at fixed parameters
# (between-chain standard errors), vs the importance-weighted FFBS maps of
# rung4-validation.R.  Run from the package root after rung4-validation.R.
# Optional arguments: iterations per chain (default 3000) and seed base
# (default 500: chain i uses seed base + i; a base other than 500 writes
# rung4-local-chains-seed<base>.rds).  The #116 runs used 20000 with bases
# 500 and 600 (data-raw/tempo-location-results/README.md).
suppressMessages(pkgload::load_all(".", quiet = TRUE)); library(parallel)
nIter <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[[1]]) else 3000L
seedBase <- if (length(commandArgs(TRUE)) > 1) as.integer(commandArgs(TRUE)[[2]]) else 500L
suffix <- if (identical(Sys.getenv("EVOTEMPO_GRID"), "default")) "-K26" else ""
r4 <- readRDS(paste0("data-raw/rung4-results/rung4", suffix, ".rds"))
cl <- makeCluster(8)
invisible(clusterCall(cl, function() suppressMessages(pkgload::load_all(".", quiet = TRUE))))
clusterExport(cl, c("r4", "nIter", "seedBase"), envir = environment())
local <- parLapply(cl, 1:8, function(i) {
  set.seed(seedBase + i)
  TempoMCMC(r4$tt, r4$model, r4$truth, list(), nIter = nIter, burnin = 100, nSlice = 20,
            moves = c(node = 2, tip = 2, branch = 2))[c("samples", "occupancy", "branchX", "acceptance")]
})
stopCluster(cl)
saveRDS(local, paste0("data-raw/rung4-results/rung4-local-chains",
                      if (seedBase != 500L) paste0("-seed", seedBase), suffix, ".rds"))
