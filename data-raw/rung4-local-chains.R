# Rung 4a(ii): eight independent chains of local path moves at fixed parameters
# (between-chain standard errors), vs the importance-weighted FFBS maps of
# rung4-validation.R.  Run from the package root after rung4-validation.R.
suppressMessages(pkgload::load_all(".", quiet = TRUE)); library(parallel)
suffix <- if (identical(Sys.getenv("EVOTEMPO_GRID"), "default")) "-K26" else ""
r4 <- readRDS(paste0("data-raw/rung4-results/rung4", suffix, ".rds"))
cl <- makeCluster(8)
invisible(clusterCall(cl, function() suppressMessages(pkgload::load_all(".", quiet = TRUE))))
clusterExport(cl, "r4", envir = environment())
local <- parLapply(cl, 1:8, function(i) {
  set.seed(500 + i)
  TempoMCMC(r4$tt, r4$model, r4$truth, list(), nIter = 3000, burnin = 100, nSlice = 20,
            moves = c(node = 2, tip = 2, branch = 2))[c("samples", "occupancy", "branchX", "acceptance")]
})
stopCluster(cl)
saveRDS(local, paste0("data-raw/rung4-results/rung4-local-chains", suffix, ".rds"))
