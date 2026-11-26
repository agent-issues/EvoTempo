# Rung 4a(ii) rerun (review of #123, B1): local path moves vs importance-weighted
# frozen FFBS at fixed parameters, with a whole-tree-proposal diagnostic arm.
# The design, seeds and decision rule are fixed in dev/plans/rung4a-local-rerun.md.
# Usage (from the package root, EvoTempo installed):
#   Rscript data-raw/rung4a-local-rerun.R <mean|bm> [nChain = 32] [nIter = 50000]
# EVOTEMPO_CORES (default 16) sets the cluster size. Writes, in
# EVOTEMPO_RERUN_OUT (default data-raw/tempo-location-results/rung4a-rerun):
#   rerun-<location>.rds          block means per chain, reference sums, provenance
#   rerun-<location>-summary.csv  z per quantity and arm, overall offsets
#   rerun-<location>-verdict.txt  the pre-registered decision
library("EvoTempo")
library("parallel")
args <- commandArgs(trailingOnly = TRUE)
location <- match.arg(args[[1]], c("mean", "bm"))
nChain <- if (length(args) > 1) as.integer(args[[2]]) else 32L
nIter <- if (length(args) > 2) as.integer(args[[3]]) else 50000L
nCore <- as.integer(Sys.getenv("EVOTEMPO_CORES", "16"))
outDir <- Sys.getenv("EVOTEMPO_RERUN_OUT", "data-raw/tempo-location-results/rung4a-rerun")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
nBlock <- 50L
nRefJob <- 16L
nRefDraw <- as.integer(Sys.getenv("EVOTEMPO_RERUN_REFDRAW", "4000")) # smoke tests only
nDiag <- 16L
nDiagIter <- as.integer(Sys.getenv("EVOTEMPO_RERUN_DIAGITER", "20000")) # smoke tests only

fixture <- readRDS("data-raw/rung4-results/rung4a-rerun-fixture.rds")
tt <- fixture$tt
truth <- fixture$truth
model <- DCETModel(rho = 0.5, location = location)
provenance <- list(
  location = model$location, scheme = model$scheme, nBin = model$nBin, zRange = model$zRange,
  package = as.character(utils::packageVersion("EvoTempo")),
  sha = Sys.getenv("EVOTEMPO_SHA", NA_character_),
  script = unname(tools::md5sum("data-raw/rung4a-local-rerun.R")),
  fixture = unname(tools::md5sum("data-raw/rung4-results/rung4a-rerun-fixture.rds")),
  nChain = nChain, nIter = nIter, started = format(Sys.time())
)
message(paste(names(provenance), unlist(lapply(provenance, paste, collapse = ",")), sep = " = ",
              collapse = "\n"))

BlockMeans <- function(m) {
  block <- cut(seq_len(nrow(m)), nBlock, labels = FALSE)
  rowsum(m, block) / as.vector(table(block))
}
jobs <- rbind(
  data.frame(arm = "local", id = seq_len(nChain), seed = 41000L + seq_len(nChain)),
  data.frame(arm = "reference", id = seq_len(nRefJob), seed = 42000L + seq_len(nRefJob)),
  data.frame(arm = "diagnostic", id = seq_len(nDiag), seed = 43000L + seq_len(nDiag))
)
cl <- makeCluster(nCore)
invisible(clusterEvalQ(cl, suppressMessages(library("EvoTempo"))))
clusterExport(cl, c("tt", "truth", "model", "nIter", "nDiagIter", "nRefDraw", "BlockMeans",
                    "nBlock", "jobs"), envir = environment())
t0 <- Sys.time()
res <- parLapplyLB(cl, seq_len(nrow(jobs)), function(j) {
  job <- jobs[j, ]
  set.seed(job$seed)
  if (job$arm == "reference") {
    maps <- SampleTempoPaths(tt, truth, model, n = nRefDraw, nSlice = 20)
    return(list(logWeight = maps$logWeight, branchX = maps$branchX,
                occupancy = maps$occupancy))
  }
  moves <- if (job$arm == "local") c(node = 2, tip = 2, branch = 2) else c(tree = 1)
  run <- TempoMCMC(tt, model, truth, list(), nIter = if (job$arm == "local") nIter else nDiagIter,
                   burnin = 100, nSlice = 20, moves = moves)
  list(branchX = BlockMeans(run$branchX), occupancy = BlockMeans(run$occupancy),
       acceptance = run$acceptance)
})
stopCluster(cl)
provenance$runtime <- format(difftime(Sys.time(), t0, units = "mins"))

# Importance-weighted reference over all draws
ref <- res[jobs$arm == "reference"]
logW <- unlist(lapply(ref, `[[`, "logWeight"))
w <- exp(logW - max(logW))
w <- w / sum(w)
essIS <- 1 / sum(w^2)
Reference <- function(part) {
  m <- do.call(rbind, lapply(ref, `[[`, part))
  mu <- colSums(w * m)
  list(mean = mu, se = sqrt(colSums(w * sweep(m, 2, mu)^2) / essIS))
}
reference <- list(branchX = Reference("branchX"), occupancy = Reference("occupancy"))

Compare <- function(arm) {
  chains <- res[jobs$arm == arm]
  rows <- lapply(c("branchX", "occupancy"), function(part) {
    means <- t(vapply(chains, function(ch) colMeans(ch[[part]]), numeric(ncol(chains[[1]][[part]]))))
    se <- apply(means, 2, stats::sd) / sqrt(nrow(means))
    r <- reference[[part]]
    data.frame(arm = arm, quantity = part, index = seq_along(r$mean), chain = colMeans(means),
               seChain = se, reference = r$mean, seReference = r$se,
               z = (colMeans(means) - r$mean) / sqrt(se^2 + r$se^2),
               tested = if (part == "branchX") TRUE else r$mean >= 1e-3)
  })
  # Overall offset: mean of the 16 branch means, per chain
  bx <- vapply(chains, function(ch) mean(colMeans(ch$branchX)), 0)
  refAvg <- mean(reference$branchX$mean)
  w16 <- rep(1 / length(reference$branchX$mean), length(reference$branchX$mean))
  refAvgSe <- sqrt(sum(w * (as.vector(do.call(rbind, lapply(ref, `[[`, "branchX")) %*% w16) -
                              refAvg)^2) / essIS)
  offset <- data.frame(arm = arm, quantity = "branchAverage", index = NA, chain = mean(bx),
                       seChain = stats::sd(bx) / sqrt(length(bx)), reference = refAvg,
                       seReference = refAvgSe,
                       z = (mean(bx) - refAvg) / sqrt(stats::var(bx) / length(bx) + refAvgSe^2),
                       tested = TRUE)
  rbind(do.call(rbind, rows), offset)
}
summary <- rbind(Compare("local"), Compare("diagnostic"))
summary$location <- location
summary$essIS <- essIS
Verdict <- function(arm) {
  s <- summary[summary$arm == arm & summary$tested, ]
  maxZ <- max(abs(s$z[s$quantity != "branchAverage"]))
  offZ <- abs(s$z[s$quantity == "branchAverage"])
  sprintf("%s arm, location %s: max |z| %.2f over %d quantities (< 3.3), offset |z| %.2f (< 3): %s",
          arm, location, maxZ, sum(s$quantity != "branchAverage"), offZ,
          if (maxZ < 3.3 && offZ < 3) "PASS" else "FAIL")
}
verdict <- c(Verdict("local"), Verdict("diagnostic"),
             sprintf("reference: %d draws, importance ESS %.1f", length(w), essIS),
             sprintf("branches with chain mean above reference (local arm): %d of %d",
                     sum(summary$arm == "local" & summary$quantity == "branchX" & summary$z > 0),
                     sum(summary$arm == "local" & summary$quantity == "branchX")))
writeLines(verdict)
writeLines(verdict, file.path(outDir, paste0("rerun-", location, "-verdict.txt")))
utils::write.csv(summary, file.path(outDir, paste0("rerun-", location, "-summary.csv")),
                 row.names = FALSE)
saveRDS(list(provenance = provenance, jobs = jobs,
             chains = res[jobs$arm != "reference"],
             reference = reference, essIS = essIS),
        file.path(outDir, paste0("rerun-", location, ".rds")))
