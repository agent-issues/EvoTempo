# Rung 4a(ii) as reported in #123 under location = "mean": 16 chains of local
# path moves (data-raw/rung4-local-chains.R 20000 500 and 20000 600, at the
# default grid, on the "mean" rung-4 fixture), pooled over every iteration,
# against the weighted-FFBS reference of that fixture (rungs/rung4a-ffbs.csv,
# columns ffbs and seFFBS).  The chains are stored as 50 block means each
# (rungs/rung4a-local-blockmeans.rds); equal blocks, so their mean is the chain
# mean.  Writes rungs/rung4a-ffbs-pooled16.csv.
# Run from the package root: Rscript data-raw/tempo-location-results/rung4a-pool.R
dir <- "data-raw/tempo-location-results/rungs"
chains <- readRDS(file.path(dir, "rung4a-local-blockmeans.rds"))
ref <- utils::read.csv(file.path(dir, "rung4a-ffbs.csv"))
out <- do.call(rbind, lapply(c("occupancy", "branchX"), function(q) {
  r <- ref[ref$quantity == q, ]
  m <- t(vapply(chains, function(ch) colMeans(ch[[q]]), numeric(nrow(r))))
  r$chain <- colMeans(m)
  r$seChain <- apply(m, 2, stats::sd) / sqrt(nrow(m))
  r$z <- (r$chain - r$ffbs) / sqrt(r$seChain^2 + r$seFFBS^2)
  r
}))
utils::write.csv(out, file.path(dir, "rung4a-ffbs-pooled16.csv"), row.names = FALSE)
b <- out[out$quantity == "branchX", ]
cat(sprintf("branch means: %d of %d above FFBS, z %.2f to %.2f\n", sum(b$z > 0), nrow(b),
            min(b$z), max(b$z)))
