# Sensitivity of the clade-size inversion: P(n | n > 0) for n < 1000 from an
# independent contour (4096 points, aliasing bound 1e-12) vs the main run.
library("EvoTempo")
outDir <- "data-raw/bm-results"
cores <- as.integer(Sys.getenv("EVOTEMPO_CORES", "13"))
Check <- function(Q, x, start, main) {
  r <- BinRates(x, 0.51, 0.5)
  p <- CladeSizeDistribution(500, Q, r$lambda, r$mu, start = start, nMax = 1000,
                             nPoint = 4096, radius = 1e-12^(1 / 4096), cores = cores)[, 1]
  ref <- readRDS(file.path(outDir, main))[1:1001, 1]
  c(maxRelDiff = max(abs(p / ref - 1)), lt1000 = sum(p[2:1000]) / (1 - p[1]),
    lt1000ref = sum(ref[2:1000]) / (1 - ref[1]))
}
bm <- BuddMannGenerator()
grid <- TempoGrid(101, c(-5, 5))
check <- rbind(
  "budd-mann" = Check(bm$Q, bm$x, 51, "pn-budd-mann.rds"),
  "central-101" = Check(TempoGenerator(grid, 0.01), grid$x, 51, "pn-central-101.rds")
)
write.csv(check, file.path(outDir, "pn-check.csv"))
print(check)
