# Clade-size distribution at T = 500 from x0 = 0 (Budd & Mann Fig 1A, 1B) by
# damped-contour inversion: 2^17 points (half solved, by conjugate symmetry),
# aliasing bounded at 1e-8. B&M's evalQ used 2^21 undamped points.
library("EvoTempo")
outDir <- "data-raw/bm-results"
nMax <- 50000
nPoint <- 2^17
radius <- 1e-8^(1 / nPoint)
cores <- as.integer(Sys.getenv("EVOTEMPO_CORES", "14"))
bm <- BuddMannGenerator()
rates <- BinRates(bm$x, 0.51, 0.5)
t0 <- Sys.time()
if (!file.exists(file.path(outDir, "pn-budd-mann.rds"))) {
pBM <- CladeSizeDistribution(500, bm$Q, rates$lambda, rates$mu, start = 51,
                             nMax = nMax, nPoint = nPoint, radius = radius,
                             cores = cores)
saveRDS(pBM, file.path(outDir, "pn-budd-mann.rds"))
message("B&M grid: ", format(Sys.time() - t0))
}
# Converged reference on z in +-5, h = 0.1 (mean clade size within 0.06% of
# K = 401 on +-10; see convergence.csv); K = 161 on +-8 costs ~15x more.
grid <- TempoGrid(101, c(-5, 5))
rates <- BinRates(grid, 0.51, 0.5)
t0 <- Sys.time()
pC <- CladeSizeDistribution(500, TempoGenerator(grid, 0.01), rates$lambda, rates$mu,
                            start = 51, nMax = nMax, nPoint = nPoint,
                            radius = radius, cores = cores)
saveRDS(pC, file.path(outDir, "pn-central-101.rds"))
message("central K = 101: ", format(Sys.time() - t0))
