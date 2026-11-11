# Check that the coda-free diagnostics in summarise.R reproduce coda.
# Run where coda is installed: Rscript check-diagnostics.R
src <- readLines("summarise.R")
eval(parse(text = src[grep("^ESS <- function", src):(grep("^logFiles <- ", src) - 1)]))
library(coda)
set.seed(1)
for (k in 1:5) {
  x <- as.numeric(arima.sim(list(ar = 0.15 * k), n = 1500))
  y <- as.numeric(arima.sim(list(ar = 0.15 * k), n = 1500)) + 0.05 * k
  z <- as.numeric(arima.sim(list(ar = 0.15 * k), n = 1500)) - 0.05 * k
  stopifnot(
    all.equal(ESS(x), unname(effectiveSize(x))),
    all.equal(PSRF(list(x, y, z)),
              gelman.diag(mcmc.list(mcmc(x), mcmc(y), mcmc(z)), autoburnin = FALSE)$psrf[[1, 1]]),
    all.equal(HPD(x), unname(HPDinterval(mcmc(x))[1, ]))
  )
}
cat("ESS, PSRF and HPD agree with coda", format(packageVersion("coda")), "\n")
