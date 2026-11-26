suppressMessages(library(EvoTempo))
t <- as.numeric(commandArgs(TRUE)[1])
states <- unlist(lapply(sort(Sys.glob("D:/claude-out/evotempo/v0/pilot-g/seg*_states.rds")), readRDS), recursive = FALSE)
first <- states[[1]]; last <- states[[length(states)]]; tree <- last$tree
pars <- c("lambda", "mu", "psi", "theta", "s")
p0 <- unlist(first$par[pars]); p1 <- unlist(last$par[pars])
p <- as.list(exp((1 - t) * log(p0) + t * log(p1)))
d <- exp((1 - t) * log(p0[["lambda"]] - p0[["mu"]]) + t * log(p1[["lambda"]] - p1[["mu"]])); p$mu <- p$lambda - d
s <- p$s; pm <- p; for (q in c("lambda", "mu", "theta")) pm[[q]] <- p[[q]] * exp(-s^2 / 2)
M <- function(K, z, loc, scheme = NULL) DCETModel(nBin = K, zRange = z, beta = 1, gamma = 0, rho = 5e-4, condition = "extant", location = loc, scheme = scheme)
L <- function(par, m) { t0 <- proc.time()[[3]]; v <- tryCatch(TreeLikelihood(tree, par, m), error = function(e) NA); c(v, proc.time()[[3]] - t0) }
r <- rbind(
  mean_K26_5   = L(pm, M(26, c(-5, 5), "mean")),
  mean_K36_7   = L(pm, M(36, c(-7, 7), "mean")),
  mean_K51_5   = L(pm, M(51, c(-5, 5), "mean")),
  mean_K31_m5p7 = L(pm, M(31, c(-5, 7), "mean")),
  bm_K26_5_central = L(p, M(26, c(-5, 5), "bm")),
  bm_K51_5_central = L(p, M(51, c(-5, 5), "bm")),
  bm_K51_5_sqra = L(p, M(51, c(-5, 5), "bm", "sqra")))
cat("t =", t, "s =", s, "\n"); print(r, digits = 9)
