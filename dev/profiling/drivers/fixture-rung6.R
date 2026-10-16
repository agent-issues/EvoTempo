# Shared fixture for the TempoMCMC() drivers: one rung-6-style replicate
# (beta = gamma = 1, K = 16, z in [-3, 3], crown 8, 20-40 tips, 30 binary
# characters).  Defines model, priors, truth, tt (TempoTree) and chars.
.Setup <- EvoTempo:::.Setup
.Sample1 <- EvoTempo:::.Sample1
meanlog <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.05), theta = log(0.2),
             s = log(0.6), rate = log(0.1))
sdlog <- c(lambda = 0.3, mu = 0.3, psi = 0.3, theta = 0.5, s = 0.25, rate = 0.4)
priors <- TempoPriors(meanlog, sdlog)
crown <- 8
model <- DCETModel(nBin = 16, zRange = c(-3, 3), rho = 0.5, beta = 1, gamma = 1)

Join <- function(h1, h2, x0) {
  n1 <- nrow(h1$edges)
  e1 <- h1$edges
  e2 <- h2$edges
  e1$parent <- ifelse(is.na(e1$parent), 1L, e1$parent + 1L)
  e2$parent <- ifelse(is.na(e2$parent), 1L, e2$parent + 1L + n1)
  f1 <- h1$fossils
  f2 <- h2$fossils
  f1$edge <- f1$edge + 1L
  f2$edge <- f2$edge + 1L + n1
  structure(list(edges = rbind(data.frame(parent = NA_integer_, start = 0, end = 0,
                                          type = "speciation", sampled = NA), e1, e2),
                 paths = c(list(list(t = 0, x = x0)), h1$paths, h2$paths),
                 fossils = rbind(f1, f2), time = h1$time),
            class = "TempoHistory", kind = "discrete")
}

set.seed(6001)
repeat {
  truth <- as.list(setNames(stats::rlnorm(6, meanlog, sdlog), names(meanlog)))
  set <- .Setup(truth, model)
  repeat {
    k <- .Sample1(set$prior)
    h <- tryCatch(lapply(1:2, function(i) {
      SimulateDCET(crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi,
                   start = k, rho = model$rho, maxEdges = 5000)
    }), error = function(e) NULL)
    if (is.null(h) ||
        all(vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0) > 0)) {
      break
    }
  }
  if (is.null(h)) next
  history <- Join(h[[1]], h[[2]], set$grid$x[k])
  tree <- ReconstructedTree(history, gamma = 1, rate = truth$rate)
  nTip <- sum(attr(tree, "tipType") != "ancestor")
  if (nTip >= 20 && nTip <= 40) break
}
chars <- SimulateCharacters(tree, 30L, 2)
tree$root.edge <- NULL
tt <- TempoTree(tree, present = crown)
cat("Tips:", tt$nTip, " fossils:", sum(tt$tipType == "fossil"),
    " ancestors:", sum(tt$tipType == "ancestor"), "\n")

