# Rung 4b at fixed parameters: importance-sampling estimate of the joint
# tree + character likelihood from data-augmented paths (FFBS draws with
# exact/frozen weights times the Mk likelihood) vs the exact product-space
# likelihood TreeLikelihood(character =), at six (psi, lambda) points.
# Run from the package root after data-raw/rung5-validation.R 4b.
suppressMessages(pkgload::load_all(".", quiet = TRUE)); library(parallel)
r <- readRDS("data-raw/rung5-results/rung4b.rds")
grid <- expand.grid(psi = c(0.06, 0.1, 0.16), lambda = c(0.3, 0.5))
cl <- makeCluster(6); invisible(clusterCall(cl, function() suppressMessages(pkgload::load_all(".", quiet = TRUE))))
clusterExport(cl, c("r", "grid"))
out <- parLapply(cl, seq_len(nrow(grid)), function(i) {
  par <- r$truth; par$psi <- grid$psi[i]; par$lambda <- grid$lambda[i]
  exact <- TreeLikelihood(r$tt, par, r$model, character = r$character)
  set <- .AddSlices(r$tt, .Precompute(r$tt, par, r$model), 20); msg <- .Messages(r$tt, set)
  frozen <- .FrozenLikelihood(r$tt, set, msg)
  set.seed(i)
  lw <- vapply(1:2500, function(k) {
    d <- .SampleTree(r$tt, set, msg)
    .LogWeight(r$tt, set, msg, d$paths, .PathNodes(r$tt)) +
      MkLikelihood(r$tt, matrix(r$character), .CharacterLengths(r$tt, set$grid$x, d$paths, par$rate, r$model$gamma))
  }, 0)
  m <- max(lw); w <- exp(lw - m)
  c(grid[i, ], exact = exact, is = frozen + m + log(mean(w)), se = sd(w) / mean(w) / sqrt(length(w)))
})
stopCluster(cl)
res <- do.call(rbind, lapply(out, as.data.frame)); res$z <- (res$is - res$exact) / res$se; print(res, digits = 7)
saveRDS(res, "data-raw/rung5-results/rung4b-likelihood-check.rds")
