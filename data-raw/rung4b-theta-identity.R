# Rung 4b localisation (#20): importance-sampling identity across theta, with
# (argument "character") or without the character: FFBS maps weighted by
# exact/frozen (times the Mk likelihood) vs TreeLikelihood(character = ).
suppressMessages(pkgload::load_all(".", quiet = TRUE)); library(parallel)
r <- readRDS("data-raw/rung5-results/rung4b.rds")
thetas <- c(0.1, 0.2, 0.35, 0.5, 0.7)
char <- if ("character" %in% commandArgs(trailingOnly = TRUE)) r$character else NULL
cl <- makeCluster(5); invisible(clusterCall(cl, function() suppressMessages(pkgload::load_all(".", quiet = TRUE))))
clusterExport(cl, c("r", "thetas", "char"))
out <- parLapply(cl, seq_along(thetas), function(i) {
  par <- r$truth; par$theta <- thetas[i]
  exact <- TreeLikelihood(r$tt, par, r$model, character = char)
  set <- .AddSlices(r$tt, .Precompute(r$tt, par, r$model), 20); msg <- .Messages(r$tt, set)
  frozen <- .FrozenLikelihood(r$tt, set, msg)
  set.seed(i)
  lw <- vapply(1:3000, function(k) {
    d <- .SampleTree(r$tt, set, msg)
    mk <- if (is.null(char)) 0 else MkLikelihood(r$tt, matrix(char),
      .CharacterLengths(r$tt, set$grid$x, d$paths, par$rate, r$model$gamma))
    c(.LogWeight(r$tt, set, msg, d$paths, .PathNodes(r$tt)) + mk,
      .LogAugmented(r$tt, set, d$paths, d$bin))
  }, c(0, 0))
  m <- max(lw[1, ]); w <- exp(lw[1, ] - m)
  # Consistency: log p_frozen(paths) = .LogAugmented - .LogWeight should equal
  # log q(paths) + frozen, i.e. the proposal density times the frozen likelihood
  data.frame(theta = thetas[i], exact = exact, is = frozen + m + log(mean(w)),
             se = sd(w) / mean(w) / sqrt(length(w)), essFrac = mean(w)^2 / mean(w^2))
})
stopCluster(cl)
res <- do.call(rbind, out); res$z <- (res$is - res$exact) / res$se
print(res, digits = 6)
saveRDS(res, paste0("data-raw/rung5-results/rung4b-theta-identity",
                   if (is.null(char)) "" else "-character", ".rds"))
