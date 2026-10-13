# Rung 3: marginal dCET tree likelihood.
#
# RevBayes was not runnable on the development machine (no `rb` executable in
# the local checkout), so the dnCDBDP comparison of spec section 10 could not
# be made.  In its place:
#   (a) beta = 0 collapse to the closed-form FBD likelihood, for random dense
#       generators, every conditioning, removal and start option;
#   (b) castor's MuSSE likelihood (independent C++ code), extant-only and with
#       Poissonian-sampled fossil tips (castor removes sampled lineages), and
#       castor's crown conditioning;
#   (c) brute-force Monte Carlo with the exact dCET simulator on tiny trees,
#       with fossil tips (no removal) and sampled ancestors: probabilities of
#       tree configurations, binned by their one free age.
# Run from the package root: Rscript data-raw/rung3-validation.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(ape)
library(parallel)
outDir <- "data-raw/rung3-results"
dir.create(outDir, showWarnings = FALSE)

# (a) beta = 0 collapse ------------------------------------------------------
tr <- read.tree(text = "((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7);")
tr$root.edge <- 0.8
tt <- TempoTree(tr)
set.seed(1)
collapse <- NULL
for (cond in c("none", "extant", "survival")) for (rem in c(0, 0.4)) {
  for (st in c("root", "origin")) {
    K <- 7
    Q <- matrix(runif(K * K), K)
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    prior <- runif(K)
    a <- .PruneTree(tt, Q, 0.5, 0.2, 0.1, 0.6, rem, prior / sum(prior), cond,
                    FALSE, st)
    b <- FBDTreeLikelihood(tt, 0.5, 0.2, 0.1, 0.6, rem, cond, st)
    collapse <- rbind(collapse, data.frame(condition = cond, removal = rem,
                                           start = st, dcet = a, fbd = b))
  }
}
collapse$diff <- collapse$dcet - collapse$fbd

# (b) castor MuSSE ------------------------------------------------------------
CastorLL <- function(tree, Q, lam, mu, rho, prior, cond = "none", psi = 0,
                     fossilTips = integer(0)) {
  n <- length(tree$tip.label)
  K <- nrow(Q)
  dist <- castor::get_all_distances_to_root(tree)
  ages <- max(dist[1:n]) - dist
  tp <- matrix(rho, n, K)
  tp[fossilTips, ] <- rep(psi, each = length(fossilTips))
  # fit_musse() fails when every parameter is fixed (`fits` undefined), so
  # call the likelihood routine that it wraps
  castor:::get_MuSSE_loglikelihood_CPP(
    Ntips = n, Nnodes = tree$Nnode, Nedges = nrow(tree$edge), Nstates = K,
    oldest_age = max(ages), tree_edge = as.vector(t(tree$edge)) - 1,
    clade_ages = ages, transition_rates = as.vector(t(Q)),
    speciation_rates = lam, extinction_rates = mu,
    sampling_rates = rep_len(psi, K), initial_D_per_tip = as.vector(t(tp)),
    initial_E_per_state = rep(1 - rho, K), root_prior_type = "custom",
    root_prior = prior, root_conditioning = cond,
    include_ancestral_likelihoods = FALSE, include_warnings = FALSE,
    max_condition_number = 1e4, relative_ODE_step = 1e-5, E_value_step = 1e-6,
    D_temporal_resolution = 10000, runtime_out_seconds = 0)$loglikelihood
}
set.seed(3)
ultra <- rcoal(12)
ultra$edge.length <- ultra$edge.length * 3
set.seed(11)
nonUltra <- rtree(10)
nonUltra$edge.length <- nonUltra$edge.length * 2
ttNonUltra <- TempoTree(nonUltra)
fossils <- which(ttNonUltra$tipType == "fossil")
set.seed(5)
points <- data.frame(lambda = runif(20, 0.3, 1), mu = runif(20, 0, 0.3),
                     theta = runif(20, 0.05, 1), s = runif(20, 0.5, 1.5))
nBin <- 5
castorCmp <- do.call(rbind, lapply(seq_len(nrow(points)), function(i) {
  p <- points[i, ]
  grid <- TempoGrid(nBin, c(-1.5, 1.5), p$s)
  Q <- TempoGenerator(grid, p$theta)
  r <- BinRates(grid, p$lambda, p$mu, psi = 0.15)
  prior <- TempoEquilibrium(Q, r$lambda, r$mu)$p
  w <- prior * r$lambda # castor includes lambda at the root bifurcation
  data.frame(
    point = i,
    extant = CastorLL(ultra, Q, r$lambda, r$mu, 0.7, prior) -
      (.PruneTree(TempoTree(ultra), Q, r$lambda, r$mu, 0, 0.7, 0, w / sum(w),
                  "none", FALSE, "root") + log(sum(w))),
    crown = CastorLL(ultra, Q, r$lambda, r$mu, 0.7, prior, "crown") -
      .PruneTree(TempoTree(ultra), Q, r$lambda, r$mu, 0, 0.7, 0, prior,
                 "extant", TRUE, "root"),
    fossil = CastorLL(nonUltra, Q, r$lambda, r$mu, 0.7, prior, psi = r$psi,
                      fossilTips = fossils) -
      (.PruneTree(ttNonUltra, Q, r$lambda, r$mu, r$psi, 0.7, 1, w / sum(w),
                  "none", FALSE, "root") + log(sum(w)))
  )
}))

# (c) Monte Carlo on tiny trees -----------------------------------------------
mcPar <- list(lambda = 0.6, mu = 0.3, psi = 0.25, theta = 0.5, s = 1)
mcModel <- DCETModel(nBin = 9, zRange = c(-2, 2), rho = 0.7, condition = "none")
crown <- 3
set <- .Setup(mcPar, mcModel)
breaks <- seq(0, crown, length.out = 7)

# Classify one daughter's history
Daughter <- function(h) {
  sampled <- which(h$edges$sampled %in% TRUE)
  nFossil <- nrow(h$fossils)
  if (nFossil == 0 && length(sampled) == 1) {
    return(list(kind = "tip"))
  }
  if (nFossil == 1 && length(sampled) == 0) {
    return(list(kind = "fossil", age = crown - h$fossils$time))
  }
  Ancestors <- function(e) {
    out <- e
    while (!is.na(h$edges$parent[e])) {
      e <- h$edges$parent[e]
      out <- c(out, e)
    }
    out
  }
  if (nFossil == 1 && length(sampled) == 1 &&
      h$fossils$edge %in% Ancestors(sampled)) {
    return(list(kind = "ancestor", age = crown - h$fossils$time))
  }
  if (nFossil == 0 && length(sampled) == 2) {
    mrca <- max(intersect(Ancestors(sampled[1]), Ancestors(sampled[2])))
    return(list(kind = "pair", age = crown - h$edges$end[mrca]))
  }
  list(kind = "other")
}

Batch <- function(seed, n) {
  set.seed(seed)
  starts <- sample.int(length(set$prior), n, replace = TRUE, prob = set$prior)
  out <- matrix(NA_real_, n, 2)
  kinds <- character(n)
  for (i in seq_len(n)) {
    d <- lapply(1:2, function(j) {
      Daughter(SimulateDCET(crown, set$grid$x, set$Q, set$lambda, set$mu,
                            set$psi, start = starts[i], rho = 0.7))
    })
    k <- vapply(d, `[[`, "", "kind")
    if (all(k == "tip")) {
      kinds[i] <- "A"
    } else if (sum(k == "tip") == 1) {
      other <- d[[which(k != "tip")]]
      kinds[i] <- c(fossil = "B", ancestor = "C", pair = "D",
                    other = "other")[[other$kind]]
      if (kinds[i] != "other") out[i, 1] <- other$age
    } else {
      kinds[i] <- "other"
    }
  }
  data.frame(kind = kinds, age = out[, 1])
}
nPerBatch <- 25000
nBatch <- 112
cl <- makeCluster(14)
invisible(clusterCall(cl, function() {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
}))
clusterExport(cl, c("Daughter", "set", "crown"))
mc <- do.call(rbind, parLapply(cl, seq_len(nBatch), Batch, n = nPerBatch))
stopCluster(cl)
nSim <- nPerBatch * nBatch

# Expected probabilities from TreeLikelihood
Lik <- function(newick, types) {
  exp(TreeLikelihood(TempoTree(read.tree(text = newick), tipType = types),
                     mcPar, mcModel))
}
Shape <- function(kind, a) {
  switch(kind,
         B = Lik(sprintf("(x:%.12f,f:%.12f);", crown, crown - a),
                 c("extant", "fossil")),
         C = Lik(sprintf("(x:%.12f,(y:%.12f,f:0):%.12f);", crown, a, crown - a),
                 c("extant", "extant", "ancestor")),
         D = Lik(sprintf("(x:%.12f,(y:%.12f,z:%.12f):%.12f);", crown, a, a,
                         crown - a), c("extant", "extant", "extant")))
}
gl <- local({ # 8-point Gauss-Legendre nodes and weights
  n <- 8
  j <- seq_len(n - 1)
  b <- j / sqrt(4 * j^2 - 1)
  e <- eigen(diag(0, n) + rbind(cbind(0, diag(b, n - 1)), 0) +
               t(rbind(cbind(0, diag(b, n - 1)), 0)), symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
})
expected <- data.frame(kind = "A", lo = NA, hi = NA,
                       p = Lik(sprintf("(x:%.12f,y:%.12f);", crown, crown),
                               c("extant", "extant")))
for (kind in c("B", "C", "D")) {
  for (b in seq_len(length(breaks) - 1)) {
    lo <- breaks[b]
    hi <- breaks[b + 1]
    a <- (hi - lo) / 2 * gl$x + (hi + lo) / 2
    p <- 2 * (hi - lo) / 2 * sum(gl$w * vapply(a, Shape, 0, kind = kind))
    expected <- rbind(expected, data.frame(kind = kind, lo = lo, hi = hi, p = p))
  }
}
expected$observed <- vapply(seq_len(nrow(expected)), function(i) {
  e <- expected[i, ]
  hit <- mc$kind == e$kind
  if (!is.na(e$lo)) hit <- hit & mc$age >= e$lo & mc$age < e$hi
  sum(hit)
}, 0)
expected$expectedCount <- expected$p * nSim
expected$z <- (expected$observed - expected$expectedCount) /
  sqrt(expected$expectedCount * (1 - expected$p))
chi2 <- sum(expected$z^2)
mcSummary <- list(nSim = nSim, chi2 = chi2, df = nrow(expected),
                  p = pchisq(chi2, nrow(expected), lower.tail = FALSE))

saveRDS(list(collapse = collapse, castor = castorCmp, mc = expected,
             mcSummary = mcSummary, mcPar = mcPar, mcModel = mcModel),
        file.path(outDir, "rung3.rds"))
print(collapse, digits = 10)
print(sapply(castorCmp[, -1], function(x) c(spread = diff(range(x)),
                                            mean = mean(x))))
print(expected, digits = 6)
print(mcSummary)
