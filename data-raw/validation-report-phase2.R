# Summaries of rungs 3-6 for the Phase 2 validation report: CSV tables in
# data-raw/validation-report-phase2/ (figures are drawn from these).
# Run from the package root after the rung scripts.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
outDir <- "data-raw/validation-report-phase2"
# EVOTEMPO_GRID=default summarises the "-K26" results of the default grid (#17)
suffix <- if (identical(Sys.getenv("EVOTEMPO_GRID"), "default")) "-K26" else ""
if (nzchar(suffix)) outDir <- "data-raw/validation-report-default-grid"
Res <- function(rung, name) {
  file.path("data-raw", paste0(rung, "-results"), paste0(name, suffix, ".rds"))
}
dir.create(outDir, showWarnings = FALSE)
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)

# Quantiles of pooled chains, with standard errors from the spread of
# per-chain quantiles (chains are independent)
ChainQuantiles <- function(chains, p) {
  per <- vapply(chains, function(x) stats::quantile(x, p, names = FALSE), 0)
  c(estimate = stats::quantile(unlist(chains), p, names = FALSE),
    se = stats::sd(per) / sqrt(length(per)))
}

# Two-sample comparison of chain sets A and B for one quantity
Compare <- function(a, b, name) {
  rows <- lapply(probs, function(p) {
    qa <- ChainQuantiles(a, p)
    qb <- if (is.list(b)) ChainQuantiles(b, p) else {
      # Independent draws: bootstrap standard error
      c(estimate = stats::quantile(b, p, names = FALSE),
        se = stats::sd(replicate(200, stats::quantile(sample(b, replace = TRUE), p,
                                                       names = FALSE))))
    }
    data.frame(quantity = name, p = p, a = qa[["estimate"]], seA = qa[["se"]],
               b = qb[["estimate"]], seB = qb[["se"]],
               z = (qa[["estimate"]] - qb[["estimate"]]) / sqrt(qa[["se"]]^2 + qb[["se"]]^2))
  })
  essA <- sum(vapply(a, function(x) coda::effectiveSize(x), 0))
  thinA <- unlist(lapply(a, function(x) {
    x[seq(1, length(x), length.out = max(2, floor(coda::effectiveSize(x))))]
  }))
  bb <- if (is.list(b)) unlist(lapply(b, function(x) {
    x[seq(1, length(x), length.out = max(2, floor(coda::effectiveSize(x))))]
  })) else b
  ks <- suppressWarnings(stats::ks.test(thinA, bb))
  list(quantiles = do.call(rbind, rows),
       summary = data.frame(quantity = name, essA = essA,
                            essB = if (is.list(b)) sum(vapply(b, coda::effectiveSize, 0)) else length(b),
                            ksD = unname(ks$statistic), ksP = ks$p.value))
}

Bind <- function(results, part) do.call(rbind, lapply(results, `[[`, part))

# Batch-means standard error of the mean of one chain
BatchSE <- function(x, nBatch = 50) {
  m <- floor(length(x) / nBatch) * nBatch
  stats::sd(colMeans(matrix(x[seq_len(m)], ncol = nBatch))) / sqrt(nBatch)
}

# Rung 3 -----------------------------------------------------------------------
if (file.exists(Res("rung3", "rung3"))) {
  r3 <- readRDS(Res("rung3", "rung3"))
  write.csv(r3$collapse, file.path(outDir, "rung3-collapse.csv"), row.names = FALSE)
  if (!is.null(r3$castor)) {
    write.csv(r3$castor, file.path(outDir, "rung3-castor.csv"), row.names = FALSE)
  }
  write.csv(r3$mc, file.path(outDir, "rung3-montecarlo.csv"), row.names = FALSE)
  write.csv(as.data.frame(r3$mcSummary), file.path(outDir, "rung3-montecarlo-summary.csv"),
            row.names = FALSE)
}

# Rung 4 -----------------------------------------------------------------------
if (file.exists(Res("rung4", "rung4"))) {
  r4 <- readRDS(Res("rung4", "rung4"))
  pars4 <- names(r4$priors)
  da <- r4$chains[r4$jobs$kind == "DA"]
  marg <- r4$chains[r4$jobs$kind == "marginal"]
  cmp4 <- lapply(pars4, function(p) {
    Compare(lapply(da, function(ch) ch$samples[[p]]),
            lapply(marg, function(ch) ch$samples[[p]]), p)
  })
  write.csv(Bind(cmp4, "quantiles"), file.path(outDir, "rung4a-quantiles.csv"), row.names = FALSE)
  write.csv(Bind(cmp4, "summary"), file.path(outDir, "rung4a-summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(seq_along(da), function(i) {
    data.frame(chain = i, t(da[[i]]$acceptance))
  })), file.path(outDir, "rung4a-acceptance.csv"), row.names = FALSE)
  # Pooled draws for figures (thinned)
  Pool <- function(chains, kind, pars, every = 10) {
    do.call(rbind, lapply(seq_along(chains), function(i) {
      s <- chains[[i]]$samples
      keep <- seq(1, nrow(s), by = every)
      data.frame(kind = kind, chain = i, s[keep, pars, drop = FALSE])
    }))
  }
  write.csv(rbind(Pool(da, "augmented", pars4), Pool(marg, "marginal", pars4)),
            file.path(outDir, "rung4a-draws.csv"), row.names = FALSE)
  # 4a(ii): local path moves at fixed parameters vs weighted FFBS.  Eight
  # independent local-move chains; standard errors from between-chain spread.
  local <- readRDS(Res("rung4", "rung4-local-chains"))
  maps <- r4$maps
  w <- exp(maps$logWeight - max(maps$logWeight))
  w <- w / sum(w)
  essIS <- 1 / sum(w^2)
  Between <- function(part) {
    means <- t(vapply(local, function(ch) colMeans(ch[[part]]), numeric(ncol(local[[1]][[part]]))))
    list(mean = colMeans(means), se = apply(means, 2, stats::sd) / sqrt(nrow(means)))
  }
  occ <- Between("occupancy")
  bx <- Between("branchX")
  occIS <- colSums(w * maps$occupancy)
  occISSE <- sqrt(colSums(w * (sweep(maps$occupancy, 2, occIS))^2) / essIS)
  bxIS <- colSums(w * maps$branchX)
  bxISSE <- sqrt(colSums(w * (sweep(maps$branchX, 2, bxIS))^2) / essIS)
  x <- .Setup(r4$truth, r4$model)$grid$x
  ffbs <- rbind(
    data.frame(quantity = "occupancy", index = seq_along(occIS), x = x, chain = occ$mean,
               seChain = occ$se, ffbs = occIS, seFFBS = occISSE),
    data.frame(quantity = "branchX", index = seq_along(bxIS), x = NA, chain = bx$mean,
               seChain = bx$se, ffbs = bxIS, seFFBS = bxISSE))
  ffbs$z <- (ffbs$chain - ffbs$ffbs) / sqrt(ffbs$seChain^2 + ffbs$seFFBS^2)
  ffbs$essIS <- essIS
  write.csv(do.call(rbind, lapply(seq_along(local), function(i) {
    data.frame(chain = i, t(local[[i]]$acceptance))
  })), file.path(outDir, "rung4a-local-acceptance.csv"), row.names = FALSE)
  write.csv(ffbs, file.path(outDir, "rung4a-ffbs.csv"), row.names = FALSE)
  write.csv(r4$laplace, file.path(outDir, "rung4c-laplace.csv"), row.names = FALSE)
  p33 <- r4$prop33
  write.csv(data.frame(from = rep(seq_len(nrow(p33$mc)), ncol(p33$mc)),
                       to = rep(seq_len(ncol(p33$mc)), each = nrow(p33$mc)),
                       mc = as.vector(p33$mc), se = as.vector(p33$se),
                       exact = as.vector(p33$exact),
                       z = as.vector((p33$mc - p33$exact) / p33$se)),
            file.path(outDir, "rung4d-sum-over-paths.csv"), row.names = FALSE)
}

# Rung 4b and 5 ----------------------------------------------------------------
if (file.exists(Res("rung5", "rung4b"))) {
  r4b <- readRDS(Res("rung5", "rung4b"))
  pars4b <- names(r4b$truth)
  pars4b <- pars4b[pars4b %in% names(r4b$chains[[1]]$samples)]
  da <- r4b$chains[r4b$jobs$kind == "DA"]
  marg <- r4b$chains[r4b$jobs$kind == "marginal"]
  cmp <- lapply(pars4b, function(p) {
    Compare(lapply(da, function(ch) ch$samples[[p]]),
            lapply(marg, function(ch) ch$samples[[p]]), p)
  })
  write.csv(Bind(cmp, "quantiles"), file.path(outDir, "rung4b-quantiles.csv"), row.names = FALSE)
  write.csv(Bind(cmp, "summary"), file.path(outDir, "rung4b-summary.csv"), row.names = FALSE)
  write.csv(rbind(Pool(da, "augmented", pars4b), Pool(marg, "marginal", pars4b)),
            file.path(outDir, "rung4b-draws.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(seq_along(da), function(i) {
    data.frame(chain = i, t(da[[i]]$acceptance))
  })), file.path(outDir, "rung4b-acceptance.csv"), row.names = FALSE)
}
if (file.exists(Res("rung5", "rung5"))) {
  r5 <- readRDS(Res("rung5", "rung5"))
  continuous <- c("lambda", "mu", "psi", "theta", "s", "treeLength", "meanNodeAge",
                  "oldestFossil", "youngestFossil", "clock", "inverseClock")
  cmp <- lapply(continuous, function(q) {
    Compare(lapply(r5$chains, function(ch) ch$samples[[q]]), r5$sims[[q]], q)
  })
  write.csv(Bind(cmp, "quantiles"), file.path(outDir, "rung5-quantiles.csv"), row.names = FALSE)
  write.csv(Bind(cmp, "summary"), file.path(outDir, "rung5-summary.csv"), row.names = FALSE)
  discrete <- do.call(rbind, lapply(c("nAncestor", "cherries", "rootBin"), function(q) {
    mcmc <- unlist(lapply(r5$chains, function(ch) ch$samples[[q]]))
    sim <- r5$sims[[q]]
    lev <- sort(unique(c(mcmc, sim)))
    pm <- tabulate(match(mcmc, lev), length(lev)) / length(mcmc)
    ps <- tabulate(match(sim, lev), length(lev)) / length(sim)
    ess <- sum(vapply(r5$chains, function(ch) coda::effectiveSize(ch$samples[[q]]), 0))
    data.frame(quantity = q, level = lev, mcmc = pm, sim = ps,
               seMCMC = sqrt(pm * (1 - pm) / ess), seSim = sqrt(ps * (1 - ps) / length(sim)),
               essMCMC = ess)
  }))
  discrete$z <- (discrete$mcmc - discrete$sim) /
    sqrt(discrete$seMCMC^2 + discrete$seSim^2)
  write.csv(discrete, file.path(outDir, "rung5-discrete.csv"), row.names = FALSE)
  write.csv(rbind(
    data.frame(kind = "MCMC", do.call(rbind, lapply(seq_along(r5$chains), function(i) {
      s <- r5$chains[[i]]$samples
      s[seq(1, nrow(s), by = 10), c(continuous, "nAncestor", "cherries")]
    }))),
    data.frame(kind = "simulation", r5$sims[, c(continuous, "nAncestor", "cherries")])),
    file.path(outDir, "rung5-draws.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(seq_along(r5$chains), function(i) {
    data.frame(chain = i, t(r5$chains[[i]]$acceptance))
  })), file.path(outDir, "rung5-acceptance.csv"), row.names = FALSE)
}

if (file.exists(Res("rung5", "rung4b-likelihood-check"))) {
  write.csv(readRDS(Res("rung5", "rung4b-likelihood-check")),
            file.path(outDir, "rung4b-likelihood-check.csv"), row.names = FALSE)
}
if (file.exists(Res("rung5", "rung5-fixed-parameters"))) {
  # Fixed-parameter diagnostic: list of 8 simulation matrices then 8 chains
  fx <- readRDS(Res("rung5", "rung5-fixed-parameters"))
  sim <- as.data.frame(do.call(rbind, fx[1:8]))
  ch <- fx[9:16]
  write.csv(do.call(rbind, lapply(c("treeLength", "meanNodeAge", "oldestFossil",
                                     "youngestFossil", "clock", "inverseClock",
                                     "nAncestor", "cherries"), function(q) {
    per <- vapply(ch, function(s) mean(s[[q]], na.rm = TRUE), 0)
    simQ <- sim[[q]][!is.na(sim[[q]])]
    data.frame(quantity = q, sim = mean(simQ), seSim = stats::sd(simQ) / sqrt(length(simQ)),
               mcmc = mean(per), seMCMC = stats::sd(per) / sqrt(length(per)),
               z = (mean(per) - mean(simQ)) /
                 sqrt(stats::var(per) / length(per) + stats::var(simQ) / length(simQ)))
  })), file.path(outDir, "rung5-fixed-parameters.csv"), row.names = FALSE)
}

# Rung 6 -----------------------------------------------------------------------
if (file.exists(Res("rung6", "rung6"))) {
  r6 <- readRDS(Res("rung6", "rung6"))
  L <- r6$nDraw
  ranks <- do.call(rbind, lapply(r6$results, function(r) {
    data.frame(cell = r$job$cell, rep = r$job$rep, quantity = names(r$ranks),
               rank = unname(r$ranks), nTip = r$nTip)
  }))
  ranks$beta <- r6$cells$beta[ranks$cell]
  ranks$gamma <- r6$cells$gamma[ranks$cell]
  write.csv(ranks, file.path(outDir, "rung6-ranks.csv"), row.names = FALSE)
  tests <- do.call(rbind, lapply(split(ranks, list(ranks$cell, ranks$quantity)), function(d) {
    bins <- cut(d$rank, breaks = seq(-0.5, L + 0.5, length.out = 6))
    counts <- table(bins)
    expected <- length(d$rank) / 5
    u <- (d$rank + stats::runif(nrow(d))) / (L + 1)
    data.frame(cell = d$cell[1], beta = d$beta[1], gamma = d$gamma[1],
               quantity = d$quantity[1], n = nrow(d),
               chi2 = sum((counts - expected)^2 / expected),
               chi2P = stats::pchisq(sum((counts - expected)^2 / expected), 4,
                                     lower.tail = FALSE),
               ksP = suppressWarnings(stats::ks.test(u, "punif")$p.value),
               meanRank = mean(d$rank) / L)
  }))
  write.csv(tests, file.path(outDir, "rung6-tests.csv"), row.names = FALSE)
  ess <- do.call(rbind, lapply(r6$results, function(r) {
    data.frame(cell = r$job$cell, rep = r$job$rep, t(r$ess), runtime = r$runtime,
               nTip = r$nTip, nAncestor = r$nAncestor, tries = r$tries)
  }))
  write.csv(ess, file.path(outDir, "rung6-ess.csv"), row.names = FALSE)
}
