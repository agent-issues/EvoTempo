# Summaries of the mixing study (issue #19): chains, probes and oracles.
#   Rscript dev/mixing/summary.R [resultsDir]
# Writes chains.csv, probes.csv and oracles.csv into the results folder.
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[[1]] else Sys.getenv("EVOTEMPO_MIXING_OUT", "dev/mixing/results")
RbindFill <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))
}
quantities <- c("lambda", "mu", "psi", "theta", "s", "rate", "rootBin", "clock",
                "charLL", "logPost", "nJumps", "meanX")

# ---- Chains -----------------------------------------------------------------
files <- list.files(dir, "^chain-.*\\.rds$", full.names = TRUE)
chains <- RbindFill(lapply(files, function(f) {
  x <- readRDS(f)
  s <- x$summary
  # Seconds of the sampling phase: total move time x nIter / (burnin + nIter)
  secs <- s$seconds * x$nIter / (x$nIter + x$burnin)
  ess <- s$ess[intersect(quantities, names(s$ess))]
  ess <- ess[quantities]
  names(ess) <- quantities
  half <- nrow(x$samples) %/% 2
  drift <- mean(x$samples$charLL[-seq_len(half)]) - mean(x$samples$charLL[seq_len(half)])
  acc <- x$acceptance
  data.frame(tree = as.integer(sub(".*-t(\\d+)-.*", "\\1", basename(f))), nChar = x$nChar,
             config = x$config, rep = x$rep, nIter = x$nIter, secs = secs,
             secPerIter = s$seconds / (x$nIter + x$burnin),
             t(setNames(ess, paste0("ess.", quantities))),
             t(setNames(ess / secs * 3600, paste0("essPerHour.", quantities))),
             essClockMin = min(s$essBranchClock, na.rm = TRUE),
             essClockMed = stats::median(s$essBranchClock, na.rm = TRUE),
             essBranchXMin = min(s$essBranchX, na.rm = TRUE),
             essBranchXMed = stats::median(s$essBranchX, na.rm = TRUE),
             charLLdrift = drift, charLLmean = mean(x$samples$charLL),
             t(setNames(vapply(intersect(c("lambda", "mu", "psi", "theta", "s", "rate"),
                                         names(x$samples)),
                               function(p) stats::sd(log(x$samples[[p]])), 0),
                        paste0("postSdLog.", intersect(c("lambda", "mu", "psi", "theta",
                                                         "s", "rate"), names(x$samples))))),
             t(setNames(acc, paste0("acc.", names(acc)))),
             t(setNames(x$seconds / sum(x$seconds), paste0("share.", names(x$seconds)))),
             check.names = FALSE)
}))
if (!is.null(chains)) {
  # Different configs have different extra-move columns
  write.csv(chains, file.path(dir, "chains.csv"), row.names = FALSE)
}

# Split R-hat (Gelman et al. 2013) over the replicate chains of each cell and
# configuration: each chain is cut in halves
SplitRhat <- function(chs) {
  halves <- unlist(lapply(chs, function(x) {
    x <- x[is.finite(x)]
    n <- length(x) %/% 2
    list(x[seq_len(n)], x[n + seq_len(n)])
  }), recursive = FALSE)
  n <- min(lengths(halves))
  m <- vapply(halves, function(h) mean(h[seq_len(n)]), 0)
  w <- mean(vapply(halves, function(h) stats::var(h[seq_len(n)]), 0))
  b <- n * stats::var(m)
  sqrt(((n - 1) / n * w + b / n) / w)
}
if (length(files)) {
  meta <- do.call(rbind, lapply(files, function(f) {
    data.frame(file = f, cell = sub("-r\\d+\\.rds$", "", basename(f)))
  }))
  rq <- c("theta", "s", "rate", "lambda", "rootBin", "clock", "charLL")
  rhat <- RbindFill(lapply(split(meta$file, meta$cell), function(fs) {
    if (length(fs) < 2) return(NULL)
    xs <- lapply(fs, function(f) readRDS(f)$samples)
    vals <- vapply(rq, function(q) {
      if (!q %in% names(xs[[1]])) return(NA_real_)
      v <- lapply(xs, function(s) if (q %in% c("theta", "s", "rate", "lambda")) log(s[[q]]) else s[[q]])
      SplitRhat(v)
    }, 0)
    data.frame(cell = sub("^chain-", "", sub("-r\\d+\\.rds$", "", basename(fs[[1]]))),
               nRep = length(fs), t(setNames(vals, paste0("rhat.", rq))), check.names = FALSE)
  }))
  if (!is.null(rhat)) write.csv(rhat, file.path(dir, "rhat.csv"), row.names = FALSE)
}

# ---- Probes -----------------------------------------------------------------
pf <- list.files(dir, "^probe-.*\\.rds$", full.names = TRUE)
probes <- do.call(rbind, lapply(pf, function(f) {
  x <- readRDS(f)
  do.call(rbind, lapply(seq_along(x$probe), function(i) {
    p <- x$probe[[i]]
    do.call(rbind, lapply(names(p$fresh), function(kind) {
      fr <- p$fresh[[kind]]
      # log acceptance of an independence proposal from the frozen law at
      # fixed parameters: (lw + char)' - (lw + char)
      d <- fr[, "lw"] + fr[, "char"] - (p$current[["lw"]] + p$current[["char"]])
      data.frame(tree = x$treeId, nChar = x$nChar, state = i, proposal = kind,
                 sdChar = stats::sd(fr[, "char"]), sdLw = stats::sd(fr[, "lw"]),
                 meanDelta = mean(d), accept = mean(pmin(1, exp(d))),
                 sdLogClock = stats::sd(log(fr[, "clock"])),
                 nJumpsFresh = mean(fr[, "nJumps"]), nJumps = p$nJumps,
                 t(setNames(p$condSdLog, paste0("condSdLog.", names(p$condSdLog)))),
                 check.names = FALSE)
    }))
  }))
}))
if (!is.null(probes)) write.csv(probes, file.path(dir, "probes.csv"), row.names = FALSE)

# ---- Oracles ----------------------------------------------------------------
BatchVar <- function(x, nBatch = 30) {
  m <- floor(length(x) / nBatch) * nBatch
  stats::var(colMeans(matrix(x[seq_len(m)], ncol = nBatch))) / nBatch
}
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
of <- list.files(dir, "^oracle-.*\\.rds$", full.names = TRUE)
# Blocks of the same variant and length (oracle-<variant>-<nIter>-bNN.rds) are pooled
key <- sub("-b\\d+$", "", sub("\\.rds$", "", basename(of)))
key[grepl("-b\\d+\\.rds$", of)] <- paste0(key[grepl("-b\\d+\\.rds$", of)], "-blocks")
oracles <- do.call(rbind, lapply(unique(key), function(k) {
  parts <- lapply(of[key == k], readRDS)
  x <- parts[[1]]
  x$chains <- do.call(c, lapply(parts, `[[`, "chains"))
  x$marginal <- do.call(c, lapply(parts, `[[`, "marginal"))
  x$variant <- paste0(x$variant, "-", x$nIter, if (length(parts) > 1) "x" else "",
                      if (length(parts) > 1) length(parts) else "")
  do.call(rbind, lapply(x$spec$free, function(p) {
    draws <- lapply(x$chains, function(ch) log(ch$samples[[p]]))
    if (!is.null(x$exact)) {
      g <- x$exact$value
      lp <- x$exact$logLik + x$exact$logPrior
      lx <- log(g)
      dens <- exp(lp - max(lp)) * g        # density of log(value)
      cdf <- c(0, cumsum((dens[-1] + dens[-length(dens)]) / 2 * diff(lx)))
      norm <- max(cdf)
      cdf <- cdf / norm
      Q <- function(pr) stats::approx(cdf, lx, pr)$y
      Dens <- function(q) stats::approx(lx, dens / norm, q)$y
      ref <- "exact"
    } else {
      mdraws <- unlist(lapply(x$marginal, function(m) log(m[[p]])))
      Q <- function(pr) stats::quantile(mdraws, pr, names = FALSE)
      dd <- stats::density(mdraws)
      Dens <- function(q) stats::approx(dd$x, dd$y, q)$y
      ref <- "marginal"
    }
    do.call(rbind, lapply(probs, function(pr) {
      qx <- Q(pr)
      ind <- lapply(draws, function(d) as.numeric(d <= qx))
      n <- vapply(ind, length, 0)
      seF <- sqrt(sum(vapply(ind, BatchVar, 0) * n^2)) / sum(n)
      # Marginal-reference SE: add the reference's own batch-means error
      if (ref == "marginal") {
        mi <- lapply(x$marginal, function(m) as.numeric(log(m[[p]]) <= qx))
        mn <- vapply(mi, length, 0)
        seF <- sqrt(seF^2 + sum(vapply(mi, BatchVar, 0) * mn^2) / sum(mn)^2)
      }
      est <- stats::quantile(unlist(draws), pr, names = FALSE)
      f0 <- Dens(qx)
      data.frame(variant = x$variant, parameter = p, reference = ref, p = pr,
                 reference.q = exp(qx), mcmc.q = exp(est), z = (est - qx) / (seF / f0))
    }))
  }))
}))
if (!is.null(oracles)) write.csv(oracles, file.path(dir, "oracles.csv"), row.names = FALSE)
cat(sprintf("chains %d, probes %d, oracle rows %d\n", length(files), length(pf),
            if (is.null(oracles)) 0L else nrow(oracles)))
