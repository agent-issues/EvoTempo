# Trilobite CET pipeline: one MCMC chain of one (beta, gamma) cell, run as a
# series of segments with a checkpoint after each (as the penguin pipeline).
#
#   Rscript data-raw/empirical/trilobites/cet/run-chain.R <cell> <rep>
#
# <cell> 1-4 = (beta, gamma) (0,0), (0,1), (1,0), (1,1), or a label such as
# b1g0.5 for a model-index ladder rung.  Run from the package root, with
# EvoTempo installed (or CET_LOAD_ALL=1 for pkgload::load_all()).  Environment:
#   CET_OUT          output root (default cet-runs); the chain writes
#                    <CET_OUT>/<label>_run<rep>/
#   CET_START_TREE   starting tree (default cet/data/start-tree-<kk>.nwk, draw
#                    k = 4 (rep - 1) + cell of prepare.R, a random sample from
#                    the pooled RevBayes posterior; a ladder rung or power chain
#                    uses draw <rep>); the starting parameter values come from
#                    the same RevBayes sample (start-values.csv), whatever
#                    CET_START_TREE says
#   CET_NITER        sampling iterations in total (default 100000)
#   CET_BURNIN       burn-in iterations, with step-size adaptation (default 2000)
#   CET_THIN         thinning (default 50)
#   CET_SEG          sampling iterations per segment (default 2000; a multiple
#                    of CET_THIN)
#   CET_HOURS        stop cleanly after the segment that passes this many hours
#                    (default 70, inside the 3-day `shared` limit); resubmit
#                    the same command to continue from the checkpoint
#   CET_POWER        power on the character likelihood (default 1)
#   CET_MKV          "pattern" (default; one Mkv correction per missing-data
#                    pattern, as RevBayes coding = "variable") or "none"; a
#                    non-default variant is part of the label
#   CET_SEED         RNG seed (default 1000 * cell + rep)
#   CET_NSLICE       frozen-model time slices (default 20)
#   CET_TREE_MOVES   scale of the tree-move counts (default 1)
#   CET_NBIN, CET_ZMAX   tempo grid (default 26 bins, z in [-5, 5])
#   CET_STATES       "true" (default) to save per-sample states for
#                    bayes-factors.R
# The checkpoint records the target (beta, gamma, power, Mkv variant, grid) and
# the thinning and segment length, and a run with any of them different stops
# rather than continue it (#74, #76).
# Outputs per segment k: seg<k>_samples.csv (TempoMCMC() samples plus charLL;
# logPrior, the log prior of the parameters; logOriginPrior, that of the origin
# age, which TempoMCMC()'s logPost also includes; the iteration), seg<k>_trees.nwk
# (ages in Myr before 485.4 Ma), seg<k>_states.rds; acceptance.csv (per segment,
# with constraint rejections per tree move); checkpoint.rds; progress.txt.
t0 <- proc.time()[["elapsed"]]
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
EnvNum <- function(name, default) as.numeric(Sys.getenv(name, as.character(default)))
cellArg <- args[[1]]
rep <- as.integer(args[[2]])
if (grepl("^[0-9]+$", cellArg)) {
  cell <- as.integer(cellArg)
  beta <- cellTable$beta[cell]
  gamma <- cellTable$gamma[cell]
} else {
  cell <- 0L
  beta <- as.numeric(sub("^b([0-9.]+)g.*$", "\\1", cellArg))
  gamma <- as.numeric(sub("^b[0-9.]+g([0-9.]+)$", "\\1", cellArg))
}
mkv <- .cet$mkv
label <- paste0(CellLabel(beta, gamma), MkvSuffix(mkv))
power <- EnvNum("CET_POWER", 1)
if (power != 1) label <- sprintf("%s_t%s", label, format(power))
outDir <- file.path(Sys.getenv("CET_OUT", "cet-runs"), sprintf("%s_run%d", label, rep))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
nIter <- EnvNum("CET_NITER", 100000)
burnin <- EnvNum("CET_BURNIN", 2000)
thin <- EnvNum("CET_THIN", 50)
seg <- EnvNum("CET_SEG", 2000)
hours <- EnvNum("CET_HOURS", 70)
nSlice <- EnvNum("CET_NSLICE", 20)
grid <- GridSettings()
saveStates <- Sys.getenv("CET_STATES", "true") == "true"
stopifnot(seg %% thin == 0, nIter %% seg == 0)
.cet$power <- power

startIndex <- StartIndex(cell, rep)
startTree <- normalizePath(Sys.getenv("CET_START_TREE", StartTreeFile(startIndex)), mustWork = TRUE)
dat <- LoadTrilobites(startTree)
model <- CetModel(beta, gamma)
priors <- CetPriors()
moves <- CetMoves(dat$tree, dat$characters, model)
Progress <- function(...) {
  cat(sprintf(...), "\n", file = file.path(outDir, "progress.txt"), append = TRUE)
}
target <- list(mkv = mkv, power = power, beta = beta, gamma = gamma, nBin = grid$nBin,
               zMax = grid$zMax, thin = thin, seg = seg)
ckpFile <- file.path(outDir, "checkpoint.rds")
if (file.exists(ckpFile)) {
  ckp <- readRDS(ckpFile)
  same <- vapply(names(target), function(n) identical(ckp[[n]], target[[n]]), TRUE)
  if (!all(same)) {
    stop(sprintf("%s holds a chain with %s; this run has %s", ckpFile,
                 paste(names(target), vapply(names(target), function(n) format(ckp[[n]]), ""),
                       collapse = ", "),
                 paste(names(target), vapply(target, format, ""), collapse = ", ")))
  }
  # The characters must be in the tip order of the chain's tree (review of #91, F1)
  if (!identical(ckp$startTree, startTree) ||
      !identical(ckp$captured$st$tree$tipLabel, rownames(dat$characters))) {
    stop(sprintf("%s was started from %s; this run loads %s (tip orders differ)", ckpFile,
                 format(ckp$startTree), startTree))
  }
  Progress("%s resuming after segment %d", format(Sys.time()), ckp$segment)
} else {
  set.seed(EnvNum("CET_SEED", 1000 * cell + rep))
  ckp <- NULL
  Progress(paste("%s start %s rep %d (beta %g, gamma %g, power %g, mkv %s, grid %d bins z +-%g):",
                 "start tree %s (draw %d); nIter %d, burnin %d, thin %d, seg %d; moves %s"),
           format(Sys.time()), label, rep, beta, gamma, power, mkv, grid$nBin, grid$zMax,
           startTree, startIndex, nIter, burnin, thin, seg,
           paste(names(moves), moves, sep = "=", collapse = " "))
}
nSeg <- nIter %/% seg
first <- if (is.null(ckp)) 1L else ckp$segment + 1L
for (k in seq.int(first, length.out = max(0L, nSeg - first + 1L))) {
  .cet$records <- list()
  .cet$injected <- FALSE
  .cet$constraintRejections <- list()
  tSeg <- proc.time()[["elapsed"]]
  stopifnot(identical(rownames(dat$characters), dat$tree$tipLabel))
  if (is.null(ckp)) {
    run <- TempoMCMC(dat$tree, model, CetStart(startIndex), priors, nIter = seg, burnin = burnin,
                     thin = thin, nSlice = nSlice, moves = moves,
                     characters = dat$characters, agePriors = CetAgePriors())
  } else {
    .cet$resume <- ckp$captured
    run <- TempoMCMC(ckp$captured$st$tree, model, ckp$captured$st$par, priors,
                     nIter = seg, burnin = 0, thin = thin, nSlice = nSlice,
                     moves = ckp$captured$moves, characters = dat$characters,
                     agePriors = CetAgePriors())
    stopifnot(.cet$injected)
  }
  secs <- proc.time()[["elapsed"]] - tSeg
  rec <- .cet$records
  stopifnot(length(rec) == nrow(run$samples))
  smp <- run$samples
  smp$charLL <- vapply(rec, `[[`, 0, "charLL")
  smp$logPrior <- vapply(seq_len(nrow(smp)), function(i) {
    sum(vapply(names(priors), function(p) priors[[p]](rec[[i]]$state$par[[p]]), 0))
  }, 0)
  smp$logOriginPrior <- vapply(rec, function(r) CetAgePriors()$origin(
    r$state$tree$age[[r$state$tree$root]] + r$state$tree$stem), 0)
  smp$iteration <- (k - 1) * seg + seq(thin, seg, by = thin)
  smp$segment <- k
  tag <- sprintf("seg%04d", k)
  utils::write.csv(smp, file.path(outDir, paste0(tag, "_samples.csv")), row.names = FALSE)
  writeLines(vapply(rec, `[[`, "", "newick"), file.path(outDir, paste0(tag, "_trees.nwk")))
  if (saveStates) {
    saveRDS(lapply(rec, `[[`, "state"), file.path(outDir, paste0(tag, "_states.rds")))
  }
  cr <- unlist(.cet$constraintRejections)
  acc <- data.frame(segment = k, move = names(run$acceptance), rate = unname(run$acceptance),
                    constraintRejections = ifelse(names(run$acceptance) %in% names(cr),
                                                  cr[names(run$acceptance)], 0),
                    secs = secs, burnin = if (is.null(ckp)) burnin else 0)
  accFile <- file.path(outDir, "acceptance.csv")
  utils::write.table(acc, accFile, sep = ",", row.names = FALSE, col.names = !file.exists(accFile),
                     append = file.exists(accFile))
  ckp <- c(list(segment = k, captured = .cet$captured, label = label, rep = rep,
                startIndex = startIndex, startTree = startTree), target)
  tmp <- paste0(ckpFile, ".tmp")
  saveRDS(ckp, tmp)
  file.rename(tmp, ckpFile)
  Progress("%s segment %d/%d done: %.1f s (%.3f s/iteration incl. %d burn-in); logPost %.2f, charLL %.2f",
           format(Sys.time()), k, nSeg, secs, secs / (seg + acc$burnin[1]), acc$burnin[1],
           utils::tail(smp$logPost, 1), utils::tail(smp$charLL, 1))
  if ((proc.time()[["elapsed"]] - t0) / 3600 > hours && k < nSeg) {
    Progress("%s stopping after %.1f h; resubmit to continue", format(Sys.time()), hours)
    break
  }
}
if (!is.null(ckp) && ckp$segment == nSeg) Progress("%s complete", format(Sys.time()))
