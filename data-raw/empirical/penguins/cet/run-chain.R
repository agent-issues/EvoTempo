# Penguin CET pipeline: one MCMC chain of one (beta, gamma) cell, run as a
# series of segments with a checkpoint after each.
#
#   Rscript data-raw/empirical/penguins/cet/run-chain.R <cell> <rep>
#
# <cell> 1-4 = (beta, gamma) (0,0), (0,1), (1,0), (1,1), or a label such as
# b1g0.5 for a model-index ladder rung.  Run from the package root, with
# EvoTempo installed (or CET_LOAD_ALL=1 for pkgload::load_all()).  Environment:
#   CET_OUT          output root (default cet-runs); the chain writes
#                    <CET_OUT>/<label>_run<rep>/
#   CET_START_TREE   starting tree (default cet/data/start-tree.nwk)
#   CET_NITER        sampling iterations in total (default 100000)
#   CET_BURNIN       burn-in iterations, with step-size adaptation (default 2000)
#   CET_THIN         thinning (default 50)
#   CET_SEG          sampling iterations per segment (default 2000; a multiple
#                    of CET_THIN)
#   CET_HOURS        stop cleanly after the segment that passes this many hours
#                    (default 70, inside the 3-day `shared` limit); resubmit
#                    the same command to continue from the checkpoint
#   CET_POWER        power on the character likelihood (default 1)
#   CET_SEED         RNG seed (default 1000 * cell + rep)
#   CET_NSLICE       frozen-model time slices (default 20)
#   CET_TREE_MOVES   scale of the tree-move counts (default 1)
#   CET_NBIN, CET_ZMAX   tempo grid (default 26 bins, z in [-5, 5])
#   CET_STATES       "true" (default) to save per-sample states for
#                    bayes-factors.R
# Outputs per segment k: seg<k>_samples.csv (TempoMCMC() samples plus the
# character log-likelihood, log prior and iteration), seg<k>_trees.nwk,
# seg<k>_states.rds; acceptance.csv (per segment); checkpoint.rds (state at the
# end of the last complete segment; written atomically); progress.txt.
t0 <- proc.time()[["elapsed"]]
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/penguins/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
EnvNum <- function(name, default) as.numeric(Sys.getenv(name, as.character(default)))
cellArg <- args[[1]]
rep <- as.integer(args[[2]])
if (grepl("^[0-9]+$", cellArg)) {
  cell <- as.integer(cellArg)
  beta <- cells$beta[cell]
  gamma <- cells$gamma[cell]
} else {
  cell <- 0L
  beta <- as.numeric(sub("^b([0-9.]+)g.*$", "\\1", cellArg))
  gamma <- as.numeric(sub("^b[0-9.]+g([0-9.]+)$", "\\1", cellArg))
}
label <- CellLabel(beta, gamma)
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
saveStates <- Sys.getenv("CET_STATES", "true") == "true"
stopifnot(seg %% thin == 0, nIter %% seg == 0)
.cet$power <- power

dat <- LoadPenguins(Sys.getenv("CET_START_TREE", file.path(cetDir, "data", "start-tree.nwk")))
model <- CetModel(beta, gamma)
priors <- CetPriors()
moves <- CetMoves(dat$tree, dat$characters, model)
Progress <- function(...) {
  cat(sprintf(...), "\n", file = file.path(outDir, "progress.txt"), append = TRUE)
}
ckpFile <- file.path(outDir, "checkpoint.rds")
if (file.exists(ckpFile)) {
  ckp <- readRDS(ckpFile)
  Progress("%s resuming after segment %d", format(Sys.time()), ckp$segment)
} else {
  set.seed(EnvNum("CET_SEED", 1000 * cell + rep))
  ckp <- NULL
  Progress("%s start %s rep %d (beta %g, gamma %g, power %g): nIter %d, burnin %d, thin %d, seg %d; moves %s",
           format(Sys.time()), label, rep, beta, gamma, power, nIter, burnin, thin, seg,
           paste(names(moves), moves, sep = "=", collapse = " "))
}
nSeg <- nIter %/% seg
first <- if (is.null(ckp)) 1L else ckp$segment + 1L
for (k in seq.int(first, length.out = max(0L, nSeg - first + 1L))) {
  .cet$records <- list()
  .cet$injected <- FALSE
  tSeg <- proc.time()[["elapsed"]]
  if (is.null(ckp)) {
    run <- TempoMCMC(dat$tree, model, CetStart(rep), priors, nIter = seg, burnin = burnin,
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
  }, 0) + vapply(rec, function(r) CetAgePriors()$origin(
    r$state$tree$age[[r$state$tree$root]] + r$state$tree$stem), 0)
  smp$iteration <- (k - 1) * seg + seq(thin, seg, by = thin)
  smp$segment <- k
  tag <- sprintf("seg%04d", k)
  utils::write.csv(smp, file.path(outDir, paste0(tag, "_samples.csv")), row.names = FALSE)
  writeLines(vapply(rec, `[[`, "", "newick"), file.path(outDir, paste0(tag, "_trees.nwk")))
  if (saveStates) {
    saveRDS(lapply(rec, `[[`, "state"), file.path(outDir, paste0(tag, "_states.rds")))
  }
  acc <- data.frame(segment = k, move = names(run$acceptance), rate = unname(run$acceptance),
                    secs = secs, burnin = if (is.null(ckp)) burnin else 0)
  accFile <- file.path(outDir, "acceptance.csv")
  utils::write.table(acc, accFile, sep = ",", row.names = FALSE, col.names = !file.exists(accFile),
                     append = file.exists(accFile))
  ckp <- list(segment = k, captured = .cet$captured, label = label, beta = beta, gamma = gamma,
              power = power, rep = rep, thin = thin, seg = seg)
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
