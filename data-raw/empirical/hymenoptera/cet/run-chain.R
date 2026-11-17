# Hymenoptera CET pipeline (as the penguin pipeline's run-chain.R): one MCMC
# chain of one (beta, gamma) cell, run as a series of segments with a
# checkpoint after each.
#
#   Rscript data-raw/empirical/hymenoptera/cet/run-chain.R <cell> <rep>
#
# <cell> 1-4 = (beta, gamma) (0,0), (0,1), (1,0), (1,1), or a label such as
# b1g0.5 for a model-index ladder rung.  Run from the package root, with
# EvoTempo installed (or CET_LOAD_ALL=1 for pkgload::load_all()).  Environment:
#   CET_OUT          output root (default cet-runs); the chain writes
#                    <CET_OUT>/<label>_run<rep>/
#   CET_START_TREE   starting tree (default: production chain k's own draw,
#                    data/start-trees/start-tree-<kk>.nwk, else the reference
#                    tree data/start-tree.nwk; common.R, StartTreeFor())
#   CET_NITER        sampling iterations in total (default 100000)
#   CET_BURNIN       burn-in iterations, with step-size adaptation (default 2000)
#   CET_THIN         thinning (default 50)
#   CET_SEG          sampling iterations per segment (default 2000; a multiple
#                    of CET_THIN)
#   CET_HOURS        stop cleanly after the segment that passes this many hours
#                    (default 70, inside the 3-day `shared` limit); resubmit
#                    the same command to continue from the checkpoint
#   CET_POWER        power on the character likelihood (default 1)
#   CET_MKV          Mkv correction: "pattern" (default; one per missing-data
#                    mask, as RevBayes coding = "variable") or "none".  A
#                    non-default variant is part of the label
#                    (<cell>_mkv<variant>[_t<power>]), and the variant is
#                    stored in the checkpoint and checked on resume (#74)
#   CET_SEED         RNG seed (default 1000 * cell + rep)
#   CET_NSLICE       frozen-model time slices (default 20)
#   CET_TREE_MOVES   scale of the tree-move counts (default 1)
#   CET_FBD_BLOCK    rounds of FBD block moves per iteration (default 1)
#   CET_NBIN, CET_ZMAX   tempo grid (default 26 bins, z in [-5, 5])
#   CET_NRATECAT     Gamma rate categories among characters (default 4; 1 is a
#                    diagnostic variant without rate variation)
#   CET_HOLO_CAL     "true" (default) or "false": Holometabola calibration on
#                    or off (constraint kept), for prior-only comparisons
#   CET_FBD_SCALE    "none" (default) or "lineage": the comparator's FBD-rate
#                    prior on the rates at x = 0 or on the lineage-averaged
#                    rates (common.R; diagnosis of cell (1, 0))
#   CET_S_MAX        upper bound of the prior on s (default 3)
#   CET_STATES       "true" (default) to save per-sample states for
#                    bayes-factors.R
# Outputs per segment k: seg<k>_samples.csv (TempoMCMC() samples, alpha
# included, plus charLL, the Mkv-corrected character log-likelihood (Gamma
# mixture included) at power 1; treeTerm, the Holometabola constraint and
# calibration term, which logPost includes; logPrior, the log prior of the
# parameters (joint FBD-rate prior and the alpha prior included), as in
# bayes-factors.R's cross.csv; logRootPrior, the log prior of the root age,
# which TempoMCMC()'s logPost also includes; holoAge, hymCrownAge,
# hymTotalAge, the ages of the comparator's monitored clades; and the
# iteration), seg<k>_trees.nwk, seg<k>_states.rds; acceptance.csv (per
# segment); checkpoint.rds (state at the end of the last complete segment;
# written atomically); progress.txt.
t0 <- proc.time()[["elapsed"]]
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/hymenoptera/cet/common.R")
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
  # Plain decimal notation only: "b1e-04g1" would give beta NA (#75)
  if (!grepl("^b[0-9.]+g[0-9.]+$", cellArg)) {
    stop("cannot parse cell label '", cellArg, "': use b<beta>g<gamma> in plain decimal notation")
  }
  beta <- suppressWarnings(as.numeric(sub("^b([0-9.]+)g.*$", "\\1", cellArg)))
  gamma <- suppressWarnings(as.numeric(sub("^b[0-9.]+g([0-9.]+)$", "\\1", cellArg)))
  if (anyNA(c(beta, gamma))) stop("cannot parse cell label '", cellArg, "'")
}
mkv <- .cet$mkv
label <- paste0(CellLabel(beta, gamma), MkvSuffix(mkv))
power <- EnvNum("CET_POWER", 1)
# FormatValue(): plain decimal notation that reads back as exactly `power`, so
# that bayes-factors.R's ParseLabel() recovers it (#75)
if (power != 1) label <- sprintf("%s_t%s", label, FormatValue(power))
parsed <- ParseLabel(sprintf("%s_run%d", label, rep))
stopifnot(parsed$beta == beta, parsed$gamma == gamma, parsed$power == power, parsed$mkv == mkv)
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

# Start tree (#89 review, R89-5): production chain k = 4 (rep - 1) + cell
# starts from its own draw from the comparator's posterior (common.R,
# StartTreeFor()).  The tree and its md5 go into the checkpoint, and a resumed
# chain must use the same file: the character rows follow the start tree's
# tip order (LoadHymenoptera()).
chainK <- ChainIndex(cell, rep)
startSrc <- StartTreeFor(chainK)
startTree <- startSrc$file
startTreeMd5 <- unname(tools::md5sum(startTree))
dat <- LoadHymenoptera(startTree)
model <- CetModel(beta, gamma)
priors <- CetPriors(beta)
charPriors <- if (.cet$nRateCat > 1L) CetCharPriors() else list()
moves <- CetMoves(dat$tree, dat$characters, model)
Progress <- function(...) {
  cat(sprintf(...), "\n", file = file.path(outDir, "progress.txt"), append = TRUE)
}
ckpFile <- file.path(outDir, "checkpoint.rds")
if (file.exists(ckpFile)) {
  ckp <- readRDS(ckpFile)
  # A checkpoint of another target must not be continued (#74)
  ckpMkv <- if (is.null(ckp$mkv)) "unrecorded (checkpoint from before #74)" else ckp$mkv
  if (!identical(ckp$mkv, mkv) || !identical(ckp$power, power) ||
      !identical(ckp$beta, beta) || !identical(ckp$gamma, gamma)) {
    stop(sprintf(paste("%s holds a chain with mkv %s, power %s, beta %s, gamma %s;",
                       "this run has mkv %s, power %s, beta %s, gamma %s"),
                 ckpFile, ckpMkv, format(ckp$power), format(ckp$beta), format(ckp$gamma),
                 mkv, format(power), format(beta), format(gamma)))
  }
  if (!identical(ckp$nRateCat, .cet$nRateCat)) {
    stop(sprintf("%s holds a chain with %s Gamma rate categories; this run has %d", ckpFile,
                 format(ckp$nRateCat), .cet$nRateCat))
  }
  if (!identical(if (is.null(ckp$fbdScale)) "none" else ckp$fbdScale, .cet$fbdScale) ||
      !identical(if (is.null(ckp$sMax)) 3 else ckp$sMax, .cet$sMax)) {
    stop(sprintf("%s holds a chain with CET_FBD_SCALE %s, CET_S_MAX %s; this run has %s, %s",
                 ckpFile, format(ckp$fbdScale), format(ckp$sMax), .cet$fbdScale, format(.cet$sMax)))
  }
  if (!identical(if (is.null(ckp$holoCal)) TRUE else ckp$holoCal, .cet$holoCal)) {
    stop(sprintf("%s holds a chain with Holometabola calibration %s; this run has %s", ckpFile,
                 format(ckp$holoCal), format(.cet$holoCal)))
  }
  if (!identical(ckp$startTree, startTree) || !identical(ckp$startTreeMd5, startTreeMd5)) {
    stop(sprintf("%s started from %s (md5 %s); this run would load %s (md5 %s)", ckpFile,
                 format(ckp$startTree), format(ckp$startTreeMd5), startTree, startTreeMd5))
  }
  if (!identical(ckp$captured$st$tree$tipLabel, rownames(dat$cells))) {
    stop("the checkpoint's tip order differs from the character rows of ", startTree)
  }
  Progress("%s resuming after segment %d", format(Sys.time()), ckp$segment)
} else {
  set.seed(EnvNum("CET_SEED", 1000 * cell + rep))
  ckp <- NULL
  if (is.null(startSrc$sample)) {
    stop("no comparator sample for the start tree ", startTree, " (CET_START_TREE given, or ",
         "data/start-tree-sample.csv missing): CetStart() needs its parameters")
  }
  # Starting parameters of the comparator sample that gave the start tree;
  # theta and s seeded by chain (production chain k) or by replicate (other
  # labels), redrawn while the start log posterior is not finite
  startIndex <- if (is.na(chainK)) rep else chainK
  startPar <- CetStart(startIndex, startSrc$sample, dat$tree, model, dat$characters, priors,
                       nSlice)
  startInfo <- attributes(startPar)[c("attempts", "logPost")]
  startPar <- as.list(stats::setNames(unlist(startPar), names(startPar)))
  # Without Gamma rates (CET_NRATECAT = 1) alpha is not a parameter
  if (.cet$nRateCat == 1L) startPar$alpha <- NULL
  Progress("%s start %s rep %d (beta %g, gamma %g, power %g, mkv %s, %d Gamma categories): nIter %d, burnin %d, thin %d, seg %d; moves %s",
           format(Sys.time()), label, rep, beta, gamma, power, mkv, .cet$nRateCat, nIter, burnin,
           thin, seg, paste(names(moves), moves, sep = "=", collapse = " "))
  Progress("%s start tree %s (md5 %s; comparator run%d iteration %d); start parameters (CetStart(%d), %d draw(s), start logPost %.2f): %s",
           format(Sys.time()), startTree, startTreeMd5, as.integer(startSrc$sample$run),
           as.integer(startSrc$sample$Iteration), startIndex, startInfo$attempts, startInfo$logPost,
           paste(names(startPar), signif(unlist(startPar), 5), sep = "=", collapse = " "))
}
nSeg <- nIter %/% seg
first <- if (is.null(ckp)) 1L else ckp$segment + 1L
for (k in seq.int(first, length.out = max(0L, nSeg - first + 1L))) {
  .cet$records <- list()
  .cet$injected <- FALSE
  .cet$blockTry[] <- 0
  .cet$blockAcc[] <- 0
  tSeg <- proc.time()[["elapsed"]]
  if (is.null(ckp)) {
    run <- TempoMCMC(dat$tree, model, startPar, priors, nIter = seg, burnin = burnin,
                     thin = thin, nSlice = nSlice, moves = moves,
                     characters = dat$characters, agePriors = CetAgePriors(),
                     charPriors = charPriors)
  } else {
    .cet$resume <- ckp$captured
    run <- TempoMCMC(ckp$captured$st$tree, model, ckp$captured$st$par, priors,
                     nIter = seg, burnin = 0, thin = thin, nSlice = nSlice,
                     moves = ckp$captured$moves, characters = dat$characters,
                     agePriors = CetAgePriors(), charPriors = charPriors)
    stopifnot(.cet$injected)
  }
  secs <- proc.time()[["elapsed"]] - tSeg
  rec <- .cet$records
  stopifnot(length(rec) == nrow(run$samples))
  smp <- run$samples
  smp$charLL <- vapply(rec, `[[`, 0, "charLL")
  smp$treeTerm <- vapply(rec, `[[`, 0, "treeTerm")
  # logPrior: the parameters only (alpha included), as cross.csv; the root age
  # prior, which the chain's logPost includes, is a separate column (#74)
  smp$logPrior <- vapply(rec, function(r) CetLogPrior(r$state$par, priors, charPriors), 0)
  smp$logRootPrior <- vapply(rec, function(r) CetAgePriors()$root(
    r$state$tree$age[[r$state$tree$root]]), 0)
  ages <- do.call(rbind, lapply(rec, `[[`, "ages"))
  smp <- cbind(smp, ages)
  smp$iteration <- (k - 1) * seg + seq(thin, seg, by = thin)
  smp$segment <- k
  tag <- sprintf("seg%04d", k)
  utils::write.csv(smp, file.path(outDir, paste0(tag, "_samples.csv")), row.names = FALSE)
  writeLines(vapply(rec, `[[`, "", "newick"), file.path(outDir, paste0(tag, "_trees.nwk")))
  if (saveStates) {
    saveRDS(lapply(rec, `[[`, "state"), file.path(outDir, paste0(tag, "_states.rds")))
  }
  # FBD block moves (common.R) run in the shift slot; their acceptance is
  # reported as moves fbdBlock_<coordinate>
  acc <- data.frame(segment = k,
                    move = c(names(run$acceptance), paste0("fbdBlock_", blockNames)),
                    rate = c(unname(run$acceptance), .cet$blockAcc / pmax(.cet$blockTry, 1)),
                    secs = secs, burnin = if (is.null(ckp)) burnin else 0)
  accFile <- file.path(outDir, "acceptance.csv")
  utils::write.table(acc, accFile, sep = ",", row.names = FALSE, col.names = !file.exists(accFile),
                     append = file.exists(accFile))
  ckp <- list(segment = k, captured = .cet$captured, label = label, beta = beta, gamma = gamma,
              power = power, mkv = mkv, nRateCat = .cet$nRateCat, holoCal = .cet$holoCal,
              fbdScale = .cet$fbdScale, sMax = .cet$sMax,
              rep = rep, thin = thin,
              seg = seg, startTree = startTree, startTreeMd5 = startTreeMd5,
              start = if (is.null(ckp)) list(par = startPar, info = startInfo,
                                             sample = startSrc$sample) else ckp$start)
  tmp <- paste0(ckpFile, ".tmp")
  saveRDS(ckp, tmp)
  file.rename(tmp, ckpFile)
  Progress("%s segment %d/%d done: %.1f s (%.3f s/iteration incl. %d burn-in); logPost %.2f, charLL %.2f, alpha %.4g",
           format(Sys.time()), k, nSeg, secs, secs / (seg + acc$burnin[1]), acc$burnin[1],
           utils::tail(smp$logPost, 1), utils::tail(smp$charLL, 1),
           if (is.null(smp$alpha)) NA else utils::tail(smp$alpha, 1))
  if ((proc.time()[["elapsed"]] - t0) / 3600 > hours && k < nSeg) {
    Progress("%s stopping after %.1f h; resubmit to continue", format(Sys.time()), hours)
    break
  }
}
if (!is.null(ckp) && ckp$segment == nSeg) Progress("%s complete", format(Sys.time()))
