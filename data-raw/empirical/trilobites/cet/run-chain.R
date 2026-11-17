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
#   CET_HOURS        hard cap: stop cleanly after the segment that passes this
#                    many hours (default 65); resubmit the same command to
#                    continue from the checkpoint
#   CET_WALL_HOURS   the job's wall-clock limit in hours (default 72, `shared`;
#                    it is not read from Slurm: set it to match a shorter
#                    --time),
#   CET_WALL_MARGIN  less this margin (default 0.25): before each segment, stop
#                    cleanly if the time elapsed plus 1.5 x the last segment's
#                    time would pass CET_WALL_HOURS - CET_WALL_MARGIN, so that
#                    Slurm never kills a chain inside a segment (review of
#                    #114, R114-1: a 10 000-iteration segment takes 2.6-4.4 h)
#   CET_ALLOW_REINSTALL  1 to resume under a reinstall of the same commit:
#                    compare EvoTempo's version and BUILD_SHA but not the md5
#                    of the installed code (do not rebuild lib-<sha> while
#                    chains await resubmission; this is the way out if it was)
#   CET_TEST_KILL_SEGMENT  testing only (smoke.sh task 10): quit with status
#                    137 after writing segment k's acceptance.csv rows and
#                    before its checkpoint, as a kill in that window would
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
# The checkpoint records the target (beta, gamma, power, Mkv variant, Gamma
# rate categories, the alpha prior, md5s of data/characters.csv and
# data/partitions.csv, the EvoTempo build: packageVersion(), the commit SHA
# recorded in the install (BUILD_SHA) and an md5 of the installed code (R114-3,
# R114-6); grid) and
# the thinning and segment length, and a run with any of them different stops
# rather than continue it (#74, #76).  It also records the start tree (md5,
# and path relative to cet/), and a resume from a tree with another md5 or tip
# order stops; the checkout path does not matter (#102).
# Outputs per segment k: seg<k>_samples.csv (TempoMCMC() samples, with the
# Gamma shape alpha, plus charLL;
# logPrior, the log prior of the parameters; logOriginPrior, that of the origin
# age, and logAlphaPrior, that of alpha, both of which TempoMCMC()'s logPost
# also includes; the iteration), seg<k>_trees.nwk
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
if (power != 1) label <- sprintf("%s_t%s", label, FormatValue(power))
outDir <- file.path(Sys.getenv("CET_OUT", "cet-runs"), sprintf("%s_run%d", label, rep))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
nIter <- EnvNum("CET_NITER", 100000)
burnin <- EnvNum("CET_BURNIN", 2000)
thin <- EnvNum("CET_THIN", 50)
seg <- EnvNum("CET_SEG", 2000)
hours <- EnvNum("CET_HOURS", 65)
wallHours <- EnvNum("CET_WALL_HOURS", 72)
wallMargin <- EnvNum("CET_WALL_MARGIN", 0.25)
killSegment <- EnvNum("CET_TEST_KILL_SEGMENT", 0)
nSlice <- EnvNum("CET_NSLICE", 20)
grid <- GridSettings()
saveStates <- Sys.getenv("CET_STATES", "true") == "true"
stopifnot(seg %% thin == 0, nIter %% seg == 0)
.cet$power <- power

startIndex <- StartIndex(cell, rep)
startTree <- normalizePath(Sys.getenv("CET_START_TREE", StartTreeFile(startIndex)), mustWork = TRUE)
# The checkpoint identifies the start tree by its content (md5) and by its path
# relative to cet/ when it lies there, so a resubmission from another checkout
# of the same commit resumes, and a different tree stops (#102)
startTreeMd5 <- unname(tools::md5sum(startTree))
startTreeRel <- local({
  base <- paste0(normalizePath(cetDir, mustWork = TRUE), "/")
  if (startsWith(startTree, base)) substring(startTree, nchar(base) + 1L) else startTree
})
dat <- LoadTrilobites(startTree)
model <- CetModel(beta, gamma)
priors <- CetPriors()
moves <- CetMoves(dat$tree, dat$characters, model)
Progress <- function(...) {
  cat(sprintf(...), "\n", file = file.path(outDir, "progress.txt"), append = TRUE)
}
# nRateCat: a checkpoint from before Gamma rates (#97), which records none, stops
# The build and the data: a resubmission must resume under the same EvoTempo
# build, characters and alpha prior, or stop (review of #114, R114-3, R114-6)
build <- BuildId()
target <- list(mkv = mkv, power = power, beta = beta, gamma = gamma, nBin = grid$nBin,
               zMax = grid$zMax, thin = thin, seg = seg, nRateCat = cetRateCats,
               alphaPrior = CetCharPriorTag(), alphaPriorCode = CetCharPriorCode(),
               charactersMd5 = DataMd5("characters.csv"), partitionsMd5 = DataMd5("partitions.csv"),
               evotempoVersion = build$version, evotempoSha = build$sha, evotempoMd5 = build$md5)
ckpFile <- file.path(outDir, "checkpoint.rds")
if (file.exists(ckpFile)) {
  ckp <- readRDS(ckpFile)
  same <- vapply(names(target), function(n) identical(ckp[[n]], target[[n]]), TRUE)
  # CET_ALLOW_REINSTALL=1: a reinstall of the same commit changes the code md5
  # but not the version or SHA; compare those only (#114 delta review, N2).
  # Refused when the SHA is unrecorded, which would leave nothing to compare
  reinstall <- Sys.getenv("CET_ALLOW_REINSTALL") == "1" && !same[["evotempoMd5"]] &&
    all(same[setdiff(names(same), "evotempoMd5")])
  if (reinstall) {
    if (identical(target$evotempoSha, "unrecorded")) {
      stop("CET_ALLOW_REINSTALL needs the build's commit in <lib>/EvoTempo/BUILD_SHA; ",
           "this install records none")
    }
    same[["evotempoMd5"]] <- TRUE
  }
  if (!all(same)) {
    # Only the fields that differ (review of #114: the build and data fields
    # made the full list hard to read)
    bad <- names(target)[!same]
    stop(sprintf("%s holds a chain with another target; it differs in: %s", ckpFile,
                 paste(sprintf("%s (checkpoint %s; this run %s)", bad,
                               vapply(bad, function(n) format(ckp[[n]] %||% "(none)"), ""),
                               vapply(target[bad], format, "")), collapse = "; ")))
  }
  # The characters must be in the tip order of the chain's tree (review of #91,
  # F1).  Same start tree = same md5 (#102; checkpoints written before #102
  # carry no md5 and fall back to the absolute path)
  sameTree <- if (is.null(ckp$startTreeMd5)) {
    identical(ckp$startTree, startTree)
  } else identical(ckp$startTreeMd5, startTreeMd5)
  if (!sameTree || !identical(ckp$captured$st$tree$tipLabel, dat$characters$tipLabel)) {
    if (is.null(ckp$startTreeMd5)) {
      # Checkpoint written before #102: no md5, so only the path can be compared
      # (#103 review, R103-11)
      stop(sprintf(paste("%s predates #102 and records no start-tree md5, so the start tree was",
                         "compared by absolute path: the chain started from %s; this run loads %s%s.",
                         "Resume from the original path, or restart the chain"), ckpFile,
                   format(ckp$startTree), startTree,
                   if (sameTree) " (same path, but the tip orders differ)" else ""))
    }
    stop(sprintf("%s was started from %s (md5 %s); this run loads %s (md5 %s)%s", ckpFile,
                 format(ckp$startTreeRel %||% ckp$startTree), ckp$startTreeMd5,
                 startTreeRel, startTreeMd5,
                 if (sameTree) ": tip orders differ" else ": a different start tree"))
  }
  if (reinstall) {
    Progress("%s CET_ALLOW_REINSTALL: same version %s and commit %s, installed code md5 %s (was %s)",
             format(Sys.time()), target$evotempoVersion, target$evotempoSha, target$evotempoMd5,
             ckp$evotempoMd5)
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
# Time of the last segment run, for the predictive stop (from the checkpoint
# when resuming)
lastSecs <- if (is.null(ckp$segSecs)) NA_real_ else ckp$segSecs
for (k in seq.int(first, length.out = max(0L, nSeg - first + 1L))) {
  elapsed <- proc.time()[["elapsed"]] - t0
  limit <- (wallHours - wallMargin) * 3600
  if (!is.na(lastSecs) && elapsed + 1.5 * lastSecs > limit) {
    Progress(paste("%s stopping before segment %d: %.3f h elapsed + 1.5 x %.3f h (last segment)",
                   "would pass %.3f h (CET_WALL_HOURS %g - CET_WALL_MARGIN %g); resubmit to continue"),
             format(Sys.time()), k, elapsed / 3600, lastSecs / 3600, limit / 3600, wallHours,
             wallMargin)
    break
  }
  .cet$records <- list()
  .cet$injected <- FALSE
  .cet$constraintRejections <- list()
  tSeg <- proc.time()[["elapsed"]]
  stopifnot(identical(dat$characters$tipLabel, dat$tree$tipLabel))
  if (is.null(ckp)) {
    run <- TempoMCMC(dat$tree, model, CetStart(startIndex), priors, nIter = seg, burnin = burnin,
                     thin = thin, nSlice = nSlice, moves = moves,
                     characters = dat$characters, agePriors = CetAgePriors(),
                     charPriors = CetCharPriors())
  } else {
    .cet$resume <- ckp$captured
    run <- TempoMCMC(ckp$captured$st$tree, model, ckp$captured$st$par, priors,
                     nIter = seg, burnin = 0, thin = thin, nSlice = nSlice,
                     moves = ckp$captured$moves, characters = dat$characters,
                     agePriors = CetAgePriors(), charPriors = CetCharPriors())
    stopifnot(.cet$injected)
  }
  secs <- proc.time()[["elapsed"]] - tSeg
  lastSecs <- secs
  rec <- .cet$records
  stopifnot(length(rec) == nrow(run$samples))
  smp <- run$samples
  smp$charLL <- vapply(rec, `[[`, 0, "charLL")
  smp$logPrior <- vapply(seq_len(nrow(smp)), function(i) {
    sum(vapply(names(priors), function(p) priors[[p]](rec[[i]]$state$par[[p]]), 0))
  }, 0)
  smp$logOriginPrior <- vapply(rec, function(r) CetAgePriors()$origin(
    r$state$tree$age[[r$state$tree$root]] + r$state$tree$stem), 0)
  smp$logAlphaPrior <- vapply(rec, function(r) CetCharPriors()$alpha(r$state$par$alpha), 0)
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
  # acceptance.csv is rewritten with any earlier rows of segment k replaced: a
  # kill after this write and before the checkpoint below reruns segment k on
  # resubmission, which used to append its rows a second time (R114-1).
  # Written before the checkpoint so that a kill between the two never loses
  # a segment's rows
  accFile <- file.path(outDir, "acceptance.csv")
  if (file.exists(accFile)) {
    old <- utils::read.csv(accFile)
    acc <- rbind(old[old$segment != k, names(acc)], acc)
  }
  accTmp <- paste0(accFile, ".tmp")
  utils::write.csv(acc, accTmp, row.names = FALSE)
  file.rename(accTmp, accFile)
  if (killSegment == k) {
    Progress("%s test kill after segment %d's acceptance.csv, before its checkpoint (CET_TEST_KILL_SEGMENT)",
             format(Sys.time()), k)
    quit(save = "no", status = 137)
  }
  acc <- acc[acc$segment == k, ]
  ckp <- c(list(segment = k, captured = .cet$captured, label = label, rep = rep,
                startIndex = startIndex, startTree = startTree, startTreeRel = startTreeRel,
                startTreeMd5 = startTreeMd5, segSecs = secs), target)
  tmp <- paste0(ckpFile, ".tmp")
  saveRDS(ckp, tmp)
  file.rename(tmp, ckpFile)
  Progress("%s segment %d/%d done: %.1f s (%.3f s/iteration incl. %d burn-in); logPost %.2f, charLL %.2f",
           format(Sys.time()), k, nSeg, secs, secs / (seg + acc$burnin[1]), acc$burnin[1],
           utils::tail(smp$logPost, 1), utils::tail(smp$charLL, 1))
  if ((proc.time()[["elapsed"]] - t0) / 3600 > hours && k < nSeg) {
    Progress("%s stopping after %.1f h (CET_HOURS); resubmit to continue", format(Sys.time()), hours)
    break
  }
}
if (!is.null(ckp) && ckp$segment == nSeg) Progress("%s complete", format(Sys.time()))
