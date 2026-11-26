# Metropolis-coupled chains of TempoMCMC() (#110): validation at chain lengths
# beyond tests/testthat/test-Mc3.R.
#
#   Rscript data-raw/mc3-validation.R [nIter] [nRep] [cores]
#
# From the package root, with EvoTempo installed (or MC3_LOAD_ALL=1 for
# pkgload::load_all()).  Writes data-raw/mc3-results/mc3-validation.txt.
#
# 1. Exact oracle (OracleMc3() in tests/testthat/helper-oracles.R): heats 1,
#    1/5, 1/9, `rate` free at gamma = 0; z-scores of the cold level against
#    the grid posterior, for `nRep` seeds.  All |z| should be below 4.
# 2. Mutants of the swap: the sign of the log ratio flipped, and every swap
#    accepted.  The oracle must reject both (|z| far above 4).
# 3. Tree moves: `nRep` single chains against `nRep` coupled runs (4 levels,
#    deltaHeat 2) with every tree move, gamma = 1 and free lambda, theta and
#    rate; two-sample z of the cold-level means of parameters, logPost and
#    tree statistics, with batch-means standard errors.  All |z| should be
#    below 4; swap rates are reported.
# 4. A mutant of `.TreeMove()` that ignores the heat, in the coupled runs of
#    part 3: the comparison must reject it.
args <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args) > 0) as.integer(args[[1]]) else 20000L
nRep <- if (length(args) > 1) as.integer(args[[2]]) else 8L
cores <- if (length(args) > 2) as.integer(args[[3]]) else 2L
loadAll <- nzchar(Sys.getenv("MC3_LOAD_ALL"))
Setup <- function(loadAll) {
  if (loadAll) pkgload::load_all(quiet = TRUE) else suppressMessages(library("EvoTempo"))
  env <- new.env(parent = asNamespace("EvoTempo"))
  sys.source("tests/testthat/helper-oracles.R", envir = env)
  env
}
cl <- parallel::makeCluster(cores)
parallel::clusterExport(cl, c("Setup", "loadAll", "nIter"))
invisible(parallel::clusterEvalQ(cl, env <- Setup(loadAll)))
out <- file.path("data-raw", "mc3-results", "mc3-validation.txt")
dir.create(dirname(out), showWarnings = FALSE)
Say <- function(...) cat(..., "\n", sep = "", file = out, append = TRUE)
Table <- function(x) utils::capture.output(print(round(x, 2)))
unlink(out)
Say("mc3-validation.R: nIter ", nIter, ", nRep ", nRep, "; EvoTempo at ",
    system("git rev-parse --short HEAD", intern = TRUE))

Oracle <- function(seed) unlist(env$OracleMc3(seed = seed, nIter = nIter, nChain = 1))
exact <- do.call(rbind, parallel::parLapply(cl, seq_len(nRep), Oracle))
Say("\n1. Exact oracle, cold level of heats 1, 1/5, 1/9 (rows: seeds)")
Say(paste(Table(exact), collapse = "\n"))
Say("max |z| ", round(max(abs(exact[, 1:4])), 2))

# Rebinds `fn` with `from` replaced by `to` in its source; returns the original
Mutant <- function(from, to, fn = ".SwapHeats") {
  pk <- asNamespace("EvoTempo")
  original <- get(fn, pk)
  src <- deparse(original, width.cutoff = 500L)
  stopifnot(any(grepl(from, src, fixed = TRUE)))
  utils::assignInNamespace(fn, eval(parse(text = sub(from, to, src, fixed = TRUE)), pk),
                           ns = "EvoTempo")
  original
}
mutants <- list(
  `sign flipped` = c("(logLik[[2]] - logLik[[1]])", "(logLik[[1]] - logLik[[2]])"),
  `always accepted` = c("!is.nan(logRatio) && log(u) < logRatio", "!is.nan(logRatio)")
)
Say("\n2. Mutants (seed 1)")
for (m in names(mutants)) {
  parallel::clusterExport(cl, "Mutant")
  original <- parallel::clusterCall(cl, function(ft) Mutant(ft[[1]], ft[[2]]),
                                    mutants[[m]])[[1]]
  z <- parallel::clusterCall(cl[1], Oracle, 1)[[1]]
  parallel::clusterCall(cl, function(f) {
    utils::assignInNamespace(".SwapHeats", f, ns = "EvoTempo")
  }, original)
  Say(m, ": ", paste(names(z), round(z, 2), collapse = ", "))
}

TreeRun <- function(job) {
  ranged <- TempoTree(
    structure(
      list(edge = rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9), c(9, 4), c(9, 5)),
           edge.length = c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3),
           tip.label = c("a", "b", "f", "c", "d"), Nnode = 4),
      class = "phylo"
    ),
    ranges = cbind(c(0.2, 0, 0.5, 0.5, 0), c(2.5, 0, 3, 3.5, 0))
  )
  model <- DCETModel(nBin = 11, zRange = c(-2, 2), rho = 0.7)
  priors <- c(TempoPriors()[c("lambda", "theta")],
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  moves <- c(TempoMoves(ranged, env$oracleChars, model), nodeAge = 2, fossilAge = 2,
             slide = 2, saJump = 1)
  set.seed(5000 + 10 * job[["rep"]] + job[["nChain"]])
  run <- TempoMCMC(ranged, model, env$oraclePar, priors, nIter = nIter, burnin = 500,
                   nSlice = 4, characters = env$oracleChars, moves = moves,
                   nChain = job[["nChain"]], deltaHeat = 2)
  stats <- c("lambda", "theta", "rate", "logPost", "meanX", "nAncestor", "treeLength",
             "cherries", "oldestFossil")
  list(mean = colMeans(run$samples[stats]),
       var = vapply(run$samples[stats], env$BatchVar, 0), swaps = run$swaps$rate,
       slide = run$acceptance[["slide"]])
}
jobs <- c(lapply(seq_len(nRep), function(r) c(rep = r, nChain = 1)),
          lapply(seq_len(nRep), function(r) c(rep = r, nChain = 4)))
runs <- parallel::parLapply(cl, jobs, TreeRun)
# 4. Mutant of a move: `.TreeMove()` accepts on the unheated character
# likelihood at every level, so heated levels sample the cold target while the
# swap assumes L^h.  The comparison of part 3 must reject it.
invisible(parallel::clusterCall(cl, function() {
  Mutant(
    'newChar - st[["charLL"]]',
    paste('CharLL(newTree, prop[["paths"]], st[["par"]], set, 1) -',
          'CharLL(tree, st[["paths"]], st[["par"]], set, 1)'),
    ".TreeMove"
  )
}))
mutantRuns <- parallel::parLapply(cl, jobs[nRep + seq_len(nRep)], TreeRun)
parallel::stopCluster(cl)
Pool <- function(sel, what) rowMeans(vapply(runs[sel], `[[`, runs[[sel[[1]]]][[what]], what))
single <- seq_len(nRep)
coupled <- nRep + seq_len(nRep)
se <- sqrt((Pool(single, "var") + Pool(coupled, "var")) / nRep)
tab <- cbind(single = Pool(single, "mean"), coupled = Pool(coupled, "mean"), se = se,
             z = (Pool(coupled, "mean") - Pool(single, "mean")) / se)
Say("\n3. Tree moves: single chains against the cold level of 4 coupled chains (deltaHeat 2: heats 1, 1/3, 1/5, 1/7)")
Say(paste(utils::capture.output(print(signif(tab, 4))), collapse = "\n"))
Say("max |z| ", round(max(abs(tab[, "z"])), 2), "; swap rates (cold to hot) ",
    paste(round(Pool(coupled, "swaps"), 3), collapse = ", "), "; slide acceptance ",
    round(mean(vapply(runs[single], `[[`, 0, "slide")), 3), " single, ",
    round(mean(vapply(runs[coupled], `[[`, 0, "slide")), 3), " cold level")
MutantPool <- function(what) rowMeans(vapply(mutantRuns, `[[`, mutantRuns[[1]][[what]], what))
seM <- sqrt((Pool(single, "var") + MutantPool("var")) / nRep)
zM <- (MutantPool("mean") - Pool(single, "mean")) / seM
Say("\n4. Mutant: tree moves accepted on the unheated likelihood at every level")
Say(paste(names(zM), round(zM, 2), collapse = ", "))
Say("max |z| ", round(max(abs(zM)), 2))
cat(readLines(out), sep = "\n")
