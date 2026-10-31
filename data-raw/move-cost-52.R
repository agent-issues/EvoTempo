# Cost and ESS per CPU hour of the local path-move count in TempoMoves() (#52),
# against tree size.  Arms (path moves per iteration; tree moves are not used):
#   old   tree 1, node 2, tip 2, branch 2, joint 1 (minAcceptance = 0), the
#         schedule before #19;
#   L25   node, tip, branch at ceiling(0.25 m) of the m internal nodes, tips
#         with paths and branches with paths, shift 1: TempoMoves() on main;
#   sqrt  ceiling(sqrt(m)) of each, shift 1;
#   fixed node 4, tip 4, branch 8, shift 1.
# Trees: the two fixtures of dev/mixing/ (17 and 75 extant tips) and two more
# simulated the same way (cell beta = gamma = 1) with about 55 and 113.
# Model, priors and characters as data-raw/mixing-schedule.R.
# Usage, from the package root with EvoTempo installed:
#   Rscript data-raw/move-cost-52.R trees
#   Rscript data-raw/move-cost-52.R time <tree> <nChar> <nRep> [nIter]
#       interleaved CPU per iteration of every arm, with the incremental
#       pruning of #44 on and off where the build has it
#   Rscript data-raw/move-cost-52.R chain <tree> <nChar> <arm> <rep> [nIter] [burnin]
#   Rscript data-raw/move-cost-52.R task <i>
#       chain i of the Hamilton design (data-raw/move-cost-52-array.sh): trees
#       1-4 x 100 / 300 characters x 4 arms x 4 replicates = 128 tasks
#   Rscript data-raw/move-cost-52.R summary
# Environment: EVOTEMPO_COST_OUT (default data-raw/move-cost-results)
source("dev/mixing/common.R")
args <- commandArgs(trailingOnly = TRUE)
outDir <- Sys.getenv("EVOTEMPO_COST_OUT", "data-raw/move-cost-results")
dir.create(file.path(outDir, "chains"), recursive = TRUE, showWarnings = FALSE)
treeFile <- file.path(outDir, "trees.rds")
model <- Model(1, 1)
priors <- TruncPriors(Meanlog("high"), sdlog)
quantities <- c("theta", "s", "rate", "lambda", "rootBin", "clock", "logPost")
logged <- c("theta", "s", "rate", "lambda", "clock")
arms <- c("old", "L25", "sqrt", "fixed")

# Elements each local move updates: internal nodes, tips with paths, branches
Elements <- function(tt) {
  tree <- TempoTree(tt)
  nodes <- .PathNodes(tree)
  c(node = length(tree$postorder), tip = sum(nodes <= tree$nTip),
    branch = length(nodes))
}
Moves <- function(arm, tt) {
  m <- Elements(tt)
  switch(arm,
         old = c(tree = 1, node = 2, tip = 2, branch = 2, joint = 1),
         L25 = c(ceiling(0.25 * m), shift = 1),
         sqrt = c(ceiling(sqrt(m)), shift = 1),
         fixed = c(node = 4, tip = 4, branch = 8, shift = 1),
         stop("Unknown arm"))
}
Run <- function(tr, nChar, arm, nIter, burnin) {
  TempoMCMC(tr$tt, model, tr$truth, priors, nIter = nIter, burnin = burnin,
            thin = max(1L, nIter %/% 1000L), nSlice = nSlice,
            moves = Moves(arm, tr$tt),
            characters = tr$chars[, seq_len(nChar), drop = FALSE],
            minAcceptance = if (arm == "old") 0 else 0.01)
}
Cpu <- function(t0) unname(sum((proc.time() - t0)[c("user.self", "sys.self")]))
# The incremental pruning of #44 (option EvoTempo.incremental) is in this build
hasCache <- exists(".CharacterLikelihood", asNamespace("EvoTempo"))

mode <- if (length(args)) args[[1]] else stop("see header for usage")

if (mode == "trees") {
  fixtures <- readRDS("dev/mixing/results/trees.rds")
  set.seed(52052L)
  sim <- list(Simulate(model, Meanlog("high"), list(crown = 10, extant = c(50, 60))),
              Simulate(model, Meanlog("high"), list(crown = 12, extant = c(105, 120))))
  trees <- list(fixtures[[1]], sim[[1]], fixtures[[2]], sim[[2]])
  for (tr in trees) {
    cat(sprintf("%d extant, %d fossil, %d SA; elements %s\n", tr$nExtant, tr$nFossil,
                tr$nAncestor, paste(Elements(tr$tt), collapse = "/")))
  }
  saveRDS(trees, treeFile)
}

if (mode == "time") {
  treeId <- as.integer(args[[2]])
  nChar <- as.integer(args[[3]])
  nRep <- as.integer(args[[4]])
  nIter <- if (length(args) > 4) as.integer(args[[5]]) else 40L
  tr <- readRDS(treeFile)[[treeId]]
  cache <- if (hasCache) c(TRUE, FALSE) else NA
  grid <- expand.grid(arm = arms, incremental = cache, rep = seq_len(nRep),
                      stringsAsFactors = FALSE)
  set.seed(5200L + treeId)
  grid <- grid[sample.int(nrow(grid)), ]
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    options(EvoTempo.incremental = if (is.na(g$incremental)) TRUE else g$incremental)
    set.seed(g$rep)
    Run(tr, nChar, g$arm, 5L, 5L)  # warm-up: caches, adaptation not timed
    set.seed(100L + g$rep)
    t0 <- proc.time()
    Run(tr, nChar, g$arm, nIter, 0L)
    moves <- Moves(g$arm, tr$tt)
    data.frame(tree = treeId, nExtant = tr$nExtant, nChar = nChar, arm = g$arm,
               incremental = g$incremental, rep = g$rep,
               localMoves = sum(moves[c("node", "tip", "branch")]),
               secPerIter = Cpu(t0) / nIter)
  }))
  file <- file.path(outDir, sprintf("time-t%d-c%03d-%s.csv", treeId, nChar,
                                    if (hasCache) "p44" else "main"))
  write.csv(out, file, row.names = FALSE)
  print(aggregate(secPerIter ~ arm + cache + localMoves,
                  within(out, cache <- format(incremental)), median), digits = 3)
}

if (mode == "task") {
  design <- expand.grid(arm = arms, rep = 1:4, nChar = c(100L, 300L), tree = 1:4,
                        stringsAsFactors = FALSE)
  d <- design[as.integer(args[[2]]), ]
  args <- c("chain", d$tree, d$nChar, d$arm, d$rep)
  mode <- "chain"
}

if (mode == "chain") {
  treeId <- as.integer(args[[2]])
  nChar <- as.integer(args[[3]])
  arm <- args[[4]]
  rep <- as.integer(args[[5]])
  nIter <- if (length(args) > 5) as.integer(args[[6]]) else 10000L
  burnin <- if (length(args) > 6) as.integer(args[[7]]) else 2000L
  file <- file.path(outDir, "chains", sprintf("t%d-c%03d-%s-r%d.rds", treeId, nChar,
                                              arm, rep))
  if (file.exists(file)) quit(save = "no")
  tr <- readRDS(treeFile)[[treeId]]
  set.seed(5200L + 100L * treeId + 10L * rep + nChar)
  t0 <- proc.time()
  run <- Run(tr, nChar, arm, nIter, burnin)
  cpu <- Cpu(t0)
  out <- list(tree = treeId, nExtant = tr$nExtant, nChar = nChar, arm = arm, rep = rep,
              nIter = nIter, burnin = burnin, cpu = cpu, moves = Moves(arm, tr$tt),
              samples = run$samples[quantities], acceptance = run$acceptance,
              ess = vapply(run$samples[quantities], Ess, 0),
              incremental = hasCache, sha = Sys.getenv("EVOTEMPO_SHA"))
  saveRDS(out, file)
  cat(sprintf("%s: %.0f s CPU, ESS %s\n", basename(file), cpu,
              paste(names(out$ess), round(out$ess), collapse = " ")))
}

if (mode == "summary") {
  files <- list.files(file.path(outDir, "chains"), "\\.rds$", full.names = TRUE)
  res <- lapply(files, readRDS)
  chains <- do.call(rbind, lapply(res, function(r) {
    hours <- r$cpu * r$nIter / (r$nIter + r$burnin) / 3600
    data.frame(nExtant = r$nExtant, nChar = r$nChar, arm = r$arm, rep = r$rep,
               secPerIter = hours * 3600 / r$nIter, quantity = names(r$ess),
               ess = unname(r$ess), essPerHour = unname(r$ess) / hours,
               essPer1000 = unname(r$ess) / r$nIter * 1000,
               sd = vapply(names(r$ess), function(q) {
                 x <- r$samples[[q]]
                 stats::sd(if (q %in% logged) log(x) else x)
               }, 0))
  }))
  write.csv(chains, file.path(outDir, "move-cost-chains.csv"), row.names = FALSE)
  # Min over the key quantities, per chain
  key <- c("theta", "s", "rate", "lambda", "rootBin", "clock")
  sub <- chains[chains$quantity %in% key, ]
  mins <- do.call(rbind, lapply(
    split(sub, paste(sub$nExtant, sub$nChar, sub$arm, sub$rep)), function(d) {
      within(d[which.min(d$essPerHour), ], quantity <- "min")
    }))
  cells <- aggregate(cbind(essPerHour, essPer1000, secPerIter) ~ nExtant + nChar +
                       arm + quantity, rbind(chains, mins), mean)
  write.csv(cells, file.path(outDir, "move-cost-cells.csv"), row.names = FALSE)
  print(cells, digits = 3, row.names = FALSE)
}
