# Penguin CET pipeline: choose the starting tree from the RevBayes
# total-evidence posterior.
#
#   Rscript data-raw/empirical/penguins/cet/select-start-tree.R <run root> <out dir> \
#     [burn-in fraction] [clock]
#
# <run root> holds one directory per RevBayes chain, <clock>_run<k>/output/
# with penguins_<clock>_run<k>.trees and .log (the layout of
# ../revbayes/run-hamilton.sh).  Defaults: burn-in 0.25, clock "shared".
#
# The post-burn-in samples (iteration >= burn-in x the chain's last iteration)
# of all chains are pooled.  Their quartet consensus (ConsTree::Quartet(),
# default settings: majority-rule start, best-improvement greedy search;
# at most 100 tips) is computed, and the sample whose topology is closest to
# it by the symmetric quartet distance 2d + r1 + r2 (the objective that
# Quartet() minimises; Quartet::QuartetStatus() counts) is the starting tree.
# Ties are broken by chain order, then iteration.  No MCC tree is used.
#
# Writes to <out dir>:
#   start-tree-sample.nex     the selected sample's tree as RevBayes wrote it
#                             (input to prepare.R)
#   start-tree-origin.txt     its origin age (origin_time at that iteration)
#   start-tree-selection.txt  chain, iteration, distances, ties, settings
#   quartet-consensus.nwk     the consensus (unrooted, possibly unresolved)
#   quartet-distances.csv     distance of every pooled sample to the consensus
#
# Packages: ape, TreeTools, ConsTree, Quartet.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
root <- args[[1]]
outDir <- args[[2]]
burnin <- if (length(args) >= 3L) as.numeric(args[[3]]) else 0.25
clock <- if (length(args) >= 4L) args[[4]] else "shared"
suppressMessages({
  library(ape)
  library(TreeTools)
  library(ConsTree)
  library(Quartet)
})
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

dirs <- sort(Sys.glob(file.path(root, paste0(clock, "_run*"))))
stopifnot(length(dirs) > 0L)
ReadChain <- function(d) {
  stem <- file.path(d, "output", paste0("penguins_", basename(d)))
  tr <- utils::read.delim(paste0(stem, ".trees"), quote = "", comment.char = "",
                          stringsAsFactors = FALSE)
  lg <- utils::read.delim(paste0(stem, ".log"), quote = "", comment.char = "",
                          check.names = FALSE)
  n <- min(nrow(tr), nrow(lg))
  tr <- tr[seq_len(n), ]
  lg <- lg[seq_len(n), ]
  stopifnot(identical(as.numeric(tr$Iteration), as.numeric(lg$Iteration)))
  keep <- tr$Iteration >= burnin * max(tr$Iteration)
  data.frame(chain = basename(d), iteration = tr$Iteration[keep],
             posterior = lg$Posterior[keep], origin = lg$origin_time[keep],
             crown = lg$crown_age[keep], nSA = lg$num_samp_anc[keep],
             newick = tr$fbd_tree[keep], stringsAsFactors = FALSE)
}
pool <- do.call(rbind, lapply(dirs, ReadChain))
cat(sprintf("%d chains, %d post-burn-in samples (%s)\n", length(dirs), nrow(pool),
            paste(sprintf("%s %d", names(table(pool$chain)), table(pool$chain)),
                  collapse = ", ")))

Strip <- function(x) gsub("\\[[^]]*\\]", "", x)
trees <- read.tree(text = Strip(pool$newick))
if (inherits(trees, "phylo")) trees <- structure(list(trees), class = "multiPhylo")
stopifnot(length(trees) == nrow(pool))
labels <- sort(trees[[1]]$tip.label)
stopifnot(all(vapply(trees, function(t) setequal(t$tip.label, labels), TRUE)))
trees <- structure(lapply(trees, RenumberTips, labels), class = "multiPhylo")
nTip <- length(labels)
cat(sprintf("%d tips\n", nTip))

t0 <- proc.time()[["elapsed"]]
cons <- ConsTree::Quartet(trees)
tCons <- proc.time()[["elapsed"]] - t0
cons <- RenumberTips(cons, labels)
cat(sprintf("Quartet consensus: %d splits (of %d possible), %.1f s\n",
            NSplits(cons), nTip - 3L, tCons))

Sym <- function(st) 2 * st[, "d"] + st[, "r1"] + st[, "r2"]
st <- QuartetStatus(trees, cf = cons)
pool$distance <- Sym(st)
pool$Q <- st[, "Q"]
# Checks: the consensus is at distance 0 from itself; tip order is irrelevant
stopifnot(Sym(QuartetStatus(structure(list(cons), class = "multiPhylo"), cf = cons)) == 0)
perm <- RenumberTips(trees[[1]], rev(labels))
stopifnot(Sym(QuartetStatus(structure(list(perm), class = "multiPhylo"), cf = cons)) ==
          pool$distance[1])

best <- order(pool$distance, match(pool$chain, basename(dirs)), pool$iteration)[1]
ties <- sum(pool$distance == pool$distance[best])
sel <- pool[best, ]
# Mean distance from the selected sample to all pooled samples, for scale
stSel <- QuartetStatus(trees, cf = trees[[best]])
meanToPool <- mean(Sym(stSel))

writeLines(c("#NEXUS", "begin trees;",
             sprintf("  tree %s_it%d = %s", sel$chain, as.integer(sel$iteration),
                     sub(";?\\s*$", ";", sel$newick)),
             "end;"),
           file.path(outDir, "start-tree-sample.nex"))
writeLines(format(sel$origin, digits = 12), file.path(outDir, "start-tree-origin.txt"))
write.tree(cons, file.path(outDir, "quartet-consensus.nwk"))
utils::write.csv(pool[, c("chain", "iteration", "posterior", "origin", "crown", "nSA",
                          "distance")],
                 file.path(outDir, "quartet-distances.csv"), row.names = FALSE)
q <- stats::quantile(pool$distance, c(0, 0.05, 0.5, 0.95, 1))
info <- c(
  sprintf("run root: %s; clock: %s; burn-in: iteration >= %.2f x last iteration", root, clock,
          burnin),
  sprintf("pooled samples: %d from %d chains (%s)", nrow(pool), length(dirs),
          paste(sprintf("%s %d", names(table(pool$chain)), table(pool$chain)), collapse = ", ")),
  sprintf("tips: %d; quartets per tree Q = %d", nTip, as.integer(pool$Q[1])),
  sprintf("quartet consensus: ConsTree %s Quartet(init = \"majority\", greedy = \"best\"); %d of %d splits; %.1f s",
          format(utils::packageVersion("ConsTree")), NSplits(cons), nTip - 3L, tCons),
  sprintf("distance: symmetric quartet distance 2d + r1 + r2 (Quartet %s QuartetStatus)",
          format(utils::packageVersion("Quartet"))),
  sprintf("distance to the consensus over the pool: min %g, 5%% %g, median %g, 95%% %g, max %g",
          q[1], q[2], q[3], q[4], q[5]),
  sprintf("sum of distances (consensus objective): %g; mean %.1f", sum(pool$distance),
          mean(pool$distance)),
  sprintf("selected: %s iteration %d; distance %g (%.5f of 2Q); ties at the minimum: %d",
          sel$chain, as.integer(sel$iteration), sel$distance, sel$distance / (2 * sel$Q), ties),
  sprintf("selected sample: posterior %.4f, origin %.6f, crown %.6f, sampled ancestors %d",
          sel$posterior, sel$origin, sel$crown, as.integer(sel$nSA)),
  sprintf("mean distance from the selected sample to the pool: %.1f (consensus: %.1f)",
          meanToPool, mean(pool$distance)))
writeLines(info, file.path(outDir, "start-tree-selection.txt"))
cat(info, sep = "\n")
