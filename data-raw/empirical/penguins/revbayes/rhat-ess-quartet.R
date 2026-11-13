#!/usr/bin/env Rscript
# Posterior summaries and convergence diagnostics for a set of RevBayes
# chains: rank-normalised split R-hat and bulk/tail ESS (Vehtari et al. 2021,
# as implemented in MkPrime), and the quartet consensus of the pooled tree
# samples (ConsTree::Quartet(), default settings, as in
# ../cet/select-start-tree.R).  No MCC tree is computed.
#
#   Rscript rhat-ess-quartet.R <root> <clock> <burnin_frac> <outdir> [max_trees] [stem]
#
# Chains: <root>/<clock>_run*/output/<stem><clock>_run<k>.{log,trees}
# (stem defaults to "penguins_"; "penguins_morph_" for the morphology-only
# comparator).  The first burnin_frac of each chain's samples is discarded and
# chains are trimmed to the shortest remaining length.  The consensus uses at
# most max_trees (default 4000) pooled samples, thinned evenly.
#
# Writes to <outdir>:
#   rhat_ess_<clock>.csv          per parameter: mean, median, 2.5% and 97.5%
#                                 quantiles, R-hat, bulk and tail ESS (pooled),
#                                 minimum per-chain bulk ESS
#   quartet_consensus_<clock>.nwk quartet consensus (unrooted); node labels are
#                                 split frequencies in the pooled samples
#   split_support_<clock>.csv     every split of the consensus with its
#                                 frequency in each chain and pooled
# Packages: ape, TreeTools, ConsTree, MkPrime (internal functions).
suppressPackageStartupMessages({
  library(ape); library(TreeTools); library(ConsTree)
})
a <- commandArgs(TRUE)
root <- a[1]; clock <- a[2]; burn <- as.numeric(a[3]); outdir <- a[4]
maxTrees <- if (length(a) >= 5) as.integer(a[5]) else 4000L
stem <- if (length(a) >= 6) a[6] else "penguins_"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
runs <- sort(Sys.glob(file.path(root, paste0(clock, "_run*"))))
stems <- file.path(runs, "output", paste0(stem, basename(runs)))
stopifnot(length(runs) > 0, file.exists(paste0(stems, ".log")))
cat("Chains:", basename(runs), "\n")

Read <- function(f) {
  x <- utils::read.delim(f, quote = "", comment.char = "", check.names = FALSE,
                         stringsAsFactors = FALSE)
  x[-seq_len(floor(burn * nrow(x))), , drop = FALSE]
}
logs <- lapply(paste0(stems, ".log"), Read)
n <- min(vapply(logs, nrow, 1L))
logs <- lapply(logs, function(x) utils::tail(x, n))
cat("Samples per chain after burn-in:", n, "; iterations",
    paste(vapply(logs, function(x) paste(range(x$Iteration), collapse = "-"), ""),
          collapse = ", "), "\n")
pars <- setdiff(Reduce(intersect, lapply(logs, names)), "Iteration")
pars <- pars[vapply(pars, function(p) is.numeric(logs[[1]][[p]]), TRUE)]
diag <- do.call(rbind, lapply(pars, function(p) {
  m <- vapply(logs, function(x) x[[p]], numeric(n))
  constant <- length(unique(as.vector(m))) == 1
  q <- stats::quantile(m, c(0.025, 0.5, 0.975), names = FALSE)
  data.frame(parameter = p, mean = mean(m), median = q[2], q025 = q[1], q975 = q[3],
             rhat = if (constant) NA else MkPrime:::.Rhat(m),
             ess_bulk = if (constant) NA else MkPrime:::.EssBulk(m),
             ess_tail = if (constant) NA else MkPrime:::.EssTail(m),
             min_chain_ess_bulk = if (constant) NA else
               min(apply(m, 2, MkPrime:::.EssBulk)))
}))
write.csv(diag, file.path(outdir, paste0("rhat_ess_", clock, ".csv")), row.names = FALSE)
ok <- diag[!is.na(diag$rhat), ]
cat(sprintf("%s: max R-hat %.4f (%s); min bulk ESS %.0f (%s); min tail ESS %.0f (%s)\n",
    clock, max(ok$rhat), ok$parameter[which.max(ok$rhat)],
    min(ok$ess_bulk), ok$parameter[which.min(ok$ess_bulk)],
    min(ok$ess_tail), ok$parameter[which.min(ok$ess_tail)]))

# Trees: same samples as the logs (matched by iteration)
Strip <- function(x) gsub("\\[[^]]*\\]", "", x)
chainTrees <- lapply(seq_along(stems), function(i) {
  tr <- utils::read.delim(paste0(stems[i], ".trees"), quote = "", comment.char = "",
                          stringsAsFactors = FALSE)
  tr <- tr[match(logs[[i]]$Iteration, tr$Iteration), ]
  stopifnot(!anyNA(tr$Iteration))
  t <- read.tree(text = Strip(tr[[ncol(tr)]]))
  labels <- sort(t[[1]]$tip.label)
  structure(lapply(t, RenumberTips, labels), class = "multiPhylo")
})
labels <- chainTrees[[1]][[1]]$tip.label
pool <- do.call(c, chainTrees)
keep <- unique(round(seq(1, length(pool), length.out = min(maxTrees, length(pool)))))
t0 <- proc.time()[["elapsed"]]
cons <- RenumberTips(ConsTree::Quartet(pool[keep]), labels)
cat(sprintf("Quartet consensus of %d trees: %d splits (of %d possible), %.0f s\n",
            length(keep), NSplits(cons), length(labels) - 3L,
            proc.time()[["elapsed"]] - t0))
cs <- as.Splits(cons)
Freq <- function(trees) SplitFrequency(cs, trees) / length(trees)
perChain <- vapply(chainTrees, Freq, numeric(length(cs)))
if (!is.matrix(perChain)) perChain <- matrix(perChain, ncol = length(chainTrees))
colnames(perChain) <- basename(runs)
pooled <- Freq(pool)
sizes <- TipsInSplits(cs)
smaller <- vapply(seq_along(cs), function(i) {
  inSplit <- as.logical(cs[[i]])
  side <- if (sum(inSplit) <= length(inSplit) / 2) inSplit else !inSplit
  paste(sort(labels[side]), collapse = ",")
}, "")
write.csv(data.frame(split = smaller, n_tips = pmin(sizes, length(labels) - sizes),
                     pooled = pooled, perChain, check.names = FALSE),
          file.path(outdir, paste0("split_support_", clock, ".csv")), row.names = FALSE)
cons$node.label <- rep("", cons$Nnode)
cons$node.label[as.integer(names(cs)) - length(labels)] <- sprintf("%.3f", pooled)
write.tree(cons, file.path(outdir, paste0("quartet_consensus_", clock, ".nwk")))
cat("Consensus written; ASDSF over consensus splits:",
    sprintf("%.4f", mean(apply(perChain, 1, stats::sd))), "\n")
