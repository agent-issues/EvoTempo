# Penguin CET pipeline: dispersed starting trees, one per production chain.
#
#   Rscript data-raw/empirical/penguins/cet/draw-start-trees.R <run root> <out dir> \
#     [n = 16] [seed = 20261111] [burn-in fraction = 0.25] [clock = "shared"]
#
# Pools the post-burn-in samples of the RevBayes chains under <run root> exactly
# as select-start-tree.R does, then draws n of them uniformly at random WITHOUT
# replacement (set.seed(seed); sample.int).  Draw k goes to production chain k,
# with k = 4 (rep - 1) + cell (cells in the order b0g0, b0g1, b1g0, b1g1).
# The quartet-nearest tree from select-start-tree.R stays the documented
# reference start tree; these draws make the chains' starts dispersed, so that
# between-chain diagnostics (R-hat, split frequencies) can detect a missed mode
# (review of PR #89, R89-5).
#
# Writes to <out dir>:
#   start-tree-<kk>.nex         draw k's tree as RevBayes wrote it (input to prepare.R)
#   start-tree-<kk>-origin.txt  its origin age
#   draws.csv                   k, cell, rep, chain, iteration, posterior, origin,
#                               crown, sampled ancestors
#
# Packages: none beyond base R and utils.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
root <- args[[1]]
outDir <- args[[2]]
n <- if (length(args) >= 3L) as.integer(args[[3]]) else 16L
seed <- if (length(args) >= 4L) as.integer(args[[4]]) else 20261111L
burnin <- if (length(args) >= 5L) as.numeric(args[[5]]) else 0.25
clock <- if (length(args) >= 6L) args[[6]] else "shared"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

dirs <- sort(Sys.glob(file.path(root, paste0(clock, "_run*"))))
stopifnot(length(dirs) > 0L)
ReadChain <- function(d) {
  stem <- file.path(d, "output", paste0("penguins_", basename(d)))
  tr <- utils::read.delim(paste0(stem, ".trees"), quote = "", comment.char = "",
                          stringsAsFactors = FALSE)
  lg <- utils::read.delim(paste0(stem, ".log"), quote = "", comment.char = "",
                          check.names = FALSE)
  m <- min(nrow(tr), nrow(lg))
  tr <- tr[seq_len(m), ]
  lg <- lg[seq_len(m), ]
  stopifnot(identical(as.numeric(tr$Iteration), as.numeric(lg$Iteration)))
  keep <- tr$Iteration >= burnin * max(tr$Iteration)
  data.frame(chain = basename(d), iteration = tr$Iteration[keep],
             posterior = lg$Posterior[keep], origin = lg$origin_time[keep],
             crown = lg$crown_age[keep], nSA = lg$num_samp_anc[keep],
             newick = tr$fbd_tree[keep], stringsAsFactors = FALSE)
}
pool <- do.call(rbind, lapply(dirs, ReadChain))
stopifnot(n <= nrow(pool))
set.seed(seed)
pick <- sample.int(nrow(pool), n)
cells <- c("b0g0", "b0g1", "b1g0", "b1g1")
out <- data.frame(k = seq_len(n), cell = cells[(seq_len(n) - 1L) %% 4L + 1L],
                  rep = (seq_len(n) - 1L) %/% 4L + 1L, pool[pick, c("chain", "iteration",
                  "posterior", "origin", "crown", "nSA")], row.names = NULL)
for (k in seq_len(n)) {
  s <- pool[pick[k], ]
  kk <- sprintf("%02d", k)
  writeLines(c("#NEXUS", "begin trees;",
               sprintf("  tree %s_it%d = %s", s$chain, as.integer(s$iteration),
                       sub(";?\\s*$", ";", s$newick)),
               "end;"),
             file.path(outDir, sprintf("start-tree-%s.nex", kk)))
  writeLines(format(s$origin, digits = 12),
             file.path(outDir, sprintf("start-tree-%s-origin.txt", kk)))
}
utils::write.csv(out, file.path(outDir, "draws.csv"), row.names = FALSE)
cat(sprintf("%d chains, %d pooled post-burn-in samples; seed %d; %d draws without replacement\n",
            length(dirs), nrow(pool), seed, n))
print(out, row.names = FALSE)
stopifnot(!anyDuplicated(pick))
