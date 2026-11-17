# Hymenoptera CET pipeline: dispersed starting trees, one per production
# chain (as the penguin pipeline's draw-start-trees.R; #89 review, R89-5).
#
#   Rscript data-raw/empirical/hymenoptera/cet/draw-start-trees.R <rbmorph root> \
#     [n = 16] [seed = 20261111] [burn-in fraction = 0.25]
#
# Run from the package root.  Pools the post-burn-in samples of the six
# chains of the morphology-only RevBayes comparator (#79) as prepare.R does
# (rb-trees.R, ReadRbPosterior()), then draws n of them uniformly at random
# WITHOUT replacement (set.seed(seed); sample.int).  Draw k goes to
# production chain k, with k = 4 (rep - 1) + cell (cells in the order b0g0,
# b0g1, b1g0, b1g1).  Each tree is converted as the reference tree
# (ConvertRbTree()).  The chain also starts from the sample's lambda, mu, psi,
# clock_rate and alpha_morph (common.R, CetStart()).  The highest-posterior
# tree of prepare.R stays the reference start tree (ladder rungs, power
# chains, checks); these draws make the production chains' starts dispersed,
# so that between-chain diagnostics (R-hat, split frequencies) can detect a
# missed mode.
#
# Writes to data/start-trees/:
#   start-tree-<kk>.nwk   draw k's starting tree
#   draws.csv             k, cell, rep, the comparator sample (run, Iteration,
#                         Posterior, Likelihood, root_age, holo_age, lambda, mu,
#                         psi, clock_rate, alpha_morph, ucld_sd, num_samp_anc),
#                         the converted tree's root age, sampled ancestors,
#                         number of age changes and largest change (Myr), and
#                         the md5 of start-tree-<kk>.nwk
#   draws-log.txt         the source, the pool and every age changed
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
suppressMessages(library(ape))
here <- "data-raw/empirical/hymenoptera"
cet <- file.path(here, "cet")
source(file.path(cet, "rb-trees.R"))
root <- args[[1]]
n <- if (length(args) >= 2L) as.integer(args[[2]]) else 16L
seed <- if (length(args) >= 3L) as.integer(args[[3]]) else 20261111L
burn <- if (length(args) >= 4L) as.numeric(args[[4]]) else 0.25
outDir <- file.path(cet, "data", "start-trees")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
tips <- read.csv(file.path(here, "hymenoptera_tips.csv"))
pool <- ReadRbPosterior(root, burn)
stopifnot(n <= nrow(pool))
set.seed(seed)
pick <- sample.int(nrow(pool), n)
stopifnot(!anyDuplicated(pick))
cellNames <- c("b0g0", "b0g1", "b1g0", "b1g1")
sampleCols <- c("run", "Iteration", "Posterior", "Likelihood", "root_age", "holo_age",
                "lambda", "mu", "psi", "clock_rate", "alpha_morph", "ucld_sd", "num_samp_anc")
log <- c(sprintf("source: %s; %d chains, %d pooled post-burn-in samples (burn-in %.2f of each chain)",
                 root, length(unique(pool$run)), nrow(pool), burn),
         sprintf("seed %d; %d draws without replacement (sample.int)", seed, n))
rows <- lapply(seq_len(n), function(k) {
  s <- pool[pick[k], ]
  tr <- ConvertRbTree(s$newick, tips)
  f <- file.path(outDir, sprintf("start-tree-%02d.nwk", k))
  write.tree(tr, f, digits = 10)
  ch <- attr(tr, "changes")
  log <<- c(log, sprintf("draw %d: run%d iteration %d; root %.6f Ma, %d sampled ancestors, %d changes, largest %.6f Myr",
                         k, s$run, s$Iteration, attr(tr, "root"), attr(tr, "nSA"), length(ch),
                         attr(tr, "maxShift")),
            if (length(ch)) paste0("  ", ch))
  data.frame(k = k, cell = cellNames[(k - 1L) %% 4L + 1L], rep = (k - 1L) %/% 4L + 1L,
             s[, sampleCols], startRoot = attr(tr, "root"), startSA = attr(tr, "nSA"),
             changes = length(ch), maxShift = attr(tr, "maxShift"),
             md5 = unname(tools::md5sum(f)), row.names = NULL)
})
out <- do.call(rbind, rows)
utils::write.csv(out, file.path(outDir, "draws.csv"), row.names = FALSE)
writeLines(log, file.path(outDir, "draws-log.txt"))
cat(log[1:2], sep = "\n")
print(out[, c("k", "cell", "rep", "run", "Iteration", "Posterior", "root_age", "holo_age",
              "clock_rate", "alpha_morph", "num_samp_anc", "startSA", "changes", "maxShift")],
      row.names = FALSE, digits = 6)
