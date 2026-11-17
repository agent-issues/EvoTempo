# Hymenoptera CET pipeline, step 1: characters and starting tree for EvoTempo.
#
#   Rscript data-raw/empirical/hymenoptera/cet/prepare.R [rbmorph root] [burn-in]
#
# Run from the package root.  Writes, under data-raw/empirical/hymenoptera/cet/:
#   data/characters.csv     taxon x character cells as in the source matrix:
#                           a state 0..K-1, an ambiguity set such as "{01}",
#                           or "" for unknown ('?' and '-')
#   data/partitions.csv     character, source column, number of states K,
#                           ordered, partition label (u<K> / o<K>)
#   data/start-tree.nwk     starting time tree (crown: no root edge), if a
#                           RevBayes run root is given
#   data/start-tree-log.txt the selected sample and every age changed
#   data/start-tree-sample.csv the selected sample's parameters (the starting
#                           parameters of chains on this tree; common.R)
#
# Characters (identical to the morphology-only RevBayes comparator,
# ../revbayes/prepare-data.R): the 343 characters not excluded in
# ../hymenoptera_characters.csv (the MrBayes EXSET of the source analysis),
# 43 of them ordered (TYPESET ORD); K = number of observed states, which are
# the contiguous 0..K-1 for every analysed character (checked).  '?' and '-'
# are unknown.  Polymorphic and ambiguous cells ((01), {01}, ...; 22 cells in
# 14 characters) are kept as state sets, as RevBayes reads them: the
# likelihood sums over the listed states (common.R).  No character is
# dropped: every analysed character varies among its unambiguous cells
# (checked).
#
# Starting tree (README.md, "Starting tree"): from the posterior of the
# merged RevBayes `morph` baseline (#79; six chains morph_run1..6 under
# <rbmorph root>/morph_run<k>/hymenoptera_morph_run<k>.{log,trees}).  Rule:
# drop the first <burn-in> (default 0.25) of each chain's samples, as the
# baseline's summaries do, and take the sample with the highest Posterior
# among the pooled remainder.  Ties: lowest chain, then earliest iteration.
# Its tree is converted as for the penguin pipeline (rb-trees.R,
# ConvertRbTree()).  This is the reference start tree; production chains start
# from random posterior draws (draw-start-trees.R).
args <- commandArgs(trailingOnly = TRUE)
here <- "data-raw/empirical/hymenoptera"
cet <- file.path(here, "cet")
dir.create(file.path(cet, "data"), showWarnings = FALSE)

# ---- Morphology -------------------------------------------------------------
nex <- readLines(file.path(here, "hymenoptera_morph.nex"), warn = FALSE)
start <- grep("^\\s*MATRIX\\s*$", nex, ignore.case = TRUE)[1]
end <- start + which(grepl("^\\s*;\\s*$", nex[-seq_len(start)]))[1]
rows <- trimws(nex[(start + 1):(end - 1)])
rows <- rows[nzchar(rows)]
taxa <- sub("\\s.*$", "", rows)
seqs <- gsub("\\s", "", sub("^\\S+\\s+", "", rows))
cells <- regmatches(seqs, gregexpr("\\{[^}]*\\}|\\([^)]*\\)|[^{}()]", seqs))
stopifnot(length(taxa) == 113L, all(lengths(cells) == 353L))
raw <- do.call(rbind, cells)
dimnames(raw) <- list(taxa, NULL)
info <- read.csv(file.path(here, "hymenoptera_characters.csv"))
stopifnot(nrow(info) == 353L, identical(info$char, seq_len(353)))
use <- info[!info$excluded, ]
stopifnot(nrow(use) == 343L, sum(use$ordered) == 43L)
x <- raw[, use$char, drop = FALSE]
x[x %in% c("?", "-")] <- ""
amb <- matrix(grepl("^[({]", x), nrow(x))
# Ambiguity sets written uniformly as {..}, digits only
x[amb] <- paste0("{", gsub("[^0-9]", "", x[amb]), "}")
k <- vapply(seq_len(ncol(x)), function(j) {
  s <- sort(unique(unlist(strsplit(gsub("[^0-9]", "", x[, j]), ""))))
  K <- length(s)
  stopifnot(identical(s, as.character(seq_len(K) - 1L)))
  K
}, 0L)
stopifnot(identical(k, as.integer(use$n_states_observed)))
# Mkv needs every character to vary among its tips; check on the unambiguous cells
varies <- vapply(seq_len(ncol(x)), function(j) {
  v <- x[!amb[, j] & nzchar(x[, j]), j]
  length(unique(v)) > 1L
}, TRUE)
stopifnot(all(varies))
cat("Characters:", ncol(x), "; ambiguous cells kept as sets:", sum(amb), "in",
    sum(colSums(amb) > 0), "characters; unknown cells:", sum(!nzchar(x)), "\n")
colnames(x) <- paste0("c", use$char)
write.csv(data.frame(taxon = rownames(x), x, check.names = FALSE),
          file.path(cet, "data", "characters.csv"), row.names = FALSE)
part <- paste0(ifelse(use$ordered, "o", "u"), k)
write.csv(data.frame(character = seq_along(k), column = use$char, nState = k,
                     ordered = use$ordered, partition = part),
          file.path(cet, "data", "partitions.csv"), row.names = FALSE)
print(table(partition = part))

tips <- read.csv(file.path(here, "hymenoptera_tips.csv"))
stopifnot(setequal(tips$taxon, rownames(x)), sum(tips$status == "fossil") == 45L)

# ---- Starting tree ------------------------------------------------------------
if (length(args) >= 1L) {
  suppressMessages(library(ape))
  source(file.path(cet, "rb-trees.R"))
  root <- args[[1]]
  burn <- if (length(args) >= 2L) as.numeric(args[[2]]) else 0.25
  post <- ReadRbPosterior(root, burn)
  cand <- do.call(rbind, lapply(split(post, post$run), function(x) x[order(-x$Posterior, x$Iteration)[1], ]))
  pick <- cand[order(-cand$Posterior, cand$run, cand$Iteration)[1], ]
  tr <- ConvertRbTree(pick$newick, tips)
  changes <- attr(tr, "changes")
  write.tree(tr, file.path(cet, "data", "start-tree.nwk"), digits = 10)
  # The sample's parameters start the chains on this tree (common.R, CetStart())
  sampleCols <- c("run", "Iteration", "Posterior", "Likelihood", "root_age", "holo_age",
                  "lambda", "mu", "psi", "clock_rate", "alpha_morph", "ucld_sd", "num_samp_anc")
  utils::write.csv(pick[, sampleCols], file.path(cet, "data", "start-tree-sample.csv"),
                   row.names = FALSE)
  writeLines(c(
    sprintf("source: %s (burn-in %.2f of each chain; highest Posterior of the pooled remainder)", root, burn),
    sprintf("candidates (best per chain): %s",
            paste(sprintf("run%d it %d Posterior %.3f (n %d)", cand$run, cand$Iteration,
                          cand$Posterior, cand$n), collapse = "; ")),
    sprintf(paste("selected: run%d iteration %d, Posterior %.3f, Likelihood %.3f, root_age %.4f,",
                  "holo_age %.4f, lambda %.6g, mu %.6g, psi %.6g, clock_rate %.6g,",
                  "alpha_morph %.6g, ucld_sd %.6g, num_samp_anc %d"),
            pick$run, pick$Iteration, pick$Posterior, pick$Likelihood, pick$root_age,
            pick$holo_age, pick$lambda, pick$mu, pick$psi, pick$clock_rate,
            pick$alpha_morph, pick$ucld_sd, as.integer(pick$num_samp_anc)),
    sprintf("start tree: root %.6f Ma, %d sampled ancestors, %d changes, largest age change %.6f Myr",
            attr(tr, "root"), attr(tr, "nSA"), length(changes), attr(tr, "maxShift")),
    changes), file.path(cet, "data", "start-tree-log.txt"))
  cat(sprintf("Start tree from run%d iteration %d: root %.3f Ma, %d sampled ancestors, %d changes (max %.4f Myr)\n",
              pick$run, pick$Iteration, attr(tr, "root"), attr(tr, "nSA"), length(changes),
              attr(tr, "maxShift")))
}
