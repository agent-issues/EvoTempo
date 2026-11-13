# Trilobite CET pipeline: data conversion and starting trees.
#
#   Rscript data-raw/empirical/trilobites/cet/prepare.R <RevBayes output dir> \
#     [n = 16] [seed = 20261111] [burn-in fraction = 0.25]
#
# Run from the package root.  <RevBayes output dir> holds the six chains of the
# merged RevBayes baseline (../revbayes, #78): trilobites_morph_run<k>.log and
# trilobites_morph_run<k>.trees (on Hamilton they were copied to
# /nobackup/pjjg18/evotempo-trilobite-cet/rbposterior/ from the baseline's
# runs-dev/trilobites/morph_run<k>/output/).  Writes cet/data/:
#
# characters.csv  107 taxa x 106 characters (character 17 excluded, as in the
#                 RevBayes baseline).  Cells: a 0-based state, "" for unknown
#                 (? or -), or the states of a polymorphism/ambiguity joined by
#                 "/" (e.g. "0/1"), which the likelihood treats as ambiguous,
#                 as RevBayes does.
# partitions.csv  per character: number of states K (highest observed state
#                 + 1, as ../revbayes/prepare-data.R), ordered (45 and 47), the
#                 RevBayes partition name (u<K> or o<K>) and the number of
#                 ambiguous cells.
# start-tree-<kk>.nwk, start-values.csv, start-tree-log.txt
#                 one starting state per production chain, kk = 01..n.
#
# Starting-state rule (review of #91, F6; as the penguin pipeline's
# draw-start-trees.R, branch penguin-start-tree-r4, 511cd4e): pool the
# post-burn-in samples of the six RevBayes chains (each chain without its first
# `burn-in fraction` of logged samples, as ../revbayes/summarise.R), and draw n
# of them uniformly at random without replacement (set.seed(seed);
# sample.int()).  Draw k starts production chain k = 4 (rep - 1) + cell (cells
# in the order b0g0, b0g1, b1g0, b1g1), i.e. Slurm array task k of
# run-array.sh, with its tree and scalar values (lambda, mu, psi, clock_morph,
# ucld_sd, origin_time).  Random draws from the posterior, unlike the earlier
# maximum-Posterior starts (3.0-4.5 sd above their chain's mean posterior, five
# of six in the youngest 8% of root ages), let between-chain diagnostics detect
# a missed mode.  Each tree's tip ages are reset to the exact point ages (the
# Newick is printed to 6 decimals), internal ages are recomputed from the tips,
# zero-length fossil tips become sampled ancestors whose node takes the
# fossil's exact age, and the origin (origin_time) becomes the root edge.
suppressMessages(library(ape))
triDir <- "data-raw/empirical/trilobites"
cetDir <- file.path(triDir, "cet")
dataDir <- file.path(cetDir, "data")
dir.create(dataDir, showWarnings = FALSE)
windowEnd <- 485.4
args <- commandArgs(trailingOnly = TRUE)
rbDir <- if (length(args)) args[[1]] else "/nobackup/pjjg18/evotempo-trilobite-cet/rbposterior"
nDraw <- if (length(args) > 1) as.integer(args[[2]]) else 16L
seed <- if (length(args) > 2) as.integer(args[[3]]) else 20261111L
burnFrac <- if (length(args) > 3) as.numeric(args[[4]]) else 0.25
logLines <- character()
Log <- function(...) {
  logLines <<- c(logLines, sprintf(...))
  cat(sprintf(...), "\n")
}

# ---- Characters ----------------------------------------------------------------
# NEXUS reader as ../revbayes/prepare-data.R: one cell per character, with
# (..) and {..} kept as one cell
ReadMatrix <- function(file) {
  lines <- readLines(file, warn = FALSE)
  start <- grep("^\\s*MATRIX\\s*$", lines, ignore.case = TRUE)[1]
  end <- start + which(grepl("^\\s*;\\s*$", lines[-seq_len(start)]))[1]
  rows <- trimws(lines[(start + 1):(end - 1)])
  rows <- rows[nzchar(rows)]
  taxa <- sub("\\s.*$", "", rows)
  seqs <- gsub("\\s", "", sub("^\\S+\\s+", "", rows))
  cells <- regmatches(seqs, gregexpr("\\{[^}]*\\}|\\([^)]*\\)|[^{}()]", seqs))
  names(cells) <- taxa
  cells
}
Symbols <- c(0:9, LETTERS)
morph <- ReadMatrix(file.path(triDir, "trilobites_morph.nex"))
stopifnot(length(morph) == 107, all(lengths(morph) == 107))
chars <- read.csv(file.path(triDir, "trilobites_characters.csv"))
stopifnot(identical(chars$char[chars$ordered_in_beast], c(45L, 47L)))
# States used in a cell (inside polymorphisms too), as ColumnStates() there
CellStates <- function(cell) {
  s <- strsplit(gsub("[^0-9A-Za-z]", "", cell), "")[[1]]
  match(s, Symbols) - 1L
}
chars$states <- lapply(chars$char, function(j) {
  sort(unique(unlist(lapply(morph, function(x) CellStates(x[[j]])))))
})
chars$k <- vapply(chars$states, function(s) max(s) + 1L, 0L)
chars$excluded <- lengths(chars$states) < 2
stopifnot(identical(chars$char[chars$excluded], 17L))
use <- chars[!chars$excluded, ]
use$part <- paste0(ifelse(use$ordered_in_beast, "o", "u"), use$k)
Code <- function(cell) {
  if (cell %in% c("?", "-")) return("")
  paste(CellStates(cell), collapse = "/")
}
taxa <- names(morph)
cells <- t(vapply(taxa, function(tx) vapply(use$char, function(j) Code(morph[[tx]][[j]]), ""),
                  character(nrow(use))))
colnames(cells) <- paste0("c", use$char)
use$nAmbiguous <- colSums(matrix(grepl("/", cells), nrow(cells)))
stopifnot(nrow(use) == 106, sum(use$nAmbiguous) == 84, sum(use$nAmbiguous > 0) == 17)
write.csv(data.frame(taxon = taxa, cells, check.names = FALSE), file.path(dataDir, "characters.csv"),
          row.names = FALSE)
write.csv(data.frame(char = use$char, k = use$k, ordered = use$ordered_in_beast, part = use$part,
                     nAmbiguous = use$nAmbiguous),
          file.path(dataDir, "partitions.csv"), row.names = FALSE)
Log("Characters: %d (character 17 excluded); partitions %s; %d ambiguous cells in %d characters (%s)",
    nrow(use), paste(names(table(use$part)), table(use$part), sep = ":", collapse = " "),
    sum(use$nAmbiguous), sum(use$nAmbiguous > 0),
    paste(use$char[use$nAmbiguous > 0], collapse = " "))

# ---- Starting trees ------------------------------------------------------------
tips <- read.csv(file.path(triDir, "trilobites_tips.csv"))
stopifnot(all(tips$min_age == tips$max_age), all(tips$min_age > windowEnd))
cons <- read.csv(file.path(triDir, "trilobites_constraints.csv"))
ingroup <- strsplit(cons$taxa[cons$constraint == "Ingroup_node"], " ")[[1]]
stopifnot(length(ingroup) == 104)
StripAnnotations <- function(x) gsub("\\[&[^]]*\\]", "", x)
pool <- do.call(rbind, lapply(1:6, function(k) {
  lg <- read.delim(file.path(rbDir, sprintf("trilobites_morph_run%d.log", k)))
  tr <- read.delim(file.path(rbDir, sprintf("trilobites_morph_run%d.trees", k)),
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  m <- min(nrow(lg), nrow(tr))
  lg <- lg[seq_len(m), ]
  tr <- tr[seq_len(m), ]
  stopifnot(identical(as.numeric(lg$Iteration), as.numeric(tr$Iteration)))
  keep <- seq_len(m) > floor(burnFrac * m)
  data.frame(chain = sprintf("trilobites_morph_run%d", k), lg[keep, ],
             newick = tr[[ncol(tr)]][keep], stringsAsFactors = FALSE)
}))
stopifnot(nDraw <= nrow(pool))
set.seed(seed)
pick <- sample.int(nrow(pool), nDraw)
stopifnot(!anyDuplicated(pick))
Log("Starting states: %d pooled post-burn-in samples from 6 RevBayes chains (first %g of each dropped); seed %d; %d draws without replacement",
    nrow(pool), burnFrac, seed, nDraw)
cellNames <- c("b0g0", "b0g1", "b1g0", "b1g1")
starts <- list()
for (k in seq_len(nDraw)) {
  s <- pool[pick[k], ]
  tr <- read.tree(text = StripAnnotations(s$newick))
  stopifnot(setequal(tr$tip.label, tips$taxon), is.binary(tr), is.rooted(tr))
  nTip <- Ntip(tr)
  tipAge <- tips$min_age[match(tr$tip.label, tips$taxon)] - windowEnd
  depth <- node.depth.edgelength(tr)
  # Internal ages from every descendant tip (mean), to average out rounding
  desc <- prop.part(tr)
  implied <- lapply(seq_len(tr$Nnode), function(i) {
    d <- desc[[i]]
    tipAge[d] + depth[d] - depth[nTip + i]
  })
  intAge <- vapply(implied, mean, 0)
  roundErr <- max(vapply(implied, function(a) diff(range(a)), 0))
  age <- c(tipAge, intAge)
  len <- numeric(nTip + tr$Nnode)
  len[tr$edge[, 2]] <- tr$edge.length
  sa <- which(len[seq_len(nTip)] < 1e-5)
  parent <- integer(nTip + tr$Nnode)
  parent[tr$edge[, 2]] <- tr$edge[, 1]
  # A sampled ancestor's node is at the fossil's exact age
  saShift <- abs(age[parent[sa]] - tipAge[sa])
  age[parent[sa]] <- tipAge[sa]
  stopifnot(!(nTip + 1L) %in% parent[sa])
  # Parents strictly older than non-ancestor children (postorder)
  changed <- 0L
  for (e in postorder(tr)) {
    p <- tr$edge[e, 1]
    ch <- tr$edge[e, 2]
    if (ch %in% sa) next
    if (age[p] < age[ch] + 1e-6) {
      if (p %in% parent[sa]) stop("a sampled-ancestor node would have to move")
      age[p] <- age[ch] + 1e-6
      changed <- changed + 1L
    }
  }
  tr$edge.length <- age[tr$edge[, 1]] - age[tr$edge[, 2]]
  tr$edge.length[tr$edge[, 2] %in% sa] <- 0
  stem <- s$origin_time - age[nTip + 1L]
  stopifnot(stem > 0, abs(s$root_age - windowEnd - age[nTip + 1L]) < 1e-3)
  tr$root.edge <- stem
  stopifnot(is.monophyletic(tr, ingroup))
  write.tree(tr, file.path(dataDir, sprintf("start-tree-%02d.nwk", k)), digits = 10)
  starts[[k]] <- data.frame(k = k, cell = cellNames[(k - 1L) %% 4L + 1L], rep = (k - 1L) %/% 4L + 1L,
                            chain = s$chain, iteration = s$Iteration,
                            posterior = s$Posterior, likelihood = s$Likelihood,
                            lambda = s$lambda, mu = s$mu, psi = s$psi,
                            clock_morph = s$clock_morph, ucld_sd = s$ucld_sd,
                            ucld_cv = s$ucld_cv, alpha_morph = s$alpha_morph,
                            origin_time = s$origin_time, root_age = s$root_age,
                            num_samp_anc = s$num_samp_anc)
  Log(paste("Draw %02d (%s rep %d): %s iteration %d, Posterior %.3f; root %.4f Ma, origin %.4f Ma,",
            "%d sampled ancestors (%s); max tip-implied age spread %.2g Myr; max SA node shift %.2g;",
            "%d parent ages raised to child + 1e-6"),
      k, cellNames[(k - 1L) %% 4L + 1L], (k - 1L) %/% 4L + 1L, s$chain, as.integer(s$Iteration),
      s$Posterior, age[nTip + 1L] + windowEnd, s$origin_time + windowEnd, length(sa),
      paste(tr$tip.label[sa], collapse = ", "), roundErr, if (length(sa)) max(saShift) else 0,
      changed)
}
sv <- do.call(rbind, starts)
# Where the draws sit in the pooled posterior (percentiles), for the README
Pct <- function(x, ref) round(100 * vapply(x, function(v) mean(ref <= v), 0), 1)
sv$posteriorPercentile <- Pct(sv$posterior, pool$Posterior)
sv$rootPercentile <- Pct(sv$root_age, pool$root_age)
write.csv(sv, file.path(dataDir, "start-values.csv"), row.names = FALSE)
writeLines(logLines, file.path(dataDir, "start-tree-log.txt"))
