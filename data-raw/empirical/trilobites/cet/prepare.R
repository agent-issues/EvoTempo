# Trilobite CET pipeline: data conversion and starting trees.
#
#   Rscript data-raw/empirical/trilobites/cet/prepare.R <RevBayes output dir> [burn-in fraction]
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
# start-tree-<k>.nwk, start-values.csv, start-tree-log.txt
#                 one starting state per RevBayes chain k = 1..6 (rule below).
#
# Starting-state rule: in each RevBayes chain, drop the first `burn-in
# fraction` (default 0.25, as ../revbayes/summarise.R) of the logged samples
# and take the sample with the highest logged Posterior; its tree (from the
# .trees file, same iteration) and scalar values start CET replicate k.  The
# tree's tip ages are reset to the exact point ages (the Newick is printed to
# 6 decimals), internal ages are recomputed from the tips, zero-length fossil
# tips become sampled ancestors whose node takes the fossil's exact age, and
# the origin (RevBayes origin_time) becomes the root edge.
suppressMessages(library(ape))
triDir <- "data-raw/empirical/trilobites"
cetDir <- file.path(triDir, "cet")
dataDir <- file.path(cetDir, "data")
dir.create(dataDir, showWarnings = FALSE)
windowEnd <- 485.4
args <- commandArgs(trailingOnly = TRUE)
rbDir <- if (length(args)) args[[1]] else "/nobackup/pjjg18/evotempo-trilobite-cet/rbposterior"
burnFrac <- if (length(args) > 1) as.numeric(args[[2]]) else 0.25
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
starts <- list()
for (k in 1:6) {
  lg <- read.delim(file.path(rbDir, sprintf("trilobites_morph_run%d.log", k)))
  keep <- seq_len(nrow(lg)) > floor(burnFrac * nrow(lg))
  best <- which(keep)[which.max(lg$Posterior[keep])]
  it <- lg$Iteration[best]
  trLines <- readLines(file.path(rbDir, sprintf("trilobites_morph_run%d.trees", k)))
  hit <- grep(sprintf("^%d\t", it), trLines)
  stopifnot(length(hit) == 1)
  nwk <- StripAnnotations(sub("^.*\t", "", trLines[hit]))
  tr <- read.tree(text = nwk)
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
  stem <- lg$origin_time[best] - age[nTip + 1L]
  stopifnot(stem > 0, abs(lg$root_age[best] - windowEnd - age[nTip + 1L]) < 1e-3)
  tr$root.edge <- stem
  stopifnot(is.monophyletic(tr, ingroup))
  write.tree(tr, file.path(dataDir, sprintf("start-tree-%d.nwk", k)), digits = 10)
  starts[[k]] <- data.frame(rep = k, chain = sprintf("trilobites_morph_run%d", k),
                            nSamples = nrow(lg), firstKept = which(keep)[1], iteration = it,
                            posterior = lg$Posterior[best], likelihood = lg$Likelihood[best],
                            lambda = lg$lambda[best], mu = lg$mu[best], psi = lg$psi[best],
                            clock_morph = lg$clock_morph[best], ucld_sd = lg$ucld_sd[best],
                            ucld_cv = lg$ucld_cv[best], alpha_morph = lg$alpha_morph[best],
                            origin_time = lg$origin_time[best], root_age = lg$root_age[best],
                            num_samp_anc = lg$num_samp_anc[best])
  Log(paste("Chain %d: %d samples, kept from row %d; best Posterior %.3f at iteration %d;",
            "root %.4f Ma, origin %.4f Ma, %d sampled ancestors (%s); max tip-implied age spread %.2g Myr;",
            "max SA node shift %.2g; %d parent ages raised to child + 1e-6"),
      k, nrow(lg), which(keep)[1], lg$Posterior[best], it, age[nTip + 1L] + windowEnd,
      lg$origin_time[best] + windowEnd, length(sa), paste(tr$tip.label[sa], collapse = ", "),
      roundErr, if (length(sa)) max(saShift) else 0, changed)
}
write.csv(do.call(rbind, starts), file.path(dataDir, "start-values.csv"), row.names = FALSE)
writeLines(logLines, file.path(dataDir, "start-tree-log.txt"))
