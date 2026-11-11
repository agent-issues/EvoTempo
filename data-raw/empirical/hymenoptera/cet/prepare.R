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
# Its tree is converted as for the penguin pipeline: extant tips at age 0;
# zero-length fossil tips are sampled ancestors; fossil ages (a sampled
# ancestor's node age) clamped into their intervals; parents made at least
# `eps` older than their children (the Newick is printed to 6 decimals); a
# sampled ancestor that cannot stay one (at the root, or above its interval)
# becomes a fossil tip.  The tree is a crown tree (start = "root"): no root
# edge.
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
  root <- args[[1]]
  burn <- if (length(args) >= 2L) as.numeric(args[[2]]) else 0.25
  eps <- 0.01
  cand <- do.call(rbind, lapply(1:6, function(r) {
    lg <- utils::read.delim(file.path(root, sprintf("morph_run%d", r),
                                      sprintf("hymenoptera_morph_run%d.log", r)))
    lg <- lg[-seq_len(floor(burn * nrow(lg))), ]
    data.frame(run = r, n = nrow(lg), lg[order(-lg$Posterior, lg$Iteration)[1], ])
  }))
  pick <- cand[order(-cand$Posterior, cand$run, cand$Iteration)[1], ]
  tl <- readLines(file.path(root, sprintf("morph_run%d", pick$run),
                            sprintf("hymenoptera_morph_run%d.trees", pick$run)))
  row <- tl[startsWith(tl, paste0(pick$Iteration, "\t"))]
  stopifnot(length(row) == 1L)
  tr <- read.tree(text = sub("^.*\t", "", row))
  nTip <- Ntip(tr)
  nAll <- nTip + tr$Nnode
  len <- numeric(nAll)
  len[tr$edge[, 2]] <- tr$edge.length
  depth <- node.depth.edgelength(tr)
  tinfo <- tips[match(tr$tip.label, tips$taxon), ]
  stopifnot(!anyNA(tinfo$taxon))
  extant <- tinfo$status == "extant"
  age <- max(depth[which(extant)]) - depth
  age[which(extant)] <- 0
  sa <- !extant & len[seq_len(nTip)] < 1e-6
  changes <- character(0)
  shifts <- numeric(0)
  Note <- function(old, new, ...) {
    changes <<- c(changes, sprintf(...))
    shifts <<- c(shifts, abs(new - old))
  }
  lo <- tinfo$min_age
  hi <- tinfo$max_age
  Clamp <- function(a, i) min(max(a, lo[i]), hi[i])
  for (i in which(!extant & !sa)) {
    new <- Clamp(age[i], i)
    if (abs(new - age[i]) > 1e-9) {
      Note(age[i], new, "fossil %s: %.6f -> %.6f", tr$tip.label[i], age[i], new)
    }
    age[i] <- new
  }
  rootNode <- nTip + 1L
  gen <- integer(nAll)
  pre <- reorder(tr, "cladewise")$edge
  for (i in seq_len(nrow(pre))) gen[pre[i, 2]] <- gen[pre[i, 1]] + 1L
  post <- (nTip + 1L):nAll
  post <- post[order(-gen[post])]
  for (pass in 1:50) {
    changed <- FALSE
    for (v in post) {
      kids <- tr$edge[tr$edge[, 1] == v, 2]
      saKid <- kids[kids <= nTip][sa[kids[kids <= nTip]]]
      other <- setdiff(kids, saKid)
      need <- max(age[other] + eps)
      target <- max(age[v], need)
      if (length(saKid)) {
        i <- saKid
        if (v == rootNode || need > hi[i]) {
          sa[i] <- FALSE
          new <- Clamp(min(target, hi[i]) - eps, i)
          target <- max(target, new + eps)
          Note(age[i], new, "sampled ancestor %s -> fossil tip at %.6f (node %d at %.6f)",
               tr$tip.label[i], new, v, target)
          age[i] <- new
        } else {
          target <- Clamp(target, i)
          age[i] <- target
        }
      }
      if (abs(target - age[v]) > 1e-9) {
        Note(age[v], target, "node %d: %.6f -> %.6f", v, age[v], target)
        age[v] <- target
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  tr$edge.length <- age[tr$edge[, 1]] - age[tr$edge[, 2]]
  tr$edge.length[sa[tr$edge[, 2]] & tr$edge[, 2] <= nTip] <- 0
  stopifnot(all(tr$edge.length >= 0))
  write.tree(tr, file.path(cet, "data", "start-tree.nwk"), digits = 10)
  maxShift <- if (length(shifts)) max(shifts) else 0
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
            age[rootNode], sum(sa), length(changes), maxShift),
    changes), file.path(cet, "data", "start-tree-log.txt"))
  cat(sprintf("Start tree from run%d iteration %d: root %.3f Ma, %d sampled ancestors, %d changes (max %.4f Myr)\n",
              pick$run, pick$Iteration, age[rootNode], sum(sa), length(changes), maxShift))
}
