# Hymenoptera CET pipeline: reading the morphology-only RevBayes comparator's
# posterior and converting its trees into CET starting trees.  Sourced by
# prepare.R (the reference start tree) and draw-start-trees.R (one start tree
# per production chain).  Needs ape.
#
# ReadRbPosterior(root, burn): the six chains morph_run1..6 under <root>
# (<root>/morph_run<k>/hymenoptera_morph_run<k>.{log,trees}), each without its
# first `burn` (default 0.25) of samples, as the baseline's summaries; one row
# per sample with the chain (run), its log columns and the Newick tree.
ReadRbPosterior <- function(root, burn = 0.25) {
  do.call(rbind, lapply(1:6, function(r) {
    stem <- file.path(root, sprintf("morph_run%d", r), sprintf("hymenoptera_morph_run%d", r))
    lg <- utils::read.delim(paste0(stem, ".log"))
    tl <- readLines(paste0(stem, ".trees"))
    tl <- tl[grepl("^[0-9]+\t", tl)]
    it <- as.numeric(sub("\t.*$", "", tl))
    m <- min(nrow(lg), length(tl))
    lg <- lg[seq_len(m), ]
    stopifnot(identical(as.numeric(lg$Iteration), it[seq_len(m)]))
    keep <- seq_len(m) > floor(burn * m)
    data.frame(run = r, n = sum(keep), lg[keep, ], newick = sub("^.*\t", "", tl[seq_len(m)][keep]),
               stringsAsFactors = FALSE)
  }))
}
# A comparator tree as a CET starting tree (crown: no root edge; start =
# "root"): extant tips at age 0; zero-length fossil tips are sampled
# ancestors; fossil ages (a sampled ancestor's node age) clamped into their
# intervals; parents made at least `eps` older than their children (the
# Newick is printed to 10 decimals); a sampled ancestor that cannot stay one
# (at the root, or above its interval) becomes a fossil tip.  `tips` is
# ../hymenoptera_tips.csv.  Returns the ape tree with attributes "root" (root
# age), "nSA" (sampled ancestors), "changes" (one line per age changed) and
# "maxShift" (largest age change, Myr).
ConvertRbTree <- function(newick, tips, eps = 0.01) {
  row <- newick
  tr <- read.tree(text = row)
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
  structure(tr, root = age[rootNode], nSA = sum(sa), changes = changes,
            maxShift = if (length(shifts)) max(shifts) else 0)
}
