# Checks every sampled tree of the Hymenoptera chains against the fossil age
# intervals (data/taxa.tsv): each fossil's age must lie in its interval, no
# fossil may be a sampled ancestor at the root or a child of the root, and the
# extant tips must be at age 0.  Every tree line of each file is checked,
# including burn-in.
#
#   Rscript check-fossil-ages.R <taxa.tsv> <out.csv> <tree file> [<tree file> ...]
#
# Writes one row per tree file.  Ages come from branch lengths, with the
# extant tips at 0.  Tolerance 1e-4 Myr (the tree files print 6 decimals).
suppressMessages(library(ape))
args <- commandArgs(TRUE)
taxa <- read.delim(args[1], stringsAsFactors = FALSE)
lo <- setNames(taxa$min_age, taxa$taxon)
hi <- setNames(taxa$max_age, taxa$taxon)
extant <- taxa$taxon[taxa$max_age == 0]
tol <- 1e-4
rows <- lapply(args[-(1:2)], function(f) {
  L <- readLines(f)[-1]
  L <- L[grepl(";\\s*$", L)]
  nw <- sub("^.*\t", "", L)
  nTree <- length(nw); nRange <- 0; nRootSA <- 0; nRootFossil <- 0; extantOff <- 0
  worst <- 0; maxRoot <- 0
  for (s in nw) {
    tr <- read.tree(text = gsub("\\[[^]]*\\]", "", s))
    nt <- length(tr$tip.label)
    d <- node.depth.edgelength(tr)
    tips <- tr$tip.label
    h <- max(d[1:nt][tips %in% extant]) - d[1:nt]
    names(h) <- tips
    extantOff <- max(extantOff, abs(h[tips %in% extant]))
    out <- pmax(lo[tips] - h, h - hi[tips], 0)
    if (any(out > tol)) nRange <- nRange + 1
    worst <- max(worst, out)
    maxRoot <- max(maxRoot, max(d[1:nt][tips %in% extant]))
    re <- which(tr$edge[, 1] == nt + 1)
    kids <- tr$edge[re, 2]
    fossilKid <- kids <= nt & !(tips[pmin(kids, nt)] %in% extant)
    if (any(fossilKid)) nRootFossil <- nRootFossil + 1
    if (any(fossilKid & tr$edge.length[re] < 1e-8)) nRootSA <- nRootSA + 1
  }
  data.frame(file = basename(f), trees = nTree, trees_fossil_outside_interval = nRange,
             max_excess_myr = worst, trees_fossil_root_child = nRootFossil,
             trees_root_sampled_ancestor = nRootSA, max_extant_tip_age = extantOff,
             max_root_age = maxRoot)
})
res <- do.call(rbind, rows)
write.csv(res, args[2], row.names = FALSE)
print(res)
