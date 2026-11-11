# Starting tree for the Hymenoptera RevBayes chains.
#
# Penguin chains started from random FBD trees stalled with crown fossils outside
# a too-young crown, so every chain here starts from a tree built from the
# source study instead:
#   1. Backbone: tree_1 of ../raw/MrBayes_Commands.nex (the 68 extant taxa, the
#      topology of Ronquist et al. 2012's fixed-topology analyses; Orthoptera
#      sister to Paraneoptera + Holometabola).  It has no branch lengths, so node
#      ages are set to 400 * ((n - 1) / 67)^0.5 Ma for a node with n extant
#      descendants (root 400 Ma, above the 315 Ma minimum of the root prior;
#      Holometabola 394 Ma, above its 302 Ma minimum).
#   2. Fossils, oldest first, each at the midpoint of its stratum interval
#      (../hymenoptera_tips.csv).  Its closest extant hymenopteran is the one
#      with the smallest proportion of differing states over the analysed
#      characters scored in both (polymorphisms match if they share a state;
#      at least 5 shared characters); the fossil joins that taxon's lineage on
#      the first branch, walking rootwards, whose parent is older than the
#      fossil, half-way between max(fossil age, child age) and the parent age.
# The result respects the holometabola_withFossils constraint (all fossils are
# Hymenoptera) and has every fossil inside its age interval.  The MCMC moves
# fossils freely from there.
#
# Run from this directory after prepare-data.R: Rscript start-tree.R  (needs ape)

# Shared NEXUS helpers (identical copies in hymenoptera/revbayes and trilobites/revbayes)
ReadMatrix <- function(file) {
  lines <- readLines(file, warn = FALSE)
  start <- grep("^\\s*MATRIX\\s*$", lines, ignore.case = TRUE)[1]
  end <- start + which(grepl("^\\s*;\\s*$", lines[-seq_len(start)]))[1]
  rows <- trimws(lines[(start + 1):(end - 1)])
  rows <- rows[nzchar(rows)]
  taxa <- sub("\\s.*$", "", rows)
  seqs <- gsub("\\s", "", sub("^\\S+\\s+", "", rows))
  # one cell per character; (..) and {..} polymorphisms/ambiguities kept as one cell
  cells <- regmatches(seqs, gregexpr("\\{[^}]*\\}|\\([^)]*\\)|[^{}()]", seqs))
  # RevBayes reads both as ambiguity; write them uniformly as {..} without commas
  cells <- lapply(cells, function(x) {
    poly <- grepl("^[({]", x)
    x[poly] <- paste0("{", gsub("[^0-9A-Za-z]", "", x[poly]), "}")
    x
  })
  names(cells) <- taxa
  cells
}

ReadCharsets <- function(file) {
  lines <- readLines(file, warn = FALSE)
  cs <- grep("^\\s*CHARSET", lines, value = TRUE, ignore.case = TRUE)
  name <- sub("^\\s*CHARSET\\s+(\\S+)\\s*=.*$", "\\1", cs, ignore.case = TRUE)
  spec <- sub("^.*=\\s*(.*);\\s*$", "\\1", cs)
  idx <- lapply(strsplit(trimws(spec), "\\s+"), function(tok) {
    unlist(lapply(tok, function(t) {
      ab <- as.integer(strsplit(t, "-")[[1]])
      if (length(ab) == 1) ab else seq(ab[1], ab[2])
    }))
  })
  setNames(idx, name)
}

WriteNexus <- function(cells, cols, file, datatype, symbols = NULL, comment) {
  taxa <- names(cells)
  width <- max(nchar(taxa)) + 2
  rows <- vapply(taxa, function(tx) {
    paste0("  ", formatC(tx, width = -width), paste(cells[[tx]][cols], collapse = ""))
  }, character(1))
  fmt <- if (is.null(symbols)) {
    sprintf("  FORMAT DATATYPE=%s MISSING=? GAP=-;", datatype)
  } else {
    sprintf("  FORMAT DATATYPE=%s MISSING=? GAP=- SYMBOLS=\"%s\";", datatype, symbols)
  }
  writeLines(c(
    "#NEXUS", sprintf("[%s]", comment), "BEGIN DATA;",
    sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", length(taxa), length(cols)),
    fmt, "  MATRIX", rows, "  ;", "END;"
  ), file)
}

# States used by one column (digits and letters; ? and - ignored)
ColumnStates <- function(cells, col) {
  s <- unlist(strsplit(gsub("[^0-9A-Za-z]", "", vapply(cells, `[`, "", col)), ""))
  sort(unique(s))
}
Symbols <- c(0:9, LETTERS)

library(ape)

src <- readLines("../raw/MrBayes_Commands.nex", warn = FALSE)
tl <- grep("^\\s*tree tree_1", src)
tr0 <- grep("^\\s*translate", src, ignore.case = TRUE)
tlines <- src[(tr0 + 1):(tl - 1)]
tlines <- trimws(sub("[,;]\\s*$", "", tlines))
tlines <- tlines[nzchar(tlines)]
translate <- setNames(sub("^\\d+\\s+", "", tlines), sub("\\s.*$", "", tlines))
newick <- sub("^.*=\\s*(\\[&R\\]\\s*)?", "", src[tl])
backbone <- read.tree(text = newick)
backbone$tip.label <- unname(translate[backbone$tip.label])

tips <- read.csv("../hymenoptera_tips.csv")
extant <- tips$taxon[tips$status == "extant"]
stopifnot(setequal(backbone$tip.label, extant), is.rooted(backbone), is.binary(backbone))

# Node table: parent pointer and age per node; tips first
nTip <- Ntip(backbone)
nNode <- backbone$Nnode
label <- c(backbone$tip.label, rep(NA_character_, nNode))
parent <- integer(nTip + nNode)
parent[backbone$edge[, 2]] <- backbone$edge[, 1]
parent[nTip + 1] <- 0L
nDesc <- vapply(seq_len(nTip + nNode), function(v) {
  if (v <= nTip) 1L else length(extract.clade(backbone, v)$tip.label)
}, 0L)
age <- ifelse(seq_along(nDesc) <= nTip, 0, 400 * ((nDesc - 1) / (nTip - 1))^0.5)

morph <- ReadMatrix("../hymenoptera_morph.nex")
chars <- read.csv("../hymenoptera_characters.csv")
cols <- chars$char[!chars$excluded]
StateSet <- function(x) strsplit(gsub("[^0-9]", "", x), "")[[1]]
Distance <- function(a, b) {
  x <- morph[[a]][cols]
  y <- morph[[b]][cols]
  ok <- !(x %in% c("?", "-")) & !(y %in% c("?", "-"))
  if (sum(ok) < 5) return(Inf)
  differ <- mapply(function(u, v) length(intersect(StateSet(u), StateSet(v))) == 0, x[ok], y[ok])
  mean(differ)
}
hymExtant <- tips$taxon[tips$status == "extant" & tips$group == "Hymenoptera"]

fossils <- tips[tips$status == "fossil", ]
fossils$start_age <- (fossils$min_age + fossils$max_age) / 2
fossils <- fossils[order(-fossils$start_age), ]
placement <- data.frame()
for (i in seq_len(nrow(fossils))) {
  f <- fossils$taxon[i]
  a <- fossils$start_age[i]
  d <- vapply(hymExtant, function(x) Distance(f, x), 0)
  nearest <- names(which.min(d))
  v <- match(nearest, label)
  while (age[parent[v]] <= a) v <- parent[v]
  p <- parent[v]
  joinAge <- (max(a, age[v]) + age[p]) / 2
  # new internal node n between p and v; new tip f below n
  n <- length(age) + 1
  label <- c(label, NA, f)
  age <- c(age, joinAge, a)
  parent <- c(parent, p, n)
  parent[v] <- n
  placement <- rbind(placement, data.frame(fossil = f, start_age = a, nearest = nearest,
                                           distance = round(min(d), 3), join_age = round(joinAge, 2)))
}

# Newick from the parent table
children <- split(seq_along(parent), parent)
Newick <- function(v) {
  kids <- children[[as.character(v)]]
  me <- if (is.null(kids)) label[v] else paste0("(", paste(vapply(kids, function(k) {
    paste0(Newick(k), ":", format(age[v] - age[k], digits = 12))
  }, ""), collapse = ","), ")")
  me
}
root <- which(parent == 0)
tr <- read.tree(text = paste0(Newick(root), ";"))
stopifnot(Ntip(tr) == 113, is.binary(tr), all(tr$edge.length > 0))
holo <- readLines("data/clade_holometabola_withFossils.txt")
stopifnot(is.monophyletic(tr, holo))
# tip ages within their intervals
depth <- node.depth.edgelength(tr)
tipAge <- max(depth) - depth[seq_len(Ntip(tr))]
m <- match(tr$tip.label, tips$taxon)
stopifnot(all(tipAge >= tips$min_age[m] - 1e-6), all(tipAge <= tips$max_age[m] + 1e-6))
write.tree(tr, "data/start_tree.nwk", digits = 10)
write.csv(placement, "data/start_tree_placement.csv", row.names = FALSE)
crown <- getMRCA(tr, readLines("data/clade_hymenoptera_crown.txt"))
cat(sprintf("start tree: root %.1f Ma, crown Hymenoptera %.1f Ma, %d fossils inside the crown\n",
            max(depth), max(depth) - depth[crown],
            sum(fossils$taxon %in% tr$tip.label[unlist(prop.part(tr)[crown - Ntip(tr)])])))
