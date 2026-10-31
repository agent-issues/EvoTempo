# Penguin CET pipeline, step 1: data for EvoTempo and for the morphology-only
# RevBayes comparator.
#
#   Rscript data-raw/empirical/penguins/cet/prepare.R [mcc.nex] [origin_age]
#
# Run from the package root.  Writes, under data-raw/empirical/penguins/cet/:
#   data/characters.csv     taxon x character states (0-based; NA = unknown)
#   data/partitions.csv     character, original column, number of states K
#   revbayes/data/morph_K.nex, taxa.tsv   the same matrix for RevBayes
#   revbayes/data/start_tree.nwk, start_origin.txt   the same starting tree
#   data/start-tree.nwk     starting tree (with an origin root edge), if an
#                           MCC tree is given; data/start-tree-log.txt lists
#                           every age that had to be changed
#
# Recoding (identical for EvoTempo and the RevBayes comparator):
# - '?' and '-' are unknown.
# - Polymorphic / ambiguous cells ({01}, {12}, {45}; 32 cells in 21
#   characters) are unknown: EvoTempo's MkLikelihood() takes one state or NA
#   per cell.  RevBayes would sum over the listed states instead.
# - Characters left with fewer than two observed states are dropped (Mkv
#   conditions on variable characters): character 81 (all 1 bar two {01}),
#   which RevBayes also drops under coding = "variable" in the total-evidence
#   baseline.
# - Partitions keep the number of states K of penguins_morph.nex (CHARSET
#   statesK, from the outgroup-inclusive matrix).
#
# Starting tree from a RevBayes MCC tree (annotated NEXUS from mccTree()):
# comments are stripped; extant tips are put at age 0; zero-length fossil
# tips are sampled ancestors; fossil ages are clamped into their intervals
# (a sampled ancestor's node age into its interval); parents are made at
# least `eps` older than their children; a sampled ancestor whose node cannot
# lie in its interval becomes a fossil tip; the root edge runs to
# `origin_age` (default: root age + 5 Myr).
args <- commandArgs(trailingOnly = TRUE)
here <- "data-raw/empirical/penguins"
cet <- file.path(here, "cet")
dir.create(file.path(cet, "data"), showWarnings = FALSE)
dir.create(file.path(cet, "revbayes", "data"), showWarnings = FALSE, recursive = TRUE)

# ---- Morphology -------------------------------------------------------------
nex <- readLines(file.path(here, "penguins_morph.nex"), warn = FALSE)
start <- grep("^\\s*MATRIX\\s*$", nex)
end <- start + grep("^\\s*;\\s*$", nex[-seq_len(start)])[1]
rows <- trimws(nex[(start + 1):(end - 1)])
rows <- rows[nzchar(rows)]
taxa <- sub("\\s.*$", "", rows)
cells <- lapply(sub("^\\S+\\s+", "", rows), function(s) {
  regmatches(s, gregexpr("\\{[^}]*\\}|\\([^)]*\\)|[^ ]", s))[[1]]
})
nChar <- unique(lengths(cells))
stopifnot(length(nChar) == 1L, nChar == 202L, length(taxa) == 55L)
raw <- do.call(rbind, cells)
dimnames(raw) <- list(taxa, NULL)
charsets <- regmatches(nex, regexpr("CHARSET states[0-9]+ = [^;]*", nex))
k <- integer(nChar)
for (cs in charsets) {
  K <- as.integer(sub("CHARSET states([0-9]+) .*", "\\1", cs))
  for (r in strsplit(sub(".*= ", "", cs), " ")[[1]]) {
    ends <- as.integer(strsplit(r, "-")[[1]])
    k[seq(ends[1], ends[length(ends)])] <- K
  }
}
stopifnot(all(k >= 2))
ambiguous <- matrix(grepl("^[{(]", raw), nrow(raw))
states <- matrix(suppressWarnings(as.integer(raw)), nrow(raw), dimnames = dimnames(raw))
states[ambiguous] <- NA
stopifnot(all(is.na(states) | (states >= 0 & states < rep(k, each = nrow(states)))))
nObs <- apply(states, 2, function(x) length(unique(x[!is.na(x)])))
keep <- which(nObs >= 2)
cat("Ambiguous cells set to unknown:", sum(ambiguous), "in",
    sum(colSums(ambiguous) > 0), "characters\n")
cat("Dropped (fewer than two observed states):", setdiff(seq_len(nChar), keep), "\n")
states <- states[, keep]
k <- k[keep]
colnames(states) <- paste0("c", keep)
write.csv(data.frame(taxon = rownames(states), states, check.names = FALSE),
          file.path(cet, "data", "characters.csv"), row.names = FALSE, na = "")
write.csv(data.frame(character = seq_along(keep), column = keep, nState = k),
          file.path(cet, "data", "partitions.csv"), row.names = FALSE)
print(table(nState = k))

# RevBayes: one NEXUS per K, as the total-evidence baseline
for (K in sort(unique(k))) {
  cols <- which(k == K)
  m <- states[, cols, drop = FALSE]
  m[is.na(m)] <- "?"
  out <- c("#NEXUS",
           paste0("[Penguin morphology for the CET comparator, Mkv partition K = ", K,
                  ": columns ", paste(keep[cols], collapse = " "),
                  " of ../../../penguins_morph.nex; ambiguous cells as missing]"),
           "BEGIN DATA;",
           sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", nrow(m), ncol(m)),
           sprintf("  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"%s\";",
                   paste(seq(0, K - 1), collapse = "")),
           "  MATRIX",
           sprintf("  %-30s%s", rownames(m), apply(m, 1, paste, collapse = "")),
           "  ;", "END;")
  writeLines(out, file.path(cet, "revbayes", "data", sprintf("morph_%d.nex", K)))
}
tips <- read.csv(file.path(here, "penguins_tips.csv"))
stopifnot(setequal(tips$taxon, rownames(states)))
write.table(tips[, c("taxon", "min_age", "max_age")],
            file.path(cet, "revbayes", "data", "taxa.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
constraints <- read.csv(file.path(here, "penguins_constraints.csv"))
stopifnot(nrow(constraints) == 0L)  # none in the source analysis

# ---- Starting tree ------------------------------------------------------------
if (length(args) >= 1L) {
  suppressMessages(library(ape))
  eps <- 0.01
  txt <- paste(readLines(args[[1]], warn = FALSE), collapse = "")
  nwk <- sub(".*?tree [^=]*=\\s*", "", txt)
  nwk <- gsub("\\[[^]]*\\]", "", nwk)
  nwk <- sub(";.*$", ";", nwk)
  tr <- read.tree(text = nwk)
  nTip <- Ntip(tr)
  nAll <- nTip + tr$Nnode
  parent <- integer(nAll)
  parent[tr$edge[, 2]] <- tr$edge[, 1]
  len <- numeric(nAll)
  len[tr$edge[, 2]] <- tr$edge.length
  depth <- node.depth.edgelength(tr)
  info <- tips[match(tr$tip.label, tips$taxon), ]
  stopifnot(!anyNA(info$taxon))
  extant <- info$status == "extant"
  age <- max(depth[which(extant)]) - depth
  age[which(extant)] <- 0
  sa <- !extant & len[seq_len(nTip)] < 1e-6
  changes <- character(0)
  Note <- function(...) changes <<- c(changes, sprintf(...))
  lo <- info$min_age
  hi <- info$max_age
  Clamp <- function(x, i) min(max(x, lo[i]), hi[i])
  for (i in which(!extant & !sa)) {
    new <- Clamp(age[i], i)
    if (abs(new - age[i]) > 1e-9) Note("fossil %s: %.4f -> %.4f", tr$tip.label[i], age[i], new)
    age[i] <- new
  }
  rootNode <- nTip + 1L
  # Internal nodes, children before parents (by number of edges from the root)
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
          Note("sampled ancestor %s -> fossil tip at %.4f (node %d at %.4f)",
               tr$tip.label[i], new, v, target)
          age[i] <- new
        } else {
          target <- Clamp(target, i)
          age[i] <- target
        }
      }
      if (abs(target - age[v]) > 1e-9) {
        Note("node %d: %.4f -> %.4f", v, age[v], target)
        age[v] <- target
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  rootAge <- age[rootNode]
  origin <- if (length(args) >= 2L) as.numeric(args[[2]]) else rootAge + 5
  if (!(origin > rootAge)) {
    Note("origin %.4f not above the root (%.4f): set to root + 1", origin, rootAge)
    origin <- rootAge + 1
  }
  tr$edge.length <- age[tr$edge[, 1]] - age[tr$edge[, 2]]
  tr$edge.length[sa[tr$edge[, 2]] & tr$edge[, 2] <= nTip] <- 0
  tr$root.edge <- origin - rootAge
  stopifnot(all(tr$edge.length >= 0))
  write.tree(tr, file.path(cet, "data", "start-tree.nwk"), digits = 10)
  # RevBayes: no root edge; the origin goes in separately
  rbTree <- tr
  rbTree$root.edge <- NULL
  write.tree(rbTree, file.path(cet, "revbayes", "data", "start_tree.nwk"), digits = 10)
  writeLines(format(origin, digits = 10), file.path(cet, "revbayes", "data", "start_origin.txt"))
  writeLines(c(sprintf("source: %s", args[[1]]),
               sprintf("root age %.4f, origin %.4f, sampled ancestors %d", rootAge,
                       origin, sum(sa)), changes),
             file.path(cet, "data", "start-tree-log.txt"))
  cat(sprintf("Start tree: root %.3f Ma, origin %.3f Ma, %d sampled ancestors, %d changes\n",
              rootAge, origin, sum(sa), length(changes)))
}
