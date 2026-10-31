# Write RevBayes input files for the trilobite tip-dating baseline.
#
# Mkv partitions by number of states K, as in the BEAST 1.8/BEAST2 analyses of
# Paterson et al. 2019.  K is the highest observed state + 1 (character 76 uses
# states 0-1, 3-8 and so K = 9, matching the BEAST 1.8 9-state partition);
# characters 45 (K = 5) and 47 (K = 4) are ordered.  Character 17 shows a single
# state in all 107 taxa: it is excluded, as it carries no information under
# Mkv and the BEAST2 analyses also used 106 characters.  States above 9 use the
# symbols A-F (character 3 has 16 states).
#
# Time: all 107 tips are Cambrian fossils (486-519 Ma).  The process is run to
# the end of the Cambrian sampling window, the Cambrian-Ordovician boundary at
# 485.4 Ma (GTS 2012), which becomes time 0, with no sampling there (rho = 0):
# a tip at time 0 would be an extant sample.  data/taxa.tsv holds ages in Myr
# before 485.4 Ma.
#
# Run from this directory: Rscript prepare-data.R

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

windowEnd <- 485.4
dir.create("data", showWarnings = FALSE)

morph <- ReadMatrix("../trilobites_morph.nex")
stopifnot(length(morph) == 107, all(lengths(morph) == 107))
chars <- read.csv("../trilobites_characters.csv")
stopifnot(nrow(chars) == 107, identical(chars$char[chars$ordered_in_beast], c(45L, 47L)))
chars$states <- lapply(chars$char, function(j) ColumnStates(morph, j))
chars$k <- vapply(chars$states, function(s) max(match(s, Symbols)), 0L)
chars$excluded <- lengths(chars$states) < 2
stopifnot(identical(chars$char[chars$excluded], 17L))
use <- chars[!chars$excluded, ]
use$part <- paste0(ifelse(use$ordered_in_beast, "o", "u"), use$k)
parts <- unique(use[order(use$ordered_in_beast, use$k), c("part", "ordered_in_beast", "k")])
names(parts)[2] <- "ordered"
for (i in seq_len(nrow(parts))) {
  cols <- use$char[use$part == parts$part[i]]
  WriteNexus(morph, cols, sprintf("data/morph_%s.nex", parts$part[i]), "STANDARD",
             symbols = paste(Symbols[seq_len(parts$k[i])], collapse = ""),
             comment = sprintf("Trilobite morphology, %s K = %d: characters %s of ../trilobites_morph.nex",
                               ifelse(parts$ordered[i], "ordered", "unordered"), parts$k[i],
                               paste(cols, collapse = " ")))
}
parts$n_char <- as.vector(table(use$part)[parts$part])
stopifnot(sum(parts$n_char) == 106)
write.table(parts, "data/morph_parts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

tips <- read.csv("../trilobites_tips.csv")
stopifnot(setequal(tips$taxon, names(morph)), all(tips$min_age == tips$max_age),
          all(tips$min_age > windowEnd))
write.table(data.frame(taxon = tips$taxon, min_age = sprintf("%.4f", tips$min_age - windowEnd),
                       max_age = sprintf("%.4f", tips$max_age - windowEnd)),
            "data/taxa.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cons <- read.csv("../trilobites_constraints.csv")
ingroup <- strsplit(cons$taxa[cons$constraint == "Ingroup_node"], " ")[[1]]
outgroup <- strsplit(cons$taxa[cons$constraint == "Fallotaspidoids"], " ")[[1]]
stopifnot(length(ingroup) == 104, length(outgroup) == 3, setequal(c(ingroup, outgroup), tips$taxon))
writeLines(ingroup, "data/clade_ingroup.txt")

RevVec <- function(x) {
  if (is.character(x)) x <- paste0('"', x, '"')
  if (is.logical(x)) x <- tolower(as.character(x))
  paste0("v(", paste(x, collapse = ", "), ")")
}
writeLines(c(
  "# Written by prepare-data.R; do not edit",
  sprintf("window_end = %s", windowEnd),
  sprintf("youngest_tip_age = %.4f", min(tips$min_age) - windowEnd),
  paste("morph_names =", RevVec(parts$part)),
  paste("morph_k =", RevVec(parts$k)),
  paste("morph_ordered =", RevVec(parts$ordered)),
  paste("morph_nchar =", RevVec(parts$n_char)),
  paste0("clade_ingroup = clade(", paste0('"', ingroup, '"', collapse = ", "), ")")
), "data/partitions.Rev")
cat("Morph partitions:", paste(parts$part, parts$n_char, collapse = "; "), "\n")
