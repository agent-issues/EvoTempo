# Write RevBayes input files for the Hymenoptera total-evidence baselines.
#
# Partitions follow Ronquist et al. 2012 / Zhang et al. 2016 (MrBayes block of
# ../raw/MrBayes_Commands.nex, partition "without_CO1_3"):
#   morphology  343 characters (10 excluded), Mkv + Gamma; the 44 TYPESET ORD
#               characters (43 analysed) ordered.  Rev cannot give a subset of a
#               STANDARD matrix its own state space, so each (ordering, K)
#               class is written as its own matrix with SYMBOLS 0..K-1
#               (K = number of observed states; every analysed character uses
#               the contiguous states 0..K-1, checked below);
#   DNA         12S+16S, 18S, 28S (SYM), CO1 positions 1-2, EF1a positions 1-2,
#               EF1a position 3; CO1 position 3 excluded.
# Tip ages (Zhang et al. 2016 stratum intervals) go to data/taxa.tsv.
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

dir.create("data", showWarnings = FALSE)

morph <- ReadMatrix("../hymenoptera_morph.nex")
stopifnot(length(morph) == 113, all(lengths(morph) == 353))
chars <- read.csv("../hymenoptera_characters.csv")
stopifnot(nrow(chars) == 353)
use <- chars[!chars$excluded, ]
stopifnot(nrow(use) == 343, sum(use$ordered) == 43)
use$k <- vapply(use$char, function(j) {
  s <- ColumnStates(morph, j)
  k <- length(s)
  stopifnot(identical(s, as.character(seq_len(k) - 1)), k == use$n_states_observed[use$char == j])
  k
}, 0L)
use$part <- paste0(ifelse(use$ordered, "o", "u"), use$k)
parts <- unique(use[order(use$ordered, use$k), c("part", "ordered", "k")])
for (i in seq_len(nrow(parts))) {
  cols <- use$char[use$part == parts$part[i]]
  WriteNexus(morph, cols, sprintf("data/morph_%s.nex", parts$part[i]), "STANDARD",
             symbols = paste(Symbols[seq_len(parts$k[i])], collapse = ""),
             comment = sprintf("Hymenoptera morphology, %s K = %d: characters %s of ../hymenoptera_morph.nex",
                               ifelse(parts$ordered[i], "ordered", "unordered"), parts$k[i],
                               paste(cols, collapse = " ")))
}
parts$n_char <- as.vector(table(use$part)[parts$part])
write.table(parts[, c("part", "k", "ordered", "n_char")], "data/morph_parts.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

dna <- ReadMatrix("../hymenoptera_dna.nex")
stopifnot(identical(names(dna), names(morph)), all(lengths(dna) == 5096))
cs <- ReadCharsets("../hymenoptera_dna.nex")
dnaParts <- list(
  p12S16S = c(cs[["12S"]], cs[["16S"]]),
  p18S = cs[["18S"]],
  p28S = cs[["28S"]],
  pCO1_12 = cs[["CO1_12"]],
  pEf1a_12 = c(cs[["Ef1aF1_12"]], cs[["Ef1aF2_12"]]),
  pEf1a_3 = c(cs[["Ef1aF1_3"]], cs[["Ef1aF2_3"]])
)
stopifnot(!anyDuplicated(unlist(dnaParts)),
          length(unlist(dnaParts)) + length(cs[["CO1_3"]]) == 5096)
for (p in names(dnaParts)) {
  cols <- sort(dnaParts[[p]])
  WriteNexus(dna, cols, sprintf("data/dna_%s.nex", p), "DNA",
             comment = sprintf("Hymenoptera DNA partition %s: %d sites of ../hymenoptera_dna.nex", p, length(cols)))
}
write.table(data.frame(part = names(dnaParts), n_sites = lengths(dnaParts),
                       model = ifelse(names(dnaParts) == "p28S", "SYM", "GTR")),
            "data/dna_parts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

tips <- read.csv("../hymenoptera_tips.csv")
stopifnot(setequal(tips$taxon, names(morph)), sum(tips$status == "fossil") == 45)
write.table(data.frame(taxon = tips$taxon, min_age = tips$min_age, max_age = tips$max_age),
            "data/taxa.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Taxon sets used by the model and the summaries
cons <- read.csv("../hymenoptera_constraints.csv")
sets <- list(
  holometabola_withFossils = strsplit(cons$taxa[cons$constraint == "holometabola_withFossils"], " ")[[1]],
  hymenoptera_crown = strsplit(cons$taxa[cons$constraint == "hymenoptera"], " ")[[1]],
  hymenoptera_all = tips$taxon[tips$group == "Hymenoptera"]
)
stopifnot(length(sets$holometabola_withFossils) == 111, length(sets$hymenoptera_crown) == 60,
          length(sets$hymenoptera_all) == 105)
for (s in names(sets)) writeLines(sets[[s]], sprintf("data/clade_%s.txt", s))

# Rev definitions of the partitions and taxon sets, sourced by hymenoptera_fbd.Rev
RevVec <- function(x) {
  if (is.character(x)) x <- paste0('"', x, '"')
  if (is.logical(x)) x <- tolower(as.character(x))
  paste0("v(", paste(x, collapse = ", "), ")")
}
writeLines(c(
  "# Written by prepare-data.R; do not edit",
  paste("morph_names =", RevVec(parts$part)),
  paste("morph_k =", RevVec(parts$k)),
  paste("morph_ordered =", RevVec(parts$ordered)),
  paste("morph_nchar =", RevVec(parts$n_char)),
  paste("dna_names =", RevVec(names(dnaParts))),
  paste("dna_sites =", RevVec(unname(lengths(dnaParts)))),
  paste("dna_sym =", RevVec(names(dnaParts) == "p28S")),
  paste0("clade_holo = clade(", paste0('"', sets$holometabola_withFossils, '"', collapse = ", "), ")"),
  paste0("clade_hym_crown = clade(", paste0('"', sets$hymenoptera_crown, '"', collapse = ", "), ")"),
  paste0("clade_hym_all = clade(", paste0('"', sets$hymenoptera_all, '"', collapse = ", "), ")")
), "data/partitions.Rev")

cat("Morph partitions:", paste(parts$part, parts$n_char, collapse = "; "), "\n")
cat("DNA partitions:", paste(names(dnaParts), lengths(dnaParts), collapse = "; "), "\n")
