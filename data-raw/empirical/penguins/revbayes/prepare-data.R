# Write RevBayes input files for the penguin total-evidence baseline.
#
# RevBayes's NEXUS reader (NCL) rejects the CHARSET ranges and CHARSTATELABELS
# in ../penguins_morph.nex, and Rev cannot give a subset of a STANDARD matrix a
# state space other than the one implied by its symbols.  So each Mkv partition
# (number of states K assumed by Gavryushkina et al. 2017, from CHARSET statesK)
# is written as its own matrix with SYMBOLS = "0..K-1", and each DNA locus as
# its own matrix.  Tip ages go to a tab-separated taxon file for readTaxonData().
#
# Run from this directory: Rscript prepare-data.R

ReadMatrix <- function(file) {
  lines <- readLines(file)
  start <- grep("^\\s*MATRIX\\s*$", lines, ignore.case = TRUE)
  end <- start + which(grepl("^\\s*;\\s*$", lines[-seq_len(start)]))[1]
  rows <- trimws(lines[(start + 1):(end - 1)])
  rows <- rows[nzchar(rows)]
  taxa <- sub("\\s.*$", "", rows)
  seqs <- sub("^\\S+\\s+", "", rows)
  # Split into cells, keeping {..} polymorphisms as one cell
  cells <- regmatches(seqs, gregexpr("\\{[^}]*\\}|[^{}]", seqs))
  names(cells) <- taxa
  cells
}

ReadCharsets <- function(file) {
  lines <- readLines(file)
  cs <- grep("^\\s*CHARSET", lines, value = TRUE, ignore.case = TRUE)
  name <- sub("^\\s*CHARSET\\s+(\\S+)\\s*=.*$", "\\1", cs, ignore.case = TRUE)
  spec <- sub("^.*=\\s*(.*);\\s*$", "\\1", cs)
  idx <- lapply(strsplit(spec, "\\s+"), function(tok) {
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

dir.create("data", showWarnings = FALSE)

morph <- ReadMatrix("../penguins_morph.nex")
stopifnot(length(morph) == 55, all(lengths(morph) == 202))
morphSets <- ReadCharsets("../penguins_morph.nex")
stopifnot(sort(unlist(morphSets)) == seq_len(202))
for (set in names(morphSets)) {
  k <- as.integer(sub("states", "", set))
  cols <- morphSets[[set]]
  used <- unique(unlist(strsplit(gsub("[^0-9]", "", unlist(lapply(morph, `[`, cols))), "")))
  stopifnot(all(as.integer(used) < k))
  WriteNexus(morph, cols, sprintf("data/morph_%d.nex", k), "STANDARD",
             symbols = paste(seq_len(k) - 1, collapse = ""),
             comment = sprintf("Penguin morphology, Mkv partition K = %d: characters %s of ../penguins_morph.nex",
                               k, paste(cols, collapse = " ")))
}

dna <- ReadMatrix("../penguins_dna.nex")
stopifnot(identical(names(dna), names(morph)), all(lengths(dna) == 8145))
dnaSets <- ReadCharsets("../penguins_dna.nex")
for (locus in names(dnaSets)) {
  WriteNexus(dna, dnaSets[[locus]], sprintf("data/dna_%s.nex", locus), "DNA",
             comment = sprintf("Penguin DNA, locus %s (sites %d-%d of ../penguins_dna.nex)", locus,
                               min(dnaSets[[locus]]), max(dnaSets[[locus]])))
}
writeLines(names(dnaSets), "data/loci.txt")

tips <- read.csv("../penguins_tips.csv")
stopifnot(setequal(tips$taxon, names(morph)))
write.table(data.frame(taxon = tips$taxon, min_age = tips$min_age, max_age = tips$max_age),
            "data/taxa.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Number of sites per DNA locus, for the relative-rate weights
nSites <- lengths(dnaSets)
writeLines(as.character(nSites), "data/dna_sites.txt")
cat("Morph partitions:", paste(names(morphSets), lengths(morphSets), collapse = "; "), "\n")
cat("DNA loci:", paste(names(nSites), nSites, collapse = "; "), "\n")
