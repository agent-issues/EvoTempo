# Rebuild the processed empirical datasets from their raw downloads.
#
# Run from the package root:  Rscript data-raw/empirical/build.R
#
# Each dataset directory holds raw/ (byte-exact downloads, listed with
# checksums in raw/MANIFEST.csv), curated/ (values transcribed from papers)
# and the processed files written here. Raw files missing from raw/ are
# downloaded; Dryad downloads need DRYAD_CLIENT_ID and DRYAD_CLIENT_SECRET.
# The script stops if a checksum fails or if matrix and tip-table taxa differ.

library("xml2")
library("readxl")

empDir <- file.path("data-raw", "empirical")

# Raw files -------------------------------------------------------------------

.DryadHandle <- function() {
  id <- Sys.getenv("DRYAD_CLIENT_ID")
  secret <- Sys.getenv("DRYAD_CLIENT_SECRET")
  if (!nzchar(id) || !nzchar(secret)) {
    stop("Set DRYAD_CLIENT_ID and DRYAD_CLIENT_SECRET to download from Dryad")
  }
  tokenHandle <- curl::new_handle()
  curl::handle_setopt(tokenHandle, postfields = paste0(
    "grant_type=client_credentials&client_id=", curl::curl_escape(id),
    "&client_secret=", curl::curl_escape(secret)
  ))
  res <- curl::curl_fetch_memory(
    "https://datadryad.org/oauth/token", handle = tokenHandle
  )
  token <- jsonlite::fromJSON(rawToChar(res[["content"]]))[["access_token"]]
  h <- curl::new_handle()
  curl::handle_setheaders(h, Authorization = paste("Bearer", token))
  # Return:
  h
}

.CheckRaw <- function(dataset) {
  rawDir <- file.path(empDir, dataset, "raw")
  manifest <- read.csv(file.path(rawDir, "MANIFEST.csv"))
  manifest <- manifest[as.logical(manifest[["committed"]]), ]
  paths <- file.path(rawDir, manifest[["file"]])
  dryadHandle <- NULL
  for (i in which(!file.exists(paths))) {
    url <- manifest[["url"]][i]
    if (grepl("datadryad.org", url, fixed = TRUE)) {
      if (is.null(dryadHandle)) dryadHandle <- .DryadHandle()
      curl::curl_download(url, paths[i], handle = dryadHandle)
    } else {
      curl::curl_download(url, paths[i])
    }
  }
  bad <- unname(tools::md5sum(paths)) != manifest[["md5"]]
  if (any(bad)) {
    stop("Checksum mismatch: ", paste(manifest[["file"]][bad], collapse = ", "))
  }
  # Return:
  invisible(paths)
}

# NEXUS helpers ----------------------------------------------------------------

# Split a character string into NEXUS tokens; polymorphic and uncertain cells
# are written compactly, e.g. "(0,1)" becomes "(01)".
.Tokens <- function(x) {
  tokens <- regmatches(x, gregexpr("\\{[^}]*\\}|\\([^)]*\\)|\\S", x))[[1]]
  # Return:
  gsub("[, ]", "", tokens)
}

.States <- function(token) {
  # Return:
  setdiff(strsplit(gsub("[{}()]", "", token), "")[[1]], c("?", "-"))
}

.ObservedStates <- function(column) {
  # Return:
  length(unique(unlist(lapply(column, .States))))
}

.Ranges <- function(x) {
  x <- sort(unique(x))
  breaks <- c(0, which(diff(x) != 1), length(x))
  # Return:
  paste(vapply(seq_len(length(breaks) - 1), function(i) {
    from <- x[breaks[i] + 1]
    to <- x[breaks[i + 1]]
    if (from == to) as.character(from) else paste0(from, "-", to)
  }, ""), collapse = " ")
}

.Quote <- function(x) {
  ifelse(grepl("^[A-Za-z0-9_.-]+$", x), x, paste0("'", gsub("'", "''", x), "'"))
}

.WriteNexus <- function(file, mat, dataType, symbols = NULL,
                        charLabels = NULL, stateLabels = NULL,
                        charSets = NULL,
                        ordered = NULL, excluded = NULL, comment = NULL) {
  taxa <- rownames(mat)
  width <- max(nchar(.Quote(taxa)))
  rows <- paste0(
    "  ", formatC(.Quote(taxa), width = -width), "  ",
    apply(mat, 1, paste, collapse = "")
  )
  format <- paste0(
    "  FORMAT DATATYPE=", dataType, " MISSING=? GAP=-",
    if (!is.null(symbols)) paste0(" SYMBOLS=\"", symbols, "\""), ";"
  )
  labels <- if (!is.null(stateLabels)) {
    states <- vapply(stateLabels, function(s) {
      paste(.Quote(s), collapse = " ")
    }, "")
    c("  CHARSTATELABELS", paste0(
      "    ", seq_along(charLabels), " ", .Quote(charLabels), " / ", states,
      c(rep(",", length(charLabels) - 1), "")
    ), "  ;")
  } else if (!is.null(charLabels)) {
    c("  CHARLABELS", paste0("    ", .Quote(charLabels)), "  ;")
  }
  sets <- if (length(charSets)) {
    c("", "BEGIN SETS;", paste0(
      "  CHARSET ", names(charSets), " = ",
      vapply(charSets, .Ranges, ""), ";"
    ), "END;")
  }
  assumptions <- if (length(ordered) || length(excluded)) {
    c("", "BEGIN ASSUMPTIONS;", "  OPTIONS DEFTYPE=UNORD;",
      if (length(ordered)) {
        paste0("  TYPESET * published = ORD: ", .Ranges(ordered), ";")
      },
      if (length(excluded)) {
        paste0("  EXSET * published = ", .Ranges(excluded), ";")
      },
      "END;")
  }
  writeLines(c(
    "#NEXUS",
    if (!is.null(comment)) paste0("[", comment, "]"),
    "", "BEGIN TAXA;",
    paste0("  DIMENSIONS NTAX=", nrow(mat), ";"),
    "  TAXLABELS", paste0("    ", .Quote(taxa)), "  ;", "END;",
    "", "BEGIN CHARACTERS;",
    paste0("  DIMENSIONS NCHAR=", ncol(mat), ";"),
    format, labels, "  MATRIX", rows, "  ;", "END;",
    sets, assumptions
  ), file, useBytes = TRUE)
}

.WriteCsv <- function(x, file) {
  write.csv(x, file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

# Read an interleaved or sequential MATRIX into a list of token vectors.
.ReadMatrix <- function(lines) {
  start <- grep("^\\s*MATRIX\\s*$", lines, ignore.case = TRUE)[1] + 1
  end <- start + grep("^\\s*;", lines[start:length(lines)])[1] - 2
  body <- trimws(lines[start:end])
  body <- body[nzchar(body) & !startsWith(body, "[")]
  quoted <- startsWith(body, "'")
  name <- ifelse(
    quoted,
    sub("^'([^']+)'.*$", "\\1", body),
    sub("^(\\S+).*$", "\\1", body)
  )
  data <- ifelse(
    quoted,
    sub("^'[^']+'\\s*", "", body),
    sub("^\\S+\\s*", "", body)
  )
  taxa <- unique(name)
  # Return:
  structure(lapply(taxa, function(tx) {
    .Tokens(paste(data[name == tx], collapse = ""))
  }), names = taxa)
}

.Summary <- function(dataset, tips, nChar, molecular) {
  fossil <- tips[["status"]] == "fossil"
  # Return:
  data.frame(
    dataset = dataset,
    n_taxa = nrow(tips),
    n_fossil = sum(fossil),
    n_extant = sum(!fossil),
    n_morph_chars = nChar,
    molecular_sites = molecular,
    oldest_tip_max_age = max(tips[["max_age"]]),
    youngest_fossil_min_age = min(tips[["min_age"]][fossil])
  )
}

.CheckTaxa <- function(nexus, tips) {
  matTaxa <- rownames(TreeTools::ReadCharacters(nexus))
  if (!setequal(matTaxa, tips[["taxon"]]) ||
      length(matTaxa) != nrow(tips)) {
    stop("Matrix and tip table taxa differ in ", nexus, ": ", paste(
      c(setdiff(matTaxa, tips[["taxon"]]), setdiff(tips[["taxon"]], matTaxa)),
      collapse = ", "
    ))
  }
  message("OK: ", length(matTaxa), " matrix taxa match tip table (", nexus, ")")
}

# Penguins: Gavryushkina et al. 2017 -------------------------------------------

.Penguins <- function() {
  dir <- file.path(empDir, "penguins")
  .CheckRaw("penguins")
  tmp <- tempfile()
  unzip(file.path(dir, "raw", "xmls.zip"), exdir = tmp)
  xmlDir <- file.path(tmp, "xmls")
  main <- read_xml(file.path(xmlDir, "8+dna_R.xml"))

  # Morphology: 245 columns from Ksepka et al. 2012, then seven dummy
  # constant columns (states 0-6) used for the Mkv correction.
  morphSeq <- xml_find_all(main, "//data[@id='penguins_morph']/sequence")
  morph <- do.call(rbind, lapply(xml_attr(morphSeq, "value"), .Tokens))
  rownames(morph) <- xml_attr(morphSeq, "taxon")
  stopifnot(ncol(morph) == 252, all(morph[, 246:252] == rep(0:6, each = 55)))
  morph <- morph[, 1:245]

  stateLab <- xml_find_all(
    main, "//data[@id='penguins_morph']//charstatelabels"
  )
  charName <- trimws(xml_attr(stateLab, "characterName"))
  charNo <- as.integer(sub("^(\\d+).*$", "\\1", charName))
  charName <- gsub("^\\d+\\s+'?|'?$", "", charName)
  stopifnot(identical(charNo, 1:245))
  stateName <- lapply(xml_attr(stateLab, "value"), function(v) {
    labs <- regmatches(v, gregexpr("'[^']*'|[^,' ]+", v))[[1]]
    # Return:
    gsub("^'|'$", "", labs)
  })

  parts <- xml_find_all(main, "//data[@spec='FilteredAlignment']")
  partId <- xml_attr(parts, "id")
  morphPart <- grepl("^penguins[2-7]$", partId)
  nState <- integer(245)
  for (i in which(morphPart)) {
    cols <- .Ranges2Int(xml_attr(parts[i], "filter"))
    nState[cols[cols <= 245]] <- as.integer(sub("penguins", "", partId[i]))
  }
  used <- which(nState > 0)
  observed <- apply(morph, 2, .ObservedStates)
  stopifnot(length(used) == 202, all(observed[-used] <= 1),
            all(observed[used] <= nState[used]))

  analysed <- morph[, used]
  partSets <- split(seq_along(used), paste0("states", nState[used]))
  .WriteNexus(
    file.path(dir, "penguins_morph.nex"), analysed, "STANDARD",
    symbols = "0123456", charLabels = charName[used],
    stateLabels = stateName[used], charSets = partSets,
    comment = paste(
      "Gavryushkina et al. 2017 penguin morphology: the 202 characters",
      "analysed (of 245 in the BEAST2 XML 8+dna_R.xml; the other 43 are",
      "constant once outgroups are removed). CHARSETs give the number of",
      "states assumed for each Mkv partition. All characters unordered."
    )
  )
  .WriteCsv(data.frame(
    char = seq_along(used), xml_column = used,
    n_states_model = nState[used], n_states_observed = observed[used],
    label = charName[used]
  ), file.path(dir, "penguins_characters.csv"))

  dnaSeq <- xml_find_all(main, "//data[@id='penguins']/sequence")
  dna <- do.call(rbind, strsplit(xml_attr(dnaSeq, "value"), ""))
  rownames(dna) <- xml_attr(dnaSeq, "taxon")
  genes <- which(grepl("^penguins_", partId) & partId != "penguins_morph")
  geneSets <- structure(
    lapply(genes, function(i) .Ranges2Int(xml_attr(parts[i], "filter"))),
    names = sub("penguins_", "", partId[genes])
  )
  .WriteNexus(
    file.path(dir, "penguins_dna.nex"), dna[rownames(analysed), ], "DNA",
    charSets = geneSets,
    comment = paste(
      "Gavryushkina et al. 2017 penguin DNA (8+dna_R.xml).",
      "CHARSET names follow the XML partition ids, whose lengths",
      "contradict the gene labels in the SI table; see README."
    )
  )

  dateTrait <- xml_attr(xml_find_first(main, "//trait"), "value")
  dateTrait <- strsplit(trimws(strsplit(dateTrait, ",")[[1]]), "\\s*=\\s*")
  taxa <- vapply(dateTrait, `[`, "", 1)
  .Interval <- function(doc) {
    sd <- xml_find_all(doc, "//*[@spec='beast.evolution.tree.SamplingDate']")
    # Return:
    data.frame(
      taxon = sub("^@", "", xml_attr(sd, "taxon")),
      min_age = as.numeric(xml_attr(sd, "lower")),
      max_age = as.numeric(xml_attr(sd, "upper"))
    )
  }
  ints <- .Interval(main)
  for (f in setdiff(list.files(xmlDir, "\\.xml$"), "8+dna_R.xml")) {
    other <- .Interval(read_xml(file.path(xmlDir, f)))
    key <- match(other[["taxon"]], ints[["taxon"]])
    stopifnot(!anyNA(key), all(other[, 2:3] == ints[key, 2:3]))
  }
  fossil <- match(taxa, ints[["taxon"]])
  tips <- data.frame(
    taxon = taxa,
    status = ifelse(is.na(fossil), "extant", "fossil"),
    min_age = ifelse(is.na(fossil), 0, ints[["min_age"]][fossil]),
    max_age = ifelse(is.na(fossil), 0, ints[["max_age"]][fossil]),
    age_source = ifelse(
      is.na(fossil), "extant (sampled at present)",
      "BEAST2 XML 8+dna_R.xml SamplingDate (Gavryushkina et al. 2017)"
    ),
    note = ""
  )
  si <- read.csv(file.path(dir, "curated", "si-table1.csv"))
  si <- si[match(tips[["taxon"]], si[["taxon"]]), ]
  differ <- !is.na(si[["taxon"]]) & (
    si[["min_age"]] != tips[["min_age"]] | si[["max_age"]] != tips[["max_age"]]
  )
  tips[["note"]][differ] <- paste0(
    "SI Table 1 gives [", si[["min_age"]][differ], ", ",
    si[["max_age"]][differ], "]; XML value used"
  )
  tips <- tips[order(tips[["status"]], tips[["taxon"]]), ]
  .WriteCsv(tips, file.path(dir, "penguins_tips.csv"))
  .WriteCsv(
    data.frame(constraint = character(0), taxa = character(0),
               type = character(0), source = character(0)),
    file.path(dir, "penguins_constraints.csv")
  )
  .CheckTaxa(file.path(dir, "penguins_morph.nex"), tips)
  .CheckTaxa(file.path(dir, "penguins_dna.nex"), tips)
  # Return:
  .Summary("penguins", tips, ncol(analysed), ncol(dna))
}

.Ranges2Int <- function(filter) {
  # Return:
  unlist(lapply(strsplit(gsub("\\s", "", filter), ",")[[1]], function(r) {
    ends <- as.integer(strsplit(r, "-")[[1]])
    ends[1]:ends[length(ends)]
  }))
}

# Hymenoptera: Ronquist et al. 2012; intervals from Zhang et al. 2016 ----------

.MbConstraints <- function(lines) {
  lines <- grep("^\\s*constraint\\s", lines, value = TRUE)
  m <- regmatches(lines, regexec(
    "constraint\\s+(\\w+)\\s*(-1)?\\s*=\\s*(\\d+)-(\\d+|\\.)", lines
  ))
  # Return:
  data.frame(
    constraint = vapply(m, `[`, "", 2),
    from = as.integer(vapply(m, `[`, "", 4)),
    to = vapply(m, `[`, "", 5)
  )
}

.Hymenoptera <- function() {
  dir <- file.path(empDir, "hymenoptera")
  .CheckRaw("hymenoptera")
  nexLines <- readLines(file.path(dir, "raw", "MrBayes_Commands.nex"),
                        encoding = "latin1")
  rows <- .ReadMatrix(nexLines)
  stopifnot(length(rows) == 113, all(lengths(rows) == 5449))
  all <- do.call(rbind, rows)
  morph <- all[, 1:353]
  dna <- all[, 354:5449]

  mb <- nexLines[grep("^begin mrbayes", nexLines, ignore.case = TRUE):
                   length(nexLines)]
  .Charset <- function(name) {
    line <- grep(paste0("charset ", name, "\\s*="), mb, value = TRUE)[1]
    spec <- trimws(sub(";.*$", "", sub("^.*=", "", line)))
    # Return:
    unlist(lapply(strsplit(spec, "\\s+")[[1]], function(r) {
      step <- if (grepl("\\\\", r)) as.integer(sub("^.*\\\\", "", r)) else 1L
      ends <- as.integer(strsplit(sub("\\\\.*$", "", r), "-")[[1]])
      seq(ends[1], ends[length(ends)], by = step)
    }))
  }
  ordered <- .Charset("morph_ordered")
  excluded <- .Charset("morph_excluded")
  .WriteNexus(
    file.path(dir, "hymenoptera_morph.nex"), morph, "STANDARD",
    symbols = "0123456",
    charSets = list(Vilhelmsen = .Charset("MV"), Schulmeister = .Charset("MS")),
    ordered = ordered, excluded = excluded,
    comment = paste(
      "Ronquist et al. 2012 Hymenoptera morphology, all 353 characters of",
      "MrBayes_Commands.nex. The analyses exclude the 10 characters in EXSET",
      "(343 analysed) and treat the TYPESET ORD characters as ordered."
    )
  )
  .WriteCsv(data.frame(
    char = 1:353,
    block = ifelse(1:353 <= 236, "Vilhelmsen", "Schulmeister"),
    ordered = 1:353 %in% ordered,
    excluded = 1:353 %in% excluded,
    n_states_observed = apply(morph, 2, .ObservedStates)
  ), file.path(dir, "hymenoptera_characters.csv"))

  geneNames <- c("12S", "16S", "18S", "28S", "CO1", "Ef1aF2", "Ef1aF1",
                 "CO1_12", "CO1_3", "Ef1aF2_12", "Ef1aF2_3",
                 "Ef1aF1_12", "Ef1aF1_3")
  geneSets <- structure(lapply(geneNames, function(g) .Charset(g) - 353L),
                        names = geneNames)
  .WriteNexus(
    file.path(dir, "hymenoptera_dna.nex"), dna, "DNA", charSets = geneSets,
    comment = paste(
      "Ronquist et al. 2012 Hymenoptera DNA: columns 354-5449 of",
      "MrBayes_Commands.nex renumbered from 1. The analyses exclude CO1_3."
    )
  )

  fixed <- regmatches(mb, regexec("(\\w+)=Fixed\\((\\d+)\\)", mb))
  fixed <- do.call(rbind, fixed[lengths(fixed) == 3])
  fixedAge <- structure(as.numeric(fixed[, 3]), names = fixed[, 2])
  fossilLine <- grep("taxset fossils", mb, value = TRUE)
  fossils <- strsplit(trimws(sub(";.*$", "", sub("^.*=", "", fossilLine))),
                      "\\s+")[[1]]
  stopifnot(length(fossils) == 45, setequal(fossils, names(fixedAge)))

  curated <- read.csv(file.path(dir, "curated", "fossil-ages.csv"))
  stopifnot(setequal(curated[["taxon"]], fossils))
  commandAge <- fixedAge[curated[["taxon"]]]
  tableDiffers <- curated[["age_ronquist2012"]] != commandAge
  curated[["note"]][tableDiffers] <- paste0(
    "Ronquist et al. 2012 Table 2 gives ",
    curated[["age_ronquist2012"]][tableDiffers],
    " Ma; MrBayes_Commands.nex fixes ", commandAge[tableDiffers],
    " Ma. ", curated[["note"]][tableDiffers]
  )

  # The numeric constraint ranges index a 114-taxon order that includes
  # Lepid_Micropterigidae, as in MrBayes's examples/hymfossil.nex; in this
  # 113-taxon file they are shifted by one.
  mfLines <- readLines(file.path(dir, "raw", "hymfossil.nex"))
  order114 <- names(.ReadMatrix(mfLines))
  stopifnot(length(order114) == 114,
            identical(setdiff(order114, "Lepid_Micropterigidae"), names(rows)))
  mfFixed <- unlist(regmatches(
    mfLines, gregexpr("\\w+=fixed\\(\\d+\\)", mfLines)
  ))
  mfAge <- structure(
    as.numeric(sub("^.*\\((\\d+)\\)$", "\\1", mfFixed)),
    names = sub("=.*$", "", mfFixed)
  )
  stopifnot(setequal(names(mfAge), fossils))

  cons <- .MbConstraints(mb)
  .Members <- function(from, to, taxa) {
    to <- if (to == ".") length(taxa) else min(as.integer(to), length(taxa))
    # Return:
    taxa[from:to]
  }
  cons[["taxa"]] <- mapply(function(f, t) {
    paste(setdiff(.Members(f, t, order114), "Lepid_Micropterigidae"),
          collapse = " ")
  }, cons[["from"]], cons[["to"]])
  cons[["n_taxa"]] <- lengths(strsplit(cons[["taxa"]], " "))
  cons[["literal_113_matches"]] <- mapply(function(f, t, tx) {
    identical(paste(.Members(f, t, names(rows)), collapse = " "), tx)
  }, cons[["from"]], cons[["to"]], cons[["taxa"]])
  cons[["used_in_ted"]] <- cons[["constraint"]] %in%
    c("root", "holometabola_withFossils")
  cons[["range_in_source"]] <- paste0(cons[["from"]], "-", cons[["to"]])
  cons <- cons[, c("constraint", "used_in_ted", "n_taxa", "range_in_source",
                   "literal_113_matches", "taxa")]
  .WriteCsv(cons, file.path(dir, "hymenoptera_constraints.csv"))
  hymenoptera <- strsplit(cons[["taxa"]][cons[["constraint"]] ==
                                           "hymenoptera"], " ")[[1]]
  stopifnot(length(hymenoptera) == 60)

  taxa <- names(rows)
  fossil <- match(taxa, curated[["taxon"]])
  tips <- data.frame(
    taxon = taxa,
    status = ifelse(is.na(fossil), "extant", "fossil"),
    group = ifelse(taxa %in% hymenoptera | !is.na(fossil),
                   "Hymenoptera", "outgroup"),
    min_age = ifelse(is.na(fossil), 0, curated[["min_age"]][fossil]),
    max_age = ifelse(is.na(fossil), 0, curated[["max_age"]][fossil]),
    point_age_ronquist2012 = ifelse(is.na(fossil), 0, fixedAge[taxa]),
    point_age_hymfossil = ifelse(is.na(fossil), 0, mfAge[taxa]),
    age_source = ifelse(
      is.na(fossil), "extant (sampled at present)",
      paste("Zhang et al. 2016 Table 3 stratum interval, matched via",
            "Ronquist et al. 2012 Table 2 stratum")
    ),
    assignment = ifelse(is.na(fossil), "", curated[["assignment"]][fossil]),
    stratum = ifelse(is.na(fossil), "",
                     curated[["stratum_provenance"]][fossil]),
    note = ifelse(is.na(fossil), "", trimws(curated[["note"]][fossil]))
  )
  midOff <- which(tips[["point_age_hymfossil"]] !=
    floor(curated[["age_class_zhang2016"]][fossil]))
  message("hymfossil.nex ages differing from floor(Zhang class): ",
          paste(taxa[midOff], collapse = ", "))
  .WriteCsv(tips, file.path(dir, "hymenoptera_tips.csv"))
  .CheckTaxa(file.path(dir, "hymenoptera_morph.nex"), tips)
  .CheckTaxa(file.path(dir, "hymenoptera_dna.nex"), tips)
  # Return:
  .Summary("hymenoptera", tips, ncol(morph) - length(excluded), ncol(dna))
}

# Trilobites: Paterson, Edgecombe & Lee 2019 -----------------------------------

.Trilobites <- function() {
  dir <- file.path(empDir, "trilobites")
  .CheckRaw("trilobites")
  nexLines <- readLines(
    file.path(dir, "raw", "Dataset_2_Discrete_characters.nex"),
    encoding = "UTF-8"
  )
  rows <- .ReadMatrix(nexLines)
  stopifnot(length(rows) == 107, all(lengths(rows) == 107))
  morph <- do.call(rbind, rows)
  rownames(morph) <- gsub(" ", "_", rownames(morph))

  labStart <- grep("^\\s*CHARLABELS", nexLines) + 1
  labLines <- grep("^\\s*\\[\\d+\\]", nexLines[labStart:length(nexLines)],
                   value = TRUE)
  charLabel <- gsub("^\\s*\\[\\d+\\]\\s*'?|'?\\s*$", "", labLines)
  stopifnot(length(charLabel) == 107)
  .WriteNexus(
    file.path(dir, "trilobites_morph.nex"), morph, "STANDARD",
    symbols = "0123456789ABCDEF", charLabels = charLabel,
    comment = paste(
      "Paterson, Edgecombe & Lee 2019 Cambrian trilobites: 107 discrete",
      "characters of Dataset_2_Discrete_characters.nex, unordered as",
      "deposited. Taxon names use underscores."
    )
  )
  # Ordered characters and state-count partitions in the BEAST 1.8 analyses
  # were identified by matching alignment columns to the NEXUS; see README.
  .WriteCsv(data.frame(
    char = 1:107, label = charLabel,
    n_states_observed = apply(morph, 2, .ObservedStates),
    ordered_in_beast = 1:107 %in% c(45, 47)
  ), file.path(dir, "trilobites_characters.csv"))

  xl <- as.data.frame(read_excel(
    file.path(dir, "raw", "Dataset_1_Taxa,Dates,CodingRefs.xlsx"),
    col_types = "text"
  ))
  names(xl) <- c("taxon", "classification", "formation", "biozone",
                 "gts2012", "age", "coding_refs")
  xl <- xl[!is.na(xl[["age"]]) & !is.na(suppressWarnings(
    as.numeric(xl[["age"]])
  )), ]
  xlName <- gsub(" ", "_", trimws(xl[["taxon"]]))
  # The spreadsheet spells one species differently from the matrix and the
  # BEAST XMLs; the matrix spelling is kept as the key.
  xlName[xlName == "Neocobboldia_chinlinica"] <- "Neocobboldia_chilinica"
  age <- as.numeric(xl[["age"]])
  tips <- data.frame(
    taxon = xlName, status = "fossil", min_age = age, max_age = age,
    age_source = "Dataset_1 Age (Ma), point age (Paterson et al. 2019)",
    name_in_spreadsheet = trimws(xl[["taxon"]]),
    order = NA_character_,
    classification = xl[["classification"]],
    gts2012_stage = xl[["gts2012"]], biozone = xl[["biozone"]],
    formation = xl[["formation"]]
  )
  orders <- c("EODISCIDA", "REDLICHIIDA", "CORYNEXOCHIDA", "AULACOPLEURIDA",
              "ASAPHIDA", "OLENIDA", "HARPIDA", "UNCERTAIN")
  allNames <- trimws(as.data.frame(read_excel(
    file.path(dir, "raw", "Dataset_1_Taxa,Dates,CodingRefs.xlsx"),
    col_types = "text"
  ))[[1]])
  header <- cumsum(allNames %in% orders)
  tips[["order"]] <- orders[header[match(tips[["name_in_spreadsheet"]],
                                         allNames)]]

  tmp <- tempfile()
  unzip(file.path(dir, "raw", "Dataset_4_BEAST_and_BEAST2.zip"), exdir = tmp)
  beastDir <- file.path(tmp, "SI_file4_BEAST_and_BEAST2",
                        "BEAST_1.8_executables_UPDATED")
  ucln <- read_xml(file.path(beastDir, paste0(
    "UCLN_CambrianTrilos_107_2_6_DisPartScaledGamUCLN_",
    "ContMuDeltaUCLN.xml"
  )))
  beastTaxa <- xml_find_all(ucln, "/beast/taxa[1]/taxon")
  beastAge <- as.numeric(xml_attr(xml_find_first(beastTaxa, "date"), "value"))
  names(beastAge) <- xml_attr(beastTaxa, "id")
  stopifnot(setequal(names(beastAge), tips[["taxon"]]),
            all(beastAge[tips[["taxon"]]] == tips[["min_age"]]))
  ingroup <- xml_attr(xml_find_all(
    ucln, "//taxa[@id='Ingroup_node']/taxon"
  ), "idref")
  outgroup <- xml_attr(xml_find_all(
    ucln, "//taxa[@id='Fallotaspidoids']/taxon"
  ), "idref")
  stopifnot(setequal(c(ingroup, outgroup), tips[["taxon"]]))
  .WriteCsv(data.frame(
    constraint = c("Ingroup_node", "Fallotaspidoids"),
    type = c("monophyletic (enforced)", "outgroup (monitored only)"),
    n_taxa = c(length(ingroup), length(outgroup)),
    taxa = c(paste(ingroup, collapse = " "), paste(outgroup, collapse = " ")),
    source = "BEAST 1.8 XMLs (Dataset_4)"
  ), file.path(dir, "trilobites_constraints.csv"))

  quant <- as.data.frame(read_excel(
    file.path(dir, "raw", "Dataset_3_Meristic_Continuous_characters.xlsx"),
    col_names = FALSE, col_types = "text", .name_repair = "minimal"
  ))
  head1 <- as.character(quant[1, ])
  head2 <- as.character(quant[2, ])
  quantNames <- trimws(ifelse(is.na(head2), head1,
                              ifelse(is.na(head1), head2,
                                     paste(head1, head2, sep = ": "))))
  quant <- quant[-(1:2), ]
  quant <- quant[!is.na(quant[[2]]), ]
  names(quant) <- c("taxon", quantNames[-1])
  quant[["taxon"]] <- gsub(" ", "_", trimws(quant[["taxon"]]))
  quant[["taxon"]][quant[["taxon"]] == "Neocobboldia_chinlinica"] <-
    "Neocobboldia_chilinica"
  stopifnot(setequal(quant[["taxon"]], tips[["taxon"]]))
  .WriteCsv(quant, file.path(dir, "trilobites_quantitative.csv"))

  tips <- tips[order(tips[["taxon"]]), ]
  .WriteCsv(tips, file.path(dir, "trilobites_tips.csv"))
  .CheckTaxa(file.path(dir, "trilobites_morph.nex"), tips)
  # Return:
  .Summary("trilobites", tips, ncol(morph), 0)
}

summary <- rbind(.Penguins(), .Hymenoptera(), .Trilobites())
.WriteCsv(summary, file.path(empDir, "summary.csv"))
print(summary)
