# Validation 1 of #97: fixed-tree character log-likelihood of MkCharacters()
# with discrete Gamma rates against RevBayes
#   dnPhyloCTMC(Q = fnJC(K) | fnOrderedRateMatrix(maxState = K),
#               branchRates = rate, siteRates = fnDiscretizeGamma(a, a, 4),
#               type = "Standard", coding = "variable")
# at 3 clock rates x 3 shapes (and without Gamma), in both engines.
#
#   Rscript data-raw/gamma-acrv-revbayes.R prepare <dir>
#     simulate the fixed tree (fossil tips, a sampled ancestor) and a matrix
#     of unordered (K = 2, 3, 4) and ordered (K = 3, 4, 5) characters with
#     missing cells, gaps and ambiguity sets, written as <dir>/tree.nwk,
#     <dir>/cells.csv, <dir>/partitions.csv and one NEXUS file per partition
#     for gamma-acrv-revbayes.Rev
#   Rscript data-raw/gamma-acrv-revbayes.R compare <dir> <rb output> <out.csv>
#     CET (both engines) against the RevBayes values printed by
#     gamma-acrv-revbayes.Rev; stops if any differs by more than 1e-6
#   Rscript data-raw/gamma-acrv-revbayes.R hymenoptera <cet dir> <rb output> <out.csv>
#     the same on the Hymenoptera matrix and start tree of the CET pipeline
#     (branch hymenoptera-cet-pipeline, data-raw/empirical/hymenoptera/cet),
#     against its RevBayes oracle (revbayes/check-mkv.Rev: coding "variable"
#     without and with Gamma(4) at alpha 1.10, and coding "all")
#
# With EvoTempo installed.  RevBayes prints 7 significant digits, so the Rev
# script prints each value as floor(lnL) and lnL - floor(lnL).
suppressMessages({
  library(EvoTempo)
  library(ape)
})
args <- commandArgs(TRUE)
mode <- args[[1]]
rates <- c(0.1, 0.5, 2)
alphas <- c(0.3, 1.1, 5)
tol <- 1e-6

# Read "RESULT i j k" blocks: floor(lnL) and the remainder on the next lines
ReadRevBayes <- function(file) {
  rb <- gsub("\"", "", trimws(readLines(file)))
  at <- grep("^RESULT ", rb)
  num <- function(s) as.numeric(trimws(s))
  do.call(rbind, lapply(at, function(i) {
    f <- as.integer(strsplit(sub("^RESULT ", "", rb[i]), " ")[[1]])
    data.frame(i1 = f[1], i2 = f[2], i3 = f[3], lnL = num(rb[i + 1]) + num(rb[i + 2]))
  }))
}

# Branch durations of a TempoTree (0 at the root)
Durations <- function(tt) {
  dur <- numeric(tt$nTip + tt$nNode)
  child <- which(!is.na(tt$parent))
  dur[child] <- tt$age[tt$parent[child]] - tt$age[child]
  dur
}

# TempoTree of a Newick time tree: tips at the maximum depth are extant,
# zero-length tips are sampled ancestors, other tips fossils
ReadTimeTree <- function(file) {
  tr <- read.tree(file)
  nTip <- Ntip(tr)
  depth <- node.depth.edgelength(tr)
  age <- max(depth) - depth
  age[abs(age) < 1e-9] <- 0
  len <- numeric(nTip + tr$Nnode)
  len[tr$edge[, 2]] <- tr$edge.length
  tip <- seq_len(nTip)
  tipType <- ifelse(age[tip] == 0, "extant", ifelse(len[tip] < 1e-12, "ancestor", "fossil"))
  list(phylo = tr, tree = TempoTree(tr, ages = age, tipType = tipType))
}

# Insert a sampled ancestor (a tip on a zero-length branch) on edge `edge`,
# at distance `pos` above its child node (ape only)
AddSampledAncestor <- function(tr, edge, pos, label) {
  n <- Ntip(tr)
  m <- tr$Nnode
  e <- tr$edge
  e[e > n] <- e[e > n] + 1L  # internal nodes move up by one
  newTip <- n + 1L
  newNode <- n + 1L + m + 1L
  child <- e[edge, 2]
  len <- tr$edge.length
  e[edge, 2] <- newNode
  len[edge] <- len[edge] - pos
  tr$edge <- rbind(e, c(newNode, child), c(newNode, newTip))
  tr$edge.length <- c(len, pos, 0)
  tr$tip.label <- c(tr$tip.label, label)
  tr$Nnode <- m + 1L
  reorder(tr, "cladewise")
}

CetLogLik <- function(tt, cells, cols, K, ordered, rate, alpha, coding = "variable",
                      engine = "cpp") {
  nCat <- if (is.na(alpha)) 1L else 4L
  mk <- MkCharacters(cells[, cols, drop = FALSE], nState = K, ordered = ordered,
                     coding = coding, nRateCat = nCat)
  MkLikelihood(tt, mk, rate * Durations(tt), alpha = if (nCat > 1) alpha,
               engine = engine)
}

if (mode == "prepare") {
  dir <- args[[2]]
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(9701)
  # Time tree: 24 extant tips (coalescent, root age 1), 4 tips cut back to
  # fossils, one sampled ancestor on an internal branch
  tr <- rcoal(24, tip.label = sprintf("t%02d", 1:24))
  tr$edge.length <- tr$edge.length / max(node.depth.edgelength(tr))
  tips <- which(tr$edge[, 2] <= 24)
  fossils <- sample(tips, 4)
  tr$edge.length[fossils] <- tr$edge.length[fossils] * runif(4, 0.2, 0.8)
  tr$tip.label[tr$edge[fossils, 2]] <- sprintf("f%d", 1:4)
  internal <- which(tr$edge[, 2] > 24)
  host <- internal[which.max(tr$edge.length[internal])]
  tr <- AddSampledAncestor(tr, host, 0.4 * tr$edge.length[host], "sa1")
  write.tree(tr, file.path(dir, "tree.nwk"), digits = 12)
  tt <- ReadTimeTree(file.path(dir, "tree.nwk"))
  tr <- tt$phylo
  stopifnot(sum(tt$tree$tipType == "fossil") == 4, sum(tt$tree$tipType == "ancestor") == 1)
  nTip <- Ntip(tr)
  parts <- data.frame(partition = c("u2", "u3", "u4", "o3", "o4", "o5"),
                      nState = c(2, 3, 4, 3, 4, 5),
                      ordered = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
                      nChar = c(30, 12, 8, 10, 8, 6))
  # Simulate on the tree in time units, rate 0.6, per-character Gamma(0.8)
  # categories
  Simulate <- function(K, ordered) {
    Q <- EvoTempo:::.MkRateMatrix(K, ordered)
    e <- eigen(Q, symmetric = TRUE)
    r <- sample(GammaRates(0.8, 4), 1) * 0.6
    state <- integer(nTip + tr$Nnode)
    state[nTip + 1] <- sample.int(K, 1) - 1L
    for (i in seq_len(nrow(tr$edge))) {  # cladewise: parents first
      P <- e$vectors %*% (exp(e$values * r * tr$edge.length[i]) * t(e$vectors))
      state[tr$edge[i, 2]] <- sample.int(K, 1, prob = pmax(P[state[tr$edge[i, 1]] + 1, ], 0)) - 1L
    }
    state[seq_len(nTip)]
  }
  cells <- NULL
  partOf <- character(0)
  for (p in seq_len(nrow(parts))) {
    K <- parts$nState[p]
    for (j in seq_len(parts$nChar[p])) {
      repeat {
        x <- as.character(Simulate(K, parts$ordered[p]))
        # Missing cells ('?' or '-'), and ambiguity sets holding the true state
        miss <- runif(nTip) < 0.15
        x[miss] <- ifelse(runif(sum(miss)) < 0.7, "?", "-")
        amb <- which(!miss & runif(nTip) < 0.08)
        for (i in amb) {
          s <- as.integer(x[i])
          extra <- sample(setdiff(0:(K - 1), s), sample.int(min(2, K - 1), 1))
          x[i] <- paste0("{", paste(sort(c(s, extra)), collapse = ""), "}")
        }
        # Variable under RevBayes's rule: two distinct non-gap cells
        full <- nchar(gsub("[^0-9]", "", x)) == K & grepl("^\\{", x)
        scoredCells <- x[!x %in% c("?", "-") & !full]
        if (length(unique(scoredCells)) >= 2) break
      }
      cells <- cbind(cells, x)
      partOf <- c(partOf, parts$partition[p])
    }
  }
  rownames(cells) <- tr$tip.label
  colnames(cells) <- sprintf("c%d", seq_len(ncol(cells)))
  write.csv(data.frame(taxon = tr$tip.label, cells, check.names = FALSE),
            file.path(dir, "cells.csv"), row.names = FALSE)
  write.csv(parts, file.path(dir, "partitions.csv"), row.names = FALSE)
  width <- max(nchar(tr$tip.label)) + 2
  for (p in seq_len(nrow(parts))) {
    cols <- which(partOf == parts$partition[p])
    rows <- paste0("  ", formatC(tr$tip.label, width = -width),
                   apply(cells[, cols, drop = FALSE], 1, paste, collapse = ""))
    writeLines(c("#NEXUS", "BEGIN DATA;",
                 sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", nTip, length(cols)),
                 sprintf("  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"%s\";",
                         paste(0:(parts$nState[p] - 1), collapse = "")),
                 "  MATRIX", rows, "  ;", "END;"),
               file.path(dir, sprintf("part_%d.nex", p)))
  }
  writeLines(c(sprintf("n_part = %d", nrow(parts)),
               sprintf("part_k = v(%s)", paste(parts$nState, collapse = ", ")),
               sprintf("part_ordered = v(%s)",
                       paste(ifelse(parts$ordered, "TRUE", "FALSE"), collapse = ", "))),
             file.path(dir, "partitions.Rev"))
  cat("prepared", ncol(cells), "characters on", nTip, "tips in", dir, "\n")
  cat("cells: ", sum(cells %in% c("?", "-")), "missing,", sum(grepl("^\\{", cells)),
      "sets\n")
}

if (mode == "compare") {
  dir <- args[[2]]
  rb <- ReadRevBayes(args[[3]])
  names(rb)[1:3] <- c("rateIndex", "partition", "alphaIndex")
  tt <- ReadTimeTree(file.path(dir, "tree.nwk"))$tree
  raw <- read.csv(file.path(dir, "cells.csv"), colClasses = "character",
                  check.names = FALSE, na.strings = character(0))
  cells <- as.matrix(raw[match(tt$tipLabel, raw$taxon), -1])
  parts <- read.csv(file.path(dir, "partitions.csv"))
  partOf <- rep(seq_len(nrow(parts)), parts$nChar)
  rows <- list()
  for (ri in seq_along(rates)) for (p in seq_len(nrow(parts))) for (ai in 1:4) {
    alpha <- if (ai <= 3) alphas[ai] else NA
    cols <- which(partOf == p)
    cpp <- CetLogLik(tt, cells, cols, parts$nState[p], parts$ordered[p], rates[ri], alpha)
    ref <- CetLogLik(tt, cells, cols, parts$nState[p], parts$ordered[p], rates[ri], alpha,
                     engine = "R")
    v <- rb$lnL[rb$rateIndex == ri & rb$partition == p & rb$alphaIndex == ai]
    stopifnot(length(v) == 1)
    rows[[length(rows) + 1]] <- data.frame(
      rate = rates[ri], partition = parts$partition[p], nChar = length(cols),
      alpha = alpha, revbayes = v, cpp = cpp, R = ref, diffCpp = cpp - v,
      diffR = ref - v, diffEngines = cpp - ref)
  }
  res <- do.call(rbind, rows)
  print(res, digits = 10, row.names = FALSE)
  byShape <- res
  byShape$alpha[is.na(byShape$alpha)] <- Inf  # Inf: no Gamma
  tot <- aggregate(cbind(revbayes, cpp, R) ~ rate + alpha, byShape, sum)
  tot$diff <- tot$cpp - tot$revbayes
  print(tot, digits = 10, row.names = FALSE)
  write.csv(res, args[[4]], row.names = FALSE)
  worst <- max(abs(c(res$diffCpp, res$diffR)) / pmax(1, abs(res$revbayes)))
  worstEng <- max(abs(res$diffEngines) / pmax(1, abs(res$R)))
  cat(sprintf("max |CET - RevBayes| / max(1, |lnL|) = %.3g; max engine difference %.3g\n",
              worst, worstEng))
  if (worst > tol || worstEng > 1e-8) stop("gamma-acrv-revbayes.R compare FAILED")
  cat("compare: all checks passed\n")
}

if (mode == "hymenoptera") {
  cetDir <- args[[2]]
  rb <- ReadRevBayes(args[[3]])
  names(rb)[1:3] <- c("rateIndex", "part", "coding")
  rbText <- gsub("\"", "", trimws(readLines(args[[3]])))
  pl <- grep("^PART ", rbText, value = TRUE)
  rbNames <- setNames(sub("^PART [0-9]+ ", "", pl), sub("^PART ([0-9]+) .*$", "\\1", pl))
  hymRates <- c(0.0007, 0.002, 0.006)
  # The pipeline's LoadHymenoptera(), without its namespace rebinding
  tips <- read.csv(file.path(cetDir, "..", "hymenoptera_tips.csv"))
  tr <- read.tree(file.path(cetDir, "data", "start-tree.nwk"))
  nTip <- Ntip(tr)
  info <- tips[match(tr$tip.label, tips$taxon), ]
  depth <- node.depth.edgelength(tr)
  extant <- info$status == "extant"
  age <- max(depth[which(extant)]) - depth
  age[which(extant)] <- 0
  len <- numeric(nTip + tr$Nnode)
  len[tr$edge[, 2]] <- tr$edge.length
  tipType <- ifelse(extant, "extant", ifelse(len[seq_len(nTip)] < 1e-9, "ancestor", "fossil"))
  ranges <- cbind(info$min_age, info$max_age)
  ranges[extant, ] <- 0
  fos <- which(tipType == "fossil")
  age[fos] <- pmin(pmax(age[fos], ranges[fos, 1]), ranges[fos, 2])
  tt <- TempoTree(tr, ages = age, tipType = tipType, ranges = ranges)
  raw <- read.csv(file.path(cetDir, "data", "characters.csv"), check.names = FALSE,
                  colClasses = "character", na.strings = character(0))
  parts <- read.csv(file.path(cetDir, "data", "partitions.csv"))
  cells <- as.matrix(raw[match(tr$tip.label, raw$taxon), -1])
  rows <- list()
  for (ri in seq_along(hymRates)) for (pi in names(rbNames)) {
    lab <- rbNames[[pi]]
    # Rows of partitions.csv index the analysed characters (columns of
    # characters.csv); `column` is the position in the source matrix
    cols <- which(parts$partition == lab)
    stopifnot(ncol(cells) == nrow(parts))
    K <- parts$nState[cols[1]]
    ord <- parts$ordered[cols[1]]
    for (cd in 1:3) {
      alpha <- if (cd == 3) 1.10 else NA
      coding <- if (cd == 2) "all" else "variable"
      cpp <- CetLogLik(tt, cells, cols, K, ord, hymRates[ri], alpha, coding)
      ref <- CetLogLik(tt, cells, cols, K, ord, hymRates[ri], alpha, coding, engine = "R")
      v <- rb$lnL[rb$rateIndex == ri & rb$part == as.integer(pi) & rb$coding == cd]
      rows[[length(rows) + 1]] <- data.frame(
        rate = hymRates[ri], partition = lab, nChar = length(cols), coding = coding,
        alpha = alpha, revbayes = v, cpp = cpp, R = ref, diffCpp = cpp - v,
        diffR = ref - v, diffEngines = cpp - ref)
    }
  }
  res <- do.call(rbind, rows)
  print(res, digits = 10, row.names = FALSE)
  write.csv(res, args[[4]], row.names = FALSE)
  worst <- max(abs(c(res$diffCpp, res$diffR)) / pmax(1, abs(res$revbayes)))
  worstEng <- max(abs(res$diffEngines) / pmax(1, abs(res$R)))
  cat(sprintf("Hymenoptera: max |CET - RevBayes| / max(1, |lnL|) = %.3g; max engine difference %.3g\n",
              worst, worstEng))
  if (worst > tol || worstEng > 1e-8) stop("gamma-acrv-revbayes.R hymenoptera FAILED")
  cat("hymenoptera: all checks passed\n")
}
