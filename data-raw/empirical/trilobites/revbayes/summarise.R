# Convergence diagnostics and posterior summaries for the RevBayes baselines
# (adapted from data-raw/empirical/penguins/revbayes/summarise.R).
#
#   Rscript summarise.R <output dir> <prefix> <config> [burn-in fraction]
#
# Reads <dir>/<prefix>_<config>_run*.log and .trees (RevBayes), discards the
# first 25 % of each chain, and writes to ./summary/:
#   convergence_<config>.csv  per parameter: ESS per chain (min), pooled ESS,
#                             PSRF, plus the topology ASDSF row
#   parameters_<config>.csv   posterior mean, median and 95 % HPD per parameter
#   published_<config>.csv    estimates vs. the source studies (published.csv)
#   fossils_<config>.csv      per fossil: P(sampled ancestor), P(in crown group)
#   clades_<config>.csv       clade frequencies per chain (for the ASDSF)
#
# Needs only ape (no coda on the cluster).  ESS is coda::effectiveSize()'s
# estimator (spectral density at zero from an AR fit chosen by AIC); PSRF is
# coda::gelman.diag()'s point estimate; HPD is the shortest interval holding
# 95 % of the samples, as coda::HPDinterval().  check-diagnostics.R checks all
# three against coda.

args <- commandArgs(trailingOnly = TRUE)
outDir <- args[1]
prefix <- args[2]
clock <- args[3]
burnin <- if (length(args) > 3) as.numeric(args[4]) else 0.25
dir.create("summary", showWarnings = FALSE)

library(ape)

ESS <- function(x) {
  n <- length(x)
  if (var(x) == 0) return(0)
  fit <- ar(x, aic = TRUE)
  spec <- fit$var.pred / (1 - sum(fit$ar))^2
  n * var(x) / spec
}

PSRF <- function(chains) {
  # coda::gelman.diag() point estimate (Brooks & Gelman 1998 df correction)
  n <- min(lengths(chains))
  chains <- lapply(chains, function(x) x[seq_len(n)])
  m <- length(chains)
  xbar <- vapply(chains, mean, 0)
  s2 <- vapply(chains, var, 0)
  w <- mean(s2)
  if (w == 0) return(NA_real_)
  b <- n * var(xbar)
  muhat <- mean(xbar)
  varW <- var(s2) / m
  varB <- 2 * b^2 / (m - 1)
  covWB <- (n / m) * (cov(s2, xbar^2) - 2 * muhat * cov(s2, xbar))
  v <- (n - 1) * w / n + (1 + 1 / m) * b / n
  varV <- ((n - 1)^2 * varW + (1 + 1 / m)^2 * varB + 2 * (n - 1) * (1 + 1 / m) * covWB) / n^2
  dfV <- 2 * v^2 / varV
  sqrt((dfV + 3) / (dfV + 1) * ((n - 1) / n + (1 + 1 / m) * b / (n * w)))
}

HPD <- function(x, prob = 0.95) {
  x <- sort(x)
  n <- length(x)
  gap <- max(1, min(n - 1, round(n * prob)))
  lower <- seq_len(n - gap)
  i <- which.min(x[lower + gap] - x[lower])
  c(x[i], x[i + gap])
}

logFiles <- sort(Sys.glob(file.path(outDir, sprintf("%s_%s_run*.log", prefix, clock))))
stopifnot(length(logFiles) >= 2)
runNames <- sub("^.*_(run[0-9]+)\\.log$", "\\1", logFiles)
Discard <- function(d) d[-seq_len(floor(nrow(d) * burnin)), , drop = FALSE]
logs <- lapply(logFiles, function(f) Discard(read.delim(f, check.names = FALSE)))
names(logs) <- runNames

params <- setdiff(Reduce(intersect, lapply(logs, names)), "Iteration")
params <- params[vapply(params, function(p) var(unlist(lapply(logs, `[[`, p))) > 0, TRUE)]

conv <- do.call(rbind, lapply(params, function(p) {
  chains <- lapply(logs, `[[`, p)
  ess <- vapply(chains, ESS, 0)
  data.frame(parameter = p, ess_min_chain = min(ess), ess_total = sum(ess),
             psrf = PSRF(chains), n_chains = length(chains),
             samples_per_chain = min(lengths(chains)))
}))

pooled <- do.call(rbind, logs)
summ <- do.call(rbind, lapply(params, function(p) {
  x <- pooled[[p]]
  h <- HPD(x)
  data.frame(parameter = p, mean = mean(x), median = median(x),
             hpd95_lower = h[1], hpd95_upper = h[2])
}))

# Trees ------------------------------------------------------------------------
taxaTable <- read.delim("data/taxa.tsv")
extant <- taxaTable$taxon[taxaTable$max_age == 0]
fossils <- if (length(extant)) taxaTable$taxon[taxaTable$max_age > 0] else taxaTable$taxon
# crown group for P(fossil in crown): data/crown_taxa.txt if present
crownTaxa <- if (file.exists("data/crown_taxa.txt")) readLines("data/crown_taxa.txt") else character(0)

ReadTrees <- function(f) {
  d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  d <- Discard(d)
  lapply(d[[ncol(d)]], function(s) read.tree(text = s))
}
treeFiles <- sub("\\.log$", ".trees", logFiles)
trees <- lapply(treeFiles, ReadTrees)
names(trees) <- runNames
allTips <- sort(trees[[1]][[1]]$tip.label)

TreeFacts <- function(tr) {
  nTip <- length(tr$tip.label)
  tipEdge <- match(seq_len(nTip), tr$edge[, 2])
  sa <- tr$edge.length[tipEdge] < 1e-8
  names(sa) <- tr$tip.label
  inCrown <- rep(NA, nTip)
  names(inCrown) <- tr$tip.label
  if (length(crownTaxa)) {
    crown <- getMRCA(tr, crownTaxa)
    inCrown[] <- tr$tip.label %in% tr$tip.label[unlist(Descendants1(tr, crown))]
  }
  # clades as sorted tip-label strings (rooted; sampled ancestors are tips)
  pp <- prop.part(tr)
  clades <- vapply(pp, function(ix) paste(sort(attr(pp, "labels")[ix]), collapse = ","), "")
  clades <- clades[lengths(pp) > 1 & lengths(pp) < nTip]
  list(sa = sa[allTips], inCrown = inCrown[allTips], clades = clades)
}
Descendants1 <- function(tr, node) {
  nTip <- length(tr$tip.label)
  out <- integer(0)
  stack <- node
  while (length(stack)) {
    v <- stack[1]
    stack <- stack[-1]
    kids <- tr$edge[tr$edge[, 1] == v, 2]
    out <- c(out, kids[kids <= nTip])
    stack <- c(stack, kids[kids > nTip])
  }
  list(out)
}

facts <- lapply(trees, function(ch) lapply(ch, TreeFacts))

cladeFreq <- lapply(facts, function(ch) table(unlist(lapply(ch, `[[`, "clades"))) / length(ch))
allClades <- unique(unlist(lapply(cladeFreq, names)))
freqMat <- vapply(cladeFreq, function(f) {
  v <- as.numeric(f[allClades])
  v[is.na(v)] <- 0
  v
}, numeric(length(allClades)))
freqMat <- matrix(freqMat, ncol = length(cladeFreq), dimnames = list(allClades, runNames))
keep <- apply(freqMat, 1, max) >= 0.1
sdFreq <- apply(freqMat[keep, , drop = FALSE], 1, sd)
asdsf <- mean(sdFreq)
conv <- rbind(conv, data.frame(parameter = "topology (ASDSF, clades >= 0.1)", ess_min_chain = NA,
                               ess_total = NA, psrf = NA, n_chains = length(trees),
                               samples_per_chain = min(lengths(trees))))
conv$asdsf <- NA
conv$asdsf[nrow(conv)] <- asdsf
conv$max_sdsf <- NA
conv$max_sdsf[nrow(conv)] <- max(sdFreq)

cladeOut <- data.frame(clade = rownames(freqMat), freqMat, mean = rowMeans(freqMat),
                       sd = apply(freqMat, 1, sd), row.names = NULL, check.names = FALSE)
cladeOut <- cladeOut[order(-cladeOut$mean), ]
cladeOut <- cladeOut[cladeOut$mean >= 0.01, ]

flat <- unlist(facts, recursive = FALSE)
pSA <- rowMeans(vapply(flat, `[[`, logical(length(allTips)), "sa"))
pCrown <- rowMeans(vapply(flat, `[[`, logical(length(allTips)), "inCrown"))
pSAchain <- vapply(facts, function(ch) rowMeans(vapply(ch, `[[`, logical(length(allTips)), "sa")),
                   numeric(length(allTips)))
fossilOut <- data.frame(taxon = allTips, p_sampled_ancestor = pSA,
                        p_sa_min_chain = apply(pSAchain, 1, min),
                        p_sa_max_chain = apply(pSAchain, 1, max),
                        p_in_crown = pCrown, row.names = NULL)
fossilOut <- merge(taxaTable, fossilOut, by = "taxon")
fossilOut <- fossilOut[fossilOut$taxon %in% fossils, ]
if (file.exists("published-fossils.csv")) {
  fossilOut <- merge(fossilOut, read.csv("published-fossils.csv"), by = "taxon", all.x = TRUE)
}
fossilOut <- fossilOut[order(-fossilOut$p_sampled_ancestor), ]

# Estimates vs. published ----------------------------------------------------------
published <- read.csv("published.csv")
published <- published[published$config %in% c(clock, "all"), ]
ageRows <- lapply(seq_len(nrow(published)), function(i) {
  p <- published$parameter[i]
  if (!p %in% names(pooled)) return(NULL)
  x <- pooled[[p]]
  h <- HPD(x)
  data.frame(published[i, ], rb_mean = mean(x), rb_median = median(x),
             rb_hpd95_lower = h[1], rb_hpd95_upper = h[2])
})
ages <- do.call(rbind, ageRows)

write.csv(conv, sprintf("summary/convergence_%s.csv", clock), row.names = FALSE)
write.csv(summ, sprintf("summary/parameters_%s.csv", clock), row.names = FALSE)
if (!is.null(ages)) write.csv(ages, sprintf("summary/published_%s.csv", clock), row.names = FALSE)
write.csv(fossilOut, sprintf("summary/fossils_%s.csv", clock), row.names = FALSE)
write.csv(cladeOut, sprintf("summary/clades_%s.csv", clock), row.names = FALSE)
writeLines(c(sprintf("config: %s", clock), sprintf("chains: %s", paste(runNames, collapse = " ")),
             sprintf("samples per chain after burn-in: %d", min(vapply(logs, nrow, 0L))),
             sprintf("ASDSF (clades >= 0.1 in any chain): %.4f over %d clades; max SDSF %.4f",
                     asdsf, sum(keep), max(sdFreq)),
             sprintf("min ESS (per chain): %.0f (%s)", min(conv$ess_min_chain, na.rm = TRUE),
                     conv$parameter[which.min(conv$ess_min_chain)]),
             sprintf("min ESS (pooled): %.0f (%s)", min(conv$ess_total, na.rm = TRUE),
                     conv$parameter[which.min(conv$ess_total)]),
             sprintf("max PSRF: %.4f (%s)", max(conv$psrf, na.rm = TRUE),
                     conv$parameter[which.max(conv$psrf)])),
           sprintf("summary/headline_%s.txt", clock))
cat(readLines(sprintf("summary/headline_%s.txt", clock)), sep = "\n")
