# Review of #98, R98-13: transition probabilities at the empirical settings.
# Both engines against an exact reference on the trilobite start trees and
# matrix (data-raw/empirical/trilobites/cet) and on the Hymenoptera ordered
# partitions (CET pipeline, branch hymenoptera-cet-pipeline), over clock
# rates 1e-9 to 30, each also scaled by e^-5 and e^5 (gamma = 1 branches at
# the edges of the default grid, z = +-5 with s = 1), with and without
# Gamma rates.
#
#   Rscript data-raw/gamma-acrv-pt-check.R <hymenoptera cet dir> <out.csv>
#
# From the package root of a checkout that includes the trilobite pipeline,
# with EvoTempo installed.  The reference repeats .MkCharactersLogLik()
# (coding = "variable": data columns, constant-pattern columns, mixture over
# categories, correction per mask) but takes P(t) from an independent
# uniformization: a direct Poisson sum without halving or squaring, to a tail
# below 1e-20, or the Taylor series with scaling and squaring once mu t > 200
# (where no transition probability of K <= 16 states is small).  Unordered
# partitions are checked between the engines and for finiteness only
# (closed-form JC in both).
suppressMessages({
  library(EvoTempo)
  library(ape)
})
ns <- asNamespace("EvoTempo")
args <- commandArgs(TRUE)
hymCet <- args[[1]]
outFile <- args[[2]]
nCore <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

MatExp <- function(A) {
  j <- max(0, ceiling(log2(max(abs(A)) * 10)))
  A <- A / 2^j
  out <- term <- diag(nrow(A))
  for (k in 1:30) {
    term <- term %*% A / k
    out <- out + term
  }
  for (i in seq_len(j)) out <- out %*% out
  out
}
RefTransition <- function(Q, t) {
  K <- nrow(Q)
  if (t == 0) return(diag(K))
  mu <- max(-diag(Q))
  x <- mu * t
  if (x > 200) return(pmax(MatExp(Q * t), 0))
  B <- diag(K) + Q / mu
  term <- diag(K)
  P <- stats::dpois(0, x) * term
  for (n in seq_len(ceiling(x + 20 * sqrt(x) + 40 + K))) {
    term <- term %*% B
    w <- stats::dpois(n, x)
    if (w > 0) P <- P + w * term
  }
  P
}
RefColumns <- function(tree, g, lengths) {
  nTip <- tree$nTip
  K <- g$nState
  tip <- array(g$tip, c(K, g$nCol, nTip))
  part <- vector("list", nTip + tree$nNode)
  for (i in seq_len(nTip)) part[[i]] <- matrix(tip[, , i], K)
  logScale <- numeric(g$nCol)
  for (v in tree$postorder) {
    a <- 1
    for (kid in tree$kids[v - nTip, ]) {
      a <- a * (RefTransition(g$Q, lengths[[kid]]) %*% part[[kid]])
    }
    mx <- apply(a, 2, max)
    part[[v]] <- a / rep(mx, each = K)
    logScale <- logScale + log(mx)
  }
  log(colSums(part[[tree$root]]) / K) + logScale
}
RefLogLik <- function(tree, g, lengths, rates) {
  logLik <- matrix(vapply(rates, function(r) RefColumns(tree, g, lengths * r),
                          numeric(g$nCol)), g$nCol)
  top <- apply(logLik, 1, max)
  siteLik <- log(rowMeans(exp(logLik - top))) + top
  data <- seq_len(g$nData)
  logConst <- siteLik[-data]
  logMask <- if (g$ordered) {
    m <- matrix(logConst, g$nConst)
    apply(m, 2, function(col) log(sum(exp(col - max(col)))) + max(col))
  } else {
    log(g$nState) + logConst
  }
  if (any(!(logMask < 0))) return(-Inf)
  sum(siteLik[data]) - sum(g$count * log(-expm1(logMask)))
}
Durations <- function(tt) {
  dur <- numeric(tt$nTip + tt$nNode)
  child <- which(!is.na(tt$parent))
  dur[child] <- tt$age[tt$parent[child]] - tt$age[child]
  dur
}

# ---- Data ----
LoadHym <- function(cetDir) {
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
  rownames(cells) <- tt$tipLabel
  ord <- which(parts$ordered)
  list(list(tree = tt, label = "start", cells = cells[, ord, drop = FALSE],
            nState = parts$nState[ord], ordered = rep(TRUE, length(ord))))
}
LoadTri <- function() {
  triDir <- "data-raw/empirical/trilobites"
  cetDir <- file.path(triDir, "cet")
  windowEnd <- 485.4  # as cet/common.R
  tips <- read.csv(file.path(triDir, "trilobites_tips.csv"))
  raw <- read.csv(file.path(cetDir, "data", "characters.csv"), check.names = FALSE,
                  colClasses = "character")
  parts <- read.csv(file.path(cetDir, "data", "partitions.csv"))
  files <- sort(list.files(file.path(cetDir, "data"), "^start-tree-[0-9]+\\.nwk$",
                           full.names = TRUE))
  lapply(files, function(f) {
    tr <- read.tree(f)
    nTip <- Ntip(tr)
    info <- tips[match(tr$tip.label, tips$taxon), ]
    tipAge <- info$min_age - windowEnd
    depth <- node.depth.edgelength(tr)
    desc <- prop.part(tr)
    age <- c(tipAge, vapply(seq_len(tr$Nnode), function(i) {
      d <- desc[[i]]
      mean(tipAge[d] + depth[d] - depth[nTip + i])
    }, 0))
    len <- numeric(nTip + tr$Nnode)
    len[tr$edge[, 2]] <- tr$edge.length
    parent <- integer(nTip + tr$Nnode)
    parent[tr$edge[, 2]] <- tr$edge[, 1]
    tipType <- ifelse(len[seq_len(nTip)] < 1e-9, "ancestor", "fossil")
    sa <- which(tipType == "ancestor")
    age[parent[sa]] <- tipAge[sa]
    tt <- TempoTree(tr, ages = age, tipType = tipType, ranges = cbind(tipAge, tipAge))
    cells <- as.matrix(raw[match(tt$tipLabel, raw$taxon), -1])
    cells[is.na(cells)] <- ""
    rownames(cells) <- tt$tipLabel
    list(tree = tt, label = basename(f), cells = cells, nState = parts$k,
         ordered = parts$ordered)
  })
}

sets <- c(lapply(LoadHym(hymCet), function(d) c(d, dataset = "hymenoptera")),
          lapply(LoadTri(), function(d) c(d, dataset = "trilobites")))
rates <- 10^seq(-9, 1, by = 1)
rates <- c(rates, 30)
scales <- c(exp(-5), 1, exp(5))
alphas <- c(NA, 1e-3, 1.1, 1e3)
jobs <- expand.grid(set = seq_along(sets), rate = rates, scale = scales, alpha = alphas)
t0 <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  d <- sets[[jobs$set[j]]]
  alpha <- jobs$alpha[j]
  nCat <- if (is.na(alpha)) 1L else 4L
  mk <- MkCharacters(d$cells, nState = d$nState, ordered = d$ordered,
                     coding = "variable", nRateCat = nCat)
  lengths <- jobs$rate[j] * jobs$scale[j] * Durations(d$tree)
  catRates <- if (nCat > 1) GammaRates(alpha, nCat) else 1
  do.call(rbind, lapply(seq_along(mk$partitions), function(p) {
    g <- mk$partitions[[p]]
    one <- mk
    one$partitions <- list(g)
    cpp <- ns$.MkCharactersLogLik(d$tree, one, lengths, alpha, engine = "cpp",
                                  rates = catRates)
    ref <- ns$.MkCharactersLogLik(d$tree, one, lengths, alpha, engine = "R",
                                  rates = catRates)
    exact <- if (g$ordered) RefLogLik(d$tree, g, lengths, catRates) else NA_real_
    data.frame(dataset = d$dataset, tree = d$label, K = g$nState, ordered = g$ordered,
               nChar = g$nData, rate = jobs$rate[j], scale = jobs$scale[j], alpha = alpha,
               cpp = cpp, R = ref, exact = exact)
  }))
}, mc.cores = nCore, mc.preschedule = FALSE))
res$relCppExact <- abs(res$cpp - res$exact) / pmax(1, abs(res$exact))
res$relRExact <- abs(res$R - res$exact) / pmax(1, abs(res$exact))
res$relCppR <- abs(res$cpp - res$R) / pmax(1, abs(res$R))
write.csv(res, outFile, row.names = FALSE)
summ <- do.call(rbind, lapply(split(res, list(res$dataset, res$ordered), drop = TRUE), function(s) {
  fin <- is.finite(s$exact) | !s$ordered
  data.frame(dataset = s$dataset[1], ordered = s$ordered[1], nEval = nrow(s),
             exactNonFinite = sum(s$ordered & !is.finite(s$exact)),
             cppNonFinite = sum(!is.finite(s$cpp)), RNonFinite = sum(!is.finite(s$R)),
             maxRelCppExact = suppressWarnings(max(s$relCppExact[fin], na.rm = TRUE)),
             maxRelRExact = suppressWarnings(max(s$relRExact[fin], na.rm = TRUE)),
             maxRelCppR = max(s$relCppR[is.finite(s$cpp) & is.finite(s$R)]))
}))
print(summ, digits = 3, row.names = FALSE)
bad <- res[!is.finite(res$cpp) | !is.finite(res$R), ]
if (nrow(bad)) {
  cat("Non-finite values (first 20):\n")
  print(utils::head(bad, 20), digits = 4, row.names = FALSE)
}
cat(sprintf("%d evaluations over %d trees in %.1f min\n", nrow(res), length(sets),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
