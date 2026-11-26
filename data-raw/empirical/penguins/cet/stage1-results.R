# Penguin CET pipeline, stage 1 results: convergence and posterior summaries of
# the four (beta, gamma) cells and of the morphology-only RevBayes comparator.
#
#   Rscript data-raw/empirical/penguins/cet/stage1-results.R <CET run root> \
#     <RevBayes root> <out dir> [burn-in = 0.25]
#
# CET chains: <CET run root>/b<beta>g<gamma>_run<rep>/ (run-chain.R).
# RevBayes chains: <RevBayes root>/<clock>_run<k>/output/penguins_morph_<clock>_run<k>.{log,trees}
# (revbayes/run-hamilton.sh).  A chain that is still running is read up to its
# last complete line, and its length is reported.
#
# Burn-in: the first [burn-in] of every chain's recorded samples (CET samples
# start after run-chain.R's own 2000 adaptive burn-in iterations).
#
# Convergence (#94; user decision 2026-11-11): rank-normalised split R-hat
# (max of bulk and tail; Vehtari et al. 2021) and bulk and tail ESS, from
# MkPrime's internal MkPrime:::.Rhat, .EssBulk and .EssTail, on matrices of
# iterations x chains (chains of unequal length are tail-equalised, as
# MkPrime's .ComputeRhat does).  Never on logPost, which is dominated by the
# tempo-path density (#94).  Flags: R-hat > 1.01 or bulk ESS < 400.
# Topology: ASDSF over rooted clades (sampled ancestors are zero-length tips,
# so they are tips here) with frequency >= 0.1 in any chain; per-chain
# posteriors of the crown age, of the extant clades whose ages the RevBayes
# baseline reported (penguins_fbd.Rev, node_2 to node_18), of sampled-ancestor
# and crown-membership probabilities of the crown fossils, and of
# "S. muizoni is a sampled ancestor of the four extant Spheniscus" (defined as
# in the baseline's summarise.R).
# Summary trees: ConsTree::Quartet() consensus of the pooled post-burn-in
# trees (no MCC trees), with the frequency of each consensus split.
#
# Ladder rungs (other b<beta>g<gamma> labels under <CET run root>) enter the
# convergence table and ASDSF only (#107 review, R107-1); meanX and rootBin
# (tempo-path state) are diagnosed in every CET cell (R107-14).
# Segment cap (R107-9): CET_MAX_SEG = m reads segments 1..m of each production
# chain, CET_MAX_SEG_RUNG (default CET_MAX_SEG) of each rung, as
# bayes-factors.R (common.R, SegCap()); the 400 000-iteration results are
# CET_MAX_SEG=40 CET_MAX_SEG_RUNG=20.  CET_MAX_RB_ITER caps the RevBayes chains.
# Comparison with RevBayes (R107-7): RevBayes allows a sampled ancestor at the
# root, EvoTempo does not, so CET is compared with RevBayes both unconditionally
# and given no root sampled ancestor (comparison.csv; chain-mean SEs).
#
# Writes to <out dir>: convergence.csv, asdsf.csv, chains.csv (per-chain
# quantities), posterior.csv (pooled summaries, with P(lambda > mu)),
# comparison.csv, consensus-<analysis>.nwk, consensus-splits.csv,
# samples.csv.gz (post-burn-in draws of the scalar quantities, for the
# figures), versions.txt (R and package versions, md5 of the MkPrime functions
# used; mkprime-functions.R, their deparsed source: R107-11), headline.txt.
# Packages: ape, TreeTools, ConsTree, MkPrime.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
cetRoot <- args[[1]]
rbRoot <- args[[2]]
outDir <- args[[3]]
burn <- if (length(args) >= 4L) as.numeric(args[[4]]) else 0.25
suppressMessages({
  library(ape)
  library(TreeTools)
  library(ConsTree)
})
Rhat <- MkPrime:::.Rhat
EssBulk <- MkPrime:::.EssBulk
EssTail <- MkPrime:::.EssTail
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
tips <- utils::read.csv("data-raw/empirical/penguins/penguins_tips.csv")
extant <- tips$taxon[tips$status == "extant"]
fossils <- tips$taxon[tips$status == "fossil"]
allTips <- sort(tips$taxon)

# Extant clades whose ages the RevBayes baseline reported (penguins_fbd.Rev)
eud <- c("Eudyptes_chrysocome", "Eudyptes_chrysolophus", "Eudyptes_filholi", "Eudyptes_moseleyi",
         "Eudyptes_pachyrhynchus", "Eudyptes_robustus", "Eudyptes_schlegeli", "Eudyptes_sclateri")
sph <- c("Spheniscus_demersus", "Spheniscus_humboldti", "Spheniscus_magellanicus",
         "Spheniscus_mendiculus")
pyg <- c("Pygoscelis_adeliae", "Pygoscelis_antarctica", "Pygoscelis_papua")
apt <- c("Aptenodytes_forsteri", "Aptenodytes_patagonicus")
clades <- list(
  crown = extant,
  node_2 = c("Megadyptes_antipodes", eud, "Eudyptula_minor", sph),
  node_3 = c("Eudyptula_minor", sph),
  node_4 = c(apt, pyg),
  node_5 = pyg,
  node_6 = c("Megadyptes_antipodes", eud),
  node_7 = c("Pygoscelis_antarctica", "Pygoscelis_papua"),
  node_8 = eud,
  node_9 = c("Eudyptes_chrysocome", "Eudyptes_filholi", "Eudyptes_moseleyi",
             "Eudyptes_chrysolophus", "Eudyptes_schlegeli"),
  node_10 = c("Eudyptes_pachyrhynchus", "Eudyptes_robustus", "Eudyptes_sclateri"),
  node_11 = sph,
  node_12 = apt,
  node_13 = c("Spheniscus_humboldti", "Spheniscus_mendiculus"),
  node_14 = c("Spheniscus_demersus", "Spheniscus_magellanicus"),
  node_15 = c("Eudyptes_chrysocome", "Eudyptes_filholi", "Eudyptes_moseleyi"),
  node_16 = c("Eudyptes_pachyrhynchus", "Eudyptes_robustus"),
  node_17 = c("Eudyptes_chrysolophus", "Eudyptes_schlegeli"),
  node_18 = c("Eudyptes_chrysocome", "Eudyptes_filholi"))
stopifnot(all(unlist(clades) %in% extant))
# Fossils that the source analysis places in or near the crown
crownFossils <- c("Spheniscus_muizoni", "Spheniscus_urbinai", "Spheniscus_megaramphus",
                  "Pygoscelis_grandis", "Madrynornis_mirandus")
stopifnot(all(crownFossils %in% fossils))

# ---- Per-tree quantities --------------------------------------------------------
TreeFacts <- function(tr) {
  nTip <- length(tr$tip.label)
  d <- node.depth.edgelength(tr)
  ext <- match(extant, tr$tip.label)
  age <- max(d[ext]) - d
  tipEdge <- match(seq_len(nTip), tr$edge[, 2])
  sa <- tr$edge.length[tipEdge] < 1e-8 & tr$tip.label %in% fossils
  # A sampled ancestor at the root (zero-length tip whose parent is the root):
  # RevBayes allows it; EvoTempo's samplers do not (.NoRootAncestor)
  rootSA <- any(sa & tr$edge[tipEdge, 1] == nTip + 1L)
  pp <- prop.part(tr)
  sizes <- lengths(pp)
  M <- matrix(FALSE, length(pp), nTip)
  for (i in seq_along(pp)) M[i, pp[[i]]] <- TRUE
  Mrca <- function(set) {
    idx <- match(set, tr$tip.label)
    ok <- which(rowSums(M[, idx, drop = FALSE]) == length(idx))
    ok[which.min(sizes[ok])]
  }
  cladeOut <- unlist(lapply(names(clades), function(nm) {
    i <- Mrca(clades[[nm]])
    below <- tr$tip.label[M[i, ]]
    stats::setNames(c(age[nTip + i], setequal(intersect(below, extant), clades[[nm]])),
                    paste0(nm, c("_age", "_mono")))
  }))
  crownNode <- Mrca(extant)
  inCrown <- M[crownNode, match(crownFossils, tr$tip.label)]
  isSA <- sa[match(crownFossils, tr$tip.label)]
  muz <- match("Spheniscus_muizoni", tr$tip.label)
  muzAnc <- FALSE
  if (sa[muz]) {
    attach <- tr$edge[tipEdge[muz], 1]
    muzAnc <- all(M[attach - nTip, match(sph, tr$tip.label)])
  }
  g <- match(tr$tip.label, allTips)
  keys <- vapply(pp[sizes > 1 & sizes < nTip], function(p) paste(sort(g[p]), collapse = "."), "")
  list(num = c(crownAge = age[nTip + crownNode], rootAge = age[nTip + 1L],
               originAge = age[nTip + 1L] + if (is.null(tr$root.edge)) NA else tr$root.edge,
               nSA = sum(sa), rootSA = rootSA, treeLengthTips = sum(tr$edge.length), muzAnc = muzAnc,
               stats::setNames(inCrown, paste0("inCrown_", crownFossils)),
               stats::setNames(isSA, paste0("SA_", crownFossils)), cladeOut),
       keys = keys)
}
Strip <- function(x) gsub("\\[[^]]*\\]", "", x)
ParseTrees <- function(nwk) {
  tr <- read.tree(text = Strip(nwk))
  if (inherits(tr, "phylo")) tr <- structure(list(tr), class = "multiPhylo")
  tr
}

# ---- Load chains ----------------------------------------------------------------
Burn <- function(n) seq_len(n) > floor(burn * n)
productionCells <- c("b0g0", "b0g1", "b1g0", "b1g1")
SegCap <- function(label) {
  m <- suppressWarnings(as.integer(Sys.getenv("CET_MAX_SEG", "0")))
  if (!sub("_run[0-9]+$", "", label) %in% productionCells && nzchar(Sys.getenv("CET_MAX_SEG_RUNG"))) {
    m <- suppressWarnings(as.integer(Sys.getenv("CET_MAX_SEG_RUNG")))
  }
  if (is.na(m)) stop("CET_MAX_SEG and CET_MAX_SEG_RUNG must be whole numbers")
  m
}
SegFiles <- function(dir, what) {
  f <- sort(Sys.glob(file.path(dir, paste0("seg*_", what))))
  m <- SegCap(basename(dir))
  if (m > 0L) {
    if (length(f) < m) stop(sprintf("%s has %d segments, fewer than the cap %d", dir, length(f), m))
    f <- utils::head(f, m)
  }
  f
}
LoadCet <- function(dir) {
  smp <- do.call(rbind, lapply(SegFiles(dir, "samples.csv"), utils::read.csv))
  nwk <- unlist(lapply(SegFiles(dir, "trees.nwk"), readLines))
  stopifnot(length(nwk) == nrow(smp))
  keep <- Burn(nrow(smp))
  list(params = smp[keep, ], nwk = nwk[keep], n = nrow(smp), lastIteration = max(smp$iteration))
}
LoadRb <- function(dir) {
  stem <- file.path(dir, "output", paste0("penguins_morph_", basename(dir)))
  Read <- function(f) {
    x <- readLines(f)
    hdr <- strsplit(x[1], "\t")[[1]]
    body <- x[-1]
    nf <- lengths(strsplit(body, "\t"))
    body <- body[nf >= length(hdr) - 1L]  # drop an incomplete last line
    utils::read.delim(text = c(x[1], body), quote = "", comment.char = "",
                      check.names = FALSE, stringsAsFactors = FALSE)
  }
  lg <- Read(paste0(stem, ".log"))
  tr <- Read(paste0(stem, ".trees"))
  tr <- tr[grepl(";\\s*$", tr$fbd_tree), ]
  it <- intersect(lg$Iteration, tr$Iteration)
  maxIt <- as.numeric(Sys.getenv("CET_MAX_RB_ITER", "Inf"))
  it <- it[it <= maxIt]
  lg <- lg[match(it, lg$Iteration), ]
  tr <- tr[match(it, tr$Iteration), ]
  keep <- Burn(nrow(lg))
  p <- lg[keep, ]
  params <- data.frame(lambda = p$lambda, mu = p$mu, psi = p$psi, rate = p$clock_morph,
                       originAge = p$origin_time, rootAge = p$root_age, crownAgeLog = p$crown_age,
                       nSAlog = p$num_samp_anc, charLL = p$Likelihood, iteration = p$Iteration)
  if ("ucld_sd" %in% names(p)) params$ucldSd <- p$ucld_sd
  list(params = params, nwk = tr$fbd_tree[keep], n = nrow(lg), lastIteration = max(lg$Iteration))
}
analyses <- list()
for (cell in c("b0g0", "b0g1", "b1g0", "b1g1")) {
  dirs <- sort(Sys.glob(file.path(cetRoot, paste0(cell, "_run*"))))
  dirs <- dirs[grepl("_run[0-9]+$", dirs)]
  analyses[[cell]] <- stats::setNames(lapply(dirs, LoadCet), basename(dirs))
}
for (clock in c("strict", "ucln")) {
  dirs <- sort(Sys.glob(file.path(rbRoot, paste0(clock, "_run*"))))
  if (length(dirs)) analyses[[paste0("revbayes_", clock)]] <- stats::setNames(lapply(dirs, LoadRb), basename(dirs))
}
# Ladder rungs: convergence and ASDSF only (R107-1)
rungs <- list()
rungDirs <- list.dirs(cetRoot, recursive = FALSE)
rungDirs <- rungDirs[grepl("^b[0-9.]+g[0-9.]+_run[0-9]+$", basename(rungDirs)) &
                       !sub("_run[0-9]+$", "", basename(rungDirs)) %in% productionCells]
for (cell in sort(unique(sub("_run[0-9]+$", "", basename(rungDirs))))) {
  dirs <- sort(rungDirs[sub("_run[0-9]+$", "", basename(rungDirs)) == cell])
  rungs[[cell]] <- stats::setNames(lapply(dirs, LoadCet), basename(dirs))
}
allSets <- c(analyses, rungs)
lengthsInfo <- do.call(rbind, lapply(names(allSets), function(a) {
  data.frame(analysis = a, chain = names(allSets[[a]]),
             samples = vapply(allSets[[a]], `[[`, 0L, "n"),
             lastIteration = vapply(allSets[[a]], `[[`, 0, "lastIteration"),
             keptAfterBurnin = vapply(allSets[[a]], function(ch) nrow(ch$params), 0L))
}))
utils::write.csv(lengthsInfo, file.path(outDir, "chain-lengths.csv"), row.names = FALSE)
print(lengthsInfo, row.names = FALSE)

# ---- Tree facts --------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
for (a in names(allSets)) {
  for (k in names(allSets[[a]])) {
    trees <- ParseTrees(allSets[[a]][[k]]$nwk)
    facts <- lapply(trees, TreeFacts)
    allSets[[a]][[k]]$facts <- as.data.frame(do.call(rbind, lapply(facts, `[[`, "num")))
    allSets[[a]][[k]]$keys <- lapply(facts, `[[`, "keys")
    # Rung trees are needed only for their facts and clade keys
    if (a %in% names(analyses)) allSets[[a]][[k]]$trees <- trees
  }
  cat(sprintf("%s: tree facts done (%.0f s)\n", a, proc.time()[["elapsed"]] - t0))
}

# Scalar quantities per draw (CET and RevBayes on the same footing)
Scalars <- function(ch, a) {
  p <- ch$params
  f <- ch$facts
  out <- data.frame(lambda = p$lambda, mu = p$mu, psi = p$psi, rate = p$rate,
                    originAge = p$originAge, rootAge = p$rootAge, crownAge = f$crownAge,
                    nSA = f$nSA, charLL = p$charLL)
  if (startsWith(a, "b")) {
    out$theta <- p$theta
    out$s <- p$s
    out$treeLength <- p$treeLength
    out$meanX <- p$meanX
    out$rootBin <- p$rootBin
  } else {
    out$treeLength <- f$treeLengthTips
    if (!is.null(p$ucldSd)) out$ucldSd <- p$ucldSd
  }
  out$diversification <- out$lambda - out$mu
  out$turnover <- out$mu / out$lambda
  out$samplingProportion <- out$psi / (out$mu + out$psi)
  out
}
for (a in names(allSets)) for (k in names(allSets[[a]])) {
  allSets[[a]][[k]]$scalars <- Scalars(allSets[[a]][[k]], a)
}
analyses <- allSets[names(analyses)]
rungs <- allSets[names(rungs)]
# Consistency checks between tree-derived and logged quantities
checks <- do.call(rbind, lapply(names(analyses), function(a) do.call(rbind, lapply(names(analyses[[a]]), function(k) {
  ch <- analyses[[a]][[k]]
  data.frame(analysis = a, chain = k,
             maxAbsRootAgeDiff = max(abs(ch$facts$rootAge - ch$params$rootAge)),
             # RevBayes tree files carry no origin (their root edge is not the stem)
             maxAbsOriginAgeDiff = if (startsWith(a, "b")) max(abs(ch$facts$originAge - ch$params$originAge)) else NA,
             nSAmismatch = if (startsWith(a, "b")) sum(ch$facts$nSA != ch$params$nAncestor) else
               sum(ch$facts$nSA != ch$params$nSAlog),
             maxAbsCrownAgeDiff = if (startsWith(a, "b")) NA else max(abs(ch$facts$crownAge - ch$params$crownAgeLog)))
}))))
print(checks, row.names = FALSE)
utils::write.csv(checks, file.path(outDir, "consistency-checks.csv"), row.names = FALSE)

# ---- Convergence ------------------------------------------------------------------
convParams <- c("lambda", "mu", "psi", "theta", "s", "rate", "originAge", "crownAge", "rootAge",
                "nSA", "charLL", "treeLength", "ucldSd", "meanX", "rootBin")
ChainMatrix <- function(a, p) {
  xs <- lapply(allSets[[a]], function(ch) ch$scalars[[p]])
  if (any(vapply(xs, is.null, TRUE))) return(NULL)
  n <- min(lengths(xs))
  do.call(cbind, lapply(xs, function(x) utils::tail(x, n)))
}
conv <- do.call(rbind, lapply(names(allSets), function(a) do.call(rbind, lapply(convParams, function(p) {
  m <- ChainMatrix(a, p)
  if (is.null(m)) return(NULL)
  data.frame(analysis = a, role = if (a %in% names(rungs)) "rung" else "cell", parameter = p,
             chains = ncol(m), drawsPerChain = nrow(m),
             rhat = Rhat(m), essBulk = EssBulk(m), essTail = EssTail(m))
}))))
conv$flag <- ifelse(conv$rhat > 1.01 | conv$essBulk < 400,
                    paste0(ifelse(conv$rhat > 1.01, "R-hat>1.01 ", ""),
                           ifelse(conv$essBulk < 400, "bulkESS<400", "")), "")
print(conv, digits = 4, row.names = FALSE)
utils::write.csv(conv, file.path(outDir, "convergence.csv"), row.names = FALSE)

# ---- ASDSF over rooted clades ------------------------------------------------------
asdsf <- do.call(rbind, lapply(names(allSets), function(a) {
  freq <- lapply(allSets[[a]], function(ch) {
    tab <- table(unlist(ch$keys))
    tab / length(ch$keys)
  })
  keys <- unique(unlist(lapply(freq, function(f) names(f)[f >= 0.1])))
  fm <- vapply(freq, function(f) { v <- as.numeric(f[keys]); v[is.na(v)] <- 0; v }, numeric(length(keys)))
  fm <- matrix(fm, nrow = length(keys))
  sds <- apply(fm, 1, stats::sd)
  data.frame(analysis = a, clades = length(keys), asdsf = mean(sds), maxSdsf = max(sds),
             cladesSdsfOver0.1 = sum(sds > 0.1))
}))
print(asdsf, digits = 4, row.names = FALSE)
utils::write.csv(asdsf, file.path(outDir, "asdsf.csv"), row.names = FALSE)

# ---- Per-chain quantities -----------------------------------------------------------
Q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
chains <- do.call(rbind, lapply(names(analyses), function(a) do.call(rbind, lapply(names(analyses[[a]]), function(k) {
  ch <- analyses[[a]][[k]]
  f <- ch$facts
  s <- ch$scalars
  rows <- list(
    data.frame(quantity = "crownAge", mean = mean(s$crownAge), median = Q(s$crownAge, 0.5),
               q025 = Q(s$crownAge, 0.025), q975 = Q(s$crownAge, 0.975)),
    data.frame(quantity = "rootAge", mean = mean(s$rootAge), median = Q(s$rootAge, 0.5),
               q025 = Q(s$rootAge, 0.025), q975 = Q(s$rootAge, 0.975)),
    data.frame(quantity = "nSA", mean = mean(s$nSA), median = Q(s$nSA, 0.5),
               q025 = Q(s$nSA, 0.025), q975 = Q(s$nSA, 0.975)),
    data.frame(quantity = "P(muizoni SA of extant Spheniscus)", mean = mean(f$muzAnc),
               median = NA, q025 = NA, q975 = NA),
    data.frame(quantity = "P(sampled ancestor at the root)", mean = mean(f$rootSA),
               median = NA, q025 = NA, q975 = NA),
    data.frame(quantity = "rootAge given no root sampled ancestor", mean = mean(s$rootAge[!f$rootSA]),
               median = Q(s$rootAge[!f$rootSA], 0.5), q025 = Q(s$rootAge[!f$rootSA], 0.025),
               q975 = Q(s$rootAge[!f$rootSA], 0.975)))
  for (nm in names(clades)[-1]) {
    rows[[length(rows) + 1]] <- data.frame(quantity = paste0("P(", nm, " monophyletic among extant)"),
                                           mean = mean(f[[paste0(nm, "_mono")]]), median = NA,
                                           q025 = NA, q975 = NA)
    x <- f[[paste0(nm, "_age")]]
    rows[[length(rows) + 1]] <- data.frame(quantity = paste0(nm, " tmrca age"), mean = mean(x),
                                           median = Q(x, 0.5), q025 = Q(x, 0.025), q975 = Q(x, 0.975))
  }
  for (fo in crownFossils) {
    rows[[length(rows) + 1]] <- data.frame(quantity = paste0("P(", fo, " in crown)"),
                                           mean = mean(f[[paste0("inCrown_", fo)]]), median = NA, q025 = NA, q975 = NA)
    rows[[length(rows) + 1]] <- data.frame(quantity = paste0("P(", fo, " sampled ancestor)"),
                                           mean = mean(f[[paste0("SA_", fo)]]), median = NA, q025 = NA, q975 = NA)
  }
  cbind(analysis = a, chain = k, do.call(rbind, rows))
}))))
utils::write.csv(chains, file.path(outDir, "chains.csv"), row.names = FALSE)

# ---- Pooled posterior summaries -----------------------------------------------------
sumParams <- c("crownAge", "rootAge", "originAge", "nSA", "lambda", "mu", "psi", "diversification",
               "turnover", "samplingProportion", "rate", "theta", "s", "ucldSd", "charLL", "treeLength")
posterior <- do.call(rbind, lapply(names(analyses), function(a) {
  pool <- do.call(rbind, lapply(names(analyses[[a]]), function(k) cbind(chain = k, analyses[[a]][[k]]$scalars)))
  do.call(rbind, lapply(intersect(sumParams, names(pool)), function(p) {
    x <- pool[[p]]
    data.frame(analysis = a, parameter = p, mean = mean(x), median = Q(x, 0.5), q025 = Q(x, 0.025),
               q975 = Q(x, 0.975), draws = length(x))
  }))
}))
# P(lambda > mu) (R107-10), as a row with the probability in mean and median
posterior <- rbind(posterior, do.call(rbind, lapply(names(analyses), function(a) {
  pool <- do.call(rbind, lapply(analyses[[a]], `[[`, "scalars"))
  p <- mean(pool$lambda > pool$mu)
  data.frame(analysis = a, parameter = "P(lambda > mu)", mean = p, median = p, q025 = NA, q975 = NA,
             draws = nrow(pool))
})))
utils::write.csv(posterior, file.path(outDir, "posterior.csv"), row.names = FALSE)

# ---- CET against RevBayes, unconditionally and given no root sampled ancestor (R107-7) ----
# Chain means (one per chain) of each quantity; SE = sd(chain means) / sqrt(chains),
# so z has few degrees of freedom (4 chains a side).  KS D on the pooled draws.
cmpQ <- c("lambda", "mu", "psi", "rate", "crownAge", "rootAge", "originAge", "nSA", "charLL", "lambdaGtMu")
ChainMeans <- function(a, q, noRootSA) {
  vapply(analyses[[a]], function(ch) {
    s <- ch$scalars
    s$lambdaGtMu <- as.numeric(s$lambda > s$mu)
    keep <- if (noRootSA) !as.logical(ch$facts$rootSA) else rep(TRUE, nrow(s))
    mean(s[[q]][keep])
  }, 0)
}
Pooled <- function(a, q, noRootSA) {
  unlist(lapply(analyses[[a]], function(ch) {
    s <- ch$scalars
    s$lambdaGtMu <- as.numeric(s$lambda > s$mu)
    keep <- if (noRootSA) !as.logical(ch$facts$rootSA) else rep(TRUE, nrow(s))
    s[[q]][keep]
  }))
}
pairs <- list(c("b0g0", "revbayes_strict"), c("b0g1", "revbayes_ucln"), c("b1g1", "revbayes_ucln"))
pairs <- Filter(function(p) all(p %in% names(analyses)), pairs)
comparison <- do.call(rbind, lapply(pairs, function(pr) do.call(rbind, lapply(c(FALSE, TRUE), function(cond) {
  do.call(rbind, lapply(cmpQ, function(q) {
    m1 <- ChainMeans(pr[1], q, FALSE)  # CET has no root sampled ancestor
    m2 <- ChainMeans(pr[2], q, cond)
    se <- sqrt(stats::var(m1) / length(m1) + stats::var(m2) / length(m2))
    data.frame(cet = pr[1], revbayes = pr[2], revbayesDraws = if (cond) "no root SA" else "all",
               quantity = if (q == "lambdaGtMu") "P(lambda > mu)" else q,
               cetMean = mean(m1), revbayesMean = mean(m2), diff = mean(m1) - mean(m2), se = se,
               z = (mean(m1) - mean(m2)) / se,
               ksD = if (q == "lambdaGtMu") NA else
                 suppressWarnings(stats::ks.test(Pooled(pr[1], q, FALSE), Pooled(pr[2], q, cond))$statistic[[1]]))
  }))
}))))
comparison <- rbind(comparison, do.call(rbind, lapply(intersect(c("revbayes_strict", "revbayes_ucln"), names(analyses)), function(a) {
  m <- vapply(analyses[[a]], function(ch) mean(as.logical(ch$facts$rootSA)), 0)
  data.frame(cet = NA, revbayes = a, revbayesDraws = "all", quantity = "P(root sampled ancestor)",
             cetMean = 0, revbayesMean = mean(m), diff = -mean(m), se = stats::sd(m) / sqrt(length(m)),
             z = NA, ksD = NA)
})))
print(comparison, digits = 4, row.names = FALSE)
utils::write.csv(comparison, file.path(outDir, "comparison.csv"), row.names = FALSE)
samples <- do.call(rbind, lapply(names(analyses), function(a) do.call(rbind, lapply(names(analyses[[a]]), function(k) {
  s <- analyses[[a]][[k]]$scalars
  for (p in setdiff(c("theta", "s", "ucldSd"), names(s))) s[[p]] <- NA_real_
  cbind(analysis = a, chain = k, draw = seq_len(nrow(s)), s[, c(sumParams)])
}))))
gz <- gzfile(file.path(outDir, "samples.csv.gz"), "w")
utils::write.csv(samples, gz, row.names = FALSE)
close(gz)

# ---- Quartet consensus ------------------------------------------------------------
splitRows <- list()
for (a in names(analyses)) {
  trees <- do.call(c, lapply(analyses[[a]], `[[`, "trees"))
  trees <- structure(lapply(trees, RenumberTips, allTips), class = "multiPhylo")
  t1 <- proc.time()[["elapsed"]]
  cons <- RenumberTips(ConsTree::Quartet(trees), allTips)
  freq <- SplitFrequency(cons, trees) / length(trees)
  lm <- as.logical(as.Splits(cons, allTips))
  lm <- matrix(lm, ncol = length(allTips))
  stopifnot(nrow(lm) == length(freq))
  side <- vapply(seq_len(nrow(lm)), function(i) {
    inSplit <- lm[i, ]
    lab <- allTips[inSplit]
    if (length(lab) > length(allTips) / 2) lab <- allTips[!inSplit]
    paste(sort(lab), collapse = " ")
  }, "")
  splitRows[[a]] <- data.frame(analysis = a, split = side, frequency = as.numeric(freq))
  write.tree(cons, file.path(outDir, sprintf("consensus-%s.nwk", a)))
  cat(sprintf("%s: quartet consensus of %d trees, %d of %d splits, %.0f s\n", a, length(trees),
              NSplits(cons), length(allTips) - 3L, proc.time()[["elapsed"]] - t1))
}
utils::write.csv(do.call(rbind, splitRows), file.path(outDir, "consensus-splits.csv"), row.names = FALSE)

# ---- Headline ----------------------------------------------------------------------------
P <- function(a, p) posterior[posterior$analysis == a & posterior$parameter == p, ]
head <- c(sprintf("burn-in: first %.0f%% of each chain's recorded samples", 100 * burn),
          sprintf("%s: chains %d, kept draws per chain %s", names(analyses),
                  vapply(analyses, length, 0L),
                  vapply(analyses, function(x) paste(range(vapply(x, function(ch) nrow(ch$params), 0L)), collapse = "-"), "")),
          vapply(names(allSets), function(a) {
            cv <- conv[conv$analysis == a, ]
            sprintf("%s: max R-hat %.4f (%s); min bulk ESS %.1f (%s); min tail ESS %.0f (%s); flagged: %s", a,
                    max(cv$rhat), cv$parameter[which.max(cv$rhat)], min(cv$essBulk),
                    cv$parameter[which.min(cv$essBulk)], min(cv$essTail), cv$parameter[which.min(cv$essTail)],
                    if (any(nzchar(cv$flag))) paste(cv$parameter[nzchar(cv$flag)], collapse = ", ") else "none")
          }, ""),
          sprintf("%s: ASDSF %.4f over %d clades; max SDSF %.4f", asdsf$analysis, asdsf$asdsf,
                  asdsf$clades, asdsf$maxSdsf),
          vapply(names(analyses), function(a) {
            x <- P(a, "crownAge")
            sprintf("%s: crown age median %.2f (95%% %.2f-%.2f)", a, x$median, x$q025, x$q975)
          }, ""))
writeLines(head, file.path(outDir, "headline.txt"))

# ---- Versions (R107-11) ------------------------------------------------------------
# MkPrime is installed without a RemoteSha, so record the md5 of the deparsed
# functions used, and their source
pkgs <- c("EvoTempo", "ape", "TreeTools", "ConsTree", "Quartet", "MkPrime")
ver <- vapply(pkgs, function(p) {
  d <- tryCatch(utils::packageDescription(p), error = function(e) NULL)
  if (is.null(d) || !is.list(d)) return(sprintf("%s: not installed", p))
  sprintf("%s %s%s (%s)", p, d$Version, if (!is.null(d$RemoteSha)) paste0(" RemoteSha ", d$RemoteSha) else "",
          dirname(dirname(attr(d, "file"))))
}, "")
fns <- c(".Rhat", ".EssBulk", ".EssTail")
src <- unlist(lapply(fns, function(f) c(paste0("MkPrime:::", f, " <- "), deparse(get(f, asNamespace("MkPrime"))), "")))
writeLines(src, file.path(outDir, "mkprime-functions.R"))
md5 <- vapply(fns, function(f) {
  tmp <- tempfile(); on.exit(unlink(tmp))
  writeLines(deparse(get(f, asNamespace("MkPrime"))), tmp)
  unname(tools::md5sum(tmp))
}, "")
evoSha <- file.path(Sys.getenv("W", "/nobackup/pjjg18/evotempo-penguin-cet"), "out", "build-sha.txt")
writeLines(c(R.version.string, ver, sprintf("MkPrime:::%s deparsed md5 %s", fns, md5),
             if (file.exists(evoSha)) paste("EvoTempo build commit", readLines(evoSha, n = 1)),
             sprintf("CET_MAX_SEG=%s CET_MAX_SEG_RUNG=%s CET_MAX_RB_ITER=%s", Sys.getenv("CET_MAX_SEG"),
                     Sys.getenv("CET_MAX_SEG_RUNG"), Sys.getenv("CET_MAX_RB_ITER"))),
           file.path(outDir, "versions.txt"))
cat(head, sep = "\n")
