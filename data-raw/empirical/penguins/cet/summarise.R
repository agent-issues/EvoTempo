# Penguin CET pipeline: convergence and posterior summaries of the CET chains
# and of the morphology-only RevBayes comparator, side by side.
#
#   Rscript data-raw/empirical/penguins/cet/summarise.R <CET run root> [<RevBayes root>] [<out dir>]
#
# CET chains: <root>/<label>_run<k>/ (run-chain.R); RevBayes chains:
# <rb root>/<clock>_run<k>/output/penguins_morph_<clock>_run<k>.{log,trees}.
# Drops the first CET_BF_BURN (default 0.25) of every chain.  Writes
# convergence.csv (mean, 95% interval, ESS per chain, total ESS, PSRF),
# clades.csv (clade frequencies, sd over chains), sampled-ancestors.csv and
# headline.txt (crown age, worst ESS and PSRF, ASDSF over clades >= 0.1).
suppressMessages(library(EvoTempo))
source("data-raw/empirical/penguins/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
cetRoot <- args[[1]]
rbRoot <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else NULL
outDir <- if (length(args) >= 3) args[[3]] else file.path(cetRoot, "summary")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
burn <- as.numeric(Sys.getenv("CET_BF_BURN", "0.25"))
tips <- read.csv("data-raw/empirical/penguins/penguins_tips.csv")
extant <- tips$taxon[tips$status == "extant"]
fossils <- tips$taxon[tips$status == "fossil"]

TreeAges <- function(tr) {
  d <- node.depth.edgelength(tr)
  ext <- match(extant, tr$tip.label)
  age <- max(d[ext]) - d
  root <- Ntip(tr) + 1L
  c(crown = age[getMRCA(tr, extant)], root = age[root],
    origin = age[root] + if (is.null(tr$root.edge)) NA else tr$root.edge)
}
SA <- function(tr) {
  len <- setNames(numeric(Ntip(tr)), tr$tip.label)
  e <- tr$edge[, 2] <= Ntip(tr)
  len[tr$edge[e, 2]] <- tr$edge.length[e]
  len[fossils] < 1e-9
}
Row <- function(name, p, x) {
  ess <- vapply(x, Ess, 0)
  u <- unlist(x)
  data.frame(analysis = name, parameter = p, mean = mean(u, na.rm = TRUE),
             q025 = unname(stats::quantile(u, 0.025, na.rm = TRUE)),
             q975 = unname(stats::quantile(u, 0.975, na.rm = TRUE)),
             minEss = min(ess), totalEss = sum(ess), psrf = Psrf(x))
}
# chains: list of list(params = data.frame, trees = multiPhylo)
Analysis <- function(name, chains) {
  conv <- do.call(rbind, lapply(names(chains[[1]]$params), function(p) {
    Row(name, p, lapply(chains, function(ch) ch$params[[p]]))
  }))
  ages <- lapply(chains, function(ch) t(vapply(ch$trees, TreeAges, numeric(3))))
  for (a in c("crown", "root", "origin")) {
    x <- lapply(ages, function(m) m[, a])
    if (!all(is.na(unlist(x)))) conv <- rbind(conv, Row(name, paste0(a, "Age"), x))
  }
  freq <- lapply(chains, function(ch) {
    pp <- prop.part(ch$trees)
    lab <- attr(pp, "labels")
    setNames(attr(pp, "number") / length(ch$trees),
             vapply(pp, function(s) paste(sort(lab[s]), collapse = ","), ""))
  })
  keys <- unique(unlist(lapply(freq, names)))
  fm <- matrix(vapply(freq, function(f) {
    v <- f[keys]
    v[is.na(v)] <- 0
    v
  }, numeric(length(keys))), length(keys))
  sel <- apply(fm, 1, max) >= 0.1
  asdsf <- if (ncol(fm) > 1) mean(apply(fm[sel, , drop = FALSE], 1, stats::sd)) else NA
  clades <- data.frame(analysis = name, clade = keys, frequency = rowMeans(fm),
                       sdOverChains = if (ncol(fm) > 1) apply(fm, 1, stats::sd) else NA)
  sa <- rowMeans(matrix(vapply(chains, function(ch) {
    rowMeans(matrix(vapply(ch$trees, SA, logical(length(fossils))), length(fossils)))
  }, numeric(length(fossils))), length(fossils)))
  list(conv = conv, clades = clades[clades$frequency >= 0.05, ],
       sa = data.frame(analysis = name, fossil = fossils, pSampledAncestor = sa),
       asdsf = asdsf, nChain = length(chains),
       nSample = sum(vapply(chains, function(ch) length(ch$trees), 0)))
}
Keep <- function(n) seq_len(n) > floor(burn * n)

out <- list()
dirs <- list.dirs(cetRoot, recursive = FALSE)
dirs <- dirs[grepl("_run[0-9]+$", basename(dirs))]
for (l in unique(sub("_run[0-9]+$", "", basename(dirs)))) {
  chains <- lapply(dirs[sub("_run[0-9]+$", "", basename(dirs)) == l], function(dd) {
    f <- sort(Sys.glob(file.path(dd, "seg*_samples.csv")))
    if (!length(f)) return(NULL)
    smp <- do.call(rbind, lapply(f, utils::read.csv))
    tr <- unlist(lapply(sort(Sys.glob(file.path(dd, "seg*_trees.nwk"))), readLines))
    keep <- Keep(nrow(smp))
    cols <- intersect(c("lambda", "mu", "psi", "theta", "s", "rate", "logPost", "charLL",
                        "rootBin", "meanX", "nAncestor", "treeLength", "clock"), names(smp))
    list(params = smp[keep, cols], trees = read.tree(text = tr[keep]))
  })
  chains <- Filter(Negate(is.null), chains)
  if (length(chains)) out[[paste("CET", l)]] <- Analysis(paste("CET", l), chains)
}
if (!is.null(rbRoot)) {
  for (clock in c("ucln", "strict")) {
    logs <- Sys.glob(file.path(rbRoot, paste0(clock, "_run*"), "output",
                               paste0("penguins_morph_", clock, "_run*.log")))
    if (!length(logs)) next
    chains <- lapply(logs, function(f) {
      lg <- utils::read.delim(f)
      tt <- utils::read.delim(sub("\\.log$", ".trees", f), stringsAsFactors = FALSE)
      n <- min(nrow(lg), nrow(tt))
      keep <- Keep(n)
      cols <- intersect(c("lambda", "mu", "psi", "clock_morph", "ucld_sd", "Posterior",
                          "Likelihood", "num_samp_anc", "origin_time"), names(lg))
      list(params = lg[seq_len(n)[keep], cols, drop = FALSE],
           trees = read.tree(text = gsub("\\[[^]]*\\]", "", tt[[ncol(tt)]][seq_len(n)[keep]])))
    })
    out[[paste("RevBayes morph", clock)]] <- Analysis(paste("RevBayes morph", clock), chains)
  }
}
utils::write.csv(do.call(rbind, lapply(out, `[[`, "conv")), file.path(outDir, "convergence.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(out, `[[`, "clades")), file.path(outDir, "clades.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(out, `[[`, "sa")),
                 file.path(outDir, "sampled-ancestors.csv"), row.names = FALSE)
headline <- vapply(names(out), function(n) {
  a <- out[[n]]
  cv <- a$conv
  ca <- cv[cv$parameter == "crownAge", ]
  Which <- function(i) if (length(i)) cv$parameter[i] else "-"
  sprintf("%s: %d chains, %d samples; crown %.2f [%.2f, %.2f] Ma; min ESS %.0f (%s); max PSRF %.3f (%s); ASDSF %.3f",
          n, a$nChain, a$nSample, ca$mean, ca$q025, ca$q975,
          suppressWarnings(min(cv$minEss, na.rm = TRUE)), Which(which.min(cv$minEss)),
          suppressWarnings(max(cv$psrf, na.rm = TRUE)), Which(which.max(cv$psrf)), a$asdsf)
}, "")
writeLines(headline, file.path(outDir, "headline.txt"))
writeLines(headline)
