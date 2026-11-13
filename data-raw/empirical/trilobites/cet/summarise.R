# Trilobite CET pipeline: convergence and posterior summaries of the CET
# chains and of the RevBayes comparator (../revbayes, #78), side by side.
#
#   Rscript data-raw/empirical/trilobites/cet/summarise.R <CET run root> [<RevBayes output dir>] [<out dir>]
#
# CET chains: <root>/<label>_run<k>/ (run-chain.R); RevBayes chains:
# <dir>/trilobites_morph_run<k>.{log,trees} (the copies prepare.R reads).
# Drops the first CET_BF_BURN (default 0.25) of every chain.  Convergence is
# MkPrime's (rank-normalised split-chain R-hat over chains and ESS of the
# pooled samples; common.R, RhatEss()), not PSRF.  Writes
# convergence.csv (per analysis and quantity: mean, 95% interval, R-hat, ESS,
# and for CET the mean and interval reweighted to the comparator's joint
# (d, nu, s) prior with the weights' MkPrime bulk ESS, which accounts for
# autocorrelation, and their iid Kish ESS), clades.csv (frequency of each
# clade with frequency >= 0.05, sd over chains), sampled-ancestors.csv and
# headline.txt (root, ingroup and origin ages in Ma, worst R-hat and ESS,
# ASDSF over clades >= 0.1).  No summary tree is made: ConsTree is not
# installed on Hamilton and its Quartet() takes at most 100 tips (107 here).
suppressMessages(library(EvoTempo))
source("data-raw/empirical/trilobites/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
cetRoot <- args[[1]]
rbDir <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else NULL
outDir <- if (length(args) >= 3) args[[3]] else file.path(cetRoot, "summary")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
burn <- as.numeric(Sys.getenv("CET_BF_BURN", "0.25"))
tips <- read.csv(file.path(triDir, "trilobites_tips.csv"))
cons <- read.csv(file.path(triDir, "trilobites_constraints.csv"))
ingroup <- strsplit(cons$taxa[cons$constraint == "Ingroup_node"], " ")[[1]]
priors <- CetPriors()

# Ages in Ma from a tree with branch lengths only: every tip has a known age
TreeAges <- function(tr) {
  d <- node.depth.edgelength(tr)
  tipMa <- tips$min_age[match(tr$tip.label, tips$taxon)]
  nTip <- Ntip(tr)
  rootMa <- mean(tipMa + d[seq_len(nTip)])
  c(root = rootMa, ingroup = rootMa - d[getMRCA(tr, ingroup)],
    origin = if (is.null(tr$root.edge)) NA else rootMa + tr$root.edge)
}
SA <- function(tr) {
  len <- setNames(numeric(Ntip(tr)), tr$tip.label)
  e <- tr$edge[, 2] <= Ntip(tr)
  len[tr$edge[e, 2]] <- tr$edge.length[e]
  len[tips$taxon] < 1e-9
}
Weighted <- function(x, lw) {
  ok <- is.finite(x) & is.finite(lw)
  if (!any(ok)) return(c(NA, NA, NA))
  w <- exp(lw[ok] - max(lw[ok]))
  o <- order(x[ok])
  cw <- cumsum(w[o]) / sum(w)
  c(sum(w * x[ok]) / sum(w), x[ok][o][which(cw >= 0.025)[1]], x[ok][o][which(cw >= 0.975)[1]])
}
# chains: list of list(params = data.frame, trees = multiPhylo, logW = numeric or NULL)
Analysis <- function(name, chains) {
  ages <- lapply(chains, function(ch) t(vapply(ch$trees, TreeAges, numeric(3))))
  vals <- lapply(seq_along(chains), function(i) {
    m <- cbind(as.matrix(chains[[i]]$params), rootMa = ages[[i]][, "root"],
               ingroupMa = ages[[i]][, "ingroup"])
    if (!all(is.na(ages[[i]][, "origin"]))) m <- cbind(m, originMa = ages[[i]][, "origin"])
    m[, colSums(is.finite(m)) == nrow(m), drop = FALSE]
  })
  cols <- Reduce(intersect, lapply(vals, colnames))
  vals <- lapply(vals, function(m) m[, cols, drop = FALSE])
  conv <- RhatEss(vals)
  pooled <- do.call(rbind, vals)
  lw <- if (!is.null(chains[[1]]$logW)) unlist(lapply(chains, `[[`, "logW")) else NULL
  we <- if (is.null(lw)) c(mkprime = NA, kish = NA) else WeightEss(lapply(chains, `[[`, "logW"))
  conv <- do.call(rbind, lapply(cols, function(p) {
    x <- pooled[, p]
    r <- conv[conv$quantity == p, ]
    wq <- if (is.null(lw)) c(NA, NA, NA) else Weighted(x, lw)
    data.frame(analysis = name, parameter = p, mean = mean(x),
               q025 = unname(stats::quantile(x, 0.025)), q975 = unname(stats::quantile(x, 0.975)),
               rhat = if (nrow(r)) r$rhat else NA, ess = if (nrow(r)) r$ess else NA,
               meanJointPrior = wq[1], q025JointPrior = wq[2], q975JointPrior = wq[3],
               weightEss = we[["mkprime"]], weightKish = we[["kish"]])
  }))
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
    rowMeans(matrix(vapply(ch$trees, SA, logical(nrow(tips))), nrow(tips)))
  }, numeric(nrow(tips))), nrow(tips)))
  list(conv = conv, clades = clades[clades$frequency >= 0.05, ],
       sa = data.frame(analysis = name, fossil = tips$taxon, pSampledAncestor = sa),
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
    logW <- vapply(which(keep), function(i) {
      JointPriorLogWeight(as.list(smp[i, c("lambda", "mu", "psi")]), priors)
    }, 0)
    list(params = smp[keep, cols], trees = read.tree(text = tr[keep]), logW = logW)
  })
  chains <- Filter(Negate(is.null), chains)
  if (length(chains)) out[[paste("CET", l)]] <- Analysis(paste("CET", l), chains)
}
if (!is.null(rbDir)) {
  logs <- sort(Sys.glob(file.path(rbDir, "trilobites_morph_run*.log")))
  chains <- lapply(logs, function(f) {
    lg <- utils::read.delim(f)
    tt <- utils::read.delim(sub("\\.log$", ".trees", f), stringsAsFactors = FALSE)
    n <- min(nrow(lg), nrow(tt))
    keep <- seq_len(n)[Keep(n)]
    cols <- intersect(c("lambda", "mu", "psi", "clock_morph", "ucld_sd", "alpha_morph",
                        "Posterior", "Likelihood", "num_samp_anc"), names(lg))
    tr <- read.tree(text = gsub("\\[[^]]*\\]", "", tt[[ncol(tt)]][keep]))
    origin <- lg$origin_time[keep] + windowEnd
    for (i in seq_along(tr)) tr[[i]]$root.edge <- origin[i] - TreeAges(tr[[i]])[["root"]]
    list(params = lg[keep, cols, drop = FALSE], trees = tr, logW = NULL)
  })
  out[["RevBayes morph UCLN"]] <- Analysis("RevBayes morph UCLN", chains)
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
  Age <- function(p) {
    r <- cv[cv$parameter == p, ]
    if (!nrow(r)) return("-")
    sprintf("%.2f [%.2f, %.2f]", r$mean, r$q025, r$q975)
  }
  Which <- function(i) if (length(i)) cv$parameter[i] else "-"
  sprintf("%s: %d chains, %d samples; root %s, ingroup %s, origin %s Ma; min ESS %.0f (%s); max R-hat %.3f (%s); ASDSF %.3f",
          n, a$nChain, a$nSample, Age("rootMa"), Age("ingroupMa"), Age("originMa"),
          suppressWarnings(min(cv$ess, na.rm = TRUE)), Which(which.min(cv$ess)),
          if (all(is.na(cv$rhat))) NA else max(cv$rhat, na.rm = TRUE), Which(which.max(cv$rhat)),
          a$asdsf)
}, "")
writeLines(headline, file.path(outDir, "headline.txt"))
writeLines(headline)
