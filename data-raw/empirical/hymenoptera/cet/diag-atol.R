# Hymenoptera CET pipeline: does the TreeLikelihood() atol issue (#90) touch
# check-prior.R?  Its oracle samples (log d, logit turnover, logit sampling)
# with TreeLikelihood() at the default atol (1e-14) as the likelihood.  An
# error that is constant across parameter values cancels in the Metropolis
# ratio; one that varies biases the oracle.
#
#   Rscript data-raw/empirical/hymenoptera/cet/diag-atol.R <prior-check-draws.rds> <out.csv> [n = 200]
#
# Rebuilds check-prior.R's target (11-tip subtree, 16 bins, theta 1/120,
# s 0.5) and evaluates, at n oracle posterior draws per beta,
# TreeLikelihood() at the default atol and at atol = 1e-60, and for beta = 0
# FBDTreeLikelihood() (the dCET likelihood at beta = 0 must equal it, up to
# the constant of the tip-age intervals).  Reports the spread of the
# differences over the draws.
if (nzchar(Sys.getenv("CET_LOAD_ALL"))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(EvoTempo))
}
source("data-raw/empirical/hymenoptera/cet/common.R")
args <- commandArgs(trailingOnly = TRUE)
draws <- readRDS(args[[1]])
outFile <- args[[2]]
n <- if (length(args) >= 3) as.integer(args[[3]]) else 200L
full <- read.tree(file.path(cetDir, "data", "start-tree.nwk"))
tips <- read.csv(file.path(hymDir, "hymenoptera_tips.csv"))
info <- tips[match(full$tip.label, tips$taxon), ]
hymExtant <- full$tip.label[info$status == "extant" & info$group == "Hymenoptera"]
fossils <- full$tip.label[info$status == "fossil"]
keep <- c("Orthoptera", "Paraneoptera", hymExtant[1:6], fossils[1:3])
sub <- drop.tip(full, setdiff(full$tip.label, keep))
si <- tips[match(sub$tip.label, tips$taxon), ]
depth <- node.depth.edgelength(sub)
extant <- si$status == "extant"
age <- max(depth[which(extant)]) - depth
age[which(extant)] <- 0
tree <- TempoTree(sub, ages = age, tipType = ifelse(extant, "extant", "fossil"))
fixed <- list(theta = 1 / 120, s = 0.5)
# check-prior.R saves list(jobs, draws): jobs give kind ("oracle", "cet",
# "control") and beta; draws hold logLambda, turnover and sampling (raw
# scale, as check-prior.R's CET_PC_REUSE path reads them)
kinds <- vapply(draws$jobs, `[[`, "", "kind")
betas <- vapply(draws$jobs, `[[`, 0, "beta")
set.seed(90)
res <- do.call(rbind, lapply(c(0, 1), function(beta) {
  x <- do.call(rbind, draws$draws[kinds == "oracle" & betas == beta])
  m <- DCETModel(nBin = 16, zRange = c(-3, 3), beta = beta, gamma = 0, rho = 5e-4,
                 condition = "extant", start = "root")
  x <- x[sample.int(nrow(x), min(n, nrow(x))), ]
  do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    lambda <- exp(x$logLambda[i])
    mu <- x$turnover[i] * lambda
    p <- c(list(lambda = lambda, mu = mu, psi = x$sampling[i] * mu / (1 - x$sampling[i])), fixed)
    a <- TreeLikelihood(tree, p, m)
    b <- TreeLikelihood(tree, p, m, atol = 1e-60)
    f <- if (beta == 0) FBDTreeLikelihood(tree, p$lambda, p$mu, p$psi, rho = 5e-4,
                                          condition = "extant", start = "root") else NA_real_
    data.frame(set = paste0("oracle, beta ", beta), beta = beta, lambda = p$lambda, mu = p$mu,
               psi = p$psi, defaultAtol = a, atol1e60 = b, fbd = f)
  }))
}))
res$defaultMinusTight <- res$defaultAtol - res$atol1e60
res$tightMinusFbd <- res$atol1e60 - res$fbd
utils::write.csv(res, outFile, row.names = FALSE)
s <- do.call(rbind, lapply(split(res, res$set), function(y) data.frame(
  set = y$set[1], beta = y$beta[1], n = nrow(y),
  defaultMinusTightMin = min(y$defaultMinusTight), defaultMinusTightMax = max(y$defaultMinusTight),
  defaultMinusTightSd = stats::sd(y$defaultMinusTight),
  tightMinusFbdSd = if (all(is.na(y$tightMinusFbd))) NA else stats::sd(y$tightMinusFbd, na.rm = TRUE))))
print(s, digits = 4, row.names = FALSE)
