# Trilobite CET pipeline: does the character-free path of TempoMCMC() (the
# path of check-window.R and check-sbc.R) give the same chains under two
# builds of EvoTempo?  If so, the window check and SBC-lite of PR #91 (run at
# 48e6bf7) stand for the build the chains now use, and need no rerun.
#
#   Rscript data-raw/empirical/trilobites/cet/check-identity.R chain <out.rds>
#   Rscript data-raw/empirical/trilobites/cet/check-identity.R compare <dirA> <dirB>
#
# `chain` runs, with fixed seeds and no characters, prior-only TempoMCMC()
# chains with free lambda, mu, psi, theta and s (parameter sweep, `joint`,
# `shift`), every path and tree move and the origin move, under condition
# "survival" and "none", on a fixed 6-tip all-fossil tree (rho = 0, start =
# "origin"), as check-sbc.R and check-window.R use them.  `compare` requires
# identical() samples and acceptance between the builds, both for these
# chains and for check-window.R's raw output (window-raw-seed*.rds, from a
# small run of check-window.R under each build; check-identity.sh).
suppressMessages(library(EvoTempo))
args <- commandArgs(trailingOnly = TRUE)
if (args[[1]] == "chain") {
  tr <- ape::read.tree(text = "(((a:2,b:3):1.5,(c:1,(d:2.5,e:0.7):1):2):0.5,f:4.2):1.5;")
  # Ages with the present 1 Myr after the youngest tip: all six tips are
  # fossils
  depth <- ape::node.depth.edgelength(tr)
  tt <- TempoTree(tr, ages = max(depth) + 1 - depth, tipType = rep("fossil", 6),
                  ranges = cbind(0, rep(1e3, 6)))
  stopifnot(!any(tt$tipType == "extant"))
  stopifnot(isTRUE(tt$stem > 0))
  out <- list(packageVersion = as.character(utils::packageVersion("EvoTempo")),
              libPath = find.package("EvoTempo"))
  for (cond in c("survival", "none")) {
    model <- DCETModel(nBin = 16, zRange = c(-3, 3), beta = 1, gamma = 1, rho = 0,
                       condition = cond, start = "origin")
    pri <- list(lambda = function(x) stats::dlnorm(x, log(0.5), 0.25, log = TRUE),
                mu = function(x) stats::dlnorm(x, log(0.4), 0.25, log = TRUE),
                psi = function(x) stats::dlnorm(x, log(0.2), 0.3, log = TRUE),
                theta = function(x) stats::dlnorm(x, log(0.2), 0.4, log = TRUE),
                s = function(x) stats::dlnorm(x, log(0.6), 0.25, log = TRUE))
    set.seed(20261116 + (cond == "none"))
    run <- TempoMCMC(tt, model, list(lambda = 0.5, mu = 0.4, psi = 0.2, theta = 0.2, s = 0.6),
                     pri, nIter = 600, burnin = 200, thin = 5, nSlice = 8,
                     moves = c(tree = 1, node = 1, tip = 1, branch = 1, joint = 1, shift = 1,
                               nodeAge = 2, fossilAge = 2, slide = 2, saJump = 2, originAge = 1),
                     agePriors = list(origin = function(a) stats::dunif(a, 3, 20, log = TRUE)))
    out[[cond]] <- list(samples = run$samples, acceptance = run$acceptance,
                        occupancy = run$occupancy)
  }
  saveRDS(out, args[[2]])
  cat(sprintf("chain: EvoTempo at %s; %d + %d samples\n", out$libPath,
              nrow(out$survival$samples), nrow(out$none$samples)))
  quit(save = "no")
}
stopifnot(args[[1]] == "compare")
A <- args[[2]]
B <- args[[3]]
res <- character()
ok <- TRUE
a <- readRDS(file.path(A, "chain.rds"))
b <- readRDS(file.path(B, "chain.rds"))
for (cond in c("survival", "none")) {
  same <- identical(a[[cond]], b[[cond]])
  ok <- ok && same
  res <- c(res, sprintf("free-parameter chain (%s, %d samples x %d columns): identical %s",
                        cond, nrow(a[[cond]]$samples), ncol(a[[cond]]$samples), same))
}
wa <- readRDS(Sys.glob(file.path(A, "window", "window-raw*.rds")))
wb <- readRDS(Sys.glob(file.path(B, "window", "window-raw*.rds")))
for (cond in c("survival", "none")) {
  same <- identical(wa$chains[[cond]], wb$chains[[cond]]) && identical(wa$sims[[cond]], wb$sims[[cond]])
  ok <- ok && same
  res <- c(res, sprintf("check-window.R (%s: %d simulated trees, %d chains x %d samples): identical %s",
                        cond, nrow(wa$sims[[cond]]), length(wa$chains[[cond]]),
                        nrow(wa$chains[[cond]][[1]]$samples), same))
}
res <- c(sprintf("builds: %s vs %s", a$libPath, b$libPath), res,
         if (ok) "check-identity.R: PASS (the character-free path is unchanged)" else
           "check-identity.R: FAIL")
writeLines(res)
if (!ok) quit(status = 1)
