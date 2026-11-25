# Validation of the likelihood-weighted regraft `wfnpr` (#110;
# dev/spec/weighted-regraft.md) at fixed node ages and parameters, where the
# topology has a known distribution.
#
# With gamma = 0 (characters independent of tempo) the topology has probability
# proportional to TreeLikelihood() times the Mk likelihood.  Every topology that
# fixed-height regrafts reach is enumerated, and independent chains are
# compared with the exact probabilities of (a) each topology of probability
# above 0.02 and (b) each clade.  `fnpr` chains run alongside as a control.
# With gamma = 1 there is no exact oracle: `wfnpr` against `fnpr`.
#
# Settings (environment):
#   WFNPR_TREE    "sa5": five tips, two fossils and a sampled ancestor (18
#                 topologies, at most two targets); "extant5": five extant tips
#                 (180 topologies, up to three targets, so that WFNPR_MAX = 3
#                 draws subsets of two candidates and WFNPR_MAX = 2 of one)
#   WFNPR_CASE    "prior" (no characters), "characters" (gamma = 0) or
#                 "tempoCharacters" (gamma = 1)
#   WFNPR_MAX     .wfnprMax, the cap on candidates including the current one
#                 (default: the package's, 8, which never truncates here)
#   WFNPR_CHAINS, WFNPR_ITER   chains per move and iterations per chain (the
#                 first tenth is discarded)
#   WFNPR_SEED    seed offset;  WFNPR_CORES  parallel workers (not on Windows)
# Run from the package root with the package installed.  Results, in
# data-raw/wfnpr-results/: <label>.csv (per statistic) and <label>-chains.csv
# (per chain, for pooling with wfnpr-pool.R), with
# <label> = <tree>-<case>-max<cap>-<seed>.
suppressMessages(library("EvoTempo"))
Env <- function(name, default) Sys.getenv(name, default)
treeName <- Env("WFNPR_TREE", "sa5")
case <- Env("WFNPR_CASE", "characters")
nChain <- as.numeric(Env("WFNPR_CHAINS", "32"))
nIter <- as.numeric(Env("WFNPR_ITER", "20000"))
offset <- as.numeric(Env("WFNPR_SEED", "0"))
nCore <- as.integer(Env("WFNPR_CORES", "1"))
ns <- asNamespace("EvoTempo")
for (f in c(".PathNodes", ".SlideNodes", ".FnprTargets", ".Descendants", ".Rebuild",
            ".AddSlices", ".Precompute", ".Messages", ".SampleTree", ".PathMove",
            ".TreeMove", ".LogSumExp", ".CharacterLengths")) {
  assign(f, get(f, ns))
}
maxCandidates <- as.integer(Env("WFNPR_MAX", as.character(get(".wfnprMax", ns))))
utils::assignInNamespace(".wfnprMax", maxCandidates, ns = "EvoTempo")
stopifnot(treeName %in% c("sa5", "extant5"),
          case %in% c("prior", "characters", "tempoCharacters"))
label <- sprintf("%s-%s-max%d-%d", treeName, case, maxCandidates, offset)

Phylo <- function(edge, length, tips) {
  structure(list(edge = edge, edge.length = length, tip.label = tips,
                 Nnode = length(tips) - 1L), class = "phylo")
}
if (treeName == "sa5") {
  # ((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7); a, c fossils; f sampled ancestor
  phylo <- Phylo(
    rbind(c(6, 7), c(7, 1), c(7, 8), c(8, 2), c(8, 3), c(6, 9), c(9, 4), c(9, 5)),
    c(1, 2, 1.5, 1.5, 0, 0.7, 1.2, 3.3), c("a", "b", "f", "c", "d")
  )
  chars <- cbind(c(0, 0, 0, 1, 1), c(0, 1, 1, 0, 0), c(1, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
} else {
  # (((a:1,b:1):1.5,c:2.5):1.5,(d:2,e:2):2);
  phylo <- Phylo(
    rbind(c(6, 7), c(7, 8), c(8, 1), c(8, 2), c(7, 3), c(6, 9), c(9, 4), c(9, 5)),
    c(1.5, 1.5, 1, 1, 2.5, 2, 2, 2), c("a", "b", "c", "d", "e")
  )
  chars <- cbind(c(0, 0, 0, 1, 1), c(1, 1, 0, 0, 0), c(1, 1, 1, 0, 0), c(0, 0, 1, 1, 1),
                 c(0, 1, 1, 0, 0))
}
start <- TempoTree(phylo)
par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.15)
gamma <- if (case == "tempoCharacters") 1 else 0
model <- DCETModel(nBin = 7, zRange = c(-2, 2), gamma = gamma, rho = 0.7)
CharLL <- function(tr, paths, p, s) {
  if (case == "prior") {
    return(0)
  }
  if (gamma == 0) {
    nodes <- .PathNodes(tr)
    len <- numeric(length(tr$parent))
    len[nodes] <- p$rate * (tr$age[tr$parent[nodes]] - tr$age[nodes])
    return(MkLikelihood(tr, chars, len))
  }
  MkLikelihood(tr, chars, .CharacterLengths(tr, s$grid$x, paths, p$rate, gamma))
}

Key <- function(tr) paste(tr$parent, collapse = " ")
trees <- list(start)
names(trees) <- Key(start)
todo <- trees
while (length(todo)) {
  tr <- todo[[1]]
  todo <- todo[-1]
  for (v in .SlideNodes(tr)) {
    p <- tr$parent[[v]]
    sib <- setdiff(tr$kids[p - tr$nTip, ], v)
    for (y in .FnprTargets(tr, p, sib, .Descendants(tr, v), Inf)) {
      new <- tr
      new$parent[sib] <- tr$parent[[p]]
      new$parent[p] <- new$parent[[y]]
      new$parent[y] <- p
      new <- .Rebuild(new)
      if (is.null(trees[[Key(new)]])) {
        trees[[Key(new)]] <- todo[[Key(new)]] <- new
      }
    }
  }
}

# Clades of each topology, as sorted tip labels joined by "+"
Clades <- function(tr) {
  nTip <- tr$nTip
  below <- vector("list", length(tr$parent))
  below[seq_len(nTip)] <- as.list(tr$tipLabel)
  for (v in tr$postorder) {
    below[[v]] <- sort(unlist(below[tr$kids[v - nTip, ]]))
  }
  vapply(below[setdiff(tr$postorder, tr$root)], paste, "", collapse = "+")
}
cladeOf <- lapply(trees, Clades)
clades <- sort(unique(unlist(cladeOf)))
# Indicator of each statistic (rows) on each topology (columns)
hasClade <- vapply(cladeOf, function(x) clades %in% x, logical(length(clades)))
rownames(hasClade) <- paste0("clade ", clades)

exact <- if (gamma == 0) {
  logP <- vapply(trees, function(tr) TreeLikelihood(tr, par, model) + CharLL(tr, NULL, par, NULL), 0)
  exp(logP - .LogSumExp(logP))
} else {
  NULL
}

Chain <- function(move, seed) {
  set.seed(seed)
  set <- .AddSlices(start, .Precompute(start, par, model), 6)
  msg <- .Messages(start, set)
  d <- .SampleTree(start, set, msg)
  st <- list(tree = start, par = par, set = set, msg = msg, paths = d$paths,
             bin = d$bin)
  st$charLL <- CharLL(start, st$paths, par, set)
  visits <- character(nIter)
  accepted <- 0
  for (i in seq_len(nIter)) {
    # gamma = 0: an exact draw of every path; gamma = 1: local path moves
    for (pathMove in if (gamma == 0) "tree" else c("node", "tip", "branch")) {
      st <- .PathMove(pathMove, st, CharLL)$state
    }
    res <- .TreeMove(move, st, CharLL)
    st <- res$state
    accepted <- accepted + res$accepted
    visits[[i]] <- Key(st$tree)
  }
  keep <- visits[-seq_len(nIter %/% 10)]
  c(table(factor(keep, names(trees))) / length(keep), acceptance = accepted / nIter)
}
Chains <- function(move) {
  seeds <- offset + 10 * seq_len(nChain) + (move == "fnpr")
  runs <- if (nCore > 1) {
    parallel::mclapply(seeds, function(seed) Chain(move, seed), mc.cores = nCore)
  } else {
    lapply(seeds, function(seed) Chain(move, seed))
  }
  do.call(cbind, runs)
}
runs <- lapply(c(wfnpr = "wfnpr", fnpr = "fnpr"), Chains)

# Statistics: topologies of exact probability above 0.02 (with gamma = 1, of
# mean frequency above 0.02 over all chains), every clade, and acceptance
topologyP <- if (is.null(exact)) {
  rowMeans(cbind(runs$wfnpr, runs$fnpr)[seq_along(trees), ])
} else {
  exact
}
shown <- which(topologyP > 0.02)
Stats <- function(freq) {
  topology <- freq[seq_along(trees), , drop = FALSE]
  out <- rbind(topology[shown, , drop = FALSE], hasClade %*% topology,
               acceptance = freq["acceptance", ])
  rownames(out)[seq_along(shown)] <- sprintf("topology %03d", shown)
  out
}
perChain <- lapply(runs, Stats)
exactStat <- if (is.null(exact)) {
  rep(NA_real_, nrow(perChain$wfnpr))
} else {
  c(exact[shown], hasClade %*% exact, NA)
}
Mean <- function(move) rowMeans(perChain[[move]])
Se <- function(move) apply(perChain[[move]], 1, stats::sd) / sqrt(nChain)
out <- data.frame(
  label = label, statistic = rownames(perChain$wfnpr), exact = exactStat,
  wfnpr = Mean("wfnpr"), wfnprSe = Se("wfnpr"), fnpr = Mean("fnpr"), fnprSe = Se("fnpr"),
  zExact = (Mean("wfnpr") - exactStat) / Se("wfnpr"),
  zFnprExact = (Mean("fnpr") - exactStat) / Se("fnpr"),
  zFnpr = (Mean("wfnpr") - Mean("fnpr")) / sqrt(Se("wfnpr")^2 + Se("fnpr")^2),
  row.names = NULL
)
out[-(1:2)] <- signif(out[-(1:2)], 4)
long <- do.call(rbind, lapply(names(perChain), function(move) {
  data.frame(label = label, move = move, chain = rep(seq_len(nChain), each = nrow(out)),
             statistic = out$statistic, exact = exactStat,
             value = signif(as.vector(perChain[[move]]), 6))
}))
dir <- "data-raw/wfnpr-results"
write.csv(out, file.path(dir, paste0(label, ".csv")), row.names = FALSE)
write.csv(long, file.path(dir, paste0(label, "-chains.csv")), row.names = FALSE)
cat(sprintf("%s: %d topologies (%d of positive probability), %d chains x %d\n", label,
            length(trees), if (is.null(exact)) NA else sum(exact > 0), nChain, nIter))
print(out[-1])
