# Validation of the spliced fixed-height regraft `sfnpr` (#110 move 3;
# dev/plans/spliced-fnpr.md) at fixed node ages and parameters, where the
# topology has a known distribution even with gamma = 1.
#
# With one character, TreeLikelihood(character = ) is the exact joint
# likelihood of tree and character with tempo integrated out (product-space
# pruning, spec section 10), so P(topology) is proportional to it over the
# time trees that fixed-height regrafts reach from the start tree.  These are
# enumerated.  Independent chains are compared with the exact probabilities of
# (a) each time tree of probability above 0.02 and (b) each clade.
# `fnpr` chains, which redraw paths from the frozen model, run alongside as a
# control.
#
# Settings (environment):
#   SFNPR_TREE    "sa5": five tips, two fossils and a sampled ancestor;
#                 "extant5": five extant tips
#   SFNPR_CASE    "prior" (no character) or "character" (one character)
#   SFNPR_GAMMA   tempo exponent of character change (default 1)
#   SFNPR_CHAINS, SFNPR_ITER   chains per move and iterations per chain (the
#                 first tenth is burn-in)
#   SFNPR_SEED    seed offset;  SFNPR_CORES  parallel workers (not on Windows)
# Run from the package root with the package installed.  Results:
# data-raw/sfnpr-results/<tree>-<case>-g<gamma>-<seed>.csv.
suppressMessages(library("EvoTempo"))
Env <- function(name, default) Sys.getenv(name, default)
treeName <- Env("SFNPR_TREE", "sa5")
case <- Env("SFNPR_CASE", "character")
gamma <- as.numeric(Env("SFNPR_GAMMA", "1"))
nChain <- as.numeric(Env("SFNPR_CHAINS", "32"))
nIter <- as.numeric(Env("SFNPR_ITER", "20000"))
offset <- as.numeric(Env("SFNPR_SEED", "0"))
nCore <- as.integer(Env("SFNPR_CORES", "1"))
stopifnot(treeName %in% c("sa5", "extant5"), case %in% c("prior", "character"))
ns <- asNamespace("EvoTempo")
for (f in c(".SlideNodes", ".FnprTargets", ".Descendants", ".Rebuild")) {
  assign(f, get(f, ns))
}
label <- sprintf("%s-%s-g%s-%d", treeName, case, format(gamma), offset)

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
  character <- c(0, 0, 1, 1, 1)
} else {
  # (((a:1,b:1):1.5,c:2.5):1.5,(d:2,e:2):2);
  phylo <- Phylo(
    rbind(c(6, 7), c(7, 8), c(8, 1), c(8, 2), c(7, 3), c(6, 9), c(9, 4), c(9, 5)),
    c(1.5, 1.5, 1, 1, 2.5, 2, 2, 2), c("a", "b", "c", "d", "e")
  )
  character <- c(0, 1, 1, 0, 0)
}
start <- TempoTree(phylo)
par <- list(lambda = 0.5, mu = 0.2, psi = 0.1, theta = 0.3, s = 1, rate = 0.4)
model <- DCETModel(nBin = 7, zRange = c(-2, 2), gamma = gamma, rho = 0.7)

# A time tree: its clades (as tip sets) with their ages
Key <- function(tr) {
  paste(sort(vapply(tr$postorder, function(u) {
    tips <- .Descendants(tr, u)
    sprintf("%s@%.8f", paste(sort(tips[tips <= tr$nTip]), collapse = ","), tr$age[[u]])
  }, "")), collapse = ";")
}
Clades <- function(key) sub("@.*$", "", strsplit(key, ";")[[1]])

# Every time tree that fixed-height regrafts reach from `start`
Regraft <- function(tr, v, y) {
  p <- tr$parent[[v]]
  sib <- setdiff(tr$kids[p - tr$nTip, ], v)
  out <- tr
  out$parent[sib] <- tr$parent[[p]]
  out$parent[p] <- tr$parent[[y]]
  out$parent[y] <- p
  .Rebuild(out)
}
trees <- list()
queue <- list(start)
while (length(queue)) {
  tr <- queue[[1]]
  queue <- queue[-1]
  key <- Key(tr)
  if (!is.null(trees[[key]])) next
  trees[[key]] <- tr
  for (v in .SlideNodes(tr)) {
    p <- tr$parent[[v]]
    sib <- setdiff(tr$kids[p - tr$nTip, ], v)
    for (y in .FnprTargets(tr, p, sib, .Descendants(tr, v), Inf)) {
      queue[[length(queue) + 1L]] <- Regraft(tr, v, y)
    }
  }
}
logL <- vapply(trees, function(tr) {
  if (case == "prior") TreeLikelihood(tr, par, model)
  else TreeLikelihood(tr, par, model, character = character)
}, 0)
exact <- exp(logL - max(logL[is.finite(logL)]))
exact[!is.finite(logL)] <- 0
exact <- exact / sum(exact)
cladeSet <- sort(unique(unlist(lapply(names(trees), Clades))))
exactClade <- vapply(cladeSet, function(cl) {
  sum(exact[vapply(names(trees), function(k) cl %in% Clades(k), TRUE)])
}, 0)
cat(sprintf("%s: %d time trees reachable, %d of probability > 0\n", label,
            length(trees), sum(exact > 0)))

# Chains: the time tree of every sample, recorded through .TreeStats()
Chain <- function(move, seed) {
  keys <- character(0)
  statsOriginal <- get(".TreeStats", ns)
  Recording <- function(tree, lengths = NULL, inverse = NULL) {
    if (!is.null(lengths)) keys[length(keys) + 1L] <<- Key(tree)
    statsOriginal(tree, lengths, inverse)
  }
  utils::assignInNamespace(".TreeStats", Recording, ns = "EvoTempo")
  on.exit(utils::assignInNamespace(".TreeStats", statsOriginal, ns = "EvoTempo"))
  set.seed(seed)
  TempoMCMC(start, model, par, list(), nIter = nIter, burnin = nIter / 10,
            nSlice = 4, moves = c(node = 1, tip = 1, branch = 1, setNames(1, move)),
            characters = if (case == "character") matrix(character),
            minAcceptance = 0)
  keys
}
jobs <- expand.grid(chain = seq_len(nChain), move = c("sfnpr", "fnpr"),
                    stringsAsFactors = FALSE)
Run <- function(i) Chain(jobs$move[[i]], 1000 * offset + i)
draws <- if (nCore > 1) {
  parallel::mclapply(seq_len(nrow(jobs)), Run, mc.cores = nCore)
} else {
  lapply(seq_len(nrow(jobs)), Run)
}

# Per-chain frequencies of each statistic against the exact value
# (paste0() of a prefix and character(0) would give the prefix alone)
statNames <- c(sprintf("tree:%s", names(exact)[exact > 0.02]), sprintf("clade:%s", cladeSet))
exactStat <- c(exact[exact > 0.02], exactClade)
Freq <- function(keys) {
  c(vapply(names(exact)[exact > 0.02], function(k) mean(keys == k), 0),
    vapply(cladeSet, function(cl) mean(vapply(keys, function(k) cl %in% Clades(k), TRUE)), 0))
}
out <- do.call(rbind, lapply(c("sfnpr", "fnpr"), function(move) {
  f <- vapply(draws[jobs$move == move], Freq, numeric(length(statNames)))
  m <- rowMeans(f)
  se <- apply(f, 1, stats::sd) / sqrt(ncol(f))
  data.frame(move = move, statistic = statNames, exact = exactStat, mean = m, se = se,
             z = ifelse(se > 0, (m - exactStat) / se, NA))
}))
dir.create("data-raw/sfnpr-results", showWarnings = FALSE)
utils::write.csv(out, sprintf("data-raw/sfnpr-results/%s.csv", label), row.names = FALSE)
for (move in c("sfnpr", "fnpr")) {
  z <- out$z[out$move == move]
  cat(sprintf("%-6s statistics %d, max |z| %.2f, rms z %.2f\n", move, sum(!is.na(z)),
              max(abs(z), na.rm = TRUE), sqrt(mean(z^2, na.rm = TRUE))))
}
