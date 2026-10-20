# Mixing study (issue #19): diagnosis, candidate moves and invariance oracles.
# Run from the package root (EvoTempo installed, or EVOTEMPO_LOAD_ALL=1):
#   Rscript dev/mixing/run.R sim
#       simulate the fixed true trees (cell beta = gamma = 1) -> trees.rds
#   Rscript dev/mixing/run.R chain <tree> <nChar> <config> <rep> [nIter] [burnin]
#       one chain on tree <tree> (1 small, 2 large) with the first nChar
#       characters, under move set <config> (see `configs`)
#   Rscript dev/mixing/run.R probe <tree> <nChar>
#       fixed-parameter states at the truth: spread of fresh path proposals
#       (whole-tree FFBS, anchored bridges) and conditional parameter spread
#   Rscript dev/mixing/run.R oracle <variant> [nIter] [nChain]
#       invariance checks of the new moves on the rung 4b tree (one character)
#   Rscript dev/mixing/run.R check
#       Driver() reproduces TempoMCMC(); anchored density identities
# Environment: EVOTEMPO_MIXING_OUT (default dev/mixing/results)
source("dev/mixing/common.R")
args <- commandArgs(trailingOnly = TRUE)
outDir <- Sys.getenv("EVOTEMPO_MIXING_OUT", "dev/mixing/results")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
Save <- function(x, file) {
  tmp <- paste0(file, ".tmp")
  saveRDS(x, tmp)
  file.rename(tmp, file)
}
treeFile <- file.path(outDir, "trees.rds")
allPar <- c("lambda", "mu", "psi", "theta", "s", "rate")
compAll <- c("rate", "lambda", "mu", "theta")

configs <- list(
  base = list(moves = baseMoves),
  fixed = list(moves = baseMoves, fixed = TRUE),
  shift = list(moves = baseMoves, extra = list(ShiftExtra(compAll))),
  shiftc = list(moves = baseMoves, extra = list(ShiftExtra("rate"))),
  ajInf = list(moves = baseMoves, extra = list(AnchoredExtra(c("theta", "s"), Inf))),
  aj4 = list(moves = baseMoves, extra = list(AnchoredExtra(c("theta", "s"), 4))),
  aj2 = list(moves = baseMoves, extra = list(AnchoredExtra(c("theta", "s"), 2))),
  aj1 = list(moves = baseMoves, extra = list(AnchoredExtra(c("theta", "s"), 1))),
  combo = list(moves = baseMoves,
               extra = list(ShiftExtra(compAll), AnchoredExtra(c("theta", "s"), 2))),
  comboNJ = list(moves = replace(baseMoves, "joint", 0),
                 extra = list(ShiftExtra(compAll),
                              AnchoredExtra(c("lambda", "mu", "psi", "theta", "s"), 2))),
  # Wave 2: no whole-tree FFBS or joint moves; local path moves scaled to the
  # tree (fraction `frac` of the internal nodes, tips and branches per iteration)
  L25 = list(local = 0.25),
  L100 = list(local = 1),
  fixedL25 = list(local = 0.25, fixed = TRUE),
  fixedL100 = list(local = 1, fixed = TRUE),
  L25S = list(local = 0.25, extra = list(ShiftExtra(compAll))),
  L25SA = list(local = 0.25, extra = list(ShiftExtra(compAll),
                                          AnchoredExtra(c("theta", "s"), 1))),
  L25SAc = list(local = 0.25, extra = list(ShiftExtra(compAll),
                                           AnchoredExtra("theta", 1),
                                           AnchoredExtra("s", 1, comp = compAll))),
  L25SN = list(local = 0.25, extra = list(ShiftExtra(compAll)),
               ncmc = list(which = c("theta", "s"), nStage = 8L, frac = 0.25)),
  L25SN16 = list(local = 0.25, extra = list(ShiftExtra(compAll)),
                 ncmc = list(which = c("theta", "s"), nStage = 16L, frac = 0.25)),
  L25SAJ = list(local = 0.25, keepJoint = 0.1,
                extra = list(ShiftExtra(compAll), AnchoredExtra(c("theta", "s"), 1)))
)

# Local-move counts per iteration for configuration field `local`
LocalMoves <- function(tt, frac) {
  tree <- TempoTree(tt)
  nodes <- .PathNodes(tree)
  c(node = ceiling(frac * length(tree$postorder)),
    tip = ceiling(frac * sum(nodes <= tree$nTip)),
    branch = ceiling(frac * length(nodes)))
}

if (!length(args)) stop("see header for usage")
mode <- args[[1]]

if (mode == "sim") {
  model <- Model(1, 1)
  set.seed(19019L)
  trees <- list(Simulate(model, Meanlog("high"), sizeLevels$small),
                Simulate(model, Meanlog("high"), sizeLevels$large))
  for (tr in trees) {
    cat(sprintf("tips: %d extant, %d fossil, %d SA; truth %s\n", tr$nExtant, tr$nFossil,
                tr$nAncestor, paste(names(tr$truth), signif(unlist(tr$truth), 3),
                                    collapse = " ")))
  }
  Save(trees, treeFile)
}

if (mode == "check") {
  model <- Model(1, 1)
  trees <- readRDS(treeFile)
  tr <- trees[[1]]
  priors <- TruncPriors(Meanlog("high"), sdlog)
  cat("Driver == TempoMCMC:", CheckDriver(tr$tt, model, tr$truth, priors,
                                          tr$chars[, 1:30], nIter = 20, burnin = 20), "\n")
  tree <- TempoTree(tr$tt)
  set <- .AddSlices(tree, .Precompute(tree, tr$truth, model), nSlice)
  msg <- .Messages(tree, set)
  d <- .SampleTree(tree, set, msg)
  nodes <- .PathNodes(tree)
  cat("anchored(Inf) - bridges:", LogAnchored(tree, set, msg, d$paths, Inf, 0) -
        .LogBridges(tree, set, nodes, d$paths), "\n")
  for (ev in c(1, 2, 4, Inf)) {
    v <- replicate(4, {
      dd <- DrawAnchored(tree, set, msg, d$paths, ev, 0)
      .LogAugmented(tree, set, dd, d$bin) - LogAnchored(tree, set, msg, dd, ev, 0) -
        .LogWeight(tree, set, msg, dd, nodes)
    })
    cat("every", ev, "constant spread", diff(range(v)), "\n")
  }
}

RunChain <- function(tr, nChar, cfgName, rep, nIter, burnin, thin) {
  cfg <- configs[[cfgName]]
  if (!is.null(cfg$local)) {
    cfg$moves <- c(LocalMoves(tr$tt, cfg$local), joint = cfg$keepJoint %||% 0)
  }
  if (!is.null(cfg$ncmc)) {
    nRelax <- ceiling(cfg$ncmc$frac * length(.PathNodes(TempoTree(tr$tt))))
    cfg$extra <- c(cfg$extra, list(NcmcExtra(cfg$ncmc$which, cfg$ncmc$nStage, nRelax)))
  }
  model <- Model(1, 1)
  priors <- TruncPriors(Meanlog("high"), sdlog)
  chars <- if (nChar > 0) tr$chars[, seq_len(nChar), drop = FALSE] else NULL
  t0 <- proc.time()[["elapsed"]]
  run <- Driver(tr$tt, model, tr$truth, priors, nIter = nIter, burnin = burnin,
                thin = thin, nSlice = nSlice, moves = cfg$moves, characters = chars,
                updateParameters = !isTRUE(cfg$fixed), extra = cfg$extra %||% list())
  wall <- proc.time()[["elapsed"]] - t0
  run$state <- NULL
  c(run, list(summary = Summarise(run), wall = wall, nIter = nIter, burnin = burnin,
              thin = thin, config = cfgName, nChar = nChar, rep = rep))
}

if (mode == "chain") {
  treeId <- as.integer(args[[2]])
  nChar <- as.integer(args[[3]])
  cfgName <- args[[4]]
  rep <- as.integer(args[[5]])
  nIter <- if (length(args) > 5) as.integer(args[[6]]) else 5000L
  burnin <- if (length(args) > 6) as.integer(args[[7]]) else 2000L
  file <- file.path(outDir, sprintf("chain-t%d-c%03d-%s-r%d.rds", treeId, nChar,
                                    cfgName, rep))
  if (!file.exists(file)) {
    trees <- readRDS(treeFile)
    set.seed(1900000L + 1000L * treeId + nChar + 100000L * rep +
               7919L * match(cfgName, names(configs)))
    out <- RunChain(trees[[treeId]], nChar, cfgName, rep, nIter, burnin,
                    thin = max(1L, nIter %/% 1000L))
    Save(out, file)
    s <- out$summary
    cat(sprintf("%s: %.0f s; ESS %s\n", basename(file), out$wall,
                paste(names(s$ess), round(s$ess), collapse = " ")))
    print(round(s$acceptance, 3))
  }
}

if (mode == "probe") {
  treeId <- as.integer(args[[2]])
  nChar <- as.integer(args[[3]])
  file <- file.path(outDir, sprintf("probe-t%d-c%03d.rds", treeId, nChar))
  trees <- readRDS(treeFile)
  tr <- trees[[treeId]]
  model <- Model(1, 1)
  priors <- TruncPriors(Meanlog("high"), sdlog)
  chars <- if (nChar > 0) tr$chars[, seq_len(nChar), drop = FALSE] else NULL
  CharLL <- CharLLFun(chars, model)
  set.seed(1910000L + 1000L * treeId + nChar)
  run <- Driver(tr$tt, model, tr$truth, priors, nIter = 1, burnin = 1500,
                nSlice = nSlice, characters = chars, updateParameters = FALSE)
  st <- run$state
  tree <- st$tree
  nodes <- .PathNodes(tree)
  nState <- 10L
  nDraw <- 40L
  everyLevels <- c(Inf, 4, 2, 1)
  probe <- vector("list", nState)
  for (i in seq_len(nState)) {
    run <- Driver(tr$tt, model, tr$truth, priors, nIter = 1, burnin = 100,
                  nSlice = nSlice, characters = chars, updateParameters = FALSE,
                  initial = st)
    st <- run$state
    cur <- c(lw = .LogWeight(tree, st$set, st$msg, st$paths, nodes), char = st$charLL)
    Spread <- function(Draw) {
      t(replicate(nDraw, {
        p <- Draw()
        c(lw = .LogWeight(tree, st$set, st$msg, p, nodes),
          char = CharLL(tree, p, st$par, st$set),
          clock = sum(.CharacterLengths(tree, st$set$grid$x, p, 1, 1)),
          nJumps = sum(vapply(nodes, function(v) length(p[[v]]$t), 0)))
      }))
    }
    fresh <- list(tree = Spread(function() .SampleTree(tree, st$set, st$msg)$paths))
    for (ev in everyLevels) {
      fresh[[paste0("anchored", ev)]] <- Spread(function() {
        DrawAnchored(tree, st$set, st$msg, st$paths, ev,
                     if (is.finite(ev)) sample.int(ev, 1L) - 1L else 0L)
      })
    }
    # Conditional spread of the parameters given the paths: parameter sweep only
    cond <- Driver(tr$tt, model, tr$truth, priors, nIter = 300, burnin = 100,
                   nSlice = nSlice, moves = c(tree = 0), characters = chars,
                   initial = st)$samples
    probe[[i]] <- list(current = cur, fresh = fresh,
                       condSdLog = vapply(cond[allPar], function(x) stats::sd(log(x)), 0),
                       condMeanLog = vapply(cond[allPar], function(x) mean(log(x)), 0),
                       nJumps = NJumps(st), clock = sum(BranchClock(st)))
  }
  Save(list(probe = probe, treeId = treeId, nChar = nChar, truth = tr$truth), file)
  cat("probe done\n")
}

if (mode == "oracle") {
  variant <- args[[2]]
  nIter <- if (length(args) > 2) as.integer(args[[3]]) else 12000L
  nChain <- if (length(args) > 3) as.integer(args[[4]]) else 8L
  block <- if (length(args) > 4) as.integer(args[[5]]) else 0L
  file <- file.path(outDir, if (block) {
    sprintf("oracle-%s-%d-b%02d.rds", variant, nIter, block)
  } else {
    sprintf("oracle-%s-%d.rds", variant, nIter)
  })
  r4b <- readRDS("data-raw/rung5-results/rung4b.rds")
  tt <- r4b$tt
  model <- r4b$model
  truth <- r4b$truth
  chars <- matrix(r4b$character)
  priors <- c(TempoPriors(meanlog = c(lambda = log(0.5), mu = log(0.2), psi = log(0.1),
                                      theta = log(0.1), s = log(0.8)),
                          sdlog = c(lambda = 0.5, mu = 0.5, psi = 0.5, theta = 0.75,
                                    s = 0.3)),
              rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  pathMoves <- c(tree = 1, node = 2, tip = 2, branch = 2)
  spec <- switch(
    variant,
    thetaAjInf = list(free = "theta", sweep = FALSE,
                      extra = list(AnchoredExtra("theta", Inf))),
    thetaAj2 = list(free = "theta", sweep = FALSE, extra = list(AnchoredExtra("theta", 2))),
    sAj2 = list(free = "s", sweep = FALSE, extra = list(AnchoredExtra("s", 2))),
    rateShift = list(free = "rate", sweep = TRUE, extra = list(ShiftExtra("rate", n = 3))),
    multi = list(free = c("lambda", "mu", "theta", "s", "rate"), sweep = TRUE,
                 extra = list(ShiftExtra(compAll, n = 2), AnchoredExtra(c("theta", "s"), 2))),
    thetaAj1 = list(free = "theta", sweep = FALSE, extra = list(AnchoredExtra("theta", 1))),
    sAj1 = list(free = "s", sweep = FALSE, extra = list(AnchoredExtra("s", 1))),
    thetaSweep = list(free = "theta", sweep = TRUE, extra = list()),
    thetaShift = list(free = "theta", sweep = TRUE, extra = list(ShiftExtra("theta", n = 3))),
    lambdaShift = list(free = "lambda", sweep = TRUE, extra = list(ShiftExtra("lambda", n = 3))),
    muShift = list(free = "mu", sweep = TRUE, extra = list(ShiftExtra("mu", n = 3))),
    thetaNcmc = list(free = "theta", sweep = FALSE, extra = list(NcmcExtra("theta", 4L, 10L))),
    sNcmc = list(free = "s", sweep = FALSE, extra = list(NcmcExtra("s", 4L, 10L))),
    sComp = list(free = c("lambda", "mu", "theta", "s", "rate"), sweep = TRUE,
                 extra = list(ShiftExtra(compAll, n = 2), AnchoredExtra("theta", 1),
                              AnchoredExtra("s", 1, comp = compAll, n = 2)))
  )
  meanlog <- c(lambda = log(0.5), mu = log(0.2), psi = log(0.1), theta = log(0.1),
               s = log(0.8), rate = log(0.3))
  sdl <- c(lambda = 0.5, mu = 0.5, psi = 0.5, theta = 0.75, s = 0.3, rate = 0.5)
  GridFor <- function(p) exp(seq(meanlog[[p]] - 3.5 * sdl[[p]], meanlog[[p]] + 3.5 * sdl[[p]],
                                 length.out = 161))
  set.seed(1920000L + sum(utf8ToInt(variant)) + 7777L * block)
  burnin <- 500L
  chains <- lapply(seq_len(nChain), function(i) {
    Driver(tt, model, truth, priors[spec$free], nIter = nIter, burnin = burnin,
           nSlice = 20L, moves = pathMoves, characters = chars, extra = spec$extra,
           sweep = spec$sweep)[c("samples", "acceptance", "seconds")]
  })
  exact <- NULL
  marginal <- NULL
  if (length(spec$free) == 1L) {
    p <- spec$free
    g <- GridFor(p)
    ll <- vapply(g, function(v) {
      par <- truth
      par[[p]] <- v
      TreeLikelihood(tt, par, model, character = r4b$character)
    }, 0)
    exact <- data.frame(value = g, logLik = ll, logPrior = priors[[p]](g))
  } else {
    marginal <- lapply(seq_len(nChain), function(i) {
      MarginalMCMC(tt, model, truth, priors[spec$free], nIter = nIter, burnin = burnin,
                   character = r4b$character)
    })
  }
  Save(list(variant = variant, spec = spec[c("free", "sweep")], chains = chains,
            exact = exact, marginal = marginal, nIter = nIter), file)
  cat("oracle", variant, "done\n")
}
