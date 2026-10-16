# C++ engine vs the R reference on the rung 3 and rung 4 fixtures, at
# K = 16 (z in [-3, 3]) and K = 26 (z in [-5, 5]); and run times.
#   Agreement: relative difference |cpp - R| / |R| of the log-likelihood
#   (TreeLikelihood(), with and without one character), of the augmented
#   density, exact/frozen weights and frozen likelihood of FFBS draws, and of
#   the character likelihood; FFBS draws under one seed must give identical
#   bins.  A tight-tolerance R run (rtol 1e-13) shows which engine is nearer
#   the exact value.
#   Benchmark: median seconds per call in each engine.
# Run from the package root, with EvoTempo installed:
#   Rscript data-raw/cpp-engine-agreement.R [nCore]
library("EvoTempo")
library("parallel")
args <- commandArgs(trailingOnly = TRUE)
nCore <- if (length(args)) as.integer(args[[1]]) else 4L
outDir <- "data-raw/cpp-engine-results"
dir.create(outDir, showWarnings = FALSE)
ns <- asNamespace("EvoTempo")
Internal <- function(name) get(name, envir = ns)
.Precompute <- Internal(".Precompute")
.AddSlices <- Internal(".AddSlices")
.Messages <- Internal(".Messages")
.SampleTree <- Internal(".SampleTree")
.LogAugmented <- Internal(".LogAugmented")
.LogWeight <- Internal(".LogWeight")
.FrozenLikelihood <- Internal(".FrozenLikelihood")
.PathNodes <- Internal(".PathNodes")
.CharacterLengths <- Internal(".CharacterLengths")

WithEngine <- function(engine, expr) {
  old <- options(EvoTempo.engine = engine)
  on.exit(options(old))
  expr
}
RelDiff <- function(a, b) ifelse(a == b, 0, abs(a - b) / abs(b))

grids <- list(K16 = list(nBin = 16, zRange = c(-3, 3)),
              K26 = list(nBin = 26, zRange = c(-5, 5)))

# Fixtures ---------------------------------------------------------------------
# Rung 3 (data-raw/rung3-validation.R): sampled-ancestor tree with a stem,
# 12-tip ultrametric tree, 10-tip tree with fossil tips; 20 parameter points
tr <- ape::read.tree(text = "((a:2,(b:1.5,f:0):1.5):1,(c:1.2,d:3.3):0.7);")
tr$root.edge <- 0.8
set.seed(3)
ultra <- ape::rcoal(12)
ultra$edge.length <- ultra$edge.length * 3
set.seed(11)
nonUltra <- ape::rtree(10)
nonUltra$edge.length <- nonUltra$edge.length * 2
rung3Trees <- list(sampledAncestor = TempoTree(tr), ultrametric = TempoTree(ultra),
                   fossilTips = TempoTree(nonUltra))
set.seed(5)
points <- data.frame(lambda = runif(20, 0.3, 1), mu = runif(20, 0, 0.3),
                     theta = runif(20, 0.05, 1), s = runif(20, 0.5, 1.5))
points$psi <- 0.15
rung3Settings <- list(
  sampledAncestor = list(condition = "survival", start = "origin", joint = FALSE,
                         removal = 0.4, rho = 0.6),
  ultrametric = list(condition = "extant", start = "root", joint = TRUE,
                     removal = 0, rho = 0.7),
  fossilTips = list(condition = "none", start = "root", joint = FALSE,
                    removal = 1, rho = 0.7)
)
# Rung 4 (data-raw/rung4-validation.R): data tree, truth, and 20 draws from
# the priors; rung 4b character on the same tree
r4 <- readRDS("data-raw/rung4-results/rung4.rds")
r4b <- readRDS("data-raw/rung5-results/rung4b.rds")
priors <- TempoPriors()
set.seed(8)
rung4Points <- rbind(
  as.data.frame(r4$truth),
  data.frame(lambda = rlnorm(19, log(0.5), 0.5), mu = rlnorm(19, log(0.2), 0.5),
             psi = rlnorm(19, log(0.1), 0.5), theta = rlnorm(19, log(0.1), 0.75),
             s = rlnorm(19, log(0.8), 0.3))
)

jobs <- rbind(
  expand.grid(fixture = names(rung3Trees), point = seq_len(20), grid = names(grids),
              stringsAsFactors = FALSE),
  expand.grid(fixture = c("rung4", "rung4b"), point = seq_len(20), grid = names(grids),
              stringsAsFactors = FALSE)
)

Job <- function(j) {
  job <- jobs[j, ]
  g <- grids[[job$grid]]
  if (job$fixture %in% names(rung3Trees)) {
    tree <- rung3Trees[[job$fixture]]
    st <- rung3Settings[[job$fixture]]
    model <- DCETModel(nBin = g$nBin, zRange = g$zRange, rho = st$rho,
                       removal = st$removal, condition = st$condition,
                       joint = st$joint, start = st$start)
    par <- as.list(points[job$point, ])
    character <- NULL
  } else {
    tree <- r4$tt
    model <- DCETModel(nBin = g$nBin, zRange = g$zRange, rho = 0.5)
    par <- as.list(rung4Points[job$point, ])
    par$rate <- 0.3
    character <- if (job$fixture == "rung4b") r4b$character else NULL
  }
  Time <- function(f) {
    t0 <- proc.time()[["elapsed"]]
    value <- f()
    list(value = value, time = proc.time()[["elapsed"]] - t0)
  }
  cpp <- Time(function() TreeLikelihood(tree, par, model, character, engine = "cpp"))
  ref <- Time(function() TreeLikelihood(tree, par, model, character, engine = "R"))
  out <- data.frame(job[rep(1, 1), ], quantity = "logLik", cpp = cpp$value,
                    R = ref$value, cppTime = cpp$time, RTime = ref$time,
                    tight = NA_real_)
  if (job$point <= 3 && is.null(character)) {
    out$tight <- TreeLikelihood(tree, par, model, rtol = 1e-13, atol = 1e-18,
                                engine = "R")
  }
  if (job$fixture == "rung4" && job$point <= 5) {
    # Augmented density of FFBS draws, identical seeds in both engines
    for (engine in c("cpp", "R")) {
      t0 <- proc.time()[["elapsed"]]
      res <- WithEngine(engine, {
        set <- .AddSlices(tree, .Precompute(tree, par, model), 20)
        msg <- .Messages(tree, set)
        tPre <- proc.time()[["elapsed"]] - t0
        set.seed(100 + job$point)
        t1 <- proc.time()[["elapsed"]]
        draws <- lapply(1:10, function(i) .SampleTree(tree, set, msg))
        tDraw <- (proc.time()[["elapsed"]] - t1) / 10
        t2 <- proc.time()[["elapsed"]]
        aug <- vapply(draws, function(d) .LogAugmented(tree, set, d$paths, d$bin), 0)
        tAug <- (proc.time()[["elapsed"]] - t2) / 10
        lw <- vapply(draws, function(d) {
          .LogWeight(tree, set, msg, d$paths, .PathNodes(tree))
        }, 0)
        t3 <- proc.time()[["elapsed"]]
        mk <- vapply(draws, function(d) {
          MkLikelihood(tree, matrix(r4b$character),
                       .CharacterLengths(tree, set$grid$x, d$paths, 0.3, 1))
        }, 0)
        tMk <- (proc.time()[["elapsed"]] - t3) / 10
        list(frozen = .FrozenLikelihood(tree, set, msg), aug = aug, lw = lw, mk = mk,
             bins = lapply(draws, `[[`, "bin"), tPre = tPre, tDraw = tDraw,
             tAug = tAug, tMk = tMk)
      })
      assign(paste0("res", engine), res)
    }
    same <- mean(mapply(identical, rescpp$bins, resR$bins))
    extra <- rbind(
      data.frame(quantity = "frozenLogLik", cpp = rescpp$frozen, R = resR$frozen,
                 cppTime = rescpp$tPre, RTime = resR$tPre),
      data.frame(quantity = "logAugmented", cpp = rescpp$aug, R = resR$aug,
                 cppTime = rescpp$tAug, RTime = resR$tAug),
      data.frame(quantity = "logWeight", cpp = rescpp$lw, R = resR$lw,
                 cppTime = NA, RTime = NA),
      data.frame(quantity = "ffbsDraw", cpp = same, R = 1,
                 cppTime = rescpp$tDraw, RTime = resR$tDraw),
      data.frame(quantity = "charLogLik", cpp = rescpp$mk, R = resR$mk,
                 cppTime = rescpp$tMk, RTime = resR$tMk)
    )
    out <- rbind(out, data.frame(job[rep(1, nrow(extra)), ], extra, tight = NA_real_))
  }
  out
}

t0 <- Sys.time()
res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), Job, mc.cores = nCore,
                               mc.preschedule = FALSE))
cat("Agreement runtime", format(difftime(Sys.time(), t0, units = "mins")), "\n")
res$relDiff <- RelDiff(res$cpp, res$R)
res$cppVsTight <- RelDiff(res$cpp, res$tight)
res$RVsTight <- RelDiff(res$R, res$tight)
saveRDS(res, file.path(outDir, "agreement.rds"))

agreement <- do.call(rbind, lapply(split(res, list(res$fixture, res$grid, res$quantity),
                                         drop = TRUE), function(d) {
  data.frame(fixture = d$fixture[[1]], grid = d$grid[[1]], quantity = d$quantity[[1]],
             n = nrow(d), maxRelDiff = max(d$relDiff),
             cppVsTight = suppressWarnings(max(d$cppVsTight, na.rm = TRUE)),
             RVsTight = suppressWarnings(max(d$RVsTight, na.rm = TRUE)))
}))
agreement$cppVsTight[!is.finite(agreement$cppVsTight)] <- NA
agreement$RVsTight[!is.finite(agreement$RVsTight)] <- NA
write.csv(agreement, file.path(outDir, "agreement.csv"), row.names = FALSE)

timed <- res[!is.na(res$cppTime) & !is.na(res$RTime), ]
benchmark <- do.call(rbind, lapply(split(timed, list(timed$fixture, timed$grid, timed$quantity),
                                         drop = TRUE), function(d) {
  data.frame(fixture = d$fixture[[1]], grid = d$grid[[1]], quantity = d$quantity[[1]],
             n = nrow(d), RSeconds = median(d$RTime), cppSeconds = median(d$cppTime),
             speedup = median(d$RTime) / median(d$cppTime))
}))
write.csv(benchmark, file.path(outDir, "benchmark.csv"), row.names = FALSE)

# One TempoMCMC run per engine, rung 4 tree with the rung 4b character
chainTime <- do.call(rbind, lapply(names(grids), function(gn) {
  g <- grids[[gn]]
  model <- DCETModel(nBin = g$nBin, zRange = g$zRange, rho = 0.5)
  start <- c(r4$truth, rate = 0.3)
  pri <- c(priors, rate = function(x) stats::dlnorm(x, log(0.3), 0.5, log = TRUE))
  do.call(rbind, lapply(c("cpp", "R"), function(engine) {
    set.seed(1)
    t0 <- proc.time()[["elapsed"]]
    run <- TempoMCMC(r4$tt, model, start, pri, nIter = 20, nSlice = 20,
                     characters = matrix(r4b$character),
                     moves = c(tree = 1, node = 1, tip = 1, branch = 1, joint = 1),
                     engine = engine)
    data.frame(grid = gn, engine = engine,
               secondsPerIteration = (proc.time()[["elapsed"]] - t0) / 20,
               finalLogPost = run$samples$logPost[20])
  }))
}))
write.csv(chainTime, file.path(outDir, "chain-benchmark.csv"), row.names = FALSE)

print(agreement, digits = 3)
print(benchmark, digits = 3)
print(chainTime, digits = 6)
cat("Max relative difference in logL:",
    format(max(res$relDiff[res$quantity == "logLik"])), "\n")
