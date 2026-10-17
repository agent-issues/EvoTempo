# Component timings for one TempoMCMC() iteration on the rung-6 fixture:
# each building block of the moves, timed in isolation.
# Usage: Rscript dev/profiling/drivers/tempo-components.R   (from the package root)
lib <- Sys.getenv("EVOTEMPO_LIB", "")
if (nzchar(lib)) library(EvoTempo, lib.loc = lib) else library(EvoTempo)
source("dev/profiling/drivers/fixture-rung6.R")
ns <- asNamespace("EvoTempo")
for (f in c(".Precompute", ".AddSlices", ".Messages", ".SampleTree", ".LogAugmented",
            ".LogWeight", ".PathNodes", ".CharacterLengths", ".PathSummary",
            ".TreeStats", ".FrozenLikelihood", ".EGrid", ".Setup", ".TreeMove",
            ".PathMove", ".ProposeTreeMove")) {
  assign(f, get(f, ns))
}
options(EvoTempo.engine = "cpp")
set.seed(1)
par <- truth
s0 <- .Precompute(tt, par, model)
set <- .AddSlices(tt, s0, 16)
msg <- .Messages(tt, set)
draw <- .SampleTree(tt, set, msg)
nodes <- .PathNodes(tt)
x <- set$grid$x
lens <- .CharacterLengths(tt, x, draw$paths, par$rate, 1)
CharLL <- function(tr, paths, p, set) {
  MkLikelihood(tr, chars, .CharacterLengths(tr, set$grid$x, paths, p$rate, 1), 2L)
}
st <- list(tree = tt, par = par, set = set, msg = msg, paths = draw$paths,
           bin = draw$bin)
st$charLL <- CharLL(tt, st$paths, par, set)
# Median over 5 batches of the time per call, in microseconds (no bench dependency,
# so that the same script runs on Hamilton).  Windows' clock ticks at ~16 ms, so
# batches there run for at least 0.5 s
Time <- function(f, minTime = if (.Platform$OS.type == "windows") 0.5 else 0.05) {
  f()
  n <- 1L
  repeat {
    t <- system.time(for (i in seq_len(n)) f())[["elapsed"]]
    if (t >= minTime) break
    n <- n * 2L
  }
  median(vapply(1:5, function(b) {
    system.time(for (i in seq_len(n)) f())[["elapsed"]] / n
  }, 0)) * 1e6
}
calls <- alist(
  Setup = .Setup(par, model),
  EGrid = .EGrid(set$Q, set$lambda, set$mu, set$psi, model$rho,
                 tt$age[tt$root], tt$age),
  Precompute = .Precompute(tt, par, model),
  AddSlices = .AddSlices(tt, s0, 16),
  Messages = .Messages(tt, set),
  SampleTree = .SampleTree(tt, set, msg),
  LogAugmented = .LogAugmented(tt, set, draw$paths, draw$bin),
  LogWeight = .LogWeight(tt, set, msg, draw$paths, nodes),
  FrozenLik = .FrozenLikelihood(tt, set, msg, FALSE),
  CharLengths = .CharacterLengths(tt, x, draw$paths, par$rate, 1),
  MkLik = MkLikelihood(tt, chars, lens, 2L),
  CharLL = CharLL(tt, draw$paths, par, set),
  PathSummary = .PathSummary(tt, st, nodes),
  TreeStats = .TreeStats(tt, lens, lens),
  move_tree = .PathMove("tree", st, CharLL),
  move_node = .PathMove("node", st, CharLL),
  move_tip = .PathMove("tip", st, CharLL),
  move_branch = .PathMove("branch", st, CharLL),
  move_nodeAge = .TreeMove("nodeAge", st, CharLL),
  move_fossilAge = .TreeMove("fossilAge", st, CharLL),
  move_slide = .TreeMove("slide", st, CharLL),
  move_saJump = .TreeMove("saJump", st, CharLL)
)
us <- vapply(calls, function(cl) Time(function() eval(cl)), 0)
print(data.frame(expr = names(calls), median_us = round(us, 1)), row.names = FALSE)
