# Seed and mutation sweep of the exact-target oracles of
# tests/testthat/helper-oracles.R (#51): psi compensation by `shift`
# (delta != 0), `shift` with `s` free, and the gamma = 0 default schedule.
# A mutation replaces .Shifted() in the loaded namespace for this process only;
# R/MCMC.R is not edited.
#   none           the package as it is;
#   psiJacobian    psi compensated, but its term dropped from the Jacobian;
#   psiSign        psi multiplied by exp(+delta d), Jacobian unchanged;
#   psiNone        psi not compensated (a symmetric bins-only shift: still exact);
#   rateJacobian   rate compensated, but its term dropped from the Jacobian;
#   staleS         d computed at s = 1 instead of the current s;
#   rateGamma0     at gamma = 0, rate compensated as if gamma were 1, while the
#                  stored character likelihood is kept.
# Usage, from the package root with EvoTempo installed:
#   Rscript data-raw/oracles-51.R <oracle> <mutation,...> <seed,...> [nIter] [nChain] [burnin]
# with <oracle> psi, s, sNoShift (OracleS() without `shift`), gamma0 or identity.
# `identity` is ShiftIdentityError(): its "z" is the largest error in the kept
# character likelihood.  Appends one row per run to
# data-raw/mixing-results/oracles-51-<oracle>.csv.
suppressMessages(library("EvoTempo"))
ns <- asNamespace("EvoTempo")
for (fn in ls(ns, all.names = TRUE)) {
  if (!exists(fn, globalenv(), inherits = FALSE)) assign(fn, get(fn, ns))
}
source("tests/testthat/helper-oracles.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("see header for usage")
oracle <- args[[1]]
nIter <- if (length(args) > 3) as.integer(args[[4]]) else 3000L
nChain <- if (length(args) > 4) as.integer(args[[5]]) else 2L
burnin <- if (length(args) > 5) as.integer(args[[6]]) else 200L
outFile <- sprintf("data-raw/mixing-results/oracles-51-%s.csv", oracle)
Original <- get(".Shifted", ns)

Mutant <- function(kind) {
  if (kind == "none") {
    return(Original)
  }
  function(st, m, model, free) {
    d <- m * st[["set"]][["grid"]][["h"]] * st[["par"]][["s"]]
    if (kind == "staleS") {
      st1 <- st
      st1[["par"]][["s"]] <- 1
      out <- Original(st1, m, model, free)
      if (!is.null(out)) out[["par"]][["s"]] <- st[["par"]][["s"]]
      return(out)
    }
    if (kind == "rateJacobian") {
      out <- Original(st, m, model, free)
      if (!is.null(out) && "rate" %in% out[["comp"]]) {
        out[["logJac"]] <- out[["logJac"]] + model[["gamma"]] * d
      }
      return(out)
    }
    if (kind == "rateGamma0") {
      out <- Original(st, m, model, free)
      if (!is.null(out) && "rate" %in% free) {
        out[["par"]][["rate"]] <- st[["par"]][["rate"]] * exp(-d)
        out[["comp"]] <- union(out[["comp"]], "rate")
        out[["logJac"]] <- out[["logJac"]] - d
      }
      return(out)
    }
    out <- Original(st, m, model, free)
    if (is.null(out) || !"psi" %in% out[["comp"]]) {
      return(out)
    }
    delta <- model[["delta"]]
    switch(kind,
           psiJacobian = out[["logJac"]] <- out[["logJac"]] + delta * d,
           psiSign = out[["par"]][["psi"]] <- st[["par"]][["psi"]] * exp(delta * d),
           psiNone = {
             out[["par"]][["psi"]] <- st[["par"]][["psi"]]
             out[["comp"]] <- setdiff(out[["comp"]], "psi")
             out[["logJac"]] <- out[["logJac"]] + delta * d
           },
           stop("Unknown mutation ", kind))
    out
  }
}

Run <- switch(oracle, psi = OraclePsi, s = OracleS, gamma0 = OracleGamma0,
              sNoShift = function(...) OracleS(..., moves = c(node = 1, tip = 1, branch = 1)),
              identity = function(seed, nIter, nChain, burnin) {
                res <- ShiftIdentityError(seed)
                list(z = c(error = res[["error"]]), shift = res[["accepted"]] / 40)
              },
              stop("Unknown oracle ", oracle))
for (kind in strsplit(args[[2]], ",")[[1]]) {
  utils::assignInNamespace(".Shifted", Mutant(kind), "EvoTempo")
  for (seed in as.integer(strsplit(args[[3]], ",")[[1]])) {
    t0 <- Sys.time()
    res <- Run(seed, nIter = nIter, nChain = nChain, burnin = burnin)
    sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    row <- data.frame(oracle = oracle, mutation = kind, seed = seed, nIter = nIter,
                      nChain = nChain, burnin = burnin, maxAbsZ = max(abs(res$z)),
                      shiftAcceptance = res$shift, seconds = sec,
                      t(round(res$z, 3)))
    cat(sprintf("%s %s seed %d: max|z| %.2f (%s), shift %.2f, %.0f s\n", oracle, kind, seed,
                row$maxAbsZ, paste(round(res$z, 2), collapse = " "), res$shift, sec))
    write.table(row, outFile, sep = ",", row.names = FALSE,
                col.names = !file.exists(outFile), append = file.exists(outFile))
  }
}
utils::assignInNamespace(".Shifted", Original, "EvoTempo")
