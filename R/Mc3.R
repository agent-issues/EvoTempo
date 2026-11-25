# Metropolis-coupled chains of TempoMCMC() (#110; dev/plans/mc3.md).
# Level r of the ladder has heat h_r = 1 / (1 + deltaHeat (r - 1)) on the
# character likelihood.  States stay in their slot (`chain`), tuning stays
# with its level (`tune`), and `slot[r]` names the state at level r, so a swap
# exchanges two entries of `slot` only.
.Mc3Init <- function(st, nChain, deltaHeat, tune, CharLL) {
  heats <- 1 / (1 + deltaHeat * (seq_len(nChain) - 1))
  chain <- c(list(st), lapply(heats[-1], function(h) {
    draw <- .SampleTree(st[["tree"]], st[["set"]], st[["msg"]])
    st[["paths"]] <- draw[["paths"]]
    st[["bin"]] <- draw[["bin"]]
    st[["charLL"]] <- CharLL(st[["tree"]], st[["paths"]], st[["par"]], st[["set"]], h)
    st
  }))
  list(heats = heats, chain = chain, slot = seq_len(nChain),
       tune = rep(list(tune), nChain),
       swapTry = numeric(nChain - 1L), swapAcc = numeric(nChain - 1L))
}

# One swap attempt between a uniformly chosen pair of adjacent levels a, a + 1,
# accepted with probability min(1, (L_B / L_A)^(h_a - h_b)), L the unheated
# character likelihoods of the states A and B at the two levels
.SwapHeats <- function(mc3, CharLL) {
  a <- 1L + floor(stats::runif(1) * length(mc3[["swapTry"]]))
  u <- stats::runif(1)
  pair <- c(a, a + 1L)
  h <- mc3[["heats"]][pair]
  id <- mc3[["slot"]][pair]
  logLik <- vapply(mc3[["chain"]][id], `[[`, 0, "charLL") / h
  logRatio <- (h[[1]] - h[[2]]) * (logLik[[2]] - logLik[[1]])
  mc3[["swapTry"]][a] <- mc3[["swapTry"]][a] + 1
  if (is.finite(logRatio) && log(u) < logRatio) {
    mc3[["swapAcc"]][a] <- mc3[["swapAcc"]][a] + 1
    mc3[["slot"]][pair] <- rev(id)
    for (k in 1:2) {
      st <- mc3[["chain"]][[id[[k]]]]
      st[["charLL"]] <- CharLL(st[["tree"]], st[["paths"]], st[["par"]], st[["set"]],
                               h[[3L - k]])
      mc3[["chain"]][[id[[k]]]] <- st
    }
  }
  # Return:
  mc3
}

.SwapTable <- function(mc3) {
  n <- length(mc3[["swapTry"]])
  data.frame(
    hotter = mc3[["heats"]][-1], colder = mc3[["heats"]][seq_len(n)],
    tries = mc3[["swapTry"]], accepted = mc3[["swapAcc"]],
    rate = mc3[["swapAcc"]] / pmax(mc3[["swapTry"]], 1)
  )
}
