# Rung 3 (#15): re-evaluate the EvoTempo side of the stored dnCDBDP
# comparisons with each available TreeLikelihood() engine ("R" and, when it
# is compiled, "cpp"), against the stored RevBayes values.
# EVOTEMPO_NO_COMPILE=true loads the package without compiling (machines
# without a toolchain); the C++ engine is then skipped.
# Run from the package root after data-raw/rung3-revbayes.R (REVBAYES_DIR
# points at the results when run from another checkout).
noCompile <- nzchar(Sys.getenv("EVOTEMPO_NO_COMPILE"))
suppressWarnings(suppressMessages(
  pkgload::load_all(".", quiet = TRUE, compile = if (noCompile) FALSE else NA)))
dir <- Sys.getenv("REVBAYES_DIR", "data-raw/rung3-results/revbayes")
hasEngine <- "engine" %in% names(formals(TreeLikelihood))
engines <- "R"
if (hasEngine) {
  cppWorks <- tryCatch({
    TreeLikelihood(ape::read.tree(text = "((a:1,b:1):1,c:2);"),
                   list(lambda = 0.5, mu = 0.2, theta = 0.3, s = 1),
                   DCETModel(nBin = 5, zRange = c(-1, 1)), engine = "cpp")
    TRUE
  }, error = function(e) FALSE)
  if (cppWorks) engines <- c("R", "cpp") else message("C++ engine unavailable: R only")
}
RbLnFactorial <- function(n) {
  ifelse(n == 0, 0, (n + 0.5) * log(n) - n + 0.918938533204672722 +
           (1 / n) * (1 / 12 + (1 / n)^2 * (-1 / 360)))
}
grids <- list(K16 = list(nBin = 16, zRange = c(-3, 3), beta = 1),
              K26 = list(nBin = 26, zRange = c(-5, 5), beta = 1),
              K16b0 = list(nBin = 16, zRange = c(-3, 3), beta = 0))
tags <- c("-full", "-point5", "-point5-survival", "-beta0")
out <- do.call(rbind, lapply(tags, function(tag) {
  x <- readRDS(file.path(dir, paste0("rung3-revbayes", tag, ".rds")))
  res <- x$res
  if (tag == "-point5") res <- res[res$nTimeSlices == 64000, ]
  for (engine in engines) {
    res[[paste0("evo_", engine)]] <- vapply(seq_len(nrow(res)), function(i) {
      cs <- res[i, ]
      g <- grids[[cs$grid]]
      par <- c(as.list(x$points[cs$point, ]), psi = 0.15)
      set <- .Setup(par, DCETModel(nBin = g$nBin, zRange = g$zRange, rho = 0.7, beta = g$beta))
      w <- set$prior * set$lambda
      model <- DCETModel(nBin = g$nBin, zRange = g$zRange, rho = 0.7, removal = 0,
                         beta = g$beta, rootPrior = w / sum(w),
                         condition = if (cs$condition == "time") "none" else "survival")
      info <- x$info[x$info$tree == cs$tree, ]
      shape <- (info$nTip - info$nSA - 1) * log(2) - RbLnFactorial(info$nExtant)
      tl <- if (!hasEngine) {
        TreeLikelihood(x$trees[[cs$tree]], par, model)
      } else {
        TreeLikelihood(x$trees[[cs$tree]], par, model, engine = engine)
      }
      tl + shape + if (cs$condition == "time") log(sum(w)) else 0
    }, 0)
  }
  res$tag <- tag
  res
}))
# Drop the 16000-slice point-5 rows superseded by their 64000-slice runs
key <- paste(out$point, out$tree, out$condition, out$grid)
p5 <- out$tag %in% c("-point5", "-point5-survival")
out <- out[!(out$tag == "-full" & key %in% key[p5]), ]
for (engine in engines) {
  out[[paste0("delta_", engine)]] <- out$rb - out[[paste0("evo_", engine)]] -
    ifelse(out$tag == "-beta0" & !is.na(out$saOffset), out$saOffset, 0)
}
write.csv(out, file.path(dir, "engines.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(engines, function(engine) {
  d <- out
  d$ad <- abs(d[[paste0("delta_", engine)]])
  s <- aggregate(ad ~ grid + tree + condition, d, max)
  s$engine <- engine
  s$sa <- ifelse(s$grid == "K16b0", "residual after predicted SA offset", "")
  s
}))
names(summ)[names(summ) == "ad"] <- "maxAbsDelta"
write.csv(summ, file.path(dir, "engines-summary.csv"), row.names = FALSE)
if (length(engines) > 1) {
  cat("max |R - cpp| over cases:", max(abs(out$evo_R - out$evo_cpp)), "\n")
}
print(summ[order(summ$engine, summ$tree, summ$grid, summ$condition), ], digits = 3, row.names = FALSE)
