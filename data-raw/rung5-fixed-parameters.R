# Rung 5 diagnostic: as rung5-validation.R, but with parameters fixed at the
# prior medians, separating tree moves from parameter moves.
# Run from the package root.
suppressMessages(pkgload::load_all(".", quiet = TRUE)); library(parallel)
args <- "5"; lines <- readLines("data-raw/rung5-validation.R")
eval(parse(text = lines[grep("^  suffix <- ", lines):(grep("^  Simulate <- function", lines) - 1)]))
fixed <- as.list(exp(meanlog5))
SimFixed <- function(seed, n) {
  set.seed(seed); out <- list(); set <- .Setup(fixed, model5)
  while (length(out) < n) {
    k <- .Sample1(set$prior)
    h <- lapply(1:2, function(i) SimulateDCET(crown, set$grid$x, set$Q, set$lambda, set$mu, set$psi, start = k, rho = model5$rho, maxEdges = 5000))
    ext <- vapply(h, function(hi) sum(hi$edges$sampled %in% TRUE), 0); fos <- vapply(h, function(hi) nrow(hi$fossils), 0)
    if (any(ext == 0) || sum(ext) != nExtant || sum(fos) != nFossil) next
    joint <- Join(h[[1]], h[[2]], set$grid$x[k])
    tree <- ReconstructedTree(joint, gamma = 1); inv <- ReconstructedTree(joint, gamma = -1)
    out[[length(out) + 1L]] <- c(rootBin = k, .TreeStats(TempoTree(tree, present = crown), attr(tree, "characterLength"), attr(inv, "characterLength")))
  }
  do.call(rbind, out)
}
r5 <- readRDS(paste0("data-raw/rung5-results/rung5", suffix, ".rds"))
tt5 <- TempoTree(r5$startTree, present = crown)
tt5 <- TempoTree(r5$startTree, present = crown, ranges = cbind(0, ifelse(tt5$tipType == "extant", 0, crown)))
cl <- makeCluster(16); invisible(clusterCall(cl, function() suppressMessages(pkgload::load_all(".", quiet = TRUE))))
clusterExport(cl, c("model5", "crown", "nExtant", "nFossil", "Join", "fixed", "SimFixed", "tt5"))
res <- parLapply(cl, 1:16, function(i) {
  if (i <= 8) return(SimFixed(900 + i, 250))
  set.seed(i)
  TempoMCMC(tt5, model5, fixed, list(), nIter = 4000, burnin = 200, nSlice = 16,
            moves = c(tree = 1, node = 1, tip = 1, branch = 1, nodeAge = 2, fossilAge = 2, slide = 2, saJump = 2))$samples
})
stopCluster(cl)
saveRDS(res, paste0("data-raw/rung5-results/rung5-fixed-parameters", suffix, ".rds"))
sim <- as.data.frame(do.call(rbind, res[1:8])); ch <- res[9:16]
for (q in c("treeLength", "meanNodeAge", "oldestFossil", "youngestFossil", "clock", "inverseClock", "nAncestor", "cherries")) {
  per <- sapply(ch, function(s) mean(s[[q]], na.rm = TRUE))
  cat(sprintf("%-15s sim %.4f (se %.4f)  mcmc %.4f (se %.4f)  z %.2f\n", q, mean(sim[[q]], na.rm = TRUE),
      sd(sim[[q]], na.rm = TRUE) / sqrt(nrow(sim)), mean(per), sd(per) / sqrt(8),
      (mean(per) - mean(sim[[q]], na.rm = TRUE)) / sqrt(var(per) / 8 + var(sim[[q]], na.rm = TRUE) / nrow(sim))))
}
