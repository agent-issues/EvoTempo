# verify43c / #72: per-chain and per-run summaries of every completed run
W <- "/nobackup/pjjg18/evotempo-bias72"
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
pc <- list(); pr <- list()
for (arm in c("A", "B", "C")) {
  for (f in list.files(file.path(W, "out", arm), pattern = "\\.rds$", full.names = TRUE)) {
    x <- readRDS(f)
    off <- as.integer(sub(".*-seed([0-9]+).*", "\\1", basename(f)))
    q <- x$out$exact[!is.na(x$out$p)]; em <- x$out$exact[is.na(x$out$p)]
    stopifnot(isTRUE(all.equal(x$out$p[!is.na(x$out$p)], probs)))
    Fm <- t(vapply(x$chains, function(ch) vapply(q, function(qq) mean(ch$age <= qq), 0), numeric(5)))
    colnames(Fm) <- paste0("F", c("05", "25", "50", "75", "95"))
    pc[[f]] <- data.frame(arm = arm, case = x$case, offset = off, chain = seq_along(x$chains),
                          seed = 118900 + off + seq_along(x$chains),
                          nSample = vapply(x$chains, function(ch) length(ch$age), 0),
                          mean = vapply(x$chains, function(ch) mean(ch$age), 0),
                          exactMean = em, Fm,
                          q05 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q95 = q[5])
    acc <- x$acceptance[grepl("Age", names(x$acceptance))]
    pr[[f]] <- data.frame(arm = arm, case = x$case, offset = off, nChain = length(x$chains),
                          nIter = x$nIter, sha = x$sha, runtimeMin = as.numeric(x$runtime),
                          zMeanScript = x$out$z[is.na(x$out$p)], ownMaxAbsZ = max(abs(x$out$z)),
                          ageMoveAcceptance = if (length(acc)) acc[[1]] else NA)
  }
}
write.csv(do.call(rbind, pc), "bias72-perchain.csv", row.names = FALSE)
write.csv(do.call(rbind, pr), "bias72-perrun.csv", row.names = FALSE)
cat("runs:", length(pr), " chains:", nrow(do.call(rbind, pc)), "\n")
