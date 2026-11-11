# Data-free (prior-only) comparison of two RevBayes builds
#   Rscript prior-summary.R <dir> <prefix> <stat1,stat2,...>
args <- commandArgs(TRUE); dir <- args[1]; prefix <- args[2]; stats <- strsplit(args[3], ",")[[1]]
ESS <- function(x) { n <- length(x); if (var(x) == 0) return(0)
  fit <- ar(x, aic = TRUE); n * var(x) / (fit$var.pred / (1 - sum(fit$ar))^2) }
out <- list()
for (b in c("v141", "dev")) {
  files <- Sys.glob(file.path(dir, b, "output", sprintf("%s_morph_run*.log", prefix)))
  logs <- lapply(files, function(f) { d <- read.delim(f); d <- d[-nrow(d), ]; d[-seq_len(floor(nrow(d) * 0.1)), ] })
  for (s in stats) {
    ch <- lapply(logs, `[[`, s); ess <- vapply(ch, ESS, 0)
    mcse <- sqrt(sum(vapply(seq_along(ch), function(i) var(ch[[i]]) / max(ess[i], 1), 0))) / length(ch)
    out[[length(out) + 1]] <- data.frame(build = b, stat = s, chains = length(ch),
      last_iteration = min(vapply(logs, function(d) max(d$Iteration), 0)), samples_per_chain = min(lengths(ch)),
      mean = mean(unlist(ch)), mcse = mcse, ess_total = sum(ess),
      chain_means = paste(signif(vapply(ch, mean, 0), 4), collapse = " "))
  }
}
out <- do.call(rbind, out); print(out, digits = 4)
for (s in stats) { a <- out[out$build == "v141" & out$stat == s, ]; d <- out[out$build == "dev" & out$stat == s, ]
  cat(sprintf("%-18s v141 %.4g +- %.2g   dev %.4g +- %.2g   z %+.1f\n", s, a$mean, a$mcse, d$mean, d$mcse,
              (a$mean - d$mean) / sqrt(a$mcse^2 + d$mcse^2))) }
write.csv(out, file.path(dir, "prior-summary.csv"), row.names = FALSE)

