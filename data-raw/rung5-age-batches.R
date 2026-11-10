# Rung 5 with a free start age (#18, #57): per-batch summaries of the
# originFree fixed-parameter replicates of section 3 of
# data-raw/rung5-age-results/README.md, from the committed per-replicate z of
# the means (rung5-age-pooled-per-run.csv).  Batches: A = seed offsets
# 6000-45000, B = 46000-125000, C = 126000-245000.  "stouffer" is the sum of the
# per-replicate z over sqrt(number of replicates).  The batch-level pooled
# quantile tests quoted in the README were computed from the per-replicate
# .rds files, which are no longer kept, and cannot be regenerated from these.
# Usage: Rscript data-raw/rung5-age-batches.R [dir]
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[[1]] else "data-raw/rung5-age-results"
perRun <- read.csv(file.path(dir, "rung5-age-pooled-per-run.csv"))
perRun <- perRun[perRun$case == "originFree-fixed", ]
offset <- as.numeric(sub(".*-seed", "", perRun$run))
perRun$batch <- cut(offset, c(5999, 45000, 125000, 245000), c("A", "B", "C"))
perRun <- perRun[!is.na(perRun$batch), ]
out <- do.call(rbind, lapply(split(perRun, list(perRun$batch, perRun$quantity), drop = TRUE),
                             function(x) {
  data.frame(batch = x$batch[[1]], quantity = x$quantity[[1]], nRep = nrow(x),
             meanZ = mean(x$z), sdZ = stats::sd(x$z), stouffer = sum(x$z) / sqrt(nrow(x)))
}))
out <- out[order(out$batch, out$quantity), ]
write.csv(out, file.path(dir, "rung5-age-pooled-per-batch.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
