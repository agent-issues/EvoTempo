# Rung 3 (#15) summary: tables from data-raw/rung3-revbayes.R runs.
# Run from the package root after those runs (tags -full at 16000 slices,
# -ladder, -beta0, -point5 and -point5-survival).
dir <- "data-raw/rung3-results/revbayes"
Read <- function(tag) readRDS(file.path(dir, paste0("rung3-revbayes", tag, ".rds")))$res
full <- Read("-full")
p5 <- rbind(Read("-point5"), Read("-point5-survival"))
# The stiffest K = 26 point converges only by 64000 slices: use that run there
p5 <- p5[p5$nTimeSlices == 64000, ]
key <- function(d) paste(d$point, d$tree, d$condition, d$grid)
full <- full[!(key(full) %in% key(p5)), ]
best <- rbind(full, p5[, names(full)])
table <- do.call(rbind, lapply(split(best, list(best$grid, best$tree, best$condition), drop = TRUE),
  function(d) {
    worst <- which.max(abs(d$delta))
    data.frame(grid = d$grid[1], tree = d$tree[1], condition = d$condition[1],
               nCase = nrow(d), nTimeSlices = paste(sort(unique(d$nTimeSlices)), collapse = "/"),
               maxAbsDelta = max(abs(d$delta)), worstPoint = d$point[worst],
               meanDelta = mean(d$delta))
  }))
table$verdict <- ifelse(table$maxAbsDelta < 1e-6, "pass",
                        ifelse(table$tree %in% c("sa", "sa2"),
                               "fail: dnCDBDP scores sampled ancestors psi * E(t)", "fail"))
ladder <- Read("-ladder")
convergence <- aggregate(abs(delta) ~ grid + tree + nTimeSlices, ladder, max)
names(convergence)[4] <- "maxAbsDelta"
b0 <- Read("-beta0")
b0$residual <- b0$delta - ifelse(is.na(b0$saOffset), 0, b0$saOffset)
saCheck <- aggregate(cbind(delta = abs(delta), residual = abs(residual)) ~ tree + nTimeSlices,
                     b0, max)
write.csv(table, file.path(dir, "table.csv"), row.names = FALSE)
write.csv(convergence, file.path(dir, "convergence.csv"), row.names = FALSE)
write.csv(saCheck, file.path(dir, "sa-diagnosis.csv"), row.names = FALSE)
print(table[order(table$tree, table$grid, table$condition), ], digits = 3, row.names = FALSE)
print(convergence[order(convergence$tree, convergence$grid, convergence$nTimeSlices), ], digits = 3, row.names = FALSE)
print(saCheck, digits = 3, row.names = FALSE)
