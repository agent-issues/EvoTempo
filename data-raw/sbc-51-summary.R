# Summary of the independent SBC of the TempoMoves() schedule (#51;
# data-raw/sbc-51-array.sh), cell beta = 0, gamma = 1.  Writes, in
# data-raw/mixing-results/sbc-51/:
#   sbc-51-tests.csv   rank uniformity (chi-square on 10 bins, KS on jittered
#                      ranks; Holm over the 13 quantities of each subset) for
#                      the 150 datasets new to the schedule (g from 202), the
#                      50 of #19 (g up to 198), all 200, and the long arm;
#   sbc-51-paired.csv  paired comparisons on the same datasets: standard arm vs
#                      #21 (200), long vs standard arm (50).  Mean rank
#                      difference with paired t, and z = difference of
#                      posterior means / combined MCSE, whose sd is about 1
#                      when the MCSEs are right (below 1 for long vs standard,
#                      whose chains share their first 20 000 iterations).
# Usage: Rscript data-raw/sbc-51-summary.R <standard dir> [long dir]
args <- commandArgs(trailingOnly = TRUE)
outDir <- "data-raw/mixing-results/sbc-51"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
nRank <- 99L

Read <- function(dir) {
  files <- list.files(file.path(dir, "replicates"), "^rep[0-9]+\\.rds$", full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    r <- readRDS(f)
    q <- names(r$ranks)
    data.frame(g = r$g, quantity = q, rank = unname(r$ranks),
               mean = unname(r$posteriorMean[q]), mcse = unname(r$posteriorMcse[q]),
               nIter = r$nIter, runtimeMin = r$runtime)
  }))
}
standard <- Read(args[[1]])
long <- if (length(args) > 1) Read(args[[2]]) else NULL
quantities <- unique(standard$quantity)
old <- utils::read.csv("data-raw/rung6-full-results/rung6-full-ranks.csv")
old <- old[old$cell == 2, ]
old <- data.frame(g = old$g, quantity = old$quantity, rank = old$rank,
                  mean = old$posteriorMean, mcse = old$posteriorSd / sqrt(old$ess))

Tests <- function(d, label) {
  out <- do.call(rbind, lapply(quantities, function(q) {
    rk <- d$rank[d$quantity == q & !is.na(d$rank)]
    counts <- tabulate(rk %/% ((nRank + 1L) / 10L) + 1L, 10L)
    u <- (rk + stats::runif(length(rk))) / (nRank + 1L)
    data.frame(subset = label, quantity = q, n = length(rk), meanRank = mean(rk),
               pChisq = suppressWarnings(stats::chisq.test(counts)$p.value),
               pKS = stats::ks.test(u, "punif")$p.value)
  }))
  out$pChisqHolm <- stats::p.adjust(out$pChisq, "holm")
  out$pKSHolm <- stats::p.adjust(out$pKS, "holm")
  out$pass <- out$pChisqHolm >= 0.05 & out$pKSHolm >= 0.05
  out
}
set.seed(5101)
tests <- rbind(Tests(standard[standard$g >= 202, ], "standard, new datasets (g >= 202)"),
               Tests(standard[standard$g <= 198, ], "standard, #19 datasets (g <= 198)"),
               Tests(standard, "standard, all"),
               if (!is.null(long)) Tests(long, "long (10x), #19 datasets"))
utils::write.csv(tests, file.path(outDir, "sbc-51-tests.csv"), row.names = FALSE)

Paired <- function(a, b, label) {
  m <- merge(a, b, by = c("g", "quantity"), suffixes = c("A", "B"))
  out <- do.call(rbind, lapply(quantities, function(q) {
    d <- m[m$quantity == q & !is.na(m$rankA) & !is.na(m$rankB), ]
    diff <- d$rankA - d$rankB
    z <- (d$meanA - d$meanB) / sqrt(d$mcseA^2 + d$mcseB^2)
    z <- z[is.finite(z)]
    data.frame(comparison = label, quantity = q, n = nrow(d),
               meanRankDiff = mean(diff),
               pPaired = if (stats::sd(diff) > 0) stats::t.test(diff)$p.value else NA_real_,
               rankCor = stats::cor(d$rankA, d$rankB, method = "spearman"),
               sdZ = stats::sd(z), fracAbsZOver2 = mean(abs(z) > 2))
  }))
  out$pPairedHolm <- stats::p.adjust(out$pPaired, "holm")
  out
}
paired <- rbind(Paired(standard, old, "standard - #21"),
                if (!is.null(long)) Paired(long, standard, "long - standard"))
utils::write.csv(paired, file.path(outDir, "sbc-51-paired.csv"), row.names = FALSE)

cost <- do.call(rbind, lapply(list(standard = standard, long = long), function(d) {
  if (is.null(d)) return(NULL)
  r <- d[!duplicated(d$g), ]
  data.frame(n = nrow(r), medianIter = stats::median(r$nIter),
             medianMinutes = stats::median(r$runtimeMin),
             maxMinutes = max(r$runtimeMin), coreHours = sum(r$runtimeMin) / 60)
}))
print(tests[!tests$pass | tests$pChisq < 0.05 | tests$pKS < 0.05, ], digits = 3)
print(paired, digits = 3)
print(cost, digits = 3)
