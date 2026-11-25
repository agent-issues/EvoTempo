# #110 pilot of Metropolis-coupled chains (mc3-pilot.sh): topology mixing per
# arm, from the cold chain's sampled trees reduced to the extant taxa.
#
#   Rscript data-raw/empirical/penguins/cet/mc3-pilot-summary.R <pilot root> [burn]
#
# From the package root.  <pilot root> holds one directory per arm, each with
# b0g0_run<r>/ from run-chain.R; `burn` is the fraction of each chain's samples
# dropped (default 0.25).  Per arm it prints:
# - samples and CPU hours (sum of the segments' seconds), likelihood
#   evaluations in units of chain-iterations;
# - clades of the extant taxa with pooled frequency in [0.1, 0.9]: the
#   largest and mean range of a clade's frequency over replicates, and the
#   ASDSF over clades reaching 0.1 in any replicate;
# - the effective sample size of each such clade's indicator (Ess() of
#   common.R, summed over replicates): median and minimum, and per CPU hour;
# - acceptance of `slide`, and swap rates per adjacent pair.
suppressMessages(library("ape"))
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
burn <- if (length(args) > 1) as.numeric(args[[2]]) else 0.25
src <- readLines("data-raw/empirical/penguins/cet/common.R")
eval(parse(text = src[seq(grep("^LogSumExp <- ", src), length(src))]))
tips <- utils::read.csv("data-raw/empirical/penguins/penguins_tips.csv")
extant <- tips$taxon[tips$status == "extant"]
Clades <- function(newick) {
  tr <- keep.tip(read.tree(text = newick), extant)
  pp <- prop.part(tr)
  lab <- attr(pp, "labels")
  vapply(pp[-1], function(s) paste(sort(lab[s]), collapse = ","), "")
}
Arm <- function(dir) {
  runs <- list.files(dir, pattern = "_run[0-9]+$", full.names = TRUE)
  chains <- lapply(runs, function(run) {
    files <- sort(list.files(run, pattern = "^seg[0-9]+_trees[.]nwk$", full.names = TRUE))
    trees <- unlist(lapply(files, readLines))
    acc <- utils::read.csv(file.path(run, "acceptance.csv"))
    swapFile <- file.path(run, "swaps.csv")
    list(clades = lapply(trees[seq_along(trees) > floor(burn * length(trees))], Clades),
         hours = sum(acc$secs[acc$move == "slide"]) / 3600,
         slide = mean(acc$rate[acc$move == "slide" & acc$burnin == 0]),
         swaps = if (file.exists(swapFile)) {
           s <- utils::read.csv(swapFile)
           tapply(s$accepted, s$colder, sum) / tapply(s$tries, s$colder, sum)
         })
  })
  keys <- unique(unlist(lapply(chains, function(ch) unique(unlist(ch$clades)))))
  ind <- lapply(chains, function(ch) {
    vapply(ch$clades, function(cl) keys %in% cl, logical(length(keys)))
  })
  freq <- vapply(ind, rowMeans, numeric(length(keys)))
  open <- rowMeans(freq) >= 0.1 & rowMeans(freq) <= 0.9
  range <- apply(freq[open, , drop = FALSE], 1, function(f) max(f) - min(f))
  ess <- rowSums(vapply(ind, function(m) {
    apply(m[open, , drop = FALSE], 1, function(x) if (length(unique(x)) > 1) Ess(as.numeric(x)) else 0)
  }, numeric(sum(open))))
  hours <- sum(vapply(chains, `[[`, 0, "hours"))
  swaps <- do.call(rbind, lapply(chains, `[[`, "swaps"))
  list(
    row = data.frame(
      arm = basename(dir), runs = length(chains), samples = sum(vapply(ind, ncol, 0)),
      cpuHours = round(hours, 2), openClades = sum(open),
      maxRange = round(max(range), 3), meanRange = round(mean(range), 3),
      asdsf = round(mean(apply(freq[apply(freq, 1, max) >= 0.1, , drop = FALSE], 1, stats::sd)), 4),
      medianEss = round(stats::median(ess)), minEss = round(min(ess)),
      medianEssPerHour = round(stats::median(ess) / hours, 1),
      minEssPerHour = round(min(ess) / hours, 1),
      slide = round(mean(vapply(chains, `[[`, 0, "slide")), 4)
    ),
    swaps = if (!is.null(swaps)) round(colMeans(swaps), 3),
    freq = data.frame(clade = keys[open], round(freq[open, , drop = FALSE], 3),
                      ess = round(ess))
  )
}
arms <- lapply(list.dirs(root, recursive = FALSE), Arm)
print(do.call(rbind, lapply(arms, `[[`, "row")), row.names = FALSE)
for (a in arms) {
  cat("\n", a$row$arm, ": swap rates by colder heat ",
      if (is.null(a$swaps)) "-" else paste(signif(as.numeric(names(a$swaps)), 3), a$swaps, sep = ": ", collapse = "; "),
      "\n", sep = "")
  a$freq$clade <- abbreviate(a$freq$clade, 40)
  print(a$freq, row.names = FALSE)
}
