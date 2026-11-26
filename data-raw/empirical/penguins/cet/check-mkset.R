# The one-call route of PartitionedMk() (MkSetFor(), .MkSetLogLik(); #110):
# 1. equals the route of one call per matrix, exactly, on random lengths;
# 2. under EVOTEMPO_MK_CACHE_CHECK = 1 stops when the handle was registered for
#    other data, or when its index map is wrong; without the check neither is
#    seen, which is why checked runs are part of validating a change here.
#
#   Rscript data-raw/empirical/penguins/cet/check-mkset.R
#
# Run from the package root with EvoTempo installed.
suppressMessages(library("EvoTempo"))
source("data-raw/empirical/penguins/cet/common.R")
stopifnot(".MkSet() is missing from this build" = .hasMkSet)
dat <- LoadPenguins()
tree <- dat$tree
chars <- dat$characters
parts <- attr(chars, "parts")
pat <- attr(chars, "mkv")
Check <- function(on) Sys.setenv(EVOTEMPO_MK_CACHE_CHECK = if (on) "1" else "0")
Failures <- function() get(".MkCache", .ns)()$checkFailures
Message <- function(lengths) {
  tryCatch({
    PartitionedMk(tree, chars, lengths)
    ""
  }, error = conditionMessage)
}

set.seed(1)
Check(FALSE)
for (i in 1:20) {
  len <- runif(tree$nTip + tree$nNode, 0.001, 0.4)
  oneCall <- PartitionedMk(tree, chars, len)
  options(EvoTempo.mkCache = FALSE)
  perMatrix <- PartitionedMk(tree, chars, len)
  options(EvoTempo.mkCache = NULL)
  stopifnot(identical(oneCall, perMatrix))
}

Check(TRUE)
stopifnot(Message(len) == "", Failures() == 0)
intact <- .cet$mkSet

other <- parts
cell <- which(!is.na(other[[1]]$x))[[1]]
other[[1]]$x[cell] <- 1L - other[[1]]$x[cell]
wrong <- MkSetFor(other, pat, tree$nTip)
wrong$parts <- parts
.cet$mkSet <- wrong
stopifnot(grepl("differ from those registered", Message(len)), Failures() == 1)
Check(FALSE)
stopifnot(Message(len) == "")

wrong <- intact
wrong$pattern <- rev(wrong$pattern)
.cet$mkSet <- wrong
stopifnot(Message(len) == "")
Check(TRUE)
stopifnot(grepl("differs from one call per matrix", Message(len)))

.cet$mkSet <- NULL
stopifnot(identical(PartitionedMk(tree, chars, len), oneCall))
cat("check-mkset: PASS\n")
