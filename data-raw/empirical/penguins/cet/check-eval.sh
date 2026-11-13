#!/bin/bash
# Penguin CET pipeline: eval's side of #104 on a real chain (one with
# seg*_states.rds and a checkpoint, e.g. a copy of a smoke chain).
#
#   bash data-raw/empirical/penguins/cet/check-eval.sh <chain dir> [out.txt]
#
# Run from the package root with EvoTempo installed.  Works on <chain dir> in
# place: use a copy.  Checks:
# 1. With CET_BF_TOL = 1e-15 (below any real agreement), eval exits non-zero,
#    removes a stale cross.csv left from an earlier run, and writes none
#    (before #104 it wrote cross.csv before its check, and never failed).
# 2. With the default tolerance, eval exits 0 and writes cross.csv with check
#    and checkChar columns, every |value| within CET_BF_TOL (default 1e-4),
#    and leaves no temporary file behind.
set -u
d=$1
out=${2:-/dev/stdout}
C=data-raw/empirical/penguins/cet
fail=0
{
  echo "check-eval.sh on $d ($(git rev-parse --short HEAD 2>/dev/null))"
  echo stale > "$d/cross.csv"
  CET_BF_TOL=1e-15 Rscript $C/bayes-factors.R eval "$d" > "$d/../check-eval-1.log" 2>&1
  s1=$?
  tail -2 "$d/../check-eval-1.log"
  if [ -e "$d/cross.csv" ]; then f1=present; else f1=absent; fi
  echo "1. CET_BF_TOL=1e-15: exit $s1, cross.csv $f1"
  if [ "$s1" -ne 0 ] && [ "$f1" = absent ]; then echo "ok   1. a failed check leaves no cross.csv"; else echo "FAIL 1. a failed check leaves no cross.csv"; fail=1; fi
  Rscript $C/bayes-factors.R eval "$d" > "$d/../check-eval-2.log" 2>&1
  s2=$?
  tail -1 "$d/../check-eval-2.log"
  Rscript -e 'x <- read.csv(file.path(commandArgs(TRUE)[1], "cross.csv")); ok <- all(c("check", "checkChar") %in% names(x)) && all(abs(as.matrix(x[, c("check", "checkChar")])) <= 1e-4); cat(sprintf("2. cross.csv: %d rows; max |check| %.3g; max |checkChar| %.3g\n", nrow(x), max(abs(x$check)), max(abs(x$checkChar)))); quit(status = if (ok) 0 else 1)' "$d"
  s3=$?
  ntmp=$(ls "$d" | grep -c '^cross-.*\.csv$')
  if [ "$s2" -eq 0 ] && [ "$s3" -eq 0 ] && [ "$ntmp" -eq 0 ]; then echo "ok   2. passing checks write cross.csv with check and checkChar, no temporary file"; else echo "FAIL 2. (eval exit $s2, column check $s3, temporary files $ntmp)"; fail=1; fi
  if [ $fail -eq 0 ]; then echo "check-eval.sh: all checks passed"; else echo "check-eval.sh: FAILED"; fi
} > "$out"
exit $fail
