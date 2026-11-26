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
# 3. The record (#105 review, M9): evalSamples equals the rows of the chain's
#    seg*_samples.csv, evalSegMd5 is "md5:..." and evalTol, evalTolChar equal
#    CET_BF_TOL, CET_BF_TOL_CHAR (defaults 1e-4, 1e-8), and evalMaxSeg is 0 (no
#    segment cap; #105 delta review, C3).
# 4. charLL + 1e-3 in one kept state of seg*_samples.csv (M8): eval exits
#    non-zero and leaves no cross.csv (checkChar, not check, catches it).
# 5. After a good eval, an eval that stops early (CET_BF_TOL=abc) leaves no
#    cross.csv (R105-1).
# 6. With cross.csv a directory, the final rename fails: eval exits non-zero
#    and leaves no temporary cross-*.csv (R105-2; #105 delta review, R105b-1).
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
  Rscript -e 'd <- commandArgs(TRUE)[1]; x <- read.csv(file.path(d, "cross.csv")); n <- sum(vapply(Sys.glob(file.path(d, "seg*_samples.csv")), function(f) nrow(read.csv(f)), 0L)); ok <- all(c("evalSamples", "evalSegMd5", "evalTol", "evalTolChar", "evalMaxSeg") %in% names(x)) && all(x$evalSamples == n) && all(grepl("^md5:[0-9a-f]{32}$", x$evalSegMd5)) && all(x$evalTol == 1e-4) && all(x$evalTolChar == 1e-8) && all(x$evalMaxSeg == 0); cat(sprintf("3. record: evalSamples %s (segment rows %d), evalSegMd5 %s, evalTol %s, evalTolChar %s\n", x$evalSamples[1], n, x$evalSegMd5[1], x$evalTol[1], x$evalTolChar[1])); cat(sprintf("   evalMaxSeg %s\n", x$evalMaxSeg[1])); quit(status = if (isTRUE(ok)) 0 else 1)' "$d"
  if [ $? -eq 0 ]; then echo "ok   3. cross.csv records what eval read and its tolerances"; else echo "FAIL 3. cross.csv record"; fail=1; fi
  last=$(ls "$d"/seg*_samples.csv | tail -1)
  cp "$last" "$d/../check-eval-backup.csv"
  Rscript -e 'f <- commandArgs(TRUE)[1]; x <- read.csv(f); x$charLL[nrow(x)] <- x$charLL[nrow(x)] + 1e-3; write.csv(x, f, row.names = FALSE)' "$last"
  Rscript $C/bayes-factors.R eval "$d" > "$d/../check-eval-4.log" 2>&1
  s4=$?
  cp "$d/../check-eval-backup.csv" "$last"; rm -f "$d/../check-eval-backup.csv"
  tail -1 "$d/../check-eval-4.log"
  if [ -e "$d/cross.csv" ]; then f4=present; else f4=absent; fi
  if [ "$s4" -ne 0 ] && [ "$f4" = absent ] && grep -q "charLL by" "$d/../check-eval-4.log"; then echo "ok   4. charLL + 1e-3 in one state: eval stops, no cross.csv"; else echo "FAIL 4. charLL + 1e-3 (exit $s4, cross.csv $f4)"; fail=1; fi
  Rscript $C/bayes-factors.R eval "$d" > "$d/../check-eval-5a.log" 2>&1
  s5a=$?
  CET_BF_TOL=abc Rscript $C/bayes-factors.R eval "$d" > "$d/../check-eval-5b.log" 2>&1
  s5b=$?
  if [ -e "$d/cross.csv" ]; then f5=present; else f5=absent; fi
  if [ "$s5a" -eq 0 ] && [ "$s5b" -ne 0 ] && [ "$f5" = absent ]; then echo "ok   5. an eval that stops early leaves no earlier cross.csv"; else echo "FAIL 5. (good eval exit $s5a, early stop exit $s5b, cross.csv $f5)"; fail=1; fi
  mkdir "$d/cross.csv"
  Rscript $C/bayes-factors.R eval "$d" > "$d/../check-eval-6.log" 2>&1
  s6=$?
  n6=$(ls "$d" | grep -c '^cross-.*\.csv$')
  rmdir "$d/cross.csv"
  tail -1 "$d/../check-eval-6.log"
  if [ "$s6" -ne 0 ] && [ "$n6" -eq 0 ]; then echo "ok   6. a failed rename leaves no temporary file"; else echo "FAIL 6. (exit $s6, temporary files $n6)"; fail=1; fi
  if [ $fail -eq 0 ]; then echo "check-eval.sh: all checks passed"; else echo "check-eval.sh: FAILED"; fi
} > "$out"
exit $fail
