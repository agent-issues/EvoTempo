#!/bin/bash
# Mutation check of the #105 tests: each mutant is a one-line change to a copy
# of the code; the suite (check-labels.R, check-eval.sh, check-segcap.R) must
# fail on it.  Usage: bash mutants.sh <worktree> <real chain dir to copy>
set -u
src=$1; chain=$2
C=data-raw/empirical/penguins/cet
declare -A FILE SED
FILE[M4]=common.R;  SED[M4]='s/^    bad <- sum(!is.finite(chk)) + sum(!is.finite(chc))$/    bad <- sum(!is.finite(chk))/; s/^    if (max(abs(chk)) > use || max(abs(chc)) > useChar) {$/    if (max(abs(chk)) > use) {/'
FILE[M5]=common.R;  SED[M5]='s/^    bad <- sum(!is.finite(chk)) + sum(!is.finite(chc))$/    bad <- 0/'
FILE[M8]=bayes-factors.R; SED[M8]='s/^  if (!isTRUE(worst <= tol) || !isTRUE(worstChar <= tolChar)) {$/  if (!isTRUE(worst <= tol)) {/'
FILE[M9]=bayes-factors.R; SED[M9]='/^  out\$evalTol <- tol$/d'
FILE[M9b]=bayes-factors.R; SED[M9b]='/^  out\$evalTolChar <- tolChar$/d'
FILE[N1]=bayes-factors.R; SED[N1]='s/^  }, finally = unlink(tmp))$/  })/'
FILE[N2]=common.R;  SED[N2]='s/^    return(isTRUE(seg >= m) \&\& length(SegFiles(dir, strict = FALSE)) == m)$/    return(length(SegFiles(dir, strict = FALSE)) == m)/'
FILE[N3]=common.R;  SED[N3]='s/^  if (strict \&\& !all(have)) {$/  if (FALSE) {/'
FILE[N4]=common.R;  SED[N4]='s/^    if (!identical(suppressWarnings(as.integer(x\$evalMaxSeg\[1\])), capNow)) {$/    if (FALSE) {/'
FILE[N5]=common.R;  SED[N5]='s/^  if (!grepl("\^\[0-9\]+\$", v) || as.numeric(v) < 1) {$/  if (FALSE) {/'
for m in head M4 M5 M8 M9 M9b N1 N2 N3 N4 N5; do
  w=$(mktemp -d); cp -r $src/. $w/; cd $w
  if [ $m != head ]; then
    f=$C/${FILE[$m]}; before=$(md5sum < $f); sed -i "${SED[$m]}" $f
    [ "$(md5sum < $f)" = "$before" ] && { echo "$m: SED DID NOT APPLY"; cd /; rm -rf $w; continue; }
  fi
  Rscript $C/check-labels.R > cl.txt 2>&1; a=$?
  rm -rf $w/chaincopy; mkdir -p $w/chaincopy/root; cp -r $chain $w/chaincopy/root/
  bash $C/check-eval.sh $w/chaincopy/root/$(basename $chain) ce.txt > /dev/null 2>&1; b=$?
  Rscript $C/check-segcap.R > cs.txt 2>&1; s=$?
  fl=$(grep -c '^FAIL' cl.txt); fe=$(grep -c '^FAIL' ce.txt)
  echo "$m: check-labels exit $a ($fl FAIL: $(grep '^FAIL' cl.txt | head -3 | cut -c1-90 | paste -sd'|')); check-eval exit $b ($fe FAIL: $(grep '^FAIL' ce.txt | cut -c1-40 | paste -sd'|')); check-segcap exit $s ($(grep -m1 FAILED cs.txt | cut -c1-90))"
  cd /; rm -rf $w
done
