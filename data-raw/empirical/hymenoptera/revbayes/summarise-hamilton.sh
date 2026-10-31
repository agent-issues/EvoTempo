#!/bin/bash
# Summaries of the hymenoptera RevBayes chains on Hamilton.
#
#   bash summarise-hamilton.sh <run root> [burn-in fraction]
#
# <run root> holds one directory per chain, <config>_run<rep>/output/.  Run from
# a copy of this directory (it needs data/, mcc.Rev, summarise.R and the
# published*.csv files).  Writes ./summary/: pooled post-burn-in tree files
# (large; not committed), MCC trees and the CSVs of summarise.R.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00

set +e; set +o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$1
BURNIN=${2:-0.25}
CONFIGS=${CONFIGS:-"te-shared te-separate morph"}
mkdir -p summary flat
for config in $CONFIGS; do
  ls -d "$ROOT"/${config}_run*/output >/dev/null 2>&1 || continue
  # flat copies of every chain's log and tree file: drop the last line (may be
  # incomplete while a chain runs) and any line with the wrong number of
  # fields (a line cut off when a job hit its time limit, followed by the
  # appended output of the job that resumed from the checkpoint); where a
  # resumed segment repeats iterations already logged after the checkpoint,
  # keep the resumed copy (the chain continues from the checkpoint state).
  for d in "$ROOT"/${config}_run*/output; do
    for f in "$d"/hymenoptera_${config}_run*.log "$d"/hymenoptera_${config}_run*.trees; do
      case "$f" in *_rates.trees) continue ;; esac
      [ -f "$f" ] || continue
      out=flat/$(basename "$f")
      head -n 1 "$f" > "$out"
      head -n -1 "$f" | awk -F'\t' 'NR == 1 { nf = NF; next } NF == nf && $1 != "Iteration"' | tac | awk -F'\t' '!seen[$1]++' | tac >> "$out"
    done
  done
  # pool post-burn-in trees (one header line per file)
  pooled=summary/pooled_${config}.trees
  first=1
  for f in flat/hymenoptera_${config}_run*.trees; do
    n=$(($(wc -l < "$f") - 1))
    skip=$(awk -v n=$n -v b=$BURNIN 'BEGIN { print int(n * b) }')
    if [ $first -eq 1 ]; then head -n 1 "$f" > "$pooled"; first=0; fi
    tail -n +$((skip + 2)) "$f" >> "$pooled"
  done
  module purge
  module load gcc/14.2 boost/1.88.0
  printf 'config = "%s"\nsource("mcc.Rev")\nq()\n' "$config" > summary/mcc_${config}_driver.Rev
  "${RB:-/nobackup/pjjg18/revbayes-1.4.1/projects/cmake/build/rb}" summary/mcc_${config}_driver.Rev > summary/mcc_${config}.out 2>&1
  module purge
  module load r/4.5.1
  R_LIBS=/nobackup/pjjg18/evotempo-dev/lib Rscript summarise.R flat hymenoptera ${config} ${BURNIN}
done
ls -lh summary
