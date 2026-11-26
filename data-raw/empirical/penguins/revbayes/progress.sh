#!/bin/bash
# Last logged state of every chain under a run root:
#   bash progress.sh /nobackup/pjjg18/evotempo-penguins-rbdev-r4
# prints <chain> n=<samples> last=<iteration> | Likelihood crown_age num_samp_anc and the clock log-SD(s)
ROOT=${1:-.}
date
for d in "$ROOT"/*_run*; do b=$(basename $d); f=$d/output/penguins_${b}.log; [ -f $f ] || continue
  h=$(head -1 $f | tr '\t' '\n' | grep -n -x -e Likelihood -e crown_age -e num_samp_anc -e ucld_sd -e ucld_sd_dna -e ucld_sd_morph | cut -d: -f1 | paste -sd,)
  echo "$b n=$(($(wc -l < $f)-1)) last=$(tail -1 $f | cut -f1) | $(tail -1 $f | cut -f$h)"
done
