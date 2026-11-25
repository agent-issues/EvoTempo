#!/bin/bash
# Run root for the CET character-power ladders (ml-ladder.sh): links to the
# production chains of the ladder cells, the t = 1 rung that ml-summary.R reads.
#
#   bash ml-root.sh <production run root> <ladder run root> <cells>
#   bash ml-root.sh $W/runs $W/ml-runs "b0g0 b1g0"
#
# Stops if a cell has no production chains, rather than leaving a dangling
# link (review of #111, R111-3); ml-summary.R also checks the links.
set -e
prod=$1; root=$2; cells=$3
mkdir -p $root
for cell in $cells; do
  n=0
  for d in $prod/${cell}_run*; do
    [ -d "$d" ] || { echo "no production chains $prod/${cell}_run*"; exit 1; }
    ln -sfn $d $root/$(basename $d); n=$(( n + 1 ))
  done
  echo "$cell: $n production chains linked"
done
