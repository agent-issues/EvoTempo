#!/bin/bash
# Run root for the CET character-power ladders (ml-ladder.sh): links to the
# production chains of the ladder cells (the t = 1 rung) and a launch record
# for `bayes-factors.R combine`.
#
#   bash ml-root.sh <production run root> <ladder run root> <cells> <ladder reps> [npow]
#   bash ml-root.sh $W/runs $W/ml-runs "b0g0 b1g0" 2
set -e
prod=$1; root=$2; cells=$3; reps=$4; npow=${5:-64}
mkdir -p $root
echo "cell,reps,power" > $root/launch.csv
for cell in $cells; do
  n=0
  for d in $prod/${cell}_run*; do ln -sfn $d $root/$(basename $d); n=$(( n + 1 )); done
  echo "$cell,$n,1" >> $root/launch.csv
  for j in $(seq 0 $(( npow - 1 ))); do
    t=$(Rscript -e "cat(format(($j / $npow)^5, digits = 17, scientific = FALSE))")
    echo "$cell,$reps,$t" >> $root/launch.csv
  done
done
cat $root/launch.csv | wc -l
