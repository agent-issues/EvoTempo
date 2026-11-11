#!/bin/bash
# Hymenoptera CET pipeline: run the fixed-tree RevBayes oracle (check-mkv.Rev)
# on the RevBayes development builds, for check-mkv.R.
#
#   bash data-raw/empirical/hymenoptera/cet/revbayes/run-check-mkv.sh <scratch dir> [rbmorph data dir]
#
# From the package root.  Regenerates the comparator's data/ with
# ../../revbayes/prepare-data.R (and, if given, checks it byte for byte
# against the data/ directory of the comparator's chains), copies it with the
# CET start tree into <scratch dir>, and writes
# cet/smoke-results/revbayes-check-mkv-<build>.txt for each build.
set +e; set +o pipefail
SCR=$1
RBDATA=$2
CET=data-raw/empirical/hymenoptera/cet
ROOT=$(pwd)
mkdir -p $SCR
(cd data-raw/empirical/hymenoptera/revbayes && Rscript prepare-data.R > $SCR/prepare-data.txt 2>&1)
if [ -n "$RBDATA" ]; then
  (cd data-raw/empirical/hymenoptera/revbayes/data && md5sum morph_*.nex partitions.Rev taxa.tsv) > $SCR/md5-regenerated.txt
  (cd $RBDATA && md5sum morph_*.nex partitions.Rev taxa.tsv) > $SCR/md5-chains.txt
  if cmp -s $SCR/md5-regenerated.txt $SCR/md5-chains.txt; then
    echo "comparator data regenerated identically ($(wc -l < $SCR/md5-chains.txt) files)" > $ROOT/$CET/smoke-results/revbayes-data-check.txt
  else
    echo "comparator data DIFFER from the chains' data/" > $ROOT/$CET/smoke-results/revbayes-data-check.txt
    diff $SCR/md5-regenerated.txt $SCR/md5-chains.txt >> $ROOT/$CET/smoke-results/revbayes-data-check.txt
  fi
fi
rm -rf $SCR/rb && mkdir -p $SCR/rb
cp -r data-raw/empirical/hymenoptera/revbayes/data $SCR/rb/data
cp $CET/data/start-tree.nwk $SCR/rb/start_tree.nwk
cp $CET/revbayes/check-mkv.Rev $SCR/rb/
for build in fd807d10 a2a8e2db; do
  RB=/nobackup/pjjg18/revbayes-dev-$build/projects/cmake/build/rb
  (source /etc/profile > /dev/null 2>&1 || true; module purge || true
   module load gcc/14.2 boost/1.88.0
   cd $SCR/rb && $RB --version > $SCR/rb-version-$build.txt 2>&1
   $RB check-mkv.Rev > $ROOT/$CET/smoke-results/revbayes-check-mkv-$build.txt 2>&1)
done
