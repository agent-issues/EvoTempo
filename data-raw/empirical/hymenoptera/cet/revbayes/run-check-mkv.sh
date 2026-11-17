#!/bin/bash
# Hymenoptera CET pipeline: run the fixed-tree RevBayes oracle (check-mkv.Rev)
# for check-mkv.R.
#
#   bash data-raw/empirical/hymenoptera/cet/revbayes/run-check-mkv.sh <scratch dir> [rbmorph data dir]
#
# From the package root.  Regenerates the comparator's data/ with
# ../../revbayes/prepare-data.R (and, if given, checks it byte for byte
# against the data/ directory of the comparator's chains), and runs the
# oracle with the RevBayes development build fd807d10 on two trees: the
# reference start tree (data/start-tree.nwk; tag "reference") and production
# chain 6's start tree (data/start-trees/start-tree-06.nwk; tag "draw06", one
# of the two draws whose conversion moved an age).  Writes
# cet/smoke-results/revbayes-check-mkv-fd807d10-<tag>.txt.
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
build=fd807d10
RB=/nobackup/pjjg18/revbayes-dev-$build/projects/cmake/build/rb
for tag in reference draw06; do
  if [ $tag = reference ]; then tree=$CET/data/start-tree.nwk; else tree=$CET/data/start-trees/start-tree-06.nwk; fi
  rm -rf $SCR/rb-$tag && mkdir -p $SCR/rb-$tag
  cp -r data-raw/empirical/hymenoptera/revbayes/data $SCR/rb-$tag/data
  cp $tree $SCR/rb-$tag/start_tree.nwk
  cp $CET/revbayes/check-mkv.Rev $SCR/rb-$tag/
  (source /etc/profile > /dev/null 2>&1 || true; module purge || true
   module load gcc/14.2 boost/1.88.0
   cd $SCR/rb-$tag && $RB --version > $SCR/rb-version-$build.txt 2>&1
   $RB check-mkv.Rev > $ROOT/$CET/smoke-results/revbayes-check-mkv-$build-$tag.txt 2>&1)
done
