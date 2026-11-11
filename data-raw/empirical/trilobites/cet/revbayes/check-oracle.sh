#!/bin/bash
# Trilobite CET pipeline: RevBayes oracle and check-mkv.R on Hamilton.
#
#   sbatch -J evotempo-trilobite-cet-oracle data-raw/empirical/trilobites/cet/revbayes/check-oracle.sh
#
# From a checkout at $SRC (default $W/src) with EvoTempo installed in $LIB and
# cet/data/ written by prepare.R.  Copies ../../revbayes/ to $W/oracle/, runs
# its prepare-data.R there (in a copy of the dataset folder) (the comparator's own matrices and taxa), writes the
# oracle inputs (check-mkv.R write), runs check-oracle.Rev with $RB, and
# compares (check-mkv.R compare), writing to $OUT (default
# $SRC/data-raw/empirical/trilobites/cet/smoke-results).
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH -o /nobackup/pjjg18/evotempo-trilobite-cet/logs/%x-%j.out
set +e; set +o pipefail
source /etc/profile || true
module purge || true
W=${W:-/nobackup/pjjg18/evotempo-trilobite-cet}
SRC=${SRC:-$W/src}
LIB=${LIB:-$W/lib-$(cat $W/out/build-sha.txt)}
RB=${RB:-/nobackup/pjjg18/revbayes-dev-fd807d10/projects/cmake/build/rb}
CET=data-raw/empirical/trilobites/cet
OUT=${OUT:-$SRC/$CET/smoke-results}
O=$W/oracle
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
module load r/4.5.1
export R_LIBS=$LIB:/nobackup/pjjg18/evotempo-dev/lib
cd $SRC
rm -rf $O; mkdir -p $O/trilobites
# prepare-data.R reads ../trilobites_*.{nex,csv}: keep the dataset's layout
cp -r data-raw/empirical/trilobites/revbayes data-raw/empirical/trilobites/trilobites_* $O/trilobites/
RBD=$O/trilobites/revbayes
(cd $RBD && Rscript prepare-data.R) > $O/prepare-data.txt 2>&1
Rscript $CET/check-mkv.R write $RBD/oracle > $O/write.txt 2>&1
cp $CET/revbayes/check-oracle.Rev $RBD/
(module purge; module load gcc/14.2 boost/1.88.0; cd $RBD && $RB --version > $O/rb-version.txt 2>&1;
 $RB check-oracle.Rev > $O/rb-oracle.txt 2>&1; echo "rb exit $?" >> $O/rb-oracle.txt)
mkdir -p $OUT
cp $RBD/oracle/revbayes-oracle.tsv $O/rb-version.txt $OUT/ 2>/dev/null
Rscript $CET/check-mkv.R compare $RBD/oracle $OUT > $OUT/check-mkv.txt 2>&1
echo "exit $?" >> $OUT/check-mkv.txt
