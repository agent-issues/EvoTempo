#!/bin/bash
# Production combine after bf-eval-40.sh (EVAL_JOB = its Slurm id): per-chain
# eval check table, then combine.  W and O are exported for the Rscript
# one-liner (#105 review, R105-12: before, O was empty and eval-checks.csv
# failed with 'cannot open file /eval-checks.csv').
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH -o /nobackup/pjjg18/evotempo-penguin-cet/logs/%x-%j.out
set +e; set +o pipefail; source /etc/profile || true; module purge || true
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export W=/nobackup/pjjg18/evotempo-penguin-cet
export R_LIBS=$W/lib-48e6bf7:$W/lib-tools:/nobackup/pjjg18/evotempo-dev/lib
cd $W/src-105
export O=${O:-$W/out/bf-400k}; mkdir -p $O
echo "src $(git rev-parse HEAD)" > $O/eval-checks.txt
sacct -j $EVAL_JOB --format=JobID%16,State,ExitCode,Elapsed -n -X >> $O/eval-checks.txt
Rscript -e 'f <- Sys.glob(file.path(Sys.getenv("W"), "runs", "*_run*", "cross.csv")); r <- do.call(rbind, lapply(f, function(x) { d <- read.csv(x); data.frame(chain = basename(dirname(x)), rows = nrow(d), maxAbsCheck = max(abs(d$check)), maxAbsCheckChar = if (is.null(d$checkChar)) NA else max(abs(d$checkChar))) })); write.csv(r, file.path(Sys.getenv("O"), "eval-checks.csv"), row.names = FALSE); cat(nrow(r), "cross.csv files; max |check|", max(r$maxAbsCheck), "; max |checkChar|", max(r$maxAbsCheckChar), "\n")' >> $O/eval-checks.txt 2>&1
/usr/bin/time -v Rscript data-raw/empirical/penguins/cet/bayes-factors.R combine $W/runs > $O/combine.txt 2>&1
echo "combine exit $?" >> $O/combine.txt
cp $W/runs/bayes-factors.csv $W/runs/closure.txt $W/runs/model-probabilities.csv $O/ 2>/dev/null
