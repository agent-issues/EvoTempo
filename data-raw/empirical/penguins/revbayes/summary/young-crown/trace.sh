#!/bin/bash
#SBATCH -J evotempo-verify-99-trace
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o %x-%j.out
set +e; set +o pipefail
cd /nobackup/pjjg18/evotempo-verify-99/trace
python3 trace99.py
