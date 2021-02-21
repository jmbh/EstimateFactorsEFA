#!/bin/bash
#SBATCH -N 1
#SBATCH -t 01:00:00

module purge
module load 2020
module load R/4.0.2-intel-2020a

cp -r "$HOME"/FactorSi4 "$TMPDIR"
cd "$TMPDIR"/FactorSi4

echo $SLURM_ARRAY_TASK_ID

Rscript --vanilla Simulation12.R $SLURM_ARRAY_TASK_ID

cp -r ./*.RDS "$HOME"/FactorSi4

