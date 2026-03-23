#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=48G
#SBATCH -p cpu
#SBATCH -c 5
#SBATCH --array=1-9
#SBATCH -o 6.features/logs/array_prom_%j.out

module load conda/latest
conda activate r

Rscript --no-restore --no-save dmaProm.R $SLURM_ARRAY_TASK_ID
