#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=24G
#SBATCH -p cpu
#SBATCH -c 6
#SBATCH --array=1-9
#SBATCH -o 6.fst/logs/array_%j.out

module load conda/latest
conda activate r-genetics

Rscript --no-restore --no-save locFst.R $SLURM_ARRAY_TASK_ID
