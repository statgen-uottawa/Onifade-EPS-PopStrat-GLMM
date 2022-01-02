#!/bin/bash
#SBATCH --job-name=PCarray_job_3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10-00:0:0
#SBATCH --mem=800GB

module load r

Rscript PC_test2.R $SLURM_ARRAY_TASK_ID

