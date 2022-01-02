#!/bin/bash
#SBATCH --job-name=LEAP_arrayjob2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-200
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5-00:0:0	# days-hours:minutes:seconds; How long you expect your job to run for (default is 3 hours).
#SBATCH --mem= 10GB	#Memory requested in megabytes (default is 1024 MB).

# commands for your job go here

module load r
#echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
Rscript LEAP_EPS.R $SLURM_ARRAY_TASK_ID



