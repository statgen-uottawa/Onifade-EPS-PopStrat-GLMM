#!/bin/bash
#SBATCH --job-name=GMMAT_sup_07_09
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5-00:0:0	# days-hours:minutes:seconds; How long you expect your job to run for (default is 3 hours).
#SBATCH --mem=400GB	#Memory requested in megabytes (default is 1024 MB).

# commands for your job go here

module load r
#echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
Rscript GMMAT_filessim.R $SLURM_ARRAY_TASK_ID

