#!/bin/bash
#SBATCH --job-name=Full_Gemmajob3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:0:0 # days-hours:minutes:seconds; How long you expect your job to run for (default is 3 hours).
#SBATCH --mem=500GB     #Memory requested in megabytes (default is 1024 MB).

# commands for your job go here

module load r

#echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
Rscript GEMMA2.R $SLURM_ARRAY_TASK_ID


gemmadir=../../../../gemma-0.98.1-linux-static

for i in {1..2000}; do
	    name=simu${i}
	    #../../../../gemma-0.98.1-linux-static -g ${name}.geno.txt -p ${name}.pheno.txt -a ${name}.anno.txt  -gk 1 -o ${name}
	    ${gemmadir} -g ${name}.geno.txt -p ${name}.pheno.txt -a ${name}.anno.txt  -gk 1 -o ${name}

done

##run R module to convert the matrices to positive semi definite

#Rscript norm_gemma.R $SLURM_ARRAY_TASK_ID
