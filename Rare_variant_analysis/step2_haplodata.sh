#!/bin/bash
#SBATCH --job-name=step_2_Haplodata
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-24:0:0
#SBATCH --mem=100GB

module load r

Rscript step2_EPS_rare_variants.R
