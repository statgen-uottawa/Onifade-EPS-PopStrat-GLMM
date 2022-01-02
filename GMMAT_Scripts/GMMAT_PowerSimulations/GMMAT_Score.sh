#!/bin/bash
#SBATCH --job-name=GMMAT_Power1Score
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-24:00:0  
#SBATCH --mem=20GB


#This is the first script to run.
#This script generates phenotype and genetic datasets.

module load r

Rscript  GMMAT_R.R


