#!/bin/bash
#SBATCH --job-name=Gemma_Matrix_job1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-24:0:0
#SBATCH --mem=400GB


##run multiple gemma line for generating grm matrix for geno, pheno and snp annotation files times

gemmadir=../../../../gemma-0.98.1-linux-static

for i in {1..10}; do
    name=simu${i}
#../../../../gemma-0.98.1-linux-static -g ${name}.geno.txt -p ${name}.pheno.txt -a ${name}.anno.txt  -gk 1 -o ${name}
${gemmadir} -g ${name}.geno.txt -p ${name}.pheno.txt -a ${name}.anno.txt  -gk 2 -o ${name}

done
