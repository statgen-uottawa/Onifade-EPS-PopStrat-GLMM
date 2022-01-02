#!/bin/bash
#SBATCH --job-name=Gemma_assoc1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:0:0
#SBATCH --mem=100GB

gemmadir=../../../gemma-0.98.1-linux-static

##Genetic relatedness matrix

for i in {1..1000}; do

	name=simu${i}       
	candfile=cand${i}.geno
	output2=assoc_test${i}

	${gemmadir} -g ${name}.geno.txt -p ${name}.pheno.txt -a ${name}.anno.txt -gk 1 -o ${name}

##fitting a linear mixed model

	${gemmadir} -g ${candfile}.txt -p ${name}.pheno.txt -n 1 -a ${name}.anno.txt  -k ./output/${name}.cXX.txt -lmm 3 -o ${output2}

done
