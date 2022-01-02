#!/bin/bash
#SBATCH --job-name=Gemma_powerasso2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:0:0
#SBATCH --mem=40GB

# Include starting and ending times
start=`date +%s`

#gemmadir= ../../../../../../global/home/hpc4298/gemma-0.98.1-linux-static 
gemmadir=../../../../gemma-0.98.1-linux-static
##Genetic relatedness matrix

for i in {1..1000}; do

	name=genomat${i}       
	candfile=geno_cand${i}.geno
	output2=Bin_assoc_test${i}

	${gemmadir} -g ${name}.geno.txt -p ${name}.pheno.txt -a ${name}.anno.txt -gk 1 -o ${name}

##fitting a linear mixed model

	${gemmadir} -g ${candfile}.txt -p ${name}.pheno.txt -n 1 -a ${name}.anno.txt  -k ./output/${name}.cXX.txt -lmm 3 -o ${output2}

done

end=`date +%s`

runtime=$((end-start))
