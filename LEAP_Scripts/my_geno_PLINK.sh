#!/bin/bash
#SBATCH --job-name=LEAP_2000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:0:0
#SBATCH --mem=200GB


cd /global/scratch/hpc4298/LEAP_1_files

mkdir -p extracts

#plinkdir=../../../../../../../plink

plinkdir=../../../../global/home/hpc4298/plink

for i in {1..1000}; do

	${plinkdir} --noweb --file my_geno${i}  --make-bed --out my_geno${i}

	awk '{print $1, $2, $6}' my_geno${i}.fam > my_pheno${i}.phe
	
	awk '{ if ( $3 == 2 ) { $3 = 1 } else if ( $3 == 1 ) { $3 = 0 }; print}' my_pheno${i}.phe > Phenotype${i}.phe
	
	awk '{print $2}' my_geno${i}.bim | grep 'SNP' > extracts/ext${i}

	awk '{print $2}' my_geno${i}.bim | grep 'Csnp' > extracts/no_ext${i}
done	

