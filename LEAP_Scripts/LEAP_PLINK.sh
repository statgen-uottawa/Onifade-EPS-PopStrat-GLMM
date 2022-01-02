#!/bin/bash
#SBATCH --job-name=LEAP_BMCJob2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:0:0
#SBATCH --mem=10GB


plinkdir=../../../../../../../plink


for i in {1..2000}; do
:q

	${plinkdir} --noweb --file EPS_geno${i} --1  --make-bed  --out EPS_geno${i}

	awk '{print $1,$2,$6}' EPS_geno${i}.fam > EPS_pheno${i}.phe
	
	${plinkdir}  --noweb --file geno${i} --no-pheno --make-bed --out geno${i}	

done
