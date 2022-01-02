#!/bin/bash
#SBATCH --job-name=LEAP_Association1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-24:0:0
#SBATCH --mem=100GB


mkdir -p results3

for i in {1..1000}; do
############ Liability Estimation Part ###############

#Find related individuals to remove
#python findRelated.py --bfilesim EPS_geno${i} --out results3/EPS_geno${i}.related

#Compute eigendecompositions for each left-out chromosome
python eigenDecompose.py --bfilesim EPS_geno${i}  --out results3/EPS_geno${i}.npz

#Compute heritability estimates for every excluded chromosome
python calc_h2.py --bfilesim EPS_geno${i}  --prev 0.001 --numRemovePCs 10 --pheno pheno${i}.phe  --related results/EPS_geno${i}.related --h2coeff 1.0 --eigen results3/EPS_geno${i}.npz | tail -1 | cut -d\" \" -f2 > results3/EPS_geno_nochr${i}.h2


#python calc_h2.py --bfilesim EPS_geno1  --prev 0.001 --numRemovePCs 10 --pheno pheno1.phe  --eigen results3/EPS_geno1.npz  | tail -1 | cut -d ' ' -f 2  > results3/EPS_geno_chr1.h2

#python ../calc_h2.py --bfilesim dataset1/dataset1 --extractSim dataset1/extracts/nochr{}_extract.txt --prev 0.001 --numRemovePCs 10 --pheno dataset1/dataset1.phe --related results/dataset1.related --h2coeff 1.0 --eigen results/dataset1_nochr{}.npz | tail -1 | cut -d\" \" -f2 > results/dataset1_nochr{}.h2

# #Estimate liabilities
python probit.py --bfilesim EPS_geno${i} --pheno pheno${i}.phe --prev 0.001  --out results3/EPS_geno_nochr${i} --related results3/EPS_geno${i}.related --h2 \`cat results3/EPS_geno_nochr${i}.h2\` --hess 0  --eigen results3/EPS_geno${i}.npz


#python probit.py --bfilesim EPS_geno1 --pheno pheno1.phe --prev 0.001  --out results3/EPS_geno_nochr1  --h2  $(cat results3/EPS_geno_chr1.h2) --hess 0  --eigen results3/EPS_geno1.npz


#seq 1 10 | xargs -i --max-procs=2 bash -c "python ../probit.py --bfilesim dataset1/dataset1 --pheno dataset1/dataset1.phe --prev 0.001 --extractSim dataset1/extracts/nochr{}_extract.txt --out results/dataset1_nochr{} --related results/dataset1.related --h2 \`cat results/dataset1_nochr{}.h2\` --hess 0  --eigen results/dataset1_nochr{}.npz"


# GWAS for each pseudo-chromosome
python leap_gwas.py --bfilesim EPS_geno${i} --bfile geno${i} --pheno results3/EPS_geno_nochr${i}.liabs  --out results3/EPS_geno${i}.gwas.out.txt --h2 \`cat results3/EPS_geno_nochr${i}.h2\`  --eigen results3/EPS_geno${i}.npz

done

