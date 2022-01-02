#!/bin/bash
#SBATCH --job-name=LEAP_gwas2000
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

mkdir -p results

for i in {1..1000}; do

###Perform eigendecomposition while specifiying SNPs to be used for GRM compututation
python eigenDecompose.py --bfilesim my_geno${i} --extractSim extracts/ext${i}  --out results/my_geno_eigen${i}.npz

###Calculate heritability: this stage outputs heritability estimates for each file. 

python calc_h2.py --bfilesim my_geno${i}  --prev 0.001 --numRemovePCs 10 --pheno Phenotype${i}.phe  --eigen results/my_geno_eigen${i}.npz | tail -1 | cut -d ' ' -f 2  > results/my_geno_heri${i}.h2

###Compute liabilities for each individual: this stage uses heritability estimate and the eigen file computed from stages 1 and 2 above
python probit.py --bfilesim my_geno${i} --pheno Phenotype${i}.phe --prev 0.001  --out results/my_geno_li${i}  --h2  $(cat results/my_geno_heri${i}.h2) --hess 0  --eigen results/my_geno_eigen${i}.npz --extractSim extracts/ext${i}

##Test for associations.

python leap_gwas.py --bfilesim my_geno${i} --bfile my_geno${i} --pheno results/my_geno_li${i}.liabs  --out results/my_geno${i}.gwas.out.txt --h2 $(cat results/my_geno_heri${i}.h2)  --eigen results/my_geno_eigen${i}.npz --extractSim extracts/ext${i} --extract extracts/no_ext${i}


done
