
python  findRelated.py --bfilesim leap_EPS/EPS_data  --out leap_EPS/results
python eigenDecompose.py --bfilesim leap_EPS/EPS_data --out leap_EPS/eigenfile
python calc_h2.py --bfilesim leap_EPS/EPS_data --prev 0.001 --numRemovePCs 10 --pheno leap_EPS/real_pheno.txt --related leap_EPS/results  --eigen leap_EPS/eigenfile.npz

python probit.py --bfilesim leap_EPS/EPS_data --pheno leap_EPS/real_pheno.txt --prev 0.001 --out leap_EPS/EPS_liab --h2 0.443679  --eigen leap_EPS/eigenfile.npz  --related leap_EPS/results

python leap_gwas.py --bfilesim leap_EPS/EPS_data  --bfile leap_EPS/EPS_data --pheno leap_EPS/EPS_liab.liabs --out leap_EPS/GWAS_EPS --h2 0.443679 --eigen leap_EPS/eigenfile.npz
