## program for the genotype PLINK files of the rare variants
##Obtain PED and MAP files, then convert to bed, fam and bim 
## via PLINK commands and then use seqarray to convert to gds


## Read in genotype file: recode from 0,1,2 to 11, 12, 22
## to be able to use in the plink files

genotypes<- EPS_geno
phenotypes<- EPS_pheno

for (i in (1: ncol(genotypes))){
  if (genotypes[i] == 0) genotypes[i] = "1 1"
  if (genotypes[i] == 1) genotypes[i] = "1 2"
  if (genotypes[i] == 2) genotypes[i] = "2 2"
}

### making PED files: 6 mandatory columns, then MAP files
Family_ID<- rep(1, nrow(genotypes))
Individual_ID<- seq(1, nrow(genotypes), 1)
Paternal_ID<- rep(0, nrow(genotypes))
Maternal_ID<- rep(0, nrow(genotypes))
Sex<- sample(c(1,2), nrow(genotypes), replace=FALSE, prob=NULL)
Phenotype<- phenotypes[,2]

PED_file<- cbind(Family_ID, Individual_ID, Paternal_ID, Maternal_ID, Sex, Phenotype, genotypes)

