## complete program to test the SMMAT software before either splitting it or 
## making intermediate programs and running it about a couple times like we did
## for common variants:


# Start by reading in the simulated haplotypes and forming genotypes
## STEP 1: Read in the haplotype data. Must specify that it is type "character". 
#filename = paste("haplodata",j, ".txt", sep="")
#filename= "haplodata.txt"

#haplodat=read.table(filename, colClasses=c("character"))

library(withr)
library(stringr)
library(SKAT)

#setwd("~/Documents/my_R_files/Rare_variants_analysis")
filename<- "haplodata2.txt"
haplo_file<- read.table(filename, colClasses = c("character"))

##obtain segsites
FirstLine = readLines(filename)[1]
FirstLine=unlist(strsplit(FirstLine,split=""))
segsites=length(FirstLine)
newhaplodat=matrix(as.numeric(unlist(strsplit(haplo_file[,1],split=""))),ncol=segsites,byrow=T)

## obtain genotype for computing the allele frequencies by combining random rows of the haplotype file
half_data = (nrow(newhaplodat)/2)
center= (nrow(newhaplodat)/2) + 1
full_data= nrow(newhaplodat)
sorted_rows<- c(sample(1:half_data), sample(center:full_data)) 
#sorted_rows<- c(sample(0:5000), sample(5001:10000))
sorted_haplodat <- newhaplodat[sorted_rows,]

genodat=matrix(nrow=nrow(haplo_file)/2,ncol=ncol(newhaplodat))
#initializing step...

for (i in 1:nrow(genodat)){
  genodat[i,]=newhaplodat[2*i-1,]+newhaplodat[2*i,]
}

### simulate phenotypes and select extreme samples. 
#rows are individuals and columns are SNP sites

X<- genodat
pheno_indep <-c()
pheno1<- rnorm(n=nrow(X)/2, mean= 0.07, sd=1)
pheno2<- rnorm(length((nrow(X)/2):(nrow(X) -1)), mean= -0.07, sd=1)
pheno_indep<- c(pheno1,pheno2)  
geno_pheno= data.frame(pheno_indep,genodat)

##sort and subset the extremes

#labelling to keep track of the individuals in the extremes 
IND<- 1:nrow(geno_pheno)
combined_genodata<- cbind(IND, geno_pheno)

##sort in ascending order to subset for EPS data
sorted_combined <- combined_genodata[order(combined_genodata[,2]),] 

##Now Obtain the EPSdata 
#K = propnExtreme #proportion in the extremes
K= 0.1
Nums = nrow(sorted_combined)
keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
epsdat<- c(rep(1, K*Nums),rep(0, K*Nums))  ##cases are 1 , controls are 0: Extreme phenotype data
EPS_pheno <- as.matrix(cbind(epsdat,sorted_combined[keep,]))


#extract genotype file corresponding to the EPS individuals selected

EPS_geno<-data.frame(EPS_pheno[,-c(1:3)]) ##This is the genotype file we are using for the rare variants association analysis

##use SKAT function to calculate Minor allele frequencies (MAF) and weights(beta weights)

MAF = colSums(EPS_geno)/(2*nrow(EPS_geno)) 

weights.beta=c(1,25)  ##default weight parameters for SKAT weights

##function to compute the weights using the beta function as obtained  from the SKAT functions page

Beta.Weights<-function(MAF,weights.beta){
  
  n<-length(MAF)
  weights<-rep(0,n)	
  IDX_0<-which(MAF == 0)
  if(length(IDX_0) == n){
    stop("No polymorphic SNPs")
  } else if( length(IDX_0) == 0){
    weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
  } else {
    weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
  }
  
  
  #print(length(IDX_0))
  #print(weights[-IDX_0])
  return(weights)
  
}

w<- Beta.Weights(MAF, weights.beta) 


## Now edit files needed for the SMMAT program
#Now form group definition file
##group definition file
nSNPs= ncol(EPS_geno)

setid<- paste("set", rep(1:nSNPs), sep="")
chr<- rep(1:nSNPs)
pos<- seq(1:nSNPs)
referenceAllele<- sample(c("G","A"), nSNPs, replace=TRUE, prob=NULL)
variantAllele<- sample(c("G","A"), nSNPs, replace=TRUE, prob=NULL)
variantweights<- w

GroupFile <- cbind(setid, chr,pos,referenceAllele,variantAllele,variantweights)
write.table(GroupFile, file= "SetID.withweights.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)


 ## First files needed for the glmmkin model fitting
## geno file, pheno file, GRM is from GEMMA (input is the geno file)

## pheno file
id<- seq(1, nrow(EPS_pheno), 1)
disease<- EPS_pheno[,1]
traits<- EPS_pheno[,3]
pheno_file<- cbind(id, disease, traits)
write.table(pheno_file, file= "rare_pheno.file", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)


##genofile of all the EPS variants; use this as file for the GRM input as well:
##rename the row names
rownames(EPS_geno)<-  seq(1, nrow(EPS_geno),1)
colnames(EPS_geno)<-  seq(1, ncol(EPS_geno),1)
write.table(EPS_geno, file= "rare_geno.file", col.names =TRUE, row.names = TRUE, sep=" ", quote = FALSE)


## program for the genotype PLINK files of the rare variants
##Obtain PED and MAP files, then convert to bed, fam and bim 
## via PLINK commands and then use seqarray to convert to gds


## Read in genotype file: recode from 0,1,2 to 11, 12, 22
## to be able to use in the plink files

genotypes<- EPS_geno
phenotypes<- pheno_file

# for (i in (1:genotypes)){
#   if (genotypes[,i] == 0) genotypes[,i] = "1 1"
#   if (genotypes[,i] == 1) genotypes[,i] = "1 2"
#   if (genotypes[,i] == 2) genotypes[,i] = "2 2"
# }


genotypes[genotypes == 0] <- "1 1"
genotypes[genotypes == 1] <- "1 2"
genotypes[genotypes == 2] <- "2 2"
 
Recoded_geno<- genotypes



### making PED files: 6 mandatory columns, then MAP files
Family_ID<- rep(1, nrow(genotypes))
Individual_ID<- seq(1, nrow(genotypes), 1)
Paternal_ID<- rep(0, nrow(genotypes))
Maternal_ID<- rep(0, nrow(genotypes))
Sex<- sample(c(1,2), nrow(genotypes), replace=TRUE, prob=NULL)
Phenotype<- phenotypes[,2]

PED_file<- cbind(Family_ID, Individual_ID, Paternal_ID, Maternal_ID, Sex, Phenotype, genotypes)
write.table(PED_file, file="rare_eps.ped", col.names=FALSE, row.names=FALSE, sep= "\t", quote=FALSE)

###MAP file: 4 columns, mostly family information
##Physical position in cm

physical_pos = as.numeric(c(0:length(2:nSNPs) %o% 10^4))
for (i in physical_pos){
  gen_posit = physical_pos * 0.01 / 1000000
  with_options(c(scipen = 999.0), str_pad(gen_posit, 3, pad = "0"))
}

Chr<- rep(1, nSNPs)
SNP<- paste("snp", 1:nSNPs, sep="")
cm<- rep(1:nSNPs)
bp<- physical_pos

MAP_file <- cbind(Chr, SNP, cm, bp)
write.table(MAP_file, file="rare_eps.map", col.names=FALSE, row.names=FALSE, sep= "\t", quote=FALSE)

