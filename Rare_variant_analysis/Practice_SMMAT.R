## complete program to test the SMMAT software before either splitting it or 
## making intermediate programs and running it about a couple times like we did
## for common variants:


# Start by reading in the simulated haplotypes and forming genotypes
## STEP 1: Read in the haplotype data. Must specify that it is type "character". 
#filename = paste("haplodata",j, ".txt", sep="")
#filename= "haplodata.txt"

#haplodat=read.table(filename, colClasses=c("character"))


filename<- "edited_ms"
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

EPS_geno<-data.frame(EPS_pheno[,-c(1:3)])

##calculate Minor allele frequencies

Minor_alele_frequencies = 

  
  
  get_rare_causals = function(UpperBound=0.05, LowerBound=0.01, EPS_geno){
    #This function randomly chooses 10 SNP site with proper MAF to be the rare causal sites.
    #By default, my decision criteria for the rare causal sites:
    #Rare causal MAF: 0.01 < MAF < 0.05; Choose the same # of rare causals for ALL simulations.
    #select all rare causal variants using this function.
    
    Minor_Allele_Frequencies = colSums(EPS_geno)/(2*nrow(EPS_geno))
    
    possible_common_rare_sites = which( (Minor_Allele_Frequencies < UpperBound) & (Minor_Allele_Frequencies > LowerBound ) )
    #Create vector of possible SNP sites to be the rare causal variants.
    
    
    # if (length(possible_common_rare_sites)< numsites){
    #  return(NULL)
    #}
    if (length(possible_common_rare_sites)< 0){
      return(NULL)
    }
    
    #rare_sites = sample(possible_common_rare_sites, numsites)
    
    return(possible_common_rare_sites)
  }


## extract genomatrix and vector of binary phenotype
## to compute weights from buRden program
##weights of variants are the MAFs. 


## Now edit files needed for the SMMAT program
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

genotypes<- rare_geno.file
phenotypes<- rare_pheno_file

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
Sex<- sample(c(1,2), nrow(genotypes), replace=TRUE, prob=NULL)
Phenotype<- phenotypes[,2]

PED_file<- cbind(Family_ID, Individual_ID, Paternal_ID, Maternal_ID, Sex, Phenotype, genotypes)
write.table(PED_file, file="rare_eps.ped", col.names=FALSE, row.names=FALSE, sep= "\t", quote=FALSE)

###MAP file: 4 columns, mostly family information
##Physical position in cm
physical_pos = as.numeric(c(1:length(2:5) %o% 10^4))
for (i in physical_pos){
  gen_posit = physical_pos * 0.01 / 1000000
  with_options(c(scipen = 999.0), str_pad(gen_posit, 3, pad = "0"))
}

Chr<- rep(1, ncol(genotypes))
SNP<- paste("snp", 1:ncol(genotypes), sep="")
cm<- rep(0, ncol(genotypes))
bp<- physical_pos

MAP_file <- cbind(Chr, SNP, cm, bp)
write.table(MAP_file, file="rare_eps.map", col.names=FALSE, row.names=FALSE, sep= "\t", quote=FALSE)

View(MAP_file)



###perform CMC test using code obtained from thea bjornland:
collapse=function(genotypes,mafs){
  uncollapsed=genotypes[,mafs>=0.05]
  collapsedgens=rowSums(genotypes[,mafs<0.05])
  collapsedgens[collapsedgens>0]=1
  return(cbind(uncollapsed,collapsedgens))
}

cmc=function(phenotypes,covariates,genotypes,mafs){
  newgen=collapse(genotypes,mafs)
  pval=scorenormal(phenotypes,covariates,as.matrix(newgen))
  return(pval)
}

####


genotypes=EPS_pheno[,-c(1:3)]
mafs= colSums(genotypes)/(2*nrow(genotypes))
uncollapsed=genotypes[,mafs>=0.05]  ### This code filters the genotypes having a mafs >=0.05
collapsedgens=rowSums(genotypes[,mafs<0.05])  ## take the sum of the rows having genotypes with mafs >=0.05
result<-cbind(uncollapsed,collapsedgens)

###to contruct pvalues from a score function using CMC. The score function is GMMAT
phenotypes<- as.numeric(EPS_pheno[,1])
