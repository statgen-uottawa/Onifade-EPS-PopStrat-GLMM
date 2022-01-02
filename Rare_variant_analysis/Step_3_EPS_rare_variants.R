#### This is step 3 in the EPS  rare variants analysis.
#### Here, we simulate the normal phenotype for the two populations we are considering and create genotypes
#### from the haplotypes by merging random rows. We then reformat the data and combining the simulated 
#### phenotype and the geneotype matrix. The last step here involves obtaining the EPS data by extracting 
#### the top and bottom 10th percentile

##Read in the genotype file from the ms program output::
##Genotypes are coded interms of the number of minor alleles. 


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)
#Obtain Slurm Task ID.  


#Now, determine indices of data files to analyse:
total_files=seq(from=1, to= 1001, by=100)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

for (j in (starting:ending)){ 
  
  setwd("/global/scratch/hpc4298/EPS_rareVariant_test/step1_haplo_data")
  
  ##Firt set the working directory to where the haplotypes data is saved
  
  ## STEP 1: Read in the haplotype data. Must specify that it is type "character". 
  filename = paste("haplodata",j, ".txt", sep="")
  
  haplodat=read.table(filename, colClasses=c("character"))
  
  ##obtain segsites
  FirstLine = readLines(filename)[1]
  FirstLine=unlist(strsplit(FirstLine,split=""))
  segsites=length(FirstLine)
  newhaplodat=matrix(as.numeric(unlist(strsplit(haplodat[,1],split=""))),ncol=segsites,byrow=T)
  
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
  
  
  ## simulate phenotypes and select extreme samples. 
  
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
  
  ##sort the phenotypes in ascending order to subset for EPS data
  sorted_combined <- combined_genodata[order(combined_genodata[,2]),] 
  
  ##Now Obtain the EPSdata 
  K= 0.1
  Nums = nrow(sorted_combined)
  keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
  epsdat<- c(rep(1, K*Nums),rep(0, K*Nums))  ##cases are 1 , controls are 0: Extreme phenotype data
  EPS_pheno <- as.matrix(cbind(epsdat,sorted_combined[keep,]))
  
  
  ### Compute the weights of the variants i.e the MAF of each variant
  f1=colSums(EPS_pheno)/nrow(EPS_pheno)
  
  
  
  
  
  ## Construct genotype and phenotype files for use in the buRden program: 
  ## this will be the genotype matrix corresponding to the extreme phenotypes.
  
  my_ccdata<- EPS_pheno[,-c(1:3)]
  
  my_genodata<- matrix(my_ccdata, nrow= nrow(EPS_pheno), ncol= ncol(EPS_pheno)-3, byrow = FALSE, dimnames = NULL)
  
  my_ccstatus<- as.numeric(EPS_pheno[,1])
  
  
  ##save these two data structures in a list as required by the program buRden
  rec_mydata <- list("ccgeno"= my_genodata, "ccpheno" = my_ccstatus)
  
  #saveRDS(rec_mydata, "MBdata")
  
  #set working directory to where we want to save the genotype and phenotype lists
  
  setwd("/global/scratch/hpc4298/EPS_rareVariant_test/step_3_data")
  
  saveRDS(rec_mydata, file=paste("MBdata", j, ".perm", sep=""))
  
}
  