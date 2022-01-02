## Simulating genotype and phenotype  files for LEAP power simulation:

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

#total_files = seq(from=1, to= 1001, by=10)
#changing this line of code to address the increase in number of simulations to 2000 as requested in the BMC manuscript review 
total_files = seq(from=1, to= 1001, by=10)

starting = total_files[task_id]
#compute starting index

ending = total_files[task_id +1] -1

for (k in (starting:ending)){
  
  #LT_sim<- function (numpop=2,npop=500,nSNP=1000,Fst=0.01,nsim=100){
  numpop=2
  N=5000 ##number in each  population 
  nSNP=5000  #number of  snps used in GRM computation.
  Fst=0.01
  omega=c(0.5,0.5) ##MAF
  propnExtreme=0.1
  nsim=100
  Fst.obs=vector(length=nSNP)
  pdiffs=vector(length=nSNP)
  genomat=matrix(nrow=N,ncol=nSNP)
  ##Simulate GRM snps
  for (i in 1:nSNP){
    p=runif(1,0.1,0.9)
    alpha=p*(1-Fst)/Fst
    beta=(1-p)*(1-Fst)/Fst
    ps=rbeta(numpop,shape1=alpha,shape2=beta)
    
    for (j in 1:numpop){
      ind1=(j-1)*N*omega[j]+1 ##index values for individuals in population 1
      ind2=j*N*omega[j]  #individuals in population 2
      freqs=c(ps[j]^2,2*ps[j]*(1-ps[j]),(1-ps[j])^2)
      genomat[ind1:ind2,i]=sample(c("1 1","1 2","2 2"),size=N*omega[j],replace=TRUE,prob=freqs)
    }
  }
  X= genomat  ### Ancestry SNPS for generating the kinhip matrix
  
  ## Simulate a causal SNP
  
  ##create a candidate SNP and make corresponding PED and MAP files
  p1=0.2; p2 = 0.2
  geno_cand <- sample(c(0:2), nSNP, prob=c(p1^2, 2*p1*p2, p2^2), replace=T)
  ## Causal SNP should be created using an additive model. LEAP files are formatted with some preceeding 
  ## columns followed by the genotypes. 
  ## This version worked for the initial simulations where we are making sure to follow the LEAP manual in terms
  ## of the files generation. 
  
  ## For this program to investigate power of the LEAP method, we need to first create an additively coded 
  ## causal SNP to be used in the phenotypes generation then convert from additive to this form of coding. 
  
  #Reformatting the causal SNP into the form required by the LEAP Program.
  New_geno<- c("2 2", "1 2", "1 1", geno_cand)[match(geno_cand, c(0, 1 ,2, geno_cand))]
  cand_geno<- t(c("0", "1", "0", "0", "1", New_geno))
  #cand_geno<- t(c("0", "1", "0", "0", "1", sample(c("1 1", "1 2", "2 2"), nSNP, prob=c(p1^2, 2*p1*p2, p2^2), replace=T), sep= " "))
  write.table(cand_geno, file = paste("geno", k, ".ped", sep = ""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)
  
  
  ##simulate phenotypes dependent on causal locus 
  beta<- 0.0352 ## Effect Size to guarantee a power of around 0.6
  pheno_indep <-c()
  pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07 + (beta * geno_cand), sd=1))
  pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07 + (beta * geno_cand), sd=1))
  pheno_indep<- c(pheno1,pheno2)
  ##individual ID's and 
  IND<- 1:N
  
  ##Full data of ID's phenotypes and genotype.
  combined_indep <- data.frame(cbind(IND, pheno_indep, X))
  
  #sort by quantitative phenotype values to obtain the EPS data.  
  sorted_combined <- data.frame(combined_indep[order(combined_indep[,2]),])
  
  ##Obtain the EPSdata 
  K = propnExtreme #proportion to be selected from the extreme and from a random proportion of the remaining population
  
  Nums = nrow(sorted_combined)
  ##obtain just the upper extremes known as the cases:
  #K1 <- 1:(K*Nums)
  #eps1<- (sorted_combined[K1,])
  #remain<- sorted_combined[setdiff(rownames(sorted_combined),rownames(eps1)),] ## substract the cases from the bigger dataset
  
  ##Sample controls from remaining
  ##controls
  #eps2<- remain[sample(nrow(remain), (K*Nums)), ]
  
  ##EPS data: combine dataframes eps1 and eps2
  #EPS<- rbind(eps1, eps2)
  
  keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
  epsdat<- c(rep("1",K*Nums),rep("0",K*Nums))
  EPS_pheno<- cbind(epsdat,sorted_combined)
  
  ##Just the EPS genotypes 
  my_eps<- EPS_pheno[,-c(1:3)]
  
  ##make PED and MAP file for the GRM: this excludes the SNP to be tested. 
  
  FID <-seq(1, nrow(EPS_pheno), 1)
  IID <-seq(1, nrow(EPS_pheno), 1)
  PID <- rep(0, nrow(EPS_pheno))
  MID<- rep(0, nrow(EPS_pheno))
  sex<- sample(c("1", "2"), nrow(EPS_pheno), prob = NULL, replace = T)
  Phenotype<- epsdat
  
  Ped_file<- cbind(FID, IID, PID, MID,sex,Phenotype, my_eps)
  
  write.table(Ped_file, file = paste("EPS_geno", k, ".ped", sep = ""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)
  
  ##corresponding MAP file:4 columns
  
  chr<- rep(1, nSNP)
  SNP_identifier<- paste("SNP", seq(1:nSNP), sep="")
  genetic_distance<- rep(0, nSNP)
  physical_distance<- seq(1, nSNP, 1)
  
  MAP_file<- cbind(chr, SNP_identifier, genetic_distance, physical_distance)
  
  write.table(MAP_file, file = paste("EPS_geno", k, ".map", sep = ""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)
  
  
  
  #geno map file:
  
  MAP_geno <-cbind("2", paste("SNP", seq(1:nSNP), sep=""), "0", "0")
  write.table(MAP_geno, file = paste("geno", k, ".map", sep = ""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)
  
}




