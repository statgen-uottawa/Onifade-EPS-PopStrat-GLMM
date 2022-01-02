##PC test 3

library("stringr")
library(withr)

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files = seq(from=1, to= 1001, by=10)

starting = total_files[task_id]
#compute starting index

ending = total_files[task_id +1] -1

for (k in (starting:ending)){

numpop=2
N=10000 ##population size 
nSNP=5000  #number of  snps used in GRM computation.
Fst=0.001
omega=c(0.5,0.5) ##MAF
propnExtreme=0.1
#nsim=1000

#for (k in (1:nsim)){
  
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
      genomat[ind1:ind2,i]=sample(c(0,1,2),size=N*omega[j],replace=TRUE,prob=freqs)
    }
    
  }
  
  X <- genomat
  
  ## perform standardization of genotype values and perform PC decomposition on the full dataset
  snpmeans=colMeans(X)
  pi=(1+colSums(X))/(2 + 2*nrow(X))
  Xmat=scale(X,center=snpmeans,scale=sqrt(pi*(1-pi)))
  pr=prcomp(Xmat,center=F,scale=F)
  #plot(pr$x,col=c(rep("red", N*omega[1]),rep("blue", N*omega[2])))
  
  ##Subset the first 5 PCS
  PC<- pr$x
  PC_five<- PC[,1:5]
  PC1<- PC[,1]
  ## Candidate SNPS 
  p1=0.5; p2 = 0.9
  geno_cand1<- sample(c(0:2), length(PC1), prob=c(p1^2, 2*p1*p2, p2^2), replace=T)
  ##simulate phenotypes dependent on population since we are testing for type 1 error. 
  pheno_indep <-c()
  pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07, sd=1))
  pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07, sd=1))
  pheno_indep<- c(pheno1,pheno2)
  ##individual ID's and Phenotype sorting
  IND<- 1:N
  #combined_indep <- data.frame(cbind(IND, pheno_indep, geno_cand1))
  combined_indep <- cbind(IND, pheno_indep, geno_cand1)
  sorted_combined <- combined_indep[order(combined_indep[,2]),] ##sort to subset for EPS data based on increasing pehnotype values
  ##Selecting the EPS rowsto keep
  ##10%in the upper and lower phenotypes distribution
  
  K = propnExtreme #proportion in the extremes  
  Nums = nrow(sorted_combined)
  keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
  epsdat<- c(rep("0",K*Nums),rep("1",K*Nums))
  EPS_pheno <- as.numeric(cbind(sorted_combined[keep,], epsdat)) ; dim(EPS_pheno)<- c(length(epsdat),4)
  colnames(EPS_pheno)<- c("ID", "Phenotype", "Genotype", "EPS")
  
  ##subset the principal components values for the EPS individuals
  PC_EPS_subset<- PC_five[sorted_combined[keep,1],]
  PC_EPS <- prcomp((genomat[EPS_pheno[,1],])) ##principal analysis of the EPS sample
  
  ##store the pvalues from the regression
  
  Fulldata<- c()
  EPS_Pvalues<- c()
  NO_correction<- c()
  
  ##Logistic regression 
  ##1. logistic regression using the PC covariates from the full dataset
  logit_full <- glm(EPS_pheno[,4]~EPS_pheno[,3]+PC_EPS_subset, family = "binomial"); Fulldata<- c(Fulldata, coef(summary(logit_full))[2,4])
  
  ##2.logistic regression with PC covariates from the EPS dataset
  
  logit_EPS <- glm(EPS_pheno[,4]~EPS_pheno[,3]+PC_EPS$x[,1:5], family = "binomial") ; EPS_Pvalues<- c(EPS_Pvalues, coef(summary(logit_EPS))[2,4])
  
  
  #3.logistic regresion with no PC correction
  logit_noPC <- glm(EPS_pheno[,4]~EPS_pheno[,3], family = "binomial"); NO_correction<- c(NO_correction, coef(summary(logit_noPC))[2,4])
  
  
  table<- round(cbind.data.frame(Fulldata, EPS_Pvalues, NO_correction),4)
 
  write.table(table, file=paste("PC_Values", k, ".csv", sep="" ), col.names = FALSE, row.names = FALSE, quote=FALSE, sep="," )
  
  
}
