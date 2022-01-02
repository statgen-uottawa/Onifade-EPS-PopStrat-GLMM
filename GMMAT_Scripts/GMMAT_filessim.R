
setwd("/global/scratch/hpc4298/GMMAT_Sup3")

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files=seq(from=1, to= 1001, by=10)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

for (i in (starting:ending)){

#for (i in 1:1000){
  #GM_sim <- function(numpop=2, N=2500, nSNP=5000, Fst=0.001, nsim=1) {
  numpop=2
  N=5000
  nSNP=5000
  Fst=0.001
  omega=c(0.5,0.5) 
  propnExtreme=0.1
  Fst.obs=vector(length=nSNP)
  pdiffs=vector(length=nSNP)
  genomat=matrix(nrow=N,ncol=nSNP)


  for (k in 1:nSNP){
    p=runif(1,0.1,0.9)
    alpha=p*(1-Fst)/Fst
    beta=(1-p)*(1-Fst)/Fst
    ps=rbeta(numpop,shape1=alpha,shape2=beta)
    vars=var(ps)
    pdiffs[k]=max(ps)-min(ps)
    Fst.obs[k]=vars/(p*(1-p))
  
    for (j in 1:numpop){
      ind1=(j-1)*N*omega[j]+1
      ind2=j*N*omega[j]
      freqs=c(ps[j]^2,2*ps[j]*(1-ps[j]),(1-ps[j])^2)
      genomat[ind1:ind2,]=sample(c(0,1,2),size=N*omega[j],replace=TRUE,prob=freqs)
      }

  }
  #X<- genomat[,i]
  X2<- genomat

##simulate phenotypes dependent on population since we are testing for type 1 error. 
  pheno_indep <-c()
  pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07, sd=1))
  pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07, sd=1))
  cont_pheno<- c(pheno1,pheno2)
##individual ID's and 
  IND<- 1:N
  combined_Ind<- data.frame(cbind(IND, cont_pheno, X2))
  sorted_contPheno<- combined_Ind[order(combined_Ind[,2]),] ##sort to subset for EPS data

#combined_indep <- data.frame(cbind(IND, pheno_indep, X))
#sorted_combined <- combined_indep[order(combined_indep[,2]),] 

##Obtain the EPSdata 
  K = propnExtreme #proportion in the extremes
  Nums = nrow(sorted_contPheno)
  keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
  epsdat<- c(rep("0",K*Nums),rep("1",K*Nums))
  EPS_pheno <- cbind(epsdat, sorted_contPheno[keep,])
  my_snp_names=paste("SNP", 1:ncol(X2), sep = "")
  colnames(EPS_pheno)<- c("Disease","IDs", "Trait", my_snp_names) ##Datafile of eps data and genotypes . Will be used for association studies. 
  EPS_phenotypes<- EPS_pheno[, c(1:3)]

  #setwd("/global/scratch/hpc4298/GMMAT_Pheno")
  #phenodata_name = paste("phenodata", i, ".txt", sep="")
  #phenodata_name = "phenodata"
  write.table(EPS_phenotypes, file =  paste("phenodata", i, ".txt", sep="") , col= T, row= FALSE, sep= "\t", quote=FALSE)
 #write.table(EPS_phenotypes, file = "phenodata.txt", col.names = TRUE, row.names=FALSE, sep= "\t", quote =FALSE) 
  
  #randomly simulate reference and effect alleles, we are frmatting the genotype file now!
  EPS_geno<- EPS_pheno[,-c(1:3)]
  EPS_2<- t(EPS_geno)
  Ref<- c(rep("A", nrow(EPS_2)))
  Eff<- sample(c("C","G","T"), nrow(EPS_2), replace=TRUE, prob= NULL)
  EP<- cbind(Ref,Eff, EPS_2)
  write.table(EP, file = "geno.txt", col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)

# One Causal Snp with MAF 0.25 and 0.85 for the two populations
  p1=0.7; p2 = 0.9
  geno_cand <- sample(c(0:2), nrow(EPS_pheno), prob=c(p1^2, 2*p1*p2, p2^2), replace=T)
  cand<- c("SNP_Cand", "A","C", geno_cand)

#Converting strings into numerics
  my_df = data.frame(ncol=length(geno_cand)+3, nrow=1)
  my_df[1,1] = cand[1]
  my_df[1,2] = cand[2]
  my_df[1,3] = cand[3]

  for (l in 1:length(geno_cand)){
    my_df[1,l+3] = as.numeric(cand[l+3])
  } 

  #setwd("/global/scratch/hpc4298/GMMAT_genos")
  #geno_name = paste("genos", i, ".txt", sep="")
  #geno_name = "genos"
  write.table(my_df, file = paste("genos", i, ".txt", sep=""), col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)
 #write.table (my_df, file = "genos.txt", col.names = FALSE, row.names=FALSE, sep = "\t", quote=FALSE)
  #}
  
}
 

