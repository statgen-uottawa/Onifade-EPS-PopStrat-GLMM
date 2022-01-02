## Author: Maryam Onifade
## Simulation to assess power of the methods:
## Our method to assess the power is to simulate a causal variant with allele frequency that is independent of ancestry. 
## We will be using a causal allele frequency of about 0.2 that will be the same for both populations. 
## Decide on the effect size for the causal allele ??? 

##Library running function for calculaing the effect size.
#library(pwr)

setwd("/global/scratch/hpc4298/GMMAT_PowerSimulations")

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files=seq(from=1, to=1001, by=10)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

for (i in (starting:ending)){
	numpop=2
	N=5000    #Total number of individuals
	nSNP=5000  #number of  snps
	Fst=0.01
	omega=c(0.5,0.5) ##MAF
	propnExtreme=0.1
	#nsim=100
	Fst.obs=vector(length=nSNP)
	pdiffs=vector(length=nSNP)
	genomat=matrix(nrow=N,ncol=nSNP)
	##Simulate snps for each population
	for (i in 1:nSNP){
 		p=runif(1,0.1,0.9) ##generating allele frequency p from uniform distribution
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
## X is the genotype matrix to be used for the GRM matrix
	X2=genomat 


## Power simulation starts here: simulate a causal locus: 

#set minor alelle frequencies in the two subpopulations to be p1=p2=0.2 since this scenario assumes there is 
##population stratificaation. 
##p1 is the causal allele frequency in population 1 and p2 is the causal allele frequency in  population 2

	p1=0.2; p2 = 0.2
	geno_SNP<- sample(c(0:2), N, prob=c(p1^2, 2*p1*p2, p2^2), replace=T)
#summary(as.data.frame(geno_SNP)) same as using table(geno_SNP). 
# The result gives us counts for each of the alelles

## Calculate power at alpha=0.05 by performing an ANOVA analysis of the continous phenotype and
## the three genotypes as the factor levels; This helps us to determine the effect size beta for the 
## phenotypes computation. 

#Effect_size= pwr.anova.test(k=3, n= 1666, f=NULL, sig.level = 0.05, power = 0.60)
	beta<- 0.0352
##An effect size of 0.0352 should give us a poer of around 0.6

##simulate phenotype values from normal distribution to agree with 
## The phenotype model should follow the model: y ~ N(mu_i + betaG_{ij}, 1 ) where 
## beta is the effect size per allele determined from the one way ANOVA analysis. 
## mu is the effect for subpopulation 1 while -mu is the effect for subpopulation 2

	pheno_indep <-c()
	pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07 + (beta * geno_SNP), sd=1))
	pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07 + (beta * geno_SNP), sd=1))
	pheno_indep<- c(pheno1,pheno2)

##Now obtain the EPS sample:
##individual ID's and 
	IND<- 1:N
	combined_Ind<- data.frame(cbind(IND, pheno_indep, X2))
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
	write.table(EPS_phenotypes, file =  paste("phenodata", i, ".txt", sep="") , col= T, row= FALSE, sep= "\t", quote=FALSE)
#write.table(EPS_phenotypes, file = "phenodata.txt", col.names = TRUE, row.names=FALSE, sep= "\t", quote =FALSE) 

#randomly simulate reference and effect alleles, we are frmatting the genotype file now!
	EPS_geno<- EPS_pheno[,-c(1:3)]
	EPS_2<- t(EPS_geno)
	Ref<- c(rep("A", nrow(EPS_2)))
	Eff<- sample(c("C","G","T"), nrow(EPS_2), replace=TRUE, prob= NULL)
	EP<- cbind(Ref,Eff, EPS_2)
	write.table(EP, file = "geno.txt", col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)

#Reformat the causal locus 
	cand<- c("SNP_Cand", "A","C", geno_SNP)

	#Converting strings into numerics
	my_df = data.frame(ncol=length(geno_SNP)+3, nrow=1)
	my_df[1,1] = cand[1]
	my_df[1,2] = cand[2]
	my_df[1,3] = cand[3]

	for (l in 1:length(geno_SNP)){
  	my_df[1,l+3] = as.numeric(cand[l+3])
	
	} 

	write.table(my_df, file = paste("genos", i, ".txt", sep=""), col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)
#write.table (my_df, file = "genos.txt", col.names = FALSE, row.names=FALSE, sep = "\t", quote=FALSE)
#}

}


