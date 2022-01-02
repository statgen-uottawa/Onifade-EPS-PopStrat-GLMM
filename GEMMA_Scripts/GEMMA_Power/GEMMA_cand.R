########### Program for simulating files for running GEMMA programs  #################
##########  In this version, we include the candidate SNP into the   ###############
########### whole genome and just focus on that when obtaining the result.#################

### Set directory to the where the files will be saved
### Why not set it to the scratch folder??

#setwd("/global/home/hpc4298/GEMMA-0.98/myexample/GEMMA_quan_3")

#GEMMA_sim<- function(N, nSNP){

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files=seq(from=1, to=1001, by=10)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

for (k in (starting:ending)){

numpop=2
N=10000    #Total number of individuals
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

X= genomat


####To simulate the putative causal locus: set minor alelle frequencies in the two subpopulations to be p1=0.25 and p2=0.85
p1=0.2; p2 = 0.2
geno_SNP<- sample(c(0:2), nSNP, prob=c(p1^2, 2*p1*p2, p2^2), replace=T)

##Effect Size chosen so that power is around 0.6
Effect_Size<- 0.0352
beta<-Effect_Size
##simulate phenotype values from normal distribution following a specified model:check README file for a desription of the model. The models includes also an effect for subpopulation as well as the causal allele effect size.  
pheno_indep <-c()
pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07+(beta*geno_SNP), sd=1))
pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean=-0.07+(beta*geno_SNP),  sd=1))
pheno_indep<- c(pheno1,pheno2)
##individual ID's and 
IND<- 1:N
combined_indep <- cbind(IND, pheno_indep, X)
sorted_combined <- combined_indep[order(combined_indep[,2]),] ##sort to subset for EPS data
##Obtain the EPSdata 
K = propnExtreme #proportion in the extremes
Nums = nrow(sorted_combined)
keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
epsdat<- c(rep(1,K*Nums),rep(0,K*Nums))
EPS_pheno <- as.matrix(cbind(epsdat,sorted_combined[keep,])) 

##isolate the genotype matrix corresponding to the EPS sampples and add the candidate SNP as the first SNP
##Not anymore: the genotype matrix here doesnt include the candidate SNP
EPS_geno<- EPS_pheno[,-c(1:3)]  
#EPS_geno_cand <- rbind(geno_SNP, EPS_geno)
my_eps<- t(EPS_geno)

##Naming the columns.
pre<- "IND"; suf<- seq(1:ncol(my_eps))
colnames(my_eps)<- paste(pre,suf, sep="")

##rs names for SNPS
## name the first SNP in a different form from the remaining since its the candidate gene 
prefix<- "rs"; suffix<- seq(from=1, to=nrow(my_eps) , by = 1)
SNP_IDs<- paste(prefix,suffix,sep="")
##alleles effect and reference 
Ref<- c(rep("A", nrow(my_eps)))
Eff<- sample(c("C","G","T"), nrow(my_eps), replace=TRUE, prob= NULL)
EP<- cbind(SNP_IDs, Ref,Eff, my_eps)

write.table(EP, file = paste("genomat", k, ".geno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\t", quote=FALSE)

Pheno<- round(EPS_pheno[,3], 4)
#Pheno<- EPS_pheno[,1]
###GEMMA phenotype file: A single line for the phenotypes: if binary, code cases 1 and controls 0
#write.table(Pheno, file = "simu.pheno.txt", col.names = FALSE, row.names = FALSE, sep= "\n", quote=FALSE)
write.table(Pheno, file = paste("genomat", k, ".pheno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\n", quote=FALSE)

###GEMMA SNP annotation file 3 columns, SNP ID, basepair position, chromosome number
POS_val= seq(1000, 5000, length.out = nrow(EP)) ##explain what lines does
POS= format(POS_val,decimal.mark = ' ')  ##what does line do??
SNP_anno<- cbind(SNP_IDs, POS, CHR=rep(1, nrow(EP)))
#write.table(SNP_anno, file="simu.anno.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(SNP_anno, file = paste("genomat", k, ".anno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\t", quote=FALSE)


##single candidate snp

cand_ID<- c("rs243556"); r<- c("A"); E<- c("C")
cand_gene<- c(cand_ID, r, E, geno_SNP)
#cg<- as.matrix(cand_gene, nrow=1, ncol=length(geno_SNP)+3)
write.table(t(cand_gene), file = paste("geno_cand", k, ".geno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= ",", quote=FALSE)

#write.table(t(cand_gene), file="cand.geno.txt", col.names = FALSE, row.names=FALSE, sep=",", quote=FALSE)



}
