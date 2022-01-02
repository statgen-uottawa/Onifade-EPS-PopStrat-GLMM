##Program for simulating files for running GEMMA programs

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

task_id <- as.numeric(slurm_arrayid)
##To obtain SLURmtaskid

##Determine the indices of the files to analyse

total_files= seq(from=1, to= 1001, by=10)

starting = total_files[task_id]
##compute starting index

ending = total_files[task_id  +1] -1
##compute ending index


for (k in (starting:ending)){

numpop=2
N=20000 ##number in each population 
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

##To simulate the putative causal locus
#set minor alelle frequencies in the two subpopulations to be p1=0.25 and p2=0.85
p1=0.25; p2 = 0.85
geno_SNP<- sample(c(0:2), N, prob=c(p1^2, 2*p1*p2, p2^2), replace=T)
##simulate phenotype values from normal distribution.
##simulate phenotypes dependent on population since we are testing for type 1 error. 
pheno_indep <-c()
pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07, sd=1))
pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07, sd=1))
pheno_indep<- c(pheno1,pheno2)
##individual ID's and 
IND<- 1:N
combined_indep <- cbind(IND, pheno_indep, X)
sorted_combined <- combined_indep[order(combined_indep[,2]),] ##sort to subset for EPS data
##Obtain the EPSdata 
K = propnExtreme #proportion in the extremes
Nums = nrow(sorted_combined)
keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
epsdat<- c(rep("1",K*Nums),rep("0",K*Nums))
EPS_pheno <- as.data.frame(cbind(epsdat,sorted_combined[keep,])) 
my_eps<-t(EPS_pheno[,-c(1:3)]) 
colnames(EPS_pheno)<- c("Bin_pheno", "Labels", "Phenotype", "Genotype")
###format GEMMA genotype file
##rs names for SNPS
n= nrow(my_eps); prefix<- "rs"; suffix<- seq(from=1, to= n , by = 1)
SNP_IDs = paste(prefix,suffix,sep="")
##alleles effect and reference 
Ref<- c(rep("A", nrow(my_eps)))
Eff<- sample(c("C","G","T"), nrow(my_eps), replace=TRUE, prob= NULL)
EP<- cbind(SNP_IDs, Ref,Eff, my_eps)
write.table(EP, file = paste("simu", k, ".geno.txt", sep=""),  col.names = FALSE, row.names = FALSE, sep= "\t", quote=FALSE)

###GEMMA phenotype file: A single line for the phenotypes: if binary, code cases 1 and controls 0
write.table(EPS_pheno[,1], file = paste("simu", k, ".pheno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep= "\n", quote=FALSE)

###GEMMA SNP annotation file 3 columns, SNP ID, basepair position, chromosome number
POS_val= seq(1000, 5000, length.out = nrow(my_eps)) ##explain what lines does
POS= format(POS_val,decimal.mark = ' ')  ##what does line do??
SNP_anno<- cbind(SNP_IDs, POS, CHR=rep(1, nrow(my_eps)))
write.table(SNP_anno, file=paste("simu", k, ".anno.txt", sep=""), col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)

}

