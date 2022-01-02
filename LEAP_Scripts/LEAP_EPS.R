## Simulating  files for LEAP:
## Using the same code used for simulating LTMLM with slight adjustments:

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files = seq(from=1, to= 3001, by=10)

starting = total_files[task_id]
#compute starting index

ending = total_files[task_id +1] -1

for (k in (starting:ending)){
 
  #LT_sim<- function (numpop=2,npop=500,nSNP=1000,Fst=0.01,nsim=100){
  numpop=2
  N=10000 ##number in each  population 
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
  X= genomat
  ##simulate phenotypes dependent on population since we are testing for type 1 error. 
  pheno_indep <-c()
  pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07, sd=1))
  pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07, sd=1))
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
  EPS_pheno<- cbind(epsdat,sorted_combined[keep,])
  
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
  
  ##create a candidate SNP and make corresponding PED and MAP files
  p1=0.25; p2 = 0.85
  cand_geno<- t(c("1", "1", "0", "0", "1", "1", sample(c("1 1", "1 2", "2 2"), nSNP, prob=c(p1^2, 2*p1*p2, p2^2), replace=T), sep= " "))
  
  write.table(cand_geno, file = paste("geno", k, ".ped", sep = ""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)
  
  
  #geno map file:
  
  MAP_geno <-cbind("2", paste("SNP", seq(1:nSNP), sep=""), "0", "0")
  write.table(MAP_geno, file = paste("geno", k, ".map", sep = ""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)

    }
  


  
  #genos2<- t(replicate(5, geno_cand1))
  #genos2<-t(genos2)
  #write.table(Data_EPS, file = paste("Sim1GRM", i, ".geno", sep=""), col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)
  #write.table(genos2, file = paste("First", k, ".geno", sep=""), col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)
  
  
  ##Create more putative SNPs: 
  # geno_cand1<- matrix(nrow=10, ncol=nrow(my_eps))
  # for (i in 1:nrow(geno_cand1)){
  #   p1=0.25; p2 = 0.85
  #   geno_cand1[i,]<- sample(c(0:2), nrow(my_eps), prob=c(p1^2, 2*p1*p2, p2^2), replace=T)
  # }
  # genos2<- t(geno_cand1)
  # write.table(genos2, file = paste("First", k, ".geno", sep=""), col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)
  # 
  # ##.ind file Thisfile has three columns. Individual ID's, sex of the individuals and phenotypes value:we are using our binary EPS phenotypes here
  # #creating ID labels
  # n= nrow(my_eps) ; prefix<- "SAMPLE"; suffix<- seq(from=0, to=n-1)
  # IDs = paste(prefix,suffix,sep="")
  # #sample randomly for the second column of gender for the sample numbers
  # col2 = sample(c("M", "F"), size= nrow(my_eps), replace=TRUE, prob= NULL)
  # indiv = data.frame(IDs, col2, col3 = as.factor(epsdat))
  # write.table(indiv, file= paste("First", k, ".ind", sep=""), col.names=FALSE, row.names = FALSE, sep=" ", quote = FALSE)
  # 
  # #.snp file: 6 columns, last two are optional, omitted here: 1st find a way to give snp names the suffix rs
  # # generate snps Id's in the form of "rsXXX"
  # # rgen=function(m=1111,to=(nSNP+1110)){
  # #   # m: this controls the number of digits
  # #   # if you want to start with n digits, then start m with (n-1) digits.
  # #   r=seq(m,to,1)
  # #   r=sprintf("rs%d",r)   #concatenate rs to numbers
  # #   return(r)
  # # }
  # 
  # p1<- "rs"; p2<- seq(100, (nSNP+99), 1)
  # SNPIds<- paste(p1,p2,sep="")
  # 
  # 
  # ##to simulate genetic and physical distances for the snp file
  # #explain what this block of code does
  # physical_pos = as.numeric(c(0:length(2:nSNP) %o% 10^4))
  # for (i in physical_pos){
  #   gen_pos1 = physical_pos * 0.01 / 1000000
  #   with_options(c(scipen = 999), str_pad(gen_pos1, 3, pad = "0"))
  # }
  # 
  # snpfile= data.frame(SNPIds, chromosome_number= rep(11, nSNP), gen_pos = rep(0, nSNP) , physical_pos)
  # 
  # write.table(snpfile, file = paste("FirstGRM", k, ".snp", sep=""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)
  # 
  # ##.snp file for 5 causal snps
  # # rgen2=function(m=1000,to=(4+1000)){
  # #   # m: this controls the number of digits
  # #   # if you want to start with n digits, then start m with (n-1) digits.
  # #   r=seq(m,to,1)
  # #   r=sprintf("rs%d",r)   #concatenate rs to numbers
  # #   return(r)
  # # }
  # 
  # rgen2<- paste(p1, seq(1, 5,1), sep="")
  # 
  # 
  # 
  # ##to simulate genetic and physical distances for the snp file
  # #explain what this block of code does
  # physical_pos2 = as.numeric(c(0:length(2:5) %o% 10^4))
  # for (i in physical_pos){
  #   gen_posit = physical_pos * 0.01 / 1000000
  #   with_options(c(scipen = 999.0), str_pad(gen_posit, 3, pad = "0"))
  # }
  # 
  # snpfile2= data.frame(SNPName1= rgen2(), chromosome_number= rep(11, 5), gen_pos2 =  rep(0, 5), physical_pos2)
  # #snps= snpfile2[-6]
  # write.table(snpfile2, file = paste("First", k, ".snp", sep=""), col.names = FALSE, row.names = FALSE, sep= " ", quote=FALSE)
  #
  



