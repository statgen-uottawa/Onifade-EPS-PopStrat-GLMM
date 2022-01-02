## R script to run GMMAT on the cluster using the output files from GMMAT_filessim.R in directory ()
##RUNNING GMMAT using output files from GM_sim()

#setwd("/global/home/hpc4298/R/x86_64-pc-linux-gnu-library/3.5/GMMAT/mydata/Testdata")
setwd("/global/scratch/hpc4298/GMMAT_Sup3")

#library(GMMAT)

#pheno <- read.table("phenodata.txt", header = TRUE)
#GRM <- matrix(scan("sim1.cXX.txt", n = 1000^2, sep=""), nrow = 1000, ncol = 1000)
#geno.file<- read.table("genos.txt", header = FALSE)
#model0 <- glmmkin(Disease ~ 1, data = pheno, kins = GRM,  family = binomial(link = "logit"))
#glmm.score(model0, infile = "genos.txt", outfile = "glmm.score.text.testout22file", infile.ncol.skip = 3, infile.ncol.print = 1:3, infile.header.print = c("SNP", "Allele1", "Allele2"))


library(GMMAT)

#GRM <- matrix(scan("sim1.cXX.txt", n = 1000^2, sep=""), nrow = 1000, ncol = 1000)

for (i in 1:1000){

	phenodata_name= paste("phenodata", i, ".txt", sep="")
	pheno <- read.table(phenodata_name, header = TRUE)
	geno_name= paste ("genos", i, ".txt", sep="")
	cov_matrices<- paste("normalized_kin", i, ".cXX.txt", sep="")
	GRM<- matrix(scan(cov_matrices, n=1000^2, sep=""), nrow=1000, ncol=1000)
	#geno.file<- read.table(geno_name, header = FALSE)
	model0 <- glmmkin(Disease ~ 1, data = pheno, kins = GRM,  family = binomial(link = "logit"))
	glmm.score(model0, infile = geno_name, outfile = paste("glmm.sup_07_09.score",i, ".txt", sep=""), infile.ncol.skip = 3)
	
	}









#library(GMMAT)
#pheno <- read.table(phenodata_name, header = TRUE)
#GRM <- matrix(scan("sim1.cXX.txt", n = 1000^2, sep=""), nrow = 1000, ncol = 1000)
#geno.file<- read.table("geno.txt", header = FALSE)
#model0 <- glmmkin(Disease ~ 1, data = pheno, kins = GRM,  family = binomial(link = "logit"))
#glmm.score(model0, infile = geno_name, outfile = paste("glmm.score.text.testoutfile",i,".txt", sep=""), infile.ncol.skip = 3, infile.ncol.print = 1:3, infile.header.print = c("SNP", "Allele1", "Allele2"))
