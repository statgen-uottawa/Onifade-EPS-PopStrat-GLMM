### Step 4 of the Extreme Phenotype sampling fror rare variants
### Here, we use the functions in the R buRden program to compute Madsen and Browning association 
### statistics by permutation. 
## We first load the library(buRden): This was installed on the cluster and its use on Rstudio was not tested

library(buRden)

for (i in (1:1000)){
  
  setwd("/global/scratch/hpc4298/EPS_rareVariant_test/step_3_data")
  
  cc_data <- paste("MBdata", i, "perm", sep="")
  rec_cc_data <- readRDS(cc_data)
  keep <- filter_sites(rec_cc_data$ccgeno, rec_cc_data$ccpheno,0,0.05,0.8)
  rec_cc_data.MB.p <- MB.p.perm(cc_data$ccgeno[,which(keep==1)], cc_data$ccpheno,100)
  
  ## obtain pvalue of the general genetic model for all 1000 simulations
  
  Pvals= write(rec_cc_data.MB.p$p.value.general, file = paste("P_value_general", i, ".txt", sep=""), ncolumns= 1, sep = "\n")
  
}