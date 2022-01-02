####Step 2 of the EPS rare variants analysis########

##The input to this program are the output from the ms program output which have been stored 
## in the same format.
##This script will take each of the .txt file saved from the ms program in each of the 1000 iterations and 
##extracts only the haplotype data which is our interest here:

for (i in 1:1000){
  setwd("/global/scratch/hpc4298/EPS_rareVariant_test/raw_data")
  #--> Set working directory to where we originally saved the ms program output.
  
  title = paste("results",i,".txt", sep="")
  data=readLines(title)
  #All ms outputs saved  in the same name format. i.e: results1.txt
 
  haplodata = data[8:length(data)]
  #--> Haplotype data is saved, beginning on the 8th line. There should be 5000 lines
  # in total, representing 5000 different haplotypes simulated.
  
  filename1 = paste("haplodata",i, ".txt", sep="")
  
  setwd("/global/scratch/hpc4298/EPS_rareVariant_test/step1_haplo_data")
  #--> Set working directory to where we want to save the extracted data.
  
  writeLines(haplodata, filename1)
  
}