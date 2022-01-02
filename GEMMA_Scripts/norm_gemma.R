# upload kinship matrix from GEMMA

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

task_id <- as.numeric(slurm_arrayid)
##To obtain SLURmtaskid

##Determine the indices of the files to analyse

total_files= seq(from=1, to= 1001, by=10)

starting = total_files[task_id]
##compute starting index

ending = total_files[task_id  +1] -1
##compute ending index


for (i in (starting:ending)){

#for (i in 1:1000){
	  
	gemma_kins<- paste("simu", i, ".cXX.txt", sep="")
	rrk <- as.matrix(read.csv(gemma_kins, sep = "", row.names = NULL, header = FALSE))

	# R function to normalize the kinship matrix, making it positive semi definite
	normalize_kinmat <- function(kinmat){
		  #normalize kinship so that Kij \in [0,1]
		  tmp=kinmat - min(kinmat)
	  	  tmp=tmp/max(tmp)
	    	  tmp[1:9,1:9]
	    #fix eigenvalues to positive
	          diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
	          tmp[1:9,1:9]  
	          return(tmp)
	}

	normk <- normalize_kinmat(rrk) 
	str(normk)
	eigen_normk <- eigen(normk) # Check the eigenvalues to see whether the properties have improved
	eigen_normk$values

	# export matrix to use in GEMMA
	setwd("/global/scratch/hpc4298/GMMATFiles2")
	write.table(normk, file= paste("normalized_kin", i, ".cXX.txt", sep="" ),  sep = "\t", col.names = FALSE, row.names = FALSE)

}
