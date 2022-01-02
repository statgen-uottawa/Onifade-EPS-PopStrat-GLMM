#!/bin/bash
#SBATCH --job-name=GMMATfiles1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=monif064@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-50:0:0  #0 days 2 hours
#SBATCH --mem=200GB
# commands for your job follow

# variable giving file name
Rfile=GMMAT_filessim.R

# print the name of the variable Rfile
echo $Rfile

# print the working directory
pwd

echo 'Running R script'

# Call my R script. This will run the R script in the background. 
date
R --vanilla < $Rfile
date
