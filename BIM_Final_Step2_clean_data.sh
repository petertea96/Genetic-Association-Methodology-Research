#!/bin/bash
#SBATCH --job-name=test1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-9:0:0  
#SBATCH --mem=8GB
# commands for your job go here
#load-sse3

#This is the second script to run.
#We take our raw data in the BIM_Final_Raw_Data directory, clean it,
#and save our clean data in the BIM_Final_Clean_Data directory.
module load r

Rscript /global/home/hpc4300/BIM_Final_RCodes/BIM_Final_step2_prelim_codes.R

echo "End of program at 'date'"