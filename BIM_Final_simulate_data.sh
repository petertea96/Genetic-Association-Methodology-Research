#!/bin/bash
#SBATCH --job-name=test1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-9:0:0  #0 days 2 hours
#SBATCH --mem=8GB
# commands for your job go here
#load-sse3

#This is the second script to run.
#First we reformat the data. Then, we simulate phenotype data given the 
#genetic data we simulated from the previous script.
module load r

Rscript /global/home/hpc4300/BIM_Final_Rcodes/BIM_Final_step1_prelim_codes.R


Rscript /global/home/hpc4300/BIM_Final_Rcodes/BIM_Final_step2_subsetdata_Prepare_Dataset.R


echo "End of program at 'date'"