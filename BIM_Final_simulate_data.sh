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
module load r

cd /global/home/hpc4300/BIM_Final_Rcodes



R CMD BATCH BIM_Final_step1_prelim_codes.R

R CMD BATCH BIM_Final_step2_subsetdata_Prepare_Dataset.R


echo "End of program at 'date'"