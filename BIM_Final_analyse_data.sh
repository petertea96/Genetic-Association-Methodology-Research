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

Rscript /global/home/hpc4300/BIM_Final_RCodes/BIM_Final_step3_analyse_data.R
