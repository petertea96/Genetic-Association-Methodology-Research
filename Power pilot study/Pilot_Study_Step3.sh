#!/bin/bash
#SBATCH --job-name=PilotStudy3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-9:0:0  #
#SBATCH --mem=8GB


#This is the third script to run.
#We simulate phenotype data given the 
#genetic data we simulated from the previous script.
module load r

Rscript /global/home/hpc4300/Pilot_Study/Pilot_Study_RCode/Pilot_Study_Part3.R $SLURM_ARRAY_TASK_ID


echo "End of program at 'date'"