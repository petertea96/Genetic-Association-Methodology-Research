#!/bin/bash
#SBATCH --job-name=PilotStudy2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-9:0:0  
#SBATCH --mem=8GB
# commands for your job go here

#This is the second script to run.
#We take our raw data in the Pilot_Study_Raw_Data directory, clean it,
#and save our clean data in the Pilot_Study_Clean_Data directory.
module load r

Rscript /global/home/hpc4300/Pilot_Study/Pilot_Study_RCode/Pilot_Study_Part2.R $SLURM_ARRAY_TASK_ID

echo "End of program at 'date'"