#!/bin/bash
#SBATCH --job-name=Job4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-9:0:0  #
#SBATCH --mem=8GB


#This is the fourth script to run. Here, we finally analyse all of our data.
module load r

Rscript /global/home/hpc4300/BIM_Final_RCodes/BIM_Final_step4_analyse_data.R $SLURM_ARRAY_TASK_ID
