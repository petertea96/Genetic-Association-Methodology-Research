#!/bin/bash
#SBATCH --job-name=Crohn3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-100:0:0  #
#SBATCH --mem=16GB



module load r

Rscript /global/project/hpcg1578/Crohn/Kernel_Analysis/Crohns_disease-analysis_part2.R $SLURM_ARRAY_TASK_ID
