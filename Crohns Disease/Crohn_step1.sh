#!/bin/bash
#SBATCH --job-name=Job1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1:0:0  
#SBATCH --mem=40GB


module load r
Rscript /global/project/hpcg1578/Crohn/Kernel_Analysis/Crohn_data_preprocessing.R