#!/bin/bash
#SBATCH --job-name=Focal_Point_Job
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-9:0:0  
#SBATCH --mem=8GB


#This is the second script to run.
#We take our raw data in the BIM_Final_Raw_Data directory, clean it,
#and save our clean data in the BIM_Final_Clean_Data directory.
module load r

Rscript /global/project/hpcg1578/Crohn/Kernel_Analysis/Crohn_focal_point_analysis.R

echo "End of program at 'date'"