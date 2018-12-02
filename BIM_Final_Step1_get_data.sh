#!/bin/bash
#SBATCH --job-name=test1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1:0:0  
#SBATCH --mem=40GB
# commands for your job go here

#This is the first script to run.
#This script generates genetic datasets.

cd /global/project/hpcg1578/Peter
for i in `seq 1 2500`;
do 
./ms 200 1 -t 4.5 -T > /global/home/hpc4300/BIM_Final_Raw_Data/results$i.txt
done
 
