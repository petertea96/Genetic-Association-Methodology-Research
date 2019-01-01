#!/bin/bash
#SBATCH --job-name=Job5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH -c 1                      # Number of CPUS requested (default is 1 CPU).
#SBATCH --mem=4                  # Memory requested in megabytes (default is 1024 MB).
#SBATCH -t 0-2:5:0      # How long you expect your job to run for (default is 3 hours).


#--> Since the previous job was an array job, we aggregate all files into one.

cd /global/home/hpc4300/BIM_Final_Results/

cat Pheno1Results_{1..100}.txt > /global/home/hpc4300/Aggregate_Pheno1Results.txt

cat Pheno2Results_{1..100}.txt > /global/home/hpc4300/Aggregate_Pheno2Results.txt

cat Pheno1Results_{1..100}_SLT.txt > /global/home/hpc4300/Aggregate_Pheno1Results_SLT.txt

cat Pheno2Results_{1..100}_SLT.txt > /global/home/hpc4300/Aggregate_Pheno2Results_SLT.txt