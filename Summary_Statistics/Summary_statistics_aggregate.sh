#!/bin/bash
#SBATCH --job-name=Summary_job_aggregate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH -c 1                      # Number of CPUS requested (default is 1 CPU).
#SBATCH --mem=4                  # Memory requested in megabytes (default is 1024 MB).
#SBATCH -t 0-2:5:0      # How long you expect your job to run for (default is 3 hours).


#--> Since the previous job was an array job, we aggregate all files into one.

cd /global/home/hpc4300/BIM_Final_PhenoAndGeno_Data

cat SummaryResults_{1..100}.txt > /global/home/hpc4300/Aggregate_Summary.txt




