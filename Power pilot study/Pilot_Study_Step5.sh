#!/bin/bash
#SBATCH --job-name=PilotStudy5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ptea035@uottawa.ca
#SBATCH -c 1                      # Number of CPUS requested (default is 1 CPU).
#SBATCH --mem=4                  # Memory requested in megabytes (default is 1024 MB).
#SBATCH -t 0-2:5:0      # How long you expect your job to run for (default is 3 hours).

cd /global/home/hpc4300/Pilot_Study/Pilot_Study_Results/

cat Actual_common_causal_vector_{1..100}.txt > /global/home/hpc4300/Pilot_Study/Pilot_Study_Pretty_Results/Aggregate_Actual_common_causal_vector.txt

cat Actual_rare_causal_table_{1..100}.txt > /global/home/hpc4300/Pilot_Study/Pilot_Study_Pretty_Results/Aggregate_Actual_rare_causal_table.txt

cat Chosen_common_causal_{1..100}.txt > /global/home/hpc4300/Pilot_Study/Pilot_Study_Pretty_Results/Aggregate_Chosen_common_causal.txt

cat Chosen_rare_causal_{1..100}.txt > /global/home/hpc4300/Pilot_Study/Pilot_Study_Pretty_Results/Aggregate_Chosen_rare_causal.txt

cat Pheno1Results_{1..100}_SLT.txt > /global/home/hpc4300/Pilot_Study/Pilot_Study_Pretty_Results/Aggregate_Pheno1Results.txt

cat Pheno2Results_{1..100}_SLT.txt > /global/home/hpc4300/Pilot_Study/Pilot_Study_Pretty_Results/Aggregate_Pheno2Results.txt