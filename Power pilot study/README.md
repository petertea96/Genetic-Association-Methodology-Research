### Pilot Study Directory Readme
This is the readme file for for all experiments ran for the Pilot Study of my project on Kernel based association methodologies research. This study was run in order to determine which Beta parameter should be used when simulating my two phenotype models. Power was determined by using a Single Locus Test (SLT). Since I simulate continuous phenotypes, the SLT just runs a basic ANOVA where each SNP site is used as the predictor variable. The minimal p-value is then saved. For fun, I also calculated the proportion of times the SLT correctly determines which SNP site is the causal site.

## Directory List:

Directory | Description
--------- | ---------
Bash Scripts | Contains all BASH shell scripts used to submit jobs to Frontenac. All scripts are broken down into 5 distinct steps
R Scripts | Contains all R scripts used in this experiment
Final Results | Contains aggregated results files.


### Bash Scripts Directory

* __Pilot_Study_Step1.sh__
  * Uses the *ms* program to simulate haplotype and gene tree data. There were 10 000 iterations (not all files will be useable though).
  * The raw results from the _ms_ program are saved in the format of results{ i }.txt;
  where i = 1, …, 10 000. 
  
  * The raw data will be saved in the directory “Pilot_Study_Raw_Data”.


* __Pilot_Study_Step2.sh__
  * This sends an Array Job.
  * This bash script calls the Pilot_Study_Step2.R R script. This R script takes the raw ms data output and cleans it to produce one set of .txt files containing the haplotype data (haplodata{ i }.txt) and another set of .txt containing the tree data (treedata{ i }.txt).
  * These cleaned datasets will be saved in the Pilot_Study_Clean_Data directory.
  * FYI This took 8 minutes to complete.


* __Pilot_Study_Step3.sh__
  * This also is an Array Job. 
  * This script submits the Pilot_Study_Step3.R script for submission. This R script takes the haplotype data and simulates 2 sets of phenotype data. The parameter for simulation is based on 16 different Beta values, hence we get 16 different phenotypes for both phenotype models. 
  * The results will be saved in the Pilot_Study_PhenoAndGeno_Data directory as NoRecomb_PhenoAndGeno{ i }.txt files. Each text file will have 32 columns depicting the 32 different phenotypes simulated (16 for each type of phenotype) and the remaining columns depicting the genetic data.
 
  * Note: the code will remove the haplodata{ i }.txt files that are unable to produce the causal SNP sites, based on the conditions of the minor allele frequencies we have set.
  
  * This R script will also save which sites were chosen as the actual common causal and which sites were chosen as the actual rare causal in the Pilot_Study_Results directory.
  
  * FYI - This took 1 minute to complete.


* __Pilot_Study_Step4.sh__
  * This is another Array job.  This shell script which submits the Pilot_Study_Step4.R script for submission
  
  * Here, we finally analyse all of our data. That is, we determine the proportion of SLT tests that were statistically significant at the level of 5%. We do this across all phenotype models under each beta value.
  
  * Additionally, this script saves the chosen and rare causal site which were the most statistically significant. we will compare these sites to the sites that were actually causal.


* __Pilot_Study_Step5.sh__
  ** Since we use Array jobs, all of our results are saved in many different files! We concatenate all of our results files into a single file. 


## Final Results Directory

File name | Description
--------- | ---------
Bash Scripts | Contains all BASH shell scripts used to submit jobs to Frontenac. All scripts are broken down into 5 distinct steps
R Scripts | Contains all R scripts used in this experiment
Final Results | Contains aggregated results files.
