### Pilot Study Directory Readme
This is the readme file for for all experiments ran for the Pilot Study of my project on Kernel based association methodologies research. This study was run in order to determine which Beta parameter should be used when simulating my two phenotype models. Power was determined by using a Single Locus Test (SLT). Since I simulate continuous phenotypes, the SLT just runs a basic ANOVA where each SNP site is used as the predictor variable. The minimal p-value is then saved. For fun, I also calculated the proportion of times the SLT correctly determines which SNP site is the causal site.

## Directory List:

Directory | Description
--------- | ---------
Bash Scripts | Contains all BASH shell scripts used to submit jobs to Frontenac. All scripts are broken down into 5 distinct steps
Final Results | Contains aggregated results files.
