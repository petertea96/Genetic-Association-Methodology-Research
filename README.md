# UOttawa-Honours-Project
This repository contains work that was completed during a summer funded research project (Natural Sciences and Engineering Research Council - Undergraduate Student Research Award)  in 2018. It is also continuing as an Honour's Project which will eventually be submitted to the faculty of Science upon completion. 


# Project Description

> Recent advancements in sequencing technologies have made it easier to identify both rare and common genetic variants in the human genome. Along with these sequencing improvements, genetic association studies have become more prominent and are used to test for association between genetic variation and a phenotype of interest such as a disease state. One prospering field of genetic association tests are kernel-based association methodologies. These methodologies first require specification of a kernel function, which outputs a map that describes the degree of genetic similarity between pairs of individuals; many kernel functions have been proposed with strategies ranging from scoring genotype similarity to tree-based approaches. Then, using these maps, a variety of different kernel statistics can be employed to measure the strength of association between genetic similarity with a trait of interest. Due to the rapid expansion of this field, there has been no study that has described and compared the performances of all these different types of kernel-based association statistics in conjunction with the different kernel functions. 


## Directory List:

Directory | Description
--------- | ---------
Bash Scripts | Contains all BASH shell scripts used to submit jobs to Frontenac. 
R Scripts | Contains all R scripts used in this Project
Final Results | Contains aggregated results files.
PSD Error files | Contains files (from a previous pilot study), that contains data files that have issues with PSD.


## BASH Scripts Directory

File name | Description
--------- | ---------
BIM_Final_Step1_get_data.sh | ...
BIM_Final_Step2_clean_data.sh	 | ...
BIM_Final_step4_analyse_data.R | ...
BIM_Final_Step3_simulate_phenotype_data.sh | ...
BIM_Final_Step4_analyse_data.sh | ...


## R Scripts Directory

File name | Description
--------- | ---------
BIM_Final_step2_prelim_codes.R | ...
BIM_Final_step3_simulate_phenotype_data.R	 | ...
BIM_Final_step4_analyse_data.R | ...
BIM_RCode_NEW_gtsm.R | Written by Zhe Gao with slight modifications made by me...
BIM_RCode_SLT.R | ...
BIM_Rcode_Calculate_all_kernels.R | ...
BIM_Rcode_Function.R | This R file was taken from the SKAT package. I needed it since it contains some functions to compute the SKAT kernel.
BIM_Rcode_MDMR_Code.R | Written by Zhe Gao. 
BIM_Rcode_ORDER.R | The distance matrix outputted from “treeSimilarity()” doesn’t order the haplotypes correctly. For example, we assume that hap1 and hap2 belong to individual 1. To order the rows and columns to go from hap1 to hap200, we use the “order()” function contained in the “ORDER.R” file.
BIM_Rcode_SKAT_Linear.R | This R file was takeen from the SKAT package. I needed it since it contains some functions to compute the SKAT kernel.
BIM_Rcode_Simulation_help.R | ...
BIM_Rcode_Solution_function.R | The kernel matrices we are studying should be individual based and not haplotype based. The tree matrices so far are haplotype based. To truncate these (200x200) matrices to (100x100) matrices.
BIM_Rcode_treesimilarityMODIFIED.R | Written by Kelly Burkett. This script takes tree data and can output any of the 5 specified tree kernels.
