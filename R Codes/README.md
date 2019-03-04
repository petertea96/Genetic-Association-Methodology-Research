
## R Codes Directory

File name | Description
--------- | ---------
BIM_Final_step2_prelim_codes.R | Processes the **ms** simulated data ouput into two separate .txt files: One for haplotype data and the other for Tree data 
BIM_Final_step3_simulate_phenotype_data.R	 | More data processing, and simulation of 2 continuous phenotypes
BIM_Final_step4_analyse_data.R | Analysis of the simulated data
BIM_RCode_NEW_gtsm.R | Written by Zhe Gao with slight modifications made by me... this script calculates the Gene Trait Similarity Regression Kernel statistic.
BIM_RCode_SLT.R | Written function returns a bonferonni corrected single locus test (since phenotype is continuous, this is an ANOVA)
BIM_Rcode_Calculate_all_kernels.R | Calculate all kernels (both allele and tree based)
BIM_Rcode_MDMR_Code.R | Written by Zhe Gao. This script calculates the MDMR Kernel statistic
BIM_Rcode_ORDER.R | The distance matrix outputted from “treeSimilarity()” doesn’t order the haplotypes correctly. For example, we assume that hap1 and hap2 belong to individual 1. To order the rows and columns to go from hap1 to hap200, we use the “order()” function contained in the “ORDER.R” file.
BIM_Rcode_Simulation_help.R | Contains functions used to assign MAFs for simulation of continuous phenotypes
BIM_Rcode_Solution_function.R | The kernel matrices we are studying should be individual based and not haplotype based. The tree matrices so far are haplotype based. To truncate these (200x200) matrices to (100x100) matrices.
BIM_Rcode_treesimilarityMODIFIED.R | Written by Kelly Burkett. This script takes tree data and can output any of the 5 specified tree kernels.



Please note that to fully replicate the work, you will require two additional R files found in the SKAT package:

SKAT R File | Description
--------- | ---------
Function.R | This R file was taken from the SKAT package. I needed it since it contains some functions to compute the SKAT kernel
SKAT_Linear.R | This R file was takeen from the SKAT package. I needed it since it contains some functions to compute the SKAT kernel.