## PSD Error Files Directory
PSD = Positive Semi-Definite

File name | Description
--------- | ---------
NoRecomb_PhenoAndGeno{i}.txt | Phenotype and genotype data of a data file that experienced PSD issues
treedata{ i }.txt	 | Gene tree data of a data file that experienced PSD issues
PSD_RScript.R | R script that attempts to compute all 16 kernel functions of interest. It prints an error message whenever we produce a kernel that is not PSD.
