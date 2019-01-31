#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||  Project Data Analysis   ||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#


#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----|| Goal: Compare the Sensitivity of each test at the           ||-----||-----#
#-----||-----|| significance level of 0.05.                                 ||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

#--> Read in the data:
setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's project/Cluster files/Final codes/Genetic-Association-Methodology-Research")

Kernel_P1_Results = read.table("Aggregate_Pheno1Results.txt")
Kernel_P2_Results = read.table("Aggregate_Pheno2Results.txt")

SLT_P1_Results = read.table("Aggregate_Pheno1Results_SLT.txt")
SLT_P2_Results = read.table("Aggregate_Pheno2Results_SLT.txt")


column_names = c("Iteration", "Assoc_Stat", "S_IBS", "S_AM", "S_AS", "S_LIN", "S_REC", "S_QUAD", "S_012", "S_123", "S_124",
                 "S_h1", "skat.kernel","tree1", "tree2", "tree3", "tree4", "tree5")
colnames(Kernel_P1_Results) = column_names; colnames(Kernel_P2_Results) = column_names


# --> Remove rows with NA
Kernel_P1_Results = Kernel_P1_Results[complete.cases(Kernel_P1_Results),]
Kernel_P2_Results = Kernel_P2_Results[complete.cases(Kernel_P2_Results),]
SLT_P1_Results = SLT_P1_Results[complete.cases(SLT_P1_Results),]
SLT_P2_Results = SLT_P2_Results[complete.cases(SLT_P2_Results),]


# --> Re-write data tables
colnames(Kernel_P1_Results) = colnames(Kernel_P2_Results) = column_names
colnames(SLT_P1_Results) =   colnames(SLT_P2_Results) = c("Iteration", "Assoc_Stat", "P-val")

write.table(x = Kernel_P1_Results, file="Clean_Pheno1Results.txt")
write.table(x = Kernel_P2_Results, file="Clean_Pheno2Results.txt")
write.table(SLT_P1_Results, file="Clean_SLT1Results.txt")
write.table(SLT_P2_Results, file="Clean_SLT2Results.txt")
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||   Kernel Analysis ||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
library(dplyr)

### Phenotype 1 Analysis ###
P1_SKAT = filter(Kernel_P1_Results, Assoc_Stat == "SKAT")
P1_MDMR = filter(Kernel_P1_Results, Assoc_Stat == "MDMR")
P1_GTSR = filter(Kernel_P1_Results, Assoc_Stat == "GTSR")

nrow(P1_SKAT);nrow(P1_MDMR); nrow(P1_GTSR)


Kernel_Sensitivity = function(mydata){
  #Create matrix that checks for each element if the p-value is statisticcally significant.
  Pvalues_which_significant = (mydata[,-c(1,2)] < 0.05) 
  
  n = nrow(mydata)
  Sensitivity = colSums(Pvalues_which_significant)/n
  
  
  Kernel_names = c("S_IBS", "S_AM", "S_AS", "S_LIN", "S_REC", "S_QUAD", "S_012", "S_123", "S_124",
                   "S_h1", "skat.kernel","tree1", "tree2", "tree3", "tree4", "tree5")

  names(Sensitivity) = Kernel_names
  
  return(Sensitivity)
  
}

round(Kernel_Sensitivity(P1_SKAT), 3)
round(Kernel_Sensitivity(P1_GTSR), 3)
round(Kernel_Sensitivity(P1_MDMR), 3)

### Phenotype 1 Analysis ###
P2_SKAT = filter(Kernel_P2_Results, Assoc_Stat == "SKAT")
P2_MDMR = filter(Kernel_P2_Results, Assoc_Stat == "MDMR")
P2_GTSR = filter(Kernel_P2_Results, Assoc_Stat == "GTSR")

nrow(P2_SKAT); nrow(P2_MDMR); nrow(P2_GTSR)

round(Kernel_Sensitivity(P2_SKAT),3); 
round(Kernel_Sensitivity(P2_MDMR), 3); 
round(Kernel_Sensitivity(P2_GTSR), 3)


#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||SLT Analysis||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

#--> Create matrix that checks for each element if the p-value is statisticcally significant.
SLT_P1_pvalues_which_significant = (SLT_P1_Results[,-c(1,2)] < 0.05) 

#--> Calculate Sensitivity
Sensitivity_SLT_P1 = sum (SLT_P1_pvalues_which_significant) / nrow(SLT_P1_Results) 
Sensitivity_SLT_P1

######### -- SLT for Phenotype 2 -- ##########

#--> Create matrix that checks for each element if the p-value is statisticcally significant.
SLT_P2_pvalues_which_significant = (SLT_P2_Results[,-c(1,2)] < 0.05) 

#--> Calculate Sensitivity
Sensitivity_SLT_P2 = sum (SLT_P2_pvalues_which_significant) / nrow(SLT_P2_Results) 
Sensitivity_SLT_P2






