#This R script is created to further investigate errors in Positive Semi-Definiteness (PSD).
#Some generated genotype data still return errors about PSD.

#First source these R files...
source("BIM_RCode_treesimilarityMODIFIED.R") # Calculate distance matrices
source("BIM_RCode_SKAT_Linear.R") #SKAT package
source("BIM_RCode_Function.R") # From SKAT package
source("BIM_RCode_ORDER.R") # Order the matrix from 1 - 200.
source("BIM_RCode_Solution_function.R") #Get the max. distance for each individual.
source("BIM_RCode_Calculate_all_kernels.R")
setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's Project/PSD error files")

#File numbers that have error codes: 264,520,678,775, 793

#I could have written for loop, but just modify the numbers below to look @ the different
#data files...
data1=read.table("NoRecomb_PhenoAndGeno793.txt")
G=data1[,3:ncol(data1)]
P1=data1$V1
P2=data1$V2
n=nrow(data1)
K=ncol(data1)-2

my.k1=get.kernels(G,P1, P2, n, K, treename="treedata793.txt")

library(SKAT)
null.model1 = SKAT_Null_Model(P1 ~ 1, out_type="C")
p.val.pheno1.SKAT= vector("list", length(my.k1))

for (j in 1:length(my.k1)){
  p.val.pheno1.SKAT[[j]]= tryCatch(SKAT(Z=as.matrix(G), obj=null.model1, kernel = my.k1[[j]])$p.value,
                                   error = function(e) 
                                     paste("PSD ERROR HERE"))
}

print(p.val.pheno1.SKAT)

