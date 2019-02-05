#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||       Crohn's Disease Data Analysis - Allele kernels only   ||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

.libPaths("/global/home/hpc4300/RPackages")
library(ape)
library(SKAT)
library(MDMR)


source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Tree/Calculate_Tree_kernels.R")
source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Tree/Tree_GTSR.R")
source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Tree/Tree_MDMR.R")



# Response variable
transmit=cbind(1:(4*129), rep(c(1,0),2*129))
status = transmit[,2]

#Read in the data:
setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Trees")




# Input file  
#treefilename=paste("Merge8Million_trees_",runnum,".out",sep="")
#trees=read.tree(treefilename)

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

treename = paste("Merge8Million_trees_", task_id, ".out", sep = "")
trees = read.tree(treename)

# Loop to compute relevant statistics
GTSR_results=matrix(ncol=5,nrow=length(trees))
MDMR_results=matrix(ncol=5,nrow=length(trees))
SKAT_results=matrix(ncol=5,nrow=length(trees))
K=floor(129*4/25)

setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Results")
G=read.table("Processed_haplodata.txt")
for (j in 1:length(trees)){
  print("This is tree #")
  print(j)
  mytree = trees[[j]]
  
  mytree_kernels = get.tree.kernels(mytree)

  #-->These are just the tree kernels...
  
  
  #-----||-----||-----|| - SKAT Code - ||-----||-----||-----#
  null.model = SKAT_Null_Model(status ~ 1, out_type="D")
  #Fit the null model
  
  SKAT.stat=vector("list", length(mytree_kernels))
  #Initialize list of SKAT p-values
  
  for (i in 1:length(mytree_kernels)){
    #Fill in list of SKAT p-values
    SKAT.stat[[i]]= tryCatch(SKAT(Z=as.matrix(G), obj=null.model, kernel = mytree_kernels[[i]])$Q,
                              error = function(e) 
                                paste("NA"))    
    #-----||-----||-----|| - What is tryCatch()?  - ||-----||-----||-----#
    #--> tryCatch() is implemented because in some datasets, we get an error.
    #--> Essentially, some product kernels seem to not produce any positive eigenvalues
    # and we can't compute p-values. This only happens a small amount of times.
    #--> With the tryCatch(), if this error occurs we assign an NA value and just move on without
    #disrupting the code.
    #-----||-----||-----|| - What is tryCatch()?  - ||-----||-----||-----#
  }
  for (i in (1:(length(SKAT.stat)))){
    #Add SKAT p-values to our results table.
    SKAT_results[j,i] = SKAT.stat[[i]]
  }
  
  
  #-----||-----||-----|| - Gene Trait Similarity Regression  - ||-----||-----||-----#
  stat.gtsm = similarity.regression.tree(P1=status,n=length(status), kernel.list=mytree_kernels)
  #Obtain list of GTSR p-values
  
  for (i in (1:length(stat.gtsm))){
    #Fill in results table with GTSR p-values
    GTSR_results[j,i] = stat.gtsm[i] 
  }
  
  
  #-----||-----||-----||-----|| - MDMR - ||-----||-----||-----||-----#
  stat.MDMR = pval_P1_MDMR_function(P1=status, kernel.list=mytree_kernels)
  #Obtain list of MDMR p-values
  
  for (i in (1:length(stat.MDMR))){
    #Fill in results table with MDMR p-values
    MDMR_results[j,i] = stat.MDMR[i]
  }
  colnames(SKAT_results) = c("Tree1","Tree2", "Tree3", "Tree4", "Tree5")
  colnames(MDMR_results) = c("Tree1","Tree2", "Tree3", "Tree4", "Tree5")
  colnames(GTSR_results) = c("Tree1","Tree2", "Tree3", "Tree4", "Tree5")
}


setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Results")


SKAT_name = paste("SKAT_part", task_id, ".txt", sep="")
GTSR_name = paste("GTSR_part", task_id, ".txt", sep="")
MDMR_name = paste("MDMR_part", task_id, ".txt", sep="")

write.table(SKAT_results, SKAT_name, quote=F, row=F, col=F)
write.table(GTSR_results, GTSR_name, quote=F, row=F, col=F)
write.table(MDMR_results, MDMR_name, quote=F, row=F, col=F)