#This R script calculate all kernel matrices...
#Today is July 26th, 2018


.libPaths("/global/home/hpc4300/RPackages")
#Set directory where all R packages are in CAC cluster

require("ape")

#Source some R codes needed to calculate the kernels...

source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Tree/Crohn_Treesimilarity.R")
source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Tree/Crohn_order.R")
source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Tree/Tree_MDMR.R")




get.tree.kernels = function(mytree){
  
  # Various (tree) similarity matrices 
  sim1=treeSimilarity(mytree) # TMRCA - time to present of pairs mrca
  sim2=treeSimilarity(mytree,propn=TRUE) # (TMRCA - time to present of pairs mrca)/TMRCA
  sim3=treeSimilarity(mytree, method="ranks") # Rank of coalescent event for each pairs mrca (oldest ranked 1)
  sim4=treeSimilarity(mytree, method="sdscaled") # Like 1, but intercoalescence times divided by mean/sd 
  sim5=treeSimilarity(mytree, method="sdscaled",propn=T) # Like 2, but intercoalescence times divided by mean/sd 
  
  #Sorted issue:
  tree1=order(sim1[[1]], sim1[[2]])
  tree2=order(sim2[[1]], sim2[[2]])
  tree3=order(sim3[[1]], sim3[[2]])
  tree4=order(sim4[[1]], sim4[[2]])
  tree5=order(sim5[[1]], sim5[[2]])
  

  
  
  
  
  #Kernel Matrices
  kernel.list = list(tree1, tree2, tree3, tree4, tree5)
  #kernel.list = list(S_IBS, S_AM, S_AS, S_h1, skat.kernel,
  #                  tree1, tree2, tree3, tree4, tree5)
  names(kernel.list) = c("tree1", "tree2", "tree3",
                         "tree4", "tree5")
  return(kernel.list)
  #names(kernel.list) = c("S_IBS", "S_AM", "S_AS",
  #                      "S_h1", "skat.kernel","tree1", "tree2", "tree3", "tree4", "tree5")
}


