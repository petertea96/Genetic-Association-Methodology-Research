#####################################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Apply the MANY different kernels to the MDMR Statistic                  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#####################################################################################################
#June 18th, 2018
#MDMR
#I'm just replicating Zhe's results, while also inputting the new tree kernels for comparison.

.libPaths("/global/home/hpc4300/RPackages")
library(MDMR)
##---MDMR -- -##
pval_P1_MDMR_function = function(P1, kernel.list){
  p.value_MDMR_p1 <- matrix(rep (0, length(kernel.list)), nrow = 1)
  
  ## phenotype 1
  for (i in 1:length(kernel.list)) {
    kernel = kernel.list[[i]]
    MDMR <-mdmr(X=P1 , D = kernel, nperm=0)
    p.value_MDMR_p1[i] <- MDMR$stat[2 ,1]
  }
  #chr=c("S_IBS", "S_AM", "S_AS", "S_h1", "skat.kernel","tree1", "tree2", "tree3", "tree4", "tree5")
  chr=c("tree1", "tree2", "tree3", "tree4", "tree5")
  colnames(p.value_MDMR_p1) <-chr
  return(as.vector(p.value_MDMR_p1))
  
}
