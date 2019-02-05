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
##---MDMR -- -##
pval_P1_MDMR_function = function(P1, kernel.list){
  p.value_MDMR_p1 <- matrix(rep (0, length(kernel.list)), nrow = 1)
  
  ## phenotype 1
  for (i in 1:length(kernel.list)) {
    kernel = kernel.list[[i]]
    MDMR <-mdmr(X=P1 , D = kernel, perm.p = TRUE)
    p.value_MDMR_p1[i] <- MDMR$pv[2 ,1]
  }
  #chr=c("S_IBS", "S_AM", "S_AS", "S_h1", "skat.kernel","tree1", "tree2", "tree3", "tree4", "tree5")
  chr=c("S_IBS", "S_AM", "S_AS", "S_h1", "skat.kernel")
  colnames(p.value_MDMR_p1) <-chr
  return(as.vector(p.value_MDMR_p1))
  
}

pval_P2_MDMR_function = function(P2,kernel.list){
  p.value_MDMR_p2 <- matrix(rep (0, length(kernel.list)), nrow = 1)
  
  ## phenotype 2
  for (i in 1: length(kernel.list)) {
    kernel = kernel.list[[i]]
    MDMR <-mdmr (X=P2 , D = kernel)
    p.value_MDMR_p2[i] <-MDMR $pv[2 ,1]
  }
  chr=c("S_IBS", "S_AM", "S_AS", "S_LIN", "S_REC", "S_QUAD", "S_012", "S_123",
        "S_124", "S_h1", "skat.kernel","tree1", "tree2", "tree3", "tree4", "tree5")
  colnames(p.value_MDMR_p2) <-chr
  return(as.vector(p.value_MDMR_p2))
  
}








