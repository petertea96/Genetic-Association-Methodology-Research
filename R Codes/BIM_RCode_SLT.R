get.min.pval = function(data, phenotype){
  #Function computes all association tests for all SNP sites, and then outputs
  #the maximum of these association tests.
  
  number.sites = ncol(data) - 2 #Compute the number of SNP sites in the data
  Pval.list = numeric(number.sites) #Initialize list of p-values
  
  for(i in 1:number.sites){
    x = as.factor(data[,i+2])
    p.val = anova(lm(phenotype~x))[1,5] #Extract p-value.
    Pval.list[i] = p.val
  }
  which_SNP = which(Pval.list==min(Pval.list))
  if(length(which_SNP) > 1){
    print("WARNING, multiple SNP sites with the minimum p-value")
  }
  corrected_pvalue = number.sites*min(Pval.list)
  myresult = list(corrected_pvalue, which_SNP[1], Pval.list)
  return(myresult) #Return correction factor
  
}
