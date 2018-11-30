#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#These functions are used to simulate the two phenotype models.
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

get_common_causal = function(UpperBound=0.35, LowerBound=0.25, genodata){
  #This function randomly chooses one SNP site with proper MAF to be the common causal site.
  
  Minor_Allele_Frequencies = colSums(genodata)/(2*nrow(genodata))
  possible_common_causal_sites = which( (Minor_Allele_Frequencies < UpperBound) & (Minor_Allele_Frequencies > LowerBound ) )
  #Create vector of possible SNP sites to be the common causal variant.
  
  if (length(possible_common_causal_sites)==1){
    Causal_site = possible_common_causal_sites[1]
    #Running sample() on a vector of length 1 will cause problems...
    
  } 
  
  else if (length(possible_common_causal_sites)==0){
    return(NULL)
  }  
  
  else{
    Causal_site = sample(possible_common_causal_sites, 1) 
  }
  
  return(Causal_site)
}

#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

get_rare_causals = function(UpperBound=0.05, LowerBound=0, genodata, numsites=10){
  #This function randomly chooses 10 SNP site with proper MAF to be the rare causal sites.
  
  Minor_Allele_Frequencies = colSums(genodata)/(2*nrow(genodata))
  
  possible_common_rare_sites = which( (Minor_Allele_Frequencies < UpperBound) & (Minor_Allele_Frequencies > LowerBound ) )
  #Create vector of possible SNP sites to be the rare causal variants.
  
  
   if (length(possible_common_rare_sites)< numsites){
    return(NULL)
  }
  
  rare_sites = sample(possible_common_rare_sites, numsites)
  
  return(rare_sites)
}

#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
sim_pheno1 = function(beta, common.causal){
  #This function simulates the phenotype 1 model.
  
  phenotype1 = rnorm(n=100, mean = genodat[,common.causal]*beta, sd=1)
  return(phenotype1)
}

#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
sim_pheno2 = function(beta, rare.causal){
  #This function simulates the phenotype 2 model.
  
  phenotype2 = rnorm(n = 100, mean = rare.causal*beta, sd=1)
  return(phenotype2)
}