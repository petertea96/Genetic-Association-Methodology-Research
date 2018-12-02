#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#These functions are used to simulate the two phenotype models.
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#This R file is important because...
#I need a way to choose which SNP site is common causal for simulating phenotype 1
#And I need to choose which SNP sites will be rare causal for phenotype 2.




get_common_causal = function(UpperBound=0.35, LowerBound=0.25, genodata){
  #This function randomly chooses one SNP site with proper MAF to be the common causal site.
  #By default, my decision criteria for the common causal site:
  #Common causal MAF: 0.25 < MAF < 0.35
  
  
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
  #By default, my decision criteria for the rare causal sites:
  #Rare causal MAF: 0 < MAF < 0.01; Choose the same # of rare causals for ALL simulations.
  
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