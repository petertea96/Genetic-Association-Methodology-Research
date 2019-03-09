#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####             Allele/Genotype Analysis at different Focal Points        #####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#Today is March 9th, 2019
.libPaths("/global/home/hpc4300/RPackages")

setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/JobFiles")



for (file_number in (1:100)){
  directory_name = paste("Run_", file_number)
  
  haplodat = read.delim("crohn5q31_haplo.dat", colClasses = "character")
  
  
  FirstLine=unlist(strsplit(haplodat[1,1],split=""))
  segsites=length(FirstLine)
  
  
  newhaplodat=matrix(as.numeric(unlist(strsplit(haplodat[,1],split=""))),
                     ncol=segsites,byrow=T)
  
  
  #STEP 2: Make sure our data is in terms of the Minor Allele.
  
  ## Get allele frequencies
  f1=colSums(newhaplodat)/nrow(newhaplodat)
  
  ## We want allele frequencies to be in terms of minor allele frequencues (MAF)
  # So.... if MAF>0.5, we ned to reverse the coding:
  tochange=which(f1>0.5)
  
  if (length(tochange) !=0){
    for (i in 1:length(tochange)){
      index=tochange[i]
      newhaplodat[,index]=1-newhaplodat[,index]
    }
  }
  
  
  #Code for phenotypes:
  #Data set up such that the binary phenotype alternates between 1 and 0 for consecutive rows.
  if (nrow(newhaplodat) %% 2 == 1){
    status_phenotype=rep(c(1,0),ceiling( nrow(newhaplodat)/2) )[-length(status_phenotype)]
  }
  else {
    status_phenotype=rep(c(1,0),nrow(newhaplodat)/2)
  }
    
  
  #Source some scripts to help calculate Kernels and to calculate assoc. statistics
  source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Crohn_Calculate_all_kernels.R")
  source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Crohn_NEW_gtsm.R")
  source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Crohn_MDMR_Code.R")
  
  
  mykernels = get.kernels(G=newhaplodat, n=nrow(newhaplodat), K=ncol(newhaplodat), treename="blah")  
  
  Phenotype_Results = data.frame(Method = character(),IBS = numeric(),
                                 AM = numeric(), AS = numeric(), h1 = numeric(),
                                 Skat = numeric(), stringsAsFactors = FALSE )
  
  #Phenotype_Results = data.frame(Method = character(),IBS = numeric(),
  #                                     AM = numeric(), AS = numeric(), h1 = numeric(),
  #                                     Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
  #                                     Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )
  
  
  
  
  #-----||-----||-----|| - SKAT Code - ||-----||-----||-----#
  null.model = SKAT_Null_Model(status_phenotype ~ 1, out_type="D")
  #Fit the null model
  
  p.val.SKAT=vector("list", length(mykernels))
  #Initialize list of SKAT p-values
  
  for (j in 1:length(mykernels)){
    #Fill in list of SKAT p-values
    p.val.SKAT[[j]]= tryCatch(SKAT(Z=as.matrix(G), obj=null.model, kernel = mykernels[[j]])$p.value,
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
  
  for (j in (1:(length(p.val.SKAT)))){
    #Add SKAT p-values to our results table.
    Phenotype_Results[1,j+1] = p.val.SKAT[[j]]
  }
  
  
  #-----||-----||-----|| - Gene Trait Similarity Regression  - ||-----||-----||-----#
  p.val.gtsm = similarity.regression.pheno1(P1=status_phenotype,n=nrow(newhaplodat), kernel.list=mykernels)
  #Obtain list of GTSR p-values
  
  for (j in (1:length(p.val.gtsm))){
    #Fill in results table with GTSR p-values
    Phenotype_Results[2,j+1] = p.val.gtsm[j] 
  }
  
  
  #-----||-----||-----||-----|| - MDMR - ||-----||-----||-----||-----#
  p.val.MDMR = pval_P1_MDMR_function(P1=status_phenotype, kernel.list=mykernels)
  #Obtain list of MDMR p-values
  
  for (j in (1:length(p.val.MDMR))){
    #Fill in results table with MDMR p-values
    Phenotype_Results[3,j+1] = p.val.MDMR[j]
  }
  
  labels = c("SKAT", "GTSM", "MDMR")
  for (i in 1:length(labels)){
    Phenotype_Results[i,1] = labels[i]
  } 
  
  write.table(Phenotype_Results,paste("New_Allele_Kernel_Data", file_number, ".txt"),quote=F,row=F,col=F)
  
  
  
  
  
  
  
  
} 









