.libPaths("/global/home/hpc4300/RPackages")

#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
#-----||-----||-----||Initialize Phenotype 1 data analysis results:||-----||-----||-----||-----#
#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
Phenotype1_Results_SKAT = data.frame(Rep. = integer(), Method = character(),IBS = numeric(),
                                     AM = numeric(), AS = numeric(), LIN = numeric(),
                                     REC = numeric(),QUAD = numeric(), S012 = numeric(),
                                     S123 = numeric(),S124 = numeric(), h1 = numeric(),
                                     Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
                                     Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )

Phenotype1_Results_Reg = data.frame(Rep. = integer(), Method = character(),IBS = numeric(),
                                    AM = numeric(), AS = numeric(), LIN = numeric(),
                                    REC = numeric(),QUAD = numeric(), S012 = numeric(),
                                    S123 = numeric(),S124 = numeric(), h1 = numeric(),
                                    Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
                                    Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )

Phenotype1_Results_MDMR = data.frame(Rep. = integer(), Method = character(),IBS = numeric(),
                                     AM = numeric(), AS = numeric(), LIN = numeric(),
                                     REC = numeric(),QUAD = numeric(), S012 = numeric(),
                                     S123 = numeric(),S124 = numeric(), h1 = numeric(),
                                     Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
                                     Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )

Phenotype1_Results_SLT = data.frame(Rep. = integer(), Method = character(),Pval = numeric(), stringsAsFactors = FALSE )


#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
#-----||-----||-----||Initialize Phenotype 2 data analysis results:||-----||-----||-----||-----#
#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
Phenotype2_Results_SKAT = data.frame(Rep. = integer(), Method = character(),IBS = numeric(),
                                     AM = numeric(), AS = numeric(), LIN = numeric(),
                                     REC = numeric(),QUAD = numeric(), S012 = numeric(),
                                     S123 = numeric(),S124 = numeric(), h1 = numeric(),
                                     Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
                                     Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )

Phenotype2_Results_Reg = data.frame(Rep. = integer(), Method =  character(),IBS = numeric(),
                                    AM = numeric(), AS = numeric(), LIN = numeric(),
                                    REC = numeric(),QUAD = numeric(), S012 = numeric(),
                                    S123 = numeric(),S124 = numeric(), h1 = numeric(),
                                    Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
                                    Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )

Phenotype2_Results_MDMR = data.frame(Rep. = integer(), Method = character(),IBS = numeric(),
                                     AM = numeric(), AS = numeric(), LIN = numeric(),
                                     REC = numeric(),QUAD = numeric(), S012 = numeric(),
                                     S123 = numeric(),S124 = numeric(), h1 = numeric(),
                                     Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
                                     Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )

Phenotype2_Results_SLT = data.frame(Rep. = integer(), Method = character(), Pval = numeric(), stringsAsFactors = FALSE )


#As double check, I want to save some results from the SLT code. I save which sites turn out
#to be the most significant, as well as its corresponding p-value.
Chosen_common_causal = list()
Chosen_rare_causal = list()
P1_pval_list=list()
P2_pval_list=list()


source("/global/home/hpc4300/BIM_Final_RCodes/BIM_Rcode_Calculate_all_kernels.R")
source("/global/home/hpc4300/BIM_Final_RCodes/BIM_RCode_NEW_gtsm.R")
source("/global/home/hpc4300/BIM_Final_RCodes/BIM_Rcode_MDMR_Code.R")
source("/global/home/hpc4300/BIM_Final_RCodes/BIM_RCode_SLT.R")


#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
#-----||-----||-----||-----||-----||Array Job Code:||-----||-----||-----||-----||-----||-----#
#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
#I set an aray job with 10 "arrays". I will split up the files to analyse into 10 chunks:
#There are 2500 files in total to potentially analyse, so I split this into 10 chunks with
#250 files belonging to each chunk.

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)
#Obtain Slurm Task ID.  


#Now, determine indices of data files to analyse:
total_files=seq(from=1, to= 2501, by=250)

starting = total_files[task_id]
  #Compute starting index

ending = total_files[task_id + 1] - 1
  #Compute ending index

for (index in (starting:ending)){
  Phenotype1_Results_SKAT[index,1] = Phenotype1_Results_Reg[index,1] = Phenotype1_Results_MDMR[index,1] = Phenotype1_Results_SLT[index,1]=Phenotype2_Results_SKAT[index,1] = Phenotype2_Results_Reg[index,1] = Phenotype2_Results_MDMR[index,1] = Phenotype2_Results_SLT[index,1]=index
  
  Phenotype1_Results_SKAT[index,2] = Phenotype2_Results_SKAT[index,2] = "SKAT"
  Phenotype1_Results_Reg[index,2] = Phenotype2_Results_Reg[index,2] = "GTSR"
  Phenotype1_Results_MDMR[index,2] = Phenotype2_Results_MDMR[index,2] = "MDMR"
  Phenotype1_Results_SLT[index,2] = Phenotype2_Results_SLT[index,2] = "SLT"
  
  table_name = paste("NoRecomb_PhenoAndGeno", index, ".txt", sep="")
  
setwd("/global/home/hpc4300/BIM_Final_PhenoAndGeno2_Data/")

  if(!file.exists(table_name)){
    print(table_name)
    next
  }
  
  
  data = read.table(table_name)
  
  G=data[,3:ncol(data)]
  #Extract genotype information.
  
  P1=data$V1
  #Extract Phenotype 1 information.
  
  P2=data$V2
  #Extract Phenotype 2 information.
  
  n=nrow(data)
  #Extract population size information.
  
  K=ncol(data)-2
  #Extract the number of segregation sites (SNP sites).
  
  treename = paste("treedata",index, ".txt", sep="")
  
setwd("/global/home/hpc4300/BIM_Final_Clean_Data/")
  mykernels = get.kernels(G=G, P1=P1, P2=P2, n=n, K=K, treename=treename)  
  #Calculate all kernel functions.
  
 
  #-----||-----||-----|| - ********************  - ||-----||-----||-----#
  #-----||-----||-----|| - Phenotype 1 Analysis - ||-----||-----||-----#
  #-----||-----||-----|| - ********************  - ||-----||-----||-----#
  
  
  #-----||-----||-----|| - SKAT Code - ||-----||-----||-----#
  null.model1 = SKAT_Null_Model(P1 ~ 1, out_type="C")
  #Fit the null model
  
  p.val.pheno1.SKAT=vector("list", length(mykernels))
  #Initialize list of SKAT p-values
  
  for (j in 1:length(mykernels)){
    #Fill in list of SKAT p-values
    p.val.pheno1.SKAT[[j]]= tryCatch(SKAT(Z=as.matrix(G), obj=null.model1, kernel = mykernels[[j]])$p.value,
                                   error = function(e) 
                                     paste("NA"))    
    #-----||-----||-----|| - What is tryCatch()?  - ||-----||-----||-----#
    #tryCatch() is implemented because in some datasets, we get an error.
    #Essentially, some product kernels seem to not produce any positive eigenvalues
    # and we can't compute p-values. This only happens a small maount of times.
    #With the tryCatch(), if this error occurs we assign an NA value and just move on without
    #disrupting the code.
    #-----||-----||-----|| - What is tryCatch()?  - ||-----||-----||-----#
  }
  for (j in (1:(length(p.val.pheno1.SKAT)))){
    #Add SKAT p-values to our results table.
    Phenotype1_Results_SKAT[index,j+2] = p.val.pheno1.SKAT[[j]]
  }
  
  
  #-----||-----||-----|| - Gene Trait Similarity Regression  - ||-----||-----||-----#
  p.val.pheno1.gtsm = similarity.regression.pheno1(P1=P1,n=n, kernel.list=mykernels)
  #Obtain list of GTSR p-values
  
  for (j in (1:length(p.val.pheno1.gtsm))){
    #Fill in results table with GTSR p-values
    Phenotype1_Results_Reg[index,j+2] = p.val.pheno1.gtsm[j] 
  }
  
  
  #-----||-----||-----||-----|| - MDMR - ||-----||-----||-----||-----#
  p.val.pheno1.MDMR = pval_P1_MDMR_function(P1=P1, kernel.list=mykernels)
  #Obtain list of MDMR p-values
  
  for (j in (1:length(p.val.pheno1.MDMR))){
    #Fill in results table with MDMR p-values
    Phenotype1_Results_MDMR[index,j+2] = p.val.pheno1.MDMR[j]
  }
  
  #-----||-----||-----|| - Single Locus Test  - ||-----||-----||-----#
  p.val.pheno1.SLT = get.min.pval(data=data, phenotype=P1)
  Phenotype1_Results_SLT[index,3]=   p.val.pheno1.SLT[[1]]
  #Populate our results table with SLT p-value
  
  Chosen_common_causal = c(Chosen_common_causal, list(p.val.pheno1.SLT[[2]]) )
  #Save which common causal variant was the most significant
  
  P1_pval_list=c(P1_pval_list, list(p.val.pheno1.SLT[[3]]))
  #Save its corresponsding pvalue as well.
 
  
  
  #-----||-----||-----|| - ********************** - ||-----||-----||------#  
  #-----||-----||-----|| - Phenotype 2 Analysis   - ||-----||-----||-----#
  #-----||-----||-----|| - ********************** - ||-----||-----||------#  
  #See comments in Phenotype 1 Analysis. The code here is pretty much the exact replica.
  
  #-----||-----||-----||-----|| - SKAT - ||-----||-----||-----||-----#
  null.model2 = SKAT_Null_Model(P2 ~ 1, out_type="C")
  p.val.pheno2.SKAT=vector("list", length(mykernels))

  
  for (j in 1:length(mykernels)){
    p.val.pheno2.SKAT[[j]]= tryCatch(SKAT(Z=as.matrix(G), obj=null.model2, kernel = mykernels[[j]])$p.value,
                                   error = function(e) 
                                     paste("NA"))   }
  
  
  for (j in (1:(length(p.val.pheno1.SKAT)))){
    Phenotype2_Results_SKAT[index,j+2] = p.val.pheno2.SKAT[[j]]
  }
  
  
  #-----||-----||-----|| - Gene Trait Similarity Regression  - ||-----||-----||-----#
  p.val.pheno2.gtsm = similarity.regression.pheno2(P2=P2,n=n, kernel.list=mykernels)
  
  for (j in (1:length(p.val.pheno2.gtsm))){
    Phenotype2_Results_Reg[index,j+2] = p.val.pheno2.gtsm[j] 
  }
  
  
  #-----||-----||-----||-----|| - MDMR - ||-----||-----||-----||-----#
  p.val.pheno2.MDMR = pval_P2_MDMR_function(P2=P2, kernel.list=mykernels)
  for (j in (1:length(p.val.pheno2.MDMR))){
    Phenotype2_Results_MDMR[index,j+2] = p.val.pheno2.MDMR[j]
  }
  
  #-----||-----||-----|| - Single Locus Test  - ||-----||-----||-----#
  p.val.pheno2.SLT = get.min.pval(data=data, phenotype=P2)
  Phenotype2_Results_SLT[index,3]=   p.val.pheno2.SLT[[1]]
  
  
  #-----||-----||-----|| - Sanity Check  - ||-----||-----||-----#
  Chosen_rare_causal = c(Chosen_rare_causal, list(p.val.pheno2.SLT[[2]]))
  P2_pval_list=c(P2_pval_list, list(p.val.pheno2.SLT[[3]]))
  
}

setwd("/global/home/hpc4300/BIM_Final_Results/")

Pheno1Results = rbind(Phenotype1_Results_SKAT, Phenotype1_Results_Reg, Phenotype1_Results_MDMR)
Pheno1Results = Pheno1Results[rowSums(is.na(Pheno1Results)) != ncol(Pheno1Results), ]

Pheno2Results = rbind(Phenotype2_Results_SKAT, Phenotype2_Results_Reg, Phenotype2_Results_MDMR)
Pheno1Results = Pheno2Results[rowSums(is.na(Pheno2Results)) != ncol(Pheno2Results), ]

P1_results_name = paste("Pheno1Results_", task_id, ".txt", sep="")
P2_results_name = paste("Pheno2Results_", task_id, ".txt", sep="")

write.table(Pheno1Results,P1_results_name,quote=F,row=F,col=F)
write.table(Pheno2Results,P2_results_name,quote=F,row=F,col=F)




Phenotype1_Results_SLT = Phenotype1_Results_SLT[rowSums(is.na(Phenotype1_Results_SLT)) != ncol(Phenotype1_Results_SLT), ]
Phenotype2_Results_SLT = Phenotype2_Results_SLT[rowSums(is.na(Phenotype2_Results_SLT)) != ncol(Phenotype2_Results_SLT), ]


P1_SLT = paste("Pheno1Results_", task_id, "_SLT.txt", sep="")
P2_SLT = paste("Pheno2Results_", task_id, "_SLT.txt", sep="")
write.table(Phenotype1_Results_SLT,P1_SLT,quote=F,row=F,col=F)
write.table(Phenotype2_Results_SLT,P2_SLT,quote=F,row=F,col=F)


write.table(Chosen_common_causal, "Chosen_common_causal.txt", quote = F, row=F, col=F)
write.table(Chosen_rare_causal, "Chosen_rare_causal.txt", quote = F, row=F, col=F)
