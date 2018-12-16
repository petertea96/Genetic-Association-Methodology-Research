.libPaths("/global/home/hpc4300/RPackages")
set.seed(55)
#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
#-----||-----||-----||Initialize Phenotype 1 data analysis results:||-----||-----||-----||-----#
#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #

Phenotype1_Results_SLT = data.frame(matrix(ncol = 18, nrow = 1))


#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #
#-----||-----||-----||Initialize Phenotype 2 data analysis results:||-----||-----||-----||-----#
#  -----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----   #

Phenotype2_Results_SLT = data.frame(matrix(ncol = 18, nrow = 1))

#As double check, I want to save some results from the SLT code. I save which sites turn out
#to be the most significant, as well as its corresponding p-value.
Chosen_common_causal = vector()
Chosen_rare_causal = vector()
P1_pval_list=list()
P2_pval_list=list()



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
total_files=seq(from=1, to= 10001, by=100)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id + 1] - 1
#Compute ending index

for (index in (starting:ending)){
  Phenotype1_Results_SLT[index,1] = Phenotype2_Results_SLT[index,1] = index
  Phenotype1_Results_SLT[index,2] = Phenotype2_Results_SLT[index,2] = "SLT"
  
  table_name = paste("NoRecomb_PhenoAndGeno", index, ".txt", sep="")
  
  setwd("/global/home/hpc4300/Pilot_Study/Pilot_Study_PhenoAndGeno_Data/")
  
  if(!file.exists(table_name)){
    print(table_name)
    next
  }
  
  
  data = read.table(table_name)
  
  G=data[,33:ncol(data)]
  #Extract genotype information.
  
  K=ncol(data)-32
  #Extract the number of segregation sites (SNP sites).


    #Extract Phenotype 2 information.
  

  
  
  #-----||-----||-----|| - ********************  - ||-----||-----||-----#
  #-----||-----||-----|| - Phenotype 1 Analysis - ||-----||-----||-----#
  #-----||-----||-----|| - ********************  - ||-----||-----||-----#


  #-----||-----||-----|| - Single Locus Test  - ||-----||-----||-----#
  #Populate our results table with SLT p-value

  
  for(h in 1:16){
    P1 = data[,h]
    p.val.pheno1.SLT = get.min.pval(data=data[,31:ncol(data)], phenotype=P1)
    Phenotype1_Results_SLT[index,2+h]=   p.val.pheno1.SLT[[1]]
    
    Chosen_common_causal = c(Chosen_common_causal, p.val.pheno1.SLT[[2]]) 
    #Save which common causal variant was the most significant
  }
  
  if (length(Chosen_common_causal)%%16 !=0){
    print(index)
    print("After adding this index, the vector is no longer a multiple of 16")
  }
  
  
  #-----||-----||-----|| - ********************** - ||-----||-----||------#  
  #-----||-----||-----|| - Phenotype 2 Analysis   - ||-----||-----||-----#
  #-----||-----||-----|| - ********************** - ||-----||-----||------#  
  
  #-----||-----||-----|| - Single Locus Test  - ||-----||-----||-----#
  
  for(h in 1:16){
    P2 = data[,h+16]
    p.val.pheno2.SLT = get.min.pval(data=data[,31:ncol(data)], phenotype=P2)
    Phenotype2_Results_SLT[index,2+h]=   p.val.pheno2.SLT[[1]]
    
    Chosen_rare_causal = c(Chosen_rare_causal, p.val.pheno2.SLT[[2]])
    #Save which common causal variant was the most significant
  }
}
  
 
setwd("/global/home/hpc4300/Pilot_Study/Pilot_Study_Results")

Phenotype1_Results_SLT = Phenotype1_Results_SLT[rowSums(is.na(Phenotype1_Results_SLT)) != ncol(Phenotype1_Results_SLT), ]
Phenotype2_Results_SLT = Phenotype2_Results_SLT[rowSums(is.na(Phenotype2_Results_SLT)) != ncol(Phenotype2_Results_SLT), ]


P1_SLT = paste("Pheno1Results_", task_id, "_SLT.txt", sep="")
P2_SLT = paste("Pheno2Results_", task_id, "_SLT.txt", sep="")
write.table(Phenotype1_Results_SLT,P1_SLT,quote=F,row=F,col=F)
write.table(Phenotype2_Results_SLT,P2_SLT,quote=F,row=F,col=F)

chosen_common_causal_name = paste("Chosen_common_causal_", task_id, ".txt", sep="")
write.table(Chosen_common_causal, chosen_common_causal_name, quote = F, row=F, col=F)

chosen_rare_causal_name = paste("Chosen_rare_causal_", task_id, ".txt", sep="")
write.table(Chosen_rare_causal, chosen_rare_causal_name, quote = F, row=F, col=F)
