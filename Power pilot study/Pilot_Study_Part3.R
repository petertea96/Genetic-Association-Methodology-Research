set.seed(5)
common_causal_vector = vector()
rare_causal_list = list()
#For precaution, I will save which SNP sites is the common causal and which ones are
#rare causal. This will be compared to when I run the Single Locus Tests.


source("/global/home/hpc4300/BIM_Final_RCodes/BIM_Rcode_Simulation_help.R")

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)
#Obtain Slurm Task ID.  


#Now, determine indices of data files to analyse:
total_files=seq(from=1, to= 10001, by=100)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

for (j in (starting:ending)){
  setwd("/global/home/hpc4300/Pilot_Study/Pilot_Study_Clean_Data")
  
  ## Read in the haplotype data. Must specify that it is type "character". 
  filename = paste("haplodata",j, ".txt", sep="")
  
  if (!file.exists(filename)){
    next
    #This is just an added sanity check. To proceed, we need to ensure the file actually exists...
  }
  
  
  haplodat=read.table(filename, colClasses=c("character"))
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  #-----||-----||-----||-----||-----||Obsolete Code||-----||-----||-----||-----||-----||-----#
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  #The following code was used to obtain the number of SNP sites. But, further research (by me)
  # has shown that some SNPs are perfectly correlated. We will only keep SNP sites that are not
  # perfectly correlated.
  #
  #
  #The first line of the haplotype data will be used to calculate the number of segregation sites,
  #otherwise known as the number of SNP sites.
  FirstLine = readLines(filename)[1]
  FirstLine=unlist(strsplit(FirstLine,split=""))
  segsites=length(FirstLine)
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  
  
  ## Note that the loci are not separated by spaces. So we must split them in order to be able to manipulate the data. 
  ## Also convert to type numeric so that we can add allele counts
  newhaplodat=matrix(as.numeric(unlist(strsplit(haplodat[,1],split=""))),ncol=segsites,byrow=T)
  
  
  ##The following code can be used to check if our data manipulation above succeeded:
  #haplodat[24,]
  #newhaplodat[24,]
  #haplodat[186,]
  #newhaplodat[186,]
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  
  
  ## Get allele frequencies
  f1=colSums(newhaplodat)/nrow(newhaplodat)
  
  ## We want allele frequencies to be in terms of minor allele frequencues (MAF)
  # So.... if MAF>0.5, we ned to reverse the coding:
  tochange=which(f1>0.5)
  
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  #***I ADDED THE FOLLOWING IF STATEMENT IN CASE tochange RETURNS logical(0)
  #I.E. We do not need to change our data format if it already in terms of the MAF...
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  if (length(tochange) !=0){
    for (i in 1:length(tochange)){
      index=tochange[i]
      newhaplodat[,index]=1-newhaplodat[,index]
    }
  }
  
  
  
  ## Now create genotype dataset. This is done by pairing consecutive rows
  genodat=matrix(nrow=nrow(haplodat)/2,ncol=ncol(newhaplodat))
  
  for (i in 1:nrow(genodat)){
    genodat[i,]=newhaplodat[2*i-1,]+newhaplodat[2*i,]
  }
  
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  #What if columns are the exact same? (I.E 100% correlated between SNPs?) Let's remove them.
  
  genodat = unique.matrix(t(genodat)) 
  #Returns a matrix with only unique ROWS (so I need the transpose)
  
  genodat = t(genodat)
  #*** If this section was too confusing, I've added a more in-depth explanation in my ReadME file.
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  
  segsites = ncol(genodat)
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  ##Again, we can use the following code to check if our data manipulation above succeeded:
  #i=99
  #newhaplodat[(2*i-1):(2*i),]
  #genodat[i,]
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  
  
  
  ## Simulate some continuous phenotype data. Both models chosen so that power
  ## Assume that minor allele of causal increases phenotype mean by ___SD
  
  
  #-----||-----||-----||-----||Simulate Phenotype 1||-----||-----||-----||-----||-----#
  common.causal=get_common_causal(genodata=genodat)
  #Obtain a random common causal site (we can cange the bounds of the MAF. It is 
  #automatically set at (0.25, 0.35)).
  
  if(is.null(common.causal)){
    file.remove(filename)
    next
    #If our dataset does not have any SNP sites that can be the common causal variant,
    #then we just skip to next iteration.
  }
  

  
  #-----||-----||-----||-----||Simulate Phenotype 2||-----||-----||-----||-----||-----#
  rare.causal=get_rare_causals(genodata=genodat)
  #Obtain 10 random sites to act as the rare causal sites.
  
  if(is.null(rare.causal)){
    file.remove(filename)
    next
  }
  
  
  common_causal_vector=c(common_causal_vector, common.causal)
  #I'm just saving a vector containing all causal variants for each data.
  
  rare_causal_list=c(rare_causal_list, list(rare.causal))
  #I'm saving a list of the chosen rare causal variants for each dataset.
  
  hascausal=rowSums(genodat[,rare.causal])
  #This computes marginal counts of rare causal variants each individual has.
  
  
  
  ##Now, we simulate the data.
  Beta = seq(from=0.25, to = 1, by = 0.05)
  
  phen1=vector()
  phen2=vector()
  for(g in 1:length(Beta)){
    y1 = sim_pheno1(beta=Beta[g], common.causal = common.causal)
    y2 = sim_pheno2(beta =Beta[g], rare.causal = hascausal )
    
    phen1 = cbind(phen1, y1)
    phen2 = cbind(phen2, y2)
  }
  
  genodat=data.frame(phen1,phen2,genodat)
  #colnames(genodat)=c(rep(paste("Pheno", 1:16)),paste("V",1:segsites,sep=""))
  
  
  #-----||-----||Save our simulated phenotypes in a file||-----||-----||-----#
  table_name = paste("NoRecomb_PhenoAndGeno", j, ".txt", sep="")
  
  setwd("/global/home/hpc4300/Pilot_Study/Pilot_Study_PhenoAndGeno_Data")
  write.table(genodat, table_name, quote=F,row=F,col=F)
  
  
}


#Save our vector and list of causal variants:
setwd("/global/home/hpc4300/Pilot_Study/Pilot_Study_Results")
common_causal_name = paste("Actual_common_causal_vector_", task_id, ".txt", sep="")
write.table(common_causal_vector, common_causal_name, quote = F, row=F, col=F)


rare_table = matrix(nrow = length(rare_causal_list), ncol=10)
for (i in 1:length(rare_causal_list)){
  rare_table[i,] = rare_causal_list[[i]]
}

rare_causal_name = paste("Actual_rare_causal_table_", task_id, ".txt", sep="")
write.table(rare_table, rare_causal_name, quote = F, row=F, col=F)

