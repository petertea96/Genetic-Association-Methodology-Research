#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#Simulate the phenotypes, then re-format our data with phenotypes.
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

common_causal_vector = vector()
rare_causal_list = list()
#For precaution, I will save which SNP sites is the common causal and which ones are
#rare causal. This will be compared to when I run the Single Locus Tests.


source("BIM_Rcode_Simulation_help.R")
for (j in (1:2500)){
  
  ## Read in the haplotype data. Must specify that it is type "character". 
  filename = paste("haplodata",j, ".txt", sep="")
  
  if (!file.exists(filename)){
    next
  #This is just an added sanity check. To proceed, we need to ensure the file actually exists...
    }
  
  
  haplodat=read.table(filename, colClasses=c("character"))
  
  #The first line of the haplotype data will be used to calculate the number of segregation sites,
  #otherwise known as the number of SNP sites.
  FirstLine = readLines(filename)[1]
  FirstLine=unlist(strsplit(FirstLine,split=""))
  segsites=length(FirstLine)
  
  ## Note that the loci are not separated by spaces. So we must split them in order to be able to manipulate the data. 
  ## Also convert to type numeric so that we can add allele counts
  newhaplodat=matrix(as.numeric(unlist(strsplit(haplodat[,1],split=""))),ncol=segsites,byrow=T)
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  #What if columns are the exact same? (I.E 100% correlated between SNPs?) Let's remove them.
  
  newhaplodat = unique.matrix(t(newhaplodat)) 
  #Returns a matrix with only unique ROWS (so I need the transpose)
  
  newhaplodat = t(newhaplodat)
  #*** If this section was too confusing, I've added a more in-depth explanation in my ReadME file.
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
  
  
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
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
  
  common_causal_vector=c(common_causal_vector, common.causal)
  #I'm just saving a vector containing all causal variants for each data.
  
  
  
  #-----||-----||-----||-----||Simulate Phenotype 2||-----||-----||-----||-----||-----#
  rare.causal=get_rare_causals(genodata=genodat)
  #Obtain 10 random sites to act as the rare causal sites.
  
  if(is.null(rare.causal)){
    file.remove(filename)
    next
  }
  
  rare_causal_list=c(rare_causal_list, list(rare.causal))
  #I'm saving a list of the chosen rare causal variants for each dataset.
  
  hascausal=rowSums(genodat[,rare.causal])
  #This computes marginal counts of rare causal variants each individual has.
  
  
  
  ##Now, we simulate the data.
  y1 = sim_pheno1(beta=0.45, common.causal = common.causal)
  y2 = sim_pheno2(beta =0.8, rare.causal = hascausal )
  
  genodat=data.frame(y1,y2,genodat)
  colnames(genodat)=c("Pheno1","Pheno2",paste("V",1:segsites,sep=""))

  
  #-----||-----||Save our simulated phenotypes in a file||-----||-----||-----#
  table_name = paste("NoRecomb_PhenoAndGeno", j, ".txt", sep="")
  write.table(genodat, table_name, quote=F,row=F,col=F)
  
  
}


#Save our vector and list of causal variants:
write.table(common_causal_vector, "common_causal_vector.txt", quote = F, row=F, col=F)


rare_table = matrix(nrow = length(rare_causal_list), ncol=10)
for (i in 1:length(rare_causal_list)){
  rare_table[i,] = rare_causal_list[[i]]
}

write.table(rare_table, "rare_causal_table.txt", quote = F, row=F, col=F)

