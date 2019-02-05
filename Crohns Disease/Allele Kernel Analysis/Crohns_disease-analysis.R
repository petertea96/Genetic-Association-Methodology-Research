#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||       Crohn's Disease Data Analysis - Allele kernels only   ||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

.libPaths("/global/home/hpc4300/RPackages")
#Read in the data:
setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Results")
# Read in the processed genotype data
dat=read.table("Processed_haplodata.txt")

# Create the binary response variable (transmitted/not)
status=rep(c(1,0),nrow(dat)/2)




source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Crohn_Calculate_all_kernels.R")
source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Crohn_NEW_gtsm.R")
source("/global/project/hpcg1578/Crohn/Kernel_Analysis/Peter_R_Code/Crohn_MDMR_Code.R")




G=dat
#Extract genotype information.



n=nrow(dat)
#Extract population size information.

K=ncol(dat)


mykernels = get.kernels(G=G, n=n, K=K, treename=treename)  
#Calculate all kernel functions.


Phenotype_Results = data.frame(Method = character(),IBS = numeric(),
                               AM = numeric(), AS = numeric(), h1 = numeric(),
                               Skat = numeric(), stringsAsFactors = FALSE )

#Phenotype_Results = data.frame(Method = character(),IBS = numeric(),
#                                     AM = numeric(), AS = numeric(), h1 = numeric(),
#                                     Skat = numeric(),Tree1 = numeric(), Tree2 = numeric(),
#                                     Tree3 = numeric(), Tree4 = numeric(), Tree5 = numeric(), stringsAsFactors = FALSE )




#-----||-----||-----|| - SKAT Code - ||-----||-----||-----#
null.model = SKAT_Null_Model(status ~ 1, out_type="D")
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
p.val.gtsm = similarity.regression.pheno1(P1=status,n=n, kernel.list=mykernels)
#Obtain list of GTSR p-values

for (j in (1:length(p.val.gtsm))){
  #Fill in results table with GTSR p-values
  Phenotype_Results[2,j+1] = p.val.gtsm[j] 
}


#-----||-----||-----||-----|| - MDMR - ||-----||-----||-----||-----#
p.val.MDMR = pval_P1_MDMR_function(P1=status, kernel.list=mykernels)
#Obtain list of MDMR p-values

for (j in (1:length(p.val.MDMR))){
  #Fill in results table with MDMR p-values
  Phenotype_Results[3,j+1] = p.val.MDMR[j]
}

labels = c("SKAT", "GTSM", "MDMR")
for (i in 1:length(labels)){
  Phenotype_Results[i,1] = labels[i]
} 

write.table(Phenotype_Results,"Data_part1.txt",quote=F,row=F,col=F)


