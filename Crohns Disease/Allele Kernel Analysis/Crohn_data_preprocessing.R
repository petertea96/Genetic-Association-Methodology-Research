#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||      Pre-process Crohn's data set      ||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Data")

## STEP 1: Read in the haplotype data. Must specify that it is type "character". 
haplodat = read.table("crohn5q31_haplo.dat", colClasses="character")

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


#-----||-----||Save our processed data||-----||-----||-----#
table_name = paste("Processed_haplodata", ".txt", sep="")

setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Results")
write.table(newhaplodat, table_name, quote=F,row=F,col=F)


