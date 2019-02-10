#This script calculates the mean tree association statistics for each tree file 
# (each tree file had around 797 to 800 trees). There are 100 tree files.

#There are 3 Association statistics...
#There are 5 Kernel functions...

setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Results")


GTSR_results = matrix(ncol = 5, nrow = 100 )
MDMR_results = matrix(ncol = 5, nrow = 100 )
SKAT_results = matrix(ncol = 5, nrow = 100 )

for (i in (1:100)){

  GTSR_file_name = paste("GTSR_part", i, ".txt", sep="")

  if (!file.exists(GTSR_file_name)){
	print(i)
    next
  #This is just an added sanity check. To proceed, we need to ensure the file actually exists...
	}
  GTSR_my_table = read.table(GTSR_file_name)
  GTSR_results[i,] = colMeans(GTSR_my_table)
  
  MDMR_file_name = paste("MDMR_part", i, ".txt", sep="")
  MDMR_my_table = read.table(MDMR_file_name)
  MDMR_results[i,] = colMeans(MDMR_my_table)
  
print(i)
  SKAT_file_name = paste("SKAT_part", i, ".txt", sep="")
  SKAT_my_table = read.table(SKAT_file_name)
  SKAT_results[i,] = colMeans(SKAT_my_table)
  
}

setwd("/global/project/hpcg1578/Crohn/Kernel_Analysis/Mean_Results")
write.table(GTSR_results, "GTSR_mean_results.txt",quote=F,row=F,col=F)
write.table(MDMR_results, "MDMR_mean_results.txt",quote=F,row=F,col=F)
write.table(SKAT_results, "SKAT_mean_results.txt",quote=F,row=F,col=F)





