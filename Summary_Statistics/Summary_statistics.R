#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####        Calculate some summary statistics in the simulated data        #####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#Date:  March 2nd, 2019

#I set an aray job with 100 "arrays". I will split up the files to analyse into 100 chunks:
#There are 3000 files in total to potentially analyse, so I split this into 100 chunks with
#30 files belonging to each chunk.

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)
#Obtain Slurm Task ID.  


#Now, determine indices of data files to analyse:
total_files=seq(from=1, to= 3001, by=30)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id +1] - 1
#Compute ending index

Summary_statistics_table = matrix(ncol=4)
setwd("/global/home/hpc4300/BIM_Final_PhenoAndGeno_Data")
for (j in starting:ending){
  table_name = paste("NoRecomb_PhenoAndGeno", j, ".txt", sep="")
  
  if (!file.exists(table_name)){
    next
    #This is just an added sanity check. To proceed, we need to ensure the file actually exists...
  }
  
  
  
  dat=read.table(table_name, colClasses=c("numeric"))
  just_geno_dat = dat[,-c(1,2)]
  
  number_variant_sites = ncol(just_geno_dat)
  
  
  f1=colSums(just_geno_dat)/ (2*nrow(just_geno_dat))
  
  number_rare_variant = sum(f1<0.05)
  number_common_variant = sum(f1 < 0.35 & f1 > 0.2)
  
  row_to_add = c(j, number_variant_sites, number_rare_variant, number_common_variant)
  
  Summary_statistics_table = rbind(Summary_statistics_table, row_to_add)
  
} 

Summary_statistics_table = Summary_statistics_table[-1,]
colnames(Summary_statistics_table) = c("Data File", "Number poly. sites", "Number rare", "Number Common")


output_name = paste("SummaryResults_", task_id, ".txt", sep="")
write.table(Summary_statistics_table,output_name,quote=F,row=F,col=F)








