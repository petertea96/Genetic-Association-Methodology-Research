#Takes the ms text results and extracts the gene tree data and the haplotype data into two separate files.
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)
#Obtain Slurm Task ID.  


#Now, determine indices of data files to analyse:
total_files=seq(from=1, to= 10001, by=100)

starting = total_files[task_id]
#Compute starting index

ending = total_files[task_id + 1] - 1
#Compute ending index
for (i in starting:ending){
  setwd("/global/home/hpc4300/Pilot_Study/Pilot_Study_Raw_Data/")
  title = paste("results",i,".txt", sep="")
  data=readLines(title)
  treedata = data[5]
  haplodata = data[8:length(data)]
  
  filename1 = paste("treedata",i, ".txt", sep="")
  filename2 = paste("haplodata",i, ".txt", sep="")
  
  setwd("/global/home/hpc4300/Pilot_Study/Pilot_Study_Clean_Data")
  writeLines(treedata, filename1)
  writeLines(haplodata, filename2)
  
}
