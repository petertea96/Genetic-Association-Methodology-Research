#Takes the ms text results and extracts the gene tree data and the haplotype data into two separate files.

for (i in 1:2500){
  setwd("/global/home/hpc4300/BIM_Final_Raw_Data")
  title = paste("results",i,".txt", sep="")
  data=readLines(title)
  treedata = data[5]
  haplodata = data[8:length(data)]
  
  filename1 = paste("treedata",i, ".txt", sep="")
  filename2 = paste("haplodata",i, ".txt", sep="")
  
  setwd("/global/home/hpc4300/BIM_Final_Clean_Data")
  writeLines(treedata, filename1)
  writeLines(haplodata, filename2)
  
}

