#Takes the ms text results and extracts the gene tree data and the haplotype data into two separate files.

for (i in 1:2500){
  title = paste("results",i,".txt", sep="")
  data=readLines(title)
  treedata = data[5]
  haplodata = data[8:length(data)]
  
  filename1 = paste("treedata",i, ".txt", sep="")
  filename2 = paste("haplodata",i, ".txt", sep="")
  writeLines(treedata, filename1)
  writeLines(haplodata, filename2)
  
}

