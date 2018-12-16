#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||    Step 2 - RCode preliminary codes    ||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#


#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----|| Goal: In step 1 of this project, we obtain genetic data     ||-----||-----#
#-----||-----||(genotype and gene tree data) using the ms program. The      ||-----||-----#
#-----||-----||output is saved in a single .txt file for each iteration...  ||-----||-----#
#-----||-----||This R script takes each .txt file and creates 2 new .txt    ||-----||-----#
#-----||-----||files that extracts and separates the haplotype data and the 
#-----||-----||gene tree data .
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

for (i in 1:3000){
  setwd("/global/home/hpc4300/BIM_Final_Raw_Data")
  #--> Set working directory to where we originally saved the ms program output.
  
  title = paste("results",i,".txt", sep="")
  data=readLines(title)
  #--> I've saved all ms outputs in the same name format. Ex: results1.txt
  
  treedata = data[5]
  #--> Treedata saved on the 5th line of the output.
  
  haplodata = data[8:length(data)]
  #--> Haplotype data is saved, beginning on the 8th line. There should be 200 lines
  # in total, representing 200 different haplotypes simulated.
  
  filename1 = paste("treedata",i, ".txt", sep="")
  filename2 = paste("haplodata",i, ".txt", sep="")
  
  setwd("/global/home/hpc4300/BIM_Final_Clean_Data")
  #--> Set working directory to where we want to save the extracted data.
  
  writeLines(treedata, filename1)
  writeLines(haplodata, filename2)
  
}

