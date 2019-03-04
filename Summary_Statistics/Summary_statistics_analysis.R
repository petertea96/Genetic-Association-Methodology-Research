#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####           Analysis of summary statistics in the simulated data        #####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#####-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#####
#Date:  March 3rd, 2019
setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's project/Cluster files/Final codes/Genetic-Association-Methodology-Research/")

summary_data = read.table("Aggregate_Summary.txt")
#View(summary_data)
colnames(summary_data) =c("Data File", "Number poly. sites", "Number rare", "Number Common")

write.table(summary_data,"Aggregate_Summary.txt")

colMeans(summary_data)
nrow(summary_data)
