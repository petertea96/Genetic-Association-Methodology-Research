#This script produces boxplots of tree association statistics....

setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's project/Cluster files/Final codes/Genetic-Association-Methodology-Research/Crohns Disease/Tree Kernel Analysis")

#Data Processing steps...
GTSR_data = read.table("GTSR_mean_results.txt" ) 
MDMR_data = read.table("MDMR_mean_results.txt")
SKAT_data = read.table("SKAT_mean_results.txt")


my_column_names= c("Tree1", "Tree2", "Tree3", "Tree4", "Tree5" )

colnames(GTSR_data) = my_column_names
colnames(MDMR_data) = my_column_names
colnames(SKAT_data) = my_column_names

library(reshape)

GTSR_reshaped = melt(na.omit(GTSR_data))
MDMR_reshaped = melt(na.omit(MDMR_data))
SKAT_reshaped = melt(na.omit(SKAT_data))

#####-----------------------------------------------------------------------------------#####
library(ggplot2)
ggplot(GTSR_reshaped, aes(x=variable, y = value, fill=variable)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  ggtitle("SimReg - Tree Kernel association statistic distributions") +
  xlab("Tree kernel") + ylab("Statistic value") +
  theme_classic()


ggplot(MDMR_reshaped, aes(x=variable, y = value, fill=variable)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  ggtitle("MDMR - Tree Kernel association statistic distributions") +
  xlab("Tree kernel") + ylab("Statistic value") +
  theme_classic()


ggplot(SKAT_reshaped, aes(x=variable, y = value, fill=variable)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  ggtitle("SKAT - Tree Kernel association statistic distributions") +
  xlab("Tree kernel") + ylab("Statistic value") +
  theme_classic()













