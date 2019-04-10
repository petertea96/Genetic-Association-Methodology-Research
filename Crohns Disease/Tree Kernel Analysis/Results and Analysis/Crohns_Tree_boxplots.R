#This script produces boxplots of tree association statistics....

setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's project/Cluster files/Final codes/Genetic-Association-Methodology-Research/Crohns Disease/Tree Kernel Analysis/Results and Analysis")

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
GTSR_reshaped$x_column = rep(c(1:96),5) 
#x_column is a variable denoting the focal points

MDMR_reshaped = melt(na.omit(MDMR_data))
MDMR_reshaped$x_column = rep(c(1:96),5)

SKAT_reshaped = melt(na.omit(SKAT_data))
SKAT_reshaped$x_column = rep(c(1:96),5)

#####-----------------------------------------------------------------------------------#####
library(ggplot2)



ggplot(GTSR_reshaped, aes(x=x_column, y = value))+ facet_wrap(~variable) +
  #geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  geom_line(color="blue") + 
  geom_point(colour = "black", fill="white", size = 1)+
  ggtitle("SimReg - Tree Kernel association statistic distributions") +
  xlab("Focal Point") + ylab("Statistic value") +
  theme_classic() + 
  theme(strip.background = element_rect(fill="lightcoral")) +
  scale_fill_discrete(name = "Kernel")

ggplot(MDMR_reshaped, aes(x=x_column, y = value))+ facet_wrap(~variable) +
  geom_line(color="blue") + 
  geom_point(colour = "black", fill="white", size = 1)+
  ggtitle("MDMR - Tree Kernel association statistic distributions") +
  xlab("Focal Point") + ylab("Statistic value") +
  theme_classic() + 
  theme(strip.background = element_rect(fill="lightgreen")) +
  scale_fill_discrete(name = "Kernel")


ggplot(SKAT_reshaped, aes(x=x_column, y = value))+ facet_wrap(~variable) +
  geom_line(color="blue") + 
  geom_point(colour = "black", fill="white", size = 1)+ 
  ggtitle("SKAT - Tree Kernel association statistic distributions") +
  xlab("Focal Point") + ylab("Statistic value") +
  theme_classic() + 
  theme(strip.background = element_rect(fill="lightblue")) +
  scale_fill_discrete(name = "Kernel")













