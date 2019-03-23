setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's project/Cluster files/Final codes/Genetic-Association-Methodology-Research/Crohns Disease/Allele Kernel Analysis/Focal Point Analysis")

focal_data = read.table("New_Allele_Kernel_Data.txt" )
View(focal_data)
colnames(focal_data) = c("Statistic", "IBS", "AM", "AS", "H1", "SKAT")


### ----- SKAT ----- ###
library(dplyr)
SKAT = filter(focal_data, Statistic == "SKAT")

library(reshape2)
Stacked_SKAT = melt(data=SKAT, id.vars = 1) 

library(ggplot2)
ggplot(Stacked_SKAT, aes(x=value, y = ..scaled.., fill=variable)) + geom_density() + facet_wrap(~variable) +
  ggtitle("SKAT - Allele/Genotype kernel association distributions") +
  xlab("Allele/Genotype kernel") + ylab("Statistic Density") +
  theme_classic() +
  theme(strip.background = element_rect(fill="lightblue")) 



### ----- MDMR ----- ###
library(dplyr)
MDMR = filter(focal_data, Statistic == "MDMR")

library(reshape2)
Stacked_MDMR = melt(data=MDMR, id.vars = 1) 

library(ggplot2)
ggplot(Stacked_MDMR, aes(x=value, y = ..scaled.., fill=variable)) + geom_density() + facet_wrap(~variable) +
  ggtitle("MDMR - Allele/Genotype kernel association distributions") +
  xlab("Allele/Genotype kernel") + ylab("Statistic Density") +
  theme_classic() +
  theme(strip.background = element_rect(fill="lightgreen")) 


### ----- GTSR ----- ###
library(dplyr)
GTSM = filter(focal_data, Statistic == "GTSM")

library(reshape2)
Stacked_GTSM = melt(data=GTSM, id.vars = 1) 

library(ggplot2)
ggplot(Stacked_GTSM, aes(x=value, y = ..scaled.., fill=variable)) + geom_density() + facet_wrap(~variable) +
  ggtitle("GTSM - Allele/Genotype kernel association distributions") +
  xlab("Allele/Genotype kernel") + ylab("Statistic Density") +
  theme_classic() +
  theme(strip.background = element_rect(fill="lightgreen")) 




