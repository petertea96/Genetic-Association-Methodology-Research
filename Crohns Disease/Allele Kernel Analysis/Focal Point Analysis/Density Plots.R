setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's project/Cluster files/Final codes/Genetic-Association-Methodology-Research/Crohns Disease/Allele Kernel Analysis/Focal Point Analysis")

focal_data = read.table("New_Allele_Kernel_Data.txt" )
View(focal_data)
colnames(focal_data) = c("Statistic", "IBS", "AM", "AS", "H1", "SKAT")

focal_data = focal_data[,-4] #Accidently added AS kernel...but we don' want to study AS
# - 3/27/2019

### ----- SKAT ----- ###
library(dplyr)
SKAT = filter(focal_data, Statistic == "SKAT")

library(reshape2)
Stacked_SKAT = melt(data=SKAT, id.vars = 1) 
Stacked_SKAT$x_column = rep(c(1:100),4)
Stacked_SKAT$y = -log10(Stacked_SKAT$value)

library(ggplot2)
ggplot(Stacked_SKAT, aes(x = x_column, y = y)) + facet_wrap(~variable) +
  geom_line(color="blue") + 
  geom_point(colour = "black", fill="white", size = 1)+
  #geom_vline(xintercept = 0.05, show.legend = TRUE, linetype = "dashed", color = "red") +
  ggtitle("SKAT - Genotype scoring kernel association distributions") +
  xlab("Focal Point") + ylab("-log10(P-value)") +
  theme_classic() +
  theme(strip.background = element_rect(fill="lightblue")) +
  scale_fill_discrete(name = "Kernel")



### ----- MDMR ----- ###
library(dplyr)
MDMR = filter(focal_data, Statistic == "MDMR")

library(reshape2)
Stacked_MDMR = melt(data=MDMR, id.vars = 1) 
Stacked_MDMR$x_column = rep(c(1:100),4)

### -- Replace p-values that are 0 to (1 / 501) -- ###
sum(Stacked_MDMR$value==0)
perm.p = 1/ 501
Stacked_MDMR$value[which(Stacked_MDMR$value==0)]
Stacked_MDMR$value[which(Stacked_MDMR$value==0)] = perm.p
Stacked_MDMR$value
### -- Replace p-values that are 0 to (1 / 501) -- ###

Stacked_MDMR$y = -log10(Stacked_MDMR$value)
library(ggplot2)

ggplot(Stacked_MDMR, aes(x = x_column, y = y)) + facet_wrap(~variable) +  ylim(0,3.5) +
  geom_line(color="blue") + 
  geom_point(colour = "black", fill="white", size = 1)+  
  #geom_vline(xintercept = 0.05, show.legend = TRUE, linetype = "dashed", color = "red") +
  ggtitle("MDMR - Genotype scoring kernel association distributions") +
  xlab("Focal Point") + ylab("-log10(P-value)") +
  theme_classic() +
  theme(strip.background = element_rect(fill="lightgreen")) +
  scale_fill_discrete(name = "Kernel")


### ----- GTSR or "SimReg"----- ###
library(dplyr)
GTSM = filter(focal_data, Statistic == "GTSM")

library(reshape2)
Stacked_GTSM = melt(data=GTSM, id.vars = 1) 
Stacked_GTSM$x_column = rep(c(1:100),4)
### -- Replace p-values that are 0 to (1 / 501) -- ###
sum(Stacked_GTSM$value==0)
perm.p = 1/ 501
Stacked_GTSM$value[which(Stacked_GTSM$value==0)]
Stacked_GTSM$value[which(Stacked_GTSM$value==0)] = perm.p
Stacked_GTSM$value
### -- Replace p-values that are 0 to (1 / 501) -- ###
Stacked_GTSM$y = -log10(Stacked_GTSM$value)

library(ggplot2)
ggplot(Stacked_GTSM, aes(x = x_column, y = y)) + facet_wrap(~variable) + ylim(0,3.5) +
  geom_line(color="blue") + 
  geom_point(colour = "black", fill="white", size = 1)+ 
  
  #geom_vline(xintercept = 0.05, show.legend = TRUE, linetype = "dashed", color = "red") +
  ggtitle("SimReg - Genotype scoring kernel association distributions") +
  xlab("Focal Point") + ylab("-log10(P-value)") +
  theme_classic() + 
  theme(strip.background = element_rect(fill="lightcoral")) +
  scale_fill_discrete(name = "Kernel")




