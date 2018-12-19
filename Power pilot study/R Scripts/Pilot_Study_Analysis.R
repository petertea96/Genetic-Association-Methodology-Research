#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----|| Pilot Study Data Analysis||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#


#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----|| Episode 1: Determine the proportion of times the actual     ||-----||-----#
#-----||-----|| causal site matches with the most statistically significant ||-----||-----#
#-----||-----|| SNP site chosen for the SLT test. Compare this proportion   ||-----||-----#
#-----||-----|| across all beta values used in phenotype simulation.        ||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

#Read in the data:

setwd("C:/Users/Peter/Documents/Uottawa/2018 - 2019 Honour's project/Cluster files/Final codes/Genetic-Association-Methodology-Research/Power pilot study/Final Results")

Actual_common_causal = read.table("Aggregate_Actual_common_causal_vector.txt")
Actual_rare_causal = read.table("Aggregate_Actual_rare_causal_table.txt")


Chosen_common_causal = scan("Aggregate_Chosen_common_causal.txt")
df1 <- matrix(Chosen_common_causal, nrow=length(Chosen_common_causal)/16, ncol=16, byrow=TRUE)


Chosen_rare_causal = scan("Aggregate_Chosen_rare_causal.txt")
df2 <- matrix(Chosen_rare_causal, nrow=length(Chosen_rare_causal)/16, ncol=16, byrow=TRUE)


Beta = seq(from=0.25, to = 1, by = 0.05)
mynames = as.character(Beta)
colnames(df1) = mynames
colnames(df2) = mynames

#All data should have the same number of rows
nrow(Actual_common_causal)
nrow(Actual_rare_causal)
nrow(df1)
nrow(df2)

n = nrow(df2)

#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----|| Goal: Count the Number of times the chosen common    ||-----||-----||-----#
#-----||-----||-----||  causal SNP matches the Actual common causal SNP.    ||-----||-----||-----#
#-----||-----||-----||  Do this across all 16 Beta values.                  ||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

#Let's create a matrix of True and Falses to see if the causal alleles match.
P1_results=matrix(nrow=1, ncol=16)
for (i in 1:n){
  real_common_causal = Actual_common_causal[i,1]
  to_test = df1[i,]
  results = (to_test == real_common_causal)
  P1_results = rbind(P1_results, results)
  
}

#First row is full of NAs (just due to the nature of how I initialized the matrix)
P1_results = P1_results[-1,]

#Check now
head(P1_results)

#Change Row names
row_name = paste("Iteration", 1:n, sep = " ")
rownames(P1_results) = row_name

common_causal_success = colSums(P1_results)/n
names(common_causal_success) = mynames



#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----|| Goal: Count the Number of times the chosen rare      ||-----||-----||-----#
#-----||-----||-----||  causal SNP matches one of the 10 Actual common      ||-----||-----||-----#
#-----||-----||-----||  causal SNP. Do this across all 16 Beta values.      ||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
P2_results=matrix(nrow=1, ncol=16)
for (i in 1:n){
  real_actual_causal_vector = Actual_rare_causal[i,]
  to_test = df2[i,]
  
  #returns a TRUE element if the chosen rare causal is in the vector of actual rare causals.
  results = (to_test %in% real_actual_causal_vector)
  P2_results = rbind(P2_results, results)
}

#First row is full of NAs (just due to the nature of how I initialized the matrix)
P2_results = P2_results[-1,]

#Checks:
head(P2_results) 
nrow(P2_results)


rownames(P2_results) = row_name

rare_causal_success = colSums(P2_results)/n
names(rare_causal_success) = mynames
rare_causal_success


#####
# --> Produce plots
common_causal_success
rare_causal_success
myresults_part1 = cbind(Beta, common_causal_success, rare_causal_success)



library(reshape2)
Melted_myresults_part1 = melt(myresults_part1, id.var='Beta')
colnames(Melted_myresults_part1) = c("Beta","Causal_model", "True_detection_rate" )
Melted_myresults_part1 = Melted_myresults_part1[-c(1:16),]

library(ggplot2)
myplot1 = ggplot(data = Melted_myresults_part1,
       aes(x = Beta, y = True_detection_rate, col=Causal_model)) + 
  geom_line(size = 1) + 
  scale_color_manual(values=c("#483D8B", "#F08080"))+
  labs(x= "Beta",
       y = "Proportion of true causal detection",
       title = "Comparison of Single Locus Test true causal variant detection",
       caption = "Peter Tea")+
  theme_bw() +
  theme(legend.title = element_text(size=10, 
                                     face="bold"),
        legend.text = element_text(size=10, 
                                   face="bold"),
    legend.position = c(0.8, 0.2)) 

myplot1

#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----|| Episode 2: Count the proportion of SLT p-values that are    ||-----||-----#
#-----||-----||  statistically significant (i.e. p-value below 0.05).       ||-----||-----#
#-----||-----||  Do this for both phenotypes across all different beta      ||-----||-----#
#-----||-----||  values.                                                    ||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#
#-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||-----#

#Phenotype 1 Results:

#First, read in the data...

P1_pvalues = read.table("Aggregate_Pheno1Results.txt")

#Remove NA values...
P1_pvalues = P1_pvalues[complete.cases(P1_pvalues),]
nrow(P1_pvalues)

#Create matrix that checks for each element if the p-value is statisticcally significant.
P1_pvalues_which_significant = (P1_pvalues[,-c(1,2)] < 0.05) 
rownames(P1_pvalues_which_significant) = row_name
colnames(P1_pvalues_which_significant) = mynames

P1_success = colSums(P1_pvalues_which_significant)/n
names(P1_success) = mynames
P1_success



#Phenotype 2 Results:

#First, read in the data...

P2_pvalues = read.table("Aggregate_Pheno2Results.txt")

#Remove NA values...
P2_pvalues = P2_pvalues[complete.cases(P2_pvalues),]
nrow(P2_pvalues)

#Create matrix that checks for each element if the p-value is statisticcally significant.
P2_pvalues_which_significant = (P2_pvalues[,-c(1,2)] < 0.05) 
rownames(P2_pvalues_which_significant) = row_name
colnames(P2_pvalues_which_significant) = mynames

P2_success = colSums(P2_pvalues_which_significant)/n
names(P2_success) = mynames
P2_success


###
# --> Let's create some plots

myresults_part2 = cbind(Beta, P1_success, P2_success)



library(reshape2)
Melted_myresults_part2 = melt(myresults_part2, id.var='Beta')
colnames(Melted_myresults_part2) = c("Beta","Phenotype model", "Sensitivity" )
Melted_myresults_part2 = Melted_myresults_part2[-c(1:16),]

library(ggplot2)
myplot2 = ggplot(data = Melted_myresults_part2,
                 aes(x = Beta, y = Sensitivity, col= Phenotype_model)) + 
  geom_line(size = 1) + 
  scale_color_manual(values=c("#483D8B", "#F08080"))+
  labs(x= "Beta",
       y = "Sensitivity",
       title = "Comparison of Sensitivities of the Single Locus Test under varying phenotype-genotype association strengths",
       caption = "Peter Tea")+
  theme_bw() +
  theme(legend.title = element_text(size=10, 
                                    face="bold"),
        legend.text = element_text(size=10, 
                                   face="bold"),
        legend.position = c(0.82, 0.18)) 

myplot2

