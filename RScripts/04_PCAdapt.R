
#PCAdapt script
library(pcadapt)
library(qvalue)
library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)
library(reshape)



#Subset region of interest for PCAdapt
#filter datasets 
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldat_maf005_NORWAY_only --chr-set 30 --chr 18 --from-bp 50000000 --to-bp 53000000 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldat_maf005_NORWAY_only_Ssa18_region_50000000_53000000 --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldat_maf005_CANADA_only --chr-set 30 --chr 18 --from-bp 50000000 --to-bp 53000000 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldat_maf005_CANADA_only_Ssa18_region_50000000_53000000 --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldat_maf005_NORWAY_only_Ssa18_region_50000000_53000000 --bmerge ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldat_maf005_CANADA_only_Ssa18_region_50000000_53000000 --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldat_maf005_CANADA_NORWAY_Ssa18_region_50000000_53000000 --make-bed")


setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/")


#Run PCAdapt with Norway and Canada together
nor_can_salmon_120pops <- read.pcadapt("Genomic_data/CIGENE_alldat_maf005_CANADA_NORWAY_Ssa18_region_50000000_53000000.bed",type="bed")

#Run with k=2
pcadapt_nor_can_salmon_2pc <- pcadapt(nor_can_salmon_120pops, method = "mahalanobis", K=2,min.maf = 0.05) #MAF already set in plink
plot(pcadapt_nor_can_salmon_2pc,option="screeplot") 

#Get percent variance explained by each axis
pcadapt_nor_can_salmon_2pc$singular.values^2
pct_pc <- round(((pcadapt_nor_can_salmon_2pc$singular.values)^2) / length(pcadapt_nor_can_salmon_2pc$maf),digits = 2)

#Read in loci info
file_info <- fread("Genomic_data/CIGENE_alldat_maf005_CANADA_NORWAY_Ssa18_region_50000000_53000000.fam", header=F)
head(file_info)

#Subset individual info
pop_list2 <- cbind(file_info$V1, file_info$V2)
pop_list2 <- as.data.frame(pop_list2)

#Combine pop info and PCA results
plot_data <- as.data.frame(cbind(pcadapt_nor_can_salmon_2pc$scores[,1:2],pop_list2))
head(plot_data)

#add column names
colnames(plot_data)=c("PC1", "PC2", "Pop", "ID")

#To add continent info - import data on location
#note updated Laks and Lakj coordinates - they were backwards before
continent <- read.csv("Locations/Locations_Norway_Can_Pops_140pops.csv", header=T)

#Merge site information with PCA data
plot_data_location<- merge(x=plot_data, y=continent[,1:5], by.x=3, by.y=1)
head(plot_data_location)
plot_data_location$Location <- plot_data_location$Continent

ggplot(data=plot_data_location, aes(x=PC1, y=PC2))+
  geom_point(aes(fill=Location),size=3, alpha=0.9, pch=21)+theme_bw()+
  ylab(paste0("PC2 (", round(pcadapt_nor_can_salmon_2pc$singular.values[2]^2, digits = 3)*100, "%)"))+
  xlab(paste0("PC1 (", round(pcadapt_nor_can_salmon_2pc$singular.values[1]^2, digits = 3)*100, "%)"))
  
##### Run this the first time - but remove if running script again - save results of this analysis#####  

#Kmeans clusters - This should be done once - as different results can be produced each time
#use Kmeans clusters to get 6 groups that fit the 6 clusters:
set.seed(7)
clusters <- kmeans(plot_data_location[,c("PC1", "PC2")], centers = 6)
plot_data_location$Kmeans <- clusters$cluster
plot_data_location$Kmeans <- as.factor(as.character(plot_data_location$Kmeans))


#Plot to see Kmeans cluster assignments
ggplot(data=plot_data_location, aes(x=PC1, y=PC2))+
geom_point(aes(fill=Kmeans),size=3, alpha=0.9, pch=21)+theme_bw()+
ylab(paste0("PC2 (", round(pcadapt_nor_can_salmon_2pc$singular.values[2]^2, digits = 3)*100, "%)"))+
xlab(paste0("PC1 (", round(pcadapt_nor_can_salmon_2pc$singular.values[1]^2, digits = 3)*100, "%)"))+
#geom_text(data=plot_data_location[which(plot_data_location$Pop=="Alta"),], aes(x=PC1, y=PC2,label=Pop))+
  theme_bw()
#levels(plot_data_location$Pop)
  
#Based on plotting - will allocate some 'names' to these groupings
plot_data_location$Genotype_names <- c("x")

table(plot_data_location$Kmeans, plot_data_location$Continent)

#This would always be different numbers for genotype class depending on how clusters are assigned in kmeans
plot_data_location$Genotype_names[which(plot_data_location$Kmeans==1)] <- "BB"
plot_data_location$Genotype_names[which(plot_data_location$Kmeans==2)] <- "CC"
plot_data_location$Genotype_names[which(plot_data_location$Kmeans==3)] <- "AA"
plot_data_location$Genotype_names[which(plot_data_location$Kmeans==4)] <- "BC"
plot_data_location$Genotype_names[which(plot_data_location$Kmeans==5)] <- "AC"
plot_data_location$Genotype_names[which(plot_data_location$Kmeans==6)] <- "AB"


#save the data to keep - and import it for future
#Remove a few columns - (site name had spaces in it)
write.table(plot_data_location[, c(1:6,9:11)], "PCADAPT/PCAdapt_results_can_nor_Ssa18region_with_clustering.txt", quote = F, row.names = F, col.names = T, sep="\t")

##### End of section - read in file produced to use here #####  


Pcadapt_full_results <- read.table("PCADAPT/PCAdapt_results_can_nor_Ssa18region_with_clustering.txt", header=T)
head(Pcadapt_full_results)

#Plot to see Kmeans cluster assignments
ggplot(data=Pcadapt_full_results, aes(x=PC1, y=PC2))+
  geom_point(aes(fill=as.factor(Genotype_names)),size=3, alpha=0.9, pch=21)+theme_bw()+
  scale_fill_manual(values=c("dodgerblue3","firebrick2", "pink", "orange", "black", "green"))+
  ylab(paste0("PC2 (", round(pcadapt_nor_can_salmon_2pc$singular.values[2]^2, digits = 3)*100, "%)"))+
  xlab(paste0("PC1 (", round(pcadapt_nor_can_salmon_2pc$singular.values[1]^2, digits = 3)*100, "%)"))


#Plot to see Results by Continent - Saved as 9x7 PDF
ggplot(data=Pcadapt_full_results, aes(x=PC1, y=PC2))+
  geom_point(aes(fill=Location),size=3, alpha=0.9, pch=21)+theme_bw()+
  scale_fill_manual(values=c("dodgerblue4","firebrick2"))+
  theme(panel.grid = element_blank(), axis.text = element_text(size=13), axis.title = element_text(size=16))+
  ylab(paste0("PC2 (", round(pcadapt_nor_can_salmon_2pc$singular.values[2]^2, digits = 3)*100, "%)"))+
  xlab(paste0("PC1 (", round(pcadapt_nor_can_salmon_2pc$singular.values[1]^2, digits = 3)*100, "%)"))


#Save R data - in case needed later - don't overwrite 


#Determine genotype proportion by continent 
counts_pops <- as.data.frame(table(Pcadapt_full_results$Pop, Pcadapt_full_results$Genotype_names))
head(counts_pops)

#Reshape data frame - wide
count_pops_tally <- reshape(counts_pops, idvar=c("Var1"), timevar="Var2", direction="wide")
head(count_pops_tally)

#Merge counts with pop info 
count_pops_tally_info <- as.data.frame(merge(x=count_pops_tally,by.x=1, y=continent[,c(1,2,3,5)], by.y=1))

#Count total number of individuals:
count_pops_tally_info$N_total <-count_pops_tally_info$Freq.AA + count_pops_tally_info$Freq.AB + count_pops_tally_info$Freq.AC + count_pops_tally_info$Freq.BB + count_pops_tally_info$Freq.BC + count_pops_tally_info$Freq.CC

#Get frequency of A allele;
count_pops_tally_info$FreqA_allele <- (((count_pops_tally_info$Freq.AA)*2) +(count_pops_tally_info$Freq.AB) +(count_pops_tally_info$Freq.AC)) / (count_pops_tally_info$N_total*2)

count_pops_tally_info$Freq_B_allele <- (((count_pops_tally_info$Freq.BB)*2) +(count_pops_tally_info$Freq.AB) +(count_pops_tally_info$Freq.BC)) / (count_pops_tally_info$N_total*2)

count_pops_tally_info$Freq_C_allele <- (((count_pops_tally_info$Freq.CC)*2) +(count_pops_tally_info$Freq.AC) +(count_pops_tally_info$Freq.BC)) / (count_pops_tally_info$N_total*2)



#Save file
write.table(count_pops_tally_info, "PCADAPT/Population_Clustering_PCA_by_6groups_Chr18_counts.txt", quote = F, row.names = F, col.names = T, sep="\t" )


#Plot some maps of the frequencies

