#Combining smolt age data and enviro data for plotting

#set diretory
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/")

library(dplyr)

#read in data
enviro <- read.table("Environmental/BioclimData_all_Sites_Present_nor_can.txt", header=T)
head(enviro)

#By continent
enviro_nor <- enviro[which(enviro$Continent=="Norway"),]
enviro_can <- enviro[which(enviro$Continent=="NorthAm"),]


#Relationship between latitude and allele frequency A
nor_results<-summary(glm(FreqA_allele~Lat, data=enviro_nor))
#r2 calculation for glm
1-(nor_results$deviance/nor_results$null.deviance)

#For canada
can_results <- summary(glm(FreqA_allele~Lat, data=enviro_can))
1-(can_results$deviance/can_results$null.deviance) #r2


##Smolt age data:
smoltage <- read.csv("Life_history/SmoltAge_by_Pop_updatedApril2022.csv")
dim(enviro)

#Norway smolt age data
nor_smoltage <- read.table("Life_history/Phenotype_Norway_Smolts.txt", header=T)
Id_info_nor<- read.table("Genomic_data/CIGENE_alldat_maf005_NORWAY_only.fam", header=F)
#remove extract character string
Id_info_nor$names <- gsub('.{4}$', '', Id_info_nor$V2)

#Merge data for norway and sum by river
combined_Nor_smoltdat<- merge(y=nor_smoltage, x=Id_info_nor, by.x="names", by.y=1)

#Read in genotype info too
geno_nor <- read.table("PCADAPT/PCAdapt_results_can_nor_Ssa18region_with_clustering.txt", header=T)
geno_with_smoltage_nor <- merge(x=combined_Nor_smoltdat, by.x="V2", y=geno_nor, by.y="ID")
write.table(geno_with_smoltage_nor, "Life_history/Geno_by_Individual_with_Phenotype_Norway.txt", quote = F, row.names = F, col.names = T, sep="\t")


#get info by river
smoltage_nor <- combined_Nor_smoltdat %>% 
   group_by(V1) %>%
   summarise(mean = mean(FWAgeC, na.rm=T), n = n())

smoltage_enviro_Norway <- merge(as.data.frame(smoltage_nor), enviro_nor, by=1)

#Not sure if this is the best way to do this? - but need to combine data for tributaries of Margaree and Restigouche to get one value for River system (only 1 smolt age value)
restigouche_combined <- as.data.frame(rbind(c("Var1"="RES", 
  (colSums(enviro[which(enviro$Var1=="PAT" | enviro$Var1=="KED" | enviro$Var1=="MAT" | enviro$Var1=="UPS"  ), 2:7])),
  (colMeans(enviro[which(enviro$Var1=="PAT" | enviro$Var1=="KED" | enviro$Var1=="MAT" | enviro$Var1=="UPS"  ), 8:9])),
  "Continent" = "NorthAm", 
  "N_total" =  (sum(enviro[which(enviro$Var1=="PAT" | enviro$Var1=="KED" | enviro$Var1=="MAT" | enviro$Var1=="UPS"  ), "N_total"])),
  (colMeans(enviro[which(enviro$Var1=="PAT" | enviro$Var1=="KED" | enviro$Var1=="MAT" | enviro$Var1=="UPS"  ), 12:33])))))

margaree_combined <- as.data.frame(rbind(c("Var1"="MRG", 
                                              (colSums(enviro[which(enviro$Var1=="MRS" | enviro$Var1=="MNE"), 2:7])),
                                              (colMeans(enviro[which(enviro$Var1=="MRS" | enviro$Var1=="MNE"), 8:9])),
                                              "Continent" = "NorthAm", 
                                              "N_total" =  (sum(enviro[which(enviro$Var1=="MRS" | enviro$Var1=="MNE" ), "N_total"])),
                                             (colMeans(enviro[which(enviro$Var1=="MRS" | enviro$Var1=="MNE" ), 12:33])))))

#Combine RES and Marg data to summarize smolt age info
to_replace <- rbind(margaree_combined, restigouche_combined)
colnames(to_replace) <- colnames(enviro)

#Change to numeric
to_replace$Freq.AA <- as.numeric(as.character(to_replace$Freq.AA))
to_replace$Freq.AB <- as.numeric(as.character(to_replace$Freq.AB))
to_replace$Freq.AC <- as.numeric(as.character(to_replace$Freq.AC))
to_replace$Freq.BB <- as.numeric(as.character(to_replace$Freq.BB))
to_replace$Freq.CC <- as.numeric(as.character(to_replace$Freq.CC))
to_replace$Freq.BC <- as.numeric(as.character(to_replace$Freq.BC))
to_replace$N_total <- as.numeric(as.character(to_replace$N_total))



#Uppdate freq based on allele counts (not average)
to_replace$FreqA_allele <- (((to_replace$Freq.AA)*2) + (to_replace$Freq.AB)+(to_replace$Freq.AC) ) / (to_replace$N_total*2)
to_replace$Freq_B_allele <- (((to_replace$Freq.BB)*2)+(to_replace$Freq.AB)+(to_replace$Freq.BC) )/ (to_replace$N_total*2)
to_replace$Freq_C_allele <- (((to_replace$Freq.CC)*2)+(to_replace$Freq.AC)+(to_replace$Freq.BC) )/ (to_replace$N_total*2)
                                                          
#Add in new MRG and RES info for combining
all_dat_env<- rbind(enviro, to_replace)

#All data to compare smolt age and enviro data for Canada
combine_smolt_age_enviro <- merge(y=all_dat_env, x=smoltage, by=1, all.x=T)


#Save data and perform plotting/analyses in new script
write.table(combine_smolt_age_enviro, "Environmental/SmoltAge_Enviro_canada.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(smoltage_enviro_Norway, "Environmental/SmoltAge_Enviro_norway.txt", quote = F, row.names = F, col.names = T, sep="\t")




