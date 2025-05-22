#Get FST between genotype groups:
#Script to run plink and update files with genotype groupings and compare FST between them

#libraries
library(ggplot2)

#set directory
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/")

#read in clustering groups - individual info 
results<-read.table("PCADAPT/PCAdapt_results_can_nor_Ssa18region_with_clustering.txt", header=T)

results[which(results$Pop=="NRH"),]


#Create sample list of AA
all_AA <- results[which(results$Genotype_names=="AA"),]
write.table(all_AA[,c("Pop", "ID")], file = "PCADAPT/List_Individuals_AA_all.txt", col.names = F, row.names = F, sep="\t", quote = F)

all_BB <- results[which(results$Genotype_names=="BB"),]
write.table(all_BB[,c("Pop", "ID")], file = "PCADAPT/List_Individuals_BB_all.txt", col.names = F, row.names = F, sep="\t", quote = F)

all_CC <- results[which(results$Genotype_names=="CC"),]
write.table(all_CC[,c("Pop", "ID")], file = "PCADAPT/List_Individuals_CC_all.txt", col.names = F, row.names = F, sep="\t", quote = F)

#For AA - compare between Europe AA and Canada AA
all_AA_northAm <- results[which(results$Genotype_names=="AA" & results$Location=="NorthAm"),]
write.table(all_AA_northAm[,c("Pop", "ID")], file = "PCADAPT/List_Individuals_AA_NorthAm.txt", col.names = F, row.names = F, sep="\t", quote = F)

all_AA_europe <- results[which(results$Genotype_names=="AA" & results$Location=="Norway"),]
write.table(all_AA_europe[,c("Pop", "ID")], file = "PCADAPT/List_Individuals_AA_Europe.txt", col.names = F, row.names = F, sep="\t", quote = F)

#Alternative genotypes in each continent
all_CC_northAm <- results[which(results$Genotype_names=="CC" & results$Location=="NorthAm"),]
write.table(all_CC_northAm[,c("Pop", "ID")], file = "PCADAPT/List_Individuals_CC_NorthAm.txt", col.names = F, row.names = F, sep="\t", quote = F)

all_BB_europe <- results[which(results$Genotype_names=="BB" & results$Location=="Norway"),]
write.table(all_BB_europe[,c("Pop", "ID")], file = "PCADAPT/List_Individuals_BB_Europe.txt", col.names = F, row.names = F, sep="\t", quote = F)




#subset for EU vs Can for AA, and alternative genoypes - FOR WHOLE CHROMOSOME 18 (will use for selection analysis - SweeD) and for FST
#Make as VCF
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005 --chr 18 --chr-set 30 --keep ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/PCADAPT/List_Individuals_AA_Europe.txt --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_wholeChr_Europe --make-bed --recode-vcf")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005 --chr 18 --chr-set 30 --keep ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/PCADAPT/List_Individuals_AA_NorthAm.txt --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_wholeChr_NorthAm --make-bed --recode-vcf")

system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005 --chr 18 --chr-set 30 --keep ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/PCADAPT/List_Individuals_CC_NorthAm.txt --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CC_genotypes_only_Ssa18_wholeChr_NorthAm --make-bed --recode-vcf")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005 --chr 18 --chr-set 30 --keep ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/PCADAPT/List_Individuals_BB_Europe.txt --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_genotypes_only_Ssa18_wholeChr_Europe --make-bed --recode-vcf")


##To compare FST will need to update IDs to do FST between groups (AA vs. BB; AA vs. CC, CC vs. BB)


#Create BB updating codes
fam_file_BB <- read.table("Genomic_data/BB_genotypes_only_Ssa18_wholeChr_Europe.fam", header=F)
fam_file_BB$New <- rep("BB")
new_ID_codes_BB <- cbind(fam_file_BB[,1:2], fam_file_BB[,c(7,2)])
write.table(new_ID_codes_BB, file = "Genomic_data/update_BB_ids_all_Europe.txt", quote = F, row.names = F, col.names = F, sep="\t")

#Create CC updating codes
fam_file_CC <- read.table("Genomic_data/CC_genotypes_only_Ssa18_wholeChr_NorthAm.fam", header=F)
fam_file_CC$New <- rep("CC")
new_ID_codes_CC <- cbind(fam_file_CC[,1:2], fam_file_CC[,c(7,2)])
write.table(new_ID_codes_CC, file = "Genomic_data/update_CC_ids_all_NorthAm.txt", quote = F, row.names = F, col.names = F, sep="\t")

##Create AA - EU and AA - Can updating ID files
fam_file_AA_NorthAM <- read.table("Genomic_data/AA_genotypes_only_Ssa18_wholeChr_NorthAm.fam", header=F)
fam_file_AA_NorthAM$New <- rep("AACAN")
new_ID_codes_NorthAM_AA <- cbind(fam_file_AA_NorthAM[,1:2], fam_file_AA_NorthAM[,c(7,2)])
write.table(new_ID_codes_NorthAM_AA, file = "Genomic_data/update_AA_ids_CANADA.txt", quote = F, row.names = F, col.names = F, sep="\t")

fam_file_AA_Europe <- read.table("Genomic_data/AA_genotypes_only_Ssa18_wholeChr_Europe.fam", header=F)
fam_file_AA_Europe$New <- rep("AAEU")
new_ID_codes_Europe_AA <- cbind(fam_file_AA_Europe[,1:2], fam_file_AA_Europe[,c(7,2)])
write.table(new_ID_codes_Europe_AA, file = "Genomic_data/update_AA_ids_Europe.txt", quote = F, row.names = F, col.names = F, sep="\t")



#Update IDs for each group
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_genotypes_only_Ssa18_wholeChr_Europe --update-ids ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/update_BB_ids_all_Europe.txt --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_genotypes_only_Ssa18_region_update_IDs_Europe --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CC_genotypes_only_Ssa18_wholeChr_NorthAm --update-ids ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/update_CC_ids_all_NorthAm.txt --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CC_genotypes_only_Ssa18_region_update_IDs_NorthAm --make-bed")

#Update IDs for Can vs Eu for AA genotype
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_wholeChr_NorthAm --update-ids ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/update_AA_ids_CANADA.txt --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_region_update_IDs_CANADA --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_wholeChr_Europe --update-ids ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/update_AA_ids_Europe.txt --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_region_update_IDs_EUROPE --make-bed")

#Files will need to be merged for comparison (FST)
#Average FST values are used from these estimates - Note the number of fixed loci in comparisons is relevant - as FST excludes those loci - but in some cases they represent 50% of loci
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --fst --family --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_region_update_IDs_CANADA --bmerge ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_region_update_IDs_EUROPE --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_Canada_Vs_EU_coded --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --fst --family --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_region_update_IDs_EUROPE --bmerge ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_genotypes_only_Ssa18_region_update_IDs_Europe --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_vs_BB_genotypes_only_Ssa18_Europe --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --fst --family --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_region_update_IDs_CANADA --bmerge ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CC_genotypes_only_Ssa18_region_update_IDs_NorthAm --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_vs_CC_genotypes_only_Ssa18_NorthAmerica --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --fst --family --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_genotypes_only_Ssa18_region_update_IDs_Europe --bmerge ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CC_genotypes_only_Ssa18_region_update_IDs_NorthAm --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_vs_CC_genotypes_only_Ssa18_BtwnContinents --make-bed")

#combined all karyotypes and make raw file to use later
#system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --fst --family --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_genotypes_only_Ssa18_region_update_IDs --bmerge ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_vs_CC_genotypes_only_Ssa18_region_all --chr-set 30 --recodeA --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_BB__CC_genotypes_only_Ssa18_region_all --make-bed")

##Plotting FST between groups
#Read in AA vs AA between continents
FST_AA_EU_vs_Canada <- read.table("Genomic_data/AA_genotypes_only_Ssa18_Canada_Vs_EU_coded.fst", header=T)
FST_AA_EU_vs_Canada$FST2 <- FST_AA_EU_vs_Canada$FST

#Change NAs to 0 for plotting
FST_AA_EU_vs_Canada$FST2[is.na(FST_AA_EU_vs_Canada$FST)] <- 0
FST_AA_EU_vs_Canada$Group<-rep("AA vs AA: Between Continents")

#read in AA vs BB
FST_AA_vs_BB <- read.table("Genomic_data/AA_vs_BB_genotypes_only_Ssa18_Europe.fst", header=T)
FST_AA_vs_BB$FST2 <- FST_AA_vs_BB$FST
FST_AA_vs_BB$FST2[is.na(FST_AA_vs_BB$FST)] <- 0
FST_AA_vs_BB$Group<-rep("AA vs BB: Within Europe")

#Read in CC vs AA
FST_AA_vs_CC <- read.table("Genomic_data/AA_vs_CC_genotypes_only_Ssa18_NorthAmerica.fst", header=T)
FST_AA_vs_CC$FST2 <- FST_AA_vs_CC$FST
FST_AA_vs_CC$FST2[is.na(FST_AA_vs_CC$FST)] <- 0
FST_AA_vs_CC$Group<-rep("AA vs CC: Within North America")


### AA vs CC
nrow(FST_AA_vs_CC[which(FST_AA_vs_CC$FST==1),])/(nrow(FST_AA_vs_CC))
nrow(FST_AA_vs_CC[is.na(FST_AA_vs_CC$FST),])/(nrow(FST_AA_vs_CC))

## AA vs BB 
nrow(FST_AA_vs_BB[is.na(FST_AA_vs_BB$FST),])/(nrow(FST_AA_vs_BB))
nrow(FST_AA_vs_BB[which(FST_AA_vs_BB$FST==1),])/(nrow(FST_AA_vs_BB))



#Combine FST results
combined_FST_results<- rbind(FST_AA_EU_vs_Canada,FST_AA_vs_BB,FST_AA_vs_CC)

#save as file
write.table(combined_FST_results[,1:6], "FST/FST_between_karyotypes_groups_Updated_Feb17_2023.txt", quote = F, row.names = F, col.names = T, sep="\t")

#Plot of FST for different comparisons - save as pdf
str(combined_FST_results)
 
combined_FST_results$Group <- as.factor(as.character(combined_FST_results$Group))
plot_wholeChr<- ggplot(data=combined_FST_results)+
     geom_vline(xintercept = 50282126, lty=2)+
    geom_vline(xintercept = 52929573, lty=2)+
     geom_point(aes(x=POS, y=FST, col=Group), size=2, pch=21)+facet_wrap(.~Group, ncol=1)+
   ylab("Locus-specific FST")+xlab("Genomic position on Ssa18 (bp)")+
   scale_colour_manual(values=c("firebrick", "orange", "dodgerblue4"))+
   theme_bw()+theme(panel.grid = element_blank(), legend.position = "none",
                    strip.text.x = element_text(size = 15),
                    axis.title = element_text(size=15),
                    axis.text =  element_text(size=13))+
    geom_point(data=combined_FST_results[is.na(combined_FST_results$FST),],
               aes(x=POS, y=FST2, fill=Group), col="black",  size=4, pch=24)

    
plot_region<- ggplot(data=combined_FST_results)+
    geom_vline(xintercept = 50282126, lty=2)+
    geom_vline(xintercept = 52929573, lty=2)+
    geom_point(data=combined_FST_results[is.na(combined_FST_results$FST),],
               aes(x=POS, y=FST2, fill=Group), col="black",  size=4, pch=24)+

    geom_point(aes(x=POS, y=FST, col=Group), size=2, pch=21)+facet_wrap(.~Group, ncol=1)+
    ylab("Locus-specific FST")+xlab("Genomic position on Ssa18 (bp)")+
    scale_colour_manual(values=c("firebrick", "orange", "dodgerblue4"))+
    theme_bw()+theme(panel.grid = element_blank(), legend.position = "none",
                     strip.text.x = element_text(size = 15),
                     axis.title = element_text(size=15),
                     axis.text =  element_text(size=13))+
    xlim(49e6, 54e6)


 
cowplot::plot_grid(plot_wholeChr,plot_region, nrow = 1, rel_widths = c(2,1)) 

#Mann Whitney - should use FST = 0 or use NAs?
 combined_FST_results$Rearrangement <- rep("Outside")
 combined_FST_results$Rearrangement[which(combined_FST_results$POS > 50282126 &combined_FST_results$POS < 52929573)] <- rep("Inside")
 combined_FST_results$Rearrangement <- as.factor(as.character(combined_FST_results$Rearrangement))
 
 str(combined_FST_results)

nrow(combined_FST_results[which(combined_FST_results$Rearrangement=="Inside"  &  combined_FST_results$Group=="AA vs AA: Between Continents" ),])
 
64/106
36/106

#fixed for same alleles
aggregate(combined_FST_results$FST[which(combined_FST_results$Rearrangement=="Inside")], 
           list(combined_FST_results$Group[which(combined_FST_results$Rearrangement=="Inside")]), function(x) sum(is.na(x)))
 
aggregate(combined_FST_results$FST[which(combined_FST_results$Rearrangement=="Outside")], 
           list(combined_FST_results$Group[which(combined_FST_results$Rearrangement=="Outside")]), function(x) sum(is.na(x)))
 

aggregate(combined_FST_results$FST, list(combined_FST_results$Group, combined_FST_results$Rearrangement), FUN=mean, na.rm=T) 

#fixed for different alleles
nrow(combined_FST_results[which(combined_FST_results$Rearrangement=="Inside"  &  combined_FST_results$Group=="AA vs AA: Between Continents" &  combined_FST_results$FST==1),])
nrow(combined_FST_results[which(combined_FST_results$Rearrangement=="Inside"  &  combined_FST_results$Group=="AA vs BB: Within Europe" &  combined_FST_results$FST==1),])
nrow(combined_FST_results[which(combined_FST_results$Rearrangement=="Inside"  &  combined_FST_results$Group=="AA vs CC: Within North America" &  combined_FST_results$FST==1),])

76/106
26/106
 
0.988593849/ 0.007997335
0.813948736/0.017545815
 
AA_comparison_MannWhit<-wilcox.test(FST ~ Rearrangement, data=combined_FST_results[which(combined_FST_results$Group=="AA vs AA: Between Continents"),], na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
 print(AA_comparison_MannWhit)
 
 AA_CC_comparison_MannWhit<-wilcox.test(FST ~ Rearrangement, data=combined_FST_results[which(combined_FST_results$Group=="AA vs CC: Within North America"),], na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
 print(AA_CC_comparison_MannWhit)
 
 AA_BB_comparison_MannWhit<-wilcox.test(FST ~ Rearrangement, data=combined_FST_results[which(combined_FST_results$Group=="AA vs BB: Within Europe"),], na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
 print(AA_BB_comparison_MannWhit)
 
 
