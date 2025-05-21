##Smolt age project
library(tidyr)
library(data.table)

#set directory
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/")

##Update Positions using new genome v3.1
newPos<-read.table("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/NewPositions_v3.1_220kArraySNPs/Ssa220k_vs_Simon3.1.txt", header=F)
newPos$Chr <- gsub(pattern = "ssa0",replacement = "",x = newPos$V3)
newPos$Chr <- gsub(pattern = "ssa",replacement = "",x = newPos$Chr)

#Write data for new map information (updated positions)
write.table(as.data.frame(newPos[,c(5,2)]), file = "New_Chr.txt", quote = F, row.names = F, col.names = F, sep="\t") #chromosome
write.table(as.data.frame(newPos[,c(4,2)]), file = "New_map.txt", quote = F, row.names = F, col.names = F, sep="\t") #genomic position

#update SNP data using plink with new chr and position names
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse --chr-set 30 --update-chr ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/New_Chr.txt 1 2 --update-map ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/New_map.txt 1 2 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse_updatedPos_v3 --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/HUNT_RIVER --chr-set 30 --update-chr ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/New_Chr.txt 1 2 --update-map ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/New_map.txt 1 2 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/HUNT_RIVER_updatedPos_v3 --make-bed")


#First need to subset 220K data for populations with available smolt age information
#Use plink to subset the file and filter for MAF>0.05

#Run throuhg plink command line here --  don't add any new lines - keep code on one line
#Note - for Margaree River and Restigouche - combine pop data (multiple tributaries) and get one allele freq value for the population (as only one smolt age value for the river system)
#to do this rename individuals UPS, KED, PAT, MAT <- All restigouche, MRS and MNE <- All Margaree

fam_ids<-read.table("CIGENE_2021_NAEU_forgrilse_updatedPos_v3.fam")

fam_ids$V1<-as.character(fam_ids$V1)
fam_ids$New <- fam_ids$V1

fam_ids$New[which(fam_ids$V1 == "UPS")] <- c("RES")
fam_ids$New[which(fam_ids$V1 == "KED")] <- c("RES")
fam_ids$New[which(fam_ids$V1 == "PAT")] <- c("RES")
fam_ids$New[which(fam_ids$V1 == "MAT")] <- c("RES")
fam_ids$New[which(fam_ids$V1 == "MRS")] <- c("MRG")
fam_ids$New[which(fam_ids$V1 == "MNE")] <- c("MRG")

#Save new IDs
write.table(fam_ids[, c(1,2,7,2)], file = "merge_MRG_RES_sample_names.txt", quote = F, row.names = F, col.names = F)
nrow(fam_ids)

#use plink to update ID names based on above info (merge HUNT file and full data file)
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse_updatedPos_v3 --chr-set 30 --update-ids ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/merge_MRG_RES_sample_names.txt --out  ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse_updatenames --make-bed")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse_updatenames --chr-set 30 -bmerge ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/HUNT_RIVER_updatedPos_v3 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse_updatenames_with_HUNT --make-bed")

#Use plink to filter and create subset data for different analyses. Get allele frequency data
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse_updatenames_with_HUNT  --keep-fam ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Life_history/Pop_List_CanadaSmoltAges_mergenames --make-bed --allow-extra-chr --maf 0.05 --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_filtered_smoltage_pops_maf005")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_2021_NAEU_forgrilse_updatenames_with_HUNT  --keep-fam ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Life_history/Pop_List_Norway_CanadaSmoltAges_mergenames --make-bed --allow-extra-chr --maf 0.05 --chr-set 30 --freq --family --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_filtered_smoltage_pops_RANGE_WIDE_maf005")
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_filtered_smoltage_pops_maf005  --freq --family --allow-extra-chr --chr-set 30 --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_filtered_smoltage_pops_maf005")


