##Extracting environmental data

setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/")
library(raster)
library(rbioclim)
library(rgdal)
library(raster)
library(data.table)
library(tidyr)
library(psych)

#Get allele frequency data
all_snps_allelefreq <- read.table("PCADAPT/Population_Clustering_PCA_by_6groups_Chr18_counts.txt", header=T)


#Get bioclim data - for present
bioclim <- recursive.getData(times="pres") # by default resolution 2.5m and the 19 bioclim variables. (but can be changed)
present_bioclim <- bioclim[["pres"]]

#Extract data for Bio 1-19
bio_present <- present_bioclim[[paste0("bio",c(1:19))]]

##Site locations
all_loc <- read.csv("Locations/Locations_Norway_Can_Pops_140pops.csv", header=T)
head(all_loc)

#Merge site locations with SNP data (include a few things for site location info lat/long)
all_data <- merge(all_loc[,c(1:3,5)], as.data.frame(all_snps_allelefreq[,c(1,11:14)]),  by=1)
nrow(all_data)
all_data[1:10,1:4]



#change ARO -- location does not appear to be on land where Bioclim data is present
all_data$Long[which(all_data$Code=="ARO")] <- -66.93775
all_data$Lat[which(all_data$Code=="ARO")] <- 50.07073

#change LSR -- almost in ocean
all_data$Lat[which(all_data$Code=="LSR")] <- all_data$Lat[which(all_data$Code=="LSR")]+0.01


#Create coordinate dataframe
coordinates <- cbind(all_data$Long, all_data$Lat)

#Get bioclim data for site coordinates
bioclim_all_sites_present <- raster::extract(bio_present, coordinates)
#Combine bioclim data with Site codes
bioclim_all_sites_present_info <- as.data.frame(cbind(as.character(all_data$Code), bioclim_all_sites_present))

#Check for NAs -- may need to edit coordinates if any not on land
sum(is.na(bioclim_all_sites_present_info))

#merge with pop info
bioclim_all_sites_present_info2 <- merge(all_snps_allelefreq,bioclim_all_sites_present_info, by=1 )

#Save data for all sites and Bioclim 1-19
write.table(bioclim_all_sites_present_info2, "Environmental/BioclimData_all_Sites_Present_nor_can.txt", quote=F, row.names = F, col.names = T, sep="\t")


