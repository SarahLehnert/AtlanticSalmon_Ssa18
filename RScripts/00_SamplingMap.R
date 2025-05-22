#require(maptools)
library(maps)
library(mapproj)
library(mapdata)
library(rgeos)
library(maptools)
library(sp)
library(raster)
library(ggspatial)
library(rgdal)
library(grid)
library(lattice)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(openxlsx)
library(cowplot)

#This script was adapted from one shared with me from Gerald Chaput (DFO)
#Some packages may not work with newer versions of R

options(stringsAsFactors = FALSE)
setwd("c:/Users/LehnertS/Desktop/Smolt_Paper/figures/Rscript_map/")
#+++++++++++++++++++


# ggplot does not work with spatial objects like shape files or raster images. You need to convert them into data.frames.
# In case you're working with SpatialPointsDataFrame, this is as simple as:
#    mapdata <- data.frame(yourshapefile)
#    ggplot() +  geom_point( data= mapdata, aes(x=long, y=lat), color="red")
# In case you're working with SpatialPolygonsDataFrame, you need to fortify the object, like this:
#   yourshapefile_df <- fortify(yourshapefile, region ="id")
# now create the map
#  ggplot() + geom_point(data= yourshapefile_df, aes(x=long, y=lat, group=group), .....color="red"...
# or geom_polygon


#############################
# geospatial data were obtained from : https://www.naturalearthdata.com/downloads/
# a zip file is downloaded from this site when you select a feature
# copy all those files in a directory, for example "natural earth spatial files" 
# you then call them with the main file name that has multiple extensions using readOGR from the rgdal package
# I keep track of which spatial files are polygons and which are points because you use a different function to plot them

# skip over these lines because the Rdata file for eastern Canada is already created

tmpcoastal <- readOGR(dsn = "C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                       layer = "ne_10m_coastline")
  
coastal <- fortify(tmpcoastal)  # this is a spatial lines dataframe

#note that this layer states_provinces matches exactly the coastline layer
tmpprov <- readOGR(dsn = "C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles/10m_cultural", 
                    layer = "ne_10m_admin_1_states_provinces")  

prov2 <- fortify(tmpprov)  # this is a spatial polygon dataframe
 
tmpland <- readOGR(dsn = "C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                  layer = "ne_10m_land_scale_rank")
land2 <- fortify(tmpland)  # this is a spatial polygon dataframe

tmprivers <- readOGR(dsn ="C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                    layer = "ne_10m_rivers_lake_centerlines_scale_rank")
riv2 <- fortify(tmprivers)  # this is a spatial lines dataframe
  
tmprivers <- readOGR(dsn ="C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                      layer = "ne_10m_rivers_north_america")
rivnac2 <- fortify(tmprivers)  # this is a spatial lines dataframe
  
tmplakes <- readOGR(dsn ="C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                      layer = "ne_10m_lakes")  
lakes2 <- fortify(tmplakes)  # this is a spatial polygons dataframe
  
tmplakes <- readOGR(dsn ="C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                   layer = "ne_10m_lakes_historic")   # this is a spatial polygons dataframe
lakes3 <- fortify(tmplakes) # this is a spatial polygons dataframe
  
tmplakes <- readOGR(dsn ="C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                   layer = "ne_10m_lakes_north_america")   # this is a spatial polygons dataframe
lakesnac2 <- fortify(tmplakes) # this is a spatial polygons dataframe
 
# bathymetry line for 200 m
tmpbath <- readOGR(dsn ="C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles", 
                  layer = "ne_10m_bathymetry_K_200")  
bath200 <- fortify(tmpbath)  # this is a spatial polygon dataframe
 
#these are big files and they take time to plot as they contain information for the entire globe
#to streamline the plotting process, I selected the data corresponding to the area for eastern Canada
#keep only geographic area of interest between latitude 40 to 70, longitude -45 to -80
#you use the pipeline %>% and functions that are in package dplyr to filter longitude and latitude

prov_nor <- prov2 %>% filter(long <= 45, long >= -10, lat <= 80, lat >= 45)
 
land_nor <- land2 %>% filter(long <= 45, long >= -10, lat <= 80, lat >= 45)
riv_nor <- riv2 %>% filter(long <= 45, long >= -10, lat <= 80, lat >= 45)
rivnac_nor <- rivnac2 %>% filter(long <= 45, long >= -10, lat <= 80, lat >= 45)
lakes_nor <- lakes2 %>% filter(long <= 45, long >= -10,lat <= 80, lat >= 45)
lakeshist_nor <- lakes3 %>% filter(long <= 45, long >= -10, lat <= 80, lat >= 45)
lakesnac_nor <- lakesnac2 %>% filter(long <= 45, long >= -10, lat <= 80, lat >= 45)

# #For eastern canada
prov <- prov2 %>% filter(long <= -45, long >= -80, lat <= 70, lat >= 40)
land <- land2 %>% filter(long <= -45, long >= -80, lat <= 70, lat >= 40)
riv <- riv2 %>% filter(long <= -45, long >= -80, lat <= 70, lat >= 40)
rivnac <- rivnac2 %>% filter(long <= -45, long >= -80, lat <= 70, lat >= 40)
lakes <- lakes2 %>% filter(long <= -45, long >= -80, lat <= 70, lat >= 40)
lakeshist <- lakes3 %>% filter(long <= -45, long >= -80, lat <= 70, lat >= 40)
lakesnac <- lakesnac2 %>% filter(long <= -45, long >= -80, lat <= 70, lat >= 40)
 
# #save these as a list and save it as Rdata
geodata <- list("coastal" = coastal, "prov" = c(prov, prov_nor), "land" = c(land, land_nor), "riv" = c(riv, riv_nor), "rivnac" = c(rivnac, rivnac_nor), "lakes" = c(lakes, lakes_nor), "lakeshist" = c(lakeshist, lakeshist_nor), "lakesnac" = c(lakesnac, lakesnac_nor),
                 "bath200" = bath200)
# 
# 
save(geodata, file = "C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles/both_geodata_dataframe.Rdata")

# start from here to make maps (if running later)
# load the reduced geodata as dataframes
#load(file = "C:/Users/LehnertS/Desktop/COSEWIC/DU_Maps/rasterfiles/both_geodata_dataframe.Rdata")
#attach(geodata)  # this breaks up the list into its dataframes that you can then call by name

# DNA sampled rivers database file for plotting
dnaall <- read.table("Population_Clustering_PCA_by_6groups_Chr18_counts.txt", header = TRUE, stringsAsFactors=FALSE)
names(dnaall)
dnarivers <- dnaall %>% dplyr::select(Var1, Lat, Long, Continent) %>% distinct()  # select river name, lat and long, and keep only unique river names, some rivers have multiple locations


#######################


#editting for plotting Norway
minLong <- 2; maxLong <- 34; meanLong <- mean(c(minLong, maxLong))
minLat <- 56; maxLat <- 72; meanLat <- mean(c(minLat, maxLat))

# figure name to save in specified directory for each DU as a .tif file


pdf("Map_GeneticSamples_Ssa18.pdf",  height=15, width=8)

# gglplot mapping commands


testmap_norway <-   ggplot() +
  theme_bw() + 
  geom_polygon(data = land_nor, 
               aes(x=long, y = lat, group = group), fill = "gray95", colour="gray50") +
#geom_polygon(data = bath200, aes(x=long, y = lat, group = group), linetype = 2, fill = NA, color="grey30") +
geom_path(data=riv_nor, aes(x = long, y = lat, group = group), color="lightblue") +
geom_path(data=rivnac_nor, aes(x = long, y = lat, group = group), color = "lightblue") +
geom_polygon(data = lakes_nor, aes(x=long, y = lat, group = group), fill = "lightblue", color="lightblue") +
geom_polygon(data = lakesnac_nor, aes(x=long, y = lat, group = group), fill = "lightblue", color="lightblue") +
geom_polygon(data = lakeshist_nor, aes(x=long, y = lat, group = group), fill = "lightblue", color="lightblue") +
geom_polygon(data = prov_nor, aes(x=long, y = lat, group = group), fill = NA, color="gray50", size = 0.1) +
#geom_polygon(data = fedlands, aes(x=long, y = lat, group = group), fill = "lightgreen", color="lightgreen") +
  geom_point(data = dnarivers[which(dnarivers$Continent=="Norway"),], aes(x=Long, y = Lat), colour = "black", fill = "firebrick3", shape = 21, size =4, alpha = 0.95) +
  labs( x = "Longitude",  y = "Latitude") + 
  theme(
    panel.grid.major = element_line(colour = "grey75", linetype = 3),
    #legend.position = "top",
    legend.position = c(0.875, 0.84),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6,"line"),
    legend.background = element_rect(colour = NA),
    legend.margin = margin(c(0,0,0,0)),
    legend.direction = "vertical",
    legend.title = element_blank(),
    axis.text.y   = element_text(size=15),
    axis.text.x   = element_text(size=15),
    axis.title.y  = element_text(size=15),
    axis.title.x  = element_text(size=15),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)) +
  geom_text(aes(x=-5, y=70), label="Europe", size=10)+
  coord_map("ortho", xlim = c(minLong*1.01, maxLong*0.98),  ylim = c(minLat, maxLat), orientation = c(meanLat, meanLong, 0))



#editting for plotting
minLong_can <- -72; maxLong_can <- -52; meanLong_can <- mean(c(minLong_can, maxLong_can))
minLat_can <- 42; maxLat_can <- 59; meanLat_can <- mean(c(minLat_can, maxLat_can))


testmap_can <-   ggplot() +
  theme_bw() + 
  geom_polygon(data = land, 
               aes(x=long, y = lat, group = group), fill = "gray95", colour="gray50") +
  #geom_polygon(data = bath200, aes(x=long, y = lat, group = group), linetype = 2, fill = NA, color="grey30") +
  geom_path(data=riv, aes(x = long, y = lat, group = group), color="lightblue") +
 geom_path(data=rivnac, aes(x = long, y = lat, group = group), color = "lightblue") +
  geom_polygon(data = lakes, aes(x=long, y = lat, group = group), fill = "lightblue", color="lightblue") +
 geom_polygon(data = lakesnac, aes(x=long, y = lat, group = group), fill = "lightblue", color="lightblue") +
 geom_polygon(data = lakeshist, aes(x=long, y = lat, group = group), fill = "lightblue", color="lightblue") +
 geom_polygon(data = prov, aes(x=long, y = lat, group = group), fill = NA, color="gray50", size = 0.1) +
 #geom_polygon(data = fedlands, aes(x=long, y = lat, group = group), fill = "lightgreen", color="lightgreen") +
  geom_point(data = dnarivers[which(dnarivers$Continent=="NorthAm"),], aes(x=Long, y = Lat), colour = "black", fill = "firebrick3", shape = 21, size =4, alpha = 0.95) +
  labs( x = "Longitude",  y = "Latitude") + 
  theme(
    panel.grid.major = element_line(colour = "grey75", linetype = 3),
    #legend.position = "top",
    legend.position = c(0.875, 0.84),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6,"line"),
    legend.background = element_rect(colour = NA),
    legend.margin = margin(c(0,0,0,0)),
    legend.direction = "vertical",
    legend.title = element_blank(),
    axis.text.y   = element_text(size=15),
    axis.text.x   = element_text(size=15),
    axis.title.y  = element_text(size=15),
    axis.title.x  = element_text(size=15),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)) +
geom_text(aes(x=-54, y=58), label="North America", size=10)+
coord_map("ortho", xlim = c(minLong_can*1.01, maxLong_can*0.98),  ylim = c(minLat_can, maxLat_can), orientation = c(meanLat_can, meanLong_can, 0))


cowplot::plot_grid(testmap_can, testmap_norway, align = "hv", nrow = 2, labels = c("A", "B"), label_size = 36)

dev.off()
