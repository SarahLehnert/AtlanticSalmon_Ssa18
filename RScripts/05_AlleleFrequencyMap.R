#R script to plot pie charts of allele frequencies

library(maps) # tool for maps
library(mapdata) # all your basemaps are here
library(marmap) # for bathymetry if needed
library(mapplots) # for add.pie
library(gplots) # for colour range
library(rworldmap)
library(maptools)
library(lattice)


#Set wd
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/PCADAPT/")

#######################
allele_freq_results1 <- read.table("Population_Clustering_PCA_by_6groups_Chr18_counts.txt", header=T)

#get adjust lat/long for mapping
mapLocs <- read.table("Maps_AlleleFreq/OLD_Results_PCAdapt_frequency_with_location_info.txt", header=T)[,c(1,15:19)]

#revise locations for Laks and Lakj 
mapLocs[which(mapLocs$Var1=="Laks"), ]
mapLocs[which(mapLocs$Var1=="Lakj"), ]
#swap adj coordinates for these two locations
mapLocs$Lat_Adj_nor[which(mapLocs$Var1=="Laks") ] <- 71
mapLocs$Lat_Adj_nor[which(mapLocs$Var1=="Lakj") ] <- 70.06

mapLocs$Long_Adj_nor[which(mapLocs$Var1=="Laks") ] <- 23
mapLocs$Long_Adj_nor[which(mapLocs$Var1=="Lakj") ] <- 27.55

#Merge allele freq with locations
allele_freq_results<- merge(x=allele_freq_results1, y=mapLocs, by=1)


#Turn frequencies as numbers
allele_freq_results$FreqA <-  allele_freq_results$FreqA_allele*100
allele_freq_results$FreqB <-  allele_freq_results$Freq_B_allele*100
allele_freq_results$FreqC <-  allele_freq_results$Freq_C_allele*100

Sample.Lat.lim=c(42,57)
Sample.Long.lim=c(-71,-49)

map("worldHires", col="gray", bg="white", 
    xlim=Sample.Long.lim, xlab= "C", ylim=Sample.Lat.lim, fill=T, border=NA,
    lwd=0.000001, mar=rep(5,4),
    resolution=0);map.axes()
mtext(c("Longitude", "Latitude"), side=c(1,2), line = 2.5, cex=1.5)

##Add line segments to plot
#Edit coordinates for pie charts
head(allele_freq_results)

#Add segents for populations
for(i in 1:nrow(allele_freq_results))
{
  segments(
    x0=allele_freq_results$Adj_Long[i],y0=allele_freq_results$Adj_Lat[i],
    x1=allele_freq_results$Long[i],y1=allele_freq_results$Lat[i], lwd =2
    
  )
}


for(i in 1:nrow(allele_freq_results))
{
  add.pie(as.integer(allele_freq_results[i, c("FreqA", "FreqB", "FreqC")]),
          x=allele_freq_results$Adj_Long[i],y=allele_freq_results$Adj_Lat[i],labels="",radius = 0.35,
          col=c("dodgerblue4", "gold", "firebrick2"), border = T)
}
#text(y= allele_freq_results$Adj_Lat+0.5, x=allele_freq_results$Adj_Long, label=allele_freq_results$Var1)

rect(-62, 53, -59, 54.5,
     col=NULL, border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)
rect(-57, 46.5, -52, 48.5,
     col=NULL, border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)


legend("topright",
       c("Shared A", "European B", "North Am C"), 
       pch=rep(15,length(3)),
       col=c("dodgerblue4", "gold", "firebrick2"))


#Save as 14 x 14 PDF

##########NL only

Sample.Lat.lim=c(46.5,48.5)
Sample.Long.lim=c(-57,-52)


map("worldHires", col="gray", bg="white", 
    xlim=Sample.Long.lim, xlab= "C", ylim=Sample.Lat.lim, fill=T, border=NA,
    lwd=0.000001, mar=rep(5,4),
    resolution=0);map.axes()
#mtext(c("Longitude", "Latitude"), side=c(1,2), line = 2.5, cex=1.5)

##Add line segments to plot
#Edit coordinates for pie charts


head(allele_freq_results)

for(i in 1:nrow(allele_freq_results))
{
  segments(
    x0=allele_freq_results$Adj_Long[i],y0=allele_freq_results$Adj_Lat[i],
    x1=allele_freq_results$Long[i],y1=allele_freq_results$Lat[i], lwd =2
    
  )
}



for(i in 1:nrow(allele_freq_results))
{
  add.pie(as.integer(allele_freq_results[i, c("FreqA", "FreqB", "FreqC")]),
          x=allele_freq_results$Adj_Long[i],y=allele_freq_results$Adj_Lat[i],labels="",radius = 0.09,
          col=c("dodgerblue4", "gold", "firebrick2"), border = T)
}
#text(y= allele_freq_results$Adj_Lat+0.1, x=allele_freq_results$Adj_Long, label=allele_freq_results$Var1)


#legend("bottomright",
 #      c("Shared A", "European B", "North Am C"), 
  #     pch=rep(15,length(3)),
   #    col=c("dodgerblue4", "gold", "firebrick2"))

#Save as 10 x 6 

#text(allele_freq_results$Long,allele_freq_results$Lat,allele_freq_results$Code , cex=0.5)

##########



Sample.Lat.lim=c(53,54.5)
Sample.Long.lim=c(-62,-59)


map("worldHires", col="gray", bg="white", 
    xlim=Sample.Long.lim, xlab= "C", ylim=Sample.Lat.lim, fill=T, border=NA,
    lwd=0.000001, mar=rep(5,4),
    resolution=0);map.axes()
#mtext(c("Longitude", "Latitude"), side=c(1,2), line = 2.5, cex=1.5)

##Add line segments to plot
#Edit coordinates for pie charts


head(allele_freq_results)

for(i in 1:nrow(allele_freq_results))
{
  segments(
    x0=allele_freq_results$Adj_Long[i],y0=allele_freq_results$Adj_Lat[i],
    x1=allele_freq_results$Long[i],y1=allele_freq_results$Lat[i], lwd =2
    
  )
}


for(i in 1:nrow(allele_freq_results))
{
  add.pie(as.integer(allele_freq_results[i, c("FreqA", "FreqB", "FreqC")]),
          x=allele_freq_results$Adj_Long[i],y=allele_freq_results$Adj_Lat[i],labels="",radius = 0.08,
          col=c("dodgerblue4", "gold", "firebrick2"), border = T)
}

#text(y= allele_freq_results$Adj_Lat+0.1, x=allele_freq_results$Adj_Long, label=allele_freq_results$Var1)

#legend("bottomright",
 #      c("Shared A", "European B", "North Am C"), 
  #     pch=rep(15,length(3)),
   #    col=c("dodgerblue4", "gold", "firebrick2"))
#Save as 10 x8 pdf


##### Norway

Sample.Lat.lim=c(57,72)
Sample.Long.lim=c(0,35)



map("worldHires", col="gray", bg="white", 
    xlim=Sample.Long.lim, xlab= "C", ylim=Sample.Lat.lim, fill=T, border=NA,
    lwd=0.000001, mar=rep(5,4),
    resolution=0);map.axes()
mtext(c("Longitude", "Latitude"), side=c(1,2), line = 2.5, cex=1.5)

##Add line segments to plot
#Edit coordinates for pie charts


head(allele_freq_results)

norway_only <- allele_freq_results[which(allele_freq_results$Continent=="Norway"),]

for(i in 1:nrow(norway_only))
{
  segments(
    x0=norway_only$Long_Adj_nor[i],y0=norway_only$Lat_Adj_nor[i],
    x1=norway_only$Long[i],y1=norway_only$Lat[i], lwd =2
    
  )
}


for(i in 1:nrow(norway_only))
{
  add.pie(as.integer(norway_only[i, c("FreqA", "FreqB", "FreqC")]),
            x=norway_only$Long_Adj_nor[i],y=norway_only$Lat_Adj_nor[i],labels="",radius = 0.35,
          col=c("dodgerblue4", "gold", "firebrick2"), border = T)
}

#text(y= allele_freq_results$Adj_Lat+0.1, x=allele_freq_results$Adj_Long, label=allele_freq_results$Var1)

legend("bottomright",
       c("Shared A", "European B", "North Am C"), 
       pch=rep(15,length(3)),
       col=c("dodgerblue4", "gold", "firebrick2"))
#Save as 14 x 14 PDF

#Figures were combined in illustrator
