#R script for running lostruct package (local PCAs)
#Note I created a function (located at end of script) to run lostruct 

#open libraries
library(lostruct)
library(ggplot2)
library(cowplot)
library(adegenet)
library(qqman)


#Set working directory
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/lostruc/")

#Need vcf.gz file
#convert to vcf in Plink then vcf.gz - see code below
system("cd ~/Desktop/Software/plink_mac_20200219/; ./plink --bfile ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005 --chr-set 30  --recode-vcf --out ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005_vcf")
system("bgzip -c  ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005_vcf.vcf >  ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005_vcf.vcf.gz")
system("tabix -p vcf ~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005_vcf.vcf.gz")

#Read in vcf - takes a while - so may be easier to import and save data in case you want to run again. 
#This is for North America only
genome220 <- read_vcf("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005_vcf.vcf.gz")
save.image("Genome220K_loaded_vcf.RData")

#If saving and opening later, can skip above steps and load .RData
#load("Genome220K_loaded_vcf.RData")

#The VCF file was read in with this function above 
#groups and IDs for individuals (for making plots with groups -- this assumes only two groups but could edit below to include more)

#Set IDs and Groups for plotting PCAs later
IDs <- data.table::fread("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005.fam") #Fam from plink for group assignments
group <- as.character(IDs$V6) #Groupings for Canada Vs. Europe (for PCA plotting) - file already had groupings here.

#Make map from VCF data
map <- data.table::fread("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CIGENE_alldata_2021_NAEU_with_HUNT_maf005_vcf.vcf", header=F, skip=35)[,1:3]
map[1:10,]
map2 <- cbind(map[,c(1,3)], rep(0), map[,2])
map2[1:10,]
nrow(map2)
map2 <- as.data.frame(map2)
colnames(map2) <- c("V1", "V2", "V3", "V4")

#My input
vcf <- genome220 #220K in loaded data (see read_vcf above)
dim(vcf)
map <- map2 #map in loaded data
chr_num="17" #try ssa08
cores=40 
window=100
k_num=2
groupings=group #groupings for PCA plotting
Inds=IDs$V2 #vector of individual IDs to add to PCA plotting


##See script at bottom for running lostruc - developed this to have quick script for running individual chromosome data
#Run function at bottom first - and then run code here.

setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/lostruc/50snp_window")

#I wrote the function below 'PCA_lostruc_genome_function' to run the lostruct script for one chromosome at a time (see function below at end of script)
# I probably should have wrote this as a loop instead - but I didnt :)
#Information includes a map object, vcf object, chr #, number of cores to use, window size (# SNPs), number of K to run, groupings, and individuals
#Not all of this is needed for lostruct - some info is used to create plots (see function at end of script for details).

#Run lostruc on each chromosome
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="1", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="2", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="3", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="4", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="5", cores=40, window=50, k_num=2, groupings=group, Inds=Inds) 
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="6", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="7", cores=40, window=50, k_num=2, groupings=group, Inds=Inds) 
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="8", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="9", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="10", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="11", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="12", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="13", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="14", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="15", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="16", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="17", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="18", cores=40, window=50, k_num=2, groupings=group, Inds=Inds) 
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="19", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="20", cores=40, window=50, k_num=2, groupings=group, Inds=Inds) 
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="21", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="22", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="23", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="24", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="25", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="26", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="27", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="28", cores=40, window=50, k_num=2, groupings=group, Inds=Inds) 
PCA_lostruc_genome_function(map=map2, vcf=vcf, chr_num="29", cores=40, window=50, k_num=2, groupings=group, Inds=Inds)



#Next step is to combine all output data (note the output file names in the fucntion have "220K_PCA.txt" so its using those files)

files <- list.files()
dbf.files <- files[grep("220K_PCA.txt", files, fixed=T)]

df_total = data.frame()

#combine files for all chromosomes
for(i in 1:29){
  tt <- read.table(dbf.files[i], header=T)
  df_total <- rbind(df_total,tt)
  rm(tt)  
}


#all results - save file - includes MDS1 and MDS2 values for all windows
write.table(df_total,file = "results_all_chr_lostruc_220k_msd_50snp_window.txt", quote = F, row.names = F, col.names = T, sep="\t")

#Read in results
all_data<-read.table("results_all_chr_lostruc_220k_msd_50snp_window.txt", header = T)
head(all_data)

#create new column with unique window names
all_data$SNP <- c(paste0("Window_",  1:nrow(all_data)))

#Outliers function - get outliers in the dataset - 
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#Outliers based on SD form mean
outlier_MSD1 <- outliers(all_data$MDS1,3) #Using 3 standard deviations from mean
outlier_MSD2 <- outliers(all_data$MDS2,3) #Using 3 standard deviations from mean

#Save outlier windows for MDS 1
write.table(all_data[which(abs(all_data$MDS1) >= min(abs(outlier_MSD1))), ],file = "OutlierWindow_50SNP_lostrcut_220K.txt",
            quote = F, row.names = F, col.names = T, sep="\t")

#Save outlier windows for MDS 2
write.table(all_data[which(abs(all_data$MDS2) >= min(abs(outlier_MSD2))), ],file = "OutlierWindow_50SNP_lostrcut_220K_MDS2.txt",
            quote = F, row.names = F, col.names = T, sep="\t")

#Plot manhattan of results
manhattan(x = all_data, chr ="Chr", bp = "MeanWindow_BP", p = "MDS1", snp = "SNP", logp = F, ylim=c(-1,1), ylab="MDS1",
          col = c('gray30', "gray70"), highlight = all_data$SNP[which(abs(all_data$MDS1) >= min(abs(outlier_MSD1))) ])

manhattan(x = all_data, chr ="Chr", bp = "MeanWindow_BP", p = "MDS2", snp = "SNP", logp = F, ylim=c(-1,1), ylab="MDS2",
          col = c('gray30', "gray70"), highlight = all_data$SNP[which(abs(all_data$MDS2) >= min(abs(outlier_MSD2))) ])



############################################################################################################################################################################################################################


##Here is a quick script for subsetting the vcf data for a specific region - and running PCA in adegenet  

#combine vcf data with map info
VCF_data_with_info=cbind(map, vcf)

#read in ID names and combine with grouping info - in case want to use any groupigns in the PCA (not needed) 
ind_names<- read.table("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/All_Genomic_Data_Norway_vs_Canada_twopops_maf005.fam", header=F)
IDs$order=c(1:nrow(IDs))
combined_names<-merge(x=IDs,y=ind_names, by.x="V2", by.y="V2",sort = F)

#note these groupings were already in the file 
groupings <- combined_names$V1.y

#Get first and last SNP in an outlier window -eg for Ssa18 here - extracted from VCF info
min_forpeak_translocation=which(VCF_data_with_info$V4 == 36027179 & VCF_data_with_info$V1 == 18) 
max_forpeak_translocation=which(VCF_data_with_info$V4 == 43424204 & VCF_data_with_info$V1 == 18) 

#Create genind data for PCA (adegenet) and run PCA
chr_matrix_highestMDS_trans <- as.genind(t(vcf[min_forpeak_translocation:max_forpeak_translocation,]))
geno_chr_peakMDS_trans <- tab(chr_matrix_highestMDS_trans, NA.method="mean")
salmon_pca_chr_highMDS_tras<-dudi.pca(df=geno_chr_peakMDS_trans,center = T, scale = F,nf = 2, scannf = F)

#Plot PCA results for outlier window
#This plot shows two groupings - may need to edit
ggplot(data=salmon_pca_chr_highMDS_tras$li, aes(salmon_pca_chr_highMDS_tras$li[,1], salmon_pca_chr_highMDS_tras$li[,2], col=groupings))+
  geom_vline(xintercept=0, col="gray")+geom_hline(yintercept=0, col="gray")+
  geom_point(size=2, alpha=0.2)+
  geom_text(aes(label=as.character(Inds)), size=3)+
  #  ylim(-5,5)+xlim(-20, 20)+
  scale_color_manual(values = c("red", "blue"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top",
        text = element_text(size=10),title = element_text(size=10), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlab("PC Axis 1")+ ylab("PC Axis 2")




####################################################################################################################################################################################



###Function to run lostruct on a chromosome - gives results (text file with MDS1 and MDS2 for all windows, and outlier windows for both coordinates (in this case based on 2 sd from mean))
#Also gives plots including one with MDS1 vs MDS2, MDS1 vs pos, MDS2 vs pos
#Other plots of PCAs for highest and lowest MDS2 peaks to show local PCA in outlier regions
#Need lostruct, ggplot2, cowplot, adegenet, data.table

#Can edit script if you don't want all this output - I thought it was useful at first - but may be too much info - creates a bunch of files
#Also - individuals and groupings are indicated here - this is just used for plotting after running lostruct - so can include so that you can view your groups in the output - or edit script to exclude this part

PCA_lostruc_genome_function <- function(map, vcf, chr_num, cores, window, k_num, groupings, Inds ) {
  #Map file (plink format)
  #vcf = use read_vcf first to read it in (see above)
  #chr_num = code for chromosome in map/vcf
  #cores = number of cores for lostruct
  #window = number of SNPs for windows
  #k_num = number of clusters for pca
  #groupings = vector of grouping value for all individuals
  #Inds =  vector of individual IDs
  
  #Combine VCF data with map information
  VCF_data_with_info=cbind(map, vcf)
  
  ##Positions map - for chromosome of interest
  Chrom_map=map[which(map$V1==chr_num),]
  
  #Chrom_map #Check number of rows for dataset
  
  #Get row in matrix for this chromosome 
  row_chr=which(VCF_data_with_info$V1 == chr_num) 
  length(row_chr)
  
  #Pick rows for this chromosome from Map file
  map[c(row_chr[1]:row_chr[length(row_chr)]),]
  
  message("Starting lostruc to run local PCAs")
  message(paste0("Number of SNPs on ", chr_num, " is ", length(row_chr)))
  
  #Run pca for VCF matrix with only rows of interest (for chromosome)
  VCF_for_Chr=vcf[c(row_chr[1]:row_chr[length(row_chr)]),] #Select genotypes for chromosomes
  
  pcs_ssa_for_chr <- eigen_windows(VCF_for_Chr,  k=k_num, win = window, mc.cores = cores) #Run local PCAs
  pcdist_ssa_for_chr <- pc_dist(pcs_ssa_for_chr, npc=k_num, mc.cores = cores) #Get PC distance
  
  message("
        Finished lostruc to run local PCAs
        ")
  missing_data_inVCF=sum(is.na(VCF_for_Chr)) / ( (sum(!is.na(VCF_for_Chr))) + (sum(is.na(VCF_for_Chr))) )
  message(paste0("
               Proportion of missing data in VCF matrix ", round(missing_data_inVCF, digits = 4), "  
               ---> NOTE: if too high for a locus or individual may not work below
               "))
  
  #Get mean position within windows for Chromosome (for windows of 100 SNPs)
  n=length(row_chr)/window
  output_map<- matrix(ncol=3, nrow=n)
  for(i in 1:(n)){
    y=(window*i)
    x=(window*i)-(window-1)
    output_map[i,]<-cbind(mean(c(Chrom_map[x,4], Chrom_map[y,4])), Chrom_map[x,4], Chrom_map[y,4])
  }
  
  
  #output_map has position information
  
  #Plot PCA of coord 1 and coord 2  (MDS1 vs MDS2) - regions that are outside of center represent potential regions of interest (outliers) (i.e., local PCA differs from rest of chromosome)
  message("Getting MDS1 and MDS2 values for results")
  message("If fails here, there are NAs in the pc matrix --> probably need to filter sites/individuals for missing data")
  
  #check missing
  missing_data=sum(is.na(pcdist_ssa_for_chr))
  message("Missing data in matrix (NAs) ---- PROBLEM!")
  
  fit2d_Chrom<- cmdscale(pcdist_ssa_for_chr, eig=TRUE, k=k_num)
  
  #Basic plots (ggplots below)
  message("Quick plots")
  plot(fit2d_Chrom$points,  xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(pcdist_ssa_for_chr)), pch=19)
  plot(cbind(output_map[,1], fit2d_Chrom$points[,1]),  xlab=paste0("Position on ", chr_num," (bp)" ), ylab="Coordinate 1", col=rainbow(1.2*nrow(pcdist_ssa_for_chr)), pch=19)
  plot(cbind(output_map[,1], fit2d_Chrom$points[,2]), xlab=paste0("Position on ", chr_num," (bp)" ), ylab="Coordinate 2", col=rainbow(1.2*nrow(pcdist_ssa_for_chr)), pch=19)
  
  
  #Create data frame of Chr results merge MDS1 and MDS2 results with Position information
  Chrom_combined_results_pca=as.data.frame(cbind(output_map[,], fit2d_Chrom$point))
  colnames(Chrom_combined_results_pca)=c("MeanWindow_BP", "Start", "Stop", "MDS1", "MDS2")
  Chrom_combined_results_pca$Chr=rep(chr_num)
  
  #Save results of Chr as text file
  message("Save data for MDS1 and MDS2 with window information")
  write.table(Chrom_combined_results_pca, file= paste0(chr_num, "_220K_PCA.txt"), quote=F, row.names = F, col.names = T)
  
  #Plot nice figures
  #Open save data file
  chrom_results=read.table(file = paste0(chr_num, "_220K_PCA.txt"), header=T)
  #head(chrom_results)
  
  #Set row names
  rownames(chrom_results)<-chrom_results$MeanWindow_BP
  
  #Get outlier values (more than 2 sd from mean) - using Forrester outlier function from her webpage
  chrom_mds1_outliers=chrom_results[which(chrom_results$MDS1 > mean(chrom_results$MDS1)+2*sd(chrom_results$MDS1) | chrom_results$MDS1 < mean(chrom_results$MDS1)+-2*sd(chrom_results$MDS1)),]
  chrom_mds2_outliers=chrom_results[which(chrom_results$MDS2 > mean(chrom_results$MDS2)+2*sd(chrom_results$MDS2) | chrom_results$MDS2 < mean(chrom_results$MDS2)+-2*sd(chrom_results$MDS2)),]
  
  #If want values close to zero
  chrom_zeros=chrom_results[which(chrom_results$MDS2 > -0.01  & chrom_results$MDS1 > -0.01 & chrom_results$MDS2 < 0.01  & chrom_results$MDS1 < 0.01 ),]
  
  #Save outliers MDS1 and MDS2
  message("Saving outliers (+/- 2 SD from mean) for MDS1 and MDS2 with window information")
  
  write.table(chrom_mds1_outliers, file= paste0(chr_num, "_220K_PCA_MDS1_outliers.txt"), quote=F, row.names = F, col.names = T)
  write.table(chrom_mds2_outliers, file= paste0(chr_num, "_220K_PCA_MDS2_outliers.txt"), quote=F, row.names = F, col.names = T)
  
  #Create plots for PCA (MDS1 vs MDS2), and MDS1 vs position, and MDS2 vs position (MDS2 seems to be the primary interst)
  message("Creating and saving plots")
  
  pca_chr=ggplot()+geom_point(data=chrom_results, aes(x=MDS1, y=MDS2))+
    ylab("MDS2")+xlab("MDS1")+
    geom_point(data=chrom_mds1_outliers, aes(x=MDS1, y=MDS2), col="purple" )+
    geom_point(data=chrom_mds2_outliers, aes(x=MDS1, y=MDS2), col="red" )+
    theme_bw()+theme(panel.grid = element_blank())
  
  mds1_chr=ggplot()+geom_point(data=chrom_results, aes(x=MeanWindow_BP, y=MDS1))+
    geom_smooth(data=chrom_results, aes(x=MeanWindow_BP, y=MDS1), span=length(row_chr)/100000 )+
    geom_point(data=chrom_mds1_outliers, aes(x=MeanWindow_BP, y=MDS1), col="purple" )+
    geom_point(data=chrom_mds2_outliers, aes(x=MeanWindow_BP, y=MDS1), col="red" )+
    ylab("MDS1")+xlab(paste0("Window position ", chr_num," (bp)"))+
    theme_bw()+theme(panel.grid = element_blank())
  
  mds2_chr=ggplot()+geom_point(data=chrom_results, aes(x=MeanWindow_BP, y=MDS2))+
    geom_smooth(data=chrom_results, aes(x=MeanWindow_BP, y=MDS2), span=0.1)+
    geom_point(data=chrom_mds1_outliers, aes(x=MeanWindow_BP, y=MDS2), col="purple" )+
    geom_point(data=chrom_mds2_outliers, aes(x=MeanWindow_BP, y=MDS2), col="red" )+
    ylab("MDS2")+xlab(paste0("Window position ", chr_num," (bp)"))+
    theme_bw()+theme(panel.grid = element_blank())
  
  
  #combine plots and save
  ggsave(plot_grid(pca_chr,mds1_chr , mds2_chr, nrow=1), file = paste0(chr_num, "_220K_PCA_Plots.pdf"), height=4, width=10)
  
  
  
  # ##Run PCA (adegenet) for peaks in MDS2 (low and high peaks)
  message("Running PCA (adegenet) and plotting results for lowest and highest MDS2 region... almost done")
  # 
  # #Get lowest peak region
  min_mds1=chrom_mds2_outliers[which.min(chrom_mds2_outliers$MDS2),]
  #Get first and last SNP in that low peak window
  min=which(VCF_data_with_info$V4 ==  min_mds1$Start & VCF_data_with_info$V1 == chr_num)
  max=which(VCF_data_with_info$V4 == min_mds1$Stop & VCF_data_with_info$V1 == chr_num)
  #
  #Create genind data for PCA (adegenet) and run PCA
  chr_matrix_lowestMDS=as.genind(t(vcf[min:max,]))
  geno_chr_lowMDS <- tab(chr_matrix_lowestMDS, NA.method="mean")
  salmon_pca_chr_lowMDS<-dudi.pca(df=geno_chr_lowMDS,center = T, scale = F,nf = 2, scannf = F)
  #
  # #Plot PCA results for lowest MDS2 peak
  #
  lowest_MDS2=ggplot(data=salmon_pca_chr_lowMDS$li, aes(salmon_pca_chr_lowMDS$li[,1], salmon_pca_chr_lowMDS$li[,2]))+
    geom_vline(xintercept=0, col="gray")+geom_hline(yintercept=0, col="gray")+
    geom_point(size=2, alpha=0.2)+ggtitle(paste0(chr_num, " (Lowest region for MDS2)"))+
    #    geom_text(aes(label=as.character(Inds)), size=3)+
    #  ylim(-5,5)+xlim(-20, 20)+
    scale_color_manual(values = c("red", "blue"))+theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top",
          text = element_text(size=10),title = element_text(size=10),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlab("PC Axis 1")+ ylab("PC Axis 2")
  
  ggsave(lowest_MDS2, file = paste0(chr_num, "_220K_PCA_Plots_lowestMDS2.pdf"), height = 12, width = 12)
  #
  # ##Do same thing but for highest MDS2
  #
  # #Get highest peak region
  max_mds1=chrom_mds2_outliers[which.max(chrom_mds2_outliers$MDS2),]
  # #Get first and last SNP in that high peak window
  min_forpeak=which(VCF_data_with_info$V4 ==  max_mds1$Start & VCF_data_with_info$V1 == chr_num)
  max_forpeak=which(VCF_data_with_info$V4 == max_mds1$Stop & VCF_data_with_info$V1 == chr_num)
  #
  # #Create genind data for PCA (adegenet) and run PCA
  chr_matrix_highestMDS=as.genind(t(vcf[min_forpeak:max_forpeak,]))
  geno_chr_peakMDS <- tab(chr_matrix_highestMDS, NA.method="mean")
  salmon_pca_chr_highMDS<-dudi.pca(df=geno_chr_peakMDS,center = T, scale = F,nf = 2, scannf = F)
  #
  # #Plot PCA results for highest MDS2 peak
  highest_MDS2=ggplot(data=salmon_pca_chr_highMDS$li, aes(salmon_pca_chr_highMDS$li[,1], salmon_pca_chr_highMDS$li[,2]))+
    geom_vline(xintercept=0, col="gray")+geom_hline(yintercept=0, col="gray")+
    geom_point(size=2, alpha=0.2)+ggtitle(paste0(chr_num, " (Highest region for MDS2)"))+
    #   geom_text(aes(label=as.character(Inds)), size=3)+
    #  #  ylim(-5,5)+xlim(-20, 20)+
    scale_color_manual(values = c("red", "blue"))+theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top",
          text = element_text(size=10),title = element_text(size=10),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlab("PC Axis 1")+ ylab("PC Axis 2")
  
  ggsave(highest_MDS2, file = paste0(chr_num, "_220K_PCA_Plots_highestMDS2.pdf"), height = 12, width = 12)
  
  #End of function
  
}
