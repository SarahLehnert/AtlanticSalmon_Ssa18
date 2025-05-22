##SweeD analysis of selection

#Run Ssa18 chromosome for European AA and BB individuals seperately
#ran files trhough plink (in previous filter script - and added --recode-vcf to make these files as vcf format for Sweed 4.0)

#run sweed
system("cd ~/Desktop/Software/sweed-master/; ./SweeD -input /Users/ianbradbury/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_wholeChr_Europe.vcf -grid 200 -name ssa18_SmoltAge_sweed_EuropeAA_geno")
system("cd ~/Desktop/Software/sweed-master/; ./SweeD -input /Users/ianbradbury/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/BB_genotypes_only_Ssa18_wholeChr_Europe.vcf -grid 200 -name ssa18_SmoltAge_sweed_EuropeBB_geno")
system("cd ~/Desktop/Software/sweed-master/; ./SweeD -input /Users/ianbradbury/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/AA_genotypes_only_Ssa18_wholeChr_NorthAm.vcf -grid 200 -name ssa18_SmoltAge_sweed_NorthamAA_geno")
system("cd ~/Desktop/Software/sweed-master/; ./SweeD -input /Users/ianbradbury/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Genomic_data/CC_genotypes_only_Ssa18_wholeChr_NorthAm.vcf -grid 200 -name ssa18_SmoltAge_sweed_NorthamCC_geno")


##Combining and plotting data:

##Plotting SweeD seelection reuslts
library(ggplot2)

#set directory
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/sweed/")

sweed_aa_nor <- read.table("SweeD_Report.ssa18_SmoltAge_sweed_EuropeAA_geno", skip=2, header=T)
sweed_bb_nor <- read.table("SweeD_Report.ssa18_SmoltAge_sweed_EuropeBB_geno", skip=2, header=T)

#Get 95% quantile
aa_quantile <- quantile(sweed_aa_nor$Likelihood, 0.95)
sweed_aa_nor[which(sweed_aa_nor$Likelihood >= aa_quantile),]

bb_quantile <- quantile(sweed_bb_nor$Likelihood, 0.95)
sweed_bb_nor[which(sweed_bb_nor$Likelihood >= bb_quantile),]

max(sweed_bb_nor$Likelihood)
max(sweed_aa_nor$Likelihood)

 
Europe_plot <- ggplot()+
  geom_rect(aes(xmin=50282126, xmax=52929573, ymin=0, ymax=max(sweed_aa_nor$Likelihood)+1), alpha=0.9, fill="lightgray")+
  geom_area(data=sweed_bb_nor, aes(y=Likelihood, x=Position, fill="BB"), size=3, alpha=0.8)+
  geom_area(data=sweed_aa_nor, aes(y=Likelihood, x=Position,fill="AA"), size=3, alpha=0.6)+
  scale_fill_manual(values=c("firebrick2", "dodgerblue3"), name="Genotype")+
  theme_classic()+
  annotate(geom = "text", x = 0, y=20, label="Europe", size=5)+
 geom_hline(yintercept = aa_quantile[1], aes(col="AA"), col="firebrick2", lty=2)+
 geom_hline(yintercept = bb_quantile[1], aes(col="BB"), col="dodgerblue3", lty=2)+
  theme(panel.background = element_rect(colour = "white"),
        plot.background = element_rect(colour = "white"),
        panel.grid = element_blank(), axis.text = element_text(size=13),
        axis.title =element_text(size=18))+ylab("Composite likelihood ratio (CLR)")+xlab("Ssa18 position (bp)")
#####
# canada
sweed_aa_can <- read.table("SweeD_Report.ssa18_SmoltAge_sweed_NorthamAA_geno", skip=2, header=T)
sweed_cc_can <- read.table("SweeD_Report.ssa18_SmoltAge_sweed_NorthamCC_geno", skip=2, header=T)

#Get 95% quantile
aa_quantile_can <- quantile(sweed_aa_can$Likelihood, 0.95)
cc_quantile_can <- quantile(sweed_cc_can$Likelihood, 0.95)


max(sweed_aa_can$Likelihood)
max(sweed_cc_can$Likelihood)

NorthAM_plot<- ggplot()+
  geom_rect(aes(xmin=50282126, xmax=52929573, ymin=0, ymax=max(sweed_cc_can$Likelihood)+1), alpha=0.9, fill="lightgray")+
  geom_area(data=sweed_cc_can, aes(y=Likelihood, x=Position, fill="CC"), size=3, alpha=0.8)+
  geom_area(data=sweed_aa_can, aes(y=Likelihood, x=Position,fill="AA"), size=3, alpha=0.6)+
  scale_fill_manual(values=c("firebrick2", "dodgerblue3"), name="Genotype")+
  theme_classic()+
  annotate(geom = "text", x = 5e6, y=51, label="North America", size=5)+
  geom_hline(yintercept = aa_quantile_can[1], aes(col="AA"), col="firebrick2", lty=2)+
  geom_hline(yintercept = cc_quantile_can[1], aes(col="CC"), col="dodgerblue3", lty=2)+
  theme(panel.background = element_rect(colour = "white"),
        plot.background = element_rect(colour = "white"),
        panel.grid = element_blank(), axis.text = element_text(size=13),
        axis.title =element_text(size=18))+ylab("Composite likelihood ratio (CLR)")+xlab("Ssa18 position (bp)")


cowplot::plot_grid(NorthAM_plot, Europe_plot, nrow = 2, labels = c("A", "B"))
#save as 16 x 10 pdf

