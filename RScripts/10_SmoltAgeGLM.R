
#Combining smolt age data and enviro data for plotting
library(ggplot2)

#set diretory
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/")

#Read in smolt age data
can_smolt <- read.table("Environmental/SmoltAge_Enviro_canada.txt", header=T)
nor_smolt <- read.table("Environmental/SmoltAge_Enviro_norway.txt", header=T)

head(can_smolt)


#Canada plots
summary(glm(can_smolt$FreqA_allele~can_smolt$SmoltAge))
1-(0.6397/3.1684)

#Set up binomial logistic regression with a pass fail based on 1000 instances (decimal points)
can_smolt$success <- round(can_smolt$FreqA_allele,3)*1000 #adults
can_smolt$fails <- (1-round(can_smolt$FreqA_allele,3))*1000
resp=cbind(can_smolt$success,can_smolt$fails)
Lat2=can_smolt$SmoltAge #model variable  #Lat2 = winter bottom temperature
mod=glm(resp~Lat2,family=binomial(logit))

#Generate new data and develop model prediction for plotting
newdata=with(can_smolt, data.frame(Lat2 = seq(min(Lat2), max(Lat2), length = 200))) # example data (200 lats)
newdata$slaggedscores <- predict(mod,newdata=newdata,type="response")
newdata$se <- predict(mod,newdata=newdata,type="response",se.fit=T)$se.fit

NorthAm_plot <- 
  ggplot()+
  geom_line(data=newdata, aes(Lat2, slaggedscores), lwd=1.5)+
  geom_point(data=can_smolt, aes(x=SmoltAge,y=FreqA_allele, col=Lat),size=3)+
  ylab("Frequency A allele")+xlab("Mean smolt age")+theme_bw()+
  scale_colour_gradientn(colors=c("red","purple", "blue", "midnightblue"), name="Lat")+theme(panel.grid = element_blank())+
 # annotate("text", x=4.2, y=0.9, label="GLM r2=0.80")+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=15))



#Europe relationship
summary(glm(nor_smolt$FreqA_allele~nor_smolt$mean))
1-(0.86272/1.24479)

nor_smolt$success <- round(nor_smolt$FreqA_allele,3)*1000 #adults
nor_smolt$fails <- (1-round(nor_smolt$FreqA_allele,3))*1000
resp_nor <- cbind(nor_smolt$success,nor_smolt$fails)
Lat2_nor <- nor_smolt$mean #model variable  #Lat2 = winter bottom temperature
mod_2_nor <- glm(resp_nor~Lat2_nor,family=binomial(logit))

#Generate new data and develop model prediction for plotting
newdata2_nor <- with(nor_smolt, data.frame(Lat2_nor = seq(min(Lat2_nor), max(Lat2_nor), length = 200))) # example data (200 lats)
newdata2_nor$slaggedscores <- predict(mod_2_nor,newdata=newdata2_nor,type="response")
newdata2_nor$se <- predict(mod_2_nor,newdata=newdata2_nor,type="response",se.fit=T)$se.fit



EU_plot <-  
  ggplot()+
  geom_line(data=newdata2_nor, aes(Lat2_nor, slaggedscores), lwd=1.5)+
  geom_point(data=nor_smolt, aes(x=mean,y=FreqA_allele, col=Lat),size=3)+
  ylab("Frequency A allele")+xlab("Mean smolt age")+theme_bw()+
  scale_colour_gradientn(colors=c("red","purple", "blue", "midnightblue"),
                         name = "Lat")+ theme(panel.grid = element_blank())+
 # annotate("text", x=4.7, y=0.9, label="GLM r2=0.31")+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=15))



cowplot::plot_grid(NorthAm_plot,EU_plot, nrow=2)
#saved as PDF 7.22 x 9.39

