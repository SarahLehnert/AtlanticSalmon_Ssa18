#Relationship between smolt age, allele frequency, and environmental data
library(ggplot2)
library(data.table)

#set workding directory
setwd("~/Desktop/Sarah/Salmon/Smolt_age/UpdatedData_April2022/Environmental/")

#read in data and scale it

#scale environmental data for analyses - apply scaling only to enviro data (bio1-19)
#For canada
canada_data <- read.table("SmoltAge_Enviro_canada.txt", header = T)
scaleenv_canada<-as.data.frame(apply(canada_data[,16:34], 2, function(x) scale(x)))
scaleenv_canada$Pop <- canada_data$PopCode
scaleenv_canada$FreqA_allele <- canada_data$FreqA_allele

#for norway
nor_data <- read.table("SmoltAge_Enviro_norway.txt", header = T)
scaleenv_norway<-as.data.frame(apply(nor_data[,17:35], 2, function(x) scale(x)))
scaleenv_norway$Pop <- nor_data$V1
scaleenv_norway$FreqA_allele <- nor_data$FreqA_allele


##Run glm between Freq A and BioClim data 1-19 using for loop - using scaled environmental data
#Run glm and get R2 value (calcualted as residual deviance divided by null deviance and subtracted from 1)
#For canada
dat<-NULL
for(i in 1:19){
  mod<-glm(scaleenv_canada$FreqA_allele~scaleenv_canada[,i])
  r2_mod= round(1-(mod$deviance/mod$null.deviance),3)
  res<-cbind(colnames(scaleenv_canada[,i]), r2_mod)
  dat<-as.data.frame(rbind(dat, res))
  rm(res)
}

Canada_R2_relationship <- as.data.frame(cbind(dat, paste0("bio",1:19)))
colnames(Canada_R2_relationship) <- c("R2_can", "Var")

##For europe
dat_nor<-NULL
for(i in 1:19){
  mod<-glm(scaleenv_norway$FreqA_allele~scaleenv_norway[,i])
  r2_mod= round(1-(mod$deviance/mod$null.deviance),3)
  res1<-cbind(colnames(scaleenv_norway[,i]), r2_mod)
  dat_nor<-as.data.frame(rbind(dat_nor, res1))
  rm(res1)
}

Nor_R2_relationship <- as.data.frame(cbind(dat_nor, paste0("bio",1:19)))
colnames(Nor_R2_relationship) <- c("R2", "Var")


#Combine Norway and Canada results for table
Results_glm <- cbind("R2_nor"=Nor_R2_relationship[,1], Canada_R2_relationship)
Results_glm[order(Results_glm$R2_nor, decreasing = T),]
Results_glm[order(Results_glm$R2_can, decreasing = T),]

#save results
write.table(Results_glm, file = "GLM_results_Bioclim_by_Allelfreq_both_continents.txt", quote = F, row.names = F, col.names = T, sep="\t")

#Check GLM results for BIO10 for both continents:
summary(glm(FreqA_allele~bio10, data=scaleenv_canada))
summary(glm(FreqA_allele~bio10, data=scaleenv_norway))


#plot results of GLM
#Create long data frame for each continent dataset to plot using facet

#Norway data long
d1<-NULL
for(i in 1:19){
temp <- (scaleenv_norway[, c("Pop", "FreqA_allele", paste0("bio", i))])
colnames(temp)<-c("pop", "FreqA_allele", "bio")
temp<- cbind(temp, "var"=rep(paste0("bio", i) ))
d1<-as.data.frame(rbind(d1, temp))
rm(temp)  
}
EU_longdata <- d1

#Canada data long
t1<-NULL
for(i in 1:19){
  temp <- (scaleenv_canada[, c("Pop", "FreqA_allele", paste0("bio", i))])
  colnames(temp)<-c("pop", "FreqA_allele", "bio")
  temp<- cbind(temp, "var"= rep(paste0("bio", i) ))
  t1<-as.data.frame(rbind(t1, temp))
  rm(temp)  
}
Can_longdata <- t1
head(Can_longdata)


##Plot results using facet 

#Plot for canada
ggplot(data=Can_longdata, aes(y=FreqA_allele, x=bio))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(.~var, scales="free_x")+xlab("Scaled environmental data (bioclim)")+ylab("A allele frequency")+
  theme_bw()+theme(axis.text = element_text(size=12),
  axis.title = element_text(size=15), strip.text = element_text(size=13))

#Plot for Europe
ggplot(data=EU_longdata, aes(y=FreqA_allele, x=bio))+
  geom_point()+
  geom_smooth(method = "lm", se = T)+
  facet_wrap(.~var, scales="free_x")+xlab("Scaled environmental data (bioclim)")+ylab("A allele frequency")+
  theme_bw()+theme(axis.text = element_text(size=12),
                   axis.title = element_text(size=15), strip.text = element_text(size=13))



###Nice plot of glm for BIO10


#Set up binomial logistic regression with a pass fail based on 1000 instances (decimal points)
can_dat_bio10<-Can_longdata[which(Can_longdata$var=="bio10"),]

can_dat_bio10$success <- round(can_dat_bio10$FreqA_allele,3)*1000 #adults
can_dat_bio10$fails <- (1-round(can_dat_bio10$FreqA_allele,3))*1000
resp=cbind(can_dat_bio10$success,can_dat_bio10$fails)
Lat2=can_dat_bio10$bio #model variable  #Lat2 = winter bottom temperature
mod=glm(resp~Lat2,family=binomial(logit))

#Generate new data and develop model prediction for plotting
newdata=with(can_dat_bio10, data.frame(Lat2 = seq(min(Lat2), max(Lat2), length = 200))) # example data (200 lats)
newdata$slaggedscores <- predict(mod,newdata=newdata,type="response")
newdata$se <- predict(mod,newdata=newdata,type="response",se.fit=T)$se.fit

#for plotting latitude
can_dat_bio10_withLat<-merge(x=can_dat_bio10,y=canada_data[,c(1,9)], by.x=1, by.y=1, sort = F, all.x=T )


NorthAm_plot_bio10 <- 
  ggplot()+
  geom_line(data=newdata, aes(Lat2, slaggedscores), lwd=1.5)+
  geom_point(data=can_dat_bio10_withLat, aes(x=bio,y=FreqA_allele, col=Lat),size=3)+
  ylab("Frequency A allele")+xlab("Summer temperature (bio10 scaled)")+theme_bw()+
  scale_colour_gradientn(colors=c("red","purple", "blue", "midnightblue"), name="Lat")+theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=15))

##Europe

EU_dat_bio10<-EU_longdata[which(EU_longdata$var=="bio10"),]

EU_dat_bio10$success <- round(EU_dat_bio10$FreqA_allele,3)*1000 #adults
EU_dat_bio10$fails <- (1-round(EU_dat_bio10$FreqA_allele,3))*1000
resp_eu=cbind(EU_dat_bio10$success,EU_dat_bio10$fails)
Lat2_eu=EU_dat_bio10$bio #model variable  #Lat2_eu = winter bottom temperature
mod2=glm(resp_eu~Lat2_eu,family=binomial(logit))

#Generate new data and develop model prediction for plotting
newdata_eu=with(EU_dat_bio10, data.frame(Lat2_eu = seq(min(Lat2_eu), max(Lat2_eu), length = 200))) # example data (200 lats)
newdata_eu$slaggedscores <- predict(mod2,newdata=newdata_eu,type="response")
newdata_eu$se <- predict(mod2,newdata=newdata_eu,type="response",se.fit=T)$se.fit

EU_dat_bio10_withLat<-merge(x=EU_dat_bio10,y=nor_data[,c(1,10)], by.x=1, by.y=1, sort = F, all.x=T )


Eu_plot_bio10 <- 
  ggplot()+
  geom_line(data=newdata_eu, aes(Lat2_eu, slaggedscores), lwd=1.5)+
  geom_point(data=EU_dat_bio10_withLat, aes(x=bio,y=FreqA_allele, col=Lat),size=3)+
  ylab("Frequency A allele")+xlab("Summer temperature (bio10 scaled)")+theme_bw()+
  scale_colour_gradientn(colors=c("red","purple", "blue", "midnightblue"), name="Lat")+theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=15))

cowplot::plot_grid(NorthAm_plot_bio10,Eu_plot_bio10, nrow=2)
#saved as PDF 7.22 x 9.39
