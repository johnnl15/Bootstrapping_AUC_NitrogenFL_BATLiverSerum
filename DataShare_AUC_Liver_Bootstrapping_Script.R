#We thank Jinsu Park for help with bootstrapping and statistical methods. Method based on 
# Bootstrapping, Using Efron, B. and Tibshirani, R. (1993). An Introduction to the Bootstrap. 

library(stringi)
library(stringr)
library(purrr)
library(ggpubr)
library(reshape2)
library(cowplot)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr,warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(gtools)
library(RColorBrewer)
library(tidyr)

mycurrentdirectory = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mycurrentdirectory)

Tissue = "BAT"   #Designate tissue type here!!!####
Label = "13C5_15N_Glutamine" #Designate Labeling here

BAT<-read.csv("Liver_TotalFractionalNitrogenLabeling.csv")
Serum<-read.csv("Serum_TotalFractionalNitrogenLabeling_ForLiver.csv")

BATmelt<-melt(BAT,value.name = "batvalue")
Serummelt<-melt(Serum,value.name = "serumvalue")

b<-merge(BATmelt,Serummelt)

time<-b$variable     ####Change Names here #####
time=stri_replace_last_regex(time,"[:alpha:]",'')
time=stri_replace_last_regex(time,"_[:digit:]$",'')
time=stri_replace_last_regex(time,"2_5",'2.5')

temp<-b$variable ####Change Names here #####
temp=stri_replace_last_regex(temp,"[:digit:]+_[:digit:]+_[:digit:]+",'')
temp=stri_replace_last_regex(temp,"[:digit:]+_[:digit:]+",'')

b$time<-time
b$temp<-temp
b$time<-as.numeric(b$time)

compoundtypes<-unique(b$Compound)

b$Compound_temp<-paste(b$Compound,b$temp,sep = "_")

compound_temp_type<-unique(b$Compound_temp)

b2<-b

bootmet<-data.frame(compound_temp=(rep(NA,length(compound_temp_type))),
                    liver_AUC_average=(rep(NA,length(compound_temp_type))),
                    serum_AUC_average=(rep(NA,length(compound_temp_type))),
                    tstat=(rep(NA,length(compound_temp_type))),
                    CILower97.5=(rep(NA,length(compound_temp_type))),
                    CIHigher2.5=(rep(NA,length(compound_temp_type))),
                    pvalue_boot=(rep(NA,length(compound_temp_type))))

for (i in 1:length(compound_temp_type)) { 
  b<-b2%>% filter(Compound_temp==compound_temp_type[i])

alpha<-0.05 # significance level

time_set<-unique(b$time)[order(unique(b$time),decreasing = FALSE)]

n<-length(time_set)

bbat_mean<-b %>% group_by(time) %>% summarise(average=mean(batvalue)) 
colnames(bbat_mean)[colnames(bbat_mean)=="average"]<-"bataverage"
bserum_mean<-b %>% group_by(time) %>% summarise(average=mean(serumvalue)) 
colnames(bserum_mean)[colnames(bserum_mean)=="average"]<-"serumaverage"

b_mean_sep<-merge(bbat_mean,bserum_mean)

bbat_sd<-b %>% group_by(time) %>% summarise(average=sd(batvalue)) 
bserum_sd<-b %>% group_by(time) %>% summarise(average=sd(serumvalue)) 

################################################################
# Setting the Time
################################################################
a_set<-rep(NA,n)

a_set[1]<-(time_set[2]-time_set[1])/2
for(j in 2:n) {
  a_set[j]<-(time_set[j+1]-time_set[j-1])/2
}
a_set[n]<-(time_set[n]-time_set[n-1])/2

timetable<-data.frame(time=time_set,
                      a_set=a_set)

################################################################

b_mean_sep_merged<-merge(b_mean_sep,timetable)

AUC_mean<-b_mean_sep %>% summarise(BAT_AUC=sum(a_set*bataverage),Serum_AUC=sum(a_set*serumaverage))

AUC_sd<-b %>% 
  group_by(time) %>% 
  summarise(BATsd=sd(batvalue)^2/n(),
            Serumsd=sd(serumvalue)^2/n())
AUC_sd<-merge(AUC_sd,timetable)
AUC_sd$a_set<-AUC_sd$a_set^2
AUC_sd <- AUC_sd %>% group_by(time) %>% mutate(BATsdaset=BATsd*a_set,
                                                   Serumsdaset=Serumsd*a_set)
AUC_sd <- AUC_sd %>% ungroup() %>%
  summarise(BAT_AUC_sd=sqrt(sum(BATsdaset)),Serum_AUC_sd=sqrt(sum(Serumsdaset)))# %>%

BatSerum_mean<-b %>% group_by(time) %>% summarise(BatSerumAve=mean(c(batvalue,serumvalue)))

Boot_B<-merge(b,BatSerum_mean,by="time")
Boot_B<-merge(Boot_B,b_mean_sep)
Boot_BSmooth<-Boot_B %>% mutate(BATSmoothValue=batvalue-bataverage+BatSerumAve)
Boot_BSmooth<-Boot_BSmooth %>% mutate(SerumSmoothValue=serumvalue-serumaverage+BatSerumAve)
 
Boot_BSmooth<-merge(Boot_BSmooth,timetable,by="time")
Boot_BSmooth$BATSmoothValueAset<-Boot_BSmooth$BATSmoothValue*Boot_BSmooth$a_set #not using the smoothened value
Boot_BSmooth$SerumSmoothValueAset<-Boot_BSmooth$SerumSmoothValue*Boot_BSmooth$a_set #not using the smoothened value

library(foreach)
library(doParallel)
library(boot)
numCores <- detectCores()
registerDoParallel(numCores)  # use multicore, set to the number of our cores

trials <- 10

system.time({
Bootstrapping_stat <- foreach(icount(trials), .combine=c) %dopar% {
    samples<-Boot_BSmooth %>% 
      group_by(time) %>%
      do(sample_n(.,n(),replace = TRUE))
    samples2<-samples %>% group_by(time) %>% 
      summarise(BATaverage=mean(BATSmoothValueAset),Serumaverage=mean(SerumSmoothValueAset)) %>%
      summarise(BATAUC=sum(BATaverage),SerumAUC=sum(Serumaverage)) %>%
      summarise(Bootstrapping_difference=BATAUC-SerumAUC)
    samples3<-samples %>% 
      group_by(time) %>% summarise(BATsd=sd(BATSmoothValueAset)^2/n(),Serumsd=sd(SerumSmoothValueAset)^2/n())
    samples3<-merge(samples3,timetable)
    samples3$a_set<-samples3$a_set^2
    samples3 <- samples3 %>% group_by(time) %>% 
      mutate(BATsdaset=BATsd*a_set,Serumsdaset=Serumsd*a_set) %>% 
      ungroup() %>%
      summarise(BAT_sd=sqrt(sum(BATsdaset)),Serum_sd=sqrt(sum(Serumsdaset))) %>%
      summarise(Bootstrapping_sd=sqrt(BAT_sd^2+Serum_sd^2))
    samples2[[1]]/samples3[[1]] 
  }
})

t_stat<-(AUC_mean$BAT_AUC-AUC_mean$Serum_AUC)/sqrt(AUC_sd$BAT_AUC_sd^2+AUC_sd$Serum_AUC_sd^2)
CI_Boot<-AUC_mean$BAT_AUC-AUC_mean$Serum_AUC-sqrt(AUC_sd$BAT_AUC_sd^2+AUC_sd$Serum_AUC_sd^2)*quantile(Bootstrapping_stat,probs=c(1-alpha/2,alpha/2),na.rm = TRUE)
pvalue_Boot<-2*min(mean(t_stat>=Bootstrapping_stat),mean(t_stat<=Bootstrapping_stat))

bootmet[i,1]<-compound_temp_type[i]
bootmet[i,2]<-AUC_mean$BAT_AUC
bootmet[i,3]<-AUC_mean$Serum_AUC
bootmet[i,4]<-t_stat
bootmet[i,5]<-CI_Boot[1]
bootmet[i,6]<-CI_Boot[2]
bootmet[i,7]<-pvalue_Boot
}


boostrapped_result<-bootmet

boostrapped_result<-separate(boostrapped_result,compound_temp,into = c("compound","temp"),sep = "_")

boostrapped_result$compound<-tolower(boostrapped_result$compound)

boostrapped_result$compound[boostrapped_result$compound=="gln-15n-ta"]<-"15N-Gln"
boostrapped_result$compound[boostrapped_result$compound=="gln-15n2-13c5-tracer"]<-"15N,13C-Gln"
boostrapped_result$compound[boostrapped_result$compound=="asparagine"]<-"Asn"
boostrapped_result$compound[boostrapped_result$compound=="leucine"]<-"Leu"
boostrapped_result$compound[boostrapped_result$compound=="isoleucine"]<-"Ile"
boostrapped_result$compound[boostrapped_result$compound=="valine"]<-"Val"
boostrapped_result$compound[boostrapped_result$compound=="serine"]<-"Ser"
boostrapped_result$compound[boostrapped_result$compound=="glutamate"]<-"Glu"
boostrapped_result$compound[boostrapped_result$compound=="alanine"]<-"Ala"
boostrapped_result$compound[boostrapped_result$compound=="aspartate"]<-"Asp"
boostrapped_result$compound[boostrapped_result$compound=="phenylalanine"]<-"Phe"
boostrapped_result$compound[boostrapped_result$compound=="proline"]<-"Pro"
boostrapped_result$compound[boostrapped_result$compound=="tryptophan"]<-"Trp"
boostrapped_result$compound[boostrapped_result$compound=="tyrosine"]<-"Tyr"
boostrapped_result$compound[boostrapped_result$compound=="glycine"]<-"Gly"
boostrapped_result$compound[boostrapped_result$compound=="citrulline"]<-"Citrulline"
boostrapped_result$compound[boostrapped_result$compound=="arginine"]<-"Arg"
boostrapped_result$compound[boostrapped_result$compound=="methionine"]<-"Met"
boostrapped_result$compound[boostrapped_result$compound=="ornithine"]<-"Ornithine"
boostrapped_result$compound[boostrapped_result$compound=="serine"]<-"Ser"

boostrapped_result$CILower97.5_foldchange<-boostrapped_result$CILower97.5+boostrapped_result$serum_AUC_average
boostrapped_result$CIHigher2.5_foldchange<-boostrapped_result$CIHigher2.5+boostrapped_result$serum_AUC_average
boostrapped_result$CILower97.5_foldchange<-boostrapped_result$CILower97.5_foldchange/boostrapped_result$serum_AUC_average
boostrapped_result$CIHigher2.5_foldchange<-boostrapped_result$CIHigher2.5_foldchange/boostrapped_result$serum_AUC_average
boostrapped_result$Fold_average_CI<-(boostrapped_result$CIHigher2.5_foldchange+boostrapped_result$CILower97.5_foldchange)/2

boostrapped_result$Tissue<-"Liver"

write.csv(boostrapped_result,file = "Bootstrapped_Result_Liver.csv",row.names = FALSE)



