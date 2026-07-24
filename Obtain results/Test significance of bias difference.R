#CLPM
remove(list=ls()) 
require(foreach)

true_value_FT <- c(.3,.3,.1,.1,1.1383929,1.1383929,.4092262,1,1,.3)
true_value_MVN <- c(.3,.3,.1,.1,1.1383929,1.1383929,.4092262,1,1,.3)
RES<-foreach (TYPE=c("MVN","FT"),.combine='rbind') %do% { 
  foreach(i=c(1,2,4,8),.combine='rbind') %do% {
    dat<-read.csv(paste0("~/Library/CloudStorage/OneDrive-UCB-O365/research/project/Simulation comparison/Power results/CLPM/",TYPE,"_Estimates_", i,".txt"))
    dat_S<-dat[,2:ncol(dat)]
    names(dat_S)<-c("CLPM_chisq","XXest", "YYest", "XYest", "YXest", "VXest","VYest", "CXYest", "VeXest", "VeYest", "CeXYest","CLPMnoXX_chisq","XXest_noXX", "YYest_noXX", "XYest_noXX", "YXest_noXX", "VXest_noXX","VYest_noXX", "CXYest_noXX", "VeXest_noXX", "VeYest_noXX", "CeXYest_noXX","CLPMnoYY_chisq","XXest_noYY", "YYest_noYY", "XYest_noYY", "YXest_noYY", "VXest_noYY","VYest_noYY", "CXYest_noYY", "VeXest_noYY", "VeYest_noYY", "CeXYest_noYY","CLPMnoXY_chisq","XXest_noXY", "YYest_noXY", "XYest_noXY", "YXest_noXY", "VXest_noXY","VYest_noXY", "CXYest_noXY", "VeXest_noXY", "VeYest_noXY", "CeXYest_noXY","CLPMnoYX_chisq","XX_noYXest", "YYest_noYX", "XYest_noYX", "YXest_noYX", "VXest_noYX","VYest_noYX", "CXYest_noYX", "VeXest_noYX", "VeYest_noYX", "CeXYest_noYX","CLPMnoCXY_chisq","XXest_noCXY", "YYest_noCXY", "XYest_noCXY", "YXest_noCXY", "VXest_noCXY","VYest_noCXY", "CXYest_noCXY", "VeXest_noCXY", "VeYest_noCXY", "CeXYest_noCXY","CLPMnoCeXY_chisq","XXest_noCeXY", "YYest_noCeXY", "XYest_noCeXY", "YXest_noCeXY", "VXest_noCeXY","VYest_noCeXY", "CXYest_noCeXY", "VeXest_noCeXY", "VeYest_noCeXY", "CeXYest_noCeXY")
    dat2 <- as.data.frame(lapply(dat_S, function(x) as.numeric(as.character(x))))
    
    dat2$TYPE<-TYPE
    dat2$Sim <- ifelse(dat2$TYPE == "MVN", 0, 1)
    
    dat2$biasXX<-dat2$XXest-get(paste0("true_value_",TYPE))[1]
    dat2$biasYY<-dat2$YYest-get(paste0("true_value_",TYPE))[2]
    dat2$biasXY<-dat2$XYest-get(paste0("true_value_",TYPE))[3]
    dat2$biasYX<-dat2$YXest-get(paste0("true_value_",TYPE))[4]
    dat2$biasVX<-dat2$VXest-get(paste0("true_value_",TYPE))[5]
    dat2$biasVY<-dat2$VYest-get(paste0("true_value_",TYPE))[6]
    dat2$biasCXY<-dat2$CXYest-get(paste0("true_value_",TYPE))[7]
    dat2$biasVeX<-dat2$VeXest-get(paste0("true_value_",TYPE))[8]
    dat2$biasVeY<-dat2$VeYest-get(paste0("true_value_",TYPE))[9]
    dat2$biasCeXY<-dat2$CeXYest-get(paste0("true_value_",TYPE))[10]
    
    return(dat2)
  }}
modelXX<-lm(biasXX~Sim, data = RES)
modelYY<-lm(biasYY~Sim, data = RES)
modelXY<-lm(biasXY~Sim, data = RES)
modelYX<-lm(biasYX~Sim, data = RES)
modelVX<-lm(biasVX~Sim, data = RES)
modelVY<-lm(biasVY~Sim, data = RES)
modelCXY<-lm(biasCXY~Sim, data = RES)
modelVeX<-lm(biasVeX~Sim, data = RES)
modelVeY<-lm(biasVeY~Sim, data = RES)
modelCeXY<-lm(biasCeXY~Sim, data = RES)
summary(modelXX)
summary(modelYY)
summary(modelYX)
summary(modelXY)
summary(modelVX)
summary(modelVY)
summary(modelCXY)
summary(modelVeX)
summary(modelVeY)
summary(modelCeXY)

#CoT
remove(list=ls()) 
require(foreach)

true_value_GE <- c(.212042, .2000046, .3151295, .2000117, .1999925, .1094478, sqrt(.2), sqrt(0.2/(2*1)))
true_value_MVN <- c(.212175, .2, .3150366, .2, .2, 0.1095111, sqrt(.2), sqrt(0.2/(2*1)))
RES<-foreach (TYPE=c("MVN","GE"),.combine='rbind') %do% { 
  foreach(i=c(1,2,4,8),.combine='rbind') %do% {
  # Conditionally set the file name string
  file_suffix <- ifelse(TYPE == "GE", "_Estimates_noabound_nomissing_", "_Estimates_noabound_")
  
  # Read the CSV using the dynamic suffix
  dat <- read.csv(paste0("~/Library/CloudStorage/OneDrive-UCB-O365/research/project/Simulation comparison/Power results/UniCoTP/", TYPE, file_suffix, i))
  dat_S<-dat[,2:ncol(dat)]
  names(dat_S)<-c(
    "CoTP_ll","CoTP_status","Aest","Sest","Eest","Mest","Dest","Test","MUest","Qest","Xest","West",
    "CoTPnoA_ll","CoTPnoA_status","Aest_noA","Sest_noA","Eest_noA","Mest_noA","Dest_noA","Test_noA","MUest_noA",
    "CoTPnoS_ll","CoTPnoS_status","Aest_noS","Sest_noS","Eest_noS","Mest_noS","Dest_noS","Test_noS","MUest_noS",
    "CoTPnoM_ll","CoTPnoM_status","Aest_noM","Sest_noM","Eest_noM","Mest_noM","Dest_noM","Test_noM","MUest_noM",
    "CoTPnoD_ll","CoTPnoD_status","Aest_noD","Sest_noD","Eest_noD","Mest_noD","Dest_noD","Test_noD","MUest_noD",
    "CoTPnoE_ll","CoTPnoE_status","Aest_noE","Sest_noE","Eest_noE","Mest_noE","Dest_noE","Test_noE","MUest_noE",
    "CoTPnoMU_ll","CoTPnoMU_status","Aest_noMU","Sest_noMU","Eest_noMU","Mest_noMU","Dest_noMU","Test_noMU","MUest_noMU"
  )
  
  dat2a<-subset(dat_S,CoTP_status<=1)
  dat2<-subset(dat2a,round(CoTP_ll,10)<=(round(CoTPnoA_ll,10)+.00001)|CoTPnoA_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoS_ll,10)+.00001)|CoTPnoS_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoM_ll,10)+.00001)|CoTPnoM_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoD_ll,10)+.00001)|CoTPnoD_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoMU_ll,10)+.00001)|CoTPnoMU_status>1)
  dat2 <- as.data.frame(lapply(dat2, function(x) as.numeric(as.character(x))))
  dat2$VAest<-dat2$Aest*dat2$Aest*dat2$Qest

  dat2$TYPE<-TYPE
  dat2$Sim <- ifelse(dat2$TYPE == "MVN", 0, 1)

  dat2$biasVA<-dat2$VAest-get(paste0("true_value_",TYPE))[1]
  dat2$biasVS<-dat2$Sest-get(paste0("true_value_",TYPE))[2]
  dat2$biasVF<-dat2$Xest-get(paste0("true_value_",TYPE))[3]
  dat2$biasVD<-dat2$Dest-get(paste0("true_value_",TYPE))[4]
  dat2$biasVE<-dat2$Eest-get(paste0("true_value_",TYPE))[5]
  dat2$biasMU<-dat2$MUest-get(paste0("true_value_",TYPE))[6]
  
  return(dat2)
  }}

modelVA<-lm(biasVA~Sim, data = RES)
modelVS<-lm(biasVS~Sim, data = RES)
modelVF<-lm(biasVF~Sim, data = RES)
modelVD<-lm(biasVD~Sim, data = RES)
modelVE<-lm(biasVE~Sim, data = RES)
modelMU<-lm(biasMU~Sim, data = RES)
summary(modelVA)
summary(modelVS)
summary(modelVF)
summary(modelVD)
summary(modelVE)
summary(modelMU)