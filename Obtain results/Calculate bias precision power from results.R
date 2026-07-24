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
    
    #cor(dat2[,2:11])
    param_cols <- c("XXest", "YYest", "XYest", "YXest", "VXest", "VYest", "CXYest", "VeXest", "VeYest", "CeXYest")
    means <- sapply(dat2[, param_cols], mean, na.rm = TRUE)
    sds   <- sapply(dat2[, param_cols], sd, na.rm = TRUE) #empirical standard error
    bias  <- means - get(paste0("true_value_",TYPE)) #absolute bias
    rel_bias <- bias/get(paste0("true_value_",TYPE)) #relative bias

    # Calculate p-values for testing if bias is significantly different from 0
    N_valid <- nrow(dat2)
    se <- sds / sqrt(N_valid)
    t_stat <- bias / se
    p_val_bias <- 2 * pt(-abs(t_stat), df = N_valid - 1)
    
    #calculate power
    datXX1<-dat2
    datXX1$dif<-datXX1$CLPMnoXX_chisq-datXX1$CLPM_chisq
    datXX2<-subset(datXX1,datXX1$dif>3.841459)
    powerXX<-nrow(datXX2)/nrow(datXX1)
    
    datYY1<-dat2
    datYY1$dif<-datYY1$CLPMnoYY_chisq-datYY1$CLPM_chisq
    datYY2<-subset(datYY1,datYY1$dif>3.841459)
    powerYY<-nrow(datYY2)/nrow(datYY1)
    
    datXY1<-dat2
    datXY1$dif<-datXY1$CLPMnoXY_chisq-datXY1$CLPM_chisq
    datXY2<-subset(datXY1,datXY1$dif>3.841459)
    powerXY<-nrow(datXY2)/nrow(datXY1)
    
    datYX1<-dat2
    datYX1$dif<-datYX1$CLPMnoYX_chisq-datYX1$CLPM_chisq
    datYX2<-subset(datYX1,datYX1$dif>3.841459)
    powerYX<-nrow(datYX2)/nrow(datYX1)
    
    datCXY1<-dat2
    datCXY1$dif<-datCXY1$CLPMnoCXY_chisq-datCXY1$CLPM_chisq
    datCXY2<-subset(datCXY1,datCXY1$dif>3.841459)
    powerCXY<-nrow(datCXY2)/nrow(datCXY1)
    
    datCeXY1<-dat2
    datCeXY1$dif<-datCeXY1$CLPMnoCeXY_chisq-datCeXY1$CLPM_chisq
    datCeXY2<-subset(datCeXY1,datCeXY1$dif>3.841459)
    powerCeXY<-nrow(datCeXY2)/nrow(datCeXY1)
    
    #plot estimates
    par(mfrow=c(3,4))  
    for (col in names(dat2)[1:11]) {
      hist(dat2[[col]], main=col, xlab="", col="lightblue",breaks=100)
    }
    
    return(c(nrow(dat2),powerXX,powerYY,powerXY,powerYX,powerCXY,powerCeXY,bias,rel_bias,sds,p_val_bias))
  }}
colnames(RES) <- c(
  "checkN","powerXX","powerYY","powerXY","powerYX","powerCXY","powerCeXY",
  "bias_XX","bias_YY","bias_XY","bias_YX","bias_VX","bias_VY","bias_CXY","bias_VeX","bias_VeY","bias_CeXY",
  "relbias_XX","relbias_YY","relbias_XY","relbias_YX","relbias_VX","relbias_VY","relbias_CXY","relbias_VeX","relbias_VeY","relbias_CeXY",
  "sds_XX","sds_YY","sds_XY","sds_YX","sds_VX","sds_VY","sds_CXY","sds_VeX","sds_VeY","sds_CeXY",
  "pval_XX","pval_YY","pval_XY","pval_YX","pval_VX","pval_VY","pval_CXY","pval_VeX","pval_VeY","pval_CeXY"
)
n <- 4
rearranged_RES <- do.call(data.frame, unlist(lapply(1:ncol(RES), function(i) {
  list(RES[1:n, i], RES[(n+1):(2*n), i])
}), recursive = FALSE))
old_names <- colnames(RES)
new_colnames <- as.vector(t(outer(old_names, c("MVN", "FT"), paste, sep = "_")))
colnames(rearranged_RES) <- new_colnames
print(rearranged_RES)

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
  #only keep status code green results
  #And, discard results if the more constrained model fits better than the less constrained model (it should not happen, and if it happens, it means the less constrained model did not find the optimal solution and the results are untrustworthy)
  dat2a<-subset(dat_S,CoTP_status<=1)
  dat2<-subset(dat2a,round(CoTP_ll,10)<=(round(CoTPnoA_ll,10)+.00001)|CoTPnoA_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoS_ll,10)+.00001)|CoTPnoS_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoM_ll,10)+.00001)|CoTPnoM_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoD_ll,10)+.00001)|CoTPnoD_status>1)
  dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoMU_ll,10)+.00001)|CoTPnoMU_status>1)
  
  dat2 <- as.data.frame(lapply(dat2, function(x) as.numeric(as.character(x))))
  dat2$VAest<-dat2$Aest*dat2$Aest*dat2$Qest #calculate VA estimate
  dat2$absAest<-abs(dat2$Aest)
  #cor(dat2[,c(3:7,9:12)])
  param_cols <- c("VAest","Sest","Xest","Dest","Eest","MUest","absAest","Mest")
  means <- sapply(dat2[, param_cols], mean, na.rm = TRUE)
  sds   <- sapply(dat2[, param_cols], sd, na.rm = TRUE) #empirical standard error
  bias  <- means - get(paste0("true_value_",TYPE)) #absolute bias
  rel_bias <- bias/get(paste0("true_value_",TYPE)) #relative bias
  
  #calculate p values for absolute bias
  N_valid <- nrow(dat2)
  se <- sds / sqrt(N_valid)
  t_stat <- bias / se
  p_val_bias <- 2 * pt(-abs(t_stat), df = N_valid - 1)
  
  #calculate statistical power
  datA1<-subset(dat2,CoTPnoA_status<=1)
  datA1$dif<-datA1$CoTPnoA_ll-datA1$CoTP_ll
  datA2<-subset(datA1,datA1$dif>3.841459)
  powerA<-nrow(datA2)/nrow(datA1)
  
  datS1 <- subset(dat2, CoTPnoS_status <= 1)
  datS1$dif <- datS1$CoTPnoS_ll -datS1$CoTP_ll
  datS2 <- subset(datS1, datS1$dif > 3.841459)
  powerS <- nrow(datS2) / nrow(datS1)
  
  datM1 <- subset(dat2, CoTPnoM_status <= 1)
  datM1$dif <- datM1$CoTPnoM_ll -datM1$CoTP_ll
  datM2 <- subset(datM1, datM1$dif > 3.841459)
  powerM <- nrow(datM2) / nrow(datM1)
  
  datD1 <- subset(dat2, CoTPnoD_status <= 1)
  datD1$dif <- datD1$CoTPnoD_ll -datD1$CoTP_ll
  datD2 <- subset(datD1, datD1$dif > 3.841459)
  powerD <- nrow(datD2) / nrow(datD1)
  
  datMU1 <- subset(dat2, CoTPnoMU_status <= 1)
  datMU1$dif <- datMU1$CoTPnoMU_ll -datMU1$CoTP_ll
  datMU2 <- subset(datMU1, datMU1$dif > 3.841459)
  powerMU <- nrow(datMU2) / nrow(datMU1)
  
  #plot estimates
  par(mfrow=c(3,4))  
  for (col in names(dat2)[1:12]) {
    hist(dat2[[col]], main=col, xlab="", col="lightblue",breaks=100)
  }
  
  return(c(min(dat_S$CoTPnoE_status),powerA,powerS,powerM,powerD,powerMU,nrow(datA1),nrow(datS1),nrow(datM1),nrow(datD1),nrow(datMU1),1000-nrow(dat2a),nrow(dat2a)-nrow(dat2),bias,rel_bias,sds,p_val_bias))
}}
colnames(RES) <- c(
  "checkE","powerA","powerS","powerM","powerD","powerMU",
  "n_datA1","n_datS1","n_datM1","n_datD1","n_datMU1","n_red","n_nonoptimal",
  "bias_VA","bias_VS","bias_VF","bias_VD","bias_VE","bias_MU","bias_absAest","bias_M",
  "relbias_VA","relbias_VS","relbias_VF","relbias_VD","relbias_VE","relbias_MU","relbias_absAest","relbias_M",
  "sds_VA","sds_VS","sds_VF","sds_VD","sds_VE","sds_MU","sds_absAest","sds_M",
  "pval_VA","pval_VS","pval_VF","pval_VD","pval_VE","pval_MU","pval_absAest","pval_M"
)
n <- 4
rearranged_RES <- do.call(data.frame, unlist(lapply(1:ncol(RES), function(i) {
  list(RES[1:n, i], RES[(n+1):(2*n), i])
}), recursive = FALSE))
old_names <- colnames(RES)
new_colnames <- as.vector(t(outer(old_names, c("MVN", "GE"), paste, sep = "_")))
colnames(rearranged_RES) <- new_colnames
print(rearranged_RES)