remove(list=ls()) 
require(foreach)
require(doMC)
registerDoMC(64)
numits <- 1000

RES <- foreach (IT = 1:numits,.combine='rbind') %do% {
    require(psych)
 
  output.directory<-paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/UniCoTP/Sim",IT)
  
  #read in datafiles
  MZM <- read.table(file.path(output.directory, "MZM"), quote="\"", comment.char="", na.strings=-999)
  MZF <- read.table(file.path(output.directory, "MZF"), quote="\"", comment.char="", na.strings=-999)
  DZM <- read.table(file.path(output.directory, "DZM"), quote="\"", comment.char="", na.strings=-999)
  DZF <- read.table(file.path(output.directory, "DZF"), quote="\"", comment.char="", na.strings=-999)
  DZO <- read.table(file.path(output.directory, "DZOS"), quote="\"", comment.char="", na.strings=-999)
  PedigreeData <- read.table(file.path(output.directory, "PedigreeData"), quote="\"", comment.char="",na.strings=-999,header = T)
    
    #Select relative phenotype data to use 
    mzm <- MZM[,c(2:19)] 
    mzf <- MZF[,c(2:19)]
    dzm <- DZM[,c(2:19)]
    dzf <- DZF[,c(2:19)]
    dzo <- DZO[,c(2:19)]
    names(mzm)<-c("Tw1","Tw2","Fa","Mo","bro1","bro2","sis1","sis2","sp_tw1","sp_tw2","son1_tw1","son2_tw1","dau1_tw1","dau2_tw1","son1_tw2","son2_tw2","dau1_tw2","dau2_tw2")
    names(mzf)<-c("Tw1","Tw2","Fa","Mo","bro1","bro2","sis1","sis2","sp_tw1","sp_tw2","son1_tw1","son2_tw1","dau1_tw1","dau2_tw1","son1_tw2","son2_tw2","dau1_tw2","dau2_tw2")
    names(dzm)<-c("Tw1","Tw2","Fa","Mo","bro1","bro2","sis1","sis2","sp_tw1","sp_tw2","son1_tw1","son2_tw1","dau1_tw1","dau2_tw1","son1_tw2","son2_tw2","dau1_tw2","dau2_tw2")
    names(dzf)<-c("Tw1","Tw2","Fa","Mo","bro1","bro2","sis1","sis2","sp_tw1","sp_tw2","son1_tw1","son2_tw1","dau1_tw1","dau2_tw1","son1_tw2","son2_tw2","dau1_tw2","dau2_tw2")
    names(dzo)<-c("Tw1","Tw2","Fa","Mo","bro1","bro2","sis1","sis2","sp_tw1","sp_tw2","son1_tw1","son2_tw1","dau1_tw1","dau2_tw1","son1_tw2","son2_tw2","dau1_tw2","dau2_tw2")
    widedata<-rbind(mzm,mzf,dzm,dzf,dzo)
    
    #get variance components
    true_VA<-var(PedigreeData$A*sqrt(.2))
    true_VS<-var(PedigreeData$S*sqrt(.2))
    true_X<-var(PedigreeData$F)
    true_VD<-var(PedigreeData$D*sqrt(.2))
    true_VE<-var(PedigreeData$U*sqrt(.2))
    true_q<-var(PedigreeData$A)
    true_W<-cov(PedigreeData$A,PedigreeData$F)
    
    #calculate mu from spouse covariance
    tau_td<-cov(PedigreeData$mating.phenotype,PedigreeData$cur.phenotype)
    CV_sp<-1/3*(cov(widedata$Tw1,widedata$sp_tw1,use="pairwise.complete.obs")+cov(widedata$Tw2,widedata$sp_tw2,use="pairwise.complete.obs")+cov(widedata$Fa,widedata$Mo,use="pairwise.complete.obs"))
    true_mu<-CV_sp/(tau_td^2)
    
    rm(MZM, MZF, DZM, DZF, DZO, PedigreeData, mzm, mzf, dzm, dzf, dzo, widedata)
    gc()
    
    return(c(true_VA,true_VS,true_X,true_VD,true_VE,true_mu,true_q,true_W))
}
write.csv(RES,paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/Power/UniCoTP/GE_true_values_100genes_100kpop"))