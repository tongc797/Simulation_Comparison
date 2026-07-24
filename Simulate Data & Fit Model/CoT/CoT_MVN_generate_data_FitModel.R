remove(list=ls()) 
require(foreach)
require(doMC)
registerDoMC(64)
numits <- 1000

foreach (j=c(1,8)) %do% {
  RES <- foreach (IT = 1:numits,.combine='rbind', .errorhandling = "pass") %dopar% {
    require(OpenMx)		
    
    require(psych)
    
    require(MASS)

######################################################################
# Part 1: SIMULATE MZ and DZ Expected Covariance Matrices:
#---------------------------------------------------------------------
# Specify the variance components
# ------------------------------

a <- sqrt(.2)
VS <- .2
VE <- .2
VD <- .2
VTw <- 0
m <- sqrt(0.2/(2*1))
f <- 1

#specify nonlinear constraints
q1<-1.060875
x1<-0.3150366
w1<-0.2711373
mu<-0.1095111

######################################################################
# Generate Dataset based on above
######################################################################

# Specify sample size and proportion of MZ and DZ pairs
TotN <- 200*j	# Total Number of pairs to be simulated
pMZ <- 0.5		# Proportion of MZ pairs
pDZ <- 1-pMZ	# Proportion of DZ pairs
######################################################################

# Shortcuts
VP <- a*q1*a + VE + VD + VS + VTw + f*x1*f + 2*a*w1*f
delta <- q1*a + w1*f
theta <- (q1 - 0.5)*a + w1*f
xiMZ <- 0.5*a*q1 + 0.5*a*delta*mu*delta + f*m*delta + f*m*VP*mu*delta
lambdaDZ <- 0.5*a*(q1 - 0.5) + 0.5*a*delta*mu*theta + f*m*theta + f*m*VP*mu*theta

# MZ covariance
CvMZ <- a*q1*a + VD + VS + VTw + f*x1*f + 2*a*w1*f

# DZ covariance
CvDZ <- a*(q1 - 0.5)*a + 0.25*VD + VS + VTw + f*x1*f + 2*a*w1*f

# Spouse
CvSps <- VP*mu*VP

# Parent-child
CvPO <- 0.5*a*delta + 0.5*a*delta*mu*VP + f*m*VP + f*m*VP*mu*VP

# MZ-Nephew/Niece (MZ twin is the Uncle/Aunt)
CvAVMZ <- 0.5*a*delta + 0.5*a*delta*mu*CvMZ + f*m*CvMZ + f*m*CvMZ*mu*VP

# DZ-Nephew/Niece (DZ twin is the Uncle/Aunt)
CvAVDZ <- 0.5*a*theta + 0.5*a*delta*mu*CvDZ + f*m*CvDZ + f*m*CvDZ*mu*VP

# Cousin-Cousin
CvCSMZ <- 0.5*a*xiMZ + 0.5*a*delta*mu*CvAVMZ + f*m*CvAVMZ + f*m*CvAVMZ*mu*VP  # Cousin via MZ
CvCSDZ <- 0.5*a*lambdaDZ + 0.5*a*delta*mu*CvAVDZ + f*m*CvAVDZ + f*m*CvAVDZ*mu*VP  # Cousin via DZ

# Siblings-in-law
CvSpMZ <- VP*mu*CvMZ         # Siblings-in-law MZ
CvSpDZ <- VP*mu*CvDZ         # Siblings-in-law DZ
CvSpSpMZ <- VP*mu*CvMZ*mu*VP # Spouse of MZ1 and MZ2
CvSpSpDZ <- VP*mu*CvDZ*mu*VP # Spouse of DZ1 and DZ2

# Uncle/Aunt-in-law
CvSpAvMZ <- VP*mu*CvAVMZ     # Uncle/Aunt-in-law MZ
CvSpAvDZ <- VP*mu*CvAVDZ     # Uncle/Aunt-in-law DZ

MZ.rel.cv =as.matrix(rbind(
  cbind(VP,      CvMZ,    CvSps,     CvSpMZ,    CvPO,      CvAVMZ),
  cbind(CvMZ,    VP,      CvSpMZ,    CvSps,     CvAVMZ,    CvPO),
  cbind(CvSps,   CvSpMZ,  VP,        CvSpSpMZ,  CvPO,      CvSpAvMZ),
  cbind(CvSpMZ,  CvSps,   CvSpSpMZ,  VP,        CvSpAvMZ,  CvPO),
  cbind(CvPO,    CvAVMZ,  CvPO,      CvSpAvMZ,  VP,        CvCSMZ),
  cbind(CvAVMZ,  CvPO,    CvSpAvMZ,  CvPO,      CvCSMZ,    VP)))

DZ.rel.cv =as.matrix(rbind(
  cbind(VP,      CvDZ,    CvSps,     CvSpDZ,    CvPO,      CvAVDZ),
  cbind(CvDZ,    VP,      CvSpDZ,    CvSps,     CvAVDZ,    CvPO),
  cbind(CvSps,   CvSpDZ,  VP,        CvSpSpDZ,  CvPO,      CvSpAvDZ),
  cbind(CvSpDZ,  CvSps,   CvSpSpDZ,  VP,        CvSpAvDZ,  CvPO),
  cbind(CvPO,    CvAVDZ,  CvPO,      CvSpAvDZ,  VP,        CvCSDZ),
  cbind(CvAVDZ,  CvPO,    CvSpAvDZ,  CvPO,      CvCSDZ,    VP)))

mzraw = as.data.frame(mvrnorm(n=(TotN*pMZ),c(0,0,0,0,0,0), MZ.rel.cv, empirical=F))
names(mzraw) = c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2")
mzraw$zyg = 1

dzraw = as.data.frame(mvrnorm(n=(TotN*pDZ),c(0,0,0,0,0,0), DZ.rel.cv, empirical=F))
names(dzraw) = c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2")
dzraw$zyg = 2

Twin_data <-rbind(mzraw,dzraw)
write.csv(Twin_data, file = paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/UniCoTP_MVN_2/TwinData_j", j, "_IT", IT, ".csv"), row.names = FALSE)

# Select Data for Analysis
mzData <- subset(Twin_data, zyg==1)
dzData <- subset(Twin_data, zyg==2)

##########################################################################################
# Part 2: Set up and run CoT model
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
mxOption(NULL, "Default optimizer", "NPSOL")
#start values
(VA.st.I <- .3)
(VA.est.I <- TRUE) #TRUE" means we allow VA to be estimated freely
(VD.st.I <- 0)
(VD.est.I <- TRUE) #TRUE" means we allow VD to be estimated freely
(VS.st.I <- 0)
(VS.est.I <- TRUE) #TRUE" means we allow VS to be estimated freely
(VT.st.I <- 0)
(VT.est.I <- FALSE) #TRUE" means we allow VT(twin environment) to be estimated freely
(VF.st.I <- 0)
(VF.est.I <- TRUE) #TRUE" means we allow VF to be estimated freely
(VE.st.I <- 0)

#Set up the CoT Sib model
a <- mxMatrix(type="Lower",nrow=1,ncol=1,free=VA.est.I,values=VA.st.I,label="AddGenPath",name="a")
VD <- mxMatrix( type="Symm",nrow=1,ncol=1,free=VD.est.I,values=VD.st.I,label="DomVar", name="VD" ) 
VS <- mxMatrix( type="Symm",nrow=1,ncol=1,free=VS.est.I,values=VS.st.I,label="SibVar", name="VS" ) 
VTw <- mxMatrix( type="Symm",nrow=1,ncol=1,free=VT.est.I,values=VT.st.I,label="TwVar", name="VTw" ) 
VE <- mxMatrix( type="Symm",nrow=1,ncol=1,free=T,values=VE.st.I,label="ENvVar", name="VE" ) 
m <- mxMatrix(type="Lower",nrow=1,ncol=1,free=VF.est.I,values=VF.st.I,label="VTMPath", name="m")
mu <- mxMatrix(type="Lower",nrow=1,ncol=1,free=T,values=0,label="AMCopath",name="mu")
f <- mxMatrix(type="Lower",nrow=1,ncol=1,free=F,values=1,label="fPath",name="f")

#Defining & equating nonlinear constraint for q (the variance of latent variable A) and x (the variance of latent variable F) and w (the covariance between A and F).
q1 <- mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=1,label="LatentVarAddGen",name="q1")
q2 <- mxAlgebra(1 + delta*mu*delta,name="q2")
qCon <- mxConstraint(q1==q2,name='qCon')   

x1 <- mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=.01,label="LatentVarF",name="x1")
x2 <- mxAlgebra(m*VP*m + m*VP*m + m*VP*mu*VP*m + m*VP*mu*VP*m,name="x2")
xCon <- mxConstraint(x1==x2,name='xCon') 

w1 <- mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=.01,label="CovAF",name="w1")
w2 <- mxAlgebra(.5*m*delta + .5*m*delta + .5*m*VP*mu*delta + .5*m*VP*mu*delta, name="w2")
wCon <- mxConstraint(w1==w2,name='wCon') 

#mxAlgebra 
#Shortcuts
delta <- mxAlgebra(q1*a + w1*f,name="delta")
theta <- mxAlgebra((q1-.5)*a + w1*f,name="theta")
xiMZ <- mxAlgebra(.5*a*q1 + .5*a*delta*mu*delta + f*m*delta + f*m*VP*mu*delta,name="xiMZ") 
lambdaDZ <- mxAlgebra(.5*a*(q1-.5) + .5*a*delta*mu*theta + f*m*theta + f*m*VP*mu*theta,name="lambdaDZ")

VP <- mxAlgebra(a*q1*a + VE + VD + VS + VTw + f*x1*f + 2*a*w1*f,name="VP")

#MZ covariance
CvMZ <- mxAlgebra(a*q1*a + VD + VS + VTw + f*x1*f + 2*a*w1*f,name="CvMZ") 

#DZ covariance
CvDZ <- mxAlgebra(a*(q1-.5)*a + .25*VD + VS + VTw + f*x1*f + 2*a*w1*f,name="CvDZ") 

#Spouse
CvSps <- mxAlgebra(VP*mu*VP,name="CvSps") 

#Parent-child
CvPO <- mxAlgebra(.5*a*delta + .5*a*delta*mu*VP + f*m*VP + f*m*VP*mu*VP,name="CvPO") 

#MZ-Nephew/Niece(MZ twin is the Uncle/Aunt)
CvAVMZ <- mxAlgebra(.5*a*delta + .5*a*delta*mu*CvMZ + f*m*CvMZ + f*m*CvMZ*mu*VP,name="CvAVMZ") 

#DZ-Nephew/Niece(DZ twin is the Uncle/Aunt)
CvAVDZ <- mxAlgebra(.5*a*theta + .5*a*delta*mu*CvDZ + f*m*CvDZ + f*m*CvDZ*mu*VP,name="CvAVDZ") 

#Cousin-Cousin
CvCSMZ <- mxAlgebra(.5*a*xiMZ + .5*a*delta*mu*CvAVMZ + f*m*CvAVMZ + f*m*CvAVMZ*mu*VP,name="CvCSMZ") #Cousin via MZ
CvCSDZ <- mxAlgebra(.5*a*lambdaDZ + .5*a*delta*mu*CvAVDZ + f*m*CvAVDZ + f*m*CvAVDZ*mu*VP,name="CvCSDZ") #Cousin via DZ

#Siblings-in-law
CvSpMZ <- mxAlgebra(VP*mu*CvMZ,name="CvSpMZ") #Siblings-in-law MZ
CvSpDZ <- mxAlgebra(VP*mu*CvDZ,name="CvSpDZ") #Siblings-in-law DZ
CvSpSpMZ <- mxAlgebra(VP*mu*CvMZ*mu*VP,name="CvSpSpMZ") #Spouse of MZ1 and MZ2
CvSpSpDZ <- mxAlgebra(VP*mu*CvDZ*mu*VP,name="CvSpSpDZ") #Spouse of DZ1 and DZ2

#Uncle/Aunt-in-law
CvSpAvMZ <- mxAlgebra(VP*mu*CvAVMZ,name="CvSpAvMZ") #Uncle/Aunt-in-law MZ
CvSpAvDZ <- mxAlgebra(VP*mu*CvAVDZ,name="CvSpAvDZ") #Uncle/Aunt-in-law DZ

#Put these relative covariances together into MZ relatives and DZ relatives matrices
MZ.rel.cv <-    mxAlgebra(rbind(
  cbind(VP,      CvMZ,    CvSps,     CvSpMZ,    CvPO,      CvAVMZ),
  cbind(CvMZ,    VP,      CvSpMZ,    CvSps,     CvAVMZ,    CvPO),
  cbind(CvSps,   CvSpMZ,  VP,        CvSpSpMZ,  CvPO,      CvSpAvMZ),
  cbind(CvSpMZ,  CvSps,   CvSpSpMZ,  VP,        CvSpAvMZ,  CvPO),
  cbind(CvPO,    CvAVMZ,  CvPO,      CvSpAvMZ,  VP,        CvCSMZ),
  cbind(CvAVMZ,  CvPO,    CvSpAvMZ,  CvPO,      CvCSMZ,    VP)),
  dimnames=list(c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2"),c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2")),name="expCovMZRels")

DZ.rel.cv <-    mxAlgebra(rbind(
  cbind(VP,      CvDZ,    CvSps,     CvSpDZ,    CvPO,      CvAVDZ),
  cbind(CvDZ,    VP,      CvSpDZ,    CvSps,     CvAVDZ,    CvPO),
  cbind(CvSps,   CvSpDZ,  VP,        CvSpSpDZ,  CvPO,      CvSpAvDZ),
  cbind(CvSpDZ,  CvSps,   CvSpSpDZ,  VP,        CvSpAvDZ,  CvPO),
  cbind(CvPO,    CvAVDZ,  CvPO,      CvSpAvDZ,  VP,        CvCSDZ),
  cbind(CvAVDZ,  CvPO,    CvSpAvDZ,  CvPO,      CvCSDZ,    VP)),
  dimnames=list(c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2"),c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2")),name="expCovDZRels")

#Put the relatives means together into means matrix (same for MZ and DZ rels)
means <-   mxMatrix(type="Full", nrow=1, ncol=6, free=TRUE, values= .01, label=c("mTw1","mTw2","msp_tw1","msp_tw2","mson1_tw1","mson1_tw2"),dimnames=list(NULL,c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2")), name="expMean")

#Put in the data
dataMZRel <- mxData(observed=mzData, type="raw")
dataDZRel <- mxData(observed=dzData, type="raw")

#Objectives for two groups
objMZRels <- mxExpectationNormal(covariance="expCovMZRels",means="expMean", dimnames=c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2"))
objDZRels <- mxExpectationNormal(covariance="expCovDZRels",means="expMean", dimnames=c("Tw1","Tw2","sp_tw1","sp_tw2","son1_tw1","son1_tw2"))

#Combine groups
myfitfun <- mxFitFunctionML()
paramsI <- list(a,VD,VE,mu,VS,VTw,m,f,VP,q1,q2,qCon,x1,x2,xCon,w1,w2,wCon,delta,theta,xiMZ,lambdaDZ,CvMZ,CvDZ,CvSps,CvPO,CvAVMZ,CvAVDZ,CvCSMZ,CvCSDZ,CvSpMZ,CvSpDZ,CvSpSpMZ,CvSpSpDZ,CvSpAvMZ,CvSpAvDZ,means,MZ.rel.cv,DZ.rel.cv,myfitfun)
modelMZ.CoTP <- mxModel("MZCoTP",paramsI,dataMZRel,objMZRels)
modelDZ.CoTP <- mxModel("DZCoTP",paramsI,dataDZRel,objDZRels)
objI <- mxFitFunctionMultigroup(c("MZCoTP","DZCoTP"))
CoTP.Model <- mxModel("CoTP",modelMZ.CoTP,modelDZ.CoTP,objI)
# -----------------------------------------------------------------------



#Run the model and get estimates
# -----------------------------------------------------------------------
#Run the model
CoTP.Fit <- mxTryHardWideSearch(CoTP.Model,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=.5, scale = .1)
Aest<-mxEval(a,CoTP.Fit$MZCoTP)
#if Aest very close to zero, re-run to see if more optimal likelihood is obtained
if(abs(Aest)<.01){CoTP.Fit <- mxTryHardWideSearch(CoTP.Model,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=1, scale = .1)}
CoTP.Fit.Summary <- summary(CoTP.Fit)
Aest<-mxEval(a,CoTP.Fit$MZCoTP)
Sest<-mxEval(VS,CoTP.Fit$MZCoTP)
Eest<-mxEval(VE,CoTP.Fit$MZCoTP)
Mest<-mxEval(m,CoTP.Fit$MZCoTP)
Dest<-mxEval(VD,CoTP.Fit$MZCoTP)
Test<-mxEval(VTw,CoTP.Fit$MZCoTP)
MUest<-mxEval(mu,CoTP.Fit$MZCoTP)
Qest<-mxEval(q1,CoTP.Fit$MZCoTP)
Xest<-mxEval(x1,CoTP.Fit$MZCoTP)
West<-mxEval(w1,CoTP.Fit$MZCoTP)

#fit submodels to get statistial power
modelnoA <- mxModel(CoTP.Fit, name="CoTPnoA")
modelnoA <- omxSetParameters(modelnoA, labels="AddGenPath", free=FALSE, values=0)
CoTPnoA.Fit <- mxTryHardWideSearch(modelnoA,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=.5, scale = .1)
CoTPnoA.Fit.Summary <- summary(CoTPnoA.Fit)
Aest_noA<-mxEval(a,CoTPnoA.Fit$MZCoTP)
Sest_noA<-mxEval(VS,CoTPnoA.Fit$MZCoTP)
Eest_noA<-mxEval(VE,CoTPnoA.Fit$MZCoTP)
Mest_noA<-mxEval(m,CoTPnoA.Fit$MZCoTP)
Dest_noA<-mxEval(VD,CoTPnoA.Fit$MZCoTP)
Test_noA<-mxEval(VTw,CoTPnoA.Fit$MZCoTP)
MUest_noA<-mxEval(mu,CoTPnoA.Fit$MZCoTP)

modelnoS <- mxModel(CoTP.Fit, name="CoTPnoS")
modelnoS <- omxSetParameters(modelnoS, labels="SibVar", free=FALSE, values=0)
CoTPnoS.Fit <- mxTryHardWideSearch(modelnoS,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=.5, scale = .1)
CoTPnoS.Fit.Summary <- summary(CoTPnoS.Fit)
Aest_noS <- mxEval(a, CoTPnoS.Fit$MZCoTP)
Sest_noS <- mxEval(VS, CoTPnoS.Fit$MZCoTP)
Eest_noS <- mxEval(VE, CoTPnoS.Fit$MZCoTP)
Mest_noS <- mxEval(m, CoTPnoS.Fit$MZCoTP)
Dest_noS <- mxEval(VD, CoTPnoS.Fit$MZCoTP)
Test_noS <- mxEval(VTw, CoTPnoS.Fit$MZCoTP)
MUest_noS <- mxEval(mu, CoTPnoS.Fit$MZCoTP)

modelnoM <- mxModel(CoTP.Fit, name="CoTPnoM")
modelnoM <- omxSetParameters(modelnoM, labels="VTMPath", free=FALSE, values=0)
CoTPnoM.Fit <- mxTryHardWideSearch(modelnoM,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=.5, scale = .1)
CoTPnoM.Fit.Summary <- summary(CoTPnoM.Fit)
Aest_noM <- mxEval(a, CoTPnoM.Fit$MZCoTP)
Sest_noM <- mxEval(VS, CoTPnoM.Fit$MZCoTP)
Eest_noM <- mxEval(VE, CoTPnoM.Fit$MZCoTP)
Mest_noM <- mxEval(m, CoTPnoM.Fit$MZCoTP)
Dest_noM <- mxEval(VD, CoTPnoM.Fit$MZCoTP)
Test_noM <- mxEval(VTw, CoTPnoM.Fit$MZCoTP)
MUest_noM <- mxEval(mu, CoTPnoM.Fit$MZCoTP)

modelnoD <- mxModel(CoTP.Fit, name="CoTPnoD")
modelnoD <- omxSetParameters(modelnoD, labels="DomVar", free=FALSE, values=0)
CoTPnoD.Fit <- mxTryHardWideSearch(modelnoD,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=.5, scale = .1)
CoTPnoD.Fit.Summary <- summary(CoTPnoD.Fit)
Aest_noD <- mxEval(a, CoTPnoD.Fit$MZCoTP)
Sest_noD <- mxEval(VS, CoTPnoD.Fit$MZCoTP)
Eest_noD <- mxEval(VE, CoTPnoD.Fit$MZCoTP)
Mest_noD <- mxEval(m, CoTPnoD.Fit$MZCoTP)
Dest_noD <- mxEval(VD, CoTPnoD.Fit$MZCoTP)
Test_noD <- mxEval(VTw, CoTPnoD.Fit$MZCoTP)
MUest_noD <- mxEval(mu, CoTPnoD.Fit$MZCoTP)

modelnoE <- mxModel(CoTP.Fit, name="CoTPnoE")
modelnoE <- omxSetParameters(modelnoE, labels="ENvVar", free=FALSE, values=0)
CoTPnoE.Fit <- mxTryHardWideSearch(modelnoE,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=.5, scale = .1)
CoTPnoE.Fit.Summary <- summary(CoTPnoE.Fit)
Aest_noE <- mxEval(a, CoTPnoE.Fit$MZCoTP)
Sest_noE <- mxEval(VS, CoTPnoE.Fit$MZCoTP)
Eest_noE <- mxEval(VE, CoTPnoE.Fit$MZCoTP)
Mest_noE <- mxEval(m, CoTPnoE.Fit$MZCoTP)
Dest_noE <- mxEval(VD, CoTPnoE.Fit$MZCoTP)
Test_noE <- mxEval(VTw, CoTPnoE.Fit$MZCoTP)
MUest_noE <- mxEval(mu, CoTPnoE.Fit$MZCoTP)

modelnoMU <- mxModel(CoTP.Fit, name="CoTPnoMU")
modelnoMU <- omxSetParameters(modelnoMU, labels="AMCopath", free=FALSE, values=0)
CoTPnoMU.Fit <- mxTryHardWideSearch(modelnoMU,intervals=FALSE,extraTries=30,exhaustive=T, OKstatuscodes = c(0,1), jitterDistrib = "rnorm", loc=.5, scale = .1)
CoTPnoMU.Fit.Summary <- summary(CoTPnoMU.Fit)
Aest_noMU <- mxEval(a, CoTPnoMU.Fit$MZCoTP)
Sest_noMU <- mxEval(VS, CoTPnoMU.Fit$MZCoTP)
Eest_noMU <- mxEval(VE, CoTPnoMU.Fit$MZCoTP)
Mest_noMU <- mxEval(m, CoTPnoMU.Fit$MZCoTP)
Dest_noMU <- mxEval(VD, CoTPnoMU.Fit$MZCoTP)
Test_noMU <- mxEval(VTw, CoTPnoMU.Fit$MZCoTP)
MUest_noMU <- mxEval(mu, CoTPnoMU.Fit$MZCoTP)

return(c(
  CoTP.Fit.Summary$Minus2LogLikelihood,
  CoTP.Fit$output$status$code,
  Aest[1,1],
  Sest[1,1],
  Eest[1,1],
  Mest[1,1],
  Dest[1,1],
  Test[1,1],
  MUest[1,1],
  Qest[1,1],
  Xest[1,1],
  West[1,1],
  
  CoTPnoA.Fit.Summary$Minus2LogLikelihood,
  CoTPnoA.Fit$output$status$code,
  Aest_noA[1,1],
  Sest_noA[1,1],
  Eest_noA[1,1],
  Mest_noA[1,1],
  Dest_noA[1,1],
  Test_noA[1,1],
  MUest_noA[1,1],
  
  CoTPnoS.Fit.Summary$Minus2LogLikelihood,
  CoTPnoS.Fit$output$status$code,
  Aest_noS[1,1],
  Sest_noS[1,1],
  Eest_noS[1,1],
  Mest_noS[1,1],
  Dest_noS[1,1],
  Test_noS[1,1],
  MUest_noS[1,1],
  
  CoTPnoM.Fit.Summary$Minus2LogLikelihood,
  CoTPnoM.Fit$output$status$code,
  Aest_noM[1,1],
  Sest_noM[1,1],
  Eest_noM[1,1],
  Mest_noM[1,1],
  Dest_noM[1,1],
  Test_noM[1,1],
  MUest_noM[1,1],
  
  CoTPnoD.Fit.Summary$Minus2LogLikelihood,
  CoTPnoD.Fit$output$status$code,
  Aest_noD[1,1],
  Sest_noD[1,1],
  Eest_noD[1,1],
  Mest_noD[1,1],
  Dest_noD[1,1],
  Test_noD[1,1],
  MUest_noD[1,1],
  
  CoTPnoE.Fit.Summary$Minus2LogLikelihood,
  CoTPnoE.Fit$output$status$code,
  Aest_noE[1,1],
  Sest_noE[1,1],
  Eest_noE[1,1],
  Mest_noE[1,1],
  Dest_noE[1,1],
  Test_noE[1,1],
  MUest_noE[1,1],
  
  CoTPnoMU.Fit.Summary$Minus2LogLikelihood,
  CoTPnoMU.Fit$output$status$code,
  Aest_noMU[1,1],
  Sest_noMU[1,1],
  Eest_noMU[1,1],
  Mest_noMU[1,1],
  Dest_noMU[1,1],
  Test_noMU[1,1],
  MUest_noMU[1,1]
))

  }
  write.csv(RES,paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/Power/UniCoTP/MVN_Estimates_noabound_",j))
}