remove(list=ls()) 
require(foreach)
require(doMC)
registerDoMC(64)
numits <- 1000

foreach (j=c(1,2,4,8)) %do% {
RES <- foreach (IT = 1:numits,.combine='rbind', .errorhandling = "pass") %dopar% {
    require(lavaan)		
    
    output.directory<-paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/CLPM_MVN/CLPM_MVN_j", j, "_IT", IT, ".csv")
    
dat <- read.csv(output.directory)
#select simulated data wave 18-20 to fit CLPM
clpm_lavaan_syntax <- function(T = 20) {
  stopifnot(T >= 2)
  
  lines <- c(
    "x18 ~~ vx0*x18",
    "y18 ~~ vy0*y18",
    "x18 ~~ cxy0*y18",
    "x18 ~ intx*1",
    "y18 ~ inty*1"
  )
  
  for (t in 19:T) {
    xm1 <- sprintf("x%d", t-1); ym1 <- sprintf("y%d", t-1)
    xt  <- sprintf("x%d", t);   yt  <- sprintf("y%d", t)
    
    # regressions (fixed coefficients) + zero intercepts
    lines <- c(lines,
               sprintf("%s ~ a_xx*%s + a_xy*%s", xt, xm1, ym1),
               sprintf("%s ~ a_yx*%s + a_yy*%s", yt, xm1, ym1),
               sprintf("%s ~ intx*1", xt),
               sprintf("%s ~ inty*1", yt),
               
               # innovation variances and covariance (same each wave)
               sprintf("%s ~~ vex*%s", xt, xt),
               sprintf("%s ~~ vey*%s", yt, yt),
               sprintf("%s ~~ cexy*%s", xt, yt)
    )
  }
  paste(lines, collapse = "\n")
}

CLPM<-clpm_lavaan_syntax()

# Run CLPM Model and save estimates
fitCLPM <- sem(CLPM, data=dat)
XXest<-unname(coef(fitCLPM)["a_xx"])
YYest<-unname(coef(fitCLPM)["a_yy"])
XYest<-unname(coef(fitCLPM)["a_xy"])
YXest<-unname(coef(fitCLPM)["a_yx"])
VXest<-unname(coef(fitCLPM)["vx0"])
VYest<-unname(coef(fitCLPM)["vy0"])
CXYest<-unname(coef(fitCLPM)["cxy0"])
VeXest<-unname(coef(fitCLPM)["vex"])
VeYest<-unname(coef(fitCLPM)["vey"])
CeXYest<-unname(coef(fitCLPM)["cexy"])

#fit submodels to obtain statistical power
CLPMnoXX <- paste(CLPM, "\na_xx == 0")
fitCLPMnoXX <- sem(CLPMnoXX, data = dat)
XXest_noXX<-unname(coef(fitCLPMnoXX)["a_xx"])
YYest_noXX<-unname(coef(fitCLPMnoXX)["a_yy"])
XYest_noXX<-unname(coef(fitCLPMnoXX)["a_xy"])
YXest_noXX<-unname(coef(fitCLPMnoXX)["a_yx"])
VXest_noXX<-unname(coef(fitCLPMnoXX)["vx0"])
VYest_noXX<-unname(coef(fitCLPMnoXX)["vy0"])
CXYest_noXX<-unname(coef(fitCLPMnoXX)["cxy0"])
VeXest_noXX<-unname(coef(fitCLPMnoXX)["vex"])
VeYest_noXX<-unname(coef(fitCLPMnoXX)["vey"])
CeXYest_noXX<-unname(coef(fitCLPMnoXX)["cexy"])

CLPMnoYY <- paste(CLPM, "\na_yy == 0")
fitCLPMnoYY <- sem(CLPMnoYY, data = dat)
XXest_noYY<-unname(coef(fitCLPMnoYY)["a_xx"])
YYest_noYY<-unname(coef(fitCLPMnoYY)["a_yy"])
XYest_noYY<-unname(coef(fitCLPMnoYY)["a_xy"])
YXest_noYY<-unname(coef(fitCLPMnoYY)["a_yx"])
VXest_noYY<-unname(coef(fitCLPMnoYY)["vx0"])
VYest_noYY<-unname(coef(fitCLPMnoYY)["vy0"])
CXYest_noYY<-unname(coef(fitCLPMnoYY)["cxy0"])
VeXest_noYY<-unname(coef(fitCLPMnoYY)["vex"])
VeYest_noYY<-unname(coef(fitCLPMnoYY)["vey"])
CeXYest_noYY<-unname(coef(fitCLPMnoYY)["cexy"])

CLPMnoXY <- paste(CLPM, "\na_xy == 0")
fitCLPMnoXY <- sem(CLPMnoXY, data = dat)
XXest_noXY<-unname(coef(fitCLPMnoXY)["a_xx"])
YYest_noXY<-unname(coef(fitCLPMnoXY)["a_yy"])
XYest_noXY<-unname(coef(fitCLPMnoXY)["a_xy"])
YXest_noXY<-unname(coef(fitCLPMnoXY)["a_yx"])
VXest_noXY<-unname(coef(fitCLPMnoXY)["vx0"])
VYest_noXY<-unname(coef(fitCLPMnoXY)["vy0"])
CXYest_noXY<-unname(coef(fitCLPMnoXY)["cxy0"])
VeXest_noXY<-unname(coef(fitCLPMnoXY)["vex"])
VeYest_noXY<-unname(coef(fitCLPMnoXY)["vey"])
CeXYest_noXY<-unname(coef(fitCLPMnoXY)["cexy"])

CLPMnoYX <- paste(CLPM, "\na_yx == 0")
fitCLPMnoYX <- sem(CLPMnoYX, data = dat)
XXest_noYX<-unname(coef(fitCLPMnoYX)["a_xx"])
YYest_noYX<-unname(coef(fitCLPMnoYX)["a_yy"])
XYest_noYX<-unname(coef(fitCLPMnoYX)["a_xy"])
YXest_noYX<-unname(coef(fitCLPMnoYX)["a_yx"])
VXest_noYX<-unname(coef(fitCLPMnoYX)["vx0"])
VYest_noYX<-unname(coef(fitCLPMnoYX)["vy0"])
CXYest_noYX<-unname(coef(fitCLPMnoYX)["cxy0"])
VeXest_noYX<-unname(coef(fitCLPMnoYX)["vex"])
VeYest_noYX<-unname(coef(fitCLPMnoYX)["vey"])
CeXYest_noYX<-unname(coef(fitCLPMnoYX)["cexy"])

CLPMnoCXY <- paste(CLPM, "\ncxy0 == 0")
fitCLPMnoCXY <- sem(CLPMnoCXY, data = dat)
XXest_noCXY<-unname(coef(fitCLPMnoCXY)["a_xx"])
YYest_noCXY<-unname(coef(fitCLPMnoCXY)["a_yy"])
XYest_noCXY<-unname(coef(fitCLPMnoCXY)["a_xy"])
YXest_noCXY<-unname(coef(fitCLPMnoCXY)["a_yx"])
VXest_noCXY<-unname(coef(fitCLPMnoCXY)["vx0"])
VYest_noCXY<-unname(coef(fitCLPMnoCXY)["vy0"])
CXYest_noCXY<-unname(coef(fitCLPMnoCXY)["cxy0"])
VeXest_noCXY<-unname(coef(fitCLPMnoCXY)["vex"])
VeYest_noCXY<-unname(coef(fitCLPMnoCXY)["vey"])
CeXYest_noCXY<-unname(coef(fitCLPMnoCXY)["cexy"])

CLPMnoCeXY <- paste(CLPM, "\ncexy == 0")
fitCLPMnoCeXY <- sem(CLPMnoCeXY, data = dat)
XXest_noCeXY<-unname(coef(fitCLPMnoCeXY)["a_xx"])
YYest_noCeXY<-unname(coef(fitCLPMnoCeXY)["a_yy"])
XYest_noCeXY<-unname(coef(fitCLPMnoCeXY)["a_xy"])
YXest_noCeXY<-unname(coef(fitCLPMnoCeXY)["a_yx"])
VXest_noCeXY<-unname(coef(fitCLPMnoCeXY)["vx0"])
VYest_noCeXY<-unname(coef(fitCLPMnoCeXY)["vy0"])
CXYest_noCeXY<-unname(coef(fitCLPMnoCeXY)["cxy0"])
VeXest_noCeXY<-unname(coef(fitCLPMnoCeXY)["vex"])
VeYest_noCeXY<-unname(coef(fitCLPMnoCeXY)["vey"])
CeXYest_noCeXY<-unname(coef(fitCLPMnoCeXY)["cexy"])

return(c(unname(fitMeasures(fitCLPM, "chisq"))[1],XXest[1],YYest[1],XYest[1],YXest[1],VXest[1],VYest[1],CXYest[1],VeXest[1],VeYest[1],CeXYest[1],unname(fitMeasures(fitCLPMnoXX, "chisq"))[1],XXest_noXX[1],YYest_noXX[1],XYest_noXX[1],YXest_noXX[1],VXest_noXX[1],VYest_noXX[1],CXYest_noXX[1],VeXest_noXX[1],VeYest_noXX[1],CeXYest_noXX[1],unname(fitMeasures(fitCLPMnoYY, "chisq"))[1],XXest_noYY[1],YYest_noYY[1],XYest_noYY[1],YXest_noYY[1],VXest_noYY[1],VYest_noYY[1],CXYest_noYY[1],VeXest_noYY[1],VeYest_noYY[1],CeXYest_noYY[1],unname(fitMeasures(fitCLPMnoXY, "chisq"))[1],XXest_noXY[1],YYest_noXY[1],XYest_noXY[1],YXest_noXY[1],VXest_noXY[1],VYest_noXY[1],CXYest_noXY[1],VeXest_noXY[1],VeYest_noXY[1],CeXYest_noXY[1],unname(fitMeasures(fitCLPMnoYX, "chisq"))[1],XXest_noYX[1],YYest_noYX[1],XYest_noYX[1],YXest_noYX[1],VXest_noYX[1],VYest_noYX[1],CXYest_noYX[1],VeXest_noYX[1],VeYest_noYX[1],CeXYest_noYX[1],unname(fitMeasures(fitCLPMnoCXY, "chisq"))[1],XXest_noCXY[1],YYest_noCXY[1],XYest_noCXY[1],YXest_noCXY[1],VXest_noCXY[1],VYest_noCXY[1],CXYest_noCXY[1],VeXest_noCXY[1],VeYest_noCXY[1],CeXYest_noCXY[1],unname(fitMeasures(fitCLPMnoCeXY, "chisq"))[1],XXest_noCeXY[1],YYest_noCeXY[1],XYest_noCeXY[1],YXest_noCeXY[1],VXest_noCeXY[1],VYest_noCeXY[1],CXYest_noCeXY[1],VeXest_noCeXY[1],VeYest_noCeXY[1],CeXYest_noCeXY[1]))
}
write.csv(RES,paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/Power/CLPM/MVN_Estimates_",j))
}