remove(list=ls()) 
require(foreach)
require(doMC)
library(lavaan)
registerDoMC(64)
numits <- 1000

# Build lavaan syntax for a 2-var CLPM with time-invariant parameters
clpm_lavaan_syntax <- function(T = 20,
                               A = matrix(c(0.3, 0.1,
                                            0.1, 0.3), 2, 2, byrow = TRUE),
                               intercept = c(0, 0),
                               Sigma_e = matrix(c(1, .3,
                                                  .3, 1), 2, 2, byrow = TRUE),
                               mu0 = c(0, 0),
                               Sigma0 = matrix(c(1, 0.3,
                                                 0.3, 1), 2, 2, byrow = TRUE)) {
  stopifnot(T >= 2)
  a_xx <- A[1,1]; a_xy <- A[1,2]
  a_yx <- A[2,1]; a_yy <- A[2,2]
  v_e_x <- Sigma_e[1,1]; v_e_y <- Sigma_e[2,2]; c_e_xy <- Sigma_e[1,2]
  v0_x  <- Sigma0[1,1]; v0_y  <- Sigma0[2,2];  c0_xy  <- Sigma0[1,2]
  int_x <- intercept[1]; int_y <- intercept[2]
  int_x0 <- mu0[1]; int_y0 <- mu0[2]
  
  # wave-1 variances/covariances and means
  lines <- c(
    sprintf("x1 ~~ %.10f*x1", v0_x),
    sprintf("y1 ~~ %.10f*y1", v0_y),
    sprintf("x1 ~~ %.10f*y1", c0_xy),
    sprintf("x1 ~ %.10f*1",  int_x0),
    sprintf("y1 ~ %.10f*1",  int_y0)
  )
  
  # waves 2..T
  for (t in 2:T) {
    xm1 <- sprintf("x%d", t-1); ym1 <- sprintf("y%d", t-1)
    xt  <- sprintf("x%d", t);   yt  <- sprintf("y%d", t)
    
    # regressions (fixed coefficients) + zero intercepts
    lines <- c(lines,
               sprintf("%s ~ %.10f*%s + %.10f*%s", xt, a_xx, xm1, a_xy, ym1),
               sprintf("%s ~ %.10f*%s + %.10f*%s", yt, a_yx, xm1, a_yy, ym1),
               sprintf("%s ~ %.10f*1", xt, int_x),
               sprintf("%s ~ %.10f*1", yt, int_y),
               
               # innovation variances and covariance (same each wave)
               sprintf("%s ~~ %.10f*%s", xt, v_e_x, xt),
               sprintf("%s ~~ %.10f*%s", yt, v_e_y, yt),
               sprintf("%s ~~ %.10f*%s", xt, c_e_xy, yt)
    )
  }
  
  paste(lines, collapse = "\n")
}
# simulate data at sample size of 200, 400, 800, 1600 individuals
foreach(j=c(1,2,4,8),.combine='rbind')%do%{
  foreach (IT = 1:numits,.combine='rbind') %dopar% {
  syntax <- clpm_lavaan_syntax()
  
  dat <- simulateData(model = syntax,
                          sample.nobs = 200*j,
                          empirical = FALSE)
  
  write.csv(dat, file = paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/CLPM_MVN/CLPM_MVN_j", j, "_IT", IT, ".csv"), row.names = FALSE)
  }
}