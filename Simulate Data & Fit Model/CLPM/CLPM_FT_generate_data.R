remove(list=ls()) 
require(foreach)
require(doMC)
registerDoMC(64)
numits <- 1000

# Forward-time simulation for a 2-variable CLPM (time-invariant parameters)
# Model:
#   x_t = c_x + a_xx*x_{t-1} + a_xy*y_{t-1} + e_xt
#   y_t = c_y + a_yx*x_{t-1} + a_yy*y_{t-1} + e_yt
#   e_t ~ N(0, Σ_e)

simulate_clpm_invariant <- function(N = 200,
                                    T = 20,
                                    A = matrix(c(0.3, 0.1,
                                                 0.1, 0.3), 2, 2, byrow = TRUE),
                                    intercept = c(0, 0),
                                    Sigma_e = matrix(c(1, .3,
                                                       .3, 1), 2, 2, byrow = TRUE),
                                    mu0 = c(0, 0),
                                    Sigma0 = matrix(c(1, 0.3,
                                                      0.3, 1), 2, 2, byrow = TRUE)) {
  stopifnot(T >= 2)
  
  # Storage matrices for X and Y (individuals × waves)
  X <- matrix(NA_real_, N, T)
  Y <- matrix(NA_real_, N, T)
  
  # Sample initial states for wave 1
  XY1 <- MASS::mvrnorm(N, mu = mu0, Sigma = Sigma0)
  X[, 1] <- XY1[, 1]
  Y[, 1] <- XY1[, 2]
  
  # Forward simulation across waves 2..T
  for (t in 2:T) {
    prev <- cbind(X[, t - 1], Y[, t - 1])           # previous time step
    mean_t <- sweep(prev %*% t(A), 2, intercept, "+")       # deterministic part
    eps <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_e)  # residuals
    
    X[, t] <- mean_t[, 1] + eps[, 1]
    Y[, t] <- mean_t[, 2] + eps[, 2]
  }
  
  # Return data in tidy long format
  # create wide-format data: one row per person
  data <- data.frame(
    id = seq_len(N),
    X,
    Y
  )
  
  # rename columns: x1...xT, y1...yT
  colnames(data) <- c(
    "id",
    paste0("x", 1:T),
    paste0("y", 1:T)
  )
  
  list(data = data, X = X, Y = Y)
}
# simulate data at sample size of 200, 400, 800, 1600 individuals
foreach(j=c(1,2,4,8),.combine='rbind')%do%{
  foreach (IT = 1:numits,.combine='rbind') %dopar% {
  sim <- simulate_clpm_invariant(
    N = 200*j
  )
  write.csv(sim$data, file = paste0("/gpfs/alpine1/scratch/toch9438/GE75/GE.Output/CLPM_FT/CLPM_FT_j", j, "_IT", IT, ".csv"), row.names = FALSE)
  }
}