# Define matrices from the simulation
A <- matrix(c(0.3, 0.1, 
              0.1, 0.3), 2, 2, byrow = TRUE)

Sigma_e <- matrix(c(1, 0.3, 
                    0.3, 1), 2, 2, byrow = TRUE)

# Wave 1 starting covariance (Sigma0)
Sigma_t <- matrix(c(1, 0.3, 
                    0.3, 1), 2, 2, byrow = TRUE)

# Iterate forward to wave 18
for (t in 2:18) {
  Sigma_t <- A %*% Sigma_t %*% t(A) + Sigma_e
}

# The resulting matrix contains the true values
print(Sigma_t)

true_var_x18 <- Sigma_t[1,1]
true_var_y18 <- Sigma_t[2,2]
true_cov_x18_y18 <- Sigma_t[1,2]