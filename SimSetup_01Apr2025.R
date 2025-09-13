library(rineq)
library(tidyverse)

# Set constants 
# n <- 1000

# 0) Be reproducible queens
set.seed(205)

beta0 = 6
beta1 = - 2
beta2 = 1
sim_ef_data = function(sigmaU, n) {
  # 1) Rural/not rural indicator 
  Z = rbinom(n = n, size = 1, prob = 0.2)
  
  # 2) True proximity | Rural / not rural 
  X = rnorm(n = n, mean = 1 + 4 * Z, sd = 1)
  
  # Calculate rank
  R <- (rank(X) - 1)/n + 1/(2 * n)
  
  # 3) Dietary inflammation score | True proximity, rural / not rural 
  eps = rnorm(n = n, mean = 0, sd = 1)
  Y = beta0 + beta1 * R + eps
  
  # Generate errors and error-prone X
  U = rnorm(n = n, mean = 0, sd = sigmaU)
  Xstar = X + U
  
  # Calculate error-prone rank
  Rstar <- (rank(Xstar) - 1)/n + 1/(2 * n)
  
  # 4) Return 
  data.frame(Y, X, R, Z, U, Xstar, Rstar)
}



simulate_rep = function(sigmaU, n) {
  
  temp <- sim_ef_data(sigmaU, n) # First step of replication
  
  # Calculating CIs steps 2a and 2b
  ci_xstar <- ci(temp$Xstar,
     temp$Y,
     type = "CI",
     method = "linreg_delta", # Check method
     df_correction = TRUE
  )
  
  ci_x <- ci(temp$X,
     temp$Y,
     type = "CI",
     method = "linreg_delta", # Check method
     df_correction = TRUE
  )
  
  data.frame(ci_xstar = ci_xstar$concentration_index, ci_x = ci_x$concentration_index, sigmaU, n)
}

########
# Testing different values for sigmaU
## Choose a range of variances 
try_sigmas <- seq(from = 0, to = 10, by = 0.15)
n <- length(try_sigmas)

## Create a data frame to store simulation results (average concentration indices and variance used)
ci_df <- data.frame(sigmaU = try_sigmas, avg_ci_x = rep(NA, n), avg_ci_xstar = rep(NA, n))

for(i in 1:n) {
  sim <- do.call(rbind, replicate(10000, simulate_rep(try_sigmas[i], 1000), simplify = FALSE))
  ci_df[i, 'avg_ci_x'] <- mean(sim$ci_x)
  ci_df[i, 'avg_ci_xstar'] <- mean(sim$ci_xstar)
}
########

# Plot some simulation results
## Pivot data frame
ci_df_long <- ci_df |>
  pivot_longer(cols = c('avg_ci_x', 'avg_ci_xstar'), names_to = 'variable', values_to = 'ci')

ggplot(ci_df_long, aes(x = sigmaU, y = ci, color = variable)) + geom_point() +
  labs(title = "Simulated Concentration Indices", x = 'Variance of Error', y = 'Concentration Index') +
  theme_bw()
########


do.call(rbind, replicate(10000, simulate_rep(0.1, 1000), simplify = FALSE))
do.call(rbind, replicate(10000, simulate_rep(0.25, 1000), simplify = FALSE))


# ALT: Y = rnorm(n = n, mean = beta0 + beta1 * X + beta2 * Z, sd = 1)

# 4a) Error-prone proximity | True proximity 
# sigmaU = 0.25 ## Low (0.1 - 0.25), Moderate (0.5), High (1)


# 4b) Error-prone proximity | True proximity, rural / not rural 
U = rnorm(n = n, mean = 0, sd = Z * (sigmaU + 0.25) + (1 - Z) * sigmaU)
Xstar = X + U

# 4c) Error-prone proximity | True proximity, dietary inflammation score 
U = rnorm(n = n, mean = 0, sd = (Y > 0) * (sigmaU + 0.25) + (Y <= 0) * sigmaU) 
Xstar = X + U

# 5) Simulations
# Number of replications
n <- 1000

# Data frame to store CIs
results <- data.frame(a = rep(NA, n), 
                      b = rep(NA, n), 
                      c = rep(NA, n), 
                      truth = rep(NA, n))

# For reference: The attenuation factor is defined 
lambda <- 1 / (1 + sigmaU ^ 2)

# 5a) Simulation using Error-prone proximity | True proximity 
sigmaU = 0.25 ## Low (0.1 - 0.25), Moderate (0.5), High (1)
for(i in 1:n) {
  temp = sim_ef_data()
  U = rnorm(n = n, mean = 0, sd = sigmaU)
  temp$Xstar = temp$X + U
  
  results[i, "a"] <- ci(temp$Xstar,
                        temp$Y,
                        type = "CI",
                        method = "linreg_delta", # Check method
                        df_correction = TRUE
  )
}

# 5b) Simulation using Error-prone proximity | True proximity, rural / not rural
for(i in 1:n) {
  U = rnorm(n = n, mean = 0, sd = Z * (sigmaU + 0.25) + (1 - Z) * sigmaU)
  Xstar = X + U
  
  results[i, "b"] <- ci(Xstar, 
                        Y,
                        type = "CI",
                        method = "linreg_delta", # Check method
                        df_correction = TRUE
  )
}

# 5c) Simulation using Error-prone proximity | True proximity, dietary inflammation score
for(i in 1:n) {
  U = rnorm(n = n, mean = 0, sd = (Y > 0) * (sigmaU + 0.25) + (Y <= 0) * sigmaU) 
  Xstar = X + U
  
  results[i, "c"] <- ci(Xstar,
                        Y,
                        type = "CI",
                        method = "linreg_delta", # Check method
                        df_correction = TRUE
  )
}

# 5d) Simulation using True proximity | Rural / not rural 
for(i in 1:n) {
  X = rnorm(n = n, mean = 1 + 4 * Z, sd = 1)
  
  results[i, "truth"] <- ci(X,
                        Y,
                        type = "CI",
                        method = "linreg_delta", # Check method
                        df_correction = TRUE
  )
}

# 6) Plot distribution of CIs
# 6a) Histograms
par(mfrow = c(2,2))
hist(results[, "a"], main = "Error-prone proximity | True proximity")
hist(results[, "b"], main = "Error-prone proximity | True proximity, rural / not rural")
hist(results[, "c"], main = "Error-prone proximity | True proximity, dietary inflammation score")
hist(results[, "truth"], main = "True proximity | Rural / not rural")

# 6b) Boxplots
par(mfrow = c(2,2), mar = c(1,1,1,1))
boxplot(results[, "a"], main = "Error-prone proximity | True proximity")
boxplot(results[, "b"], main = "Error-prone proximity | True proximity, rural / not rural")
boxplot(results[, "c"], main = "Error-prone proximity | True proximity, dietary inflammation score")
boxplot(results[, "truth"], main = "True proximity | Rural / not rural")


