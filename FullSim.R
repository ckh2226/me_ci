library(rineq)
library(tidyverse)

# Set constants 
# n <- 1000

# 0) Be reproducible queens
set.seed(205)

# Function to simulate data- choose variance of errors, sample size, and coefficients for model (to simulate CI magnitude)
## Generate true proximities (X) and rank them
## Generate response/health score (Y) using true ranks
## Generate additive errors (U); create error-prone Xstar by adding U to X
## Re-compute ranks using error-prone Xstar, save difference between true and error-prone ranks (w)
## Return variables: response (Y), true X, true R, indicator used to generate X (Z), error (U), 
## error-prone Xstar, error-prone Rstar, difference between true and error-prone ranks (w), mean of response (Y_bar)

sim_ef_data = function(sigmaU, n, alpha1, beta1) {
  # 1) Rural/not rural indicator 
  Z = rbinom(n = n, size = 1, prob = 0.2)
  
  # 2) True proximity | Rural / not rural 
  X = rnorm(n = n, mean = 1 + 4 * Z, sd = 1)
  
  # Calculate rank
  R <- (rank(X) - 1)/n + 1/(2 * n)
  
  # 3) Dietary inflammation score | Fractional rank of true proximity
  eps = rnorm(n = n, mean = 0, sd = 1)
  Y = alpha1 + beta1 * R + eps
  
  # Generate errors and error-prone X
  U = rnorm(n = n, mean = 0, sd = sigmaU)
  Xstar = X + U
  
  # Calculate error-prone rank
  Rstar <- (rank(Xstar) - 1)/n + 1/(2 * n)
  w <- Rstar - R
  
  # 4) Return 
  data.frame(Y, X, R, Z, U, Xstar, Rstar, w, Y_bar = mean(Y))
}

data <- sim_ef_data(0.5, 10000, -0.5, 3)
mu_hat <- data$Y_bar[1]
var_r <- var(data$R)
fit1 <- lm(Y ~ R, data = data)
beta1_hat <- fit1$coefficients[2]

our_ci <- 2 * var_r/mu_hat * beta1_hat

their_ci <- ci(data$X,
           data$Y,
           type = "CI",
           method = "linreg_delta", # Check method
           df_correction = TRUE)


# Function to compute concentration induces from above simulations
simulate_rep = function(sigmaU, n, alpha1, beta1) {
  
  temp <- sim_ef_data(sigmaU, n, alpha1, beta1) # First step of replication
  
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
  
  var_w <- var(temp$w)
  var_R <- var(temp$R)
  lambda <- var_R/(var_R + var_w)
  
  data.frame(ci_xstar = ci_xstar$concentration_index, ci_x = ci_x$concentration_index, sigmaU, var_w, var_R, lambda, n)
}

# Run multiple simulations at different levels of alpha1 and beta1
df_low_ci <- do.call(rbind, replicate(10000, simulate_rep(0.5, 1000, 2.5, -3), simplify = FALSE)) 
df_zero_ci <- do.call(rbind, replicate(10000, simulate_rep(0.5, 1000, 3, 0), simplify = FALSE)) 
df_high_ci <- do.call(rbind, replicate(10000, simulate_rep(0.5, 1000, -0.5, 3), simplify = FALSE)) 

df_low_small_n <- do.call(rbind, replicate(10000, simulate_rep(0.5, 1000, 2.5, -3), simplify = FALSE)) 
df_low_mid_n <- do.call(rbind, replicate(10000, simulate_rep(0.5, 5000, 2.5, -3), simplify = FALSE)) 
df_low_large_n <- do.call(rbind, replicate(10000, simulate_rep(0.5, 10000, 2.5, -3), simplify = FALSE)) 

df_low_small_var <- do.call(rbind, replicate(10000, simulate_rep(0.1, 5000, 2.5, -3), simplify = FALSE)) 
df_low_mid_var <- do.call(rbind, replicate(10000, simulate_rep(0.5, 5000, 2.5, -3), simplify = FALSE)) 
df_low_large_var <- do.call(rbind, replicate(10000, simulate_rep(1, 5000, 2.5, -3), simplify = FALSE)) 

ggplot(data = df_high_ci, aes(x = ci_xstar/ci_x, y = lambda)) + geom_point() + geom_abline(intercept = 0, slope = 1) +
  geom_smooth()


df_low_ci <- do.call(rbind, replicate(10000, simulate_rep(0.5, 5000, 2.5, -3), simplify = FALSE)) 
df_zero_ci <- do.call(rbind, replicate(10000, simulate_rep(0.5, 5000, 3, 0), simplify = FALSE)) 
df_high_ci <- do.call(rbind, replicate(10000, simulate_rep(0.5, 5000, -0.5, 3), simplify = FALSE)) 
                                                                            


df2 <- do.call(rbind, replicate(10000, simulate_rep(0.25, 1000), simplify = FALSE))
df3 <- do.call(rbind, replicate(10000, simulate_rep(3, 1000), simplify = FALSE))

g1 <- ggplot(data = df1, aes(x = ci_x)) + geom_histogram() + theme_bw() +
  labs(title = "sigmaU= 0.1")
g2 <- ggplot(data = df1, aes(x = ci_xstar)) + geom_histogram() + theme_bw() + 
  labs(title = "sigmaU= 0.1")

g3 <- ggplot(data = df2, aes(x = ci_x)) + geom_histogram() + theme_bw() + 
  labs(title = "sigmaU= 0.25")
g4 <- ggplot(data = df2, aes(x = ci_xstar)) + geom_histogram() + theme_bw() +
  labs(title = "sigmaU= 0.25")

g5 <- ggplot(data = df3, aes(x = ci_x)) + geom_histogram() + 
  theme_bw() + labs(title = "sigmaU= 3")
g6 <- ggplot(data = df3, aes(x = ci_xstar)) + geom_histogram() + theme_bw() + 
  labs(title = "sigmaU= 3")

gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 2)






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


