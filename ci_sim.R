sim_ci = function(alpha1_beta1, n = 1000) {
  alpha1 = alpha1_beta1[1]
  beta1= alpha1_beta1[2]

  # 1) Rural/not rural indicator 
  Z = rbinom(n = n, size = 1, prob = 0.2)
  
  # 2) True proximity | Rural / not rural 
  X = rnorm(n = n, mean = 1 + 4 * Z, sd = 1)
  
  # Calculate rank
  R <- (rank(X) - 1)/n + 1/(2 * n)
  
  # 3) Dietary inflammation score | True proximity, rural / not rural 
  eps = rnorm(n = n, mean = 0, sd = 1)
  Y = alpha1 + beta1 * R + eps
  
  # 4) Return 
  ci_x <- ci(X,
             Y,
             type = "CI",
             method = "linreg_delta", # Check method
             df_correction = TRUE)
  ci <- ci_x$concentration_index
  data.frame(alpha1, beta1, ci, Y_bar = mean(Y))
  # data.frame(Y, X, R, Z, U, Xstar, Rstar)
}

params <- expand.grid(alpha1 = seq(from = -3, to = 3, by = 0.5),
                      beta1 = seq(from = -3, to = 3, by = 0.5))

apply(
  X = params,
  MARGIN = 1,
  FUN = sim_ci,
  n = 1500
)

do.call(rbind,
        apply(
          X = params,
          MARGIN = 1,
          FUN = sim_ci,
          n = 10000
        )
        ) |>
  filter(Y_bar >= 0) |>
  ggplot(aes(x = ci)) + geom_histogram() + xlim(c(-2, 2))

df1 <- do.call(rbind,
               apply(
                 X = params,
                 MARGIN = 1,
                 FUN = sim_ci,
                 n = 10000
               )
               )

low_ci <- df1 |>
  filter(ci > -0.55 & ci < -0.45) # try alpha1 = 2.5 and beta1 = -3.0, y_bar > 0

zero_ci <- df1 |>
  filter(ci > -0.01 & ci < 0.01 & Y_bar > 0) # try alpha1 = 3.0 and beta1 = 0.0

high_ci <- df1 |>
  filter(ci > 0.45 & ci < 0.55 & Y_bar > 0) # try alpha1 = -0.5 and beta = 3.0




## Re-run simulations using different values of alpha1, beta1
sim_ef_data = function(sigmaU, n, alpha1, beta1) {
  # 1) Rural/not rural indicator 
  Z = rbinom(n = n, size = 1, prob = 0.2)
  
  # 2) True proximity | Rural / not rural 
  X = rnorm(n = n, mean = 1 + 4 * Z, sd = 1)
  
  # Calculate rank
  R <- (rank(X) - 1)/n + 1/(2 * n)
  
  # 3) Dietary inflammation score | True proximity, rural / not rural 
  eps = rnorm(n = n, mean = 0, sd = 1)
  Y = alpha1 + beta1 * R + eps
  
  # Generate errors and error-prone X
  U = rnorm(n = n, mean = 0, sd = sigmaU)
  Xstar = X + U
  
  # Calculate error-prone rank
  Rstar <- (rank(Xstar) - 1)/n + 1/(2 * n)
  
  # 4) Return 
  data.frame(Y, X, R, Z, U, Xstar, Rstar)
}



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
  
  data.frame(ci_xstar = ci_xstar$concentration_index, ci_x = ci_x$concentration_index, sigmaU, n)
}

# Testing functions, running multiple simulations
low_ci_sim <- do.call(rbind, replicate(10000, simulate_rep(0.25, 1000, 2.5, -3), simplify = FALSE))
mid_ci_sim <- do.call(rbind, replicate(10000, simulate_rep(0.25, 1000, 3, 0), simplify = FALSE))
high_ci_sim <- do.call(rbind, replicate(10000, simulate_rep(0.25, 1000, -0.5, 3), simplify = FALSE))

g1 <- ggplot(data = low_ci_sim, aes(x = ci_x)) + geom_histogram() + theme_bw() +
  labs(title = "CI ~ -0.5")
g2 <- ggplot(data = low_ci_sim, aes(x = ci_xstar)) + geom_histogram() + theme_bw() + 
  labs(title = "CI ~ -0.5")

g3 <- ggplot(data = mid_ci_sim, aes(x = ci_x)) + geom_histogram() + theme_bw() + 
  labs(title = "CI ~ 0")
g4 <- ggplot(data = mid_ci_sim, aes(x = ci_xstar)) + geom_histogram() + theme_bw() +
  labs(title = "CI ~ 0")

g5 <- ggplot(data = high_ci_sim , aes(x = ci_x)) + geom_histogram() + 
  theme_bw() + labs(title = "CI ~ 0.5")
g6 <- ggplot(data = high_ci_sim, aes(x = ci_xstar)) + geom_histogram() + theme_bw() + 
  labs(title = "CI ~ 0.5")

gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 2)


########
# Testing different values for sigmaU
## Choose a range of variances 
try_sigmas <- seq(from = 0, to = 10, by = 0.15)
n <- length(try_sigmas)

## Create a data frame to store simulation results (average concentration indices and variance used)
low_ci_df <- data.frame(sigmaU = try_sigmas, avg_ci_x = rep(NA, n), avg_ci_xstar = rep(NA, n))
mid_ci_df <- data.frame(sigmaU = try_sigmas, avg_ci_x = rep(NA, n), avg_ci_xstar = rep(NA, n))
high_ci_df <- data.frame(sigmaU = try_sigmas, avg_ci_x = rep(NA, n), avg_ci_xstar = rep(NA, n))

for(i in 1:n) {
  sim <- do.call(rbind, replicate(1000, simulate_rep(try_sigmas[i], 1000, 2.5, -3), simplify = FALSE))
  low_ci_df[i, 'avg_ci_x'] <- mean(sim$ci_x)
  low_ci_df[i, 'avg_ci_xstar'] <- mean(sim$ci_xstar)
}

for(i in 1:n) {
  sim <- do.call(rbind, replicate(1000, simulate_rep(try_sigmas[i], 1000, 3, 0), simplify = FALSE))
  mid_ci_df[i, 'avg_ci_x'] <- mean(sim$ci_x)
  mid_ci_df[i, 'avg_ci_xstar'] <- mean(sim$ci_xstar)
}

for(i in 1:n) {
  sim <- do.call(rbind, replicate(1000, simulate_rep(try_sigmas[i], 1000, -0.5, 3), simplify = FALSE))
  high_ci_df[i, 'avg_ci_x'] <- mean(sim$ci_x)
  high_ci_df[i, 'avg_ci_xstar'] <- mean(sim$ci_xstar)
}
########

# Plot some simulation results
## Pivot data frame
low_ci_df_long <- low_ci_df |>
  pivot_longer(cols = c('avg_ci_x', 'avg_ci_xstar'), names_to = 'variable', values_to = 'ci')

ggplot(low_ci_df_long, aes(x = sigmaU, y = ci, color = variable)) + geom_point() +
  labs(title = "Simulated Concentration Indices", x = "Variance of Error Term", y = "Concentration Index") +
  theme_bw()


mid_ci_df_long <- mid_ci_df |>
  pivot_longer(cols = c('avg_ci_x', 'avg_ci_xstar'), names_to = 'variable', values_to = 'ci')

ggplot(mid_ci_df_long, aes(x = sigmaU, y = ci, color = variable)) + geom_point() +
  labs(title = "Simulated Concentration Indices", x = "Variance of Error Term", y = "Concentration Index") +
  theme_bw()



high_ci_df_long <- high_ci_df |>
  pivot_longer(cols = c('avg_ci_x', 'avg_ci_xstar'), names_to = 'variable', values_to = 'ci')

ggplot(high_ci_df_long, aes(x = sigmaU, y = ci, color = variable)) + geom_point() +
  labs(title = "Simulated Concentration Indices", x = "Variance of Error Term", y = "Concentration Index") +
  theme_bw()
########
