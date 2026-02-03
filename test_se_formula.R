# Load packages
library(rineq) ## for CI calculations
library(tidyverse) ## for data wrangling, etc

# Be reproducible queens
set.seed(205)

# Function to simulate data- choose variance of errors, sample size, and coefficients for model (to simulate CI magnitude)
sim_data = function(sigmaU, n, alpha1, beta1) {
  # 1) Rural/not rural indicator
  Z = rbinom(n = n, size = 1, prob = 0.2)

  # 2) True proximity | Rural / not rural
  X = rnorm(n = n, mean = 1 + 4 * Z, sd = 1)

  # Calculate rank
  R = (rank(X) - 1)/n + 1/(2 * n)

  # 3) Dietary inflammation score | Fractional rank of true proximity
  eps = rnorm(n = n, mean = 0, sd = 1)
  Y = alpha1 + beta1 * R + eps

  # Generate errors and error-prone X
  U = rnorm(n = n, mean = 0, sd = sigmaU)
  Xstar = X + U

  # Calculate error-prone rank
  Rstar = (rank(Xstar) - 1) / n + 1 / (2 * n)
  w = Rstar - R

  # 4) Return
  data.frame(Y, X, R, Z, U, Xstar, Rstar, w, Y_bar = mean(Y))
}

# Simulation to calculate the beta1/beta1* coefficients
simulate_betas = function(sigmaU = 0.5, n = 1000, alpha1 = 2.5, beta1 = -3) {
  data = sim_data(sigmaU = sigmaU,
                     n = n,
                     alpha1 = alpha1,
                     beta1 = beta1)

  rineq1 = ci(ineqvar = data$X, outcome = data$Y, type = "CI", )
  rineq2 = ci(ineqvar = data$Xstar, outcome = data$Y, type = "CI")

  mu_hat = data$Y_bar[1]
  var_r = var(data$R)
  var_w = var(data$Rstar - data$R)

  fit1 = lm(Y ~ R, data = data)
  alpha1_hat = fit1$coefficients[1]
  beta1_hat = fit1$coefficients[2]
  cov_hat = vcov(fit1)

  fit2 = lm(Y ~ Rstar, data = data)
  alpha1star_hat = fit2$coefficients[1]
  beta1star_hat = fit2$coefficients[2]
  covstar_hat = vcov(fit2)

  return(data.frame(alpha1_hat,
                    beta1_hat,
                    alpha1star_hat,
                    beta1star_hat,
                    c_hat = 2 * var_r / mu_hat * beta1_hat,
                    se_c_hat = (2 * var_r / ((alpha1_hat + beta1_hat / 2) ^ 2)) *
                      sqrt(beta1_hat ^ 2 * cov_hat[1, 1] - 2 * alpha1_hat * beta1_hat * cov_hat[1, 2] + alpha1_hat ^ 2 * cov_hat[2, 2]),
                    cstar_hat = 2 * var_r / mu_hat * beta1star_hat,
                    se_cstar_hat = (2 * var_r / ((alpha1star_hat + beta1star_hat / 2) ^ 2)) *
                      sqrt(beta1star_hat ^ 2 * covstar_hat[1, 1] - 2 * alpha1star_hat * beta1star_hat * covstar_hat[1, 2] + alpha1star_hat ^ 2 * covstar_hat[2, 2])
                    )
         )
}
