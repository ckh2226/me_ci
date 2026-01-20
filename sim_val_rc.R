source("~/Documents/me_ci/sim_data.R")

# Be reproducible queens. 
set.seed(120)

# Load packages
library(ggplot2)
library(dplyr)

# Write function to partially validate and try regression calibration
sim_val_rc = function(sigmaU, n, alpha1, beta1, pv = 0.1) {
  # Simulate data 
  dat = sim_data(sigmaU, n, alpha1, beta1)
  
  # Partially validate 
  dat$V = sample(x = c(FALSE, TRUE), 
                 size = 1000, 
                 replace = TRUE, 
                 prob = c(1 - pv, pv))
  dat$Xval = dat$X ## initialize Xval = X 
  dat$Xval[!dat$V] = NA ## but then redact Xval if V = FALSE (unvalidated)

  # Fit a calibration model to the subset with V = TRUE 
  ## X | X* (validated | unvalidated + other things that might be helpful)
  fit_cal = lm(formula = Xval ~ Xstar, ## X | X* 
               data = dat)
  ## Predict X from X* for everyone 
  dat$Xcal = predict(object = fit_cal, 
                     newdata = dat)

  # But... we need to get ranks now. 
  dat$Rcal = (rank(dat$Xcal) - 1) / 1000 + 1 / (2 * 1000)
  
  # Now, we can try to fit the concentration index 
  mu_hat <- mean(dat$Y)
  var_R <- var(dat$R) ## = var(R*) = var(Rcal) -- Nice! 
  fit_ci_cal <- lm(Y ~ Rcal, data = dat)
  fit_ci_naive <- lm(Y ~ Rstar, data = dat)
  fit_ci_oracle <- lm(Y ~ R, data = dat)
  beta1cal_hat <- fit_ci_cal$coefficients[2] 
  beta1star_hat <- fit_ci_naive$coefficients[2] 
  beta1_hat <- fit_ci_oracle$coefficients[2] 
  
  ci_premult = 2 * var_R/mu_hat 
  ci_xcal <- ci_premult * beta1cal_hat # Regression calibration CI
  ci_xstar <- ci_premult * beta1star_hat # Error-prone CI
  ci_x <- ci_premult * beta1_hat # Error-free CI
  return(c(ci_x, ci_xcal, ci_xstar))
}

# Run multiple simulations at different levels of alpha1 and beta1
df_low_ci <- do.call(rbind, replicate(10000, sim_val_rc(0.5, 1000, 2.5, -3), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = -0.5) # CI ~ -0.5
df_zero_ci <- do.call(rbind, replicate(10000, sim_val_rc(0.5, 1000, 3, 0), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = 0.0) # CI ~ 0.0
df_high_ci <- do.call(rbind, replicate(10000, sim_val_rc(0.5, 1000, -0.5, 3), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = 0.5) #CI ~ 0.5

library(ggplot2)
dat |> 
  ggplot() + 
  geom_point(aes(x = Xstar, y = Xval), color = "cornflowerblue", alpha = 0.5) + 
  geom_point(aes(x = Xstar, y = Xcal), color = "magenta", alpha = 0.5)

# Normally, you would then fit Y ~ Xcal 
fit_cal_outcome1 = lm(formula = Y ~ Xcal, 
                      data = dat)
summary(fit_cal_outcome1)
fit_naive_outcome1 = lm(formula = Y ~ Xstar, 
                        data = dat)
summary(fit_naive_outcome1)
fit_oracle_outcome1 = lm(formula = Y ~ X, 
                         data = dat)
summary(fit_oracle_outcome1)

dat |> 
  ggplot() + 
  geom_point(aes(x = Rstar, y = R), color = "cornflowerblue", alpha = 0.5) + 
  geom_point(aes(x = Rstar, y = Rcal), color = "magenta", alpha = 0.5)