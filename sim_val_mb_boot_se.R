source("~/Documents/me_ci/sim_data.R")

# Be reproducible queens. 
set.seed(120)

# Load packages
library(ggplot2) ## for pretty plots
library(dplyr) ## for data manipulation
library(tidyr) ## for data transformation

# Write function to partially validate and try regression calibration
sim_val_mb = function(sigmaU, n, approx_ci, pv = 0.1) {
  # Simulate data 
  dat = sim_data(sigmaU, n, approx_ci, pv)
  
  # Use validation subset to estimate quantities in the bias factor
  varRval = var(dat$Rval, na.rm = TRUE) ## Var(R)
  varWval = var(dat$Wval, na.rm = TRUE) ## Var(W)
  covRWval = cov(dat$Rval, dat$Wval, use = "complete.obs") ## Cov(R,W)
  lambdahat_varRval = (varRval + covRWval) / 
    (varRval + varWval + 2 * covRWval) ## Estimated bias factor
  
  # Now, we can try to fit the concentration index 
  ## Full-sample methods (oracle and naive)
  mu_hat <- mean(dat$Y)
  varR <- var(dat$R) ## = var(R*) -- Nice! 
  fit_ci_naive <- lm(Y ~ Rstar, data = dat)
  fit_ci_oracle <- lm(Y ~ R, data = dat)
  beta1star_hat <- fit_ci_naive$coefficients[2] 
  beta1_hat <- fit_ci_oracle$coefficients[2] 
  ci_premult = 2 * varR / mu_hat 
  ci_xstar <- ci_premult * beta1star_hat # Error-prone CI
  ci_x <- ci_premult * beta1_hat # Error-free CI
  
  ## Partially validated methods 
  ### Using Var(R) from validation subsample
  ci_xmb_varRval <-  ci_xstar / lambdahat_varRval # Moment-corrected CI using subsample Var(R)
  varRstar <- var(dat$Rstar) ## Var(R*)
  ### Using Var(R*) from the the sample
  lambdahat_varRstar = (varRstar + covRWval) / 
    (varRstar + varWval + 2 * covRWval) ## Estimated bias factor
  ci_xmb_varRstar <- ci_xstar / lambdahat_varRstar # Moment-corrected CI using subsample Var(R)
  
  ## Add oracle for reference 
  covRW = cov(dat$R, dat$W)
  varW = var(dat$W)
  lambdahat = (varR + covRW) / 
    (varR + varW + 2 * covRW) ## Estimated bias factor
  
  ## Bootstrap resample and re-estimate moment-based CI 
  for (b in 1:B) { ### let B be a large number (like 10000) 
    ### Resampled with replacement (ignores validation status)
    boot_rows = sample(x = 1:nrow(dat), size = nrow(dat), replace = TRUE) ### rows to sample
    boot_dat = dat[boot_rows, ] #### ensure you have the same number validated/unvalidated 
    ### Resampled with replacement (holds validation percentage equal)
    boot_rows = c(sample(x = which(!is.na(dat$Rval)), size = sum(!is.na(dat$Rval)), replace = TRUE), ### rows to sample (from validated)
                  sample()) ### rows to sample (from unvalidated)
    boot_dat = dat[boot_rows, ] #### ensure you have the same number validated/unvalidated 
    ### Copy + paste and re-estimate ci_xmb_varRstar 
  }
  ## Summarize bootstrapped estimates 
  boot_se = ## sd(boot strapped estimates)
  boot_lb = ## 2.5th percentile of bootstrapped estimates (for a 95% conf interval)
  boot_ub = ## 97.5th percentile of bootstrapped estimates (for a 95% conf interval)
  
  # Formula for SEs for the oracle / naive   
  # se_c_hat = (2 * var_r / ((beta0_hat + beta1_hat / 2) ^ 2)) *
  #   sqrt(beta1_hat ^ 2 * cov_hat[1, 1] - 2 * beta0_hat * beta1_hat * cov_hat[1, 2] + beta0_hat ^ 2 * cov_hat[2, 2]),
  # se_cstar_hat = (2 * var_r / ((beta0star_hat + beta1star_hat / 2) ^ 2)) *
  #   sqrt(beta1star_hat ^ 2 * covstar_hat[1, 1] - 2 * beta0star_hat * beta1star_hat * covstar_hat[1, 2] + beta0star_hat ^ 2 * covstar_hat[2, 2])
  return(c(ci_x = as.numeric(ci_x), 
           ci_xmb_varRval = as.numeric(ci_xmb_varRval), 
           ci_xmb_varRstar = as.numeric(ci_xmb_varRstar), 
           ci_xstar = as.numeric(ci_xstar), 
           lambdahat, lambdahat_varRval, lambdahat_varRstar, 
           varRval = varRval, varWval = varWval, covRWval = covRWval,
           varR = varR, varW = varW, covRW = covRW))
}

psi_theta_mat = function(data) {
  function(theta) {
    ## Separate elements of parameter vector to align with equations 
    beta0star = theta[1]
    beta1star = theta[2]
    muX = theta[3]
    sigma2X = theta[4]
    sigma2U = theta[5]
    
    with(data, 
         c(Y - beta0star - beta1star * Xstar, ### equation 1 (naive intercept)
           (Y - beta0star - beta1star * Xstar) * Xstar, ### equation 2 (naive slope)
           (X - muX) * V, ### equation 3 (mean of X)
           ((X - muX) ^ 2 - sigma2X) * V, ### equation 4 (variance of X)
           ((U - 0) ^ 2 - sigma2U) * V)) ### equation 5 (variance of U), assume E(U) = 0 
    
    # ## First two elements (estimating equations for naive OLS model)
    # ### Uses all observations
    # r = data$Y - beta0star - beta1star * data$Xstar ### residuals 
    # g1 = r ### equation 1 (naive intercept)
    # g2 = r * data$Xstar ### equation 2 (naive slope)
    # 
    # ## Remaining elements (estimating equations based on moments of X/U)
    # ### Uses only validated observations
    # g3 = (data$X - muX) * data$V ### equation 3 (mean of X)
    # g4 = ((data$X - muX) ^ 2 - sigma2X) * data$V ### equation 4 (variance of X)
    # g5 = ((data$U - 0) ^ 2 - sigma2U) * data$V ### equation 5 (variance of U), assume E(U) = 0
    # 
    # ## Replace NA elements (for unvalidated observations) with 0s
    # g3[is.na(g3)] = 0
    # g4[is.na(g4)] = 0
    # g5[is.na(g5)] = 0
    # 
    # ## Return as a vector 
    # return(c(g1, g2, g3, g4, g5))
  }
}

library(geex)
dat = sim_data(0.5, 1000, -0.5, 3)
m_estimate(
  estFUN = psi_theta_mat,
  data  = temp,
  root_control = setup_root_control(start = c(1,1,1,1,1)))

# Run multiple simulations at different levels of alpha1 and beta1
df_low_ci <- do.call(rbind, replicate(10000, sim_val_mb(0.5, 1000, approx_ci = -0.5), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = -0.5) # CI ~ -0.5
df_zero_ci <- do.call(rbind, replicate(10000, sim_val_mb(0.5, 1000, approx_ci = 0), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = 0.0) # CI ~ 0.0
df_high_ci <- do.call(rbind, replicate(10000, sim_val_mb(0.5, 1000, approx_ci = 0.5), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = 0.5) #CI ~ 0.5

# Combine simulations from all three settings 
all_df = df_low_ci |> 
  bind_rows(df_zero_ci) |> 
  bind_rows(df_high_ci) 

# Make a boxplot of CI estimates by method x setting
all_df |> 
  select(starts_with("ci_"), approx_ci) |> ## only pull Ci vals and settings
  gather(key = "Method", value = "Estimate", -5) |> ## pivot from wide --> long 
  ggplot(aes(x = Method, y = Estimate, fill = Method)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = approx_ci), 
             linetype = "dashed") + 
  facet_wrap(~approx_ci, scales = "free")

# Make a boxplot of var/covar estimates by method (ignored setting because they shouldn't vary)
all_df |> 
  select(starts_with(c("var", "cov"))) |> ## only pull Var(R), Var(W), Cov(R,w) and settings
  gather(key = "Method", value = "Estimate") |> ## pivot from wide --> long 
  mutate(validated = grepl(pattern = "val", x = Method), 
         Method = sub(pattern = "val", replacement = "", x = Method)) |> 
  ggplot(aes(x = validated, y = Estimate, fill = validated)) + 
  geom_boxplot() + 
  facet_wrap(~Method, scales = "free")