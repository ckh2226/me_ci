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
  
  return(c(ci_x = as.numeric(ci_x), 
           ci_xmb_varRval = as.numeric(ci_xmb_varRval), 
           ci_xmb_varRstar = as.numeric(ci_xmb_varRstar), 
           ci_xstar = as.numeric(ci_xstar), 
           lambdahat, lambdahat_varRval, lambdahat_varRstar, 
           varRval = varRval, varWval = varWval, covRWval = covRWval,
           varR = varR, varW = varW, covRW = covRW))
}

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