source("~/Documents/me_ci/sim_data.R")

# Be reproducible queens. 
set.seed(120)

# Load packages
library(ggplot2)
library(dplyr)

# Write function to partially validate and try regression calibration
sim_val_mb = function(sigmaU, n, alpha1, beta1, pv = 0.1) {
  # Simulate data 
  dat = sim_data(sigmaU, n, alpha1, beta1)
  
  # Partially validate 
  dat$V = sample(x = c(FALSE, TRUE), 
                 size = 1000, 
                 replace = TRUE, 
                 prob = c(1 - pv, pv))
  dat$Xval = dat$X ## initialize Xval = X 
  dat$Xval[!dat$V] = NA ## but then redact Xval if V = FALSE (unvalidated)
  dat$Wval = dat$w ## initialize Wval = w 
  dat$Wval[!dat$V] = NA ## but then redact Xval if V = FALSE (unvalidated)
  
  # Calculate ranks based on X in the validation subsample
  nv = sum(dat$V) ## sample size validated
  dat$Rval = NA ## initialize as NA
  dat$Rval[dat$V] = (rank(dat$Xval[dat$V]) - 1) / nv + 1 / (2 * nv)
  
  # Use validation subset to estimate quantities in the bias factor
  varRval = var(dat$Rval, na.rm = TRUE) ## Var(R)
  varWval = var(dat$Wval, na.rm = TRUE) ## Var(W)
  covRWval = cov(dat$Rval, dat$Wval, use = "complete.obs") ## Cov(R,W)
  lambdahat_varRval = (varRval + covRWval) / 
    (varRval + varWval + 2 * covRWval) ## Estimated bias factor
  
  # Now, we can try to fit the concentration index 
  ## Full-sample methods (oracle and naive)
  mu_hat <- mean(dat$Y)
  var_R <- var(dat$R) ## = var(R*) -- Nice! 
  fit_ci_naive <- lm(Y ~ Rstar, data = dat)
  fit_ci_oracle <- lm(Y ~ R, data = dat)
  beta1star_hat <- fit_ci_naive$coefficients[2] 
  beta1_hat <- fit_ci_oracle$coefficients[2] 
  ci_premult = 2 * var_R / mu_hat 
  ci_xstar <- ci_premult * beta1star_hat # Error-prone CI
  ci_x <- ci_premult * beta1_hat # Error-free CI
  
  ## Partially validated methods 
  ### Using Var(R) from validation subsample
  ci_xmb_varRval <- lambdahat_varRval * ci_xstar # Moment-corrected CI using subsample Var(R)
  varRstar <- var(dat$Rstar) 
  ### Using Var(R*) from the the sample
  lambdahat_varRstar = (varRstar + covRWval) / 
    (varRstar + varWval + 2 * covRWval) ## Estimated bias factor
  ci_xmb_varRstar <- lambdahat_varRstar * ci_xstar # Moment-corrected CI using subsample Var(R)
  
  return(c(ci_x = as.numeric(ci_x), 
           ci_xmb_varRval = as.numeric(ci_xmb_varRval), 
           ci_xmb_varRstar = as.numeric(ci_xmb_varRstar), 
           ci_xstar = as.numeric(ci_xstar), 
           varRval = varRval, varWval = varWval, covRWval = covRWval,
           varR = var_R, varW = var(dat$w), covRW = cov(dat$R, dat$w)))
}

# Run multiple simulations at different levels of alpha1 and beta1
df_low_ci <- do.call(rbind, replicate(1000, sim_val_mb(0.5, 1000, 2.5, -3), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = -0.5) # CI ~ -0.5
df_zero_ci <- do.call(rbind, replicate(1000, sim_val_mb(0.5, 1000, 3, 0), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = 0.0) # CI ~ 0.0
df_high_ci <- do.call(rbind, replicate(1000, sim_val_mb(0.5, 1000, -0.5, 3), simplify = FALSE)) |>
  data.frame() |> 
  mutate(approx_ci = 0.5) #CI ~ 0.5

library(ggplot2)
all_df = df_low_ci |> 
  bind_rows(df_zero_ci) |> 
  bind_rows(df_high_ci) 
all_df |> 
  dplyr::select(dplyr::starts_with("ci_"), approx_ci) |>
  tidyr::gather(key = "Method", value = "Estimate", -5) |> 
  ggplot(aes(x = Method, y = Estimate, fill = Method)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = approx_ci), 
             linetype = "dashed") + 
  facet_wrap(~approx_ci, scales = "free")

all_df |> 
  dplyr::select(dplyr::starts_with(c("var", "cov"))) |>
  tidyr::gather(key = "Method", value = "Estimate") |> 
  dplyr::mutate(validated = grepl(pattern = "val", x = Method), 
                Method = sub(pattern = "val", replacement = "", x = Method)) |> 
  ggplot(aes(x = validated, y = Estimate, fill = validated)) + 
  geom_boxplot() + 
  facet_wrap(~Method, scales = "free")

# dat |> 
#   ggplot() + 
#   geom_point(aes(x = Xstar, y = Xval), color = "cornflowerblue", alpha = 0.5) + 
#   geom_point(aes(x = Xstar, y = Xcal), color = "magenta", alpha = 0.5)