## Load libraries
library(LorenzRegression) ### for quick concentration curves

## Set sample size 
N = 100 

## Be reproducible
set.seed(212)

## Define variance of the fractional rank 
sigma2R = var(1:N / N) ### treated as a constant, only dictated by N

## Function to simulate data 
sim = function(alpha1 = -1, beta1 = 0.1, return_data = FALSE) {
  ### Simulate METRO status (yes/no)
  METRO = rbinom(n = N, 
                 size = 1, 
                 prob = 0.15)
  
  ### Simulate INCOME | METRO
  INCOME = rgamma(n = N, shape = 1.75 + 0.25 * METRO, scale = 5)
  
  ### Combine them 
  dat = data.frame(METRO, INCOME) |> 
    dplyr::arrange(INCOME) |> 
    dplyr::mutate(R = 1:dplyr::n() / dplyr::n())
  
  ### Define fractional ranks based on INCOME 
  # R = order(INCOME) / N #### I THINK THIS IS THE PROBLEM
  
  ### Define intercept in model based on CI and beta1
  # alpha1 = (2 * sigma2R * beta1 / CI) - (beta1 / 2)
  
  ### Simulate ACCESS | RANK under (8.10) in Health Equity Ch 8
  dat = dat |> 
    dplyr::mutate(ACCESS = rlnorm(n = N, 
                                  meanlog = alpha1 + beta1 * INCOME, 
                                  sdlog = 1 / 8))
  
  ### Calculate from the covariance formula 
  ci_eq8.3 = with(dat, (2 / mean(ACCESS)) * cov(ACCESS, R)) #### concentration index
  
  ### Fit the transformed model 
  dat = dat |> 
    dplyr::mutate(T_ACCESS = 2 * sigma2R * ACCESS / mean(ACCESS))
  fit_eq8_7 = lm(formula = T_ACCESS ~ R, 
                 data = dat)
  ci_eq8.7 = fit_eq8_7$coefficients[2]
  
  ### Fit the untransformed model 
  fit_eq8_10 = lm(formula = ACCESS ~ R, 
                  data = dat)
  alpha1_hat = fit_eq8_10$coefficients[1] #### intercept
  beta1_hat = fit_eq8_10$coefficients[2] #### slope on rank
  ci_eq8.11 = beta1_hat * (2 * sigma2R) / (alpha1_hat + beta1_hat / 2) #### concentration index
  
  ### Return results 
  if (return_data) {
    list(ci_est = data.frame(ci_eq8.3, ci_eq8.7, ci_eq8.11), 
         data = dat)
  } else {
    data.frame(ci_eq8.3, ci_eq8.7, ci_eq8.11)  
  }
}

## Try all the different combinations of alpha1, beta1
try_params = function(param_combos) {
  for (r in 1:nrow(param_combos)) {
    sim_ci = ## Try it out
      colMeans(do.call(what = rbind, 
                       args = replicate(n = 1000, 
                                        expr = sim(alpha1 = param_combos$alpha1[r],
                                                   beta1 = param_combos$beta1[r]), 
                                        simplify = FALSE)))
    param_combos$CI[r] = sim_ci[1]
    
    ## Plot the corresponding CI for one 
    temp = sim(alpha1 = param_combos$alpha1[r],
               beta1 = param_combos$beta1[r], 
               return_data = TRUE)$data
    param_combos$minACCESS[r] = min(temp$ACCESS)
    param_combos$medACCESS[r] = median(temp$ACCESS)
    param_combos$maxACCESS[r] = max(temp$ACCESS)
  }
  return(param_combos)
}

## Pick an alpha1, beta1 for a heavy disparity among low-income
heavy_neg_params = expand.grid(alpha1 = seq(2.5, 3.5, by = 0.1), 
                               beta1 = seq(-0.5, -0.4, by = 0.1)) |>
  dplyr::mutate(CI = NA, minACCESS = NA, medACCESS = NA, maxACCESS = NA) |> 
  try_params()

### CHOSEN: alpha1 = 3.5, beta1 = -0.5 --> CI = -0.72, on avg, with 0 < ACCESS < 27 (1.8 avg)
heavy_neg_params[11, ]

### Plot the corresponding CI for one 
temp = sim(alpha1 = 3.5, beta1 = -0.5, return_data = TRUE)$data
Lorenz.curve(y = temp$ACCESS, ## farthest --> nearest proximity 
             x = temp$INCOME, ## worse --> better income
             graph = TRUE)

## Pick an alpha1, beta1 for a moderate disparity among low-income
mod_neg_params = expand.grid(alpha1 = seq(1, 5, by = 0.25), 
                             beta1 = seq(-0.3, -0.2, by = 0.1)) |>
  dplyr::mutate(CI = NA, minACCESS = NA, medACCESS = NA, maxACCESS = NA) |> 
  try_params()

### CHOSEN: alpha1 = 3.25, beta1 = -0.25 --> CI = -0.46, on avg, with 0 < ACCESS < 27 (5.7 avg)
mod_neg_params[27, ]

### Plot the corresponding CI for one 
temp = sim(alpha1 = 3.25, beta1 = -0.2, return_data = TRUE)$data
Lorenz.curve(y = temp$ACCESS, ## farthest --> nearest proximity 
             x = temp$INCOME, ## worse --> better income
             graph = TRUE)

## Pick an alpha1, beta1 for a light disparity among low-income
light_neg_params = expand.grid(alpha1 = seq(1, 5, by = 0.25), 
                             beta1 = seq(-0.1, 0, by = 0.05)) |>
  dplyr::mutate(CI = NA, minACCESS = NA, medACCESS = NA, maxACCESS = NA) |> 
  try_params()

### CHOSEN: alpha1 = 3.25, beta1 = -0.25 --> CI = -0.46, on avg, with 0 < ACCESS < 27 (5.7 avg)
mod_neg_params[27, ]

### Plot the corresponding CI for one 
temp = sim(alpha1 = 3.25, beta1 = -0.2, return_data = TRUE)$data
Lorenz.curve(y = temp$ACCESS, ## farthest --> nearest proximity 
             x = temp$INCOME, ## worse --> better income
             graph = TRUE)

