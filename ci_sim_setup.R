## Set sample size 
N = 100 

## Define variance of the fractional rank 
sigma2R = var(1:N / N) ### treated as a constant, only dictated by N

## Function to simulate data 
sim = function(CI = -0.5, beta1 = -0.1, use_rineq = FALSE, return_data = FALSE) {
  ### Simulate METRO status (yes/no)
  METRO = rbinom(n = N, 
                 size = 1, 
                 prob = 0.15)
  
  ### Simulate INCOME | METRO
  INCOME = rgamma(n = N, 
                  shape = 1.75 + 0.25 * METRO, 
                  scale = 5)
  
  ### Combine them 
  dat = data.frame(METRO, INCOME) |> 
    dplyr::arrange(INCOME) #### Order by INCOME (ascending)
  
  ### Define fractional ranks based on INCOME 
  #### Adapted from the rineq::rank_wt() function
  dat$R = c(0, 
            cumsum(x = rep(x = 1 / N, 
                           times = (N - 1)))) + 1 / (2 * N)

  ### Define intercept in model based on CI and beta1
  alpha1 = (2 * sigma2R * beta1 / CI) - (beta1 / 2)
  
  ### Simulate ACCESS | RANK under (8.10) in Health Equity Ch 8
  dat = dat |> 
    dplyr::mutate(ACCESS = rnorm(n = N, 
                                 mean = alpha1 + beta1 * R, 
                                 sd = 1 / 4))
  
  if (use_rineq) {
    ### Calculate from the covariance formula 
    ci_eq8.3 = ci(ineqvar = dat$INCOME, 
                  outcome = dat$ACCESS, 
                  method = "cov_convenience", 
                  df_correction = FALSE)$concentration_index
    
    ### Fit the transformed model 
    ci_eq8.7 = ci(ineqvar = dat$INCOME, 
                  outcome = dat$ACCESS, 
                  method = "linreg_convenience", 
                  df_correction = FALSE)$concentration_index
    
    ### Fit the untransformed model 
    ci_eq8.11 = ci(ineqvar = dat$INCOME, 
                   outcome = dat$ACCESS, 
                   method = "linreg_delta", 
                   df_correction = FALSE)$concentration_index
  } else {
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
  }
  
  ### Return results 
  if (return_data) {
    list(ci_est = data.frame(ci_eq8.3, ci_eq8.7, ci_eq8.11), 
         data = dat)
  } else {
    data.frame(ci_eq8.3, ci_eq8.7, ci_eq8.11)  
  }
}

## Try it out
### My code 
set.seed(212) #### Be reproducible
colMeans(do.call(what = rbind, 
                 args = replicate(n = 1000, 
                                  expr = sim(CI = 0.9, beta1 = 0.2), 
                                  simplify = FALSE)))
### The rineq package
set.seed(212) #### Be reproducible
rineq_reps = do.call(what = rbind, 
                     args = replicate(n = 1000, 
                                      expr = sim(CI = 0.9, beta1 = 0.2, use_rineq = TRUE), 
                                      simplify = FALSE))
colMeans(rineq_reps)

library(LorenzRegression)
temp = sim(CI = -0.9, beta1 = -10, return_data = TRUE)$data
Lorenz.curve(y = temp$ACCESS, ## farthest --> nearest proximity 
             x = temp$INCOME, ## worse --> better income
             graph = TRUE)

## Try different alpha1, beta1 
heavy_neg_params = expand.grid(alpha1 = seq(1, 2, by = 0.1), 
                               beta1 = -0.5) |>
  dplyr::mutate(CI = NA, minACCESS = NA, medACCESS = NA, maxACCESS = NA)

# try_params = expand.grid(alpha1 = seq(-3, 3), 
#                          beta1 = seq(-0.5, 0.5, by = 0.1)) |> 
#   dplyr::mutate(CI = NA, minACCESS = NA, medACCESS = NA, maxACCESS = NA)
for (r in 1:nrow(heavy_neg_params)) {
  sim_ci = ## Try it out
    colMeans(do.call(what = rbind, 
                     args = replicate(n = 1000, 
                                      expr = sim(alpha1 = heavy_neg_params$alpha1[r],
                                                 beta1 = heavy_neg_params$beta1[r]), 
                                      simplify = FALSE)))
  heavy_neg_params$CI[r] = sim_ci[1]
  
  ## Plot the corresponding CI for one 
  temp = sim(alpha1 = heavy_neg_params$alpha1[r],
             beta1 = heavy_neg_params$beta1[r], 
             return_data = TRUE)$data
  heavy_neg_params$minACCESS[r] = min(temp$ACCESS)
  heavy_neg_params$medACCESS[r] = median(temp$ACCESS)
  heavy_neg_params$maxACCESS[r] = max(temp$ACCESS)
}

## Plot the corresponding CI for one 
temp = sim(alpha1 = 1.6, beta1 = -0.5, return_data = TRUE)$data |> 
  dplyr::mutate(c.ACCESS = cumsum(ACCESS) / sum(ACCESS))

# #calculate cumulative proportion of distance traveled to healthy foods 
# ## --> all routine appts/number of individuals (not sure if this is right)
# conc$c.PROXIMITY <- cumsum(conc$PROXIMITY) / sum(conc$PROXIMITY)
# 
# # cumulative distribution of median income (percentile rank from 0 to 1)
# library(magrittr) ## needed for %$%
# conc$c.INCOME <- conc %$% cume_dist(INCOME)
# 
# 

library(ggplot2)
data.frame(x = 1:N/N) |> 
  ggplot(aes(x = x)) + 
  stat_function(fun = function(x) x ^ (1 / 4)) +
  stat_function(fun = function(x) x ^ (1 / 3)) +
  stat_function(fun = function(x) x ^ (1 / 2)) +
  stat_function(fun = function(x) x ^ 2) +
  stat_function(fun = function(x) x ^ 3) +
  stat_function(fun = function(x) x ^ 4) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  labs(x = "Fractional Rank by Outcome", 
       y = "Cumulative Proportion of Distance to Healthy Foods")

beta1 * (2 * sigma2R) / (alpha1 + (beta1 / 2)) ## CHECK: = CI (YES)
ci(ineqvar = INCOME, outcome = ACCESS, method = "linreg_delta")
ci(ineqvar = INCOME, outcome = ACCESS, method = "cov_convenience")