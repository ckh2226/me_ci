source("~/Documents/GitHub/me_ci/sim_data.R")
source("~/Documents/GitHub/mochi/mochi_estimate.R")

set.seed(24)

pak::pak("sarahlotspeich/mochi")
library(cinch)
library(dplyr)
library(tidyr)
library(auditDesignR)
library(ggplot2)

# Simulate data
# Run multiple simulations at different levels of alpha1 and beta1
df_low_ci <- do.call(rbind, replicate(1000, sim_data(0.5, 1000, approx_ci = -0.5), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = -0.5) # CI ~ -0.5
df_zero_ci <- do.call(rbind, replicate(1000, sim_data(0.5, 1000, approx_ci = 0), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = 0.0) # CI ~ 0.0
df_high_ci <- do.call(rbind, replicate(1000, sim_data(0.5, 1000, approx_ci = 0.5), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = 0.5) #CI ~ 0.5

# Combine simulations from all three settings
# all_df = df_low_ci |>
#   bind_rows(df_zero_ci) |>
#   bind_rows(df_high_ci)

# Divide data into individual samples
groups <- rep(1:1000, each = 1000)
low_ci_samples <- split(df_low_ci, groups)
zero_ci_samples <- split(df_zero_ci, groups)
high_ci_samples <- split(df_high_ci, groups)

all_df <- list(low_ci_samples, zero_ci_samples, high_ci_samples)

# results <- data.frame(ci_naive = rep(NA, 1000),
#                       ci_mb = rep(NA, 1000),
#                       ci_mb_se = rep(NA, 1000),
#                       ci_oracle = rep(NA, 1000),
#                       approx_ci = rep(NA, 1000))

all_results <- list()

# Compute moment-based, naive, oracle, and jackknife standard errors for moment estimate
for(l in 1:3) {
  list <- all_df[[l]]
  results <- data.frame(ci_naive = rep(NA, 1000),
                        ci_mb = rep(NA, 1000),
                        ci_mb_se = rep(NA, 1000),
                        ci_oracle = rep(NA, 1000),
                        approx_ci = rep(NA, 1000))
  
  for(i in 1:length(list)) {
    samp <- list[[i]]
    
    # Calculate moment-based and naive CI using mochi
    est <- mochi(health = samp$Y, unval_exposure = samp$Xstar, val_exposure = samp$Xval, include_se = TRUE, return_naive = TRUE)
    results[i, "ci_naive"] <- est$ci.naive
    results[i, "ci_mb"] <- est$ci.moment
    results[i, "ci_mb_se"] <- est$ci.moment.se
    
    # Calculate true CI
    mu_hat <- mean(samp$Y)
    varR <- var(samp$R)
    fit_ci_oracle <- lm(Y ~ R, data = samp)
    beta1_hat <- fit_ci_oracle$coefficients[2] 
    ci_premult = 2 * varR / mu_hat 
    results[i, "ci_oracle"] <- ci_premult * beta1_hat # Error-free CI
    
    # Add approximate true CI group
    results[i, "approx_ci"] <- unique(samp$approx_ci)
  }
  
  all_results[[l]] <- results
}

# Save as RDS
# saveRDS(all_results, file = "~/Desktop/Thesis/all_results.rds")

# for(i in 1:length(low_ci_samples)) {
#   samp <- low_ci_samples[[i]]
#   
#   # Calculate moment-based and naive CI using mochi
#   est <- mochi(health = samp$Y, unval_exposure = samp$Xstar, val_exposure = samp$Xval, include_se = TRUE, return_naive = TRUE)
#   results[i, "ci_naive"] <- est$ci.naive
#   results[i, "ci_mb"] <- est$ci.moment
#   results[i, "ci_mb_se"] <- est$ci.moment.se
#   
#   # Calculate true CI
#   mu_hat <- mean(samp$Y)
#   varR <- var(samp$R)
#   fit_ci_oracle <- lm(Y ~ R, data = df_low_ci)
#   beta1_hat <- fit_ci_oracle$coefficients[2] 
#   ci_premult = 2 * varR / mu_hat 
#   results[i, "ci_oracle"] <- ci_premult * beta1_hat # Error-free CI
# }

# write.csv(results, "~/Desktop/Thesis/sim_results.csv")

# results <- read.csv("~/Desktop/Thesis/sim_results.csv")

# Build table for bias, ESE, and ASE
sum_tab <- data.frame(bias_naive = rep(NA, 3),
                      bias_perc_naive = rep(NA, 3),
                      bias_mb = rep(NA, 3),
                      bias_perc_mb = rep(NA, 3),
                      ese = rep(NA, 3),
                      ase = rep(NA, 3),
                      truth = rep(NA, 3))

for(l in 1:3) {
  sum_tab[l,"bias_naive"] <- mean(all_results[[l]]$ci_naive - all_results[[l]]$ci_oracle)
  sum_tab[l, "bias_perc_naive"] <- mean((all_results[[l]]$ci_naive - all_results[[l]]$ci_oracle)/all_results[[l]]$ci_oracle) 
  sum_tab[l,"bias_mb"] <- mean(all_results[[l]]$ci_mb - all_results[[l]]$ci_oracle)
  sum_tab[l, "bias_perc_mb"] <- mean((all_results[[l]]$ci_mb - all_results[[l]]$ci_oracle)/all_results[[l]]$ci_oracle)
  sum_tab[l, "ese"] <- sd(all_results[[l]]$ci_mb)
  sum_tab[l, "ase"] <- mean(all_results[[l]]$ci_mb_se)
  sum_tab[l, "truth"] <- unique(all_results[[l]]$approx_ci)
}

library(gt)
sum_tab |>
  gt() |>
  cols_label(
    bias_naive = md("**Bias (Naive)**"),
    bias_perc_naive = md("**Bias (Naive) %**"),
    bias_mb = md("**Bias (MB)**"),
    bias_perc_mb = md("**Bias (%)**"),
    ese = md("**ESE (MB)**"),
    ase = md("**ASE (MB)**"),
    truth = md("**True CI**")) |>
  fmt_number(
    columns = c(bias_naive, bias_mb, ese, ase),
    decimals = 4
  ) |>
  fmt_percent(
    columns = c(bias_perc_naive, bias_perc_mb),
    decimals = 1) |>
  fmt_number(
    columns = c(truth),
    drop_trailing_zeros = TRUE
  )
    

gt::fmt_number(sum_tab, 
           decimals = 5,
           drop_trailing_zeros = TRUE)


bias <- c(bias_mb, bias_perc, ese, ase)

knitr::kable(t(bias),
             digits = 5,
             col.names = c("Moment-Based Bias", "Moment-Based Bias (%)", "ESE", "ASE"))


  # | Moment-Based Bias| Moment-Based Bias (%)|    ESE|     ASE|
  # |-----------------:|---------------------:|------:|-------:|
  # |          -0.00057|               0.11394| 0.0324| 0.03513|

# Repeat for other CI levels and add to table

# For each samples, compute naive, oracle, and mb CI
results <- data.frame(ci_oracle = rep(NA, 1000),
                      ci_naive = rep(NA, 1000),
                      ci_mb = rep(NA, 1000))
for(i in 1:length(low_ci_samples)) {
  # Select i-th sample
  samp <- low_ci_samples[[i]]
  
  # Compute constants
  mu_hat <- mean(samp$Y)
  varR <- var(df_low_ci$R)
  
  # Fit models
  fit_ci_naive <- lm(Y ~ Rstar, data = df_low_ci)
  fit_ci_oracle <- lm(Y ~ R, data = df_low_ci)
  
  # Store coefficients
  beta1star_hat <- fit_ci_naive$coefficients[2] 
  beta1_hat <- fit_ci_oracle$coefficients[2] 
  
  # Compute concentration indeicies
  ci_premult = 2 * varR / mu_hat 
  results[i, "ci_oracle"] <- ci_premult * beta1_hat # Error-free CI
  results[i, "ci_naive"] <- ci_premult * beta1star_hat # Error-prone CI
  results[i, "ci_mb"] <- mochi(health = samp$Y, unval_exposure = samp$Xstar, val_exposure = samp$Xval)
}

# Pivot for plotting
results_new <- results |>
  pivot_longer(cols = c(ci_oracle, ci_naive, ci_mb),
               names_to = "Estimate")

ggplot(data = results_new, aes(x = "", y = value, fill = Estimate)) +
  geom_boxplot() +
  theme_bw() + 
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkgrey") +
  labs(title = "Distribution of Simulated Concentration Indices",
       x = "n = 1000",
       y = "Concentration Index") +
  scale_fill_manual(values = c("#725def", "#dd227d", "#ffb10e"), 
                               breaks = c("ci_mb", "ci_naive", "ci_oracle"),
                               labels = c("Moment-Based CI", "Naive CI", "Oracle CI")) +
  annotate("text", x = 0.5, y = -0.495, label = "True CI = -0.5", size = 3)
  # scale_fill_manual(values = c("#82A641", "#9ABC59", "#A8CC64", "#B4D17D", 
                               # "#F2D68F", "#FFD061", "#EBB940", "#C49525", 
                               # "#386092", "#3D6CA6", "#5981B5", "#6892C4")





