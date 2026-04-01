source("~/Documents/GitHub/me_ci/sim_data.R")

set.seed(24)

pak::pak("sarahlotspeich/cinch")
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
df_zero_ci <- do.call(rbind, replicate(10000, sim_data(0.5, 1000, approx_ci = 0), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = 0.0) # CI ~ 0.0
df_high_ci <- do.call(rbind, replicate(10000, sim_data(0.5, 1000, approx_ci = 0.5), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = 0.5) #CI ~ 0.5

# Combine simulations from all three settings
all_df = df_low_ci |>
  bind_rows(df_zero_ci) |>
  bind_rows(df_high_ci)

# Divide data into individual samples
groups <- rep(1:1000, each = 1000)
low_ci_samples <- split(df_low_ci, groups)

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
  scale_fill_manual(values = c("#82A641", "#EBB940", "#5981B5")) +
  annotate("text", x = 0.5, y = -0.495, label = "True CI = -0.5", size = 3)
  # scale_fill_manual(values = c("#82A641", "#9ABC59", "#A8CC64", "#B4D17D", 
                               # "#F2D68F", "#FFD061", "#EBB940", "#C49525", 
                               # "#386092", "#3D6CA6", "#5981B5", "#6892C4"))





