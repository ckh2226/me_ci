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
  varR <- var(samp$R)
  
  # Fit models
  fit_ci_naive <- lm(Y ~ Rstar, data = samp)
  fit_ci_oracle <- lm(Y ~ R, data = samp)
  
  # Store coefficients
  beta1star_hat <- fit_ci_naive$coefficients[2] 
  beta1_hat <- fit_ci_oracle$coefficients[2] 
  
  # Compute concentration indices
  ci_premult = 2 * varR / mu_hat 
  results[i, "ci_oracle"] <- ci_premult * beta1_hat # Error-free CI
  results[i, "ci_naive"] <- ci_premult * beta1star_hat # Error-prone CI
  results[i, "ci_mb"] <- mochi(health = samp$Y, unval_exposure = samp$Xstar, val_exposure = samp$Xval)
}

# Pivot for plotting
results_new <- results |>
  pivot_longer(cols = c(ci_oracle, ci_naive, ci_mb),
               names_to = "Estimator")

ggplot(data = results_new, aes(x = "", y = value, fill = Estimator)) +
  geom_boxplot() +
  theme_bw() + 
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkgrey") +
  labs(title = "Distribution of Simulated Concentration Indices",
       x = "n = 1000",
       y = "Concentration Index",
       fill = "CI Estimator") +
  scale_fill_manual(values = c("#82A641", "#EBB940", "#5981B5"), 
                    labels = c("Moment-Based", "Naive", "Oracle")) +
  annotate("text", x = 0.515, y = -0.495, label = "True CI = -0.5", size = 3) 
  # scale_fill_manual(values = c("#82A641", "#9ABC59", "#A8CC64", "#B4D17D", 
                               # "#F2D68F", "#FFD061", "#EBB940", "#C49525", 
                               # "#386092", "#3D6CA6", "#5981B5", "#6892C4"))

# Compute jackknife sample and standard error using sample of n = 1000 observations
one_samp <- low_ci_samples[[1]]
n <- nrow(one_samp)[1]
jack_ci_mb <- rep(NA, n)
for(i in 1:n) {
  temp_samp <- one_samp[-i, ]
  # Compute constants
  mu_hat <- mean(temp_samp$Y)
  varR <- var(temp_samp$R)
  
  # Fit models
  fit_ci_naive <- lm(Y ~ Rstar, data = temp_samp)
  fit_ci_oracle <- lm(Y ~ R, data = temp_samp)
  
  # Store coefficients
  beta1star_hat <- fit_ci_naive$coefficients[2] 
  beta1_hat <- fit_ci_oracle$coefficients[2] 
  
  # Compute moment-based concentration index
  jack_ci_mb[i] <- mochi(health = temp_samp$Y, 
                                    unval_exposure = temp_samp$Xstar, 
                                    val_exposure = temp_samp$Xval)
}

# Distribution of jackknife sample
ggplot(aes(x = jack_ci_mb), data = as.data.frame(jack_ci_mb)) +
  geom_histogram(fill = "#5981B5") +
  labs(title = "Distribution of n = 1000 Jackknife Moment-Based Estimators",
       x = "Jackknife Moment-Based Concentration Indices",
       y = "Count") +
  theme_bw()

# Compute jackknife standard error
jack_ci_mb <- as.numeric(jack_ci_mb)
se_jack <- sqrt((n - 1)/n * sum((jack_ci_mb - mean(jack_ci_mb))^2))

# Compute jackknife confidence interval
c(quantile(jack_ci_mb, 0.025), quantile(jack_ci_mb, 0.975))





