# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)

# Source the function sim_data() to simulate data
#source("https://raw.githubusercontent.com/ckh2226/me_ci/refs/heads/main/sim_data.R")
source("~/Documents/me_ci/sim_data.R")

# Simulation to calculate the beta1/beta1* coefficients
sim_rep = function(sigmaU = 1, N = 1000, beta0 = 2.5, beta1 = -3) {
  data = sim_data(sigmaU = sigmaU,
                  n = N,
                  beta0 = beta0,
                  beta1 = beta1)
  data$W = data$Rstar - data$R

  mu_hat = data$Y_bar[1]
  var_r = var(data$R)
  var_w = var(data$W)
  cov_rw = cov(x = data$R, y = data$W)

  fit1 = lm(Y ~ R, data = data)
  beta0_hat = fit1$coefficients[1]
  beta1_hat = fit1$coefficients[2]
  cov_hat = vcov(fit1)

  fit2 = lm(Y ~ Rstar, data = data)
  beta0star_hat = fit2$coefficients[1]
  beta1star_hat = fit2$coefficients[2]
  covstar_hat = vcov(fit2)

  return(data.frame(beta0_hat,
                    beta1_hat,
                    beta0star_hat,
                    beta1star_hat,
                    var_r,
                    var_w,
                    cov_rw,
                    c_hat = 2 * var_r / mu_hat * beta1_hat,
                    se_c_hat = (2 * var_r / ((beta0_hat + beta1_hat / 2) ^ 2)) *
                      sqrt(beta1_hat ^ 2 * cov_hat[1, 1] - 2 * beta0_hat * beta1_hat * cov_hat[1, 2] + beta0_hat ^ 2 * cov_hat[2, 2]),
                    cstar_hat = 2 * var_r / mu_hat * beta1star_hat,
                    se_cstar_hat = (2 * var_r / ((beta0star_hat + beta1star_hat / 2) ^ 2)) *
                      sqrt(beta1star_hat ^ 2 * covstar_hat[1, 1] - 2 * beta0star_hat * beta1star_hat * covstar_hat[1, 2] + beta0star_hat ^ 2 * covstar_hat[2, 2])
                    )
         )
}

# Sett 1: C = -0.5 --------------------------------------------------------------
set.seed(11182025) ## For reproducibility
res_sigmaU_100 = do.call(what = rbind,
                         args = replicate(n = 1000,
                                          expr = sim_rep(sigmaU = 1),
                                          simplify = FALSE))
res_sigmaU_50 = do.call(what = rbind,
                       args = replicate(n = 1000,
                                        expr = sim_rep(sigmaU = 0.5),
                                        simplify = FALSE))
res_sigmaU_25 = do.call(what = rbind,
                       args = replicate(n = 1000,
                                        expr = sim_rep(sigmaU = 0.25),
                                        simplify = FALSE))
res_sigmaU_10 = do.call(what = rbind,
                        args = replicate(n = 1000,
                                         expr = sim_rep(sigmaU = 0.1),
                                         simplify = FALSE))
## Combine results and save
all_res_negC = res_sigmaU_100 |>
  mutate(var_u = 1) |>
  bind_rows(
    res_sigmaU_50 |>
      mutate(var_u = 0.5 ^ 2)
  ) |>
  bind_rows(
    res_sigmaU_25 |>
      mutate(var_u = 0.25 ^ 2)
  ) |>
  bind_rows(
    res_sigmaU_10 |>
      mutate(var_u = 0.1 ^ 2)
  )
all_res_negC |>
  write.csv("~/Documents/me_ci/sim_res_vary_sigmaU_negC.csv",
            row.names = FALSE)
## Calculate average Var(W), Cov(R,W) per Var(U)
all_res_negC |>
  group_by(var_u) |>
  summarize(avg_var_w = mean(var_w),
            avg_cov_rw = mean(cov_rw))
# Plot Var(W) vs Var(U) and Cov(R,W) vs. Var(U) --------------------------------
plot_var_w = all_res_negC |>
  ggplot(aes(x = factor(var_u), y = var_w)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_boxplot(fill = "#39558CFF") +
  xlab(label = "Error Variance in Exposure: Var(U)") +
  ylab(label = "Error Variance in Ranks: Var(W)") +
  ylim(c(-0.04, 0.07)) +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_text(face = "bold"))
plot_cov_rw = all_res_negC |>
  ggplot(aes(x = factor(var_u), y = cov_rw)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_boxplot(fill = "#AB337CFF") +
  xlab(label = "Error Variance in Exposure: Var(U)") +
  ylab(label = "Covariance Between Ranks and Errors: Cov(R,W)") +
  ylim(c(-0.04, 0.07)) +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_text(face = "bold"))

# Combine and save them
ggpubr::ggarrange(plot_var_w, plot_cov_rw, labels = "AUTO")
ggsave("~/Documents/me_ci/rank_error_var_covar_boxplot_negC.png",
       device = "png", width = 10, height = 7, units = "in")

# Sett 2: C = 0.5 -------------------------------------------------------------
set.seed(11182025) ## For reproducibility
res_sigmaU_100 = do.call(what = rbind,
                         args = replicate(n = 1000,
                                          expr = sim_rep(sigmaU = 1,
                                                         beta0 = -0.5,
                                                         beta1 = 3),
                                          simplify = FALSE))
res_sigmaU_50 = do.call(what = rbind,
                        args = replicate(n = 1000,
                                         expr = sim_rep(sigmaU = 0.5,
                                                        beta0 = -0.5,
                                                        beta1 = 3),
                                         simplify = FALSE))
res_sigmaU_25 = do.call(what = rbind,
                        args = replicate(n = 1000,
                                         expr = sim_rep(sigmaU = 0.25,
                                                        beta0 = -0.5,
                                                        beta1 = 3),
                                         simplify = FALSE))
res_sigmaU_10 = do.call(what = rbind,
                        args = replicate(n = 1000,
                                         expr = sim_rep(sigmaU = 0.1,
                                                        beta0 = -0.5,
                                                        beta1 = 3),
                                         simplify = FALSE))

# Combine results
all_res_posC = res_sigmaU_100 |>
  mutate(var_u = 1) |>
  bind_rows(
    res_sigmaU_50 |>
      mutate(var_u = 0.5 ^ 2)
  ) |>
  bind_rows(
    res_sigmaU_25 |>
      mutate(var_u = 0.25 ^ 2)
  ) |>
  bind_rows(
    res_sigmaU_10 |>
      mutate(var_u = 0.1 ^ 2)
  )
all_res_posC |>
  write.csv("~/Documents/me_ci/sim_res_vary_sigmaU_posC.csv",
            row.names = FALSE)

# Sett 3: C = 0 ----------------------------------------------------------------
set.seed(11182025) ## For reproducibility
res_sigmaU_100 = do.call(what = rbind,
                         args = replicate(n = 1000,
                                          expr = sim_rep(sigmaU = 1,
                                                         beta0 = 3,
                                                         beta1 = 0),
                                          simplify = FALSE))
res_sigmaU_50 = do.call(what = rbind,
                        args = replicate(n = 1000,
                                         expr = sim_rep(sigmaU = 0.5,
                                                        beta0 = 3,
                                                        beta1 = 0),
                                         simplify = FALSE))
res_sigmaU_25 = do.call(what = rbind,
                        args = replicate(n = 1000,
                                         expr = sim_rep(sigmaU = 0.25,
                                                        beta0 = 3,
                                                        beta1 = 0),
                                         simplify = FALSE))
res_sigmaU_10 = do.call(what = rbind,
                        args = replicate(n = 1000,
                                         expr = sim_rep(sigmaU = 0.1,
                                                        beta0 = 3,
                                                        beta1 = 0),
                                         simplify = FALSE))

# Combine results
all_res_nullC = res_sigmaU_100 |>
  mutate(var_u = 1) |>
  bind_rows(
    res_sigmaU_50 |>
      mutate(var_u = 0.5 ^ 2)
  ) |>
  bind_rows(
    res_sigmaU_25 |>
      mutate(var_u = 0.25 ^ 2)
  ) |>
  bind_rows(
    res_sigmaU_10 |>
      mutate(var_u = 0.1 ^ 2)
  )
all_res_nullC |>
  write.csv("~/Documents/me_ci/sim_res_vary_sigmaU_nullC.csv",
            row.names = FALSE)
## Average Var(W), Cov(R,W) per Var(U)
all_res_nullC |>
  group_by(var_u) |>
  summarize(avg_var_w = mean(var_w),
            avg_cov_rw = mean(cov_rw))
# Table of all settings --------------------------------------------------------
all_res = all_res_posC |>
  mutate(C = 0.5) |>
  bind_rows(
    all_res_negC |>
      mutate(C = -0.5)
  ) |>
  bind_rows(
    all_res_nullC |>
      mutate(C = 0)
  )
all_res = all_res |>
  # group_by(C, var_u) |>
  # summarize(avg_c_hat = mean(c_hat),
  #           ese_c_hat = sd(c_hat),
  #           ase_c_hat = mean(se_c_hat),
  #           avg_cstar_hat = mean(cstar_hat),
  #           ese_cstar_hat = sd(cstar_hat),
  #           ase_cstar_hat = mean(se_cstar_hat)) |>
  mutate(var_u = factor(x = var_u,
                        levels = c(0.01, 0.0625, 0.25, 1),
                        labels = paste0("Var(U)=",
                                        c(0.01, 0.0625, 0.25, 1))),
         Cf = factor(x = C,
                    levels = c(-0.5, 0, 0.5),
                    labels = paste0("True C=", c(-0.5, 0, 0.5))))

all_res |>
  ggplot(aes(x = cstar_hat, fill = factor(var_u))) +
  geom_vline(aes(xintercept = C), linetype = "dashed") +
  geom_density(alpha = 0.7) +
  facet_wrap(~Cf, scales = "free") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white", face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_fill_manual(name = "Exposure Error Variance",
                    values = c("#39558CFF", "#AB337CFF", "#287D8EFF",  "#E85362FF",  "#481568FF")) +
  xlab("Naive Concentration Index Estimate") +
  ylab("Density of Replicates")
ggsave("~/Documents/me_ci/density_vary_sigmaU.png",
       device = "png", width = 10, height = 7, units = "in")

## Average Var(W), Cov(R,W) per Var(U)
all_res |>
  ggplot(aes(x = var_w, y = -cov_rw, color = factor(var_u))) +
  geom_point() +
  facet_wrap(~C) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Exposure Error Variance",
                     values = c("#39558CFF", "#AB337CFF", "#287D8EFF",  "#E85362FF",  "#481568FF")) +
  facet_wrap(~Cf, scales = "free") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white", face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  xlab("Error Variance in Ranks: Var(W)") +
  ylab("Negative Covariance Between Ranks and Errors: -Cov(R,W)")
ggsave("~/Documents/me_ci/show_sett2b_vary_sigmaU.png",
       device = "png", width = 10, height = 7, units = "in")

## Average Var(W), Cov(R,W) per Var(U)
all_res |>
  group_by(C, var_u) |>
  summarize(avg_var_r = mean(var_r),
            avg_var_w = mean(var_w),
            avg_cov_rw = mean(cov_rw),
            avg_lambda = mean(cstar_hat / c_hat),
            avg_c_hat = mean(c_hat),
            # ese_c_hat = sd(c_hat),
            # ase_c_hat = mean(se_c_hat),
            avg_cstar_hat = mean(cstar_hat),
            # ese_cstar_hat = sd(cstar_hat),
            # ase_cstar_hat = mean(se_cstar_hat)
  ) |>
  # group_by(C, var_u) |>
  # summarize() |>
  mutate(nu = avg_var_r + avg_cov_rw,
         delta = avg_var_r + avg_var_w + 2 * avg_cov_rw,
         avg_exp_lambda = nu / delta) |>
  ggplot(aes(x = avg_exp_lambda, y = avg_lambda)) +
  geom_point() +
  facet_wrap(~C) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
