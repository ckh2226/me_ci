n = 10000

# 1) Rural/not rural indicator
Z = rbinom(n = n, size = 1, prob = 0.2)

# 2) True proximity | Rural / not rural
X = rnorm(n = n, mean = 1 + 4 * Z, sd = 1)

# Calculate rank
R = (rank(X) - 1)/n + 1/(2 * n)
sigma2R = var(R)

# Generate errors and error-prone X
U = rnorm(n = n, mean = 0, sd = 1)
Xstar = X + U

# Calculate error-prone rank
Rstar = (rank(Xstar) - 1)/n + 1/(2 * n)
sigma2Rstar = var(Rstar)
sigma2W = var(Rstar - R)

try_sigmaRW = seq(-0.15, 0.15, by = 0.0001)
try_sigma2W = seq(0.01, 0.15, by = 0.0001)
try_vals = expand.grid(try_sigmaRW = try_sigmaRW,
                       try_sigma2W = try_sigma2W)
try_lambda = with(try_vals,
                  (sigma2R + try_sigmaRW) /
                    (sigma2R + try_sigma2W + 2 * try_sigmaRW))
try_vals$sigma2R = sigma2R
try_vals$lambda = try_lambda
cauchy_schwartz = with(try_vals,
                       try_sigmaRW ^ 2 <= (sigma2R * try_sigma2W))
try_vals = try_vals[cauchy_schwartz, ]

library(ggplot2)
summary(try_vals$lambda[try_vals$lambda < 0]) ## -45, -1, -0.4, -0.14
summary(try_vals$lambda[try_vals$lambda >= 0 & try_vals$lambda <= 1]) ## 0, 0.4, 0.47, 0.6
summary(try_vals$lambda[try_vals$lambda > 1]) ## 1, 1.1, 1.3, 1.7, 45
contour_plot = try_vals |>
  ggplot(aes(x = try_sigma2W, y = try_sigmaRW, z = lambda)) +
  # geom_contour_filled(breaks = c(-40, -5, -1,
  #                                seq(0, 1, by = 0.25),
  #                                5, 10, 40)) +
  geom_contour_filled(breaks = c(-45, -1, -0.4, -0.14,
                                 0, 0.4, 0.47, 0.6,
                                 1, 1.1, 1.3, 1.7, 45)) +
  theme_minimal() +
  labs(x = "Var(W)",
       y = "Cov(R,W)",
       fill = latex2exp::TeX("Bias\nFactor: $\\lambda$",
                             bold = TRUE)) +
  theme(axis.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "azure"),
        panel.grid = element_blank())
contour_plot
ggsave(filename = "~/Documents/me_ci/lambda_contour_plot_plain.png",
       device = "png",
       width = 7,
       height = 5,
       units = "in")
ggsave(filename = "~/Documents/me_ci/lambda_contour_plot_plain.svg",
       device = "svg",
       width = 7,
       height = 5,
       units = "in")

contour_plot +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "azure") +
  annotate(geom = "text",
           label = "Cov(R,W)=0",
           color = "azure",
           x = 0.026,
           y = 0.005) +
  geom_abline(slope = -1,
              intercept = 0,
              linetype = "dashed",
              color = "azure") +
  annotate(geom = "text",
           label = "-Cov(R,W) <= Var(W)",
           color = "azure",
           x = 0.026,
           y = 0.005,
           angle = -40)
