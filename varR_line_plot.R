library(dplyr) ## for data wrangling
library(ggplot2) ## for plots

set.seed(1)
n = 10000

rank_var = function(x) {
  sapply(X = x,
         FUN = function(x) var((rank(1:x) - 1)/x + 1/(2 * x)))
}

rank_plot_dat = data.frame(n = seq(25, 10000, by = 25))
rank_plot_dat$sigma2R = rank_var(x = rank_plot_dat$n)

rank_plot_dat |>
  ggplot(aes(x = n, y = sigma2R)) +
  geom_line(color = "#E85362FF",
            size = 1.2) +
  theme_minimal() +
  labs(x = "N",
       y = "Var(R)") +
  theme(axis.title = element_text(face = "bold")) +
  scale_x_continuous(labels = scales::label_comma())
ggsave(filename = "~/Documents/me_ci/varR_line_plot.png",
       device = "png",
       width = 7,
       height = 5,
       units = "in")
ggsave(filename = "~/Documents/me_ci/varR_line_plot.svg",
       device = "svg",
       width = 7,
       height = 5,
       units = "in")

expand_rank_plot_dat = rank_plot_dat |>
  mutate(sigma2W = 1) |>
  bind_rows(
    rank_plot_dat |>
      mutate(sigma2W = 0.5)
  ) |>
  bind_rows(
    rank_plot_dat |>
      mutate(sigma2W = 0.1)
  ) |>
  mutate(ub_lambda = (sigma2R + sqrt(sigma2R * sigma2W)) /
           (sigma2R + sigma2W + 2 * sqrt(sigma2R * sigma2W)),
         alt_ub_lambda = sqrt(sigma2R) * (sqrt(sigma2R) + sqrt(sigma2W)) /
           (sigma2R + sigma2W + 2 * sqrt(sigma2R * sigma2W)))
with(expand_rank_plot_dat, all.equal(ub_lambda, alt_ub_lambda))

expand_rank_plot_dat |>
  ggplot(aes(x = sigma2R, y = ub_lambda)) +
  geom_point(color = "#482677FF",
             size = 1.2) +
  theme_minimal() +
  labs(x = "Var(R)",
       y = "lambda") +
  theme(axis.title = element_text(face = "bold")) +
  scale_x_continuous(labels = scales::label_comma())
ggsave(filename = "~/Documents/me_ci/varR_line_plot.png",
       device = "png",
       width = 7,
       height = 5,
       units = "in")
ggsave(filename = "~/Documents/me_ci/varR_line_plot.svg",
       device = "svg",
       width = 7,
       height = 5,
       units = "in")
