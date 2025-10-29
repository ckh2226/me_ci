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
