source("~/Documents/GitHub/me_ci/sim_data.R")
source("~/Documents/GitHub/me_ci/sim_val_mb.R")

set.seed(214)

val_srs <- do.call(rbind, replicate(10000, sim_val_mb(0.5, 1000, approx_ci = -0.5, pv = 0.1, design = "SRS"), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = -0.5, design = "SRS")

val_cc <- do.call(rbind, replicate(10000, sim_val_mb(0.5, 1000, approx_ci = -0.5, pv = 0.1, design = "CC"), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = -0.5, design = "CC")

val_bcc <- do.call(rbind, replicate(10000, sim_val_mb(0.5, 1000, approx_ci = -0.5, pv = 0.1, design = "BCC"), simplify = FALSE)) |>
  data.frame() |>
  mutate(approx_ci = -0.5, design = "BCC")

all_df <- val_srs |>
  rbind(val_cc) |>
  rbind(val_bcc)

all_df |>
  select(starts_with(c("cov")), design) |> ## only pull Cov(R,w) and study design
  pivot_longer(cols = -design, 
               names_to = "Method", 
               values_to = "Estimate") |> ## pivot from wide --> long
  mutate(validated = grepl(pattern = "val", x = Method),
              Method = sub(pattern = "val", replacement = "", x = Method)) |>
  ggplot(aes(x = validated, y = Estimate, fill = validated)) +
  geom_boxplot() +
  facet_wrap(~ Method + design, scales = "free") +
  theme_bw()
