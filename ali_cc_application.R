# Load packages
library(dplyr)
library(tidyr)
library(latex2exp)
library(ggplot2)
library(rineq)
library(ggconc)

# Read in data
## Demographics + ALI components
all_ali_dat <- read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv", header=TRUE)
## Separate hospitalizations + ED visits outcomes
encounters_dat <- read.csv("~/Documents/Allostatic_load_audits/Raw_Data/hospital_ed_separate.csv")
# Merge demographics + ALI components + Hospitalizations/ED visits
all_ali_dat <- all_ali_dat |>
  left_join(y = encounters_dat) |>
  replace_na(replace = list(NUM_ADMIT = 0, NUM_ED = 0)) |> ## replace NA with 0
  mutate(ANY_ADMIT = as.numeric(NUM_ADMIT > 0), ## Define binary indicators
         ANY_ED = as.numeric(NUM_ED > 0), ## Define binary indicators
         ANY_ENCOUNTERS = ANY_ADMIT, ## Override ANY_ENCOUNTER with ANY_ADMIT
         SEX = factor(x = SEX, ## Convert SEX to factor variable
                      levels = c("Male","Female")))
# Create separate dataframe for EHR (Before Validation)
unvalid_dat <- all_ali_dat |>
  filter(DATA == "EHR (Before Validation)")

## More demographics
more_demo = read.csv("~/Documents/Allostatic_load_audits/summary_data.csv", header=TRUE) |>
  select(PAT_MRN_ID, RACE, ETHNICITY) |>
  mutate(RACE3 = if_else(RACE %in% c("American Indian or Alaska Native",
                                     "Asian Indian",
                                     "Other"), "Other Race", RACE))
unvalid_dat <- unvalid_dat |>
  left_join(more_demo)

everyone_cc = unvalid_dat |>
  make_conc_data(outcome_var = NUM_ADMIT,
                 rank_var = ALI,
                 rank_ascend = FALSE) |>
  make_conc_curve() +
  ggtitle(label = "All Patients")

sex_cc = unvalid_dat |>
  make_conc_data(outcome_var = NUM_ADMIT,
                 rank_var = ALI,
                 rank_ascend = FALSE,
                 group_var = SEX) |>
  make_conc_curve() +
  scale_color_manual(values = c("#39558CFF", "#AB337CFF", "#287D8EFF",  "#E85362FF",  "#481568FF"),
                     name = "Sex") +
  labs(x = "Cumulative Proportion of Patients (Ranked by Error-Prone ALI)",
       y = "Cumulative Proportion of Hospitalizations") +
  facet_wrap(~Group) +
  geom_text(data = data.frame(x_pos = 0,
                              y_pos = 1,
                              label = factor(x = c(-0.31, -0.08),
                                             labels = c(TeX("$\\widehat{C} = -0.31 (-0.44, -0.18)$"),
                                                        TeX("$\\widehat{C} = -0.08 (-0.24, -0.09)$"))),
                              Group = c("Female",
                                        "Male")),
            aes(x = x_pos, y = y_pos, label = label),
            hjust = 0, vjust = 1, parse = TRUE) +
  theme(strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"))
sex_cc
ci(ineqvar = -unvalid_dat$ALI[unvalid_dat$SEX == "Female"],
   outcome = unvalid_dat$NUM_ADMIT[unvalid_dat$SEX == "Female"]) |>
  summary()
## -0.3097126 (-0.4432892 -0.1761361)
ci(ineqvar = -unvalid_dat$ALI[unvalid_dat$SEX == "Male"],
   outcome = unvalid_dat$NUM_ADMIT[unvalid_dat$SEX == "Male"]) |>
  summary()
## -0.07934105 (-0.2488523 0.09017021 )

race_cc = unvalid_dat |>
  make_conc_data(outcome_var = NUM_ADMIT,
                 rank_var = ALI,
                 rank_ascend = FALSE,
                 group_var = RACE3) |>
  make_conc_curve() +
  scale_color_manual(values = rev(c("#39558CFF", "#AB337CFF", "#287D8EFF",  "#E85362FF",  "#481568FF")),
                     name = "Race",
                     guide = "none") +
  labs(x = "Cumulative Proportion of Patients (Ranked by Error-Prone ALI)",
       y = "Cumulative Proportion of Hospitalizations") +
  facet_wrap(~Group) +
  geom_text(data = data.frame(x_pos = 0,
                              y_pos = 1,
                              label = factor(x = c(-0.20, -0.22, -0.18),
                                             labels = c(TeX("$\\widehat{C} = -0.20 (-0.43, 0.03)$"),
                                                        TeX("$\\widehat{C} = -0.22 (-0.35, -0.09)$"),
                                                        TeX("$\\widehat{C} = -0.18 (-0.57, 0.21)$"))),
                              Group = c("Black or African American",
                                        "White or Caucasian",
                                        "Other Race")),
            aes(x = x_pos, y = y_pos, label = label),
            hjust = 0, vjust = 1, parse = TRUE) +
  theme(strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"))
race_cc
ci(ineqvar = -unvalid_dat$ALI[unvalid_dat$RACE3 == "White or Caucasian"],
   outcome = unvalid_dat$NUM_ADMIT[unvalid_dat$RACE3 == "White or Caucasian"]) |>
  summary()
## -0.2200671 (-0.3469516 -0.0931826)
ci(ineqvar = -unvalid_dat$ALI[unvalid_dat$RACE3 == "Black or African American"],
   outcome = unvalid_dat$NUM_ADMIT[unvalid_dat$RACE3 == "Black or African American"]) |>
  summary()
## -0.2001342 (-0.4321174 0.03184889)
ci(ineqvar = -unvalid_dat$ALI[unvalid_dat$RACE3 == "Other Race"],
   outcome = unvalid_dat$NUM_ADMIT[unvalid_dat$RACE3 == "Other Race"]) |>
  summary()
## -0.1802373 (-0.5699532 0.2094786 )
