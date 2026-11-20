# ==============================================================================
# 14_age_model_checks.R
# Two checks of the age model (ESM Methods S3):
#   1. Elapsed time: for colonies that grew between the surveys, the modelled
#      time between their 2015 and 2022 sizes should match the true 7 years.
#   2. Detection floor: how much do modelled ages change if A0 is doubled?
#      Compared for colonies above and below A_thresh (2022 size), with the
#      same 5,000 parameter draws, so only A0 changes.
#
# Reads:  results/colony_ages.csv, models/age_simulation_inputs.rds (script 08),
#         models/athresh_sizes.rds (script 07)
# Writes: results/age_elapsed_time_check.csv   check 1, by species
#         results/a0_sensitivity_ages.csv      check 2, per colony
#         figures/a0_sensitivity.png           check 2, baseline vs doubled A0
# ==============================================================================

source("setup.R")
source("age_model.R")     # species_params(), age_at_size()
suppressPackageStartupMessages(library(tidyverse))

A0_SHIFT_PCT <- 1.00      # +100%: A0 doubled

inputs      <- readRDS(in_data("models", "age_simulation_inputs.rds"))
colonies    <- inputs$stand_age_audited %>% mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id))
A0_shifted  <- inputs$A0_species * (1 + A0_SHIFT_PCT)
athresh     <- readRDS(in_data("models", "athresh_sizes.rds"))
colony_ages <- read_csv(in_data("results", "colony_ages.csv"), show_col_types = FALSE) %>%
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id))

cat("A0 (cm2), baseline and doubled:\n")
print(round(rbind(baseline = inputs$A0_species, doubled = A0_shifted), 3))

# ------------------------------------------------------------------------------
# Median age over the draws of every colony, with A0 doubled (about 4 minutes)
# ------------------------------------------------------------------------------
draw_params <- lapply(setNames(SPECIES, SPECIES), function(sp) species_params(inputs$param_draws, sp))
median_age <- function(size, sp) {
  if (is.na(size) || size <= 0 || !sp %in% SPECIES) return(NA_real_)
  P <- draw_params[[sp]]
  median(sapply(seq_len(nrow(P)), function(k) age_at_size(size, P[k, ], A0_shifted[sp])), na.rm = TRUE)
}
shifted <- colonies %>%
  mutate(Age_2015_Median_shifted = map2_dbl(A1, Species, median_age),
         Age_2022_Median_shifted = map2_dbl(A2, Species, median_age))

comparison <- colony_ages %>%
  select(TagLab.Genet.Id, Species, A1, A2, tracked_both,
         Age_2015_Median, Age_2022_Median, Delta_Age_Median) %>%
  left_join(select(shifted, TagLab.Genet.Id, Age_2015_Median_shifted, Age_2022_Median_shifted),
            by = "TagLab.Genet.Id") %>%
  mutate(Delta_Age_Median_shifted = ifelse(
    tracked_both & !is.na(Age_2015_Median_shifted) & !is.na(Age_2022_Median_shifted) & A2 > A1,
    Age_2022_Median_shifted - Age_2015_Median_shifted, NA_real_))

# ------------------------------------------------------------------------------
# Check 1: modelled elapsed time vs the true 7 years
# ------------------------------------------------------------------------------
elapsed_check <- comparison %>%
  filter(tracked_both, !is.na(Delta_Age_Median)) %>%
  group_by(Species) %>%
  summarise(n                         = n(),
            True_Interval_Yrs         = SURVEY_INTERVAL,
            Median_Delta_Age_Baseline = round(median(Delta_Age_Median, na.rm = TRUE), 1),
            Median_Delta_Age_Shifted  = round(median(Delta_Age_Median_shifted, na.rm = TRUE), 1),
            .groups = "drop") %>%
  mutate(Baseline_minus_True = round(Median_Delta_Age_Baseline - True_Interval_Yrs, 1),
         Shifted_minus_True  = round(Median_Delta_Age_Shifted  - True_Interval_Yrs, 1))

cat("\nCheck 1: modelled time between 2015 and 2022 size vs the true", SURVEY_INTERVAL, "years:\n")
print(as.data.frame(elapsed_check), row.names = FALSE)
write_csv(elapsed_check, in_data("results", "age_elapsed_time_check.csv"))

# ------------------------------------------------------------------------------
# Check 2: change in 2022 age with A0 doubled, above and below A_thresh
# ------------------------------------------------------------------------------
shift_check_colonies <- comparison %>%
  filter(Species %in% SPECIES, !is.na(A2), A2 > 0,
         !is.na(Age_2022_Median), !is.na(Age_2022_Median_shifted)) %>%
  mutate(Athresh_cm2   = athresh[Species],
         Size_Class    = ifelse(A2 > Athresh_cm2, "Adult (A > Athresh)", "Sub-threshold (A < Athresh)"),
         Abs_Shift_Yrs = Age_2022_Median_shifted - Age_2022_Median,
         Pct_Shift     = ifelse(Age_2022_Median > 0, 100 * Abs_Shift_Yrs / Age_2022_Median, NA_real_))

shift_summary <- shift_check_colonies %>%
  group_by(Size_Class) %>%
  summarise(n                    = n(),
            Median_Pct_Shift     = round(median(Pct_Shift, na.rm = TRUE), 2),
            Mean_Pct_Shift       = round(mean(Pct_Shift, na.rm = TRUE), 2),
            Median_Abs_Shift_Yrs = round(median(Abs_Shift_Yrs, na.rm = TRUE), 2),
            .groups = "drop")
shift_by_species <- shift_check_colonies %>%
  group_by(Species, Size_Class) %>%
  summarise(n                    = n(),
            Median_Pct_Shift     = round(median(Pct_Shift, na.rm = TRUE), 2),
            Median_Abs_Shift_Yrs = round(median(Abs_Shift_Yrs, na.rm = TRUE), 2),
            .groups = "drop")

cat("\nCheck 2: change in 2022 age with A0 doubled, all species:\n")
print(as.data.frame(shift_summary), row.names = FALSE)
cat("\nBy species:\n")
print(as.data.frame(shift_by_species), row.names = FALSE)

write_csv(shift_check_colonies %>%
            select(TagLab.Genet.Id, Species, A2, Athresh_cm2, Size_Class,
                   Age_2022_Median, Age_2022_Median_shifted, Abs_Shift_Yrs, Pct_Shift),
          in_data("results", "a0_sensitivity_ages.csv"))

p_shift <- ggplot(shift_check_colonies %>% mutate(Species = factor(Species, levels = SPECIES)),
                  aes(x = Age_2022_Median, y = Age_2022_Median_shifted, color = Species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.5, size = 1.4) +
  facet_wrap(~ Size_Class) +
  scale_color_manual(values = SPECIES_COLOURS, name = "Species") +
  labs(title = paste0("Colony age under baseline vs. A0 +", A0_SHIFT_PCT * 100, "% detection floor"),
       x = "Baseline A0 -- 2022 median age (yr)",
       y = "Shifted A0 -- 2022 median age (yr)") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom")
save_figure(p_shift, "a0_sensitivity.png", width = 178, height = 100)
