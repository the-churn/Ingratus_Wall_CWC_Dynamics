# ==============================================================================
# 18_manuscript_numbers.R
# Numbers quoted in the Results and Discussion that are not printed by the
# other scripts. Prints to its log only.
#   1. age at A_thresh (median and 95% interval over the 5,000 draws)
#   2. ages of the colonies present in 2022
#   3. colonies tracked in both surveys
#   4. share of colonies older than the age at A_thresh
#   5. pairwise Wald tests of the growth-model parameters (all species pairs)
#   6. persistent colonies that did not grow between the surveys
#   7. modelled ages of the colonies that recruited between the surveys
#   8. could differences in survival explain the colonisation order?
#
# Reads:  results/colony_ages.csv, models/age_simulation_inputs.rds (script 08),
#         models/athresh_sizes.rds (script 07), models/growth_model.rds (script 05)
# ==============================================================================

source("setup.R")
source("age_model.R")     # species_params(), age_at_size()
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
})
options(width = 160)

final_df <- read_csv(in_data("results", "colony_ages.csv"), show_col_types = FALSE)
inputs   <- readRDS(in_data("models", "age_simulation_inputs.rds"))
athresh  <- readRDS(in_data("models", "athresh_sizes.rds"))
m_gnls   <- readRDS(in_data("models", "growth_model.rds"))

section <- function(title) cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")

# ------------------------------------------------------------------------------
section("1. Age at A_thresh (median and 95% interval over the parameter draws)")
# ------------------------------------------------------------------------------
athresh_summary <- map_dfr(names(athresh), function(sp) {
  P    <- species_params(inputs$param_draws, sp)
  ages <- sapply(seq_len(nrow(P)), function(k) age_at_size(athresh[[sp]], P[k, ], inputs$A0_species[sp]))
  data.frame(Species     = sp,
             Athresh_cm2 = athresh[[sp]],
             Median_Age  = round(median(ages, na.rm = TRUE), 1),
             Lower_95    = round(quantile(ages, 0.025, na.rm = TRUE), 2),
             Upper_95    = round(quantile(ages, 0.975, na.rm = TRUE), 2),
             row.names   = NULL)
})
print(athresh_summary, row.names = FALSE)

# ------------------------------------------------------------------------------
section("2. Modelled ages of the colonies present in 2022")
# ------------------------------------------------------------------------------
pop_2022 <- final_df %>%
  filter(!is.na(Age_2022_Median), Species %in% SPECIES) %>%
  group_by(Species) %>%
  summarise(Total_Colonies_2022 = n(),
            Median_Age_2022     = round(median(Age_2022_Median), 1),
            IQR_Age_2022        = round(IQR(Age_2022_Median), 1),
            Min_Age             = round(min(Age_2022_Median), 1),
            Max_Age             = round(max(Age_2022_Median), 1),
            .groups = "drop")
print(as.data.frame(pop_2022), row.names = FALSE)
cat("Colonies with a 2022 age:", sum(pop_2022$Total_Colonies_2022), "\n")

# ------------------------------------------------------------------------------
section("3. Colonies tracked in both surveys")
# ------------------------------------------------------------------------------
print(as.data.frame(final_df %>% filter(tracked_both, Species %in% SPECIES) %>%
                      count(Species, name = "Tracked_Growth_Pairs")), row.names = FALSE)

# ------------------------------------------------------------------------------
section("4. Colonies older in 2022 than the age at A_thresh")
# ------------------------------------------------------------------------------
maturation <- final_df %>%
  filter(!is.na(Age_2022_Median), Species %in% SPECIES) %>%
  left_join(select(athresh_summary, Species, Athresh_Age = Median_Age), by = "Species") %>%
  group_by(Species) %>%
  summarise(Total_n             = n(),
            Athresh_Age_Yrs     = first(Athresh_Age),
            n_Above_Threshold   = sum(Age_2022_Median >= Athresh_Age),
            Pct_Above_Threshold = round(100 * mean(Age_2022_Median >= Athresh_Age), 1),
            .groups = "drop")
print(as.data.frame(maturation), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Pairwise Wald tests: the model summary compares each species with
# M. oculata only; here every pair is tested from the coefficients and their
# covariance (difference, SE, t on the residual df)
# ------------------------------------------------------------------------------
b_hat  <- coef(m_gnls)
V_hat  <- vcov(m_gnls)
df_res <- m_gnls$dims$N - m_gnls$dims$p
section(paste0("5. Pairwise Wald tests of the growth-model parameters (", df_res, " df)"))

contrast_vec <- function(par, sp_a, sp_b) {   # species a minus species b
  v <- setNames(numeric(length(b_hat)), names(b_hat))
  if (sp_a != SPECIES[1]) v[paste0(par, ".fSpecies", sp_a)] <- v[paste0(par, ".fSpecies", sp_a)] + 1
  if (sp_b != SPECIES[1]) v[paste0(par, ".fSpecies", sp_b)] <- v[paste0(par, ".fSpecies", sp_b)] - 1
  v
}
species_pairs <- list(c("Desmophyllum pertusum", "Madrepora oculata"),
                      c("Primnoa msp.",          "Madrepora oculata"),
                      c("Desmophyllum pertusum", "Primnoa msp."))
par_labels <- c(Asym = "RGR_inf", R0 = "RGR0", lrc = "lrc")
wald_tab <- map_dfr(names(par_labels), function(par) map_dfr(species_pairs, function(pr) {
  l   <- contrast_vec(par, pr[1], pr[2])
  est <- sum(l * b_hat)
  se  <- sqrt(drop(t(l) %*% V_hat %*% l))
  data.frame(Parameter = par_labels[[par]], Contrast = paste(pr[1], "-", pr[2]),
             Estimate = signif(est, 4), SE = signif(se, 3), t = round(est / se, 2),
             p = signif(2 * pt(-abs(est / se), df_res), 3))
}))
print(wald_tab, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. "Did not grow" = 2022 area no larger than in 2015; "shrank detectably" =
# area loss larger than the colony's own MDC95. Also split by 2015 size
# relative to A_thresh.
# ------------------------------------------------------------------------------
section("6. Persistent colonies that did not grow between the surveys")
persistent <- final_df %>%
  filter(Species %in% SPECIES, state_transition == "persistence") %>%
  mutate(no_growth         = change_area <= 0,
         shrank_detectably = as.logical(detectable_change) %in% TRUE & change_area < 0,
         above_thresh      = A1 > unname(athresh[Species]))
print(as.data.frame(persistent %>% group_by(Species) %>%
  summarise(n_persistent = n(), n_no_growth = sum(no_growth), Pct_no_growth = round(100 * mean(no_growth), 1),
            n_shrank_detectably = sum(shrank_detectably), .groups = "drop")), row.names = FALSE)
cat("\nBy 2015 size relative to A_thresh:\n")
print(as.data.frame(persistent %>%
  group_by(Species, Size_2015 = ifelse(above_thresh, "Larger than A_thresh", "Smaller than A_thresh")) %>%
  summarise(n = n(), n_no_growth = sum(no_growth), Pct_no_growth = round(100 * mean(no_growth), 1),
            .groups = "drop")), row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Colonies absent in 2015 and present in 2022 reached a detectable size
# within the 7 years, so a modelled age above 7 years places them too early
# ------------------------------------------------------------------------------
section("7. Modelled ages of the colonies that recruited between the surveys")
recruit_tab <- final_df %>%
  filter(Species %in% SPECIES, state_transition == "recruitment", !is.na(Age_2022_Median)) %>%
  group_by(Species) %>%
  summarise(n_recruits            = n(),
            Median_age            = round(median(Age_2022_Median), 1),
            n_age_over_interval   = sum(Age_2022_Median > SURVEY_INTERVAL),
            Pct_age_over_interval = round(100 * mean(Age_2022_Median > SURVEY_INTERVAL), 1),
            Max_age               = round(max(Age_2022_Median), 1),
            .groups = "drop") %>%
  left_join(final_df %>% filter(Species %in% SPECIES, !is.na(Age_2022_Median)) %>%
              group_by(Species) %>%
              summarise(n_all_age_le_interval = sum(Age_2022_Median <= SURVEY_INTERVAL), .groups = "drop"),
            by = "Species")
print(as.data.frame(recruit_tab), row.names = FALSE)
cat("(n_all_age_le_interval: all colonies the model places in the last", SURVEY_INTERVAL,
    "years, recruits or not)\n")

# ------------------------------------------------------------------------------
# 8. With per-capita mortality h (per yr, from the 2015-2022 mortality) and
# steady recruitment for T years, the share of living colonies older than A is
# (exp(-h A) - exp(-h T)) / (1 - exp(-h T)). Expected counts use each species'
# 2022 sample size; observed counts use the modelled ages.
# ------------------------------------------------------------------------------
SURVIVAL_AGE <- 40
STEADY_YEARS <- c(60, 100)
section(paste("8. Survival check: colonies older than", SURVIVAL_AGE,
              "yr, observed vs expected under steady recruitment"))
survival_tab <- final_df %>%
  filter(Species %in% SPECIES) %>%
  group_by(Species) %>%
  summarise(n_2015           = sum(state_transition %in% c("persistence", "mortality")),
            n_died           = sum(state_transition == "mortality"),
            n_2022           = sum(!is.na(Age_2022_Median)),
            n_older_observed = sum(Age_2022_Median > SURVIVAL_AGE, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Mortality_7yr_pct = round(100 * n_died / n_2015, 1),
         Hazard_per_yr     = -log(1 - n_died / n_2015) / SURVEY_INTERVAL)
for (T_yr in STEADY_YEARS) {
  survival_tab[[paste0("n_older_expected_", T_yr, "yr")]] <- with(survival_tab,
    round(n_2022 * (exp(-Hazard_per_yr * SURVIVAL_AGE) - exp(-Hazard_per_yr * T_yr)) /
            (1 - exp(-Hazard_per_yr * T_yr)), 1))
}
print(as.data.frame(survival_tab %>% mutate(Hazard_per_yr = signif(Hazard_per_yr, 3))), row.names = FALSE)
