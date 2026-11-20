# ==============================================================================
# 08_colony_ages.R
# Colony ages from colony size (ESM Methods S3). The time a colony needs to
# grow from the detection floor A0 to its observed area A is
#
#   t(A) = integral from A0 to A of  da / (a * RGR(a))
#
# with RGR(a) the fitted growth curve of its species. Uncertainty: the
# integral is repeated for 5,000 parameter sets drawn from the model's
# estimates and their covariance (multivariate normal); ages are the median
# and 2.5-97.5% quantiles over the draws. A0 is the 3rd percentile of all
# observed areas of the species (both surveys).
#
# Colonies seen only in 2022 whose modelled age was 10 years or more were
# checked by eye in TagLab (TagLab_QuickJump_Audit.xlsx): "growth" = already
# present in 2015 (its 2015 area is taken from the audit), "occluded" = hidden
# in 2015, "artifact" = not a colony (dropped).
#
# Reads:  processed/colonies_with_mdc.xlsx, models/growth_model.rds (script 05),
#         TagLab_QuickJump_Audit.xlsx
# Writes: results/colony_ages.csv             one row per colony, with ages
#         models/age_simulation_inputs.rds    data, parameter draws and A0 (later scripts)
# ==============================================================================

source("setup.R")
source("age_model.R")     # species_params(), age_at_size()
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(nlme)    # vcov() for the gnls model
})

N_SIMS      <- 5000   # parameter draws
P_FLOOR     <- 0.03   # percentile of observed areas that sets A0
SURVEY_YEAR <- 2022

set.seed(42)

# ------------------------------------------------------------------------------
# Data and detection floor
# ------------------------------------------------------------------------------
flag_tracked <- function(df) {
  df %>%
    mutate(tracked_both = !is.na(A1) & !is.na(A2) & A1 > 0 & A2 > 0,
           delta_area   = ifelse(tracked_both, A2 - A1, NA_real_),
           is_growing   = ifelse(tracked_both, delta_area > 0, NA))
}

colonies <- openxlsx::read.xlsx(in_data("processed", "colonies_with_mdc.xlsx")) %>%
  mutate(Species = standard_species(trimws(Species))) %>%
  flag_tracked()

# A0: 3rd percentile of all areas of the species, 2015 and 2022 together
A0_species <- data.frame(Species = c(colonies$Species, colonies$Species),
                         Area    = c(colonies$A1, colonies$A2)) %>%
  filter(Species %in% SPECIES, !is.na(Area), Area > 0) %>%
  group_by(Species) %>%
  summarise(min_size = quantile(Area, probs = P_FLOOR), .groups = "drop") %>%
  with(setNames(min_size, Species))
cat("Detection floor A0 (cm2, 3rd percentile of observed areas):\n")
print(round(A0_species, 3))

# ------------------------------------------------------------------------------
# Parameter draws
# ------------------------------------------------------------------------------
m_growth_gnls <- readRDS(in_data("models", "growth_model.rds"))
param_draws <- MASS::mvrnorm(n = N_SIMS, mu = coef(m_growth_gnls), Sigma = vcov(m_growth_gnls))

# Median and 95% interval of the ages over all draws, for 2015 and 2022 size
draw_params <- lapply(setNames(SPECIES, SPECIES), function(sp) species_params(param_draws, sp))

simulate_ages <- function(df) {
  out <- matrix(NA_real_, nrow(df), 9)
  for (i in seq_len(nrow(df))) {
    sp <- df$Species[i]
    if (!sp %in% SPECIES) next
    ages <- function(size) {
      if (is.na(size) || size <= 0) return(NULL)
      sapply(seq_len(N_SIMS), function(k) age_at_size(size, draw_params[[sp]][k, ], A0_species[sp]))
    }
    summary_3 <- function(x) c(median(x, na.rm = TRUE),
                               quantile(x, probs = 0.025, na.rm = TRUE),
                               quantile(x, probs = 0.975, na.rm = TRUE))
    sim_2015 <- ages(df$A1[i])
    sim_2022 <- ages(df$A2[i])
    if (!is.null(sim_2015)) out[i, 1:3] <- summary_3(sim_2015)
    if (!is.null(sim_2022)) out[i, 4:6] <- summary_3(sim_2022)
    # Modelled time between the two sizes, for colonies that grew
    if (!is.null(sim_2015) && !is.null(sim_2022) && df$A2[i] > df$A1[i]) {
      out[i, 7:9] <- summary_3(sim_2022 - sim_2015)
    }
  }
  df %>%
    mutate(
      Age_2015_Median   = out[, 1], Age_2015_Lower95  = out[, 2], Age_2015_Upper95  = out[, 3],
      Age_2022_Median   = out[, 4], Age_2022_Lower95  = out[, 5], Age_2022_Upper95  = out[, 6],
      Delta_Age_Median  = out[, 7], Delta_Age_Lower95 = out[, 8], Delta_Age_Upper95 = out[, 9],
      Colonization_Year_Est = ifelse(!is.na(Age_2022_Median), SURVEY_YEAR - Age_2022_Median, NA_real_)
    )
}

# ------------------------------------------------------------------------------
# Apply the recruit audit
# ------------------------------------------------------------------------------
audit <- read_excel(in_data("TagLab_QuickJump_Audit.xlsx")) %>%
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id), A1_2015 = as.numeric(A1_2015)) %>%
  select(TagLab.Genet.Id, Audit_Status, Match_2015, A1_2015)

colonies_audited <- colonies %>%
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id)) %>%
  left_join(audit, by = "TagLab.Genet.Id")
stopifnot(nrow(colonies_audited) == nrow(colonies))   # IDs must be unique

colonies_audited <- colonies_audited %>%
  mutate(
    A1 = case_when(!is.na(Audit_Status) & Audit_Status == "growth" & !is.na(A1_2015) ~ A1_2015,
                   TRUE ~ A1),
    Known_Non_Recruit = (!is.na(Audit_Status) & Audit_Status == "occluded") |
                        (!is.na(Audit_Status) & Audit_Status == "growth" & !is.na(A1_2015))
  ) %>%
  filter(is.na(Audit_Status) | Audit_Status != "artifact") %>%
  flag_tracked()

# ------------------------------------------------------------------------------
# Ages of every colony (about 5 minutes)
# ------------------------------------------------------------------------------
colony_ages <- simulate_ages(colonies_audited) %>%
  mutate(Annual_Growth_cm2_yr = ifelse(!is.na(A1) & !is.na(A2) & A1 > 0,
                                       round((A2 - A1) / SURVEY_INTERVAL, 2), NA_real_))

cat("\nBy species:\n")
colony_ages %>%
  group_by(Species) %>%
  summarise(Total_Colonies     = n(),
            Tracked_Colonies   = sum(tracked_both, na.rm = TRUE),
            Median_Age_2022    = round(median(Age_2022_Median, na.rm = TRUE), 1),
            Mean_Growth_cm2_yr = round(mean(Annual_Growth_cm2_yr, na.rm = TRUE), 2),
            .groups = "drop") %>%
  as.data.frame() %>%
  print(row.names = FALSE)

write_csv(colony_ages, in_data("results", "colony_ages.csv"))
saveRDS(list(stand_age_audited = colonies_audited, param_draws = param_draws,
             A0_species = A0_species, valid_species = SPECIES),
        in_data("models", "age_simulation_inputs.rds"))
