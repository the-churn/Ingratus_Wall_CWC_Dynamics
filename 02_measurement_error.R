# ==============================================================================
# 02_measurement_error.R
# Models the measurement error of colony area from repeat CNN segmentations,
# and gives each colony its own 95% minimum detectable change (MDC95).
#
# Reads:  CNN Error.xlsx                    repeat segmentations of the same colonies
#         processed/colonies_tracked.csv    (script 01)
# Writes: processed/colonies_with_mdc.xlsx  the colony table used by all later scripts
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
})

# ------------------------------------------------------------------------------
# Error model: SD of a colony's area across repeat segmentations, as a function
# of its mean area and species
# ------------------------------------------------------------------------------
# Species names as in the paper: the file labels D. pertusum both
# "Desmophyllum pertusum" (timepoint t1) and "D. pertusum" (t2)
repeats <- read_excel(in_data("CNN Error.xlsx")) %>%
  mutate(species = standard_species(species)) %>%
  group_by(colony_id, species, timepoint) %>%
  summarise(mean_area = mean(area, na.rm = TRUE),
            sd_area   = sd(area, na.rm = TRUE),
            .groups   = "drop")

error_model <- lm(sd_area ~ mean_area + species, data = repeats)
cat("Measurement error model:\n")
print(summary(error_model))

# ------------------------------------------------------------------------------
# State of each colony between the surveys
# ------------------------------------------------------------------------------
colonies <- read_csv(in_data("processed", "colonies_tracked.csv"), show_col_types = FALSE) %>%
  mutate(
    A1 = Area_2015,                       # planar area in 2015 (cm2)
    A2 = Area_2022,                       # planar area in 2022 (cm2)
    present_2015 = !is.na(A1) & A1 > 0,
    present_2022 = !is.na(A2) & A2 > 0,
    state_transition = case_when(
      fate == "Occluded"                                     ~ "occluded",
      fate == "Occluded 2015"                                ~ "occluded_2015",
      fate == "Total Mortality"                              ~ "mortality",
      fate == "Recruit"                                      ~ "recruitment",
      fate %in% c("Survivor", "Fission", "Fusion", "Complex") ~ "persistence",
      TRUE                                                   ~ "unknown"
    )
  )

# ------------------------------------------------------------------------------
# MDC95 = 1.96 x the SD of the change, from the predicted error of each area
# ------------------------------------------------------------------------------
# Error-model species of each colony; unidentified recruits ("Coral Recruit")
# get no MDC95
error_species <- factor(ifelse(colonies$Species == "Coral Recruit", NA, standard_species(colonies$Species)),
                        levels = levels(factor(repeats$species)))
predicted_sd <- function(area) {
  predict(error_model, newdata = data.frame(mean_area = area, species = error_species))
}

# No MDC95 for unidentified recruits or for colonies hidden in either survey
no_mdc <- colonies$Species == "Coral Recruit" |
  colonies$state_transition %in% c("occluded", "occluded_2015")

colonies <- colonies %>%
  mutate(
    sigma1       = predicted_sd(A1),
    sigma2       = predicted_sd(A2),
    change_area  = A2 - A1,
    sigma_change = ifelse(no_mdc, NA, sqrt(sigma1^2 + sigma2^2)),
    MDC95        = ifelse(no_mdc, NA, 1.96 * sigma_change),
    # Only colonies present in both surveys can show a detectable change
    detectable_change = ifelse(state_transition == "persistence" & !no_mdc,
                               abs(change_area) > MDC95, NA),
    # Recomputed at full precision (the copy read from the csv is rounded)
    Annual_Area_Change = ifelse(state_transition %in% c("occluded", "occluded_2015"),
                                NA_real_, change_area / SURVEY_INTERVAL),
    MDC95_Annual = MDC95 / SURVEY_INTERVAL
  )

write_xlsx(colonies, in_data("processed", "colonies_with_mdc.xlsx"))

cat("\nPersistent colonies whose change exceeds their MDC95:\n")
colonies %>%
  filter(state_transition == "persistence", !no_mdc) %>%
  count(Species = standard_species(Species), detectable_change) %>%
  print()
