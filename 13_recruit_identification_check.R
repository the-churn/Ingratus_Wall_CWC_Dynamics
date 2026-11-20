# ==============================================================================
# 13_recruit_identification_check.R
# Is the 2015-2022 shift in colony size an artefact of how small recruits were
# identified? (ESM Methods S1, Fig. S1). Small scleractinian colonies are hard
# to identify to species and some were labelled "Coral Recruit". All
# scleractinian colonies up to a size threshold (and every "Coral Recruit",
# whatever its size) are pooled into one cohort; larger colonies stay in their
# species. Colony sizes in 2015 and 2022 are then compared within each cohort
# (two-sample Kolmogorov-Smirnov tests on ln(area)).
#
# Threshold: the 90th percentile of the areas labelled "Coral Recruit", in
# each year, then the larger of the two (robust to a single mislabelled colony).
#
# Reads:  processed/colony_outlines.csv, processed/colonies_tracked.csv (script 01)
# Writes: results/recruit_threshold.csv
#         results/recruit_identification_ks.csv
#         figures/recruit_identification_check.png   (Fig. S1)
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages(library(tidyverse))

POOLED <- "Pooled scleractinian (\u2264 threshold)"

# ------------------------------------------------------------------------------
# Class of each colony in each year: a colony labelled "Coral Recruit" in 2015
# may have been identified to species by 2022, so the years are kept apart
# ------------------------------------------------------------------------------
species_by_year <- read_csv(in_data("processed", "colony_outlines.csv"), show_col_types = FALSE,
                            col_types = cols(Join_ID = col_character())) %>%
  group_by(Join_ID, Year) %>%
  summarise(Species = first(TagLab.Class.name), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = Species, names_prefix = "Species_")

tracked <- read_csv(in_data("processed", "colonies_tracked.csv"), show_col_types = FALSE,
                    col_types = cols(TagLab.Genet.Id = col_character())) %>%
  left_join(species_by_year, by = c("TagLab.Genet.Id" = "Join_ID")) %>%
  mutate(Species_2015_clean = standard_species(Species_2015),
         Species_2022_clean = standard_species(Species_2022))
stopifnot(!any(is.na(tracked$Species_2015) & is.na(tracked$Species_2022)))

# ------------------------------------------------------------------------------
# Size threshold
# ------------------------------------------------------------------------------
recruit_areas_2015 <- tracked$Area_2015[tracked$Species_2015_clean %in% "Coral Recruit"]
recruit_areas_2022 <- tracked$Area_2022[tracked$Species_2022_clean %in% "Coral Recruit"]

threshold_summary <- tibble(
  n_recruit_2015   = sum(!is.na(recruit_areas_2015)),
  n_recruit_2022   = sum(!is.na(recruit_areas_2022)),
  max_recruit_2015 = max(recruit_areas_2015, na.rm = TRUE),
  max_recruit_2022 = max(recruit_areas_2022, na.rm = TRUE),
  p90_recruit_2015 = quantile(recruit_areas_2015, 0.90, na.rm = TRUE, names = FALSE),
  p90_recruit_2022 = quantile(recruit_areas_2022, 0.90, na.rm = TRUE, names = FALSE)
) %>%
  mutate(recruit_threshold = max(p90_recruit_2015, p90_recruit_2022))
recruit_threshold <- threshold_summary$recruit_threshold

cat("Colonies labelled 'Coral Recruit' (cm2):\n")
with(threshold_summary, {
  cat("  2015: n =", n_recruit_2015, "| p90 =", round(p90_recruit_2015, 2), "| max =", round(max_recruit_2015, 2), "\n")
  cat("  2022: n =", n_recruit_2022, "| p90 =", round(p90_recruit_2022, 2), "| max =", round(max_recruit_2022, 2), "\n")
})
cat("Pooling threshold (larger p90):", round(recruit_threshold, 2), "cm2\n\n")
write_excel_csv(threshold_summary, in_data("results", "recruit_threshold.csv"))

# ------------------------------------------------------------------------------
# Cohorts and K-S tests (ln area, 2015 vs 2022)
# ------------------------------------------------------------------------------
scleractinians <- c("Desmophyllum pertusum", "Madrepora oculata", "Coral Recruit")
demo_pop <- bind_rows(
  tracked %>% filter(!is.na(Area_2015), Area_2015 > 0) %>%
    transmute(TagLab.Genet.Id, Year = 2015L, Area = Area_2015, Species = Species_2015_clean),
  tracked %>% filter(!is.na(Area_2022), Area_2022 > 0) %>%
    transmute(TagLab.Genet.Id, Year = 2022L, Area = Area_2022, Species = Species_2022_clean)
) %>%
  mutate(cohort = case_when(
           Species == "Coral Recruit"                                 ~ POOLED,   # any size
           Species %in% scleractinians & Area <= recruit_threshold   ~ POOLED,
           Species %in% scleractinians                               ~ Species,
           TRUE                                                      ~ NA_character_),   # Primnoa
         ln_Area = log(Area)) %>%
  filter(!is.na(cohort))

cat("Colonies per cohort and year:\n")
print(demo_pop %>% count(cohort, Year) %>% pivot_wider(names_from = Year, values_from = n, values_fill = 0))

run_ks <- function(df) {
  x2015 <- df$ln_Area[df$Year == 2015]
  x2022 <- df$ln_Area[df$Year == 2022]
  has_ties <- length(unique(x2015)) < length(x2015) || length(unique(x2022)) < length(x2022) ||
    length(intersect(x2015, x2022)) > 0
  res <- suppressWarnings(ks.test(x2015, x2022))
  tibble(n_2015 = length(x2015), n_2022 = length(x2022),
         D = unname(res$statistic), p_value = res$p.value,
         note = if (has_ties) "Ties present -- asymptotic p-value (see ?ks.test)" else "")
}

ks_results <- demo_pop %>%
  group_by(cohort) %>%
  group_modify(~ run_ks(.x)) %>%
  ungroup() %>%
  mutate(across(c(D, p_value), ~ round(.x, 4)))

cat("\nK-S tests, ln(area) 2015 vs 2022, by cohort:\n")
print(ks_results, width = Inf)
write_excel_csv(ks_results, in_data("results", "recruit_identification_ks.csv"))   # BOM keeps the \u2264 readable in Excel

# ------------------------------------------------------------------------------
# Fig. S1: size-frequency distributions by cohort and year
# ------------------------------------------------------------------------------
fmt3 <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.3f", x))
cohort_levels <- c(POOLED, setdiff(unique(demo_pop$cohort), POOLED))

plot_df <- demo_pop %>% mutate(cohort = factor(cohort, levels = cohort_levels), Year_f = factor(Year))
ann_df  <- ks_results %>%
  mutate(cohort = factor(cohort, levels = cohort_levels),
         label  = paste0("D = ", fmt3(D), "\np ", ifelse(p_value < 0.001, "< 0.001", paste0("= ", fmt3(p_value)))))

bw <- diff(range(plot_df$ln_Area)) / 20      # one bin width for all panels
year_colors <- c("2015" = "#3B0F4B", "2022" = "#F69336")

fig_s1 <- ggplot(plot_df, aes(x = ln_Area, fill = Year_f, color = Year_f)) +
  geom_histogram(position = "identity", alpha = 0.55, binwidth = bw, linewidth = 0.2) +
  geom_text(data = ann_df, aes(label = label), x = -Inf, y = Inf,
            hjust = -0.08, vjust = 1.15, inherit.aes = FALSE, size = 2.6, lineheight = 0.9) +
  facet_wrap(~cohort, scales = "free_y", ncol = length(cohort_levels)) +
  scale_fill_manual(values = year_colors, name = "Survey year") +
  scale_color_manual(values = year_colors, name = "Survey year") +
  labs(x = expression(paste("ln(Colony Planar Area, cm"^2, ")")), y = "Number of colonies") +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 7, face = "italic"),
    axis.title       = element_text(size = 8, color = "black"),
    axis.text        = element_text(size = 7, color = "black"),
    axis.line        = element_line(linewidth = 0.4, color = "black"),
    axis.ticks       = element_line(linewidth = 0.4, color = "black"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    plot.margin      = margin(t = 5, r = 8, b = 5, l = 5, unit = "pt")
  )
save_figure(fig_s1, "recruit_identification_check.png", width = 178, height = 70)
