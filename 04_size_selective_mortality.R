# ==============================================================================
# 04_size_selective_mortality.R
# Fig. 3: does mortality depend on colony size? Compares the 2015 size of
# colonies that died, shrank (lost more area than their MDC95) or survived,
# with Wilcoxon rank-sum tests.
#
# Reads:  processed/colonies_with_mdc.xlsx   (script 02)
# Writes: figures/mortality_by_size.png       (Fig. 3)
#         results/mortality_size_tests.csv    (Wilcoxon tests, died vs survived)
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

FATES <- c("Total Mortality", "Partial Mortality (Shrinkage)", "Survived (Stable/Grew)")
fate_colors <- c("Survived (Stable/Grew)"        = "steelblue",
                 "Partial Mortality (Shrinkage)" = "darkorange",
                 "Total Mortality"               = "darkred")
fate_labels <- c("Survived (Stable/Grew)"        = "Survived (stable or grew)",
                 "Partial Mortality (Shrinkage)" = "Partial mortality",
                 "Total Mortality"               = "Whole-colony mortality")

# Cells with fewer colonies than this are drawn as points only (no violin)
MIN_N_FOR_VIOLIN <- 5

# ------------------------------------------------------------------------------
# Fate of every identified colony present in 2015. Colonies hidden in 2022
# are censored rather than counted as dead.
# ------------------------------------------------------------------------------
colonies <- read_excel(in_data("processed", "colonies_with_mdc.xlsx"))
cat("Censored", sum(colonies$state_transition == "occluded"), "colonies hidden in 2022.\n")

fates <- colonies %>%
  filter(state_transition != "occluded", present_2015, Species != "Coral Recruit") %>%
  mutate(
    Species = factor(standard_species(Species), levels = SPECIES),
    is_dead = if_else(present_2022, 0, 1),
    mortality_type = case_when(
      !present_2022                              ~ "Total Mortality",
      change_area < 0 & detectable_change %in% TRUE ~ "Partial Mortality (Shrinkage)",
      TRUE                                       ~ "Survived (Stable/Grew)"
    ),
    mortality_type = factor(mortality_type, levels = FATES)
  ) %>%
  arrange(desc(mortality_type))      # survived, shrank, died

# ------------------------------------------------------------------------------
# Wilcoxon rank-sum tests on 2015 area
# ------------------------------------------------------------------------------
p_to_stars <- function(p) {
  as.character(cut(p, breaks = c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "ns")))
}

# Died vs survived (shrinking colonies count as survivors), per species
mortality_tests <- fates %>%
  group_by(Species) %>%
  summarise(
    n_died               = sum(is_dead == 1),
    n_survived           = sum(is_dead == 0),
    median_area_died     = median(Area_2015[is_dead == 1]),
    median_area_survived = median(Area_2015[is_dead == 0]),
    W_stat               = wilcox.test(Area_2015 ~ is_dead)$statistic,
    p_value              = wilcox.test(Area_2015 ~ is_dead)$p.value,
    .groups              = "drop"
  ) %>%
  mutate(sig_stars = p_to_stars(p_value))

cat("\nWhole-colony mortality vs 2015 area (Wilcoxon rank-sum):\n")
print(as.data.frame(mortality_tests), row.names = FALSE)
write_csv(mortality_tests, in_data("results", "mortality_size_tests.csv"))

# Shrank vs stable or grew, M. oculata only (too few shrinking colonies otherwise)
shrink_test <- wilcox.test(Area_2015 ~ mortality_type,
                           data = filter(fates, Species == "Madrepora oculata",
                                         mortality_type != "Total Mortality"))
cat("\nM. oculata partial mortality vs 2015 area: W =", shrink_test$statistic,
    "| p =", shrink_test$p.value, "\n")
shrink_stars <- p_to_stars(shrink_test$p.value)

# ------------------------------------------------------------------------------
# Figure: violin + inner boxplot + every colony as a point
# ------------------------------------------------------------------------------
# n labels (with significance stars) in one column right of the largest colony
label_x <- max(fates$Area_2015) * 1.7
n_labels <- fates %>%
  group_by(Species, mortality_type, .drop = FALSE) %>%
  summarise(n_count = n(), .groups = "drop") %>%
  left_join(select(mortality_tests, Species, sig_stars), by = "Species") %>%
  mutate(
    sig = case_when(
      mortality_type == "Total Mortality" & sig_stars != "ns" ~ sig_stars,
      Species == "Madrepora oculata" &
        mortality_type == "Partial Mortality (Shrinkage)" & shrink_stars != "ns" ~ shrink_stars,
      TRUE ~ ""
    ),
    label_text = if_else(sig == "", paste0("n = ", n_count), paste0("n = ", n_count, " ", sig))
  )

cell_n <- count(fates, Species, mortality_type, name = "n_count")
cat("\nDrawn as points only (n <", MIN_N_FOR_VIOLIN, "):\n")
print(as.data.frame(filter(cell_n, n_count < MIN_N_FOR_VIOLIN)), row.names = FALSE)
violin_data <- fates %>%
  left_join(cell_n, by = c("Species", "mortality_type")) %>%
  filter(n_count >= MIN_N_FOR_VIOLIN)

strip_labeller <- as_labeller(
  c("Madrepora oculata"     = "bolditalic('Madrepora oculata')",
    "Desmophyllum pertusum" = "bolditalic('Desmophyllum pertusum')",
    "Primnoa msp."          = "bolditalic('Primnoa')~bold('msp.')"),
  label_parsed
)

fig3 <- ggplot(fates, aes(x = Area_2015, y = mortality_type)) +
  geom_violin(data = violin_data, aes(fill = mortality_type),
              alpha = 0.55, color = "black", linewidth = 0.4, width = 0.7, trim = TRUE) +
  geom_boxplot(data = violin_data, width = 0.12, fill = "white", outlier.shape = NA,
               alpha = 0.95, color = "black", linewidth = 0.35) +
  geom_point(aes(color = mortality_type), alpha = 0.4, size = 1.1,
             position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  geom_text(data = n_labels, aes(x = label_x, y = mortality_type, label = label_text),
            inherit.aes = FALSE, hjust = 0, size = 8 / .pt, color = "black") +
  facet_wrap(~ Species, ncol = 1, labeller = strip_labeller) +
  scale_x_log10(breaks = c(1, 10, 100, 1000), labels = scales::label_comma(),
                expand = expansion(mult = c(0.04, 0.16))) +
  scale_y_discrete(drop = FALSE, labels = fate_labels) +
  scale_fill_manual(values = fate_colors) +
  scale_color_manual(values = fate_colors) +
  theme_bw(base_size = 9) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "plain", size = 9, color = "black"),
    axis.text          = element_text(size = 8.5, color = "black"),
    axis.title         = element_text(size = 9.5, face = "bold", color = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing      = unit(0.8, "lines"),
    plot.margin        = margin(t = 6, r = 8, b = 6, l = 6, unit = "pt"),
    axis.title.y       = element_blank()
  ) +
  labs(x = expression(bold(paste("2015 planar area (cm"^2, ")"))))   # log axis

save_figure(fig3, "mortality_by_size.png", width = 174, height = 180)
