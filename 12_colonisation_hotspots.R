# ==============================================================================
# 12_colonisation_hotspots.R
# Where on the wall did each species arrive early or late? Getis-Ord Gi* of
# colonisation year per species (8 nearest neighbours, binary weights
# including the colony itself): hotspots are clusters of late colonisation
# years, coldspots clusters of early ones.
#
# Reads:  results/colonisation_years.csv (script 10)
# Writes: results/colonisation_hotspots.csv   Gi* z-score and class per colony
#         figures/colonisation_hotspots.png   map per species
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(spdep)
  library(patchwork)
})

K_NEIGHBORS <- 8
PX_TO_M     <- 0.00216   # orthomosaic pixel size (2.16 mm)

spatial_df <- read_csv(in_data("results", "colonisation_years.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Colonization_Year_Est), !is.na(Centroid_x), !is.na(Centroid_y_flipped)) %>%
  filter(Species %in% SPECIES) %>%
  mutate(x_m = Centroid_x * PX_TO_M, y_m = Centroid_y_flipped * PX_TO_M)

gi_classes <- c("Hotspot (99% CI)", "Hotspot (95% CI)", "Hotspot (90% CI)", "Not Significant",
                "Coldspot (90% CI)", "Coldspot (95% CI)", "Coldspot (99% CI)")

compute_species_gi <- function(sp) {
  sub_df <- spatial_df %>% filter(Species == sp)
  nb     <- include.self(knn2nb(knearneigh(as.matrix(sub_df[, c("x_m", "y_m")]), k = K_NEIGHBORS)))
  gi_z   <- as.numeric(localG(sub_df$Colonization_Year_Est,
                              nb2listw(nb, style = "B", zero.policy = TRUE), zero.policy = TRUE))
  sub_df %>%
    mutate(Gi_Zscore = gi_z,
           Gi_Class  = case_when(
             gi_z >=  2.58 ~ "Hotspot (99% CI)",
             gi_z >=  1.96 ~ "Hotspot (95% CI)",
             gi_z >=  1.65 ~ "Hotspot (90% CI)",
             gi_z <= -2.58 ~ "Coldspot (99% CI)",
             gi_z <= -1.96 ~ "Coldspot (95% CI)",
             gi_z <= -1.65 ~ "Coldspot (90% CI)",
             TRUE          ~ "Not Significant"),
           Gi_Class = factor(Gi_Class, levels = gi_classes))
}

final_gi_df <- bind_rows(lapply(SPECIES, compute_species_gi))
write_csv(final_gi_df, in_data("results", "colonisation_hotspots.csv"))
cat("Getis-Ord Gi* classes of colonisation year:\n")
print(table(final_gi_df$Species, final_gi_df$Gi_Class)[SPECIES, ])

# ------------------------------------------------------------------------------
# Figure: one panel per species, same axes
# ------------------------------------------------------------------------------
x_min <- 0
x_max <- max(final_gi_df$x_m) * 1.02
y_min <- min(final_gi_df$y_m) * 0.90
y_max <- max(final_gi_df$y_m) * 1.08

hotspot_colors <- c("Hotspot (99% CI)"  = "#D73027", "Hotspot (95% CI)"  = "#FC8D59",
                    "Hotspot (90% CI)"  = "#FEE090", "Not Significant"   = "#E5E7EB",
                    "Coldspot (90% CI)" = "#E0F3F8", "Coldspot (95% CI)" = "#91BFDB",
                    "Coldspot (99% CI)" = "#4575B4")

make_species_panel <- function(sp, show_x = FALSE) {
  p_df <- final_gi_df %>% filter(Species == sp)
  ggplot() +
    geom_point(data = filter(p_df, Gi_Class == "Not Significant"), aes(x = x_m, y = y_m),
               color = "#D1D5DB", alpha = 0.4, size = 1.2, stroke = 0) +
    geom_point(data = filter(p_df, Gi_Class != "Not Significant"), aes(x = x_m, y = y_m, fill = Gi_Class),
               shape = 21, color = "#111827", stroke = 0.35, size = 2.2, alpha = 0.95,
               show.legend = TRUE) +   # a key for every class, also those absent from this panel
    scale_fill_manual(values = hotspot_colors, drop = FALSE,
                      name = "Getis\u2013Ord Gi* cluster of colonisation year (confidence)") +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
    annotate("text", x = x_min + 0.01 * (x_max - x_min), y = y_max - 0.08 * (y_max - y_min),
             label = sp, fontface = "bold.italic", size = 3.3, hjust = 0, color = "black") +
    labs(x = if (show_x) "Distance along wall (m)" else NULL, y = "Height on wall (m)") +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(
      axis.title       = element_text(size = 8.5, face = "bold", color = "black"),
      axis.text        = element_text(size = 7.5, color = "black"),
      axis.line        = element_line(linewidth = 0.4, color = "black"),
      axis.ticks       = element_line(linewidth = 0.4, color = "black"),
      panel.background = element_rect(fill = "#FAFAFA", color = "#D1D5DB"),
      plot.margin      = margin(t = 3, r = 8, b = 3, l = 5, unit = "pt")
    )
}

figure <- (make_species_panel(SPECIES[1]) / make_species_panel(SPECIES[2]) /
             make_species_panel(SPECIES[3], show_x = TRUE)) +
  plot_layout(guides = "collect") &
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.title     = element_text(size = 8, face = "bold"),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.35, "cm"),
        legend.margin    = margin(t = 5, b = 2, unit = "pt"),
        plot.tag         = element_text(size = 11, face = "bold", family = "sans")) &
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5,
                             override.aes = list(size = 3)))
figure <- figure + plot_annotation(tag_levels = "A")

save_figure(figure, "colonisation_hotspots.png", width = 178, height = 160)
