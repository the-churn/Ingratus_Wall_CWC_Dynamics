# ==============================================================================
# 07_growth_curves.R
# Fig. 4: relative growth rate against 2015 colony size, with the fitted
# growth curve of each species, its asymptote (RGR-infinity) and A_thresh, the
# size at which the curve has covered 95% of the way from R0 to the asymptote:
#   A_thresh = -log(0.05) / exp(lrc) = ln(20) / exp(lrc)
#
# Reads:  processed/colonies_with_mdc.xlsx, models/growth_model.rds (script 05)
# Writes: models/athresh_sizes.rds     A_thresh per species (scripts 09, 14, 18)
#         figures/growth_curves.png    (Fig. 4)
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
})

FONT_FAMILY <- "sans"

df_growth     <- read_growth_data()
m_growth_gnls <- readRDS(in_data("models", "growth_model.rds"))

# ------------------------------------------------------------------------------
# Species-level parameters (M. oculata is the model's reference level)
# ------------------------------------------------------------------------------
cf <- coef(m_growth_gnls)
species_coef <- function(par) {
  ref <- cf[[paste0(par, ".(Intercept)")]]
  c(ref,
    ref + cf[[paste0(par, ".fSpecies", SPECIES[2])]],
    ref + cf[[paste0(par, ".fSpecies", SPECIES[3])]])
}
asym    <- species_coef("Asym")
athresh <- setNames(-log(0.05) / exp(species_coef("lrc")), SPECIES)

saveRDS(athresh, in_data("models", "athresh_sizes.rds"))
cat("A_thresh (cm2):\n")
print(round(athresh, 1))

# Asymptote labels (plotmath), above (vjust < 0) or below (> 1) their own
# dashed line, so that no label crosses another species' line
asym_df <- data.frame(
  fSpecies    = factor(SPECIES, levels = SPECIES),
  asym_val    = asym,
  label_text  = sprintf("bold(RGR[infinity]~'= %.3f yr'^-1)", asym),
  label_vjust = c(1.3, -0.4, 1.3)
)
size_thresholds <- data.frame(
  fSpecies   = factor(SPECIES, levels = SPECIES),
  asym_size  = unname(athresh),
  label_text = paste0(round(unname(athresh), 1), " cm\u00B2")
)

# Fitted curves across the observed size range
new_data <- expand.grid(A1 = seq(0, max(df_growth$A1), length.out = 1000),
                        fSpecies = factor(SPECIES, levels = SPECIES)) %>%
  mutate(pred = predict(m_growth_gnls, newdata = .))

# ------------------------------------------------------------------------------
# Figure. Axis limits: x 0-300 cm2, y -0.22 to 0.52 (6 Primnoa colonies of
# 332-718 cm2 and 2 colonies with extreme RGR fall outside; see caption)
# ------------------------------------------------------------------------------
species_labels_expr <- c("Madrepora oculata"     = expression(italic("M. oculata")),
                         "Desmophyllum pertusum" = expression(italic("D. pertusum")),
                         "Primnoa msp."          = expression(italic("Primnoa") ~ "msp."))
y_breaks <- seq(-0.2, 0.5, by = 0.1)

fig4 <- ggplot(df_growth, aes(x = A1, y = RGR_Planar, color = fSpecies)) +
  geom_point(alpha = 0.30, size = 1.2, stroke = 0) +
  geom_line(data = new_data, aes(y = pred), linewidth = 0.9, show.legend = TRUE) +
  geom_hline(data = asym_df, aes(yintercept = asym_val, color = fSpecies),
             linetype = "dashed", linewidth = 0.4, alpha = 0.8, show.legend = FALSE) +
  geom_vline(data = size_thresholds, aes(xintercept = asym_size, color = fSpecies),
             linetype = "dotted", linewidth = 0.5, alpha = 0.85, show.legend = FALSE) +
  geom_text(data = size_thresholds,
            aes(x = asym_size, y = 0.48, label = label_text, color = fSpecies),
            angle = 90, vjust = -0.4, hjust = 1, size = 2.5, fontface = "bold",
            show.legend = FALSE) +
  geom_text(data = asym_df,
            aes(x = 295, y = asym_val, label = label_text, color = fSpecies, vjust = label_vjust),
            hjust = 1, size = 2.5, parse = TRUE, show.legend = FALSE) +
  scale_color_manual(values = SPECIES_COLOURS, labels = species_labels_expr) +
  scale_x_continuous(breaks = seq(0, 300, by = 50)) +
  scale_y_continuous(breaks = y_breaks, labels = sub("^-", "\u2212", sprintf("%.1f", y_breaks))) +
  coord_cartesian(xlim = c(0, 300), ylim = c(-0.22, 0.52)) +
  labs(x = expression(bold(paste("2015 planar area, ", bolditalic("A"), " (cm"^2, ")"))),
       y = expression(bold(paste("RGR (yr"^-1, ")"))),
       color = NULL) +
  theme_classic(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    legend.position  = "top",
    legend.key       = element_blank(),
    legend.key.width = unit(8, "mm"),
    legend.text      = element_text(size = 8, face = "plain"),
    axis.title       = element_text(face = "bold", size = 8.5),
    axis.text        = element_text(color = "black", size = 7.5),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.4)
  )

save_figure(fig4, "growth_curves.png", width = 129, height = 100)
