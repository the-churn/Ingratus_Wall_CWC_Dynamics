# ==============================================================================
# 09_age_figure.R
# Fig. 5: how the age model works and how uncertain it is.
#   a-c  modelled age against planar area for 200 of the 5,000 parameter draws
#        (coloured) and the fitted curve (black); dashed line at A_thresh
#   d-f  modelled age at A_thresh across all 5,000 draws (median, 95% interval)
#   g    size at age of the three species over the first 30 years (fitted
#        curve, band with 95% of the draws, point at the age of A_thresh)
# Also prints the size reached at 5, 10 and 20 years and the age at 100 cm2.
#
# Reads:  models/growth_model.rds (script 05), models/athresh_sizes.rds (script 07),
#         models/age_simulation_inputs.rds (script 08)
# Writes: figures/age_model.png (Fig. 5)
# ==============================================================================

source("setup.R")
source("age_model.R")     # species_params(), age_at_size()
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(nlme)
})

FONT_FAMILY      <- "sans"
N_SIMS_SPAGHETTI <- 200            # draws shown in panels a-c
BAND_GRID_N      <- 800            # sizes at which the panel g bands are computed
G_MAX_AGE        <- 30             # panel g x-axis (yr)
G_MAX_AREA       <- 200            # panel g y-axis (cm2)
REPORT_AGES      <- c(5, 10, 20)   # size at these ages is printed
REPORT_AREA      <- 100            # age at this size is printed

param_means <- coef(readRDS(in_data("models", "growth_model.rds")))
athresh     <- readRDS(in_data("models", "athresh_sizes.rds"))[SPECIES]
inputs      <- readRDS(in_data("models", "age_simulation_inputs.rds"))
draws_full  <- inputs$param_draws

# Detection floor and largest observed colony of each species
A0_list       <- list()
max_size_list <- list()
for (sp in SPECIES) {
  A0_list[[sp]] <- unname(inputs$A0_species[sp])
  d <- inputs$stand_age_audited
  sizes <- c(d$A1[d$Species == sp], d$A2[d$Species == sp])
  max_size_list[[sp]] <- max(sizes[!is.na(sizes) & sizes > 0])
}

species_titles <- list(
  "Madrepora oculata"     = quote(bolditalic("Madrepora oculata")),
  "Desmophyllum pertusum" = quote(bolditalic("Desmophyllum pertusum")),
  "Primnoa msp."          = quote(bolditalic("Primnoa") * bold(" msp."))
)

# ------------------------------------------------------------------------------
# Ages for all draws at once on a log-spaced size grid (for panel g): the age
# integral by the trapezoid rule in ln(area), t = integral of d(ln a) / RGR(a).
# Agrees with integrate() to about 0.01 yr. Returns the grid and a
# draws x sizes matrix of ages (NA where RGR <= 0).
# ------------------------------------------------------------------------------
age_grid_all_draws <- function(P, A0, A_max, n = BAND_GRID_N) {
  A   <- exp(seq(log(A0), log(A_max), length.out = n))
  du  <- diff(log(A))[1]
  rgr <- P[, "asym"] + (P[, "r0"] - P[, "asym"]) * exp(-outer(exp(P[, "lrc"]), A))
  f   <- 1 / rgr
  f[rgr <= 0] <- NA
  ages <- t(apply(f, 1, function(v) c(0, cumsum((v[-1] + v[-length(v)]) / 2 * du))))
  list(A = A, ages = ages)
}

panel_theme <- theme_classic(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    plot.tag            = element_text(face = "bold", size = 10, vjust = 1),
    plot.title.position = "panel",
    axis.title          = element_text(size = 8.5, face = "bold", color = "black"),
    axis.text           = element_text(size = 8, color = "black"),
    axis.line           = element_line(linewidth = 0.3, color = "black"),
    axis.ticks          = element_line(linewidth = 0.3, color = "black"),
    plot.margin         = margin(t = 4, r = 4, b = 4, l = 4)
  )

# ------------------------------------------------------------------------------
# Panels a-f, one species at a time
# ------------------------------------------------------------------------------
plots_row1 <- list()
plots_row2 <- list()
band_store <- list()
thresh_pts <- list()

for (i in seq_along(SPECIES)) {
  sp       <- SPECIES[i]
  A_thresh <- athresh[[sp]]
  A0       <- A0_list[[sp]]
  sp_color <- SPECIES_COLOURS[[sp]]

  # a-c: age against area, fitted curve and 200 draws
  size_grid  <- exp(seq(log(A0), log(max_size_list[[sp]]), length.out = 60))
  P_fit  <- species_params(param_means, sp)[1, ]
  P_all  <- species_params(draws_full, sp)
  curve_mean <- data.frame(Size = size_grid,
                           Age  = sapply(size_grid, age_at_size, p = P_fit, A0 = A0))
  spaghetti_df <- bind_rows(lapply(seq_len(N_SIMS_SPAGHETTI), function(k) {
    data.frame(Draw = k, Size = size_grid,
               Age  = sapply(size_grid, age_at_size, p = P_all[k, ], A0 = A0))
  }))
  band_store[[sp]] <- age_grid_all_draws(P_all, A0, max_size_list[[sp]])

  plots_row1[[sp]] <- ggplot(spaghetti_df, aes(Size, Age)) +
    geom_line(aes(group = Draw), color = sp_color, alpha = 0.08, linewidth = 0.35) +
    geom_line(data = curve_mean, color = "black", linewidth = 0.7) +
    geom_vline(xintercept = A_thresh, linetype = "dashed", color = "grey20", linewidth = 0.4) +
    labs(tag = letters[i], title = species_titles[[sp]],
         x = expression(bold(paste("Planar area (cm"^2, ")"))),
         y = expression(bold("Modelled age (yr)"))) +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4), labels = scales::label_comma()) +
    panel_theme +
    theme(plot.title = element_text(size = 9, hjust = 0))

  # d-f: age at A_thresh over all draws
  ages_at_thresh <- sapply(seq_len(nrow(P_all)), function(k) age_at_size(A_thresh, P_all[k, ], A0))
  ages_at_thresh <- ages_at_thresh[!is.na(ages_at_thresh)]
  age_med <- median(ages_at_thresh)
  age_lo  <- unname(quantile(ages_at_thresh, 0.025))
  age_hi  <- unname(quantile(ages_at_thresh, 0.975))
  cat(sprintf("%-22s age at A_thresh (%.1f cm2): %.1f yr [%.1f-%.1f], %d draws\n",
              sp, A_thresh, age_med, age_lo, age_hi, length(ages_at_thresh)))

  plots_row2[[sp]] <- ggplot(data.frame(Age = ages_at_thresh), aes(Age)) +
    geom_histogram(fill = sp_color, alpha = 0.75, bins = 35, color = "white", linewidth = 0.15) +
    geom_vline(xintercept = age_med, color = "black", linewidth = 0.6) +
    geom_vline(xintercept = c(age_lo, age_hi), color = "black", linetype = "dotted", linewidth = 0.5) +
    labs(tag = letters[i + 3],
         title    = sprintf("Age at %.1f cm\u00B2", A_thresh),
         subtitle = sprintf("%.1f yr [%.1f\u2013%.1f]", age_med, age_lo, age_hi),
         x = expression(bold("Modelled age (yr)")),
         y = expression(bold("Number of draws"))) +
    panel_theme +
    theme(plot.title    = element_text(face = "bold", size = 9, hjust = 0),
          plot.subtitle = element_text(color = "grey30", size = 8.5, face = "bold", hjust = 0))

  thresh_pts[[sp]] <- data.frame(Species = sp, Age = age_med, Area = A_thresh)
}

# ------------------------------------------------------------------------------
# Panel g: size at age, first 30 years, all species on common axes
# ------------------------------------------------------------------------------
age_axis <- seq(0, G_MAX_AGE, by = 0.1)
size_at_age <- bind_rows(lapply(SPECIES, function(sp) {
  g_all <- band_store[[sp]]
  sizes <- apply(g_all$ages, 1, function(tt) approx(tt, g_all$A, xout = age_axis)$y)   # ages x draws
  fit_A <- exp(seq(log(A0_list[[sp]]), log(max_size_list[[sp]]), length.out = BAND_GRID_N))
  fit_t <- age_grid_all_draws(species_params(param_means, sp), A0_list[[sp]],
                              max_size_list[[sp]])$ages[1, ]
  data.frame(Species = sp, Age = age_axis,
             fit = approx(fit_t, fit_A, xout = age_axis)$y,
             lo  = apply(sizes, 1, quantile, probs = 0.025, na.rm = TRUE),
             hi  = apply(sizes, 1, quantile, probs = 0.975, na.rm = TRUE))
}))
size_at_age$Species <- factor(size_at_age$Species, levels = SPECIES)
thresh_df <- bind_rows(thresh_pts)
thresh_df$Species <- factor(thresh_df$Species, levels = SPECIES)

cat("\nSize (cm2) reached at selected ages, fitted curve [95% range across draws]:\n")
for (sp in SPECIES) {
  r <- size_at_age[size_at_age$Species == sp & round(size_at_age$Age, 1) %in% REPORT_AGES, ]
  cat(sprintf("  %-22s %s\n", sp, paste(sprintf("%g yr: %.1f [%.1f-%.1f]", r$Age, r$fit, r$lo, r$hi),
                                       collapse = "; ")))
}
cat(sprintf("\nAge (yr) at %g cm2, fitted curve [95%% range across draws]:\n", REPORT_AREA))
for (sp in SPECIES) {
  g_all <- band_store[[sp]]
  a_all <- apply(g_all$ages, 1, function(tt) approx(g_all$A, tt, xout = REPORT_AREA)$y)
  cat(sprintf("  %-22s %.1f [%.1f-%.1f]\n", sp,
              age_at_size(REPORT_AREA, species_params(param_means, sp)[1, ], A0_list[[sp]]),
              quantile(a_all, 0.025, na.rm = TRUE), quantile(a_all, 0.975, na.rm = TRUE)))
}

species_labels <- c("Madrepora oculata"     = expression(italic("M. oculata")),
                    "Desmophyllum pertusum" = expression(italic("D. pertusum")),
                    "Primnoa msp."          = expression(italic("Primnoa") * " msp."))

p_g <- ggplot(size_at_age, aes(Age)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = Species), alpha = 0.25) +
  geom_line(aes(y = fit, color = Species), linewidth = 0.8) +
  geom_point(data = thresh_df, aes(Age, Area, color = Species, shape = "thresh"), size = 1.8) +
  scale_color_manual(values = SPECIES_COLOURS, labels = species_labels, name = NULL) +
  scale_fill_manual(values = SPECIES_COLOURS, labels = species_labels, name = NULL) +
  scale_shape_manual(values = c(thresh = 16), labels = c(thresh = expression("Age at " * A[thresh])),
                     name = NULL) +
  guides(color = guide_legend(order = 1, override.aes = list(shape = NA)),
         fill  = guide_legend(order = 1),
         shape = guide_legend(order = 2, override.aes = list(color = "grey20"))) +
  coord_cartesian(xlim = c(0, G_MAX_AGE), ylim = c(0, G_MAX_AREA), expand = TRUE) +
  scale_x_continuous(breaks = seq(0, G_MAX_AGE, by = 5)) +
  labs(tag = "g", title = "All three species, first 30 years",
       x = expression(bold("Modelled age (yr)")),
       y = expression(bold(paste("Planar area (cm"^2, ")")))) +
  panel_theme +
  theme(plot.title       = element_text(face = "bold", size = 9, hjust = 0),
        legend.text      = element_text(size = 8.5, hjust = 0),
        legend.key.width = unit(7, "mm"),
        legend.spacing.y = unit(1, "mm"))

# ------------------------------------------------------------------------------
# Layout: a-c, d-f, then g with its legend to the right
# ------------------------------------------------------------------------------
fig5 <- wrap_plots(plots_row1[[1]], plots_row1[[2]], plots_row1[[3]],
                   plots_row2[[1]], plots_row2[[2]], plots_row2[[3]],
                   p_g, guide_area(),
                   design = "ABC\nDEF\nGGH") +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(plot.margin = margin(t = 3, r = 8, b = 3, l = 3),
        legend.position = "right", legend.justification = "left")

save_figure(fig5, "age_model.png", width = 174, height = 185)
