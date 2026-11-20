# ==============================================================================
# 06_growth_model_spatial_check.R
# Are the growth model's residuals spatially autocorrelated? (ESM Methods S2,
# Fig. S3). Neighbouring colonies share an environment, so their growth may be
# correlated, which gnls does not allow for. Global Moran's I of the normalised
# residuals, with the five-nearest-neighbour weights used for the colonisation
# test (script 10), for all colonies and for each species. If any test has
# p < 0.05, the model is refitted with an exponential spatial correlation
# structure and the two fits are compared. The refit is a check only: later
# scripts keep using the script 05 model.
#
# Reads:  processed/colonies_with_mdc.xlsx, models/growth_model.rds (script 05)
# Writes: results/growth_residuals_morans_i.csv
#         and, only when the model is refitted:
#         results/growth_spatial_refit_comparison.csv  fit, parameters, A_thresh
#         results/growth_spatial_refit_shifts.csv      coefficient shifts in SEs
#         figures/growth_spatial_refit_curves.png      both growth curves
#         figures/growth_residual_semivariograms.png   (Fig. S3)
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(spdep)
  library(patchwork)
})
options(width = 160)          # keep the printed tables on one line

K_NEIGHBOURS <- 5             # as in the colonisation test (script 10)
PX_TO_M      <- 0.00216       # orthomosaic pixel size (2.16 mm)

refit_outputs <- c(in_data("results", "growth_spatial_refit_comparison.csv"),
                   in_data("results", "growth_spatial_refit_shifts.csv"),
                   in_data("figures", "growth_spatial_refit_curves.png"),
                   in_data("figures", "growth_residual_semivariograms.png"))
# Remove refit outputs from an earlier run, so they can't be mistaken for this run's
invisible(file.remove(Filter(file.exists, refit_outputs)))

# ------------------------------------------------------------------------------
# Residuals of the growth model
# ------------------------------------------------------------------------------
df_growth     <- read_growth_data()
m_growth_gnls <- readRDS(in_data("models", "growth_model.rds"))

# The rows must be the ones the model was fitted to, in the same order
resid_response <- as.numeric(residuals(m_growth_gnls, type = "response"))
if (length(resid_response) != nrow(df_growth) ||
    !isTRUE(all.equal(as.numeric(fitted(m_growth_gnls)) + resid_response, df_growth$RGR_Planar))) {
  stop("The data don't match the saved growth model. Run script 05 again first.", call. = FALSE)
}
# corExp needs a position for every colony, and no two colonies at the same one
stopifnot(!anyNA(df_growth$Centroid_x), !anyDuplicated(df_growth[, c("Centroid_x", "Centroid_y")]))

# Normalised residuals allow for the model's variance function
df_growth <- df_growth %>%
  mutate(Resid = as.numeric(residuals(m_growth_gnls, type = "normalized")))

# ------------------------------------------------------------------------------
# Global Moran's I: all species, then each species
# ------------------------------------------------------------------------------
moran_residuals <- function(resid, x, y, scope) {
  lw <- nb2listw(knn2nb(knearneigh(cbind(x, y), k = K_NEIGHBOURS)), style = "W")
  mt <- moran.test(resid, lw)   # randomisation; one-sided (positive autocorrelation)
  data.frame(Scope = scope, n = length(resid), k = K_NEIGHBOURS,
             Moran_I    = unname(mt$estimate["Moran I statistic"]),
             Expected_I = unname(mt$estimate["Expectation"]),
             Variance   = unname(mt$estimate["Variance"]),
             Z          = unname(mt$statistic),
             p_value    = mt$p.value)
}

moran_tab <- bind_rows(
  moran_residuals(df_growth$Resid, df_growth$Centroid_x, df_growth$Centroid_y, "All species"),
  lapply(SPECIES, function(sp) {
    d <- df_growth %>% filter(fSpecies == sp)
    moran_residuals(d$Resid, d$Centroid_x, d$Centroid_y, sp)
  })
)

p_text <- function(p) ifelse(p < 0.001, "< 0.001", paste("=", sprintf("%.3f", p)))
cat("Global Moran's I of the normalised growth-model residuals (k = ", K_NEIGHBOURS, "):\n", sep = "")
for (i in seq_len(nrow(moran_tab))) {
  r <- moran_tab[i, ]
  cat(sprintf("  %-23s n = %4d   I = %7.4f   E(I) = %7.4f   z = %6.2f   p %s\n",
              r$Scope, r$n, r$Moran_I, r$Expected_I, r$Z, p_text(r$p_value)))
}

write.csv(moran_tab %>% mutate(across(c(Moran_I, Expected_I, Variance, Z), ~ signif(.x, 4)),
                               p_value = signif(p_value, 3)),
          in_data("results", "growth_residuals_morans_i.csv"), row.names = FALSE)

if (!any(moran_tab$p_value < 0.05)) {
  cat("\nNo spatial autocorrelation in the residuals (all p >= 0.05); no refit needed.\n")
  quit(save = "no")
}
cat("\nAutocorrelated (p < 0.05) in:",
    paste(moran_tab$Scope[moran_tab$p_value < 0.05], collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# Refit with an exponential spatial correlation structure
# ------------------------------------------------------------------------------
# Same variance function as the model used, starting at its estimate
w_expr <- m_growth_gnls$call$weights
vs     <- coef(m_growth_gnls$modelStruct$varStruct, unconstrained = FALSE)
w_fun  <- if (is.call(w_expr)) as.character(w_expr[[1]]) else ""
if (w_fun %in% c("varPower", "varExp") && length(vs) == 1) w_expr$value <- unname(vs)
if (w_fun == "varConstPower" && all(c("const", "power") %in% names(vs))) {
  w_expr$const <- vs[["const"]]; w_expr$power <- vs[["power"]]
}
cat("\nRefitting with corExp; variance function:", deparse(w_expr), "\n")

# Starting range: median distance to the 5 nearest neighbours; nugget 0.5
xy     <- cbind(df_growth$Centroid_x, df_growth$Centroid_y)
range0 <- median(unlist(nbdists(knn2nb(knearneigh(xy, k = K_NEIGHBOURS)), xy)))

m_spatial <- eval(bquote(gnls(
  RGR_Planar ~ SSasymp(A1, Asym, R0, lrc),
  params      = list(Asym ~ fSpecies, R0 ~ fSpecies, lrc ~ fSpecies),
  data        = df_growth,
  weights     = .(w_expr),
  correlation = corExp(value = c(range0, 0.5), form = ~ Centroid_x + Centroid_y, nugget = TRUE),
  start       = unname(coef(m_growth_gnls)),
  control     = gnlsControl(maxIter = 500, minScale = 1e-7, tolerance = 1e-4),
  na.action   = na.omit
)))

m_base   <- m_growth_gnls
cs       <- coef(m_spatial$modelStruct$corStruct, unconstrained = FALSE)
mor_base <- moran_tab[1, ]
mor_sp   <- moran_residuals(as.numeric(residuals(m_spatial, type = "normalized")),
                            df_growth$Centroid_x, df_growth$Centroid_y, "All species (spatial model)")

# Species-level parameters from the treatment contrasts, and A_thresh
species_parameters <- function(m) {
  cf <- coef(m)
  term <- function(par, sp) {
    ref <- cf[[paste0(par, ".(Intercept)")]]
    if (sp == SPECIES[1]) ref else ref + cf[[paste0(par, ".fSpecies", sp)]]
  }
  data.frame(Species = SPECIES,
             Asym    = sapply(SPECIES, function(sp) term("Asym", sp)),
             R0      = sapply(SPECIES, function(sp) term("R0", sp)),
             lrc     = sapply(SPECIES, function(sp) term("lrc", sp)),
             row.names = NULL) %>%
    mutate(A_thresh = log(20) / exp(lrc))
}
par_base <- species_parameters(m_base)
par_sp   <- species_parameters(m_spatial)

model_rows <- data.frame(
  Quantity       = c("AIC", "Log-likelihood", "Parameters (df)", "Correlation range (pixels)",
                     "Correlation range (m)", "Nugget", "Residual Moran's I (all species)",
                     "Residual Moran's I p-value"),
  Species        = "",
  Script04_model = c(AIC(m_base), as.numeric(logLik(m_base)), attr(logLik(m_base), "df"),
                     NA, NA, NA, mor_base$Moran_I, mor_base$p_value),
  Spatial_model  = c(AIC(m_spatial), as.numeric(logLik(m_spatial)), attr(logLik(m_spatial), "df"),
                     cs[["range"]], cs[["range"]] * PX_TO_M, cs[["nugget"]],
                     mor_sp$Moran_I, mor_sp$p_value)
)
par_rows <- bind_rows(lapply(c("Asym", "R0", "lrc", "A_thresh"), function(q) {
  data.frame(Quantity = if (q == "A_thresh") "A_thresh (cm2)" else q,
             Species = SPECIES, Script04_model = par_base[[q]], Spatial_model = par_sp[[q]])
}))
comparison <- bind_rows(model_rows, par_rows) %>%
  mutate(Percent_change = ifelse(Species == "", NA,
                                 100 * (Spatial_model - Script04_model) / abs(Script04_model)),
         across(c(Script04_model, Spatial_model, Percent_change), ~ signif(.x, 5))) %>%
  rename(Model_used = Script04_model)

cat("\nModel used vs model allowing spatial correlation:\n")
print(comparison %>%
        mutate(across(c(Model_used, Spatial_model, Percent_change),
                      ~ ifelse(is.na(.x), "", formatC(.x, digits = 4, format = "g")))),
      row.names = FALSE)
d_aic <- AIC(m_spatial) - AIC(m_base)
write.csv(comparison, refit_outputs[1], row.names = FALSE)

# Shift of each coefficient in units of its standard error (a percentage
# change blows up for contrasts near zero)
tt_base <- summary(m_base)$tTable
tt_sp   <- summary(m_spatial)$tTable
shift_tab <- data.frame(
  Parameter         = rownames(tt_base),
  Estimate_Script04 = tt_base[, "Value"],
  SE_Script04       = tt_base[, "Std.Error"],
  Estimate_spatial  = tt_sp[rownames(tt_base), "Value"],
  SE_spatial        = tt_sp[rownames(tt_base), "Std.Error"],
  row.names         = NULL
) %>%
  mutate(Shift_in_SEs  = (Estimate_spatial - Estimate_Script04) / SE_Script04,
         SE_change_pct = 100 * (SE_spatial - SE_Script04) / SE_Script04) %>%
  rename(Estimate_used = Estimate_Script04, SE_used = SE_Script04)
cat("\nCoefficient shifts (spatial - model used), in standard errors of the model used:\n")
print(shift_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
write.csv(shift_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 5))),
          refit_outputs[2], row.names = FALSE)

# ------------------------------------------------------------------------------
# Figure: the two growth curves on top of each other (axes as in Fig. 4)
# ------------------------------------------------------------------------------
growth_curve <- function(par, model) {
  bind_rows(lapply(seq_len(nrow(par)), function(i) {
    A <- seq(0, 300, length.out = 601)
    data.frame(Species = par$Species[i], A1 = A, Model = model,
               RGR = par$Asym[i] + (par$R0[i] - par$Asym[i]) * exp(-exp(par$lrc[i]) * A))
  }))
}
model_levels <- c("Model used (independent residuals)", "Model allowing spatial correlation")
curve_base <- growth_curve(par_base, model_levels[1])
curve_sp   <- growth_curve(par_sp, model_levels[2])

gap_tab <- curve_base %>%
  mutate(Gap = abs(RGR - curve_sp$RGR)) %>%
  group_by(Species) %>%
  summarise(Max_gap = max(Gap), .groups = "drop")
cat("\nLargest gap between the two growth curves, 0-300 cm2 (per yr):\n")
print(as.data.frame(gap_tab), row.names = FALSE)

species_expr <- c("Madrepora oculata"     = 'italic("M. oculata")',
                  "Desmophyllum pertusum" = 'italic("D. pertusum")',
                  "Primnoa msp."          = 'italic("Primnoa")~"msp."')
panel_of <- function(x) factor(species_expr[as.character(x)], levels = species_expr[SPECIES])
as_model <- function(x) factor(x, levels = model_levels)
curve_base <- curve_base %>% mutate(Panel = panel_of(Species), Model = as_model(Model))
curve_sp   <- curve_sp   %>% mutate(Panel = panel_of(Species), Model = as_model(Model))
points_df  <- df_growth %>% transmute(A1, RGR = RGR_Planar, Panel = panel_of(fSpecies))
gap_labels <- gap_tab %>%
  mutate(Panel = panel_of(Species), Label = sprintf('"Largest gap %.3f "*yr^-1', Max_gap))

y_breaks <- seq(-0.2, 0.5, by = 0.1)
p_curves <- ggplot(mapping = aes(x = A1, y = RGR)) +
  geom_point(data = points_df, colour = "grey70", alpha = 0.35, size = 0.9, stroke = 0) +
  geom_line(data = curve_base, aes(colour = Species, linetype = Model), linewidth = 1.3) +
  geom_line(data = curve_sp, aes(linetype = Model), colour = "black", linewidth = 0.45) +
  geom_text(data = gap_labels, aes(x = 292, y = 0.47, label = Label), parse = TRUE,
            hjust = 1, size = 2.8, inherit.aes = FALSE) +
  scale_colour_manual(values = SPECIES_COLOURS, guide = "none") +
  scale_linetype_manual(values = setNames(c("solid", "22"), model_levels), name = NULL,
                        guide = guide_legend(override.aes = list(colour = c("grey45", "black"),
                                                                 linewidth = c(1.3, 0.45)))) +
  scale_x_continuous(breaks = seq(0, 300, 100)) +
  scale_y_continuous(breaks = y_breaks, labels = sub("^-", "\u2212", sprintf("%.1f", y_breaks))) +
  coord_cartesian(xlim = c(0, 300), ylim = c(-0.22, 0.52)) +
  facet_wrap(~ Panel, nrow = 1, labeller = label_parsed) +
  labs(x = expression(paste("2015 planar area, ", italic("A"), " (cm"^2, ")")),
       y = expression(paste("RGR (yr"^-1, ")"))) +
  theme_classic(base_family = "sans", base_size = 9) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.line        = element_blank(),
    axis.text        = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0, size = 9),
    legend.position  = "bottom",
    legend.key.width = unit(11, "mm"),
    legend.text      = element_text(size = 8.5),
    legend.margin    = margin(t = -2, b = 0)
  )
save_figure(p_curves, "growth_spatial_refit_curves.png", width = 174, height = 72)

# ------------------------------------------------------------------------------
# Fig. S3: semivariograms, residuals vs fitted, and maps of the normalised
# residuals of both models. Without spatial structure the semivariance sits
# near 1 at every distance. Distances in metres, up to half the largest.
# ------------------------------------------------------------------------------
s3_levels <- c("Without spatial correlation", "With spatial correlation (corExp)")
df_vg <- df_growth %>%
  mutate(x_m = (Centroid_x - min(Centroid_x)) * PX_TO_M,
         y_m = (max(Centroid_y) - Centroid_y) * PX_TO_M)   # image rows run downwards, so flip
max_d <- 0.5 * max(dist(df_vg[, c("x_m", "y_m")]))
vg <- bind_rows(
  as.data.frame(Variogram(m_base, form = ~ x_m + y_m, resType = "normalized",
                          data = df_vg, maxDist = max_d)) %>% mutate(Model = s3_levels[1]),
  as.data.frame(Variogram(m_spatial, form = ~ x_m + y_m, resType = "normalized",
                          data = df_vg, maxDist = max_d)) %>% mutate(Model = s3_levels[2])
) %>%
  mutate(Model = factor(Model, levels = s3_levels))

s3_theme <- theme_classic(base_family = "sans", base_size = 9) +
  theme(panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.line        = element_blank(),
        axis.text        = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text       = element_text(hjust = 0, size = 9),
        plot.title       = element_text(face = "bold", size = 9.5))

p_vg <- ggplot(vg, aes(x = dist, y = variog)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_point(colour = "grey30", alpha = 0.75, size = 1.3) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#B63679", linewidth = 0.7) +
  expand_limits(y = 0) +
  facet_wrap(~ Model, ncol = 2) +
  labs(title = "a) Semivariograms of normalised residuals",
       x = "Distance between colonies (m)", y = "Semivariance") +
  s3_theme

res_s3 <- bind_rows(
  df_vg %>% mutate(Resid = as.numeric(residuals(m_base, type = "normalized")),
                   Fitted = as.numeric(fitted(m_base)), Model = s3_levels[1]),
  df_vg %>% mutate(Resid = as.numeric(residuals(m_spatial, type = "normalized")),
                   Fitted = as.numeric(fitted(m_spatial)), Model = s3_levels[2])
) %>%
  mutate(Model = factor(Model, levels = s3_levels), Resid_capped = pmax(pmin(Resid, 3), -3)) %>%
  arrange(abs(Resid_capped))   # strongest departures drawn last, on top

p_res_map <- ggplot(res_s3, aes(x = x_m, y = y_m, colour = Resid_capped)) +
  geom_point(size = 0.7, alpha = 0.9) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0,
                         limits = c(-3, 3), breaks = -3:3,
                         labels = sub("^-", "\u2212", as.character(-3:3)),
                         name = "Normalised\nresidual") +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  scale_y_continuous(breaks = seq(0, 100, 2)) +
  coord_equal() +
  facet_wrap(~ Model, ncol = 1) +
  labs(title = "c) Normalised residuals across the wall",
       x = "Distance along wall (m)", y = "Height on wall (m)") +
  s3_theme +
  theme(legend.position = "right", legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8))

p_res_fit <- ggplot(res_s3, aes(x = Fitted, y = Resid)) +
  geom_point(aes(colour = fSpecies), alpha = 0.5, size = 0.9, stroke = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", colour = "red", linewidth = 0.4) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "black",
              fill = "grey70", linewidth = 0.6) +
  scale_colour_manual(values = SPECIES_COLOURS, labels = function(x) x, name = NULL) +
  scale_y_continuous(labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE))) +
  facet_wrap(~ Model, ncol = 2) +
  labs(title = "b) Residuals vs fitted values",
       x = "Fitted relative growth rate (RGR)", y = "Normalised residual") +
  s3_theme +
  theme(legend.position = "bottom", legend.text = element_text(size = 8, face = "italic"),
        legend.margin = margin(t = -4, b = 0)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

s3_h <- 125 + 2 * (174 * diff(range(df_vg$y_m)) / max(diff(range(df_vg$x_m)), 1) + 14)
save_figure(p_vg / p_res_fit / p_res_map + plot_layout(heights = c(0.8, 0.9, 1)),
            "growth_residual_semivariograms.png", width = 174, height = min(max(s3_h, 190), 240))

cat("\nNumbers for ESM Methods S2:\n")
cat(sprintf("  Moran's I, all colonies (n = %d): I = %.3f, p = %.3f; with the spatial structure I = %.3f, p = %.3f\n",
            mor_base$n, mor_base$Moran_I, mor_base$p_value, mor_sp$Moran_I, mor_sp$p_value))
cat(sprintf("  Delta AIC (spatial - model used) = %+.1f\n", d_aic))
cat(sprintf("  Nugget = %.2f (largest correlation between neighbouring colonies = %.2f); range = %.2f m\n",
            cs[["nugget"]], 1 - cs[["nugget"]], cs[["range"]] * PX_TO_M))
cat(sprintf("  Largest coefficient shift = %.2f SE; largest A_thresh change = %.1f%%\n",
            max(abs(shift_tab$Shift_in_SEs)),
            max(abs(comparison$Percent_change[comparison$Quantity == "A_thresh (cm2)"]))))
cat(sprintf("  Semivariograms run from %.2f to %.1f m\n", min(vg$dist), max(vg$dist)))
