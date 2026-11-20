# ==============================================================================
# 05_growth_model.R
# Growth model (ESM Methods S2, Table S3, Fig. S2). Relative growth rate (RGR)
# of persistent colonies as an asymptotic function of 2015 size, with
# species-specific parameters, fitted by generalised nonlinear least squares:
#
#   RGR = Asym + (R0 - Asym) * exp(-exp(lrc) * A1)        (nlme::SSasymp)
#
# The residual variance structure is chosen following Zuur et al. (2009):
# candidate structures are fitted to this full model and compared by AIC; the
# simplest structure within 2 AIC of the best is used.
#
# Reads:  processed/colonies_with_mdc.xlsx      (script 02)
# Writes: models/growth_model.rds               the selected fit, used by later scripts
#         models/growth_model_candidates.rds    all candidate fits (script 16)
#         tables/growth_model.xlsx              Table S3, coefficients, residual checks
#         figures/growth_model_diagnostics.png  (Fig. S2)
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
  library(nlme)
  library(patchwork)
})

FONT_FAMILY <- "sans"

# Persistent colonies with a growth rate. NOTE: the saved models refer to the
# data as `df_growth` and the starting values as `start_vec`; scripts that
# refit them with update() recreate both under these names.
df_growth <- read_growth_data()
cat("Persistent colonies in the growth model:", nrow(df_growth), "\n")

# ------------------------------------------------------------------------------
# Starting values from separate fits per species, on the model's contrast
# scale (M. oculata = intercept; other species = difference from it)
# ------------------------------------------------------------------------------
per_species <- coef(nlsList(RGR_Planar ~ SSasymp(A1, Asym, R0, lrc) | fSpecies,
                            data = df_growth, na.action = na.omit))[SPECIES, ]
to_contrast <- function(v) c(v[1], v[-1] - v[1])
start_vec <- c(to_contrast(per_species$Asym), to_contrast(per_species$R0),
               to_contrast(per_species$lrc))

# ------------------------------------------------------------------------------
# Candidate variance structures. Each is fitted from one or more starting
# values for its variance parameters and the best fit is kept (varConstPower
# can otherwise stop where its constant shrinks to zero). The structures are
# unevaluated expressions so each fit stores the literal structure in its call.
# ------------------------------------------------------------------------------
cp_starts <- function(f) {
  grid <- expand.grid(const = c(0.01, 0.1, 1), power = c(-0.25, -0.5))
  lapply(seq_len(nrow(grid)), function(i)
    bquote(varConstPower(const = .(grid$const[i]), power = .(grid$power[i]), form = .(f))))
}

variance_candidates <- list(
  "Baseline (Homoscedastic)" = list(NULL),
  "VarIdent"             = list(quote(varIdent(form = ~ 1 | fSpecies))),
  "VarPower_Global"      = list(quote(varPower(form = ~ A1)),
                                quote(varPower(value = -0.5, form = ~ A1))),
  "VarPower_Sp"          = list(quote(varPower(form = ~ A1 | fSpecies)),
                                quote(varPower(value = -0.5, form = ~ A1 | fSpecies))),
  "VarExp_Global"        = list(quote(varExp(form = ~ A1)),
                                quote(varExp(value = -0.005, form = ~ A1))),
  "VarExp_Sp"            = list(quote(varExp(form = ~ A1 | fSpecies)),
                                quote(varExp(value = -0.005, form = ~ A1 | fSpecies))),
  "VarConstPower_Global" = c(list(quote(varConstPower(form = ~ A1))), cp_starts(~ A1)),
  "VarConstPower_Sp"     = c(list(quote(varConstPower(form = ~ A1 | fSpecies))),
                             cp_starts(~ A1 | fSpecies))
)

fit_growth_gnls <- function(weights_expr) {
  eval(bquote(gnls(
    RGR_Planar ~ SSasymp(A1, Asym, R0, lrc),
    params    = list(Asym ~ fSpecies, R0 ~ fSpecies, lrc ~ fSpecies),
    data      = df_growth,
    weights   = .(weights_expr),
    start     = start_vec,
    control   = gnlsControl(maxIter = 500, minScale = 1e-7, tolerance = 1e-4),
    na.action = na.omit
  )))
}

# Best fit (highest likelihood) over the starting values of one structure
fit_best_start <- function(starts) {
  fits <- lapply(starts, function(w) tryCatch(fit_growth_gnls(w), error = function(e) NULL))
  ok   <- !vapply(fits, is.null, logical(1))
  if (!any(ok)) return(list(fit = NULL, n_ok = 0L, n_tried = length(starts)))
  ll   <- vapply(fits[ok], function(m) as.numeric(logLik(m)), numeric(1))
  list(fit = fits[ok][[which.max(ll)]], n_ok = sum(ok), n_tried = length(starts))
}

variance_runs <- lapply(variance_candidates, fit_best_start)
variance_fits <- lapply(variance_runs, `[[`, "fit")

# ------------------------------------------------------------------------------
# AIC comparison (Table S3) and selection
# ------------------------------------------------------------------------------
from_fit <- function(f) sapply(variance_fits, function(m) if (is.null(m)) NA else f(m))
model_comparison <- data.frame(
  Model            = names(variance_fits),
  df               = from_fit(function(m) attr(logLik(m), "df")),
  logLik           = from_fit(function(m) as.numeric(logLik(m))),
  AIC              = from_fit(AIC),
  Starts_converged = sapply(variance_runs, function(r) paste0(r$n_ok, "/", r$n_tried)),
  row.names        = NULL
) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  arrange(AIC)

# Simplest structure (fewest parameters) within 2 AIC of the best
selected_structure <- model_comparison %>%
  filter(!is.na(dAIC), dAIC < 2) %>%
  arrange(df, AIC) %>%
  slice(1) %>%
  pull(Model)
model_comparison <- model_comparison %>% mutate(Selected = Model == selected_structure)

cat("\nVariance structures compared by AIC (Table S3):\n")
print(model_comparison)
cat("\nSelected:", selected_structure, "\n")

m_growth_gnls <- variance_fits[[selected_structure]]
m_0           <- variance_fits[["Baseline (Homoscedastic)"]]

cat("\nGrowth model (", selected_structure, "):\n", sep = "")
print(summary(m_growth_gnls))

var_params <- if (is.null(m_growth_gnls$modelStruct$varStruct)) {
  data.frame(Parameter = character(0), Value = numeric(0))
} else {
  vp <- coef(m_growth_gnls$modelStruct$varStruct, unconstrained = FALSE)
  data.frame(Parameter = names(vp), Value = unname(vp))
}
cat("\nVariance-function parameters:\n")
print(var_params)

# Species-level estimates; they should be close to the separate per-species fits
cf <- coef(m_growth_gnls)
species_value <- function(par, i) {
  ref <- cf[[paste0(par, ".(Intercept)")]]
  if (i == 1) ref else ref + cf[[paste0(par, ".fSpecies", SPECIES[i])]]
}
cat("\nSpecies-level estimates, this model vs separate fits per species:\n")
print(data.frame(
  Species   = SPECIES,
  Asym_gnls = sapply(1:3, function(i) species_value("Asym", i)), Asym_nls = per_species$Asym,
  R0_gnls   = sapply(1:3, function(i) species_value("R0", i)),   R0_nls   = per_species$R0,
  lrc_gnls  = sapply(1:3, function(i) species_value("lrc", i)),  lrc_nls  = per_species$lrc
), digits = 4, row.names = FALSE)

saveRDS(m_growth_gnls, in_data("models", "growth_model.rds"))
saveRDS(list(fits = variance_fits, table = model_comparison, selected = selected_structure),
        in_data("models", "growth_model_candidates.rds"))

# ------------------------------------------------------------------------------
# Residual diagnostics (Fig. S2): equal variance vs the model used.
# Standardised (Pearson) residuals = residual / its modelled SD.
# ------------------------------------------------------------------------------
baseline_label <- "Equal variance assumed"
winning_label  <- switch(selected_structure,
  "Baseline (Homoscedastic)" = "Equal variance (model used)",
  "VarIdent"                 = "Variance differs by species (model used)",
  "Variance changes with colony size (model used)"
)

df_diag_m0 <- df_growth %>%
  mutate(Resid = residuals(m_0, type = "pearson"), Fitted = fitted(m_0), Model = baseline_label)
df_diag_win <- df_growth %>%
  mutate(Resid = residuals(m_growth_gnls, type = "pearson"), Fitted = fitted(m_growth_gnls),
         Model = winning_label)
df_diag_all <- bind_rows(df_diag_m0, df_diag_win) %>%
  mutate(Model = factor(Model, levels = c(baseline_label, winning_label)))

# Does the residual spread change with size? SD of the residuals in 8
# equal-count size classes (should be 1 throughout), and Spearman's rho
# between |residual| and log area (should be 0)
N_SIZE_BINS <- 8
spread_df <- df_diag_all %>%
  mutate(Size_bin = ntile(A1, N_SIZE_BINS)) %>%
  group_by(Model, Size_bin) %>%
  summarise(A1_median = median(A1), SD = sd(Resid), n = n(), .groups = "drop") %>%
  mutate(Lower95 = SD - 1.96 * SD / sqrt(2 * (n - 1)),
         Upper95 = SD + 1.96 * SD / sqrt(2 * (n - 1)))

spread_rho <- df_diag_all %>%
  group_by(Model) %>%
  summarise(
    Spearman_rho    = cor(abs(Resid), log(A1), method = "spearman"),
    p_value         = cor.test(abs(Resid), log(A1), method = "spearman", exact = FALSE)$p.value,
    SD_smallest_bin = sd(Resid[ntile(A1, N_SIZE_BINS) == 1]),
    SD_largest_bin  = sd(Resid[ntile(A1, N_SIZE_BINS) == N_SIZE_BINS]),
    .groups = "drop"
  )
cat("\nResidual spread vs colony size (|residual| ~ log area):\n")
print(spread_rho)

minus <- function(x) sub("^-", "\u2212", x)            # true minus sign
axis_minus <- function(b) minus(format(b, trim = TRUE))

species_lab <- c("Madrepora oculata"     = 'italic("Madrepora oculata")',
                 "Desmophyllum pertusum" = 'italic("Desmophyllum pertusum")',
                 "Primnoa msp."          = 'italic("Primnoa")~"msp."')

fig_theme <- theme_classic(base_family = FONT_FAMILY, base_size = 9) +
  theme(
    plot.title       = element_text(face = "bold", size = 9.5),
    axis.text        = element_text(colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.line        = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0, size = 9)
  )

# Shared layers of panels a and b: points, zero line, +-2 lines, loess mean
resid_layers <- list(
  geom_point(aes(colour = fSpecies), alpha = 0.5, size = 0.9, stroke = 0),
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5),
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", colour = "grey40", linewidth = 0.4),
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "black",
              fill = "grey70", linewidth = 0.6),
  scale_colour_manual(values = SPECIES_COLOURS, labels = parse(text = species_lab[SPECIES]),
                      name = NULL),
  scale_y_continuous(labels = axis_minus),
  facet_wrap(~ Model, ncol = 2),
  fig_theme
)

p_fitted <- ggplot(df_diag_all, aes(x = Fitted, y = Resid)) +
  resid_layers +
  labs(title = "a) Residuals vs fitted values",
       x = expression(paste("Fitted RGR (yr"^-1, ")")), y = "Standardised residual") +
  theme(legend.position = "none")

p_area <- ggplot(df_diag_all, aes(x = A1, y = Resid)) +
  resid_layers +
  scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1,000")) +
  labs(title = "b) Residuals vs colony size",
       x = expression(paste("2015 planar area, ", italic("A"), " (cm"^2, ", log scale)")),
       y = "Standardised residual") +
  theme(legend.position = "bottom", legend.text = element_text(size = 8.5),
        legend.margin = margin(t = -3, b = 0)) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))

short_label <- setNames(c("Equal variance assumed", "Model used"), c(baseline_label, winning_label))
spread_lab <- with(spread_rho, setNames(
  sprintf("%s (\u03C1 = %s, %s)", short_label[as.character(Model)],
          minus(sprintf("%.2f", Spearman_rho)),
          ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.2f", p_value))),
  as.character(Model)))
spread_plot_df <- spread_df %>%
  mutate(x_plot = A1_median * ifelse(Model == baseline_label, 0.96, 1.04))  # side by side

p_spread <- ggplot(spread_plot_df, aes(x = x_plot, y = SD, colour = Model, shape = Model)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, linewidth = 0.5) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.9, fill = "white") +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100, 300)) +
  scale_colour_manual(values = setNames(c("grey55", "black"), levels(df_diag_all$Model)),
                      labels = spread_lab, name = NULL) +
  scale_shape_manual(values = setNames(c(21, 16), levels(df_diag_all$Model)),
                     labels = spread_lab, name = NULL) +
  labs(title = "c) Spread of residuals by colony size",
       x = expression(paste("Median 2015 planar area of size class (cm"^2, ", log scale)")),
       y = "Residual SD") +
  fig_theme +
  theme(legend.position = "bottom", legend.text = element_text(size = 8.5),
        legend.margin = margin(t = -3, b = 0))

qq_df <- df_diag_win %>%
  mutate(Panel = factor(species_lab[as.character(fSpecies)], levels = species_lab[SPECIES]))
p_qq <- ggplot(qq_df, aes(sample = Resid, colour = fSpecies)) +
  stat_qq_line(colour = "black", linetype = "dashed", linewidth = 0.5) +
  stat_qq(alpha = 0.6, size = 0.9, stroke = 0) +
  scale_colour_manual(values = SPECIES_COLOURS, guide = "none") +
  scale_x_continuous(labels = axis_minus) +
  scale_y_continuous(labels = axis_minus) +
  facet_wrap(~ Panel, nrow = 1, scales = "free_y", labeller = label_parsed) +
  labs(title = "d) Normal Q\u2013Q plots of standardised residuals (model used)",
       x = "Theoretical quantiles", y = "Sample quantiles") +
  fig_theme

fig_s2 <- p_fitted / p_area / p_spread / p_qq + plot_layout(heights = c(1, 1, 0.7, 0.75))
save_figure(fig_s2, "growth_model_diagnostics.png", width = 174, height = 232)

write_xlsx(list(
  variance_selection      = model_comparison,
  coefficients            = as.data.frame(summary(m_growth_gnls)$tTable) %>%
                              rownames_to_column("Parameter"),
  variance_parameters     = var_params,
  residual_spread         = spread_rho,
  residual_spread_by_size = spread_df
), in_data("tables", "growth_model.xlsx"))
