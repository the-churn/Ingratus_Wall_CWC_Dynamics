# ==============================================================================
# 16_variance_structure_sensitivity.R
# Does the choice of variance structure change the results? (ESM Methods S2)
# The best variance structures lie within about 2 AIC of each other (Table S3).
# The model used and the next three best-supported structures (fitted in
# script 05) are compared on what the paper reports:
#   - species parameters (RGR-infinity, RGR0, rc, A_thresh) with 95% CIs
#   - Wald tests between species
#   - colony ages in 2015 and 2022 (fitted curve)
#   - age at A_thresh (median and 95% interval over 2,000 parameter draws)
#
# Reads:  models/growth_model.rds, models/growth_model_candidates.rds (script 05),
#         models/age_simulation_inputs.rds, results/colony_ages.csv (script 08)
# Writes: tables/variance_structure_sensitivity.xlsx
#         figures/variance_structure_sensitivity.png
# ==============================================================================

source("setup.R")
source("age_model.R")     # age_at_size()
suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
  library(nlme)
  library(patchwork)
})

FONT_FAMILY    <- "sans"
N_SIMS         <- 2000    # parameter draws for the age at A_thresh
N_ALTERNATIVES <- 3       # next best-supported structures to compare with
set.seed(42)

m_saved <- readRDS(in_data("models", "growth_model.rds"))
cand    <- readRDS(in_data("models", "growth_model_candidates.rds"))
A0      <- readRDS(in_data("models", "age_simulation_inputs.rds"))$A0_species
current_name <- paste0(cand$selected, " (current)")

# Colony sizes (after the recruit audit) for the age comparison
df_age <- read.csv(in_data("results", "colony_ages.csv"), check.names = FALSE) %>%
  mutate(A1 = as.numeric(A1), A2 = as.numeric(A2)) %>%
  filter(Species %in% SPECIES) %>%
  select(Species, A1, A2)

# ------------------------------------------------------------------------------
# The model used and the next best-supported structures
# ------------------------------------------------------------------------------
alt_names <- cand$table %>%
  filter(!is.na(AIC), Model != cand$selected) %>%
  arrange(AIC) %>%
  slice_head(n = N_ALTERNATIVES) %>%
  pull(Model)
models <- c(setNames(list(m_saved), current_name), cand$fits[alt_names])
struct_levels <- names(models)

aic_tab <- data.frame(
  Structure = struct_levels,
  df        = vapply(models, function(m) attr(logLik(m), "df"), numeric(1)),
  AIC       = vapply(models, AIC, numeric(1)),
  row.names = NULL
) %>%
  mutate(dAIC = AIC - min(AIC)) %>%
  arrange(AIC)

# ------------------------------------------------------------------------------
# Species parameters and Wald tests between species
# ------------------------------------------------------------------------------
# Weights that pick one species' value of a parameter from the coefficients
L_vec <- function(b, par, sp) {
  L <- setNames(numeric(length(b)), names(b))
  L[paste0(par, ".(Intercept)")] <- 1
  if (sp != SPECIES[1]) L[paste0(par, ".fSpecies", sp)] <- 1
  L
}
lin_comb <- function(m, L) c(est = sum(L * coef(m)), se = sqrt(drop(t(L) %*% vcov(m) %*% L)))

param_tab <- bind_rows(lapply(struct_levels, function(st) {
  m <- models[[st]]
  bind_rows(lapply(SPECIES, function(sp) {
    a  <- lin_comb(m, L_vec(coef(m), "Asym", sp))
    r0 <- lin_comb(m, L_vec(coef(m), "R0",   sp))
    lr <- lin_comb(m, L_vec(coef(m), "lrc",  sp))
    data.frame(
      Structure = st, Species = sp,
      Parameter = c("RGR_inf", "RGR_0", "rc", "A_thresh (cm2)"),
      Estimate  = c(a["est"], r0["est"], exp(lr["est"]), log(20) / exp(lr["est"])),
      Lower95   = c(a["est"] - 1.96 * a["se"], r0["est"] - 1.96 * r0["se"],
                    exp(lr["est"] - 1.96 * lr["se"]), log(20) / exp(lr["est"] + 1.96 * lr["se"])),
      Upper95   = c(a["est"] + 1.96 * a["se"], r0["est"] + 1.96 * r0["se"],
                    exp(lr["est"] + 1.96 * lr["se"]), log(20) / exp(lr["est"] - 1.96 * lr["se"])),
      row.names = NULL
    )
  }))
}))

pairs <- list(SPECIES[c(1, 2)], SPECIES[c(1, 3)], SPECIES[c(2, 3)])
wald_tab <- bind_rows(lapply(struct_levels, function(st) {
  m <- models[[st]]
  res_df <- length(resid(m)) - length(coef(m))
  bind_rows(lapply(c("Asym", "R0", "lrc"), function(par) {
    bind_rows(lapply(pairs, function(pr) {
      lc <- lin_comb(m, L_vec(coef(m), par, pr[2]) - L_vec(coef(m), par, pr[1]))
      t  <- lc["est"] / lc["se"]
      data.frame(Structure = st, Parameter = par, Comparison = paste(pr[2], "vs", pr[1]),
                 Difference = lc["est"], SE = lc["se"], t = t, df = res_df,
                 p = 2 * pt(-abs(t), res_df), Significant = 2 * pt(-abs(t), res_df) < 0.05,
                 row.names = NULL)
    }))
  }))
}))

# ------------------------------------------------------------------------------
# Colony ages from the fitted curve of each structure
# ------------------------------------------------------------------------------
pars_for <- function(b) {
  setNames(lapply(SPECIES, function(sp) {
    g <- function(par) sum(L_vec(b, par, sp) * b)
    list(asym = g("Asym"), r0 = g("R0"), lrc = g("lrc"))
  }), SPECIES)
}
age_or_na <- function(size, sp, P) {
  if (is.na(size) || size <= 0) NA_real_ else age_at_size(size, P[[sp]], A0[[sp]])
}

ages <- bind_rows(lapply(struct_levels, function(st) {
  P <- pars_for(coef(models[[st]]))
  df_age %>%
    mutate(Colony    = row_number(),
           Structure = st,
           Age_2015  = mapply(age_or_na, A1, Species, MoreArgs = list(P = P)),
           Age_2022  = mapply(age_or_na, A2, Species, MoreArgs = list(P = P)))
}))
ages_long <- ages %>%
  pivot_longer(c(Age_2015, Age_2022), names_to = "Year", values_to = "Age") %>%
  filter(!is.na(Age))

age_summary <- ages_long %>%
  group_by(Structure, Species, Year) %>%
  summarise(n = n(), Median = median(Age), P90 = quantile(Age, 0.9), Max = max(Age), .groups = "drop")

age_diff <- ages_long %>%
  filter(Structure != current_name) %>%
  inner_join(ages_long %>% filter(Structure == current_name) %>% select(Colony, Year, Age_current = Age),
             by = c("Colony", "Year")) %>%
  mutate(Diff = Age - Age_current,
         Rel_pct = ifelse(Age_current > 0, 100 * Diff / Age_current, NA_real_)) %>%
  group_by(Structure, Species, Year) %>%
  summarise(n = n(),
            Median_abs_diff_yrs = median(abs(Diff)),
            Max_abs_diff_yrs    = max(abs(Diff)),
            Median_rel_diff_pct = median(Rel_pct, na.rm = TRUE),
            Correlation         = cor(Age, Age_current),
            .groups = "drop")

# ------------------------------------------------------------------------------
# Age at A_thresh over parameter draws
# ------------------------------------------------------------------------------
athresh_age <- bind_rows(lapply(struct_levels, function(st) {
  m <- models[[st]]
  draws <- MASS::mvrnorm(N_SIMS, mu = coef(m), Sigma = vcov(m))
  bind_rows(lapply(SPECIES, function(sp) {
    a <- vapply(seq_len(N_SIMS), function(k) {
      P <- pars_for(draws[k, ])[[sp]]
      age_at_size(log(20) / exp(P$lrc), P, A0[[sp]])
    }, numeric(1))
    data.frame(Structure = st, Species = sp,
               Median = median(a, na.rm = TRUE),
               Lower95 = quantile(a, 0.025, na.rm = TRUE),
               Upper95 = quantile(a, 0.975, na.rm = TRUE),
               Failed_draws = sum(is.na(a)), row.names = NULL)
  }))
}))

# ------------------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------------------
cur_par <- param_tab %>%
  filter(Structure == current_name) %>%
  select(Species, Parameter, Est_cur = Estimate, Lo_cur = Lower95, Hi_cur = Upper95)
par_shift <- param_tab %>%
  filter(Structure != current_name) %>%
  inner_join(cur_par, by = c("Species", "Parameter")) %>%
  mutate(Inside_current_CI = Estimate >= Lo_cur & Estimate <= Hi_cur,
         Change_pct = 100 * (Estimate - Est_cur) / abs(Est_cur))

sig_cur <- wald_tab %>%
  filter(Structure == current_name) %>%
  select(Parameter, Comparison, Sig_cur = Significant)
sig_flips <- wald_tab %>%
  filter(Structure != current_name) %>%
  inner_join(sig_cur, by = c("Parameter", "Comparison")) %>%
  filter(Significant != Sig_cur)

cat("Variance structures compared:\n")
print(aic_tab)
cat("\nSpecies parameters (estimate, 95% CI):\n")
print(param_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))
cat("\nAge at A_thresh (yr):\n")
print(athresh_age %>% mutate(across(where(is.numeric), ~ round(.x, 1))))
cat("\nColony ages vs the model used:\n")
print(age_diff %>% mutate(across(where(is.numeric), ~ round(.x, 2))), n = Inf)

cat("\nAlternative estimates inside the 95% CI of the model used:",
    sum(par_shift$Inside_current_CI), "of", nrow(par_shift), "\n")
cat("Largest change in any species parameter:", round(max(abs(par_shift$Change_pct)), 1), "%\n")
cat("Wald tests that change significance (p < 0.05):", nrow(sig_flips), "\n")
if (nrow(sig_flips) > 0) print(sig_flips %>% select(Structure, Parameter, Comparison, p, Sig_cur))
cat("Largest median change in colony age:", round(max(age_diff$Median_abs_diff_yrs), 2), "years (",
    round(max(abs(age_diff$Median_rel_diff_pct), na.rm = TRUE), 1), "% )\n")

# ------------------------------------------------------------------------------
# Figure: a) species parameters, b) 2022 ages, model used vs alternatives
# ------------------------------------------------------------------------------
theme_chk <- theme_classic(base_family = FONT_FAMILY) +
  theme(plot.title      = element_text(face = "bold", size = 11),
        plot.subtitle   = element_text(size = 8.5, color = "grey30"),
        axis.text       = element_text(size = 8.5, color = "black"),
        axis.title      = element_text(size = 9.5, face = "bold"),
        strip.text      = element_text(face = "bold", size = 8.5),
        panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "bottom", legend.title = element_blank())

p_par <- param_tab %>%
  filter(Parameter != "rc") %>%
  mutate(Structure = factor(Structure, levels = struct_levels),
         Species   = factor(Species, levels = SPECIES)) %>%
  ggplot(aes(x = Structure, y = Estimate, ymin = Lower95, ymax = Upper95, color = Structure)) +
  geom_pointrange(size = 0.3) +
  facet_wrap(Parameter ~ Species, scales = "free_y", ncol = 3) +
  labs(title = "a) Species parameters (95% CI) under each variance structure", x = NULL, y = NULL) +
  theme_chk +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_age <- ages_long %>%
  filter(Year == "Age_2022") %>%
  select(Colony, Species, Structure, Age) %>%
  pivot_wider(names_from = Structure, values_from = Age) %>%
  pivot_longer(-c(Colony, Species, all_of(current_name)), names_to = "Alternative", values_to = "Age_alt") %>%
  ggplot(aes(x = .data[[current_name]], y = Age_alt, color = Species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5, size = 0.9) +
  facet_wrap(~ Alternative, nrow = 1) +
  scale_color_manual(values = SPECIES_COLOURS) +
  labs(title = "b) Colony ages in 2022: current model vs alternatives",
       subtitle = "Points on the dashed line = identical age",
       x = paste0("Age with ", cand$selected, " (years)"), y = "Age with alternative (years)") +
  theme_chk +
  theme(legend.text = element_text(face = "italic"))

save_figure(p_par / p_age + plot_layout(heights = c(1.3, 1)),
            "variance_structure_sensitivity.png", width = 260, height = 240, dpi = 400)

write_xlsx(list(AIC = aic_tab, Species_Parameters = param_tab, Parameter_Shifts = par_shift,
                Wald_Tests = wald_tab, Age_Summary = age_summary, Age_Differences = age_diff,
                Age_at_Athresh = athresh_age),
           in_data("tables", "variance_structure_sensitivity.xlsx"))
