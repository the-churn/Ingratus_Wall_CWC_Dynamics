# ==============================================================================
# 15_age_bias_sensitivity.R
# How far could a bias in the age model change the colonisation comparisons?
# (ESM Methods S3, Fig. S4). Script 14 (check 1) finds that the modelled time
# between the 2015 and 2022 sizes of colonies that grew is longer than the true
# 7 years (about 1.2x for D. pertusum and Primnoa, 2x for M. oculata): the
# model describes average net growth, including colonies that lost tissue,
# while the check can only use colonies that grew.
#
# Worst case: every colony grew all its life as fast as those that grew
# between the surveys, which shortens all ages of a species by the factor
# 7 / (median modelled interval). The colonisation-year comparisons of script
# 11 are repeated with model ages, halfway, and worst-case ages: median age,
# Kruskal-Wallis, and quantile-regression contrasts at tau = 0.05, 0.10 and
# 0.50 (bootstrapped SEs, 1,000 resamples), M. oculata against each species.
# Then: how much of the gap must be real bias for the conclusions to change?
# A share w of it shortens ages by 1 - w * (1 - factor) (w = 0: model ages,
# w = 0.5: halfway, w = 1: worst case).
#
# Reads:  results/colonisation_years.csv (script 10),
#         results/age_elapsed_time_check.csv (script 14)
# Writes: results/age_bias_colonisation.csv   model, halfway and worst case
#         results/age_bias_grid.csv           the same for w = 0 to 1
#         figures/age_bias_sensitivity.png    (Fig. S4)
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(quantreg)
})
options(width = 160)

SURVEY_YEAR <- 2022
TAUS        <- c(0.05, 0.10, 0.50)
N_BOOT      <- 1000
species_lv  <- c("Desmophyllum pertusum", "Madrepora oculata", "Primnoa msp.")

# Age-shortening factor per species: true interval / median modelled interval
elapsed <- read_csv(in_data("results", "age_elapsed_time_check.csv"), show_col_types = FALSE)
shorten <- setNames(SURVEY_INTERVAL / elapsed$Median_Delta_Age_Baseline, elapsed$Species)
cat("Median modelled interval (yr) and age-shortening factor:\n")
print(data.frame(Species = names(shorten), Predicted_interval = elapsed$Median_Delta_Age_Baseline,
                 Factor = round(unname(shorten), 3)), row.names = FALSE)

d <- read_csv(in_data("results", "colonisation_years.csv"), show_col_types = FALSE) %>%
  filter(Species %in% species_lv, !is.na(Age_2022_Median)) %>%
  mutate(Species  = factor(Species, levels = species_lv),
         Year_obs = SURVEY_YEAR - Age_2022_Median,
         Age_wc   = Age_2022_Median * shorten[as.character(Species)],
         Year_wc  = SURVEY_YEAR - Age_wc,
         Age_mid  = (Age_2022_Median + Age_wc) / 2,
         Year_mid = SURVEY_YEAR - Age_mid)

# ------------------------------------------------------------------------------
# Model ages, halfway and worst case
# ------------------------------------------------------------------------------
med <- d %>%
  group_by(Species) %>%
  summarise(n = n(), Median_age = median(Age_2022_Median), Median_age_halfway = median(Age_mid),
            Median_age_worst_case = median(Age_wc), .groups = "drop")
kw_obs <- kruskal.test(Year_obs ~ Species, data = d)$p.value
kw_mid <- kruskal.test(Year_mid ~ Species, data = d)$p.value
kw_wc  <- kruskal.test(Year_wc ~ Species, data = d)$p.value

set.seed(2026)
contrast <- function(yvar, tau, ref) {
  dd  <- d %>% mutate(Species = relevel(Species, ref = ref))
  fit <- rq(as.formula(paste(yvar, "~ Species")), data = dd, tau = tau, method = "fn")
  s   <- summary(fit, se = "boot", R = N_BOOT)$coefficients
  r   <- grep("Madrepora", rownames(s))
  data.frame(Ages = c(Year_obs = "Model", Year_mid = "Halfway", Year_wc = "Worst case")[[yvar]],
             tau = tau, Contrast = paste("M. oculata vs", ref),
             Difference_yr = s[r, 1], SE = s[r, 2], p_value = s[r, 4])
}
qr_tab <- bind_rows(lapply(c("Year_obs", "Year_mid", "Year_wc"), function(y)
  bind_rows(lapply(TAUS, function(t)
    bind_rows(lapply(c("Desmophyllum pertusum", "Primnoa msp."), function(r) contrast(y, t, r)))))))

cat("\nMedian age in 2022 (yr):\n")
print(as.data.frame(med), row.names = FALSE, digits = 3)
cat(sprintf("\nKruskal-Wallis on colonisation year: model ages p = %.2g; halfway p = %.2g; worst case p = %.2g\n",
            kw_obs, kw_mid, kw_wc))
cat("\nM. oculata minus other species, colonisation year (negative = earlier):\n")
print(qr_tab %>% mutate(across(c(Difference_yr, SE), ~ round(.x, 1)), p_value = signif(p_value, 2)),
      row.names = FALSE)

out <- bind_rows(
  med %>% transmute(Quantity = "Median age in 2022 (yr)", Group = as.character(Species),
                    Model = Median_age, Halfway = Median_age_halfway, Worst_case = Median_age_worst_case),
  data.frame(Quantity = "Kruskal-Wallis p", Group = "All species", Model = kw_obs, Halfway = kw_mid,
             Worst_case = kw_wc),
  qr_tab %>%
    mutate(Quantity = paste0("Quantile contrast tau = ", tau, " (yr; SE; p)"), Group = Contrast) %>%
    select(Quantity, Group, Ages, Difference_yr, SE, p_value) %>%
    pivot_wider(names_from = Ages, values_from = c(Difference_yr, SE, p_value))
)
write_csv(out, in_data("results", "age_bias_colonisation.csv"))

# ------------------------------------------------------------------------------
# How much of the gap would have to be bias? Ages x (1 - w * (1 - factor))
# ------------------------------------------------------------------------------
SHARES <- sort(unique(round(c(seq(0, 1, by = 0.1), 0.25, 0.75), 2)))
grid <- bind_rows(lapply(SHARES, function(w) {
  dd  <- d %>% mutate(Age  = Age_2022_Median * (1 - w * (1 - shorten[as.character(Species)])),
                      Year = SURVEY_YEAR - Age)
  med <- tapply(dd$Age, dd$Species, median)
  qr_est <- function(tau, ref) {
    cf <- coef(rq(Year ~ Species, data = dd %>% mutate(Species = relevel(Species, ref = ref)),
                  tau = tau, method = "fn"))
    unname(cf["SpeciesMadrepora oculata"])
  }
  mo <- dd$Year[dd$Species == "Madrepora oculata"]
  data.frame(
    Share_of_bias         = w,
    Median_age_Moc        = med[["Madrepora oculata"]],
    Median_age_Dpe        = med[["Desmophyllum pertusum"]],
    Median_age_Pri        = med[["Primnoa msp."]],
    KW_p                  = kruskal.test(Year ~ Species, data = dd)$p.value,
    MannWhitney_p_Moc_Dpe = wilcox.test(mo, dd$Year[dd$Species == "Desmophyllum pertusum"])$p.value,
    MannWhitney_p_Moc_Pri = wilcox.test(mo, dd$Year[dd$Species == "Primnoa msp."])$p.value,
    P05_Moc_minus_Dpe_yr  = qr_est(0.05, "Desmophyllum pertusum"),
    P05_Moc_minus_Pri_yr  = qr_est(0.05, "Primnoa msp."),
    P10_Moc_minus_Dpe_yr  = qr_est(0.10, "Desmophyllum pertusum"),
    P10_Moc_minus_Pri_yr  = qr_est(0.10, "Primnoa msp."),
    P50_Moc_minus_Dpe_yr  = qr_est(0.50, "Desmophyllum pertusum"),
    P50_Moc_minus_Pri_yr  = qr_est(0.50, "Primnoa msp.")
  )
}))
cat("\nHow much of the gap would have to be bias? (0 = model ages, 1 = worst case)\n")
print(grid %>% mutate(across(starts_with(c("Median", "P0", "P1", "P5")), ~ round(.x, 1)),
                      across(ends_with("_p") | contains("_p_"), ~ signif(.x, 2))), row.names = FALSE)
write_csv(grid, in_data("results", "age_bias_grid.csv"))

# ------------------------------------------------------------------------------
# Fig. S4: cumulative share of colonies by colonisation year, per species.
# Where a curve crosses a dashed line gives that species' 5th, 10th or 50th
# percentile; the leftmost curve colonised earliest.
# ------------------------------------------------------------------------------
species_lab <- c("Madrepora oculata"     = 'italic("M. oculata")',
                 "Desmophyllum pertusum" = 'italic("D. pertusum")',
                 "Primnoa msp."          = 'italic("Primnoa")~"msp."')
ages_lv <- c("a) Ages from the\ngrowth model", "b) Halfway between\na and c",
             "c) Ages shortened to match\nthe validation check")
ecdf_df <- bind_rows(
  d %>% transmute(Species, Year = Year_obs, Ages = ages_lv[1]),
  d %>% transmute(Species, Year = Year_mid, Ages = ages_lv[2]),
  d %>% transmute(Species, Year = Year_wc,  Ages = ages_lv[3])
) %>% mutate(Ages = factor(Ages, levels = ages_lv))
p_lab <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.2f", p))
kw_text <- data.frame(Ages = factor(ages_lv, levels = ages_lv),
                      label = paste("Kruskal\u2013Wallis", p_lab(c(kw_obs, kw_mid, kw_wc))))
# 5th and 10th labelled at the right edge, median at the left (no curves there)
ref_lines <- data.frame(y = c(0.05, 0.10, 0.50), label = c("5th", "10th", "Median"),
                        x = c(Inf, Inf, -Inf), hjust = c(1.1, 1.1, -0.1))

fig_s4 <- ggplot(ecdf_df, aes(x = Year, colour = Species)) +
  geom_hline(data = ref_lines, aes(yintercept = y), linetype = "dashed", colour = "grey60",
             linewidth = 0.35) +
  geom_text(data = ref_lines, aes(x = x, y = y, label = label, hjust = hjust), inherit.aes = FALSE,
            vjust = -0.35, size = 2.5, colour = "grey35") +
  geom_text(data = kw_text, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.6, size = 2.5, colour = "grey20") +
  stat_ecdf(geom = "step", linewidth = 0.7) +
  scale_x_continuous(breaks = c(1900, 1950, 2000)) +
  scale_colour_manual(values = SPECIES_COLOURS, breaks = names(species_lab),
                      labels = parse(text = species_lab), name = NULL) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = 0.02)) +
  facet_wrap(~ Ages, nrow = 1) +
  labs(x = "Colonisation year", y = "Cumulative proportion of colonies") +
  theme_classic(base_family = "sans", base_size = 9) +
  theme(panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.line        = element_blank(),
        axis.text        = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text       = element_text(hjust = 0, size = 9, face = "bold"),
        legend.position  = "bottom",
        legend.text      = element_text(size = 8.5),
        legend.margin    = margin(t = -3, b = 0),
        panel.spacing.x  = unit(5, "mm"))
save_figure(fig_s4, "age_bias_sensitivity.png", width = 174, height = 92)
