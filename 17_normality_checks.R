# ==============================================================================
# 17_normality_checks.R
# The distributions behind the choice of non-parametric tests. For each
# species: Shapiro-Wilk test, skewness and a normal Q-Q plot of planar area
# (2015 and 2022), AGR and RGR of persistent colonies, and colonisation year.
# Also the dotplots and boxplots of the growth model's variables (2015 area
# and RGR; ESM Methods S2). Shapiro-Wilk flags even slight departures in large
# samples, so read the tests together with the Q-Q plots and skewness.
#
# Reads:  processed/colonies_with_mdc.xlsx (script 02), results/colony_ages.csv (script 08)
# Writes: results/normality_shapiro_wilk.csv
#         figures/normality_qq.png
#         figures/growth_variables_dotplots.png
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
})
options(width = 160)

short_names <- c("Madrepora oculata" = "M. oculata", "Desmophyllum pertusum" = "D. pertusum",
                 "Primnoa msp." = "Primnoa msp.")
p_text <- function(p) ifelse(is.na(p), "p = NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))

# Each variable uses the colonies of the test it informs: areas, every colony
# present that year; AGR and RGR, persistent colonies; colonisation year,
# colonies present in 2022
df <- read_excel(in_data("processed", "colonies_with_mdc.xlsx")) %>%
  mutate(Species = standard_species(Species)) %>%
  filter(Species %in% SPECIES) %>%
  mutate(Species = factor(Species, levels = SPECIES))
mc <- read_csv(in_data("results", "colony_ages.csv"), show_col_types = FALSE) %>%
  filter(Species %in% SPECIES, !is.na(Colonization_Year_Est)) %>%
  mutate(Species = factor(Species, levels = SPECIES))

var_levels <- c("Planar area 2015 (cm2)", "Planar area 2022 (cm2)", "AGR (cm2 per yr)",
                "RGR (per yr)", "Colonisation year")
var_short  <- setNames(c("Area 2015", "Area 2022", "AGR", "RGR", "Colonisation\nyear"), var_levels)

long <- bind_rows(
  df %>% filter(!is.na(Area_2015), Area_2015 > 0) %>%
    transmute(Species, Variable = var_levels[1], Value = Area_2015),
  df %>% filter(!is.na(Area_2022), Area_2022 > 0) %>%
    transmute(Species, Variable = var_levels[2], Value = Area_2022),
  df %>% filter(state_transition == "persistence", !is.na(Annual_Area_Change)) %>%
    transmute(Species, Variable = var_levels[3], Value = Annual_Area_Change),
  df %>% filter(state_transition == "persistence", !is.na(RGR_Planar), !is.na(A1)) %>%
    transmute(Species, Variable = var_levels[4], Value = RGR_Planar),
  mc %>% transmute(Species, Variable = var_levels[5], Value = as.numeric(Colonization_Year_Est))
) %>%
  filter(is.finite(Value)) %>%
  mutate(Variable = factor(Variable, levels = var_levels))

# ------------------------------------------------------------------------------
# Shapiro-Wilk tests and skewness
# ------------------------------------------------------------------------------
shapiro_wp <- function(x) {   # shapiro.test() needs 3-5000 values, not all identical
  if (length(x) < 3 || length(x) > 5000 || length(unique(x)) < 3) return(c(W = NA_real_, p = NA_real_))
  t <- shapiro.test(x)
  c(W = unname(t$statistic), p = t$p.value)
}
skewness <- function(x) {     # Fisher-Pearson coefficient, g1
  m2 <- mean((x - mean(x))^2)
  if (m2 == 0) return(NA_real_)
  mean((x - mean(x))^3) / m2^1.5
}

normality <- long %>%
  group_by(Variable, Species) %>%
  summarise(n        = n(),
            Median   = median(Value),
            IQR      = IQR(Value),
            Skewness = skewness(Value),
            W        = shapiro_wp(Value)[["W"]],
            p_value  = shapiro_wp(Value)[["p"]],
            .groups  = "drop") %>%
  mutate(Normal_at_0.05 = ifelse(is.na(p_value), NA, p_value >= 0.05))

cat("Shapiro-Wilk tests by species:\n")
print(as.data.frame(normality %>%
  mutate(across(c(Median, IQR), ~ formatC(.x, digits = 4, format = "g")),
         Skewness = round(Skewness, 2), W = round(W, 4),
         p_value = formatC(p_value, digits = 3, format = "g"))), row.names = FALSE)
cat(sprintf("\nDeparts from normality (p < 0.05): %d of %d variable-species combinations.\n",
            sum(normality$p_value < 0.05, na.rm = TRUE), sum(!is.na(normality$p_value))))

write.csv(normality %>%
            mutate(across(c(Median, IQR), ~ signif(.x, 5)), Skewness = round(Skewness, 3),
                   W = round(W, 4), p_value = signif(p_value, 3)),
          in_data("results", "normality_shapiro_wilk.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Normal Q-Q plots. Values are standardised within each panel (z-scores):
# the shapes are unchanged, but the species in a row share one axis.
# ------------------------------------------------------------------------------
long_z <- long %>%
  group_by(Variable, Species) %>%
  mutate(z = (Value - mean(Value)) / sd(Value)) %>%
  ungroup() %>%
  mutate(Row = factor(var_short[as.character(Variable)], levels = var_short),
         Col = factor(short_names[as.character(Species)], levels = short_names))
qq_labels <- normality %>%
  mutate(Row   = factor(var_short[as.character(Variable)], levels = var_short),
         Col   = factor(short_names[as.character(Species)], levels = short_names),
         label = sprintf("n = %d, W = %.3f\n%s", n, W, p_text(p_value)))

p_qq <- ggplot(long_z, aes(sample = z, colour = Species)) +
  stat_qq(size = 0.7, alpha = 0.6) +
  stat_qq_line(colour = "black", linetype = "dashed", linewidth = 0.4) +
  geom_text(data = qq_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.2, size = 2.3, lineheight = 0.9) +
  facet_grid(Row ~ Col, scales = "free_y") +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  scale_colour_manual(values = SPECIES_COLOURS, guide = "none") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles (standardised)") +
  theme_classic(base_family = "sans", base_size = 9) +
  theme(panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        strip.background = element_blank(),
        strip.text.x     = element_text(face = "italic", size = 9),
        strip.text.y     = element_text(size = 8.5),
        axis.text        = element_text(colour = "black", size = 7.5))
save_figure(p_qq, "normality_qq.png", width = 174, height = 220)

# ------------------------------------------------------------------------------
# Dotplots (each value against its row in the data: isolated points are
# outliers) and boxplots of the growth model's variables
# ------------------------------------------------------------------------------
df_model <- df %>%
  filter(state_transition == "persistence", !is.na(RGR_Planar), !is.na(A1)) %>%
  mutate(Col = factor(short_names[as.character(Species)], levels = short_names)) %>%
  group_by(Species) %>%
  mutate(Order = row_number()) %>%
  ungroup()

explore_theme <- theme_classic(base_family = "sans", base_size = 9) +
  theme(panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        strip.background = element_blank(),
        strip.text.y     = element_text(face = "italic", angle = 0, hjust = 0, size = 8.5),
        axis.text        = element_text(colour = "black", size = 7.5),
        legend.position  = "none",
        plot.tag         = element_text(face = "bold", size = 10))

dotplot <- function(var, xlab, tag) {
  ggplot(df_model, aes(x = .data[[var]], y = Order, colour = Species)) +
    geom_point(size = 0.8, alpha = 0.7) +
    facet_grid(Col ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = SPECIES_COLOURS) +
    labs(x = xlab, y = "Order of the data", tag = tag) +
    explore_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
boxplot_sp <- function(var, ylab, tag) {
  ggplot(df_model, aes(x = Col, y = .data[[var]], fill = Species)) +
    geom_boxplot(width = 0.6, linewidth = 0.3, outlier.size = 0.8, outlier.alpha = 0.7) +
    scale_fill_manual(values = SPECIES_COLOURS) +
    labs(x = NULL, y = ylab, tag = tag) +
    explore_theme +
    theme(axis.text.x = element_text(face = "italic"))
}

lab_area <- expression("2015 planar area" ~ (cm^2))
lab_rgr  <- expression("RGR" ~ (yr^-1))
p_explore <- (dotplot("A1", lab_area, "a") | boxplot_sp("A1", lab_area, "b")) /
  (dotplot("RGR_Planar", lab_rgr, "c") | boxplot_sp("RGR_Planar", lab_rgr, "d")) +
  plot_layout(widths = c(1.5, 1))
save_figure(p_explore, "growth_variables_dotplots.png", width = 174, height = 150)
