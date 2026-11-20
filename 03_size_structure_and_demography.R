# ==============================================================================
# 03_size_structure_and_demography.R
# Fig. 2 and Table 2: colony size distributions in 2015 and 2022 (KS tests),
# annual growth of persistent colonies (Kruskal-Wallis and Dunn tests), and
# recruitment and mortality rates with Wilson 95% confidence intervals.
#
# Reads:  processed/colonies_with_mdc.xlsx         (script 02)
# Writes: figures/size_structure_and_demography.png (Fig. 2)
#         tables/size_structure_and_demography.xlsx (the tables printed below)
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(patchwork)
  library(ggsignif)
  library(dunn.test)
})

FONT_FAMILY <- "sans"
c_2015_fill <- "#E5E5E5"
c_2015_line <- "#404040"

colonies <- read_excel(in_data("processed", "colonies_with_mdc.xlsx")) %>%
  mutate(Species = standard_species(Species)) %>%
  filter(Species %in% SPECIES) %>%                 # drops unidentified recruits
  mutate(Species = factor(Species, levels = SPECIES))

# ------------------------------------------------------------------------------
# Panels a-c: size distributions, 2015 vs 2022
# ------------------------------------------------------------------------------
sizes <- colonies %>%
  arrange(Species) %>%
  select(TagLab.Genet.Id, Group = Species, Area_2015, Area_2022) %>%
  pivot_longer(c(Area_2015, Area_2022), names_to = "Year", values_to = "Area") %>%
  mutate(Year = str_remove(Year, "Area_")) %>%
  filter(Area > 0) %>%                             # absent (or hidden) in that survey
  mutate(Log_Area = log(Area))

size_stats <- sizes %>%
  group_by(Group) %>%
  summarise(
    n_15    = sum(Year == "2015"),
    n_22    = sum(Year == "2022"),
    med_15  = median(Area[Year == "2015"]),
    med_22  = median(Area[Year == "2022"]),
    ks_stat = unname(ks.test(Area[Year == "2015"], Area[Year == "2022"])$statistic),
    ks_p    = ks.test(Area[Year == "2015"], Area[Year == "2022"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_label = if_else(ks_p < 0.001, "p < 0.001", paste0("p = ", round(ks_p, 3))))

# ------------------------------------------------------------------------------
# Panels d-f: growth, recruitment and mortality
# ------------------------------------------------------------------------------
persistent <- colonies %>% filter(state_transition == "persistence")

# Colonies hidden in either survey are left out of both denominators
rates <- colonies %>%
  group_by(Species) %>%
  summarise(
    n_2015          = sum(state_transition %in% c("persistence", "mortality")),
    n_2022          = sum(state_transition %in% c("persistence", "recruitment")),
    n_persistence   = sum(state_transition == "persistence"),
    n_mortality     = sum(state_transition == "mortality"),
    n_recruitment   = sum(state_transition == "recruitment"),
    n_occluded      = sum(state_transition == "occluded"),
    n_occluded_2015 = sum(state_transition == "occluded_2015"),
    mean_mdc        = mean(MDC95_Annual[state_transition == "persistence"], na.rm = TRUE),
    .groups         = "drop"
  )

wilson <- function(x, n) Hmisc::binconf(x, n, method = "wilson") * 100   # % with 95% CI
mort_ci <- wilson(rates$n_mortality, rates$n_2015)
rec_ci  <- wilson(rates$n_recruitment, rates$n_2022)
rates <- rates %>%
  mutate(
    mortality_rate   = mort_ci[, "PointEst"],
    mort_ci_lower    = mort_ci[, "Lower"],
    mort_ci_upper    = mort_ci[, "Upper"],
    recruitment_rate = rec_ci[, "PointEst"],
    rec_ci_lower     = rec_ci[, "Lower"],
    rec_ci_upper     = rec_ci[, "Upper"],
    occlusion_rate   = 100 * n_occluded / (n_2015 + n_occluded)   # reported only
  )

n_labels <- function(n) setNames(paste0("(n = ", n, ")"), rates$Species)

# Growth differences between species: Kruskal-Wallis, then Dunn (BH-adjusted)
kw_test  <- kruskal.test(Annual_Area_Change ~ Species, data = persistent)
dunn_res <- dunn.test(x = persistent$Annual_Area_Change, g = persistent$Species,
                      method = "bh", kw = TRUE)
dunn_p <- setNames(dunn_res$P.adjusted, dunn_res$comparisons)

comparison_list <- list(c("Madrepora oculata", "Desmophyllum pertusum"),
                        c("Desmophyllum pertusum", "Primnoa msp."),
                        c("Madrepora oculata", "Primnoa msp."))
pair_p <- function(pair) {
  keys <- c(paste(pair[1], "-", pair[2]), paste(pair[2], "-", pair[1]))
  unname(dunn_p[keys[keys %in% names(dunn_p)]])
}
p_to_stars <- function(p) {
  as.character(cut(p, breaks = c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "ns")))
}
dunn_labels <- p_to_stars(sapply(comparison_list, pair_p))

# ------------------------------------------------------------------------------
# Themes
# ------------------------------------------------------------------------------
strip_labeller <- as_labeller(
  c("Madrepora oculata"     = "bolditalic('Madrepora oculata')",
    "Desmophyllum pertusum" = "bolditalic('Desmophyllum pertusum')",
    "Primnoa msp."          = "bolditalic('Primnoa')~bold('msp.')"),
  label_parsed
)

theme_horizontal_panel <- function() {
  theme_minimal(base_size = 10, base_family = FONT_FAMILY) +
    theme(
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.text.x      = element_text(color = "black", size = 8.5),
      axis.title.x     = element_text(face = "bold", size = 9.5, margin = margin(t = 8, b = 2)),
      axis.text.y      = element_text(color = "black", size = 8.5),
      axis.title.y     = element_text(face = "bold", size = 9.5, margin = margin(r = 6, b = 2)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      panel.spacing    = unit(8, "pt"),
      strip.background = element_rect(fill = "white", color = "grey80", linewidth = 0.4),
      strip.text       = element_text(face = "plain", size = 9, margin = margin(t = 4, b = 4)),
      legend.position  = "none",
      plot.tag         = element_text(face = "bold", size = 12),
      plot.margin      = margin(t = 6, r = 6, b = 8, l = 6)
    )
}

# Base theme for panels d-f and the legends
theme_set(theme_classic(base_size = 8, base_family = FONT_FAMILY))

publication_theme <- theme(
  plot.title         = element_blank(),
  axis.title.x       = element_blank(),
  axis.title.y       = element_text(margin = margin(r = 6, b = 2), face = "bold", color = "black", size = 8.5),
  axis.text.x        = element_text(color = "black", size = 7.5, hjust = 0.5, vjust = 0.5),
  axis.text.y        = element_text(color = "black", size = 7.5),
  axis.line          = element_line(linewidth = 0.3, color = "black"),
  axis.ticks         = element_line(linewidth = 0.3, color = "black"),
  panel.grid.major.y = element_line(color = "gray96", linetype = "solid"),
  legend.position    = "none",
  plot.tag           = element_text(face = "bold", size = 11, vjust = 1),
  plot.margin        = margin(t = 8, r = 6, b = 8, l = 6)
)

# ------------------------------------------------------------------------------
# Panels a-c
# ------------------------------------------------------------------------------
area_axis_title <- expression(bold(paste("Colony planar area (cm"^2, ")")))   # log axis

make_density_plot <- function(species, tag) {
  sub_data  <- sizes %>% filter(Group == species)
  sub_stats <- size_stats %>% filter(Group == species)
  col_22    <- SPECIES_COLOURS[[species]]
  stats_block <- paste(paste0("n: ", sub_stats$n_15, " \u2192 ", sub_stats$n_22),
                       paste0("KS: D = ", round(sub_stats$ks_stat, 2)),
                       sub_stats$p_label, sep = "\n")

  ggplot(sub_data, aes(x = Log_Area)) +
    geom_density(data = filter(sub_data, Year == "2015"),
                 fill = c_2015_fill, color = c_2015_line, linetype = "dashed", alpha = 0.85, linewidth = 0.5) +
    geom_density(data = filter(sub_data, Year == "2022"),
                 fill = col_22, color = col_22, linetype = "solid", alpha = 0.4, linewidth = 0.5) +
    geom_vline(xintercept = log(sub_stats$med_15), color = c_2015_line, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = log(sub_stats$med_22), color = col_22,      linetype = "solid",  linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf, label = stats_block,
             size = 2.3, hjust = -0.1, vjust = 1.15, lineheight = 0.9,
             color = "black", family = FONT_FAMILY) +
    facet_wrap(~Group, labeller = strip_labeller) +
    scale_x_continuous(breaks = log(c(1, 10, 100, 1000)), labels = scales::comma(c(1, 10, 100, 1000)),
                       limits = log(c(0.8, 2500))) +
    coord_cartesian(clip = "off") +
    labs(x = area_axis_title, y = "Probability density", tag = tag) +
    theme_horizontal_panel()
}

p_a <- make_density_plot("Madrepora oculata", "a") + theme(axis.title.x = element_blank())
p_b <- make_density_plot("Desmophyllum pertusum", "b") + theme(axis.title.y = element_blank())
p_c <- make_density_plot("Primnoa msp.", "c") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

# ------------------------------------------------------------------------------
# Panel d: annual growth. Filled points: change larger than the colony's own
# MDC95; open circles: smaller (drawn on top so they stay visible).
# ------------------------------------------------------------------------------
p_d <- ggplot(persistent, aes(x = Species, y = Annual_Area_Change)) +
  geom_violin(aes(fill = Species), alpha = 0.08, color = NA) +
  geom_point(data = ~ dplyr::filter(.x, detectable_change %in% TRUE), aes(color = Species),
             position = position_jitter(width = 0.18, height = 0, seed = 1), size = 0.6, alpha = 0.8) +
  geom_point(data = ~ dplyr::filter(.x, detectable_change %in% FALSE),
             position = position_jitter(width = 0.18, height = 0, seed = 2),
             shape = 21, fill = "white", color = "grey15", size = 0.9, stroke = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.95, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.3) +
  geom_signif(comparisons = comparison_list, annotations = dunn_labels,
              y_position = c(58, 66, 74), textsize = 3.0, vjust = 0.3,
              tip_length = 0.02, linewidth = 0.35, color = "black") +
  scale_fill_manual(values = SPECIES_COLOURS) +
  scale_color_manual(values = SPECIES_COLOURS) +
  scale_x_discrete(labels = n_labels(rates$n_persistence)) +
  scale_y_continuous(breaks = seq(0, 75, by = 25), limits = c(NA, 80),
                     expand = expansion(mult = c(0.08, 0))) +
  coord_cartesian(clip = "off") +
  labs(tag = "d", y = expression(bold(paste("AGR (", cm^2, " ", yr^-1, ")")))) +
  publication_theme

# Panels e and f: recruitment and mortality rates with Wilson 95% CIs
rate_plot <- function(rate, lower, upper, n, tag, y_title) {
  ggplot(rates, aes(x = Species, y = {{ rate }}, fill = Species)) +
    geom_col(width = 0.45, alpha = 0.85, color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = {{ lower }}, ymax = {{ upper }}), width = 0.10, color = "black", linewidth = 0.4) +
    geom_text(aes(y = {{ upper }}, label = sprintf("%.1f%%", {{ rate }})), vjust = -0.6, size = 2.4, fontface = "bold") +
    scale_fill_manual(values = SPECIES_COLOURS) +
    scale_x_discrete(labels = n_labels(n)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    coord_cartesian(clip = "off") +
    labs(tag = tag, y = y_title) +
    publication_theme
}
p_e <- rate_plot(recruitment_rate, rec_ci_lower, rec_ci_upper, rates$n_2022, "e", "Recruitment rate (%)")
p_f <- rate_plot(mortality_rate, mort_ci_lower, mort_ci_upper, rates$n_2015, "f", "Mortality rate (%)")

# ------------------------------------------------------------------------------
# Legends and layout
# ------------------------------------------------------------------------------
year_key <- ggplot() +
  geom_rect(aes(xmin = 1, xmax = 1.5, ymin = 1, ymax = 1.4), fill = c_2015_fill, color = c_2015_line, linetype = "dashed") +
  annotate("text", x = 1.7, y = 1.2, label = "2015", size = 3.2, fontface = "bold", hjust = 0, family = FONT_FAMILY) +
  geom_rect(aes(xmin = 3.0, xmax = 3.5, ymin = 1, ymax = 1.4), fill = "white", color = "grey30", linetype = "solid") +
  annotate("text", x = 3.7, y = 1.2, label = "2022", size = 3.2, fontface = "bold", hjust = 0, family = FONT_FAMILY) +
  coord_cartesian(xlim = c(0, 6), ylim = c(0.8, 1.6), clip = "off") +
  labs(title = "Survey year") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.title  = element_text(size = 9.5, face = "bold", hjust = 0.5, margin = margin(b = 2)),
        plot.margin = margin(t = 6, r = 0, b = 6, l = 0))

species_key <- ggplot() +
  geom_point(aes(x = 1.0, y = 1), color = SPECIES_COLOURS[["Madrepora oculata"]], size = 4.0) +
  geom_point(aes(x = 2.7, y = 1), color = SPECIES_COLOURS[["Desmophyllum pertusum"]], size = 4.0) +
  geom_point(aes(x = 5.0, y = 1), color = SPECIES_COLOURS[["Primnoa msp."]], size = 4.0) +
  annotate("text", x = 1.15, y = 1, label = "Madrepora oculata",
           fontface = "italic", size = 3.2, hjust = 0, family = FONT_FAMILY) +
  annotate("text", x = 2.85, y = 1, label = "Desmophyllum pertusum",
           fontface = "italic", size = 3.2, hjust = 0, family = FONT_FAMILY) +
  annotate("text", x = 5.15, y = 1, label = "italic('Primnoa')~'msp.'", parse = TRUE,
           size = 3.2, hjust = 0, family = FONT_FAMILY) +
  coord_cartesian(xlim = c(0.8, 6.2), ylim = c(0.9, 1.1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 6, r = 0, b = 4, l = 0))

fig2 <- wrap_plots(year_key, p_a, p_b, p_c, p_d, p_e, p_f, species_key,
                   design = "AAA\nBCD\nEFG\nHHH") +
  plot_layout(heights = c(0.12, 1, 1, 0.10))

save_figure(fig2, "size_structure_and_demography.png", width = 174, height = 135)

# ------------------------------------------------------------------------------
# Table 2
# ------------------------------------------------------------------------------
size_table <- sizes %>%
  group_by(Group, Year) %>%
  summarise(n = n(), p10 = quantile(Area, 0.10), p50 = median(Area), p90 = quantile(Area, 0.90),
            .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = c(n, p10, p50, p90)) %>%
  left_join(select(size_stats, Group, ks_stat, ks_p), by = "Group") %>%
  transmute(Species = Group, n_2015, n_2022, p10_2015, p10_2022, p50_2015, p50_2022,
            p90_2015, p90_2022, ks_D = round(ks_stat, 3),
            ks_p_val = if_else(ks_p < 0.001, "< 0.001", as.character(round(ks_p, 4))))

growth_table <- persistent %>%
  group_by(Species) %>%
  summarise(
    n_persistent  = n(),
    n_below_mdc   = sum(!detectable_change, na.rm = TRUE),   # open circles in panel d
    mean_growth   = mean(Annual_Area_Change, na.rm = TRUE),
    median_growth = median(Annual_Area_Change, na.rm = TRUE),
    p25_growth    = quantile(Annual_Area_Change, 0.25, na.rm = TRUE),
    p75_growth    = quantile(Annual_Area_Change, 0.75, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  left_join(select(rates, Species, mean_mdc), by = "Species")

pairwise_table <- data.frame(Comparison   = names(dunn_p),
                             P_Adjusted   = round(unname(dunn_p), 5),
                             Significance = p_to_stars(unname(dunn_p)))

demography_table <- rates %>%
  transmute(Species,
            n_2015_resolved = n_2015, n_mortality, mortality_rate, mort_ci_lower, mort_ci_upper,
            n_2022_total = n_2022, n_recruitment, recruitment_rate, rec_ci_lower, rec_ci_upper,
            n_occluded, occlusion_rate, n_occluded_2015,
            net_change_pct = recruitment_rate - mortality_rate)

cat("\nSize structure, 2015 vs 2022 (cm2; panels a-c):\n")
print(as.data.frame(size_table), row.names = FALSE)
cat("\nAnnual growth of persistent colonies (cm2 per yr) and mean annual MDC95 (panel d):\n")
print(as.data.frame(growth_table), row.names = FALSE)
cat(sprintf("\nKruskal-Wallis: H = %.2f, df = %d, p = %.4e\n",
            unname(kw_test$statistic), unname(kw_test$parameter), unname(kw_test$p.value)))
cat("\nDunn's pairwise tests (BH-adjusted):\n")
print(pairwise_table, row.names = FALSE)
cat("\nMortality and recruitment rates (%) with Wilson 95% CIs (panels e-f):\n")
print(as.data.frame(demography_table), row.names = FALSE)
cat("\nLeft out of the denominators:", sum(rates$n_occluded), "colonies hidden in 2022 and",
    sum(rates$n_occluded_2015), "hidden in 2015.\n")

write_xlsx(list(size_structure    = size_table,
                growth            = growth_table,
                growth_pairwise   = pairwise_table,
                demographic_rates = demography_table),
           in_data("tables", "size_structure_and_demography.xlsx"))
