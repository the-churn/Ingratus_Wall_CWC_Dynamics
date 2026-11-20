# ==============================================================================
# 10_colonisation_sequence.R
# Fig. 6: modelled colonisation years (2022 minus modelled age) of the three
# species, and the tests of whether the species arrived at different times:
# Kruskal-Wallis with Dunn's tests (BH-adjusted), quantile regressions of the
# earliest arrivals (5th and 10th percentiles, bootstrapped SEs), and global
# Moran's I of colonisation year per species (5 nearest neighbours).
#
# Reads:  results/colony_ages.csv (script 08)
# Writes: figures/colonisation_years.png   (Fig. 6)
#         results/colonisation_years.csv   colonies with a colonisation year,
#                                          y-axis flipped to "up the wall" and
#                                          the 95% interval of the year
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(quantreg)
  library(dunn.test)
  library(spdep)
})

SURVEY_YEAR <- 2022
BIN_YEARS   <- 5      # bar width (yr); bins end at 2022
row_h       <- 0.88   # height of each row's tallest bar or curve point (rows are 1 apart)
fig_levels  <- c("Primnoa msp.", "Desmophyllum pertusum", "Madrepora oculata")   # bottom to top

colony_ages <- read_csv(in_data("results", "colony_ages.csv"), show_col_types = FALSE)

# D. pertusum first: the reference level of the quantile regressions below
df_clean <- colony_ages %>%
  filter(Species %in% SPECIES, !is.na(Colonization_Year_Est)) %>%
  mutate(Species = factor(Species, levels = c("Desmophyllum pertusum", "Primnoa msp.",
                                              "Madrepora oculata")))

# ------------------------------------------------------------------------------
# Fig. 6: one row per species, each scaled to its own maximum
#   bars:  share of colonies per 5-yr bin (the data)
#   curve: kernel density reflected at 2022 (no colony arrived after the survey)
#   lines: median (labelled), 10th (dashed) and 5th (dotted) percentiles
# ------------------------------------------------------------------------------
yrs_by_sp <- split(df_clean$Colonization_Year_Est, as.character(df_clean$Species))[fig_levels]
bw_joint  <- mean(sapply(yrs_by_sp, bw.nrd0))   # one bandwidth for all species
cat("Kernel bandwidth (yr):", round(bw_joint, 2), "\n")
x_from    <- floor((min(df_clean$Colonization_Year_Est) - 3 * bw_joint) / 10) * 10

dens_df <- bind_rows(lapply(fig_levels, function(sp) {
  v  <- yrs_by_sp[[sp]]
  dd <- density(c(v, 2 * SURVEY_YEAR - v), bw = bw_joint, from = x_from, to = SURVEY_YEAR, n = 1024)
  data.frame(Species = sp, x = dd$x, dens = 2 * dd$y)   # the reflected sample has 2n points
}))

bin_breaks <- seq(SURVEY_YEAR - BIN_YEARS * ceiling((SURVEY_YEAR - x_from) / BIN_YEARS), SURVEY_YEAR,
                  by = BIN_YEARS)
bars_df <- bind_rows(lapply(fig_levels, function(sp) {
  h <- hist(yrs_by_sp[[sp]], breaks = bin_breaks, plot = FALSE)
  data.frame(Species = sp, xmin = head(bin_breaks, -1), xmax = bin_breaks[-1],
             dens = h$density, count = h$counts)
})) %>% filter(count > 0)

row_max <- sapply(fig_levels, function(sp)
  max(c(dens_df$dens[dens_df$Species == sp], bars_df$dens[bars_df$Species == sp])))
dens_df <- dens_df %>%
  mutate(y0 = match(Species, fig_levels), top = y0 + row_h * dens / row_max[Species]) %>%
  group_by(Species) %>%
  filter(x >= min(x[dens >= 0.005 * max(dens)])) %>%   # trim the empty far-left tail
  ungroup()
bars_df <- bars_df %>% mutate(y0 = match(Species, fig_levels), top = y0 + row_h * dens / row_max[Species])

curve_top <- function(sp, xv) {
  d <- dens_df[dens_df$Species == sp, ]
  approx(d$x, d$top, xout = xv, rule = 2)$y
}
fig_stats <- df_clean %>%
  mutate(Species = as.character(Species)) %>%
  group_by(Species) %>%
  summarise(n   = n(),
            p05 = quantile(Colonization_Year_Est, 0.05),
            p10 = quantile(Colonization_Year_Est, 0.10),
            p50 = median(Colonization_Year_Est),
            .groups = "drop") %>%
  mutate(y0 = match(Species, fig_levels)) %>%
  rowwise() %>%
  mutate(p05_top = curve_top(Species, p05), p10_top = curve_top(Species, p10),
         p50_top = curve_top(Species, p50)) %>%
  ungroup()
cat("\nMedian, 10th and 5th percentile colonisation years:\n")
print(as.data.frame(select(fig_stats, Species, n, p50, p10, p05)), row.names = FALSE)

n_of <- function(sp) fig_stats$n[fig_stats$Species == sp]
y_labels <- c(
  parse(text = sprintf("italic('Primnoa')*' msp. ('*italic(n)*' = %d)'", n_of("Primnoa msp."))),
  parse(text = sprintf("italic('Desmophyllum pertusum')*' ('*italic(n)*' = %d)'", n_of("Desmophyllum pertusum"))),
  parse(text = sprintf("italic('Madrepora oculata')*' ('*italic(n)*' = %d)'", n_of("Madrepora oculata")))
)

row_layers <- lapply(fig_levels, function(sp) list(
  geom_rect(data = bars_df[bars_df$Species == sp, ], aes(xmin = xmin, xmax = xmax, ymin = y0, ymax = top),
            fill = SPECIES_COLOURS[[sp]], alpha = 0.35, colour = "white", linewidth = 0.2),
  geom_line(data = dens_df[dens_df$Species == sp, ], aes(x = x, y = top),
            colour = SPECIES_COLOURS[[sp]], linewidth = 0.8)
))

fig6 <- ggplot() +
  row_layers +
  geom_hline(yintercept = seq_along(fig_levels), linewidth = 0.3, colour = "grey30") +
  geom_segment(data = fig_stats, aes(x = p05, xend = p05, y = y0, yend = p05_top), linetype = "dotted", linewidth = 0.5) +
  geom_segment(data = fig_stats, aes(x = p10, xend = p10, y = y0, yend = p10_top), linetype = "dashed", linewidth = 0.5) +
  geom_segment(data = fig_stats, aes(x = p50, xend = p50, y = y0, yend = p50_top), linewidth = 0.7) +
  geom_label(data = fig_stats, aes(x = p50, y = y0 + (p50_top - y0) * 0.55, label = round(p50)),
             fill = "white", colour = "black", fontface = "bold", size = 8 / .pt,
             label.padding = unit(0.15, "lines")) +
  scale_x_continuous(
    breaks   = seq(1900, 2020, by = 20),
    expand   = expansion(mult = c(0.01, 0.02)),
    sec.axis = sec_axis(~ SURVEY_YEAR - ., name = expression(bold("Modelled age in 2022 (yr)")),
                        breaks = seq(0, 140, by = 20))
  ) +
  scale_y_continuous(breaks = seq_along(fig_levels), labels = y_labels,
                     limits = c(0.95, length(fig_levels) + row_h + 0.05), expand = c(0, 0)) +
  labs(x = expression(bold("Colonisation year")), y = NULL) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.title   = element_text(size = 8.5, face = "bold", colour = "black"),
    axis.text.x  = element_text(size = 8, colour = "black"),
    axis.text.y  = element_text(size = 8.5, colour = "black", vjust = 0),
    axis.ticks.y = element_blank(),
    axis.line    = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.3, colour = "black"),
    plot.margin  = margin(t = 6, r = 10, b = 6, l = 6)
  )
save_figure(fig6, "colonisation_years.png", width = 174, height = 105)

# ------------------------------------------------------------------------------
# Did the species arrive at different times?
# ------------------------------------------------------------------------------
cat("\nKruskal-Wallis test of colonisation year:\n")
print(kruskal.test(Colonization_Year_Est ~ Species, data = df_clean))

cat("Dunn's tests (BH-adjusted):\n")
dunn_res <- dunn.test(x = df_clean$Colonization_Year_Est, g = df_clean$Species,
                      method = "bh", altp = TRUE)

# Earliest arrivals: quantile regressions at the 5th and 10th percentiles,
# against D. pertusum and then against Primnoa. Bootstrapped standard errors
# are random, so the seed keeps them identical between runs.
set.seed(2026)
quantile_fit <- function(data, tau, reference) {
  q <- rq(Colonization_Year_Est ~ Species, data = data, tau = tau, method = "fn")
  cat(sprintf("\nQuantile regression, tau = %.2f, reference %s:\n", tau, reference))
  print(summary(q, se = "boot", R = 1000)$coefficients)
}
df_primnoa_ref <- df_clean %>% mutate(Species = relevel(factor(Species), ref = "Primnoa msp."))
quantile_fit(df_clean, 0.05, "D. pertusum")
quantile_fit(df_clean, 0.10, "D. pertusum")
quantile_fit(df_primnoa_ref, 0.05, "Primnoa msp.")
quantile_fit(df_primnoa_ref, 0.10, "Primnoa msp.")

# ------------------------------------------------------------------------------
# Colonisation years with their 95% interval, y flipped to run up the wall
# (used by scripts 12 and 15), and Moran's I per species
# ------------------------------------------------------------------------------
max_y <- max(colony_ages$Centroid_y, na.rm = TRUE)
min_y <- min(colony_ages$Centroid_y, na.rm = TRUE)

colonisation <- colony_ages %>%
  filter(!is.na(Colonization_Year_Est), !is.na(Centroid_x), !is.na(Centroid_y)) %>%
  mutate(
    Centroid_y_flipped         = max_y - Centroid_y + min_y,
    Colonization_Year_Est      = as.numeric(Colonization_Year_Est),
    Age_2022_Median            = as.numeric(Age_2022_Median),
    Colonization_Year_Lower95  = 2022 - as.numeric(Age_2022_Upper95),
    Colonization_Year_Upper95  = 2022 - as.numeric(Age_2022_Lower95),
    Colonization_Year_CI_Range = Colonization_Year_Upper95 - Colonization_Year_Lower95
  )
write_csv(colonisation, in_data("results", "colonisation_years.csv"))

for (sp in SPECIES) {
  d  <- colonisation %>% filter(Species == sp)
  lw <- nb2listw(knn2nb(knearneigh(cbind(d$Centroid_x, d$Centroid_y_flipped), k = 5)), style = "W")
  cat("\nGlobal Moran's I of colonisation year -", sp, ":\n")
  print(moran.test(d$Colonization_Year_Est, lw))
}
