# ==============================================================================
# 11_colonisation_spatial_robustness.R
# Do the species differences in colonisation year survive spatial
# autocorrelation? (ESM Methods S4). Colonisation year is autocorrelated
# within species (script 10), so tests that treat colonies as independent can
# overstate significance. Two checks:
#   (i)  spatial block bootstrap: the quantile-regression contrasts at
#        tau = 0.05, 0.10 and 0.50 are recomputed after resampling square
#        blocks of the wall (1, 2 and 4 m) instead of single colonies (999 times)
#   (ii) within-block permutation: species labels are shuffled only among
#        colonies in the same 2-m block (9,999 times), so species are compared
#        within the same parts of the wall (Kruskal-Wallis H and differences in
#        the 5th, 10th and 50th percentiles)
# A Moran's I correlogram for M. oculata shows the distance over which
# colonisation years stay similar, which sets the block sizes.
#
# Reads:  results/colony_ages.csv (script 08)
# Writes: results/colonisation_spatial_robustness.csv
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(quantreg)
  library(spdep)
})
options(width = 160)

PX_TO_M     <- 0.00216     # orthomosaic pixel size (2.16 mm)
BLOCK_SIZES <- c(1, 2, 4)  # m
PERM_BLOCK  <- 2           # m
TAUS        <- c(0.05, 0.10, 0.50)
N_BOOT      <- 999
N_PERM      <- 9999
set.seed(2026)

# D. pertusum is the reference level, as in script 10
sp_levels <- c("Desmophyllum pertusum", "Primnoa msp.", "Madrepora oculata")
x <- read_csv(in_data("results", "colony_ages.csv"), show_col_types = FALSE) %>%
  filter(Species %in% sp_levels, !is.na(Colonization_Year_Est), !is.na(Centroid_x), !is.na(Centroid_y)) %>%
  mutate(sp = factor(Species, levels = sp_levels), yr = Colonization_Year_Est,
         xm = Centroid_x * PX_TO_M, ym = Centroid_y * PX_TO_M)
cat("Colonies:", nrow(x), "| wall extent:", round(diff(range(x$xm)), 1), "x",
    round(diff(range(x$ym)), 1), "m\n")
print(table(x$sp))

block_id <- function(size) {
  interaction(floor((x$xm - min(x$xm)) / size), floor((x$ym - min(x$ym)) / size), drop = TRUE)
}

# ------------------------------------------------------------------------------
# Correlogram (M. oculata): over what distance do colonisation years stay similar?
# ------------------------------------------------------------------------------
mo  <- x %>% filter(sp == "Madrepora oculata")
brk <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6)
cat("\nMoran's I correlogram, M. oculata (distance band in m):\n")
for (i in seq_len(length(brk) - 1)) {
  nb <- dnearneigh(cbind(mo$xm, mo$ym), brk[i], brk[i + 1])
  mt <- moran.test(mo$yr, nb2listw(nb, style = "W", zero.policy = TRUE), zero.policy = TRUE,
                   alternative = "greater")
  cat(sprintf("  %4.2f-%4.2f m: I = %5.3f, p = %.2g\n", brk[i], brk[i + 1], mt$estimate[1], mt$p.value))
}

# ------------------------------------------------------------------------------
# (i) Spatial block bootstrap of the quantile-regression contrasts
# ------------------------------------------------------------------------------
contrasts_rq <- function(dd, tau) {
  b <- coef(rq(yr ~ sp, data = dd, tau = tau, method = "fn"))
  c(Moc_vs_Dpe = unname(b["spMadrepora oculata"]),
    Pri_vs_Dpe = unname(b["spPrimnoa msp."]),
    Moc_vs_Pri = unname(b["spMadrepora oculata"] - b["spPrimnoa msp."]))
}
# Two-sided bootstrap p-value: how often the resampled contrast crosses zero
boot_p <- function(e) {
  e <- e[!is.na(e)]
  min(1, 2 * min((1 + sum(e >= 0)) / (length(e) + 1), (1 + sum(e <= 0)) / (length(e) + 1)))
}

boot_rows <- list()
for (size in BLOCK_SIZES) {
  idx <- split(seq_len(nrow(x)), block_id(size))
  ids <- names(idx)
  cat(sprintf("\nBlock bootstrap, %g-m blocks (%d occupied), %d resamples\n", size, length(ids), N_BOOT))
  for (tau in TAUS) {
    est  <- contrasts_rq(x, tau)
    reps <- t(replicate(N_BOOT, {
      dd <- x[unlist(idx[sample(ids, length(ids), replace = TRUE)], use.names = FALSE), ]
      # rq warns "solution may be nonunique" for tied quantiles in some resamples
      if (any(table(dd$sp) == 0)) rep(NA_real_, 3) else suppressWarnings(contrasts_rq(dd, tau))
    }))
    boot_rows[[length(boot_rows) + 1]] <- data.frame(
      check = "block bootstrap", block_m = size, tau = tau, contrast = names(est), estimate = est,
      se    = apply(reps, 2, sd, na.rm = TRUE),
      ci_lo = apply(reps, 2, quantile, 0.025, na.rm = TRUE),
      ci_hi = apply(reps, 2, quantile, 0.975, na.rm = TRUE),
      p     = apply(reps, 2, boot_p), row.names = NULL)
  }
}
boot_tab <- bind_rows(boot_rows)

# ------------------------------------------------------------------------------
# (ii) Within-block permutation
# ------------------------------------------------------------------------------
perm_stats <- function(sp_lab) {
  q <- sapply(TAUS, function(t) tapply(x$yr, sp_lab, quantile, probs = t))   # species x taus
  dq <- c(q["Madrepora oculata", ] - q["Desmophyllum pertusum", ],
          q["Primnoa msp.", ]      - q["Desmophyllum pertusum", ],
          q["Madrepora oculata", ] - q["Primnoa msp.", ])
  names(dq) <- paste(rep(c("Moc_vs_Dpe", "Pri_vs_Dpe", "Moc_vs_Pri"), each = length(TAUS)),
                     rep(TAUS, 3), sep = "@")
  c(KW_H = unname(kruskal.test(x$yr, sp_lab)$statistic), dq)
}
cat(sprintf("\nWithin-block permutation, %g-m blocks, %d permutations\n", PERM_BLOCK, N_PERM))
blk_perm <- split(seq_len(nrow(x)), block_id(PERM_BLOCK))
obs  <- perm_stats(x$sp)
perm <- replicate(N_PERM, {
  s <- x$sp
  for (g in blk_perm) if (length(g) > 1) s[g] <- s[g][sample.int(length(g))]
  perm_stats(s)
})
p_perm <- sapply(names(obs), function(k) {
  if (k == "KW_H") (1 + sum(perm[k, ] >= obs[k])) / (N_PERM + 1)
  else (1 + sum(abs(perm[k, ]) >= abs(obs[k]))) / (N_PERM + 1)
})
perm_tab <- data.frame(check = "within-block permutation", block_m = PERM_BLOCK,
                       tau = c(NA, rep(TAUS, 3)),
                       contrast = c("Kruskal-Wallis H", sub("@.*", "", names(obs)[-1])),
                       estimate = unname(obs), se = NA, ci_lo = NA, ci_hi = NA, p = unname(p_perm),
                       row.names = NULL)

res <- bind_rows(boot_tab, perm_tab)
write_csv(res, in_data("results", "colonisation_spatial_robustness.csv"))
cat("\nSpatial robustness of the colonisation-year comparisons (yr; negative = earlier):\n")
print(res %>% mutate(across(c(estimate, se, ci_lo, ci_hi), ~ round(.x, 1)), p = signif(p, 2)),
      row.names = FALSE)
