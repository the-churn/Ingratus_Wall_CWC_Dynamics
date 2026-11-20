# ==============================================================================
# setup.R
# Shared settings, sourced at the top of every numbered script: the data
# folder (from config.R), the output folders, the species names and colours,
# and a helper that saves figures.
# ==============================================================================

if (!file.exists("config.R")) {
  stop("config.R not found. Run the scripts from the code folder (the one ",
       "that holds config.R and setup.R).", call. = FALSE)
}
source("config.R")                       # sets CWC_DATA_DIR

# Every input and output lives in the data folder
in_data <- function(...) file.path(CWC_DATA_DIR, ...)
for (folder in c("processed", "models", "results", "tables", "figures")) {
  dir.create(in_data(folder), showWarnings = FALSE, recursive = TRUE)
}

# Years between the two surveys (2015 and 2022)
SURVEY_INTERVAL <- 7

# The three species, in the order used in every table and figure
SPECIES <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")
SPECIES_COLOURS <- c("Madrepora oculata"     = "#722082",
                     "Desmophyllum pertusum" = "#B63679",
                     "Primnoa msp."          = "#fb9f3a")

# TagLab class names -> species names used in the paper. "Coral Recruit"
# (colonies too small to identify) is kept as it is.
standard_species <- function(x) {
  dplyr::case_when(
    x %in% c("D. pertusum", "Desmophyllum pertusum") ~ "Desmophyllum pertusum",
    grepl("^Primnoa", x)                            ~ "Primnoa msp.",
    TRUE                                            ~ x
  )
}

# The colonies the growth model is fitted to: identified colonies present in
# both surveys, with their 2015 area (A1) and relative growth rate
read_growth_data <- function() {
  readxl::read_excel(in_data("processed", "colonies_with_mdc.xlsx")) |>
    dplyr::mutate(Group = standard_species(Species)) |>
    dplyr::filter(Group %in% SPECIES) |>
    dplyr::mutate(fSpecies = factor(Group, levels = SPECIES)) |>
    dplyr::filter(state_transition == "persistence", !is.na(RGR_Planar), !is.na(A1))
}

# Save a ggplot as a PNG in figures/. Stops if the file was not rewritten
# (on Windows this happens when the PNG is open in an image viewer).
save_figure <- function(plot, file, width, height, units = "mm", dpi = 600) {
  out <- in_data("figures", file)
  started <- Sys.time()
  ggplot2::ggsave(out, plot = plot, width = width, height = height, units = units,
                  dpi = dpi, device = ragg::agg_png, bg = "white")
  if (!file.exists(out) || file.mtime(out) < started - 2) {
    stop("Could not write ", out, ". Close it if it is open in another program.",
         call. = FALSE)
  }
  cat("Saved figures/", file, "\n", sep = "")
}
