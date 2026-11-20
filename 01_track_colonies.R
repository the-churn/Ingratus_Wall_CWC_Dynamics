# ==============================================================================
# 01_track_colonies.R
# Joins the 2015 and 2022 TagLab annotations, follows each colony between the
# two surveys, assigns its fate and computes its annual growth.
#
# Reads:  OW15_1fps_full.csv, OW22_1fps_coreg_full.csv   (TagLab exports)
# Writes: processed/colony_outlines.csv   one row per outline, both years
#         processed/colonies_tracked.csv  one row per colony
# ==============================================================================

source("setup.R")
suppressPackageStartupMessages(library(tidyverse))

# ------------------------------------------------------------------------------
# Colony lists from the visual audit of the imagery and 3D models
# ------------------------------------------------------------------------------
# Outside the area covered by both surveys: dropped
outside_overlap_ids <- c(941, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1008, 1005,
                         1004, 1006, 1007, 984, 971, 975, 973, 974, 966, 965, 1086,
                         989, 993, 994, 995, 1115)

# Annotation artefacts, not colonies (recruit audit, TagLab_QuickJump_Audit.xlsx): dropped
artefact_ids <- c(1025)

# Seen in 2015 but hidden in 2022 (overgrown, buried, collapsed or in shadow).
# Fate "Occluded": censored, not counted as dead.
occluded_2022_ids <- c(46, 67, 292, 371, 447, 469, 470, 497, 498, 499, 500, 501, 502,
                       507, 599, 828)

# Seen in 2022 but hidden in 2015 (recruit audit). Fate "Occluded 2015":
# censored from recruitment, but kept in the 2022 size distributions.
occluded_2015_ids <- c(940, 943, 947, 956, 957, 969, 987, 992, 1001, 1009, 1030, 1031,
                       1058, 1059, 1060, 1062, 1105, 1106, 1108)

# ------------------------------------------------------------------------------
# Colony outlines from both surveys
# ------------------------------------------------------------------------------
read_export <- function(file, year) {
  read_csv(in_data(file), col_types = cols(`TagLab Note` = col_character())) %>%
    mutate(Year = year)
}
annotations <- bind_rows(read_export("OW15_1fps_full.csv", 2015),
                         read_export("OW22_1fps_coreg_full.csv", 2022)) %>%
  relocate(Year)

# Coral classes are drawn as outlines ("regions"); everything else is a point
coral_classes <- c("Madrepora oculata", "D. pertusum", "Coral Recruit")

outlines <- annotations %>%
  filter(TagLab.Class.name %in% coral_classes | str_detect(TagLab.Class.name, "^Primnoa")) %>%
  filter(!is.na(Image.name) & Image.name != "") %>%
  filter(!TagLab.Genet.Id %in% c(outside_overlap_ids, artefact_ids)) %>%
  # A colony keeps its TagLab genet ID across years. Outlines without one
  # (ID 0 or missing) were seen in one survey only and get their own ID.
  mutate(Join_ID = ifelse(TagLab.Genet.Id == 0 | is.na(TagLab.Genet.Id),
                          paste0(Year, "_Single_", row_number()),
                          as.character(TagLab.Genet.Id)))

write_csv(outlines, in_data("processed", "colony_outlines.csv"))

# ------------------------------------------------------------------------------
# One row per colony and year (a colony can have several outlines)
# ------------------------------------------------------------------------------
summarise_year <- function(year) {
  outlines %>%
    filter(Year == year) %>%
    group_by(Join_ID) %>%
    summarise(Species    = first(TagLab.Class.name),
              Area       = sum(TagLab.Area, na.rm = TRUE),
              Surf_Area  = sum(TagLab.Surf.area, na.rm = TRUE),
              Perimeter  = sum(TagLab.Perimeter, na.rm = TRUE),
              Centroid_x = mean(TagLab.Centroid.x, na.rm = TRUE),
              Centroid_y = mean(TagLab.Centroid.y, na.rm = TRUE),
              Poly_Count = n(),
              .groups    = "drop") %>%
    rename_with(~ paste0(.x, "_", year), -Join_ID)
}

# ------------------------------------------------------------------------------
# Fate and growth of each colony between 2015 and 2022
# ------------------------------------------------------------------------------
colonies <- full_join(summarise_year(2015), summarise_year(2022), by = "Join_ID") %>%
  mutate(
    Species = coalesce(Species_2015, Species_2022),
    n_2015  = replace_na(Poly_Count_2015, 0),
    n_2022  = replace_na(Poly_Count_2022, 0),

    # The occlusion lists come first: case_when() keeps the first match
    fate = case_when(
      Join_ID %in% occluded_2022_ids & n_2015 > 0 & n_2022 == 0 ~ "Occluded",
      Join_ID %in% occluded_2015_ids & n_2015 == 0 & n_2022 > 0 ~ "Occluded 2015",
      n_2015 == 0 & n_2022 > 0  ~ "Recruit",
      n_2015 > 0  & n_2022 == 0 ~ "Total Mortality",
      n_2015 == 1 & n_2022 > 1  ~ "Fission",
      n_2015 > 1  & n_2022 == 1 ~ "Fusion",
      n_2015 > 1  & n_2022 > 1  ~ "Complex",
      n_2015 == 1 & n_2022 == 1 ~ "Survivor",
      TRUE                      ~ "Unknown"
    ),

    # A colony absent in one survey has size 0 there. Colonies hidden in
    # 2022 were not measured, so their 2022 sizes stay NA (not 0).
    Area_2015      = replace_na(Area_2015, 0),
    Area_2022      = ifelse(fate == "Occluded", NA_real_, replace_na(Area_2022, 0)),
    Surf_Area_2015 = replace_na(Surf_Area_2015, 0),
    Surf_Area_2022 = ifelse(fate == "Occluded", NA_real_, replace_na(Surf_Area_2022, 0)),
    Perimeter_2015 = replace_na(Perimeter_2015, 0),
    Perimeter_2022 = ifelse(fate == "Occluded", NA_real_, replace_na(Perimeter_2022, 0)),
    Centroid_x     = coalesce(Centroid_x_2015, Centroid_x_2022),
    Centroid_y     = coalesce(Centroid_y_2015, Centroid_y_2022),

    # Annual change and relative growth rate over the 7 years. None for
    # colonies hidden in either survey.
    hidden = fate %in% c("Occluded", "Occluded 2015"),
    Annual_Area_Change      = ifelse(hidden, NA_real_, (Area_2022 - Area_2015) / SURVEY_INTERVAL),
    Annual_Surf_Area_Change = ifelse(hidden, NA_real_, (Surf_Area_2022 - Surf_Area_2015) / SURVEY_INTERVAL),
    RGR_Planar  = ifelse(!hidden & Area_2015 > 0 & Area_2022 > 0,
                         (log(Area_2022) - log(Area_2015)) / SURVEY_INTERVAL, NA),
    RGR_Surface = ifelse(!hidden & Surf_Area_2015 > 0 & Surf_Area_2022 > 0,
                         (log(Surf_Area_2022) - log(Surf_Area_2015)) / SURVEY_INTERVAL, NA)
  )

# Every colony on an occlusion list should have received that fate
check_list <- function(ids, expected) {
  listed <- colonies %>% filter(Join_ID %in% ids)
  wrong   <- listed$Join_ID[listed$fate != expected]
  missing <- setdiff(as.character(ids), listed$Join_ID)
  if (length(wrong) > 0) {
    warning("Not given the fate '", expected, "' (check the list): ",
            paste(wrong, collapse = ", "), call. = FALSE)
  }
  if (length(missing) > 0) {
    warning("On the '", expected, "' list but not in the data: ",
            paste(missing, collapse = ", "), call. = FALSE)
  }
  sum(listed$fate == expected)
}
n_hidden_2022 <- check_list(occluded_2022_ids, "Occluded")
n_hidden_2015 <- check_list(occluded_2015_ids, "Occluded 2015")

colonies_tracked <- colonies %>%
  select(TagLab.Genet.Id = Join_ID, Species, fate, Centroid_x, Centroid_y,
         Area_2015, Area_2022, Annual_Area_Change, RGR_Planar,
         Surf_Area_2015, Surf_Area_2022, Annual_Surf_Area_Change, RGR_Surface,
         Perimeter_2015, Perimeter_2022)

write_csv(colonies_tracked, in_data("processed", "colonies_tracked.csv"))

cat(nrow(colonies_tracked), "colonies tracked;", n_hidden_2022, "hidden in 2022 and",
    n_hidden_2015, "hidden in 2015;", length(artefact_ids), "artefact removed.\n\nFates:\n")
print(table(colonies_tracked$fate))
