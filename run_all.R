# ==============================================================================
# run_all.R
# Runs the whole analysis in order, from the four raw inputs in the data folder
# set in config.R. Each script runs in its own fresh R session, exactly as if
# run by hand, and writes its log to <data folder>/run_logs/. The run stops at
# the first error. At the end it checks that every expected output was written
# by this run and records the R and package versions (sessionInfo.txt).
#
# Usage: 1. Put the four raw inputs in a folder and set CWC_DATA_DIR in config.R.
#        2. Set R's working directory to this folder, then
#             source("run_all.R")     (RStudio or the R console), or
#             Rscript run_all.R       (a terminal)
#        To resume after an error, set START_AT to the script that failed.
#        To run some scripts only, list them in RUN_ONLY.
# ==============================================================================

START_AT <- "01_track_colonies.R"
RUN_ONLY <- character(0)          # e.g. c("10_colonisation_sequence.R")

if (!file.exists("config.R")) stop("Set the working directory to the code folder first.", call. = FALSE)
source("config.R")
data_dir <- CWC_DATA_DIR
log_dir  <- file.path(data_dir, "run_logs")

raw_inputs <- c("OW15_1fps_full.csv", "OW22_1fps_coreg_full.csv", "CNN Error.xlsx",
                "TagLab_QuickJump_Audit.xlsx")

# Run order, with the outputs each script writes (relative to the data folder)
pipeline <- list(
  "01_track_colonies.R"                   = c("processed/colony_outlines.csv", "processed/colonies_tracked.csv"),
  "02_measurement_error.R"                = c("processed/colonies_with_mdc.xlsx"),
  "03_size_structure_and_demography.R"    = c("figures/size_structure_and_demography.png",
                                              "tables/size_structure_and_demography.xlsx"),
  "04_size_selective_mortality.R"         = c("figures/mortality_by_size.png", "results/mortality_size_tests.csv"),
  "05_growth_model.R"                     = c("models/growth_model.rds", "models/growth_model_candidates.rds",
                                              "tables/growth_model.xlsx", "figures/growth_model_diagnostics.png"),
  # the last four are written only when the residuals are autocorrelated (they are)
  "06_growth_model_spatial_check.R"       = c("results/growth_residuals_morans_i.csv",
                                              "results/growth_spatial_refit_comparison.csv",
                                              "results/growth_spatial_refit_shifts.csv",
                                              "figures/growth_spatial_refit_curves.png",
                                              "figures/growth_residual_semivariograms.png"),
  "07_growth_curves.R"                    = c("models/athresh_sizes.rds", "figures/growth_curves.png"),
  "08_colony_ages.R"                      = c("results/colony_ages.csv", "models/age_simulation_inputs.rds"),
  "09_age_figure.R"                       = c("figures/age_model.png"),
  "10_colonisation_sequence.R"            = c("figures/colonisation_years.png", "results/colonisation_years.csv"),
  "11_colonisation_spatial_robustness.R"  = c("results/colonisation_spatial_robustness.csv"),
  "12_colonisation_hotspots.R"            = c("results/colonisation_hotspots.csv",
                                              "figures/colonisation_hotspots.png"),
  "13_recruit_identification_check.R"     = c("results/recruit_threshold.csv",
                                              "results/recruit_identification_ks.csv",
                                              "figures/recruit_identification_check.png"),
  "14_age_model_checks.R"                 = c("results/age_elapsed_time_check.csv",
                                              "results/a0_sensitivity_ages.csv", "figures/a0_sensitivity.png"),
  "15_age_bias_sensitivity.R"             = c("results/age_bias_colonisation.csv", "results/age_bias_grid.csv",
                                              "figures/age_bias_sensitivity.png"),
  "16_variance_structure_sensitivity.R"   = c("tables/variance_structure_sensitivity.xlsx",
                                              "figures/variance_structure_sensitivity.png"),
  "17_normality_checks.R"                 = c("results/normality_shapiro_wilk.csv", "figures/normality_qq.png",
                                              "figures/growth_variables_dotplots.png"),
  "18_manuscript_numbers.R"               = character(0)   # prints to its log only
)

# Packages used by the scripts (versions recorded at the end, for ESM Table S4)
packages <- c("tidyverse", "readxl", "writexl", "openxlsx", "nlme", "MASS", "dunn.test",
              "Hmisc", "quantreg", "spdep", "patchwork", "ggsignif", "scales", "ragg")

# ------------------------------------------------------------------------------
# Checks before running
# ------------------------------------------------------------------------------
scripts <- names(pipeline)
missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install these packages first: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
missing_inputs <- raw_inputs[!file.exists(file.path(data_dir, raw_inputs))]
if (length(missing_inputs) > 0) {
  stop("Copy these raw inputs into ", data_dir, " first:\n  ",
       paste(missing_inputs, collapse = "\n  "), call. = FALSE)
}
if (length(RUN_ONLY) > 0) {
  unknown <- setdiff(RUN_ONLY, scripts)
  if (length(unknown) > 0) stop("Not in the pipeline: ", paste(unknown, collapse = ", "), call. = FALSE)
  to_run <- scripts[scripts %in% RUN_ONLY]
} else {
  if (!START_AT %in% scripts) stop("START_AT is not one of the scripts: ", START_AT, call. = FALSE)
  to_run <- scripts[match(START_AT, scripts):length(scripts)]
}

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
rscript   <- file.path(R.home("bin"), "Rscript")
run_start <- Sys.time()
cat("Data folder:", data_dir, "\nRunning", length(to_run), "scripts (about 15 minutes for all)\n\n")

# ------------------------------------------------------------------------------
# Run each script in a fresh R session
# ------------------------------------------------------------------------------
timings <- data.frame(script = to_run, minutes = NA_real_, status = "not run")
for (i in seq_along(to_run)) {
  s        <- to_run[i]
  log_file <- file.path(log_dir, sub("\\.R$", ".log", s))
  cat(sprintf("[%2d/%d] %-42s ", i, length(to_run), s))
  t0 <- Sys.time()
  exit_code <- system2(rscript, args = shQuote(s), stdout = log_file, stderr = log_file)
  timings$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  if (!identical(as.integer(exit_code), 0L)) {
    timings$status[i] <- "FAILED"
    write.csv(timings, file.path(log_dir, "timings.csv"), row.names = FALSE)
    cat("FAILED\n\n---- last lines of", basename(log_file), "----\n")
    cat(tail(readLines(log_file, warn = FALSE), 25), sep = "\n")
    stop("Stopped at ", s, ". Fix it, then set START_AT <- \"", s, "\" and run again.", call. = FALSE)
  }
  timings$status[i] <- "ok"
  cat(sprintf("ok (%.1f min)\n", timings$minutes[i]))
}
write.csv(timings, file.path(log_dir, "timings.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Every expected output must be newer than the start of this run (on Windows,
# a figure open in an image viewer cannot be overwritten)
# ------------------------------------------------------------------------------
expected <- unlist(pipeline[to_run], use.names = FALSE)
mtime    <- file.mtime(file.path(data_dir, expected))
stale    <- expected[is.na(mtime) | mtime < run_start]
if (length(stale) == 0) {
  cat("\nAll", length(expected), "expected outputs were written by this run.\n")
} else {
  cat("\nWARNING: missing or not rewritten by this run:\n  ", paste(stale, collapse = "\n  "), "\n", sep = "")
}

# R and package versions, recorded in a separate session so the packages are
# not left attached in this one
si_file   <- normalizePath(file.path(log_dir, "sessionInfo.txt"), winslash = "/", mustWork = FALSE)
si_script <- tempfile(fileext = ".R")
writeLines(c(
  sprintf("pkgs <- c(%s)", paste(sprintf('"%s"', packages), collapse = ", ")),
  "invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(require(p, character.only = TRUE, quietly = TRUE))))",
  sprintf("writeLines(capture.output(sessionInfo()), %s)", deparse(si_file))
), si_script)
system2(rscript, args = shQuote(si_script), stdout = FALSE, stderr = FALSE)

cat(sprintf("\nDone in %.1f min. Logs, timings and package versions: %s\n",
            as.numeric(difftime(Sys.time(), run_start, units = "mins")), log_dir))
