# Cold-water coral demography on a Whittard Canyon wall, 2015–2022: data and R code

Data and R code for the paper on the growth, survival and colonisation history of *Madrepora oculata*, *Desmophyllum pertusum* and *Primnoa* msp. on a vertical wall at 1,080 m depth in Whittard Canyon (NE Atlantic), from repeat photogrammetry surveys in 2015 and 2022. Running the code reproduces, from the raw annotations, every result in the paper and its Electronic Supplementary Material (ESM), and every figure except Fig. 1.

**Citation:** [authors, year, title, *Coral Reefs*, DOI — added on publication]

## Contents

```
README.md
config.R          the data folder (one line)
setup.R           shared settings, sourced by every script
age_model.R       the age model, sourced by the scripts that compute ages
run_all.R         runs every script in order
01_ ... 18_       the analysis, one script per step
data/             the raw inputs and the modelled colony ages; all outputs are written here
```

### Data (`data/`)

| File | What it is |
|---|---|
| `OW15_1fps_full.csv`, `OW22_1fps_coreg_full.csv` | TagLab exports of the 2015 and 2022 orthomosaics (2.16 mm per pixel): one row per annotation, with class, genet ID (the same colony keeps its ID in both years), area, 3D surface area, perimeter and centroid (pixels) |
| `CNN Error.xlsx` | Repeat CNN segmentations of the same colonies, used for the measurement-error model and each colony's 95% minimum detectable change (MDC95) |
| `TagLab_QuickJump_Audit.xlsx` | Visual check in TagLab of the colonies seen only in 2022 whose modelled age was 10 years or more: `growth` (present in 2015), `occluded` (hidden in 2015) or `artifact` (not a colony) |
| `colony_ages.csv` | Not an input: a copy of `results/colony_ages.csv` from script 08, with every colony's modelled age in 2015 and 2022 (median and 95% interval) and its colonisation year, for anyone who wants the ages without running the model |

## How to run

1. Install R (4.3 or later) and the packages:
   ```r
   install.packages(c("tidyverse", "readxl", "writexl", "openxlsx", "nlme", "MASS", "dunn.test",
                      "Hmisc", "quantreg", "spdep", "patchwork", "ggsignif", "scales", "ragg"))
   ```
2. Open R in this folder (in RStudio: *Session > Set Working Directory > To Source File Location* with `run_all.R` open).
3. Run `source("run_all.R")`, or `Rscript run_all.R` from a terminal.

The whole run takes about 15 minutes; most of it is the Monte Carlo age model (scripts 08 and 14) and the spatial refit of the growth model (script 06). Each script runs in a fresh R session and can also be run on its own once the scripts before it have run. Every random step has a fixed seed, so the results are identical on every run.

To keep the data elsewhere, change the path in `config.R`.

## The scripts

| Script | What it does | In the paper |
|---|---|---|
| `01_track_colonies.R` | Joins the 2015 and 2022 annotations, follows each colony, assigns its fate (survived, died, recruited, fission, fusion, or hidden in one survey) and computes annual growth | Methods |
| `02_measurement_error.R` | Measurement-error model from the repeat segmentations; each colony's MDC95 | Methods |
| `03_size_structure_and_demography.R` | Size distributions in 2015 and 2022 (KS tests), annual growth (Kruskal–Wallis, Dunn), recruitment and mortality rates (Wilson CIs) | Fig. 2 |
| `04_size_selective_mortality.R` | 2015 size of colonies that died, shrank or survived (Wilcoxon tests) | Fig. 3 |
| `05_growth_model.R` | Asymptotic growth model (gnls) and the choice of its variance structure by AIC | Table S3, Fig. S2, Methods S2 |
| `06_growth_model_spatial_check.R` | Moran's I of the growth-model residuals; refit with a spatial correlation structure | Fig. S3, Methods S2 |
| `07_growth_curves.R` | Growth curves, asymptotes and A<sub>thresh</sub> | Fig. 4 |
| `08_colony_ages.R` | Colony ages from size, with 5,000 parameter draws | Methods S3 |
| `09_age_figure.R` | How the age model works; size at age | Fig. 5 |
| `10_colonisation_sequence.R` | Colonisation years and the tests of arrival order (Kruskal–Wallis, Dunn, quantile regression, Moran's I) | Fig. 6 |
| `11_colonisation_spatial_robustness.R` | Spatial block bootstrap and within-block permutation of the colonisation comparisons | Methods S4 |
| `12_colonisation_hotspots.R` | Getis–Ord Gi* clusters of early and late colonisation | |
| `13_recruit_identification_check.R` | Whether the size shift is an artefact of how small recruits were identified | Fig. S1, Methods S1 |
| `14_age_model_checks.R` | Modelled vs true time between surveys; ages with the detection floor doubled | Methods S3 |
| `15_age_bias_sensitivity.R` | The colonisation comparisons if the age model overestimates ages | Fig. S4, Methods S3 |
| `16_variance_structure_sensitivity.R` | Parameters and ages under the next best variance structures | Methods S2 |
| `17_normality_checks.R` | Shapiro–Wilk tests, Q–Q plots and dotplots behind the choice of tests | Methods S2 |
| `18_manuscript_numbers.R` | Other numbers quoted in the Results and Discussion | Results, Discussion |

## Outputs

Everything is written inside the data folder:

| Folder | Contents |
|---|---|
| `processed/` | The colony tables: `colonies_tracked.csv` (script 01) and `colonies_with_mdc.xlsx` (script 02, used by all later scripts) |
| `models/` | The fitted growth model, A<sub>thresh</sub> and the parameter draws of the age model |
| `results/` | Tables as csv, including `colony_ages.csv` (every colony's modelled age with its 95% interval) |
| `tables/` | Tables as xlsx (one sheet per table) |
| `figures/` | Every figure, as PNG |
| `run_logs/` | One log per script with every printed result, the run times, and `sessionInfo.txt` (R and package versions, ESM Table S4) |

The logs hold all the statistics quoted in the paper; script 18 collects those not printed elsewhere.
