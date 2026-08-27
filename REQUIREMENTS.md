# Requirements

This repository is an analysis-code publication snapshot. It includes source
code, selected result summaries, and figure sources, but it does not include
participant-level input data.

## R

The current scripts were developed for R 4.x.

Required CRAN packages for the active workflow:

- `dplyr`
- `ggplot2`
- `MASS`
- `mvtnorm`
- `pacman`
- `patchwork`
- `purrr`
- `readr`
- `scales`
- `tibble`
- `tidyr`
- `tidyverse`

`HM_Data_Application_Conditional.R` uses `pacman::p_load()` for the data
application packages. Other scripts use ordinary `library()` calls or explicit
`pkg::function()` references.

## External Data Inputs

The participant-level HM data files are not included in the publication
snapshot. To rerun the data application or refresh data-like simulation
calibration files, provide these inputs locally:

- Weekly data file: `McoachWeeklyForKellyFinal.csv`
- Outcome data file: `imp_total_1.csv`

The active data-application and calibration scripts read these paths from
environment variables:

```bash
export MCOACH_WEEKLY_PATH="/path/to/McoachWeeklyForKellyFinal.csv"
export MCOACH_OUTCOME_PATH="/path/to/imp_total_1.csv"
```

These variables are required for rerunning the data-dependent scripts because
the publication snapshot does not include participant-level data.

## Validation

To check the publication file set and parse the active R workflow:

```bash
Rscript scripts/validate_publication_manifest.R
```
