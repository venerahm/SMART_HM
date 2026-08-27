# HM GEE Code

This folder contains the publication workflow for the hurdle-model GEE analysis
of zero-inflated SMART data. The outcome model has two pieces:

- a binary model for whether the outcome is zero or positive;
- a positive-count model for the outcome among nonzero observations.

The publication workflow uses the `Y0_Baseline` formulation, where baseline
outcome `Y0` enters as a covariate and DTR summaries are marginalized over the
baseline `Y0` distribution.

## Folder Map

- `Y0_Baseline/`: current conditional-first manuscript workflow.
- `Y0_Baseline/Data_Results_Cond_New_Surrogate/`: selected data-application
  summaries, calibration files, and final figure sources.
- `Y0_Baseline/New_Surrogate_Results/`: selected simulation summaries and final
  figure sources.

Most scripts assume they are run from the repository root.
