# Publication Cleanup Inventory

This inventory separates the current manuscript workflow from legacy,
diagnostic, and generated files. It is intended to guide creation of a smaller
publication-ready repository.

## Current Manuscript Workflow

The current workflow is the conditional-first `Y0_Baseline` analysis with `Y0`
used as a baseline covariate.

Keep these source files:

- `README.md`
- `GEE_Code/HM/README.md`
- `GEE_Code/HM/Y0_Baseline/README.md`
- `GEE_Code/HM/Y0_Baseline/HM_Data_Application_Conditional.R`
- `GEE_Code/HM/Y0_Baseline/Derive_Data_Like_Simulation_Values.R`
- `GEE_Code/HM/Y0_Baseline/Manuscript_Simulation.R`
- `GEE_Code/HM/Y0_Baseline/Simulation_Stress_Tests.R`
- `GEE_Code/HM/Y0_Baseline/Simulation_Stress_Tests_README.md`
- `GEE_Code/HM/Y0_Baseline/functions_HM_conditional.R`
- `GEE_Code/HM/Y0_Baseline/generateSMART_GEE_conditional.R`
- `GEE_Code/HM/Y0_Baseline/fit_gcomp_conditional.R`
- `GEE_Code/HM/Y0_Baseline/plot_conditional_results.R`
- `GEE_Code/HM/Y0_Baseline/simulation_wrapper_utils.R`

Current data-application output root:

- `GEE_Code/HM/Y0_Baseline/Data_Results_Cond_New_Surrogate/`

Current simulation output root:

- `GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/`

## Development Reference

These files or folders are useful for methods history, checks, or development
reference, but they are not the current manuscript driver:

- `GEE_Code/HM/Y0_Baseline/Tests_conditional.R`
- `GEE_Code/HM/Y0_Baseline/Projection/`
- `GEE_Code/HM/Y0_Baseline/Archive/`
- `GEE_Code/HM/Y0_Baseline/Data_Results/`
- `GEE_Code/HM/Y0_Baseline/Sim_Results/`
- `GEE_Code/HM/Y0_Baseline/Sim_Results_Cond/`
- `GEE_Code/HM/Y0_Outcome/`
- `GEE_Code/HM/NB/`
- `GEE_Code/HM/Code_Checks/`
- `GEE_Code/HM/Archive/`
- `GLMM_Code/`
- `Batch/`
- `logs/`

## Generated Or Local Files To Exclude

Exclude these from the publication clone unless there is a specific archival
reason to include them:

- `.DS_Store`
- `.Rhistory`
- `.RData`
- `.RDataTmp`
- `.RDataTmp1`
- `.Rproj.user/`
- `*.aux`
- `*.log`
- large simulation `.rds` files
- temporary run logs

## Notes Before Cloning

- The active manuscript simulation writes to
  `GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/Sim_Results_Cond/`.
- A clean publication repo can be created from the manifest with:

  ```bash
  Rscript scripts/create_publication_snapshot.R /path/to/new/repo --init-git
  ```

  Add `--remote-url=https://github.com/<user>/<repo>.git` when the target
  GitHub repository already exists. Add
  `--commit-message="Initial publication-ready repository"` to create the
  initial local commit in the clean snapshot.
- Several older scripts in the reference folders source legacy files that have
  been moved or deleted from the active top-level workflow. Treat those scripts
  as archival unless they are deliberately restored as a separate reproducible
  chapter workflow.
- For a publication repository, prefer including source code, documentation,
  small summary CSV files, and final manuscript figures. Keep large raw
  simulation outputs outside git or attach them through an external archive.
