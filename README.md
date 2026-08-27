# Zero Inflated SMART

Publication repository for the conditional-first hurdle-model GEE analysis of
zero-inflated SMART data.

The active manuscript workflow is the `Y0_Baseline` analysis, where baseline
outcome `Y0` is used as an ANCOVA-like baseline covariate and dynamic treatment
regime summaries are marginalized over the baseline `Y0` distribution.

Most scripts assume they are run from the repository root so paths beginning
with `GEE_Code/HM/...` resolve correctly.

## Repository Structure

- `GEE_Code/HM/`: hurdle-model GEE code for the HM application.
- `GEE_Code/HM/Y0_Baseline/`: current conditional-first manuscript workflow.
- `GEE_Code/HM/Y0_Baseline/Data_Results_Cond_New_Surrogate/`: selected current
  data-application summaries, figures, and simulation calibration files.
- `GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/`: selected current simulation
  summaries and manuscript figures.
- `scripts/`: manifest validation and publication snapshot utilities.

## Current Workflow

Run these scripts from the repository root:

1. `GEE_Code/HM/Y0_Baseline/HM_Data_Application_Conditional.R`
   constructs the corrected surrogate monthly alcohol outcomes, fits the
   conditional hurdle model to the data application, and writes current tables
   and figures.
2. `GEE_Code/HM/Y0_Baseline/Derive_Data_Like_Simulation_Values.R`
   derives data-like calibration values used by the simulations.
3. `GEE_Code/HM/Y0_Baseline/Manuscript_Simulation.R`
   runs the primary manuscript simulation.
4. `GEE_Code/HM/Y0_Baseline/Simulation_Stress_Tests.R`
   runs stress-test, zero-prediction, and response-misspecification
   simulations.

The main helper files are:

- `functions_HM_conditional.R`
- `generateSMART_GEE_conditional.R`
- `fit_gcomp_conditional.R`
- `plot_conditional_results.R`
- `simulation_wrapper_utils.R`

## Publication Contents

This repository keeps source code, documentation, selected summary CSV files,
and final EPS figure sources. Large raw simulation outputs, local workspace
files, and development archives are intentionally excluded from the publication
snapshot.

See `REQUIREMENTS.md` for R package requirements and the external data input
paths needed to rerun the data application.

Use the manifest to check the exact publication file set:

```bash
Rscript scripts/validate_publication_manifest.R
```

To create a clean publication repository snapshot:

```bash
Rscript scripts/create_publication_snapshot.R /path/to/new/repo --init-git
```

Add `--include-optional` to include development-reference files listed as
optional in `PUBLICATION_MANIFEST.tsv`. The snapshot script initializes git on
branch `main` by default. To attach a GitHub remote at creation time, add:

```bash
--remote-url=https://github.com/<user>/<repo>.git
```

To also create the initial local commit:

```bash
--commit-message="Initial publication-ready repository"
```

## Cleanup Metadata

- `PUBLICATION_CLEANUP_INVENTORY.md` documents which files and folders are
  current, legacy/reference, or generated/local.
- `PUBLICATION_MANIFEST.tsv` is the machine-checkable keep/archive/exclude
  list used by the validation and snapshot scripts.
- `RELEASE_CHECKLIST.md` gives the final validation, snapshot, commit, and push
  steps for the publication repository.

## Code Origins

Simulation code adapted from Nick Seewald's
[rmSMARTSize](https://github.com/nickseewald/rmSMARTSize/tree/8484e1f56b551b8f9188aaacdd3b0da1fd978593).
