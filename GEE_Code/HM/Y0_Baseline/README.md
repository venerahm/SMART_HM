# Y0 Baseline Conditional Hurdle Workflow

This folder contains the current publication workflow for the conditional-first
hurdle-model GEE analysis of zero-inflated SMART data.

In this formulation, baseline outcome `Y0` is used as an ANCOVA-like baseline
covariate. The binary and positive-count models condition on `Y0`, and final
dynamic treatment regime summaries are marginalized over the baseline `Y0`
distribution. This separates baseline adjustment from the post-baseline outcome
trajectory.

## Source Files

- `HM_Data_Application_Conditional.R`: data-application driver. It constructs
  corrected surrogate monthly alcohol outcomes, fits the conditional hurdle
  model, and writes selected tables and figures.
- `Derive_Data_Like_Simulation_Values.R`: computes data-like calibration values
  used by the simulation scripts.
- `Manuscript_Simulation.R`: runs the primary manuscript simulation.
- `Simulation_Stress_Tests.R`: runs response-rate, zero-inflation,
  zero-prediction, and response-model misspecification simulations.
- `Simulation_Stress_Tests_README.md`: documents stress-test toggles and output
  conventions.
- `functions_HM_conditional.R`: conditional hurdle-model helper functions.
- `generateSMART_GEE_conditional.R`: conditional SMART data generator.
- `fit_gcomp_conditional.R`: model fitting, g-computation, and delta-method
  estimation utilities.
- `plot_conditional_results.R`: plotting helpers for simulation summaries.
- `simulation_wrapper_utils.R`: shared simulation wrapper utilities.

## Current Workflow

Run scripts from the repository root so paths like
`GEE_Code/HM/Y0_Baseline/...` resolve correctly.

1. Run `HM_Data_Application_Conditional.R` to produce current data-application
   estimates, plots, and tables.
2. Run `Derive_Data_Like_Simulation_Values.R` after the data application to
   refresh calibration values for the data-like simulation scenarios.
3. Run `Manuscript_Simulation.R` for the primary manuscript simulation.
4. Run `Simulation_Stress_Tests.R` for stress-test simulations,
   zero-prediction diagnostics, and the response-model misspecification check.

## Output Folders

`Data_Results_Cond_New_Surrogate/` contains selected current data-application
outputs and calibration files:

- `Outcome_Proportions_By_Time_Conditional.csv`
- `DTR_Estimates_Conditional_Month_4.csv`
- `DTR_Estimates_By_Model_Conditional_Month_4.csv`
- `DTR_Trajectories_Conditional_PostSplit.csv`
- `DTR_Positive_Outcome_Summary_Month_4.csv`
- `DTR_Positive_Outcome_Summary_By_Branch_Month_4.csv`
- `DTR_Sample_Support_Month_4.csv`
- `dtr_number_key.csv`
- `figure_alt_text.txt`
- `figures/`: selected EPS figure sources for the manuscript
- `Data_Like_Calibration/`: data-like simulation calibration summaries

`New_Surrogate_Results/Sim_Results_Cond/Manuscript_Primary/niter_5000/`
contains selected primary manuscript simulation summaries and the final
simulation figure source. Large raw iteration-level simulation outputs are
excluded from the publication snapshot.

## Corrected Surrogate Outcome

The current data application constructs month-one and month-two surrogate
monthly alcohol consumption outcomes from weekly reports. For each week, weekly
total drinks are calculated as drinking days multiplied by typical drinks per
drinking day. Weeks with zero drinking days are assigned zero drinks,
consistent with the skip pattern for the typical-drinks item. Month-one and
month-two surrogates are formed by averaging available weekly total-drink values
in weeks 1-4 and 5-8, multiplying by `30 / 7`, and rounding to the nearest
integer. Participants missing all weekly values in either window are excluded,
and the current data application filters baseline or month-four values greater
than or equal to 600 drinks.

## Conditional-First GLM Formulation

The observed-data models are specified directly:

```text
R | A1, Y0
B = 1(Y > 0) | R, A1, A2, Y0, time
Y | Y > 0, R, A1, A2, Y0, time
```

Responder and nonresponder branches have separate hurdle and positive-count
GLMs. DTR truth is computed by g-computation:

```text
E[ r(a1, Y0) q_R(a1, a2R, Y0) m_R(a1, a2R, Y0)
 + {1 - r(a1, Y0)} q_NR(a1, a2NR, Y0) m_NR(a1, a2NR, Y0) ]
```

The simulation estimates two DTR-level hurdle quantities:

- `p_zero`: marginal probability of a zero outcome under the DTR, `P(Y = 0)`.
- `m_plus`: marginal mean among positives, `E[Y | Y > 0]`.

The code may compute the positive-probability complement and a positive-count
numerator internally to form `m_plus`, but those are not reported as final
simulation estimands.

## Publication Cleanup

The publication file set is controlled by `PUBLICATION_MANIFEST.tsv` at the
repository root. Run:

```bash
Rscript scripts/validate_publication_manifest.R
```

to check that required files exist, active R scripts parse, and large generated
simulation outputs remain outside git.
