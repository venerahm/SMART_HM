# SMART_HM

Simulation code for evaluating conditional-first hurdle-model estimators for
zero-inflated sequential multiple assignment randomized trial (SMART) outcomes.

This repository is intentionally simulation-only. The original data application
used participant-level M-COACH data that cannot be publicly shared, so private
data-analysis scripts and generated application results are not included here.
The public code instead defines a reproducible simulation design in
`config/simulation_design.R`.

## Repository Structure

- `R/`: reusable functions for the conditional hurdle model, SMART data
  generation, g-computation estimation, simulation summaries, and plotting.
- `scripts/`: executable simulation drivers.
- `config/`: public, non-identifying simulation design inputs.
- `docs/`: notes for secondary simulation diagnostics.
- `results/`: default location for generated simulation outputs. This folder
  is ignored by git.

## Main Scripts

Run scripts from the repository root so the relative `source()` calls resolve.

```bash
Rscript scripts/run_primary_simulation.R
Rscript scripts/run_stress_tests.R
```

The primary simulation writes to:

```text
results/primary_simulation/niter_<N>/
```

The stress-test script writes to:

```text
results/stress_tests/
```

## Quick Smoke Test

The manuscript-scale simulations use many iterations. For a quick local check,
override the iteration counts with environment variables:

```bash
SMART_HM_NITER=1 Rscript scripts/run_primary_simulation.R
SMART_HM_ZERO_PREDICTION_NITER=1 Rscript scripts/run_stress_tests.R
```

You can redirect all generated outputs with:

```bash
SMART_HM_RESULTS_DIR=/path/to/results Rscript scripts/run_primary_simulation.R
```

## Simulation Design

`config/simulation_design.R` defines the public design inputs used to translate
interpretable response, zero-inflation, and positive-count settings into the
conditional hurdle GLM parameters. These values are not participant-level data.

The primary simulation evaluates dynamic treatment regimes of the form
`d = (A1, A2R, A2NR)`, where `A2R` is the second-stage treatment for responders
and `A2NR` is the second-stage treatment for nonresponders.

## Code Map

- `R/conditional_hurdle_model.R`: link functions, zero-truncated Poisson
  utilities, design matrices, and conversion from design inputs to conditional
  hurdle-model parameters.
- `R/simulate_smart_data.R`: SMART data generation under the conditional hurdle
  model.
- `R/fit_g_computation.R`: fitted conditional hurdle models and standardized
  g-computation estimators for dynamic treatment regimes.
- `R/plot_simulation_results.R`: manuscript-style summary plots.
- `R/simulation_helpers.R`: shared wrappers for iteration fitting, diagnostics,
  and summary calculations.
- `scripts/run_primary_simulation.R`: primary manuscript simulation grid.
- `scripts/run_stress_tests.R`: response-rate, zero-prediction, and
  response-misspecification diagnostics.

## Requirements

See `REQUIREMENTS.md` for R package requirements.

## Code Origins

Simulation code adapted from Nick Seewald's
[rmSMARTSize](https://github.com/nickseewald/rmSMARTSize/tree/8484e1f56b551b8f9188aaacdd3b0da1fd978593).
