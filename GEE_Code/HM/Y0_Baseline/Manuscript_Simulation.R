#####################################################################
# Manuscript Conditional-First Hurdle SMART Simulation: Primary Grid #
#####################################################################
#
# Primary manuscript scenarios for the exact conditional hurdle model with
# Y0 handled as a baseline adjustment covariate.
#
# Scenario families:
#   A. Data-like primary
#   B. Correlation sensitivity
#   C. Zero-inflation sensitivity
#
# This script follows Tests_conditional.R:
#   1. Generate observed SMART data from exact conditional hurdle GLMs.
#   2. Compute DTR truth by g-computation.
#   3. Fit the same observed-data conditional GLMs.
#   4. Estimate DTR means by standardized g-computation.
#   5. Compute bias and Wald coverage using delta-method SEs.

suppressPackageStartupMessages({
  library(stats)
})

source("GEE_Code/HM/Y0_Baseline/functions_HM_conditional.R")
source("GEE_Code/HM/Y0_Baseline/generateSMART_GEE_conditional.R")
source("GEE_Code/HM/Y0_Baseline/fit_gcomp_conditional.R")
source("GEE_Code/HM/Y0_Baseline/plot_conditional_results.R")
source("GEE_Code/HM/Y0_Baseline/simulation_wrapper_utils.R")

# -----------------------
# User-configurable block
# -----------------------

n <- 400
niter <- 5000
progress_every <- max(1, floor(niter / 10))
times <- c(1, 2, 4)
spltime <- 1
t_eval <- max(times)
seed <- 20260706
use_delta_se <- TRUE

output_subdir <- paste0("niter_", niter)
output_dir <- file.path(
  "GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/Sim_Results_Cond",
  "Manuscript_Primary",
  output_subdir
)

calibration_dir <- "GEE_Code/HM/Y0_Baseline/Data_Results_Cond_New_Surrogate/Data_Like_Calibration"
design_params <- load_data_like_design_params(calibration_dir)
Y0_ref <- load_data_like_Y0_ref(calibration_dir, default = 12)

# Dynamic treatment regimes d = (A1, A2R, A2NR).
dtr_grid <- expand.grid(A1 = c(-1, 1), A2R = c(-1, 1), A2NR = c(-1, 1))

# Large reference population for exact g-computation truth.
set.seed(seed)
Y0_truth <- gen_Y0_conditional(100000)

scenario_grid <- data.frame(
  scenario_label = c("A1", "A2", "B1", "B2", "C1"),
  scenario_family = c(
    "A. Correlated outcomes, rho 0.3",
    "A. Correlated outcomes, rho 0.3",
    "B. Correlated outcomes, rho 0.5",
    "B. Correlated outcomes, rho 0.5",
    "C. Independent outcomes"
  ),
  n = n,
  ZI = c(37, 37, 37, 37, 37),
  corstr = c("exchangeable", "ar1", "exchangeable", "ar1", "independence"),
  rho = c(0.3, 0.3, 0.5, 0.5, 0),
  response_truth = "A1 + Y0",
  positive_count_truth = "ZTP",
  primary_fit = "Correct HM",
  why_include = c(
    "Application-like zero inflation with moderate exchangeable within-person correlation.",
    "Application-like zero inflation with moderate AR(1) within-person correlation.",
    "Application-like zero inflation with stronger exchangeable within-person correlation.",
    "Application-like zero inflation with stronger AR(1) within-person correlation.",
    "Application-like zero inflation without within-person outcome correlation."
  ),
  stringsAsFactors = FALSE
)
scenario_grid$scenario_id <- paste0(
  scenario_grid$scenario_label,
  "_ZI_", scenario_grid$ZI,
  "_", scenario_grid$corstr,
  "_rho_", scenario_grid$rho
)

truth_list <- vector("list", nrow(scenario_grid))
design_list <- vector("list", nrow(scenario_grid))
design_check_list <- vector("list", nrow(scenario_grid))
parameter_list <- vector("list", nrow(scenario_grid))
estimate_list <- vector("list", nrow(scenario_grid))
cell_diagnostic_list <- vector("list", nrow(scenario_grid))
error_list <- vector("list", nrow(scenario_grid))

immediate_cat("============================================================")
immediate_cat("Manuscript primary conditional-first hurdle SMART simulation")
immediate_cat("Scenarios: ", nrow(scenario_grid), " | iterations per scenario: ", niter)
immediate_cat("n: ", n, " | times: ", paste(times, collapse = ", "), " | t_eval: ", t_eval)
immediate_cat("Y0_ref for calibrated design inputs: ", Y0_ref)
immediate_cat("Output directory: ", output_dir)
immediate_cat("Progress printed every ", progress_every, " iteration(s).")
immediate_cat("============================================================")

for (ss in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[ss, , drop = FALSE]
  scenario_start <- Sys.time()
  
  immediate_cat("------------------------------------------------------------")
  immediate_cat("Scenario ", ss, "/", nrow(scenario_grid), ": ", scenario$scenario_id)
  immediate_cat("Family: ", scenario$scenario_family)
  immediate_cat("ZI: ", scenario$ZI, "% | corstr: ", scenario$corstr,
                " | rho: ", scenario$rho)
  immediate_cat("Converting design inputs to conditional GLM parameters and truth...")
  
  params <- conditional_params_from_design(
    design_params = design_params,
    zero_probability = scenario$ZI / 100,
    Y0_ref = Y0_ref,
    time_ref = t_eval,
    spltime = spltime
  )
  
  truth_s <- conditional_dtr_truth(
    params = params,
    dtr_grid = dtr_grid,
    Y0_vec = Y0_truth,
    times = times,
    spltime = spltime
  )
  truth_s$scenario_id <- scenario$scenario_id
  truth_s$scenario_label <- scenario$scenario_label
  truth_s$scenario_family <- scenario$scenario_family
  truth_s$ZI <- scenario$ZI
  truth_s$corstr <- scenario$corstr
  truth_s$rho <- scenario$rho
  truth_list[[ss]] <- truth_s
  
  design_s <- conditional_design_table(design_params)
  design_s$scenario_id <- scenario$scenario_id
  design_s$scenario_label <- scenario$scenario_label
  design_s$scenario_family <- scenario$scenario_family
  design_s$ZI <- scenario$ZI
  design_s$corstr <- scenario$corstr
  design_s$rho <- scenario$rho
  design_s$zero_probability <- scenario$ZI / 100
  design_s$Y0_ref <- Y0_ref
  design_s$time_ref <- t_eval
  design_list[[ss]] <- design_s
  
  design_check_s <- conditional_design_check_table(
    params = params,
    Y0_ref = Y0_ref,
    time_ref = t_eval,
    spltime = spltime
  )
  design_check_s$scenario_id <- scenario$scenario_id
  design_check_s$scenario_label <- scenario$scenario_label
  design_check_s$scenario_family <- scenario$scenario_family
  design_check_s$ZI <- scenario$ZI
  design_check_s$corstr <- scenario$corstr
  design_check_s$rho <- scenario$rho
  design_check_s$zero_probability <- scenario$ZI / 100
  design_check_s$Y0_ref <- Y0_ref
  design_check_s$time_ref <- t_eval
  design_check_list[[ss]] <- design_check_s
  
  params_s <- conditional_parameter_table(params)
  params_s$scenario_id <- scenario$scenario_id
  params_s$scenario_label <- scenario$scenario_label
  params_s$scenario_family <- scenario$scenario_family
  params_s$ZI <- scenario$ZI
  params_s$corstr <- scenario$corstr
  params_s$rho <- scenario$rho
  parameter_list[[ss]] <- params_s
  
  iteration_estimates <- vector("list", niter)
  iteration_cell_diagnostics <- vector("list", niter)
  iteration_errors <- vector("list", niter)
  
  for (iter in seq_len(niter)) {
    if (iter == 1 || iter == niter || iter %% progress_every == 0) {
      immediate_cat(sprintf(
        "  Scenario %d/%d (%s: ZI=%s%%, corstr=%s, rho=%s) - iteration %d/%d",
        ss, nrow(scenario_grid), scenario$scenario_label,
        scenario$ZI, scenario$corstr, scenario$rho, iter, niter
      ))
    }
    
    set.seed(seed + ss * 100000 + iter)
    generated <- generateSMART_gee_conditional(
      n = scenario$n,
      times = times,
      spltime = spltime,
      params = params,
      corstr = scenario$corstr,
      rho = scenario$rho
    )
    
    cell_diag <- positive_cell_summary(generated$long_data, spltime = spltime)
    cell_diag$iteration <- iter
    cell_diag$scenario_id <- scenario$scenario_id
    cell_diag$scenario_label <- scenario$scenario_label
    cell_diag$scenario_family <- scenario$scenario_family
    cell_diag$ZI <- scenario$ZI
    cell_diag$corstr <- scenario$corstr
    cell_diag$rho <- scenario$rho
    iteration_cell_diagnostics[[iter]] <- cell_diag
    
    iter_result <- fit_conditional_iteration_for_wrapper(
      generated = generated,
      spltime = spltime,
      dtr_grid = dtr_grid,
      t_eval = t_eval,
      use_delta_se = use_delta_se,
      scenario = scenario,
      iter = iter
    )
    iteration_estimates[[iter]] <- iter_result$estimate
    iteration_errors[[iter]] <- iter_result$error
  }
  
  estimate_list[[ss]] <- do.call(rbind, iteration_estimates)
  cell_diagnostic_list[[ss]] <- do.call(rbind, iteration_cell_diagnostics)
  error_list[[ss]] <- do.call(rbind, iteration_errors)
  
  scenario_elapsed <- round(as.numeric(difftime(Sys.time(), scenario_start, units = "mins")), 2)
  immediate_cat("Scenario ", ss, "/", nrow(scenario_grid), " completed in ",
                scenario_elapsed, " minutes.")
}

truth <- do.call(rbind, truth_list)
design_inputs <- do.call(rbind, design_list)
design_checks <- do.call(rbind, design_check_list)
parameters <- do.call(rbind, parameter_list)
estimates <- do.call(rbind, estimate_list)
cell_diagnostics <- do.call(rbind, cell_diagnostic_list)
errors <- do.call(rbind, error_list)

truth_eval <- truth[truth$time == t_eval, c(
  "scenario_id", "scenario_label", "scenario_family",
  "ZI", "corstr", "rho", "DTR",
  "p_zero_true", "m_plus_true"
), drop = FALSE]

results <- merge(
  estimates,
  truth_eval,
  by = c("scenario_id", "scenario_label", "scenario_family", "ZI", "corstr", "rho", "DTR"),
  all.x = TRUE
)
results$p_zero_bias <- results$p_zero - results$p_zero_true
results$m_plus_bias <- results$m_plus - results$m_plus_true

if (all(c("p_zero_se", "m_plus_se") %in% names(results))) {
  results$p_zero_lower <- results$p_zero - 1.96 * results$p_zero_se
  results$p_zero_upper <- results$p_zero + 1.96 * results$p_zero_se
  results$m_plus_lower <- results$m_plus - 1.96 * results$m_plus_se
  results$m_plus_upper <- results$m_plus + 1.96 * results$m_plus_se
  
  results$p_zero_covered <- as.numeric(
    results$p_zero_lower <= results$p_zero_true & results$p_zero_true <= results$p_zero_upper
  )
  results$m_plus_covered <- as.numeric(
    results$m_plus_lower <= results$m_plus_true & results$m_plus_true <= results$m_plus_upper
  )
}

summary_columns <- c(
  "p_zero", "m_plus",
  "p_zero_bias", "m_plus_bias",
  intersect(c("p_zero_se", "m_plus_se",
              "p_zero_covered", "m_plus_covered"),
            names(results))
)

summary_by_dtr <- aggregate(
  results[, summary_columns, drop = FALSE],
  by = list(
    scenario_id = results$scenario_id,
    scenario_label = results$scenario_label,
    scenario_family = results$scenario_family,
    ZI = results$ZI,
    corstr = results$corstr,
    rho = results$rho,
    DTR = results$DTR
  ),
  FUN = mean,
  na.rm = TRUE
)

empirical_sd_by_dtr <- aggregate(
  results[, c("p_zero", "m_plus"), drop = FALSE],
  by = list(
    scenario_id = results$scenario_id,
    scenario_label = results$scenario_label,
    scenario_family = results$scenario_family,
    ZI = results$ZI,
    corstr = results$corstr,
    rho = results$rho,
    DTR = results$DTR
  ),
  FUN = stats::sd,
  na.rm = TRUE
)
names(empirical_sd_by_dtr)[-(1:7)] <- paste0(names(empirical_sd_by_dtr)[-(1:7)], "_emp_sd")
summary_by_dtr <- merge(
  summary_by_dtr,
  empirical_sd_by_dtr,
  by = c("scenario_id", "scenario_label", "scenario_family", "ZI", "corstr", "rho", "DTR"),
  all.x = TRUE
)

cell_diagnostic_summary <- aggregate(
  cell_diagnostics[, c(
    "min_postsplit_positive_cell_n",
    "median_postsplit_positive_cell_n",
    "cells_below_10_positive",
    "cells_below_20_positive",
    "cells_below_30_positive",
    "cells_below_50_positive"
  ), drop = FALSE],
  by = list(
    scenario_id = cell_diagnostics$scenario_id,
    scenario_label = cell_diagnostics$scenario_label,
    scenario_family = cell_diagnostics$scenario_family,
    ZI = cell_diagnostics$ZI,
    corstr = cell_diagnostics$corstr,
    rho = cell_diagnostics$rho
  ),
  FUN = mean,
  na.rm = TRUE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(scenario_grid, file.path(output_dir, "scenario_grid_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(design_inputs, file.path(output_dir, "scenario_design_inputs_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(design_checks, file.path(output_dir, "scenario_design_checks_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(parameters, file.path(output_dir, "scenario_parameters_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(truth, file.path(output_dir, "truth_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(results, file.path(output_dir, "iteration_estimates_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(summary_by_dtr, file.path(output_dir, "bias_coverage_summary_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(cell_diagnostics, file.path(output_dir, "positive_cell_diagnostics_manuscript_primary.csv"), row.names = FALSE)
utils::write.csv(cell_diagnostic_summary, file.path(output_dir, "positive_cell_diagnostic_summary_manuscript_primary.csv"), row.names = FALSE)
if (!is.null(errors) && nrow(errors) > 0) {
  utils::write.csv(errors, file.path(output_dir, "iteration_errors_manuscript_primary.csv"), row.names = FALSE)
}

plot_conditional_bias_coverage(
  summary_by_dtr = summary_by_dtr,
  output_dir = output_dir,
  n = n,
  niter = niter,
  month = t_eval,
  message_fun = immediate_cat
)

immediate_cat("============================================================")
immediate_cat("All manuscript primary scenarios completed.")
immediate_cat("Output directory: ", output_dir)
if (!is.null(errors) && nrow(errors) > 0) {
  immediate_cat("Iteration errors recorded: ", nrow(errors))
} else {
  immediate_cat("Iteration errors recorded: 0")
}
immediate_cat("============================================================")
