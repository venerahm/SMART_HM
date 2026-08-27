#####################################################################
# Conditional-First Hurdle SMART Simulation: Stress Tests            #
#####################################################################
#
# Stress-test families:
#   1. Correct-HM zero-inflation/response-rate stress at n = 400.
#   2. Correct-HM zero-inflation/response-rate stress at n = 1000.
#   3. All-count zero prediction: compare observed zero proportions to zeros
#      implied by ordinary one-part Poisson and NB count models fit to all
#      generated counts, including zeros.
#   4. Response-model misspecification: generate response using
#      A1 + Y0 + BLDEPPOS, then fit the working response model A1 + Y0.
#
# These diagnostics use the same public design configuration as the primary
# simulation and write generated outputs under results/ by default.

suppressPackageStartupMessages({
  library(stats)
})

source("R/conditional_hurdle_model.R")
source("R/simulate_smart_data.R")
source("R/fit_g_computation.R")
source("R/plot_simulation_results.R")
source("R/simulation_helpers.R")
source("config/simulation_design.R")

# -----------------------
# User-configurable block
# -----------------------

# Each section runs and writes output independently. Toggle only the section(s)
# needed for the current run to avoid rerunning expensive diagnostics.
run_response_stress_n400 <- FALSE
run_response_stress_n1000 <- FALSE
run_zero_prediction <- TRUE
run_response_misspecification <- FALSE

# Section-specific grids and iteration counts.
response_sample_size_values <- c(400, 1000)
zero_prediction_sample_size_values <- c(400, 1000)
response_misspecification_sample_size_values <- c(400)

response_stress_niter <- as.integer(Sys.getenv("SMART_HM_RESPONSE_STRESS_NITER", "5000"))
zero_prediction_niter <- as.integer(Sys.getenv("SMART_HM_ZERO_PREDICTION_NITER", "1000"))
response_misspecification_niter <- as.integer(Sys.getenv("SMART_HM_MISSPECIFICATION_NITER", "5000"))
response_progress_every <- max(1, floor(response_stress_niter / 10))
zero_prediction_progress_every <- max(1, floor(zero_prediction_niter / 10))
response_misspecification_progress_every <- max(1, floor(response_misspecification_niter / 10))
times <- c(1, 2, 4)
spltime <- 1
t_eval <- max(times)
seed <- 20260706

stress_corstr <- "ar1"
stress_rho <- 0.5

stress_ZI_values <- c(15, 30, 60)

# Data-like response-misspecification settings. BLDEPPOS prevalence is based
# on the approximate observed baseline depression prevalence in the data
# application; the response OR controls misspecification strength.
response_misspecification_ZI_values <- c(37)
response_misspecification_corstr <- "ar1"
response_misspecification_rho <- 0.5
BLDEPPOS_prevalence <- 0.43
BLDEPPOS_response_log_or <- log(2)

stress_output_root <- file.path(
  Sys.getenv("SMART_HM_RESULTS_DIR", "results"),
  "stress_tests"
)
zero_prediction_output_dir <- file.path(
  stress_output_root,
  "Zero_Prediction",
  paste0("niter_", zero_prediction_niter)
)
response_misspecification_output_dir <- file.path(
  stress_output_root,
  "Response_Misspecification",
  paste0("niter_", response_misspecification_niter)
)

design_params <- simulation_design_params
Y0_ref <- simulation_Y0_ref

dtr_grid <- expand.grid(A1 = c(-1, 1), A2R = c(-1, 1), A2NR = c(-1, 1))

set.seed(seed)
Y0_truth <- gen_Y0_conditional(100000)

response_profiles <- data.frame(
  response_profile = c(
    "r0_0.5_r1_0.5",
    "r0_0.4_r1_0.6"
  ),
  response_prob_A1_minus = c(0.50, 0.40),
  response_prob_A1_plus = c(0.50, 0.60),
  response_profile_label = c(
    "r0 = 0.5, r1 = 0.5",
    "r0 = 0.4, r1 = 0.6"
  ),
  stringsAsFactors = FALSE
)

response_stress_grid <- merge(
  merge(
    data.frame(sample_size = response_sample_size_values),
    data.frame(ZI = stress_ZI_values),
    all = TRUE
  ),
  response_profiles,
  all = TRUE
)
response_stress_grid$scenario_label <- paste0(
  "R_", response_stress_grid$response_profile,
  "_n", response_stress_grid$sample_size,
  "_ZI", response_stress_grid$ZI
)
response_stress_grid$scenario_family <- "Response-rate and zero-inflation stress"
response_stress_grid$corstr <- stress_corstr
response_stress_grid$rho <- stress_rho
response_stress_grid$response_truth <- "A1 + Y0"
response_stress_grid$positive_count_truth <- "ZTP"
response_stress_grid$primary_fit <- "Correct HM"
response_stress_grid$why_include <- paste0(
  "Correct-HM performance with ",
  response_stress_grid$response_profile_label,
  " and ZI ", response_stress_grid$ZI, "%."
)
response_stress_grid$scenario_id <- paste0(
  response_stress_grid$scenario_label,
  "_", response_stress_grid$corstr,
  "_rho_", response_stress_grid$rho
)
response_stress_grid <- response_stress_grid[order(
  response_stress_grid$sample_size,
  response_stress_grid$ZI,
  response_stress_grid$response_profile
), ]
rownames(response_stress_grid) <- NULL
response_stress_grid_all <- response_stress_grid

zero_prediction_grid <- expand.grid(
  sample_size = zero_prediction_sample_size_values,
  ZI = stress_ZI_values
)
zero_prediction_grid$scenario_label <- paste0(
  "n", zero_prediction_grid$sample_size,
  "_ZI", zero_prediction_grid$ZI
)
zero_prediction_grid$scenario_family <- "All-count zero prediction"
zero_prediction_grid$corstr <- stress_corstr
zero_prediction_grid$rho <- stress_rho
zero_prediction_grid$response_prob_A1_minus <- 0.50
zero_prediction_grid$response_prob_A1_plus <- 0.50
zero_prediction_grid$response_profile <- "r0_0.5_r1_0.5"
zero_prediction_grid$response_profile_label <- "r0 = 0.5, r1 = 0.5"
zero_prediction_grid$response_truth <- "A1 + Y0"
zero_prediction_grid$positive_count_truth <- "ZTP"
zero_prediction_grid$comparator_fit <- "HM binary part, all-count Poisson, all-count NB"
zero_prediction_grid$why_include <- paste0(
  "Tests whether one-part all-count models recover observed zeros at n ",
  zero_prediction_grid$sample_size,
  " and ZI ", zero_prediction_grid$ZI, "%."
)
zero_prediction_grid <- zero_prediction_grid[, c(
  "scenario_label",
  "scenario_family",
  "sample_size", "ZI", "corstr", "rho",
  "response_prob_A1_minus", "response_prob_A1_plus",
  "response_profile", "response_profile_label",
  "response_truth", "positive_count_truth",
  "comparator_fit", "why_include"
)]
zero_prediction_grid$scenario_id <- paste0(
  zero_prediction_grid$scenario_label,
  "_", zero_prediction_grid$corstr,
  "_rho_", zero_prediction_grid$rho
)

response_misspecification_grid <- expand.grid(
  sample_size = response_misspecification_sample_size_values,
  ZI = response_misspecification_ZI_values
)
response_misspecification_grid$scenario_label <- paste0(
  "response_misspec_n",
  response_misspecification_grid$sample_size,
  "_ZI",
  response_misspecification_grid$ZI
)
response_misspecification_grid$scenario_family <- "Response model misspecification"
response_misspecification_grid$corstr <- response_misspecification_corstr
response_misspecification_grid$rho <- response_misspecification_rho
response_misspecification_grid$response_prob_A1_minus <- 0.50
response_misspecification_grid$response_prob_A1_plus <- 0.50
response_misspecification_grid$response_profile <- "truth_A1_Y0_BLDEPPOS_fit_A1_Y0"
response_misspecification_grid$response_profile_label <- "truth A1 + Y0 + BLDEPPOS; fit A1 + Y0"
response_misspecification_grid$response_truth <- "A1 + Y0 + BLDEPPOS"
response_misspecification_grid$response_fit <- "A1 + Y0"
response_misspecification_grid$positive_count_truth <- "ZTP"
response_misspecification_grid$primary_fit <- "HM with misspecified response model"
response_misspecification_grid$BLDEPPOS_prevalence <- BLDEPPOS_prevalence
response_misspecification_grid$BLDEPPOS_response_log_or <- BLDEPPOS_response_log_or
response_misspecification_grid$BLDEPPOS_response_or <- exp(BLDEPPOS_response_log_or)
response_misspecification_grid$why_include <- paste0(
  "Tests omitted baseline response predictor BLDEPPOS with prevalence ",
  round(BLDEPPOS_prevalence, 2),
  " and response OR ",
  round(exp(BLDEPPOS_response_log_or), 2),
  "."
)
response_misspecification_grid$scenario_id <- paste0(
  response_misspecification_grid$scenario_label,
  "_",
  response_misspecification_grid$corstr,
  "_rho_",
  response_misspecification_grid$rho
)

immediate_cat("============================================================")
immediate_cat("Conditional-first hurdle SMART stress tests")
immediate_cat("Run response stress n=400: ", run_response_stress_n400,
              " | run response stress n=1000: ", run_response_stress_n1000,
              " | niter: ", response_stress_niter,
              " | sample sizes: ", paste(response_sample_size_values, collapse = ", "))
immediate_cat("Run zero-prediction diagnostics: ", run_zero_prediction,
              " | niter: ", zero_prediction_niter,
              " | sample sizes: ", paste(zero_prediction_sample_size_values, collapse = ", "))
immediate_cat("Run response misspecification: ", run_response_misspecification,
              " | niter: ", response_misspecification_niter,
              " | sample sizes: ", paste(response_misspecification_sample_size_values, collapse = ", "))
immediate_cat("Response misspecification data-like settings: ZI ",
              paste(response_misspecification_ZI_values, collapse = ", "),
              "% | corstr: ", response_misspecification_corstr,
              " | rho: ", response_misspecification_rho)
immediate_cat("Stress ZI values: ", paste(stress_ZI_values, collapse = ", "), "%")
immediate_cat("times: ", paste(times, collapse = ", "), " | t_eval: ", t_eval)
immediate_cat("Stress corstr: ", stress_corstr, " | stress rho: ", stress_rho)
immediate_cat("Y0_ref for calibrated design inputs: ", Y0_ref)
immediate_cat("Zero-prediction output directory: ", zero_prediction_output_dir)
immediate_cat("Response misspecification output directory: ", response_misspecification_output_dir)
immediate_cat("============================================================")

#####################################################################
# 1. Correct-HM Performance Under Response And ZI Stress
#####################################################################

errors <- NULL
response_stress_sections <- data.frame(
  sample_size = c(400, 1000),
  run_section = c(run_response_stress_n400, run_response_stress_n1000)
)
response_stress_sections <- response_stress_sections[response_stress_sections$run_section, , drop = FALSE]

if (nrow(response_stress_sections) > 0) {

for (response_section_idx in seq_len(nrow(response_stress_sections))) {
response_section_n <- response_stress_sections$sample_size[response_section_idx]
response_output_dir <- file.path(
  stress_output_root,
  paste0("Response_Stress_n", response_section_n),
  paste0("niter_", response_stress_niter)
)
response_stress_grid <- response_stress_grid_all[
  response_stress_grid_all$sample_size == response_section_n,
  ,
  drop = FALSE
]
rownames(response_stress_grid) <- NULL

truth_list <- vector("list", nrow(response_stress_grid))
design_list <- vector("list", nrow(response_stress_grid))
design_check_list <- vector("list", nrow(response_stress_grid))
parameter_list <- vector("list", nrow(response_stress_grid))
estimate_list <- vector("list", nrow(response_stress_grid))
cell_diagnostic_list <- vector("list", nrow(response_stress_grid))
error_list <- vector("list", nrow(response_stress_grid))

immediate_cat("------------------------------------------------------------")
immediate_cat("Response-rate stress scenarios for n=", response_section_n, ": ",
              nrow(response_stress_grid))
immediate_cat("Response stress output directory: ", response_output_dir)

for (ss in seq_len(nrow(response_stress_grid))) {
  scenario <- response_stress_grid[ss, , drop = FALSE]
  scenario_start <- Sys.time()
  
  immediate_cat("------------------------------------------------------------")
  immediate_cat("Response stress scenario ", ss, "/", nrow(response_stress_grid), ": ",
                scenario$scenario_id)
  immediate_cat("ZI: ", scenario$ZI, "% | r0: ", scenario$response_prob_A1_minus,
                " | r1: ", scenario$response_prob_A1_plus)
  
  design_params_s <- design_params
  design_params_s$response_prob_A1_minus <- scenario$response_prob_A1_minus
  design_params_s$response_prob_A1_plus <- scenario$response_prob_A1_plus
  
  params <- conditional_params_from_design(
    design_params = design_params_s,
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
  truth_s$sample_size <- scenario$sample_size
  truth_s$ZI <- scenario$ZI
  truth_s$corstr <- scenario$corstr
  truth_s$rho <- scenario$rho
  truth_s$response_profile <- scenario$response_profile
  truth_s$response_prob_A1_minus <- scenario$response_prob_A1_minus
  truth_s$response_prob_A1_plus <- scenario$response_prob_A1_plus
  truth_list[[ss]] <- truth_s
  
  design_s <- conditional_design_table(design_params_s)
  design_s$scenario_id <- scenario$scenario_id
  design_s$scenario_label <- scenario$scenario_label
  design_s$scenario_family <- scenario$scenario_family
  design_s$sample_size <- scenario$sample_size
  design_s$ZI <- scenario$ZI
  design_s$corstr <- scenario$corstr
  design_s$rho <- scenario$rho
  design_s$response_profile <- scenario$response_profile
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
  design_check_s$sample_size <- scenario$sample_size
  design_check_s$ZI <- scenario$ZI
  design_check_s$corstr <- scenario$corstr
  design_check_s$rho <- scenario$rho
  design_check_s$response_profile <- scenario$response_profile
  design_check_s$zero_probability <- scenario$ZI / 100
  design_check_s$Y0_ref <- Y0_ref
  design_check_s$time_ref <- t_eval
  design_check_list[[ss]] <- design_check_s
  
  params_s <- conditional_parameter_table(params)
  params_s$scenario_id <- scenario$scenario_id
  params_s$scenario_label <- scenario$scenario_label
  params_s$scenario_family <- scenario$scenario_family
  params_s$sample_size <- scenario$sample_size
  params_s$ZI <- scenario$ZI
  params_s$corstr <- scenario$corstr
  params_s$rho <- scenario$rho
  params_s$response_profile <- scenario$response_profile
  parameter_list[[ss]] <- params_s
  
  iteration_estimates <- vector("list", response_stress_niter)
  iteration_cell_diagnostics <- vector("list", response_stress_niter)
  iteration_errors <- vector("list", response_stress_niter)
  
  for (iter in seq_len(response_stress_niter)) {
    if (iter == 1 || iter == response_stress_niter || iter %% response_progress_every == 0) {
      immediate_cat(sprintf(
        "  Response stress %d/%d (%s) - iteration %d/%d",
        ss, nrow(response_stress_grid), scenario$scenario_label, iter, response_stress_niter
      ))
    }
    
    set.seed(seed + ss * 100000 + iter)
    generated <- generateSMART_gee_conditional(
      n = scenario$sample_size,
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
    cell_diag$sample_size <- scenario$sample_size
    cell_diag$ZI <- scenario$ZI
    cell_diag$corstr <- scenario$corstr
    cell_diag$rho <- scenario$rho
    cell_diag$response_profile <- scenario$response_profile
    iteration_cell_diagnostics[[iter]] <- cell_diag
    
    iter_result <- fit_conditional_iteration_for_wrapper(
      generated = generated,
      spltime = spltime,
      dtr_grid = dtr_grid,
      t_eval = t_eval,
      use_delta_se = TRUE,
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
  immediate_cat("Response stress scenario ", ss, "/", nrow(response_stress_grid),
                " completed in ", scenario_elapsed, " minutes.")
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
  "sample_size", "ZI", "corstr", "rho", "response_profile", "DTR",
  "p_zero_true", "m_plus_true"
), drop = FALSE]

response_results <- merge(
  estimates,
  truth_eval,
  by = c("scenario_id", "scenario_label", "scenario_family",
         "sample_size", "ZI", "corstr", "rho", "response_profile", "DTR"),
  all.x = TRUE
)
response_results$p_zero_bias <- response_results$p_zero - response_results$p_zero_true
response_results$m_plus_bias <- response_results$m_plus - response_results$m_plus_true
response_results$p_zero_lower <- response_results$p_zero - 1.96 * response_results$p_zero_se
response_results$p_zero_upper <- response_results$p_zero + 1.96 * response_results$p_zero_se
response_results$m_plus_lower <- response_results$m_plus - 1.96 * response_results$m_plus_se
response_results$m_plus_upper <- response_results$m_plus + 1.96 * response_results$m_plus_se
response_results$p_zero_covered <- as.numeric(
  response_results$p_zero_lower <= response_results$p_zero_true &
    response_results$p_zero_true <= response_results$p_zero_upper
)
response_results$m_plus_covered <- as.numeric(
  response_results$m_plus_lower <= response_results$m_plus_true &
    response_results$m_plus_true <= response_results$m_plus_upper
)

response_summary_columns <- c(
  "p_zero", "m_plus",
  "p_zero_bias", "m_plus_bias",
  "p_zero_se", "m_plus_se",
  "p_zero_covered", "m_plus_covered"
)
response_summary_by_dtr <- stats::aggregate(
  response_results[, response_summary_columns, drop = FALSE],
  by = list(
    scenario_id = response_results$scenario_id,
    scenario_label = response_results$scenario_label,
    scenario_family = response_results$scenario_family,
    sample_size = response_results$sample_size,
    ZI = response_results$ZI,
    corstr = response_results$corstr,
    rho = response_results$rho,
    response_profile = response_results$response_profile,
    DTR = response_results$DTR
  ),
  FUN = mean,
  na.rm = TRUE
)

response_empirical_sd_by_dtr <- stats::aggregate(
  response_results[, c("p_zero", "m_plus"), drop = FALSE],
  by = list(
    scenario_id = response_results$scenario_id,
    scenario_label = response_results$scenario_label,
    scenario_family = response_results$scenario_family,
    sample_size = response_results$sample_size,
    ZI = response_results$ZI,
    corstr = response_results$corstr,
    rho = response_results$rho,
    response_profile = response_results$response_profile,
    DTR = response_results$DTR
  ),
  FUN = stats::sd,
  na.rm = TRUE
)
names(response_empirical_sd_by_dtr)[-(1:9)] <- paste0(
  names(response_empirical_sd_by_dtr)[-(1:9)], "_emp_sd"
)
response_summary_by_dtr <- merge(
  response_summary_by_dtr,
  response_empirical_sd_by_dtr,
  by = c("scenario_id", "scenario_label", "scenario_family",
         "sample_size", "ZI", "corstr", "rho", "response_profile", "DTR"),
  all.x = TRUE
)

response_cell_diagnostic_summary <- stats::aggregate(
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
    sample_size = cell_diagnostics$sample_size,
    ZI = cell_diagnostics$ZI,
    corstr = cell_diagnostics$corstr,
    rho = cell_diagnostics$rho,
    response_profile = cell_diagnostics$response_profile
  ),
  FUN = mean,
  na.rm = TRUE
)

dir.create(response_output_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(response_stress_grid, file.path(response_output_dir, "scenario_grid_response_stress.csv"), row.names = FALSE)
utils::write.csv(design_inputs, file.path(response_output_dir, "scenario_design_inputs_response_stress.csv"), row.names = FALSE)
utils::write.csv(design_checks, file.path(response_output_dir, "scenario_design_checks_response_stress.csv"), row.names = FALSE)
utils::write.csv(parameters, file.path(response_output_dir, "scenario_parameters_response_stress.csv"), row.names = FALSE)
utils::write.csv(truth, file.path(response_output_dir, "truth_response_stress.csv"), row.names = FALSE)
utils::write.csv(response_results, file.path(response_output_dir, "iteration_estimates_response_stress.csv"), row.names = FALSE)
utils::write.csv(response_summary_by_dtr, file.path(response_output_dir, "bias_coverage_summary_response_stress.csv"), row.names = FALSE)
utils::write.csv(cell_diagnostics, file.path(response_output_dir, "positive_cell_diagnostics_response_stress.csv"), row.names = FALSE)
utils::write.csv(response_cell_diagnostic_summary, file.path(response_output_dir, "positive_cell_diagnostic_summary_response_stress.csv"), row.names = FALSE)
if (!is.null(errors) && nrow(errors) > 0) {
  utils::write.csv(errors, file.path(response_output_dir, "iteration_errors_response_stress.csv"), row.names = FALSE)
}

plot_response_stress_coverage(
  summary_by_dtr = response_summary_by_dtr,
  output_dir = response_output_dir,
  niter = response_stress_niter,
  month = t_eval,
  message_fun = immediate_cat
)

immediate_cat("Completed response stress section for n=", response_section_n)
}

} else {
  immediate_cat("Skipping response-rate/zero-inflation bias and coverage stress tests.")
}

#####################################################################
# 2. All-Count Zero Prediction Stress
#####################################################################

zero_errors <- NULL
if (run_zero_prediction) {

diagnostic_list <- vector("list", nrow(zero_prediction_grid))
cell_diagnostic_list <- vector("list", nrow(zero_prediction_grid))
zero_design_list <- vector("list", nrow(zero_prediction_grid))
zero_design_check_list <- vector("list", nrow(zero_prediction_grid))
zero_parameter_list <- vector("list", nrow(zero_prediction_grid))
zero_error_list <- vector("list", nrow(zero_prediction_grid))

immediate_cat("------------------------------------------------------------")
immediate_cat("All-count zero-prediction scenarios: ", nrow(zero_prediction_grid))

for (ss in seq_len(nrow(zero_prediction_grid))) {
  scenario <- zero_prediction_grid[ss, , drop = FALSE]
  scenario_start <- Sys.time()
  
  immediate_cat("------------------------------------------------------------")
  immediate_cat("Zero-prediction scenario ", ss, "/", nrow(zero_prediction_grid), ": ",
                scenario$scenario_id)
  immediate_cat("ZI: ", scenario$ZI, "% | corstr: ", scenario$corstr,
                " | rho: ", scenario$rho)
  
  design_params_s <- design_params
  design_params_s$response_prob_A1_minus <- scenario$response_prob_A1_minus
  design_params_s$response_prob_A1_plus <- scenario$response_prob_A1_plus
  
  params <- conditional_params_from_design(
    design_params = design_params_s,
    zero_probability = scenario$ZI / 100,
    Y0_ref = Y0_ref,
    time_ref = t_eval,
    spltime = spltime
  )
  
  design_s <- conditional_design_table(design_params_s)
  design_s$scenario_id <- scenario$scenario_id
  design_s$scenario_label <- scenario$scenario_label
  design_s$scenario_family <- scenario$scenario_family
  design_s$sample_size <- scenario$sample_size
  design_s$ZI <- scenario$ZI
  design_s$corstr <- scenario$corstr
  design_s$rho <- scenario$rho
  design_s$zero_probability <- scenario$ZI / 100
  design_s$Y0_ref <- Y0_ref
  design_s$time_ref <- t_eval
  zero_design_list[[ss]] <- design_s
  
  design_check_s <- conditional_design_check_table(
    params = params,
    Y0_ref = Y0_ref,
    time_ref = t_eval,
    spltime = spltime
  )
  design_check_s$scenario_id <- scenario$scenario_id
  design_check_s$scenario_label <- scenario$scenario_label
  design_check_s$scenario_family <- scenario$scenario_family
  design_check_s$sample_size <- scenario$sample_size
  design_check_s$ZI <- scenario$ZI
  design_check_s$corstr <- scenario$corstr
  design_check_s$rho <- scenario$rho
  design_check_s$zero_probability <- scenario$ZI / 100
  design_check_s$Y0_ref <- Y0_ref
  design_check_s$time_ref <- t_eval
  zero_design_check_list[[ss]] <- design_check_s
  
  params_s <- conditional_parameter_table(params)
  params_s$scenario_id <- scenario$scenario_id
  params_s$scenario_label <- scenario$scenario_label
  params_s$scenario_family <- scenario$scenario_family
  params_s$sample_size <- scenario$sample_size
  params_s$ZI <- scenario$ZI
  params_s$corstr <- scenario$corstr
  params_s$rho <- scenario$rho
  zero_parameter_list[[ss]] <- params_s
  
  iteration_diagnostics <- vector("list", zero_prediction_niter)
  iteration_cell_diagnostics <- vector("list", zero_prediction_niter)
  iteration_errors <- vector("list", zero_prediction_niter)
  
  for (iter in seq_len(zero_prediction_niter)) {
    if (iter == 1 || iter == zero_prediction_niter ||
        iter %% zero_prediction_progress_every == 0) {
      immediate_cat(sprintf(
        "  Zero-prediction %d/%d (%s: ZI=%s%%, corstr=%s, rho=%s) - iteration %d/%d",
        ss, nrow(zero_prediction_grid), scenario$scenario_label,
        scenario$ZI, scenario$corstr, scenario$rho, iter, zero_prediction_niter
      ))
    }
    
    set.seed(seed + 10000000 + ss * 100000 + iter)
    generated <- generateSMART_gee_conditional(
      n = scenario$sample_size,
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
    cell_diag$sample_size <- scenario$sample_size
    cell_diag$ZI <- scenario$ZI
    cell_diag$corstr <- scenario$corstr
    cell_diag$rho <- scenario$rho
    iteration_cell_diagnostics[[iter]] <- cell_diag
    
    iter_result <- fit_positive_count_zero_diagnostics(
      long_data = generated$long_data,
      spltime = spltime,
      scenario = scenario,
      iter = iter
    )
    iteration_diagnostics[[iter]] <- iter_result$summary
    iteration_errors[[iter]] <- iter_result$error
  }
  
  diagnostic_list[[ss]] <- do.call(rbind, iteration_diagnostics)
  cell_diagnostic_list[[ss]] <- do.call(rbind, iteration_cell_diagnostics)
  zero_error_list[[ss]] <- do.call(rbind, iteration_errors)
  
  scenario_elapsed <- round(as.numeric(difftime(Sys.time(), scenario_start, units = "mins")), 2)
  immediate_cat("Zero-prediction scenario ", ss, "/", nrow(zero_prediction_grid),
                " completed in ", scenario_elapsed, " minutes.")
}

zero_diagnostics <- do.call(rbind, diagnostic_list)
zero_cell_diagnostics <- do.call(rbind, cell_diagnostic_list)
zero_design_inputs <- do.call(rbind, zero_design_list)
zero_design_checks <- do.call(rbind, zero_design_check_list)
zero_parameters <- do.call(rbind, zero_parameter_list)
zero_errors <- do.call(rbind, zero_error_list)

zero_diagnostic_summary <- stats::aggregate(
  zero_diagnostics[, c(
    "observed_zero",
    "predicted_zero_hm",
    "predicted_zero_all_poisson",
    "predicted_zero_all_nb",
    "excess_zero_ratio_all_poisson",
    "excess_zero_ratio_all_nb",
    "n",
    "n_positive"
  ), drop = FALSE],
  by = list(
    scenario_id = zero_diagnostics$scenario_id,
    scenario_label = zero_diagnostics$scenario_label,
    scenario_family = zero_diagnostics$scenario_family,
    sample_size = zero_diagnostics$sample_size,
    ZI = zero_diagnostics$ZI,
    corstr = zero_diagnostics$corstr,
    rho = zero_diagnostics$rho,
    time = zero_diagnostics$time
  ),
  FUN = mean,
  na.rm = TRUE
)

zero_cell_diagnostic_summary <- stats::aggregate(
  zero_cell_diagnostics[, c(
    "min_postsplit_positive_cell_n",
    "median_postsplit_positive_cell_n",
    "cells_below_10_positive",
    "cells_below_20_positive",
    "cells_below_30_positive",
    "cells_below_50_positive"
  ), drop = FALSE],
  by = list(
    scenario_id = zero_cell_diagnostics$scenario_id,
    scenario_label = zero_cell_diagnostics$scenario_label,
    scenario_family = zero_cell_diagnostics$scenario_family,
    sample_size = zero_cell_diagnostics$sample_size,
    ZI = zero_cell_diagnostics$ZI,
    corstr = zero_cell_diagnostics$corstr,
    rho = zero_cell_diagnostics$rho
  ),
  FUN = mean,
  na.rm = TRUE
)

} else {
  immediate_cat("Skipping all-count zero-prediction diagnostics.")
}

#####################################################################
# 3. Response Model Misspecification: Omitted BLDEPPOS
#####################################################################

misspec_errors <- NULL
if (run_response_misspecification) {

misspec_truth_list <- vector("list", nrow(response_misspecification_grid))
misspec_design_list <- vector("list", nrow(response_misspecification_grid))
misspec_design_check_list <- vector("list", nrow(response_misspecification_grid))
misspec_parameter_list <- vector("list", nrow(response_misspecification_grid))
misspec_estimate_list <- vector("list", nrow(response_misspecification_grid))
misspec_cell_diagnostic_list <- vector("list", nrow(response_misspecification_grid))
misspec_error_list <- vector("list", nrow(response_misspecification_grid))

immediate_cat("------------------------------------------------------------")
immediate_cat("Response misspecification scenarios: ", nrow(response_misspecification_grid))
immediate_cat("Truth response model: A1 + Y0 + BLDEPPOS")
immediate_cat("Fitted response model: A1 + Y0")
immediate_cat("BLDEPPOS prevalence: ", BLDEPPOS_prevalence,
              " | response OR: ", exp(BLDEPPOS_response_log_or))

set.seed(seed + 30000000)
BLDEPPOS_truth <- stats::rbinom(length(Y0_truth), 1, BLDEPPOS_prevalence)
BLDEPPOS_truth_data <- data.frame(BLDEPPOS = BLDEPPOS_truth)

for (ss in seq_len(nrow(response_misspecification_grid))) {
  scenario <- response_misspecification_grid[ss, , drop = FALSE]
  scenario_start <- Sys.time()
  
  immediate_cat("------------------------------------------------------------")
  immediate_cat("Response misspecification scenario ", ss, "/",
                nrow(response_misspecification_grid), ": ",
                scenario$scenario_id)
  immediate_cat("ZI: ", scenario$ZI, "% | n: ", scenario$sample_size,
                " | BLDEPPOS OR: ", scenario$BLDEPPOS_response_or)
  
  design_params_s <- design_params
  design_params_s$response_prob_A1_minus <- scenario$response_prob_A1_minus
  design_params_s$response_prob_A1_plus <- scenario$response_prob_A1_plus
  
  params <- conditional_params_from_design(
    design_params = design_params_s,
    zero_probability = scenario$ZI / 100,
    Y0_ref = Y0_ref,
    time_ref = t_eval,
    spltime = spltime
  )
  params$alpha_response["(Intercept)"] <-
    params$alpha_response["(Intercept)"] -
    BLDEPPOS_response_log_or * BLDEPPOS_prevalence
  params$alpha_response["BLDEPPOS"] <- BLDEPPOS_response_log_or
  
  truth_s <- conditional_dtr_truth(
    params = params,
    dtr_grid = dtr_grid,
    Y0_vec = Y0_truth,
    times = times,
    spltime = spltime,
    baseline_response_covariates = BLDEPPOS_truth_data
  )
  truth_s$scenario_id <- scenario$scenario_id
  truth_s$scenario_label <- scenario$scenario_label
  truth_s$scenario_family <- scenario$scenario_family
  truth_s$sample_size <- scenario$sample_size
  truth_s$ZI <- scenario$ZI
  truth_s$corstr <- scenario$corstr
  truth_s$rho <- scenario$rho
  truth_s$response_profile <- scenario$response_profile
  truth_s$response_prob_A1_minus <- scenario$response_prob_A1_minus
  truth_s$response_prob_A1_plus <- scenario$response_prob_A1_plus
  truth_s$BLDEPPOS_prevalence <- BLDEPPOS_prevalence
  truth_s$BLDEPPOS_response_log_or <- BLDEPPOS_response_log_or
  misspec_truth_list[[ss]] <- truth_s
  
  design_s <- conditional_design_table(design_params_s)
  design_s$scenario_id <- scenario$scenario_id
  design_s$scenario_label <- scenario$scenario_label
  design_s$scenario_family <- scenario$scenario_family
  design_s$sample_size <- scenario$sample_size
  design_s$ZI <- scenario$ZI
  design_s$corstr <- scenario$corstr
  design_s$rho <- scenario$rho
  design_s$response_profile <- scenario$response_profile
  design_s$zero_probability <- scenario$ZI / 100
  design_s$Y0_ref <- Y0_ref
  design_s$time_ref <- t_eval
  design_s$BLDEPPOS_prevalence <- BLDEPPOS_prevalence
  design_s$BLDEPPOS_response_log_or <- BLDEPPOS_response_log_or
  misspec_design_list[[ss]] <- design_s
  
  design_check_s <- conditional_design_check_table(
    params = params,
    Y0_ref = Y0_ref,
    time_ref = t_eval,
    spltime = spltime,
    baseline_response_covariates_ref = data.frame(BLDEPPOS = BLDEPPOS_prevalence)
  )
  design_check_s$scenario_id <- scenario$scenario_id
  design_check_s$scenario_label <- scenario$scenario_label
  design_check_s$scenario_family <- scenario$scenario_family
  design_check_s$sample_size <- scenario$sample_size
  design_check_s$ZI <- scenario$ZI
  design_check_s$corstr <- scenario$corstr
  design_check_s$rho <- scenario$rho
  design_check_s$response_profile <- scenario$response_profile
  design_check_s$zero_probability <- scenario$ZI / 100
  design_check_s$Y0_ref <- Y0_ref
  design_check_s$time_ref <- t_eval
  design_check_s$BLDEPPOS_prevalence <- BLDEPPOS_prevalence
  design_check_s$BLDEPPOS_response_log_or <- BLDEPPOS_response_log_or
  misspec_design_check_list[[ss]] <- design_check_s
  
  params_s <- conditional_parameter_table(params)
  params_s$scenario_id <- scenario$scenario_id
  params_s$scenario_label <- scenario$scenario_label
  params_s$scenario_family <- scenario$scenario_family
  params_s$sample_size <- scenario$sample_size
  params_s$ZI <- scenario$ZI
  params_s$corstr <- scenario$corstr
  params_s$rho <- scenario$rho
  params_s$response_profile <- scenario$response_profile
  misspec_parameter_list[[ss]] <- params_s
  
  iteration_estimates <- vector("list", response_misspecification_niter)
  iteration_cell_diagnostics <- vector("list", response_misspecification_niter)
  iteration_errors <- vector("list", response_misspecification_niter)
  
  for (iter in seq_len(response_misspecification_niter)) {
    if (iter == 1 || iter == response_misspecification_niter ||
        iter %% response_misspecification_progress_every == 0) {
      immediate_cat(sprintf(
        "  Response misspecification %d/%d (%s) - iteration %d/%d",
        ss, nrow(response_misspecification_grid), scenario$scenario_label,
        iter, response_misspecification_niter
      ))
    }
    
    set.seed(seed + 40000000 + ss * 100000 + iter)
    BLDEPPOS_iter <- stats::rbinom(scenario$sample_size, 1, BLDEPPOS_prevalence)
    generated <- generateSMART_gee_conditional(
      n = scenario$sample_size,
      times = times,
      spltime = spltime,
      params = params,
      corstr = scenario$corstr,
      rho = scenario$rho,
      baseline_response_covariates = data.frame(BLDEPPOS = BLDEPPOS_iter)
    )
    
    cell_diag <- positive_cell_summary(generated$long_data, spltime = spltime)
    cell_diag$iteration <- iter
    cell_diag$scenario_id <- scenario$scenario_id
    cell_diag$scenario_label <- scenario$scenario_label
    cell_diag$scenario_family <- scenario$scenario_family
    cell_diag$sample_size <- scenario$sample_size
    cell_diag$ZI <- scenario$ZI
    cell_diag$corstr <- scenario$corstr
    cell_diag$rho <- scenario$rho
    cell_diag$response_profile <- scenario$response_profile
    cell_diag$BLDEPPOS_prevalence_observed <- mean(BLDEPPOS_iter)
    iteration_cell_diagnostics[[iter]] <- cell_diag
    
    iter_result <- fit_conditional_iteration_for_wrapper(
      generated = generated,
      spltime = spltime,
      dtr_grid = dtr_grid,
      t_eval = t_eval,
      use_delta_se = TRUE,
      scenario = scenario,
      iter = iter
    )
    iteration_estimates[[iter]] <- iter_result$estimate
    iteration_errors[[iter]] <- iter_result$error
  }
  
  misspec_estimate_list[[ss]] <- do.call(rbind, iteration_estimates)
  misspec_cell_diagnostic_list[[ss]] <- do.call(rbind, iteration_cell_diagnostics)
  misspec_error_list[[ss]] <- do.call(rbind, iteration_errors)
  
  scenario_elapsed <- round(as.numeric(difftime(Sys.time(), scenario_start, units = "mins")), 2)
  immediate_cat("Response misspecification scenario ", ss, "/",
                nrow(response_misspecification_grid),
                " completed in ", scenario_elapsed, " minutes.")
}

misspec_truth <- do.call(rbind, misspec_truth_list)
misspec_design_inputs <- do.call(rbind, misspec_design_list)
misspec_design_checks <- do.call(rbind, misspec_design_check_list)
misspec_parameters <- do.call(rbind, misspec_parameter_list)
misspec_estimates <- do.call(rbind, misspec_estimate_list)
misspec_cell_diagnostics <- do.call(rbind, misspec_cell_diagnostic_list)
misspec_errors <- do.call(rbind, misspec_error_list)

misspec_truth_eval <- misspec_truth[misspec_truth$time == t_eval, c(
  "scenario_id", "scenario_label", "scenario_family",
  "sample_size", "ZI", "corstr", "rho", "response_profile", "DTR",
  "p_zero_true", "m_plus_true"
), drop = FALSE]

misspec_results <- merge(
  misspec_estimates,
  misspec_truth_eval,
  by = c("scenario_id", "scenario_label", "scenario_family",
         "sample_size", "ZI", "corstr", "rho", "response_profile", "DTR"),
  all.x = TRUE
)
misspec_results$p_zero_bias <- misspec_results$p_zero - misspec_results$p_zero_true
misspec_results$m_plus_bias <- misspec_results$m_plus - misspec_results$m_plus_true
misspec_results$p_zero_lower <- misspec_results$p_zero - 1.96 * misspec_results$p_zero_se
misspec_results$p_zero_upper <- misspec_results$p_zero + 1.96 * misspec_results$p_zero_se
misspec_results$m_plus_lower <- misspec_results$m_plus - 1.96 * misspec_results$m_plus_se
misspec_results$m_plus_upper <- misspec_results$m_plus + 1.96 * misspec_results$m_plus_se
misspec_results$p_zero_covered <- as.numeric(
  misspec_results$p_zero_lower <= misspec_results$p_zero_true &
    misspec_results$p_zero_true <= misspec_results$p_zero_upper
)
misspec_results$m_plus_covered <- as.numeric(
  misspec_results$m_plus_lower <= misspec_results$m_plus_true &
    misspec_results$m_plus_true <= misspec_results$m_plus_upper
)

misspec_summary_columns <- c(
  "p_zero", "m_plus",
  "p_zero_bias", "m_plus_bias",
  "p_zero_se", "m_plus_se",
  "p_zero_covered", "m_plus_covered"
)
misspec_summary_by_dtr <- stats::aggregate(
  misspec_results[, misspec_summary_columns, drop = FALSE],
  by = list(
    scenario_id = misspec_results$scenario_id,
    scenario_label = misspec_results$scenario_label,
    scenario_family = misspec_results$scenario_family,
    sample_size = misspec_results$sample_size,
    ZI = misspec_results$ZI,
    corstr = misspec_results$corstr,
    rho = misspec_results$rho,
    response_profile = misspec_results$response_profile,
    DTR = misspec_results$DTR
  ),
  FUN = mean,
  na.rm = TRUE
)

misspec_empirical_sd_by_dtr <- stats::aggregate(
  misspec_results[, c("p_zero", "m_plus"), drop = FALSE],
  by = list(
    scenario_id = misspec_results$scenario_id,
    scenario_label = misspec_results$scenario_label,
    scenario_family = misspec_results$scenario_family,
    sample_size = misspec_results$sample_size,
    ZI = misspec_results$ZI,
    corstr = misspec_results$corstr,
    rho = misspec_results$rho,
    response_profile = misspec_results$response_profile,
    DTR = misspec_results$DTR
  ),
  FUN = stats::sd,
  na.rm = TRUE
)
names(misspec_empirical_sd_by_dtr)[-(1:9)] <- paste0(
  names(misspec_empirical_sd_by_dtr)[-(1:9)], "_emp_sd"
)
misspec_summary_by_dtr <- merge(
  misspec_summary_by_dtr,
  misspec_empirical_sd_by_dtr,
  by = c("scenario_id", "scenario_label", "scenario_family",
         "sample_size", "ZI", "corstr", "rho", "response_profile", "DTR"),
  all.x = TRUE
)

misspec_cell_diagnostic_summary <- stats::aggregate(
  misspec_cell_diagnostics[, c(
    "min_postsplit_positive_cell_n",
    "median_postsplit_positive_cell_n",
    "cells_below_10_positive",
    "cells_below_20_positive",
    "cells_below_30_positive",
    "cells_below_50_positive",
    "BLDEPPOS_prevalence_observed"
  ), drop = FALSE],
  by = list(
    scenario_id = misspec_cell_diagnostics$scenario_id,
    scenario_label = misspec_cell_diagnostics$scenario_label,
    scenario_family = misspec_cell_diagnostics$scenario_family,
    sample_size = misspec_cell_diagnostics$sample_size,
    ZI = misspec_cell_diagnostics$ZI,
    corstr = misspec_cell_diagnostics$corstr,
    rho = misspec_cell_diagnostics$rho,
    response_profile = misspec_cell_diagnostics$response_profile
  ),
  FUN = mean,
  na.rm = TRUE
)

dir.create(response_misspecification_output_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(response_misspecification_grid, file.path(response_misspecification_output_dir, "scenario_grid_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_design_inputs, file.path(response_misspecification_output_dir, "scenario_design_inputs_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_design_checks, file.path(response_misspecification_output_dir, "scenario_design_checks_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_parameters, file.path(response_misspecification_output_dir, "scenario_parameters_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_truth, file.path(response_misspecification_output_dir, "truth_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_results, file.path(response_misspecification_output_dir, "iteration_estimates_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_summary_by_dtr, file.path(response_misspecification_output_dir, "bias_coverage_summary_response_stress.csv"), row.names = FALSE)
utils::write.csv(misspec_summary_by_dtr, file.path(response_misspecification_output_dir, "bias_coverage_summary_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_cell_diagnostics, file.path(response_misspecification_output_dir, "positive_cell_diagnostics_response_misspecification.csv"), row.names = FALSE)
utils::write.csv(misspec_cell_diagnostic_summary, file.path(response_misspecification_output_dir, "positive_cell_diagnostic_summary_response_misspecification.csv"), row.names = FALSE)
if (!is.null(misspec_errors) && nrow(misspec_errors) > 0) {
  utils::write.csv(misspec_errors, file.path(response_misspecification_output_dir, "iteration_errors_response_misspecification.csv"), row.names = FALSE)
}

plot_response_stress_coverage(
  summary_by_dtr = misspec_summary_by_dtr,
  output_dir = response_misspecification_output_dir,
  niter = response_misspecification_niter,
  month = t_eval,
  message_fun = immediate_cat
)

immediate_cat("Response misspecification output directory: ",
              response_misspecification_output_dir)

} else {
  immediate_cat("Skipping response model misspecification scenarios.")
}

#####################################################################
# Write outputs
#####################################################################

if (run_zero_prediction) {
  dir.create(zero_prediction_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  utils::write.csv(zero_prediction_grid, file.path(zero_prediction_output_dir, "scenario_grid_zero_prediction.csv"), row.names = FALSE)
  utils::write.csv(zero_design_inputs, file.path(zero_prediction_output_dir, "scenario_design_inputs_zero_prediction.csv"), row.names = FALSE)
  utils::write.csv(zero_design_checks, file.path(zero_prediction_output_dir, "scenario_design_checks_zero_prediction.csv"), row.names = FALSE)
  utils::write.csv(zero_parameters, file.path(zero_prediction_output_dir, "scenario_parameters_zero_prediction.csv"), row.names = FALSE)
  utils::write.csv(zero_diagnostics, file.path(zero_prediction_output_dir, "iteration_zero_prediction_diagnostics.csv"), row.names = FALSE)
  utils::write.csv(zero_diagnostic_summary, file.path(zero_prediction_output_dir, "zero_prediction_summary_stress_tests.csv"), row.names = FALSE)
  utils::write.csv(zero_cell_diagnostics, file.path(zero_prediction_output_dir, "positive_cell_diagnostics_zero_prediction.csv"), row.names = FALSE)
  utils::write.csv(zero_cell_diagnostic_summary, file.path(zero_prediction_output_dir, "positive_cell_diagnostic_summary_zero_prediction.csv"), row.names = FALSE)
  if (!is.null(zero_errors) && nrow(zero_errors) > 0) {
    utils::write.csv(zero_errors, file.path(zero_prediction_output_dir, "iteration_errors_zero_prediction.csv"), row.names = FALSE)
  }
  
  plot_zero_prediction_diagnostics(
    zero_summary = zero_diagnostic_summary,
    output_dir = zero_prediction_output_dir,
    month = t_eval,
    message_fun = immediate_cat
  )
}

immediate_cat("============================================================")
immediate_cat("All stress-test scenarios completed.")
if (run_response_stress_n400) {
  immediate_cat("Response n=400 output directory: ",
                file.path(stress_output_root, "Response_Stress_n400",
                          paste0("niter_", response_stress_niter)))
}
if (run_response_stress_n1000) {
  immediate_cat("Response n=1000 output directory: ",
                file.path(stress_output_root, "Response_Stress_n1000",
                          paste0("niter_", response_stress_niter)))
}
if (run_zero_prediction) immediate_cat("Zero-prediction output directory: ", zero_prediction_output_dir)
if (run_response_misspecification) {
  immediate_cat("Response misspecification output directory: ",
                response_misspecification_output_dir)
}
immediate_cat("Response stress errors: ", ifelse(is.null(errors), 0, nrow(errors)))
immediate_cat("Zero-prediction errors: ", ifelse(is.null(zero_errors), 0, nrow(zero_errors)))
immediate_cat("Response misspecification errors: ",
              ifelse(is.null(misspec_errors), 0, nrow(misspec_errors)))
immediate_cat("============================================================")
