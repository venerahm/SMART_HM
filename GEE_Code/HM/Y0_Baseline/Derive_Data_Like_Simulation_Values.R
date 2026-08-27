#####################################################################
# Derive Data-Like Simulation Values For Conditional Hurdle SMARTs   #
#####################################################################
#
# This script computes aggregate, application-calibrated quantities for the
# conditional/exact Y0-baseline simulation design. It does not write subject-
# level data; outputs are summary CSV/RDS files only.
#
# Optional environment overrides:
#   MCOACH_WEEKLY_PATH
#   MCOACH_OUTCOME_PATH
#
# Main outputs:
#   Data_Results_Cond_New_Surrogate/Data_Like_Calibration/scenario_calibration_values.csv
#   Data_Results_Cond_New_Surrogate/Data_Like_Calibration/data_like_design_params.csv
#   Data_Results_Cond_New_Surrogate/Data_Like_Calibration/rho_recommendations.csv
#   Data_Results_Cond_New_Surrogate/Data_Like_Calibration/one_part_zero_diagnostics.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

source("GEE_Code/HM/Y0_Baseline/functions_HM_conditional.R")

weekly_path <- Sys.getenv(
  "MCOACH_WEEKLY_PATH",
  unset = ""
)
outcome_path <- Sys.getenv(
  "MCOACH_OUTCOME_PATH",
  unset = ""
)

output_dir <- "GEE_Code/HM/Y0_Baseline/Data_Results_Cond_New_Surrogate/Data_Like_Calibration"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spltime <- 1
times <- c(1, 2, 4)
t_eval <- 4
candidate_rhos <- c(0, 0.1, 0.3, 0.5)

clamp_prob <- function(x, eps = 1e-4) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

odds <- function(p) {
  p <- clamp_prob(p)
  p / (1 - p)
}

odds_ratio <- function(p_plus, p_minus) {
  odds(p_plus) / odds(p_minus)
}

safe_coef <- function(fit, term, default = 0) {
  coefs <- stats::coef(fit)
  if (!term %in% names(coefs) || !is.finite(coefs[[term]])) {
    return(default)
  }
  unname(coefs[[term]])
}

closest_candidate <- function(x, candidates = candidate_rhos) {
  if (!is.finite(x)) return(NA_real_)
  candidates[which.min(abs(candidates - x))]
}

predict_response <- function(fit, A1, Y0_ref) {
  clamp_prob(stats::predict(
    fit,
    newdata = data.frame(A1 = A1, Y0 = Y0_ref),
    type = "response"
  ))
}

predict_positive <- function(fit, A1, A2 = 0, R = NULL, Y0_ref, time_ref) {
  newdata <- data.frame(A1 = A1, A2_observed = A2, Y0 = Y0_ref, time = time_ref)
  if (!is.null(R)) newdata$R <- R
  clamp_prob(stats::predict(fit, newdata = newdata, type = "response"))
}

predict_count_mean <- function(fit, A1, A2 = 0, Y0_ref, time_ref) {
  as.numeric(stats::predict(
    fit,
    newdata = data.frame(A1 = A1, A2_observed = A2, Y0 = Y0_ref, time = time_ref),
    type = "response"
  ))
}

stop_if_missing <- function(path, label) {
  if (!nzchar(path) || !file.exists(path)) {
    stop(label, " file not found: ", path,
         "\nSet the corresponding environment variable if the path changed.")
  }
}

stop_if_missing(weekly_path, "Weekly")
stop_if_missing(outcome_path, "Outcome")

data_weekly <- read.csv(weekly_path, header = TRUE)
data_outcome <- read.csv(outcome_path, header = TRUE)

data_outcome_clean <- data_outcome %>%
  mutate(
    A1 = FIRSTRANDOM,
    A2 = T2,
    Y0 = BLALCCONSUMPTION,
    Y4 = FL4MALCCONSUMPTION,
    R = RESPONDENT,
    id = PID_SC
  ) %>%
  select(id, A1, A2, R, Y0, Y4) %>%
  filter(Y0 < 600, Y4 < 600)

data_weekly_1 <- data_weekly %>%
  select(
    PID_SC,
    AUD1_WK1, AUD2_WK1,
    AUD1_WK2, AUD2_WK2,
    AUD1_WK3, AUD2_WK3,
    AUD1_WK4, AUD2_WK4
  ) %>%
  mutate(
    id = PID_SC,
    wk1_drinks = if_else(AUD1_WK1 == 0, 0, AUD1_WK1 * AUD2_WK1),
    wk2_drinks = if_else(AUD1_WK2 == 0, 0, AUD1_WK2 * AUD2_WK2),
    wk3_drinks = if_else(AUD1_WK3 == 0, 0, AUD1_WK3 * AUD2_WK3),
    wk4_drinks = if_else(AUD1_WK4 == 0, 0, AUD1_WK4 * AUD2_WK4)
  ) %>%
  filter(!if_all(starts_with("wk"), is.na)) %>%
  mutate(
    Y1_weekly_mean = rowMeans(across(wk1_drinks:wk4_drinks), na.rm = TRUE),
    Y1 = round(Y1_weekly_mean * 30 / 7, 0)
  ) %>%
  select(id, Y1)

data_weekly_2 <- data_weekly %>%
  select(
    PID_SC,
    AUD1_WK5, AUD2_WK5,
    AUD1_WK6, AUD2_WK6,
    AUD1_WK7, AUD2_WK7,
    AUD1_WK8, AUD2_WK8
  ) %>%
  mutate(
    id = PID_SC,
    wk5_drinks = if_else(AUD1_WK5 == 0, 0, AUD1_WK5 * AUD2_WK5),
    wk6_drinks = if_else(AUD1_WK6 == 0, 0, AUD1_WK6 * AUD2_WK6),
    wk7_drinks = if_else(AUD1_WK7 == 0, 0, AUD1_WK7 * AUD2_WK7),
    wk8_drinks = if_else(AUD1_WK8 == 0, 0, AUD1_WK8 * AUD2_WK8)
  ) %>%
  filter(!if_all(starts_with("wk"), is.na)) %>%
  mutate(
    Y2_weekly_mean = rowMeans(across(wk5_drinks:wk8_drinks), na.rm = TRUE),
    Y2 = round(Y2_weekly_mean * 30 / 7, 0)
  ) %>%
  select(id, Y2)

data_wide <- merge(data_weekly_1, data_weekly_2, by = "id") %>%
  merge(data_outcome_clean, by = "id") %>%
  select(id, A1, A2, R, Y0, Y1, Y2, Y4)

if (!all(stats::na.omit(data_wide$A1) %in% c(-1, 1))) {
  stop("A1 must be coded as -1/1.")
}
if (!all(stats::na.omit(data_wide$A2) %in% c(-1, 1))) {
  stop("A2 must be coded as -1/1.")
}
if (!all(stats::na.omit(data_wide$R) %in% c(0, 1))) {
  stop("R must be coded as 0/1.")
}
if (anyNA(data_wide[, c("Y0", "Y1", "Y2", "Y4")])) {
  stop("Outcome data contain missing Y0/Y1/Y2/Y4 values after merging.")
}

data_long_all <- data_wide %>%
  pivot_longer(
    cols = c(Y0, Y1, Y2, Y4),
    names_to = "time",
    names_prefix = "Y",
    values_to = "Y"
  ) %>%
  mutate(
    time = as.numeric(time),
    is_zero = Y == 0,
    is_positive = Y > 0
  ) %>%
  arrange(id, time)

data_long_conditional <- data_wide %>%
  pivot_longer(
    cols = c(Y1, Y2, Y4),
    names_to = "time",
    names_prefix = "Y",
    values_to = "Y"
  ) %>%
  mutate(
    time = as.numeric(time),
    branch = case_when(
      time <= spltime ~ "pre",
      R == 1 ~ "responder",
      R == 0 ~ "nonresponder",
      TRUE ~ NA_character_
    ),
    A2_observed = ifelse(time <= spltime, 0, A2),
    is_zero = Y == 0,
    is_positive = Y > 0
  ) %>%
  arrange(id, time)

y0_summary <- data_wide %>%
  summarise(
    n = n(),
    Y0_mean = mean(Y0),
    Y0_sd = sd(Y0),
    Y0_median = median(Y0),
    Y0_prop_zero = mean(Y0 == 0),
    Y0_positive_mean = mean(Y0[Y0 > 0]),
    Y0_positive_sd = sd(Y0[Y0 > 0])
  )

Y0_ref <- y0_summary$Y0_median[[1]]

zero_by_time <- data_long_all %>%
  group_by(time) %>%
  summarise(
    n = n(),
    prop_zero = mean(is_zero),
    prop_positive = mean(is_positive),
    positive_count_mean = mean(Y[Y > 0]),
    positive_count_sd = sd(Y[Y > 0]),
    .groups = "drop"
  )

zero_by_branch_time <- data_long_conditional %>%
  group_by(branch, R, time) %>%
  summarise(
    n = n(),
    prop_zero = mean(is_zero),
    prop_positive = mean(is_positive),
    positive_count_mean = mean(Y[Y > 0]),
    positive_count_sd = sd(Y[Y > 0]),
    .groups = "drop"
  )

positive_cell_sizes <- data_long_conditional %>%
  filter(time > spltime) %>%
  group_by(time, A1, R, A2 = A2_observed) %>%
  summarise(
    n = n(),
    n_positive = sum(is_positive),
    prop_positive = mean(is_positive),
    .groups = "drop"
  )

response_by_A1 <- data_wide %>%
  group_by(A1) %>%
  summarise(
    n = n(),
    response_probability = mean(R == 1),
    Y0_mean = mean(Y0),
    Y0_median = median(Y0),
    .groups = "drop"
  )

response_fit <- stats::glm(R ~ A1 + Y0, family = stats::binomial(), data = data_wide)
response_fit_interaction <- stats::glm(R ~ A1 * Y0, family = stats::binomial(), data = data_wide)

response_prob_A1_minus <- predict_response(response_fit, A1 = -1, Y0_ref = Y0_ref)
response_prob_A1_plus <- predict_response(response_fit, A1 = 1, Y0_ref = Y0_ref)
response_Y0_log_or <- safe_coef(response_fit, "Y0")
response_A1_Y0_interaction_log_or <- safe_coef(response_fit_interaction, "A1:Y0")

pre_data <- data_long_conditional %>% filter(time <= spltime)
post_data <- data_long_conditional %>% filter(time > spltime)
responder_data <- post_data %>% filter(R == 1)
nonresponder_data <- post_data %>% filter(R == 0)

positive_pre_fit <- stats::glm(
  is_positive ~ A1 + Y0 + time,
  family = stats::binomial(),
  data = pre_data
)
positive_responder_fit <- stats::glm(
  is_positive ~ A1 + A2_observed + Y0 + time,
  family = stats::binomial(),
  data = responder_data
)
positive_nonresponder_fit <- stats::glm(
  is_positive ~ A1 + A2_observed + Y0 + time,
  family = stats::binomial(),
  data = nonresponder_data
)
positive_post_fit <- stats::glm(
  is_positive ~ R + A1 + A2_observed + Y0 + time,
  family = stats::binomial(),
  data = post_data
)

positive_Y0_log_or <- mean(c(
  safe_coef(positive_pre_fit, "Y0"),
  safe_coef(positive_responder_fit, "Y0"),
  safe_coef(positive_nonresponder_fit, "Y0")
), na.rm = TRUE)
positive_time_log_or <- mean(c(
  safe_coef(positive_pre_fit, "time"),
  safe_coef(positive_responder_fit, "time"),
  safe_coef(positive_nonresponder_fit, "time")
), na.rm = TRUE)

pre_pos_minus <- predict_positive(
  positive_pre_fit, A1 = -1, Y0_ref = Y0_ref, time_ref = spltime
)
pre_pos_plus <- predict_positive(
  positive_pre_fit, A1 = 1, Y0_ref = Y0_ref, time_ref = spltime
)
responder_pos_A2_minus <- predict_positive(
  positive_responder_fit, A1 = 0, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
responder_pos_A2_plus <- predict_positive(
  positive_responder_fit, A1 = 0, A2 = 1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_pos_A2_minus <- predict_positive(
  positive_nonresponder_fit, A1 = 0, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_pos_A2_plus <- predict_positive(
  positive_nonresponder_fit, A1 = 0, A2 = 1, Y0_ref = Y0_ref, time_ref = t_eval
)
responder_pos_A1_minus <- predict_positive(
  positive_responder_fit, A1 = -1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
responder_pos_A1_plus <- predict_positive(
  positive_responder_fit, A1 = 1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_pos_A1_minus <- predict_positive(
  positive_nonresponder_fit, A1 = -1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_pos_A1_plus <- predict_positive(
  positive_nonresponder_fit, A1 = 1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
post_responder_reference_pos <- predict_positive(
  positive_post_fit, A1 = 0, A2 = -1, R = 1, Y0_ref = Y0_ref, time_ref = t_eval
)
post_nonresponder_reference_pos <- predict_positive(
  positive_post_fit, A1 = 0, A2 = -1, R = 0, Y0_ref = Y0_ref, time_ref = t_eval
)

positive_data <- data_long_conditional %>% filter(Y > 0)
count_pre_fit <- stats::glm(
  Y ~ A1 + Y0 + time,
  family = stats::poisson(),
  data = positive_data %>% filter(time <= spltime)
)
count_responder_fit <- stats::glm(
  Y ~ A1 + A2_observed + Y0 + time,
  family = stats::poisson(),
  data = positive_data %>% filter(time > spltime, R == 1)
)
count_nonresponder_fit <- stats::glm(
  Y ~ A1 + A2_observed + Y0 + time,
  family = stats::poisson(),
  data = positive_data %>% filter(time > spltime, R == 0)
)

pre_count_mean_A1_minus <- predict_count_mean(
  count_pre_fit, A1 = -1, Y0_ref = Y0_ref, time_ref = spltime
)
pre_count_mean_A1_plus <- predict_count_mean(
  count_pre_fit, A1 = 1, Y0_ref = Y0_ref, time_ref = spltime
)
responder_count_mean_A2_minus <- predict_count_mean(
  count_responder_fit, A1 = 0, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
responder_count_mean_A2_plus <- predict_count_mean(
  count_responder_fit, A1 = 0, A2 = 1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_count_mean_A2_minus <- predict_count_mean(
  count_nonresponder_fit, A1 = 0, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_count_mean_A2_plus <- predict_count_mean(
  count_nonresponder_fit, A1 = 0, A2 = 1, Y0_ref = Y0_ref, time_ref = t_eval
)
responder_count_mean_A1_minus <- predict_count_mean(
  count_responder_fit, A1 = -1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
responder_count_mean_A1_plus <- predict_count_mean(
  count_responder_fit, A1 = 1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_count_mean_A1_minus <- predict_count_mean(
  count_nonresponder_fit, A1 = -1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)
nonresponder_count_mean_A1_plus <- predict_count_mean(
  count_nonresponder_fit, A1 = 1, A2 = -1, Y0_ref = Y0_ref, time_ref = t_eval
)

count_Y0_log_ratio <- mean(c(
  safe_coef(count_pre_fit, "Y0"),
  safe_coef(count_responder_fit, "Y0"),
  safe_coef(count_nonresponder_fit, "Y0")
), na.rm = TRUE)
count_time_log_ratio <- mean(c(
  safe_coef(count_pre_fit, "time"),
  safe_coef(count_responder_fit, "time"),
  safe_coef(count_nonresponder_fit, "time")
), na.rm = TRUE)

data_like_design_params <- tibble(
  design_input = c(
    "response_prob_A1_minus",
    "response_prob_A1_plus",
    "response_Y0_log_or",
    "positive_Y0_log_or",
    "positive_time_log_or",
    "first_stage_positive_or_A1_plus_vs_minus",
    "responder_positive_or_A1_plus_vs_minus",
    "nonresponder_positive_or_A1_plus_vs_minus",
    "responder_positive_or_A2_plus_vs_minus",
    "nonresponder_positive_or_vs_responder",
    "nonresponder_positive_or_A2_plus_vs_minus",
    "count_Y0_log_ratio",
    "count_time_log_ratio",
    "pre_count_mean_A1_minus",
    "pre_count_mean_A1_plus",
    "responder_count_mean_ratio_A1_plus_vs_minus",
    "nonresponder_count_mean_ratio_A1_plus_vs_minus",
    "responder_count_mean_A2_minus",
    "responder_count_mean_A2_plus",
    "nonresponder_count_mean_A2_minus",
    "nonresponder_count_mean_A2_plus"
  ),
  value = c(
    response_prob_A1_minus,
    response_prob_A1_plus,
    response_Y0_log_or,
    positive_Y0_log_or,
    positive_time_log_or,
    odds_ratio(pre_pos_plus, pre_pos_minus),
    odds_ratio(responder_pos_A1_plus, responder_pos_A1_minus),
    odds_ratio(nonresponder_pos_A1_plus, nonresponder_pos_A1_minus),
    odds_ratio(responder_pos_A2_plus, responder_pos_A2_minus),
    odds_ratio(post_nonresponder_reference_pos, post_responder_reference_pos),
    odds_ratio(nonresponder_pos_A2_plus, nonresponder_pos_A2_minus),
    count_Y0_log_ratio,
    count_time_log_ratio,
    pre_count_mean_A1_minus,
    pre_count_mean_A1_plus,
    responder_count_mean_A1_plus / responder_count_mean_A1_minus,
    nonresponder_count_mean_A1_plus / nonresponder_count_mean_A1_minus,
    responder_count_mean_A2_minus,
    responder_count_mean_A2_plus,
    nonresponder_count_mean_A2_minus,
    nonresponder_count_mean_A2_plus
  )
)

data_like_design_list <- stats::setNames(
  as.list(data_like_design_params$value),
  data_like_design_params$design_input
)

design_params_check <- conditional_params_from_design(
  design_params = data_like_design_list,
  zero_probability = zero_by_time$prop_zero[zero_by_time$time == t_eval],
  Y0_ref = Y0_ref,
  time_ref = t_eval,
  spltime = spltime
)
design_check_table <- conditional_design_check_table(
  params = design_params_check,
  Y0_ref = Y0_ref,
  time_ref = t_eval,
  spltime = spltime
)

wide_outcomes <- data_wide %>%
  select(id, Y0, Y1, Y2, Y4)

lag_pairs <- wide_outcomes %>%
  pivot_longer(cols = c(Y0, Y1, Y2, Y4), names_to = "time", values_to = "Y") %>%
  mutate(time = as.numeric(sub("^Y", "", time))) %>%
  arrange(id, time) %>%
  group_by(id) %>%
  summarise(
    time_pairs = list(t(utils::combn(time, 2))),
    Y_pairs = list(t(utils::combn(Y, 2))),
    .groups = "drop"
  ) %>%
  tidyr::unnest(cols = c(time_pairs, Y_pairs)) %>%
  mutate(
    time_1 = time_pairs[, 1],
    time_2 = time_pairs[, 2],
    lag = abs(time_2 - time_1),
    Y_1 = Y_pairs[, 1],
    Y_2 = Y_pairs[, 2],
    zero_1 = as.integer(Y_1 == 0),
    zero_2 = as.integer(Y_2 == 0)
  ) %>%
  select(-time_pairs, -Y_pairs)

lag_correlations <- lag_pairs %>%
  group_by(lag) %>%
  summarise(
    n_pairs = n(),
    raw_count_cor = stats::cor(Y_1, Y_2, use = "complete.obs"),
    zero_indicator_cor = stats::cor(zero_1, zero_2, use = "complete.obs"),
    positive_count_cor = {
      idx <- Y_1 > 0 & Y_2 > 0
      if (sum(idx) > 2) stats::cor(Y_1[idx], Y_2[idx], use = "complete.obs") else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    raw_count_ar1_monthly_rho = ifelse(raw_count_cor > 0, raw_count_cor^(1 / lag), NA_real_),
    zero_indicator_ar1_monthly_rho = ifelse(zero_indicator_cor > 0, zero_indicator_cor^(1 / lag), NA_real_),
    positive_count_ar1_monthly_rho = ifelse(positive_count_cor > 0, positive_count_cor^(1 / lag), NA_real_)
  )

rho_targets <- tibble(
  quantity = c(
    "raw_count_exchangeable_rho",
    "raw_count_ar1_monthly_rho",
    "zero_indicator_exchangeable_rho",
    "zero_indicator_ar1_monthly_rho",
    "positive_count_exchangeable_rho",
    "positive_count_ar1_monthly_rho"
  ),
  empirical_value = c(
    mean(lag_correlations$raw_count_cor, na.rm = TRUE),
    mean(lag_correlations$raw_count_ar1_monthly_rho, na.rm = TRUE),
    mean(lag_correlations$zero_indicator_cor, na.rm = TRUE),
    mean(lag_correlations$zero_indicator_ar1_monthly_rho, na.rm = TRUE),
    mean(lag_correlations$positive_count_cor, na.rm = TRUE),
    mean(lag_correlations$positive_count_ar1_monthly_rho, na.rm = TRUE)
  )
) %>%
  mutate(
    closest_candidate_rho = vapply(empirical_value, closest_candidate, numeric(1)),
    note = case_when(
      grepl("exchangeable", quantity) ~ "Use as an exchangeable data-like rho candidate.",
      grepl("ar1", quantity) ~ "Use as an AR1 data-like monthly rho candidate.",
      TRUE ~ NA_character_
    )
  )

one_part_formula <- Y ~ Y0 + time + A1 + R + A2_observed
poisson_fit <- stats::glm(one_part_formula, family = stats::poisson(), data = data_long_conditional)
data_long_conditional$pred_zero_poisson <- exp(-stats::predict(poisson_fit, type = "response"))

nb_available <- requireNamespace("MASS", quietly = TRUE)
if (nb_available) {
  nb_fit <- MASS::glm.nb(one_part_formula, data = data_long_conditional)
  nb_mu <- stats::predict(nb_fit, type = "response")
  nb_theta <- nb_fit$theta
  data_long_conditional$pred_zero_nb <- (nb_theta / (nb_theta + nb_mu))^nb_theta
} else {
  data_long_conditional$pred_zero_nb <- NA_real_
}

one_part_zero_diagnostics <- data_long_conditional %>%
  group_by(time) %>%
  summarise(
    observed_zero = mean(is_zero),
    predicted_zero_poisson = mean(pred_zero_poisson),
    predicted_zero_nb = mean(pred_zero_nb, na.rm = TRUE),
    excess_zero_ratio_poisson = observed_zero / predicted_zero_poisson,
    excess_zero_ratio_nb = observed_zero / predicted_zero_nb,
    n = n(),
    .groups = "drop"
  )

scenario_calibration_values <- tibble(
  quantity = c(
    "application_n",
    "Y0_ref_median",
    "Y0_mean",
    "Y0_sd",
    "Y0_prop_zero",
    "data_like_zero_percent_month4",
    "data_like_zero_percent_post_months_2_4_mean",
    "data_like_positive_count_mean_month4",
    "response_prob_A1_minus_at_Y0_ref",
    "response_prob_A1_plus_at_Y0_ref",
    "response_Y0_log_or",
    "response_A1_Y0_interaction_log_or",
    "raw_count_exchangeable_rho_closest",
    "raw_count_ar1_rho_closest",
    "zero_indicator_exchangeable_rho_closest",
    "zero_indicator_ar1_rho_closest",
    "min_postsplit_positive_cell_n",
    "median_postsplit_positive_cell_n",
    "month4_excess_zero_ratio_poisson",
    "month4_excess_zero_ratio_nb"
  ),
  value = c(
    nrow(data_wide),
    Y0_ref,
    y0_summary$Y0_mean,
    y0_summary$Y0_sd,
    y0_summary$Y0_prop_zero,
    100 * zero_by_time$prop_zero[zero_by_time$time == t_eval],
    100 * mean(zero_by_time$prop_zero[zero_by_time$time %in% c(2, 4)]),
    zero_by_time$positive_count_mean[zero_by_time$time == t_eval],
    response_prob_A1_minus,
    response_prob_A1_plus,
    response_Y0_log_or,
    response_A1_Y0_interaction_log_or,
    rho_targets$closest_candidate_rho[rho_targets$quantity == "raw_count_exchangeable_rho"],
    rho_targets$closest_candidate_rho[rho_targets$quantity == "raw_count_ar1_monthly_rho"],
    rho_targets$closest_candidate_rho[rho_targets$quantity == "zero_indicator_exchangeable_rho"],
    rho_targets$closest_candidate_rho[rho_targets$quantity == "zero_indicator_ar1_monthly_rho"],
    min(positive_cell_sizes$n_positive),
    stats::median(positive_cell_sizes$n_positive),
    one_part_zero_diagnostics$excess_zero_ratio_poisson[one_part_zero_diagnostics$time == t_eval],
    one_part_zero_diagnostics$excess_zero_ratio_nb[one_part_zero_diagnostics$time == t_eval]
  )
)

readr::write_csv(y0_summary, file.path(output_dir, "y0_summary.csv"))
readr::write_csv(zero_by_time, file.path(output_dir, "zero_by_time.csv"))
readr::write_csv(zero_by_branch_time, file.path(output_dir, "zero_by_branch_time.csv"))
readr::write_csv(positive_cell_sizes, file.path(output_dir, "positive_cell_sizes.csv"))
readr::write_csv(response_by_A1, file.path(output_dir, "response_by_A1.csv"))
readr::write_csv(lag_correlations, file.path(output_dir, "lag_correlations.csv"))
readr::write_csv(rho_targets, file.path(output_dir, "rho_recommendations.csv"))
readr::write_csv(one_part_zero_diagnostics, file.path(output_dir, "one_part_zero_diagnostics.csv"))
readr::write_csv(data_like_design_params, file.path(output_dir, "data_like_design_params.csv"))
readr::write_csv(design_check_table, file.path(output_dir, "data_like_design_check.csv"))
readr::write_csv(scenario_calibration_values, file.path(output_dir, "scenario_calibration_values.csv"))

saveRDS(
  list(
    data_like_design_params = data_like_design_list,
    scenario_calibration_values = scenario_calibration_values,
    zero_by_time = zero_by_time,
    rho_targets = rho_targets,
    one_part_zero_diagnostics = one_part_zero_diagnostics
  ),
  file.path(output_dir, "data_like_calibration.rds")
)

cat("Data-like calibration complete.\n")
cat("Output directory:", output_dir, "\n")
print(scenario_calibration_values, n = Inf)
