#######################################
# Conditional-First Hurdle SMART Utils #
#######################################
#
# Standalone helpers for the exact conditional-GLM simulation described in
# conditional_hurdle_glm_dtr.pdf. This file intentionally does not source or
# depend on the marginal-first projection utilities.

expit_conditional <- function(x) {
  stats::plogis(x)
}

ztp_mean_conditional <- function(lambda) {
  lambda <- pmax(as.numeric(lambda), .Machine$double.eps)
  lambda / pmax(1 - exp(-lambda), .Machine$double.eps)
}

ztp_lambda_from_mean_conditional <- function(mu, tol = 1e-10, maxit = 50) {
  mu <- pmax(as.numeric(mu), 1 + sqrt(.Machine$double.eps))
  lambda <- ifelse(mu < 1.05, pmax(2 * (mu - 1), 1e-8), pmax(mu - 0.5, 1e-8))
  for (iter in seq_len(maxit)) {
    exp_neg <- exp(-lambda)
    denom <- pmax(1 - exp_neg, .Machine$double.eps)
    f <- lambda / denom - mu
    fp <- (1 - exp_neg - lambda * exp_neg) / (denom^2)
    lambda_new <- pmax(lambda - f / pmax(fp, .Machine$double.eps), 1e-10)
    if (max(abs(lambda_new - lambda), na.rm = TRUE) < tol) {
      lambda <- lambda_new
      break
    }
    lambda <- lambda_new
  }
  lambda
}

ztp_var_conditional <- function(lambda) {
  lambda <- pmax(as.numeric(lambda), .Machine$double.eps)
  denom <- pmax(1 - exp(-lambda), .Machine$double.eps)
  second <- (lambda + lambda^2) / denom
  second - ztp_mean_conditional(lambda)^2
}

rtpois_conditional <- function(n, lambda) {
  lambda <- rep(lambda, length.out = n)
  p0 <- exp(-lambda)
  target <- stats::runif(n) * (1 - p0) + p0
  y <- stats::qpois(target, lambda = lambda)
  as.integer(pmax(y, 1L))
}

gen_Y0_conditional <- function(n) {
  y0_positive <- stats::rbinom(n, 1, expit_conditional(-0.3))
  y0_count <- rtpois_conditional(n, lambda = exp(1.7) + 2)
  ifelse(y0_positive == 1, y0_count, 0)
}

make_correlation_matrix_conditional <- function(times, corstr = c("independence", "exchangeable", "ar1"),
                                                rho = 0) {
  corstr <- match.arg(corstr)
  n_time <- length(times)
  if (corstr == "independence" || rho == 0) {
    return(diag(1, n_time))
  }
  if (corstr == "exchangeable") {
    return(diag(rep(1 - rho, n_time)) + matrix(rho, n_time, n_time))
  }
  lag_mat <- abs(outer(times, times, "-"))
  rho^lag_mat
}

correlated_uniforms_conditional <- function(times, corstr = c("independence", "exchangeable", "ar1"),
                                            rho = 0) {
  corstr <- match.arg(corstr)
  if (corstr == "independence" || rho == 0 || length(times) == 1) {
    return(stats::runif(length(times)))
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for correlated conditional simulation.")
  }
  cor_mat <- make_correlation_matrix_conditional(times = times, corstr = corstr, rho = rho)
  z <- as.numeric(mvtnorm::rmvnorm(1, mean = rep(0, length(times)), sigma = cor_mat, method = "svd"))
  stats::pnorm(z)
}

response_design_conditional <- function(A1, Y0, baseline_covariates = NULL) {
  design <- data.frame(
    `(Intercept)` = 1,
    A1 = A1,
    Y0 = Y0,
    check.names = FALSE
  )
  if (!is.null(baseline_covariates)) {
    baseline_covariates <- as.data.frame(baseline_covariates)
    if (nrow(baseline_covariates) != length(Y0)) {
      stop("baseline_covariates must have one row per subject.")
    }
    design <- cbind(design, baseline_covariates)
  }
  design
}

pre_design_conditional <- function(A1, Y0, time) {
  data.frame(
    `(Intercept)` = 1,
    Y0 = Y0,
    time = time,
    A1_time = A1 * time,
    check.names = FALSE
  )
}

unified_design_conditional <- function(A1, A2, R, Y0, time, spltime) {
  post <- as.numeric(time > spltime)
  R_ind <- as.numeric(R == 1)
  NR_ind <- as.numeric(R == 0)
  data.frame(
    `(Intercept)` = 1,
    Y0 = Y0,
    time = time,
    A1_time = A1 * time,
    R_time = post * R_ind * time,
    R_A1_time = post * R_ind * A1 * time,
    R_A2_time = post * R_ind * A2 * time,
    R_A1_A2_time = post * R_ind * A1 * A2 * time,
    NR_time = post * NR_ind * time,
    NR_A1_time = post * NR_ind * A1 * time,
    NR_A2_time = post * NR_ind * A2 * time,
    NR_A1_A2_time = post * NR_ind * A1 * A2 * time,
    check.names = FALSE
  )
}

add_unified_terms_conditional <- function(data, spltime) {
  data$post <- as.numeric(data$time > spltime)
  data$R_ind <- as.numeric(data$R == 1)
  data$NR_ind <- as.numeric(data$R == 0)
  data$A1_time <- data$A1 * data$time
  data$R_time <- data$post * data$R_ind * data$time
  data$R_A1_time <- data$post * data$R_ind * data$A1 * data$time
  data$R_A2_time <- data$post * data$R_ind * data$A2 * data$time
  data$R_A1_A2_time <- data$post * data$R_ind * data$A1 * data$A2 * data$time
  data$NR_time <- data$post * data$NR_ind * data$time
  data$NR_A1_time <- data$post * data$NR_ind * data$A1 * data$time
  data$NR_A2_time <- data$post * data$NR_ind * data$A2 * data$time
  data$NR_A1_A2_time <- data$post * data$NR_ind * data$A1 * data$A2 * data$time
  data
}

linear_predictor_conditional <- function(design, coefs) {
  coef_names <- names(coefs)
  coefs <- as.numeric(coefs)
  if (!is.null(coef_names)) {
    names(coefs) <- coef_names
    missing <- setdiff(coef_names, colnames(design))
    if (length(missing) > 0) {
      stop("Coefficient names not found in design: ", paste(missing, collapse = ", "))
    }
    design <- design[, coef_names, drop = FALSE]
  } else if (length(coefs) != ncol(design)) {
    stop("Unnamed coefficient vector length does not match design width.")
  }
  as.numeric(as.matrix(design) %*% coefs)
}

prob_from_odds_multiplier <- function(prob, odds_multiplier) {
  odds <- prob / (1 - prob)
  odds_new <- odds * odds_multiplier
  odds_new / (1 + odds_new)
}

make_two_level_logit_coefs <- function(prob_minus, prob_plus,
                                       varying_term_name,
                                       term_multiplier,
                                       y0_coef, time_coef,
                                       Y0_ref, time_ref) {
  eta_minus <- stats::qlogis(prob_minus)
  eta_plus <- stats::qlogis(prob_plus)
  varying_coef <- (eta_plus - eta_minus) / (2 * term_multiplier)
  intercept <- (eta_minus + eta_plus) / 2 -
    y0_coef * Y0_ref -
    time_coef * time_ref
  out <- c(`(Intercept)` = unname(intercept), Y0 = y0_coef, time = time_coef)
  out[varying_term_name] <- unname(varying_coef)
  out
}

make_two_level_log_lambda_coefs <- function(mean_minus, mean_plus,
                                            varying_term_name,
                                            term_multiplier,
                                            y0_coef, time_coef,
                                            Y0_ref, time_ref) {
  eta_minus <- log(ztp_lambda_from_mean_conditional(mean_minus))
  eta_plus <- log(ztp_lambda_from_mean_conditional(mean_plus))
  varying_coef <- (eta_plus - eta_minus) / (2 * term_multiplier)
  intercept <- (eta_minus + eta_plus) / 2 -
    y0_coef * Y0_ref -
    time_coef * time_ref
  out <- c(`(Intercept)` = unname(intercept), Y0 = y0_coef, time = time_coef)
  out[varying_term_name] <- unname(varying_coef)
  out
}

make_positive_or_logit_coef <- function(reference_positive_prob,
                                        positive_or_plus_vs_minus,
                                        term_multiplier) {
  positive_minus <- reference_positive_prob
  positive_plus <- prob_from_odds_multiplier(
    positive_minus,
    positive_or_plus_vs_minus
  )
  zero_minus <- 1 - positive_minus
  zero_plus <- 1 - positive_plus
  (stats::qlogis(zero_plus) - stats::qlogis(zero_minus)) / (2 * term_multiplier)
}

make_ztp_log_lambda_ratio_coef <- function(reference_count_mean,
                                           count_mean_ratio_plus_vs_minus,
                                           term_multiplier) {
  mean_minus <- reference_count_mean
  mean_plus <- reference_count_mean * count_mean_ratio_plus_vs_minus
  eta_minus <- log(ztp_lambda_from_mean_conditional(mean_minus))
  eta_plus <- log(ztp_lambda_from_mean_conditional(mean_plus))
  (eta_plus - eta_minus) / (2 * term_multiplier)
}

conditional_design_defaults <- function() {
  list(
    response_prob_A1_minus = 0.40,
    response_prob_A1_plus = 0.60,
    response_Y0_log_or = 0.02,
    
    positive_Y0_log_or = 0.02,
    positive_time_log_or = 0.00,
    first_stage_positive_or_A1_plus_vs_minus = 0.90,
    responder_positive_or_A1_plus_vs_minus = 0.75,
    nonresponder_positive_or_A1_plus_vs_minus = 0.95,
    responder_positive_or_A2_plus_vs_minus = 0.80,
    nonresponder_positive_or_vs_responder = 1.25,
    nonresponder_positive_or_A2_plus_vs_minus = 0.70,
    
    count_Y0_log_ratio = 0.02,
    count_time_log_ratio = 0.00,
    pre_count_mean_A1_minus = 12,
    pre_count_mean_A1_plus = 10,
    responder_count_mean_ratio_A1_plus_vs_minus = 0.75,
    nonresponder_count_mean_ratio_A1_plus_vs_minus = 0.95,
    responder_count_mean_A2_minus = 12,
    responder_count_mean_A2_plus = 10,
    nonresponder_count_mean_A2_minus = 15,
    nonresponder_count_mean_A2_plus = 12
  )
}

conditional_params_from_design <- function(design_params, zero_probability,
                                           Y0_ref = 0, time_ref = 4,
                                           spltime = 1) {
  if (!is.finite(zero_probability) || zero_probability <= 0 || zero_probability >= 1) {
    stop("zero_probability must be strictly between 0 and 1.")
  }
  if (!is.finite(time_ref) || time_ref <= spltime) {
    stop("time_ref must be a post-split time greater than spltime.")
  }
  positive_reference_prob <- 1 - zero_probability
  
  response_eta_minus <- stats::qlogis(design_params$response_prob_A1_minus)
  response_eta_plus <- stats::qlogis(design_params$response_prob_A1_plus)
  response_A1 <- (response_eta_plus - response_eta_minus) / 2
  response_intercept <- (response_eta_minus + response_eta_plus) / 2 -
    design_params$response_Y0_log_or * Y0_ref
  
  pre_positive_minus <- positive_reference_prob
  pre_positive_plus <- prob_from_odds_multiplier(
    pre_positive_minus,
    design_params$first_stage_positive_or_A1_plus_vs_minus
  )
  
  responder_positive_minus <- positive_reference_prob
  responder_positive_plus <- prob_from_odds_multiplier(
    responder_positive_minus,
    design_params$responder_positive_or_A2_plus_vs_minus
  )
  
  nonresponder_positive_minus <- prob_from_odds_multiplier(
    responder_positive_minus,
    design_params$nonresponder_positive_or_vs_responder
  )
  nonresponder_positive_plus <- prob_from_odds_multiplier(
    nonresponder_positive_minus,
    design_params$nonresponder_positive_or_A2_plus_vs_minus
  )
  
  pre_zero_minus <- 1 - pre_positive_minus
  pre_zero_plus <- 1 - pre_positive_plus
  responder_zero_minus <- 1 - responder_positive_minus
  responder_zero_plus <- 1 - responder_positive_plus
  nonresponder_zero_minus <- 1 - nonresponder_positive_minus
  nonresponder_zero_plus <- 1 - nonresponder_positive_plus
  
  beta_pre <- make_two_level_logit_coefs(
    prob_minus = pre_zero_minus,
    prob_plus = pre_zero_plus,
    varying_term_name = "A1_time",
    term_multiplier = spltime,
    y0_coef = -design_params$positive_Y0_log_or,
    time_coef = -design_params$positive_time_log_or,
    Y0_ref = Y0_ref,
    time_ref = spltime
  )
  
  beta_R_base <- make_two_level_logit_coefs(
    prob_minus = responder_zero_minus,
    prob_plus = responder_zero_plus,
    varying_term_name = "R_A2_time",
    term_multiplier = time_ref,
    y0_coef = -design_params$positive_Y0_log_or,
    time_coef = -design_params$positive_time_log_or,
    Y0_ref = Y0_ref,
    time_ref = time_ref
  )
  
  beta_NR_base <- make_two_level_logit_coefs(
    prob_minus = nonresponder_zero_minus,
    prob_plus = nonresponder_zero_plus,
    varying_term_name = "NR_A2_time",
    term_multiplier = time_ref,
    y0_coef = -design_params$positive_Y0_log_or,
    time_coef = -design_params$positive_time_log_or,
    Y0_ref = Y0_ref,
    time_ref = time_ref
  )

  beta_R_A1_total <- make_positive_or_logit_coef(
    reference_positive_prob = responder_positive_minus,
    positive_or_plus_vs_minus = design_params$responder_positive_or_A1_plus_vs_minus,
    term_multiplier = time_ref
  )
  beta_NR_A1_total <- make_positive_or_logit_coef(
    reference_positive_prob = nonresponder_positive_minus,
    positive_or_plus_vs_minus = design_params$nonresponder_positive_or_A1_plus_vs_minus,
    term_multiplier = time_ref
  )
  
  gamma_pre <- make_two_level_log_lambda_coefs(
    mean_minus = design_params$pre_count_mean_A1_minus,
    mean_plus = design_params$pre_count_mean_A1_plus,
    varying_term_name = "A1_time",
    term_multiplier = spltime,
    y0_coef = design_params$count_Y0_log_ratio,
    time_coef = design_params$count_time_log_ratio,
    Y0_ref = Y0_ref,
    time_ref = spltime
  )
  
  gamma_R_base <- make_two_level_log_lambda_coefs(
    mean_minus = design_params$responder_count_mean_A2_minus,
    mean_plus = design_params$responder_count_mean_A2_plus,
    varying_term_name = "R_A2_time",
    term_multiplier = time_ref,
    y0_coef = design_params$count_Y0_log_ratio,
    time_coef = design_params$count_time_log_ratio,
    Y0_ref = Y0_ref,
    time_ref = time_ref
  )
  
  gamma_NR_base <- make_two_level_log_lambda_coefs(
    mean_minus = design_params$nonresponder_count_mean_A2_minus,
    mean_plus = design_params$nonresponder_count_mean_A2_plus,
    varying_term_name = "NR_A2_time",
    term_multiplier = time_ref,
    y0_coef = design_params$count_Y0_log_ratio,
    time_coef = design_params$count_time_log_ratio,
    Y0_ref = Y0_ref,
    time_ref = time_ref
  )

  gamma_R_A1_total <- make_ztp_log_lambda_ratio_coef(
    reference_count_mean = design_params$responder_count_mean_A2_minus,
    count_mean_ratio_plus_vs_minus = design_params$responder_count_mean_ratio_A1_plus_vs_minus,
    term_multiplier = time_ref
  )
  gamma_NR_A1_total <- make_ztp_log_lambda_ratio_coef(
    reference_count_mean = design_params$nonresponder_count_mean_A2_minus,
    count_mean_ratio_plus_vs_minus = design_params$nonresponder_count_mean_ratio_A1_plus_vs_minus,
    term_multiplier = time_ref
  )
  
  shared_beta_names <- c("(Intercept)", "Y0", "time", "A1_time")
  beta <- c(
    beta_pre[shared_beta_names],
    R_time = unname((beta_R_base["(Intercept)"] - beta_pre["(Intercept)"]) / time_ref),
    R_A1_time = unname(beta_R_A1_total - beta_pre["A1_time"]),
    R_A2_time = unname(beta_R_base["R_A2_time"]),
    R_A1_A2_time = 0,
    NR_time = unname((beta_NR_base["(Intercept)"] - beta_pre["(Intercept)"]) / time_ref),
    NR_A1_time = unname(beta_NR_A1_total - beta_pre["A1_time"]),
    NR_A2_time = unname(beta_NR_base["NR_A2_time"]),
    NR_A1_A2_time = 0
  )
  
  shared_gamma_names <- c("(Intercept)", "Y0", "time", "A1_time")
  gamma <- c(
    gamma_pre[shared_gamma_names],
    R_time = unname((gamma_R_base["(Intercept)"] - gamma_pre["(Intercept)"]) / time_ref),
    R_A1_time = unname(gamma_R_A1_total - gamma_pre["A1_time"]),
    R_A2_time = unname(gamma_R_base["R_A2_time"]),
    R_A1_A2_time = 0,
    NR_time = unname((gamma_NR_base["(Intercept)"] - gamma_pre["(Intercept)"]) / time_ref),
    NR_A1_time = unname(gamma_NR_A1_total - gamma_pre["A1_time"]),
    NR_A2_time = unname(gamma_NR_base["NR_A2_time"]),
    NR_A1_A2_time = 0
  )
  
  list(
    alpha_response = c(
      `(Intercept)` = unname(response_intercept),
      A1 = unname(response_A1),
      Y0 = design_params$response_Y0_log_or
    ),
    beta = beta,
    gamma = gamma
  )
}

conditional_default_params <- function() {
  conditional_params_from_design(
    design_params = conditional_design_defaults(),
    zero_probability = 0.30,
    Y0_ref = 0,
    time_ref = 4,
    spltime = 1
  )
}

conditional_parameter_table <- function(params) {
  pieces <- names(params)
  out <- do.call(rbind, lapply(pieces, function(piece) {
    data.frame(
      model_piece = piece,
      parameter = names(params[[piece]]),
      value = as.numeric(params[[piece]]),
      row.names = NULL
    )
  }))
  rownames(out) <- NULL
  out
}

conditional_design_table <- function(design_params) {
  data.frame(
    design_input = names(design_params),
    value = as.numeric(unlist(design_params, use.names = FALSE)),
    row.names = NULL
  )
}

conditional_design_check_table <- function(params, Y0_ref = 0, time_ref = 4,
                                           spltime = 1,
                                           baseline_response_covariates_ref = NULL) {
  response_minus <- response_probability_conditional(
    -1,
    Y0_ref,
    params$alpha_response,
    baseline_covariates = baseline_response_covariates_ref
  )
  response_plus <- response_probability_conditional(
    1,
    Y0_ref,
    params$alpha_response,
    baseline_covariates = baseline_response_covariates_ref
  )
  
  pre_zero_minus <- zero_probability_unified_conditional(-1, 0, 1, Y0_ref, spltime, spltime, params$beta)
  pre_zero_plus <- zero_probability_unified_conditional(1, 0, 1, Y0_ref, spltime, spltime, params$beta)
  
  responder_zero_minus <- zero_probability_unified_conditional(0, -1, 1, Y0_ref, time_ref, spltime, params$beta)
  responder_zero_plus <- zero_probability_unified_conditional(0, 1, 1, Y0_ref, time_ref, spltime, params$beta)
  nonresponder_zero_minus <- zero_probability_unified_conditional(0, -1, 0, Y0_ref, time_ref, spltime, params$beta)
  nonresponder_zero_plus <- zero_probability_unified_conditional(0, 1, 0, Y0_ref, time_ref, spltime, params$beta)
  responder_zero_A1_minus <- zero_probability_unified_conditional(-1, -1, 1, Y0_ref, time_ref, spltime, params$beta)
  responder_zero_A1_plus <- zero_probability_unified_conditional(1, -1, 1, Y0_ref, time_ref, spltime, params$beta)
  nonresponder_zero_A1_minus <- zero_probability_unified_conditional(-1, -1, 0, Y0_ref, time_ref, spltime, params$beta)
  nonresponder_zero_A1_plus <- zero_probability_unified_conditional(1, -1, 0, Y0_ref, time_ref, spltime, params$beta)
  
  pre_mean_minus <- ztp_mean_conditional(lambda_unified_conditional(-1, 0, 1, Y0_ref, spltime, spltime, params$gamma))
  pre_mean_plus <- ztp_mean_conditional(lambda_unified_conditional(1, 0, 1, Y0_ref, spltime, spltime, params$gamma))
  responder_mean_minus <- ztp_mean_conditional(lambda_unified_conditional(0, -1, 1, Y0_ref, time_ref, spltime, params$gamma))
  responder_mean_plus <- ztp_mean_conditional(lambda_unified_conditional(0, 1, 1, Y0_ref, time_ref, spltime, params$gamma))
  nonresponder_mean_minus <- ztp_mean_conditional(lambda_unified_conditional(0, -1, 0, Y0_ref, time_ref, spltime, params$gamma))
  nonresponder_mean_plus <- ztp_mean_conditional(lambda_unified_conditional(0, 1, 0, Y0_ref, time_ref, spltime, params$gamma))
  responder_mean_A1_minus <- ztp_mean_conditional(lambda_unified_conditional(-1, -1, 1, Y0_ref, time_ref, spltime, params$gamma))
  responder_mean_A1_plus <- ztp_mean_conditional(lambda_unified_conditional(1, -1, 1, Y0_ref, time_ref, spltime, params$gamma))
  nonresponder_mean_A1_minus <- ztp_mean_conditional(lambda_unified_conditional(-1, -1, 0, Y0_ref, time_ref, spltime, params$gamma))
  nonresponder_mean_A1_plus <- ztp_mean_conditional(lambda_unified_conditional(1, -1, 0, Y0_ref, time_ref, spltime, params$gamma))
  
  data.frame(
    quantity = c(
      "response_prob_A1_minus",
      "response_prob_A1_plus",
      "pre_positive_prob_A1_minus",
      "pre_positive_prob_A1_plus",
      "pre_zero_prob_A1_minus",
      "pre_zero_prob_A1_plus",
      "responder_positive_prob_A2_minus",
      "responder_positive_prob_A2_plus",
      "responder_zero_prob_A2_minus",
      "responder_zero_prob_A2_plus",
      "nonresponder_positive_prob_A2_minus",
      "nonresponder_positive_prob_A2_plus",
      "nonresponder_zero_prob_A2_minus",
      "nonresponder_zero_prob_A2_plus",
      "responder_positive_prob_A1_minus_at_A2_minus",
      "responder_positive_prob_A1_plus_at_A2_minus",
      "responder_zero_prob_A1_minus_at_A2_minus",
      "responder_zero_prob_A1_plus_at_A2_minus",
      "nonresponder_positive_prob_A1_minus_at_A2_minus",
      "nonresponder_positive_prob_A1_plus_at_A2_minus",
      "nonresponder_zero_prob_A1_minus_at_A2_minus",
      "nonresponder_zero_prob_A1_plus_at_A2_minus",
      "pre_positive_count_mean_A1_minus",
      "pre_positive_count_mean_A1_plus",
      "responder_positive_count_mean_A2_minus",
      "responder_positive_count_mean_A2_plus",
      "nonresponder_positive_count_mean_A2_minus",
      "nonresponder_positive_count_mean_A2_plus",
      "responder_positive_count_mean_A1_minus_at_A2_minus",
      "responder_positive_count_mean_A1_plus_at_A2_minus",
      "nonresponder_positive_count_mean_A1_minus_at_A2_minus",
      "nonresponder_positive_count_mean_A1_plus_at_A2_minus"
    ),
    value = c(
      response_minus,
      response_plus,
      1 - pre_zero_minus,
      1 - pre_zero_plus,
      pre_zero_minus,
      pre_zero_plus,
      1 - responder_zero_minus,
      1 - responder_zero_plus,
      responder_zero_minus,
      responder_zero_plus,
      1 - nonresponder_zero_minus,
      1 - nonresponder_zero_plus,
      nonresponder_zero_minus,
      nonresponder_zero_plus,
      1 - responder_zero_A1_minus,
      1 - responder_zero_A1_plus,
      responder_zero_A1_minus,
      responder_zero_A1_plus,
      1 - nonresponder_zero_A1_minus,
      1 - nonresponder_zero_A1_plus,
      nonresponder_zero_A1_minus,
      nonresponder_zero_A1_plus,
      pre_mean_minus,
      pre_mean_plus,
      responder_mean_minus,
      responder_mean_plus,
      nonresponder_mean_minus,
      nonresponder_mean_plus,
      responder_mean_A1_minus,
      responder_mean_A1_plus,
      nonresponder_mean_A1_minus,
      nonresponder_mean_A1_plus
    ),
    row.names = NULL
  )
}

response_probability_conditional <- function(A1, Y0, alpha_response,
                                             baseline_covariates = NULL) {
  eta <- linear_predictor_conditional(
    response_design_conditional(A1, Y0, baseline_covariates = baseline_covariates),
    alpha_response
  )
  expit_conditional(eta)
}

zero_probability_pre_conditional <- function(A1, Y0, time, beta_pre) {
  eta <- linear_predictor_conditional(pre_design_conditional(A1, Y0, time), beta_pre)
  expit_conditional(eta)
}

positive_probability_pre_conditional <- function(A1, Y0, time, beta_pre) {
  1 - zero_probability_pre_conditional(A1, Y0, time, beta_pre)
}

lambda_pre_conditional <- function(A1, Y0, time, gamma_pre) {
  eta <- linear_predictor_conditional(pre_design_conditional(A1, Y0, time), gamma_pre)
  exp(eta)
}

zero_probability_unified_conditional <- function(A1, A2, R, Y0, time, spltime, beta) {
  eta <- linear_predictor_conditional(
    unified_design_conditional(A1, A2, R, Y0, time, spltime),
    beta
  )
  expit_conditional(eta)
}

positive_probability_unified_conditional <- function(A1, A2, R, Y0, time, spltime, beta) {
  1 - zero_probability_unified_conditional(A1, A2, R, Y0, time, spltime, beta)
}

lambda_unified_conditional <- function(A1, A2, R, Y0, time, spltime, gamma) {
  eta <- linear_predictor_conditional(
    unified_design_conditional(A1, A2, R, Y0, time, spltime),
    gamma
  )
  exp(eta)
}

conditional_dtr_truth <- function(params, dtr_grid, Y0_vec, times, spltime,
                                  baseline_response_covariates = NULL) {
  if (!is.null(baseline_response_covariates)) {
    baseline_response_covariates <- as.data.frame(baseline_response_covariates)
    if (nrow(baseline_response_covariates) != length(Y0_vec)) {
      stop("baseline_response_covariates must have one row per Y0 value.")
    }
  }
  rows <- list()
  idx <- 1L
  for (tt in times) {
    for (ii in seq_len(nrow(dtr_grid))) {
      A1 <- dtr_grid$A1[ii]
      A2R <- dtr_grid$A2R[ii]
      A2NR <- dtr_grid$A2NR[ii]
      dtr_label <- paste(A1, A2R, A2NR, sep = ",")
      
      if (tt <= spltime) {
        p_zero <- zero_probability_unified_conditional(A1, 0, 1, Y0_vec, tt, spltime, params$beta)
        q <- 1 - p_zero
        lambda <- lambda_unified_conditional(A1, 0, 1, Y0_vec, tt, spltime, params$gamma)
        pi_plus <- mean(q, na.rm = TRUE)
        mu_total <- mean(q * ztp_mean_conditional(lambda), na.rm = TRUE)
        rows[[idx]] <- data.frame(
          time = tt, A1 = A1, A2R = A2R, A2NR = A2NR, DTR = dtr_label,
          response_probability = NA_real_,
          p_zero_true = mean(p_zero, na.rm = TRUE),
          m_plus_true = ifelse(pi_plus > 0, mu_total / pi_plus, NA_real_),
          q_R_mean = NA_real_,
          q_NR_mean = NA_real_,
          m_R_mean = NA_real_,
          m_NR_mean = NA_real_
        )
      } else {
        r <- response_probability_conditional(
          A1,
          Y0_vec,
          params$alpha_response,
          baseline_covariates = baseline_response_covariates
        )
        p_R <- zero_probability_unified_conditional(A1, A2R, 1, Y0_vec, tt, spltime, params$beta)
        p_NR <- zero_probability_unified_conditional(A1, A2NR, 0, Y0_vec, tt, spltime, params$beta)
        q_R <- 1 - p_R
        q_NR <- 1 - p_NR
        lambda_R <- lambda_unified_conditional(A1, A2R, 1, Y0_vec, tt, spltime, params$gamma)
        lambda_NR <- lambda_unified_conditional(A1, A2NR, 0, Y0_vec, tt, spltime, params$gamma)
        m_R <- ztp_mean_conditional(lambda_R)
        m_NR <- ztp_mean_conditional(lambda_NR)
        
        pi_plus_i <- r * q_R + (1 - r) * q_NR
        mu_i <- r * q_R * m_R + (1 - r) * q_NR * m_NR
        pi_plus <- mean(pi_plus_i, na.rm = TRUE)
        mu_total <- mean(mu_i, na.rm = TRUE)
        p_zero <- mean(r * p_R + (1 - r) * p_NR, na.rm = TRUE)
        
        rows[[idx]] <- data.frame(
          time = tt, A1 = A1, A2R = A2R, A2NR = A2NR, DTR = dtr_label,
          response_probability = mean(r, na.rm = TRUE),
          p_zero_true = p_zero,
          m_plus_true = ifelse(pi_plus > 0, mu_total / pi_plus, NA_real_),
          q_R_mean = mean(q_R, na.rm = TRUE),
          q_NR_mean = mean(q_NR, na.rm = TRUE),
          m_R_mean = mean(m_R, na.rm = TRUE),
          m_NR_mean = mean(m_NR, na.rm = TRUE)
        )
      }
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}
