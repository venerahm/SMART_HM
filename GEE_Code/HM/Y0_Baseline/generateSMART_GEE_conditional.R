#########################################
# Conditional-First SMART Data Generator #
#########################################
#
# Generates observed SMART data from exact conditional hurdle GLMs. This file
# intentionally does not source or depend on the marginal-first generator.

generateSMART_gee_conditional <- function(n, times, spltime, params,
                                          balanceRand = FALSE,
                                          corstr = c("independence", "exchangeable", "ar1"),
                                          rho = 0,
                                          baseline_response_covariates = NULL) {
  required_helpers <- c(
    "gen_Y0_conditional",
    "response_probability_conditional",
    "zero_probability_unified_conditional",
    "lambda_unified_conditional",
    "ztp_mean_conditional",
    "rtpois_conditional",
    "correlated_uniforms_conditional"
  )
  missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]
  if (length(missing_helpers) > 0) {
    stop("Source functions_HM_conditional.R before using generateSMART_gee_conditional(). Missing: ",
         paste(missing_helpers, collapse = ", "))
  }
  
  corstr <- match.arg(corstr)
  if (sum(times > spltime) == 0) {
    stop("At least one post-split time is required.")
  }
  
  subject_data <- data.frame(
    id = seq_len(n),
    Y0 = gen_Y0_conditional(n),
    A1 = 0,
    R = NA_integer_,
    A2R = 0,
    A2NR = 0
  )
  if (!is.null(baseline_response_covariates)) {
    baseline_response_covariates <- as.data.frame(baseline_response_covariates)
    if (nrow(baseline_response_covariates) != n) {
      stop("baseline_response_covariates must have n rows.")
    }
    subject_data <- cbind(subject_data, baseline_response_covariates)
  }
  
  if (balanceRand) {
    treated <- sample(seq_len(n), size = floor(n / 2), replace = FALSE)
    subject_data$A1[treated] <- 1
    subject_data$A1[-treated] <- -1
  } else {
    subject_data$A1 <- 2 * stats::rbinom(n, 1, 0.5) - 1
  }
  
  r_prob <- response_probability_conditional(
    A1 = subject_data$A1,
    Y0 = subject_data$Y0,
    alpha_response = params$alpha_response,
    baseline_covariates = baseline_response_covariates
  )
  subject_data$R <- stats::rbinom(n, 1, r_prob)
  
  if (balanceRand) {
    strata <- expand.grid(A1 = c(-1, 1), R = c(0, 1))
    for (ss in seq_len(nrow(strata))) {
      idx <- which(subject_data$A1 == strata$A1[ss] & subject_data$R == strata$R[ss])
      if (length(idx) == 0) next
      plus_idx <- sample(idx, size = floor(length(idx) / 2), replace = FALSE)
      minus_idx <- setdiff(idx, plus_idx)
      if (strata$R[ss] == 1) {
        subject_data$A2R[plus_idx] <- 1
        subject_data$A2R[minus_idx] <- -1
      } else {
        subject_data$A2NR[plus_idx] <- 1
        subject_data$A2NR[minus_idx] <- -1
      }
    }
  } else {
    resp_idx <- which(subject_data$R == 1)
    nonresp_idx <- which(subject_data$R == 0)
    subject_data$A2R[resp_idx] <- 2 * stats::rbinom(length(resp_idx), 1, 0.5) - 1
    subject_data$A2NR[nonresp_idx] <- 2 * stats::rbinom(length(nonresp_idx), 1, 0.5) - 1
  }
  
  long_rows <- vector("list", n)
  y_wide <- matrix(NA_integer_, nrow = n, ncol = length(times))
  colnames(y_wide) <- paste0("Y", times)
  
  for (ii in seq_len(n)) {
    p_zero <- numeric(length(times))
    lambda <- numeric(length(times))
    branch <- character(length(times))
    A2_observed <- numeric(length(times))
    
    for (jj in seq_along(times)) {
      tt <- times[jj]
      if (tt <= spltime) {
        branch[jj] <- "pre"
        A2_observed[jj] <- 0
        p_zero[jj] <- zero_probability_unified_conditional(
          A1 = subject_data$A1[ii],
          A2 = 0,
          R = 1,
          Y0 = subject_data$Y0[ii],
          time = tt,
          spltime = spltime,
          beta = params$beta
        )
        lambda[jj] <- lambda_unified_conditional(
          A1 = subject_data$A1[ii],
          A2 = 0,
          R = 1,
          Y0 = subject_data$Y0[ii],
          time = tt,
          spltime = spltime,
          gamma = params$gamma
        )
      } else if (subject_data$R[ii] == 1) {
        branch[jj] <- "responder"
        A2_observed[jj] <- subject_data$A2R[ii]
        p_zero[jj] <- zero_probability_unified_conditional(
          A1 = subject_data$A1[ii],
          A2 = subject_data$A2R[ii],
          R = 1,
          Y0 = subject_data$Y0[ii],
          time = tt,
          spltime = spltime,
          beta = params$beta
        )
        lambda[jj] <- lambda_unified_conditional(
          A1 = subject_data$A1[ii],
          A2 = subject_data$A2R[ii],
          R = 1,
          Y0 = subject_data$Y0[ii],
          time = tt,
          spltime = spltime,
          gamma = params$gamma
        )
      } else {
        branch[jj] <- "nonresponder"
        A2_observed[jj] <- subject_data$A2NR[ii]
        p_zero[jj] <- zero_probability_unified_conditional(
          A1 = subject_data$A1[ii],
          A2 = subject_data$A2NR[ii],
          R = 0,
          Y0 = subject_data$Y0[ii],
          time = tt,
          spltime = spltime,
          beta = params$beta
        )
        lambda[jj] <- lambda_unified_conditional(
          A1 = subject_data$A1[ii],
          A2 = subject_data$A2NR[ii],
          R = 0,
          Y0 = subject_data$Y0[ii],
          time = tt,
          spltime = spltime,
          gamma = params$gamma
        )
      }
    }
    
    u_binary <- correlated_uniforms_conditional(times = times, corstr = corstr, rho = rho)
    B <- as.integer(u_binary < p_zero)
    Y <- integer(length(times))
    pos_idx <- which(B == 0)
    if (length(pos_idx) > 0) {
      u_count <- correlated_uniforms_conditional(times = times, corstr = corstr, rho = rho)
      p0 <- exp(-lambda[pos_idx])
      target <- u_count[pos_idx] * (1 - p0) + p0
      Y[pos_idx] <- as.integer(pmax(stats::qpois(target, lambda = lambda[pos_idx]), 1L))
    }
    y_wide[ii, ] <- Y
    
    row_ii <- data.frame(
      id = subject_data$id[ii],
      Y0 = subject_data$Y0[ii],
      A1 = subject_data$A1[ii],
      R = subject_data$R[ii],
      A2R = subject_data$A2R[ii],
      A2NR = subject_data$A2NR[ii],
      time = times,
      branch = branch,
      A2 = A2_observed,
      p_zero_true_row = p_zero,
      lambda_true_row = lambda,
      B = B,
      Y = Y
    )
    if (!is.null(baseline_response_covariates)) {
      cov_row <- subject_data[
        rep(ii, length(times)),
        names(baseline_response_covariates),
        drop = FALSE
      ]
      rownames(cov_row) <- NULL
      row_ii <- cbind(row_ii, cov_row)
    }
    long_rows[[ii]] <- row_ii
  }
  
  wide_data <- cbind(subject_data, as.data.frame(y_wide))
  long_data <- do.call(rbind, long_rows)
  rownames(long_data) <- NULL
  
  output <- list(
    data = wide_data,
    long_data = long_data,
    params = list(
      times = times,
      spltime = spltime,
      conditional_params = params,
      corstr = corstr,
      rho = rho
    )
  )
  class(output) <- c("generateSMART_conditional", class(output))
  output
}
