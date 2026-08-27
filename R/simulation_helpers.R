########################################################
# Shared Utilities For Conditional Simulation Wrappers #
########################################################

immediate_cat <- function(...) {
  cat(..., "\n", file = stdout())
  try(flush.console(), silent = TRUE)
  invisible(NULL)
}

positive_cell_summary <- function(long_data, spltime) {
  post_data <- long_data[long_data$time > spltime, , drop = FALSE]
  if (nrow(post_data) == 0) {
    return(data.frame(
      min_postsplit_positive_cell_n = NA_integer_,
      median_postsplit_positive_cell_n = NA_real_,
      cells_below_10_positive = NA_integer_,
      cells_below_20_positive = NA_integer_,
      cells_below_30_positive = NA_integer_,
      cells_below_50_positive = NA_integer_
    ))
  }
  
  post_data$cell <- interaction(
    post_data$time,
    post_data$A1,
    post_data$R,
    post_data$A2,
    drop = TRUE
  )
  n_positive <- stats::aggregate(
    as.integer(post_data$Y > 0),
    by = list(cell = post_data$cell),
    FUN = sum
  )$x
  
  data.frame(
    min_postsplit_positive_cell_n = min(n_positive),
    median_postsplit_positive_cell_n = stats::median(n_positive),
    cells_below_10_positive = sum(n_positive < 10),
    cells_below_20_positive = sum(n_positive < 20),
    cells_below_30_positive = sum(n_positive < 30),
    cells_below_50_positive = sum(n_positive < 50)
  )
}

fit_conditional_iteration_for_wrapper <- function(generated, spltime, dtr_grid,
                                                  t_eval, use_delta_se,
                                                  scenario, iter) {
  out <- tryCatch({
    fits <- fit_conditional_hurdle_models(generated$long_data, spltime = spltime)
    est <- if (isTRUE(use_delta_se)) {
      predict_conditional_dtr_means_delta(
        fits = fits,
        dtr_grid = dtr_grid,
        Y0_vec = generated$data$Y0,
        t_eval = t_eval
      )
    } else {
      predict_conditional_dtr_means(
        fits = fits,
        dtr_grid = dtr_grid,
        Y0_vec = generated$data$Y0,
        t_eval = t_eval
      )
    }
    
    est$iteration <- iter
    est$scenario_id <- scenario$scenario_id
    est$scenario_label <- scenario$scenario_label
    est$scenario_family <- scenario$scenario_family
    est$ZI <- scenario$ZI
    est$corstr <- scenario$corstr
    est$rho <- scenario$rho
    if ("sample_size" %in% names(scenario)) {
      est$sample_size <- scenario$sample_size
    }
    if ("response_profile" %in% names(scenario)) {
      est$response_profile <- scenario$response_profile
    }
    if ("response_prob_A1_minus" %in% names(scenario)) {
      est$response_prob_A1_minus <- scenario$response_prob_A1_minus
    }
    if ("response_prob_A1_plus" %in% names(scenario)) {
      est$response_prob_A1_plus <- scenario$response_prob_A1_plus
    }
    
    list(estimate = est, error = NULL)
  }, error = function(e) {
    immediate_cat("    ERROR in scenario ", scenario$scenario_id,
                  ", iteration ", iter, ": ", conditionMessage(e))
    error <- data.frame(
      scenario_id = scenario$scenario_id,
      scenario_label = scenario$scenario_label,
      scenario_family = scenario$scenario_family,
      ZI = scenario$ZI,
      corstr = scenario$corstr,
      rho = scenario$rho,
      iteration = iter,
      message = conditionMessage(e)
    )
    if ("sample_size" %in% names(scenario)) {
      error$sample_size <- scenario$sample_size
    }
    if ("response_profile" %in% names(scenario)) {
      error$response_profile <- scenario$response_profile
    }
    list(
      estimate = NULL,
      error = error
    )
  })
  out
}

nb_zero_probability <- function(mu, theta) {
  exp(stats::dnbinom(0, size = theta, mu = mu, log = TRUE))
}

poisson_zero_probability <- function(mu) {
  exp(stats::dpois(0, lambda = mu, log = TRUE))
}

fit_nb_regression_mle <- function(formula, data) {
  model_frame <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- stats::model.response(model_frame)
  x <- stats::model.matrix(attr(model_frame, "terms"), data = model_frame)
  
  poisson_fit <- stats::glm.fit(
    x = x,
    y = y,
    family = stats::poisson()
  )
  beta_start <- stats::coef(poisson_fit)
  beta_start[!is.finite(beta_start)] <- 0
  
  y_mean <- mean(y)
  y_var <- stats::var(y)
  theta_start <- if (is.finite(y_var) && y_var > y_mean) {
    y_mean^2 / (y_var - y_mean)
  } else {
    10
  }
  theta_start <- min(max(theta_start, 0.01), 100)
  
  objective <- function(par) {
    beta <- par[-length(par)]
    theta <- exp(par[length(par)])
    eta <- as.vector(x %*% beta)
    mu <- pmax(exp(pmin(eta, 700)), .Machine$double.eps)
    -sum(stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
  }
  
  fit <- stats::optim(
    par = c(beta_start, log(theta_start)),
    fn = objective,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  if (!is.finite(fit$value) || fit$convergence != 0) {
    fit <- stats::optim(
      par = c(beta_start, log(theta_start)),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 2000)
    )
  }
  if (!is.finite(fit$value)) {
    stop("Negative binomial likelihood optimization failed.")
  }
  
  beta <- fit$par[-length(fit$par)]
  names(beta) <- colnames(x)
  list(
    coefficients = beta,
    theta = exp(fit$par[length(fit$par)]),
    terms = stats::terms(formula),
    xlevels = .getXlevels(stats::terms(formula), model_frame),
    contrasts = attr(x, "contrasts"),
    convergence = fit$convergence,
    value = fit$value
  )
}

predict_nb_regression_mle <- function(fit, newdata) {
  model_frame <- stats::model.frame(
    fit$terms,
    data = newdata,
    na.action = stats::na.pass,
    xlev = fit$xlevels
  )
  x <- stats::model.matrix(fit$terms, data = model_frame, contrasts.arg = fit$contrasts)
  eta <- as.vector(x %*% fit$coefficients)
  pmax(exp(pmin(eta, 700)), .Machine$double.eps)
}

fit_positive_count_zero_diagnostics <- function(long_data, spltime,
                                                scenario, iter) {
  count_formula <- Y ~ Y0 + time + A1_time +
    R_time + R_A1_time + R_A2_time + R_A1_A2_time +
    NR_time + NR_A1_time + NR_A2_time + NR_A1_A2_time
  
  out <- tryCatch({
    model_data <- add_unified_terms_conditional(long_data, spltime = spltime)
    
    fits_hm <- fit_conditional_hurdle_models(
      long_data = long_data,
      spltime = spltime,
      robust = FALSE
    )
    model_data$pred_zero_hm <- stats::predict(
      fits_hm$binary,
      newdata = model_data,
      type = "response"
    )

    poisson_all_fit <- stats::glm(
      count_formula,
      family = stats::poisson(),
      data = model_data
    )
    mu_poisson_all <- stats::predict(
      poisson_all_fit,
      newdata = model_data,
      type = "response"
    )
    model_data$pred_zero_all_poisson <- poisson_zero_probability(mu_poisson_all)
    
    model_data$pred_zero_all_nb <- NA_real_
    nb_all_fit <- tryCatch(
      fit_nb_regression_mle(count_formula, data = model_data),
      error = function(e) NULL
    )
    if (!is.null(nb_all_fit)) {
      mu_nb_all <- predict_nb_regression_mle(nb_all_fit, newdata = model_data)
      theta_nb_all <- nb_all_fit$theta
      model_data$pred_zero_all_nb <- nb_zero_probability(
        mu = mu_nb_all,
        theta = theta_nb_all
      )
    }
    
    overall_summary <- data.frame(
      time = "overall",
      observed_zero = mean(model_data$Y == 0),
      predicted_zero_hm = mean(model_data$pred_zero_hm),
      predicted_zero_all_poisson = mean(model_data$pred_zero_all_poisson),
      predicted_zero_all_nb = mean(model_data$pred_zero_all_nb, na.rm = TRUE),
      n = nrow(model_data),
      n_positive = sum(model_data$Y > 0)
    )
    
    time_summary <- stats::aggregate(
      cbind(
        observed_zero = as.numeric(model_data$Y == 0),
        predicted_zero_hm = model_data$pred_zero_hm,
        predicted_zero_all_poisson = model_data$pred_zero_all_poisson,
        predicted_zero_all_nb = model_data$pred_zero_all_nb,
        n_positive = as.numeric(model_data$Y > 0)
      ),
      by = list(time = model_data$time),
      FUN = mean,
      na.rm = TRUE
    )
    time_n <- stats::aggregate(
      model_data$Y,
      by = list(time = model_data$time),
      FUN = length
    )
    names(time_n)[2] <- "n"
    time_summary <- merge(time_summary, time_n, by = "time", all.x = TRUE)
    time_summary$n_positive <- time_summary$n_positive * time_summary$n
    time_summary$time <- as.character(time_summary$time)
    
    summary <- rbind(overall_summary, time_summary)
    summary$excess_zero_ratio_all_poisson <- summary$observed_zero /
      summary$predicted_zero_all_poisson
    summary$excess_zero_ratio_all_nb <- summary$observed_zero /
      summary$predicted_zero_all_nb
    
    summary$iteration <- iter
    summary$scenario_id <- scenario$scenario_id
    summary$scenario_label <- scenario$scenario_label
    summary$scenario_family <- scenario$scenario_family
    summary$ZI <- scenario$ZI
    summary$corstr <- scenario$corstr
    summary$rho <- scenario$rho
    if ("sample_size" %in% names(scenario)) {
      summary$sample_size <- scenario$sample_size
    }
    if ("response_profile" %in% names(scenario)) {
      summary$response_profile <- scenario$response_profile
    }
    
    list(summary = summary, error = NULL)
  }, error = function(e) {
    immediate_cat("    ERROR in stress diagnostic scenario ", scenario$scenario_id,
                  ", iteration ", iter, ": ", conditionMessage(e))
    error <- data.frame(
      scenario_id = scenario$scenario_id,
      scenario_label = scenario$scenario_label,
      scenario_family = scenario$scenario_family,
      ZI = scenario$ZI,
      corstr = scenario$corstr,
      rho = scenario$rho,
      iteration = iter,
      message = conditionMessage(e)
    )
    if ("sample_size" %in% names(scenario)) {
      error$sample_size <- scenario$sample_size
    }
    if ("response_profile" %in% names(scenario)) {
      error$response_profile <- scenario$response_profile
    }
    list(
      summary = NULL,
      error = error
    )
  })
  out
}
