###############################################
# Conditional Hurdle GLM G-Computation Fitting #
###############################################
#
# Fits the exact observed-data model used by generateSMART_gee_conditional()
# and estimates DTR means by standardized g-computation.

log1mexp_conditional <- function(x) {
  x <- pmax(as.numeric(x), .Machine$double.eps)
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

cluster_sandwich_conditional <- function(score_mat, bread_inv, cluster) {
  cluster <- as.factor(cluster)
  cluster_levels <- levels(cluster)
  p <- ncol(score_mat)
  meat <- matrix(0, p, p)
  
  for (cl in cluster_levels) {
    idx <- which(cluster == cl)
    score_sum <- colSums(score_mat[idx, , drop = FALSE])
    meat <- meat + tcrossprod(score_sum)
  }
  
  n <- nrow(score_mat)
  g <- length(cluster_levels)
  correction <- if (g > 1 && n > p) {
    (g / (g - 1)) * ((n - 1) / (n - p))
  } else {
    1
  }
  
  vcov <- correction * bread_inv %*% meat %*% bread_inv
  rownames(vcov) <- colnames(vcov) <- colnames(score_mat)
  vcov
}

robust_vcov_glm_conditional <- function(fit, data, cluster_id = "id") {
  if (!cluster_id %in% names(data)) {
    stop("cluster_id column not found in data: ", cluster_id)
  }
  X <- stats::model.matrix(fit)
  y <- stats::model.response(stats::model.frame(fit))
  mu <- stats::fitted(fit)
  score_mat <- X * as.numeric(y - mu)
  bread_inv <- stats::vcov(fit)
  cluster_sandwich_conditional(score_mat, bread_inv, data[[cluster_id]])
}

fit_ztp_glm_conditional <- function(formula, data, maxit = 500,
                                    cluster_id = "id",
                                    robust = TRUE) {
  mf <- stats::model.frame(formula, data = data)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  y <- as.numeric(stats::model.response(mf))
  if (any(y <= 0, na.rm = TRUE)) {
    stop("Zero-truncated Poisson GLM requires positive counts only.")
  }
  
  start_fit <- stats::glm(formula, family = stats::poisson(), data = data)
  start <- stats::coef(start_fit)
  start[!is.finite(start)] <- 0
  
  nll <- function(beta) {
    eta <- as.numeric(X %*% beta)
    lambda <- pmin(exp(eta), 1e8)
    -sum(-lambda + y * eta - lgamma(y + 1) - log1mexp_conditional(lambda))
  }
  
  opt <- stats::optim(start, nll, method = "BFGS", hessian = TRUE,
                      control = list(maxit = maxit))
  vcov_model <- tryCatch(solve(opt$hessian), error = function(e) {
    matrix(NA_real_, nrow = length(opt$par), ncol = length(opt$par))
  })
  rownames(vcov_model) <- colnames(vcov_model) <- names(opt$par)
  
  eta <- as.numeric(X %*% opt$par)
  lambda <- exp(eta)
  m <- ztp_mean_conditional(lambda)
  score_mat <- X * as.numeric(y - m)
  colnames(score_mat) <- names(opt$par)
  
  vcov <- if (isTRUE(robust) && cluster_id %in% names(data)) {
    cluster_sandwich_conditional(score_mat, vcov_model, data[[cluster_id]])
  } else {
    vcov_model
  }
  
  fit <- list(
    coefficients = opt$par,
    vcov = vcov,
    vcov_model = vcov_model,
    formula = formula,
    terms = stats::terms(formula, data = data),
    converged = opt$convergence == 0,
    nll = opt$value,
    data = data,
    cluster_id = cluster_id
  )
  class(fit) <- "ztp_glm_conditional"
  fit
}

predict_ztp_glm_conditional <- function(object, newdata,
                                        type = c("lambda", "m_plus", "link")) {
  type <- match.arg(type)
  X <- stats::model.matrix(delete.response(object$terms), newdata)
  eta <- as.numeric(X %*% object$coefficients)
  if (type == "link") return(eta)
  lambda <- exp(eta)
  if (type == "lambda") return(lambda)
  ztp_mean_conditional(lambda)
}

fit_conditional_hurdle_models <- function(long_data, spltime,
                                          response_formula = R ~ A1 + Y0,
                                          binary_formula = B ~ Y0 + time + A1_time +
                                            R_time + R_A1_time + R_A2_time + R_A1_A2_time +
                                            NR_time + NR_A1_time + NR_A2_time + NR_A1_A2_time,
                                          count_formula = Y ~ Y0 + time + A1_time +
                                            R_time + R_A1_time + R_A2_time + R_A1_A2_time +
                                            NR_time + NR_A1_time + NR_A2_time + NR_A1_A2_time,
                                          cluster_id = "id",
                                          robust = TRUE) {
  required_helpers <- c("ztp_mean_conditional", "add_unified_terms_conditional")
  missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]
  if (length(missing_helpers) > 0) {
    stop("Source functions_HM_conditional.R before fitting. Missing: ",
         paste(missing_helpers, collapse = ", "))
  }
  
  model_data <- add_unified_terms_conditional(long_data, spltime = spltime)
  subject_data <- model_data[!duplicated(model_data$id), c("id", "Y0", "A1", "R")]
  positive_data <- model_data[model_data$Y > 0, , drop = FALSE]
  
  if (nrow(positive_data) == 0) {
    stop("Positive counts are required to fit the ZTP model.")
  }
  
  response_fit <- stats::glm(response_formula, family = stats::binomial(), data = subject_data)
  binary_fit <- stats::glm(binary_formula, family = stats::binomial(), data = model_data)
  
  if (isTRUE(robust)) {
    response_fit$vcov_robust <- robust_vcov_glm_conditional(response_fit, subject_data, cluster_id)
    binary_fit$vcov_robust <- robust_vcov_glm_conditional(binary_fit, model_data, cluster_id)
  }
  
  list(
    response = response_fit,
    binary = binary_fit,
    count = fit_ztp_glm_conditional(count_formula, data = positive_data,
                                    cluster_id = cluster_id, robust = robust),
    formulas = list(
      response = response_formula,
      binary = binary_formula,
      count = count_formula
    ),
    model_data = list(
      subject = subject_data,
      outcome = model_data,
      positive = positive_data,
      cluster_id = cluster_id,
      robust = robust,
      spltime = spltime
    )
  )
}

make_prediction_rows_conditional <- function(A1, A2R, A2NR, Y0_vec, t_eval, spltime) {
  response_new <- data.frame(A1 = A1, Y0 = Y0_vec)
  responder_new <- data.frame(
    id = seq_along(Y0_vec), Y0 = Y0_vec, A1 = A1, R = 1,
    A2 = ifelse(t_eval > spltime, A2R, 0), time = t_eval
  )
  nonresponder_new <- data.frame(
    id = seq_along(Y0_vec), Y0 = Y0_vec, A1 = A1, R = 0,
    A2 = ifelse(t_eval > spltime, A2NR, 0), time = t_eval
  )
  list(
    response = response_new,
    responder = add_unified_terms_conditional(responder_new, spltime = spltime),
    nonresponder = add_unified_terms_conditional(nonresponder_new, spltime = spltime)
  )
}

predict_conditional_dtr_means <- function(fits, dtr_grid, Y0_vec, t_eval) {
  rows <- vector("list", nrow(dtr_grid))
  
  for (ii in seq_len(nrow(dtr_grid))) {
    A1 <- dtr_grid$A1[ii]
    A2R <- dtr_grid$A2R[ii]
    A2NR <- dtr_grid$A2NR[ii]
    
    pred_rows <- make_prediction_rows_conditional(A1, A2R, A2NR, Y0_vec, t_eval, fits$model_data$spltime)
    r_hat <- stats::predict(fits$response, newdata = pred_rows$response, type = "response")
    
    p_R_hat <- stats::predict(fits$binary, newdata = pred_rows$responder, type = "response")
    p_NR_hat <- stats::predict(fits$binary, newdata = pred_rows$nonresponder, type = "response")
    q_R_hat <- 1 - p_R_hat
    q_NR_hat <- 1 - p_NR_hat
    m_R_hat <- predict_ztp_glm_conditional(fits$count, newdata = pred_rows$responder, type = "m_plus")
    m_NR_hat <- predict_ztp_glm_conditional(fits$count, newdata = pred_rows$nonresponder, type = "m_plus")
    
    pi_plus_i <- r_hat * q_R_hat + (1 - r_hat) * q_NR_hat
    p_zero_i <- r_hat * p_R_hat + (1 - r_hat) * p_NR_hat
    mu_i <- r_hat * q_R_hat * m_R_hat + (1 - r_hat) * q_NR_hat * m_NR_hat
    pi_plus <- mean(pi_plus_i, na.rm = TRUE)
    p_zero <- mean(p_zero_i, na.rm = TRUE)
    mu_total <- mean(mu_i, na.rm = TRUE)
    
    rows[[ii]] <- data.frame(
      A1 = A1,
      A2R = A2R,
      A2NR = A2NR,
      DTR = paste(A1, A2R, A2NR, sep = ","),
      time = t_eval,
      p_zero = p_zero,
      m_plus = ifelse(pi_plus > 0, mu_total / pi_plus, NA_real_)
    )
  }
  
  do.call(rbind, rows)
}

block_diag_conditional <- function(mats) {
  nr <- sum(vapply(mats, nrow, integer(1)))
  nc <- sum(vapply(mats, ncol, integer(1)))
  out <- matrix(0, nr, nc)
  r_start <- 1L
  c_start <- 1L
  for (mat in mats) {
    rr <- seq.int(r_start, length.out = nrow(mat))
    cc <- seq.int(c_start, length.out = ncol(mat))
    out[rr, cc] <- mat
    r_start <- r_start + nrow(mat)
    c_start <- c_start + ncol(mat)
  }
  out
}

cluster_score_matrix_conditional <- function(score_mat, cluster, all_clusters) {
  out <- matrix(0, nrow = length(all_clusters), ncol = ncol(score_mat))
  colnames(out) <- colnames(score_mat)
  rownames(out) <- as.character(all_clusters)
  
  cluster <- as.character(cluster)
  for (cl in unique(cluster)) {
    idx <- which(cluster == cl)
    out[as.character(cl), ] <- colSums(score_mat[idx, , drop = FALSE])
  }
  out
}

glm_score_matrix_conditional <- function(fit) {
  X <- stats::model.matrix(fit)
  y <- stats::model.response(stats::model.frame(fit))
  mu <- stats::fitted(fit)
  score <- X * as.numeric(y - mu)
  colnames(score) <- names(stats::coef(fit))
  score
}

ztp_score_matrix_conditional <- function(fit) {
  data <- fit$data
  mf <- stats::model.frame(fit$formula, data = data)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  y <- as.numeric(stats::model.response(mf))
  eta <- as.numeric(X %*% fit$coefficients)
  m <- ztp_mean_conditional(exp(eta))
  score <- X * as.numeric(y - m)
  colnames(score) <- names(fit$coefficients)
  score
}

joint_coef_vcov_conditional <- function(fits) {
  if (is.null(fits$model_data)) return(NULL)
  if (!isTRUE(fits$model_data$robust)) return(NULL)
  cluster_id <- fits$model_data$cluster_id
  all_clusters <- sort(unique(as.character(fits$model_data$subject[[cluster_id]])))
  
  coefs <- list(
    response = stats::coef(fits$response),
    binary = stats::coef(fits$binary),
    count = fits$count$coefficients
  )
  bread_inv <- block_diag_conditional(list(
    response = stats::vcov(fits$response),
    binary = stats::vcov(fits$binary),
    count = fits$count$vcov_model
  ))
  
  score_blocks <- list(
    response = cluster_score_matrix_conditional(
      glm_score_matrix_conditional(fits$response),
      fits$model_data$subject[[cluster_id]],
      all_clusters
    ),
    binary = cluster_score_matrix_conditional(
      glm_score_matrix_conditional(fits$binary),
      fits$model_data$outcome[[cluster_id]],
      all_clusters
    ),
    count = cluster_score_matrix_conditional(
      ztp_score_matrix_conditional(fits$count),
      fits$model_data$positive[[cluster_id]],
      all_clusters
    )
  )
  
  coef_vec <- unlist(coefs, use.names = TRUE)
  score_stack <- do.call(cbind, score_blocks)
  colnames(score_stack) <- names(coef_vec)
  dimnames(bread_inv) <- list(names(coef_vec), names(coef_vec))
  
  meat <- t(score_stack) %*% score_stack
  g <- nrow(score_stack)
  p <- ncol(score_stack)
  correction <- if (g > 1 && g > p) {
    (g / (g - 1)) * ((g - 1) / (g - p))
  } else if (g > 1) {
    g / (g - 1)
  } else {
    1
  }
  
  vcov_joint <- correction * bread_inv %*% meat %*% bread_inv
  dimnames(vcov_joint) <- list(names(coef_vec), names(coef_vec))
  
  list(coefs = coefs, coef_vec = coef_vec, vcov = vcov_joint)
}

coef_vcov_conditional <- function(fits) {
  joint <- joint_coef_vcov_conditional(fits)
  if (!is.null(joint)) return(joint)
  
  coefs <- list(
    response = stats::coef(fits$response),
    binary = stats::coef(fits$binary),
    count = fits$count$coefficients
  )
  vcovs <- list(
    response = if (!is.null(fits$response$vcov_robust)) fits$response$vcov_robust else stats::vcov(fits$response),
    binary = if (!is.null(fits$binary$vcov_robust)) fits$binary$vcov_robust else stats::vcov(fits$binary),
    count = fits$count$vcov
  )
  
  coef_vec <- unlist(coefs, use.names = TRUE)
  vcov_mat <- block_diag_conditional(vcovs)
  dimnames(vcov_mat) <- list(names(coef_vec), names(coef_vec))
  
  list(coefs = coefs, coef_vec = coef_vec, vcov = vcov_mat)
}

conditional_dtr_from_coef_vector <- function(coef_vec, templates, sizes) {
  idx <- cumsum(c(0L, sizes))
  b_response <- coef_vec[(idx[1] + 1L):idx[2]]
  b_binary <- coef_vec[(idx[2] + 1L):idx[3]]
  b_count <- coef_vec[(idx[3] + 1L):idx[4]]
  
  r_hat <- stats::plogis(as.numeric(templates$X_response %*% b_response))
  p_R_hat <- stats::plogis(as.numeric(templates$X_binary_R %*% b_binary))
  p_NR_hat <- stats::plogis(as.numeric(templates$X_binary_NR %*% b_binary))
  q_R_hat <- 1 - p_R_hat
  q_NR_hat <- 1 - p_NR_hat
  m_R_hat <- ztp_mean_conditional(exp(as.numeric(templates$X_count_R %*% b_count)))
  m_NR_hat <- ztp_mean_conditional(exp(as.numeric(templates$X_count_NR %*% b_count)))
  
  pi_plus_i <- r_hat * q_R_hat + (1 - r_hat) * q_NR_hat
  p_zero_i <- r_hat * p_R_hat + (1 - r_hat) * p_NR_hat
  mu_i <- r_hat * q_R_hat * m_R_hat + (1 - r_hat) * q_NR_hat * m_NR_hat
  pi_plus <- mean(pi_plus_i, na.rm = TRUE)
  p_zero <- mean(p_zero_i, na.rm = TRUE)
  mu_total <- mean(mu_i, na.rm = TRUE)
  
  c(
    p_zero = p_zero,
    m_plus = ifelse(pi_plus > 0, mu_total / pi_plus, NA_real_)
  )
}

numeric_gradient_conditional <- function(fn, x, eps = 1e-5) {
  base_scale <- pmax(abs(x), 1)
  grad <- matrix(NA_real_, nrow = length(fn(x)), ncol = length(x))
  for (jj in seq_along(x)) {
    step <- eps * base_scale[jj]
    x_hi <- x
    x_lo <- x
    x_hi[jj] <- x_hi[jj] + step
    x_lo[jj] <- x_lo[jj] - step
    grad[, jj] <- (fn(x_hi) - fn(x_lo)) / (2 * step)
  }
  rownames(grad) <- names(fn(x))
  colnames(grad) <- names(x)
  grad
}

predict_conditional_dtr_means_delta <- function(fits, dtr_grid, Y0_vec, t_eval) {
  coef_info <- coef_vcov_conditional(fits)
  coef_vec <- coef_info$coef_vec
  vcov_mat <- coef_info$vcov
  sizes <- vapply(coef_info$coefs, length, integer(1))
  
  rows <- vector("list", nrow(dtr_grid))
  
  for (ii in seq_len(nrow(dtr_grid))) {
    A1 <- dtr_grid$A1[ii]
    A2R <- dtr_grid$A2R[ii]
    A2NR <- dtr_grid$A2NR[ii]
    
    pred_rows <- make_prediction_rows_conditional(A1, A2R, A2NR, Y0_vec, t_eval, fits$model_data$spltime)
    
    templates <- list(
      X_response = stats::model.matrix(delete.response(stats::terms(fits$response)), pred_rows$response),
      X_binary_R = stats::model.matrix(delete.response(stats::terms(fits$binary)), pred_rows$responder),
      X_binary_NR = stats::model.matrix(delete.response(stats::terms(fits$binary)), pred_rows$nonresponder),
      X_count_R = stats::model.matrix(delete.response(fits$count$terms), pred_rows$responder),
      X_count_NR = stats::model.matrix(delete.response(fits$count$terms), pred_rows$nonresponder)
    )
    
    estimate_fn <- function(beta) {
      conditional_dtr_from_coef_vector(beta, templates = templates, sizes = sizes)
    }
    theta_hat <- estimate_fn(coef_vec)
    grad <- numeric_gradient_conditional(estimate_fn, coef_vec)
    var_hat <- diag(grad %*% vcov_mat %*% t(grad))
    se_hat <- sqrt(pmax(var_hat, 0))
    names(se_hat) <- paste0(names(theta_hat), "_se")
    
    rows[[ii]] <- data.frame(
      A1 = A1,
      A2R = A2R,
      A2NR = A2NR,
      DTR = paste(A1, A2R, A2NR, sep = ","),
      time = t_eval,
      p_zero = unname(theta_hat["p_zero"]),
      m_plus = unname(theta_hat["m_plus"]),
      p_zero_se = unname(se_hat["p_zero_se"]),
      m_plus_se = unname(se_hat["m_plus_se"])
    )
  }
  
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
