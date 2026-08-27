##########################################
# Conditional Simulation Performance Plots #
##########################################

biom_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.55),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.65),
      plot.margin = ggplot2::margin(6, 8, 6, 8)
    )
}

save_biom_figure <- function(plot, filename, output_dir, width, height,
                             formats = c("png", "eps"), dpi = 600) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (fmt in formats) {
    path <- file.path(output_dir, paste0(filename, ".", fmt))
    if (fmt == "eps") {
      ggplot2::ggsave(
        path, plot, width = width, height = height, units = "in",
        device = grDevices::cairo_ps, fallback_resolution = dpi
      )
    } else {
      ggplot2::ggsave(path, plot, width = width, height = height,
                      units = "in", dpi = dpi)
    }
  }
  invisible(plot)
}

pretty_corstr <- function(corstr) {
  labels <- c(
    independence = "Independence",
    exchangeable = "Exchangeable",
    ar1 = "AR(1)"
  )
  out <- unname(labels[as.character(corstr)])
  out[is.na(out)] <- tools::toTitleCase(as.character(corstr)[is.na(out)])
  out
}

format_dtr_data <- function(df) {
  dtr_levels <- c(
    "-1,-1,-1", "-1,-1,1", "-1,1,-1", "-1,1,1",
    "1,-1,-1", "1,-1,1", "1,1,-1", "1,1,1"
  )
  dtr_numbers <- paste0("DTR ", seq_along(dtr_levels))
  dtr_number_lookup <- setNames(dtr_numbers, dtr_levels)
  df$DTR <- factor(df$DTR, levels = dtr_levels)
  df$DTR_Number <- factor(
    dtr_number_lookup[as.character(df$DTR)],
    levels = dtr_numbers
  )
  df
}

make_greyscale_values <- function(levels, start = 0.90, end = 0.10) {
  setNames(grDevices::gray.colors(length(levels), start = start, end = end), levels)
}

performance_y_breaks <- function(limits) {
  if (all(is.finite(limits)) && mean(limits) > 50) {
    return(c(92.5, 93.75, 95, 96.25, 97.5))
  }
  pretty(limits, n = 5)
}

plot_conditional_bias_coverage <- function(summary_by_dtr, output_dir,
                                           n = NA_integer_,
                                           niter = NA_integer_,
                                           month = NA_integer_,
                                           message_fun = NULL) {
  if (is.null(message_fun)) {
    message_fun <- function(...) {
      cat(..., "\n", file = stdout())
      invisible(NULL)
    }
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    message_fun("Skipping conditional plots: ggplot2 and/or patchwork unavailable.")
    return(invisible(NULL))
  }
  
  required_columns <- c(
    "scenario_id", "ZI", "corstr", "rho", "DTR",
    "p_zero", "m_plus",
    "p_zero_bias", "m_plus_bias",
    "p_zero_se", "m_plus_se",
    "p_zero_covered", "m_plus_covered",
    "p_zero_emp_sd", "m_plus_emp_sd"
  )
  if (!"p_zero" %in% names(summary_by_dtr) &&
      all(c("pi_plus", "pi_plus_bias", "pi_plus_se",
            "pi_plus_covered", "pi_plus_emp_sd") %in% names(summary_by_dtr))) {
    summary_by_dtr$p_zero <- 1 - summary_by_dtr$pi_plus
    summary_by_dtr$p_zero_bias <- -summary_by_dtr$pi_plus_bias
    summary_by_dtr$p_zero_se <- summary_by_dtr$pi_plus_se
    summary_by_dtr$p_zero_covered <- summary_by_dtr$pi_plus_covered
    summary_by_dtr$p_zero_emp_sd <- summary_by_dtr$pi_plus_emp_sd
  }
  missing_columns <- setdiff(required_columns, names(summary_by_dtr))
  if (length(missing_columns) > 0) {
    stop("Cannot plot conditional results. Missing columns: ",
         paste(missing_columns, collapse = ", "))
  }
  
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  make_metric_data <- function(metric, label) {
    out <- data.frame(
      scenario_id = summary_by_dtr$scenario_id,
      ZI = summary_by_dtr$ZI,
      corstr = summary_by_dtr$corstr,
      rho = summary_by_dtr$rho,
      DTR = summary_by_dtr$DTR,
      estimand = label,
      estimate = summary_by_dtr[[metric]],
      bias = summary_by_dtr[[paste0(metric, "_bias")]],
      se_mean = summary_by_dtr[[paste0(metric, "_se")]],
      coverage = summary_by_dtr[[paste0(metric, "_covered")]] * 100,
      emp_sd = summary_by_dtr[[paste0(metric, "_emp_sd")]],
      row.names = NULL
    )
    if ("scenario_label" %in% names(summary_by_dtr)) {
      out$scenario_label_short <- summary_by_dtr$scenario_label
    }
    if ("scenario_family" %in% names(summary_by_dtr)) {
      out$scenario_family <- summary_by_dtr$scenario_family
    }
    out
  }
  
  plot_data <- do.call(rbind, list(
    make_metric_data("p_zero", "Binary Model"),
    make_metric_data("m_plus", "Truncated Count Model")
  ))
  plot_data <- format_dtr_data(plot_data)
  has_manuscript_labels <- all(c("scenario_label_short", "scenario_family") %in% names(plot_data))
  if (has_manuscript_labels) {
    plot_data$scenario_label <- paste0(
      pretty_corstr(plot_data$corstr),
      ", rho ", plot_data$rho,
      ", ZI ", plot_data$ZI, "%"
    )
    scenario_levels <- unique(
      plot_data$scenario_label[
        order(plot_data$scenario_label_short, plot_data$corstr, plot_data$rho, plot_data$ZI)
      ]
    )
  } else {
    plot_data$scenario_label <- paste0(
      pretty_corstr(plot_data$corstr), ", rho ", plot_data$rho, ", ZI ", plot_data$ZI, "%"
    )
    corstr_order <- c("independence", "exchangeable", "ar1")
    scenario_levels <- unique(
      plot_data$scenario_label[
        order(match(plot_data$corstr, corstr_order), plot_data$rho, plot_data$ZI)
      ]
    )
  }
  plot_data$scenario_label <- factor(plot_data$scenario_label,
                                     levels = scenario_levels)
  if (has_manuscript_labels) {
    scenario_colors <- c(
      "Exchangeable, rho 0.3, ZI 37%" = "grey80",
      "AR(1), rho 0.3, ZI 37%" = "grey60",
      "Exchangeable, rho 0.5, ZI 37%" = "grey35",
      "AR(1), rho 0.5, ZI 37%" = "black",
      "Independence, rho 0, ZI 37%" = "white"
    )
  } else {
    scenario_colors <- make_greyscale_values(levels(plot_data$scenario_label), start = 0.95, end = 0.05)
  }
  missing_color_levels <- setdiff(levels(plot_data$scenario_label), names(scenario_colors))
  if (length(missing_color_levels) > 0) {
    scenario_colors <- c(
      scenario_colors,
      setNames(
        grDevices::gray.colors(length(missing_color_levels), start = 0.95, end = 0.05),
        missing_color_levels
      )
    )
  }
  scenario_colors <- scenario_colors[levels(plot_data$scenario_label)]
  plot_data$estimand <- factor(
    plot_data$estimand,
    levels = c("Binary Model", "Truncated Count Model")
  )
  utils::write.csv(
    plot_data,
    file.path(plot_dir, "conditional_bias_coverage_plot_data.csv"),
    row.names = FALSE
  )
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
  })
  
  fill_levels <- levels(plot_data$scenario_label)
  if (length(fill_levels) == 5) {
    scenario_colors <- setNames(c("white", "grey75", "grey50", "grey25", "black"), fill_levels)
  } else {
    scenario_colors <- make_greyscale_values(fill_levels, start = 1, end = 0)
  }
  
  sample_size <- if (is.finite(n)) n else NA_integer_
  bias_limit <- if (is.finite(sample_size) && isTRUE(sample_size == 400)) {
    0.035
  } else if (is.finite(sample_size) && isTRUE(sample_size == 1000)) {
    0.03
  } else {
    max_abs_bias <- max(abs(plot_data$bias), na.rm = TRUE)
    if (!is.finite(max_abs_bias) || max_abs_bias == 0) {
      max_abs_bias <- .Machine$double.eps
    }
    ceiling(max_abs_bias * 1000) / 1000
  }
  bias_limits <- c(-bias_limit, bias_limit)
  coverage_limits <- c(92.5, 97.5)
  
  performance_data <- rbind(
    data.frame(
      DTR_Number = plot_data$DTR_Number,
      estimand = plot_data$estimand,
      scenario_label = plot_data$scenario_label,
      metric = "Mean Bias",
      value = plot_data$bias,
      row.names = NULL
    ),
    data.frame(
      DTR_Number = plot_data$DTR_Number,
      estimand = plot_data$estimand,
      scenario_label = plot_data$scenario_label,
      metric = "Coverage (%)",
      value = plot_data$coverage,
      row.names = NULL
    )
  )
  performance_data$metric <- factor(
    performance_data$metric,
    levels = c("Mean Bias", "Coverage (%)")
  )
  performance_data$DTR_Index <- as.numeric(performance_data$DTR_Number)
  performance_data$Scenario_Index <- as.numeric(performance_data$scenario_label)
  performance_data$Bar_Width <- 0.12
  performance_data$Bar_X <- performance_data$DTR_Index +
    (performance_data$Scenario_Index - mean(seq_along(fill_levels))) * performance_data$Bar_Width
  performance_data$Bar_Xmin <- performance_data$Bar_X - performance_data$Bar_Width / 2
  performance_data$Bar_Xmax <- performance_data$Bar_X + performance_data$Bar_Width / 2
  blank_data <- rbind(
    expand.grid(
      DTR_Index = seq_along(levels(plot_data$DTR_Number)),
      estimand = levels(plot_data$estimand),
      metric = factor("Mean Bias", levels = levels(performance_data$metric)),
      value = bias_limits
    ),
    expand.grid(
      DTR_Index = seq_along(levels(plot_data$DTR_Number)),
      estimand = levels(plot_data$estimand),
      metric = factor("Coverage (%)", levels = levels(performance_data$metric)),
      value = coverage_limits
    )
  )
  
  month_label <- if (is.finite(month)) as.character(month) else "__"
  n_title <- if (is.finite(n)) as.character(n) else "__"
  niter_title <- if (is.finite(niter)) as.character(niter) else "__"
  plot_title <- paste0(
    "DTR performance at Month ", month_label,
    " | n = ", n_title,
    " | niter = ", niter_title
  )
  
  combined_plot <- ggplot2::ggplot() +
    ggplot2::geom_blank(
      data = blank_data,
      ggplot2::aes(x = DTR_Index, y = value)
    ) +
    ggplot2::geom_hline(
      data = data.frame(metric = factor("Mean Bias", levels = levels(performance_data$metric))),
      ggplot2::aes(yintercept = 0),
      inherit.aes = FALSE,
      linewidth = 0.55,
      color = "black"
    ) +
    ggplot2::geom_hline(
      data = data.frame(metric = factor("Coverage (%)", levels = levels(performance_data$metric))),
      ggplot2::aes(yintercept = 95),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.55,
      color = "black"
    ) +
    ggplot2::geom_rect(
      data = performance_data[performance_data$metric == "Mean Bias", ],
      ggplot2::aes(
        xmin = Bar_Xmin,
        xmax = Bar_Xmax,
        ymin = 0,
        ymax = value,
        fill = scenario_label
      ),
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::geom_rect(
      data = performance_data[performance_data$metric == "Coverage (%)", ],
      ggplot2::aes(
        xmin = Bar_Xmin,
        xmax = Bar_Xmax,
        ymin = coverage_limits[1],
        ymax = value,
        fill = scenario_label
      ),
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::facet_grid(metric ~ estimand, scales = "free_y", switch = "y") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(levels(plot_data$DTR_Number)),
      labels = levels(plot_data$DTR_Number)
    ) +
    ggplot2::scale_y_continuous(breaks = performance_y_breaks) +
    ggplot2::scale_fill_manual(values = scenario_colors, drop = FALSE, name = "Simulation Scenario") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
    biom_theme(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black", size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5),
      axis.text.y = ggplot2::element_text(size = 9.5),
      strip.text = ggplot2::element_text(face = "bold", color = "black", size = 10),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 90, size = 10),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 9),
      legend.text = ggplot2::element_text(size = 7.5),
      legend.key.size = ggplot2::unit(0.18, "in")
    )
  
  save_biom_figure(
    combined_plot,
    filename = "conditional_bias_coverage",
    output_dir = plot_dir,
    width = 10,
    height = 6.8
  )
  
  message_fun("Conditional performance plots saved to: ", plot_dir)
  invisible(plot_data)
}

plot_response_stress_coverage <- function(summary_by_dtr, output_dir,
                                          niter = NA_integer_,
                                          month = NA_integer_,
                                          message_fun = NULL) {
  if (is.null(message_fun)) {
    message_fun <- function(...) {
      cat(..., "\n", file = stdout())
      invisible(NULL)
    }
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    message_fun("Skipping stress bias/coverage plot: ggplot2 and/or patchwork unavailable.")
    return(invisible(NULL))
  }
  
  required_columns <- c(
    "sample_size", "ZI", "response_profile", "DTR",
    "p_zero_bias", "m_plus_bias",
    "p_zero_covered", "m_plus_covered"
  )
  missing_columns <- setdiff(required_columns, names(summary_by_dtr))
  if (length(missing_columns) > 0) {
    stop("Cannot plot response stress bias/coverage. Missing columns: ",
         paste(missing_columns, collapse = ", "))
  }
  
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  make_component_data <- function(bias_metric, coverage_metric, label) {
    data.frame(
      sample_size = summary_by_dtr$sample_size,
      ZI = summary_by_dtr$ZI,
      response_profile = summary_by_dtr$response_profile,
      DTR = summary_by_dtr$DTR,
      estimand = label,
      bias = summary_by_dtr[[bias_metric]],
      coverage = 100 * summary_by_dtr[[coverage_metric]],
      row.names = NULL
    )
  }
  plot_data <- rbind(
    make_component_data("p_zero_bias", "p_zero_covered", "Binary Model"),
    make_component_data("m_plus_bias", "m_plus_covered", "Truncated Count Model")
  )
  
  response_labels <- c(
    "r0_0.5_r1_0.5" = "r0 = 0.5, r1 = 0.5",
    "r0_0.6_r1_0.4" = "r0 = 0.6, r1 = 0.4",
    "r0_0.4_r1_0.6" = "r0 = 0.4, r1 = 0.6",
    "truth_A1_Y0_BLDEPPOS_fit_A1_Y0" = "Response misspecification"
  )
  plot_data$response_label <- unname(response_labels[plot_data$response_profile])
  plot_data$response_label[is.na(plot_data$response_label)] <-
    plot_data$response_profile[is.na(plot_data$response_label)]
  
  response_order <- c(
    "r0 = 0.5, r1 = 0.5",
    "r0 = 0.6, r1 = 0.4",
    "r0 = 0.4, r1 = 0.6",
    "Response misspecification"
  )
  response_order <- response_order[response_order %in% unique(plot_data$response_label)]
  response_order <- c(response_order, setdiff(unique(plot_data$response_label), response_order))
  plot_data$response_label <- factor(plot_data$response_label, levels = response_order)
  plot_data$estimand <- factor(
    plot_data$estimand,
    levels = c("Binary Model", "Truncated Count Model")
  )
  plot_data <- format_dtr_data(plot_data)
  plot_data$ZI_label <- factor(
    paste0(plot_data$ZI, "%"),
    levels = paste0(sort(unique(plot_data$ZI)), "%")
  )
  plot_data$sample_size_label <- factor(
    paste0("n = ", plot_data$sample_size),
    levels = paste0("n = ", sort(unique(plot_data$sample_size)))
  )
  plot_data$scenario_label <- paste0("ZI ", plot_data$ZI_label, "; ", plot_data$response_label)
  scenario_order <- unique(plot_data$scenario_label[
    order(plot_data$response_label, plot_data$ZI)
  ])
  plot_data$scenario_label <- factor(plot_data$scenario_label, levels = scenario_order)
  fill_levels <- levels(plot_data$scenario_label)
  utils::write.csv(
    plot_data,
    file.path(plot_dir, "response_stress_bias_coverage_plot_data.csv"),
    row.names = FALSE
  )
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
  })
  
  if (length(fill_levels) == 1) {
    fill_values <- "black"
  } else if (length(fill_levels) == 6) {
    fill_values <- c("white", grDevices::gray(c(0.6, 0.2, 0.8, 0.4)), "black")
  } else {
    fill_values <- grDevices::gray.colors(length(fill_levels), start = 1, end = 0)
  }
  fill_values <- setNames(fill_values, fill_levels)
  
  sample_sizes <- sort(unique(plot_data$sample_size))
  bias_limit <- if (length(fill_levels) == 1 &&
                    identical(as.character(fill_levels), "ZI 37%; Response misspecification")) {
    0.02
  } else if (length(sample_sizes) == 1 && isTRUE(sample_sizes == 400)) {
    0.035
  } else if (length(sample_sizes) == 1 && isTRUE(sample_sizes == 1000)) {
    0.03
  } else {
    max_abs_bias <- max(abs(plot_data$bias), na.rm = TRUE)
    if (!is.finite(max_abs_bias) || max_abs_bias == 0) {
      max_abs_bias <- .Machine$double.eps
    }
    ceiling(max_abs_bias * 1000) / 1000
  }
  bias_limits <- c(-bias_limit, bias_limit)
  coverage_limits <- c(92.5, 97.5)
  
  performance_data <- rbind(
    data.frame(
      DTR_Number = plot_data$DTR_Number,
      estimand = plot_data$estimand,
      scenario_label = plot_data$scenario_label,
      metric = "Mean Bias",
      value = plot_data$bias,
      row.names = NULL
    ),
    data.frame(
      DTR_Number = plot_data$DTR_Number,
      estimand = plot_data$estimand,
      scenario_label = plot_data$scenario_label,
      metric = "Coverage (%)",
      value = plot_data$coverage,
      row.names = NULL
    )
  )
  performance_data$metric <- factor(
    performance_data$metric,
    levels = c("Mean Bias", "Coverage (%)")
  )
  performance_data$DTR_Index <- as.numeric(performance_data$DTR_Number)
  performance_data$Scenario_Index <- as.numeric(performance_data$scenario_label)
  performance_data$Bar_Width <- 0.12
  performance_data$Bar_X <- performance_data$DTR_Index +
    (performance_data$Scenario_Index - mean(seq_along(fill_levels))) * performance_data$Bar_Width
  performance_data$Bar_Xmin <- performance_data$Bar_X - performance_data$Bar_Width / 2
  performance_data$Bar_Xmax <- performance_data$Bar_X + performance_data$Bar_Width / 2
  blank_data <- rbind(
    expand.grid(
      DTR_Index = seq_along(levels(plot_data$DTR_Number)),
      estimand = levels(plot_data$estimand),
      metric = factor("Mean Bias", levels = levels(performance_data$metric)),
      value = bias_limits
    ),
    expand.grid(
      DTR_Index = seq_along(levels(plot_data$DTR_Number)),
      estimand = levels(plot_data$estimand),
      metric = factor("Coverage (%)", levels = levels(performance_data$metric)),
      value = coverage_limits
    )
  )
  
  month_label <- if (is.finite(month)) as.character(month) else "__"
  sample_size_title <- if (length(sample_sizes) == 1) {
    as.character(sample_sizes)
  } else {
    paste(sample_sizes, collapse = ", ")
  }
  niter_title <- if (is.finite(niter)) as.character(niter) else "__"
  plot_title <- paste0(
    "DTR performance at Month ", month_label,
    " | n = ", sample_size_title,
    " | niter = ", niter_title
  )
  
  combined_plot <- ggplot2::ggplot() +
    ggplot2::geom_blank(
      data = blank_data,
      ggplot2::aes(x = DTR_Index, y = value)
    ) +
    ggplot2::geom_hline(
      data = data.frame(metric = factor("Mean Bias", levels = levels(performance_data$metric))),
      ggplot2::aes(yintercept = 0),
      inherit.aes = FALSE,
      linewidth = 0.55,
      color = "black"
    ) +
    ggplot2::geom_hline(
      data = data.frame(metric = factor("Coverage (%)", levels = levels(performance_data$metric))),
      ggplot2::aes(yintercept = 95),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.55,
      color = "black"
    ) +
    ggplot2::geom_rect(
      data = performance_data[performance_data$metric == "Mean Bias", ],
      ggplot2::aes(
        xmin = Bar_Xmin,
        xmax = Bar_Xmax,
        ymin = 0,
        ymax = value,
        fill = scenario_label
      ),
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::geom_rect(
      data = performance_data[performance_data$metric == "Coverage (%)", ],
      ggplot2::aes(
        xmin = Bar_Xmin,
        xmax = Bar_Xmax,
        ymin = coverage_limits[1],
        ymax = value,
        fill = scenario_label
      ),
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::facet_grid(metric ~ estimand, scales = "free_y", switch = "y") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(levels(plot_data$DTR_Number)),
      labels = levels(plot_data$DTR_Number)
    ) +
    ggplot2::scale_y_continuous(breaks = performance_y_breaks) +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE, name = "Simulation Scenario") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
    biom_theme(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black", size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5),
      axis.text.y = ggplot2::element_text(size = 9.5),
      strip.text = ggplot2::element_text(face = "bold", color = "black", size = 10),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 90, size = 10),
      legend.position = if (length(fill_levels) == 1) "none" else "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 9),
      legend.text = ggplot2::element_text(size = 7.5),
      legend.key.size = ggplot2::unit(0.18, "in")
    )
  
  save_biom_figure(
    combined_plot,
    filename = "response_stress_bias_coverage",
    output_dir = plot_dir,
    width = 10,
    height = 6.9
  )
  
  message_fun("Response stress bias/coverage plot saved to: ", plot_dir)
  invisible(plot_data)
}

plot_zero_prediction_diagnostics <- function(zero_summary, output_dir,
                                             month = 4,
                                             message_fun = NULL) {
  if (is.null(message_fun)) {
    message_fun <- function(...) {
      cat(..., "\n", file = stdout())
      invisible(NULL)
    }
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message_fun("Skipping zero-prediction plot: ggplot2 unavailable.")
    return(invisible(NULL))
  }
  
  required_columns <- c(
    "sample_size", "ZI", "time", "observed_zero", "predicted_zero_hm",
    "predicted_zero_all_poisson", "predicted_zero_all_nb"
  )
  missing_columns <- setdiff(required_columns, names(zero_summary))
  if (length(missing_columns) > 0) {
    stop("Cannot plot zero-prediction diagnostics. Missing columns: ",
         paste(missing_columns, collapse = ", "))
  }
  
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_time <- as.character(month)
  plot_data_wide <- zero_summary[as.character(zero_summary$time) == plot_time, , drop = FALSE]
  if (nrow(plot_data_wide) == 0) {
    stop("Cannot plot zero-prediction diagnostics. No rows found for time = ", plot_time)
  }
  
  model_columns <- c(
    observed_zero = "Observed",
    predicted_zero_hm = "Hurdle Model",
    predicted_zero_all_poisson = "All-Count Poisson",
    predicted_zero_all_nb = "All-Count NB"
  )
  plot_data <- do.call(rbind, lapply(names(model_columns), function(column_name) {
    data.frame(
      sample_size = plot_data_wide$sample_size,
      ZI = plot_data_wide$ZI,
      model = unname(model_columns[column_name]),
      zero_probability = plot_data_wide[[column_name]],
      row.names = NULL
    )
  }))
  model_levels <- unname(model_columns)
  plot_data$model <- factor(plot_data$model, levels = model_levels)
  plot_data$sample_size_label <- factor(
    paste0("n = ", plot_data$sample_size),
    levels = paste0("n = ", sort(unique(plot_data$sample_size)))
  )
  utils::write.csv(
    plot_data,
    file.path(plot_dir, "zero_prediction_diagnostics_plot_data.csv"),
    row.names = FALSE
  )
  
  suppressPackageStartupMessages(library(ggplot2))
  
  model_values <- setNames(
    c("black", "grey35", "grey60", "grey80")[seq_along(model_levels)],
    model_levels
  )
  model_shapes <- setNames(c(16, 1, 17, 2)[seq_along(model_levels)], model_levels)
  model_linetypes <- setNames(c("solid", "dashed", "dotdash", "dotted")[
    seq_along(model_levels)
  ], model_levels)
  
  zero_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = ZI,
      y = zero_probability,
      group = model,
      color = model,
      shape = model,
      linetype = model
    )
  ) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(size = 2.4, fill = "white", stroke = 0.8) +
    ggplot2::facet_wrap(~ sample_size_label, nrow = 1) +
    ggplot2::scale_color_manual(values = model_values, drop = FALSE) +
    ggplot2::scale_shape_manual(values = model_shapes, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = model_linetypes, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(plot_data$ZI))) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Zero Inflation (%)",
      y = "Predicted/Observed Zero Proportion, Log Scale",
      color = "Model",
      shape = "Model",
      linetype = "Model"
    ) +
    biom_theme(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
  
  save_biom_figure(
    zero_plot,
    filename = "zero_prediction_diagnostics",
    output_dir = plot_dir,
    width = 7.25,
    height = 4.6
  )
  
  message_fun("Zero-prediction diagnostic plot saved to: ", plot_dir)
  invisible(plot_data)
}

plot_conditional_results_from_dir <- function(
    output_dir = file.path("results", "primary_simulation", "niter_100"),
    n = NA_integer_,
    niter = NA_integer_,
    month = NA_integer_,
    summary_filename = NULL) {
  if (is.null(summary_filename)) {
    candidates <- c(
      "bias_coverage_summary_conditional.csv",
      "bias_coverage_summary_manuscript_primary.csv"
    )
    summary_path <- candidates[file.exists(file.path(output_dir, candidates))][1]
    summary_path <- file.path(output_dir, summary_path)
  } else {
    summary_path <- file.path(output_dir, summary_filename)
  }
  if (!file.exists(summary_path)) {
    stop("Cannot find conditional summary file: ", summary_path)
  }
  
  summary_by_dtr <- utils::read.csv(summary_path)
  plot_conditional_bias_coverage(
    summary_by_dtr = summary_by_dtr,
    output_dir = output_dir,
    n = n,
    niter = niter,
    month = month
  )
}

plot_stress_results_from_dir <- function(output_dir,
                                         niter = NA_integer_,
                                         month = 4) {
  response_summary_path <- file.path(output_dir, "bias_coverage_summary_response_stress.csv")
  zero_summary_path <- file.path(output_dir, "zero_prediction_summary_stress_tests.csv")
  
  out <- list(response_stress_coverage = NULL, zero_prediction = NULL)
  if (file.exists(response_summary_path)) {
    response_summary <- utils::read.csv(response_summary_path)
    out$response_stress_coverage <- plot_response_stress_coverage(
      summary_by_dtr = response_summary,
      output_dir = output_dir,
      niter = niter,
      month = month
    )
  }
  if (file.exists(zero_summary_path)) {
    zero_summary <- utils::read.csv(zero_summary_path)
    out$zero_prediction <- plot_zero_prediction_diagnostics(
      zero_summary = zero_summary,
      output_dir = output_dir,
      month = month
    )
  }
  invisible(out)
}
