###############################################################
# Conditional Hurdle Model Data Application: Y0 Baseline Case #
###############################################################
#
# Parallel to HM_Data_Application.R, but using the conditional-first analysis:
#   1. Fit R | A1, Y0.
#   2. Fit a binary model for B = 1(Y = 0) with shared pre-split terms and
#      branch-specific post-split terms.
#   3. Fit a zero-truncated count model for Y | Y > 0 with the same structure.
#   4. Estimate DTR means by g-computation with cluster-robust delta-method SEs.

if (!("pacman" %in% installed.packages()[, "Package"])) install.packages("pacman")
library(pacman)
p_load(tidyverse, patchwork, scales)

source("GEE_Code/HM/Y0_Baseline/functions_HM_conditional.R")
source("GEE_Code/HM/Y0_Baseline/fit_gcomp_conditional.R")

immediate_cat <- function(...) {
  cat(..., "\n", file = stdout())
  try(flush.console(), silent = TRUE)
  invisible(NULL)
}

output_dir_data <- "GEE_Code/HM/Y0_Baseline/Data_Results_Cond_New_Surrogate"
dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir_data)) {
  stop("Failed to create output directory: ", output_dir_data)
}

output_dir_figures <- file.path(output_dir_data, "figures")
dir.create(output_dir_figures, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir_figures)) {
  stop("Failed to create output figure directory: ", output_dir_figures)
}

weekly_path <- Sys.getenv(
  "MCOACH_WEEKLY_PATH",
  unset = ""
)
outcome_path <- Sys.getenv(
  "MCOACH_OUTCOME_PATH",
  unset = ""
)

stop_if_missing <- function(path, label) {
  if (!nzchar(path) || !file.exists(path)) {
    stop(label, " file not found: ", path,
         "\nSet the corresponding environment variable if the path changed.")
  }
}

figure_dpi <- 600

save_data_app_figure <- function(plot, filename, width, height, dpi = figure_dpi) {
  ggsave(
    file.path(output_dir_figures, paste0(filename, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  ggsave(
    file.path(output_dir_figures, paste0(filename, ".eps")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_ps,
    fallback_resolution = dpi
  )
  invisible(plot)
}

biom_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.55),
      strip.text = element_text(face = "bold", color = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.65),
      plot.margin = margin(6, 8, 6, 8)
    )
}

##### Data Cleaning #####
stop_if_missing(weekly_path, "Weekly")
stop_if_missing(outcome_path, "Outcome")

data_weekly <- read.csv(weekly_path, header = TRUE)
data_outcome <- read.csv(outcome_path, header = TRUE)

data_outcome1 <- data_outcome %>%
  mutate(
    A1 = FIRSTRANDOM,
    A2 = T2,
    Y0 = BLALCCONSUMPTION,
    Y4 = FL4MALCCONSUMPTION,
    R = RESPONDENT,
    id = PID_SC
  ) %>%
  dplyr::select(id, A1, A2, R, Y0, Y4) %>%
  filter(Y0 < 600, Y4 < 600)

data_weekly1 <- data_weekly %>%
  dplyr::select(
    PID_SC,
    AUD1_WK1, AUD2_WK1,
    AUD1_WK2, AUD2_WK2,
    AUD1_WK3, AUD2_WK3,
    AUD1_WK4, AUD2_WK4
  ) %>%
  mutate(
    id = PID_SC,
    wk1_drinks = dplyr::if_else(AUD1_WK1 == 0, 0, AUD1_WK1 * AUD2_WK1),
    wk2_drinks = dplyr::if_else(AUD1_WK2 == 0, 0, AUD1_WK2 * AUD2_WK2),
    wk3_drinks = dplyr::if_else(AUD1_WK3 == 0, 0, AUD1_WK3 * AUD2_WK3),
    wk4_drinks = dplyr::if_else(AUD1_WK4 == 0, 0, AUD1_WK4 * AUD2_WK4)
  ) %>%
  filter(!if_all(starts_with("wk"), is.na)) %>%
  mutate(
    Y1_weekly_mean = rowMeans(across(wk1_drinks:wk4_drinks), na.rm = TRUE),
    Y1_mo = Y1_weekly_mean * 30 / 7,
    Y1 = round(Y1_mo, 0)
  ) %>%
  dplyr::select(id, Y1)

data_weekly2 <- data_weekly %>%
  dplyr::select(
    PID_SC,
    AUD1_WK5, AUD2_WK5,
    AUD1_WK6, AUD2_WK6,
    AUD1_WK7, AUD2_WK7,
    AUD1_WK8, AUD2_WK8
  ) %>%
  mutate(
    id = PID_SC,
    wk5_drinks = dplyr::if_else(AUD1_WK5 == 0, 0, AUD1_WK5 * AUD2_WK5),
    wk6_drinks = dplyr::if_else(AUD1_WK6 == 0, 0, AUD1_WK6 * AUD2_WK6),
    wk7_drinks = dplyr::if_else(AUD1_WK7 == 0, 0, AUD1_WK7 * AUD2_WK7),
    wk8_drinks = dplyr::if_else(AUD1_WK8 == 0, 0, AUD1_WK8 * AUD2_WK8)
  ) %>%
  filter(!if_all(starts_with("wk"), is.na)) %>%
  mutate(
    Y2_weekly_mean = rowMeans(across(wk5_drinks:wk8_drinks), na.rm = TRUE),
    Y2_mo = Y2_weekly_mean * 30 / 7,
    Y2 = round(Y2_mo, 0)
  ) %>%
  dplyr::select(id, Y2)

DataWeeklyFull <- merge(data_weekly1, data_weekly2, by = "id")
# DataWideFormat <- merge(data_outcome2, DataWeeklyFull, by = "id")
DataWideFormat <- merge(data_outcome1, DataWeeklyFull, by = "id")

if (!all(stats::na.omit(DataWideFormat$A1) %in% c(-1, 1))) {
  stop("A1 must be coded as -1/1 before running the conditional HM data application.")
}
if (!all(stats::na.omit(DataWideFormat$A2) %in% c(-1, 1))) {
  stop("A2 must be coded as -1/1 before running the conditional HM data application.")
}
if (!all(stats::na.omit(DataWideFormat$R) %in% c(0, 1))) {
  stop("R must be coded as 0/1 before running the conditional HM data application.")
}
if (anyNA(DataWideFormat[, c("Y0", "Y1", "Y2", "Y4")])) {
  stop("Outcome data contain missing Y0/Y1/Y2/Y4 values after merging.")
}

##### Diagnostics #####
DataLongAll <- DataWideFormat %>%
  pivot_longer(
    cols = starts_with("Y"),
    names_to = "time",
    names_prefix = "Y",
    values_to = "Y"
  ) %>%
  mutate(time = as.numeric(time)) %>%
  arrange(id, time) %>%
  filter(Y < 600)

proportions_by_time <- DataLongAll %>%
  mutate(Ybin_zero = as.integer(Y == 0)) %>%
  group_by(time) %>%
  summarise(
    prop_zero = mean(Ybin_zero == 1),
    prop_positive = mean(Ybin_zero == 0),
    .groups = "drop"
  )

boxplot_outcome <- ggplot(DataLongAll, aes(x = factor(time), y = Y)) +
  geom_boxplot(fill = "gray90", outlier.size = 0.8) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 2) +
  labs(
    x = "Time (Months)",
    y = "Number of Drinks Consumed Over the Last 30 Days",
    title = "Distribution of Drinks Over Time"
  ) +
  theme_minimal()

time_colors <- c(
  "Month 0" = "#A2E1A6",
  "Month 1" = "#FCD366",
  "Month 2" = "#7FB2E5",
  "Month 4" = "#0A2540"
)

density_plot <- DataLongAll %>%
  mutate(month_label = paste0("Month ", time)) %>%
  ggplot(aes(x = Y, color = month_label, fill = month_label)) +
  geom_density(alpha = 0.3, position = "identity", linewidth = 0.6) + 
  scale_fill_manual(values = time_colors) +
  scale_color_manual(values = time_colors) + 
  labs(
    title = "Distribution of Alcoholic Drinks Consumed Over the Last 30 Days",
    x = "Number of Alcoholic Drinks Consumed",
    y = "Density",
    color = "Time",
    fill = "Time"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15)
  )

# density_plot <- DataLongAll %>%
#   mutate(month_label = paste0("Month ", time)) %>%
#   ggplot(aes(x = Y, color = month_label, fill = month_label)) +
#   geom_density(alpha = 0.4, position = "identity") +
#   labs(
#     title = "Distribution of Alcoholic Drinks Consumed Over the Last 30 Days",
#     x = "Number of Alcoholic Drinks Consumed",
#     y = "Density",
#     color = "Time",
#     fill = "Time"
#   ) +
#   theme_minimal(base_size = 16)

##### Conditional Analysis Dataset #####
spltime <- 1
times <- c(1, 2, 4)
t_eval <- 4

DataLongConditional <- DataWideFormat %>%
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
    A2 = ifelse(time <= spltime, 0, A2),
    A2R = ifelse(R == 1, A2, 0),
    A2NR = ifelse(R == 0, A2, 0),
    B = as.integer(Y == 0)
  ) %>%
  arrange(id, time)

##### Fit Conditional Hurdle Models #####
fits <- tryCatch({
  fit_conditional_hurdle_models(
    long_data = DataLongConditional,
    spltime = spltime,
    response_formula = R ~ A1 + Y0,
    cluster_id = "id",
    robust = TRUE
  )
}, error = function(e) {
  immediate_cat("ERROR: conditional hurdle fit failed: ", e$message)
  NULL
})
if (is.null(fits)) {
  stop("Conditional hurdle fit failed; cannot compute DTR estimates.")
}

##### DTR Settings #####
dtr_levels <- c(
  "-1,-1,-1", "-1,-1,1", "-1,1,-1", "-1,1,1",
  "1,-1,-1", "1,-1,1", "1,1,-1", "1,1,1"
)
# Intervention labels follow Table 1: A2R = -1 means continue the
# first-stage treatment; A2R = 1 means step down to resource brochure.
dtr_labels_intervention <- c(
  "TM,TM,TM",
  "TM,TM,HC",
  "TM,RB,TM",
  "TM,RB,HC",
  "HC,HC,HC",
  "HC,HC,HC+",
  "HC,RB,HC",
  "HC,RB,HC+"
)
dtr_labels_number <- paste0("DTR ", seq_along(dtr_levels))

dtr_grid <- tibble(DTR = dtr_levels) %>%
  separate(DTR, into = c("A1", "A2R", "A2NR"), sep = ",",
           convert = TRUE, remove = FALSE)

dtr_lookup <- tibble(
  DTR = dtr_levels,
  DTR_label_int = dtr_labels_intervention,
  dtr_labels_number = dtr_labels_number,
  DTR_label = paste(dtr_labels_number, dtr_labels_intervention, sep = ": ")
)

expected_dtr_lookup <- tibble(
  DTR = c("-1,-1,-1", "-1,-1,1", "-1,1,-1", "-1,1,1",
          "1,-1,-1", "1,-1,1", "1,1,-1", "1,1,1"),
  DTR_label_int = c("TM,TM,TM", "TM,TM,HC", "TM,RB,TM", "TM,RB,HC",
                    "HC,HC,HC", "HC,HC,HC+", "HC,RB,HC", "HC,RB,HC+")
)
if (!identical(
  dtr_lookup %>% dplyr::select(DTR, DTR_label_int),
  expected_dtr_lookup
)) {
  stop("DTR labels do not match Table 1 definitions.")
}

##### DTR Estimates #####
dtr_results <- predict_conditional_dtr_means_delta(
  fits = fits,
  dtr_grid = dtr_grid %>% dplyr::select(A1, A2R, A2NR),
  Y0_vec = DataWideFormat$Y0,
  t_eval = t_eval
) %>%
  mutate(
    p_zero_lower = p_zero - 1.96 * p_zero_se,
    p_zero_upper = p_zero + 1.96 * p_zero_se,
    m_plus_lower = m_plus - 1.96 * m_plus_se,
    m_plus_upper = m_plus + 1.96 * m_plus_se
  ) %>%
  left_join(dtr_lookup, by = c("DTR"))

dtr_model_results <- bind_rows(
  dtr_results %>%
    transmute(
      model = "Binary Model: P(Y=0)",
      outcome_label = "P(Y = 0)",
      DTR, A1, A2R, A2NR, dtr_labels_number, DTR_label_int, DTR_label,
      Estimate = p_zero,
      Lower = p_zero_lower,
      Upper = p_zero_upper
    ),
  dtr_results %>%
    transmute(
      model = "Truncated Count Model: E[Y|Y>0]",
      outcome_label = "E[Y | Y > 0]",
      DTR, A1, A2R, A2NR, dtr_labels_number, DTR_label_int, DTR_label,
      Estimate = m_plus,
      Lower = m_plus_lower,
      Upper = m_plus_upper
    )
) %>%
  mutate(
    model = factor(model, levels = c("Binary Model: P(Y=0)", "Truncated Count Model: E[Y|Y>0]")),
    DTR_label = factor(DTR_label, levels = dtr_lookup$DTR_label),
    dtr_plot_label = factor(DTR_label, levels = rev(dtr_lookup$DTR_label))
  )

forest_plot <- ggplot(
  dtr_model_results,
  aes(x = Estimate, y = dtr_plot_label, color = model)
) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), orientation = "y", width = 0.25) +
  facet_wrap(~ model, scales = "free_x", nrow = 1) +
  scale_color_manual(
    values = c(
      "Binary Model: P(Y=0)" = "#002B5B",        # Deep Navy
      "Truncated Count Model: E[Y|Y>0]" = "#B30006" # Crimson
    ),
    guide = "none"
  ) +
  labs(
    title = paste0("Conditional G-Computation DTR Estimates at Month ", t_eval),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "#000000", margin = margin(b = 12)),
    strip.text = element_text(size = 13, face = "bold", color = "#000000", margin = margin(b = 8)),
    axis.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    axis.text.x = element_text(size = 10, color = "#000000"),
    panel.grid.major.x = element_line(color = "#EAEAEA", linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.spacing = unit(1.5, "lines")
  )

# presentation version
forest_plot_ppt <- ggplot(
  dtr_model_results,
  aes(x = Estimate, y = dtr_plot_label, color = model)
) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), orientation = "y", width = 0.3, linewidth = 0.8) +
  geom_point(size = 4.5) +
  facet_wrap(~ model, scales = "free_x", nrow = 1) +
    scale_color_manual(
    values = c(
      "Binary Model: P(Y=0)" = "#00274C",
      "Truncated Count Model: E[Y|Y>0]" = "#D9A100" # Deep presentation gold
    ),
    guide = "none"
  ) +
  labs(
    title = paste0("Conditional G-Computation DTR Estimates at Month ", t_eval),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 16) + 
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    strip.text = element_text(size = 16, face = "bold", margin = margin(b = 10)),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 13),
    panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.6),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.5), 
    panel.spacing = unit(2, "lines")
  )

##### Post-Split Trajectory Estimates #####
dtr_time_results <- purrr::map_dfr(c(2, 4), function(tt) {
  predict_conditional_dtr_means_delta(
    fits = fits,
    dtr_grid = dtr_grid %>% dplyr::select(A1, A2R, A2NR),
    Y0_vec = DataWideFormat$Y0,
    t_eval = tt
  ) %>%
    mutate(
      time = tt,
      p_zero_lower = p_zero - 1.96 * p_zero_se,
      p_zero_upper = p_zero + 1.96 * p_zero_se,
      m_plus_lower = m_plus - 1.96 * m_plus_se,
      m_plus_upper = m_plus + 1.96 * m_plus_se
    )
}) %>%
  left_join(dtr_lookup, by = "DTR")

##### Manuscript-Style Data Application Figures #####
DataLongAll <- DataLongAll %>%
  mutate(
    Time = factor(
      time,
      levels = c(0, 1, 2, 4),
      labels = paste0("Month ", c(0, 1, 2, 4))
    )
  )

fig_data_density <- ggplot(
  DataLongAll,
  aes(x = Y, linetype = Time, color = Time)
) +
  geom_density(linewidth = 0.85, adjust = 1.1) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "grey25", "grey50", "grey75")) +
  coord_cartesian(xlim = quantile(DataLongAll$Y, probs = c(0, 0.99), na.rm = TRUE)) +
  labs(
    x = "Number of Alcoholic Drinks Consumed Over the Last 30 Days",
    y = "Density",
    linetype = "Time",
    color = "Time"
  ) +
  biom_theme(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0))

outcome_proportions_long <- proportions_by_time %>%
  mutate(
    Time = factor(
      time,
      levels = c(0, 1, 2, 4),
      labels = paste0("Month ", c(0, 1, 2, 4))
    )
  ) %>%
  pivot_longer(
    cols = c(prop_zero, prop_positive),
    names_to = "Outcome",
    values_to = "Proportion"
  ) %>%
  mutate(
    Outcome = recode(
      Outcome,
      prop_zero = "Zero Outcomes",
      prop_positive = "Positive Outcomes"
    ),
    Outcome = factor(Outcome, levels = c("Zero Outcomes", "Positive Outcomes"))
  )

fig_observed_proportions <- ggplot(
  outcome_proportions_long,
  aes(x = Time, y = Proportion, group = Outcome)
) +
  geom_line(aes(linetype = Outcome), color = "black", linewidth = 0.85) +
  geom_point(aes(shape = Outcome), color = "black", fill = "white", size = 3.2, stroke = 0.9) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_shape_manual(values = c(16, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = NULL, y = "Observed Proportion", linetype = NULL, shape = NULL) +
  biom_theme(base_size = 12)

dtr_estimates_plot_data <- dtr_model_results %>%
  mutate(
    DTR_Number = factor(dtr_labels_number, levels = dtr_lookup$dtr_labels_number),
    Estimand = recode(
      as.character(model),
      "Binary Model: P(Y=0)" = "Probability of Zero Outcome",
      "Truncated Count Model: E[Y|Y>0]" = "Mean Count Among Positive Outcomes"
    ),
    Estimand = factor(
      Estimand,
      levels = c("Probability of Zero Outcome", "Mean Count Among Positive Outcomes")
    ),
    DTR_Forest_Label = factor(DTR_label, levels = unique(DTR_label))
  )

fig_dtr_estimates <- ggplot(
  dtr_estimates_plot_data,
  aes(x = Estimate, y = fct_rev(DTR_Forest_Label))
) +
  geom_errorbar(
    aes(xmin = Lower, xmax = Upper),
    orientation = "y",
    width = 0.22,
    linewidth = 0.7,
    color = "black"
  ) +
  geom_point(shape = 21, size = 2.8, stroke = 0.75, color = "black", fill = "black") +
  facet_wrap(~ Estimand, scales = "free_x", nrow = 1) +
  labs(x = "Estimate with 95% Confidence Interval", y = NULL) +
  biom_theme(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(size = 8.8)
  )

dtr_key <- dtr_estimates_plot_data %>%
  distinct(DTR_Number, DTR_label_int, DTR) %>%
  arrange(DTR_Number)

##### Save Outputs #####
save_data_app_figure(density_plot, "Density_Plot_Conditional_ppt", width = 13, height = 7)
save_data_app_figure(boxplot_outcome, "Boxplot_Outcome_Conditional_ppt", width = 13, height = 7)
save_data_app_figure(forest_plot, "DTR_ForestPlots_Conditional_ppt", width = 15, height = 8)
save_data_app_figure(forest_plot_ppt, "DTR_ForestPlots_Conditional_presentation_ppt", width = 15, height = 8)
save_data_app_figure(fig_data_density, "data_density_plot_conditional_biometrics", width = 6.5, height = 4.2)
save_data_app_figure(fig_observed_proportions, "figure1_observed_zero_positive_proportions", width = 6.5, height = 4.2)
save_data_app_figure(fig_dtr_estimates, "figure2_dtr_estimates_month4", width = 7.25, height = 5.4)

write.csv(proportions_by_time, file.path(output_dir_data, "Outcome_Proportions_By_Time_Conditional.csv"), row.names = FALSE)
write.csv(dtr_results, file.path(output_dir_data, paste0("DTR_Estimates_Conditional_Month_", t_eval, ".csv")), row.names = FALSE)
write.csv(dtr_model_results, file.path(output_dir_data, paste0("DTR_Estimates_By_Model_Conditional_Month_", t_eval, ".csv")), row.names = FALSE)
write.csv(dtr_time_results, file.path(output_dir_data, "DTR_Trajectories_Conditional_PostSplit.csv"), row.names = FALSE)
write.csv(dtr_key, file.path(output_dir_data, "dtr_number_key.csv"), row.names = FALSE)
readr::write_lines(
  c(
    "Data Density Plot. Overlaid density curves show the distribution of alcoholic drinks consumed over the last 30 days at months 0, 1, 2, and 4. Time points are distinguished by line type.",
    "Figure 1. Line graph of observed zero and positive outcome proportions at months 0, 1, 2, and 4.",
    "Figure 2. Two-panel forest plot of conditional g-computation DTR estimates at month 4. The first panel shows probability of zero outcome and the second shows mean count among positive outcomes; each DTR has a point estimate and 95% confidence interval."
  ),
  file.path(output_dir_data, "figure_alt_text.txt")
)

# print(dtr_model_results)
