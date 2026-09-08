#!/usr/bin/env Rscript

# Generate ROC and operating-point plots from the completed 32-model,
# patient-grouped repeated nested cross-validation analysis.
#
# Manuscript role
#   The manuscript uses four outputs from the August 14 run of this script:
#     * Figure3A_patient_grouped_repeated_nested_cv -> Extended Data Figure 5A
#     * Figure4A_patient_grouped_repeated_nested_cv -> Extended Data Figure 7A
#     * ExtendedDataFigure9A_patient_grouped_repeated_nested_cv
#       -> Extended Data Figure 9A
#     * ExtendedDataFigure9B_patient_grouped_repeated_nested_cv
#       -> Extended Data Figure 9B
#   The first two source filenames reflect an earlier panel plan. Script
#   6_19_Sync_ROC_Main_Operating_Extended_Data.R copies them to their final
#   Extended Data locations. The final main Figure 3A and Figure 4A ROC plots
#   come from 6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R.
#
# Inputs
#   The combined output directory created by
#   6_13_Assemble_All_Model_Grouped_CV_Results.R, including RUN_COMPLETE,
#   outer-held-out predictions, outer-fold metrics, and the 32-model performance
#   summary. The default is the 50-repeat result used in the manuscript.
#
# Analysis steps
#   1. Validate the completed combined run and its expected outer folds.
#   2. Reconstruct a mean ROC curve across outer repeats for each model.
#   3. Summarize fold-wise sensitivity and specificity for each model.
#   4. Draw the BM, blood, and fragmentomics ROC and operating-point plots.
#   5. Export the plotted source tables and a RUN_COMPLETE marker.
#
# Outputs
#   A new directory under
#   Output_figures_2025/patient_grouped_repeated_nested_cv/. It contains PNG and
#   PDF plots plus repeated_grouped_cv_mean_roc_source_data.csv and
#   grouped_cv_fold_operating_point_source_data.csv.
#
# R packages
#   dplyr, ggplot2, patchwork, readr, scales, and tidyr.
#
# Run from the repository root
#   Rscript 6_14_Generate_Grouped_CV_Manuscript_Replacement_Panels.R \
#     --input-run-id=<completed-6_13-run-id> \
#     --output-run-id=<new-figure-run-id>
#
# This script requires --name=value argument syntax. The output run ID must be
# new because the script refuses to overwrite an existing output directory.
#
# Historical figures and source files are read-only and are never overwritten.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

parse_named_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax: ", arg, call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[pieces[[1]]]] <- paste(pieces[-1], collapse = "=")
  }
  out
}

args <- parse_named_args(commandArgs(trailingOnly = TRUE))
arg_or_default <- function(name, default) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) default else value
}
input_run_id <- arg_or_default(
  "input-run-id", "2026-08-05_all_models_50repeats_combined_v3"
)
output_run_id <- arg_or_default(
  "output-run-id", "2026-08-14_no_grid_larger_titles_50repeats_v2"
)

input_dir <- file.path(
  "Output_tables_2025", "patient_grouped_repeated_nested_cv", input_run_id
)
output_dir <- file.path(
  "Output_figures_2025", "patient_grouped_repeated_nested_cv", output_run_id
)

if (!file.exists(file.path(input_dir, "RUN_COMPLETE"))) {
  stop("Definitive grouped-CV input is incomplete: ", input_dir, call. = FALSE)
}
if (dir.exists(output_dir)) {
  stop("Refusing to overwrite existing output directory: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

predictions <- read_csv(
  file.path(input_dir, "outer_heldout_predictions.csv"),
  show_col_types = FALSE
)
fold_metrics <- read_csv(
  file.path(input_dir, "outer_fold_metrics.csv"),
  show_col_types = FALSE
)
summary_tbl <- read_csv(
  file.path(input_dir, "publication_model_performance_and_legacy_comparison.csv"),
  show_col_types = FALSE
)

required_prediction_columns <- c(
  "model", "outer_repeat", "outer_fold", "Patient", "Sample_Code", "obs", "prob"
)
if (!all(required_prediction_columns %in% names(predictions))) {
  stop("Held-out prediction columns are incomplete.", call. = FALSE)
}
if (!all(c("model", "outer_repeat", "outer_fold", "sensitivity", "specificity") %in%
         names(fold_metrics))) {
  stop("Outer-fold metric columns are incomplete.", call. = FALSE)
}
if (any(!is.finite(predictions$prob))) {
  stop("Non-finite held-out probabilities detected.", call. = FALSE)
}

outer_repeats <- sort(unique(predictions$outer_repeat))
outer_folds <- sort(unique(predictions$outer_fold))
if (!setequal(outer_repeats, seq_len(max(outer_repeats))) ||
    length(outer_repeats) != max(outer_repeats)) {
  stop("Outer-repeat identifiers are incomplete or non-sequential.", call. = FALSE)
}
if (!setequal(outer_folds, seq_len(5L)) || length(outer_folds) != 5L) {
  stop("Expected exactly five outer folds numbered 1-5.", call. = FALSE)
}
n_outer_repeats <- length(outer_repeats)
expected_fold_rows <- n_outer_repeats * length(outer_folds)
fold_count_qc <- fold_metrics |>
  count(model, name = "n_rows")
if (any(fold_count_qc$n_rows != expected_fold_rows)) {
  stop("Each model must contain one metric row per outer assessment fold.",
       call. = FALSE)
}

model_metadata <- tribble(
  ~model, ~cohort, ~label, ~short_label, ~color,
  "BM_Sites", "BM", "Sites Model", "Sites", "#2C7FB8",
  "BM_cVAF", "BM", "cVAF Model", "cVAF z-score", "#41AB5D",
  "BM_Raw_cVAF", "BM", "Raw cVAF", "Raw cVAF", "#F28E2B",
  "BM_Combined_Mutation_Zscores", "BM", "Combined Model", "Combined", "#6A3D9A",
  "BM_All_Mutation_Features", "BM", "All mutation features", "All mutation", "#1B9E77",
  "BM_Mutation_Fragmentomics_Full", "BM", "Mutation + fragmentomics (full)", "Mut. + frag. full", "#D95F02",
  "BM_Mutation_Fragmentomics_Min", "BM", "Mutation + fragmentomics (minimal)", "Mut. + frag. min", "#7570B3",
  "Blood_Sites", "Blood", "Sites Model", "Sites", "#2C7FB8",
  "Blood_cVAF", "Blood", "cVAF Model", "cVAF z-score", "#41AB5D",
  "Blood_Raw_cVAF", "Blood", "Raw cVAF", "Raw cVAF", "#F28E2B",
  "Blood_Combined_Mutation_Zscores", "Blood", "Combined Model", "Combined", "#6A3D9A",
  "Blood_All_Mutation_Features", "Blood", "All mutation features", "All mutation", "#1B9E77",
  "Blood_Mutation_Fragmentomics_Full", "Blood", "Mutation + fragmentomics (full)", "Mut. + frag. full", "#D95F02",
  "Blood_Mutation_Fragmentomics_Min", "Blood", "Mutation + fragmentomics (minimal)", "Mut. + frag. min", "#7570B3",
  "Fragmentomics_FullCohort_Full", "FullFrag", "Fragmentomics (full)", "Full", "#000000",
  "Fragmentomics_FullCohort_Min", "FullFrag", "Fragmentomics (minimal)", "FS + mean coverage", "#E66101",
  "Fragmentomics_FullCohort_FS", "FullFrag", "Fragment size score", "Fragment size", "#1F78B4",
  "Fragmentomics_FullCohort_Mean_Coverage", "FullFrag", "Mean coverage", "Mean coverage", "#009E73",
  "Fragmentomics_FullCohort_Proportion_Short", "FullFrag", "Proportion short fragments", "Proportion short", "#E6AB02",
  "Fragmentomics_FullCohort_Tumor_Fraction", "FullFrag", "Tumor fraction", "Tumor fraction", "#D95F02"
)

if (anyDuplicated(model_metadata$model)) {
  stop("Model metadata names must be unique.", call. = FALSE)
}

auc_direct <- function(obs, prob) {
  y <- obs == "pos"
  n_pos <- sum(y)
  n_neg <- sum(!y)
  ranks <- rank(prob, ties.method = "average")
  (sum(ranks[y]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

roc_points <- function(obs, prob) {
  ord <- order(prob, decreasing = TRUE)
  y <- obs[ord] == "pos"
  tp <- c(0, cumsum(y))
  fp <- c(0, cumsum(!y))
  tibble(fpr = fp / sum(!y), tpr = tp / sum(y)) |>
    group_by(fpr) |>
    summarise(tpr = max(tpr), .groups = "drop") |>
    arrange(fpr)
}

fpr_grid <- seq(0, 1, length.out = 201)
roc_by_repeat <- predictions |>
  group_by(model, outer_repeat) |>
  group_modify(~ {
    curve <- roc_points(.x$obs, .x$prob)
    tibble(
      fpr = fpr_grid,
      tpr = approx(curve$fpr, curve$tpr, xout = fpr_grid,
                   method = "constant", f = 0, rule = 2)$y,
      auc = auc_direct(.x$obs, .x$prob)
    )
  }) |>
  ungroup()

mean_roc <- roc_by_repeat |>
  group_by(model, fpr) |>
  summarise(
    mean_tpr = mean(tpr),
    tpr_q025 = quantile(tpr, 0.025),
    tpr_q975 = quantile(tpr, 0.975),
    .groups = "drop"
  ) |>
  left_join(model_metadata, by = "model") |>
  left_join(
    summary_tbl |>
      select(model, repeat_pooled_auc_mean, auc_cluster_boot_q025,
             auc_cluster_boot_q975),
    by = "model"
  ) |>
  mutate(
    legend_label = sprintf(
      "%s: AUC %.2f (95%% CI %.2f-%.2f)",
      short_label, repeat_pooled_auc_mean,
      auc_cluster_boot_q025, auc_cluster_boot_q975
    )
  )

fold_operating_summary <- fold_metrics |>
  group_by(model) |>
  summarise(
    sensitivity_mean = mean(sensitivity, na.rm = TRUE),
    sensitivity_sd = sd(sensitivity, na.rm = TRUE),
    specificity_mean = mean(specificity, na.rm = TRUE),
    specificity_sd = sd(specificity, na.rm = TRUE),
    accuracy_mean = mean(accuracy, na.rm = TRUE),
    n_outer_folds = n(),
    .groups = "drop"
  ) |>
  left_join(model_metadata, by = "model") |>
  left_join(
    summary_tbl |>
      select(model, repeat_pooled_auc_mean),
    by = "model"
  ) |>
  mutate(
    legend_label = sprintf(
      "%s (%d%% sens, %d%% spec)", short_label,
      round(100 * sensitivity_mean), round(100 * specificity_mean)
    )
  )

write_csv(
  mean_roc,
  file.path(output_dir, "repeated_grouped_cv_mean_roc_source_data.csv")
)
write_csv(
  fold_operating_summary,
  file.path(output_dir, "grouped_cv_fold_operating_point_source_data.csv")
)

# Match the established manuscript visual language used by the original
# 3_1-generated ROC and operating-point panels. Only presentation changes here;
# all values continue to come from the definitive grouped-CV analysis above.
base_theme <- theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 8.2, hjust = 0.5,
                              lineheight = 0.95),
    axis.title = element_text(size = 7.4),
    axis.text = element_text(size = 6.5, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey55", linewidth = 0.35),
    legend.title = element_text(size = 5.8, face = "plain"),
    legend.text = element_text(size = 5.15, lineheight = 0.85),
    legend.key = element_blank(),
    legend.key.height = unit(0.17, "cm"),
    legend.key.width = unit(0.30, "cm"),
    legend.spacing.y = unit(0.015, "cm"),
    legend.margin = margin(1, 2, 1, 2),
    legend.background = element_rect(
      fill = alpha("white", 0.80), color = NA
    ),
    plot.margin = margin(4, 4, 4, 4)
  )

make_operating_plot <- function(cohort_name, title_text, title_size = 8.2) {
  dat <- fold_operating_summary |>
    filter(cohort == cohort_name) |>
    arrange(desc(repeat_pooled_auc_mean)) |>
    mutate(legend_label = factor(legend_label, levels = legend_label))
  colors <- setNames(dat$color, dat$legend_label)
  ggplot(dat, aes(sensitivity_mean, specificity_mean, color = legend_label)) +
    # Retain only the 50% quadrant guides; all regular panel gridlines are
    # removed by base_theme.
    geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.35,
               color = "grey65") +
    geom_hline(yintercept = 0.5, linetype = "dotted", linewidth = 0.35,
               color = "grey65") +
    geom_errorbar(
      aes(ymin = pmax(0, specificity_mean - specificity_sd),
          ymax = pmin(1, specificity_mean + specificity_sd)),
      width = 0, linewidth = 0.35
    ) +
    geom_errorbar(
      aes(xmin = pmax(0, sensitivity_mean - sensitivity_sd),
          xmax = pmin(1, sensitivity_mean + sensitivity_sd)),
      width = 0, linewidth = 0.35, orientation = "y"
    ) +
    geom_point(size = 2) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = label_number(accuracy = 0.01)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = label_number(accuracy = 0.01)) +
    labs(
      x = "Mean sensitivity (CV)",
      y = "Mean specificity (CV)",
      color = "Model",
      title = title_text
    ) +
    base_theme +
    coord_fixed(ratio = 1) +
    theme(
      plot.title = element_text(
        face = "bold", size = title_size, hjust = 0.5, lineheight = 0.95
      ),
      legend.position = c(0.015, 0.015),
      legend.justification = c(0, 0)
    ) +
    guides(color = guide_legend(
      ncol = 1, byrow = TRUE,
      override.aes = list(size = 1.7, linewidth = 0.5)
    ))
}

make_roc_plot <- function(cohort_name, title_text, title_size = 8.2) {
  dat <- mean_roc |>
    filter(cohort == cohort_name) |>
    arrange(desc(repeat_pooled_auc_mean), fpr) |>
    mutate(legend_label = factor(legend_label, levels = unique(legend_label)))
  colors <- setNames(unique(dat[, c("legend_label", "color")])$color,
                     unique(dat[, c("legend_label", "color")])$legend_label)
  ggplot(dat, aes(fpr, mean_tpr, color = legend_label)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                linewidth = 0.35, color = "grey60") +
    geom_line(linewidth = 0.65) +
    scale_color_manual(values = colors, drop = FALSE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    labs(
      x = "False-positive rate (1 - specificity)",
      y = "True-positive rate (sensitivity)",
      color = NULL,
      title = title_text
    ) +
    base_theme +
    theme(
      plot.title = element_text(
        face = "bold", size = title_size, hjust = 0.5, lineheight = 0.95
      ),
      legend.position = c(0.985, 0.015),
      legend.justification = c(1, 0)
    ) +
    guides(color = guide_legend(
      ncol = 1, byrow = TRUE,
      override.aes = list(linewidth = 0.7)
    ))
}

save_panel <- function(plot, basename, width, height) {
  ggsave(file.path(output_dir, paste0(basename, ".png")), plot,
         width = width, height = height, units = "in", dpi = 600,
         bg = "white")
  ggsave(file.path(output_dir, paste0(basename, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf,
         bg = "white")
}

f3a <- make_operating_plot(
  "BM", "Fold-Wise Sensitivity & Specificity for\nBM-Derived Features (mean ± SD)",
  title_size = 9.5
)
f4a <- make_operating_plot(
  "Blood", "Fold-Wise Sensitivity & Specificity for\ncfDNA-Derived Features (mean ± SD)",
  title_size = 9.5
)
ed5a <- make_roc_plot(
  "BM", sprintf(
    "Mean Outer-Held-Out ROC Curve on\nBM-Derived Features (%d × 5 grouped CV)",
    n_outer_repeats
  ),
  title_size = 9.5
)
ed7a <- make_roc_plot(
  "Blood", sprintf(
    "Mean Outer-Held-Out ROC Curve on\ncfDNA-Derived Features (%d × 5 grouped CV)",
    n_outer_repeats
  ),
  title_size = 9.5
)
ed9a <- make_roc_plot(
  "FullFrag", sprintf(
    "Mean Outer-Held-Out ROC Curve for\nFragmentomics-Only Models (%d × 5 grouped CV)",
    n_outer_repeats
  ),
  title_size = 11.5
)
ed9b <- make_operating_plot(
  "FullFrag", "Fold-Wise Sensitivity & Specificity for\nFragmentomics-Only Models (mean ± SD)",
  title_size = 11.5
)

save_panel(f3a, "Figure3A_patient_grouped_repeated_nested_cv", 3.6, 3.6)
save_panel(f4a, "Figure4A_patient_grouped_repeated_nested_cv", 3.6, 3.6)
save_panel(ed5a, "ExtendedDataFigure5A_patient_grouped_repeated_nested_cv", 4.3, 4.3)
save_panel(ed7a, "ExtendedDataFigure7A_patient_grouped_repeated_nested_cv", 4.3, 4.3)
save_panel(ed9a, "ExtendedDataFigure9A_patient_grouped_repeated_nested_cv", 4.5, 4.3)
save_panel(ed9b, "ExtendedDataFigure9B_patient_grouped_repeated_nested_cv", 4.5, 4.3)
save_panel(
  ed9a + ed9b + plot_layout(widths = c(1, 1)),
  "ExtendedDataFigure9AB_patient_grouped_repeated_nested_cv", 9.2, 4.3
)

expected_files <- c(
  "Figure3A_patient_grouped_repeated_nested_cv.png",
  "Figure4A_patient_grouped_repeated_nested_cv.png",
  "ExtendedDataFigure5A_patient_grouped_repeated_nested_cv.png",
  "ExtendedDataFigure7A_patient_grouped_repeated_nested_cv.png",
  "ExtendedDataFigure9A_patient_grouped_repeated_nested_cv.png",
  "ExtendedDataFigure9B_patient_grouped_repeated_nested_cv.png",
  "ExtendedDataFigure9AB_patient_grouped_repeated_nested_cv.png"
)
if (!all(file.exists(file.path(output_dir, expected_files)))) {
  stop("One or more replacement panels were not created.", call. = FALSE)
}

writeLines(
  c(
    "Definitive grouped-CV manuscript replacement panels",
    paste0("Input: ", input_run_id),
    sprintf(
      "Outer resampling: %d repeats of patient-grouped five-fold CV",
      n_outer_repeats
    ),
    "Inner resampling: 5 repeats of patient-grouped five-fold CV",
    "Thresholds: derived only from grouped inner out-of-fold predictions",
    "Display order: decreasing repeat-pooled AUC in every multi-model panel",
    paste(
      "Visual style: no regular panel gridlines; operating-point panels retain",
      paste0(
        "dotted 50% horizontal and vertical quadrant guides; panel-specific ",
        "title sizing; values unchanged"
      )
    ),
    "No historical figure or result file was overwritten."
  ),
  file.path(output_dir, "README.txt")
)
writeLines("PASS", file.path(output_dir, "RUN_COMPLETE"))

message("Replacement panels written to: ", output_dir)
