#!/usr/bin/env Rscript

# Generate compact all-model ROC panels from the patient-grouped repeated
# nested cross-validation analysis. Curves connect empirical mean ROC
# coordinates; no statistical smoothing is applied.
#
# Manuscript role
#   This script generates the 50-repeat ROC plots used in Figure 3A (BM-informed
#   models) and Figure 4A (baseline-plasma-informed models). The PNG/PDF files
#   written by this script are the panel files used in those figures.
#
# Inputs
#   The combined 32-model output from
#   6_13_Assemble_All_Model_Grouped_CV_Results.R, specifically RUN_COMPLETE,
#   outer_heldout_predictions.csv, and
#   publication_model_performance_and_legacy_comparison.csv. The default is the
#   50-repeat combined result used in the manuscript.
#
# Analysis steps
#   1. Validate the requested BM and blood model sets.
#   2. Calculate an empirical ROC curve within every outer repeat.
#   3. Interpolate each curve on a common false-positive-rate grid and average
#      sensitivity across repeats.
#   4. Join the repeat-pooled AUCs and patient-clustered confidence intervals.
#   5. Draw and export the BM and blood all-model ROC panels and source tables.
#
# Outputs
#   A new versioned directory under
#   Output_figures_2025/patient_grouped_repeated_nested_cv/ containing Figure 3A
#   and Figure 4A as PNG/PDF files, the mean-ROC source table, panel statistics,
#   README.txt, and RUN_COMPLETE.
#
# R packages
#   dplyr, ggplot2, readr, scales, and tibble.
#
# Run from the analysis project directory that contains Output_tables_2025/
#   Rscript Scripts_2025/Final_Scripts/6_17_Generate_Compact_All_Model_Grouped_CV_ROC_Panels.R \
#     --input-run-id <completed-6_13-run-id> \
#     --output-run-id <new-figure-run-id>
#
# The output run ID must be new because the script refuses to overwrite an
# existing output directory.
#
# Historical outputs are not overwritten.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) && length(args) > hit[1]) args[hit[1] + 1] else default
}

input_run_id <- get_arg(
  "--input-run-id", "2026-08-05_all_models_50repeats_combined_v3"
)
output_run_id <- get_arg(
  "--output-run-id", "2026-08-05_compact_all_model_roc_50repeats_v3"
)
valid_run_id <- "^[A-Za-z0-9._-]+$"
if (!grepl(valid_run_id, input_run_id) || !grepl(valid_run_id, output_run_id)) {
  stop("Run IDs may contain only letters, numbers, period, underscore, and hyphen.",
       call. = FALSE)
}

input_dir <- file.path(
  "Output_tables_2025", "patient_grouped_repeated_nested_cv",
  input_run_id
)
output_dir <- file.path(
  "Output_figures_2025", "patient_grouped_repeated_nested_cv",
  output_run_id
)

if (!file.exists(file.path(input_dir, "RUN_COMPLETE"))) {
  stop("Completed grouped-CV input is missing RUN_COMPLETE: ", input_dir,
       call. = FALSE)
}
if (dir.exists(output_dir)) {
  stop("Refusing to overwrite existing output directory: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

predictions <- read_csv(
  file.path(input_dir, "outer_heldout_predictions.csv"), show_col_types = FALSE
)
performance <- read_csv(
  file.path(input_dir, "publication_model_performance_and_legacy_comparison.csv"),
  show_col_types = FALSE
)

metadata <- tribble(
  ~model, ~cohort, ~display_name, ~colour, ~selected, ~order,
  "BM_Sites", "BM", "BM Sites Z-score", "#D55E00", TRUE, 1,
  "BM_Raw_cVAF", "BM", "BM cVAF", "#CC79A7", FALSE, 2,
  "BM_cVAF", "BM", "BM cVAF z-score", "#0072B2", FALSE, 3,
  "BM_Combined_Mutation_Zscores", "BM", "BM Sites + cVAF z-score", "#E69F00", FALSE, 4,
  "BM_All_Mutation_Features", "BM", "All Mut. Features", "#000000", FALSE, 5,
  "BM_Mutation_Fragmentomics_Full", "BM", "BM Mut. + Fragmentomics", "#56B4E9", FALSE, 6,
  "BM_Mutation_Fragmentomics_Min", "BM", "BM Mut. + Fragments (min)", "#009E73", FALSE, 7,
  "Blood_Mutation_Fragmentomics_Full", "Blood", "Blood Mut. + Fragmentomics", "#56B4E9", TRUE, 1,
  "Blood_All_Mutation_Features", "Blood", "All Mut. Features", "#000000", FALSE, 2,
  "Blood_Mutation_Fragmentomics_Min", "Blood", "Blood Mut. + Fragments (min)", "#009E73", FALSE, 3,
  "Blood_Sites", "Blood", "Blood Sites Z-score", "#D55E00", FALSE, 4,
  "Blood_Raw_cVAF", "Blood", "Blood cVAF", "#CC79A7", FALSE, 5,
  "Blood_Combined_Mutation_Zscores", "Blood", "Blood Sites + cVAF Z-score", "#E69F00", FALSE, 6,
  "Blood_cVAF", "Blood", "Blood cVAF Z-score", "#0072B2", FALSE, 7
)

required_columns <- c("model", "outer_repeat", "outer_fold", "obs", "prob")
if (!all(required_columns %in% names(predictions))) {
  stop("Outer held-out predictions are missing required columns.", call. = FALSE)
}
if (!all(metadata$model %in% predictions$model) ||
    !all(metadata$model %in% performance$model)) {
  stop("One or more requested models are absent from the final run.", call. = FALSE)
}

empirical_roc <- function(obs, prob) {
  ordering <- order(prob, decreasing = TRUE)
  positive <- obs[ordering] == "pos"
  if (sum(positive) == 0L || sum(!positive) == 0L) {
    stop("Both classes are required for an ROC curve.", call. = FALSE)
  }
  tibble(
    fpr = c(0, cumsum(!positive)) / sum(!positive),
    tpr = c(0, cumsum(positive)) / sum(positive)
  ) |>
    group_by(fpr) |>
    summarise(tpr = max(tpr), .groups = "drop") |>
    arrange(fpr)
}

fpr_grid <- seq(0, 1, length.out = 401)
repeat_rocs <- predictions |>
  filter(model %in% metadata$model) |>
  group_by(model, outer_repeat) |>
  group_modify(~ {
    curve <- empirical_roc(.x$obs, .x$prob)
    tibble(
      fpr = fpr_grid,
      tpr = approx(
        curve$fpr, curve$tpr, xout = fpr_grid,
        method = "constant", f = 0, rule = 2
      )$y
    )
  }) |>
  ungroup()

mean_rocs <- repeat_rocs |>
  summarise(mean_tpr = mean(tpr), .by = c(model, fpr)) |>
  left_join(metadata, by = "model") |>
  left_join(
    performance |>
      select(model, repeat_pooled_auc_mean, auc_cluster_boot_q025,
             auc_cluster_boot_q975, sensitivity_estimate, specificity_estimate,
             outer_repeats, outer_folds, inner_repeats, inner_folds),
    by = "model"
  ) |>
  mutate(legend_label = sprintf("%s, AUC = %.2f", display_name,
                                repeat_pooled_auc_mean))

design_qc <- mean_rocs |>
  distinct(model, outer_repeats, outer_folds, inner_repeats, inner_folds)
if (length(unique(design_qc$outer_repeats)) != 1L ||
    any(design_qc$outer_folds != 5L) ||
    any(design_qc$inner_repeats != 5L) || any(design_qc$inner_folds != 5L)) {
  stop("Final resampling parameters do not match the expected design.", call. = FALSE)
}

panel_statistics <- mean_rocs |>
  distinct(model, cohort, display_name, selected, order,
           repeat_pooled_auc_mean, auc_cluster_boot_q025,
           auc_cluster_boot_q975, sensitivity_estimate, specificity_estimate) |>
  arrange(cohort, order)
write_csv(mean_rocs, file.path(output_dir, "compact_all_model_mean_roc_source_data.csv"))
write_csv(panel_statistics, file.path(output_dir, "compact_all_model_panel_statistics.csv"))

build_panel <- function(cohort_name, feature_title) {
  dat <- mean_rocs |>
    filter(cohort == cohort_name) |>
    arrange(desc(repeat_pooled_auc_mean)) |>
    mutate(legend_label = factor(legend_label, levels = unique(legend_label)))
  selected_curve <- dat |> filter(selected)
  colour_table <- dat |> distinct(legend_label, colour)
  colours <- setNames(colour_table$colour, colour_table$legend_label)

  ggplot(dat, aes(fpr, mean_tpr, colour = legend_label, group = model)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                linewidth = 0.35, colour = "grey70") +
    geom_line(data = dat |> filter(!selected), linewidth = 1) +
    geom_line(data = selected_curve, linewidth = 1) +
    scale_colour_manual(values = colours, drop = FALSE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25))) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25))) +
    labs(
      title = paste0("Mean Outer-Fold ROC Curves on\n", feature_title,
                     "\n(patient-grouped nested CV)"),
      x = "False-positive rate (1 - specificity)",
      y = "True-positive rate (sensitivity)",
      colour = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = c(0.63, 0.15),
      legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                       colour = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(colour = guide_legend(
      ncol = 1, byrow = TRUE,
      override.aes = list(linewidth = 1)
    ))
}

figure_3a <- build_panel("BM", "BM-Derived Features")
figure_4a <- build_panel("Blood", "cfDNA-Derived Features")

save_panel <- function(plot, stem) {
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
         width = 6, height = 6, units = "in", dpi = 500, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
         width = 6, height = 6, units = "in", device = cairo_pdf,
         bg = "white")
}

save_panel(figure_3a, "Figure3A_BM_all_models_compact_updated_ROC")
save_panel(figure_4a, "Figure4A_Blood_all_models_compact_updated_ROC")

expected <- c(
  "Figure3A_BM_all_models_compact_updated_ROC.png",
  "Figure3A_BM_all_models_compact_updated_ROC.pdf",
  "Figure4A_Blood_all_models_compact_updated_ROC.png",
  "Figure4A_Blood_all_models_compact_updated_ROC.pdf",
  "compact_all_model_mean_roc_source_data.csv",
  "compact_all_model_panel_statistics.csv"
)
if (!all(file.exists(file.path(output_dir, expected)))) {
  stop("One or more compact ROC artifacts were not generated.", call. = FALSE)
}

writeLines(c(
  "Compact all-model ROC panels completed successfully.",
  paste0("Input: ", input_run_id),
  paste0(
    "Outer resampling: ", unique(design_qc$outer_repeats),
    " repeats of patient-grouped five-fold CV"
  ),
  "Inner resampling: 5 repeats of patient-grouped five-fold CV",
  paste0(
    "Curves: mean empirical ROC across the ",
    unique(design_qc$outer_repeats), " outer-held-out prediction sets"
  ),
  "Rendering: empirical mean coordinates connected by straight line segments; no smoothing",
  "No operating-point marker is displayed in these all-model panels",
  "Historical results and figures were not overwritten."
), file.path(output_dir, "README.txt"))
writeLines("PASS", file.path(output_dir, "RUN_COMPLETE"))

message("Compact all-model ROC panels written to: ", output_dir)
