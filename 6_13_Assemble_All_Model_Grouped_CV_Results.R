#!/usr/bin/env Rscript

# =============================================================================
# 6_13_Assemble_All_Model_Grouped_CV_Results.R
#
# Goal
#   Assemble independently executed blocks from the complete 32-model grouped
#   nested-CV sensitivity analysis into one result set used for the manuscript.
#
# Why this is separate
#   The full model library was run in smaller BM, blood, and fragmentomics
#   blocks for runtime management. Combining them requires more than appending
#   tables: every model must be present once, source QC must pass, probabilities
#   must be finite, and patient-fold assignments must remain consistent.
#
# Manuscript role
#   This script combines the three 50-repeat source runs used for Figure 3A,
#   Figure 4A, Extended Data Figure 5A, Extended Data Figure 7A, Extended Data
#   Figure 9A-B, and Supplementary Table 4. It writes the 32-model result set
#   read by the downstream grouped-CV plotting and table scripts.
#
# Inputs
#   Three completed output directories from
#   6_12_Patient_Grouped_Repeated_Nested_CV.R. Each directory must contain its
#   RUN_COMPLETE marker, QC table, summaries, held-out predictions, patient-fold
#   assignments, and captured model-fitting warnings. The preserved legacy
#   validation RDS files are read only to report the grouped-versus-row-level
#   AUC differences.
#
# Outputs
#   A new combined directory under
#   Output_tables_2025/patient_grouped_repeated_nested_cv/. It contains the
#   merged detail tables, the 32-model summary, source-run manifest, combined QC
#   table, warning summary, and RUN_COMPLETE marker.
#
# R packages
#   dplyr, purrr, readr, tibble, and tidyr.
#
# Analysis steps
#   1. Resolve the explicitly named completed source runs.
#   2. Reject missing files or incomplete source runs.
#   3. Merge detailed folds, predictions, tuning results, and bootstraps.
#   4. Validate model completeness, uniqueness, and fold integrity.
#   5. Compare grouped-CV AUCs with preserved legacy row-level estimates.
#   6. Write a combined manifest, QC table, and RUN_COMPLETE marker.
#
# Scientific boundary
#   This script does not fit models, select features, change thresholds,
#   overwrite source runs, or rescore test samples. Source and output run IDs
#   are deliberately explicit below and must be reviewed for a new release.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) && length(args) > hit[1]) args[hit[1] + 1] else default
}

# ----------------------------------------------------------------------------
# 1. Locked source-run selection and expected model inventory
# ----------------------------------------------------------------------------
project_root <- normalizePath(getwd(), mustWork = TRUE)
results_root <- file.path(
  project_root,
  "Output_tables_2025",
  "patient_grouped_repeated_nested_cv"
)

default_source_run_ids <- c(
  "2026-08-05_all_models_50repeats_bm_v4",
  "2026-08-05_all_models_50repeats_blood_v4",
  "2026-08-05_all_models_50repeats_fullfrag_v4"
)
source_run_ids <- trimws(strsplit(
  get_arg("--source-runs", paste(default_source_run_ids, collapse = ",")),
  ",", fixed = TRUE
)[[1]])
output_run_id <- get_arg(
  "--output-run-id",
  "2026-08-05_all_models_50repeats_combined_v3"
)
if (length(source_run_ids) != 3L || any(!nzchar(source_run_ids))) {
  stop("--source-runs must contain exactly three comma-separated run IDs.",
       call. = FALSE)
}
valid_run_id <- "^[A-Za-z0-9._-]+$"
if (any(!grepl(valid_run_id, source_run_ids)) ||
    !grepl(valid_run_id, output_run_id)) {
  stop("Run IDs may contain only letters, numbers, period, underscore, and hyphen.",
       call. = FALSE)
}
source_dirs <- file.path(results_root, source_run_ids)
output_dir <- file.path(results_root, output_run_id)

expected_models <- c(
  "BM_Sites", "BM_cVAF", "BM_Raw_cVAF", "BM_All_Mutation_Features",
  "BM_Combined_Mutation_Zscores", "BM_Mutation_Fragmentomics_Full",
  "BM_Mutation_Fragmentomics_Min",
  "Blood_Sites", "Blood_cVAF", "Blood_Raw_cVAF",
  "Blood_All_Mutation_Features", "Blood_Combined_Mutation_Zscores",
  "Blood_Mutation_Fragmentomics_Full", "Blood_Mutation_Fragmentomics_Min",
  paste0("Fragmentomics_FullCohort_", c(
    "Full", "Min", "FS", "Mean_Coverage", "Proportion_Short", "Tumor_Fraction"
  )),
  paste0("Fragmentomics_BMRestricted_", c(
    "Full", "Min", "FS", "Mean_Coverage", "Proportion_Short", "Tumor_Fraction"
  )),
  paste0("Fragmentomics_BloodRestricted_", c(
    "Full", "Min", "FS", "Mean_Coverage", "Proportion_Short", "Tumor_Fraction"
  ))
)

required_files <- c(
  "RUN_COMPLETE", "quality_control_checks.csv", "primary_performance_summary.csv",
  "outer_heldout_predictions.csv", "outer_patient_fold_assignments.csv",
  "inner_patient_fold_assignments.csv", "model_fitting_warnings.csv"
)

# ----------------------------------------------------------------------------
# 2. Reject incomplete blocks before creating a combined result
# ----------------------------------------------------------------------------
missing_inputs <- map2(source_dirs, source_run_ids, function(path, run_id) {
  absent <- required_files[!file.exists(file.path(path, required_files))]
  if (length(absent) == 0) return(NULL)
  tibble(source_run_id = run_id, missing_file = absent)
}) |> bind_rows()

if (nrow(missing_inputs) > 0) {
  stop(
    "Cannot assemble incomplete source runs:\n",
    paste0(missing_inputs$source_run_id, ": ", missing_inputs$missing_file, collapse = "\n"),
    call. = FALSE
  )
}

if (dir.exists(output_dir)) {
  stop("Refusing to overwrite existing combined result directory: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 3. Merge row-level records without modifying source runs
# ----------------------------------------------------------------------------
read_blocks <- function(filename) {
  map2_dfr(source_dirs, source_run_ids, function(path, run_id) {
    read_csv(file.path(path, filename), show_col_types = FALSE) |>
      mutate(source_run_id = run_id, .before = 1)
  })
}

merge_files <- c(
  "training_frame_manifest.csv",
  "outer_patient_fold_assignments.csv",
  "outer_fold_balance.csv",
  "outer_heldout_predictions.csv",
  "outer_fold_metrics.csv",
  "outer_repeat_pooled_metrics.csv",
  "inner_tuning_and_thresholds.csv",
  "inner_patient_fold_assignments.csv",
  "model_fitting_warnings.csv",
  "patient_cluster_bootstrap_metrics.csv",
  "primary_performance_summary.csv",
  "current_rows_excluded_to_match_frozen_training.csv",
  "input_and_legacy_artifact_manifest.csv"
)

walk(merge_files, function(filename) {
  write_csv(read_blocks(filename), file.path(output_dir, filename), na = "")
})

summary_tbl <- read_blocks("primary_performance_summary.csv")
fold_metrics_tbl <- read_blocks("outer_fold_metrics.csv")
pred_tbl <- read_blocks("outer_heldout_predictions.csv")
outer_assignments <- read_blocks("outer_patient_fold_assignments.csv")
inner_assignments <- read_blocks("inner_patient_fold_assignments.csv")
warnings_tbl <- read_blocks("model_fitting_warnings.csv")
source_qc <- read_blocks("quality_control_checks.csv")
observed_outer_repeats <- sort(unique(summary_tbl$outer_repeats))
observed_inner_repeats <- sort(unique(summary_tbl$inner_repeats))

# ----------------------------------------------------------------------------
# 4. Validate completeness, uniqueness, fold integrity, and finite predictions
# ----------------------------------------------------------------------------
observed_models <- sort(unique(summary_tbl$model))

prediction_key <- intersect(
  c("model", "outer_repeat", "Sample_Code", "sample_id", "sample_key", "Patient", "Patient_ID"),
  names(pred_tbl)
)
if (!all(c("model", "outer_repeat") %in% prediction_key)) {
  stop("Held-out prediction table lacks model/repeat identifiers.", call. = FALSE)
}

probability_column <- intersect(c("probability", "predicted_probability", "prob"), names(pred_tbl))
if (length(probability_column) != 1) {
  stop("Could not identify exactly one probability column.", call. = FALSE)
}

duplicate_predictions <- pred_tbl |>
  count(across(all_of(prediction_key)), name = "n") |>
  filter(n != 1)

fold_consistency <- outer_assignments |>
  distinct(model, outer_repeat, Patient, fold) |>
  mutate(
    resampling_cohort = case_when(
      grepl("^BM_|^Fragmentomics_BMRestricted", model) ~ "BM",
      grepl("^Blood_|^Fragmentomics_BloodRestricted", model) ~ "Blood",
      grepl("^Fragmentomics_FullCohort", model) ~ "FullCohort",
      TRUE ~ NA_character_
    )
  ) |>
  group_by(resampling_cohort, outer_repeat, Patient) |>
  summarise(n_folds = n_distinct(fold), .groups = "drop") |>
  filter(is.na(resampling_cohort) | n_folds != 1)

qc <- tribble(
  ~check, ~status, ~detail,
  "all_source_runs_complete", "PASS", paste(source_run_ids, collapse = "; "),
  "all_source_qc_checks_pass", if_else(all(source_qc$status == "PASS"), "PASS", "FAIL"), paste0(nrow(source_qc), " checks inspected"),
  "exact_expected_model_library", if_else(identical(observed_models, sort(expected_models)), "PASS", "FAIL"), paste0(length(observed_models), " models observed; 32 expected"),
  "one_primary_summary_row_per_model", if_else(nrow(summary_tbl) == 32L && !anyDuplicated(summary_tbl$model), "PASS", "FAIL"), paste0(nrow(summary_tbl), " rows"),
  "consistent_outer_repeats_for_every_model", if_else(length(observed_outer_repeats) == 1L, "PASS", "FAIL"), paste0("Observed outer repeats: ", paste(observed_outer_repeats, collapse = ", ")),
  "five_inner_repeats_for_every_model", if_else(identical(observed_inner_repeats, 5), "PASS", "FAIL"), paste0("Observed inner repeats: ", paste(observed_inner_repeats, collapse = ", ")),
  "finite_primary_metrics", if_else(all(is.finite(summary_tbl$auc_estimate)), "PASS", "FAIL"), "AUC checked for every model",
  "finite_heldout_probabilities", if_else(all(is.finite(pred_tbl[[probability_column]])), "PASS", "FAIL"), paste0(nrow(pred_tbl), " held-out predictions inspected"),
  "unique_model_repeat_sample_predictions", if_else(nrow(duplicate_predictions) == 0, "PASS", "FAIL"), paste0(nrow(duplicate_predictions), " duplicate/non-unique keys"),
  "identical_patient_folds_within_each_cohort", if_else(nrow(fold_consistency) == 0, "PASS", "FAIL"), paste0(nrow(fold_consistency), " inconsistent patient/repeat assignments"),
  "inner_assignments_present", if_else(nrow(inner_assignments) > 0, "PASS", "FAIL"), paste0(nrow(inner_assignments), " inner patient-fold assignments"),
  "warnings_retained_for_audit", "PASS", paste0(nrow(warnings_tbl), " captured warning events")
)
write_csv(qc, file.path(output_dir, "quality_control_checks.csv"))

if (any(qc$status != "PASS")) {
  write_csv(duplicate_predictions, file.path(output_dir, "failed_duplicate_prediction_keys.csv"))
  write_csv(fold_consistency, file.path(output_dir, "failed_fold_consistency.csv"))
  stop("Combined-result quality control failed; inspect ", output_dir, call. = FALSE)
}

warning_summary <- warnings_tbl |>
  mutate(
    warning_category = case_when(
      grepl("fitted probabilities numerically 0 or 1", warning) ~ "logistic_separation",
      grepl("algorithm did not converge", warning) ~ "logistic_nonconvergence",
      grepl("glmnet C\\+\\+ code", warning) ~ "glmnet_late_lambda_nonconvergence",
      TRUE ~ "other"
    )
  ) |>
  count(model, warning_category, name = "n_warning_events") |>
  arrange(model, warning_category)
write_csv(warning_summary, file.path(output_dir, "model_fitting_warning_summary.csv"))

# ----------------------------------------------------------------------------
# 5. Compare grouped results with preserved legacy validation estimates
# ----------------------------------------------------------------------------
# The difference is descriptive: it quantifies the impact of patient grouping.
# It is not external validation and does not change the frozen deployed model.

legacy_sources <- tribble(
  ~cohort_group, ~legacy_file,
  "BM", "nested_bm_validation_updated5.rds",
  "Blood", "nested_blood_validation_updated5.rds",
  "FullCohort", "nested_fragmentomics_validation_updated3_original.rds",
  "BMRestricted", "nested_fragmentomics_bm_validation_updated5.rds",
  "BloodRestricted", "nested_fragmentomics_blood_validation_updated5.rds"
)

legacy_metrics <- pmap_dfr(legacy_sources, function(cohort_group, legacy_file) {
  legacy_path <- file.path(project_root, legacy_file)
  if (!file.exists(legacy_path)) stop("Missing legacy object: ", legacy_path, call. = FALSE)
  obj <- readRDS(legacy_path)
  obj$nested_metrics |>
    transmute(
      cohort_group = cohort_group,
      legacy_combo = combo,
      legacy_row_level_auc_mean = auc_mean,
      legacy_row_level_auc_sd = auc_sd
    )
})

model_map <- tribble(
  ~model, ~cohort_group, ~legacy_combo, ~paper_label,
  "BM_Sites", "BM", "BM_zscore_only_sites", "Sites Model",
  "BM_cVAF", "BM", "BM_zscore_only_detection_rate", "cVAF Model",
  "BM_Raw_cVAF", "BM", "BM_rate_only", "Raw cVAF",
  "BM_All_Mutation_Features", "BM", "BM_base", "All mutation features",
  "BM_Combined_Mutation_Zscores", "BM", "BM_base_zscore", "Combined Model",
  "BM_Mutation_Fragmentomics_Full", "BM", "BM_plus_fragment", "Mutation + fragmentomics (full)",
  "BM_Mutation_Fragmentomics_Min", "BM", "BM_plus_fragment_min", "Mutation + fragmentomics (minimal)",
  "Blood_Sites", "Blood", "Blood_zscore_only_sites", "Sites Model",
  "Blood_cVAF", "Blood", "Blood_zscore_only_detection_rate", "cVAF Model",
  "Blood_Raw_cVAF", "Blood", "Blood_rate_only", "Raw cVAF",
  "Blood_All_Mutation_Features", "Blood", "Blood_base", "All mutation features",
  "Blood_Combined_Mutation_Zscores", "Blood", "Blood_base_zscore", "Combined Model",
  "Blood_Mutation_Fragmentomics_Full", "Blood", "Blood_plus_fragment", "Mutation + fragmentomics (full)",
  "Blood_Mutation_Fragmentomics_Min", "Blood", "Blood_plus_fragment_min", "Mutation + fragmentomics (minimal)"
)

frag_suffixes <- tribble(
  ~suffix, ~legacy_combo, ~paper_label,
  "Full", "Fragmentomics_full", "Fragmentomics (full)",
  "Min", "Fragmentomics_min", "Fragmentomics (minimal)",
  "FS", "Fragmentomics_FS_only", "Fragment size score",
  "Mean_Coverage", "Fragmentomics_mean_coverage_only", "Mean coverage",
  "Proportion_Short", "Fragmentomics_prop_short_only", "Proportion short fragments",
  "Tumor_Fraction", "Fragmentomics_tumor_fraction_only", "Tumor fraction"
)

frag_groups <- tribble(
  ~model_prefix, ~cohort_group,
  "Fragmentomics_FullCohort_", "FullCohort",
  "Fragmentomics_BMRestricted_", "BMRestricted",
  "Fragmentomics_BloodRestricted_", "BloodRestricted"
)

frag_map <- tidyr::crossing(frag_groups, frag_suffixes) |>
  transmute(
    model = paste0(model_prefix, suffix),
    cohort_group,
    legacy_combo,
    paper_label
  )
model_map <- bind_rows(model_map, frag_map)

fold_auc_summary <- fold_metrics_tbl |>
  group_by(model) |>
  summarise(
    outer_fold_auc_mean = mean(auc, na.rm = TRUE),
    outer_fold_auc_sd = sd(auc, na.rm = TRUE),
    n_outer_fold_estimates = sum(is.finite(auc)),
    .groups = "drop"
  )

publication_summary <- summary_tbl |>
  select(-source_run_id) |>
  rename(
    repeat_pooled_auc_mean = auc_estimate,
    repeat_pooled_auc_sd = auc_split_sd,
    repeat_pooled_auc_q025 = auc_split_q025,
    repeat_pooled_auc_q975 = auc_split_q975
  ) |>
  left_join(fold_auc_summary, by = "model") |>
  left_join(model_map, by = "model") |>
  left_join(legacy_metrics, by = c("cohort_group", "legacy_combo")) |>
  mutate(
    auc_change_grouped_minus_legacy = repeat_pooled_auc_mean - legacy_row_level_auc_mean,
    grouped_auc_display = sprintf(
      "%.3f (patient-clustered 95%% CI %.3f-%.3f)",
      repeat_pooled_auc_mean, auc_cluster_boot_q025, auc_cluster_boot_q975
    )
  ) |>
  relocate(cohort_group, model, paper_label, legacy_combo) |>
  arrange(factor(cohort_group, levels = c("BM", "Blood", "FullCohort", "BMRestricted", "BloodRestricted")), desc(repeat_pooled_auc_mean))
write_csv(publication_summary, file.path(output_dir, "publication_model_performance_and_legacy_comparison.csv"), na = "")

# ----------------------------------------------------------------------------
# 6. Record the source runs and mark the combined run complete
# ----------------------------------------------------------------------------

run_manifest <- tibble(
  combined_run_id = output_run_id,
  source_run_id = source_run_ids,
  source_directory = source_dirs,
  source_complete = file.exists(file.path(source_dirs, "RUN_COMPLETE")),
  assembled_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)
write_csv(run_manifest, file.path(output_dir, "combined_run_manifest.csv"))

writeLines(
  c(
    "Complete patient-grouped repeated nested cross-validation model library",
    "32 model/cohort combinations",
    paste0(observed_outer_repeats, " repeated outer five-fold splits; 5 repeated inner five-fold splits"),
    "2,000 patient-clustered bootstrap replicates",
    "All combined quality-control checks passed.",
    paste0("Source runs: ", paste(source_run_ids, collapse = ", "))
  ),
  file.path(output_dir, "README.txt")
)
writeLines("PASS", file.path(output_dir, "RUN_COMPLETE"))

message("Assembled complete 32-model result set: ", output_dir)
print(
  publication_summary |>
    select(cohort_group, model, repeat_pooled_auc_mean, repeat_pooled_auc_sd,
           outer_fold_auc_mean, outer_fold_auc_sd, auc_cluster_boot_q025,
           auc_cluster_boot_q975, legacy_row_level_auc_mean,
           auc_change_grouped_minus_legacy),
  n = Inf
)
