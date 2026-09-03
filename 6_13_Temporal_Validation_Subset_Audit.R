#!/usr/bin/env Rscript

# =============================================================================
# 6_13_Temporal_Validation_Subset_Audit.R
#
# Goal
#   Separate the original seven-patient hold-out cohort, whose results informed
#   selection of displayed feature specifications, from the subsequently
#   accrued temporal validation patients. Quantify fixed-model performance in
#   each phase without changing model coefficients or decision thresholds.
#
# Inputs
#   * The current expanded-test sample manifest produced by
#     3_1C_Expanded_test_clustered_sensitivity.R.
#   * The pre-expansion patient cohort assignment. Its seven Non-frontline
#     patients define the original selection-informed hold-out set.
#
# Outputs
#   A new versioned directory under:
#     Output_tables_2025/temporal_validation_subset/<run-id>/
#   Existing expanded-test outputs and manuscript packages are never modified.
#
# Inference
#   Point estimates remain sample-level because the clinical comparison is
#   sample-level. Percentile intervals resample patients and retain all samples
#   belonging to every sampled patient. Earliest-evaluable-sample-per-patient
#   results are reported separately as a deterministic sensitivity analysis.
#
# Script roadmap
#   1. Lock the scored-manifest and pre-expansion cohort-assignment inputs.
#   2. Identify the original seven selection-informed hold-out patients.
#   3. Assign every current evaluable sample to original or later-accrued phase.
#   4. Compute fixed-model point estimates separately for each phase.
#   5. Bootstrap patients within phase and retain earliest-sample sensitivity.
#   6. Export exact membership, denominators, QC, and completion provenance.
#
# Why this matters
#   Results from the original hold-out influenced which feature specifications
#   were emphasized, so they are not equivalent to data accrued afterward.
#   This audit keeps both phases visible without changing frozen models or
#   thresholds and without discarding valid evaluable samples.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

# ----------------------------------------------------------------------------
# 1. Versioned settings and explicit input selection
# ----------------------------------------------------------------------------
# For a release run, pass --scored-manifest explicitly. Automatic latest-file
# discovery is convenient for local review but is weaker release provenance.
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) && length(args) > hit[1]) args[hit[1] + 1] else default
}
as_positive_integer <- function(value, label) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(label, " must be a positive integer; received: ", value, call. = FALSE)
  }
  parsed
}

RUN_ID <- get_arg("--run-id", format(Sys.Date(), "%Y-%m-%d_v1"))
BOOTSTRAP_REPS <- as_positive_integer(get_arg("--bootstrap-reps", "10000"),
                                      "--bootstrap-reps")
BASE_SEED <- as_positive_integer(get_arg("--seed", "20260731"), "--seed")
OUT_ROOT <- get_arg(
  "--output-root", file.path("Output_tables_2025", "temporal_validation_subset")
)
ORIGINAL_ASSIGNMENT_PATH <- get_arg(
  "--original-assignment", "Output_tables_2025/patient_cohort_assignment.csv"
)

default_manifests <- list.files(
  file.path("Output_tables_2025", "expanded_test_clustered_sensitivity"),
  pattern = "^expanded_test_repeated_measures_sample_manifest_[0-9-]+\\.csv$",
  full.names = TRUE
)
if (!length(default_manifests)) {
  stop("No expanded-test repeated-measures sample manifest was found.", call. = FALSE)
}
default_manifest <- default_manifests[which.max(file.info(default_manifests)$mtime)]
SCORED_MANIFEST_PATH <- get_arg("--scored-manifest", default_manifest)

if (!grepl("^[A-Za-z0-9._-]+$", RUN_ID)) {
  stop("--run-id contains unsupported characters.", call. = FALSE)
}
OUT_DIR <- file.path(OUT_ROOT, RUN_ID)
if (dir.exists(OUT_DIR) || file.exists(OUT_DIR)) {
  stop("Refusing to overwrite existing run directory: ", OUT_DIR, call. = FALSE)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(SCORED_MANIFEST_PATH, ORIGINAL_ASSIGNMENT_PATH)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing required input(s): ", paste(missing_paths, collapse = ", "),
       call. = FALSE)
}

file_manifest <- tibble(
  role = c("current_scored_manifest", "pre_expansion_cohort_assignment"),
  path = required_paths,
  size_bytes = as.numeric(file.info(required_paths)$size),
  modified_time = format(file.info(required_paths)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  md5 = unname(tools::md5sum(required_paths))
)
write_csv(file_manifest, file.path(OUT_DIR, "input_file_manifest.csv"))

# ----------------------------------------------------------------------------
# 2. Recover the original selection-informed seven-patient hold-out
# ----------------------------------------------------------------------------
# Failing unless exactly seven are recovered prevents a newer cohort assignment
# from silently redefining the historical comparison group.
original_assignment <- read_csv(ORIGINAL_ASSIGNMENT_PATH, show_col_types = FALSE)
required_assignment_columns <- c("Patient", "Cohort")
if (!all(required_assignment_columns %in% names(original_assignment))) {
  stop("Original cohort assignment lacks Patient or Cohort.", call. = FALSE)
}
original_holdout <- original_assignment %>%
  filter(Cohort == "Non-frontline") %>%
  distinct(Patient) %>%
  arrange(Patient)
if (nrow(original_holdout) != 7L) {
  stop("Expected exactly seven original Non-frontline hold-out patients; found ",
       nrow(original_holdout), ".", call. = FALSE)
}
write_csv(original_holdout, file.path(OUT_DIR, "original_seven_patient_ids.csv"))

scored <- read_csv(SCORED_MANIFEST_PATH, show_col_types = FALSE)
required_scored_columns <- c(
  "analysis_family", "model_key", "model", "threshold", "Patient",
  "Sample_Code", "truth", "probability", "prediction", "selected_earliest"
)
missing_columns <- setdiff(required_scored_columns, names(scored))
if (length(missing_columns)) {
  stop("Scored manifest lacks required columns: ",
       paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(scored[c("analysis_family", "model_key", "Patient", "Sample_Code")])) {
  stop("Scored manifest contains duplicate model/patient/sample rows.", call. = FALSE)
}
if (any(!scored$truth %in% c(0, 1)) || any(!scored$prediction %in% c(0, 1))) {
  stop("truth and prediction must be binary 0/1 values.", call. = FALSE)
}

# ----------------------------------------------------------------------------
# 3. Label validation phase while retaining all evaluable samples
# ----------------------------------------------------------------------------
scored <- scored %>%
  mutate(
    validation_phase = if_else(
      Patient %in% original_holdout$Patient,
      "original_seven_selection_informed",
      "newly_accrued_temporal_validation"
    ),
    truth_label = factor(truth, levels = c(0, 1), labels = c("neg", "pos")),
    prediction_label = factor(prediction, levels = c(0, 1), labels = c("neg", "pos"))
  )

if (!all(original_holdout$Patient %in% scored$Patient)) {
  missing_original <- setdiff(original_holdout$Patient, scored$Patient)
  # Not every original patient must be evaluable for every workflow, but every
  # original patient should appear somewhere in the combined scored manifest.
  stop("Original hold-out patient(s) absent from all scored workflows: ",
       paste(missing_original, collapse = ", "), call. = FALSE)
}

# ----------------------------------------------------------------------------
# 4. Compute fixed-model performance with sample and patient denominators
# ----------------------------------------------------------------------------
# Point estimates remain sample-level because that is the stated estimand.
# Both denominator types are exported so repeated samples cannot be mistaken
# for independent patients.
auc_fast <- function(truth, probability) {
  y <- as.integer(truth == 1)
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (!n_pos || !n_neg) return(NA_real_)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

metric_row <- function(data) {
  tp <- sum(data$truth == 1 & data$prediction == 1)
  fn <- sum(data$truth == 1 & data$prediction == 0)
  tn <- sum(data$truth == 0 & data$prediction == 0)
  fp <- sum(data$truth == 0 & data$prediction == 1)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  tibble(
    n_patients = n_distinct(data$Patient),
    n_samples = nrow(data),
    n_positive = tp + fn,
    n_negative = tn + fp,
    tp = tp, fn = fn, tn = tn, fp = fp,
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    accuracy = (tp + tn) / nrow(data),
    ppv = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    npv = if ((tn + fn) > 0) tn / (tn + fn) else NA_real_,
    auc = auc_fast(data$truth, data$probability)
  )
}

analysis_sets <- bind_rows(
  scored %>% mutate(analysis_set = "all_expanded_evaluable"),
  scored %>%
    filter(validation_phase == "original_seven_selection_informed") %>%
    mutate(analysis_set = "original_seven_selection_informed"),
  scored %>%
    filter(validation_phase == "newly_accrued_temporal_validation") %>%
    mutate(analysis_set = "newly_accrued_temporal_validation")
) %>%
  mutate(selection = "all_evaluable_samples")

earliest_sets <- analysis_sets %>%
  filter(selected_earliest) %>%
  mutate(selection = "earliest_evaluable_sample_per_patient")

all_sets <- bind_rows(analysis_sets, earliest_sets)

point_metrics <- all_sets %>%
  group_by(analysis_family, model_key, model, threshold, analysis_set, selection) %>%
  group_modify(~ metric_row(.x)) %>%
  ungroup()

# ----------------------------------------------------------------------------
# 5. Bootstrap patients within phase
# ----------------------------------------------------------------------------
# All samples for a resampled patient travel together. Earliest-evaluable
# results remain a separate sensitivity analysis rather than replacing the
# all-evaluable-sample estimand.
bootstrap_one_group <- function(data, reps, seed) {
  patients <- unique(as.character(data$Patient))
  rows_by_patient <- split(seq_len(nrow(data)), as.character(data$Patient))
  set.seed(seed)
  map_dfr(seq_len(reps), function(b) {
    sampled <- sample(patients, length(patients), replace = TRUE)
    idx <- unlist(rows_by_patient[sampled], use.names = FALSE)
    metric_row(data[idx, , drop = FALSE]) %>%
      select(sensitivity, specificity, balanced_accuracy, accuracy, ppv, npv, auc) %>%
      mutate(bootstrap_rep = b, .before = 1)
  })
}

group_keys <- all_sets %>%
  distinct(analysis_family, model_key, model, threshold, analysis_set, selection) %>%
  arrange(analysis_family, model_key, analysis_set, selection) %>%
  mutate(group_seed = BASE_SEED + row_number() * 100003L)

bootstrap_long <- pmap_dfr(group_keys, function(
    analysis_family, model_key, model, threshold, analysis_set, selection,
    group_seed) {
  data <- all_sets %>%
    filter(.data$analysis_family == .env$analysis_family,
           .data$model_key == .env$model_key,
           .data$analysis_set == .env$analysis_set,
           .data$selection == .env$selection)
  bootstrap_one_group(data, BOOTSTRAP_REPS, group_seed) %>%
    mutate(
      analysis_family = analysis_family, model_key = model_key, model = model,
      threshold = threshold, analysis_set = analysis_set, selection = selection,
      .before = 1
    )
})

metric_names <- c("sensitivity", "specificity", "balanced_accuracy", "accuracy",
                  "ppv", "npv", "auc")
bootstrap_intervals <- bootstrap_long %>%
  pivot_longer(cols = all_of(metric_names), names_to = "metric", values_to = "value") %>%
  group_by(analysis_family, model_key, model, threshold, analysis_set, selection,
           metric) %>%
  summarise(
    ci_lower = if (any(is.finite(value))) quantile(value, 0.025, na.rm = TRUE) else NA_real_,
    ci_upper = if (any(is.finite(value))) quantile(value, 0.975, na.rm = TRUE) else NA_real_,
    valid_replicates = sum(is.finite(value)),
    total_replicates = n(),
    .groups = "drop"
  )

# ----------------------------------------------------------------------------
# 6. Export phase membership, exact denominators, estimates, and hard QC
# ----------------------------------------------------------------------------
phase_manifest <- scored %>%
  select(analysis_family, model_key, model, threshold, validation_phase,
         Patient, Sample_Code, truth, probability, prediction, selected_earliest,
         everything())

cohort_counts <- scored %>%
  distinct(validation_phase, Patient) %>%
  count(validation_phase, name = "n_evaluable_patients_in_any_workflow")

qc <- tibble(
  check = c(
    "original_holdout_contains_exactly_seven_patients",
    "all_original_holdout_patients_appear_in_scored_manifest",
    "no_duplicate_model_patient_sample_rows",
    "fixed_threshold_is_constant_within_model",
    "temporal_validation_subset_is_nonempty"
  ),
  status = c(
    ifelse(nrow(original_holdout) == 7L, "PASS", "FAIL"),
    ifelse(all(original_holdout$Patient %in% scored$Patient), "PASS", "FAIL"),
    ifelse(!anyDuplicated(scored[c("analysis_family", "model_key", "Patient", "Sample_Code")]),
           "PASS", "FAIL"),
    ifelse(all(scored %>% group_by(model_key) %>%
                 summarise(n_thresholds = n_distinct(threshold), .groups = "drop") %>%
                 pull(n_thresholds) == 1L), "PASS", "FAIL"),
    ifelse(any(scored$validation_phase == "newly_accrued_temporal_validation"),
           "PASS", "FAIL")
  )
)
if (any(qc$status != "PASS")) stop("Temporal validation QC failed.", call. = FALSE)

write_csv(phase_manifest, file.path(OUT_DIR, "scored_sample_phase_manifest.csv"))
write_csv(cohort_counts, file.path(OUT_DIR, "evaluable_patient_counts_by_phase.csv"))
write_csv(point_metrics, file.path(OUT_DIR, "fixed_model_performance_by_phase.csv"))
write_csv(bootstrap_intervals,
          file.path(OUT_DIR, "patient_clustered_intervals_by_phase.csv"))
write_csv(bootstrap_long,
          file.path(OUT_DIR, "patient_clustered_bootstrap_replicates.csv"))
write_csv(qc, file.path(OUT_DIR, "quality_control_checks.csv"))
write_csv(
  tibble(
    parameter = c("run_id", "scored_manifest", "original_assignment",
                  "bootstrap_reps", "base_seed"),
    value = c(RUN_ID, SCORED_MANIFEST_PATH, ORIGINAL_ASSIGNMENT_PATH,
              BOOTSTRAP_REPS, BASE_SEED)
  ),
  file.path(OUT_DIR, "run_parameters.csv")
)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "session_info.txt"))
writeLines(
  c(
    "Temporal validation subset audit completed successfully.",
    paste("Run ID:", RUN_ID),
    paste("Output directory:", OUT_DIR),
    "Model probabilities and fixed thresholds were read from the current scored manifest and were not changed."
  ),
  file.path(OUT_DIR, "RUN_COMPLETE")
)

message("Temporal validation subset audit completed. Newly accrued results:")
print(point_metrics %>%
        filter(analysis_set == "newly_accrued_temporal_validation",
               selection == "all_evaluable_samples") %>%
        select(analysis_family, model_key, n_patients, n_samples, n_positive,
               n_negative, sensitivity, specificity, accuracy, auc))
message("Outputs written to: ", OUT_DIR)
