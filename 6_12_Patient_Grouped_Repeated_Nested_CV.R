#!/usr/bin/env Rscript

# =============================================================================
# 6_12_Patient_Grouped_Repeated_Nested_CV.R
#
# Goal
#   Estimate internal discrimination and operating-point performance for either
#   the frozen headline feature specifications or the complete original model
#   library while respecting repeated samples from the same patient in every
#   resampling layer.
#
# Scientific estimand
#   The estimand is sample-level MRD classification for a previously unseen
#   patient, conditional on each already-selected feature specification. This
#   script does not reselect the manuscript's feature combinations and does not
#   refit or overwrite the frozen models used for test-cohort scoring.
#
# Design
#   * Repeated five-fold outer cross-validation grouped by patient.
#   * Repeated five-fold inner cross-validation grouped by patient.
#   * Inner folds tune glmnet hyperparameters and generate out-of-fold
#     predictions used to select the outer-training Youden threshold.
#   * The threshold is then applied unchanged to the outer-held-out patients.
#   * Primary summaries pool held-out predictions within each outer repeat and
#     average the resulting repeat-level metrics.
#   * Percentile confidence intervals use a patient-clustered bootstrap that
#     resamples patients and retains all their samples and repeated-CV
#     predictions.
#
# Script roadmap
#   1. Parse and validate run settings.
#      Why: every result must record an explicit, reproducible resampling design.
#   2. Protect and hash frozen artifacts.
#      Why: this correction estimates validation performance; it must not alter
#      models or thresholds used for independent test-cohort scoring.
#   3. Reconstruct each model's historical training frame.
#      Why: grouped and legacy estimates are comparable only when they use the
#      same eligible rows, outcome, and preselected predictors.
#   4. Create outer and inner folds grouped by patient.
#      Why: all samples from one patient must remain on one side of every split.
#   5. Tune glmnet and select a threshold inside each outer-training set.
#      Why: held-out samples cannot contribute to tuning or threshold selection.
#   6. Score outer-held-out patients and pool predictions within each repeat.
#      Why: the primary AUC uses all held-out predictions in a repeat rather
#      than an average of potentially unstable fold-level AUCs.
#   7. Bootstrap patients and export QC/provenance tables.
#      Why: uncertainty must preserve within-patient dependence and every
#      estimate must remain traceable to its inputs, folds, and predictions.
#
# Outputs
#   A new, explicitly named run directory under:
#     Output_tables_2025/patient_grouped_repeated_nested_cv/<run-id>/
#   The script refuses to overwrite an existing run directory. Historical
#   row-level CV outputs, the July 29 grouped diagnostic, frozen models, and
#   frozen thresholds are read-only inputs and remain unchanged.
#
# Example
#   Rscript Scripts_2025/Final_Scripts/6_12_Patient_Grouped_Repeated_Nested_CV.R \
#     --run-id 2026-07-31_v1 --outer-repeats 100 --inner-repeats 5 \
#     --bootstrap-reps 2000
#
#   Complete original model library (32 model/cohort combinations):
#     Rscript Scripts_2025/Final_Scripts/6_12_Patient_Grouped_Repeated_Nested_CV.R \
#       --run-id 2026-07-31_all_models_v1 --model-library all \
#       --outer-repeats 50 --inner-repeats 5 --bootstrap-reps 2000
#
# Expected failure modes
#   The script stops if frozen artifacts are writable, canonical frame sizes
#   drift, a patient appears on both sides of a split, a fold lacks either MRD
#   class, inner out-of-fold predictions are incomplete, or an output directory
#   already exists.
# =============================================================================

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(glmnet)
  library(pROC)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

# ----------------------------------------------------------------------------
# 1. Run configuration: explicit arguments and versioned output directory
# ----------------------------------------------------------------------------
# Existing run directories are never overwritten because results generated
# with different seeds or resampling settings are different scientific objects.
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
OUTER_REPEATS <- as_positive_integer(get_arg("--outer-repeats", "100"),
                                     "--outer-repeats")
INNER_REPEATS <- as_positive_integer(get_arg("--inner-repeats", "5"),
                                     "--inner-repeats")
BOOTSTRAP_REPS <- as_positive_integer(get_arg("--bootstrap-reps", "2000"),
                                      "--bootstrap-reps")
BASE_SEED <- as_positive_integer(get_arg("--seed", "20260731"), "--seed")
MODEL_FILTER <- get_arg("--models", NULL)
MODEL_LIBRARY <- tolower(get_arg("--model-library", "headline"))
ALLOW_DIM_MISMATCH <- "--allow-dim-mismatch" %in% args
MANIFEST_ONLY <- "--manifest-only" %in% args
SHOW_WARNINGS <- "--show-warnings" %in% args
AGG_PATH <- get_arg(
  "--data",
  "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds"
)

if (!grepl("^[A-Za-z0-9._-]+$", RUN_ID)) {
  stop("--run-id may contain only letters, numbers, period, underscore, and hyphen.",
       call. = FALSE)
}
if (!MODEL_LIBRARY %in% c("headline", "all")) {
  stop("--model-library must be either 'headline' or 'all'.", call. = FALSE)
}

OUT_ROOT <- get_arg(
  "--output-root",
  file.path("Output_tables_2025", "patient_grouped_repeated_nested_cv")
)
OUT_DIR <- file.path(OUT_ROOT, RUN_ID)
if (dir.exists(OUT_DIR) || file.exists(OUT_DIR)) {
  stop("Refusing to overwrite existing run directory: ", OUT_DIR, call. = FALSE)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cleanup_incomplete_run <- function() {
  if (!file.exists(file.path(OUT_DIR, "RUN_COMPLETE"))) {
    message("Run did not complete. Partial diagnostic files remain in: ", OUT_DIR)
  }
}
on.exit(cleanup_incomplete_run(), add = TRUE)

# ----------------------------------------------------------------------------
# 2. Frozen-artifact guardrail and input provenance
# ----------------------------------------------------------------------------
# These files define the historical model boundary. They are hashed and must be
# read-only so validation cannot silently become retraining or change test calls.
frozen_artifacts <- c(
  "Output_tables_2025/selected_combo_models_2026-02-16.rds",
  "Output_tables_2025/selected_combo_thresholds_2026-02-16.rds",
  "nested_bm_validation_updated5.rds",
  "nested_blood_validation_updated5.rds",
  "nested_fragmentomics_validation_updated3_original.rds"
)

legacy_diagnostic_files <- list.files(
  file.path("Output_tables_2025", "patient_grouped_cv_diagnostic"),
  full.names = TRUE,
  recursive = TRUE
)

assert_frozen_artifacts_protected <- function(paths) {
  present <- paths[file.exists(paths)]
  if (!length(present)) {
    stop("No frozen model artifacts were found from the project root.", call. = FALSE)
  }
  writable <- present[file.access(present, mode = 2) == 0]
  if (length(writable)) {
    stop(
      "Refusing to run because frozen model artifacts are writable:\n",
      paste0("  ", writable, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(present)
}

hash_manifest <- function(paths, role) {
  paths <- unique(paths[file.exists(paths) & !dir.exists(paths)])
  if (!length(paths)) return(tibble())
  info <- file.info(paths)
  tibble(
    role = role,
    path = paths,
    size_bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    md5 = unname(tools::md5sum(paths))
  )
}

present_frozen <- assert_frozen_artifacts_protected(frozen_artifacts)
input_paths <- unique(c(AGG_PATH, "baseline_high_quality_patients_updated.csv",
                        present_frozen, legacy_diagnostic_files))
write_csv(
  bind_rows(
    hash_manifest(c(AGG_PATH, "baseline_high_quality_patients_updated.csv"),
                  "analysis_input"),
    hash_manifest(present_frozen, "frozen_artifact_preserved"),
    hash_manifest(legacy_diagnostic_files, "legacy_diagnostic_preserved")
  ),
  file.path(OUT_DIR, "input_and_legacy_artifact_manifest.csv")
)

if (!file.exists(AGG_PATH)) stop("Missing aggregate input: ", AGG_PATH, call. = FALSE)

# ----------------------------------------------------------------------------
# 3. Reconstruct the canonical training population and model specifications
# ----------------------------------------------------------------------------
# The headline library contains preselected manuscript models. The optional
# all-model library is a sensitivity expansion; it does not repeat feature
# selection or redefine the independent test cohort.

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) .helpers_path <- "helpers.R"
if (!file.exists(.helpers_path)) stop("Missing helpers.R", call. = FALSE)
source(.helpers_path)
rm(.helpers_path)

dat <- readRDS(AGG_PATH)
cohort_df <- load_final_cohort_assignment()
if (anyDuplicated(cohort_df$Patient)) {
  stop("Canonical cohort assignment contains duplicate Patient values.", call. = FALSE)
}

dat <- dat %>%
  left_join(cohort_df %>% select(Patient, Cohort), by = "Patient") %>%
  mutate(
    MRD_truth = case_when(
      !is.na(Adaptive_Binary) ~ Adaptive_Binary,
      is.na(Adaptive_Binary) & !is.na(Flow_Binary) ~ Flow_Binary,
      TRUE ~ NA_real_
    )
  )

dat_mrd <- dat %>% filter(!timepoint_info %in% c("Diagnosis", "Baseline"))
.demoted_path <- file.path(
  "Output_tables_2025", "clinical_support",
  "nonbaseline_timepoints_with_baseline_or_diagnosis_label_demoted_audit.csv"
)
if (file.exists(.demoted_path)) {
  frozen_non_mrd_rows <- read_csv(.demoted_path, show_col_types = FALSE) %>%
    filter(.data$timepoint_info %in% c("Diagnosis", "Baseline")) %>%
    distinct(.data$Patient, .data$Sample_Code)
  dat_mrd <- dat_mrd %>%
    anti_join(frozen_non_mrd_rows, by = c("Patient", "Sample_Code"))
}

train_df <- dat_mrd %>%
  filter(Cohort == "Frontline", !is.na(MRD_truth)) %>%
  mutate(MRD_truth = factor(MRD_truth, levels = c(0, 1), labels = c("neg", "pos")))

if (any(is.na(train_df$MRD_truth))) {
  stop("MRD_truth contains values other than 0/1 after filtering.", call. = FALSE)
}

Good_pts <- read.csv("baseline_high_quality_patients_updated.csv",
                     stringsAsFactors = FALSE)
bm_good <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_BM_cells == 1) %>%
  distinct(Patient)
blood_good <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
  distinct(Patient)

bm_eligibility_preds <- c(
  "zscore_BM", "z_score_detection_rate_BM", "detect_rate_BM",
  "FS", "Mean.Coverage", "Proportion.Short",
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
)
blood_eligibility_preds <- c(
  "zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood",
  "FS", "Mean.Coverage", "Proportion.Short",
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
)

fragmentomics_eligibility_preds <- c(
  "FS", "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
  "Proportion.Short"
)

SPECS <- list(
  BM_Sites = list(
    label = "BM mutant-sites z-score",
    predictors = "zscore_BM",
    eligibility_predictors = bm_eligibility_preds,
    eligible_patients = bm_good$Patient,
    frozen_reference_path = "nested_bm_validation_updated5.rds",
    frozen_reference_model = "BM_plus_fragment",
    membership_key_predictors = "zscore_BM",
    expected_samples = 42L,
    expected_patients = 26L
  ),
  BM_cVAF = list(
    label = "BM cVAF z-score",
    predictors = "z_score_detection_rate_BM",
    eligibility_predictors = bm_eligibility_preds,
    eligible_patients = bm_good$Patient,
    frozen_reference_path = "nested_bm_validation_updated5.rds",
    frozen_reference_model = "BM_plus_fragment",
    membership_key_predictors = "zscore_BM",
    expected_samples = 42L,
    expected_patients = 26L
  ),
  Blood_Combined = list(
    label = "Baseline-plasma mutation plus fragmentomics",
    predictors = c(
      "zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood",
      "FS", "Mean.Coverage", "Proportion.Short",
      "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
    ),
    eligibility_predictors = blood_eligibility_preds,
    eligible_patients = blood_good$Patient,
    frozen_reference_path = "nested_blood_validation_updated5.rds",
    frozen_reference_model = "Blood_plus_fragment",
    membership_key_predictors = "zscore_blood",
    expected_samples = 44L,
    expected_patients = 28L
  ),
  Blood_Sites = list(
    label = "Baseline-plasma mutant-sites z-score",
    predictors = "zscore_blood",
    eligibility_predictors = blood_eligibility_preds,
    eligible_patients = blood_good$Patient,
    frozen_reference_path = "nested_blood_validation_updated5.rds",
    frozen_reference_model = "Blood_plus_fragment",
    membership_key_predictors = "zscore_blood",
    expected_samples = 44L,
    expected_patients = 28L
  )
)

if (MODEL_LIBRARY == "all") {
  bm_common <- list(
    eligibility_predictors = bm_eligibility_preds,
    eligible_patients = bm_good$Patient,
    frozen_reference_path = "nested_bm_validation_updated5.rds",
    frozen_reference_model = "BM_plus_fragment",
    membership_key_predictors = "zscore_BM",
    expected_samples = 42L,
    expected_patients = 26L
  )
  blood_common <- list(
    eligibility_predictors = blood_eligibility_preds,
    eligible_patients = blood_good$Patient,
    frozen_reference_path = "nested_blood_validation_updated5.rds",
    frozen_reference_model = "Blood_plus_fragment",
    membership_key_predictors = "zscore_blood",
    expected_samples = 44L,
    expected_patients = 28L
  )

  make_spec <- function(label, predictors, common) {
    c(list(label = label, predictors = predictors), common)
  }

  bm_specs <- list(
    BM_Sites = make_spec(
      "BM mutant-sites z-score", "zscore_BM", bm_common
    ),
    BM_cVAF = make_spec(
      "BM cVAF z-score", "z_score_detection_rate_BM", bm_common
    ),
    BM_Raw_cVAF = make_spec(
      "BM raw cVAF", "detect_rate_BM", bm_common
    ),
    BM_All_Mutation_Features = make_spec(
      "BM sites z-score plus raw and standardized cVAF",
      c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM"),
      bm_common
    ),
    BM_Combined_Mutation_Zscores = make_spec(
      "BM sites and cVAF z-scores",
      c("zscore_BM", "z_score_detection_rate_BM"), bm_common
    ),
    BM_Mutation_Fragmentomics_Full = make_spec(
      "BM mutation plus full fragmentomics",
      c(
        "zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM",
        "FS", "Mean.Coverage", "Proportion.Short",
        "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
      ),
      bm_common
    ),
    BM_Mutation_Fragmentomics_Min = make_spec(
      "BM mutation plus minimal fragmentomics",
      c(
        "zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM",
        "FS", "Mean.Coverage"
      ),
      bm_common
    )
  )

  blood_specs <- list(
    Blood_Sites = make_spec(
      "Baseline-plasma mutant-sites z-score", "zscore_blood", blood_common
    ),
    Blood_cVAF = make_spec(
      "Baseline-plasma cVAF z-score", "z_score_detection_rate_blood",
      blood_common
    ),
    Blood_Raw_cVAF = make_spec(
      "Baseline-plasma raw cVAF", "detect_rate_blood", blood_common
    ),
    Blood_All_Mutation_Features = make_spec(
      "Baseline-plasma sites z-score plus raw and standardized cVAF",
      c("zscore_blood", "detect_rate_blood", "z_score_detection_rate_blood"),
      blood_common
    ),
    Blood_Combined_Mutation_Zscores = make_spec(
      "Baseline-plasma sites and cVAF z-scores",
      c("zscore_blood", "z_score_detection_rate_blood"), blood_common
    ),
    Blood_Mutation_Fragmentomics_Full = make_spec(
      "Baseline-plasma mutation plus full fragmentomics",
      c(
        "zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood",
        "FS", "Mean.Coverage", "Proportion.Short",
        "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
      ),
      blood_common
    ),
    Blood_Mutation_Fragmentomics_Min = make_spec(
      "Baseline-plasma mutation plus minimal fragmentomics",
      c(
        "zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood",
        "FS", "Mean.Coverage"
      ),
      blood_common
    )
  )

  fragmentomics_predictors <- list(
    Full = fragmentomics_eligibility_preds,
    Min = c("FS", "Mean.Coverage"),
    FS = "FS",
    Mean_Coverage = "Mean.Coverage",
    Proportion_Short = "Proportion.Short",
    Tumor_Fraction = "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
  )

  make_fragmentomics_specs <- function(prefix, cohort_label, reference_path,
                                        eligible_patients, expected_samples,
                                        expected_patients) {
    common <- list(
      eligibility_predictors = fragmentomics_eligibility_preds,
      eligible_patients = eligible_patients,
      frozen_reference_path = reference_path,
      frozen_reference_model = "Fragmentomics_full",
      # Tumour-fraction values for one historical row were recalibrated after
      # the frozen run. Match membership on the three unchanged fragmentomic
      # measurements, then restore every predictor (including the historical
      # tumour fraction) from the preserved trainingData below.
      membership_key_predictors = c(
        "FS", "Mean.Coverage", "Proportion.Short"
      ),
      expected_samples = expected_samples,
      expected_patients = expected_patients
    )
    imap(fragmentomics_predictors, function(predictors, model_suffix) {
      make_spec(
        paste0("Fragmentomics-only ", model_suffix, "; ", cohort_label),
        predictors, common
      )
    }) %>% set_names(paste0(prefix, "_", names(fragmentomics_predictors)))
  }

  fragmentomics_full_specs <- make_fragmentomics_specs(
    prefix = "Fragmentomics_FullCohort",
    cohort_label = "full training cohort",
    reference_path = "nested_fragmentomics_validation_updated3_original.rds",
    eligible_patients = unique(train_df$Patient),
    expected_samples = 65L,
    expected_patients = 41L
  )
  fragmentomics_bm_specs <- make_fragmentomics_specs(
    prefix = "Fragmentomics_BMRestricted",
    cohort_label = "BM-restricted training cohort",
    reference_path = "nested_fragmentomics_bm_validation_updated5.rds",
    eligible_patients = bm_good$Patient,
    expected_samples = 42L,
    expected_patients = 26L
  )
  fragmentomics_blood_specs <- make_fragmentomics_specs(
    prefix = "Fragmentomics_BloodRestricted",
    cohort_label = "baseline-plasma-restricted training cohort",
    reference_path = "nested_fragmentomics_blood_validation_updated5.rds",
    eligible_patients = blood_good$Patient,
    expected_samples = 44L,
    expected_patients = 28L
  )

  SPECS <- c(
    bm_specs, blood_specs,
    fragmentomics_full_specs,
    fragmentomics_bm_specs,
    fragmentomics_blood_specs
  )
  if (length(SPECS) != 32L || anyDuplicated(names(SPECS))) {
    stop("Complete model library must contain 32 uniquely named specifications.",
         call. = FALSE)
  }
}

if (!is.null(MODEL_FILTER)) {
  requested_models <- trimws(strsplit(MODEL_FILTER, ",", fixed = TRUE)[[1]])
  unknown_models <- setdiff(requested_models, names(SPECS))
  if (length(unknown_models)) {
    stop("Unknown --models value(s): ", paste(unknown_models, collapse = ", "),
         ". Available models: ", paste(names(SPECS), collapse = ", "),
         call. = FALSE)
  }
  SPECS <- SPECS[requested_models]
}

# ----------------------------------------------------------------------------
# 4. Match each specification to its frozen historical training frame
# ----------------------------------------------------------------------------
# Matching uses the outcome and unchanged predictor values, with occurrence
# indices for duplicate keys. This makes every excluded or unmatched row
# reviewable and holds the analysis population constant across CV designs.
make_membership_key <- function(outcome, data, predictors) {
  formatted <- lapply(data[predictors], function(x) sprintf("%.10g", as.numeric(x)))
  do.call(paste, c(list(as.character(outcome)), formatted, sep = "\r"))
}

occurrence_within_key <- function(key) {
  ave(seq_along(key), key, FUN = seq_along)
}

build_frame <- function(spec) {
  required <- unique(c("Patient", "Sample_Code", "MRD_truth",
                       spec$predictors, spec$eligibility_predictors))
  missing <- setdiff(required, names(train_df))
  if (length(missing)) {
    stop("Missing required training-frame columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  reconstructed <- train_df %>%
    drop_na(all_of(c("MRD_truth", spec$eligibility_predictors))) %>%
    filter(Patient %in% spec$eligible_patients) %>%
    select(all_of(required)) %>%
    mutate(.row_id = row_number())
  if (anyDuplicated(reconstructed[c("Patient", "Sample_Code")])) {
    stop("Training frame contains duplicate Patient/Sample_Code rows.", call. = FALSE)
  }

  # Recover the exact sample membership of the frozen development run. The
  # current aggregate table can acquire later longitudinal samples from a
  # patient who was already in the training cohort. Such samples must not enter
  # a retrospective correction of the frozen training analysis. The preserved
  # full-feature caret object supplies the exact outcome/predictor rows; the
  # current table supplies their patient and sample identifiers.
  frozen_object <- readRDS(spec$frozen_reference_path)
  frozen_training <- frozen_object$models[[spec$frozen_reference_model]]$trainingData
  if (is.null(frozen_training)) {
    stop("Frozen reference model lacks trainingData: ",
         spec$frozen_reference_path, " / ", spec$frozen_reference_model,
         call. = FALSE)
  }
  reference_predictors <- setdiff(names(frozen_training), ".outcome")
  if (!setequal(reference_predictors, spec$eligibility_predictors)) {
    stop("Frozen reference predictors do not match the canonical eligibility set.",
         call. = FALSE)
  }

  if (!all(spec$membership_key_predictors %in% reference_predictors)) {
    stop("Membership-key predictor is absent from the frozen reference model.",
         call. = FALSE)
  }
  current_key <- make_membership_key(
    reconstructed$MRD_truth, reconstructed, spec$membership_key_predictors
  )
  frozen_key <- make_membership_key(
    frozen_training$.outcome, frozen_training, spec$membership_key_predictors
  )
  current_composite <- paste(current_key, occurrence_within_key(current_key), sep = "\f")
  frozen_composite <- paste(frozen_key, occurrence_within_key(frozen_key), sep = "\f")
  matched_idx <- match(frozen_composite, current_composite)
  if (anyNA(matched_idx)) {
    stop(
      "The current aggregate table cannot reproduce every row in the frozen ",
      "training frame for ", spec$frozen_reference_model,
      ". Predictor values or eligibility have drifted; review before continuing.",
      call. = FALSE
    )
  }

  frame <- reconstructed[matched_idx, , drop = FALSE]
  if (!all(as.character(frame$MRD_truth) == as.character(frozen_training$.outcome))) {
    stop("Frozen membership matching produced an outcome mismatch.", call. = FALSE)
  }
  # Use the preserved predictor values themselves, not recalculated versions
  # from the current aggregate table. This makes the corrected resampling an
  # exact reanalysis of the frozen development data while retaining current
  # patient/sample identifiers solely for grouping and provenance.
  for (predictor in reference_predictors) {
    frame[[predictor]] <- frozen_training[[predictor]]
  }
  frame <- frame %>% mutate(.row_id = row_number())
  excluded_idx <- setdiff(seq_len(nrow(reconstructed)), matched_idx)
  exclusions <- reconstructed[excluded_idx, , drop = FALSE] %>%
    transmute(
      Patient, Sample_Code, MRD_truth,
      exclusion_reason = "present in current reconstructed frame but absent from frozen trainingData"
    )
  attr(frame, "membership_exclusions") <- exclusions
  frame
}

FRAME_CACHE <- imap(SPECS, function(spec, model_name) {
  frame <- build_frame(spec)
  list(
    frame = frame,
    exclusions = attr(frame, "membership_exclusions") %>%
      mutate(model = model_name, .before = 1)
  )
})

membership_exclusions <- map_dfr(FRAME_CACHE, "exclusions")
write_csv(membership_exclusions,
          file.path(OUT_DIR, "current_rows_excluded_to_match_frozen_training.csv"))

if (MANIFEST_ONLY) {
  manifest_only <- imap_dfr(SPECS, function(spec, model_name) {
    FRAME_CACHE[[model_name]]$frame %>%
      transmute(model = model_name, Patient, Sample_Code, MRD_truth,
                across(all_of(spec$predictors)))
  })
  write_csv(manifest_only, file.path(OUT_DIR, "training_frame_manifest_only.csv"))
  writeLines(
    c("Training-frame manifest completed without model fitting.",
      paste("Run ID:", RUN_ID), paste("Output directory:", OUT_DIR)),
    file.path(OUT_DIR, "RUN_COMPLETE")
  )
  message("Manifest-only run completed: ", OUT_DIR)
  quit(save = "no", status = 0)
}

# Randomized greedy assignment minimizes imbalance in positive samples,
# negative samples, total samples, and patient counts while retaining whole
# patients. Multiple candidate assignments are generated and the best valid
# one is returned deterministically for the supplied seed.
#
# ----------------------------------------------------------------------------
# 5. Build leakage-free outer and inner patient folds
# ----------------------------------------------------------------------------
make_grouped_folds <- function(frame, k = 5L, seed, attempts = 500L) {
  group_summary <- frame %>%
    transmute(Patient = as.character(Patient),
              is_pos = MRD_truth == "pos", is_neg = MRD_truth == "neg") %>%
    group_by(Patient) %>%
    summarise(n_pos = sum(is_pos), n_neg = sum(is_neg), n = n(), .groups = "drop")

  if (nrow(group_summary) < k) {
    stop("Cannot create ", k, " patient-grouped folds from ", nrow(group_summary),
         " patients.", call. = FALSE)
  }

  target <- c(
    pos = sum(group_summary$n_pos) / k,
    neg = sum(group_summary$n_neg) / k,
    n = sum(group_summary$n) / k,
    patients = nrow(group_summary) / k
  )

  best <- NULL
  best_score <- Inf
  set.seed(seed)

  for (attempt in seq_len(attempts)) {
    difficulty <- pmax(
      group_summary$n_pos / max(target["pos"], 1),
      group_summary$n_neg / max(target["neg"], 1),
      group_summary$n / max(target["n"], 1)
    )
    ord <- order(-(difficulty + runif(nrow(group_summary), 0, 0.05)))
    gs <- group_summary[ord, , drop = FALSE]
    tallies <- matrix(0, nrow = k, ncol = 4,
                      dimnames = list(NULL, c("pos", "neg", "n", "patients")))
    assignment <- integer(nrow(gs))

    for (i in seq_len(nrow(gs))) {
      addition <- c(gs$n_pos[i], gs$n_neg[i], gs$n[i], 1)
      candidate_scores <- vapply(seq_len(k), function(fold_id) {
        proposed <- tallies
        proposed[fold_id, ] <- proposed[fold_id, ] + addition
        sum(((sweep(proposed, 2, target, "-") /
                pmax(matrix(target, nrow = k, ncol = 4, byrow = TRUE), 1))^2) %*%
              c(1.5, 1.5, 0.4, 0.2))
      }, numeric(1))
      choices <- which(candidate_scores == min(candidate_scores))
      chosen <- sample(choices, 1)
      assignment[i] <- chosen
      tallies[chosen, ] <- tallies[chosen, ] + addition
    }

    gs$fold <- assignment
    valid <- all(tallies[, "pos"] > 0) && all(tallies[, "neg"] > 0) &&
      all(tallies[, "patients"] > 0)
    if (!valid) next

    score <- sum(((sweep(tallies, 2, target, "-") /
                     pmax(matrix(target, nrow = k, ncol = 4, byrow = TRUE), 1))^2) %*%
                   c(1.5, 1.5, 0.4, 0.2))
    if (score < best_score) {
      best <- gs %>% select(Patient, fold)
      best_score <- score
    }
  }

  if (is.null(best)) {
    stop("Unable to construct grouped folds containing both MRD classes. ",
         "Consider fewer folds only after reviewing the patient-level outcome table.",
         call. = FALSE)
  }

  fold_id <- best$fold[match(as.character(frame$Patient), best$Patient)]
  if (anyNA(fold_id)) stop("Internal fold assignment failure.", call. = FALSE)

  train_indices <- lapply(seq_len(k), function(i) which(fold_id != i))
  test_indices <- lapply(seq_len(k), function(i) which(fold_id == i))
  names(train_indices) <- names(test_indices) <- sprintf("Fold%02d", seq_len(k))

  for (i in seq_len(k)) {
    leaked <- intersect(unique(frame$Patient[train_indices[[i]]]),
                        unique(frame$Patient[test_indices[[i]]]))
    if (length(leaked)) stop("Patient leakage detected in grouped folds.", call. = FALSE)
    if (nlevels(droplevels(frame$MRD_truth[test_indices[[i]]])) < 2L) {
      stop("A grouped held-out fold lacks one MRD class.", call. = FALSE)
    }
  }

  list(index = train_indices, indexOut = test_indices,
       patient_assignment = best, balance_score = best_score)
}

make_repeated_inner_indices <- function(frame, repeats, seed) {
  index <- list()
  index_out <- list()
  assignments <- list()
  for (r in seq_len(repeats)) {
    folds <- make_grouped_folds(frame, k = 5L,
                                seed = seed + r * 1009L, attempts = 300L)
    for (f in seq_along(folds$index)) {
      nm <- sprintf("Repeat%02d.Fold%02d", r, f)
      index[[nm]] <- folds$index[[f]]
      index_out[[nm]] <- folds$indexOut[[f]]
    }
    assignments[[r]] <- folds$patient_assignment %>% mutate(inner_repeat = r)
  }
  list(index = index, indexOut = index_out,
       assignments = bind_rows(assignments))
}

# ----------------------------------------------------------------------------
# 6. Tune inside outer training data and select an inner-OOF threshold
# ----------------------------------------------------------------------------
# Tuning uses only outer-training patients. The Youden threshold comes from
# inner out-of-fold predictions and is then applied unchanged to outer-held-out
# patients, keeping the outer evaluation independent.
filter_final_tuning_predictions <- function(fit) {
  pred <- fit$pred
  if (is.null(pred) || !nrow(pred)) {
    stop("caret did not retain inner out-of-fold predictions.", call. = FALSE)
  }
  tune_cols <- intersect(names(fit$bestTune), names(pred))
  if (length(tune_cols)) {
    for (nm in tune_cols) {
      target <- fit$bestTune[[nm]][1]
      if (is.numeric(target)) {
        pred <- pred[abs(pred[[nm]] - target) < sqrt(.Machine$double.eps), , drop = FALSE]
      } else {
        pred <- pred[pred[[nm]] == target, , drop = FALSE]
      }
    }
  }
  pred
}

finite_youden_threshold <- function(obs, prob) {
  if (any(!is.finite(prob)) || any(prob < 0 | prob > 1)) {
    stop("Youden threshold requires finite probabilities in [0, 1].",
         call. = FALSE)
  }
  if (nlevels(droplevels(obs)) < 2L) {
    stop("Youden threshold requires both outcome classes.", call. = FALSE)
  }
  values <- sort(unique(prob))
  candidates <- if (length(values) == 1L) {
    c(0, 1)
  } else {
    unique(c(0, (values[-1] + values[-length(values)]) / 2, 1))
  }
  youden <- vapply(candidates, function(threshold) {
    predicted <- prob >= threshold
    positive <- obs == "pos"
    sensitivity <- sum(predicted & positive) / sum(positive)
    specificity <- sum(!predicted & !positive) / sum(!positive)
    sensitivity + specificity - 1
  }, numeric(1))
  best <- candidates[youden == max(youden)]
  # Resolve exact ties deterministically without privileging an extreme cutoff.
  threshold <- best[which.min(abs(best - stats::median(values)))]
  if (!is.finite(threshold)) {
    stop("Finite Youden threshold calculation failed.", call. = FALSE)
  }
  threshold
}

fit_outer_training_model <- function(frame, predictors, inner_repeats, seed) {
  inner <- make_repeated_inner_indices(frame, repeats = inner_repeats, seed = seed)
  ctrl <- trainControl(
    method = "cv",
    index = inner$index,
    indexOut = inner$indexOut,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    allowParallel = FALSE
  )
  form <- reformulate(predictors, response = "MRD_truth")
  set.seed(seed + 900001L)
  captured_warnings <- character()
  fit <- withCallingHandlers(
    {
      if (length(predictors) == 1L) {
        train(
          form, data = frame, method = "glm", family = binomial(), metric = "ROC",
          trControl = ctrl, preProcess = c("center", "scale")
        )
      } else {
        train(
          form, data = frame, method = "glmnet", metric = "ROC",
          tuneGrid = expand.grid(
            alpha = seq(0, 1, by = 0.25),
            lambda = 10^seq(-3, 1, length.out = 50)
          ),
          trControl = ctrl, preProcess = c("center", "scale")
        )
      }
    },
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      if (SHOW_WARNINGS) message("Training warning: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  inner_pred <- filter_final_tuning_predictions(fit) %>%
    group_by(rowIndex) %>%
    summarise(obs = first(obs), prob = mean(pos), n_predictions = n(), .groups = "drop")

  if (nrow(inner_pred) != nrow(frame) ||
      any(inner_pred$n_predictions != inner_repeats)) {
    stop("Inner grouped out-of-fold predictions are incomplete: expected ",
         inner_repeats, " per outer-training row.", call. = FALSE)
  }

  roc_obj <- pROC::roc(inner_pred$obs, inner_pred$prob, quiet = TRUE,
                       levels = c("neg", "pos"), direction = "<")
  threshold <- finite_youden_threshold(inner_pred$obs, inner_pred$prob)

  list(fit = fit, threshold = threshold,
       inner_auc = as.numeric(pROC::auc(roc_obj)),
       inner_predictions = inner_pred,
       inner_assignments = inner$assignments,
       warnings = captured_warnings)
}

auc_fast <- function(obs, prob) {
  y <- as.integer(obs == "pos")
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (!n_pos || !n_neg) return(NA_real_)
  ranks <- rank(prob, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

metric_row <- function(obs, prob, call) {
  pos <- obs == "pos"
  pred_pos <- call == "pos"
  tp <- sum(pos & pred_pos)
  fn <- sum(pos & !pred_pos)
  tn <- sum(!pos & !pred_pos)
  fp <- sum(!pos & pred_pos)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  tibble(
    auc = auc_fast(obs, prob),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    accuracy = (tp + tn) / (tp + tn + fp + fn),
    brier = mean((as.integer(obs == "pos") - prob)^2),
    tp = tp, fn = fn, tn = tn, fp = fp
  )
}

cluster_bootstrap_metrics <- function(predictions, reps, seed) {
  patients <- unique(as.character(predictions$Patient))
  rows_by_patient <- split(seq_len(nrow(predictions)), as.character(predictions$Patient))
  repeats <- sort(unique(predictions$outer_repeat))
  set.seed(seed)

  one_bootstrap <- function(b) {
    sampled <- sample(patients, length(patients), replace = TRUE)
    idx <- unlist(rows_by_patient[sampled], use.names = FALSE)
    boot <- predictions[idx, , drop = FALSE]
    by_repeat <- lapply(repeats, function(r) {
      x <- boot[boot$outer_repeat == r, , drop = FALSE]
      metric_row(x$obs, x$prob, x$call)
    })
    bind_rows(by_repeat) %>%
      summarise(across(c(auc, sensitivity, specificity, balanced_accuracy,
                         accuracy, brier), ~ mean(.x, na.rm = TRUE))) %>%
      mutate(bootstrap_rep = b, .before = 1)
  }

  map_dfr(seq_len(reps), one_bootstrap)
}

# ----------------------------------------------------------------------------
# 7. Execute repeated outer validation and retain row-level provenance
# ----------------------------------------------------------------------------
# Every prediction is held out at the patient level. Fold assignments, tuning
# results, warnings, and probabilities are retained so summary estimates can be
# independently reconstructed rather than accepted as opaque scalar outputs.
all_predictions <- list()
all_folds <- list()
all_inner_tuning <- list()
all_inner_assignments <- list()
all_training_warnings <- list()
all_outer_assignments <- list()
frame_manifests <- list()

resampling_group_index <- function(model_name) {
  if (startsWith(model_name, "BM_") ||
      startsWith(model_name, "Fragmentomics_BMRestricted_")) return(1L)
  if (startsWith(model_name, "Blood_") ||
      startsWith(model_name, "Fragmentomics_BloodRestricted_")) return(2L)
  if (startsWith(model_name, "Fragmentomics_FullCohort_")) return(3L)
  stop("Unable to assign resampling group for model: ", model_name,
       call. = FALSE)
}

for (model_name in names(SPECS)) {
  spec <- SPECS[[model_name]]
  frame <- FRAME_CACHE[[model_name]]$frame
  n_samples <- nrow(frame)
  n_patients <- n_distinct(frame$Patient)
  if (!ALLOW_DIM_MISMATCH &&
      (n_samples != spec$expected_samples || n_patients != spec$expected_patients)) {
    stop(
      model_name, " frame has ", n_samples, " samples from ", n_patients,
      " patients; expected ", spec$expected_samples, "/", spec$expected_patients,
      ". Resolve the drift or use --allow-dim-mismatch with a documented reason.",
      call. = FALSE
    )
  }

  frame_manifests[[model_name]] <- frame %>%
    transmute(model = model_name, Patient, Sample_Code, MRD_truth,
              repeated_patient = ave(.row_id, Patient, FUN = length) > 1L)
  message("\n", model_name, ": ", n_samples, " samples from ", n_patients,
          " patients; starting ", OUTER_REPEATS, " grouped outer repeats.")

  for (outer_repeat in seq_len(OUTER_REPEATS)) {
    repeat_seed <- BASE_SEED + resampling_group_index(model_name) * 10000000L +
      outer_repeat * 10000L
    outer <- make_grouped_folds(frame, k = 5L, seed = repeat_seed, attempts = 500L)
    outer_assignment <- outer$patient_assignment %>%
      mutate(model = model_name, outer_repeat = outer_repeat,
             balance_score = outer$balance_score, .before = 1)
    all_outer_assignments[[length(all_outer_assignments) + 1L]] <- outer_assignment

    for (outer_fold in seq_along(outer$index)) {
      train_idx <- outer$index[[outer_fold]]
      test_idx <- outer$indexOut[[outer_fold]]
      tr <- frame[train_idx, , drop = FALSE]
      te <- frame[test_idx, , drop = FALSE]
      leaked <- intersect(unique(tr$Patient), unique(te$Patient))
      if (length(leaked)) stop("Outer patient leakage detected.", call. = FALSE)

      fit_result <- fit_outer_training_model(
        tr, predictors = spec$predictors, inner_repeats = INNER_REPEATS,
        seed = repeat_seed + outer_fold * 101L
      )
      prob <- predict(fit_result$fit, newdata = te, type = "prob")[["pos"]]
      call <- factor(ifelse(prob >= fit_result$threshold, "pos", "neg"),
                     levels = c("neg", "pos"))
      fold_metrics <- metric_row(te$MRD_truth, prob, call) %>%
        mutate(
          model = model_name, outer_repeat = outer_repeat, outer_fold = outer_fold,
          threshold = fit_result$threshold, inner_auc = fit_result$inner_auc,
          n_train_samples = nrow(tr), n_train_patients = n_distinct(tr$Patient),
          n_test_samples = nrow(te), n_test_patients = n_distinct(te$Patient),
          leaked_patients = length(leaked), .before = 1
        )
      all_folds[[length(all_folds) + 1L]] <- fold_metrics

      pred_tbl <- tibble(
        model = model_name,
        outer_repeat = outer_repeat,
        outer_fold = outer_fold,
        Patient = te$Patient,
        Sample_Code = te$Sample_Code,
        obs = te$MRD_truth,
        prob = prob,
        threshold = fit_result$threshold,
        call = call
      )
      all_predictions[[length(all_predictions) + 1L]] <- pred_tbl

      tuning <- as_tibble(fit_result$fit$bestTune) %>%
        mutate(model = model_name, outer_repeat = outer_repeat,
               outer_fold = outer_fold, inner_auc = fit_result$inner_auc,
               threshold = fit_result$threshold, .before = 1)
      all_inner_tuning[[length(all_inner_tuning) + 1L]] <- tuning
      inner_assignment <- fit_result$inner_assignments %>%
        mutate(model = model_name, outer_repeat = outer_repeat,
               outer_fold = outer_fold, .before = 1)
      all_inner_assignments[[length(all_inner_assignments) + 1L]] <- inner_assignment
      if (length(fit_result$warnings)) {
        all_training_warnings[[length(all_training_warnings) + 1L]] <- tibble(
          model = model_name,
          outer_repeat = outer_repeat,
          outer_fold = outer_fold,
          warning = fit_result$warnings
        )
      }
    }

    if (outer_repeat %% max(1L, min(5L, OUTER_REPEATS)) == 0L ||
        outer_repeat == OUTER_REPEATS) {
      message(model_name, ": completed outer repeat ", outer_repeat, "/",
              OUTER_REPEATS)
    }
  }
}

predictions <- bind_rows(all_predictions)
fold_metrics <- bind_rows(all_folds)
inner_tuning <- bind_rows(all_inner_tuning)
inner_assignments <- bind_rows(all_inner_assignments)
training_warnings <- bind_rows(all_training_warnings)
outer_assignments <- bind_rows(all_outer_assignments)
frame_manifest <- bind_rows(frame_manifests)

expected_predictions <- frame_manifest %>% count(model, name = "n") %>%
  mutate(expected = n * OUTER_REPEATS)
observed_predictions <- predictions %>% count(model, name = "observed")
prediction_check <- left_join(expected_predictions, observed_predictions, by = "model")
if (any(prediction_check$expected != prediction_check$observed)) {
  stop("Each sample must have exactly one held-out prediction per outer repeat.",
       call. = FALSE)
}
if (any(fold_metrics$leaked_patients != 0L)) {
  stop("Patient leakage was detected in final fold metrics.", call. = FALSE)
}

repeat_metrics <- predictions %>%
  group_by(model, outer_repeat) %>%
  group_modify(~ metric_row(.x$obs, .x$prob, .x$call)) %>%
  ungroup()

# ----------------------------------------------------------------------------
# 8. Summarize repeat-level performance and patient-clustered uncertainty
# ----------------------------------------------------------------------------
# Metrics are pooled within outer repeat. Bootstrap resampling uses patients,
# not rows, so all longitudinal samples from a patient remain a dependent unit.
bootstrap_results <- list()
for (model_name in names(SPECS)) {
  message("Patient-clustered bootstrap for ", model_name, " (", BOOTSTRAP_REPS,
          " replicates).")
  bootstrap_results[[model_name]] <- cluster_bootstrap_metrics(
    predictions %>% filter(model == model_name),
    reps = BOOTSTRAP_REPS,
    seed = BASE_SEED + match(model_name, names(SPECS)) * 5000003L
  ) %>% mutate(model = model_name, .before = 1)
}
bootstrap_metrics <- bind_rows(bootstrap_results)

metric_names <- c("auc", "sensitivity", "specificity", "balanced_accuracy",
                  "accuracy", "brier")
primary_summary <- repeat_metrics %>%
  group_by(model) %>%
  summarise(
    across(all_of(metric_names),
           list(estimate = ~ mean(.x, na.rm = TRUE),
                split_sd = ~ sd(.x, na.rm = TRUE),
                split_q025 = ~ quantile(.x, 0.025, na.rm = TRUE),
                split_q975 = ~ quantile(.x, 0.975, na.rm = TRUE))),
    .groups = "drop"
  )

bootstrap_ci <- bootstrap_metrics %>%
  group_by(model) %>%
  summarise(
    across(all_of(metric_names),
           list(cluster_boot_q025 = ~ quantile(.x, 0.025, na.rm = TRUE),
                cluster_boot_q975 = ~ quantile(.x, 0.975, na.rm = TRUE))),
    .groups = "drop"
  )

primary_summary <- primary_summary %>%
  left_join(bootstrap_ci, by = "model") %>%
  left_join(
    frame_manifest %>%
      group_by(model) %>%
      summarise(
        n_samples = n(), n_patients = n_distinct(Patient),
        n_positive_samples = sum(MRD_truth == "pos"),
        n_negative_samples = sum(MRD_truth == "neg"),
        n_patients_with_repeated_samples = n_distinct(Patient[repeated_patient]),
        .groups = "drop"
      ),
    by = "model"
  ) %>%
  mutate(
    outer_repeats = OUTER_REPEATS,
    outer_folds = 5L,
    inner_repeats = INNER_REPEATS,
    inner_folds = 5L,
    bootstrap_reps = BOOTSTRAP_REPS,
    estimand = "sample-level MRD classification in a previously unseen patient, conditional on the frozen feature specification"
  )

fold_balance <- predictions %>%
  group_by(model, outer_repeat, outer_fold) %>%
  summarise(
    n_samples = n(), n_patients = n_distinct(Patient),
    n_positive = sum(obs == "pos"), n_negative = sum(obs == "neg"),
    .groups = "drop"
  )

# ----------------------------------------------------------------------------
# 9. Export auditable results and enforce completion QC
# ----------------------------------------------------------------------------
# RUN_COMPLETE is written only after the fold, prediction, frame-size, and
# finite-metric invariants below pass. Assemblers reject runs without it.
write_csv(frame_manifest, file.path(OUT_DIR, "training_frame_manifest.csv"))
write_csv(outer_assignments, file.path(OUT_DIR, "outer_patient_fold_assignments.csv"))
write_csv(fold_balance, file.path(OUT_DIR, "outer_fold_balance.csv"))
write_csv(predictions, file.path(OUT_DIR, "outer_heldout_predictions.csv"))
write_csv(fold_metrics, file.path(OUT_DIR, "outer_fold_metrics.csv"))
write_csv(repeat_metrics, file.path(OUT_DIR, "outer_repeat_pooled_metrics.csv"))
write_csv(inner_tuning, file.path(OUT_DIR, "inner_tuning_and_thresholds.csv"))
write_csv(inner_assignments, file.path(OUT_DIR, "inner_patient_fold_assignments.csv"))
write_csv(training_warnings, file.path(OUT_DIR, "model_fitting_warnings.csv"))
write_csv(bootstrap_metrics, file.path(OUT_DIR, "patient_cluster_bootstrap_metrics.csv"))
write_csv(primary_summary, file.path(OUT_DIR, "primary_performance_summary.csv"))

run_parameters <- tibble(
  parameter = c("run_id", "output_root", "aggregate_input", "outer_repeats", "outer_folds",
                "inner_repeats", "inner_folds", "bootstrap_reps", "base_seed",
                "fold_assignment_attempts", "frozen_feature_selection", "models"),
  value = c(RUN_ID, OUT_ROOT, AGG_PATH, OUTER_REPEATS, 5L, INNER_REPEATS, 5L,
            BOOTSTRAP_REPS, BASE_SEED, 500L, "TRUE", paste(names(SPECS), collapse = ","))
)
write_csv(run_parameters, file.path(OUT_DIR, "run_parameters.csv"))

capture.output(sessionInfo(), file = file.path(OUT_DIR, "session_info.txt"))

qc <- tibble(
  check = c(
    "all_outer_splits_have_zero_patient_leakage",
    "each_sample_predicted_once_per_outer_repeat",
    "all_outer_folds_have_both_classes",
    "each_inner_patient_has_one_fold_per_inner_repeat",
    "frozen_artifacts_remained_read_only",
    "legacy_diagnostic_files_preserved"
  ),
  status = c(
    ifelse(all(fold_metrics$leaked_patients == 0L), "PASS", "FAIL"),
    ifelse(all(prediction_check$expected == prediction_check$observed), "PASS", "FAIL"),
    ifelse(all(fold_balance$n_positive > 0 & fold_balance$n_negative > 0), "PASS", "FAIL"),
    ifelse(
      all(inner_assignments %>%
            count(model, outer_repeat, outer_fold, inner_repeat, Patient) %>%
            pull(n) == 1L),
      "PASS", "FAIL"
    ),
    ifelse(all(file.access(present_frozen, mode = 2) != 0), "PASS", "FAIL"),
    ifelse(all(file.exists(legacy_diagnostic_files)), "PASS", "FAIL")
  )
)
write_csv(qc, file.path(OUT_DIR, "quality_control_checks.csv"))
if (any(qc$status != "PASS")) stop("One or more final QC checks failed.", call. = FALSE)

writeLines(
  c(
    "Patient-grouped repeated nested cross-validation completed successfully.",
    paste("Run ID:", RUN_ID),
    paste("Output directory:", OUT_DIR),
    "Historical results and frozen model artifacts were not overwritten.",
    "See primary_performance_summary.csv for the primary estimates and patient-clustered intervals."
  ),
  file.path(OUT_DIR, "RUN_COMPLETE")
)

message("\nCompleted successfully. Primary summary:")
print(primary_summary %>%
        select(model, n_samples, n_patients, auc_estimate,
               auc_cluster_boot_q025, auc_cluster_boot_q975,
               sensitivity_estimate, specificity_estimate))
message("Outputs written to: ", OUT_DIR)
