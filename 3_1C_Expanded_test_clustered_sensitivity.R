#!/usr/bin/env Rscript

# Expanded-test repeated-measures sensitivity analysis
# ====================================================
#
# Goal
# ----
# Quantify how repeated samples from the same patient affect the expanded-test
# confusion-matrix estimates used in Main Figures 3B and 4B.
#
# Analyses
# --------
# 1. Patient-clustered nonparametric bootstrap:
#    patients are sampled with replacement and every evaluable sample from a
#    selected patient is retained. This preserves within-patient correlation
#    while keeping the published sample-level estimand.
# 2. Deterministic one-sample-per-patient check:
#    the earliest evaluable post-baseline sample is selected without reference
#    to truth or prediction, with ties broken by numeric timepoint and then
#    Sample_Code.
#
# Inputs
# ------
# - Current masked scored table from 3_1_Optimize_cfWGS_thresholds.R
# - Current Figure 3B and Figure 4B test-confusion source tables, which provide
#   the frozen model thresholds and manuscript-facing model labels
#
# Outputs
# -------
# Stable CSVs and an audit text file are written to:
#   Output_tables_2025/expanded_test_clustered_sensitivity/
#
# Reproducibility
# ---------------
# The bootstrap and within-patient sampling use fixed, model-specific seeds.
# The script fails on cohort-size, truth-label, threshold, probability-range,
# or identifier drift rather than silently changing the analysis population.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve the current script path from commandArgs().", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

input_scored <- file.path(
  project_root,
  "Output_tables_2025",
  "all_patients_with_BM_and_blood_calls_updated6.csv"
)
input_bm_confusion <- file.path(
  script_dir,
  "final_manuscript_objects",
  "01_main_figures",
  "Figure_3",
  "Figure_3B",
  "F3B_supporting_data_csv_test_confusion.csv"
)
input_blood_confusion <- file.path(
  script_dir,
  "final_manuscript_objects",
  "01_main_figures",
  "Figure_4",
  "Figure_4B",
  "F4B_supporting_data_csv_test_confusion.csv"
)
output_dir <- file.path(
  project_root,
  "Output_tables_2025",
  "expanded_test_clustered_sensitivity"
)

required_files <- c(input_scored, input_bm_confusion, input_blood_confusion)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) are missing:\n- ",
    paste(missing_files, collapse = "\n- "),
    call. = FALSE
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_bootstrap <- 10000L
ci_level <- 0.95
ci_probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
analysis_date <- "2026-07-23"

model_specs <- tribble(
  ~analysis_family, ~model_key, ~probability_column, ~expected_samples, ~expected_patients, ~seed,
  "BM-informed", "BM_zscore_only_detection_rate", "BM_zscore_only_detection_rate_prob", 28L, 19L, 202607231L,
  "BM-informed", "BM_base_zscore", "BM_base_zscore_prob", 28L, 19L, 202607231L,
  "Baseline-plasma-informed", "Blood_zscore_only_sites", "Blood_zscore_only_sites_prob", 24L, 15L, 202607232L,
  "Baseline-plasma-informed", "Blood_plus_fragment", "Blood_plus_fragment_prob", 24L, 15L, 202607232L
)

scored <- read_csv(input_scored, show_col_types = FALSE)
bm_thresholds <- read_csv(input_bm_confusion, show_col_types = FALSE)
blood_thresholds <- read_csv(input_blood_confusion, show_col_types = FALSE)
threshold_table <- bind_rows(bm_thresholds, blood_thresholds) %>%
  distinct(.data$model_key, .data$model, .data$threshold)

required_scored_columns <- unique(c(
  "Patient", "Timepoint", "Sample_Code", "Date", "timepoint_info", "Cohort",
  "Study", "MRD_truth", model_specs$probability_column
))
missing_columns <- setdiff(required_scored_columns, names(scored))
if (length(missing_columns) > 0L) {
  stop(
    "The scored table is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (!inherits(scored$Date, "Date")) {
  scored <- scored %>% mutate(Date = as.Date(.data$Date))
}

threshold_counts <- threshold_table %>% count(.data$model_key, name = "n_threshold_rows")
invalid_threshold_keys <- threshold_counts %>%
  filter(.data$n_threshold_rows != 1L) %>%
  pull(.data$model_key)
if (length(invalid_threshold_keys) > 0L) {
  stop(
    "Each model must have exactly one frozen threshold. Invalid model(s): ",
    paste(invalid_threshold_keys, collapse = ", "),
    call. = FALSE
  )
}

model_specs <- model_specs %>%
  left_join(threshold_table, by = "model_key")
if (anyNA(model_specs$threshold) || anyNA(model_specs$model)) {
  stop("Frozen threshold or manuscript model label could not be resolved.", call. = FALSE)
}

# Both models within a specimen family must use the same evaluable rows. This
# preserves paired comparisons and ensures the same patient-cluster draw is
# applied to both family models.
walk(unique(model_specs$analysis_family), function(family_name) {
  family_specs <- model_specs %>% filter(.data$analysis_family == family_name)
  family_probability_columns <- family_specs$probability_column
  family_rows <- scored %>%
    filter(
      .data$Cohort == "Non-frontline",
      !.data$timepoint_info %in% c("Baseline", "Diagnosis"),
      !is.na(.data$MRD_truth),
      if_all(all_of(family_probability_columns), ~ !is.na(.x))
    ) %>%
    pull(.data$Sample_Code)

  walk(family_probability_columns, function(probability_column) {
    model_rows <- scored %>%
      filter(
        .data$Cohort == "Non-frontline",
        !.data$timepoint_info %in% c("Baseline", "Diagnosis"),
        !is.na(.data$MRD_truth),
        !is.na(.data[[probability_column]])
      ) %>%
      pull(.data$Sample_Code)
    if (!setequal(family_rows, model_rows)) {
      stop(
        "Models within ", family_name,
        " do not have identical evaluable Sample_Code sets.",
        call. = FALSE
      )
    }
  })
})

safe_divide <- function(numerator, denominator) {
  if (is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

auc_rank <- function(truth, probability) {
  positive <- probability[truth == 1L]
  negative <- probability[truth == 0L]
  if (length(positive) == 0L || length(negative) == 0L) {
    return(NA_real_)
  }
  comparisons <- outer(positive, negative, FUN = "-")
  mean(comparisons > 0) + 0.5 * mean(comparisons == 0)
}

calculate_metrics <- function(truth, probability, threshold) {
  prediction <- as.integer(probability >= threshold)
  tp <- sum(prediction == 1L & truth == 1L)
  fn <- sum(prediction == 0L & truth == 1L)
  tn <- sum(prediction == 0L & truth == 0L)
  fp <- sum(prediction == 1L & truth == 0L)
  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  ppv <- safe_divide(tp, tp + fp)
  npv <- safe_divide(tn, tn + fn)

  c(
    n_samples = length(truth),
    n_positive = sum(truth == 1L),
    n_negative = sum(truth == 0L),
    tp = tp,
    fn = fn,
    tn = tn,
    fp = fp,
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = safe_divide(tp + tn, tp + tn + fp + fn),
    balanced_accuracy = if (is.na(sensitivity) || is.na(specificity)) {
      NA_real_
    } else {
      mean(c(sensitivity, specificity))
    },
    ppv = ppv,
    npv = npv,
    auc = auc_rank(truth, probability)
  )
}

prepare_model_data <- function(spec_row) {
  probability_column <- spec_row$probability_column[[1]]
  threshold <- spec_row$threshold[[1]]

  model_data <- scored %>%
    filter(
      .data$Cohort == "Non-frontline",
      !.data$timepoint_info %in% c("Baseline", "Diagnosis"),
      !is.na(.data$MRD_truth),
      !is.na(.data[[probability_column]])
    ) %>%
    transmute(
      analysis_family = spec_row$analysis_family[[1]],
      model_key = spec_row$model_key[[1]],
      model = spec_row$model[[1]],
      threshold = threshold,
      Patient = as.character(.data$Patient),
      Timepoint = as.character(.data$Timepoint),
      timepoint_order = readr::parse_number(as.character(.data$Timepoint)),
      Sample_Code = as.character(.data$Sample_Code),
      Date = as.Date(.data$Date),
      timepoint_info = as.character(.data$timepoint_info),
      Study = as.character(.data$Study),
      truth = as.integer(.data$MRD_truth),
      probability = as.numeric(.data[[probability_column]]),
      prediction = as.integer(.data[[probability_column]] >= threshold)
    ) %>%
    arrange(.data$Patient, .data$Date, .data$Sample_Code)

  if (anyNA(model_data$Patient) || any(model_data$Patient == "")) {
    stop("Patient identifiers are missing for ", spec_row$model_key[[1]], ".", call. = FALSE)
  }
  if (any(!model_data$truth %in% c(0L, 1L))) {
    stop("MRD_truth must contain only 0/1 for ", spec_row$model_key[[1]], ".", call. = FALSE)
  }
  if (any(model_data$probability < 0 | model_data$probability > 1)) {
    stop("Model probabilities fall outside [0,1] for ", spec_row$model_key[[1]], ".", call. = FALSE)
  }
  if (anyNA(model_data$Date)) {
    stop(
      "Date is required for the earliest-visit sensitivity analysis: ",
      spec_row$model_key[[1]], ".",
      call. = FALSE
    )
  }
  duplicate_sample_codes <- model_data %>%
    count(.data$Sample_Code) %>%
    filter(.data$n != 1L)
  if (nrow(duplicate_sample_codes) > 0L) {
    stop(
      "Sample_Code must be unique within each model analysis: ",
      spec_row$model_key[[1]], ".",
      call. = FALSE
    )
  }

  observed_samples <- nrow(model_data)
  observed_patients <- n_distinct(model_data$Patient)
  if (
    observed_samples != spec_row$expected_samples[[1]] ||
    observed_patients != spec_row$expected_patients[[1]]
  ) {
    stop(
      "Expanded-test cohort drift for ", spec_row$model_key[[1]], ": expected ",
      spec_row$expected_samples[[1]], " samples from ",
      spec_row$expected_patients[[1]], " patients; observed ",
      observed_samples, " samples from ", observed_patients, " patients.",
      call. = FALSE
    )
  }

  model_data
}

select_earliest_per_patient <- function(model_data) {
  model_data %>%
    arrange(
      .data$Patient,
      .data$Date,
      is.na(.data$timepoint_order),
      .data$timepoint_order,
      .data$Sample_Code
    ) %>%
    group_by(.data$Patient) %>%
    slice(1L) %>%
    ungroup()
}

bootstrap_by_patient <- function(model_data, threshold, n_replicates, seed) {
  patient_rows <- split(seq_len(nrow(model_data)), model_data$Patient)
  patient_ids <- names(patient_rows)
  set.seed(seed)
  replicate(
    n_replicates,
    {
      sampled_patients <- sample(patient_ids, length(patient_ids), replace = TRUE)
      sampled_rows <- unlist(patient_rows[sampled_patients], use.names = FALSE)
      calculate_metrics(
        model_data$truth[sampled_rows],
        model_data$probability[sampled_rows],
        threshold
      )
    }
  ) %>%
    t() %>%
    as_tibble()
}

metric_names <- c(
  "sensitivity", "specificity", "accuracy", "balanced_accuracy", "ppv", "npv", "auc"
)

analysis_results <- pmap(
  model_specs,
  function(
      analysis_family, model_key, probability_column, expected_samples,
      expected_patients, seed, model, threshold) {
    spec_row <- tibble(
      analysis_family = analysis_family,
      model_key = model_key,
      probability_column = probability_column,
      expected_samples = expected_samples,
      expected_patients = expected_patients,
      seed = seed,
      model = model,
      threshold = threshold
    )
    model_data <- prepare_model_data(spec_row)
    earliest <- select_earliest_per_patient(model_data)

    full_metrics <- calculate_metrics(
      model_data$truth, model_data$probability, threshold
    )
    earliest_metrics <- calculate_metrics(
      earliest$truth, earliest$probability, threshold
    )
    bootstrap_distribution <- bootstrap_by_patient(
      model_data, threshold, n_bootstrap, seed
    )
    bootstrap_summary <- map_dfr(metric_names, function(metric_name) {
      values <- bootstrap_distribution[[metric_name]]
      valid_values <- values[is.finite(values)]
      tibble(
        analysis_family = analysis_family,
        model_key = model_key,
        model = model,
        threshold = threshold,
        metric = metric_name,
        estimate = unname(full_metrics[[metric_name]]),
        ci_level = ci_level,
        ci_lower = if (length(valid_values)) {
          unname(quantile(valid_values, ci_probs[[1]], names = FALSE, type = 6))
        } else {
          NA_real_
        },
        ci_upper = if (length(valid_values)) {
          unname(quantile(valid_values, ci_probs[[2]], names = FALSE, type = 6))
        } else {
          NA_real_
        },
        valid_replicates = length(valid_values),
        total_replicates = n_bootstrap,
        resampling_unit = "Patient; all evaluable samples retained within resampled patient"
      )
    })

    one_sample_summary <- bind_rows(
      as_tibble_row(full_metrics) %>% mutate(selection = "All evaluable samples"),
      as_tibble_row(earliest_metrics) %>% mutate(selection = "Earliest evaluable sample per patient")
    ) %>%
      mutate(
        analysis_family = analysis_family,
        model_key = model_key,
        model = model,
        threshold = threshold,
        .before = 1
      )

    manifest <- model_data %>%
      mutate(
        selected_earliest = .data$Sample_Code %in% earliest$Sample_Code,
        samples_for_patient = ave(
          rep(1L, n()),
          .data$Patient,
          FUN = length
        )
      )

    list(
      bootstrap_summary = bootstrap_summary,
      one_sample_summary = one_sample_summary,
      manifest = manifest
    )
  }
)

bootstrap_summary <- map_dfr(analysis_results, "bootstrap_summary")
one_sample_summary <- map_dfr(analysis_results, "one_sample_summary")
sample_manifest <- map_dfr(analysis_results, "manifest")

round_numeric <- function(data, digits = 4L) {
  data %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
}

bootstrap_output <- round_numeric(bootstrap_summary)
one_sample_output <- round_numeric(one_sample_summary)
write_csv(
  bootstrap_output,
  file.path(output_dir, paste0(
    "expanded_test_patient_clustered_bootstrap_metrics_", analysis_date, ".csv"
  ))
)
write_csv(
  one_sample_output,
  file.path(output_dir, paste0(
    "expanded_test_one_sample_per_patient_metrics_", analysis_date, ".csv"
  ))
)
write_csv(
  sample_manifest,
  file.path(output_dir, paste0(
    "expanded_test_repeated_measures_sample_manifest_", analysis_date, ".csv"
  ))
)

audit_lines <- c(
  "Expanded-test repeated-measures sensitivity analysis",
  paste0("Analysis date: ", analysis_date),
  paste0("Script: ", script_path),
  paste0("Scored input: ", input_scored),
  paste0("BM threshold input: ", input_bm_confusion),
  paste0("Blood threshold input: ", input_blood_confusion),
  paste0("Patient-clustered bootstrap replicates: ", n_bootstrap),
  paste0("Confidence level: ", ci_level),
  "",
  "Cohort checks:",
  paste0(
    "- ", model_specs$analysis_family, " / ", model_specs$model,
    ": ", model_specs$expected_samples, " samples from ",
    model_specs$expected_patients, " patients; threshold ",
    format(model_specs$threshold, digits = 8)
  ),
  "",
  "Interpretation:",
  "- Cluster-bootstrap intervals account for patient-level clustering while preserving the sample-level estimand.",
  "- The earliest-visit check changes the estimand to one selected sample per patient."
)
write_lines(
  audit_lines,
  file.path(output_dir, paste0(
    "expanded_test_repeated_measures_analysis_audit_", analysis_date, ".txt"
  ))
)

message("Wrote repeated-measures sensitivity outputs to: ", output_dir)
