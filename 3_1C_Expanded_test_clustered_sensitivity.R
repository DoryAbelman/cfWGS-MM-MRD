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
# Stable CSVs, an audit text file, and the final Supplementary Table 6 workbook
# are written to:
#   Output_tables_2025/expanded_test_clustered_sensitivity/
# The workbook is also copied to the final manuscript-object tree as STABLE6.
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
  library(openxlsx)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve the current script path from commandArgs().", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

helper_path <- file.path(script_dir, "manuscript_output_helpers.R")
if (!file.exists(helper_path)) {
  stop("Missing manuscript output helper: ", helper_path, call. = FALSE)
}
source(helper_path)

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
input_thresholds <- file.path(
  project_root,
  "Output_tables_2025",
  "selected_combo_thresholds_2026-02-16.rds"
)
output_dir <- file.path(
  project_root,
  "Output_tables_2025",
  "expanded_test_clustered_sensitivity"
)

required_files <- c(
  input_scored, input_bm_confusion, input_blood_confusion, input_thresholds
)
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
  "Baseline-plasma-informed", "Blood_zscore_only_sites", "Blood_zscore_only_sites_prob", 25L, 16L, 202607232L,
  "Baseline-plasma-informed", "Blood_plus_fragment", "Blood_plus_fragment_prob", 25L, 16L, 202607232L
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

# Same-sample comparison requested during peer review. Every model below is
# evaluated on the identical 28-sample/19-patient BM-informed expanded-test
# stratum. Fragmentomics probabilities are evaluated at the fixed thresholds
# selected by the original training workflow; no test-set refitting or
# threshold selection is performed here.
matched_model_specs <- tribble(
  ~workflow, ~model_key, ~probability_column,
  "BM-informed mutation tracking", "BM_zscore_only_detection_rate", "BM_zscore_only_detection_rate_prob",
  "Fragmentomics only; full-cohort training", "Fragmentomics_full_Full", "Fragmentomics_full_Full_prob",
  "Fragmentomics only; full-cohort training", "Fragmentomics_mean_coverage_only_Full", "Fragmentomics_mean_coverage_only_Full_prob",
  "Fragmentomics only; BM-restricted training", "Fragmentomics_full_BM_restricted", "Fragmentomics_full_BM_restricted_prob",
  "Fragmentomics only; BM-restricted training", "Fragmentomics_mean_coverage_only_BM_restricted", "Fragmentomics_mean_coverage_only_BM_restricted_prob"
)

fixed_thresholds <- readRDS(input_thresholds)
if (!is.numeric(fixed_thresholds) || is.null(names(fixed_thresholds))) {
  stop("Frozen-threshold input must be a named numeric vector.", call. = FALSE)
}
missing_matched_columns <- setdiff(
  matched_model_specs$probability_column,
  names(scored)
)
if (length(missing_matched_columns) > 0L) {
  stop(
    "Scored table is missing matched-comparison probability column(s): ",
    paste(missing_matched_columns, collapse = ", "),
    call. = FALSE
  )
}
missing_matched_thresholds <- setdiff(
  matched_model_specs$model_key,
  names(fixed_thresholds)
)
if (length(missing_matched_thresholds) > 0L) {
  stop(
    "Frozen threshold(s) missing for matched comparison: ",
    paste(missing_matched_thresholds, collapse = ", "),
    call. = FALSE
  )
}

matched_rows <- scored %>%
  filter(
    .data$Cohort == "Non-frontline",
    !.data$timepoint_info %in% c("Baseline", "Diagnosis"),
    !is.na(.data$MRD_truth),
    if_all(all_of(matched_model_specs$probability_column), ~ !is.na(.x))
  )

if (nrow(matched_rows) != 28L || n_distinct(matched_rows$Patient) != 19L) {
  stop(
    "Matched workflow comparison cohort drift: expected 28 samples from 19 patients; observed ",
    nrow(matched_rows), " samples from ", n_distinct(matched_rows$Patient), " patients.",
    call. = FALSE
  )
}

matched_workflow_comparison <- pmap_dfr(
  matched_model_specs,
  function(workflow, model_key, probability_column) {
    threshold <- unname(fixed_thresholds[[model_key]])
    metrics <- calculate_metrics(
      truth = as.integer(matched_rows$MRD_truth),
      probability = as.numeric(matched_rows[[probability_column]]),
      threshold = threshold
    )
    as_tibble_row(metrics) %>%
      mutate(
        analysis_set = "BM-informed expanded-test evaluable stratum",
        workflow = workflow,
        model_key = model_key,
        probability_column = probability_column,
        threshold = threshold,
        n_patients = n_distinct(matched_rows$Patient),
        .before = 1
      )
  }
)

round_numeric <- function(data, digits = 4L) {
  data %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
}

# openxlsx can leave empty drawing relationships and an A1-only worksheet
# dimension in otherwise valid workbooks. Excel usually repairs these silently,
# but strict OOXML readers do not. Remove only relationships whose drawing
# targets are absent and set each worksheet dimension from its actual cell refs.
sanitize_openxlsx_package <- function(path) {
  staging <- tempfile("st6_ooxml_")
  dir.create(staging, recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, exdir = staging)

  read_xml_text <- function(xml_path) {
    paste(readLines(xml_path, warn = FALSE, encoding = "UTF-8"), collapse = "")
  }
  write_xml_text <- function(xml_path, value) {
    writeLines(value, xml_path, useBytes = TRUE)
  }
  column_number <- function(label) {
    values <- utf8ToInt(label) - utf8ToInt("A") + 1L
    as.integer(sum(values * 26L ^ rev(seq_along(values) - 1L)))
  }
  column_label <- function(number) {
    output <- character()
    while (number > 0L) {
      number <- number - 1L
      output <- c(intToUtf8(number %% 26L + utf8ToInt("A")), output)
      number <- number %/% 26L
    }
    paste(output, collapse = "")
  }

  rel_paths <- list.files(
    file.path(staging, "xl", "worksheets", "_rels"),
    pattern = "[.]rels$", full.names = TRUE
  )
  walk(rel_paths, function(rel_path) {
    xml <- read_xml_text(rel_path)
    xml <- gsub(
      '<Relationship[^>]+Type="[^"]+/relationships/(drawing|vmlDrawing)"[^>]*/>',
      "", xml, perl = TRUE
    )
    write_xml_text(rel_path, xml)
  })

  sheet_paths <- list.files(
    file.path(staging, "xl", "worksheets"),
    pattern = "^sheet[0-9]+[.]xml$", full.names = TRUE
  )
  walk(sheet_paths, function(sheet_path) {
    xml <- read_xml_text(sheet_path)
    refs <- regmatches(
      xml,
      gregexpr('(?<= r=")[A-Z]+[0-9]+(?=")', xml, perl = TRUE)
    )[[1]]
    refs <- refs[nzchar(refs)]
    if (length(refs) == 0L) return(invisible(NULL))
    rows <- as.integer(gsub("[A-Z]+", "", refs))
    cols <- vapply(
      gsub("[0-9]+", "", refs), column_number, integer(1)
    )
    dimension <- paste0("A1:", column_label(max(cols)), max(rows))
    xml <- sub(
      '<dimension ref="[^"]*"/>',
      paste0('<dimension ref="', dimension, '"/>'),
      xml
    )
    write_xml_text(sheet_path, xml)
  })

  rebuilt <- tempfile(
    "Supplementary_Table_6_sanitized_", tmpdir = dirname(path), fileext = ".xlsx"
  )
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(staging)
  archive_files <- list.files(
    ".", recursive = TRUE, all.files = TRUE, no.. = TRUE,
    include.dirs = FALSE
  )
  utils::zip(rebuilt, files = archive_files, flags = "-q -r -X")
  setwd(old_wd)
  if (!file.exists(rebuilt) || file.info(rebuilt)$size <= 0L) {
    stop("Failed to rebuild sanitized Supplementary Table 6 OOXML package.", call. = FALSE)
  }
  if (!file.copy(rebuilt, path, overwrite = TRUE)) {
    stop("Failed to replace Supplementary Table 6 with sanitized package.", call. = FALSE)
  }
  unlink(rebuilt)
  invisible(path)
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
write_csv(
  round_numeric(matched_workflow_comparison),
  file.path(output_dir, paste0(
    "expanded_test_matched_workflow_comparison_", analysis_date, ".csv"
  ))
)

audit_lines <- c(
  "Expanded-test repeated-measures sensitivity analysis",
  paste0("Analysis date: ", analysis_date),
  paste0("Script: ", script_path),
  paste0("Scored input: ", input_scored),
  paste0("BM threshold input: ", input_bm_confusion),
  paste0("Blood threshold input: ", input_blood_confusion),
  paste0("Frozen threshold vector: ", input_thresholds),
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
  ,"- The matched workflow comparison evaluates mutation-informed and fragmentomics-only models on the identical 28-sample/19-patient BM-informed stratum."
)
write_lines(
  audit_lines,
  file.path(output_dir, paste0(
    "expanded_test_repeated_measures_analysis_audit_", analysis_date, ".txt"
  ))
)

message("Wrote repeated-measures sensitivity outputs to: ", output_dir)

# -------------------------------------------------------------------------
# Build final revision-inclusive Supplementary Table 6.
#
# The original classifier-performance rows remain owned by
# 3_1_Optimize_cfWGS_thresholds.R. This block combines that preserved-model
# export with the repeated-measures results generated above. It does not
# retrain models or reselect thresholds.
# -------------------------------------------------------------------------

performance_csv <- file.path(
  project_root,
  "Final Tables and Figures",
  "Supplementary_Table_5_All_Model_performance_testing_cohort_v3_Feb2026_with_restricted_fragmentomics.csv"
)
if (!file.exists(performance_csv)) {
  stop("Missing expanded-test classifier-performance source: ", performance_csv, call. = FALSE)
}

performance <- read_csv(performance_csv, show_col_types = FALSE)
expected_performance_columns <- c(
  "combo", "auc_mean", "sens_mean", "spec_mean", "sens95_mean",
  "spec95_mean", "balacc_mean", "ppv_mean", "npv_mean", "f1_mean",
  "brier_mean", "pAUC90_mean", "acc_mean", "sample_group", "eval_cohort"
)
missing_performance_columns <- setdiff(expected_performance_columns, names(performance))
if (length(missing_performance_columns) > 0L) {
  stop(
    "Classifier-performance input is missing column(s): ",
    paste(missing_performance_columns, collapse = ", "),
    call. = FALSE
  )
}
if (nrow(performance) != 32L) {
  stop(
    "Expected 32 revision-inclusive classifier-performance rows; observed ",
    nrow(performance), ".",
    call. = FALSE
  )
}

expected_primary <- tribble(
  ~model_key, ~sensitivity, ~specificity, ~accuracy,
  "BM_zscore_only_detection_rate", 0.6667, 0.6842, 0.6786,
  "BM_base_zscore", 0.6667, 0.3684, 0.4643,
  "Blood_zscore_only_sites", 0.2857, 0.7778, 0.6400,
  "Blood_plus_fragment", 0.2857, 0.8889, 0.7200
)
metric_tolerance <- 0.0011
walk(seq_len(nrow(expected_primary)), function(i) {
  expected_row <- expected_primary[i, ]
  observed_row <- performance %>% filter(.data$combo == expected_row$model_key)
  if (nrow(observed_row) != 1L) {
    stop(
      "Expected exactly one classifier-performance row for ",
      expected_row$model_key, ".",
      call. = FALSE
    )
  }
  observed <- c(
    sensitivity = observed_row$sens_mean[[1]],
    specificity = observed_row$spec_mean[[1]],
    accuracy = observed_row$acc_mean[[1]]
  )
  expected <- c(
    sensitivity = expected_row$sensitivity[[1]],
    specificity = expected_row$specificity[[1]],
    accuracy = expected_row$accuracy[[1]]
  )
  if (any(abs(observed - expected) > metric_tolerance)) {
    stop(
      "Primary expanded-test metrics drifted for ", expected_row$model_key,
      ". Observed: ", paste(names(observed), round(observed, 4), collapse = "; "),
      ". Expected: ", paste(names(expected), round(expected, 4), collapse = "; "),
      ".",
      call. = FALSE
    )
  }
})

if (nrow(bootstrap_output) != 28L) {
  stop("Expected 28 clustered-bootstrap rows; observed ", nrow(bootstrap_output), ".", call. = FALSE)
}
if (!all(bootstrap_output$total_replicates == n_bootstrap)) {
  stop("All clustered-bootstrap rows must report 10,000 total replicates.", call. = FALSE)
}
if (nrow(one_sample_output) != 8L) {
  stop("Expected 8 one-sample-per-patient rows; observed ", nrow(one_sample_output), ".", call. = FALSE)
}

expected_family_counts <- tribble(
  ~analysis_family, ~all_samples, ~patient_samples,
  "BM-informed", 28L, 19L,
  "Baseline-plasma-informed", 25L, 16L
)
walk(seq_len(nrow(expected_family_counts)), function(i) {
  family_row <- expected_family_counts[i, ]
  family_results <- one_sample_output %>%
    filter(.data$analysis_family == family_row$analysis_family)
  all_rows <- family_results %>% filter(.data$selection == "All evaluable samples")
  patient_rows <- family_results %>%
    filter(.data$selection == "Earliest evaluable sample per patient")
  if (
    nrow(all_rows) != 2L ||
      nrow(patient_rows) != 2L ||
      any(all_rows$n_samples != family_row$all_samples) ||
      any(patient_rows$n_samples != family_row$patient_samples)
  ) {
    stop(
      "Revision-inclusive denominator check failed for ",
      family_row$analysis_family, ".",
      call. = FALSE
    )
  }
})

workbook <- createWorkbook(creator = "Abelman et al.")
sheet_names <- c(
  "README",
  "Classifier performance",
  "Clustered bootstrap",
  "One sample per patient",
  "Matched workflow comparison",
  "Sample manifest"
)
walk(sheet_names, ~ addWorksheet(workbook, .x, gridLines = FALSE))

header_style <- createStyle(
  fontName = "Arial", fontSize = 9, fontColour = "#FFFFFF",
  fgFill = "#335C67", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
body_style <- createStyle(
  fontName = "Arial", fontSize = 9, fontColour = "#222222", valign = "center"
)
label_style <- createStyle(
  fontName = "Arial", fontSize = 9, fontColour = "#22333B",
  fgFill = "#EAF0F2", textDecoration = "bold", valign = "top"
)
title_style <- createStyle(
  fontName = "Arial", fontSize = 15, fontColour = "#FFFFFF",
  fgFill = "#335C67", textDecoration = "bold", valign = "center"
)
rate_style <- createStyle(numFmt = "0.0000")
count_style <- createStyle(numFmt = "0")

readme <- tibble(
  Section = c(
    "Scope", "Classifier performance", "Clustered bootstrap",
    "One sample per patient", "Matched workflow comparison", "Sample manifest", "BM-informed denominator",
    "Baseline-plasma denominator", "Fixed random seeds", "Canonical analysis scripts"
  ),
  Description = c(
    "Revision-inclusive expanded independent test cohort; frozen training model thresholds.",
    "All model performance metrics in the current expanded test cohort.",
    paste0(
      "Patient-clustered nonparametric bootstrap; 10,000 replicates; all evaluable ",
      "samples from each resampled patient retained; percentile 95% confidence intervals."
    ),
    paste0(
      "Deterministic sensitivity analysis using the earliest evaluable post-baseline ",
      "sample per patient, selected without reference to truth or prediction."
    ),
    paste0(
      "Mutation-informed and fragmentomics-only models evaluated on the identical ",
      "28-sample/19-patient BM-informed expanded-test stratum using frozen thresholds."
    ),
    "Exact evaluable sample/patient rows used by the repeated-measures sensitivity analyses.",
    "28 samples from 19 patients.",
    "25 samples from 16 patients (including revision sample IMG-257-T2).",
    "BM-informed: 202607231; baseline-plasma-informed: 202607232.",
    paste(
      "Scripts_2025/Final_Scripts/3_1_Optimize_cfWGS_thresholds.R",
      "Scripts_2025/Final_Scripts/3_1C_Expanded_test_clustered_sensitivity.R",
      sep = "\n"
    )
  )
)

mergeCells(workbook, "README", cols = 1:2, rows = 1)
writeData(
  workbook, "README",
  "Supplementary Table 6 — Expanded-test classifier performance",
  startCol = 1, startRow = 1
)
addStyle(workbook, "README", title_style, rows = 1, cols = 1:2, gridExpand = TRUE)
writeData(workbook, "README", readme, startRow = 3, headerStyle = header_style)
addStyle(workbook, "README", label_style, rows = 4:(nrow(readme) + 3), cols = 1, gridExpand = TRUE)
addStyle(workbook, "README", body_style, rows = 4:(nrow(readme) + 3), cols = 2, gridExpand = TRUE)
setColWidths(workbook, "README", cols = 1, widths = 28)
setColWidths(workbook, "README", cols = 2, widths = 58)
setRowHeights(workbook, "README", rows = 1, heights = 28)
setRowHeights(workbook, "README", rows = 4:(nrow(readme) + 3), heights = 38)
freezePane(workbook, "README", firstActiveRow = 4)
pageSetup(
  workbook, "README", orientation = "landscape",
  paperSize = 9, fitToWidth = 1, fitToHeight = 1
)

data_sheets <- list(
  "Classifier performance" = performance,
  "Clustered bootstrap" = bootstrap_output,
  "One sample per patient" = one_sample_output,
  "Matched workflow comparison" = round_numeric(matched_workflow_comparison),
  "Sample manifest" = sample_manifest
)
rate_columns <- c(
  "auc_mean", "sens_mean", "spec_mean", "sens95_mean", "spec95_mean",
  "balacc_mean", "ppv_mean", "npv_mean", "f1_mean", "brier_mean",
  "pAUC90_mean", "acc_mean", "threshold", "estimate", "ci_level",
  "ci_lower", "ci_upper", "sensitivity", "specificity", "accuracy",
  "balanced_accuracy", "ppv", "npv", "auc"
)
count_columns <- c(
  "n_samples", "n_positive", "n_negative", "tp", "fn", "tn", "fp",
  "valid_replicates", "total_replicates", "samples_for_patient", "n_patients"
)

walk(names(data_sheets), function(sheet) {
  data <- data_sheets[[sheet]]
  writeData(workbook, sheet, data, startRow = 1, headerStyle = header_style)
  addFilter(workbook, sheet, rows = 1, cols = seq_len(ncol(data)))
  addStyle(
    workbook, sheet, body_style, rows = 2:(nrow(data) + 1),
    cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE
  )
  rate_idx <- which(names(data) %in% rate_columns)
  if (length(rate_idx)) {
    addStyle(
      workbook, sheet, rate_style, rows = 2:(nrow(data) + 1),
      cols = rate_idx, gridExpand = TRUE, stack = TRUE
    )
  }
  count_idx <- which(names(data) %in% count_columns)
  if (length(count_idx)) {
    addStyle(
      workbook, sheet, count_style, rows = 2:(nrow(data) + 1),
      cols = count_idx, gridExpand = TRUE, stack = TRUE
    )
  }
  setColWidths(workbook, sheet, cols = seq_len(ncol(data)), widths = "auto")
  wide_text <- which(vapply(
    data,
    function(x) {
      widths <- nchar(as.character(x))
      widths <- widths[!is.na(widths)]
      length(widths) > 0L && max(widths) > 28L
    },
    logical(1)
  ))
  if (length(wide_text)) {
    setColWidths(workbook, sheet, cols = wide_text, widths = 30)
  }
  freezePane(workbook, sheet, firstActiveRow = 2)
  pageSetup(
    workbook, sheet, orientation = "landscape",
    paperSize = 9, fitToWidth = 1, fitToHeight = 0
  )
})

workbook_path <- file.path(output_dir, "Supplementary_Table_6_FINAL.xlsx")
saveWorkbook(workbook, workbook_path, overwrite = TRUE)
sanitize_openxlsx_package(workbook_path)
if (!file.exists(workbook_path) || file.info(workbook_path)$size <= 0L) {
  stop("Supplementary Table 6 workbook was not written successfully.", call. = FALSE)
}

destination <- ms_copy_artifact(
  source_path = workbook_path,
  artifact_id = "STABLE6",
  role = "submission_workbook_xlsx",
  description = paste0(
    "Revision-inclusive expanded-test classifier performance with patient-clustered ",
    "bootstrap confidence intervals, one-sample-per-patient sensitivity analysis, ",
    "and exact sample manifest."
  ),
  script_name = basename(script_path),
  project_root = project_root,
  overwrite = TRUE
)

selection_readme <- file.path(
  dirname(destination),
  "SUPPLEMENTARY_TABLE_6_SUBMISSION_SELECTION.txt"
)
writeLines(
  c(
    "SUBMISSION FILE",
    basename(destination),
    "",
    "This XLSX is the authoritative revision-inclusive Supplementary Table 6.",
    "Older CSV files in this provenance folder are retained only to preserve analysis history.",
    "Do not submit an older CSV in place of this workbook.",
    "",
    "Rebuild:",
    "Rscript Scripts_2025/Final_Scripts/3_1C_Expanded_test_clustered_sensitivity.R"
  ),
  selection_readme
)

message("Supplementary Table 6 workbook complete: ", destination)
