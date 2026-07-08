################################################################################
## All-Evaluable Longitudinal and Maintenance Landmark KM Analysis
##
## Purpose:
##   Build companion Kaplan-Meier curves that include evaluable patients from
##   both the training/frontline and test/non-frontline cohorts. Because the two
##   cohorts have different sampling schedules, this script uses a patient-level
##   first prospective non-baseline assessment as a broad descriptive landmark,
##   and a first maintenance-phase assessment as the biologically aligned
##   expansion of the original one-year maintenance analysis. The script exports
##   timing/QC tables plus Cox models adjusted for cohort and months from
##   baseline.
##
## Scientific estimand:
##   Among patients who are progression-free and under observation at a selected
##   post-baseline, non-progression-labeled MRD assessment, describe subsequent
##   PFS by MRD status at that assessment.
##
## Key safeguards:
##   - One row per patient per assay/anchor.
##   - Excludes baseline/diagnosis/germline-normal samples from the landmark.
##   - Excludes relapse/progression-labeled samples from the primary anchor to
##     reduce outcome-triggered sampling leakage.
##   - Requires non-negative follow-up from assessment to event/censor date.
##   - Exports event-free-30-day sensitivity anchors to identify same-episode
##     progression concerns.
##   - Exports fixed-window availability summaries for 6 and 12 months, but
##     does not treat them as primary because counts are smaller.
##
## Outputs:
##   Scripts_2025/Final_Scripts/final_manuscript_objects/
##     01_main_figures/ or 02_extended_data_figures/<normal panel>/
##       <panel>__all_evaluable_first_nonbaseline*_figure_panel_png__*.png
##     additional_all_evaluable_first_nonbaseline_km/
##       KM_*_first_nonbaseline_primary.png
##       source_data/*_source_data.csv
##       all_evaluable_first_nonbaseline_anchor_counts.csv
##       all_evaluable_first_nonbaseline_cox_models.csv
##       all_evaluable_first_nonbaseline_delayed_entry_cox_models.csv
##       all_evaluable_time_varying_cox_models.csv
##       all_evaluable_model_comparison_summary.csv
##       all_evaluable_endpoint_denominator_audit.csv
##       all_evaluable_patient_sampling_pattern_summary.csv
##       all_evaluable_first_nonbaseline_method_note.md
##
## Manuscript outputs created/updated:
##   - Additional all-evaluable first-nonbaseline KM review outputs under
##     final_manuscript_objects/additional_all_evaluable_first_nonbaseline_km/.
##
## How to run:
##   Rscript Scripts_2025/Final_Scripts/4_1B_Build_all_evaluable_first_nonbaseline_KM.R
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(lubridate)
  library(purrr)
  library(readr)
  library(rlang)
  library(scales)
  library(stringr)
  library(survival)
  library(survminer)
  library(tibble)
  library(tidyr)
})

.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)

.endpoint_helper <- file.path("Scripts_2025", "Final_Scripts", "next_event_endpoint_helpers.R")
if (!file.exists(.endpoint_helper)) {
  .endpoint_helper <- "next_event_endpoint_helpers.R"
}
source(.endpoint_helper)
rm(.endpoint_helper)

input_calls_rds <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6.rds"
input_pfs_rds <- "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"
input_relapse_dates_rds <- "Exported_data_tables_clinical/Relapse_dates_full_updated.rds"
input_followup_rds <- "Exported_data_tables_clinical/patient_followup_dates_updated.rds"
input_followup_csv <- "Exported_data_tables_clinical/patient_followup_dates_updated.csv"
input_latest_dates_csv <- "Exported_data_tables_clinical/latest_dates_per_patient.csv"
input_latest_dates_updated_csv <- "Exported_data_tables_clinical/latest_dates_per_patient_updated.csv"

output_dir <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects",
  "additional_all_evaluable_first_nonbaseline_km"
)
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

date_tag <- format(Sys.Date(), "%Y-%m-%d")
cohorts_to_include <- c("Frontline", "Non-frontline")
min_group_n_for_km <- 3L
min_events_for_cox <- 5L
days_per_month <- 30.44

required_columns <- function(data, cols, data_name) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      data_name, " is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

read_pfs_table <- function(path) {
  if (!file.exists(path)) {
    stop("Missing PFS/censor input: ", path, call. = FALSE)
  }

  readRDS(path) %>%
    rename_with(
      tolower,
      any_of(c("Baseline_Date", "Censor_date", "Relapsed"))
    ) %>%
    transmute(
      Patient = as.character(Patient),
      baseline_date = as.Date(baseline_date),
      censor_date = as.Date(censor_date),
      relapsed = as.integer(relapsed)
    )
}

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

format_model_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

fit_adjusted_cox <- function(df, assay_label, anchor_label) {
  n_events <- sum(df$Relapsed_Binary == 1, na.rm = TRUE)
  n_groups <- n_distinct(df$Group)
  model_label <- if (n_distinct(df$Cohort) >= 2) {
    "coxph: MRD group + cohort + months from baseline"
  } else {
    "coxph: MRD group + months from baseline"
  }

  if (nrow(df) < 2 || n_groups < 2 || n_events < min_events_for_cox) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      n_patients = n_distinct(df$Patient),
      n_events = n_events,
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = "not_fit_insufficient_groups_or_events"
    ))
  }

  model_df <- df %>%
    mutate(
      Cohort = factor(Cohort),
      Group = factor(Group, levels = c("Negative", "Positive"))
  )

  fit <- tryCatch(
    {
      model_formula <- if (n_distinct(model_df$Cohort) >= 2) {
        survival::Surv(Time_to_event_months, Relapsed_Binary) ~
          Group + Cohort + months_from_baseline
      } else {
        survival::Surv(Time_to_event_months, Relapsed_Binary) ~
          Group + months_from_baseline
      }
      survival::coxph(
        model_formula,
        data = model_df,
        ties = "breslow"
      )
    },
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      n_patients = n_distinct(df$Patient),
      n_events = n_events,
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = paste0("fit_failed: ", conditionMessage(fit))
    ))
  }

  coefs <- summary(fit)$coefficients
  conf <- suppressMessages(confint(fit))
  term_name <- grep("^GroupPositive$", rownames(coefs), value = TRUE)
  if (length(term_name) != 1) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      n_patients = n_distinct(df$Patient),
      n_events = n_events,
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = "fit_no_group_positive_term"
    ))
  }

  p_value <- unname(coefs[term_name, "Pr(>|z|)"])
  tibble(
    anchor = anchor_label,
    assay_label = assay_label,
    n_patients = n_distinct(df$Patient),
    n_events = n_events,
    model = model_label,
    term = "MRD positive vs negative",
    hazard_ratio = unname(exp(coefs[term_name, "coef"])),
    conf_low = unname(exp(conf[term_name, 1])),
    conf_high = unname(exp(conf[term_name, 2])),
    p_value = p_value,
    p_value_label = format_model_p(p_value),
    status = "fit"
  )
}

fit_delayed_entry_cox <- function(df, assay_label, anchor_label) {
  n_events <- sum(df$Relapsed_Binary == 1, na.rm = TRUE)
  n_groups <- n_distinct(df$Group)
  model_label <- if (n_distinct(df$Cohort) >= 2) {
    "coxph delayed entry: MRD group + cohort; baseline-time scale"
  } else {
    "coxph delayed entry: MRD group; baseline-time scale"
  }

  model_df <- df %>%
    mutate(
      entry_month = months_from_baseline,
      exit_month = months_from_baseline + Time_to_event_months,
      Cohort = factor(Cohort),
      Group = factor(Group, levels = c("Negative", "Positive"))
    ) %>%
    filter(
      is.finite(entry_month),
      is.finite(exit_month),
      exit_month > entry_month,
      !is.na(Relapsed_Binary),
      !is.na(Group)
    )

  if (nrow(model_df) < 2 || n_groups < 2 || n_events < min_events_for_cox) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      n_patients = n_distinct(model_df$Patient),
      n_events = sum(model_df$Relapsed_Binary == 1, na.rm = TRUE),
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = "not_fit_insufficient_groups_or_events"
    ))
  }

  fit <- tryCatch(
    {
      model_formula <- if (n_distinct(model_df$Cohort) >= 2) {
        survival::Surv(entry_month, exit_month, Relapsed_Binary) ~ Group + Cohort
      } else {
        survival::Surv(entry_month, exit_month, Relapsed_Binary) ~ Group
      }
      survival::coxph(model_formula, data = model_df, ties = "breslow")
    },
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      n_patients = n_distinct(model_df$Patient),
      n_events = sum(model_df$Relapsed_Binary == 1, na.rm = TRUE),
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = paste0("fit_failed: ", conditionMessage(fit))
    ))
  }

  coefs <- summary(fit)$coefficients
  conf <- suppressMessages(confint(fit))
  term_name <- grep("^GroupPositive$", rownames(coefs), value = TRUE)
  if (length(term_name) != 1) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      n_patients = n_distinct(model_df$Patient),
      n_events = sum(model_df$Relapsed_Binary == 1, na.rm = TRUE),
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = "fit_no_group_positive_term"
    ))
  }

  p_value <- unname(coefs[term_name, "Pr(>|z|)"])
  tibble(
    anchor = anchor_label,
    assay_label = assay_label,
    n_patients = n_distinct(model_df$Patient),
    n_events = sum(model_df$Relapsed_Binary == 1, na.rm = TRUE),
    model = model_label,
    term = "MRD positive vs negative",
    hazard_ratio = unname(exp(coefs[term_name, "coef"])),
    conf_low = unname(exp(conf[term_name, 1])),
    conf_high = unname(exp(conf[term_name, 2])),
    p_value = p_value,
    p_value_label = format_model_p(p_value),
    status = "fit"
  )
}

fit_time_varying_cox <- function(data, assay_key, assay_col, assay_label,
                                 analysis_type = "time_varying_locf",
                                 event_free_days = 0L,
                                 timepoint_regex = NULL) {
  required_columns(data, assay_col, "analysis input")

  model_label <- if (event_free_days > 0) {
    paste0(
      "coxph time-varying LOCF, event-free >= ", event_free_days,
      " days at assessment: MRD group + cohort + cluster(patient); baseline-time scale"
    )
  } else {
    "coxph time-varying LOCF: MRD group + cohort + cluster(patient); baseline-time scale"
  }

  tv_candidates <- data %>%
    filter(
      !is.na(.data[[assay_col]]),
      !is.na(sample_date),
      !is.na(baseline_date),
      !is.na(censor_date),
      !is.na(relapsed),
      days_from_baseline >= 0,
      sample_date <= censor_date,
      days_to_event_from_sample >= event_free_days,
      !timepoint_info_normalized %in% c(
        "baseline", "diagnosis", "germline_normal", "relapse", "progression"
      )
    ) %>%
    mutate(
      assay_value = as.integer(.data[[assay_col]]),
      event_or_censor_month = days_to_event_from_baseline / days_per_month,
      sample_month = months_from_baseline
    ) %>%
    filter(
      assay_value %in% c(0L, 1L),
      is.finite(sample_month),
      is.finite(event_or_censor_month),
      event_or_censor_month > sample_month
    )

  if (!is.null(timepoint_regex)) {
    tv_candidates <- tv_candidates %>%
      filter(str_detect(timepoint_info_normalized, timepoint_regex))
  }

  interval_df <- tv_candidates %>%
    group_by(Patient, sample_date) %>%
    summarise(
      Cohort = first(Cohort),
      sample_month = first(sample_month),
      event_or_censor_month = first(event_or_censor_month),
      Relapsed_Binary = first(Relapsed_Binary),
      assay_value = as.integer(any(assay_value == 1L, na.rm = TRUE)),
      n_same_date_rows = n(),
      n_same_date_positive_rows = sum(assay_value == 1L, na.rm = TRUE),
      timepoint_info_normalized = paste(sort(unique(timepoint_info_normalized)), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(Patient, sample_month, sample_date) %>%
    group_by(Patient) %>%
    mutate(
      next_sample_month = lead(sample_month),
      start_month = sample_month,
      stop_month = pmin(coalesce(next_sample_month, event_or_censor_month), event_or_censor_month),
      interval_event = as.integer(
        Relapsed_Binary == 1L &
          is.na(next_sample_month) &
          abs(stop_month - event_or_censor_month) < 1e-8
      )
    ) %>%
    ungroup() %>%
    filter(stop_month > start_month) %>%
    mutate(
      Group = factor(
        if_else(assay_value == 1L, "Positive", "Negative"),
        levels = c("Negative", "Positive")
      ),
      Cohort = factor(
        Cohort,
        levels = c("Frontline", "Non-frontline"),
        labels = c("Training Cohort", "Test Cohort")
      ),
      assay_key = assay_key,
      assay_col = assay_col,
      assay_label = assay_label
    )

  interval_path <- file.path(
    source_dir,
    paste0(
      "time_varying_", assay_key,
      if_else(is.null(timepoint_regex), "", "_maintenance_only"),
      if_else(event_free_days > 0L, paste0("_event_free_", event_free_days, "d"), ""),
      "_intervals_source_data.csv"
    )
  )
  write_csv(interval_df, interval_path)

  n_events <- sum(interval_df$interval_event == 1, na.rm = TRUE)
  n_groups <- n_distinct(interval_df$Group)

  if (nrow(interval_df) < 2 || n_groups < 2 || n_events < min_events_for_cox) {
    return(tibble(
      analysis_type = analysis_type,
      assay_label = assay_label,
      n_patients = n_distinct(interval_df$Patient),
      n_intervals = nrow(interval_df),
      n_events = n_events,
      n_positive_intervals = sum(interval_df$Group == "Positive", na.rm = TRUE),
      n_negative_intervals = sum(interval_df$Group == "Negative", na.rm = TRUE),
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = "not_fit_insufficient_groups_or_events",
      interval_source_path = interval_path
    ))
  }

  fit <- tryCatch(
    {
      model_formula <- if (n_distinct(interval_df$Cohort) >= 2) {
        survival::Surv(start_month, stop_month, interval_event) ~
          Group + Cohort + cluster(Patient)
      } else {
        survival::Surv(start_month, stop_month, interval_event) ~
          Group + cluster(Patient)
      }
      survival::coxph(model_formula, data = interval_df, ties = "breslow")
    },
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(tibble(
      analysis_type = analysis_type,
      assay_label = assay_label,
      n_patients = n_distinct(interval_df$Patient),
      n_intervals = nrow(interval_df),
      n_events = n_events,
      n_positive_intervals = sum(interval_df$Group == "Positive", na.rm = TRUE),
      n_negative_intervals = sum(interval_df$Group == "Negative", na.rm = TRUE),
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = paste0("fit_failed: ", conditionMessage(fit)),
      interval_source_path = interval_path
    ))
  }

  coefs <- summary(fit)$coefficients
  conf <- suppressMessages(confint(fit))
  term_name <- grep("^GroupPositive$", rownames(coefs), value = TRUE)
  if (length(term_name) != 1) {
    return(tibble(
      analysis_type = analysis_type,
      assay_label = assay_label,
      n_patients = n_distinct(interval_df$Patient),
      n_intervals = nrow(interval_df),
      n_events = n_events,
      n_positive_intervals = sum(interval_df$Group == "Positive", na.rm = TRUE),
      n_negative_intervals = sum(interval_df$Group == "Negative", na.rm = TRUE),
      model = model_label,
      term = "MRD positive vs negative",
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      p_value_label = NA_character_,
      status = "fit_no_group_positive_term",
      interval_source_path = interval_path
    ))
  }

  p_value <- unname(coefs[term_name, "Pr(>|z|)"])
  tibble(
    analysis_type = analysis_type,
    assay_label = assay_label,
    n_patients = n_distinct(interval_df$Patient),
    n_intervals = nrow(interval_df),
    n_events = n_events,
    n_positive_intervals = sum(interval_df$Group == "Positive", na.rm = TRUE),
    n_negative_intervals = sum(interval_df$Group == "Negative", na.rm = TRUE),
    model = model_label,
    term = "MRD positive vs negative",
    hazard_ratio = unname(exp(coefs[term_name, "coef"])),
    conf_low = unname(exp(conf[term_name, 1])),
    conf_high = unname(exp(conf[term_name, 2])),
    p_value = p_value,
    p_value_label = format_model_p(p_value),
    status = "fit",
    interval_source_path = interval_path
  )
}

make_anchor_df <- function(data, assay_col, anchor_type, target_month = NA_real_,
                           window_months = NA_real_, event_free_days = 0L) {
  required_columns(data, assay_col, "analysis input")

  anchor_candidates <- data %>%
    filter(
      !is.na(.data[[assay_col]]),
      !is.na(sample_date),
      !is.na(baseline_date),
      !is.na(censor_date),
      !is.na(relapsed),
      days_from_baseline >= 0,
      days_to_event_from_sample >= event_free_days
    )

  if (anchor_type == "first_nonbaseline") {
    anchor_candidates <- anchor_candidates %>%
      filter(
        !timepoint_info_normalized %in% c(
          "baseline", "diagnosis", "germline_normal", "relapse", "progression"
        )
      ) %>%
      arrange(Patient, sample_date)
  } else if (anchor_type == "closest_window") {
    if (!is.finite(target_month) || !is.finite(window_months)) {
      stop("closest_window anchors require target_month and window_months.", call. = FALSE)
    }
    anchor_candidates <- anchor_candidates %>%
      filter(abs(months_from_baseline - target_month) <= window_months) %>%
      arrange(Patient, abs(months_from_baseline - target_month), sample_date)
  } else if (anchor_type == "first_maintenance") {
    anchor_candidates <- anchor_candidates %>%
      filter(str_detect(timepoint_info_normalized, "maintenance")) %>%
      arrange(Patient, sample_date)
  } else if (anchor_type == "exact_one_year_maintenance") {
    anchor_candidates <- anchor_candidates %>%
      filter(timepoint_info_normalized == "1yr maintenance") %>%
      arrange(Patient, sample_date)
  } else {
    stop("Unknown anchor_type: ", anchor_type, call. = FALSE)
  }

  anchor_candidates %>%
    group_by(Patient) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      assay_value = as.integer(.data[[assay_col]]),
      Group = factor(
        if_else(assay_value == 1L, "Positive", "Negative"),
        levels = c("Negative", "Positive")
      ),
      Time_to_event_months = days_to_event_from_sample / days_per_month,
      Cohort = factor(
        Cohort,
        levels = c("Frontline", "Non-frontline"),
        labels = c("Training Cohort", "Test Cohort")
      )
    ) %>%
    filter(!is.na(Group), Time_to_event_months >= 0)
}

summarise_anchor_df <- function(df, assay_label, anchor_label) {
  if (nrow(df) == 0) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      Cohort = character(),
      Group = character(),
      n_patients = integer(),
      n_events = integer(),
      median_months_from_baseline = numeric()
    ))
  }

  df %>%
    group_by(Cohort, Group) %>%
    summarise(
      anchor = .env$anchor_label,
      assay_label = .env$assay_label,
      n_patients = n_distinct(Patient),
      n_events = sum(Relapsed_Binary == 1, na.rm = TRUE),
      median_months_from_baseline = median(months_from_baseline, na.rm = TRUE),
      q25_months_from_baseline = quantile(months_from_baseline, 0.25, na.rm = TRUE),
      q75_months_from_baseline = quantile(months_from_baseline, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    relocate(anchor, assay_label)
}

make_assay_patient_eligibility <- function(assay_key, assay_col, assay_label, data) {
  required_columns(data, assay_col, "analysis input")

  excluded_primary_timepoints <- c(
    "baseline", "diagnosis", "germline_normal", "relapse", "progression"
  )

  data %>%
    mutate(
      assay_has_call = !is.na(.data[[assay_col]]),
      assay_is_positive = .data[[assay_col]] == 1L,
      has_complete_pfs = !is.na(baseline_date) & !is.na(censor_date) & !is.na(relapsed),
      pre_event_assay_row = assay_has_call &
        has_complete_pfs &
        !is.na(sample_date) &
        days_from_baseline >= 0 &
        days_to_event_from_sample >= 0,
      first_nonbaseline_eligible_row = pre_event_assay_row &
        !timepoint_info_normalized %in% excluded_primary_timepoints,
      first_nonbaseline_event_free_30d_row = first_nonbaseline_eligible_row &
        days_to_event_from_sample >= 30,
      maintenance_eligible_row = pre_event_assay_row &
        str_detect(timepoint_info_normalized, "maintenance"),
      maintenance_event_free_30d_row = maintenance_eligible_row &
        days_to_event_from_sample >= 30,
      exact_1yr_maintenance_row = pre_event_assay_row &
        timepoint_info_normalized == "1yr maintenance",
      baseline_or_diagnosis_call_row = assay_has_call &
        timepoint_info_normalized %in% c("baseline", "diagnosis"),
      relapse_or_progression_call_row = assay_has_call &
        timepoint_info_normalized %in% c("relapse", "progression"),
      post_pfs_maintenance_call_row = assay_has_call &
        has_complete_pfs &
        !is.na(sample_date) &
        str_detect(timepoint_info_normalized, "maintenance") &
        days_to_event_from_sample < 0
    ) %>%
    group_by(Patient, Cohort) %>%
    summarise(
      assay_key = .env$assay_key,
      assay_col = .env$assay_col,
      assay_label = .env$assay_label,
      has_complete_pfs = any(has_complete_pfs, na.rm = TRUE),
      has_any_assay_call = any(assay_has_call, na.rm = TRUE),
      has_any_positive_assay_call = any(assay_is_positive, na.rm = TRUE),
      has_pre_event_assay_call = any(pre_event_assay_row, na.rm = TRUE),
      has_first_nonbaseline_eligible_call = any(first_nonbaseline_eligible_row, na.rm = TRUE),
      has_first_nonbaseline_event_free_30d_call = any(first_nonbaseline_event_free_30d_row, na.rm = TRUE),
      has_first_maintenance_eligible_call = any(maintenance_eligible_row, na.rm = TRUE),
      has_first_maintenance_event_free_30d_call = any(maintenance_event_free_30d_row, na.rm = TRUE),
      has_exact_1yr_maintenance_call = any(exact_1yr_maintenance_row, na.rm = TRUE),
      has_baseline_or_diagnosis_call = any(baseline_or_diagnosis_call_row, na.rm = TRUE),
      has_relapse_or_progression_call = any(relapse_or_progression_call_row, na.rm = TRUE),
      has_post_pfs_maintenance_call = any(post_pfs_maintenance_call_row, na.rm = TRUE),
      n_rows_with_assay_call = sum(assay_has_call, na.rm = TRUE),
      n_pre_event_assay_rows = sum(pre_event_assay_row, na.rm = TRUE),
      n_first_nonbaseline_eligible_rows = sum(first_nonbaseline_eligible_row, na.rm = TRUE),
      n_maintenance_eligible_rows = sum(maintenance_eligible_row, na.rm = TRUE),
      n_post_pfs_maintenance_rows = sum(post_pfs_maintenance_call_row, na.rm = TRUE),
      .groups = "drop"
    )
}

make_km_plot <- function(df, assay_label, anchor_label, output_stub) {
  group_counts <- df %>%
    count(Group, name = "n")

  can_plot <- nrow(df) >= 2 * min_group_n_for_km &&
    n_distinct(df$Group) == 2 &&
    all(group_counts$n >= min_group_n_for_km)

  if (!can_plot) {
    return(tibble(
      anchor = anchor_label,
      assay_label = assay_label,
      figure_path = NA_character_,
      status = "not_plotted_insufficient_group_size"
    ))
  }

  fit <- survival::survfit(
    survival::Surv(Time_to_event_months, Relapsed_Binary) ~ Group,
    data = df
  )

  pval_display <- tryCatch(
    {
      lr <- survival::survdiff(
        survival::Surv(Time_to_event_months, Relapsed_Binary) ~ Group,
        data = df
      )
      p_value <- stats::pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
      if (is.finite(p_value)) {
        paste0("Log-rank p = ", format_model_p(p_value))
      } else {
        FALSE
      }
    },
    error = function(e) FALSE
  )

  title_text <- paste0(
    "PFS by ", assay_label,
    "\nAll evaluable patients - ", anchor_label
  )

  km <- survminer::ggsurvplot(
    fit,
    data = df,
    pval = pval_display,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.title = "Number at risk",
    risk.table.title.theme = element_text(hjust = 0),
    break.time.by = 12,
    palette = c("black", "red"),
    legend.title = "MRD status",
    legend.labs = c("MRD-", "MRD+"),
    xlab = "Time since MRD assessment (months)",
    ylab = "Progression-free survival",
    title = stringr::str_wrap(title_text, width = 58),
    risk.table.height = 0.25,
    ggtheme = theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13)
      )
  )

  km$table <- km$table +
    theme(
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0, face = "plain")
    )

  km$plot <- km$plot +
    theme(axis.title.x = element_blank())

  combined <- ggpubr::ggarrange(km$plot, km$table, ncol = 1, heights = c(3, 1))
  figure_path <- file.path(output_dir, paste0(output_stub, ".png"))
  figure_path_dated <- file.path(output_dir, paste0(output_stub, "_", date_tag, ".png"))

  ggsave(figure_path, combined, width = 7.5, height = 7, dpi = 500)
  ggsave(figure_path_dated, combined, width = 7.5, height = 7, dpi = 500)

  tibble(
    anchor = anchor_label,
    assay_label = assay_label,
    figure_path = figure_path,
    status = "plotted"
  )
}

if (!file.exists(input_calls_rds)) {
  stop("Missing cfWGS calls input: ", input_calls_rds, call. = FALSE)
}

calls_df <- readRDS(input_calls_rds)
pfs_df <- read_pfs_table(input_pfs_rds)

next_event_endpoint_resources <- load_next_event_endpoint_resources(
  pfs_path = input_pfs_rds,
  relapse_dates_path = input_relapse_dates_rds,
  followup_rds_path = input_followup_rds,
  followup_csv_path = input_followup_csv,
  latest_dates_paths = c(input_latest_dates_updated_csv, input_latest_dates_csv),
  sample_data = calls_df,
  patient_col = "Patient",
  sample_date_col = "Date"
)

assay_defs <- tibble::tribble(
  ~assay_key, ~assay_col, ~assay_label,
  "BM_cfWGS_cVAF", "BM_zscore_only_detection_rate_call", "BM-derived cfWGS (cVAF model)",
  "Blood_cfWGS_sites", "Blood_zscore_only_sites_call", "cfDNA-derived cfWGS (sites model)",
  "Blood_cfWGS_combined", "Blood_plus_fragment_call", "cfDNA-derived cfWGS (combined model)",
  "MFC", "Flow_Binary", "MFC",
  "clonoSEQ", "Adaptive_Binary", "clonoSEQ"
)

first_nonbaseline_artifact_map <- tibble::tribble(
  ~assay_key, ~artifact_id, ~description,
  "BM_cfWGS_cVAF", "FIG3F", "All-evaluable first-nonbaseline companion for the BM-derived cfWGS KM panel.",
  "Blood_cfWGS_sites", "FIG4E", "All-evaluable first-nonbaseline companion for the cfDNA-derived sites-model cfWGS KM panel.",
  "MFC", "EDFIG6C", "All-evaluable first-nonbaseline companion for the MFC KM panel.",
  "clonoSEQ", "EDFIG6D", "All-evaluable first-nonbaseline companion for the clonoSEQ KM panel."
)

first_maintenance_artifact_map <- tibble::tribble(
  ~assay_key, ~artifact_id, ~description,
  "BM_cfWGS_cVAF", "FIG3F", "All-evaluable first-maintenance companion for the BM-derived cfWGS KM panel.",
  "Blood_cfWGS_sites", "FIG4E", "All-evaluable first-maintenance companion for the cfDNA-derived sites-model cfWGS KM panel.",
  "MFC", "EDFIG6C", "All-evaluable first-maintenance companion for the MFC KM panel.",
  "clonoSEQ", "EDFIG6D", "All-evaluable first-maintenance companion for the clonoSEQ KM panel."
)

required_columns(
  calls_df,
  c("Patient", "Date", "Cohort", "timepoint_info", assay_defs$assay_col),
  "cfWGS calls table"
)

analysis_df <- calls_df %>%
  mutate(
    Patient = as.character(Patient),
    sample_date = as.Date(Date),
    Cohort = as.character(Cohort),
    timepoint_info_normalized = tolower(as.character(timepoint_info))
  ) %>%
  filter(Cohort %in% cohorts_to_include) %>%
  add_next_event_endpoint(
    endpoint_resources = next_event_endpoint_resources,
    sample_date_col = "sample_date",
    event_grace_days = 30L
  ) %>%
  mutate(
    days_from_baseline = as.numeric(sample_date - baseline_date),
    months_from_baseline = days_from_baseline / days_per_month,
    censor_date = endpoint_date,
    days_to_event_from_baseline = as.numeric(endpoint_date - baseline_date),
    days_to_event_from_sample = endpoint_days_from_sample,
    Relapsed_Binary = as.integer(endpoint_status),
    relapsed = as.integer(endpoint_status)
  )

input_qc <- analysis_df %>%
  summarise(
    n_rows = n(),
    n_patients = n_distinct(Patient),
    n_patients_missing_baseline = n_distinct(Patient[is.na(baseline_date)]),
    n_patients_missing_censor = n_distinct(Patient[is.na(censor_date)]),
    n_rows_after_event_or_censor = sum(days_to_event_from_sample < 0, na.rm = TRUE),
    n_rows_after_first_pfs = sum(first_pfs_days_from_sample < 0, na.rm = TRUE),
    n_rows_ignoring_prior_progression = sum(endpoint_ignores_prior_progression, na.rm = TRUE),
    n_rows_using_later_progression_after_first_pfs =
      sum(endpoint_uses_later_progression_after_first_pfs, na.rm = TRUE),
    n_rows_missing_sample_date = sum(is.na(sample_date))
  )

write_csv(input_qc, file.path(source_dir, "all_evaluable_first_nonbaseline_input_qc.csv"))

next_event_endpoint_audit <- analysis_df %>%
  filter(endpoint_ignores_prior_progression |
           endpoint_uses_later_progression_after_first_pfs |
           force_relapse_sample_day0) %>%
  select(
    Patient, Cohort, Timepoint, Sample_Code, sample_date, timepoint_info,
    first_pfs_date, first_pfs_days_from_sample,
    latest_prior_progression_date, n_prior_progressions_before_sample,
    next_progression_date, endpoint_date, endpoint_status,
    endpoint_type, endpoint_source, endpoint_days_from_sample,
    endpoint_uses_later_progression_after_first_pfs, force_relapse_sample_day0
  ) %>%
  arrange(Patient, sample_date)

write_csv(
  next_event_endpoint_audit,
  file.path(source_dir, "all_evaluable_next_event_endpoint_audit.csv")
)

next_event_endpoint_audit_compact <- next_event_endpoint_audit %>%
  mutate(
    endpoint_status_label = if_else(endpoint_status == 1L, "event", "censor"),
    ignored_reason = case_when(
      force_relapse_sample_day0 ~ paste0(
        "relapse-labeled sample assigned event at day 0; original first PFS was ",
        abs(first_pfs_days_from_sample), " days before sample"
      ),
      n_prior_progressions_before_sample > 0 ~ paste0(
        n_prior_progressions_before_sample,
        " prior progression(s) before sample; latest prior=",
        latest_prior_progression_date,
        "; first PFS was ", abs(first_pfs_days_from_sample),
        " days before sample"
      ),
      TRUE ~ NA_character_
    ),
    replacement_endpoint = case_when(
      force_relapse_sample_day0 ~ paste0(
        "relapse sample event on ", endpoint_date, " (0 days from sample)"
      ),
      endpoint_status == 1L ~ paste0(
        "next progression on ", endpoint_date, " (",
        endpoint_days_from_sample, " days after sample)"
      ),
      TRUE ~ paste0(
        "censor on ", endpoint_date, " from ", endpoint_source, " (",
        endpoint_days_from_sample, " days after sample)"
      )
    )
  ) %>%
  transmute(
    Patient, Cohort, Timepoint, Sample_Code, sample_date, timepoint_info,
    ignored_reason,
    old_first_pfs_date = first_pfs_date,
    old_first_pfs_days_from_sample = first_pfs_days_from_sample,
    replacement_endpoint, endpoint_status_label, next_progression_date,
    endpoint_date, endpoint_days_from_sample,
    endpoint_uses_later_progression_after_first_pfs,
    force_relapse_sample_day0
  )

write_csv(
  next_event_endpoint_audit_compact,
  file.path(source_dir, "all_evaluable_next_event_endpoint_audit_compact.csv")
)

next_event_endpoint_patient_summary <- next_event_endpoint_audit_compact %>%
  mutate(
    Cohort = if_else(Cohort == "Frontline", "Training", "Test"),
    affected_sample = paste0(Timepoint, " ", sample_date),
    endpoint_summary = paste0(
      endpoint_status_label, " ", endpoint_date, " (",
      endpoint_days_from_sample, "d)"
    )
  ) %>%
  group_by(Patient, Cohort) %>%
  summarise(
    n_affected_samples = n(),
    affected_samples = paste(affected_sample, collapse = "; "),
    original_first_pfs = paste(
      unique(paste0(
        old_first_pfs_date, " (", old_first_pfs_days_from_sample, "d)"
      )),
      collapse = "; "
    ),
    replacement_endpoints = paste(endpoint_summary, collapse = "; "),
    replacement_status = paste(unique(endpoint_status_label), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Cohort, Patient)

write_csv(
  next_event_endpoint_patient_summary,
  file.path(source_dir, "all_evaluable_next_event_endpoint_patient_summary.csv")
)

patient_sampling_patterns <- analysis_df %>%
  group_by(Patient, Cohort) %>%
  summarise(
    has_complete_pfs = any(!is.na(baseline_date) & !is.na(censor_date) & !is.na(relapsed), na.rm = TRUE),
    n_rows = n(),
    n_rows_with_sample_date = sum(!is.na(sample_date)),
    n_pre_event_postbaseline_nonprogression_rows = sum(
      !is.na(sample_date) &
        !is.na(baseline_date) &
        !is.na(censor_date) &
        !is.na(relapsed) &
        days_from_baseline >= 0 &
        days_to_event_from_sample >= 0 &
        !timepoint_info_normalized %in% c(
          "baseline", "diagnosis", "germline_normal", "relapse", "progression"
        ),
      na.rm = TRUE
    ),
    n_pre_event_maintenance_rows = sum(
      !is.na(sample_date) &
        !is.na(baseline_date) &
        !is.na(censor_date) &
        !is.na(relapsed) &
        days_from_baseline >= 0 &
        days_to_event_from_sample >= 0 &
        str_detect(timepoint_info_normalized, "maintenance"),
      na.rm = TRUE
    ),
    n_post_pfs_maintenance_rows = sum(
      !is.na(sample_date) &
        !is.na(baseline_date) &
        !is.na(censor_date) &
        !is.na(relapsed) &
        days_to_event_from_sample < 0 &
        str_detect(timepoint_info_normalized, "maintenance"),
      na.rm = TRUE
    ),
    timepoint_pattern = paste(
      sort(unique(coalesce(timepoint_info_normalized, "missing"))),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(Cohort, Patient)

patient_sampling_pattern_summary <- patient_sampling_patterns %>%
  mutate(
    has_pre_event_postbaseline_nonprogression_row =
      n_pre_event_postbaseline_nonprogression_rows > 0,
    has_pre_event_maintenance_row = n_pre_event_maintenance_rows > 0,
    has_post_pfs_maintenance_row = n_post_pfs_maintenance_rows > 0
  ) %>%
  count(
    Cohort,
    timepoint_pattern,
    has_complete_pfs,
    has_pre_event_postbaseline_nonprogression_row,
    has_pre_event_maintenance_row,
    has_post_pfs_maintenance_row,
    name = "n_patients"
  ) %>%
  arrange(Cohort, desc(n_patients), timepoint_pattern)

patient_assay_eligibility_flags <- purrr::pmap_dfr(
  assay_defs,
  make_assay_patient_eligibility,
  data = analysis_df
)

endpoint_denominator_audit <- patient_assay_eligibility_flags %>%
  group_by(assay_key, assay_label, Cohort) %>%
  summarise(
    total_patients = n_distinct(Patient),
    pfs_evaluable_patients = sum(has_complete_pfs, na.rm = TRUE),
    patients_with_any_assay_call = sum(has_any_assay_call, na.rm = TRUE),
    patients_with_any_pre_event_assay_call = sum(has_pre_event_assay_call, na.rm = TRUE),
    first_nonbaseline_eligible_patients = sum(has_first_nonbaseline_eligible_call, na.rm = TRUE),
    first_nonbaseline_event_free_30d_patients =
      sum(has_first_nonbaseline_event_free_30d_call, na.rm = TRUE),
    first_maintenance_eligible_patients = sum(has_first_maintenance_eligible_call, na.rm = TRUE),
    first_maintenance_event_free_30d_patients =
      sum(has_first_maintenance_event_free_30d_call, na.rm = TRUE),
    exact_1yr_maintenance_eligible_patients = sum(has_exact_1yr_maintenance_call, na.rm = TRUE),
    patients_with_baseline_or_diagnosis_call = sum(has_baseline_or_diagnosis_call, na.rm = TRUE),
    patients_with_relapse_or_progression_call = sum(has_relapse_or_progression_call, na.rm = TRUE),
    patients_with_post_pfs_maintenance_call = sum(has_post_pfs_maintenance_call, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(assay_key, Cohort)

post_pfs_maintenance_samples <- analysis_df %>%
  filter(
    !is.na(sample_date),
    !is.na(censor_date),
    str_detect(timepoint_info_normalized, "maintenance"),
    days_to_event_from_sample < 0
  ) %>%
  transmute(
    Patient,
    Cohort,
    Sample_Code = as.character(Sample_Code),
    sample_date,
    censor_date,
    days_after_pfs_event = -days_to_event_from_sample,
    timepoint_info,
    timepoint_info_normalized,
    !!!syms(assay_defs$assay_col)
  ) %>%
  arrange(Cohort, Patient, sample_date)

write_csv(
  patient_sampling_patterns,
  file.path(output_dir, "all_evaluable_patient_sampling_patterns.csv")
)
write_csv(
  patient_sampling_pattern_summary,
  file.path(output_dir, "all_evaluable_patient_sampling_pattern_summary.csv")
)
write_csv(
  patient_assay_eligibility_flags,
  file.path(output_dir, "all_evaluable_patient_assay_eligibility_flags.csv")
)
write_csv(
  endpoint_denominator_audit,
  file.path(output_dir, "all_evaluable_endpoint_denominator_audit.csv")
)
write_csv(
  post_pfs_maintenance_samples,
  file.path(source_dir, "all_evaluable_post_pfs_maintenance_samples.csv")
)

anchor_defs <- tibble::tribble(
  ~anchor_key, ~anchor_label, ~anchor_type, ~target_month, ~window_months, ~event_free_days, ~plot_primary,
  "first_nonbaseline_primary",
  "first prospective non-baseline assessment",
  "first_nonbaseline", NA_real_, NA_real_, 0L, TRUE,
  "first_nonbaseline_event_free_30d",
  "first prospective non-baseline assessment, event-free >=30 days",
  "first_nonbaseline", NA_real_, NA_real_, 30L, TRUE,
  "first_maintenance_primary",
  "first evaluable maintenance-phase assessment",
  "first_maintenance", NA_real_, NA_real_, 0L, TRUE,
  "first_maintenance_event_free_30d",
  "first evaluable maintenance-phase assessment, event-free >=30 days",
  "first_maintenance", NA_real_, NA_real_, 30L, TRUE,
  "exact_1yr_maintenance_reference",
  "exact clinically labeled one-year maintenance assessment",
  "exact_one_year_maintenance", NA_real_, NA_real_, 0L, FALSE,
  "closest_6mo_pm3",
  "closest assessment to 6 months (+/-3 months)",
  "closest_window", 6, 3, 0L, FALSE,
  "closest_12mo_pm3",
  "closest assessment to 12 months (+/-3 months)",
  "closest_window", 12, 3, 0L, FALSE
)

analysis_outputs <- purrr::pmap(
  assay_defs,
  function(assay_key, assay_col, assay_label) {
    purrr::pmap(
      anchor_defs,
      function(anchor_key, anchor_label, anchor_type, target_month,
               window_months, event_free_days, plot_primary) {
        anchor_df <- make_anchor_df(
          data = analysis_df,
          assay_col = assay_col,
          anchor_type = anchor_type,
          target_month = target_month,
          window_months = window_months,
          event_free_days = event_free_days
        ) %>%
          mutate(
            assay_key = assay_key,
            assay_col = assay_col,
            assay_label = assay_label,
            anchor_key = anchor_key,
            anchor_label = anchor_label
          )

        source_path <- file.path(
          source_dir,
          paste0(anchor_key, "_", assay_key, "_source_data.csv")
        )
        write_csv(anchor_df, source_path)

        counts <- summarise_anchor_df(anchor_df, assay_label, anchor_label) %>%
          mutate(
            assay_key = assay_key,
            anchor_key = anchor_key,
            source_path = source_path
          )

        model <- fit_adjusted_cox(anchor_df, assay_label, anchor_label) %>%
          mutate(assay_key = assay_key, anchor_key = anchor_key)

        delayed_entry_model <- fit_delayed_entry_cox(anchor_df, assay_label, anchor_label) %>%
          mutate(assay_key = assay_key, anchor_key = anchor_key)

        figure <- if (isTRUE(plot_primary)) {
          make_km_plot(
            anchor_df,
            assay_label = assay_label,
            anchor_label = anchor_label,
            output_stub = paste0("KM_", assay_key, "_", anchor_key)
          )
        } else {
          tibble(
            anchor = anchor_label,
            assay_label = assay_label,
            figure_path = NA_character_,
            status = "availability_summary_only"
          )
        }

        list(
          counts = counts,
          model = model,
          delayed_entry_model = delayed_entry_model,
          figure = figure %>%
            mutate(assay_key = assay_key, anchor_key = anchor_key)
        )
      }
    )
  }
) %>%
  purrr::flatten()

anchor_counts <- analysis_outputs %>%
  purrr::map("counts") %>%
  dplyr::bind_rows()

cox_models <- analysis_outputs %>%
  purrr::map("model") %>%
  dplyr::bind_rows()

delayed_entry_models <- analysis_outputs %>%
  purrr::map("delayed_entry_model") %>%
  dplyr::bind_rows()

time_varying_models <- bind_rows(
  purrr::pmap_dfr(
    assay_defs,
    fit_time_varying_cox,
    data = analysis_df,
    analysis_type = "time_varying_locf",
    event_free_days = 0L
  ),
  purrr::pmap_dfr(
    assay_defs,
    fit_time_varying_cox,
    data = analysis_df,
    analysis_type = "time_varying_locf_event_free_30d",
    event_free_days = 30L
  ),
  purrr::pmap_dfr(
    assay_defs,
    fit_time_varying_cox,
    data = analysis_df,
    analysis_type = "time_varying_maintenance_locf",
    event_free_days = 0L,
    timepoint_regex = "maintenance"
  ),
  purrr::pmap_dfr(
    assay_defs,
    fit_time_varying_cox,
    data = analysis_df,
    analysis_type = "time_varying_maintenance_locf_event_free_30d",
    event_free_days = 30L,
    timepoint_regex = "maintenance"
  )
)

model_comparison_summary <- bind_rows(
  cox_models %>%
    transmute(
      analysis_type = "sample_time_anchor_cox",
      anchor = anchor,
      assay_label = assay_label,
      n_patients = n_patients,
      n_intervals = NA_integer_,
      n_events = n_events,
      model = model,
      term = term,
      hazard_ratio = hazard_ratio,
      conf_low = conf_low,
      conf_high = conf_high,
      p_value = p_value,
      p_value_label = p_value_label,
      status = status,
      assay_key = assay_key,
      anchor_key = anchor_key
    ),
  delayed_entry_models %>%
    transmute(
      analysis_type = "delayed_entry_anchor_cox",
      anchor = anchor,
      assay_label = assay_label,
      n_patients = n_patients,
      n_intervals = NA_integer_,
      n_events = n_events,
      model = model,
      term = term,
      hazard_ratio = hazard_ratio,
      conf_low = conf_low,
      conf_high = conf_high,
      p_value = p_value,
      p_value_label = p_value_label,
      status = status,
      assay_key = assay_key,
      anchor_key = anchor_key
    ),
  time_varying_models %>%
    mutate(
      anchor = if_else(
        str_detect(analysis_type, "maintenance"),
        "all eligible serial maintenance-phase assessments",
        "all eligible serial non-baseline assessments"
      ),
      assay_key = assay_defs$assay_key[match(assay_label, assay_defs$assay_label)],
      anchor_key = analysis_type
    ) %>%
    transmute(
      analysis_type = analysis_type,
      anchor = anchor,
      assay_label = assay_label,
      n_patients = n_patients,
      n_intervals = n_intervals,
      n_events = n_events,
      model = model,
      term = term,
      hazard_ratio = hazard_ratio,
      conf_low = conf_low,
      conf_high = conf_high,
      p_value = p_value,
      p_value_label = p_value_label,
      status = status,
      assay_key = assay_key,
      anchor_key = anchor_key
    )
)

figure_manifest <- analysis_outputs %>%
  purrr::map("figure") %>%
  dplyr::bind_rows()

write_csv(
  anchor_counts,
  file.path(output_dir, "all_evaluable_first_nonbaseline_anchor_counts.csv")
)
write_csv(
  cox_models,
  file.path(output_dir, "all_evaluable_first_nonbaseline_cox_models.csv")
)
write_csv(
  delayed_entry_models,
  file.path(output_dir, "all_evaluable_first_nonbaseline_delayed_entry_cox_models.csv")
)
write_csv(
  time_varying_models,
  file.path(output_dir, "all_evaluable_time_varying_cox_models.csv")
)
write_csv(
  model_comparison_summary,
  file.path(output_dir, "all_evaluable_model_comparison_summary.csv")
)
write_csv(
  figure_manifest,
  file.path(output_dir, "all_evaluable_first_nonbaseline_figure_manifest.csv")
)

figure_manifest %>%
  filter(
    status == "plotted",
    str_starts(anchor_key, "first_nonbaseline"),
    !is.na(figure_path),
    file.exists(figure_path)
  ) %>%
  inner_join(first_nonbaseline_artifact_map, by = "assay_key") %>%
  mutate(
    role = case_when(
      anchor_key == "first_nonbaseline_event_free_30d" ~ "all_evaluable_first_nonbaseline_event_free_30d_figure_panel_png",
      TRUE ~ "all_evaluable_first_nonbaseline_figure_panel_png"
    )
  ) %>%
  pwalk(function(anchor, assay_label, figure_path, status, assay_key,
                 anchor_key, artifact_id, description, role) {
    ms_copy_artifact(
      source_path = figure_path,
      artifact_id = artifact_id,
      role = role,
      description = paste0(description, " Anchor: ", anchor, "."),
      script_name = "4_1B_Build_all_evaluable_first_nonbaseline_KM.R"
    )
  })

figure_manifest %>%
  filter(
    status == "plotted",
    str_starts(anchor_key, "first_maintenance"),
    !is.na(figure_path),
    file.exists(figure_path)
  ) %>%
  inner_join(first_maintenance_artifact_map, by = "assay_key") %>%
  mutate(
    role = case_when(
      anchor_key == "first_maintenance_event_free_30d" ~ "all_evaluable_first_maintenance_event_free_30d_figure_panel_png",
      TRUE ~ "all_evaluable_first_maintenance_figure_panel_png"
    )
  ) %>%
  pwalk(function(anchor, assay_label, figure_path, status, assay_key,
                 anchor_key, artifact_id, description, role) {
    ms_copy_artifact(
      source_path = figure_path,
      artifact_id = artifact_id,
      role = role,
      description = paste0(description, " Anchor: ", anchor, "."),
      script_name = "4_1B_Build_all_evaluable_first_nonbaseline_KM.R"
    )
  })

cohort_counts_text <- patient_sampling_patterns %>%
  count(Cohort, name = "n_patients") %>%
  mutate(label = paste0(Cohort, " n=", n_patients)) %>%
  pull(label) %>%
  paste(collapse = "; ")

test_pre_event_maintenance_n <- patient_sampling_patterns %>%
  filter(Cohort == "Non-frontline") %>%
  summarise(n = sum(n_pre_event_maintenance_rows > 0, na.rm = TRUE)) %>%
  pull(n)

test_post_pfs_maintenance_n <- patient_sampling_patterns %>%
  filter(Cohort == "Non-frontline") %>%
  summarise(
    n = sum(
      n_post_pfs_maintenance_rows > 0 & n_pre_event_maintenance_rows == 0,
      na.rm = TRUE
    )
  ) %>%
  pull(n)

test_first_maintenance_denominator_text <- endpoint_denominator_audit %>%
  filter(Cohort == "Non-frontline") %>%
  arrange(match(assay_key, assay_defs$assay_key)) %>%
  transmute(
    label = paste0(
      assay_label, " n=", first_maintenance_eligible_patients,
      " (exact 1yr n=", exact_1yr_maintenance_eligible_patients, ")"
    )
  ) %>%
  pull(label) %>%
  paste(collapse = "; ")

test_first_nonbaseline_denominator_text <- endpoint_denominator_audit %>%
  filter(Cohort == "Non-frontline") %>%
  arrange(match(assay_key, assay_defs$assay_key)) %>%
  transmute(label = paste0(assay_label, " n=", first_nonbaseline_eligible_patients)) %>%
  pull(label) %>%
  paste(collapse = "; ")

method_note <- c(
  "# All-evaluable longitudinal and maintenance KM companion analysis",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Recommended interpretation",
  "",
  "Use the first evaluable maintenance-phase assessment as the biologically aligned expansion of the original one-year maintenance analysis. This preserves the successful maintenance estimand while allowing the test/non-frontline cohort's broader Maintenance label to contribute when sampled before the first PFS event.",
  "",
  "Use the first prospective non-baseline assessment as a broad descriptive sensitivity analysis, not as the primary replacement for the maintenance result. It maximizes patient inclusion across the training/frontline and test/non-frontline cohorts but mixes early response, post-transplant, maintenance, and irregular follow-up timing.",
  "",
  "Because sampling intervals differ across cohorts, every KM curve is descriptive and conditional on reaching the assessment. The adjusted Cox table should be shown or cited alongside the KM because it adjusts the MRD association for cohort and months from baseline.",
  "",
  "Maintenance is defined from the clinical `timepoint_info` label (case-insensitive match to `maintenance`), not from a post hoc calendar window. This includes explicit yearly maintenance labels in the frontline cohort and the broader Maintenance label used in the non-frontline cohort. Cohort-specific counts and sampling-time distributions must accompany this analysis because label granularity differs between cohorts.",
  "",
  "For timing-robust inference, use the model comparison table. The delayed-entry Cox keeps each patient-level anchor but uses baseline as the time scale and enters each patient at the assessment date. The general time-varying Cox uses all eligible serial post-baseline, non-progression-labeled assessments with last-observation-carried-forward MRD status, baseline as the time scale, and patient-clustered robust standard errors. Maintenance-only serial models apply the same method only to clinically labeled maintenance assessments. Parallel event-free-30-day sensitivities exclude assessments collected within 30 days of event/censor to reduce immediate pre-progression sampling concerns.",
  "",
  "The KM x-axis is time since the selected MRD assessment. The delayed-entry and time-varying Cox models use time since baseline as the analysis time scale. For samples collected after an earlier recorded progression, the endpoint is the next known progression after the sample; if no later progression is known, the row is censored at the latest available follow-up/sample date after the sample.",
  "",
  "Do not present the fixed 6-month or 12-month windows as the primary all-evaluable result unless the manuscript needs a fixed-time landmark. Those windows are cleaner with respect to calendar time but have materially smaller denominators in the current data. The maintenance-phase analysis preserves the successful biological estimand more directly.",
  "",
  "## Current denominator interpretation",
  "",
  paste0("- Total patients in the included cohorts with source rows: ", n_distinct(analysis_df$Patient), " (", cohort_counts_text, ")."),
  paste0("- Test/non-frontline patients with any pre-PFS maintenance-labeled row: ", test_pre_event_maintenance_n, "."),
  paste0("- Test/non-frontline patients with maintenance-labeled rows only after the recorded first PFS event: ", test_post_pfs_maintenance_n, ". These rows are now evaluated against the next known progression after the sample, or censored if no later progression is known."),
  paste0("- Test/non-frontline first-maintenance assay denominators: ", test_first_maintenance_denominator_text, "."),
  paste0("- Test/non-frontline broad first-nonbaseline assay denominators: ", test_first_nonbaseline_denominator_text, "."),
  "",
  "The apparent gap between the total cohort size and the KM denominators is therefore expected: the landmark analysis requires a post-baseline sample with the relevant assay call and a non-negative next-event/censor interval after that sample. Baseline-only, relapse-only, baseline-plus-relapse without later follow-up, or assay-missing patients are not evaluable for that specific endpoint.",
  "",
  "## Primary exclusions",
  "",
  "- Baseline, diagnosis, and germline-normal rows are not eligible as non-baseline landmarks.",
  "- Relapse- and progression-labeled rows are not eligible for the primary anchor because they may be outcome-triggered.",
  "- Rows after an earlier recorded PFS event are retained only when a valid next-event/censor endpoint can be assigned after the sample date.",
  "- The event-free-30-day anchor is exported as a leakage sensitivity analysis.",
  "",
  "## Key outputs",
  "",
  "- `all_evaluable_first_nonbaseline_anchor_counts.csv`: patient counts, events, MRD status, and assessment timing by cohort.",
  "- `all_evaluable_first_nonbaseline_cox_models.csv`: MRD-positive hazard ratios adjusted for cohort and months from baseline.",
  "- `all_evaluable_first_nonbaseline_delayed_entry_cox_models.csv`: first-anchor hazard ratios using baseline time with delayed entry at the assessment date.",
  "- `all_evaluable_time_varying_cox_models.csv`: serial-assessment LOCF time-varying Cox models, including the event-free-30-day sensitivity.",
  "- `all_evaluable_model_comparison_summary.csv`: combined anchor, delayed-entry, and time-varying model summary.",
  "- `all_evaluable_endpoint_denominator_audit.csv`: assay-by-cohort endpoint eligibility counts, including first-maintenance and broad first-nonbaseline denominators.",
  "- `all_evaluable_patient_sampling_patterns.csv`: patient-level source-row and timepoint-pattern inventory.",
  "- `all_evaluable_patient_sampling_pattern_summary.csv`: compressed cohort-level inventory explaining why total patients do not equal endpoint-evaluable patients.",
  "- `all_evaluable_patient_assay_eligibility_flags.csv`: patient-by-assay eligibility flags for denominator review.",
  "- `source_data/*_source_data.csv`: patient-level source rows for each assay/anchor.",
  "- `source_data/time_varying_*_intervals_source_data.csv`: patient-interval source rows for the serial time-varying models.",
  "- `source_data/all_evaluable_post_pfs_maintenance_samples.csv`: maintenance-labeled rows occurring after the first PFS event; these require next-event/censor handling.",
  "- `source_data/all_evaluable_next_event_endpoint_audit.csv`: rows whose endpoint ignores at least one prior progression and uses a later progression or censoring date.",
  "- `KM_*_first_maintenance_primary.png`: expanded maintenance-phase KM figures.",
  "- `KM_*_first_maintenance_event_free_30d.png`: event-free-30-day sensitivity KM figures for the expanded maintenance endpoint.",
  "- `KM_*_first_nonbaseline_primary.png`: primary all-evaluable KM figures.",
  "- `KM_*_first_nonbaseline_event_free_30d.png`: event-free-30-day sensitivity KM figures."
)

writeLines(
  method_note,
  file.path(output_dir, "all_evaluable_first_nonbaseline_method_note.md")
)

cat("\nAll-evaluable longitudinal and maintenance KM companion analysis complete.\n")
cat("Output directory:", output_dir, "\n")
cat("Primary figure count:", sum(figure_manifest$status == "plotted"), "\n")
cat("Primary/source table rows:", nrow(anchor_counts), "\n")
cat("Time-varying models fit:", sum(time_varying_models$status == "fit"), "\n")
