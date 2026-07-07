################################################################################
## All-Evaluable First Non-Baseline Landmark KM Analysis
##
## Purpose:
##   Build companion Kaplan-Meier curves that include evaluable patients from
##   both the training/frontline and test/non-frontline cohorts. Because the two
##   cohorts have different sampling schedules, this script uses a patient-level
##   first prospective non-baseline assessment as the primary descriptive
##   landmark, and exports timing/QC tables plus Cox models adjusted for cohort
##   and months from baseline.
##
## Scientific estimand:
##   Among patients who are progression-free and under observation at their first
##   non-baseline, non-progression-labeled cfWGS assessment, describe subsequent
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
##     additional_all_evaluable_first_nonbaseline_km/
##       KM_*_first_nonbaseline_primary.png
##       source_data/*_source_data.csv
##       all_evaluable_first_nonbaseline_anchor_counts.csv
##       all_evaluable_first_nonbaseline_cox_models.csv
##       all_evaluable_first_nonbaseline_delayed_entry_cox_models.csv
##       all_evaluable_time_varying_cox_models.csv
##       all_evaluable_model_comparison_summary.csv
##       all_evaluable_first_nonbaseline_method_note.md
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

input_calls_rds <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6.rds"
input_pfs_rds <- "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"

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
                                 event_free_days = 0L) {
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

assay_defs <- tibble::tribble(
  ~assay_key, ~assay_col, ~assay_label,
  "BM_cfWGS_cVAF", "BM_zscore_only_detection_rate_call", "BM-derived cfWGS (cVAF model)",
  "Blood_cfWGS_sites", "Blood_zscore_only_sites_call", "cfDNA-derived cfWGS (sites model)",
  "Blood_cfWGS_combined", "Blood_plus_fragment_call", "cfDNA-derived cfWGS (combined model)",
  "MFC", "Flow_Binary", "MFC",
  "clonoSEQ", "Adaptive_Binary", "clonoSEQ"
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
  left_join(pfs_df, by = "Patient") %>%
  mutate(
    days_from_baseline = as.numeric(sample_date - baseline_date),
    months_from_baseline = days_from_baseline / days_per_month,
    days_to_event_from_baseline = as.numeric(censor_date - baseline_date),
    days_to_event_from_sample = as.numeric(censor_date - sample_date),
    Relapsed_Binary = as.integer(relapsed)
  )

input_qc <- analysis_df %>%
  summarise(
    n_rows = n(),
    n_patients = n_distinct(Patient),
    n_patients_missing_baseline = n_distinct(Patient[is.na(baseline_date)]),
    n_patients_missing_censor = n_distinct(Patient[is.na(censor_date)]),
    n_rows_after_event_or_censor = sum(days_to_event_from_sample < 0, na.rm = TRUE),
    n_rows_missing_sample_date = sum(is.na(sample_date))
  )

write_csv(input_qc, file.path(source_dir, "all_evaluable_first_nonbaseline_input_qc.csv"))

anchor_defs <- tibble::tribble(
  ~anchor_key, ~anchor_label, ~anchor_type, ~target_month, ~window_months, ~event_free_days, ~plot_primary,
  "first_nonbaseline_primary",
  "first prospective non-baseline assessment",
  "first_nonbaseline", NA_real_, NA_real_, 0L, TRUE,
  "first_nonbaseline_event_free_30d",
  "first prospective non-baseline assessment, event-free >=30 days",
  "first_nonbaseline", NA_real_, NA_real_, 30L, TRUE,
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
      anchor = "all eligible serial non-baseline assessments",
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

method_note <- c(
  "# All-evaluable first non-baseline KM companion analysis",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Recommended interpretation",
  "",
  "Use the first prospective non-baseline assessment as the main broad all-evaluable KM companion. This maximizes patient inclusion across the training/frontline and test/non-frontline cohorts while preserving one contribution per patient.",
  "",
  "Because sampling intervals differ across cohorts, the KM curve is descriptive and conditional on reaching the assessment. The adjusted Cox table should be shown or cited alongside the KM because it adjusts the MRD association for cohort and months from baseline.",
  "",
  "For timing-robust inference, use the model comparison table. The delayed-entry Cox keeps the first non-baseline patient-level anchor but uses baseline as the time scale and enters each patient at the assessment date. The time-varying Cox uses all eligible serial post-baseline, non-progression-labeled assessments with last-observation-carried-forward MRD status, baseline as the time scale, and patient-clustered robust standard errors. A parallel event-free-30-day serial sensitivity excludes assessments collected within 30 days of event/censor to reduce immediate pre-progression sampling concerns.",
  "",
  "The KM x-axis is time since the selected MRD assessment. The delayed-entry and time-varying Cox models use time since baseline as the analysis time scale.",
  "",
  "Do not present the fixed 6-month or 12-month windows as the primary all-evaluable result unless the manuscript needs a fixed-time landmark. Those windows are cleaner with respect to sampling time but have materially smaller denominators in the current data.",
  "",
  "## Primary exclusions",
  "",
  "- Baseline, diagnosis, and germline-normal rows are not eligible as non-baseline landmarks.",
  "- Relapse- and progression-labeled rows are not eligible for the primary anchor because they may be outcome-triggered.",
  "- Rows after the PFS event/censor date are excluded.",
  "- The event-free-30-day anchor is exported as a leakage sensitivity analysis.",
  "",
  "## Key outputs",
  "",
  "- `all_evaluable_first_nonbaseline_anchor_counts.csv`: patient counts, events, MRD status, and assessment timing by cohort.",
  "- `all_evaluable_first_nonbaseline_cox_models.csv`: MRD-positive hazard ratios adjusted for cohort and months from baseline.",
  "- `all_evaluable_first_nonbaseline_delayed_entry_cox_models.csv`: first-anchor hazard ratios using baseline time with delayed entry at the assessment date.",
  "- `all_evaluable_time_varying_cox_models.csv`: serial-assessment LOCF time-varying Cox models, including the event-free-30-day sensitivity.",
  "- `all_evaluable_model_comparison_summary.csv`: combined anchor, delayed-entry, and time-varying model summary.",
  "- `source_data/*_source_data.csv`: patient-level source rows for each assay/anchor.",
  "- `source_data/time_varying_*_intervals_source_data.csv`: patient-interval source rows for the serial time-varying models.",
  "- `KM_*_first_nonbaseline_primary.png`: primary all-evaluable KM figures.",
  "- `KM_*_first_nonbaseline_event_free_30d.png`: event-free-30-day sensitivity KM figures."
)

writeLines(
  method_note,
  file.path(output_dir, "all_evaluable_first_nonbaseline_method_note.md")
)

cat("\nAll-evaluable first non-baseline KM companion analysis complete.\n")
cat("Output directory:", output_dir, "\n")
cat("Primary figure count:", sum(figure_manifest$status == "plotted"), "\n")
cat("Primary/source table rows:", nrow(anchor_counts), "\n")
cat("Time-varying models fit:", sum(time_varying_models$status == "fit"), "\n")
