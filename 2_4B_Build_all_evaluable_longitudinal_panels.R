# ===============================================================================
# Script: 2_4B_Build_all_evaluable_longitudinal_panels.R
# ===============================================================================
#
# Purpose:
#   Build separate all-evaluable versions of the locked Figure 2B/2C/2D
#   longitudinal panels. The original manuscript panels in
#   2_4_Longitudinal_features_analysis.R intentionally filter to
#   Cohort == "Frontline"; this companion keeps both Frontline/train and
#   Non-frontline/test patients after the same evidence-of-disease masking.
#
# Outputs:
#   Scripts_2025/Final_Scripts/final_manuscript_objects/
#     01_main_figures/Figure_2/<panel>/
#       <panel>__all_evaluable_figure_panel_png__*.png
#     additional_all_evaluable_longitudinal_panels/
#       Figure_2B_all_evaluable_BM_zscore_longitudinal.png
#       Figure_2B_all_evaluable_BM_zscore_longitudinal_pseudolog.png
#       Figure_2C_all_evaluable_blood_zscore_longitudinal_pseudolog.png
#       Figure_2C_all_evaluable_blood_zscore_longitudinal.png
#       Figure_2C_all_evaluable_blood_zscore_longitudinal_capped300.png
#       Figure_2D_all_evaluable_fragmentomics_longitudinal.png
#       Figure_2D_all_evaluable_fragmentomics_longitudinal_pseudolog.png
#       Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal_pseudolog.png
#       Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal_pseudolog.png
#       Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal_pseudolog.png
#       *_source_data.csv
#       all_evaluable_longitudinal_panel_qc.csv
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/2_4B_Build_all_evaluable_longitudinal_panels.R
# ===============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
})

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

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

output_dir <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects",
  "additional_all_evaluable_longitudinal_panels"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure2_dir <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects",
  "01_main_figures", "Figure_2"
)

required_files <- c(
  "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds",
  final_cohort_assignment_path("rds"),
  "baseline_high_quality_patients_updated.csv",
  "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds",
  "Exported_data_tables_clinical/Relapse_dates_full_updated.rds"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste0("  - ", missing_files, collapse = "\n"),
    call. = FALSE
  )
}

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

cohort_df <- load_final_cohort_assignment()
feature_df <- readRDS(
  "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds"
)

coverage_fallback_path <- "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated8.rds"
coverage_fallback_audit <- tibble(
  coverage_fallback_path = coverage_fallback_path,
  fallback_file_exists = file.exists(coverage_fallback_path),
  missing_coverage_before = sum(is.na(feature_df$Mean.Coverage)),
  missing_coverage_after = NA_integer_,
  n_restored_patient_date_values = NA_integer_
)
if (file.exists(coverage_fallback_path) && "Mean.Coverage" %in% names(feature_df)) {
  coverage_fallback <- readRDS(coverage_fallback_path) %>%
    select(Patient, Date, Mean.Coverage_fallback = Mean.Coverage) %>%
    mutate(Date = as.Date(Date))

  duplicate_coverage_keys <- coverage_fallback %>%
    count(Patient, Date, name = "n") %>%
    filter(n > 1)

  if (nrow(duplicate_coverage_keys) > 0) {
    stop(
      "Cannot restore Mean.Coverage: duplicate Patient-Date keys in ",
      coverage_fallback_path,
      call. = FALSE
    )
  }

  missing_coverage_before <- sum(is.na(feature_df$Mean.Coverage))
  feature_df <- feature_df %>%
    mutate(Date = as.Date(Date)) %>%
    left_join(coverage_fallback, by = c("Patient", "Date")) %>%
    mutate(Mean.Coverage = coalesce(Mean.Coverage, Mean.Coverage_fallback)) %>%
    select(-Mean.Coverage_fallback)

  missing_coverage_after <- sum(is.na(feature_df$Mean.Coverage))
  coverage_fallback_audit <- coverage_fallback_audit %>%
    mutate(
      missing_coverage_after = .env$missing_coverage_after,
      n_restored_patient_date_values = .env$missing_coverage_before - .env$missing_coverage_after
    )
  message(
    "Mean.Coverage fallback restored ",
    missing_coverage_before - missing_coverage_after,
    " missing Patient-Date values from ",
    coverage_fallback_path,
    "."
  )
}

good_pts <- read.csv(
  "baseline_high_quality_patients_updated.csv",
  stringsAsFactors = FALSE
)
baseline_dates <- readRDS(
  "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"
)

next_event_endpoint_resources <- load_next_event_endpoint_resources(
  pfs_path = "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds",
  relapse_dates_path = "Exported_data_tables_clinical/Relapse_dates_full_updated.rds",
  followup_rds_path = "Exported_data_tables_clinical/patient_followup_dates_updated.rds",
  followup_csv_path = "Exported_data_tables_clinical/patient_followup_dates_updated.csv",
  latest_dates_paths = c(
    "Exported_data_tables_clinical/latest_dates_per_patient_updated.csv",
    "Exported_data_tables_clinical/latest_dates_per_patient.csv"
  ),
  sample_data = feature_df,
  patient_col = "Patient",
  sample_date_col = "Date"
)

required_columns(cohort_df, c("Patient", "Cohort"), "cohort_df")
required_columns(
  good_pts,
  c(
    "Patient",
    "Date",
    "WGS_Evidence_of_Disease_BM_cells",
    "WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed"
  ),
  "good_pts"
)
required_columns(baseline_dates, c("Patient", "baseline_date"), "baseline_dates")
required_columns(
  feature_df,
  c(
    "Patient", "Date", "Num_days_to_closest_relapse",
    "detect_rate_BM", "z_score_detection_rate_BM", "sites_rate_BM",
    "zscore_BM", "detect_rate_blood", "z_score_detection_rate_blood",
    "sites_rate_blood", "zscore_blood", "FS", "Mean.Coverage",
    "Proportion.Short", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
  ),
  "feature_df"
)

bm_good_pts <- good_pts %>%
  filter(WGS_Evidence_of_Disease_BM_cells == 1) %>%
  pull(Patient) %>%
  unique()

cfdna_good_pts <- good_pts %>%
  filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
  pull(Patient) %>%
  unique()

mutation_source_anchor_dates <- good_pts %>%
  transmute(
    Patient,
    Mutation_Source_Baseline_Date = as.Date(Date),
    Mutation_Source_Baseline_Timepoint = Timepoint,
    Mutation_Source_Baseline_Sample_Code = Sample_Code,
    Mutation_Source_Baseline_Timepoint_Info = timepoint_info
  )

duplicate_anchor_patients <- mutation_source_anchor_dates %>%
  count(Patient, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_anchor_patients) > 0) {
  stop(
    "Mutation-source baseline anchor table has duplicate Patient rows: ",
    paste(duplicate_anchor_patients$Patient, collapse = ", "),
    call. = FALSE
  )
}

bm_feats <- c("zscore_BM", "z_score_detection_rate_BM", "detect_rate_BM", "sites_rate_BM")
blood_feats <- c(
  "zscore_blood", "z_score_detection_rate_blood",
  "detect_rate_blood", "sites_rate_blood"
)

dat <- feature_df %>%
  left_join(cohort_df, by = "Patient") %>%
  mutate(Date = as.Date(Date)) %>%
  filter(!is.na(Cohort)) %>%
  filter(Cohort %in% c("Frontline", "Non-frontline")) %>%
  filter(
    !is.na(zscore_BM) |
      !is.na(zscore_blood) |
      !is.na(FS) |
      !is.na(Mean.Coverage) |
      !is.na(Proportion.Short) |
      !is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA)
  ) %>%
  mutate(
    across(all_of(bm_feats), ~ if_else(Patient %in% bm_good_pts, .x, NA_real_)),
    across(all_of(blood_feats), ~ if_else(Patient %in% cfdna_good_pts, .x, NA_real_))
  ) %>%
  apply_manual_longitudinal_feature_masks("img127_blood_tf_mask") %>%
  left_join(
    baseline_dates %>%
      transmute(Patient, Clinical_Baseline_Date = as.Date(baseline_date)),
    by = "Patient"
  ) %>%
  left_join(mutation_source_anchor_dates, by = "Patient") %>%
  mutate(
    Baseline_Date = Mutation_Source_Baseline_Date,
    Weeks_Since_Baseline = as.numeric(difftime(Date, Baseline_Date, units = "weeks")),
    Weeks_Since_Clinical_Baseline = as.numeric(difftime(Date, Clinical_Baseline_Date, units = "weeks")),
    Weeks_Since_Baseline = case_when(
      Weeks_Since_Baseline >= -2 & Weeks_Since_Baseline < 0 ~ 0,
      TRUE ~ Weeks_Since_Baseline
    )
  ) %>%
  add_next_event_endpoint(
    endpoint_resources = next_event_endpoint_resources,
    sample_date_col = "Date",
    event_grace_days = 30L
  ) %>%
  mutate(
    Original_Num_days_to_closest_relapse = Num_days_to_closest_relapse,
    Num_days_to_closest_relapse = endpoint_days_from_sample,
    Relapsed_Binary = as.integer(endpoint_status)
  )

longitudinal_next_event_endpoint_audit <- dat %>%
  filter(endpoint_ignores_prior_progression |
           endpoint_uses_later_progression_after_first_pfs |
           force_relapse_sample_day0) %>%
  select(
    Patient, Cohort, Date, timepoint_info,
    Original_Num_days_to_closest_relapse, Num_days_to_closest_relapse,
    first_pfs_date, first_pfs_days_from_sample,
    latest_prior_progression_date, n_prior_progressions_before_sample,
    next_progression_date, endpoint_date, endpoint_status,
    endpoint_type, endpoint_source, endpoint_days_from_sample,
    endpoint_uses_later_progression_after_first_pfs, force_relapse_sample_day0
  ) %>%
  arrange(Patient, Date)

write_csv(
  longitudinal_next_event_endpoint_audit,
  file.path(output_dir, "all_evaluable_longitudinal_next_event_endpoint_audit.csv")
)

time_anchor_exclusions <- dat %>%
  filter(is.na(Weeks_Since_Baseline)) %>%
  select(
    Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
    Mutation_Source_Baseline_Timepoint, Mutation_Source_Baseline_Sample_Code,
    Mutation_Source_Baseline_Timepoint_Info, Num_days_to_closest_relapse,
    all_of(c(bm_feats, blood_feats, "FS", "Mean.Coverage"))
  )

if (nrow(time_anchor_exclusions) > 0) {
  write_csv(
    time_anchor_exclusions,
    file.path(output_dir, "all_evaluable_rows_excluded_missing_weeks_since_baseline.csv")
  )
}

dat <- dat %>%
  filter(!is.na(Weeks_Since_Baseline))

time_anchor_audit <- dat %>%
  distinct(
    Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
    Weeks_Since_Baseline, Weeks_Since_Clinical_Baseline,
    Mutation_Source_Baseline_Timepoint, Mutation_Source_Baseline_Sample_Code,
    Mutation_Source_Baseline_Timepoint_Info
  ) %>%
  mutate(
    anchor_shift_weeks = Weeks_Since_Clinical_Baseline - Weeks_Since_Baseline,
    anchor_shift_days = as.numeric(Baseline_Date - Clinical_Baseline_Date)
  ) %>%
  arrange(desc(abs(anchor_shift_weeks)), Patient, Date)

if (nrow(dat) == 0) {
  stop("No evaluable longitudinal rows remain after all-evaluable filters.", call. = FALSE)
}

relapse_legend_labels <- c(
  `FALSE` = "No relapse within 180 days of blood draw",
  `TRUE` = "Relapse within 180 days of blood draw"
)
cohort_linetypes <- c("Frontline" = "solid", "Non-frontline" = "22")
cohort_display_labels <- c("Frontline" = "Training Cohort", "Non-frontline" = "Test Cohort")

make_relapse_facet_labels <- function(plot_data) {
  required_cols <- c("Patient", "Date", "patient_relapse180")
  missing_cols <- setdiff(required_cols, names(plot_data))
  if (length(missing_cols) > 0) {
    stop(
      "Cannot build relapse facet labels; missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  patient_counts <- plot_data %>%
    filter(!is.na(patient_relapse180)) %>%
    distinct(Patient, patient_relapse180) %>%
    count(patient_relapse180, name = "n_patients")

  sample_counts <- plot_data %>%
    filter(!is.na(patient_relapse180), !is.na(Date)) %>%
    distinct(Patient, Date, patient_relapse180) %>%
    count(patient_relapse180, name = "n_samples")

  count_for <- function(status) {
    value <- patient_counts$n_patients[
      as.character(patient_counts$patient_relapse180) == as.character(status)
    ]
    if (length(value) == 0) 0L else as.integer(value[[1]])
  }

  sample_count_for <- function(status) {
    value <- sample_counts$n_samples[
      as.character(sample_counts$patient_relapse180) == as.character(status)
    ]
    if (length(value) == 0) 0L else as.integer(value[[1]])
  }

  c(
    `FALSE` = paste0(
      "No relapse within 180 days of blood draw\n(n = ",
      count_for(FALSE), " patients; ", sample_count_for(FALSE), " samples)"
    ),
    `TRUE` = paste0(
      "Relapse within 180 days of blood draw\n(n = ",
      count_for(TRUE), " patients; ", sample_count_for(TRUE), " samples)"
    )
  )
}

add_patient_relapse_flag <- function(plot_df) {
  plot_df %>%
    group_by(Patient) %>%
    mutate(patient_relapse180 = any(replace_na(relapse_within_180, FALSE))) %>%
    ungroup() %>%
    mutate(
      relapse_within_180 = replace_na(relapse_within_180, FALSE),
      patient_relapse180 = factor(patient_relapse180, levels = c(FALSE, TRUE)),
      Cohort = factor(Cohort, levels = c("Frontline", "Non-frontline"))
    )
}

flag_followup_relapse_points <- function(plot_df) {
  plot_df %>%
    arrange(Patient, Metric, Weeks_Since_Baseline, Date) %>%
    group_by(Patient, Metric) %>%
    mutate(
      # Baseline is defined by the mutation-source time anchor (week 0), not
      # by being the first assay-evaluable observation for a patient. Thus, a
      # patient's only available observation can still be a relapse-associated
      # follow-up point when it occurs after week 0.
      is_effective_baseline_point = near(Weeks_Since_Baseline, 0),
      relapse_within_180 = if_else(
        # True week-0 baseline observations are always shown in black. A point
        # is red only when it is after week 0 and was collected 0--180 days
        # before the patient's next documented relapse/progression.
        !is_effective_baseline_point &
          Relapsed_Binary == 1L &
          Num_days_to_closest_relapse >= 0 &
          Num_days_to_closest_relapse <= 180,
        TRUE,
        FALSE,
        missing = FALSE
      )
    ) %>%
    ungroup()
}

make_segment_df <- function(df, y_col = "Value") {
  df %>%
    arrange(Patient, Weeks_Since_Baseline) %>%
    group_by(Patient) %>%
    mutate(
      x = Weeks_Since_Baseline,
      y = .data[[y_col]],
      xend = lead(Weeks_Since_Baseline),
      yend = lead(.data[[y_col]]),
      seg_relapse = lead(relapse_within_180, default = FALSE),
      seg_cohort = first(Cohort)
    ) %>%
    filter(!is.na(x), !is.na(xend), !is.na(y), !is.na(yend)) %>%
    ungroup()
}

make_metric_panel <- function(df, metric, panel_title, y_label, cap_at = NA_real_,
                              reverse_display = FALSE, show_cap_labels = TRUE,
                              x_cap_at = NA_real_,
                              y_transform = c("identity", "pseudo_log"),
                              y_transform_sigma = 1,
                              y_transform_breaks = NULL,
                              display_as_percent = FALSE,
                              overcap_shape = 17,
                              prefix_cap_axis_label = TRUE) {
  y_transform <- match.arg(y_transform)
  metric_df <- df %>% filter(Metric == metric)
  if (nrow(metric_df) == 0) {
    stop("No rows available for metric: ", metric, call. = FALSE)
  }
  facet_labels <- make_relapse_facet_labels(metric_df)
  has_x_cap <- is.finite(x_cap_at)
  percent_axis_labels <- scales::label_percent(accuracy = 1)
  percent_value_labels <- function(x) {
    sub("\\.0%$", "%", scales::percent(x, accuracy = 0.1))
  }

  if (is.finite(cap_at) && y_transform != "identity" && isTRUE(reverse_display)) {
    stop(
      "A capped transformed y-axis cannot be combined with a reversed display.",
      call. = FALSE
    )
  }

  if (is.finite(cap_at)) {
    cap_label_accuracy <- if (cap_at < 1) 0.01 else 1
    cap_y_transform <- if (y_transform == "pseudo_log") {
      scales::pseudo_log_trans(sigma = y_transform_sigma, base = 10)
    } else {
      "identity"
    }
    cap_y_breaks <- if (y_transform == "pseudo_log") {
      y_transform_breaks
    } else {
      pretty(c(0, cap_at))
    }
    cap_y_labels <- function(x) {
      formatted_x <- if (isTRUE(display_as_percent)) {
        percent_axis_labels(x)
      } else {
        scales::number(x, accuracy = cap_label_accuracy)
      }
      ifelse(
        abs(x - cap_at) < 1e-8,
        paste0(if (isTRUE(prefix_cap_axis_label)) ">" else "", formatted_x),
        formatted_x
      )
    }
    metric_plot_df <- metric_df %>%
      mutate(
        overcap = Value > cap_at,
        Value_plot = pmin(Value, cap_at),
        x_overcap = has_x_cap & Weeks_Since_Baseline > x_cap_at,
        Weeks_Since_Baseline_plot = if (has_x_cap) pmin(Weeks_Since_Baseline, x_cap_at) else Weeks_Since_Baseline,
        label_val = if_else(
          overcap,
          if (isTRUE(display_as_percent)) {
            percent_value_labels(Value)
          } else {
            scales::number(Value, accuracy = cap_label_accuracy)
          },
          NA_character_
        )
      )
    seg_df <- metric_df %>%
      arrange(Patient, Weeks_Since_Baseline) %>%
      group_by(Patient) %>%
      mutate(
        x = if (has_x_cap) pmin(Weeks_Since_Baseline, x_cap_at) else Weeks_Since_Baseline,
        y = pmin(pmax(Value, 0), cap_at),
        xend = lead(if (has_x_cap) pmin(Weeks_Since_Baseline, x_cap_at) else Weeks_Since_Baseline),
        yend = lead(pmin(pmax(Value, 0), cap_at)),
        seg_relapse = lead(relapse_within_180, default = FALSE),
        seg_cohort = first(Cohort)
      ) %>%
      filter(!is.na(x), !is.na(xend), !is.na(y), !is.na(yend), xend >= x) %>%
      ungroup()

    return(
      ggplot(metric_plot_df, aes(Weeks_Since_Baseline_plot, Value_plot, group = Patient)) +
        geom_segment(
          data = seg_df,
          aes(x = x, y = y, xend = xend, yend = yend,
              colour = seg_relapse, linetype = seg_cohort),
          linewidth = 0.35,
          alpha = 0.6
        ) +
        geom_point(aes(color = relapse_within_180, shape = overcap), size = 1.6, alpha = 0.8) +
        geom_text(
          data = ~ filter(.x, overcap, show_cap_labels),
          aes(label = label_val),
          vjust = -0.75,
          size = 2
        ) +
        facet_wrap(
          ~ patient_relapse180,
          nrow = 1,
          labeller = labeller(patient_relapse180 = facet_labels)
        ) +
        scale_color_manual(
          values = c(`FALSE` = "black", `TRUE` = "red"),
          labels = relapse_legend_labels,
          name = NULL
        ) +
        scale_linetype_manual(
          values = cohort_linetypes,
          labels = cohort_display_labels,
          name = "Cohort"
        ) +
        scale_shape_manual(
          values = c(`FALSE` = 16, `TRUE` = overcap_shape),
          labels = c(
            `FALSE` = paste0(
              "<=",
              if (isTRUE(display_as_percent)) percent_axis_labels(cap_at) else cap_at
            ),
            `TRUE` = paste0(
              ">",
              if (isTRUE(display_as_percent)) percent_axis_labels(cap_at) else cap_at,
              " (capped)"
            )
          ),
          name = NULL
        ) +
        guides(shape = "none") +
        scale_y_continuous(
          trans = cap_y_transform,
          limits = c(0, cap_at * 1.04),
          breaks = cap_y_breaks,
          labels = cap_y_labels
        ) +
        (if (has_x_cap) {
          scale_x_continuous(
            limits = c(0, x_cap_at),
            breaks = sort(unique(c(pretty(c(0, x_cap_at)), x_cap_at))),
            labels = function(x) ifelse(abs(x - x_cap_at) < 1e-8, paste0(">", x_cap_at), x)
          )
        }) +
        labs(title = panel_title, x = "Weeks Since Baseline", y = y_label) +
        theme_classic(base_size = 11) +
        theme(strip.background = element_rect(fill = "grey95", colour = NA))
    )
  }

  metric_plot_df <- metric_df %>%
    mutate(
      Weeks_Since_Baseline_plot = if (has_x_cap) pmin(Weeks_Since_Baseline, x_cap_at) else Weeks_Since_Baseline
    )
  seg_df <- metric_df %>%
    arrange(Patient, Weeks_Since_Baseline) %>%
    group_by(Patient) %>%
    mutate(
      x = if (has_x_cap) pmin(Weeks_Since_Baseline, x_cap_at) else Weeks_Since_Baseline,
      y = Value,
      xend = lead(if (has_x_cap) pmin(Weeks_Since_Baseline, x_cap_at) else Weeks_Since_Baseline),
      yend = lead(Value),
      seg_relapse = lead(relapse_within_180, default = FALSE),
      seg_cohort = first(Cohort)
    ) %>%
    filter(!is.na(x), !is.na(xend), !is.na(y), !is.na(yend), xend >= x) %>%
    ungroup()
  panel <- ggplot(metric_plot_df, aes(Weeks_Since_Baseline_plot, Value, group = Patient)) +
    geom_segment(
      data = seg_df,
      aes(x = x, y = y, xend = xend, yend = yend,
          colour = seg_relapse, linetype = seg_cohort),
      linewidth = 0.35,
      alpha = 0.6
    ) +
    geom_point(aes(color = relapse_within_180), size = 1.6, alpha = 0.8) +
    facet_wrap(
      ~ patient_relapse180,
      nrow = 1,
      labeller = labeller(patient_relapse180 = facet_labels)
    ) +
    scale_color_manual(
      values = c(`FALSE` = "black", `TRUE` = "red"),
      labels = relapse_legend_labels,
      name = NULL
    ) +
    scale_linetype_manual(
      values = cohort_linetypes,
      labels = cohort_display_labels,
      name = "Cohort"
    ) +
    (if (has_x_cap) {
      scale_x_continuous(
        limits = c(0, x_cap_at),
        breaks = sort(unique(c(pretty(c(0, x_cap_at)), x_cap_at))),
        labels = function(x) ifelse(abs(x - x_cap_at) < 1e-8, paste0(">", x_cap_at), x)
      )
    }) +
    labs(title = panel_title, x = "Weeks Since Baseline", y = y_label) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA))

  if (y_transform == "identity") {
    if (isTRUE(reverse_display)) {
      panel <- panel + scale_y_reverse(
        labels = if (isTRUE(display_as_percent)) percent_axis_labels else waiver()
      )
    } else if (isTRUE(display_as_percent)) {
      panel <- panel + scale_y_continuous(labels = percent_axis_labels)
    }
  }
  if (y_transform == "pseudo_log") {
    pseudo_log_transform <- scales::pseudo_log_trans(
      sigma = y_transform_sigma,
      base = 10
    )
    display_transform <- if (isTRUE(reverse_display)) {
      scales::transform_compose(pseudo_log_transform, scales::reverse_trans())
    } else {
      pseudo_log_transform
    }
    panel <- panel + scale_y_continuous(
      trans = display_transform,
      breaks = y_transform_breaks,
      labels = if (isTRUE(display_as_percent)) {
        percent_axis_labels
      } else {
        scales::label_number(big.mark = ",")
      }
    )
  }
  panel
}

make_mutation_plot_df <- function(dat, assay) {
  if (assay == "BM") {
    out <- dat %>%
      select(
        Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
        Weeks_Since_Baseline, Weeks_Since_Clinical_Baseline,
        Mutation_Source_Baseline_Timepoint,
        Mutation_Source_Baseline_Sample_Code,
        Mutation_Source_Baseline_Timepoint_Info,
        Num_days_to_closest_relapse, Relapsed_Binary,
        cVAF_z = z_score_detection_rate_BM,
        sites_z = zscore_BM
      )
  } else if (assay == "blood") {
    out <- dat %>%
      select(
        Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
        Weeks_Since_Baseline, Weeks_Since_Clinical_Baseline,
        Mutation_Source_Baseline_Timepoint,
        Mutation_Source_Baseline_Sample_Code,
        Mutation_Source_Baseline_Timepoint_Info,
        Num_days_to_closest_relapse, Relapsed_Binary,
        cVAF_z = z_score_detection_rate_blood,
        sites_z = zscore_blood
      )
  } else {
    stop("Unknown assay: ", assay, call. = FALSE)
  }

  out %>%
    pivot_longer(
      cols = c(cVAF_z, sites_z),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    drop_na(Value) %>%
    flag_followup_relapse_points() %>%
    add_patient_relapse_flag()
}

make_mutation_panel <- function(plot_df, title, cvaf_cap_at = NA_real_,
                                sites_cap_at = NA_real_,
                                show_cap_labels = TRUE,
                                x_cap_at = NA_real_,
                                cvaf_y_transform = c("identity", "pseudo_log"),
                                sites_y_transform = c("identity", "pseudo_log")) {
  cvaf_y_transform <- match.arg(cvaf_y_transform)
  sites_y_transform <- match.arg(sites_y_transform)
  p_cvaf <- make_metric_panel(
    plot_df,
    metric = "cVAF_z",
    panel_title = "Cumulative VAF Z-score",
    y_label = if (cvaf_y_transform == "pseudo_log") {
      "Cumulative VAF (Z; pseudo-log)"
    } else {
      "Cumulative VAF (Z)"
    },
    cap_at = cvaf_cap_at,
    show_cap_labels = show_cap_labels,
    x_cap_at = x_cap_at,
    y_transform = cvaf_y_transform,
    y_transform_sigma = 1,
    y_transform_breaks = if (cvaf_y_transform == "pseudo_log") {
      c(-10, -1, 0, 1, 10, 100, 1000, 10000)
    } else {
      NULL
    }
  )
  p_sites <- make_metric_panel(
    plot_df,
    metric = "sites_z",
    panel_title = "Proportion of Sites Detected Z-score",
    y_label = if (sites_y_transform == "pseudo_log") {
      "Sites detected (Z; pseudo-log)"
    } else {
      "Prop. Mutant Sites Detected (Z)"
    },
    cap_at = sites_cap_at,
    show_cap_labels = show_cap_labels,
    x_cap_at = x_cap_at,
    y_transform = sites_y_transform,
    y_transform_sigma = 1,
    y_transform_breaks = if (sites_y_transform == "pseudo_log") {
      c(-10, -1, 0, 1, 10, 100, 1000, 10000)
    } else {
      NULL
    }
  )

  p_cvaf + p_sites +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    ) &
    theme(legend.position = "bottom")
}

make_fragmentomics_plot_df <- function(dat) {
  dat %>%
    select(
      Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
      Weeks_Since_Baseline, Weeks_Since_Clinical_Baseline,
      Mutation_Source_Baseline_Timepoint,
      Mutation_Source_Baseline_Sample_Code,
      Mutation_Source_Baseline_Timepoint_Info,
      Num_days_to_closest_relapse, Relapsed_Binary,
      FS, Mean.Coverage
    ) %>%
    pivot_longer(
      cols = c(FS, Mean.Coverage),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    drop_na(Value) %>%
    flag_followup_relapse_points() %>%
    add_patient_relapse_flag()
}

make_fragmentomics_panel <- function(plot_df, x_cap_at = NA_real_,
                                     y_transform = c("identity", "pseudo_log")) {
  y_transform <- match.arg(y_transform)
  p_fs <- make_metric_panel(
    plot_df,
    metric = "FS",
    panel_title = "Fragment-size score",
    y_label = if (y_transform == "pseudo_log") {
      "Fragment-size score (pseudo-log)"
    } else {
      "Fragment-size score"
    },
    x_cap_at = x_cap_at,
    y_transform = y_transform,
    y_transform_sigma = 0.05,
    y_transform_breaks = if (y_transform == "pseudo_log") {
      c(-1, -0.5, -0.1, 0, 0.1)
    } else {
      NULL
    }
  )
  p_coverage <- make_metric_panel(
    plot_df,
    metric = "Mean.Coverage",
    panel_title = "cfDNA coverage at MM active regulatory sites",
    y_label = if (y_transform == "pseudo_log") {
      "Mean cfDNA coverage (pseudo-log)"
    } else {
      "Mean cfDNA coverage (MM regs)"
    },
    reverse_display = TRUE,
    x_cap_at = x_cap_at,
    y_transform = y_transform,
    y_transform_sigma = 0.01,
    y_transform_breaks = if (y_transform == "pseudo_log") {
      c(0.9, 0.95, 1, 1.05)
    } else {
      NULL
    }
  )

  p_fs + p_coverage +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Fragmentomic tumour-agnostic longitudinal tracking: Relapsed vs Non-relapsed patients",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    ) &
    theme(legend.position = "bottom")
}

make_raw_mutation_plot_df <- function(dat, assay) {
  if (assay == "BM") {
    out <- dat %>%
      select(
        Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
        Weeks_Since_Baseline, Weeks_Since_Clinical_Baseline,
        Mutation_Source_Baseline_Timepoint,
        Mutation_Source_Baseline_Sample_Code,
        Mutation_Source_Baseline_Timepoint_Info,
        Num_days_to_closest_relapse, Relapsed_Binary,
        cVAF = detect_rate_BM,
        sites = sites_rate_BM
      )
  } else if (assay == "blood") {
    out <- dat %>%
      select(
        Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
        Weeks_Since_Baseline, Weeks_Since_Clinical_Baseline,
        Mutation_Source_Baseline_Timepoint,
        Mutation_Source_Baseline_Sample_Code,
        Mutation_Source_Baseline_Timepoint_Info,
        Num_days_to_closest_relapse, Relapsed_Binary,
        cVAF = detect_rate_blood,
        sites = sites_rate_blood
      )
  } else {
    stop("Unknown assay: ", assay, call. = FALSE)
  }

  out %>%
    pivot_longer(
      cols = c(cVAF, sites),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    drop_na(Value) %>%
    flag_followup_relapse_points() %>%
    add_patient_relapse_flag()
}

make_raw_mutation_panel <- function(plot_df, title, cvaf_cap_at = NA_real_,
                                    sites_cap_at = NA_real_,
                                    x_cap_at = NA_real_,
                                    y_transform = c("identity", "pseudo_log")) {
  y_transform <- match.arg(y_transform)
  p_cvaf <- make_metric_panel(
    plot_df,
    metric = "cVAF",
    panel_title = "Cumulative VAF",
    y_label = "Cumulative VAF (%)",
    cap_at = cvaf_cap_at,
    x_cap_at = x_cap_at,
    y_transform = y_transform,
    y_transform_sigma = 0.001,
    y_transform_breaks = if (y_transform == "pseudo_log") {
      c(0, 0.001, 0.01, 0.1, 1)
    } else {
      NULL
    },
    display_as_percent = TRUE,
    overcap_shape = 18,
    prefix_cap_axis_label = FALSE
  )
  p_sites <- make_metric_panel(
    plot_df,
    metric = "sites",
    panel_title = "Proportion of Sites Detected",
    y_label = "Proportion of Mutant Sites Detected (%)",
    cap_at = sites_cap_at,
    x_cap_at = x_cap_at,
    y_transform = y_transform,
    y_transform_sigma = 0.01,
    y_transform_breaks = if (y_transform == "pseudo_log") {
      c(0, 0.01, 0.1, 1)
    } else {
      NULL
    },
    display_as_percent = TRUE,
    overcap_shape = 18,
    prefix_cap_axis_label = FALSE
  )

  p_cvaf + p_sites +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    ) &
    theme(legend.position = "bottom")
}

make_extended_fragmentomics_plot_df <- function(dat) {
  dat %>%
    select(
      Patient, Cohort, Date, Baseline_Date, Clinical_Baseline_Date,
      Weeks_Since_Baseline, Weeks_Since_Clinical_Baseline,
      Mutation_Source_Baseline_Timepoint,
      Mutation_Source_Baseline_Sample_Code,
      Mutation_Source_Baseline_Timepoint_Info,
      Num_days_to_closest_relapse, Relapsed_Binary,
      Proportion.Short,
      TF_ichorCNA = WGS_Tumor_Fraction_Blood_plasma_cfDNA
    ) %>%
    pivot_longer(
      cols = c(Proportion.Short, TF_ichorCNA),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    drop_na(Value) %>%
    flag_followup_relapse_points() %>%
    add_patient_relapse_flag()
}

make_extended_fragmentomics_panel <- function(plot_df, x_cap_at = NA_real_,
                                              y_transform = c("identity", "pseudo_log")) {
  y_transform <- match.arg(y_transform)
  p_short <- make_metric_panel(
    plot_df,
    metric = "Proportion.Short",
    panel_title = "Short-fragment proportion",
    y_label = "Short cfDNA Fragments (%)",
    x_cap_at = x_cap_at,
    y_transform = y_transform,
    y_transform_sigma = 0.01,
    y_transform_breaks = if (y_transform == "pseudo_log") {
      c(0.05, 0.1, 0.2, 0.3)
    } else {
      NULL
    },
    display_as_percent = TRUE,
    overcap_shape = 18,
    prefix_cap_axis_label = FALSE
  )
  p_tf <- make_metric_panel(
    plot_df,
    metric = "TF_ichorCNA",
    panel_title = "cfDNA tumor fraction",
    y_label = "cfDNA Tumor Fraction (%)",
    cap_at = 0.5,
    x_cap_at = x_cap_at,
    y_transform = y_transform,
    y_transform_sigma = 0.01,
    y_transform_breaks = if (y_transform == "pseudo_log") {
      c(0, 0.01, 0.1, 0.5)
    } else {
      NULL
    },
    display_as_percent = TRUE,
    overcap_shape = 18,
    prefix_cap_axis_label = FALSE
  )

  p_short + p_tf +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Fragmentomic tumour-agnostic longitudinal tracking: Relapsed vs Non-relapsed patients",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    ) &
    theme(legend.position = "bottom")
}

bm_plot_df <- make_mutation_plot_df(dat, "BM")
blood_plot_df <- make_mutation_plot_df(dat, "blood")
fragmentomics_plot_df <- make_fragmentomics_plot_df(dat)
ed3a_plot_df <- make_raw_mutation_plot_df(dat, "BM")
ed3b_plot_df <- make_raw_mutation_plot_df(dat, "blood")
ed3c_plot_df <- make_extended_fragmentomics_plot_df(dat)

if (nrow(bm_plot_df) == 0) stop("No BM-evaluable rows available.", call. = FALSE)
if (nrow(blood_plot_df) == 0) stop("No blood-evaluable rows available.", call. = FALSE)
if (nrow(fragmentomics_plot_df) == 0) stop("No fragmentomics-evaluable rows available.", call. = FALSE)
if (nrow(ed3a_plot_df) == 0) stop("No ED3A all-evaluable rows available.", call. = FALSE)
if (nrow(ed3b_plot_df) == 0) stop("No ED3B all-evaluable rows available.", call. = FALSE)
if (nrow(ed3c_plot_df) == 0) stop("No ED3C all-evaluable rows available.", call. = FALSE)

bm_panel <- make_mutation_panel(
  bm_plot_df,
  title = "Baseline BM-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  cvaf_cap_at = 1500,
  sites_cap_at = 1000,
  show_cap_labels = FALSE,
  x_cap_at = 250
)
bm_panel_pseudolog <- make_mutation_panel(
  bm_plot_df,
  title = "Baseline BM-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  cvaf_cap_at = NA_real_,
  sites_cap_at = NA_real_,
  show_cap_labels = FALSE,
  x_cap_at = 250,
  cvaf_y_transform = "pseudo_log",
  sites_y_transform = "pseudo_log"
)
blood_panel <- make_mutation_panel(
  blood_plot_df,
  title = "Baseline cfDNA-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  cvaf_cap_at = NA_real_,
  sites_cap_at = NA_real_,
  x_cap_at = 250
)
blood_panel_pseudolog <- make_mutation_panel(
  blood_plot_df,
  title = "Baseline cfDNA-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  cvaf_cap_at = NA_real_,
  sites_cap_at = NA_real_,
  x_cap_at = 250,
  cvaf_y_transform = "pseudo_log",
  sites_y_transform = "pseudo_log"
)
blood_panel_capped300 <- make_mutation_panel(
  blood_plot_df,
  title = "Baseline cfDNA-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  cvaf_cap_at = 300,
  sites_cap_at = 300,
  show_cap_labels = FALSE,
  x_cap_at = 250
)
fragmentomics_panel <- make_fragmentomics_panel(
  fragmentomics_plot_df,
  x_cap_at = 250
)
fragmentomics_panel_pseudolog <- make_fragmentomics_panel(
  fragmentomics_plot_df,
  x_cap_at = 250,
  y_transform = "pseudo_log"
)
ed3a_panel <- make_raw_mutation_panel(
  ed3a_plot_df,
  title = "Baseline BM-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  cvaf_cap_at = 0.25,
  x_cap_at = 250
)
ed3a_panel_pseudolog <- make_raw_mutation_panel(
  ed3a_plot_df,
  title = "Baseline BM-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  x_cap_at = 250,
  y_transform = "pseudo_log"
)
ed3b_panel <- make_raw_mutation_panel(
  ed3b_plot_df,
  title = "Baseline cfDNA-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  cvaf_cap_at = 0.25,
  x_cap_at = 250
)
ed3b_panel_pseudolog <- make_raw_mutation_panel(
  ed3b_plot_df,
  title = "Baseline cfDNA-informed longitudinal tracking: Relapsed vs Non-relapsed patients",
  x_cap_at = 250,
  y_transform = "pseudo_log"
)
ed3c_panel <- make_extended_fragmentomics_panel(ed3c_plot_df, x_cap_at = 250)
ed3c_panel_pseudolog <- make_extended_fragmentomics_panel(
  ed3c_plot_df,
  x_cap_at = 250,
  y_transform = "pseudo_log"
)

ggsave(
  file.path(output_dir, "Figure_2B_all_evaluable_BM_zscore_longitudinal.png"),
  plot = bm_panel,
  width = 12,
  height = 4,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Figure_2B_all_evaluable_BM_zscore_longitudinal_pseudolog.png"),
  plot = bm_panel_pseudolog,
  width = 12,
  height = 4.5,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal.png"),
  plot = blood_panel,
  width = 12,
  height = 4,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal_pseudolog.png"),
  plot = blood_panel_pseudolog,
  width = 12,
  height = 4.5,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal_capped300.png"),
  plot = blood_panel_capped300,
  width = 12,
  height = 4,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Figure_2D_all_evaluable_fragmentomics_longitudinal.png"),
  plot = fragmentomics_panel,
  width = 12,
  height = 4,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Figure_2D_all_evaluable_fragmentomics_longitudinal_pseudolog.png"),
  plot = fragmentomics_panel_pseudolog,
  width = 12,
  height = 4.5,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal.png"),
  plot = ed3a_panel,
  width = 12,
  height = 4,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal_pseudolog.png"),
  plot = ed3a_panel_pseudolog,
  width = 12,
  height = 4.5,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal.png"),
  plot = ed3b_panel,
  width = 12,
  height = 4,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal_pseudolog.png"),
  plot = ed3b_panel_pseudolog,
  width = 12,
  height = 4.5,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal.png"),
  plot = ed3c_panel,
  width = 12,
  height = 4,
  dpi = 600
)
ggsave(
  file.path(output_dir, "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal_pseudolog.png"),
  plot = ed3c_panel_pseudolog,
  width = 12,
  height = 4.5,
  dpi = 600
)

ms_copy_artifact(
  source_path = file.path(output_dir, "Figure_2B_all_evaluable_BM_zscore_longitudinal.png"),
  artifact_id = "FIG2B",
  role = "all_evaluable_figure_panel_png",
  description = "All-evaluable training/test companion version of Main Figure 2B longitudinal BM mutation z-score panel; red follow-up points use next progression within 180 days after the sample.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Figure_2B_all_evaluable_BM_zscore_longitudinal_pseudolog.png"),
  artifact_id = "FIG2B",
  role = "all_evaluable_pseudolog_figure_panel_png",
  description = "All-evaluable training/test companion version of Main Figure 2B with a signed pseudo-log cumulative VAF z-score axis that retains zero and negative values.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal.png"),
  artifact_id = "FIG2C",
  role = "all_evaluable_figure_panel_png",
  description = "All-evaluable training/test companion version of Main Figure 2C longitudinal cfDNA mutation z-score panel; red follow-up points use next progression within 180 days after the sample.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal_pseudolog.png"),
  artifact_id = "FIG2C",
  role = "all_evaluable_pseudolog_figure_panel_png",
  description = "All-evaluable training/test companion version of Main Figure 2C with signed pseudo-log cumulative VAF and site-detection z-score axes.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal_capped300.png"),
  artifact_id = "FIG2C",
  role = "all_evaluable_capped300_figure_panel_png",
  description = "All-evaluable training/test companion version of Main Figure 2C with cfDNA z-score display capped at 300; red follow-up points use next progression within 180 days after the sample.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Figure_2D_all_evaluable_fragmentomics_longitudinal.png"),
  artifact_id = "FIG2D",
  role = "all_evaluable_figure_panel_png",
  description = "All-evaluable training/test companion version of Main Figure 2D longitudinal fragmentomics panel; red follow-up points use next progression within 180 days after the sample.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Figure_2D_all_evaluable_fragmentomics_longitudinal_pseudolog.png"),
  artifact_id = "FIG2D",
  role = "all_evaluable_pseudolog_figure_panel_png",
  description = "All-evaluable training/test companion version of Main Figure 2D with signed pseudo-log fragment-size and reversed pseudo-log coverage axes.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal.png"),
  artifact_id = "EDFIG3A",
  role = "all_evaluable_figure_panel_png",
  description = "All-evaluable training/test companion version of Extended Data Figure 3A raw BM mutation trajectory panel; red follow-up points use next progression within 180 days after the sample.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal_pseudolog.png"),
  artifact_id = "EDFIG3A",
  role = "all_evaluable_pseudolog_figure_panel_png",
  description = "All-evaluable training/test companion version of Extended Data Figure 3A with pseudo-log raw BM mutation axes.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal.png"),
  artifact_id = "EDFIG3B",
  role = "all_evaluable_figure_panel_png",
  description = "All-evaluable training/test companion version of Extended Data Figure 3B raw cfDNA mutation trajectory panel; red follow-up points use next progression within 180 days after the sample.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal_pseudolog.png"),
  artifact_id = "EDFIG3B",
  role = "all_evaluable_pseudolog_figure_panel_png",
  description = "All-evaluable training/test companion version of Extended Data Figure 3B with pseudo-log raw cfDNA mutation axes.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal.png"),
  artifact_id = "EDFIG3C",
  role = "all_evaluable_figure_panel_png",
  description = "All-evaluable training/test companion version of Extended Data Figure 3C fragmentomics trajectory panel; red follow-up points use next progression within 180 days after the sample.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)
ms_copy_artifact(
  source_path = file.path(output_dir, "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal_pseudolog.png"),
  artifact_id = "EDFIG3C",
  role = "all_evaluable_pseudolog_figure_panel_png",
  description = "All-evaluable training/test companion version of Extended Data Figure 3C with pseudo-log short-fragment and tumour-fraction axes.",
  script_name = "2_4B_Build_all_evaluable_longitudinal_panels.R"
)

write_csv(
  bm_plot_df,
  file.path(output_dir, "Figure_2B_all_evaluable_BM_zscore_longitudinal_source_data.csv")
)
write_csv(
  blood_plot_df,
  file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal_source_data.csv")
)
write_csv(
  fragmentomics_plot_df,
  file.path(output_dir, "Figure_2D_all_evaluable_fragmentomics_longitudinal_source_data.csv")
)
write_csv(
  ed3a_plot_df,
  file.path(output_dir, "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal_source_data.csv")
)
write_csv(
  ed3b_plot_df,
  file.path(output_dir, "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal_source_data.csv")
)
write_csv(
  ed3c_plot_df,
  file.path(output_dir, "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal_source_data.csv")
)

all_evaluable_support_manifest <- tribble(
  ~file_name, ~description,
  "Figure_2B_all_evaluable_BM_zscore_longitudinal.png",
  "All-evaluable support version of the BM longitudinal z-score panel.",
  "Figure_2B_all_evaluable_BM_zscore_longitudinal_pseudolog.png",
  "All-evaluable support version of the BM longitudinal z-score panel with a signed pseudo-log cumulative VAF axis.",
  "Figure_2C_all_evaluable_blood_zscore_longitudinal.png",
  "All-evaluable support version of the cfDNA longitudinal z-score panel.",
  "Figure_2C_all_evaluable_blood_zscore_longitudinal_pseudolog.png",
  "All-evaluable support version of the cfDNA longitudinal z-score panel with signed pseudo-log y-axes.",
  "Figure_2C_all_evaluable_blood_zscore_longitudinal_capped300.png",
  "All-evaluable support version of the cfDNA longitudinal z-score panel with the plotting y-axis capped at 300.",
  "Figure_2D_all_evaluable_fragmentomics_longitudinal.png",
  "All-evaluable support version of the fragmentomics longitudinal z-score panel.",
  "Figure_2D_all_evaluable_fragmentomics_longitudinal_pseudolog.png",
  "All-evaluable support version of the fragmentomics panel with pseudo-log y-axes.",
  "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal.png",
  "All-evaluable support version of the raw BM mutation longitudinal panel.",
  "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal_pseudolog.png",
  "All-evaluable support version of the raw BM mutation panel with pseudo-log y-axes.",
  "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal.png",
  "All-evaluable support version of the raw cfDNA mutation longitudinal panel.",
  "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal_pseudolog.png",
  "All-evaluable support version of the raw cfDNA mutation panel with pseudo-log y-axes.",
  "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal.png",
  "All-evaluable support version of the complementary fragmentomics longitudinal panel.",
  "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal_pseudolog.png",
  "All-evaluable support version of the complementary fragmentomics panel with pseudo-log y-axes."
) %>%
  mutate(
    path = file.path(output_dir, file_name),
    exists = file.exists(path),
    md5 = if_else(exists, map_chr(path, ~ unname(tools::md5sum(.x))), NA_character_),
    note = paste(
      "All-evaluable companion panel staged beside the normal Figure 2 panel",
      "using an all_evaluable role; the normal manuscript panel remains frontline-only."
    )
  )

write_csv(
  all_evaluable_support_manifest,
  file.path(output_dir, "all_evaluable_support_manifest.csv")
)

qc_tbl <- bind_rows(
  bm_plot_df %>% mutate(panel = "Figure_2B_all_evaluable_BM"),
  blood_plot_df %>% mutate(panel = "Figure_2C_all_evaluable_blood"),
  fragmentomics_plot_df %>% mutate(panel = "Figure_2D_all_evaluable_fragmentomics"),
  ed3a_plot_df %>% mutate(panel = "Extended_Data_Figure_3A_all_evaluable_BM_raw"),
  ed3b_plot_df %>% mutate(panel = "Extended_Data_Figure_3B_all_evaluable_blood_raw"),
  ed3c_plot_df %>% mutate(panel = "Extended_Data_Figure_3C_all_evaluable_fragmentomics")
) %>%
  distinct(panel, Metric, Cohort, Patient, Weeks_Since_Baseline, relapse_within_180) %>%
  group_by(panel, Metric, Cohort) %>%
  summarise(
    n_patients = n_distinct(Patient),
    n_timepoints = n(),
    n_red_followup_points = sum(relapse_within_180, na.rm = TRUE),
    n_red_followup_patients = n_distinct(Patient[relapse_within_180]),
    min_week = min(Weeks_Since_Baseline, na.rm = TRUE),
    max_week = max(Weeks_Since_Baseline, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(panel, Metric, Cohort)

write_csv(qc_tbl, file.path(output_dir, "all_evaluable_longitudinal_panel_qc.csv"))
write_csv(
  time_anchor_audit,
  file.path(output_dir, "all_evaluable_longitudinal_time_anchor_audit.csv")
)
capped300_audit <- blood_plot_df %>%
  mutate(capped_at_300 = Value > 300) %>%
  group_by(Metric) %>%
  summarise(
    n_rows = n(),
    n_patients = n_distinct(Patient),
    n_capped_rows = sum(capped_at_300, na.rm = TRUE),
    n_capped_patients = n_distinct(Patient[capped_at_300]),
    max_observed_value = max(Value, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(
  capped300_audit,
  file.path(output_dir, "Figure_2C_all_evaluable_blood_zscore_longitudinal_capped300_audit.csv")
)
pseudolog_audit <- bind_rows(
  bm_plot_df %>% mutate(panel = "Figure_2B", sigma = 1),
  blood_plot_df %>% mutate(panel = "Figure_2C", sigma = 1),
  fragmentomics_plot_df %>% mutate(
    panel = "Figure_2D",
    sigma = if_else(Metric == "FS", 0.05, 0.01)
  )
) %>%
  group_by(panel, Metric, sigma) %>%
  summarise(
    transform = "scales::pseudo_log_trans",
    base = 10,
    n_rows = n(),
    n_patients = n_distinct(Patient),
    n_negative_values = sum(Value < 0, na.rm = TRUE),
    n_zero_values = sum(Value == 0, na.rm = TRUE),
    min_observed_value = min(Value, na.rm = TRUE),
    max_observed_value = max(Value, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(
  pseudolog_audit,
  file.path(output_dir, "Figure_2BCD_all_evaluable_pseudolog_audit.csv")
)
write_csv(
  coverage_fallback_audit,
  file.path(output_dir, "all_evaluable_mean_coverage_fallback_audit.csv")
)

message("Wrote all-evaluable longitudinal panels to: ", normalizePath(output_dir))
print(qc_tbl, n = Inf)
