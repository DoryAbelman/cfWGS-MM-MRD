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
#     additional_all_evaluable_longitudinal_panels/
#       Figure_2B_all_evaluable_BM_zscore_longitudinal.png
#       Figure_2C_all_evaluable_blood_zscore_longitudinal.png
#       Figure_2C_all_evaluable_blood_zscore_longitudinal_capped300.png
#       Figure_2D_all_evaluable_fragmentomics_longitudinal.png
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
  "cohort_assignment_table_updated.rds",
  "baseline_high_quality_patients_updated.csv",
  "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"
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

cohort_df <- readRDS("cohort_assignment_table_updated.rds") %>%
  augment_cohort_assignment_with_spring2026_revision()
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

required_columns(cohort_df, c("Patient", "Cohort"), "cohort_df")
required_columns(
  good_pts,
  c(
    "Patient",
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
    across(all_of(blood_feats), ~ if_else(Patient %in% cfdna_good_pts, .x, NA_real_)),
    WGS_Tumor_Fraction_Blood_plasma_cfDNA = if_else(
      Patient == "IMG-127" & Date == as.Date("2022-08-15"),
      NA_real_,
      WGS_Tumor_Fraction_Blood_plasma_cfDNA
    )
  ) %>%
  left_join(
    baseline_dates %>%
      transmute(Patient, Baseline_Date = as.Date(baseline_date)),
    by = "Patient"
  ) %>%
  mutate(
    Weeks_Since_Baseline = as.numeric(difftime(Date, Baseline_Date, units = "weeks")),
    Weeks_Since_Baseline = case_when(
      Weeks_Since_Baseline >= -2 & Weeks_Since_Baseline < 0 ~ 0,
      TRUE ~ Weeks_Since_Baseline
    )
  )

time_anchor_exclusions <- dat %>%
  filter(is.na(Weeks_Since_Baseline)) %>%
  select(
    Patient, Cohort, Date, Baseline_Date, Num_days_to_closest_relapse,
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

if (nrow(dat) == 0) {
  stop("No evaluable longitudinal rows remain after all-evaluable filters.", call. = FALSE)
}

relapse_labels <- c(
  `FALSE` = "No relapse <=180d",
  `TRUE` = "Relapse <=180d"
)
cohort_linetypes <- c("Frontline" = "solid", "Non-frontline" = "22")
cohort_display_labels <- c("Frontline" = "Training Cohort", "Non-frontline" = "Test Cohort")

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
                              reverse_display = FALSE, show_cap_labels = TRUE) {
  metric_df <- df %>% filter(Metric == metric)
  if (nrow(metric_df) == 0) {
    stop("No rows available for metric: ", metric, call. = FALSE)
  }

  if (is.finite(cap_at)) {
    metric_plot_df <- metric_df %>%
      mutate(
        overcap = Value > cap_at,
        Value_plot = pmin(Value, cap_at),
        label_val = if_else(overcap, as.character(round(Value)), NA_character_)
      )
    seg_df <- metric_df %>%
      arrange(Patient, Weeks_Since_Baseline) %>%
      group_by(Patient) %>%
      mutate(
        x = Weeks_Since_Baseline,
        y = pmin(pmax(Value, 0), cap_at),
        xend = lead(Weeks_Since_Baseline),
        yend = lead(pmin(pmax(Value, 0), cap_at)),
        seg_relapse = lead(relapse_within_180, default = FALSE),
        seg_cohort = first(Cohort)
      ) %>%
      filter(!is.na(x), !is.na(xend), !is.na(y), !is.na(yend)) %>%
      ungroup()

    return(
      ggplot(metric_plot_df, aes(Weeks_Since_Baseline, Value_plot, group = Patient)) +
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
          labeller = labeller(patient_relapse180 = relapse_labels)
        ) +
        scale_color_manual(
          values = c(`FALSE` = "black", `TRUE` = "red"),
          labels = relapse_labels,
          name = NULL
        ) +
        scale_linetype_manual(
          values = cohort_linetypes,
          labels = cohort_display_labels,
          name = "Cohort"
        ) +
        scale_shape_manual(
          values = c(`FALSE` = 16, `TRUE` = 17),
          labels = c(`FALSE` = paste0("<=", cap_at), `TRUE` = paste0(">", cap_at, " (capped)")),
          name = NULL
        ) +
        guides(shape = "none") +
        coord_cartesian(ylim = c(0, cap_at * 1.04), clip = "on") +
        scale_y_continuous(
          breaks = pretty(c(0, cap_at)),
          labels = function(x) ifelse(x == cap_at, paste0(">", cap_at), x)
        ) +
        labs(title = panel_title, x = "Weeks Since Baseline", y = y_label) +
        theme_classic(base_size = 11) +
        theme(strip.background = element_rect(fill = "grey95", colour = NA))
    )
  }

  seg_df <- make_segment_df(metric_df)
  panel <- ggplot(metric_df, aes(Weeks_Since_Baseline, Value, group = Patient)) +
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
      labeller = labeller(patient_relapse180 = relapse_labels)
    ) +
    scale_color_manual(
      values = c(`FALSE` = "black", `TRUE` = "red"),
      labels = relapse_labels,
      name = NULL
    ) +
    scale_linetype_manual(
      values = cohort_linetypes,
      labels = cohort_display_labels,
      name = "Cohort"
    ) +
    labs(title = panel_title, x = "Weeks Since Baseline", y = y_label) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA))

  if (isTRUE(reverse_display)) {
    panel <- panel + scale_y_reverse()
  }
  panel
}

make_mutation_plot_df <- function(dat, assay) {
  if (assay == "BM") {
    out <- dat %>%
      select(
        Patient, Cohort, Weeks_Since_Baseline, Num_days_to_closest_relapse,
        cVAF_z = z_score_detection_rate_BM,
        sites_z = zscore_BM
      )
  } else if (assay == "blood") {
    out <- dat %>%
      select(
        Patient, Cohort, Weeks_Since_Baseline, Num_days_to_closest_relapse,
        cVAF_z = z_score_detection_rate_blood,
        sites_z = zscore_blood
      )
  } else {
    stop("Unknown assay: ", assay, call. = FALSE)
  }

  out %>%
    mutate(
      relapse_within_180 = if_else(
        Num_days_to_closest_relapse <= 180,
        TRUE,
        FALSE,
        missing = FALSE
      )
    ) %>%
    pivot_longer(
      cols = c(cVAF_z, sites_z),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    drop_na(Value) %>%
    add_patient_relapse_flag()
}

make_mutation_panel <- function(plot_df, title, cvaf_cap_at = NA_real_,
                                sites_cap_at = NA_real_,
                                show_cap_labels = TRUE) {
  p_cvaf <- make_metric_panel(
    plot_df,
    metric = "cVAF_z",
    panel_title = "Cumulative VAF Z-score",
    y_label = "Cumulative VAF (Z)",
    cap_at = cvaf_cap_at,
    show_cap_labels = show_cap_labels
  )
  p_sites <- make_metric_panel(
    plot_df,
    metric = "sites_z",
    panel_title = "Proportion of Sites Detected Z-score",
    y_label = "Prop. Mutant Sites Detected (Z)",
    cap_at = sites_cap_at,
    show_cap_labels = show_cap_labels
  )

  p_cvaf + p_sites +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    ) &
    theme(legend.position = "bottom")
}

make_fragmentomics_plot_df <- function(dat) {
  dat %>%
    select(
      Patient, Cohort, Weeks_Since_Baseline, Num_days_to_closest_relapse,
      FS, Mean.Coverage
    ) %>%
    mutate(
      relapse_within_180 = if_else(
        Num_days_to_closest_relapse <= 180,
        TRUE,
        FALSE,
        missing = FALSE
      )
    ) %>%
    pivot_longer(
      cols = c(FS, Mean.Coverage),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    drop_na(Value) %>%
    add_patient_relapse_flag()
}

make_fragmentomics_panel <- function(plot_df) {
  p_fs <- make_metric_panel(
    plot_df,
    metric = "FS",
    panel_title = "Fragment-size score",
    y_label = "Fragment-size score"
  )
  p_coverage <- make_metric_panel(
    plot_df,
    metric = "Mean.Coverage",
    panel_title = "cfDNA coverage at MM active regulatory sites",
    y_label = "Mean cfDNA coverage (MM regs)",
    reverse_display = TRUE
  )

  p_fs + p_coverage +
    plot_annotation(
      title = "Longitudinal trajectories of fragmentomic features\nAll evaluable patients",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    ) &
    theme(legend.position = "bottom")
}

bm_plot_df <- make_mutation_plot_df(dat, "BM")
blood_plot_df <- make_mutation_plot_df(dat, "blood")
fragmentomics_plot_df <- make_fragmentomics_plot_df(dat)

if (nrow(bm_plot_df) == 0) stop("No BM-evaluable rows available.", call. = FALSE)
if (nrow(blood_plot_df) == 0) stop("No blood-evaluable rows available.", call. = FALSE)
if (nrow(fragmentomics_plot_df) == 0) stop("No fragmentomics-evaluable rows available.", call. = FALSE)

bm_panel <- make_mutation_panel(
  bm_plot_df,
  title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles\nAll evaluable patients",
  cvaf_cap_at = 2500
)
blood_panel <- make_mutation_panel(
  blood_plot_df,
  title = "Longitudinal trajectories of MRD metrics from baseline PB cfDNA mutation profiles\nAll evaluable patients",
  cvaf_cap_at = NA_real_,
  sites_cap_at = NA_real_
)
blood_panel_capped300 <- make_mutation_panel(
  blood_plot_df,
  title = "Longitudinal trajectories of MRD metrics from baseline PB cfDNA mutation profiles\nAll evaluable patients",
  cvaf_cap_at = 300,
  sites_cap_at = 300,
  show_cap_labels = FALSE
)
fragmentomics_panel <- make_fragmentomics_panel(fragmentomics_plot_df)

ggsave(
  file.path(output_dir, "Figure_2B_all_evaluable_BM_zscore_longitudinal.png"),
  plot = bm_panel,
  width = 12,
  height = 4,
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

figure_panel_targets <- tribble(
  ~source_file, ~panel_subdir, ~target_file,
  "Figure_2B_all_evaluable_BM_zscore_longitudinal.png",
  "Figure_2B",
  "Figure_2B_ALL_EVALUABLE__figure_panel_png__BM_zscore_longitudinal.png",
  "Figure_2C_all_evaluable_blood_zscore_longitudinal.png",
  "Figure_2C",
  "Figure_2C_ALL_EVALUABLE__figure_panel_png__blood_zscore_longitudinal.png",
  "Figure_2C_all_evaluable_blood_zscore_longitudinal_capped300.png",
  "Figure_2C",
  "Figure_2C_ALL_EVALUABLE__figure_panel_png__blood_zscore_longitudinal_capped300.png",
  "Figure_2D_all_evaluable_fragmentomics_longitudinal.png",
  "Figure_2D",
  "Figure_2D_ALL_EVALUABLE__figure_panel_png__fragmentomics_longitudinal.png"
)

walk2(
  file.path(output_dir, figure_panel_targets$source_file),
  file.path(figure2_dir, figure_panel_targets$panel_subdir, figure_panel_targets$target_file),
  function(source_path, target_path) {
    dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source_path, target_path, overwrite = TRUE)) {
      stop("Failed to copy all-evaluable panel to: ", target_path, call. = FALSE)
    }
  }
)

qc_tbl <- bind_rows(
  bm_plot_df %>% mutate(panel = "Figure_2B_all_evaluable_BM"),
  blood_plot_df %>% mutate(panel = "Figure_2C_all_evaluable_blood"),
  fragmentomics_plot_df %>% mutate(panel = "Figure_2D_all_evaluable_fragmentomics")
) %>%
  distinct(panel, Metric, Cohort, Patient, Weeks_Since_Baseline) %>%
  group_by(panel, Metric, Cohort) %>%
  summarise(
    n_patients = n_distinct(Patient),
    n_timepoints = n(),
    min_week = min(Weeks_Since_Baseline, na.rm = TRUE),
    max_week = max(Weeks_Since_Baseline, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(panel, Metric, Cohort)

write_csv(qc_tbl, file.path(output_dir, "all_evaluable_longitudinal_panel_qc.csv"))
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
write_csv(
  coverage_fallback_audit,
  file.path(output_dir, "all_evaluable_mean_coverage_fallback_audit.csv")
)

message("Wrote all-evaluable longitudinal panels to: ", normalizePath(output_dir))
print(qc_tbl, n = Inf)
