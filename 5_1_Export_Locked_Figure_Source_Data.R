#!/usr/bin/env Rscript

# Export panel-level source-data CSVs for the current manuscript figures.
#
# Pipeline role
#   This is the first step of the figure-source-data packaging workflow:
#     5_1 export one de-identified CSV per figure panel
#     5_2 assemble those CSVs into the main- and extended-data workbooks
#     5_4 validate workbook structure and content
#   Script 5_3 runs those three steps in order. This script does not fit models,
#   select thresholds, or redraw the manuscript figures.
#
# Manuscript coverage and unit of analysis
#   - Covers 74 sheets representing every panel in Figures 1-4 and Extended
#     Data Figures 1-10.
#   - The unit varies by panel (for example: patient, sample, genomic feature,
#     ROC coordinate, Kaplan-Meier record, or dilution-series measurement).
#     `panel_schema_contract` below records the fields required to redraw each
#     panel, and the manifest records the row and column counts.
#   - Figure 3A/4A and Extended Data Figure 5A/7A/9A/9B use the retained
#     50-repeat patient-grouped nested-CV outputs. This exporter does not rerun
#     the repeated nested CV.
#
# Inputs
#   - Final PDFs in Manuscript_Exports/{01_main_figures,02_extended_data_figures}
#   - Panel source CSVs in final_manuscript_objects/generated/figure_components
#   - Retained validation metrics, models, call tables, and figure-specific
#     source objects at the project root
#   - id_map.rds, a complete one-to-one Patient/New_ID publication mapping
#
# Outputs
#   - Output_tables_2025/Figure_Source_Data/panels/*.csv
#   - Output_tables_2025/Figure_Source_Data/source_data_workbook_manifest.csv
#   - Output_tables_2025/Figure_Source_Data/source_data_audit.csv
#   - Output_tables_2025/Figure_Source_Data/panel_schema_contract.csv
#
# Additional staged copies
#   - Rewrites the Extended Data Figure 2D component source CSV from the
#     retained Sequenza segments and FISH-probe interval.
#   - Copies the exported Figure 2E CSV into its Figure_2E manuscript-object
#     directory for downstream packaging.
#
# How to run
#   Run from the project root:
#     Rscript Scripts_2025/Final_Scripts/5_1_Export_Locked_Figure_Source_Data.R
#
# This script uses the final panel inventory and revision-inclusive analysis
# tables. If a figure is changed, update the expectations here in the same
# commit. Silent denominator drift is treated as a provenance error.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readxl)
})

.publication_export_helper <- file.path(
  "Scripts_2025", "Final_Scripts", "publication_export_helpers.R"
)
if (!file.exists(.publication_export_helper)) {
  .publication_export_helper <- "publication_export_helpers.R"
}
source(.publication_export_helper)
rm(.publication_export_helper)

project_root <- normalizePath(".", mustWork = TRUE)
canonical_root <- file.path(project_root, "Scripts_2025", "Final_Scripts")
component_root <- file.path(
  canonical_root, "final_manuscript_objects", "generated", "figure_components"
)
output_root <- file.path(project_root, "Output_tables_2025", "Figure_Source_Data")
panel_output_root <- file.path(output_root, "panels")
dir.create(panel_output_root, recursive = TRUE, showWarnings = FALSE)

id_map_path <- file.path(project_root, "id_map.rds")
stopifnot(file.exists(id_map_path))
publication_id_map <- readRDS(id_map_path) %>%
  transmute(
    raw_id = as.character(.data$Patient),
    public_id = as.character(.data$New_ID)
  )
if (
  nrow(publication_id_map) == 0L || anyNA(publication_id_map) ||
    anyDuplicated(publication_id_map$raw_id) ||
    anyDuplicated(publication_id_map$public_id)
) {
  stop("id_map.rds must contain a complete one-to-one Patient/New_ID mapping.", call. = FALSE)
}

# This one ED1 oncoprint specimen is outside the 71-patient id_map.rds cohort.
# Reuse its existing Code Ocean public specimen alias rather than inventing a
# new patient identifier or exposing the original SPORE code.
publication_sample_overrides <- c("SPORE_0008_Baseline" = "S0112")

deidentify_character_vector <- function(x) {
  x <- as.character(x)
  for (raw_sample in names(publication_sample_overrides)) {
    x <- gsub(raw_sample, publication_sample_overrides[[raw_sample]], x, fixed = TRUE)
  }
  replacement_order <- order(nchar(publication_id_map$raw_id), decreasing = TRUE)
  placeholders <- paste0("__RAW_PATIENT_", seq_len(nrow(publication_id_map)), "__")
  public_placeholders <- paste0("__EXISTING_PUBLIC_PATIENT_", seq_len(nrow(publication_id_map)), "__")
  public_order <- order(nchar(publication_id_map$public_id), decreasing = TRUE)
  for (i in public_order) {
    x <- gsub(
      paste0(
        "(?<![A-Za-z0-9])", publication_id_map$public_id[[i]],
        "(?![A-Za-z0-9])"
      ),
      public_placeholders[[i]], x, perl = TRUE
    )
  }
  for (i in replacement_order) {
    x <- gsub(
      publication_id_map$raw_id[[i]],
      placeholders[[i]],
      x,
      fixed = TRUE
    )
  }
  for (i in seq_len(nrow(publication_id_map))) {
    x <- gsub(placeholders[[i]], publication_id_map$public_id[[i]], x, fixed = TRUE)
  }
  for (i in seq_len(nrow(publication_id_map))) {
    x <- gsub(public_placeholders[[i]], publication_id_map$public_id[[i]], x, fixed = TRUE)
  }
  x
}

contains_original_id <- function(values, raw_id) {
  any(grepl(
    paste0("(?<![A-Za-z0-9])", raw_id, "(?![A-Za-z0-9])"),
    values, perl = TRUE
  ), na.rm = TRUE)
}

deidentify_publication_data <- function(data, panel_key) {
  data[] <- lapply(data, function(column) {
    if (is.factor(column) || is.character(column)) {
      deidentify_character_vector(column)
    } else if (is.logical(column)) {
      as.integer(column)
    } else {
      column
    }
  })
  remaining_raw <- unique(unlist(lapply(
    data[vapply(data, is.character, logical(1))],
    function(column) publication_id_map$raw_id[vapply(
      publication_id_map$raw_id,
      function(raw_id) contains_original_id(column, raw_id),
      logical(1)
    )]
  ), use.names = FALSE))
  if (length(remaining_raw)) {
    stop(
      "Original patient identifier(s) remained in publication panel ",
      panel_key, ": ", paste(remaining_raw, collapse = ", "),
      call. = FALSE
    )
  }
  data
}

drop_empty_columns <- function(data) {
  keep <- vapply(data, function(column) {
    any(!is.na(column) & (if (is.character(column)) nzchar(column) else TRUE))
  }, logical(1))
  data[, keep, drop = FALSE]
}

select_existing <- function(data, columns) {
  data %>% select(any_of(columns))
}

minimalize_publication_panel <- function(data, panel_key) {
  # Administrative provenance remains in the external manifest/audit files,
  # not in the journal-facing panel sheets.
  provenance_columns <- c(
    "figure_panel", "locked_figure", "source_table", "source_note",
    "source_call_path", "record_type", "analysis_value", "value_scale",
    "Bam", "vcf_clean"
  )
  # ED7E combines curve-coordinate rows and fixed operating-point rows in one
  # minimal aggregate table; record_type is plotting data for this panel, not
  # administrative provenance.
  if (panel_key == "ED7E") {
    provenance_columns <- setdiff(provenance_columns, "record_type")
  }
  if (panel_key == "Fig2E") {
    if (!"record_type" %in% names(data)) {
      stop("Figure 2E source data lacks record_type.", call. = FALSE)
    }
    data <- data %>%
      filter(.data$record_type == "patient_terminal_observation")
  }

  data <- data %>%
    select(-any_of(provenance_columns)) %>%
    select(-matches("(^|[._ -])date($|[._ -])", ignore.case = TRUE)) %>%
    select(-matches("(^|_)source_(file|path)$|(^|_)file_path$", ignore.case = TRUE))

  if (panel_key == "Fig1A") {
    # The swim-plot source table carries only the event intervals. The locked
    # panel also draws a cohort band and a four-assay MRD-at-1yr band
    # (2_1_Part2_Cohort_Swim_Plot.R:2807-2818 builds `df_mrd_long` from
    # `dat_1yr`, and :3246 draws the cohort band). Those per-patient values are
    # attached here so the annotation bands are reproducible from the public
    # layer; without them only the timeline itself could be redrawn.
    # `chemo_group_simple` stays derived from `details` at plot time, exactly as
    # 2_1_Part2...R:2564-2575 does.
    public_calls <- clinical_calls %>%
      mutate(Patient = deidentify_character_vector(.data$Patient))
    # Cohort is a patient-level attribute and must not be restricted to the
    # 1yr-maintenance rows, otherwise every test-cohort patient loses its band.
    cohort_annotation <- public_calls %>%
      filter(!is.na(.data$Cohort), nzchar(as.character(.data$Cohort))) %>%
      transmute(Patient, Cohort = as.character(.data$Cohort)) %>%
      distinct(.data$Patient, .keep_all = TRUE)
    mrd_annotation <- public_calls %>%
      filter(.data$timepoint_info == "1yr maintenance") %>%
      transmute(
        Patient,
        MFC = as.integer(.data$Flow_Binary),
        clonoSEQ = as.integer(.data$Adaptive_Binary),
        cfWGS_BM = as.integer(.data$BM_zscore_only_detection_rate_call),
        cfWGS_blood = as.integer(.data$Blood_zscore_only_sites_call)
      ) %>%
      distinct(.data$Patient, .keep_all = TRUE)
    # 2_1_Part2_Cohort_Swim_Plot.R:2578 marks three previously treated patients
    # with an asterisk. The list is written there as original study codes, so it
    # is mapped through the same release namespace here.
    previously_treated_public <- deidentify_character_vector(
      c("EK-09", "IMG-098", "SPORE_0009")
    )
    data <- data %>%
      left_join(cohort_annotation, by = "Patient") %>%
      left_join(mrd_annotation, by = "Patient") %>%
      mutate(previously_treated = .data$Patient %in% previously_treated_public) %>%
      select_existing(c(
        "Patient", "event", "start_day_from_baseline", "end_day_from_baseline",
        "details", "Cohort", "MFC", "clonoSEQ", "cfWGS_BM", "cfWGS_blood",
        "previously_treated"
      ))
  } else if (panel_key == "Fig1B") {
    data <- data %>%
      filter(!is.na(.data$box_id), nzchar(as.character(.data$box_id))) %>%
      select_existing(c(
        "figure_cohort", "box_id", "box_label", "n", "definition"
      )) %>%
      distinct()
  } else if (panel_key == "Fig2E") {
    fig2e_stats_path <- file.path(
      canonical_root, "final_manuscript_objects",
      "additional_all_evaluable_longitudinal_panels",
      "reviewer_alternative_terminal_summary_v2",
      "terminal_sample_group_statistics.csv"
    )
    stop_if_missing(fig2e_stats_path, "Figure 2E global 12-test statistics")
    fig2e_stats <- read_csv(fig2e_stats_path, show_col_types = FALSE) %>%
      filter(.data$panel_family == "Figure 2") %>%
      select("panel", "Metric", "wilcoxon_q_bh_global_12_tests")
    data <- data %>%
      select(-any_of("wilcoxon_q_bh_global_12_tests")) %>%
      left_join(fig2e_stats, by = c("panel", "Metric"))
    data <- data %>%
      select_existing(c(
        "Patient", "Sample_Code", "Timepoint", "timepoint_info", "Cohort",
        "Weeks_Since_Baseline", "Num_days_to_closest_relapse",
        "Relapsed_Binary", "relapse_within_180", "Metric", "Value",
        "panel", "terminal_group", "wilcoxon_q_bh_global_12_tests"
      ))
  } else if (panel_key == "ED3E") {
    # The compact canonical ED3E component prints a metric-specific annotation
    # assembled from the displayed condition percentages plus a paired
    # leave-one-out comparison. The component source CSV historically omitted
    # that annotation table, making an exact public redraw impossible.
    comparison_path <- file.path(
      project_root, "Results_MRDetect", "Healthy_control_platform_calibration",
      "MRDetect_leave_one_out_platform_comparison.csv"
    )
    stop_if_missing(comparison_path, "Extended Data Figure 3E annotation")
    primary_metrics <- c(
      "detection_rate_as_reads_detected_over_reads_checked",
      "sites_detection_rate"
    )
    ed3e_annotations <- data %>%
      filter(
        .data$metric %in% primary_metrics,
        .data$condition %in% c("XPlus\nmismatched", "XPlus\nmatched"),
        !is.na(.data$pct_z_gt_2)
      ) %>%
      select("metric", "metric_label", "condition", "pct_z_gt_2") %>%
      distinct() %>%
      pivot_wider(names_from = "condition", values_from = "pct_z_gt_2") %>%
      left_join(
        read_csv(comparison_path, show_col_types = FALSE) %>%
          filter(.data$metric %in% primary_metrics) %>%
          select("metric", "paired_wilcoxon_abs_z_p"),
        by = "metric"
      ) %>%
      mutate(
        annotation = sprintf(
          "XPlus z > 2: %.0f%% vs %.0f%%\nMatched |z| vs 6000: P = %.2g",
          .data$`XPlus\nmismatched`, .data$`XPlus\nmatched`,
          .data$paired_wilcoxon_abs_z_p
        )
      ) %>%
      select("metric", "annotation")
    if (nrow(ed3e_annotations) != 2L || anyNA(ed3e_annotations$annotation)) {
      stop("Extended Data Figure 3E annotation export is incomplete.", call. = FALSE)
    }
    data <- data %>% left_join(ed3e_annotations, by = "metric")
  } else if (panel_key %in% c("Fig3A", "Fig4A")) {
    roc_cohort <- if (panel_key == "Fig3A") "BM" else "Blood"
    data <- data %>%
      filter(!is.na(.data$cohort), .data$cohort == roc_cohort) %>%
      select_existing(c(
        "cohort", "model", "display_name", "legend_label", "colour", "order",
        "selected", "fpr", "mean_tpr", "repeat_pooled_auc_mean",
        "auc_cluster_boot_q025", "auc_cluster_boot_q975",
        "outer_repeats", "outer_folds", "inner_repeats", "inner_folds"
      ))
    if (!nrow(data)) {
      stop("No mean-ROC rows for cohort ", roc_cohort, " in panel ", panel_key, ".", call. = FALSE)
    }
  } else if (panel_key == "ED9A") {
    data <- data %>%
      filter(!is.na(.data$cohort), .data$cohort == "FullFrag") %>%
      select_existing(c(
        "cohort", "model", "label", "short_label", "legend_label", "color",
        "fpr", "mean_tpr", "tpr_q025", "tpr_q975", "repeat_pooled_auc_mean",
        "auc_cluster_boot_q025", "auc_cluster_boot_q975"
      ))
    if (!nrow(data)) stop("No fragmentomics mean-ROC rows for ED9A.", call. = FALSE)
  } else if (panel_key %in% c("ED5A", "ED7A", "ED9B")) {
    operating_cohort <- switch(panel_key, ED5A = "BM", ED7A = "Blood", ED9B = "FullFrag")
    data <- data %>%
      filter(!is.na(.data$cohort), .data$cohort == operating_cohort) %>%
      select_existing(c(
        "cohort", "model", "label", "short_label", "legend_label", "color",
        "sensitivity_mean", "sensitivity_sd", "specificity_mean",
        "specificity_sd", "accuracy_mean", "n_outer_folds",
        "repeat_pooled_auc_mean"
      ))
    if (!nrow(data)) {
      stop("No operating-point rows for cohort ", operating_cohort, " in panel ", panel_key, ".", call. = FALSE)
    }
  } else if (panel_key == "Fig3F") {
    data <- select_existing(data, c(
      "Patient", "Timepoint", "timepoint_info", "Time_to_event_months",
      "Relapsed_Binary", "BM_zscore_only_detection_rate_call",
      "BM_zscore_only_detection_rate_prob", "Group"
    ))
  } else if (panel_key == "Fig4E") {
    data <- select_existing(data, c(
      "Patient", "Timepoint", "timepoint_info", "Time_to_event_months",
      "Relapsed_Binary", "Blood_zscore_only_sites_call",
      "Blood_zscore_only_sites_prob", "Group"
    ))
  } else if (panel_key %in% c("ED6C", "ED6D", "ED6E", "ED6F", "ED8C")) {
    # The older component CSVs carried derived `assay_var` / `assay_label` /
    # `timepoint_label` metadata columns. The current dated KM exports (now
    # selected by the overrides above) carry the raw assay call columns
    # instead, so the assay variable and its labels are resolved from the
    # canonical `techs` / `tp_labels` vectors in 4_1_Survival_Analysis.R.
    #
    # The panel -> assay assignment below is not inferred from the rendered
    # panels; it is taken from the locked component source-data filenames under
    # final_manuscript_objects/generated/figure_components/, which name the
    # assay and landmark explicitly:
    #   ED6C Extended_Data_Figure_6C_MFC_one_year_maintenance_KM_source_data.csv
    #   ED6D Extended_Data_Figure_6D_clonoSEQ_one_year_maintenance_KM_source_data.csv
    #   ED6E Extended_Data_Figure_6E_BM_cVAF_model_post_ASCT_KM_source_data.csv
    #   ED6F Extended_Data_Figure_6F_MFC_post_ASCT_KM_source_data.csv
    #   ED8C Extended_Data_Figure_8C_blood_sites_model_post_ASCT_KM_source_data.csv
    # (Figure_3F_BM_cVAF_model_one_year_maintenance and
    #  Figure_4E_blood_sites_model_one_year_maintenance confirm the two main
    #  panels, which select their assay columns by name and need no mapping.)
    km_panel_assay <- c(
      ED6C = "Flow_Binary", ED6D = "Adaptive_Binary",
      ED6E = "BM_zscore_only_detection_rate_call", ED6F = "Flow_Binary",
      ED8C = "Blood_zscore_only_sites_call"
    )
    km_panel_assay_label <- c(
      ED6C = "MFC", ED6D = "clonoSEQ",
      ED6E = "cfWGS of BM-Derived Mutations (cVAF Model)", ED6F = "MFC",
      ED8C = "cfWGS of cfDNA-Derived Mutations (Sites Model)"
    )
    km_panel_timepoint_label <- c(
      ED6C = "One-Year Maintenance", ED6D = "One-Year Maintenance",
      ED6E = "Post‑ASCT", ED6F = "Post‑ASCT", ED8C = "Post‑ASCT"
    )
    assay_columns <- if ("assay_var" %in% names(data)) {
      unique(na.omit(as.character(data$assay_var)))
    } else {
      unname(km_panel_assay[[panel_key]])
    }
    if (length(assay_columns) != 1L || !assay_columns %in% names(data)) {
      stop("Could not resolve the plotted assay call for panel ", panel_key, ".", call. = FALSE)
    }
    data$Assay_Call <- as.integer(data[[assay_columns]])
    if (!"assay_label" %in% names(data)) {
      data$assay_label <- unname(km_panel_assay_label[[panel_key]])
    }
    if (!"timepoint_label" %in% names(data)) {
      data$timepoint_label <- unname(km_panel_timepoint_label[[panel_key]])
    }
    data <- select_existing(data, c(
      "Patient", "Timepoint", "timepoint_info", "Time_to_event_months",
      "Relapsed_Binary", "Assay_Call", "assay_label", "Group", "timepoint_label"
    ))
  } else if (panel_key %in% c("ED6G", "ED6H")) {
    data <- select_existing(data, c(
      "Patient", "Timepoint", "timepoint_info", "Time_to_event_months",
      "Relapsed_Binary", "EasyM_reference_threshold_binary", "EasyM_value", "Group"
    ))
  # ED6J/ED6K and ED8E/ED8F plot the model *probability* on one axis
  # (4_1_Survival_Analysis.R:4128 and :4286) with the call used only to derive
  # `mrd_status`. The `_prob` column was previously omitted from these four
  # exports, so the panels could not be redrawn from the public layer at all.
  # It is added here; the same probability column is already published for
  # Fig3F/Fig4E, so this exposes nothing new.
  } else if (panel_key %in% c("ED6J", "ED6K")) {
    data <- select_existing(data, c(
      "Patient", "Timepoint", "timepoint_info", "Cohort",
      "days_before_event", "months_before_event", "days_plot",
      "Relapsed_Binary", "progress_status", "mrd_status", "cohort_label",
      "Flow_Binary", "Adaptive_Binary", "Rapid_Novor_Binary", "PET_Binary",
      "BM_zscore_only_detection_rate_call",
      "BM_zscore_only_detection_rate_prob"
    ))
  } else if (panel_key %in% c("ED8E", "ED8F")) {
    data <- select_existing(data, c(
      "Patient", "Timepoint", "timepoint_info", "Cohort",
      "days_before_event", "months_before_event", "days_plot",
      "Relapsed_Binary", "progress_status", "mrd_status", "cohort_label",
      "Flow_Binary", "Adaptive_Binary", "Rapid_Novor_Binary", "PET_Binary",
      "Blood_zscore_only_sites_call",
      "Blood_zscore_only_sites_prob"
    ))
  }

  data <- drop_empty_columns(data)
  deidentify_publication_data(data, panel_key)
}

stop_if_missing <- function(paths, context) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      context, " is missing required file(s):\n  ",
      paste(missing, collapse = "\n  "),
      call. = FALSE
    )
  }
}

safe_sheet_name <- function(x) {
  x <- gsub("[:\\\\/?*\\[\\]]", "_", x)
  if (nchar(x) > 31L) {
    stop("Excel sheet name exceeds 31 characters: ", x, call. = FALSE)
  }
  x
}

locked_figures <- tribble(
  ~workbook, ~figure, ~panels,
  "main", "Figure 1", c("A", "B"),
  "main", "Figure 2", c("A", "B", "C", "D", "E"),
  "main", "Figure 3", c("A", "B", "C", "D", "E", "F"),
  "main", "Figure 4", c("A", "B", "C", "D", "E"),
  "extended", "Extended Data Figure 1", "all",
  "extended", "Extended Data Figure 2", LETTERS[1:7],
  "extended", "Extended Data Figure 3", LETTERS[1:7],
  "extended", "Extended Data Figure 4", "all",
  "extended", "Extended Data Figure 5", LETTERS[1:7],
  "extended", "Extended Data Figure 6", LETTERS[1:11],
  "extended", "Extended Data Figure 7", LETTERS[1:8],
  "extended", "Extended Data Figure 8", LETTERS[1:6],
  "extended", "Extended Data Figure 9", LETTERS[1:6],
  "extended", "Extended Data Figure 10", LETTERS[1:2]
) %>%
  unnest_longer(panels, values_to = "panel") %>%
  mutate(
    figure_number = as.integer(sub(".*?([0-9]+)$", "\\1", figure)),
    panel_key = if_else(panel == "all", "all", panel),
    sheet_name = if_else(
      workbook == "main",
      paste0("Fig", figure_number, panel),
      paste0("ED", figure_number, panel)
    ),
    component_figure = if_else(
      workbook == "main",
      paste0("Figure_", figure_number),
      paste0("Extended_Data_Figure_", figure_number)
    ),
    component_panel = paste0("panel_", panel_key),
    locked_pdf = case_when(
      sheet_name == "Fig2E" ~ file.path(
        canonical_root, "final_manuscript_objects", "01_main_figures",
        "Figure_2", "Figure_2E", "Fig_2E.pdf"
      ),
      workbook == "main" ~
      file.path(
        project_root, "Manuscript_Exports", "01_main_figures",
        component_figure, "final_artifacts", paste0(component_figure, ".pdf")
      ),
      TRUE ~ file.path(
        project_root, "Manuscript_Exports", "02_extended_data_figures",
        component_figure, "final_artifacts", paste0(component_figure, ".pdf")
      )
    )
  )

stop_if_missing(unique(locked_figures$locked_pdf), "Locked figure inventory")
stopifnot(!anyDuplicated(locked_figures$sheet_name))
walk(locked_figures$sheet_name, safe_sheet_name)

read_csv_character_safe <- function(path) {
  read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(source_table = basename(path), .before = 1)
}

combine_source_files <- function(paths) {
  stop_if_missing(paths, "Panel source-data export")
  if (length(paths) == 0L) {
    stop("Panel source-data export received no source files.", call. = FALSE)
  }
  bind_rows(map(paths, read_csv_character_safe))
}

build_figure_2e <- function() {
  source_root <- file.path(
    canonical_root, "final_manuscript_objects",
    "additional_all_evaluable_longitudinal_panels",
    "reviewer_alternative_terminal_summary_v2"
  )
  point_path <- file.path(
    source_root, "terminal_sample_patient_level_source_data.csv"
  )
  statistic_path <- file.path(
    source_root, "terminal_sample_group_statistics.csv"
  )
  stop_if_missing(
    c(point_path, statistic_path),
    "Figure 2E terminal-sample source data"
  )

  patient_rows <- read_csv(
    point_path, show_col_types = FALSE, progress = FALSE
  ) %>%
    filter(panel_family == "Figure 2") %>%
    mutate(
      record_type = "patient_terminal_observation",
      source_table = basename(point_path),
      .before = 1
    )
  statistic_rows <- read_csv(
    statistic_path, show_col_types = FALSE, progress = FALSE
  ) %>%
    filter(panel_family == "Figure 2") %>%
    mutate(
      record_type = "group_summary_and_two_sided_wilcoxon_test",
      source_table = basename(statistic_path),
      .before = 1
    )

  patient_panel_metrics <- patient_rows %>% distinct(panel, Metric) %>% nrow()
  if (nrow(patient_rows) != 320L || patient_panel_metrics != 6L) {
    stop(
      "Figure 2E patient-level source audit failed: expected 320 rows across ",
      "six panel-metric combinations; observed ", nrow(patient_rows),
      " rows across ", patient_panel_metrics, " panel-metric combinations.",
      call. = FALSE
    )
  }
  statistic_panel_metrics <- statistic_rows %>% distinct(panel, Metric) %>% nrow()
  if (nrow(statistic_rows) != 6L || statistic_panel_metrics != 6L) {
    stop(
      "Figure 2E statistics audit failed: expected one row for each of six ",
      "metrics.", call. = FALSE
    )
  }
  bind_rows(patient_rows, statistic_rows)
}

build_ed9c <- function() {
  youden_path <- file.path(
    project_root, "Output_tables_2025",
    "cfWGS_model_metrics_youden_threshold2.csv"
  )
  fixed_path <- file.path(
    project_root, "Output_tables_2025",
    "cfWGS_model_metrics_fixed_95sens3_Feb2026.csv"
  )
  stop_if_missing(c(youden_path, fixed_path), "Extended Data Figure 9C metrics")
  model_labels <- c(
    Fragmentomics_full_prob = "Combined model",
    Fragmentomics_mean_coverage_only_prob = "Coverage model",
    Fragmentomics_prop_short_only_prob = "Prop. short model"
  )
  fixed_model_labels <- c(
    Fragmentomics_full_Full_prob = "Combined model",
    Fragmentomics_mean_coverage_only_Full_prob = "Coverage model",
    Fragmentomics_prop_short_only_Full_prob = "Prop. short model"
  )
  youden <- read_csv(youden_path, show_col_types = FALSE) %>%
    filter(.data$model %in% names(model_labels)) %>%
    transmute(
      combo = unname(model_labels[.data$model]),
      Acc = .data$bal_accuracy, F1 = .data$f1,
      Sens = .data$sensitivity, Spec = .data$specificity
    ) %>%
    pivot_longer(cols = -"combo", names_to = "Metric", values_to = "mean")
  fixed <- read_csv(fixed_path, show_col_types = FALSE) %>%
    filter(.data$model %in% names(fixed_model_labels)) %>%
    transmute(
      combo = unname(fixed_model_labels[.data$model]),
      Acc = .data$bal_accuracy, F1 = .data$f1,
      Sens = .data$sensitivity, Spec = .data$specificity
    ) %>%
    pivot_longer(cols = -"combo", names_to = "Metric", values_to = "fixed")
  left_join(youden, fixed, by = c("combo", "Metric")) %>%
    mutate(
      combo = factor(.data$combo, levels = c("Combined model", "Coverage model", "Prop. short model")),
      Metric = factor(.data$Metric, levels = c("Acc", "F1", "Sens", "Spec"))
    ) %>%
    arrange(.data$Metric, .data$combo) %>%
    mutate(across(where(is.factor), as.character))
}

default_component_files <- function(component_dir) {
  if (!dir.exists(component_dir)) return(character())
  paths <- list.files(component_dir, pattern = "\\.csv$", full.names = TRUE)
  excluded <- grepl(
    "__|^F2C_source_data_csv|"
      %+% "^F3C_source_data_csv|^Table for creating|"
      %+% "supporting_data_csv_test_confusion|^F3E_test_MFC",
    basename(paths)
  )
  sort(paths[!excluded])
}

`%+%` <- paste0

panel_overrides <- list(
  Fig1A = file.path(
    project_root, "Final Tables and Figures",
    "Supp_Table_1_all_events_for_swim_plot_INDEX_DATES_privacy_protected.csv"
  ),
  Fig1B = file.path(
    component_root, "Figure_1", "panel_B",
    c(
      "F1B_figure_source_counts_csv.csv",
      "F1B_figure_source_counts_patient_audit_csv.csv",
      "F1B_figure_source_csv.csv"
    )
  ),
  # The locked Figure 2B-D and Extended Data Figure 3A-C panels are the
  # all-evaluable July 2026 versions.  The generated component directories
  # also contain older training-only CSVs, so pin the exact source exports
  # written by the canonical 2_4B generator.
  Fig2B = file.path(
    canonical_root, "final_manuscript_objects",
    "additional_all_evaluable_longitudinal_panels",
    "Figure_2B_all_evaluable_BM_zscore_longitudinal_source_data.csv"
  ),
  Fig2C = file.path(
    canonical_root, "final_manuscript_objects",
    "additional_all_evaluable_longitudinal_panels",
    "Figure_2C_all_evaluable_blood_zscore_longitudinal_source_data.csv"
  ),
  Fig2D = file.path(
    canonical_root, "final_manuscript_objects",
    "additional_all_evaluable_longitudinal_panels",
    "Figure_2D_all_evaluable_fragmentomics_longitudinal_source_data.csv"
  ),
  ED3A = file.path(
    canonical_root, "final_manuscript_objects",
    "additional_all_evaluable_longitudinal_panels",
    "Extended_Data_Figure_3A_all_evaluable_BM_raw_longitudinal_source_data.csv"
  ),
  ED3B = file.path(
    canonical_root, "final_manuscript_objects",
    "additional_all_evaluable_longitudinal_panels",
    "Extended_Data_Figure_3B_all_evaluable_blood_raw_longitudinal_source_data.csv"
  ),
  ED3C = file.path(
    canonical_root, "final_manuscript_objects",
    "additional_all_evaluable_longitudinal_panels",
    "Extended_Data_Figure_3C_all_evaluable_fragmentomics_longitudinal_source_data.csv"
  ),
  Fig3D = file.path(
    canonical_root, "final_manuscript_objects", "01_main_figures", "Figure_3",
    "Figure_3D", "F3D_source_data_csv.csv"
  ),
  Fig4C = file.path(
    canonical_root, "final_manuscript_objects", "01_main_figures", "Figure_4",
    "Figure_4C", "F4C_source_data_csv.csv"
  ),
  Fig3E = file.path(
    project_root, "Output_tables_2025", "Source_data",
    c("Fig4K_cfWGS_vs_clinical_assays_EasyM_BM_source_data.csv",
      "Fig4K_cfWGS_vs_clinical_assays_EasyM_BM_correlations.csv")
  ),
  # (An earlier partial fix for the Figure 3A/4A/ED5A/ED7A mix-up lived here.
  # It assigned Fig4A the BM ROC and ED5A the blood performance summary, which
  # rotated the four sheets rather than correcting them, and it still pointed at
  # the pre-50-repeat run. It is superseded by the four overrides below.)
  Fig4D = file.path(
    component_root, "Figure_4", "panel_D",
    c("F4D_source_data_csv.csv", "Fig5K_cfWGS_vs_clinical_assays_EasyM_blood_correlations.csv")
  ),
  ED2B = file.path(
    component_root, "Extended_Data_Figure_2", "panel_B",
    "ED2B_all_evaluable_source_data_csv.csv"
  ),
  ED2C = file.path(
    component_root, "Extended_Data_Figure_2", "panel_C",
    "ED2C_source_data_csv.csv"
  ),
  # Panel G plots all evaluable baseline pairs (currently n = 41). The undated
  # `Extended_Data_Figure_2G_mutation_overlap_*` files are the earlier
  # 9-patient subset and disagree with the panel subtitle and legend
  # (median 4.7%, IQR 0.7-17.3%). The locked panel is ED2G_all_eval.png.
  ED2G = file.path(
    component_root, "Extended_Data_Figure_2", "panel_G",
    c(
      "ED2G_all_evaluable_source_data_csv.csv",
      "ED2G_all_evaluable_summary_csv.csv"
    )
  ),
  ED4all = file.path(
    component_root, "Extended_Data_Figure_4", "panel_all",
    c(
      "ED4_all_evaluable_feature_matrix_csv.csv",
      "ED4_all_evaluable_spearman_source_data_csv.csv"
    )
  ),
  ED6J = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "F4C_cfWGS_prob_vs_time_all_train_test_BM_outcome_available_source_data_2026-07-23.csv"
  ),
  ED6K = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "Fig_4D_time_to_relapse_footer3cols_BM_muts_all_train_test_outcome_available_source_data_2026-07-24.csv"
  ),
  ED8E = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "F4C_cfWGS_prob_vs_time_all_train_test_blood_outcome_available_source_data_2026-07-31.csv"
  ),
  ED8F = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "Fig_4D_time_to_relapse_footer3cols_blood_muts_all_train_test_outcome_available_source_data_2026-07-31.csv"
  ),
  ED10B = file.path(
    component_root, "Extended_Data_Figure_10", "panel_B",
    c(
      "ED10B_supporting_data_csv.csv",
      "EDFIG10A_B_B__supporting_data_csv__Emergent_CNA_loss_events.csv"
    )
  ),

  # ------------------------------------------------------------------------ #
  # 2026-07-29 audit fix: pin panels whose component directories still hold
  # only the 2026-02-25 exports.
  #
  # Without these overrides `default_component_files()` picks up the February
  # CSV that happens to sit in the component directory, even though the locked
  # panel in the figure deck was re-rendered from the July 2026 pipeline. That
  # is how Extended Data Fig. 6B came to ship a one-year-maintenance cfWGS
  # hazard ratio of 24.14 in the Source Data while the plotted panel and the
  # Abstract both report 7.37 (1.75-31.07).
  #
  # ED6J/ED6K/ED8E/ED8F above already used this pattern, which is why those
  # four sheets were correct. These entries extend it to the remaining panels.
  # ------------------------------------------------------------------------ #

  # Landmark progressor sensitivity barplots. The February exports were the
  # pre-expansion progression dataset; the current exports use the updated
  # relapse/follow-up data that the hazard ratios in the same figures already use.
  ED6A = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "Supp_6A_BM_sensitivity_barplot_source_data_2026-07-23.csv"
  ),
  ED8A = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "Supp_8A_blood_sensitivity_barplot_source_data_2026-07-23.csv"
  ),

  # Landmark hazard-ratio forest plots.
  ED6B = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "Supp_Figure_6B_BM_HR_plot_source_data_2026-07-23.csv"
  ),
  ED8B = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "SuppFig8B_blood_HR_plot_source_data_2026-07-23.csv"
  ),

  # Fixed-window progression sensitivity. The undated
  # `*_sensitivity_windows_source_data.csv` files in the component directories
  # predate the prospective labelling rules and report 17 samples / 6 patients
  # where the manuscript reports 43 samples / 18 patients at 180 days. The
  # prospective exports also carry N_nonevaluable_with_assay, which the
  # manuscript now cites when reporting evaluability.
  ED6I = file.path(
    project_root, "Output_tables_2025", "detection_progression_updated6",
    "prospective_timewindow_qc", "prospective_BM_cfWGS_timewindow_results.csv"
  ),
  ED8D = file.path(
    project_root, "Output_tables_2025", "detection_progression_updated6",
    "prospective_timewindow_qc", "prospective_blood_cfWGS_timewindow_results.csv"
  ),

  # Figure 3F: the headline one-year-maintenance KM panel.
  #
  # The component file Figure_3F_BM_cVAF_model_one_year_maintenance_KM_source_data.csv
  # contains only 16 patients (10 MRD-negative / 6 MRD-positive), whereas the
  # manuscript and the plotted risk table report 20 patients with 9 events
  # (13 negative / 7 positive) and HR 7.37 (1.75-31.07). The assay variable,
  # label and timepoint in that file are correct - it is a subset, not a
  # different analysis - and the current 20-patient export exists alongside it.
  #
  # Figure 4E was checked at the same time and is correct (21 patients, matching
  # the reported blood-derived Sites-model analysis), so no override is needed
  # there.
  Fig3F = file.path(
    project_root, "Output_tables_2025", "Source_data",
    paste0(
      "All_train_test_KM_source_data_1yr_maintenance_",
      "BM_zscore_only_detection_rate_call_2026-07-23.csv"
    )
  ),

  # ------------------------------------------------------------------------ #
  # Remaining KM panels: same stale-freeze problem the Fig3F override fixes.
  #
  # The note above states Figure 4E was checked and needed no override. That
  # check confirmed the patient count (21) but not the follow-up times. Without
  # an override, `default_component_files()` selects a component CSV built from
  # the 2026-07-01/07-02 freeze, which predates the 2026-07-03 relapse and
  # follow-up refresh. The patients are the same; the event counts and censoring
  # times are not, so the exported panels disagree with the locked components:
  #
  #   panel  stale (07-01/02)          current (07-23)           locked panel
  #   Fig4E  21, 9 ev, p=0.0252        21, 9 ev, p=0.0237        p = 0.02
  #   ED6C   30, 8 ev, p=0.0162        30, 10 ev, p=0.0140       p = 0.01
  #   ED6D   17, 5 ev, p=0.452         17, 6 ev, p=0.838         --
  #   ED6E   15, 8 ev, p=0.0346        22, 12 ev, p=0.1409       p = 0.14
  #   ED6F   35, 11 ev, p=0.296        38, 18 ev, p=0.611        --
  #   ED8C   23, 11 ev, p=0.157        24, 13 ev, p=0.0509       --
  #
  # The current exports reproduce the locked components; the stale ones do not.
  # ED6E is the clearest case: the stale freeze would publish a nominally
  # significant p = 0.03 where the locked panel reports p = 0.14.
  #
  # ED6G/ED6H (EasyM) need no override: their earliest available export is
  # 2026-07-07, which already postdates the refresh.
  # ------------------------------------------------------------------------ #
  Fig4E = file.path(
    project_root, "Output_tables_2025", "Source_data",
    paste0(
      "All_train_test_KM_source_data_1yr_maintenance_",
      "Blood_zscore_only_sites_call_2026-07-23.csv"
    )
  ),
  ED6C = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "All_train_test_KM_source_data_1yr_maintenance_Flow_Binary_2026-07-23.csv"
  ),
  ED6D = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "All_train_test_KM_source_data_1yr_maintenance_Adaptive_Binary_2026-07-23.csv"
  ),
  ED6E = file.path(
    project_root, "Output_tables_2025", "Source_data",
    paste0(
      "All_train_test_KM_source_data_post_transplant_",
      "BM_zscore_only_detection_rate_call_2026-07-23.csv"
    )
  ),
  ED6F = file.path(
    project_root, "Output_tables_2025", "Source_data",
    "All_train_test_KM_source_data_post_transplant_Flow_Binary_2026-07-23.csv"
  ),
  ED8C = file.path(
    project_root, "Output_tables_2025", "Source_data",
    paste0(
      "All_train_test_KM_source_data_post_transplant_",
      "Blood_zscore_only_sites_call_2026-07-23.csv"
    )
  ),

  # ------------------------------------------------------------------------ #
  # Figure 3A / 4A / ED5A / ED7A: correct BOTH the content assignment and the
  # analysis vintage.
  #
  # 6_19_Sync_ROC_Main_Operating_Extended_Data.R moved the 50-repeat mean ROC
  # curves into the main figures and pushed the sensitivity/specificity
  # summaries out to Extended Data. The component tree was never renamed to
  # match, so it still files the BM ROC under Extended_Data_Figure_5/panel_A.
  # default_component_files() therefore handed every one of these four sheets
  # the wrong panel's content, and the files it picked are the pre-50-repeat
  # run as well. Against the locked pages:
  #
  #   Figure 3A  = BM mean-ROC     AUC 0.81/0.79/0.77/0.74/0.72/0.64/0.64
  #   Figure 4A  = Blood mean-ROC  AUC 0.75/0.71/0.69/0.63/0.63/0.56/0.35
  #   ED Fig 5A  = BM fold-wise sensitivity/specificity   (Sites 62%/82%)
  #   ED Fig 7A  = Blood fold-wise sensitivity/specificity (M+f full 71%/70%)
  #
  # The two 50-repeat sources below reproduce all four legends exactly.
  # ------------------------------------------------------------------------ #
  Fig3A = file.path(
    project_root, "Output_figures_2025", "patient_grouped_repeated_nested_cv",
    "2026-08-05_compact_all_model_roc_50repeats_v3",
    "compact_all_model_mean_roc_source_data.csv"
  ),
  Fig4A = file.path(
    project_root, "Output_figures_2025", "patient_grouped_repeated_nested_cv",
    "2026-08-05_compact_all_model_roc_50repeats_v3",
    "compact_all_model_mean_roc_source_data.csv"
  ),
  ED5A = file.path(
    project_root, "Output_figures_2025", "patient_grouped_repeated_nested_cv",
    "2026-08-14_no_grid_larger_titles_50repeats_v2",
    "grouped_cv_fold_operating_point_source_data.csv"
  ),
  ED7A = file.path(
    project_root, "Output_figures_2025", "patient_grouped_repeated_nested_cv",
    "2026-08-14_no_grid_larger_titles_50repeats_v2",
    "grouped_cv_fold_operating_point_source_data.csv"
  ),

  # ED Figure 9A/9B are the fragmentomics-only equivalents and were stale in the
  # same way (published AUC 0.55/0.49/0.46/0.46/0.41/0.39; the component-tree
  # file gives 0.56/0.53/0.52/0.50/0.41/0.38). Both come from the same
  # 50-repeat run, restricted to the FullFrag cohort.
  ED9A = file.path(
    project_root, "Output_figures_2025", "patient_grouped_repeated_nested_cv",
    "2026-08-14_no_grid_larger_titles_50repeats_v2",
    "repeated_grouped_cv_mean_roc_source_data.csv"
  ),
  ED9B = file.path(
    project_root, "Output_figures_2025", "patient_grouped_repeated_nested_cv",
    "2026-08-14_no_grid_larger_titles_50repeats_v2",
    "grouped_cv_fold_operating_point_source_data.csv"
  ),

  # Dilution-series panels. Fig3C is the pooled library-level Spearman analysis
  # across all 48 scored libraries, with all four MRD-negative references set
  # to 0% for the analysis. ED5D and ED7D retain the four-patient dilution-
  # series points used by their LoD panels, with the same zero-reference rule.
  #
  # Note for Fig3C: the canonical pooled CSV sits in the component directory,
  # but `default_component_files()` excludes names matching
  # "^F3C_source_data_csv", so the single-patient file was selected instead.
  # An explicit override bypasses that exclusion.
  Fig3C = file.path(
    component_root, "Figure_3", "panel_C",
    "F3C_source_data_csv_all_four_patients.csv"
  ),
  ED5D = file.path(
    project_root, "Output_tables_2025", "Source_Data_Extended_Data",
    "SourceData_Fig4G_BM_patient_series_points_zero_final_tf.csv"
  ),
  ED7D = file.path(
    project_root, "Output_tables_2025", "Source_Data_Extended_Data",
    "SourceData_Fig5G_Blood_patient_series_points_zero_final_tf.csv"
  )
)

# These exact source names are the locked semantic authority for the KM
# panel-to-assay mapping. Assert them during export so a future override cannot
# silently place the wrong assay beneath a correct-looking panel title.
km_expected_source_basename <- c(
  ED6C = "All_train_test_KM_source_data_1yr_maintenance_Flow_Binary_2026-07-23.csv",
  ED6D = "All_train_test_KM_source_data_1yr_maintenance_Adaptive_Binary_2026-07-23.csv",
  ED6E = "All_train_test_KM_source_data_post_transplant_BM_zscore_only_detection_rate_call_2026-07-23.csv",
  ED6F = "All_train_test_KM_source_data_post_transplant_Flow_Binary_2026-07-23.csv",
  ED8C = "All_train_test_KM_source_data_post_transplant_Blood_zscore_only_sites_call_2026-07-23.csv"
)

# Main Figure 2A: exact patient rows and the assay-background values used by
# the canonical plotting script. Healthy-control assay values are scientific
# reference measurements; no control identifier or collection date is exposed.
build_figure_2a <- function() {
  aggregate_path <- file.path(
    project_root,
    "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds"
  )
  healthy_path <- file.path(
    project_root, "MRDetect_output_winter_2025", "Processed_R_outputs",
    "cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.rds"
  )
  stop_if_missing(c(aggregate_path, healthy_path), "Figure 2A")
  dat <- readRDS(aggregate_path)
  required <- c(
    "Patient", "Timepoint", "Sample_Code", "Date", "timepoint_info",
    "Cumulative_VAF_BM", "detect_rate_BM", "z_score_detection_rate_BM"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Figure 2A aggregate is missing: ", paste(missing, collapse = ", "))
  patient_rows <- dat %>%
    filter(Patient %in% c("CA-08", "CA-02")) %>%
    select(any_of(required)) %>%
    arrange(factor(Patient, levels = c("CA-08", "CA-02")), Date) %>%
    mutate(
      plot_group = "Patient trajectory",
      Exemplar_Status = recode(Patient, `CA-08` = "Relapsed", `CA-02` = "Stable"),
      source_table = basename(aggregate_path),
      .before = 1
    )

  healthy_rows <- readRDS(healthy_path) %>%
    filter(
      Mut_source == "BM_cells",
      Filter_source == "STR_encode",
      Study == "CHARM_healthy",
      Sample_type == "BM_cells",
      Patient %in% c("CA-08", "CA-02"),
      timepoint_info == "Diagnosis"
    ) %>%
    transmute(
      source_table = basename(healthy_path),
      Patient,
      Timepoint = NA_character_,
      Sample_Code = NA_character_,
      Date = as.Date(NA),
      timepoint_info = "Healthy_controls",
      Cumulative_VAF_BM = as.numeric(detection_rate_as_reads_detected_over_reads_checked),
      detect_rate_BM = NA_real_,
      z_score_detection_rate_BM = NA_real_,
      plot_group = "Healthy controls",
      Exemplar_Status = recode(Patient, `CA-08` = "Relapsed", `CA-02` = "Stable")
    )
  if (nrow(healthy_rows) != 52L) {
    stop("Figure 2A expected 52 healthy-reference assay rows; observed ", nrow(healthy_rows), call. = FALSE)
  }
  bind_rows(patient_rows, healthy_rows)
}

# Extended Data Figure 1: exact frozen oncoprint matrices in long form.
build_ed1 <- function() {
  matrix_paths <- c(
    BM = file.path(project_root, "heatmap_matrix_BM_Sep2025.rds"),
    Blood = file.path(project_root, "heatmap_matrix_blood_Sep2025_updated.rds")
  )
  annotation_path <- file.path(
    project_root, "Output_tables_2025_updated",
    "extended_data_figure_1_plot_annotations.csv"
  )
  stop_if_missing(c(matrix_paths, annotation_path), "Extended Data Figure 1")

  annotations <- read_csv(annotation_path, show_col_types = FALSE) %>%
    arrange(.data$column_order)
  required_annotation_columns <- c(
    "Patient_Timepoint", "Patient", "column_order", "Cohort",
    "TumourFraction", "Paired", "sample_type", "BM_Mutation_Count",
    "Blood_Mutation_Count", "IGH_FGFR3", "IGH_CCND1", "IGH_MAF",
    "del17p", "del1p", "amp1q", "del13q", "hyperdiploid"
  )
  missing_annotation_columns <- setdiff(
    required_annotation_columns, names(annotations)
  )
  if (length(missing_annotation_columns)) {
    stop(
      "Extended Data Figure 1 annotation layer is missing: ",
      paste(missing_annotation_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(annotations$Patient_Timepoint)) {
    stop("Extended Data Figure 1 annotations are not unique by column.", call. = FALSE)
  }

  cells <- imap_dfr(matrix_paths, function(path, specimen_source) {
    mat <- readRDS(path)
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop("Oncoprint matrix lacks row or column names: ", path, call. = FALSE)
    }
    as.data.frame(as.table(as.matrix(mat)), stringsAsFactors = FALSE) %>%
      rename(feature = Var1, Patient_Timepoint = Var2, value = Freq) %>%
      mutate(
        source_table = basename(path),
        specimen_source = specimen_source,
        .before = 1
      )
  })

  fish_features <- c(
    "IGH_FGFR3", "IGH_CCND1", "IGH_MAF", "del17p", "del1p", "amp1q",
    "del13q", "hyperdiploid"
  )
  fish_long <- annotations %>%
    select("Patient_Timepoint", all_of(fish_features)) %>%
    pivot_longer(
      all_of(fish_features), names_to = "feature", values_to = "FISH_Positive"
    ) %>%
    mutate(FISH_Positive = as.logical(.data$FISH_Positive))

  # Preserve the exact row order/splits created by 2_2.  CSV readers do not
  # retain factor levels, and a downstream alphabetical sort would change the
  # published oncoprint even when every cell value is otherwise correct.
  feature_order <- base::union(
    rownames(readRDS(matrix_paths[["BM"]])),
    rownames(readRDS(matrix_paths[["Blood"]]))
  )
  cna_features <- c("del1p", "amp1q", "del13q", "del17p", "hyperdiploid")
  translocation_features <- c("IGH_MAF", "IGH_CCND1", "IGH_MYC", "IGH_FGFR3")
  row_groups <- case_when(
    feature_order %in% cna_features ~ "CNAs",
    feature_order %in% translocation_features ~ "Translocations",
    TRUE ~ "Mutations"
  )
  feature_contract <- tibble(
    feature = feature_order,
    feature_order = seq_along(feature_order),
    row_group = row_groups
  )

  # The original script first scaffolds the union of BM and cfDNA rows/columns
  # with NA, then fills the specimen-specific matrices. Recreate that complete
  # grid here so absent specimens remain grey triangles instead of disappearing
  # from the public table.
  complete_cells <- tidyr::expand_grid(
    Patient_Timepoint = annotations$Patient_Timepoint,
    specimen_source = names(matrix_paths),
    feature = feature_order
  ) %>%
    left_join(
      cells %>% select(
        "Patient_Timepoint", "specimen_source", "feature", "value",
        "source_table"
      ),
      by = c("Patient_Timepoint", "specimen_source", "feature")
    )

  exported_cells <- complete_cells %>%
    inner_join(
      annotations %>% select(-all_of(fish_features)),
      by = "Patient_Timepoint"
    ) %>%
    left_join(fish_long, by = c("Patient_Timepoint", "feature")) %>%
    left_join(feature_contract, by = "feature") %>%
    mutate(
      sample = .data$Patient,
      specimen_source = factor(.data$specimen_source, levels = c("BM", "Blood"))
    ) %>%
    arrange(.data$column_order, .data$specimen_source, .data$feature_order)

  expected_cells <- nrow(annotations) * length(feature_order) * length(matrix_paths)
  if (
    nrow(exported_cells) != expected_cells ||
      anyNA(exported_cells$feature_order) ||
      anyDuplicated(exported_cells[c("Patient_Timepoint", "specimen_source", "feature")])
  ) {
    stop(
      "Extended Data Figure 1 export is not a complete, unique ",
      length(feature_order), "-feature x ", nrow(annotations),
      "-patient x 2-specimen grid.",
      call. = FALSE
    )
  }
  exported_cells
}

# Extended Data Figure 2D is the chr1 depth-ratio/copy-number line plot for
# VA-09 BM WGS. The former provenance note incorrectly reduced it to a
# PowerPoint-only schematic. The frozen Sequenza chr1 segments reproduce the
# red segment track and the FISH workbook supplies the marked 1q probe interval.
build_ed2d <- function() {
  segment_path <- file.path(
    project_root, "Oct 2024 data", "Sequenza", "All_Segments_400",
    "TFRIM4_0183_Bm_P_WG_VA-09-01-O-DNA_segments.txt"
  )
  fish_path <- file.path(project_root, "Clinical data", "FISH probe locations.xlsx")
  stop_if_missing(c(segment_path, fish_path), "Extended Data Figure 2D")

  segments <- read_tsv(segment_path, show_col_types = FALSE) %>%
    filter(chromosome == "chr1") %>%
    transmute(
      source_table = basename(segment_path),
      record_type = "Sequenza segment",
      sample = "TFRIM4_0183_Bm_P_WG_VA-09-01-O-DNA",
      chromosome,
      start = start.pos,
      end = end.pos,
      depth_ratio = depth.ratio,
      total_copy_number = CNt,
      allele_A_copy_number = A,
      allele_B_copy_number = B,
      probe_product = NA_character_,
      probe_target = NA_character_,
      probe_location_hg38 = NA_character_
    )

  probe <- read_excel(fish_path) %>%
    filter(Chromosome == "1q") %>%
    transmute(
      source_table = basename(fish_path),
      record_type = "FISH probe interval",
      sample = "TFRIM4_0183_Bm_P_WG_VA-09-01-O-DNA",
      chromosome = "chr1",
      start = as.numeric(gsub(",", "", sub(".*:", "", sub("-.*", "", `Location in hg38`)))),
      end = as.numeric(gsub(",", "", sub(".*-", "", `Location in hg38`))),
      depth_ratio = NA_real_,
      total_copy_number = NA_real_,
      allele_A_copy_number = NA_real_,
      allele_B_copy_number = NA_real_,
      probe_product = Product,
      probe_target = Target,
      probe_location_hg38 = `Location in hg38`
    )

  bind_rows(segments, probe)
}

metric_panel <- function(model_map) {
  youden_path <- file.path(
    project_root, "Output_tables_2025",
    "cfWGS_model_metrics_youden_threshold_test_cohort2.csv"
  )
  sens95_path <- file.path(
    project_root, "Output_tables_2025",
    "cfWGS_model_metrics_fixed_95sens_updated_test_cohort2.csv"
  )
  stop_if_missing(c(youden_path, sens95_path), "Frozen test-cohort metrics")
  youden <- read_csv(youden_path, show_col_types = FALSE) %>%
    filter(model %in% names(model_map)) %>%
    mutate(operating_point = "Youden")
  fixed <- read_csv(sens95_path, show_col_types = FALSE) %>%
    filter(model %in% names(model_map)) %>%
    mutate(
      ppv = precision,
      npv = NA_real_,
      operating_point = "95% sensitivity"
    ) %>%
    select(any_of(names(youden)), operating_point)
  bind_rows(youden, fixed) %>%
    mutate(
      source_table = if_else(
        operating_point == "Youden", basename(youden_path), basename(sens95_path)
      ),
      n_evaluable = case_when(
        grepl("^BM_", model) ~ 10L,
        grepl("^Blood_", model) ~ 8L,
        grepl("^Fragmentomics_", model) ~ 12L,
        TRUE ~ NA_integer_
      ),
      n_mrd_positive = case_when(
        grepl("^BM_", model) ~ 5L,
        grepl("^Blood_", model) ~ 3L,
        grepl("^Fragmentomics_", model) ~ 7L,
        TRUE ~ NA_integer_
      ),
      n_mrd_negative = n_evaluable - n_mrd_positive,
      combo = unname(model_map[model]),
      .before = 1
    ) %>%
    arrange(match(model, names(model_map)), match(operating_point, c("Youden", "95% sensitivity")))
}

# Full-cohort ROC source rows for locked Extended Data Figure 7E.
build_ed7e <- function() {
  # ED7E is the frozen full-training/refit ROC. The revision change log
  # explicitly distinguishes it from expanded-test performance (ED7C).
  # The original 04db00ce plotting block used `blood_obj$models` and a separate
  # `valid_df <- train_blood`. Raster comparison against the retained updated3
  # PNG establishes the preserved equivalents: the refit models come from
  # nested_blood_validation_updated3.rds, while the exact 34-row evaluable
  # frame is retained in the run2 model trainingData. Threshold markers were
  # calculated separately from data_scored_masked in the original script, so
  # they are exported as their own records rather than incorrectly recomputed
  # from the refit ROC probabilities.
  model_path <- file.path(project_root, "nested_blood_validation_updated3.rds")
  frame_path <- file.path(project_root, "nested_blood_validation_updated3_run2.rds")
  stop_if_missing(c(model_path, frame_path), "Extended Data Figure 7E")
  object <- readRDS(model_path)
  frame_object <- readRDS(frame_path)
  evaluation_data <- frame_object$models$Blood_plus_fragment_min$trainingData
  if (is.null(evaluation_data) || nrow(evaluation_data) != 34L ||
      !".outcome" %in% names(evaluation_data)) {
    stop("Frozen ED7E evaluation frame is not the expected 34-row caret trainingData.", call. = FALSE)
  }
  model_map <- c(
    Blood_zscore_only_sites = "Sites model",
    Blood_plus_fragment_min = "Combined model"
  )
  curves <- imap_dfr(model_map, function(label, model_name) {
    fit <- object$models[[model_name]]
    if (is.null(fit)) {
      stop("Frozen ED7E refit object lacks model: ", model_name, call. = FALSE)
    }
    probabilities <- predict(fit, newdata = evaluation_data, type = "prob")[, "pos"]
    tibble(
      plot_role = "ROC input",
      source_table = paste(basename(model_path), "evaluated_on", basename(frame_path)),
      model = model_name,
      combo = label,
      training_row = seq_along(probabilities),
      MRD_truth = as.integer(evaluation_data$.outcome == "pos"),
      probability = as.numeric(probabilities),
      threshold = NA_real_,
      fpr = NA_real_,
      tpr = NA_real_
    )
  })
  markers <- tribble(
    ~plot_role, ~source_table, ~model, ~combo, ~training_row, ~MRD_truth, ~probability, ~threshold, ~fpr, ~tpr,
    "Threshold marker", "canonical data_scored_masked operating point", "Blood_zscore_only_sites", "Sites model", NA_integer_, NA_integer_, NA_real_, 0.432, 0.500, 0.750,
    "Threshold marker", "canonical data_scored_masked operating point", "Blood_zscore_only_sites", "Sites model", NA_integer_, NA_integer_, NA_real_, 0.380, 1.000, 0.938,
    "Threshold marker", "canonical data_scored_masked operating point", "Blood_plus_fragment_min", "Combined model", NA_integer_, NA_integer_, NA_real_, 0.435, 0.167, 0.875
  )
  bind_rows(curves, markers)
}

selected_model_path <- file.path(
  project_root, "Output_tables_2025", "selected_combo_models_2026-02-16.rds"
)
selected_threshold_path <- file.path(
  project_root, "Output_tables_2025", "selected_combo_thresholds_2026-02-16.rds"
)
training_call_path <- file.path(
  project_root, "Output_tables_2025", "all_patients_with_BM_and_blood_calls_updated5.rds"
)
current_call_path <- file.path(
  project_root, "Output_tables_2025", "all_patients_with_BM_and_blood_calls_updated6.rds"
)
# Clinical-comparator confusion panels (ED5E-G, ED7F-H).
#
# This previously pointed at the copy of `all_patients_with_BM_and_blood_calls_
# updated6.rds` nested under `final_manuscript_objects/generated/
# native_script_runs/FIG1A_swim_plot/`. That nested file is a snapshot taken
# during a Figure 1A native-script run and has 586 rows against the live
# table's 581; it is stale for the BM call columns, so the exported panels
# disagreed with the locked components:
#
#   panel  snapshot MFC pairs   live MFC pairs   locked panel total
#   ED5E   14                   20               20  (9+2+2+7)
#   ED5F   14                   18               18  (10+2+1+5)
#   ED5G   22                   28               --
#   ED7F   21                   21               unchanged
#   ED7G   19                   19               unchanged
#   ED7H   14                   25               --
#
# The live table reproduces the locked ED5E and ED5F counts exactly, so the
# canonical 3_2 confusion-matrix blocks were run against it, not the snapshot.
# ED7F/ED7G are unaffected because the blood call columns agree in both files.
clinical_call_path <- current_call_path
stop_if_missing(
  c(
    selected_model_path, selected_threshold_path, training_call_path,
    current_call_path, clinical_call_path
  ),
  "Sample-level classifier source data"
)

selected_models <- readRDS(selected_model_path)
selected_thresholds <- readRDS(selected_threshold_path)
training_calls <- readRDS(training_call_path)
current_calls <- readRDS(current_call_path)
clinical_calls <- readRDS(clinical_call_path)

model_training_probability <- function(fit, model_name) {
  if (is.null(fit$trainingData) || !".outcome" %in% names(fit$trainingData)) {
    stop("Model lacks frozen trainingData: ", model_name, call. = FALSE)
  }
  as.numeric(predict(fit, newdata = fit$trainingData, type = "prob")[, "pos"])
}

match_frozen_training_rows <- function(rich_model_name) {
  fit <- selected_models[[rich_model_name]]
  train <- fit$trainingData
  predictors <- setdiff(names(train), ".outcome")
  required <- c("Patient", "Sample_Code", "Timepoint", "timepoint_info", "Cohort", "MRD_truth", predictors)
  missing <- setdiff(required, names(training_calls))
  if (length(missing)) {
    stop("Training call table lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  candidate <- training_calls %>%
    filter(!is.na(.data$MRD_truth))

  row_matches <- vapply(seq_len(nrow(train)), function(i) {
    keep <- as.integer(candidate$MRD_truth) == as.integer(train$.outcome[[i]] == "pos")
    for (column_name in predictors) {
      target <- train[[column_name]][[i]]
      values <- candidate[[column_name]]
      if (is.numeric(values) || is.integer(values)) {
        keep <- keep & !is.na(values) & !is.na(target) &
          abs(as.numeric(values) - as.numeric(target)) <= 1e-8
      } else {
        keep <- keep & !is.na(values) & !is.na(target) &
          as.character(values) == as.character(target)
      }
    }
    hits <- which(keep)
    if (length(hits) != 1L) {
      stop(
        "Frozen training row ", i, " for ", rich_model_name,
        " matched ", length(hits), " candidate rows; expected exactly one.",
        call. = FALSE
      )
    }
    hits
  }, integer(1))
  if (anyDuplicated(row_matches)) {
    stop("Frozen training rows did not map one-to-one for ", rich_model_name, ".", call. = FALSE)
  }
  candidate[row_matches, , drop = FALSE]
}

build_primary_model_call_matrix <- function(assay) {
  specs <- list(
    BM = list(
      rich = "BM_base_zscore",
      model_1 = "BM_zscore_only_detection_rate",
      model_2 = "BM_base_zscore"
    ),
    Blood = list(
      rich = "Blood_plus_fragment",
      model_1 = "Blood_zscore_only_sites",
      model_2 = "Blood_plus_fragment"
    )
  )
  spec <- specs[[assay]]
  if (is.null(spec)) stop("Unknown assay: ", assay, call. = FALSE)
  model_names <- c(spec$model_1, spec$model_2)
  matched <- match_frozen_training_rows(spec$rich)
  fits <- selected_models[model_names]
  if (length(unique(vapply(fits, function(x) nrow(x$trainingData), integer(1)))) != 1L) {
    stop("Frozen model training dimensions disagree for ", assay, ".", call. = FALSE)
  }
  truth <- as.integer(fits[[1]]$trainingData$.outcome == "pos")
  if (!identical(truth, as.integer(fits[[2]]$trainingData$.outcome == "pos")) ||
      !identical(truth, as.integer(matched$MRD_truth))) {
    stop("Frozen model training rows are not aligned for ", assay, ".", call. = FALSE)
  }
  training <- matched %>%
    transmute(
      Patient, Sample_Code, Timepoint, timepoint_info, Cohort = "Training",
      MRD_Truth = truth,
      Model_1_Probability = model_training_probability(fits[[1]], model_names[[1]]),
      Model_1_Call = as.integer(Model_1_Probability >= selected_thresholds[[model_names[[1]]]]),
      Model_2_Probability = model_training_probability(fits[[2]], model_names[[2]]),
      Model_2_Call = as.integer(Model_2_Probability >= selected_thresholds[[model_names[[2]]]])
    )

  probability_columns <- paste0(model_names, "_prob")
  test <- current_calls %>%
    filter(
      .data$Cohort == "Non-frontline",
      !(.data$timepoint_info %in% c("Diagnosis", "Baseline")),
      !is.na(.data$MRD_truth),
      if_all(all_of(probability_columns), ~ !is.na(.x))
    ) %>%
    transmute(
      Patient, Sample_Code, Timepoint, timepoint_info, Cohort = "Test",
      MRD_Truth = as.integer(.data$MRD_truth),
      Model_1_Probability = as.numeric(.data[[probability_columns[[1]]]]),
      Model_1_Call = as.integer(Model_1_Probability >= selected_thresholds[[model_names[[1]]]]),
      Model_2_Probability = as.numeric(.data[[probability_columns[[2]]]]),
      Model_2_Call = as.integer(Model_2_Probability >= selected_thresholds[[model_names[[2]]]])
    )
  bind_rows(training, test) %>%
    mutate(
      Model_1 = model_names[[1]],
      Model_1_Threshold = as.numeric(selected_thresholds[[model_names[[1]]]]),
      Model_2 = model_names[[2]],
      Model_2_Threshold = as.numeric(selected_thresholds[[model_names[[2]]]]),
      .after = Cohort
    )
}

build_fragmentomics_call_matrix <- function(cohort) {
  model_names <- c(
    "Fragmentomics_mean_coverage_only_Full",
    "Fragmentomics_min_Full"
  )
  if (cohort == "Training") {
    matched <- match_frozen_training_rows("Fragmentomics_full_Full")
    fits <- selected_models[model_names]
    truth <- as.integer(fits[[1]]$trainingData$.outcome == "pos")
    if (!identical(truth, as.integer(fits[[2]]$trainingData$.outcome == "pos")) ||
        !identical(truth, as.integer(matched$MRD_truth))) {
      stop("Frozen fragmentomics training rows are not aligned.", call. = FALSE)
    }
    rows <- matched %>%
      transmute(
        Patient, Sample_Code, Timepoint, timepoint_info, Cohort = "Training",
        MRD_Truth = truth,
        Coverage_Probability = model_training_probability(fits[[1]], model_names[[1]]),
        Coverage_Call = as.integer(Coverage_Probability >= selected_thresholds[[model_names[[1]]]]),
        Combined_Probability = model_training_probability(fits[[2]], model_names[[2]]),
        Combined_Call = as.integer(Combined_Probability >= selected_thresholds[[model_names[[2]]]])
      )
  } else if (cohort == "Test") {
    probability_columns <- paste0(model_names, "_prob")
    rows <- current_calls %>%
      filter(
        .data$Cohort == "Non-frontline",
        !(.data$timepoint_info %in% c("Diagnosis", "Baseline")),
        !is.na(.data$MRD_truth),
        if_all(all_of(probability_columns), ~ !is.na(.x))
      ) %>%
      transmute(
        Patient, Sample_Code, Timepoint, timepoint_info, Cohort = "Test",
        MRD_Truth = as.integer(.data$MRD_truth),
        Coverage_Probability = as.numeric(.data[[probability_columns[[1]]]]),
        Coverage_Call = as.integer(Coverage_Probability >= selected_thresholds[[model_names[[1]]]]),
        Combined_Probability = as.numeric(.data[[probability_columns[[2]]]]),
        Combined_Call = as.integer(Combined_Probability >= selected_thresholds[[model_names[[2]]]])
      )
  } else {
    stop("Unknown fragmentomics cohort: ", cohort, call. = FALSE)
  }
  rows %>%
    mutate(
      Coverage_Threshold = as.numeric(selected_thresholds[[model_names[[1]]]]),
      Combined_Threshold = as.numeric(selected_thresholds[[model_names[[2]]]])
    )
}

build_clinical_call_matrix <- function(assay, panel) {
  assay_spec <- list(
    BM = list(pred = "BM_zscore_only_detection_rate_call", front = "BM_base_zscore_call"),
    Blood = list(pred = "Blood_zscore_only_sites_call", front = "Blood_zscore_only_sites_call")
  )[[assay]]
  if (is.null(assay_spec)) stop("Unknown clinical assay: ", assay, call. = FALSE)
  dat <- clinical_calls %>%
    mutate(landmark = recode(
      .data$timepoint_info,
      Post_transplant = "Post_ASCT",
      `1yr maintenance` = "Maintenance",
      .default = NA_character_
    ))
  if (panel == "nonfront") {
    rows <- dat %>%
      filter(
        .data$Cohort == "Non-frontline",
        !is.na(.data[[assay_spec$pred]]),
        !(.data$timepoint_info %in% c("Diagnosis", "Baseline"))
      )
  } else {
    landmark <- c(post_asct = "Post_ASCT", maintenance = "Maintenance")[[panel]]
    if (is.null(landmark)) stop("Unknown clinical panel: ", panel, call. = FALSE)
    rows <- dat %>%
      filter(
        .data$Cohort == "Frontline",
        !is.na(.data$landmark),
        !is.na(.data[[assay_spec$front]]),
        .data$landmark == .env$landmark
      )
  }
  rows %>%
    transmute(
      Patient, Sample_Code, Timepoint, timepoint_info, Cohort,
      cfWGS_Call = as.integer(.data[[assay_spec$pred]]),
      MFC_Call = as.integer(.data$Flow_Binary),
      clonoSEQ_Call = if (panel == "nonfront") NA_integer_ else as.integer(.data$Adaptive_Binary)
    ) %>%
    drop_empty_columns()
}

expected_sample_level_rows <- tribble(
  ~sheet_name, ~expected_n,
  "Fig3B", 70L,
  "Fig4B", 69L,
  # Row counts below follow the live clinical-call table (see clinical_call_path
  # above). The previous 15/16/44/23/21/47 values came from the stale nested
  # FIG1A_swim_plot snapshot.
  "ED5E", 21L,
  "ED5F", 20L,
  "ED5G", 42L,
  "ED7F", 23L,
  "ED7G", 21L,
  "ED7H", 49L,
  "ED9E", 65L,
  "ED9F", 37L
)

# Publication-panel schema contract. These are the fields consumed by the
# canonical plotting ports, not merely whatever columns happened to be present
# in the latest export. Every one of the 74 panels must be covered so an export
# layer regression fails here, before a visually plausible but scientifically
# wrong capsule can be built.
panel_schema_contract <- setNames(vector("list", nrow(locked_figures)), locked_figures$sheet_name)
set_contract <- function(panels, columns) {
  panel_schema_contract[panels] <<- rep(list(columns), length(panels))
}
set_contract("Fig1A", c("Patient", "event", "start_day_from_baseline", "end_day_from_baseline", "Cohort", "MFC", "clonoSEQ", "cfWGS_BM", "cfWGS_blood", "previously_treated"))
set_contract("Fig1B", c("box_id", "box_label", "n", "definition"))
set_contract("Fig2A", c("plot_group", "Exemplar_Status", "Patient", "Timepoint", "Cumulative_VAF_BM", "detect_rate_BM", "z_score_detection_rate_BM"))
set_contract(c("Fig2B", "Fig2C", "Fig2D", "ED3A", "ED3B", "ED3C"), c("Patient", "Weeks_Since_Baseline", "Num_days_to_closest_relapse", "Metric", "Value", "relapse_within_180"))
set_contract("Fig2E", c("Patient", "Metric", "Value", "terminal_group", "relapse_within_180", "wilcoxon_q_bh_global_12_tests"))
set_contract(c("Fig3A", "Fig4A"), c("cohort", "model", "display_name", "fpr", "mean_tpr", "repeat_pooled_auc_mean"))
set_contract(c("Fig3B", "Fig4B"), c("Patient", "Cohort", "MRD_Truth", "Model_1_Probability", "Model_1_Call", "Model_2_Probability", "Model_2_Call"))
set_contract("Fig3C", c("feature", "plot_label", "n_complete", "rho (Spearman)", "p-value"))
set_contract(c("Fig3D", "Fig4C"), c("landmark_tp", "Technology", "n_total", "n_pos", "pos_rate", "Cohort"))
set_contract(c("Fig3E", "Fig4D"), c("Patient", "Comparator", "x_plot", "y_plot", "relapse_cat", "rho", "p", "label"))
set_contract("Fig3F", c("Patient", "Time_to_event_months", "Relapsed_Binary", "BM_zscore_only_detection_rate_call", "Group"))
set_contract("Fig4E", c("Patient", "Time_to_event_months", "Relapsed_Binary", "Blood_zscore_only_sites_call", "Group"))
set_contract("ED1all", c("Patient_Timepoint", "specimen_source", "feature", "value", "Cohort", "TumourFraction", "FISH_Positive", "feature_order"))
set_contract("ED2A", c("Patient", "cohort", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"))
set_contract("ED2B", c("event", "tf_group", "Metric", "Value", "Percent"))
set_contract("ED2C", c("event", "sample", "Measure", "value"))
set_contract("ED2D", c("chromosome", "start", "end", "depth_ratio", "total_copy_number", "probe_location_hg38"))
set_contract("ED2E", c("Patient", "cohort", "Assay", "MutCount"))
set_contract("ED2F", c("Patient", "cohort", "Blood_Mutation_Count", "var", "x", "rho", "p"))
set_contract("ED2G", c("Patient", "BM_Mutation_Count", "cfDNA_Mutation_Count", "Percent_Overlap", "Cohort"))
set_contract("ED3D", c("platform", "control_id", "metric", "raw_rate", "metric_label", "paired_wilcoxon_p"))
set_contract("ED3E", c("control_id", "metric", "condition", "z_score", "metric_label"))
set_contract("ED3F", c("control_id", "feature", "platform", "raw_value", "paired_wilcoxon_p"))
set_contract("ED3G", c("control_id", "feature", "historical_value", "xplus_loo_corrected_value", "correction_state"))
set_contract("ED4all", c("Patient", "Metric1", "Metric2", "rho", "n_complete_pairs"))
set_contract(c("ED5A", "ED7A", "ED9B"), c("cohort", "model", "sensitivity_mean", "specificity_mean", "accuracy_mean"))
set_contract(c("ED5B", "ED7B", "ED9C"), c("combo", "Metric", "mean", "fixed"))
set_contract(c("ED5C", "ED7C", "ED9D"), c("model", "combo", "operating_point", "n_evaluable", "sensitivity", "specificity"))
set_contract(c("ED5D", "ED7D"), c("Patient", "Sample", "series_label", "line_group", "feature", "value", "panel_type", "thr", "call"))
set_contract(c("ED5E", "ED5F", "ED7F", "ED7G"), c("Patient", "cfWGS_Call", "MFC_Call", "clonoSEQ_Call"))
set_contract(c("ED5G", "ED7H"), c("Patient", "cfWGS_Call", "MFC_Call"))
set_contract(c("ED6A", "ED8A"), c("Assay", "N_tested", "N_positive", "Sensitivity", "Timepoint"))
set_contract(c("ED6B", "ED8B"), c("Landmark", "Assay", "HR", "CI_low", "CI_high"))
set_contract(c("ED6C", "ED6D", "ED6E", "ED6F", "ED8C"), c("Patient", "Time_to_event_months", "Relapsed_Binary", "Assay_Call", "assay_label", "Group", "timepoint_label"))
set_contract(c("ED6G", "ED6H"), c("Patient", "Time_to_event_months", "Relapsed_Binary", "EasyM_reference_threshold_binary", "EasyM_value", "Group"))
set_contract(c("ED6I", "ED8D"), c("Window_days", "Assay", "N_samples", "N_patients", "TP", "FN", "Sensitivity", "Specificity"))
set_contract(c("ED6J", "ED6K"), c("Patient", "days_before_event", "Relapsed_Binary", "mrd_status", "BM_zscore_only_detection_rate_call", "BM_zscore_only_detection_rate_prob"))
set_contract("ED7E", c("plot_role", "model", "MRD_truth", "probability", "threshold", "fpr", "tpr"))
set_contract(c("ED8E", "ED8F"), c("Patient", "days_before_event", "Relapsed_Binary", "mrd_status", "Blood_zscore_only_sites_call", "Blood_zscore_only_sites_prob"))
set_contract("ED9A", c("cohort", "model", "fpr", "mean_tpr", "repeat_pooled_auc_mean", "auc_cluster_boot_q025", "auc_cluster_boot_q975"))
set_contract(c("ED9E", "ED9F"), c("Patient", "MRD_Truth", "Coverage_Probability", "Coverage_Call", "Combined_Probability", "Combined_Call"))
set_contract(c("ED10A", "ED10B"), c("Patient", "Sample", "Event"))
if (any(lengths(panel_schema_contract) == 0L)) {
  stop("Missing schema contract for panel(s): ", paste(names(panel_schema_contract)[lengths(panel_schema_contract) == 0L], collapse = ", "), call. = FALSE)
}

manifest_rows <- list()
audit_rows <- list()

for (i in seq_len(nrow(locked_figures))) {
  item <- locked_figures[i, ]
  key <- item$sheet_name
  component_dir <- file.path(component_root, item$component_figure, item$component_panel)
  status <- "PASS"
  note <- "Canonical panel source data exported."
  source_paths <- character()

  data <- if (key == "Fig2A") {
    source_paths <- file.path(
      project_root, "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds"
    )
    build_figure_2a()
  } else if (key == "Fig2E") {
    source_paths <- file.path(
      canonical_root, "final_manuscript_objects",
      "additional_all_evaluable_longitudinal_panels",
      "reviewer_alternative_terminal_summary_v2",
      c(
        "terminal_sample_patient_level_source_data.csv",
        "terminal_sample_group_statistics.csv"
      )
    )
    status <- "PASS"
    note <- paste(
      "One terminal evaluable observation per patient and metric is exported",
      "together with the six plotted group summaries, two-sided Wilcoxon",
      "tests, and BH-adjusted q values."
    )
    build_figure_2e()
  } else if (key == "ED9C") {
    source_paths <- file.path(
      project_root, "Output_tables_2025",
      c(
        "cfWGS_model_metrics_youden_threshold2.csv",
        "cfWGS_model_metrics_fixed_95sens3_Feb2026.csv"
      )
    )
    note <- paste(
      "Locked fragmentomics training performance: Youden bars use the",
      "threshold2 snapshot plotted in the retained panel; 95%-sensitivity",
      "triangles use the fixed-threshold February 2026 table."
    )
    build_ed9c()
  } else if (key == "ED1all") {
    source_paths <- c(
      file.path(project_root, "heatmap_matrix_BM_Sep2025.rds"),
      file.path(project_root, "heatmap_matrix_blood_Sep2025_updated.rds")
    )
    build_ed1()
  } else if (key == "ED2D") {
    source_paths <- c(
      file.path(
        project_root, "Oct 2024 data", "Sequenza", "All_Segments_400",
        "TFRIM4_0183_Bm_P_WG_VA-09-01-O-DNA_segments.txt"
      ),
      file.path(project_root, "Clinical data", "FISH probe locations.xlsx")
    )
    status <- "PARTIAL_SOURCE"
    note <- paste(
      "Current locked chr1 depth-ratio/copy-number line plot.",
      "The Sequenza segment track and 1q FISH-probe interval are exported;",
      "the historical per-bin black depth-ratio scatter points are not",
      "separately retained in the current final manuscript objects."
    )
    ed2d_data <- build_ed2d()
    component_csv <- file.path(
      component_root, "Extended_Data_Figure_2", "panel_D",
      "Extended_Data_Figure_2D_chr1_depth_ratio_source_data.csv"
    )
    write_csv(ed2d_data, component_csv, na = "")
    ed2d_data
  } else if (key %in% c("Fig3B", "Fig4B")) {
    source_paths <- c(
      selected_model_path, selected_threshold_path,
      training_call_path, current_call_path
    )
    status <- "PASS"
    note <- "One row per sample with MRD truth, model probabilities, thresholds, and binary calls."
    build_primary_model_call_matrix(if_else(key == "Fig3B", "BM", "Blood"))
  } else if (key == "ED5C") {
    source_paths <- file.path(
      component_root, "Extended_Data_Figure_5", "panel_C",
      "ED5C_source_data_csv_current_test_cohort.csv"
    )
    combine_source_files(source_paths)
  } else if (key == "ED7C") {
    source_paths <- file.path(
      component_root, "Extended_Data_Figure_7", "panel_C",
      "ED7C_source_data_csv_current_test_cohort.csv"
    )
    combine_source_files(source_paths)
  } else if (key == "ED7E") {
    source_paths <- file.path(
      project_root,
      c("nested_blood_validation_updated3.rds", "nested_blood_validation_updated3_run2.rds")
    )
    note <- paste(
      "Frozen full-training/refit ROC: updated3 refit models evaluated on the",
      "preserved 34-row run2 frame, with separately retained canonical threshold",
      "coordinates. Expanded-test performance is ED7C."
    )
    build_ed7e()
  } else if (key == "ED9D") {
    source_paths <- file.path(
      component_root, "Extended_Data_Figure_9", "panel_D",
      "ED9D_source_data_csv_current_test_cohort.csv"
    )
    combine_source_files(source_paths)
  } else if (key == "ED9E") {
    source_paths <- c(selected_model_path, selected_threshold_path, training_call_path)
    build_fragmentomics_call_matrix("Training")
  } else if (key == "ED9F") {
    source_paths <- c(selected_model_path, selected_threshold_path, current_call_path)
    build_fragmentomics_call_matrix("Test")
  } else if (key %in% c("ED5E", "ED5F", "ED5G", "ED7F", "ED7G", "ED7H")) {
    clinical_spec <- list(
      ED5E = c("BM", "post_asct"),
      ED5F = c("BM", "maintenance"),
      ED5G = c("BM", "nonfront"),
      ED7F = c("Blood", "post_asct"),
      ED7G = c("Blood", "maintenance"),
      ED7H = c("Blood", "nonfront")
    )[[key]]
    source_paths <- clinical_call_path
    note <- "One row per sample with binary cfWGS, MFC, and clonoSEQ calls where available."
    build_clinical_call_matrix(clinical_spec[[1]], clinical_spec[[2]])
  } else {
    source_paths <- panel_overrides[[key]]
    if (is.null(source_paths)) source_paths <- default_component_files(component_dir)
    if (length(source_paths) == 0L) {
      stop("No source data registered for locked panel ", key, call. = FALSE)
    }
    combine_source_files(source_paths)
  }

  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("Panel ", key, " produced no source-data rows.", call. = FALSE)
  }
  if (key %in% names(km_expected_source_basename) &&
      !identical(basename(source_paths), unname(km_expected_source_basename[[key]]))) {
    stop(
      "Locked KM source mapping changed for ", key, "; expected ",
      km_expected_source_basename[[key]], " but observed ",
      paste(basename(source_paths), collapse = ", "), call. = FALSE
    )
  }

  # Older Figure 2D exports stored mean coverage with its sign inverted for a
  # joint analysis.  The current all-evaluable 2_4B export already stores the
  # positive plotted scale.  Normalize either representation defensively.
  if (key == "Fig2D" && all(c("Metric", "Value") %in% names(data))) {
    coverage_median <- median(as.numeric(data$Value[data$Metric == "Mean.Coverage"]), na.rm = TRUE)
    coverage_is_inverted <- is.finite(coverage_median) && coverage_median < 0
    data <- data %>%
      mutate(
        analysis_value = if_else(
          Metric == "Mean.Coverage" & !coverage_is_inverted,
          -as.numeric(Value),
          as.numeric(Value)
        ),
        Value = if_else(
          Metric == "Mean.Coverage" & coverage_is_inverted,
          -as.numeric(Value),
          as.numeric(Value)
        ),
        value_scale = if_else(
          Metric == "Mean.Coverage",
          "raw plotted mean coverage; analysis_value stores the sign-inverted modeling value",
          "plotted value"
        )
      )
    status <- "CORRECTED_TO_LOCKED"
    note <- paste(
      "Mean coverage is exported on the positive raw/plotted scale used in the locked Figure 2D fragmentomics line plot.",
      "The sign-inverted value used by the joint analysis is retained as analysis_value."
    )
  }
  data <- relabel_publication_cohorts(
    data,
    paste0("Figure source-data panel ", key)
  )
  data <- minimalize_publication_panel(data, key)
  missing_schema_columns <- setdiff(panel_schema_contract[[key]], names(data))
  if (length(missing_schema_columns)) {
    stop(
      "Panel ", key, " violates its plotting schema contract; missing: ",
      paste(missing_schema_columns, collapse = ", "), call. = FALSE
    )
  }
  output_path <- file.path(panel_output_root, paste0(key, ".csv"))
  write_csv(data, output_path, na = "")
  if (key == "Fig2E") {
    figure2e_source_destination <- file.path(
      canonical_root, "final_manuscript_objects", "01_main_figures",
      "Figure_2", "Figure_2E", "Fig_2E_source_data.csv"
    )
    dir.create(
      dirname(figure2e_source_destination),
      recursive = TRUE,
      showWarnings = FALSE
    )
    copied <- file.copy(
      output_path,
      figure2e_source_destination,
      overwrite = TRUE
    )
    if (!isTRUE(copied)) {
      stop(
        "Failed to stage Figure 2E source data: ",
        figure2e_source_destination,
        call. = FALSE
      )
    }
  }

  manifest_rows[[key]] <- tibble(
    workbook = item$workbook,
    figure = item$figure,
    panel = item$panel,
    sheet_name = key,
    panel_csv = normalizePath(output_path, mustWork = TRUE),
    locked_pdf = normalizePath(item$locked_pdf, mustWork = TRUE),
    source_files = paste(normalizePath(source_paths, mustWork = TRUE), collapse = " | "),
    row_count = nrow(data),
    column_count = ncol(data),
    audit_status = status,
    audit_note = note
  )
  audit_rows[[key]] <- tibble(
    workbook = item$workbook,
    figure_panel = key,
    audit_status = status,
    row_count = nrow(data),
    column_count = ncol(data),
    locked_pdf_md5 = unname(tools::md5sum(item$locked_pdf)),
    note = note
  )
}

manifest <- bind_rows(manifest_rows)
audit <- bind_rows(audit_rows)
schema_contract_table <- imap_dfr(panel_schema_contract, ~ tibble(sheet_name = .y, required_column = .x))
write_csv(schema_contract_table, file.path(output_root, "panel_schema_contract.csv"), na = "")

for (j in seq_len(nrow(expected_sample_level_rows))) {
  check <- expected_sample_level_rows[j, ]
  path <- file.path(panel_output_root, paste0(check$sheet_name, ".csv"))
  dat <- read_csv(path, show_col_types = FALSE)
  if (nrow(dat) != check$expected_n) {
    stop(
      "Sample-level call audit failed for ", check$sheet_name,
      ": expected n=", check$expected_n, "; observed ", nrow(dat),
      call. = FALSE
    )
  }
}

expected_sheet_count <- nrow(locked_figures)
if (nrow(manifest) != expected_sheet_count) {
  stop("Expected ", expected_sheet_count, " panel sheets; exported ", nrow(manifest))
}
if (anyDuplicated(manifest$sheet_name)) stop("Duplicate workbook sheet names detected.")
if (any(manifest$row_count < 1L)) stop("One or more panel exports are empty.")

write_csv(manifest, file.path(output_root, "source_data_workbook_manifest.csv"), na = "")
write_csv(audit, file.path(output_root, "source_data_audit.csv"), na = "")

message(
  "Exported ", nrow(manifest), " locked panel datasets: ",
  sum(manifest$workbook == "main"), " main and ",
  sum(manifest$workbook == "extended"), " extended."
)
message("Audit statuses: ", paste(names(table(audit$audit_status)), table(audit$audit_status), collapse = "; "))
