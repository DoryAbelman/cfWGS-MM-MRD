#!/usr/bin/env Rscript

################################################################################
##  5_0_Build_Manuscript_Text_Number_Exports.R
##
##  Script purpose:
##    Build manuscript-writing helper exports: section-organized numerical
##    values, short prose snippets, and source pointers for updating manuscript
##    text after rerunning the pipeline or adding new samples.
##
##  Why this exists:
##    Figure/table exports tell us what panels and supplementary tables were
##    generated, but manuscript revision also needs paragraph-level numbers:
##    cohort sizes, sample counts, assay denominators, model-performance values,
##    follow-up dates, and sensitivity/concordance summaries. This script creates
##    a single auditable "where do I find the number for this sentence?" export.
##
##  Scientific and reproducibility notes:
##    - This script does not refit models, select thresholds, or change any
##      analysis result.
##    - Values are derived from already-generated pipeline outputs and retained
##      manuscript source tables.
##    - Each row records the source file and, where known, the source script or
##      final manuscript artifact that owns the value.
##    - Rows are intentionally broad at first: they cover the manuscript sections
##      most likely to need manual text updates after adding samples. Additional
##      rows can be added safely as more manuscript paragraphs are audited.
##
##  Inputs:
##    - cfWGS Manuscript Drafts/*.docx
##    - docs/manuscript_artifact_source_map.tsv
##    - Scripts_2025/Final_Scripts/script_index.tsv
##    - combined_clinical_data_updated_April2025.csv
##    - cohort_assignment_table_updated.rds
##    - Output_tables_2025/all_patients_with_BM_and_blood_calls_updated5.rds
##    - Exported_data_tables_clinical/*PFS*/relapse/follow-up tables
##    - Final Tables and Figures / Output_tables_2025 model and time-window CSVs
##
##  Outputs:
##    Scripts_2025/Final_Scripts/manuscript_writing/
##      manuscript_numbers_by_section.tsv
##      manuscript_numbers_by_section.xlsx
##      manuscript_numbers_by_section.md
##      manuscript_draft_paragraph_index.tsv
##      manuscript_numeric_paragraph_review.tsv
##      manuscript_paragraph_metric_audit.tsv
##      by_section/
##      README.md
##
##  How to run:
##    Rscript Scripts_2025/Final_Scripts/5_0_Build_Manuscript_Text_Number_Exports.R
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(writexl)
  library(xml2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

get_this_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  fallback <- file.path(
    "Scripts_2025",
    "Final_Scripts",
    "5_0_Build_Manuscript_Text_Number_Exports.R"
  )
  normalizePath(fallback, mustWork = FALSE)
}

script_path <- get_this_script_path()
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
out_dir <- file.path(script_dir, "manuscript_writing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

project_path <- function(...) file.path(project_root, ...)

file_status <- function(path) {
  if (file.exists(path)) "available" else "missing"
}

read_csv_if_exists <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, name_repair = "unique_quiet", ...)
}

read_delim_if_exists <- function(path, delim = "\t", ...) {
  if (!file.exists(path)) return(NULL)
  readr::read_delim(path, delim = delim, show_col_types = FALSE, ...)
}

read_rds_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

safe_n_distinct <- function(x) {
  if (is.null(x)) return(NA_integer_)
  dplyr::n_distinct(x[!is.na(x)])
}

fmt_int <- function(x) {
  if (is.na(x)) return(NA_character_)
  format(as.integer(x), big.mark = ",", scientific = FALSE)
}

fmt_pct <- function(x, digits = 0) {
  if (is.na(x)) return(NA_character_)
  paste0(format(round(100 * x, digits), nsmall = digits, scientific = FALSE), "%")
}

fmt_num <- function(x, digits = 2) {
  if (is.na(x)) return(NA_character_)
  format(round(as.numeric(x), digits), nsmall = digits, scientific = FALSE)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[1] else paths[1]
}

slugify <- function(x, max_len = 90) {
  x <- str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- str_replace_all(x, "^_+|_+$", "")
  x <- ifelse(nzchar(x), x, "unspecified")
  str_sub(x, 1, max_len)
}

new_row <- function(
  manuscript_section,
  manuscript_subsection,
  metric_id,
  statistic_label,
  value,
  formatted_value = as.character(value),
  numerator = NA,
  denominator = NA,
  units = "",
  suggested_text = "",
  source_file = "",
  source_script = "",
  related_artifacts = "",
  update_trigger = "rerun pipeline or add samples",
  caveat = "",
  review_status = "ready_for_review"
) {
  tibble::tibble(
    manuscript_section = manuscript_section,
    manuscript_subsection = manuscript_subsection,
    metric_id = metric_id,
    statistic_label = statistic_label,
    value = as.character(value),
    formatted_value = as.character(formatted_value),
    numerator = as.character(numerator),
    denominator = as.character(denominator),
    units = units,
    suggested_text = suggested_text,
    source_file = source_file,
    source_script = source_script,
    related_artifacts = related_artifacts,
    update_trigger = update_trigger,
    caveat = caveat,
    review_status = review_status
  )
}

script_index <- read_delim_if_exists(file.path(script_dir, "script_index.tsv"))
artifact_map <- read_delim_if_exists(project_path("docs", "manuscript_artifact_source_map.tsv"))

script_for_artifact <- function(artifact_ids) {
  if (is.null(artifact_map) || !"artifact_id" %in% names(artifact_map)) return("")
  out <- artifact_map %>%
    filter(.data$artifact_id %in% artifact_ids) %>%
    distinct(.data$source_script) %>%
    pull(.data$source_script)
  paste(out[!is.na(out) & nzchar(out)], collapse = "; ")
}

artifact_label <- function(ids) paste(ids, collapse = "; ")

rows <- list()
add_rows <- function(x) {
  rows[[length(rows) + 1]] <<- x
}

################################################################################
## Draft paragraph index
################################################################################

extract_docx_paragraphs <- function(docx_path) {
  tmp_dir <- tempfile("docx_extract_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(docx_path, files = "word/document.xml", exdir = tmp_dir)
  xml_path <- file.path(tmp_dir, "word", "document.xml")
  if (!file.exists(xml_path)) {
    return(tibble::tibble(
      draft_file = basename(docx_path),
      paragraph_index = integer(),
      style = character(),
      text = character()
    ))
  }

  doc <- xml2::read_xml(xml_path)
  ns <- xml2::xml_ns(doc)
  paras <- xml2::xml_find_all(doc, ".//w:p", ns)
  purrr::map_dfr(seq_along(paras), function(i) {
    p <- paras[[i]]
    text <- xml2::xml_find_all(p, ".//w:t", ns) %>%
      xml2::xml_text() %>%
      paste(collapse = "")
    style <- xml2::xml_find_first(p, ".//w:pStyle", ns) %>%
      xml2::xml_attr("val") %||% ""
    tibble::tibble(
      draft_file = basename(docx_path),
      paragraph_index = i,
      style = style,
      text = str_squish(text)
    )
  }) %>%
    filter(nzchar(.data$text))
}

draft_paths <- list.files(project_path("cfWGS Manuscript Drafts"), pattern = "\\.docx$", full.names = TRUE)
draft_index <- purrr::map_dfr(draft_paths, extract_docx_paragraphs)

annotate_draft_sections <- function(x) {
  if (!nrow(x)) return(x)
  x %>%
    group_by(.data$draft_file) %>%
    arrange(.data$paragraph_index, .by_group = TRUE) %>%
    mutate(
      is_major_heading = replace_na(.data$style %in% c("Title", "Heading1"), FALSE),
      is_minor_heading = replace_na(.data$style %in% "Heading2", FALSE),
      major_heading_raw = if_else(.data$is_major_heading, .data$text, NA_character_),
      major_block = cumsum(!is.na(.data$major_heading_raw))
    ) %>%
    tidyr::fill("major_heading_raw", .direction = "down") %>%
    group_by(.data$draft_file, .data$major_block) %>%
    mutate(
      minor_heading_raw = if_else(.data$is_minor_heading, .data$text, NA_character_)
    ) %>%
    tidyr::fill("minor_heading_raw", .direction = "down") %>%
    ungroup() %>%
    mutate(
      manuscript_section = coalesce(.data$major_heading_raw, "Front matter"),
      manuscript_subsection = coalesce(.data$minor_heading_raw, .data$manuscript_section),
      paragraph_role = case_when(
        .data$is_major_heading ~ "major_heading",
        .data$is_minor_heading ~ "minor_heading",
        replace_na(.data$style == "Bibliography", FALSE) ~ "reference_or_bibliography",
        TRUE ~ "body_or_caption"
      )
    ) %>%
    select(
      "draft_file",
      "paragraph_index",
      "style",
      "manuscript_section",
      "manuscript_subsection",
      "paragraph_role",
      "text"
    )
}

extract_numeric_tokens <- function(text) {
  hits <- stringr::str_extract_all(
    text,
    "\\b\\d+(?:[.,]\\d+)*(?:\\.\\d+)?(?:/\\d+)?%?\\b|\\b\\d+\\s*-\\s*\\d+\\b"
  )[[1]]
  hits <- hits[!is.na(hits) & nzchar(hits)]
  paste(unique(hits), collapse = "; ")
}

infer_number_source_hint <- function(text) {
  txt <- str_to_lower(text)
  hints <- character()
  if (str_detect(txt, "patient|sample|cohort|training|test|baseline|diagnosis")) {
    hints <- c(hints, "cohort/sample-count outputs")
  }
  if (str_detect(txt, "auc|sensitivity|specificity|ppv|npv|youden|threshold|classifier|model")) {
    hints <- c(hints, "model-performance tables")
  }
  if (str_detect(txt, "relapse|progression|follow-up|pfs|survival|maintenance|post-asct")) {
    hints <- c(hints, "survival/time-window outputs")
  }
  if (str_detect(txt, "concordance|discordance|mfc|clonoseq|easym|mrd")) {
    hints <- c(hints, "clinical concordance outputs")
  }
  if (str_detect(txt, "mutation|cna|copy number|translocation|fish|genomic|feature")) {
    hints <- c(hints, "genomic feature/concordance outputs")
  }
  if (str_detect(txt, "dilution|lod|limit of detection|tumor fraction|10")) {
    hints <- c(hints, "dilution-series outputs")
  }
  if (!length(hints)) hints <- "manual review"
  paste(unique(hints), collapse = "; ")
}

draft_index <- annotate_draft_sections(draft_index)
numeric_paragraph_review <- draft_index %>%
  mutate(
    numeric_tokens = vapply(.data$text, extract_numeric_tokens, character(1)),
    candidate_source_hint = vapply(.data$text, infer_number_source_hint, character(1)),
    review_priority = case_when(
      .data$paragraph_role %in% c("major_heading", "minor_heading", "reference_or_bibliography") ~ "low",
      str_detect(.data$manuscript_section, regex("front matter|references|acknowledgements|author contributions|competing interests", ignore_case = TRUE)) ~ "low",
      str_detect(.data$manuscript_section, regex("figures|tables|extended data|supplementary table legends", ignore_case = TRUE)) ~ "medium",
      str_detect(.data$text, regex("figure|table|extended data|supplementary", ignore_case = TRUE)) ~ "medium",
      TRUE ~ "high"
    ),
    review_note = case_when(
      .data$review_priority == "high" ~ "Manuscript-body paragraph with numeric text; check against the linked pipeline number exports before updating prose.",
      .data$review_priority == "medium" ~ "Caption/table-legend numeric text; usually checked against final figure/table artifacts.",
      TRUE ~ "Low-priority numeric text such as affiliation, heading, reference, or administrative material."
    )
  ) %>%
  filter(nzchar(.data$numeric_tokens)) %>%
  select(
    "draft_file",
    "paragraph_index",
    "manuscript_section",
    "manuscript_subsection",
    "paragraph_role",
    "review_priority",
    "numeric_tokens",
    "candidate_source_hint",
    "review_note",
    "text"
  )

write_tsv(draft_index, file.path(out_dir, "manuscript_draft_paragraph_index.tsv"))
write_tsv(numeric_paragraph_review, file.path(out_dir, "manuscript_numeric_paragraph_review.tsv"))

################################################################################
## Cohort and sample availability numbers
################################################################################

clinical_path <- project_path("combined_clinical_data_updated_April2025.csv")
clinical <- read_csv_if_exists(clinical_path)
if (!is.null(clinical)) {
  add_rows(new_row(
    "Results",
    "Cohort and sample overview",
    "clinical_unique_patients",
    "Unique patients in harmonized clinical metadata",
    safe_n_distinct(clinical$Patient),
    fmt_int(safe_n_distinct(clinical$Patient)),
    source_file = clinical_path,
    source_script = "Scripts_2025/Final_Scripts/1_0_Process_clinical_metadata.R",
    related_artifacts = "Figure_1A; Supplementary_Table_1",
    suggested_text = paste0("The harmonized clinical metadata currently includes ", fmt_int(safe_n_distinct(clinical$Patient)), " unique patients.")
  ))
  add_rows(new_row(
    "Results",
    "Cohort and sample overview",
    "clinical_total_samples",
    "Total samples in harmonized clinical metadata",
    nrow(clinical),
    fmt_int(nrow(clinical)),
    source_file = clinical_path,
    source_script = "Scripts_2025/Final_Scripts/1_0_Process_clinical_metadata.R",
    related_artifacts = "Figure_1A; Figure_1B; Supplementary_Table_1"
  ))

  if ("Study" %in% names(clinical)) {
    add_rows(clinical %>%
      filter(!is.na(.data$Study), nzchar(.data$Study)) %>%
      group_by(.data$Study) %>%
      summarise(
        patients = n_distinct(.data$Patient),
        samples = n(),
        .groups = "drop"
      ) %>%
      transmute(row = purrr::pmap(
        list(.data$Study, .data$patients, .data$samples),
        ~ new_row(
          "Results",
          "Cohort and sample overview",
          paste0("clinical_study_", make.names(..1)),
          paste0("Patients and samples from ", ..1),
          paste0(..2, " patients; ", ..3, " samples"),
          paste0(fmt_int(..2), " patients; ", fmt_int(..3), " samples"),
          numerator = ..2,
          denominator = ..3,
          units = "patients; samples",
          source_file = clinical_path,
          source_script = "Scripts_2025/Final_Scripts/1_0_Process_clinical_metadata.R",
          related_artifacts = "Figure_1A; Figure_1B"
        )
      )) %>%
      pull(.data$row) %>%
      bind_rows())
  }

  if ("Sample_type" %in% names(clinical)) {
    add_rows(clinical %>%
      filter(!is.na(.data$Sample_type), nzchar(.data$Sample_type)) %>%
      count(.data$Sample_type, name = "samples") %>%
      transmute(row = purrr::pmap(
        list(.data$Sample_type, .data$samples),
        ~ new_row(
          "Results",
          "Cohort and sample overview",
          paste0("clinical_sample_type_", make.names(..1)),
          paste0("Samples with type ", ..1),
          ..2,
          fmt_int(..2),
          denominator = nrow(clinical),
          units = "samples",
          source_file = clinical_path,
          source_script = "Scripts_2025/Final_Scripts/1_0_Process_clinical_metadata.R",
          related_artifacts = "Figure_1A; Figure_1B"
        )
      )) %>%
      pull(.data$row) %>%
      bind_rows())
  }
}

cohort_path <- project_path("cohort_assignment_table_updated.rds")
cohort_df <- read_rds_if_exists(cohort_path)
if (!is.null(cohort_df) && all(c("Patient", "Cohort") %in% names(cohort_df))) {
  add_rows(cohort_df %>%
    filter(!is.na(.data$Cohort), nzchar(.data$Cohort)) %>%
    count(.data$Cohort, name = "patients") %>%
    transmute(row = purrr::pmap(
      list(.data$Cohort, .data$patients),
      ~ new_row(
        "Results",
        "Cohort definition",
        paste0("assigned_cohort_", make.names(..1)),
        paste0("Patients assigned to ", ..1, " cohort"),
        ..2,
        fmt_int(..2),
        units = "patients",
        source_file = cohort_path,
        source_script = "Scripts_2025/Final_Scripts/1_6_Identify_High_Quality_Patient_Pairs.R",
        related_artifacts = "Figure_1B"
      )
    )) %>%
    pull(.data$row) %>%
    bind_rows())
}

################################################################################
## cfWGS/clinical call table numbers
################################################################################

calls_path <- project_path("Output_tables_2025", "all_patients_with_BM_and_blood_calls_updated5.rds")
calls <- read_rds_if_exists(calls_path)
if (!is.null(calls)) {
  add_rows(new_row(
    "Results",
    "Integrated assay call table",
    "call_table_unique_patients",
    "Unique patients in integrated cfWGS/clinical call table",
    safe_n_distinct(calls$Patient),
    fmt_int(safe_n_distinct(calls$Patient)),
    source_file = calls_path,
    source_script = "Scripts_2025/Final_Scripts/3_2_Plot_optimal_cutoff_and_clinical_concordance.R",
    related_artifacts = "Figure_3D; Figure_3E; Figure_4C; Figure_4D; Supplementary_Table_8; Supplementary_Table_10"
  ))
  add_rows(new_row(
    "Results",
    "Integrated assay call table",
    "call_table_total_samples",
    "Rows/samples in integrated cfWGS/clinical call table",
    nrow(calls),
    fmt_int(nrow(calls)),
    source_file = calls_path,
    source_script = "Scripts_2025/Final_Scripts/3_2_Plot_optimal_cutoff_and_clinical_concordance.R",
    related_artifacts = "Supplementary_Table_8; Supplementary_Table_10"
  ))

  assay_cols <- c(
    MFC = "Flow_Binary",
    clonoSEQ = "Adaptive_Binary",
    EasyM = "EasyM_optimized_binary",
    BM_cfWGS = "BM_zscore_only_detection_rate_call",
    Blood_cfWGS_sites = "Blood_zscore_only_sites_call",
    Blood_cfWGS_combined = "Blood_plus_fragment_call"
  )
  assay_cols <- assay_cols[assay_cols %in% names(calls)]
  add_rows(purrr::imap_dfr(assay_cols, function(col, label) {
    n_samples <- sum(!is.na(calls[[col]]))
    n_patients <- safe_n_distinct(calls$Patient[!is.na(calls[[col]])])
    n_positive <- sum(calls[[col]] == 1, na.rm = TRUE)
    bind_rows(
      new_row(
        "Results",
        "Integrated assay call table",
        paste0("assay_available_", label),
        paste0(label, " available calls"),
        paste0(n_samples, " samples; ", n_patients, " patients"),
        paste0(fmt_int(n_samples), " samples; ", fmt_int(n_patients), " patients"),
        numerator = n_samples,
        denominator = n_patients,
        units = "samples; patients",
        source_file = calls_path,
        source_script = "Scripts_2025/Final_Scripts/3_2_Plot_optimal_cutoff_and_clinical_concordance.R",
        related_artifacts = "Figure_3D; Figure_3E; Figure_4C; Figure_4D; Supplementary_Table_8; Supplementary_Table_10"
      ),
      new_row(
        "Results",
        "Integrated assay call table",
        paste0("assay_positive_", label),
        paste0(label, " positive calls among available calls"),
        n_positive,
        paste0(fmt_int(n_positive), " / ", fmt_int(n_samples), " (", fmt_pct(n_positive / n_samples, 1), ")"),
        numerator = n_positive,
        denominator = n_samples,
        units = "positive calls / available calls",
        source_file = calls_path,
        source_script = "Scripts_2025/Final_Scripts/3_2_Plot_optimal_cutoff_and_clinical_concordance.R",
        related_artifacts = "Figure_3D; Figure_3E; Figure_4C; Figure_4D; Supplementary_Table_8; Supplementary_Table_10"
      )
    )
  }))

  if ("Cohort" %in% names(calls)) {
    add_rows(calls %>%
      filter(!is.na(.data$Cohort), nzchar(.data$Cohort)) %>%
      group_by(.data$Cohort) %>%
      summarise(samples = n(), patients = n_distinct(.data$Patient), .groups = "drop") %>%
      transmute(row = purrr::pmap(
        list(.data$Cohort, .data$samples, .data$patients),
        ~ new_row(
          "Results",
          "Integrated assay call table",
          paste0("call_table_cohort_", make.names(..1)),
          paste0("Integrated call-table samples in ", ..1, " cohort"),
          paste0(..2, " samples; ", ..3, " patients"),
          paste0(fmt_int(..2), " samples; ", fmt_int(..3), " patients"),
          numerator = ..2,
          denominator = ..3,
          units = "samples; patients",
          source_file = calls_path,
          source_script = "Scripts_2025/Final_Scripts/3_2_Plot_optimal_cutoff_and_clinical_concordance.R",
          related_artifacts = "Figure_3; Figure_4; Extended_Data_Figure_5; Extended_Data_Figure_7"
        )
      )) %>%
      pull(.data$row) %>%
      bind_rows())
  }
}

################################################################################
## Model performance tables
################################################################################

model_metric_files <- tibble::tribble(
  ~section, ~subsection, ~metric_group, ~path, ~source_script, ~artifact,
  "Results", "BM cfWGS model performance", "Supplementary Table 4 nested-CV metrics",
  project_path("Figures_Exported", "Tables_exported", "Renamed", "Supplementary_Table_4_All_Model_performance_nested_CV.csv"),
  "Scripts_2025/Final_Scripts/3_1_Optimize_cfWGS_thresholds.R", "Figure_3A; Supplementary_Table_4",
  "Results", "BM/blood refit model metrics", "Supplementary Table 5 refit metrics",
  project_path("Final Tables and Figures", "Supplementary_Table_4_All_Model_Metrics_Refit6_Feb2026.csv"),
  "Scripts_2025/Final_Scripts/3_1_Optimize_cfWGS_thresholds.R", "Figure_3B; Figure_4B; Supplementary_Table_5",
  "Results", "Blood cfWGS model performance", "Supplementary Table 6 test-cohort metrics",
  project_path("Figures_Exported", "Tables_exported", "Renamed", "Supplementary_Table_6_All_Model_performance_nested_CV_on_test_cohort.csv"),
  "Scripts_2025/Final_Scripts/3_1_Optimize_cfWGS_thresholds.R", "Figure_4A; Supplementary_Table_6"
)

for (i in seq_len(nrow(model_metric_files))) {
  info <- model_metric_files[i, ]
  tbl <- read_csv_if_exists(info$path)
  if (is.null(tbl)) next
  add_rows(new_row(
    info$section,
    info$subsection,
    paste0("model_metric_table_rows_", make.names(info$metric_group)),
    paste0(info$metric_group, " row count"),
    nrow(tbl),
    fmt_int(nrow(tbl)),
    units = "rows",
    source_file = info$path,
    source_script = info$source_script,
    related_artifacts = info$artifact,
    caveat = "Use the underlying table for exact model-specific AUC, sensitivity, specificity, PPV, and NPV values."
  ))

  auc_col <- names(tbl)[str_detect(tolower(names(tbl)), "^auc$|auc")]
  model_col <- names(tbl)[str_detect(tolower(names(tbl)), "model|assay|feature")]
  if (length(auc_col) && length(model_col)) {
    auc_col <- auc_col[1]
    model_col <- model_col[1]
    tbl_auc <- tbl %>%
      mutate(.auc = suppressWarnings(as.numeric(.data[[auc_col]]))) %>%
      filter(!is.na(.data$.auc)) %>%
      arrange(desc(.data$.auc)) %>%
      slice_head(n = 5)
    if (nrow(tbl_auc)) {
      add_rows(tbl_auc %>%
        transmute(row = purrr::pmap(
          list(.data[[model_col]], .data$.auc),
          ~ new_row(
            info$section,
            info$subsection,
            paste0("top_auc_", make.names(info$metric_group), "_", make.names(as.character(..1))),
            paste0("AUC for ", as.character(..1)),
            ..2,
            fmt_num(..2, 3),
            units = "AUC",
            source_file = info$path,
            source_script = info$source_script,
            related_artifacts = info$artifact
          )
        )) %>%
        pull(.data$row) %>%
        bind_rows())
    }
  }
}

################################################################################
## Time-window relapse/progression performance
################################################################################

timewindow_sources <- tibble::tribble(
  ~subsection, ~path, ~source_script, ~artifact, ~caveat,
  "Submitted BM time-window analysis",
  project_path("Output_tables_2025", "detection_progression_updated5", "BM_cfWGS_timewindow_results2.csv"),
  "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R", "Extended_Data_Figure_6I; Supplementary_Table_9",
  "Submitted manuscript values; preserved cached source.",
  "Submitted blood time-window analysis",
  project_path("Output_tables_2025", "detection_progression_updated5", "blood_cfWGS_timewindow_results2.csv"),
  "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R", "Extended_Data_Figure_8D; Supplementary_Table_9",
  "Submitted manuscript values; preserved cached source.",
  "Prospective BM time-window QC",
  project_path("Output_tables_2025", "detection_progression_updated6", "prospective_timewindow_qc", "prospective_BM_cfWGS_timewindow_results.csv"),
  "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R", "prospective_QC",
  "Prospective update rule for newly added test-cohort samples.",
  "Prospective blood time-window QC",
  project_path("Output_tables_2025", "detection_progression_updated6", "prospective_timewindow_qc", "prospective_blood_cfWGS_timewindow_results.csv"),
  "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R", "prospective_QC",
  "Prospective update rule for newly added test-cohort samples."
)

for (i in seq_len(nrow(timewindow_sources))) {
  info <- timewindow_sources[i, ]
  tbl <- read_csv_if_exists(info$path)
  if (is.null(tbl) || !all(c("Window_days", "Assay", "N_samples", "N_patients", "Sensitivity") %in% names(tbl))) next
  add_rows(tbl %>%
    transmute(row = purrr::pmap(
      list(.data$Window_days, .data$Assay, .data$N_samples, .data$N_patients, .data$Sensitivity),
      ~ new_row(
        "Results",
        info$subsection,
        paste0("timewindow_", make.names(info$subsection), "_", ..1, "_", make.names(as.character(..2))),
        paste0(as.character(..2), " sensitivity within ", ..1, " days"),
        ..5,
        paste0(fmt_pct(..5, 0), " sensitivity; n=", fmt_int(..3), " samples / ", fmt_int(..4), " patients"),
        numerator = ..3,
        denominator = ..4,
        units = "sensitivity; samples; patients",
        source_file = info$path,
        source_script = info$source_script,
        related_artifacts = info$artifact,
        caveat = info$caveat
      )
    )) %>%
    pull(.data$row) %>%
    bind_rows())
}

prospective_labels_path <- project_path("Output_tables_2025", "detection_progression_updated6", "prospective_timewindow_qc", "prospective_timewindow_sample_level_labels.csv")
prospective_labels <- read_csv_if_exists(prospective_labels_path)
if (!is.null(prospective_labels) && all(c("Patient", "sample_row_id", "Window_days", "evaluable_in_window") %in% names(prospective_labels))) {
  add_rows(prospective_labels %>%
    group_by(.data$Window_days) %>%
    summarise(
      samples_total = n_distinct(.data$sample_row_id),
      samples_evaluable = n_distinct(.data$sample_row_id[.data$evaluable_in_window]),
      patients_total = n_distinct(.data$Patient),
      .groups = "drop"
    ) %>%
    transmute(row = purrr::pmap(
      list(.data$Window_days, .data$samples_total, .data$samples_evaluable, .data$patients_total),
      ~ new_row(
        "Results",
        "Prospective time-window QC",
        paste0("prospective_evaluable_", ..1, "d"),
        paste0("Samples evaluable for prospective ", ..1, "-day window"),
        ..3,
        paste0(fmt_int(..3), " / ", fmt_int(..2), " samples evaluable across ", fmt_int(..4), " patients"),
        numerator = ..3,
        denominator = ..2,
        units = "evaluable samples / total samples",
        source_file = prospective_labels_path,
        source_script = "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R",
        related_artifacts = "prospective_QC",
        caveat = "Prospective QC output, not the submitted cached manuscript table."
      )
    )) %>%
    pull(.data$row) %>%
    bind_rows())
}

followup_path <- project_path("Exported_data_tables_clinical", "patient_followup_dates_updated.csv")
followup <- read_csv_if_exists(followup_path)
if (!is.null(followup) && all(c("Patient", "followup_end_date") %in% names(followup))) {
  add_rows(new_row(
    "Methods",
    "Clinical follow-up date handling",
    "patient_followup_table_rows",
    "Patients with explicit last-known follow-up dates",
    nrow(followup),
    fmt_int(nrow(followup)),
    units = "patients",
    source_file = followup_path,
    source_script = "Scripts_2025/Final_Scripts/1_0_Process_clinical_metadata.R",
    related_artifacts = "prospective_QC",
    suggested_text = paste0("Prospective time-window evaluability uses a separate patient-level last-follow-up table with ", fmt_int(nrow(followup)), " patients, rather than overwriting PFS event dates.")
  ))
}

################################################################################
## Final artifact/provenance helper rows
################################################################################

if (!is.null(artifact_map) && nrow(artifact_map)) {
  add_rows(new_row(
    "Supplementary information",
    "Manuscript output organization",
    "mapped_manuscript_artifacts",
    "Mapped final manuscript figure/table artifacts",
    nrow(artifact_map),
    fmt_int(nrow(artifact_map)),
    units = "artifact-map rows",
    source_file = project_path("docs", "manuscript_artifact_source_map.tsv"),
    source_script = "Scripts_2025/Final_Scripts/build_stage_artifact_map.R",
    related_artifacts = "all mapped figures/tables"
  ))
}

################################################################################
## Dynamic manuscript-claim helper rows
################################################################################

sample_flow_path <- project_path("Final Tables and Figures", "Table for creating sample flowchart updated3.csv")
sample_flow <- read_csv_if_exists(sample_flow_path)
if (!is.null(sample_flow) && all(c("Patient", "total_cfDNA_samples") %in% names(sample_flow))) {
  add_rows(bind_rows(
    new_row(
      "Abstract",
      "Cohort summary claims",
      "current_sample_flow_unique_patients",
      "Unique patients in Figure 1B sample-flow source table",
      safe_n_distinct(sample_flow$Patient),
      fmt_int(safe_n_distinct(sample_flow$Patient)),
      units = "patients",
      source_file = sample_flow_path,
      source_script = "Scripts_2025/Final_Scripts/1_6_Identify_High_Quality_Patient_Pairs.R",
      related_artifacts = "Figure_1B",
      update_trigger = "rerun sample-flow/cohort scripts or add samples",
      suggested_text = paste0("Current Figure 1B sample-flow source table includes ", fmt_int(safe_n_distinct(sample_flow$Patient)), " patients.")
    ),
    new_row(
      "Abstract",
      "Cohort summary claims",
      "current_sample_flow_total_cfdna_samples",
      "Total cfDNA samples in Figure 1B sample-flow source table",
      sum(sample_flow$total_cfDNA_samples, na.rm = TRUE),
      fmt_int(sum(sample_flow$total_cfDNA_samples, na.rm = TRUE)),
      units = "cfDNA samples",
      source_file = sample_flow_path,
      source_script = "Scripts_2025/Final_Scripts/1_6_Identify_High_Quality_Patient_Pairs.R",
      related_artifacts = "Figure_1B",
      update_trigger = "rerun sample-flow/cohort scripts or add samples",
      caveat = "This is the current source-table value and may not match older prose if sample inclusion has changed."
    ),
    new_row(
      "Abstract",
      "Cohort summary claims",
      "current_sample_flow_high_quality_bm_patients",
      "Patients with high-quality BM in Figure 1B sample-flow source table",
      sum(sample_flow$high_quality_BM, na.rm = TRUE),
      fmt_int(sum(sample_flow$high_quality_BM, na.rm = TRUE)),
      units = "patients",
      source_file = sample_flow_path,
      source_script = "Scripts_2025/Final_Scripts/1_6_Identify_High_Quality_Patient_Pairs.R",
      related_artifacts = "Figure_1B",
      update_trigger = "rerun sample-flow/cohort scripts or add samples"
    ),
    new_row(
      "Abstract",
      "Cohort summary claims",
      "current_sample_flow_high_quality_cfdna_patients",
      "Patients with high-quality cfDNA in Figure 1B sample-flow source table",
      sum(sample_flow$high_quality_cfDNA, na.rm = TRUE),
      fmt_int(sum(sample_flow$high_quality_cfDNA, na.rm = TRUE)),
      units = "patients",
      source_file = sample_flow_path,
      source_script = "Scripts_2025/Final_Scripts/1_6_Identify_High_Quality_Patient_Pairs.R",
      related_artifacts = "Figure_1B",
      update_trigger = "rerun sample-flow/cohort scripts or add samples"
    )
  ))
}

bm_maf_path <- project_path("combined_maf_bm_dx.rds")
bm_maf <- read_rds_if_exists(bm_maf_path)
if (!is.null(bm_maf) && "Patient" %in% names(bm_maf)) {
  bm_mut_counts <- bm_maf %>%
    count(.data$Patient, name = "n_mutations")
  add_rows(bind_rows(
    new_row(
      "Abstract",
      "Genomic-feature summary claims",
      "bm_dx_mutation_count_median",
      "Median BM diagnostic mutations per patient in available BM MAF",
      median(bm_mut_counts$n_mutations, na.rm = TRUE),
      fmt_int(median(bm_mut_counts$n_mutations, na.rm = TRUE)),
      units = "mutations per patient",
      source_file = bm_maf_path,
      source_script = "Scripts_2025/Final_Scripts/1_2_Process_Mutation_Data.R",
      related_artifacts = "Extended_Data_Figure_2; Supplementary_Table_2",
      update_trigger = "rerun mutation-processing scripts or add BM mutation data",
      caveat = "Computed from the current BM diagnostic MAF object; compare to manuscript prose before replacing older retained text values."
    ),
    new_row(
      "Abstract",
      "Genomic-feature summary claims",
      "bm_dx_mutation_count_patients",
      "Patients contributing BM diagnostic mutation counts",
      nrow(bm_mut_counts),
      fmt_int(nrow(bm_mut_counts)),
      units = "patients",
      source_file = bm_maf_path,
      source_script = "Scripts_2025/Final_Scripts/1_2_Process_Mutation_Data.R",
      related_artifacts = "Extended_Data_Figure_2; Supplementary_Table_2",
      update_trigger = "rerun mutation-processing scripts or add BM mutation data"
    )
  ))
}

bm_nested_path <- project_path("Figures_Exported", "Tables_exported", "Renamed", "Supplementary_Table_4_All_Model_performance_nested_CV.csv")
bm_nested <- read_csv_if_exists(bm_nested_path)
if (!is.null(bm_nested) && all(c("combo", "auc_mean") %in% names(bm_nested))) {
  selected <- bm_nested %>%
    filter(.data$combo %in% c("BM_zscore_only_sites", "Blood_plus_fragment", "Blood_rate_only", "Blood_base")) %>%
    arrange(match(.data$combo, c("BM_zscore_only_sites", "Blood_plus_fragment", "Blood_rate_only", "Blood_base")))
  if (nrow(selected)) {
    add_rows(selected %>%
      transmute(row = purrr::pmap(
        list(.data$combo, .data$auc_mean),
        ~ new_row(
          "Abstract",
          "Model-performance summary claims",
          paste0("nested_cv_auc_mean_", make.names(..1)),
          paste0("Mean nested-CV AUC for ", ..1),
          ..2,
          fmt_num(..2, 3),
          units = "AUC",
          source_file = bm_nested_path,
          source_script = "Scripts_2025/Final_Scripts/3_1_Optimize_cfWGS_thresholds.R",
          related_artifacts = "Figure_3A; Figure_4A; Supplementary_Table_4",
          update_trigger = "rerun cache-sensitive model outputs only when explicitly requested",
          caveat = "Model-training output is cache-sensitive; routine test-cohort expansion should not refit this value."
        )
      )) %>%
      pull(.data$row) %>%
      bind_rows())
  }
}

bm_surv_path <- project_path("Output_tables_2025", "detection_progression_updated6", "cfWGS_vs_flow_progression_summary_2026-06-22.csv")
bm_surv <- read_csv_if_exists(bm_surv_path)
if (!is.null(bm_surv) && all(c("Landmark", "HR_cf", "MedRFS_cf_pos", "Patients", "Events") %in% names(bm_surv))) {
  bm_1yr <- bm_surv %>% filter(.data$Landmark == "1yr_maintenance") %>% slice_head(n = 1)
  if (nrow(bm_1yr)) {
    add_rows(bind_rows(
      new_row(
        "Abstract",
        "Survival/progression summary claims",
        "bm_cfwgs_1yr_maintenance_hr",
        "BM-informed cfWGS hazard ratio at 1-year maintenance",
        bm_1yr$HR_cf,
        fmt_num(bm_1yr$HR_cf, 2),
        units = "hazard ratio",
        source_file = bm_surv_path,
        source_script = "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R",
        related_artifacts = "Figure_3F; Extended_Data_Figure_6",
        update_trigger = "rerun survival outputs or add outcome/sample data"
      ),
      new_row(
        "Abstract",
        "Survival/progression summary claims",
        "bm_cfwgs_1yr_maintenance_median_rfs_positive_months",
        "Median RFS among BM-informed cfWGS-positive patients at 1-year maintenance",
        bm_1yr$MedRFS_cf_pos,
        fmt_num(bm_1yr$MedRFS_cf_pos, 1),
        units = "months",
        source_file = bm_surv_path,
        source_script = "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R",
        related_artifacts = "Figure_3F; Extended_Data_Figure_6",
        update_trigger = "rerun survival outputs or add outcome/sample data"
      ),
      new_row(
        "Abstract",
        "Survival/progression summary claims",
        "bm_cfwgs_1yr_maintenance_patients_events",
        "Patients and events in BM-informed 1-year maintenance survival summary",
        paste0(bm_1yr$Patients, " patients; ", bm_1yr$Events, " events"),
        paste0(fmt_int(bm_1yr$Patients), " patients; ", fmt_int(bm_1yr$Events), " events"),
        numerator = bm_1yr$Events,
        denominator = bm_1yr$Patients,
        units = "events / patients",
        source_file = bm_surv_path,
        source_script = "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R",
        related_artifacts = "Figure_3F; Extended_Data_Figure_6",
        update_trigger = "rerun survival outputs or add outcome/sample data"
      )
    ))
  }
}

blood_surv_path <- project_path("Output_tables_2025", "detection_progression_updated6", "cfWGS_vs_flow_progression_summary_blood_muts_2026-06-22.csv")
blood_surv <- read_csv_if_exists(blood_surv_path)
if (!is.null(blood_surv) && all(c("Landmark", "HR_cf", "MedRFS_cf_pos", "Patients", "Events") %in% names(blood_surv))) {
  blood_1yr <- blood_surv %>% filter(.data$Landmark == "1yr_maintenance") %>% slice_head(n = 1)
  if (nrow(blood_1yr)) {
    add_rows(bind_rows(
      new_row(
        "Abstract",
        "Survival/progression summary claims",
        "blood_cfwgs_1yr_maintenance_hr",
        "Blood-derived cfWGS hazard ratio at 1-year maintenance",
        blood_1yr$HR_cf,
        fmt_num(blood_1yr$HR_cf, 2),
        units = "hazard ratio",
        source_file = blood_surv_path,
        source_script = "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R",
        related_artifacts = "Figure_4E; Extended_Data_Figure_8",
        update_trigger = "rerun survival outputs or add outcome/sample data"
      ),
      new_row(
        "Abstract",
        "Survival/progression summary claims",
        "blood_cfwgs_1yr_maintenance_median_rfs_positive_months",
        "Median RFS among blood-derived cfWGS-positive patients at 1-year maintenance",
        blood_1yr$MedRFS_cf_pos,
        fmt_num(blood_1yr$MedRFS_cf_pos, 1),
        units = "months",
        source_file = blood_surv_path,
        source_script = "Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R",
        related_artifacts = "Figure_4E; Extended_Data_Figure_8",
        update_trigger = "rerun survival outputs or add outcome/sample data"
      )
    ))
  }
}

revision_docx_path <- first_existing(project_path("cfWGS Manuscript Drafts", c(
  "cfWGS_MRD_MM_Manuscript_Revision1.docx",
  "Manuscript_ Cell-free DNA Whole Genome Sequencing for Non-Invasive MRD Detection in Multiple Myeloma.docx"
)))
add_rows(bind_rows(
  new_row(
    "Summary Paragraph",
    "Retained lead-time prose checks",
    "retained_blood_derived_lead_time_months_prose",
    "Retained prose value for blood-derived cfWGS lead time before progression",
    19,
    "19 months",
    units = "months",
    source_file = revision_docx_path,
    source_script = "Scripts_2025/Final_Scripts/5_0_Build_Manuscript_Text_Number_Exports.R",
    related_artifacts = "manuscript prose audit",
    update_trigger = "check against updated progression/survival outputs before manuscript revision",
    caveat = "Retained manuscript prose value, not recalculated by this helper. Use source_file_hints in manuscript_paragraph_metric_audit.tsv to confirm whether this should be updated after adding samples.",
    review_status = "requires_manual_confirmation"
  ),
  new_row(
    "Nature Summary Paragraph",
    "Retained lead-time prose checks",
    "retained_blood_derived_lead_time_days_prose",
    "Retained prose value for blood-derived cfWGS lead time before progression",
    385,
    "385 days",
    units = "days",
    source_file = revision_docx_path,
    source_script = "Scripts_2025/Final_Scripts/5_0_Build_Manuscript_Text_Number_Exports.R",
    related_artifacts = "manuscript prose audit",
    update_trigger = "check against updated progression/survival outputs before manuscript revision",
    caveat = "Retained manuscript prose value, not recalculated by this helper. Use source_file_hints in manuscript_paragraph_metric_audit.tsv to confirm whether this should be updated after adding samples.",
    review_status = "requires_manual_confirmation"
  )
))

dilution_path <- project_path("Output_tables_2025", "Supp_Table_7_dilution_series_scored_updated3.csv")
dilution <- read_csv_if_exists(dilution_path)
if (!is.null(dilution) && "LOD_updated" %in% names(dilution)) {
  dilution_rows <- list(new_row(
      "Discussion",
      "Analytical sensitivity claims",
      "dilution_series_lowest_lod_updated",
      "Lowest non-missing updated LOD in dilution-series scored table",
      min(dilution$LOD_updated, na.rm = TRUE),
      paste0(fmt_num(min(dilution$LOD_updated, na.rm = TRUE), 4), "%"),
      units = "aggregate allele fraction percent",
      source_file = dilution_path,
      source_script = "Scripts_2025/Final_Scripts/3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R",
      related_artifacts = "Figure_3C; Extended_Data_Figure_5D; Extended_Data_Figure_7D; Supplementary_Table_7",
      update_trigger = "rerun dilution-series scoring",
      caveat = "Current source-table value; compare with manuscript prose if retained text uses the lowest tested dilution rather than the minimum LOD_updated field."
    ))
  if ("LOD_original" %in% names(dilution)) {
    tested_dilution <- dilution %>%
      filter(!is.na(.data$LOD_updated), !is.na(.data$LOD_original), .data$LOD_original > 0) %>%
      arrange(.data$LOD_updated) %>%
      slice_head(n = 1)
    if (nrow(tested_dilution)) {
      dilution_rows[[length(dilution_rows) + 1]] <- new_row(
        "Discussion",
        "Analytical sensitivity claims",
        "dilution_series_lowest_tested_dilution_lod_updated",
        "Lowest non-zero dilution-series updated LOD",
        tested_dilution$LOD_updated,
        paste0(fmt_num(tested_dilution$LOD_updated, 3), "%"),
        units = "aggregate allele fraction percent",
        source_file = dilution_path,
        source_script = "Scripts_2025/Final_Scripts/3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R",
        related_artifacts = "Figure_3C; Extended_Data_Figure_5D; Extended_Data_Figure_7D; Supplementary_Table_7",
        update_trigger = "rerun dilution-series scoring",
        caveat = "Matches prose describing the lowest dilution tested, excluding the undiluted/neat positive sample."
      )
    }
  }
  add_rows(bind_rows(dilution_rows))
}

################################################################################
## Write outputs
################################################################################

numbers <- bind_rows(rows) %>%
  mutate(
    manuscript_section = replace_na(.data$manuscript_section, ""),
    manuscript_subsection = replace_na(.data$manuscript_subsection, ""),
    metric_id = replace_na(.data$metric_id, ""),
    source_exists = ifelse(nzchar(.data$source_file), file.exists(.data$source_file), NA),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ) %>%
  arrange(.data$manuscript_section, .data$manuscript_subsection, .data$metric_id)

infer_metric_source_category <- function(manuscript_subsection, source_script, related_artifacts) {
  txt <- str_to_lower(paste(manuscript_subsection, source_script, related_artifacts, collapse = " "))
  categories <- character()
  if (str_detect(txt, "clinical_metadata|cohort|sample|figure_1|table_1|supplementary_table_1")) {
    categories <- c(categories, "cohort/sample-count outputs")
  }
  if (str_detect(txt, "optimize_cfwgs|model|nested|auc|threshold|supplementary_table_4|supplementary_table_5|supplementary_table_6|figure_3a|figure_3b|figure_4a|figure_4b")) {
    categories <- c(categories, "model-performance tables")
  }
  if (str_detect(txt, "survival|time-window|progression|follow-up|supplementary_table_9|extended_data_figure_6|extended_data_figure_8")) {
    categories <- c(categories, "survival/time-window outputs")
  }
  if (str_detect(txt, "retained lead-time prose|lead.time|lead_time|manuscript prose audit")) {
    categories <- c(categories, "survival/time-window outputs")
  }
  if (str_detect(txt, "concordance|clinical_concordance|supplementary_table_8|supplementary_table_10|figure_3d|figure_3e|figure_4c|figure_4d")) {
    categories <- c(categories, "clinical concordance outputs")
  }
  if (str_detect(txt, "feature|mutation|cna|fish|translocation|extended_data_figure_1|extended_data_figure_2|supplementary_table_2|supplementary_table_3")) {
    categories <- c(categories, "genomic feature/concordance outputs")
  }
  if (str_detect(txt, "dilution|lod|figure_3c|supplementary_table_7")) {
    categories <- c(categories, "dilution-series outputs")
  }
  if (!length(categories)) categories <- "manual review"
  paste(unique(categories), collapse = "; ")
}

source_file_hints_for_category <- function(candidate_source_hint) {
  categories <- split_semicolon(candidate_source_hint)
  hints <- character()
  if ("cohort/sample-count outputs" %in% categories) {
    hints <- c(
      hints,
      "Final Tables and Figures/Table for creating sample flowchart updated3.csv",
      "combined_clinical_data_updated_April2025.csv",
      "cohort_assignment_table_updated.rds"
    )
  }
  if ("model-performance tables" %in% categories) {
    hints <- c(
      hints,
      "Figures_Exported/Tables_exported/Renamed/Supplementary_Table_4_All_Model_performance_nested_CV.csv",
      "Final Tables and Figures/Supplementary_Table_4_All_Model_Metrics_Refit6_Feb2026.csv",
      "Figures_Exported/Tables_exported/Renamed/Supplementary_Table_6_All_Model_performance_nested_CV_on_test_cohort.csv"
    )
  }
  if ("survival/time-window outputs" %in% categories) {
    hints <- c(
      hints,
      "Output_tables_2025/detection_progression_updated6/cfWGS_vs_flow_progression_summary_2026-06-22.csv",
      "Output_tables_2025/detection_progression_updated6/cfWGS_vs_flow_progression_summary_blood_muts_2026-06-22.csv",
      "Output_tables_2025/detection_progression_updated5/BM_cfWGS_timewindow_results2.csv",
      "Output_tables_2025/detection_progression_updated5/blood_cfWGS_timewindow_results2.csv"
    )
  }
  if ("clinical concordance outputs" %in% categories) {
    hints <- c(
      hints,
      "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated5.rds",
      "Final Tables and Figures/Supplementary_Table_8_model_comparisons_to_clinical_metrics3.xlsx",
      "Figures_Exported/Tables_exported/Renamed/Supplementary_Table_10_All_call_metrics_against_clinical_metrics.xlsx"
    )
  }
  if ("genomic feature/concordance outputs" %in% categories) {
    hints <- c(
      hints,
      "Final Tables and Figures/Extended_Data_Figure_2G_mutation_overlap_summary.csv",
      "Final Tables and Figures/Extended_Data_Figure_2G_mutation_overlap_source_data.csv",
      "Figures_Exported/Tables_exported/Renamed/Supplementary_Table_2_SV_CNA_performance_summary_updated.xlsx",
      "Figures_Exported/Tables_exported/Renamed/Supplementary_Table_3_All_Feature_Correlations.csv"
    )
  }
  if ("dilution-series outputs" %in% categories) {
    hints <- c(
      hints,
      "Output_tables_2025/Supp_Table_7_dilution_series_scored_updated3.csv",
      "Figures_Exported/Tables_exported/Renamed/Supplementary_Table_7_Dilution_Series_Combined.xlsx"
    )
  }
  if (!length(hints)) hints <- "Static manuscript/background/methods text; verify manually against protocol, citation, kit, or software documentation if edited."
  paste(unique(hints), collapse = "; ")
}

infer_update_relevance <- function(manuscript_section, manuscript_subsection, text, candidate_source_hint) {
  txt <- str_to_lower(paste(manuscript_section, manuscript_subsection, text))
  citation_context <- str_detect(
    txt,
    "recent prospective studies|genomic studies at relapse|references for supplementary notes"
  )
  dynamic_hint <- str_detect(
    candidate_source_hint,
    "cohort/sample-count|model-performance|survival/time-window|clinical concordance|genomic feature|dilution"
  )
  static_method <- str_detect(
    txt,
    "category #|cat#|kit|illumina|novaseq|bwa|samtools|gatk|picard|mutect|vep|bcftools|bedtools|tabix|vcf2maf|ichorcna|sequenza|igcaller|griffin|hg38|\\bng\\b|\\bbp\\b|reb #|declaration of helsinki|github|version| v[0-9]|minimum input|dual indexes|run configuration|controls were recruited"
  )
  case_when(
    citation_context ~ "static_or_manual_context",
    static_method ~ "static_methods_or_background_number",
    str_detect(txt, "abstract|summary paragraph|statement of translational significance|results|discussion") & dynamic_hint ~ "dynamic_manuscript_claim_check_on_rerun",
    dynamic_hint ~ "potentially_dynamic_check_on_rerun",
    TRUE ~ "static_or_manual_context"
  )
}

split_semicolon <- function(x) {
  out <- unlist(strsplit(x %||% "", "\\s*;\\s*"))
  out[nzchar(out)]
}

token_candidates_for_matching <- function(tokens) {
  out <- split_semicolon(tokens)
  out <- str_replace_all(out, ",", "")
  out <- out[nzchar(out)]
  keep <- str_detect(out, "[/%\\.]") | suppressWarnings(as.numeric(out) >= 10)
  unique(out[replace_na(keep, FALSE)])
}

metric_lookup <- numbers %>%
  mutate(
    metric_source_category = purrr::pmap_chr(
      list(.data$manuscript_subsection, .data$source_script, .data$related_artifacts),
      infer_metric_source_category
    ),
    metric_search_text = str_to_lower(paste(
      .data$value,
      .data$formatted_value,
      .data$numerator,
      .data$denominator,
      .data$statistic_label,
      .data$suggested_text,
      sep = " | "
    ))
  )

match_paragraph_to_metrics <- function(
  paragraph_section,
  paragraph_subsection,
  numeric_tokens,
  candidate_source_hint
) {
  hint_categories <- split_semicolon(candidate_source_hint)
  tokens <- token_candidates_for_matching(numeric_tokens)

  candidate_metrics <- metric_lookup %>%
    filter(
      purrr::map_lgl(
        .data$metric_source_category,
        ~ any(split_semicolon(.x) %in% hint_categories)
      )
    )

  exact_matches <- candidate_metrics
  if (length(tokens)) {
    exact_matches <- candidate_metrics %>%
      filter(purrr::map_lgl(
        .data$metric_search_text,
        ~ any(str_detect(.x, fixed(str_to_lower(tokens))))
      ))
  } else {
    exact_matches <- exact_matches[0, ]
  }

  status <- case_when(
    nrow(exact_matches) > 0 ~ "exact_value_match_found",
    nrow(candidate_metrics) > 0 ~ "candidate_source_category_found",
    TRUE ~ "needs_manual_metric_row_or_source_table"
  )

  tibble::tibble(
    match_status = status,
    candidate_metric_count = nrow(candidate_metrics),
    exact_value_match_count = nrow(exact_matches),
    candidate_metric_ids = paste(head(candidate_metrics$metric_id, 20), collapse = "; "),
    exact_value_metric_ids = paste(head(exact_matches$metric_id, 20), collapse = "; "),
    exact_value_metric_labels = paste(head(exact_matches$statistic_label, 20), collapse = "; ")
  )
}

paragraph_metric_audit <- numeric_paragraph_review %>%
  mutate(
    update_relevance = purrr::pmap_chr(
      list(
        .data$manuscript_section,
        .data$manuscript_subsection,
        .data$text,
        .data$candidate_source_hint
      ),
      infer_update_relevance
    ),
    source_file_hints = vapply(.data$candidate_source_hint, source_file_hints_for_category, character(1))
  ) %>%
  mutate(
    audit = purrr::pmap(
      list(
        .data$manuscript_section,
        .data$manuscript_subsection,
        .data$numeric_tokens,
        .data$candidate_source_hint
      ),
      match_paragraph_to_metrics
    )
  ) %>%
  tidyr::unnest("audit") %>%
  mutate(
    audit_note = case_when(
      .data$match_status == "exact_value_match_found" ~ "At least one extracted manuscript number appears in a pipeline-derived metric row; inspect exact_value_metric_ids before editing.",
      .data$match_status == "candidate_source_category_found" ~ "A likely source category exists, but no exact value match was found in the current metric rows; inspect candidate_metric_ids and source tables.",
      TRUE ~ "No current pipeline-derived metric row matched this paragraph category; add a metric row or document the source table before updating this paragraph."
    )
  )

numbers_tsv <- file.path(out_dir, "manuscript_numbers_by_section.tsv")
numbers_xlsx <- file.path(out_dir, "manuscript_numbers_by_section.xlsx")
numbers_md <- file.path(out_dir, "manuscript_numbers_by_section.md")
paragraph_audit_tsv <- file.path(out_dir, "manuscript_paragraph_metric_audit.tsv")
readme_md <- file.path(out_dir, "README.md")

write_tsv(numbers, numbers_tsv)
write_tsv(paragraph_metric_audit, paragraph_audit_tsv)

section_dir <- file.path(out_dir, "by_section")
if (dir.exists(section_dir)) unlink(section_dir, recursive = TRUE, force = TRUE)
dir.create(section_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(section_dir, "numbers"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(section_dir, "paragraph_audits"), recursive = TRUE, showWarnings = FALSE)

write_section_files <- function(tbl, subdir, prefix) {
  if (!nrow(tbl)) return(invisible(character()))
  tbl <- tbl %>%
    mutate(
      .section_file = paste0(
        prefix,
        "__",
        slugify(.data$manuscript_section),
        "__",
        slugify(.data$manuscript_subsection),
        ".tsv"
      )
    )
  paths <- character()
  for (file_name in unique(tbl$.section_file)) {
    section_tbl <- tbl %>%
      filter(.data$.section_file == file_name) %>%
      select(-".section_file")
    out_path <- file.path(section_dir, subdir, file_name)
    write_tsv(section_tbl, out_path)
    paths <- c(paths, out_path)
  }
  invisible(paths)
}

number_section_files <- write_section_files(numbers, "numbers", "numbers")
audit_section_files <- write_section_files(paragraph_metric_audit, "paragraph_audits", "paragraph_audit")

section_file_index <- bind_rows(
  tibble::tibble(
    file_type = "pipeline_numbers",
    path = number_section_files
  ),
  tibble::tibble(
    file_type = "paragraph_metric_audit",
    path = audit_section_files
  )
) %>%
  mutate(file_name = basename(.data$path)) %>%
  select("file_type", "file_name", "path")
write_tsv(section_file_index, file.path(section_dir, "section_file_index.tsv"))

writexl::write_xlsx(
  list(
    numbers_by_section = numbers,
    paragraph_metric_audit = paragraph_metric_audit,
    draft_numeric_review = numeric_paragraph_review,
    draft_paragraph_index = draft_index
  ),
  path = numbers_xlsx
)

section_summary <- numbers %>%
  count(.data$manuscript_section, .data$manuscript_subsection, name = "n_metrics") %>%
  arrange(.data$manuscript_section, .data$manuscript_subsection)

md_lines <- c(
  "# Manuscript Numbers By Section",
  "",
  "This file is generated by `Scripts_2025/Final_Scripts/5_0_Build_Manuscript_Text_Number_Exports.R`.",
  "",
  "Use the TSV or XLSX when updating manuscript text after adding samples or rerunning the pipeline. Each row records the manuscript section, statistic, formatted value, source file, source script, related figure/table artifact, and caveat.",
  "",
  paste0("- TSV: `", numbers_tsv, "`"),
  paste0("- XLSX: `", numbers_xlsx, "`"),
  paste0("- Draft paragraph index: `", file.path(out_dir, "manuscript_draft_paragraph_index.tsv"), "`"),
  paste0("- Numeric paragraph review: `", file.path(out_dir, "manuscript_numeric_paragraph_review.tsv"), "`"),
  paste0("- Paragraph-to-metric audit: `", paragraph_audit_tsv, "`"),
  paste0("- Per-section files: `", section_dir, "`"),
  paste0("- Generated rows: ", nrow(numbers)),
  paste0("- Draft paragraphs indexed: ", nrow(draft_index)),
  paste0("- Number-containing draft paragraphs: ", nrow(numeric_paragraph_review)),
  paste0("- Paragraphs with exact metric-value matches: ", sum(paragraph_metric_audit$match_status == "exact_value_match_found")),
  paste0("- Paragraphs with source-category candidates but no exact value match: ", sum(paragraph_metric_audit$match_status == "candidate_source_category_found")),
  paste0("- Paragraphs needing manual metric rows/source-table review: ", sum(paragraph_metric_audit$match_status == "needs_manual_metric_row_or_source_table")),
  paste0("- Dynamic manuscript claims to check on rerun: ", sum(paragraph_metric_audit$update_relevance == "dynamic_manuscript_claim_check_on_rerun")),
  paste0("- Static methods/background numbers: ", sum(paragraph_metric_audit$update_relevance == "static_methods_or_background_number")),
  "",
  "## Section Summary",
  ""
)

for (i in seq_len(nrow(section_summary))) {
  md_lines <- c(
    md_lines,
    paste0(
      "- ",
      section_summary$manuscript_section[i],
      " / ",
      section_summary$manuscript_subsection[i],
      ": ",
      section_summary$n_metrics[i],
      " metrics"
    )
  )
}

md_lines <- c(
  md_lines,
  "",
  "## Notes",
  "",
  "- This is a manuscript-writing helper, not a new analysis layer.",
  "- `manuscript_numbers_by_section.*` contains pipeline-derived values organized by manuscript section.",
  "- `manuscript_numeric_paragraph_review.tsv` contains every number-containing paragraph extracted from the working DOCX drafts and assigns a review priority/source hint.",
  "- `manuscript_paragraph_metric_audit.tsv` links number-containing paragraphs to candidate metric rows where possible and explicitly flags gaps.",
  "- `update_relevance` separates sample/model/outcome-sensitive manuscript claims from static Methods/background numbers.",
  "- `source_file_hints` gives concrete source files to inspect when an exact metric-row value is not available.",
  "- `by_section/` contains smaller TSV files for section-by-section manuscript editing.",
  "- Submitted manuscript time-window values remain separated from prospective QC values.",
  "- Rows marked with caveats should be checked against the linked source table before changing final prose.",
  "- Add new rows here when a manuscript paragraph contains a number that is not yet represented by the pipeline-derived table."
)

writeLines(md_lines, numbers_md)
writeLines(md_lines, readme_md)

cat("Manuscript writing exports written to:\n")
cat("  ", numbers_tsv, "\n", sep = "")
cat("  ", numbers_xlsx, "\n", sep = "")
cat("  ", numbers_md, "\n", sep = "")
cat("  ", paragraph_audit_tsv, "\n", sep = "")
cat("  ", section_dir, "\n", sep = "")
cat("Rows: ", nrow(numbers), "\n", sep = "")
cat("Draft paragraphs indexed: ", nrow(draft_index), "\n", sep = "")
cat("Paragraph audit rows: ", nrow(paragraph_metric_audit), "\n", sep = "")
