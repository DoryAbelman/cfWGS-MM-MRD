# =============================================================================
# Script: prepare_local_inputs.R
# Project: cfWGS MRD detection (M4 / SPORE / IMMAGINE)
#
# How to run from the project root:
#   Rscript Scripts_2025/Final_Scripts/prepare_local_inputs.R --check
#   Rscript Scripts_2025/Final_Scripts/prepare_local_inputs.R --source-dir /path/to/staged_inputs --copy-missing
#
# Script purpose:
#   Check that local input folders/files needed by the command-line manuscript
#   pipeline are present after removing private workstation/cloud-drive fallbacks.
#   When a local source mirror is supplied, this script can copy missing inputs
#   into the project while preserving the relative paths expected by the
#   numbered analysis scripts.
#
# Pipeline role:
#   This script does not perform scientific analysis and does not create
#   manuscript figures or tables. It is a local setup utility; scientific
#   calculations remain in the numbered scripts.
#   `--check` is read-only unless `--report <path>` is supplied.
#
# Reproducibility note:
#   The repository should not depend on user-specific absolute paths. Required
#   protected/de-identified inputs must be staged into the project-relative
#   paths listed below. If files are missing, this script fails loudly with a
#   concise report so the missing local inputs can be copied in before running
#   the analysis.
# =============================================================================

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) y else x
}

message_log <- function(...) {
  message(timestamp(), " | ", paste0(..., collapse = ""))
}

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

parse_flag_value <- function(args, flag, default = NULL) {
  inline <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(inline)) return(sub(paste0("^", flag, "="), "", inline[1]))
  idx <- which(args == flag)
  if (!length(idx)) return(default)
  if (idx[1] == length(args)) stop(flag, " requires a value.", call. = FALSE)
  args[idx[1] + 1]
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  list(
    # `check` is retained as a parsed status flag, but main() always performs
    # the same checks; only copy_missing changes whether inputs are copied.
    check = "--check" %in% args || !("--copy-missing" %in% args),
    copy_missing = "--copy-missing" %in% args,
    source_dir = parse_flag_value(args, "--source-dir"),
    report = parse_flag_value(args, "--report", NULL)
  )
}

expected_input_roots <- function() {
  data.frame(
    relative_path = c(
      "Clinical data",
      "M4_CMRG_Data",
      "Exported_data_tables_clinical",
      "Jan2025_exported_data",
      "MRDetect_output_winter_2025",
      "Results_Fragmentomics",
      "Fragmentomics_data",
      "Aimee additional data",
      "Output_EasyM_MRD_analysis_2025",
      "Output_tables_2025",
      "Final Tables and Figures"
    ),
    kind = "directory",
    required = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    category = c(
      "clinical_inputs",
      "wgs_processing_inputs",
      "clinical_processed_inputs",
      "wgs_processed_inputs",
      "mrdetect_inputs",
      "fragmentomics_inputs",
      "fragmentomics_raw_inputs",
      "easym_raw_inputs",
      "easym_inputs",
      "preserved_model_and_analysis_inputs",
      "mixed_inputs_and_manuscript_working_outputs"
    ),
    reason = c(
      "Clinical metadata, treatment, relapse, FISH, and cohort files used by Stage 0-2 scripts.",
      "WGS processing logs and source tables used by early feature-processing scripts.",
      "Processed clinical event/censor-date tables used by swimplot and survival scripts.",
      "Processed mutation/CNA/translocation/feature objects used by Stage 2-4 scripts.",
      "MRDetect processed outputs and healthy-reference objects used by MRD call processing.",
      "Fragmentomics metric tables used by feature integration and discordance summaries.",
      "Raw/intermediate fragmentomics files used when recomputing fragmentomics from upstream inputs.",
      "Raw EasyM/proteomic MRD collaborator CSVs used by 3_1_A to rebuild processed EasyM calls.",
      "EasyM/proteomic MRD processed inputs used by Stage 3-4 comparison scripts.",
      "Preserved model calls, metrics, thresholds, and downstream analysis objects used by default manuscript runs.",
      "Historical working outputs that are also read as inputs by later manuscript scripts."
    ),
    stringsAsFactors = FALSE
  )
}

critical_input_files <- function() {
  data.frame(
    relative_path = c(
      "combined_clinical_data_updated_April2025.csv",
      "cohort_assignment_table_updated.rds",
      "id_map.rds",
      "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds",
      "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6.rds",
      "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6_full.rds",
      "Aimee additional data/RAPID NOVOR VALUES with values with relapse.csv",
      "Aimee additional data/RAPID NOVOR pos-neg with relapse.csv",
      "Output_EasyM_MRD_analysis_2025/EasyM_all_samples_with_optimized_calls.csv",
      "Output_EasyM_MRD_analysis_2025/EasyM_threshold_values_by_timepoint.csv",
      "Output_EasyM_MRD_analysis_2025/tables/EasyM_threshold_values_by_timepoint.csv",
      "Final Tables and Figures/Baseline dates for samples.csv",
      "nested_bm_validation_updated5.rds",
      "nested_blood_validation_updated5.rds",
      "nested_fragmentomics_validation_updated3_original.rds",
      "nested_fragmentomics_bm_validation_updated5.rds",
      "nested_fragmentomics_blood_validation_updated5.rds"
    ),
    kind = "file",
    required = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    category = c(
      "clinical_processed_inputs",
      "cohort_assignment",
      "deidentified_id_map",
      "master_feature_table",
      "cfwgs_scored_calls",
      "cfwgs_scored_calls",
      "easym_raw_inputs",
      "easym_raw_inputs",
      "easym_inputs",
      "easym_threshold_reference",
      "easym_threshold_reference",
      "clinical_date_reference",
      "cache_sensitive_model_artifact",
      "cache_sensitive_model_artifact",
      "cache_sensitive_model_artifact",
      "cache_sensitive_model_artifact",
      "cache_sensitive_model_artifact"
    ),
    reason = c(
      "Master clinical metadata table used by multiple Stage 1-3 scripts.",
      "Cohort assignment object used across figure/table scripts.",
      "De-identified patient ID map used for exported manuscript tables.",
      "Latest master feature table read by multiple downstream scripts.",
      "Primary scored cfWGS calls/probabilities used by concordance and survival scripts.",
      "Full scored cfWGS calls/probabilities used by swimplot/survival support sections.",
      "Raw quantitative EasyM table used by 3_1_A to rebuild processed EasyM calls.",
      "Raw binary EasyM table used by 3_1_A to rebuild processed EasyM calls.",
      "Processed EasyM calls used by concordance and survival comparisons.",
      "Optional EasyM threshold reference read by some sections.",
      "Optional EasyM threshold reference read by 3_2.",
      "Baseline date reference used when building Supplementary Table 8.",
      "Preserved BM nested-CV/model object used for reproducibility; needed for default scoring/reference runs if not recomputing.",
      "Preserved blood nested-CV/model object used for reproducibility; needed for default scoring/reference runs if not recomputing.",
      "Preserved fragmentomics nested-CV/model object used for reproducibility; needed for default scoring/reference runs if not recomputing.",
      "Preserved BM fragmentomics model object used for reproducibility; needed for default scoring/reference runs if not recomputing.",
      "Preserved blood fragmentomics model object used for reproducibility; needed for default scoring/reference runs if not recomputing."
    ),
    stringsAsFactors = FALSE
  )
}

copy_missing_path <- function(source_root, project_root, relative_path, kind) {
  source_path <- file.path(source_root, relative_path)
  destination_path <- file.path(project_root, relative_path)

  if (!file.exists(source_path)) {
    return(list(copied = FALSE, copy_status = "missing_from_source"))
  }

  dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)

  if (kind == "directory") {
    # Copy all top-level entries recursively without overwriting existing files.
    # This is a directory mirror operation, not a selective file-manifest copy.
    dir.create(destination_path, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(
      from = list.files(source_path, all.files = TRUE, no.. = TRUE, full.names = TRUE),
      to = destination_path,
      recursive = TRUE,
      copy.date = TRUE,
      overwrite = FALSE
    )
    copied <- length(ok) > 0 && all(ok)
  } else {
    copied <- file.copy(source_path, destination_path, overwrite = FALSE, copy.date = TRUE)
  }

  list(
    copied = isTRUE(copied),
    copy_status = if (isTRUE(copied)) "copied" else "copy_failed"
  )
}

check_one <- function(row, project_root, source_root = NULL, copy_missing = FALSE) {
  relative_path <- row[["relative_path"]]
  kind <- row[["kind"]]
  destination_path <- file.path(project_root, relative_path)
  exists_before <- file.exists(destination_path)
  copied <- FALSE
  copy_status <- if (exists_before) "already_present" else "not_requested"

  if (!exists_before && isTRUE(copy_missing)) {
    if (is.null(source_root)) {
      copy_status <- "source_dir_not_supplied"
    } else {
      copy_result <- copy_missing_path(source_root, project_root, relative_path, kind)
      copied <- copy_result$copied
      copy_status <- copy_result$copy_status
    }
  }

  exists_after <- file.exists(destination_path)
  status <- if (exists_after) {
    "present"
  } else if (isTRUE(row[["required"]])) {
    "missing_required"
  } else {
    "missing_optional"
  }

  data.frame(
    relative_path = relative_path,
    kind = kind,
    required = row[["required"]],
    category = row[["category"]],
    status = status,
    exists_before = exists_before,
    copied = copied,
    copy_status = copy_status,
    reason = row[["reason"]],
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args()
  script_dir <- get_script_dir()
  legacy_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  project_root <- if (
    dir.exists(file.path(legacy_root, "Clinical data")) ||
      dir.exists(file.path(legacy_root, "M4_CMRG_Data"))
  ) legacy_root else script_dir

  source_root <- NULL
  if (!is.null(args$source_dir)) {
    source_root <- normalizePath(args$source_dir, mustWork = TRUE)
  }

  message_log("Project root: ", project_root)
  message_log("Script directory: ", script_dir)
  if (!is.null(source_root)) message_log("Source mirror: ", source_root)
  if (isTRUE(args$copy_missing)) message_log("Copy missing inputs: enabled")

  manifest <- rbind(expected_input_roots(), critical_input_files())
  checked <- do.call(
    rbind,
    lapply(seq_len(nrow(manifest)), function(i) {
      check_one(manifest[i, , drop = FALSE], project_root, source_root, args$copy_missing)
    })
  )

  report_path <- NULL
  if (!is.null(args$report)) {
    report_path <- file.path(project_root, args$report)
    dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      checked,
      report_path,
      sep = "\t",
      row.names = FALSE,
      quote = TRUE,
      na = ""
    )
  }

  n_required_missing <- sum(checked$status == "missing_required")
  n_optional_missing <- sum(checked$status == "missing_optional")
  n_copied <- sum(checked$copied)

  if (!is.null(report_path)) {
    message_log("Local input report: ", report_path)
  } else {
    message_log("Read-only check; no report file written.")
  }
  message_log("Present: ", sum(checked$status == "present"))
  message_log("Copied: ", n_copied)
  message_log("Missing required: ", n_required_missing)
  message_log("Missing optional: ", n_optional_missing)

  if (n_required_missing > 0) {
    missing_required <- checked[checked$status == "missing_required", c("relative_path", "reason")]
    message_log("Required local inputs are missing:")
    for (i in seq_len(nrow(missing_required))) {
      message("  - ", missing_required$relative_path[i], " | ", missing_required$reason[i])
    }
    stop(
      "Missing required local inputs. Stage them under the project root or rerun with --source-dir and --copy-missing.",
      call. = FALSE
    )
  }

  invisible(checked)
}

main()
