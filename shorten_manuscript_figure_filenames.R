# =============================================================================
# shorten_manuscript_figure_filenames.R
#
# Purpose:
#   Optional internal cleanup that renames existing figure image/PDF files in
#   final_manuscript_objects/ to shorter filenames. This script is not used to
#   generate a manuscript figure, table, statistic, or source-data workbook.
#
# IMPORTANT:
#   This script changes filenames in place and rewrites matching absolute paths
#   in every TSV under final_manuscript_objects/. It does not change image pixels
#   or scientific values, but it can invalidate links held outside that folder.
#   There is no dry-run or confirmation step: sourcing or running the file starts
#   the rename immediately. Preserve a backup before intentional use.
#
# Known recovery risk:
#   The rename manifest is written inside final_manuscript_objects/logs/ before
#   all TSVs under that tree are rewritten. Because that search includes the
#   manifest itself, its old_path values can be replaced with new_path values.
#   The current script therefore does not reliably preserve a reversal map.
#
# Scope:
#   - Organized manuscript figure component folders.
#   - Additional exploratory and sensitivity figure folders used in the manuscript.
#   - Top-level stray figure image exports in final_manuscript_objects/.
#
# This file-maintenance step is not required to reproduce the manuscript.
# =============================================================================

source("Scripts_2025/Final_Scripts/manuscript_output_helpers.R")

project_root <- ms_find_project_root()
output_root <- ms_output_root(project_root)

figure_ext_pattern <- "[.](png|pdf|jpg|jpeg|tif|tiff|svg)$"

target_roots <- file.path(
  output_root,
  c(
    "01_main_figures",
    "02_extended_data_figures",
    "additional_all_evaluable_first_nonbaseline_km",
    "additional_all_evaluable_longitudinal_panels",
    "generated/figure_components",
    "generated/manual_assembly_reference",
    "generated/native_script_runs"
  )
)
target_roots <- target_roots[dir.exists(target_roots)]

target_files <- unlist(
  lapply(
    target_roots,
    function(root) {
      files <- list.files(root, pattern = figure_ext_pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
      files[!grepl("/source_data/", files, fixed = FALSE)]
    }
  ),
  use.names = FALSE
)

top_level_files <- list.files(output_root, pattern = figure_ext_pattern, full.names = TRUE, ignore.case = TRUE)
target_files <- sort(unique(c(target_files, top_level_files)))

compact_label_from_text <- function(text) {
  text <- basename(text)
  hit <- regmatches(text, regexpr("Extended_Data_Figure_[0-9]+[A-Za-z]?", text))
  if (length(hit) && nzchar(hit)) return(ms_compact_artifact_label(hit))

  hit <- regmatches(text, regexpr("Figure_[0-9]+[A-Za-z]?", text))
  if (length(hit) && nzchar(hit)) return(ms_compact_artifact_label(hit))

  hit <- regmatches(text, regexpr("EDFIG[0-9]+[A-Za-z]?", text))
  if (length(hit) && nzchar(hit)) return(sub("^EDFIG", "ED", hit))

  hit <- regmatches(text, regexpr("FIG[0-9]+[A-Za-z]?", text))
  if (length(hit) && nzchar(hit)) return(sub("^FIG", "F", hit))

  ""
}

compact_label_from_path <- function(path) {
  parts <- strsplit(normalizePath(path, mustWork = FALSE), .Platform$file.sep, fixed = TRUE)[[1]]
  parts <- rev(parts[nzchar(parts)])
  for (part in parts) {
    label <- compact_label_from_text(part)
    if (nzchar(label)) return(label)
  }

  lower <- tolower(tools::file_path_sans_ext(basename(path)))
  if (grepl("^km_blood_cfwgs_combined", lower)) return("KM_blood_combined")
  if (grepl("^km_blood_cfwgs_sites", lower)) return("KM_blood_sites")
  if (grepl("^km_bm_cfwgs_cvaf", lower)) return("KM_BM_cvaf")
  if (grepl("^km_clonoseq", lower)) return("KM_clonoSEQ")
  if (grepl("^km_mfc", lower)) return("KM_MFC")
  if (grepl("^km_", lower)) return("KM")

  ""
}

compact_role_from_filename <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  lower <- tolower(stem)

  patterns <- c(
    "first_nonbaseline_event_free_30d" = "first_nonbase_EFS",
    "first_maintenance_event_free_30d" = "first_maint_EFS",
    "first_nonbaseline" = "first_nonbase_PFS",
    "first_maintenance" = "first_maint_PFS",
    "sustained_first_two_maintenance" = "sustained_maint2",
    "cross_platform_patient_lines_with_novaseq6000_zero_controls" = "alt_novaseq_zero",
    "patient_lines_zero_final_tf" = "alt_zero_TF",
    "dilution_patient_lines_zero_final_tf" = "dilution_zero_TF",
    "cvaf_zscore_model" = "alt_cvaf_z",
    "combined_model" = "alt_combined",
    "cvaf_model" = "alt_cvaf",
    "all_evaluable_capped300" = "all_eval_cap300",
    "all_evaluable" = "all_eval",
    "all_samples_mfc_only" = "all_samples_MFC",
    "all_samples" = "all_samples",
    "manual_final" = "manual_final",
    "external_final" = "external_final"
  )
  for (pattern in names(patterns)) {
    if (grepl(pattern, lower, fixed = TRUE)) return(unname(patterns[[pattern]]))
  }

  hint <- ms_compact_source_hint(path)
  if (nzchar(hint)) return(hint)

  if (grepl("figure_component", lower, fixed = TRUE)) return("component")
  if (grepl("figure_panel", lower, fixed = TRUE)) return("panel")
  "figure"
}

compact_extra_from_filename <- function(path) {
  lower <- tolower(tools::file_path_sans_ext(basename(path)))

  if (grepl("ca_08|ca-08", lower)) return("CA08")
  if (grepl("ca_02|ca-02", lower)) return("CA02")
  if (grepl("training", lower)) return("training")
  if (grepl("test", lower)) return("test")
  if (grepl("fragmentomics|fragment", lower)) return("frag")
  if (grepl("blood", lower)) return("blood")
  if (grepl("bm_", lower) || grepl("_bm", lower) || grepl("bm-derived", lower)) return("BM")
  if (grepl("easym", lower)) return("EasyM")
  if (grepl("mfc", lower)) return("MFC")
  if (grepl("clonoseq", lower)) return("clonoSEQ")
  if (grepl("sites", lower)) return("sites")
  if (grepl("combined", lower)) return("combined")
  if (grepl("cvaf_z", lower) || grepl("zscore", lower)) return("cvaf_z")
  if (grepl("cvaf", lower)) return("cvaf")
  ""
}

compact_name_for_path <- function(path) {
  current_name <- basename(path)
  if (
    nchar(current_name) <= 55 &&
      !grepl("__|[[:space:]()]|[+]", current_name)
  ) {
    return(current_name)
  }

  label <- compact_label_from_path(path)
  if (!nzchar(label)) {
    hint <- ms_compact_source_hint(path)
    if (nzchar(hint)) {
      label <- hint
    } else {
      label <- ms_slug(tools::file_path_sans_ext(basename(path)))
      if (nchar(label) > 28) label <- substr(label, 1, 28)
    }
  }

  role <- compact_role_from_filename(path)
  extra <- compact_extra_from_filename(path)
  ext <- tools::file_ext(path)
  ext <- if (nzchar(ext)) paste0(".", ext) else ""

  tokens <- unique(c(label, role, extra))
  tokens <- tokens[nzchar(tokens)]
  filename <- paste0(paste(tokens, collapse = "_"), ext)

  if (nchar(filename) > 55) {
    filename <- paste0(paste(tokens[seq_len(min(2, length(tokens)))], collapse = "_"), ext)
  }
  filename
}

planned <- data.frame(
  old_path = normalizePath(target_files, mustWork = TRUE),
  new_path = file.path(dirname(target_files), vapply(target_files, compact_name_for_path, character(1))),
  stringsAsFactors = FALSE
)
planned$new_path <- normalizePath(planned$new_path, mustWork = FALSE)
planned <- planned[planned$old_path != planned$new_path, , drop = FALSE]

if (nrow(planned)) {
  for (dir_name in unique(dirname(planned$new_path))) {
    in_dir <- planned[dirname(planned$new_path) == dir_name, , drop = FALSE]
    duplicated_targets <- duplicated(in_dir$new_path) | duplicated(in_dir$new_path, fromLast = TRUE)
    if (any(duplicated_targets)) {
      dup_targets <- unique(in_dir$new_path[duplicated_targets])
      for (target in dup_targets) {
        idx <- which(planned$new_path == target)
        for (j in seq_along(idx)) {
          stem <- tools::file_path_sans_ext(planned$new_path[idx[j]])
          ext <- tools::file_ext(planned$new_path[idx[j]])
          ext <- if (nzchar(ext)) paste0(".", ext) else ""
          planned$new_path[idx[j]] <- paste0(stem, "_", j, ext)
        }
      }
    }
  }

  existing_collision <- file.exists(planned$new_path) &
    !(normalizePath(planned$new_path, mustWork = FALSE) %in% planned$old_path)
  if (any(existing_collision)) {
    for (idx in which(existing_collision)) {
      stem <- tools::file_path_sans_ext(planned$new_path[[idx]])
      ext <- tools::file_ext(planned$new_path[[idx]])
      ext <- if (nzchar(ext)) paste0(".", ext) else ""
      suffix <- 2L
      candidate <- paste0(stem, "_", suffix, ext)
      while (file.exists(candidate) || candidate %in% planned$new_path[-idx]) {
        suffix <- suffix + 1L
        candidate <- paste0(stem, "_", suffix, ext)
      }
      planned$new_path[[idx]] <- candidate
    }
  }

  rename_manifest <- file.path(output_root, "logs", "figure_filename_shortening_manifest.tsv")
  dir.create(dirname(rename_manifest), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(planned, rename_manifest, sep = "\t", row.names = FALSE, quote = TRUE)

  # Renames occur sequentially with no transaction or rollback. A failure after
  # an earlier successful rename leaves a partially changed output tree.
  for (i in seq_len(nrow(planned))) {
    ok <- file.rename(planned$old_path[[i]], planned$new_path[[i]])
    if (!ok) stop("Failed to rename: ", planned$old_path[[i]], call. = FALSE)
  }

  # This list currently includes rename_manifest itself; see the recovery-risk
  # warning in the file header.
  tsv_files <- list.files(output_root, pattern = "[.]tsv$", recursive = TRUE, full.names = TRUE)
  for (tsv in tsv_files) {
    text <- readLines(tsv, warn = FALSE)
    if (!length(text)) next
    updated <- text
    for (i in seq_len(nrow(planned))) {
      updated <- gsub(planned$old_path[[i]], planned$new_path[[i]], updated, fixed = TRUE)
    }
    if (!identical(text, updated)) writeLines(updated, tsv, useBytes = TRUE)
  }
}

remaining <- target_files[file.exists(target_files)]
remaining <- c(remaining, planned$new_path[file.exists(planned$new_path)])
remaining <- unique(remaining[file.exists(remaining)])
remaining_lengths <- nchar(basename(remaining))

summary <- data.frame(
  renamed_files = nrow(planned),
  max_filename_length_after = if (length(remaining_lengths)) max(remaining_lengths) else 0,
  files_over_60_chars_after = sum(remaining_lengths > 60),
  stringsAsFactors = FALSE
)
utils::write.table(
  summary,
  file.path(output_root, "logs", "figure_filename_shortening_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = TRUE
)

print(summary)
