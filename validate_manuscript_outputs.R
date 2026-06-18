#!/usr/bin/env Rscript

# =============================================================================
# validate_manuscript_outputs.R
#
# Purpose:
#   Validate the manuscript-facing output tree written directly by the numbered
#   scripts in Scripts_2025/Final_Scripts.
#
# What this checks:
#   - docs/manuscript_artifact_source_map.tsv is readable.
#   - final_manuscript_objects/ exists.
#   - manuscript_direct_output_manifest.tsv exists and covers every mapped
#     artifact ID.
#   - manuscript_output_index.tsv and script_output_index.tsv exist and are
#     complete.
#   - final_manuscript_objects/README.md exists and summarizes manuscript-facing
#     output folders, regeneration status, validation status, and caveats.
#   - manuscript_output_index.tsv contains regeneration/validation status
#     columns, so known cached/manual/value-difference items are explicit.
#   - Explicit per-panel paths and generated mirror paths exist.
#   - Easy-use final folders contain assembled main figures, extended figures,
#     main tables, and supplementary tables.
#   - Each mapped generating script contains a direct manuscript-output helper
#     call, including the supporting external-ichorCNA script.
#   - Each mapped generating script references every artifact ID that the
#     manuscript artifact source map assigns to that script.
#   - Numbered source scripts have manuscript-output headers, and the README
#     lists every script that owns a mapped final artifact.
#
# This script does not generate scientific figures/tables. It refreshes the
# non-scientific output indexes/README from the current manifest and artifact
# source map, then writes a validation report for developer/release checks.
#
# Usage from project root:
#   Rscript Scripts_2025/Final_Scripts/validate_manuscript_outputs.R
#
# Outputs:
#   Scripts_2025/Final_Scripts/final_manuscript_objects/
#     manuscript_output_validation_report.tsv
# =============================================================================

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

count_files <- function(path) {
  if (!dir.exists(path)) return(0L)
  entries <- list.files(path, recursive = FALSE, full.names = TRUE, no.. = TRUE)
  sum(file.exists(entries) & !dir.exists(entries))
}

split_paths <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}

all_listed_paths_exist <- function(values) {
  paths <- unlist(lapply(values, split_paths), use.names = FALSE)
  paths <- paths[nzchar(paths)]
  if (!length(paths)) return(FALSE)
  all(file.exists(paths))
}

missing_listed_paths <- function(values) {
  paths <- unlist(lapply(values, split_paths), use.names = FALSE)
  paths <- paths[nzchar(paths)]
  paths[!file.exists(paths)]
}

strip_external_script_prefix <- function(script_name) {
  sub(
    "^External ichorCNA workflow outside this repo; supporting repo script: ",
    "",
    script_name
  )
}

script_references_artifact_ids <- function(script_path, artifact_ids) {
  if (!file.exists(script_path)) {
    return(artifact_ids)
  }
  script_text <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
  artifact_ids[!vapply(artifact_ids, grepl, logical(1), x = script_text, fixed = TRUE)]
}

record_check <- function(checks, check_name, passed, details = "") {
  rbind(
    checks,
    data.frame(
      check = check_name,
      status = if (isTRUE(passed)) "PASS" else "FAIL",
      details = details,
      stringsAsFactors = FALSE
    )
  )
}

main <- function() {
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  helper_path <- file.path(script_dir, "manuscript_output_helpers.R")
  if (!file.exists(helper_path)) {
    stop("Missing helper: ", helper_path, call. = FALSE)
  }
  source(helper_path)

  source_map_path <- file.path(project_root, "docs", "manuscript_artifact_source_map.tsv")
  source_pipeline_path <- file.path(project_root, "reproducible_workflow", "config", "source_pipeline.tsv")
  readme_path <- file.path(script_dir, "README.md")
  output_root <- ms_output_root(project_root)
  manifest_path <- ms_manifest_path(project_root)
  output_index_path <- ms_output_index_path(project_root)
  script_index_path <- ms_script_index_path(project_root)
  output_readme_path <- ms_output_readme_path(project_root)
  validation_path <- file.path(output_root, "manuscript_output_validation_report.tsv")
  logs_validation_path <- file.path(ms_output_logs_dir(project_root), "validation_report.tsv")
  logs_validation_summary_path <- file.path(ms_output_logs_dir(project_root), "validation_summary.md")

  checks <- data.frame(check = character(), status = character(), details = character())

  source_map_exists <- file.exists(source_map_path)
  checks <- record_check(checks, "source_map_exists", source_map_exists, source_map_path)
  if (!source_map_exists) stop("Cannot validate without source map: ", source_map_path, call. = FALSE)

  source_pipeline_exists <- file.exists(source_pipeline_path)
  checks <- record_check(checks, "source_pipeline_exists", source_pipeline_exists, source_pipeline_path)
  if (!source_pipeline_exists) {
    stop("Cannot validate source-script headers without source pipeline: ", source_pipeline_path, call. = FALSE)
  }

  readme_exists <- file.exists(readme_path)
  checks <- record_check(checks, "final_scripts_readme_exists", readme_exists, readme_path)

  source_map <- utils::read.delim(source_map_path, stringsAsFactors = FALSE, check.names = FALSE)
  source_pipeline <- utils::read.delim(source_pipeline_path, stringsAsFactors = FALSE, check.names = FALSE)
  mapped_artifacts <- unique(source_map$artifact_id)
  mapped_scripts <- unique(strip_external_script_prefix(source_map$generating_script))

  ms_write_output_index(project_root)

  checks <- record_check(
    checks,
    "output_root_exists",
    dir.exists(output_root),
    output_root
  )

  checks <- record_check(
    checks,
    "manifest_exists",
    file.exists(manifest_path),
    manifest_path
  )

  checks <- record_check(
    checks,
    "output_index_exists",
    file.exists(output_index_path),
    output_index_path
  )

  checks <- record_check(
    checks,
    "script_index_exists",
    file.exists(script_index_path),
    script_index_path
  )

  checks <- record_check(
    checks,
    "output_readme_exists",
    file.exists(output_readme_path),
    output_readme_path
  )

  manifest <- if (file.exists(manifest_path)) {
    utils::read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }

  output_index <- if (file.exists(output_index_path)) {
    utils::read.delim(output_index_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }

  script_index <- if (file.exists(script_index_path)) {
    utils::read.delim(script_index_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }

  manifest_artifacts <- if ("artifact_id" %in% names(manifest)) unique(manifest$artifact_id) else character()
  index_artifacts <- if ("artifact_id" %in% names(output_index)) unique(output_index$artifact_id) else character()
  index_scripts <- if ("generating_script" %in% names(script_index)) unique(script_index$generating_script) else character()

  missing_manifest <- setdiff(mapped_artifacts, manifest_artifacts)
  checks <- record_check(
    checks,
    "manifest_covers_all_mapped_artifacts",
    length(missing_manifest) == 0,
    if (length(missing_manifest)) paste(missing_manifest, collapse = "; ") else paste("n =", length(mapped_artifacts))
  )

  missing_index <- setdiff(mapped_artifacts, index_artifacts)
  checks <- record_check(
    checks,
    "output_index_covers_all_mapped_artifacts",
    length(missing_index) == 0,
    if (length(missing_index)) paste(missing_index, collapse = "; ") else paste("n =", length(mapped_artifacts))
  )

  missing_script_index <- setdiff(mapped_scripts[file.exists(mapped_scripts)], index_scripts)
  checks <- record_check(
    checks,
    "script_index_covers_all_mapped_scripts",
    length(missing_script_index) == 0,
    if (length(missing_script_index)) paste(missing_script_index, collapse = "; ") else paste("n =", length(unique(index_scripts)))
  )

  numbered_script_paths <- file.path(script_dir, source_pipeline$script)
  missing_manuscript_headers <- numbered_script_paths[
    file.exists(numbered_script_paths) &
      !vapply(numbered_script_paths, function(path) {
        any(grepl("Manuscript outputs created/updated", readLines(path, warn = FALSE), fixed = TRUE))
      }, logical(1))
  ]
  checks <- record_check(
    checks,
    "numbered_scripts_have_manuscript_output_headers",
    length(missing_manuscript_headers) == 0,
    if (length(missing_manuscript_headers)) {
      paste(basename(missing_manuscript_headers), collapse = "; ")
    } else {
      paste("n =", sum(file.exists(numbered_script_paths)))
    }
  )

  readme_text <- if (readme_exists) paste(readLines(readme_path, warn = FALSE), collapse = "\n") else ""
  missing_readme_mapped_scripts <- basename(mapped_scripts[file.exists(mapped_scripts)])
  missing_readme_mapped_scripts <- missing_readme_mapped_scripts[
    !vapply(missing_readme_mapped_scripts, grepl, logical(1), x = readme_text, fixed = TRUE)
  ]
  checks <- record_check(
    checks,
    "readme_lists_all_mapped_output_scripts",
    length(missing_readme_mapped_scripts) == 0,
    if (length(missing_readme_mapped_scripts)) {
      paste(missing_readme_mapped_scripts, collapse = "; ")
    } else {
      paste("n =", length(unique(basename(mapped_scripts[file.exists(mapped_scripts)]))))
    }
  )

  script_artifact_groups <- split(source_map$artifact_id, strip_external_script_prefix(source_map$generating_script))
  script_reference_failures <- unlist(lapply(names(script_artifact_groups), function(script_path) {
    if (!file.exists(script_path)) return(character())
    missing_ids <- script_references_artifact_ids(script_path, unique(script_artifact_groups[[script_path]]))
    if (!length(missing_ids)) return(character())
    paste0(script_path, ": ", paste(missing_ids, collapse = ", "))
  }), use.names = FALSE)
  checks <- record_check(
    checks,
    "mapped_scripts_reference_all_assigned_artifact_ids",
    length(script_reference_failures) == 0,
    if (length(script_reference_failures)) {
      paste(script_reference_failures, collapse = "; ")
    } else {
      paste("scripts =", length(names(script_artifact_groups)))
    }
  )

  if (nrow(output_index)) {
    required_status_columns <- c(
      "regeneration_status",
      "validation_status",
      "validation_scope",
      "source_confidence"
    )
    missing_status_columns <- setdiff(required_status_columns, names(output_index))
    checks <- record_check(
      checks,
      "output_index_has_status_columns",
      length(missing_status_columns) == 0,
      if (length(missing_status_columns)) {
        paste(missing_status_columns, collapse = "; ")
      } else {
        paste(required_status_columns, collapse = "; ")
      }
    )

    if (length(missing_status_columns) == 0) {
      missing_status_rows <- output_index$artifact_id[
        !nzchar(output_index$regeneration_status) |
          !nzchar(output_index$validation_status) |
          !nzchar(output_index$validation_scope) |
          !nzchar(output_index$source_confidence)
      ]
      checks <- record_check(
        checks,
        "output_index_status_values_complete",
        length(missing_status_rows) == 0,
        if (length(missing_status_rows)) {
          paste(missing_status_rows, collapse = "; ")
        } else {
          paste("n =", nrow(output_index))
        }
      )
    }

    non_organized <- output_index$artifact_id[output_index$organization_status != "organized"]
    checks <- record_check(
      checks,
      "all_output_index_rows_organized",
      length(non_organized) == 0,
      if (length(non_organized)) paste(non_organized, collapse = "; ") else paste("n =", nrow(output_index))
    )

    checks <- record_check(
      checks,
      "explicit_traceability_paths_exist",
      all_listed_paths_exist(output_index$explicit_traceability_paths),
      paste("rows =", nrow(output_index))
    )

    checks <- record_check(
      checks,
      "generated_mirror_paths_exist",
      all_listed_paths_exist(output_index$generated_mirror_paths),
      paste("rows =", nrow(output_index))
    )
  }

  reference_style_dirs <- file.path(
    output_root,
    c(
      "generated/main_figures",
      "generated/extended_figures",
      "generated/figure_components",
      "generated/main_tables",
      "generated/main_tables/source_data",
      "generated/manual_assembly_reference",
      "generated/native_script_runs",
      "generated/supplementary_tables",
      "logs/rendered_final_figures",
      "frozen/main_figures",
      "frozen/extended_figures",
      "frozen/main_tables",
      "frozen/supplementary_tables"
    )
  )
  missing_reference_style_dirs <- reference_style_dirs[!dir.exists(reference_style_dirs)]
  checks <- record_check(
    checks,
    "reference_style_output_dirs_exist",
    length(missing_reference_style_dirs) == 0,
    if (length(missing_reference_style_dirs)) {
      paste(missing_reference_style_dirs, collapse = "; ")
    } else {
      paste("n =", length(reference_style_dirs))
    }
  )

  reference_style_log_files <- file.path(
    output_root,
    "logs",
    c(
      "output_index.tsv",
      "script_output_index.tsv",
      "generation_manifest.tsv",
      "manual_assembly_reference_manifest.tsv",
      "native_script_runs_manifest.tsv",
      "source_availability_manifest.tsv",
      "artifact_dependency_manifest.tsv",
      "artifact_regeneration_status.tsv",
      "cache_rebuild_matrix.tsv",
      "cached_rds_rebuild_manifest.tsv",
      "raw_to_final_pipeline_status.tsv",
      "source_pipeline_plan.tsv",
      "test_cohort_update_plan.tsv",
      "rendered_final_figures_manifest.tsv"
    )
  )
  missing_reference_style_log_files <- reference_style_log_files[!file.exists(reference_style_log_files)]
  checks <- record_check(
    checks,
    "reference_style_log_files_exist",
    length(missing_reference_style_log_files) == 0,
    if (length(missing_reference_style_log_files)) {
      paste(missing_reference_style_log_files, collapse = "; ")
    } else {
      paste("n =", length(reference_style_log_files))
    }
  )

  if (nrow(manifest)) {
    missing_destinations <- manifest$destination_path[!file.exists(manifest$destination_path)]
    checks <- record_check(
      checks,
      "manifest_destinations_exist",
      length(missing_destinations) == 0,
      if (length(missing_destinations)) paste(missing_destinations, collapse = "; ") else paste("rows =", nrow(manifest))
    )

    if ("generated_mirror_path" %in% names(manifest)) {
      missing_generated <- missing_listed_paths(manifest$generated_mirror_path)
      checks <- record_check(
        checks,
        "manifest_generated_mirror_paths_exist",
        length(missing_generated) == 0,
        if (length(missing_generated)) paste(missing_generated, collapse = "; ") else paste("rows =", nrow(manifest))
      )
    }
  }

  rendered_manifest_path <- file.path(output_root, "logs", "rendered_final_figures_manifest.tsv")
  rendered_manifest <- if (file.exists(rendered_manifest_path)) {
    utils::read.delim(rendered_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }

  manual_reference_manifest_path <- file.path(output_root, "logs", "manual_assembly_reference_manifest.tsv")
  manual_reference_manifest <- if (file.exists(manual_reference_manifest_path)) {
    utils::read.delim(manual_reference_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }

  native_script_runs_manifest_path <- file.path(output_root, "logs", "native_script_runs_manifest.tsv")
  native_script_runs_manifest <- if (file.exists(native_script_runs_manifest_path)) {
    utils::read.delim(native_script_runs_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }

  expected_final_counts <- aggregate(
    final_filename ~ artifact_category,
    data = unique(source_map[, c("artifact_category", "final_filename")]),
    FUN = function(x) length(unique(x[!is.na(x) & nzchar(x)]))
  )
  expected_counts <- setNames(expected_final_counts$final_filename, expected_final_counts$artifact_category)

  final_folder_checks <- list(
    main_figure = file.path(output_root, "main_figures"),
    extended_data_figure = file.path(output_root, "extended_figures"),
    main_table = file.path(output_root, "main_tables"),
    supplementary_table = file.path(output_root, "supplementary_tables")
  )

  for (category in names(final_folder_checks)) {
    expected <- unname(expected_counts[[category]])
    if (is.null(expected)) expected <- 0L
    observed <- count_files(final_folder_checks[[category]])
    checks <- record_check(
      checks,
      paste0("final_folder_count_", category),
      observed >= expected,
      paste("observed =", observed, "expected >=", expected, "folder =", final_folder_checks[[category]])
    )
  }

  generated_final_folder_checks <- list(
    main_figure = file.path(output_root, "generated", "main_figures"),
    extended_data_figure = file.path(output_root, "generated", "extended_figures"),
    main_table = file.path(output_root, "generated", "main_tables"),
    supplementary_table = file.path(output_root, "generated", "supplementary_tables")
  )

  for (category in names(generated_final_folder_checks)) {
    expected <- unname(expected_counts[[category]])
    if (is.null(expected)) expected <- 0L
    observed <- count_files(generated_final_folder_checks[[category]])
    checks <- record_check(
      checks,
      paste0("generated_final_folder_count_", category),
      observed >= expected,
      paste("observed =", observed, "expected >=", expected, "folder =", generated_final_folder_checks[[category]])
    )
  }

  table1_source_data_dir <- file.path(output_root, "generated", "main_tables", "source_data")
  table1_source_files <- list.files(table1_source_data_dir, full.names = FALSE)
  expected_table1_source_patterns <- c(
    "Table_1_clinical_demographics_computed_source_counts[.]csv$",
    "Table_1_clinical_demographics_computed_source_qc[.]tsv$"
  )
  table1_source_hits <- vapply(
    expected_table1_source_patterns,
    function(pattern) any(grepl(pattern, table1_source_files)),
    logical(1)
  )
  checks <- record_check(
    checks,
    "table1_main_table_source_data_available",
    all(table1_source_hits),
    if (all(table1_source_hits)) {
      paste("n =", sum(table1_source_hits))
    } else {
      paste(expected_table1_source_patterns[!table1_source_hits], collapse = "; ")
    }
  )

  expected_figure_pdf_count <- sum(expected_counts[c("main_figure", "extended_data_figure")], na.rm = TRUE)
  if (nrow(rendered_manifest) && all(c("status", "renderer", "rendered_png") %in% names(rendered_manifest))) {
    renderers <- unique(rendered_manifest$renderer)
    rendered_pngs <- rendered_manifest$rendered_png[
      rendered_manifest$status == "rendered" & file.exists(rendered_manifest$rendered_png)
    ]
    render_ok <- if (all(rendered_manifest$status == "skipped") && all(renderers == "unavailable")) {
      file.exists(file.path(output_root, "logs", "rendered_final_figures", "README.md"))
    } else {
      length(rendered_pngs) >= expected_figure_pdf_count &&
        !any(rendered_manifest$status %in% c("failed", "pending"))
    }
    checks <- record_check(
      checks,
      "rendered_final_figure_previews_available_or_documented",
      render_ok,
      paste(
        "rendered =",
        length(rendered_pngs),
        "expected >=",
        expected_figure_pdf_count,
        "renderer =",
        paste(renderers, collapse = ",")
      )
    )
  } else {
    checks <- record_check(
      checks,
      "rendered_final_figure_previews_available_or_documented",
      FALSE,
      rendered_manifest_path
    )
  }

  manual_count <- if (nrow(output_index) && "manual_or_external_final_assembly" %in% names(output_index)) {
    sum(!is.na(output_index$manual_or_external_final_assembly) & output_index$manual_or_external_final_assembly)
  } else {
    0L
  }
  copied_manual_references <- if (
    nrow(manual_reference_manifest) &&
      all(c("artifact_id", "reference_path", "status") %in% names(manual_reference_manifest))
  ) {
    unique(manual_reference_manifest$artifact_id[
      manual_reference_manifest$status == "copied" &
        file.exists(manual_reference_manifest$reference_path)
    ])
  } else {
    character()
  }
  checks <- record_check(
    checks,
    "manual_assembly_references_available",
    if (manual_count == 0L) TRUE else length(copied_manual_references) >= manual_count,
    paste("manual artifacts =", manual_count, "with references =", length(copied_manual_references))
  )

  native_fig1a_artifacts <- if (
    nrow(native_script_runs_manifest) &&
      all(c("snapshot", "artifact_id", "snapshot_path", "status") %in% names(native_script_runs_manifest))
  ) {
    unique(native_script_runs_manifest$artifact_id[
      native_script_runs_manifest$snapshot == "FIG1A_swim_plot" &
        native_script_runs_manifest$status == "copied" &
        file.exists(native_script_runs_manifest$snapshot_path)
    ])
  } else {
    character()
  }
  checks <- record_check(
    checks,
    "native_script_run_snapshot_fig1a_available",
    all(c("FIG1A", "STABLE1") %in% native_fig1a_artifacts),
    if (length(native_fig1a_artifacts)) {
      paste(native_fig1a_artifacts, collapse = "; ")
    } else {
      native_script_runs_manifest_path
    }
  )

  script_helper_presence <- vapply(mapped_scripts[file.exists(mapped_scripts)], function(path) {
    any(grepl("ms_copy_artifact|ms_save_plot", readLines(path, warn = FALSE)))
  }, logical(1))

  missing_helper <- names(script_helper_presence)[!script_helper_presence]
  checks <- record_check(
    checks,
    "mapped_scripts_call_manuscript_helpers",
    length(missing_helper) == 0,
    if (length(missing_helper)) paste(missing_helper, collapse = "; ") else paste("n =", length(script_helper_presence))
  )

  utils::write.table(
    checks,
    file = validation_path,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )
  file.copy(validation_path, logs_validation_path, overwrite = TRUE)

  failed <- checks[checks$status != "PASS", , drop = FALSE]
  summary_lines <- c(
    "# Final_Scripts Manuscript Output Validation",
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- Checks: ", nrow(checks)),
    paste0("- Passed: ", sum(checks$status == "PASS")),
    paste0("- Failed: ", nrow(failed)),
    paste0("- Report: `", basename(logs_validation_path), "`"),
    "",
    "## Scope",
    "",
    "- Validates the direct `Final_Scripts/final_manuscript_objects` output tree.",
    "- Checks manuscript artifact coverage, script/header ownership, generated mirrors, final easy-use folders, rendered final-figure previews, and reference-style log files.",
    "- Does not rerun scientific analyses."
  )
  if (nrow(failed)) {
    summary_lines <- c(
      summary_lines,
      "",
      "## Failed Checks",
      "",
      paste0("- ", failed$check, ": ", failed$details)
    )
  }
  writeLines(summary_lines, logs_validation_summary_path)

  cat("Validation report:", validation_path, "\n")
  cat("Checks:", nrow(checks), "\n")
  cat("Passed:", sum(checks$status == "PASS"), "\n")
  cat("Failed:", nrow(failed), "\n")

  if (nrow(failed)) {
    cat("\nFailed checks:\n")
    utils::write.table(failed, sep = "\t", row.names = FALSE, quote = FALSE)
    quit(save = "no", status = 1)
  }

  cat("MANUSCRIPT_OUTPUT_VALIDATION_OK\n")
}

main()
