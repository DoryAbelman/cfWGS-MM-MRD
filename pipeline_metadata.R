# =============================================================================
# pipeline_metadata.R
#
# Purpose:
# Shared metadata helpers for the numbered cfWGS-MM-MRD scripts.
#
# These helpers intentionally do not perform any scientific analysis. They read
# the source pipeline plan and manuscript artifact source map so command-line
# tools can report what each numbered script does and which manuscript outputs
# it affects.
#
# Manuscript outputs created/updated:
#   - None directly. This support/provenance file is used by workflow runners
#     and artifact-map builders to keep script order and manuscript-output
#     mappings auditable.
#
# Usage:
#   source("Scripts_2025/Final_Scripts/pipeline_metadata.R")
#   # Usually sourced by run_pipeline.R, run_manuscript_workflow.R, and
#   # build_stage_artifact_map.R rather than run directly.
# =============================================================================

fs_require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

fs_read_tsv <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path, call. = FALSE)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

fs_project_root_from_script_dir <- function(script_dir) {
  normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
}

fs_read_source_pipeline <- function(project_root) {
  plan_path <- file.path(project_root, "reproducible_workflow", "config", "source_pipeline.tsv")
  plan <- fs_read_tsv(plan_path, "source pipeline plan")
  fs_require_columns(
    plan,
    c("order", "script_id", "stage", "script", "run_policy", "notes"),
    "source pipeline plan"
  )
  plan
}

fs_read_artifact_source_map <- function(project_root, required = FALSE) {
  map_path <- file.path(project_root, "docs", "manuscript_artifact_source_map.tsv")
  if (!file.exists(map_path)) {
    if (isTRUE(required)) stop("Missing manuscript artifact source map: ", map_path, call. = FALSE)
    return(data.frame())
  }
  source_map <- utils::read.delim(map_path, check.names = FALSE, stringsAsFactors = FALSE)
  fs_require_columns(
    source_map,
    c(
      "artifact_id", "artifact_category", "artifact", "panel_or_sheet",
      "final_filename", "current_final_artifact_path", "generating_script",
      "source_component_or_original_output", "confidence", "notes"
    ),
    "manuscript artifact source map"
  )
  source_map
}

fs_basename_script <- function(script_path) {
  basename(trimws(script_path))
}

fs_nonempty_unique <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

fs_collapse_unique <- function(x, limit = 8) {
  values <- fs_nonempty_unique(x)
  if (!length(values)) return("")
  if (length(values) > limit) {
    values <- c(values[seq_len(limit)], paste0("... +", length(values) - limit, " more"))
  }
  paste(values, collapse = "; ")
}

fs_slug_label <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

fs_format_manuscript_item <- function(artifact, panel_or_sheet) {
  artifact <- trimws(as.character(artifact))
  panel_or_sheet <- trimws(as.character(panel_or_sheet))
  no_panel <- is.na(panel_or_sheet) | !nzchar(panel_or_sheet) |
    panel_or_sheet %in% c("all", "all_sheets")
  out <- artifact
  one_letter_panel <- !no_panel & grepl("^[A-Z]$", panel_or_sheet)
  out[one_letter_panel] <- paste0(artifact[one_letter_panel], panel_or_sheet[one_letter_panel])
  other_panel <- !no_panel & !one_letter_panel
  out[other_panel] <- paste0(artifact[other_panel], "_", fs_slug_label(panel_or_sheet[other_panel]))
  out
}

fs_optional_tsv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

fs_extract_generator_function <- function(notes) {
  notes <- as.character(notes)
  out <- rep("", length(notes))
  matched <- grepl("Generator:\\s*[A-Za-z0-9_]+\\(", notes)
  out[matched] <- sub(".*Generator:\\s*([A-Za-z0-9_]+)\\(.*", "\\1", notes[matched])
  out
}

fs_collapse_rows_by_artifact <- function(data, value_col, artifact_col = "artifact_id") {
  if (!nrow(data) || !all(c(artifact_col, value_col) %in% names(data))) {
    return(data.frame(artifact_id = character(), value = character(), stringsAsFactors = FALSE))
  }
  artifacts <- fs_nonempty_unique(data[[artifact_col]])
  rows <- lapply(artifacts, function(artifact_id) {
    data.frame(
      artifact_id = artifact_id,
      value = fs_collapse_unique(data[data[[artifact_col]] == artifact_id, value_col], limit = 10),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fs_collapse_generation_by_artifact <- function(generation) {
  if (!nrow(generation) || !"artifact_id" %in% names(generation)) return(data.frame())
  generation$clean_generator_function <- fs_extract_generator_function(generation$notes)
  artifacts <- fs_nonempty_unique(generation$artifact_id)
  rows <- lapply(artifacts, function(artifact_id) {
    x <- generation[generation$artifact_id == artifact_id, , drop = FALSE]
    data.frame(
      artifact_id = artifact_id,
      artifact_type = fs_collapse_unique(x$artifact_type, limit = 5),
      output_path = fs_collapse_unique(x$output_path, limit = 10),
      source_path = fs_collapse_unique(x$source_path, limit = 10),
      clean_generator_function = fs_collapse_unique(x$clean_generator_function, limit = 5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fs_build_stage_artifact_map <- function(project_root) {
  plan <- fs_read_source_pipeline(project_root)
  source_map <- fs_read_artifact_source_map(project_root, required = TRUE)
  status_path <- file.path(project_root, "docs", "artifact_regeneration_status.tsv")
  generation_path <- file.path(project_root, "reproducible_workflow", "outputs", "logs", "generation_manifest.tsv")
  validation_path <- file.path(project_root, "reproducible_workflow", "outputs", "logs", "validation_report.tsv")

  status <- fs_optional_tsv(status_path)
  generation <- fs_optional_tsv(generation_path)
  validation <- fs_optional_tsv(validation_path)

  source_map$script <- fs_basename_script(source_map$generating_script)
  source_map$source_lines <- ifelse(
    nzchar(as.character(source_map$line_start)) | nzchar(as.character(source_map$line_end)),
    paste0(source_map$line_start, "-", source_map$line_end),
    ""
  )

  if (nrow(status)) {
    keep <- c("artifact_id", "regeneration_status", "generated_output_path", "reference_artifact_path", "key_assumptions_or_caveats")
    status <- status[, intersect(keep, names(status)), drop = FALSE]
  }
  if (nrow(generation)) {
    generation <- fs_collapse_generation_by_artifact(generation)
  }
  if (nrow(validation)) {
    validation_col <- if ("status" %in% names(validation)) "status" else "validation_status"
    val <- fs_collapse_rows_by_artifact(validation, validation_col)
    names(val)[names(val) == "value"] <- "validation_status"
    validation <- val
  }

  artifact_rows <- merge(
    source_map,
    plan,
    by.x = "script",
    by.y = "script",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_artifact", "_script")
  )
  if ("notes_script" %in% names(artifact_rows)) {
    artifact_rows$notes <- artifact_rows$notes_script
  } else if (!"notes" %in% names(artifact_rows)) {
    artifact_rows$notes <- ""
  }
  if ("notes_artifact" %in% names(artifact_rows)) {
    artifact_rows$artifact_notes <- artifact_rows$notes_artifact
  } else {
    artifact_rows$artifact_notes <- ""
  }
  if (nrow(status)) artifact_rows <- merge(artifact_rows, status, by = "artifact_id", all.x = TRUE, sort = FALSE)
  if (nrow(generation)) artifact_rows <- merge(artifact_rows, generation, by = "artifact_id", all.x = TRUE, sort = FALSE)
  if (nrow(validation)) artifact_rows <- merge(artifact_rows, validation, by = "artifact_id", all.x = TRUE, sort = FALSE)

  artifact_rows$manuscript_item <- fs_format_manuscript_item(
    artifact_rows$artifact,
    artifact_rows$panel_or_sheet
  )
  artifact_rows$generator_or_source <- ifelse(
    "clean_generator_function" %in% names(artifact_rows) & nzchar(as.character(artifact_rows$clean_generator_function)),
    artifact_rows$clean_generator_function,
    "original_script_or_manual_component"
  )
  artifact_rows$row_type <- "direct_manuscript_artifact"

  artifact_cols <- c(
    "order", "script_id", "stage", "script", "run_policy", "notes",
    "row_type", "artifact_id", "artifact_category", "manuscript_item",
    "final_filename", "panel_or_sheet", "regeneration_status", "validation_status",
    "generator_or_source", "generated_output_path", "output_path", "source_path",
    "reference_artifact_path", "current_final_artifact_path", "source_component_or_original_output",
    "source_lines", "function_or_section", "export_call", "key_objects", "input_files",
    "intermediate_files", "confidence", "artifact_notes", "key_assumptions_or_caveats"
  )
  artifact_cols <- intersect(artifact_cols, names(artifact_rows))
  artifact_rows <- artifact_rows[, artifact_cols, drop = FALSE]

  direct_scripts <- fs_nonempty_unique(artifact_rows$script)
  upstream <- plan[!plan$script %in% direct_scripts, , drop = FALSE]
  if (nrow(upstream)) {
    upstream$row_type <- "upstream_dependency_or_intermediate"
    upstream$artifact_id <- ""
    upstream$artifact_category <- ""
    upstream$manuscript_item <- ""
    upstream$final_filename <- ""
    upstream$panel_or_sheet <- ""
    upstream$regeneration_status <- ""
    upstream$validation_status <- ""
    upstream$generator_or_source <- ""
    upstream$generated_output_path <- ""
    upstream$output_path <- ""
    upstream$source_path <- ""
    upstream$reference_artifact_path <- ""
    upstream$current_final_artifact_path <- ""
    upstream$source_component_or_original_output <- ""
    upstream$source_lines <- ""
    upstream$function_or_section <- ""
    upstream$export_call <- ""
    upstream$key_objects <- ""
    upstream$input_files <- ""
    upstream$intermediate_files <- ""
    upstream$confidence <- ""
    upstream$artifact_notes <- ""
    upstream$key_assumptions_or_caveats <- "No direct manuscript figure/table export is mapped yet; this script feeds downstream artifacts."
    upstream <- upstream[, artifact_cols, drop = FALSE]
    artifact_rows <- rbind(artifact_rows, upstream)
  }

  artifact_rows <- artifact_rows[order(as.numeric(artifact_rows$order), artifact_rows$artifact_id), , drop = FALSE]
  row.names(artifact_rows) <- NULL
  artifact_rows
}

fs_write_stage_artifact_map <- function(project_root, output_tsv = NULL, output_md = NULL) {
  stage_map <- fs_build_stage_artifact_map(project_root)
  if (is.null(output_tsv)) output_tsv <- file.path(project_root, "docs", "stage_ordered_script_artifact_map.tsv")
  if (is.null(output_md)) output_md <- file.path(project_root, "docs", "stage_ordered_script_artifact_map.md")

  utils::write.table(stage_map, output_tsv, sep = "\t", row.names = FALSE, quote = TRUE, na = "")

  con <- file(output_md, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("# Stage-Ordered Script To Manuscript Artifact Map", con)
  writeLines("", con)
  writeLines("This file is generated by `Scripts_2025/Final_Scripts/build_stage_artifact_map.R`.", con)
  writeLines("It keeps the original stage/question organization while showing the clean manuscript artifact outputs, generator functions, regeneration status, and validation status.", con)
  writeLines("", con)

  script_rows <- stage_map[!duplicated(stage_map$script), c("order", "script_id", "stage", "script", "notes", "run_policy"), drop = FALSE]
  script_rows <- script_rows[order(as.numeric(script_rows$order)), , drop = FALSE]
  previous_stage <- ""

  for (i in seq_len(nrow(script_rows))) {
    stage <- script_rows$stage[i]
    if (!identical(stage, previous_stage)) {
      writeLines(paste0("## ", sprintf("%02d", as.integer(script_rows$order[i])), " ", stage), con)
      writeLines("", con)
      previous_stage <- stage
    }

    direct <- stage_map[
      stage_map$script == script_rows$script[i] &
        stage_map$row_type == "direct_manuscript_artifact",
      ,
      drop = FALSE
    ]

    writeLines(paste0("### ", script_rows$script_id[i], " - `", script_rows$script[i], "`"), con)
    writeLines("", con)
    writeLines(paste0("- Purpose: ", script_rows$notes[i]), con)
    writeLines(paste0("- Run policy: `", script_rows$run_policy[i], "`"), con)

    if (!nrow(direct)) {
      writeLines("- Manuscript outputs: none directly; upstream/intermediate dependency.", con)
    } else {
      writeLines("- Manuscript outputs:", con)
      for (j in seq_len(nrow(direct))) {
        status <- if ("regeneration_status" %in% names(direct)) direct$regeneration_status[j] else ""
        validation <- if ("validation_status" %in% names(direct)) direct$validation_status[j] else ""
        generator <- if ("generator_or_source" %in% names(direct)) direct$generator_or_source[j] else ""
        output <- if ("generated_output_path" %in% names(direct) && nzchar(as.character(direct$generated_output_path[j]))) {
          direct$generated_output_path[j]
        } else if ("output_path" %in% names(direct)) {
          direct$output_path[j]
        } else {
          ""
        }
        output_text <- if (nzchar(as.character(output))) paste0(" | output: `", output, "`") else ""
        writeLines(
          paste0(
            "  - `", direct$artifact_id[j], "`: ", direct$manuscript_item[j],
            " | status: `", status,
            "` | validation: `", validation,
            "` | generator/source: `", generator, "`",
            output_text
          ),
          con
        )
      }
    }
    writeLines("", con)
  }

  invisible(list(tsv = output_tsv, md = output_md))
}

fs_artifact_display_name <- function(source_map) {
  fs_format_manuscript_item(source_map$artifact, source_map$panel_or_sheet)
}

fs_summarize_artifacts_by_script <- function(source_map) {
  if (!nrow(source_map)) {
    return(data.frame(
      script = character(),
      artifact_ids = character(),
      manuscript_outputs = character(),
      final_filenames = character(),
      source_confidence = character(),
      stringsAsFactors = FALSE
    ))
  }

  source_map$script <- fs_basename_script(source_map$generating_script)
  source_map$artifact_display <- fs_artifact_display_name(source_map)
  scripts <- sort(unique(source_map$script))

  rows <- lapply(scripts, function(script) {
    script_rows <- source_map[source_map$script == script, , drop = FALSE]
    data.frame(
      script = script,
      artifact_ids = fs_collapse_unique(script_rows$artifact_id, limit = 20),
      manuscript_outputs = fs_collapse_unique(script_rows$artifact_display, limit = 20),
      final_filenames = fs_collapse_unique(script_rows$final_filename, limit = 10),
      source_confidence = fs_collapse_unique(script_rows$confidence, limit = 5),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

fs_annotate_plan_with_outputs <- function(plan, source_map) {
  artifact_summary <- fs_summarize_artifacts_by_script(source_map)
  plan$artifact_ids <- ""
  plan$manuscript_outputs <- ""
  plan$final_filenames <- ""
  plan$source_confidence <- ""
  plan$manuscript_role <- "upstream_dependency_or_intermediate"

  if (nrow(artifact_summary)) {
    idx <- match(plan$script, artifact_summary$script)
    matched <- !is.na(idx)
    plan$artifact_ids[matched] <- artifact_summary$artifact_ids[idx[matched]]
    plan$manuscript_outputs[matched] <- artifact_summary$manuscript_outputs[idx[matched]]
    plan$final_filenames[matched] <- artifact_summary$final_filenames[idx[matched]]
    plan$source_confidence[matched] <- artifact_summary$source_confidence[idx[matched]]
    plan$manuscript_role[matched] <- "direct_manuscript_output_or_panel"
  }

  plan
}

fs_write_script_index <- function(project_root, script_dir, output_path = NULL) {
  plan <- fs_read_source_pipeline(project_root)
  source_map <- fs_read_artifact_source_map(project_root, required = TRUE)
  script_index <- fs_annotate_plan_with_outputs(plan, source_map)

  if (is.null(output_path)) {
    output_path <- file.path(script_dir, "script_index.tsv")
  }

  utils::write.table(
    script_index,
    output_path,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = ""
  )

  invisible(output_path)
}

fs_sanitize_path_component <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- "unnamed"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}
