# =============================================================================
# manuscript_output_helpers.R
#
# Purpose:
#   Shared output utilities for the original numbered Final_Scripts pipeline.
#   These helpers do not perform any scientific analysis. They only make it easy
#   for each original script to write or copy its own final manuscript figure,
#   table, source-data, or helper object into one standardized output tree.
#
# Why this exists:
#   The legacy scripts already create the correct scientific objects, but many
#   outputs have historical filenames that do not clearly say which manuscript
#   figure/table they support. The final pipeline should make that mapping
#   explicit inside the script that creates each object, rather than relying on a
#   separate post-hoc collector.
#
# Output root:
#   Scripts_2025/Final_Scripts/final_manuscript_objects/
#
# Directory layout written by these helpers:
#   01_main_figures/, 02_extended_data_figures/, 03_main_tables/,
#   04_supplementary_tables/
#     Primary per-artifact folders. These are the most explicit location for
#     reviewers because every copied file sits under its final figure/table and
#     panel or sheet label.
#
#   generated/figure_components/
#     Reference-workflow-style component folders, e.g.
#     generated/figure_components/Figure_3/panel_F/
#
#   generated/main_figures/, generated/extended_figures/, generated/main_tables/,
#   and generated/supplementary_tables/
#     Reference-workflow-style folders for assembled manuscript artifacts.
#     Table source-data companions are additionally mirrored under
#     generated/main_tables/source_data/.
#
#   generated/manual_assembly_reference/
#     Source files for artifacts whose final visual assembly is documented as
#     manual/external in docs/manuscript_artifact_source_map.tsv.
#
#   main_figures/, extended_figures/, main_tables/, supplementary_tables/
#     Frozen/current final assembled artifacts from
#     docs/manuscript_artifact_source_map.tsv when they exist locally. These
#     preserve the final manuscript PDFs/DOCX/XLSX/CSV outputs in the same
#     easy-to-browse style as reproducible_workflow/outputs.
#
# Expected use inside a numbered script:
#   source("Scripts_2025/Final_Scripts/manuscript_output_helpers.R")
#   ms_copy_artifact(
#     source_path = "Final Tables and Figures/example.png",
#     artifact_id = "FIG1A",
#     role = "figure_panel_png",
#     description = "Final plotted component used for Main Figure 1A"
#   )
#
# Relationship to manuscript_artifact_source_map.tsv:
#   The helper reads docs/manuscript_artifact_source_map.tsv when available so
#   the directory names and final figure/table labels come from the same audited
#   provenance map used during validation. Each script still calls the helper
#   explicitly at the point where it creates the figure/table.
# =============================================================================

ms_find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    has_final_scripts <- dir.exists(file.path(current, "Scripts_2025", "Final_Scripts"))
    has_source_map <- file.exists(file.path(current, "docs", "manuscript_artifact_source_map.tsv"))
    if (has_final_scripts && has_source_map) return(current)

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not locate the project root. Run from the project root or a ",
        "subdirectory containing Scripts_2025/Final_Scripts and docs/.",
        call. = FALSE
      )
    }
    current <- parent
  }
}

ms_final_scripts_dir <- function(project_root = ms_find_project_root()) {
  file.path(project_root, "Scripts_2025", "Final_Scripts")
}

ms_output_root <- function(project_root = ms_find_project_root()) {
  root <- file.path(ms_final_scripts_dir(project_root), "final_manuscript_objects")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  ms_initialize_output_tree(root)
  root
}

ms_initialize_output_tree <- function(root) {
  dirs <- c(
    "01_main_figures",
    "02_extended_data_figures",
    "03_main_tables",
    "04_supplementary_tables",
    "05_supplementary_figures",
    "06_helpers_and_source_data",
    "main_figures",
    "extended_figures",
    "main_tables",
    "supplementary_tables",
    "supplementary_figures",
    "frozen/main_figures",
    "frozen/extended_figures",
    "frozen/main_tables",
    "frozen/supplementary_tables",
    "frozen/supplementary_figures",
    "generated/main_figures",
    "generated/extended_figures",
    "generated/figure_components",
    "generated/main_tables",
    "generated/main_tables/source_data",
    "generated/manual_assembly_reference",
    "generated/native_script_runs",
    "generated/supplementary_tables",
    "generated/supplementary_tables/audit",
    "logs",
    "logs/rendered_final_figures"
  )
  invisible(lapply(file.path(root, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))
}

ms_manifest_path <- function(project_root = ms_find_project_root()) {
  file.path(ms_output_root(project_root), "manuscript_direct_output_manifest.tsv")
}

ms_output_index_path <- function(project_root = ms_find_project_root()) {
  file.path(ms_output_root(project_root), "manuscript_output_index.tsv")
}

ms_script_index_path <- function(project_root = ms_find_project_root()) {
  file.path(ms_output_root(project_root), "script_output_index.tsv")
}

ms_output_readme_path <- function(project_root = ms_find_project_root()) {
  file.path(ms_output_root(project_root), "README.md")
}

ms_output_logs_dir <- function(project_root = ms_find_project_root()) {
  logs_dir <- file.path(ms_output_root(project_root), "logs")
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  logs_dir
}

ms_read_source_map <- function(project_root = ms_find_project_root()) {
  path <- file.path(project_root, "docs", "manuscript_artifact_source_map.tsv")
  if (!file.exists(path)) {
    stop("Missing manuscript artifact source map: ", path, call. = FALSE)
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

ms_artifact_metadata <- function(artifact_id, project_root = ms_find_project_root()) {
  source_map <- ms_read_source_map(project_root)
  hit <- source_map[source_map$artifact_id == artifact_id, , drop = FALSE]
  if (!nrow(hit)) {
    stop(
      "artifact_id '", artifact_id,
      "' is not present in docs/manuscript_artifact_source_map.tsv.",
      call. = FALSE
    )
  }
  hit[1, , drop = FALSE]
}

ms_category_dirname <- function(artifact_category) {
  category_map <- c(
    main_figure = "01_main_figures",
    extended_data_figure = "02_extended_data_figures",
    main_table = "03_main_tables",
    supplementary_table = "04_supplementary_tables",
    supplementary_figure = "05_supplementary_figures",
    helper = "06_helpers_and_source_data"
  )

  if (!artifact_category %in% names(category_map)) {
    return(file.path("99_other", ms_slug(artifact_category)))
  }
  unname(category_map[[artifact_category]])
}

ms_slug <- function(x) {
  x <- ifelse(is.na(x) | !nzchar(x), "unspecified", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unspecified")
}

ms_artifact_label <- function(meta) {
  artifact <- meta$artifact[[1]]
  panel <- meta$panel_or_sheet[[1]]
  if (!is.na(panel) && nzchar(panel) && !panel %in% c("all", "all_sheets")) {
    if (nchar(panel) == 1) {
      paste0(artifact, panel)
    } else {
      paste(artifact, ms_slug(panel), sep = "_")
    }
  } else {
    artifact
  }
}

ms_artifact_dir <- function(artifact_id, project_root = ms_find_project_root()) {
  meta <- ms_artifact_metadata(artifact_id, project_root)
  file.path(
    ms_output_root(project_root),
    ms_category_dirname(meta$artifact_category[[1]]),
    ms_slug(meta$artifact[[1]]),
    ms_slug(ms_artifact_label(meta))
  )
}

ms_normalize_source_path <- function(source_path, project_root = ms_find_project_root()) {
  if (file.exists(source_path)) return(normalizePath(source_path, mustWork = TRUE))

  project_relative <- file.path(project_root, source_path)
  if (file.exists(project_relative)) return(normalizePath(project_relative, mustWork = TRUE))

  stop(
    "Manuscript output source does not exist: ", source_path,
    "\nChecked as written and relative to project root: ", project_relative,
    call. = FALSE
  )
}

ms_destination_path <- function(source_path, artifact_id, role, project_root = ms_find_project_root()) {
  meta <- ms_artifact_metadata(artifact_id, project_root)
  destination_dir <- ms_artifact_dir(artifact_id, project_root)
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)

  ext <- tools::file_ext(source_path)
  ext <- if (nzchar(ext)) paste0(".", ext) else ""
  source_stem <- tools::file_path_sans_ext(basename(source_path))
  filename <- paste0(
    ms_slug(ms_artifact_label(meta)),
    "__",
    ms_slug(role),
    "__",
    ms_slug(source_stem),
    ext
  )
  file.path(destination_dir, filename)
}

ms_panel_dirname <- function(panel_or_sheet) {
  panel <- ifelse(is.na(panel_or_sheet) | !nzchar(panel_or_sheet), "all", panel_or_sheet)
  if (panel %in% c("all", "all_sheets")) return("panel_all")
  paste0("panel_", ms_slug(panel))
}

ms_reproducible_style_generated_path <- function(source_path,
                                                 artifact_id,
                                                 role,
                                                 project_root = ms_find_project_root()) {
  meta <- ms_artifact_metadata(artifact_id, project_root)
  generated_root <- file.path(ms_output_root(project_root), "generated")
  ext <- tools::file_ext(source_path)
  ext <- if (nzchar(ext)) paste0(".", ext) else ""
  source_stem <- tools::file_path_sans_ext(basename(source_path))
  filename <- paste0(artifact_id, "__", ms_slug(role), "__", ms_slug(source_stem), ext)

  if (meta$artifact_category[[1]] %in% c("main_figure", "extended_data_figure", "supplementary_figure")) {
    return(file.path(
      generated_root,
      "figure_components",
      ms_slug(meta$artifact[[1]]),
      ms_panel_dirname(meta$panel_or_sheet[[1]]),
      filename
    ))
  }

  if (identical(meta$artifact_category[[1]], "main_table")) {
    if (grepl("source|companion", role, ignore.case = TRUE)) {
      return(file.path(generated_root, "main_tables", "source_data", filename))
    }
    return(file.path(generated_root, "main_tables", filename))
  }

  if (identical(meta$artifact_category[[1]], "supplementary_table")) {
    return(file.path(generated_root, "supplementary_tables", filename))
  }

  file.path(generated_root, "helpers_and_source_data", filename)
}

ms_is_manual_or_external <- function(meta) {
  notes <- if ("notes" %in% names(meta)) meta$notes[[1]] else ""
  source_component <- if ("source_component_or_original_output" %in% names(meta)) {
    meta$source_component_or_original_output[[1]]
  } else {
    ""
  }

  grepl(
    paste(
      c(
        "manual_assembly_required",
        "manual final",
        "manual Word",
        "manual PDF",
        "PowerPoint-only",
        "PowerPoint assembly",
        "external ichorCNA",
        "Lucidchart"
      ),
      collapse = "|"
    ),
    notes,
    ignore.case = TRUE
  ) ||
    grepl(
      paste(
        c(
          "final visual panel created",
          "manual final",
          "manual Word",
          "manual PDF",
          "PowerPoint-only",
          "PowerPoint",
          "external ichorCNA",
          "Lucidchart"
        ),
        collapse = "|"
      ),
      source_component,
      ignore.case = TRUE
    )
}

ms_manual_assembly_reference_path <- function(source_path,
                                              artifact_id,
                                              role,
                                              project_root = ms_find_project_root()) {
  meta <- ms_artifact_metadata(artifact_id, project_root)
  ext <- tools::file_ext(source_path)
  ext <- if (nzchar(ext)) paste0(".", ext) else ""
  source_stem <- tools::file_path_sans_ext(basename(source_path))
  filename <- paste0(
    artifact_id,
    "__",
    ms_slug(ms_artifact_label(meta)),
    "__",
    ms_slug(role),
    "__",
    ms_slug(source_stem),
    ext
  )
  file.path(
    ms_output_root(project_root),
    "generated",
    "manual_assembly_reference",
    ms_slug(ms_artifact_label(meta)),
    filename
  )
}

ms_copy_file_quietly <- function(source_path, destination_path, overwrite = TRUE) {
  dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(source_path, destination_path, overwrite = overwrite)
  if (!ok) stop("Failed to copy manuscript artifact to: ", destination_path, call. = FALSE)
  normalizePath(destination_path, mustWork = TRUE)
}

ms_render_final_figure_previews <- function(project_root = ms_find_project_root()) {
  output_root <- ms_output_root(project_root)
  preview_dir <- file.path(ms_output_logs_dir(project_root), "rendered_final_figures")
  dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

  figure_pdfs <- c(
    list.files(file.path(output_root, "main_figures"), pattern = "[.]pdf$", full.names = TRUE),
    list.files(file.path(output_root, "extended_figures"), pattern = "[.]pdf$", full.names = TRUE)
  )
  figure_pdfs <- sort(unique(figure_pdfs[file.exists(figure_pdfs)]))

  renderer <- if (requireNamespace("pdftools", quietly = TRUE)) {
    "pdftools"
  } else if (nzchar(Sys.which("sips"))) {
    "sips"
  } else {
    "unavailable"
  }

  if (!length(figure_pdfs)) {
    manifest <- data.frame(
      source_pdf = character(),
      rendered_png = character(),
      renderer = character(),
      status = character(),
      message = character(),
      stringsAsFactors = FALSE
    )
  } else {
    manifest_rows <- lapply(figure_pdfs, function(pdf_path) {
      png_path <- file.path(
        preview_dir,
        paste0(tools::file_path_sans_ext(basename(pdf_path)), ".png")
      )
      status <- "pending"
      message <- ""

      if (identical(renderer, "pdftools")) {
        tryCatch(
          {
            pdftools::pdf_convert(
              pdf = pdf_path,
              format = "png",
              pages = 1,
              dpi = 150,
              filenames = png_path,
              verbose = FALSE
            )
            status <- if (file.exists(png_path)) "rendered" else "failed"
          },
          error = function(e) {
            status <<- "failed"
            message <<- conditionMessage(e)
          }
        )
      } else if (identical(renderer, "sips")) {
        result <- tryCatch(
          system2(
            "sips",
            args = c("-s", "format", "png", pdf_path, "--out", png_path),
            stdout = TRUE,
            stderr = TRUE
          ),
          error = function(e) {
            message <<- conditionMessage(e)
            integer(1)
          }
        )
        status <- if (file.exists(png_path)) "rendered" else "failed"
        if (!length(message)) message <- paste(result, collapse = " | ")
      } else {
        status <- "skipped"
        message <- "No PDF preview renderer available. Install pdftools or run on macOS with sips."
      }

      data.frame(
        source_pdf = normalizePath(pdf_path, mustWork = TRUE),
        rendered_png = if (file.exists(png_path)) normalizePath(png_path, mustWork = TRUE) else png_path,
        renderer = renderer,
        status = status,
        message = message,
        stringsAsFactors = FALSE
      )
    })
    manifest <- do.call(rbind, manifest_rows)
  }

  utils::write.table(
    manifest,
    file = file.path(ms_output_logs_dir(project_root), "rendered_final_figures_manifest.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )

  if (identical(renderer, "unavailable")) {
    writeLines(
      c(
        "# Rendered Final Figures",
        "",
        "PNG previews were not generated because no PDF renderer was available.",
        "Install the R package `pdftools` or run on macOS with `sips` to populate this folder."
      ),
      con = file.path(preview_dir, "README.md")
    )
  }

  invisible(manifest)
}

ms_copy_generated_mirror <- function(copied_artifact_path,
                                     artifact_id,
                                     role,
                                     naming_source_path = copied_artifact_path,
                                     project_root = ms_find_project_root(),
                                     overwrite = TRUE) {
  meta <- ms_artifact_metadata(artifact_id, project_root)
  generated_path <- ms_reproducible_style_generated_path(
    source_path = naming_source_path,
    artifact_id = artifact_id,
    role = role,
    project_root = project_root
  )
  copied_paths <- ms_copy_file_quietly(copied_artifact_path, generated_path, overwrite = overwrite)

  if (ms_is_manual_or_external(meta)) {
    manual_reference_path <- ms_manual_assembly_reference_path(
      source_path = naming_source_path,
      artifact_id = artifact_id,
      role = role,
      project_root = project_root
    )
    copied_paths <- c(
      copied_paths,
      ms_copy_file_quietly(copied_artifact_path, manual_reference_path, overwrite = overwrite)
    )
  }

  paste(copied_paths, collapse = "; ")
}

ms_final_category_dirname <- function(artifact_category) {
  final_map <- c(
    main_figure = "main_figures",
    extended_data_figure = "extended_figures",
    main_table = "main_tables",
    supplementary_table = "supplementary_tables",
    supplementary_figure = "supplementary_figures"
  )
  if (!artifact_category %in% names(final_map)) return(NA_character_)
  unname(final_map[[artifact_category]])
}

ms_existing_path_from_map <- function(path_value, project_root = ms_find_project_root()) {
  if (is.na(path_value) || !nzchar(path_value)) return(NA_character_)
  if (file.exists(path_value)) return(normalizePath(path_value, mustWork = TRUE))
  project_relative <- file.path(project_root, path_value)
  if (file.exists(project_relative)) return(normalizePath(project_relative, mustWork = TRUE))
  NA_character_
}

ms_split_map_paths <- function(path_value) {
  if (is.na(path_value) || !nzchar(path_value)) return(character())
  trimws(unlist(strsplit(path_value, ";", fixed = TRUE), use.names = FALSE))
}

ms_resolve_candidate_sources <- function(values, project_root = ms_find_project_root()) {
  candidates <- unique(unlist(lapply(values, ms_split_map_paths), use.names = FALSE))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  resolved <- vapply(candidates, ms_existing_path_from_map, character(1), project_root = project_root)
  unique(resolved[!is.na(resolved) & nzchar(resolved)])
}

ms_write_source_availability_manifest <- function(source_map, project_root = ms_find_project_root()) {
  rows <- lapply(seq_len(nrow(source_map)), function(i) {
    meta <- source_map[i, , drop = FALSE]
    source_values <- c(
      meta$current_final_artifact_path[[1]],
      meta$manuscript_source_path[[1]],
      meta$source_component_or_original_output[[1]],
      meta$input_files[[1]],
      meta$intermediate_files[[1]]
    )
    resolved <- ms_resolve_candidate_sources(source_values, project_root)
    data.frame(
      artifact_id = meta$artifact_id[[1]],
      artifact_category = meta$artifact_category[[1]],
      artifact = meta$artifact[[1]],
      panel_or_sheet = meta$panel_or_sheet[[1]],
      confidence = meta$confidence[[1]],
      expected_source = paste(source_values[!is.na(source_values) & nzchar(source_values)], collapse = "; "),
      resolved_source_count = length(resolved),
      resolved_sources = paste(resolved, collapse = "; "),
      inferred_status_if_not_generated = meta$notes[[1]],
      stringsAsFactors = FALSE
    )
  })
  source_availability <- do.call(rbind, rows)
  utils::write.table(
    source_availability,
    file = file.path(ms_output_logs_dir(project_root), "source_availability_manifest.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )
  invisible(source_availability)
}

ms_copy_reference_log_if_present <- function(source_path, destination_name, project_root = ms_find_project_root()) {
  if (!file.exists(source_path)) return(FALSE)
  file.copy(
    source_path,
    file.path(ms_output_logs_dir(project_root), destination_name),
    overwrite = TRUE
  )
}

ms_mirror_project_provenance_logs <- function(project_root = ms_find_project_root()) {
  docs_dir <- file.path(project_root, "docs")
  provenance_files <- c(
    "artifact_dependency_manifest.tsv",
    "artifact_regeneration_status.tsv",
    "cache_rebuild_matrix.tsv",
    "cached_rds_rebuild_manifest.tsv",
    "frozen_model_artifact_manifest.tsv",
    "frozen_model_scoring_checks.tsv",
    "frozen_model_scoring_model_schema.tsv",
    "frozen_model_scoring_readiness.tsv",
    "frozen_model_scoring_replay_check.tsv",
    "raw_to_final_pipeline_status.tsv",
    "source_pipeline_plan.tsv",
    "test_cohort_update_frozen_training_check.tsv",
    "test_cohort_update_plan.tsv"
  )
  copied <- vapply(
    provenance_files,
    function(filename) {
      ms_copy_reference_log_if_present(file.path(docs_dir, filename), filename, project_root)
    },
    logical(1)
  )
  data.frame(
    file = names(copied),
    copied = unname(copied),
    stringsAsFactors = FALSE
  )
}

ms_copy_current_final_artifact <- function(artifact_id,
                                           project_root = ms_find_project_root(),
                                           overwrite = TRUE) {
  meta <- ms_artifact_metadata(artifact_id, project_root)
  final_dir <- ms_final_category_dirname(meta$artifact_category[[1]])
  if (is.na(final_dir)) return(invisible(NA_character_))

  source <- ms_existing_path_from_map(meta$current_final_artifact_path[[1]], project_root)
  if (is.na(source)) {
    source <- ms_existing_path_from_map(meta$manuscript_source_path[[1]], project_root)
  }
  if (is.na(source)) return(invisible(NA_character_))

  final_filename <- meta$final_filename[[1]]
  if (is.na(final_filename) || !nzchar(final_filename)) {
    final_filename <- basename(source)
  }

  destination <- file.path(ms_output_root(project_root), final_dir, basename(final_filename))
  frozen_destination <- file.path(ms_output_root(project_root), "frozen", final_dir, basename(final_filename))
  generated_destination <- file.path(ms_output_root(project_root), "generated", final_dir, basename(final_filename))
  ms_copy_file_quietly(source, destination, overwrite = overwrite)
  ms_copy_file_quietly(source, frozen_destination, overwrite = overwrite)
  ms_copy_file_quietly(source, generated_destination, overwrite = overwrite)
  invisible(normalizePath(destination, mustWork = TRUE))
}

ms_record_manifest <- function(row, project_root = ms_find_project_root()) {
  manifest <- ms_manifest_path(project_root)
  dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
  row <- as.data.frame(row, stringsAsFactors = FALSE)

  if (file.exists(manifest)) {
    previous <- utils::read.delim(manifest, stringsAsFactors = FALSE, check.names = FALSE)
    all_columns <- unique(c(names(previous), names(row)))
    for (col in setdiff(all_columns, names(previous))) previous[[col]] <- ""
    for (col in setdiff(all_columns, names(row))) row[[col]] <- ""
    previous <- previous[, all_columns, drop = FALSE]
    row <- row[, all_columns, drop = FALSE]
    utils::write.table(
      rbind(previous, row),
      file = manifest,
      sep = "\t",
      row.names = FALSE,
      quote = TRUE,
      col.names = TRUE
    )
    return(invisible(manifest))
  }

  utils::write.table(
    row,
    file = manifest,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    append = FALSE,
    col.names = TRUE
  )
  invisible(manifest)
}

ms_strip_external_script_prefix <- function(script_name) {
  sub(
    "^External ichorCNA workflow outside this repo; supporting repo script: ",
    "",
    script_name
  )
}

ms_collapse_existing_paths <- function(paths) {
  split_paths <- unlist(strsplit(paths[!is.na(paths) & nzchar(paths)], ";", fixed = TRUE), use.names = FALSE)
  split_paths <- trimws(split_paths)
  paths <- unique(split_paths[nzchar(split_paths) & file.exists(split_paths)])
  if (!length(paths)) return("")
  paste(normalizePath(paths, mustWork = TRUE), collapse = "; ")
}

ms_refresh_manual_assembly_references <- function(output_index, project_root = ms_find_project_root()) {
  reference_root <- file.path(ms_output_root(project_root), "generated", "manual_assembly_reference")
  dir.create(reference_root, recursive = TRUE, showWarnings = FALSE)

  if (!nrow(output_index) || !"manual_or_external_final_assembly" %in% names(output_index)) {
    manifest <- data.frame(
      artifact_id = character(),
      final_label = character(),
      source_path = character(),
      reference_path = character(),
      status = character(),
      stringsAsFactors = FALSE
    )
  } else {
    manual_rows <- output_index[
      !is.na(output_index$manual_or_external_final_assembly) &
        output_index$manual_or_external_final_assembly,
      ,
      drop = FALSE
    ]

    manifest_rows <- lapply(seq_len(nrow(manual_rows)), function(i) {
      row <- manual_rows[i, , drop = FALSE]
      paths <- ms_split_map_paths(row$explicit_traceability_paths[[1]])
      paths <- unique(paths[file.exists(paths)])
      if (!length(paths)) {
        return(data.frame(
          artifact_id = row$artifact_id[[1]],
          final_label = row$final_label[[1]],
          source_path = "",
          reference_path = "",
          status = "missing_traceability_source",
          stringsAsFactors = FALSE
        ))
      }

      do.call(rbind, lapply(paths, function(source_path) {
        destination <- file.path(
          reference_root,
          ms_slug(row$final_label[[1]]),
          paste0(row$artifact_id[[1]], "__", ms_slug(tools::file_path_sans_ext(basename(source_path))), ".", tools::file_ext(source_path))
        )
        copied <- tryCatch(
          {
            ms_copy_file_quietly(source_path, destination, overwrite = TRUE)
            TRUE
          },
          error = function(e) FALSE
        )
        data.frame(
          artifact_id = row$artifact_id[[1]],
          final_label = row$final_label[[1]],
          source_path = normalizePath(source_path, mustWork = TRUE),
          reference_path = if (file.exists(destination)) normalizePath(destination, mustWork = TRUE) else destination,
          status = if (copied && file.exists(destination)) "copied" else "failed",
          stringsAsFactors = FALSE
        )
      }))
    })
    manifest <- do.call(rbind, manifest_rows)
  }

  utils::write.table(
    manifest,
    file = file.path(ms_output_logs_dir(project_root), "manual_assembly_reference_manifest.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )
  invisible(manifest)
}

ms_refresh_native_script_run_snapshots <- function(output_index, project_root = ms_find_project_root()) {
  native_root <- file.path(ms_output_root(project_root), "generated", "native_script_runs")
  dir.create(native_root, recursive = TRUE, showWarnings = FALSE)

  snapshot_specs <- list(
    FIG1A_swim_plot = c("FIG1A", "STABLE1")
  )

  manifest_rows <- list()
  for (snapshot_name in names(snapshot_specs)) {
    snapshot_dir <- file.path(native_root, snapshot_name)
    dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
    artifact_ids <- snapshot_specs[[snapshot_name]]
    rows <- output_index[output_index$artifact_id %in% artifact_ids, , drop = FALSE]

    for (i in seq_len(nrow(rows))) {
      row <- rows[i, , drop = FALSE]
      paths <- ms_split_map_paths(row$explicit_traceability_paths[[1]])
      paths <- unique(paths[file.exists(paths)])
      if (!length(paths)) {
        manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
          snapshot = snapshot_name,
          artifact_id = row$artifact_id[[1]],
          final_label = row$final_label[[1]],
          source_path = "",
          snapshot_path = "",
          status = "missing_traceability_source",
          stringsAsFactors = FALSE
        )
        next
      }

      for (source_path in paths) {
        ext <- tools::file_ext(source_path)
        ext <- if (nzchar(ext)) paste0(".", ext) else ""
        source_stem <- tools::file_path_sans_ext(basename(source_path))
        destination <- file.path(
          snapshot_dir,
          paste0(row$artifact_id[[1]], "__", ms_slug(source_stem), ext)
        )
        copied <- tryCatch(
          {
            ms_copy_file_quietly(source_path, destination, overwrite = TRUE)
            TRUE
          },
          error = function(e) FALSE
        )
        manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
          snapshot = snapshot_name,
          artifact_id = row$artifact_id[[1]],
          final_label = row$final_label[[1]],
          source_path = normalizePath(source_path, mustWork = TRUE),
          snapshot_path = if (file.exists(destination)) normalizePath(destination, mustWork = TRUE) else destination,
          status = if (copied && file.exists(destination)) "copied" else "failed",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  manifest <- if (length(manifest_rows)) {
    do.call(rbind, manifest_rows)
  } else {
    data.frame(
      snapshot = character(),
      artifact_id = character(),
      final_label = character(),
      source_path = character(),
      snapshot_path = character(),
      status = character(),
      stringsAsFactors = FALSE
    )
  }

  utils::write.table(
    manifest,
    file = file.path(ms_output_logs_dir(project_root), "native_script_runs_manifest.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )
  invisible(manifest)
}

ms_lookup_current_final_path <- function(meta_row, project_root = ms_find_project_root()) {
  final_dir <- ms_final_category_dirname(meta_row$artifact_category[[1]])
  final_filename <- meta_row$final_filename[[1]]

  if (!is.na(final_dir) && !is.na(final_filename) && nzchar(final_filename)) {
    organized_final <- file.path(ms_output_root(project_root), final_dir, basename(final_filename))
    if (file.exists(organized_final)) {
      return(normalizePath(organized_final, mustWork = TRUE))
    }
  }

  source <- ms_existing_path_from_map(meta_row$current_final_artifact_path[[1]], project_root)
  if (!is.na(source)) return(source)

  source <- ms_existing_path_from_map(meta_row$manuscript_source_path[[1]], project_root)
  if (!is.na(source)) return(source)

  ""
}

ms_optional_tsv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

ms_first_value <- function(data, column, default = "") {
  if (!nrow(data) || !column %in% names(data)) return(default)
  value <- data[[column]][[1]]
  if (is.na(value)) default else as.character(value)
}

ms_count_lines <- function(values) {
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return("- none")
  counts <- sort(table(values), decreasing = TRUE)
  paste0("- ", names(counts), ": ", as.integer(counts))
}

ms_write_output_readme <- function(output_index, script_index, project_root = ms_find_project_root()) {
  readme <- ms_output_readme_path(project_root)
  run_stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

  direct_count <- sum(output_index$organization_status == "organized", na.rm = TRUE)
  total_count <- nrow(output_index)
  direct_line <- paste0("- Organized mapped artifacts: ", direct_count, " / ", total_count)

  status_column <- if ("regeneration_status" %in% names(output_index)) output_index$regeneration_status else character()
  validation_column <- if ("validation_status" %in% names(output_index)) output_index$validation_status else character()

  caveat_statuses <- c(
    "manual_assembly_required",
    "requires_cached_intermediate",
    "generated_from_code_value_difference"
  )
  caveats <- output_index[
    (!is.na(output_index$manual_or_external_final_assembly) & output_index$manual_or_external_final_assembly) |
      (!is.na(status_column) & status_column %in% caveat_statuses),
    ,
    drop = FALSE
  ]

  lines <- c(
    "# Final manuscript objects",
    "",
    paste0("Generated by `manuscript_output_helpers.R` on ", run_stamp, "."),
    "",
    "This directory is written directly by the numbered scripts in `Scripts_2025/Final_Scripts`.",
    "It is intended to mirror the easy manuscript-order output structure from `reproducible_workflow/outputs`, while keeping the original scripts as the scientific source of truth.",
    "",
    "## Quick-use folders",
    "",
    "- `main_figures/`: current final assembled main figure PDFs.",
    "- `extended_figures/`: current final assembled Extended Data figure PDFs.",
    "- `main_tables/`: current final main table DOCX/PDF artifacts.",
    "- `supplementary_tables/`: current final supplementary table CSV/XLSX artifacts.",
    "",
    "## Script-owned generated components",
    "",
    "- `01_main_figures/`, `02_extended_data_figures/`, `03_main_tables/`, and `04_supplementary_tables/` contain explicit per-artifact traceability folders.",
    "- `generated/main_figures/`, `generated/extended_figures/`, `generated/figure_components/`, `generated/main_tables/`, and `generated/supplementary_tables/` mirror the reference workflow style for fast manuscript updates.",
    "- `generated/main_tables/source_data/` and `generated/manual_assembly_reference/` preserve source-data companions and manual/external assembly inputs when applicable.",
    "- `logs/manual_assembly_reference_manifest.tsv` lists the manual/external artifact references copied into `generated/manual_assembly_reference/`.",
    "- `generated/native_script_runs/` contains native-run-style snapshots for scripts where that layout helps bridge to the reference workflow output.",
    "",
    "## Index files",
    "",
    "- `manuscript_output_index.tsv`: one row per mapped final manuscript artifact.",
    "- `script_output_index.tsv`: one row per numbered script that owns a mapped final output.",
    "- `manuscript_direct_output_manifest.tsv`: every script-recorded copy/save event with source, destination, checksum, and script name.",
    "- `manuscript_output_validation_report.tsv`: validation results from `validate_manuscript_outputs.R`.",
    "- `logs/`: reference-workflow-style copies of the output index, generation manifest, source availability manifest, dependency/cache manifests, raw-to-final status, rendered final-figure previews, validation report, and validation summary.",
    "",
    "## Output status",
    "",
    direct_line,
    "",
    "## Regeneration status counts",
    "",
    ms_count_lines(status_column),
    "",
    "## Validation status counts",
    "",
    ms_count_lines(validation_column),
    "",
    "## Manual, External, Or Cached Items",
    ""
  )

  if (nrow(caveats)) {
    caveat_lines <- paste0(
      "- ",
      caveats$artifact_id,
      " (",
      caveats$final_label,
      "): ",
      if ("regeneration_status" %in% names(caveats)) caveats$regeneration_status else "manual_or_external",
      ifelse(
        caveats$manual_or_external_final_assembly,
        "; manual/external final assembly is documented",
        ""
      )
    )
  } else {
    caveat_lines <- "- none"
  }

  lines <- c(
    lines,
    caveat_lines,
    "",
    "## New Test-Cohort Data Policy",
    "",
    "New test-cohort samples should be processed through the upstream feature scripts and scored with preserved training-derived models and thresholds.",
    "Do not rerun cache-sensitive model training during routine test-cohort expansion unless the training cohort is intentionally redefined and manuscript values are being revalidated.",
    "",
    "## Source Of Truth",
    "",
    "`docs/manuscript_artifact_source_map.tsv` remains the authoritative crosswalk from final manuscript item to original script section, historical filename, final filename, and caveats."
  )

  writeLines(lines, con = readme)
  invisible(readme)
}

ms_write_output_index <- function(project_root = ms_find_project_root()) {
  source_map <- ms_read_source_map(project_root)
  invisible(lapply(unique(source_map$artifact_id), function(artifact_id) {
    ms_copy_current_final_artifact(
      artifact_id = artifact_id,
      project_root = project_root,
      overwrite = TRUE
    )
  }))

  regeneration <- ms_optional_tsv(file.path(project_root, "docs", "artifact_regeneration_status.tsv"))
  validation <- ms_optional_tsv(file.path(project_root, "reproducible_workflow", "outputs", "logs", "validation_report.tsv"))
  manifest <- if (file.exists(ms_manifest_path(project_root))) {
    utils::read.delim(ms_manifest_path(project_root), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }

  index_rows <- lapply(seq_len(nrow(source_map)), function(i) {
    meta <- source_map[i, , drop = FALSE]
    artifact_id <- meta$artifact_id[[1]]
    manifest_hit <- if (nrow(manifest) && "artifact_id" %in% names(manifest)) {
      manifest[manifest$artifact_id == artifact_id, , drop = FALSE]
    } else {
      data.frame()
    }

    explicit_paths <- if (nrow(manifest_hit) && "destination_path" %in% names(manifest_hit)) {
      ms_collapse_existing_paths(manifest_hit$destination_path)
    } else {
      ""
    }

    generated_paths <- if (nrow(manifest_hit) && "generated_mirror_path" %in% names(manifest_hit)) {
      ms_collapse_existing_paths(manifest_hit$generated_mirror_path)
    } else {
      ""
    }

    final_path <- ms_lookup_current_final_path(meta, project_root)
    regeneration_hit <- if (nrow(regeneration) && "artifact_id" %in% names(regeneration)) {
      regeneration[regeneration$artifact_id == artifact_id, , drop = FALSE]
    } else {
      data.frame()
    }
    validation_hit <- if (nrow(validation) && "artifact_id" %in% names(validation)) {
      validation[validation$artifact_id == artifact_id, , drop = FALSE]
    } else {
      data.frame()
    }
    validation_status_column <- if ("validation_status" %in% names(validation_hit)) {
      "validation_status"
    } else if ("status" %in% names(validation_hit)) {
      "status"
    } else {
      ""
    }
    generated_script <- ms_strip_external_script_prefix(meta$generating_script[[1]])
    notes <- meta$notes[[1]]
    manual_or_external <- ms_is_manual_or_external(meta)

    status <- if (nzchar(explicit_paths) && nzchar(generated_paths)) {
      "organized"
    } else if (nzchar(final_path)) {
      "final_artifact_preserved_only"
    } else {
      "missing_organized_output"
    }

    data.frame(
      artifact_id = artifact_id,
      artifact_category = meta$artifact_category[[1]],
      artifact = meta$artifact[[1]],
      panel_or_sheet = meta$panel_or_sheet[[1]],
      final_label = ms_artifact_label(meta),
      final_filename = meta$final_filename[[1]],
      regeneration_status = ms_first_value(regeneration_hit, "regeneration_status"),
      validation_status = if (nzchar(validation_status_column)) {
        ms_first_value(validation_hit, validation_status_column)
      } else {
        ""
      },
      validation_scope = ms_first_value(validation_hit, "validation_scope"),
      source_confidence = ms_first_value(regeneration_hit, "confidence", meta$confidence[[1]]),
      generating_script = generated_script,
      source_component_or_original_output = meta$source_component_or_original_output[[1]],
      final_assembled_artifact_path = final_path,
      explicit_traceability_paths = explicit_paths,
      generated_mirror_paths = generated_paths,
      manuscript_notes = notes,
      manual_or_external_final_assembly = manual_or_external,
      organization_status = status,
      stringsAsFactors = FALSE
    )
  })

  output_index <- do.call(rbind, index_rows)
  utils::write.table(
    output_index,
    file = ms_output_index_path(project_root),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )

  script_groups <- split(output_index, output_index$generating_script)
  script_index <- do.call(rbind, lapply(names(script_groups), function(script_name) {
    rows <- script_groups[[script_name]]
    data.frame(
      generating_script = script_name,
      n_artifacts = nrow(rows),
      artifacts = paste(rows$artifact_id, collapse = "; "),
      final_labels = paste(rows$final_label, collapse = "; "),
      has_manual_or_external_outputs = any(rows$manual_or_external_final_assembly),
      stringsAsFactors = FALSE
    )
  }))

  script_index <- script_index[order(script_index$generating_script), , drop = FALSE]
  utils::write.table(
    script_index,
    file = ms_script_index_path(project_root),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )

  ms_write_output_readme(output_index, script_index, project_root)

  logs_dir <- ms_output_logs_dir(project_root)
  utils::write.table(
    output_index,
    file = file.path(logs_dir, "output_index.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )
  utils::write.table(
    script_index,
    file = file.path(logs_dir, "script_output_index.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    col.names = TRUE
  )
  if (file.exists(ms_manifest_path(project_root))) {
    file.copy(
      ms_manifest_path(project_root),
      file.path(logs_dir, "generation_manifest.tsv"),
      overwrite = TRUE
    )
  }
  ms_write_source_availability_manifest(source_map, project_root)
  ms_mirror_project_provenance_logs(project_root)
  ms_refresh_manual_assembly_references(output_index, project_root)
  ms_refresh_native_script_run_snapshots(output_index, project_root)
  ms_render_final_figure_previews(project_root)

  invisible(list(output_index = ms_output_index_path(project_root), script_index = ms_script_index_path(project_root)))
}

ms_copy_artifact <- function(source_path,
                             artifact_id,
                             role,
                             description,
                             script_name = NULL,
                             project_root = ms_find_project_root(),
                             overwrite = TRUE) {
  source_abs <- ms_normalize_source_path(source_path, project_root)
  destination <- ms_destination_path(source_abs, artifact_id, role, project_root)
  if (file.exists(destination) && !overwrite) {
    stop("Destination already exists and overwrite = FALSE: ", destination, call. = FALSE)
  }

  ok <- file.copy(source_abs, destination, overwrite = overwrite)
  if (!ok) stop("Failed to copy manuscript artifact to: ", destination, call. = FALSE)
  generated_destination <- ms_copy_generated_mirror(
    copied_artifact_path = destination,
    artifact_id = artifact_id,
    role = role,
    naming_source_path = source_abs,
    project_root = project_root,
    overwrite = overwrite
  )
  current_final_destination <- ms_copy_current_final_artifact(
    artifact_id = artifact_id,
    project_root = project_root,
    overwrite = TRUE
  )

  meta <- ms_artifact_metadata(artifact_id, project_root)
  md5 <- unname(tools::md5sum(destination))
  ms_record_manifest(
    list(
      recorded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      artifact_id = artifact_id,
      artifact_category = meta$artifact_category[[1]],
      artifact = meta$artifact[[1]],
      panel_or_sheet = meta$panel_or_sheet[[1]],
      role = role,
      description = description,
      source_path = source_abs,
      destination_path = normalizePath(destination, mustWork = TRUE),
      generated_mirror_path = generated_destination,
      current_final_artifact_path = ifelse(
        is.na(current_final_destination),
        "",
        current_final_destination
      ),
      md5 = md5,
      script_name = ifelse(is.null(script_name), "", script_name)
    ),
    project_root = project_root
  )
  ms_write_output_index(project_root)

  message("Manuscript output recorded: ", artifact_id, " -> ", destination)
  invisible(destination)
}

ms_save_plot <- function(plot,
                         filename,
                         artifact_id,
                         role,
                         description,
                         width,
                         height,
                         dpi = 500,
                         units = "in",
                         script_name = NULL,
                         project_root = ms_find_project_root(),
                         ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for ms_save_plot().", call. = FALSE)
  }

  destination <- ms_destination_path(filename, artifact_id, role, project_root)
  ggplot2::ggsave(
    filename = destination,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = units,
    ...
  )

  meta <- ms_artifact_metadata(artifact_id, project_root)
  generated_destination <- ms_copy_generated_mirror(
    copied_artifact_path = destination,
    artifact_id = artifact_id,
    role = role,
    naming_source_path = filename,
    project_root = project_root,
    overwrite = TRUE
  )
  current_final_destination <- ms_copy_current_final_artifact(
    artifact_id = artifact_id,
    project_root = project_root,
    overwrite = TRUE
  )
  md5 <- unname(tools::md5sum(destination))
  ms_record_manifest(
    list(
      recorded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      artifact_id = artifact_id,
      artifact_category = meta$artifact_category[[1]],
      artifact = meta$artifact[[1]],
      panel_or_sheet = meta$panel_or_sheet[[1]],
      role = role,
      description = description,
      source_path = "",
      destination_path = normalizePath(destination, mustWork = TRUE),
      generated_mirror_path = generated_destination,
      current_final_artifact_path = ifelse(
        is.na(current_final_destination),
        "",
        current_final_destination
      ),
      md5 = md5,
      script_name = ifelse(is.null(script_name), "", script_name)
    ),
    project_root = project_root
  )
  ms_write_output_index(project_root)

  message("Manuscript plot saved: ", artifact_id, " -> ", destination)
  invisible(destination)
}
