#!/usr/bin/env Rscript

###############################################################################
# 5_2_Build_Figure_Source_Data_Workbooks.R
#
# Purpose and pipeline role
#   Build the two Excel source-data workbooks from the panel CSVs and manifests
#   produced by 5_1_Export_Locked_Figure_Source_Data.R. This script does not
#   recalculate figure statistics, refit models, or redraw panels. It does
#   prepare the CSV fields for Excel as described below, so the written cells
#   are not necessarily byte-for-byte character copies of the CSV text.
#   Script 5_3_Rebuild_Figure_Source_Data_Workbooks.R runs 5_1, this script,
#   and 5_4_Validate_Figure_Source_Data_Workbooks.R in sequence.
#
# Inputs
#   - --manifest: source_data_workbook_manifest.csv from 5_1
#   - --audit: source_data_audit.csv from 5_1
#   - Panel CSVs referenced by the manifest
#   - id_map.rds for a final check that original patient identifiers are absent
#
# Outputs
#   - <output-dir>/Source_Data_Main_Figures.xlsx
#   - <output-dir>/Source_Data_Extended_Data_Figures.xlsx
#   - Copies of the extended-data workbook under final_manuscript_objects/
#     02_extended_data_figures and 06_helpers_and_source_data
#   Existing files at all four destinations are overwritten.
#
# Workbook preparation and validation
#   - Requires exactly 74 unique panel sheets: 18 main and 56 extended-data.
#   - Converts complete TRUE/FALSE character columns to logical values and
#     numeric-looking columns to numeric values, while preserving identifiers
#     with leading zeroes as character values.
#   - Applies the manuscript cohort labels through
#     relabel_publication_cohorts().
#   - Removes calendar-date columns, then stops if date-like values or original
#     patient identifiers remain.
#   - Repairs dangling OOXML drawing relationships emitted by some openxlsx
#     versions so the files open in strict spreadsheet readers.
#
# How to run
#   Normally run through 5_3. Direct use requires:
#     Rscript Scripts_2025/Final_Scripts/5_2_Build_Figure_Source_Data_Workbooks.R \
#       --manifest <manifest.csv> --audit <audit.csv> --output-dir <directory>
#
# Required packages: readr, dplyr, openxlsx, xml2, and zip.
###############################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(openxlsx)
})

.publication_export_helper <- file.path(
  "Scripts_2025", "Final_Scripts", "publication_export_helpers.R"
)
if (!file.exists(.publication_export_helper)) {
  .publication_export_helper <- "publication_export_helpers.R"
}
source(.publication_export_helper)
rm(.publication_export_helper)

id_map_path <- "id_map.rds"
if (!file.exists(id_map_path)) stop("Missing deidentification authority: id_map.rds", call. = FALSE)
publication_id_map <- readRDS(id_map_path) %>%
  transmute(raw_id = as.character(.data$Patient), public_id = as.character(.data$New_ID))

sanitize_unused_drawing_relationships <- function(xlsx) {
  # openxlsx 4.2.8.1 can emit worksheet relationships to drawing/VML parts
  # that were never created when a workbook contains no images or comments.
  # Excel may repair these silently, while stricter OOXML readers fail. Remove
  # only drawing relationships whose target part is absent from the archive.
  if (!requireNamespace("xml2", quietly = TRUE) ||
      !requireNamespace("zip", quietly = TRUE)) {
    stop("Workbook OOXML cleanup requires the xml2 and zip packages.", call. = FALSE)
  }

  unpack_dir <- tempfile("source_data_xlsx_")
  dir.create(unpack_dir)
  on.exit(unlink(unpack_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(xlsx, exdir = unpack_dir)

  rel_dir <- file.path(unpack_dir, "xl", "worksheets", "_rels")
  rel_files <- list.files(rel_dir, pattern = "[.]rels$", full.names = TRUE)
  removed <- 0L
  for (rel_file in rel_files) {
    doc <- xml2::read_xml(rel_file)
    relationships <- xml2::xml_find_all(doc, "//*[local-name()='Relationship']")
    for (relationship in relationships) {
      type <- xml2::xml_attr(relationship, "Type")
      target <- xml2::xml_attr(relationship, "Target")
      if (!grepl("/(?:drawing|vmlDrawing)$", type, perl = TRUE)) next
      target_path <- normalizePath(
        file.path(unpack_dir, "xl", "worksheets", target),
        mustWork = FALSE
      )
      if (!file.exists(target_path)) {
        xml2::xml_remove(relationship)
        removed <- removed + 1L
      }
    }
    xml2::write_xml(doc, rel_file, options = "format")
  }

  repacked <- tempfile(fileext = ".xlsx")
  files <- list.files(unpack_dir, recursive = TRUE, all.files = TRUE,
                      no.. = TRUE, include.dirs = FALSE)
  zip::zip(
    repacked, files = files, root = unpack_dir, mode = "mirror",
    include_directories = FALSE
  )
  if (!file.rename(repacked, xlsx)) {
    if (!file.copy(repacked, xlsx, overwrite = TRUE)) {
      stop("Could not replace workbook after OOXML relationship cleanup: ", xlsx,
           call. = FALSE)
    }
    unlink(repacked)
  }
  message("Removed ", removed, " dangling drawing relationship(s) from ", basename(xlsx))
}

parse_args <- function(args) {
  if (!length(args)) {
    stop(
      "Missing arguments: --manifest, --audit, --output-dir.",
      call. = FALSE
    )
  }
  if (length(args) %% 2L != 0L) stop("Arguments must be --key value pairs.", call. = FALSE)
  keys <- sub("^--", "", args[seq(1L, length(args), by = 2L)])
  values <- args[seq(2L, length(args), by = 2L)]
  out <- setNames(as.list(values), keys)
  required <- c("manifest", "audit", "output-dir")
  missing <- required[!required %in% names(out)]
  if (length(missing)) stop("Missing argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)
  out
}

coerce_column <- function(x) {
  x <- as.character(x)
  nonmissing <- !is.na(x) & nzchar(x)
  if (!any(nonmissing)) return(x)
  values <- x[nonmissing]
  if (all(values %in% c("TRUE", "FALSE"))) return(as.logical(x))
  numeric_pattern <- "^-?(?:0|[1-9][0-9]*)(?:[.][0-9]+)?(?:[eE][+-]?[0-9]+)?$"
  protected_zero <- "^-?0[0-9]+"
  if (all(grepl(numeric_pattern, values, perl = TRUE)) && !any(grepl(protected_zero, values, perl = TRUE))) {
    return(suppressWarnings(as.numeric(x)))
  }
  x
}

read_panel_for_excel <- function(path) {
  if (!file.exists(path)) stop("Panel CSV is missing: ", path, call. = FALSE)
  dat <- read.csv(
    path, check.names = FALSE, colClasses = "character",
    na.strings = c(""), stringsAsFactors = FALSE
  )
  if (!nrow(dat) || !ncol(dat)) stop("Panel CSV is empty: ", path, call. = FALSE)
  dat[] <- lapply(dat, coerce_column)
  dat <- relabel_publication_cohorts(dat, paste0("Workbook panel ", basename(path)))

  # Calendar dates are required upstream to derive relative-time endpoints but
  # must never be distributed in supplementary or figure source-data files.
  date_columns <- grepl("(^|[._ -])date($|[._ -])", names(dat), ignore.case = TRUE)
  if (any(date_columns)) {
    message(
      "Removing calendar-date column(s) from ", basename(path), ": ",
      paste(names(dat)[date_columns], collapse = ", ")
    )
    dat <- dat[, !date_columns, drop = FALSE]
  }
  if (any(grepl("(^|[._ -])date($|[._ -])", names(dat), ignore.case = TRUE))) {
    stop("Calendar-date column remained after publication sanitization: ", path, call. = FALSE)
  }
  forbidden_headers <- c(
    "figure_panel", "locked_figure", "source_table", "source_note",
    "source_call_path", "record_type"
  )
  present_forbidden <- intersect(names(dat), forbidden_headers)
  if (length(present_forbidden)) {
    stop(
      "Administrative column(s) remained in ", path, ": ",
      paste(present_forbidden, collapse = ", "), call. = FALSE
    )
  }
  character_columns <- dat[vapply(dat, is.character, logical(1))]
  remaining_raw <- publication_id_map$raw_id[vapply(
    publication_id_map$raw_id,
    function(raw_id) any(vapply(
      character_columns,
      function(column) any(grepl(
        paste0("(?<![A-Za-z0-9])", raw_id, "(?![A-Za-z0-9])"),
        column, perl = TRUE
      ), na.rm = TRUE),
      logical(1)
    )),
    logical(1)
  )]
  if (length(remaining_raw)) {
    stop(
      "Original patient identifier(s) remained in ", path, ": ",
      paste(remaining_raw, collapse = ", "), call. = FALSE
    )
  }
  date_pattern <- paste0(
    "(^|[^0-9])(?:[12][0-9]{3}[-/][01]?[0-9][-/][0-3]?[0-9]",
    "|[01]?[0-9][-/][0-3]?[0-9][-/][12][0-9]{3})([^0-9]|$)"
  )
  if (any(vapply(
    character_columns,
    function(column) any(grepl(date_pattern, column, perl = TRUE), na.rm = TRUE),
    logical(1)
  ))) {
    stop("Calendar-date value remained in publication panel: ", path, call. = FALSE)
  }
  dat
}

header_style <- createStyle(
  fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold",
  halign = "left", valign = "center", wrapText = TRUE
)
label_style <- createStyle(
  fgFill = "#D9EAF7", fontColour = "#17365D", textDecoration = "bold",
  valign = "top", wrapText = TRUE
)
title_style <- createStyle(
  fgFill = "#17365D", fontColour = "#FFFFFF", textDecoration = "bold",
  fontSize = 16, valign = "center"
)

style_data_sheet <- function(wb, sheet, dat) {
  addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(dat)), gridExpand = TRUE)
  setRowHeights(wb, sheet, rows = 1, heights = 32)
  freezePane(wb, sheet, firstRow = TRUE)
  pageSetup(
    wb, sheet, orientation = "landscape", scale = 100,
    fitToWidth = 1, fitToHeight = FALSE, paperSize = 9,
    printTitleRows = 1
  )
  for (column in seq_len(ncol(dat))) {
    header <- names(dat)[column]
    width <- case_when(
      grepl("definition|description", header, ignore.case = TRUE) ~ 60,
      grepl("box_label|label$", header, ignore.case = TRUE) ~ 46,
      grepl("box_id", header, ignore.case = TRUE) ~ 40,
      grepl("model_[12]$", header, ignore.case = TRUE) ~ 28,
      grepl("sample|timepoint_info", header, ignore.case = TRUE) ~ 20,
      grepl("source|note", header, ignore.case = TRUE) ~ 28,
      TRUE ~ 16
    )
    setColWidths(wb, sheet, cols = column, widths = width)
  }
}

add_readme <- function(wb, workbook_type, manifest_rows) {
  addWorksheet(wb, "README", gridLines = FALSE)
  title <- if (workbook_type == "main") "Source Data - Main Figures" else "Source Data - Extended Data Figures"
  status_counts <- table(manifest_rows$audit_status)
  values <- data.frame(
    field = c(title, "Purpose", "Current authority", "Rebuild", "Panel sheets", "Audit status", "Provenance", "Important"),
    value = c(
      "",
      "One sheet per locked figure panel, or a shared panel group when the locked figure uses one common dataset.",
      "Scripts_2025/Final_Scripts/final_manuscript_objects/ and the locked Revision 1 figure PDF.",
      "Run Scripts_2025/Final_Scripts/5_3_Rebuild_Figure_Source_Data_Workbooks.R from the project root.",
      as.character(nrow(manifest_rows)),
      paste(names(status_counts), as.integer(status_counts), sep = ": ", collapse = "; "),
      "Each panel sheet retains figure-panel and source-table provenance fields where distributable.",
      "Test-cohort panels use the populations shown in the current locked figures; assay-specific denominators are retained in each panel sheet."
    ),
    stringsAsFactors = FALSE
  )
  writeData(wb, "README", values, colNames = FALSE)
  mergeCells(wb, "README", cols = 1:2, rows = 1)
  addStyle(wb, "README", title_style, rows = 1, cols = 1:2, gridExpand = TRUE)
  addStyle(wb, "README", label_style, rows = 2:nrow(values), cols = 1, gridExpand = TRUE)
  setColWidths(wb, "README", cols = 1, widths = 22)
  setColWidths(wb, "README", cols = 2, widths = 95)
  setRowHeights(wb, "README", rows = seq_len(nrow(values)), heights = "auto")
  freezePane(wb, "README", firstRow = TRUE)
  pageSetup(wb, "README", orientation = "landscape", scale = 90, paperSize = 9)
}

build_workbook <- function(workbook_type, manifest_rows, audit_rows, output_dir, manifest_dir) {
  wb <- createWorkbook(creator = "Figure source-data rebuild workflow")

  for (i in seq_len(nrow(manifest_rows))) {
    item <- manifest_rows[i, ]
    sheet <- item$sheet_name[[1]]
    if (!nzchar(sheet) || nchar(sheet) > 31L) stop("Invalid worksheet name: ", sheet, call. = FALSE)
    panel_path <- item$panel_csv[[1]]
    if (!grepl("^/", panel_path)) panel_path <- file.path(manifest_dir, panel_path)
    dat <- read_panel_for_excel(panel_path)
    addWorksheet(wb, sheet, gridLines = FALSE)
    writeData(wb, sheet, dat, keepNA = FALSE)
    style_data_sheet(wb, sheet, dat)
  }

  expected <- nrow(manifest_rows)
  if (length(names(wb)) != expected) stop("Workbook sheet-count validation failed.", call. = FALSE)
  filename <- if (workbook_type == "main") "Source_Data_Main_Figures.xlsx" else "Source_Data_Extended_Data_Figures.xlsx"
  output <- file.path(output_dir, filename)
  saveWorkbook(wb, output, overwrite = TRUE)
  sanitize_unused_drawing_relationships(output)
  output
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
manifest_path <- normalizePath(args$manifest, mustWork = TRUE)
audit_path <- normalizePath(args$audit, mustWork = TRUE)
output_dir <- args$`output-dir`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read_csv(manifest_path, show_col_types = FALSE)
audit <- read_csv(audit_path, show_col_types = FALSE)
expected_manifest_rows <- 74L # 18 main (including Fig. 2E) + 56 extended.
if (nrow(manifest) != expected_manifest_rows || anyDuplicated(manifest$sheet_name)) {
  stop(
    "Manifest must contain ", expected_manifest_rows,
    " unique locked panel rows (18 main, including Fig. 2E, + 56 extended).",
    call. = FALSE
  )
}

outputs <- vapply(c("main", "extended"), function(workbook_type) {
  build_workbook(
    workbook_type,
    manifest %>% filter(workbook == workbook_type),
    audit %>% filter(workbook == workbook_type),
    output_dir,
    dirname(manifest_path)
  )
}, character(1))

# Keep the canonical manuscript-object tree synchronized on every rebuild.
# The panel-specific ED2 CSVs/PNGs are exported by their native generators;
# this full workbook belongs at the extended-data category root because it
# contains source sheets for all extended-data figures, not only ED Figure 2.
final_objects_root <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects"
)
extended_source_destinations <- c(
  file.path(
    final_objects_root,
    "02_extended_data_figures",
    "Source_Data_Extended_Data_Figures.xlsx"
  ),
  file.path(
    final_objects_root,
    "06_helpers_and_source_data",
    "Source_Data_Extended_Data_Figures.xlsx"
  )
)

for (destination in extended_source_destinations) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(outputs[["extended"]], destination, overwrite = TRUE)
  if (!copied) {
    stop(
      "Failed to copy the extended-data source workbook to final manuscript objects: ",
      destination,
      call. = FALSE
    )
  }
}

cat(paste(normalizePath(outputs, mustWork = TRUE), collapse = "\n"), "\n", sep = "")
cat(
  paste(normalizePath(extended_source_destinations, mustWork = TRUE), collapse = "\n"),
  "\n",
  sep = ""
)
