#!/usr/bin/env Rscript

###############################################################################
# 5_4_Validate_Figure_Source_Data_Workbooks.R
#
# Purpose and pipeline role
#   Validate the Excel workbooks built by 5_2. This is the final step called by
#   5_3_Rebuild_Figure_Source_Data_Workbooks.R and does not modify the workbook
#   or any scientific value.
#
# Validation
#   - Requires 18 main-figure or 56 extended-data sheets, with no extra README
#     or AUDIT sheets. This script checks the counts, not the complete expected
#     sheet-name list or the per-panel schema contract written by 5_1.
#   - Stops on spreadsheet formula errors, original patient identifiers,
#     calendar dates, administrative provenance columns, empty sheets, or
#     dangling worksheet drawing relationships.
#   - If LibreOffice is available on PATH, renders the complete workbook to PDF
#     for later visual inspection; it does not inspect the rendered pages.
#     Otherwise structural/content validation still runs.
#
# Inputs and outputs
#   - Input: --xlsx <Source_Data_*.xlsx>; id_map.rds must be available from the
#     project root for the identifier check.
#   - Output: console validation status and, when possible, a PDF under
#     --preview-dir. The input workbook is not rewritten, but the preview
#     directory is created and an existing preview PDF may be replaced.
#
# How to run
#   From the project root with no arguments, validate both workbooks under
#   outputs/figure_source_data. The 5_3 runner instead calls this script once
#   per workbook with explicit --xlsx and --preview-dir arguments.
###############################################################################

suppressPackageStartupMessages({
  library(openxlsx)
})

validate_no_dangling_worksheet_drawing_rels <- function(xlsx) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("OOXML relationship validation requires the xml2 package.", call. = FALSE)
  }
  unpack_dir <- tempfile("source_data_validate_")
  dir.create(unpack_dir)
  on.exit(unlink(unpack_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(xlsx, exdir = unpack_dir)
  rel_dir <- file.path(unpack_dir, "xl", "worksheets", "_rels")
  rel_files <- list.files(rel_dir, pattern = "[.]rels$", full.names = TRUE)
  dangling <- character()
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
        dangling <- c(dangling, paste0(basename(rel_file), " -> ", target))
      }
    }
  }
  if (length(dangling)) {
    stop("Dangling worksheet drawing relationship(s): ",
         paste(dangling, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

# When run directly from the project root, validate both standard source-data
# workbooks. The explicit --xlsx/--preview-dir interface remains available for
# callers such as 5_3_Rebuild_Figure_Source_Data_Workbooks.R.
cli_args <- commandArgs(trailingOnly = TRUE)
if (!length(cli_args)) {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
  output_root <- file.path(project_root, "outputs", "figure_source_data")
  workbook_specs <- list(
    main = list(
      xlsx = file.path(output_root, "Source_Data_Main_Figures.xlsx"),
      preview_dir = file.path(output_root, "source_data_workbook_previews", "main")
    ),
    extended = list(
      xlsx = file.path(output_root, "Source_Data_Extended_Data_Figures.xlsx"),
      preview_dir = file.path(output_root, "source_data_workbook_previews", "extended")
    )
  )
  script_path <- file.path(
    project_root, "Scripts_2025", "Final_Scripts",
    "5_4_Validate_Figure_Source_Data_Workbooks.R"
  )
  for (spec in workbook_specs) {
    status <- system2(
      file.path(R.home("bin"), "Rscript"),
      c(
        shQuote(script_path),
        "--xlsx", shQuote(spec$xlsx),
        "--preview-dir", shQuote(spec$preview_dir)
      )
    )
    if (!identical(status, 0L)) {
      stop("Workbook validation failed with status ", status, ".", call. = FALSE)
    }
  }
  quit(save = "no", status = 0L)
}

parse_args <- function(args) {
  if (length(args) %% 2L != 0L) stop("Arguments must be --key value pairs.", call. = FALSE)
  keys <- sub("^--", "", args[seq(1L, length(args), by = 2L)])
  values <- args[seq(2L, length(args), by = 2L)]
  out <- setNames(as.list(values), keys)
  required <- c("xlsx", "preview-dir")
  missing <- required[!required %in% names(out)]
  if (length(missing)) stop("Missing argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)
  out
}

args <- parse_args(cli_args)
xlsx <- normalizePath(args$xlsx, mustWork = TRUE)
preview_dir <- args$`preview-dir`
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

validate_no_dangling_worksheet_drawing_rels(xlsx)

sheet_names <- getSheetNames(xlsx)
expected_sheet_count <- if (basename(xlsx) == "Source_Data_Main_Figures.xlsx") {
  18L
} else if (basename(xlsx) == "Source_Data_Extended_Data_Figures.xlsx") {
  56L
} else {
  stop("Unrecognized source-data workbook name: ", basename(xlsx), call. = FALSE)
}
if (length(sheet_names) != expected_sheet_count || any(sheet_names %in% c("README", "AUDIT"))) {
  stop("Unexpected workbook sheet inventory: ", paste(sheet_names, collapse = ", "), call. = FALSE)
}

id_map_path <- file.path(normalizePath(getwd(), mustWork = TRUE), "id_map.rds")
if (!file.exists(id_map_path)) stop("Missing deidentification authority: id_map.rds", call. = FALSE)
id_map <- readRDS(id_map_path)
raw_ids <- as.character(id_map$Patient)
forbidden_headers <- c(
  "figure_panel", "locked_figure", "source_table", "source_note",
  "source_call_path", "record_type"
)

formula_errors <- character()
calendar_date_violations <- character()
identifier_violations <- character()
administrative_column_violations <- character()
for (sheet in sheet_names) {
  dat <- read.xlsx(xlsx, sheet = sheet, colNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  if (!nrow(dat) || !ncol(dat)) stop("Worksheet is empty: ", sheet, call. = FALSE)
  values <- unlist(dat, use.names = FALSE)
  hits <- grep("#REF!|#DIV/0!|#VALUE!|#NAME[?]|#N/A", as.character(values), value = TRUE)
  if (length(hits)) formula_errors <- c(formula_errors, paste(sheet, unique(hits), sep = ": "))
  headers <- as.character(dat[1, , drop = TRUE])
  forbidden_present <- intersect(headers, forbidden_headers)
  if (length(forbidden_present)) {
    administrative_column_violations <- c(
      administrative_column_violations,
      paste0(sheet, ": ", paste(forbidden_present, collapse = ", "))
    )
  }
  raw_present <- raw_ids[vapply(
    raw_ids,
    function(raw_id) any(grepl(
      paste0("(?<![A-Za-z0-9])", raw_id, "(?![A-Za-z0-9])"),
      as.character(values), perl = TRUE
    ), na.rm = TRUE),
    logical(1)
  )]
  if (length(raw_present)) {
    identifier_violations <- c(
      identifier_violations,
      paste0(sheet, ": ", paste(raw_present, collapse = ", "))
    )
  }
  date_headers <- headers[grepl("(^|[._ -])date($|[._ -])", headers, ignore.case = TRUE)]
  date_values <- grep(
    "\\b(?:19|20)[0-9]{2}[-/.][0-9]{1,2}[-/.][0-9]{1,2}\\b",
    as.character(values), value = TRUE, perl = TRUE
  )
  if (length(date_headers) || length(date_values)) {
    calendar_date_violations <- c(
      calendar_date_violations,
      paste0(
        sheet, ": headers=", paste(unique(date_headers), collapse = ","),
        "; values=", paste(utils::head(unique(date_values), 5L), collapse = ",")
      )
    )
  }
}
if (length(formula_errors)) {
  stop("Potential spreadsheet error(s): ", paste(formula_errors, collapse = "; "), call. = FALSE)
}
if (length(calendar_date_violations)) {
  stop(
    "Calendar dates are prohibited in publication source-data workbooks: ",
    paste(calendar_date_violations, collapse = "; "),
    call. = FALSE
  )
}
if (length(identifier_violations)) {
  stop(
    "Original patient identifiers are prohibited in publication source-data workbooks: ",
    paste(identifier_violations, collapse = "; "), call. = FALSE
  )
}
if (length(administrative_column_violations)) {
  stop(
    "Administrative provenance columns are prohibited in panel sheets: ",
    paste(administrative_column_violations, collapse = "; "), call. = FALSE
  )
}

soffice <- Sys.which("soffice")
rendered <- FALSE
render_failed <- FALSE
if (nzchar(soffice)) {
  office_profile <- tempfile("source_workbook_libreoffice_profile_")
  dir.create(office_profile, recursive = TRUE)
  on.exit(unlink(office_profile, recursive = TRUE, force = TRUE), add = TRUE)
  status <- system2(
    soffice,
    c(
      paste0("-env:UserInstallation=file://", normalizePath(office_profile, mustWork = TRUE)),
      "--headless", "--convert-to", "pdf", "--outdir", shQuote(preview_dir), shQuote(xlsx)
    ),
    stdout = TRUE, stderr = TRUE
  )
  expected_pdf <- file.path(preview_dir, paste0(tools::file_path_sans_ext(basename(xlsx)), ".pdf"))
  render_status <- attr(status, "status")
  if (is.null(render_status)) render_status <- 0L
  rendered <- identical(as.integer(render_status), 0L) && file.exists(expected_pdf)
  render_failed <- !rendered
  if (render_failed) {
    warning(
      "LibreOffice workbook rendering failed; structural/content validation still passed: ",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }
}

cat(
  xlsx, ": ", length(sheet_names), " sheets structurally validated",
  if (rendered) {
    " and rendered to PDF"
  } else if (render_failed) {
    " (LibreOffice rendering failed; visual rendering unverified)"
  } else {
    " (LibreOffice not available; rendering skipped)"
  },
  "\n", sep = ""
)
