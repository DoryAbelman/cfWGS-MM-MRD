#!/usr/bin/env Rscript

###############################################################################
# Rebuild the main-figure and Extended Data source-data workbooks.
#
# Pipeline role
#   Runner for the figure-source-data workflow. It calls:
#     5_1  export and de-identify the 74 panel CSVs
#     5_2  assemble the main- and extended-data Excel workbooks
#     5_4  validate both workbooks and render QA previews
#   Scientific values come from the retained panel/model/source objects selected
#   by 5_1; this runner does not fit models or change thresholds.
#
# Run from the project root:
#   Rscript Scripts_2025/Final_Scripts/5_3_Rebuild_Figure_Source_Data_Workbooks.R
#
# Optional:
#   FIGURE_SOURCE_OUTPUT_DIR  Destination for the two workbooks and QA PDFs.
#
# Write scope
#   In addition to FIGURE_SOURCE_OUTPUT_DIR, 5_1 refreshes
#   Output_tables_2025/Figure_Source_Data and two panel source-data copies under
#   final_manuscript_objects; 5_2 refreshes the staged extended-data workbook
#   copies under final_manuscript_objects.
#   This is not a dry-run command. Existing outputs at those paths may be
#   overwritten by the three scripts it calls.
###############################################################################

project_root <- normalizePath(".", mustWork = TRUE)
script_root <- file.path(project_root, "Scripts_2025", "Final_Scripts")
default_output <- file.path(project_root, "outputs", "figure_source_data")
output_dir <- Sys.getenv("FIGURE_SOURCE_OUTPUT_DIR", unset = default_output)

exporter <- file.path(script_root, "5_1_Export_Locked_Figure_Source_Data.R")
builder <- file.path(script_root, "5_2_Build_Figure_Source_Data_Workbooks.R")
validator <- file.path(script_root, "5_4_Validate_Figure_Source_Data_Workbooks.R")
required <- c(exporter, builder, validator)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Rebuild dependency is missing:\n- ", paste(missing, collapse = "\n- "), call. = FALSE)

run_rscript <- function(script, args = character()) {
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(args)))
  if (status != 0L) stop(basename(script), " failed with status ", status, call. = FALSE)
}

run_rscript(exporter)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- file.path(project_root, "Output_tables_2025", "Figure_Source_Data", "source_data_workbook_manifest.csv")
audit <- file.path(project_root, "Output_tables_2025", "Figure_Source_Data", "source_data_audit.csv")
run_rscript(builder, c("--manifest", manifest, "--audit", audit, "--output-dir", output_dir))

workbooks <- c(
  main = file.path(output_dir, "Source_Data_Main_Figures.xlsx"),
  extended = file.path(output_dir, "Source_Data_Extended_Data_Figures.xlsx")
)
for (workbook_type in names(workbooks)) {
  run_rscript(
    validator,
    c(
      "--xlsx", workbooks[[workbook_type]],
      "--preview-dir", file.path(output_dir, "source_data_workbook_previews", workbook_type)
    )
  )
}

message("Rebuild complete: ", normalizePath(output_dir, mustWork = TRUE))
