
# ==============================================================================
# Script: 1_10_Estimate_MRDetect_LOD_proxy.R
# ==============================================================================
# 1_10_Estimate_MRDetect_LOD_proxy.R
#
# ## Goal
# Add transparent, read-denominator-based limit-of-detection (LOD) proxy columns
# to an MRDetect detection-rate table. The output is intended for QC,
# interpretability, and manuscript-supporting sensitivity analyses; it does not
# re-run MRDetect and it does not change detection calls.
#
# ## Scientific rationale
# MRDetect tracks a fixed set of patient-specific mutation sites. For a sample
# with `N` informative reads checked at those sites, the smallest observable
# non-zero read fraction is approximately 1/N. A 95% Poisson one-read detection
# proxy is -log(0.05)/N, the expected mutant fraction needed so that the
# probability of observing at least one mutant read is 95% under an idealized
# Poisson sampling model.
#
# This is a theoretical read-count proxy, not a validated clinical LOD. It does
# not model sequencing error, variant-specific mappability, molecule collapsing,
# strand bias, site-level dropout, tumor fraction heterogeneity, or MRDetect's
# full statistical decision rule. Interpret it as an auditable denominator-based
# sensitivity descriptor.
#
# ## How to run
# Run from the project root with defaults:
#   Rscript Scripts_2025/Final_Scripts/1_10_Estimate_MRDetect_LOD_proxy.R
#
# Optional arguments:
#   Rscript Scripts_2025/Final_Scripts/1_10_Estimate_MRDetect_LOD_proxy.R \
#     --input path/to/detection_rates.csv \
#     --output path/to/detection_rates_with_lod.csv \
#     --ir-column reads_checked
#
# ## Role in manuscript workflow
# This is an upstream QC/annotation script. It augments MRDetect detection-rate
# outputs with denominator-derived LOD proxy fields that can be used to audit
# analytical sensitivity and explain why low-informative-read samples have
# weaker detection limits.
#
# ## Main input
# - MRDetect_output_winter_2025/Processed_R_outputs/
#   All_detection_rates_baseline_and_controls_Feb2026.csv
#
# Required column:
# - reads_checked, unless another informative-read denominator is supplied with
#   --ir-column.
#
# Optional columns used when present:
# - sites_checked
# - reads_detected
# - total_reads
#
# ## Main output
# - MRDetect_output_winter_2025/Processed_R_outputs/
#   All_detection_rates_baseline_and_controls_Feb2026_with_mrdetect_lod.csv
#
# ## Added columns
# - informative_reads_mrdetect
# - informative_reads_source
# - lod_fraction_one_read and lod_ppm_one_read
# - lod95_fraction_one_read and lod95_ppm_one_read
# - poisson_reads_required_for_lod95_one_read
# - mean_reads_per_checked_site
# - observed_mutant_fraction_over_informative_reads
# - observed_mutant_fraction_over_total_reads
# - lod_interpretation
# - denominator_note
#
# ## Assumptions and audit notes
# - Informative-read denominators must be positive numeric values; non-numeric or
#   non-positive values produce NA LOD fields rather than silent zeroes.
# - CSV inputs are written back as CSV; other extensions are treated as tabular
#   delimited text.
# - The default denominator is `reads_checked`, matching the MRDetect output
#   convention used in this project.
# - This script is deterministic and has no random component.
# ==============================================================================


args <- commandArgs(trailingOnly = TRUE)

defaults <- list(
  input = "MRDetect_output_winter_2025/Processed_R_outputs/All_detection_rates_baseline_and_controls_Feb2026.csv",
  output = "MRDetect_output_winter_2025/Processed_R_outputs/All_detection_rates_baseline_and_controls_Feb2026_with_mrdetect_lod.csv",
  ir_column = "reads_checked"
)

parse_args <- function(args, defaults) {
  options <- defaults

  if (length(args) == 0) {
    return(options)
  }

  if (length(args) %% 2 != 0) {
    stop("Arguments must be supplied as --key value pairs.")
  }

  key_map <- c(
    input = "input",
    output = "output",
    "ir-column" = "ir_column",
    ir_column = "ir_column"
  )

  for (index in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[index]])
    value <- args[[index + 1]]

    if (!key %in% names(key_map)) {
      stop(sprintf("Unknown argument '%s'. Expected --input, --output, or --ir-column.", args[[index]]))
    }

    options[[key_map[[key]]]] <- value
  }

  options
}

read_table_auto <- function(path) {
  extension <- tolower(tools::file_ext(path))

  if (extension == "csv") {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }

  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
}

strip_wrapping_quotes <- function(x) {
  if (!is.character(x)) {
    return(x)
  }

  sub('^"(.*)"$', '\\1', x)
}

write_table_auto <- function(data, path) {
  extension <- tolower(tools::file_ext(path))

  if (extension == "csv") {
    write.csv(data, path, row.names = FALSE, quote = TRUE)
    return(invisible(NULL))
  }

  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

as_numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

options <- parse_args(args, defaults)

if (!file.exists(options$input)) {
  stop(sprintf("Input file not found: %s", options$input))
}

df <- read_table_auto(options$input)
names(df) <- strip_wrapping_quotes(names(df))
df[] <- lapply(df, strip_wrapping_quotes)

if (!options$ir_column %in% names(df)) {
  stop(sprintf("Column '%s' was not found in %s", options$ir_column, options$input))
}

informative_reads <- as_numeric_or_na(df[[options$ir_column]])
sites_checked <- if ("sites_checked" %in% names(df)) as_numeric_or_na(df[["sites_checked"]]) else rep(NA_real_, nrow(df))
reads_detected <- if ("reads_detected" %in% names(df)) as_numeric_or_na(df[["reads_detected"]]) else rep(NA_real_, nrow(df))
total_reads <- if ("total_reads" %in% names(df)) as_numeric_or_na(df[["total_reads"]]) else rep(NA_real_, nrow(df))

lambda_95_one_read <- -log(0.05)

df$informative_reads_mrdetect <- informative_reads
df$informative_reads_source <- options$ir_column
df$lod_fraction_one_read <- ifelse(!is.na(informative_reads) & informative_reads > 0, 1 / informative_reads, NA_real_)
df$lod_ppm_one_read <- df$lod_fraction_one_read * 1e6
df$lod95_fraction_one_read <- ifelse(!is.na(informative_reads) & informative_reads > 0, lambda_95_one_read / informative_reads, NA_real_)
df$lod95_ppm_one_read <- df$lod95_fraction_one_read * 1e6
df$poisson_reads_required_for_lod95_one_read <- lambda_95_one_read
df$mean_reads_per_checked_site <- ifelse(!is.na(sites_checked) & sites_checked > 0, informative_reads / sites_checked, NA_real_)
df$observed_mutant_fraction_over_informative_reads <- ifelse(!is.na(reads_detected) & informative_reads > 0, reads_detected / informative_reads, NA_real_)
df$observed_mutant_fraction_over_total_reads <- ifelse(!is.na(reads_detected) & !is.na(total_reads) & total_reads > 0, reads_detected / total_reads, NA_real_)
df$lod_interpretation <- if (identical(options$ir_column, "reads_checked")) {
  paste(
    "Theoretical Poisson LOD using MRDetect reads_checked as the informative-read denominator.",
    "This is exact under the MRDetect read-based counting scheme when reads_checked is defined as the number of reads checked at tracked mutation sites."
  )
} else {
  sprintf("Theoretical Poisson LOD using %s as the informative-read denominator", options$ir_column)
}
df$denominator_note <- if (identical(options$ir_column, "reads_checked")) {
  paste(
    "reads_checked = number of reads checked at sites overlapping tracked mutations;",
    "reads_detected = reads carrying the mutation;",
    "sites_checked = tracked variant sites checked;",
    "sites_detected = tracked mutant sites with detection."
  )
} else {
  sprintf("Denominator supplied from column %s", options$ir_column)
}

dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
write_table_auto(df, options$output)

cat(sprintf("Wrote %d rows to %s\n", nrow(df), options$output))
cat(sprintf("Denominator column: %s\n", options$ir_column))

if (identical(options$ir_column, "reads_checked")) {
  cat("Interpretation: read-based theoretical LOD using MRDetect reads_checked as the informative-read denominator.\n")
}
