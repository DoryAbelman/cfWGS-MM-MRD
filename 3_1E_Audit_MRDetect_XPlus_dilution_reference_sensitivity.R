#!/usr/bin/env Rscript

# Compare dilution-series scoring produced with two XPlus healthy-control
# reference panels. This script does not fit models or alter thresholds.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript 3_1E_Audit_MRDetect_XPlus_dilution_reference_sensitivity.R ",
    "<paired19.csv> <primary22.csv> <output_directory>",
    call. = FALSE
  )
}

old_path <- args[[1]]
new_path <- args[[2]]
output_dir <- args[[3]]
if (!file.exists(old_path)) stop("Missing paired-19 dilution table: ", old_path, call. = FALSE)
if (!file.exists(new_path)) stop("Missing primary-22 dilution table: ", new_path, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

old <- read_csv(old_path, show_col_types = FALSE)
new <- read_csv(new_path, show_col_types = FALSE)
key_cols <- c("Patient", "Sample_ID", "Sample", "LOD")
missing_keys <- setdiff(key_cols, intersect(names(old), names(new)))
if (length(missing_keys)) {
  stop("Missing comparison key(s): ", paste(missing_keys, collapse = ", "), call. = FALSE)
}

old <- old %>% arrange(across(all_of(key_cols)))
new <- new %>% arrange(across(all_of(key_cols)))
if (nrow(old) != nrow(new) || !identical(old[key_cols], new[key_cols])) {
  stop("Paired-19 and primary-22 dilution tables have different sample keys.", call. = FALSE)
}

value_changed <- function(x, y) {
  xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
}

feature_cols <- intersect(
  c("zscore_BM", "zscore_blood", "z_score_detection_rate_BM",
    "z_score_detection_rate_blood", "detect_rate_BM", "detect_rate_blood"),
  intersect(names(old), names(new))
)
call_cols <- intersect(grep("_call$", names(new), value = TRUE), names(old))
changed_call_cols <- call_cols[vapply(
  call_cols,
  function(column) any(value_changed(old[[column]], new[[column]])),
  logical(1)
)]
prob_cols <- sub("_call$", "_prob", changed_call_cols)
prob_cols <- prob_cols[prob_cols %in% intersect(names(old), names(new))]
audit_cols <- unique(c(feature_cols, prob_cols, changed_call_cols))

changed_matrix <- vapply(
  audit_cols,
  function(column) value_changed(old[[column]], new[[column]]),
  logical(nrow(new))
)
affected <- if (length(audit_cols) == 1L) changed_matrix else apply(changed_matrix, 1L, any)

old_values <- as_tibble(old[affected, audit_cols, drop = FALSE])
new_values <- as_tibble(new[affected, audit_cols, drop = FALSE])
names(old_values) <- paste0(audit_cols, "_paired19")
names(new_values) <- paste0(audit_cols, "_primary22")

impact <- bind_cols(
  new[affected, key_cols, drop = FALSE],
  old_values,
  new_values
)

call_summary <- bind_rows(lapply(changed_call_cols, function(column) {
  tibble(
    call_column = column,
    n_flips = sum(value_changed(old[[column]], new[[column]])),
    n_zero_to_one = sum(old[[column]] == 0 & new[[column]] == 1, na.rm = TRUE),
    n_one_to_zero = sum(old[[column]] == 1 & new[[column]] == 0, na.rm = TRUE),
    n_missingness_changes = sum(xor(is.na(old[[column]]), is.na(new[[column]])))
  )
}))

correlation_summary <- bind_rows(lapply(feature_cols, function(column) {
  paired19_complete <- complete.cases(old[[column]], old$LOD)
  primary22_complete <- complete.cases(new[[column]], new$LOD)
  tibble(
    feature = column,
    rho_paired19 = cor(old[[column]][paired19_complete], old$LOD[paired19_complete], method = "spearman"),
    rho_primary22 = cor(new[[column]][primary22_complete], new$LOD[primary22_complete], method = "spearman"),
    n_paired19 = sum(paired19_complete),
    n_primary22 = sum(primary22_complete)
  )
}))

write_csv(impact, file.path(output_dir, "dilution_reference_22_vs_19_sample_impact.csv"), na = "")
write_csv(call_summary, file.path(output_dir, "dilution_reference_22_vs_19_call_flip_summary.csv"), na = "")
write_csv(correlation_summary, file.path(output_dir, "dilution_reference_22_vs_19_correlation_summary.csv"), na = "")

message("Affected dilution rows: ", nrow(impact))
message("Call columns with changes: ", nrow(call_summary))
