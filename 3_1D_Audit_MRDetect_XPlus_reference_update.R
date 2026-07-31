#!/usr/bin/env Rscript

# Compare a preserved pre-correction scored table with the current frozen-model
# scored table. This script never fits models or changes thresholds.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript 3_1D_Audit_MRDetect_XPlus_reference_update.R ",
    "<pre_correction.rds> <corrected.rds> <output_directory>",
    call. = FALSE
  )
}

old_path <- args[[1]]
new_path <- args[[2]]
output_dir <- args[[3]]
if (!file.exists(old_path)) stop("Missing pre-correction RDS: ", old_path, call. = FALSE)
if (!file.exists(new_path)) stop("Missing corrected RDS: ", new_path, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

old <- readRDS(old_path)
new <- readRDS(new_path)
key_cols <- c("Patient", "Sample_Code", "Date", "Timepoint")
missing_keys <- setdiff(key_cols, intersect(names(old), names(new)))
if (length(missing_keys) > 0L) {
  stop("Missing comparison key(s): ", paste(missing_keys, collapse = ", "), call. = FALSE)
}
if (nrow(old) != nrow(new) || !identical(old[key_cols], new[key_cols])) {
  stop("Pre-correction and corrected tables do not have identical ordered sample keys.", call. = FALSE)
}

value_changed <- function(x, y) {
  xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
}

feature_cols <- intersect(
  c(
    "zscore_BM", "zscore_blood",
    "z_score_detection_rate_BM", "z_score_detection_rate_blood",
    "detect_rate_BM", "detect_rate_blood"
  ),
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
names(old_values) <- paste0(audit_cols, "_old")
names(new_values) <- paste0(audit_cols, "_corrected")

zscore_cols <- intersect(
  feature_cols,
  c(
    "zscore_BM", "zscore_blood",
    "z_score_detection_rate_BM", "z_score_detection_rate_blood"
  )
)
any_zscore_change <- Reduce(
  `|`,
  lapply(zscore_cols, function(column) {
    value_changed(old[[column]][affected], new[[column]][affected])
  })
)
any_call_flip <- Reduce(
  `|`,
  lapply(changed_call_cols, function(column) {
    value_changed(old[[column]][affected], new[[column]][affected])
  })
)

audit <- bind_cols(
  new[affected, key_cols, drop = FALSE],
  tibble(
    cohort = if ("Cohort" %in% names(new)) new$Cohort[affected] else NA_character_,
    img159_primary_measurement = if_else(
      new$Patient[affected] == "IMG-159",
      "Original NovaSeq 6000 plasma measurement retained; revision repeat not averaged",
      NA_character_
    )
  ),
  old_values,
  new_values
) %>%
  mutate(
    any_zscore_change = any_zscore_change,
    any_call_flip = any_call_flip
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

write_csv(audit, file.path(output_dir, "mrdetect_xplus_reference_update_sample_impact.csv"), na = "")
write_csv(call_summary, file.path(output_dir, "mrdetect_xplus_reference_update_call_flip_summary.csv"), na = "")

message("Affected scored rows: ", nrow(audit))
message("Call columns with changes: ", nrow(call_summary))
