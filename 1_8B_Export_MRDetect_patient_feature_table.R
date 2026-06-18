# =============================================================================
# Export combined MRDetect patient feature table
# Project: cfWGS MRDetect
#
# Purpose:
#   Create a non-destructive, analysis-ready table of matched-patient MRDetect
#   values using both tumor/BM-derived mutations and cfDNA-derived mutations as
#   the patient-specific baseline mutation source.
#
# Input:
#   MRDetect_output_winter_2025/Processed_R_outputs/
#     cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025_2.rds
#
# Output:
#   MRDetect_output_winter_2025/Processed_R_outputs/Derived_exports/
#     MRDetect_patient_features_BM_and_cfDNA_baselines_Feb2026.csv
#     MRDetect_patient_features_BM_and_cfDNA_baselines_Feb2026.rds
#
# Notes:
#   This script only reads existing processed MRDetect output from script 1_8
#   and writes a new derived export. It does not modify any original outputs.
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_8B_Export_MRDetect_patient_feature_table.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script exports the combined MRDetect patient
#     feature table used by downstream model and concordance scripts.
#
# =============================================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

input_file <- file.path(
  "MRDetect_output_winter_2025",
  "Processed_R_outputs",
  "cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025_2.rds"
)

output_dir <- file.path(
  "MRDetect_output_winter_2025",
  "Processed_R_outputs",
  "Derived_exports"
)

output_csv <- file.path(
  output_dir,
  "MRDetect_patient_features_BM_and_cfDNA_baselines_Feb2026.csv"
)

output_rds <- file.path(
  output_dir,
  "MRDetect_patient_features_BM_and_cfDNA_baselines_Feb2026.rds"
)

if (!file.exists(input_file)) {
  stop("Input MRDetect RDS not found: ", input_file, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mrdetect <- readRDS(input_file)

required_cols <- c(
  "Patient", "Patient_Bam", "Sample_ID", "Sample_ID_Bam", "BAM", "VCF",
  "VCF_clean", "Sample_type", "Sample_type_Bam", "timepoint_info",
  "timepoint_info_Bam", "Mut_source", "Filter_source", "plotting_type",
  "sites_checked", "reads_checked", "sites_detected", "reads_detected",
  "total_reads", "detection_rate",
  "detection_rate_as_reads_detected_over_reads_checked",
  "detection_rate_as_reads_detected_over_total_reads",
  "sites_detection_rate", "detection_rate_zscore_charm",
  "detection_rate_zscore_reads_checked_charm",
  "detection_rate_zscore_total_reads_charm", "sites_rate_zscore_charm"
)

missing_cols <- setdiff(required_cols, names(mrdetect))
if (length(missing_cols) > 0) {
  stop(
    "Input MRDetect table is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

patient_features <- mrdetect %>%
  filter(
    plotting_type == "Matched_plasma",
    Filter_source == "STR_encode",
    Mut_source %in% c("BM_cells", "Blood"),
    Sample_type_Bam == "Blood_plasma_cfDNA",
    !is.na(Patient),
    !is.na(Patient_Bam),
    Patient == Patient_Bam
  ) %>%
  mutate(
    baseline_source = case_when(
      Mut_source == "BM_cells" ~ "Tumor_BM_baseline",
      Mut_source == "Blood" ~ "cfDNA_baseline",
      TRUE ~ NA_character_
    ),
    mrdetect_positive_sites_z4_5 = sites_rate_zscore_charm > 4.5,
    cumulative_vaf_if_positive = if_else(
      mrdetect_positive_sites_z4_5,
      detection_rate_as_reads_detected_over_reads_checked,
      0
    )
  ) %>%
  group_by(
    baseline_source,
    Patient,
    Sample_ID_Bam,
    Sample_ID,
    VCF
  ) %>%
  # Match script 1_8's duplicate-run handling: keep the deepest run.
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    baseline_source,
    Study,
    Patient,
    sample_id_tested = Sample_ID_Bam,
    bam_tested = BAM,
    sample_type_tested = Sample_type_Bam,
    timepoint_tested = timepoint_info_Bam,
    mutation_source_sample_id = Sample_ID,
    mutation_source_type = Sample_type,
    mutation_source_timepoint = timepoint_info,
    mutation_source_vcf = VCF,
    mutation_source_vcf_clean = VCF_clean,
    Mut_source,
    Filter_source,
    sites_checked,
    reads_checked,
    sites_detected,
    reads_detected,
    total_reads,
    detection_rate,
    detection_rate_as_reads_detected_over_reads_checked,
    detection_rate_as_reads_detected_over_total_reads,
    sites_detection_rate,
    detection_rate_zscore_charm,
    detection_rate_zscore_reads_checked_charm,
    detection_rate_zscore_total_reads_charm,
    sites_rate_zscore_charm,
    mrdetect_positive_sites_z4_5,
    cumulative_vaf_if_positive,
    mean_det_charm,
    sd_det_charm,
    mean_det_checked_charm,
    sd_det_checked_charm,
    mean_det_total_charm,
    sd_det_total_charm,
    mean_sites_charm,
    sd_sites_charm,
    Relapsed,
    Num_days_to_closest_relapse
  ) %>%
  arrange(Patient, baseline_source, sample_id_tested, mutation_source_sample_id)

if (nrow(patient_features) == 0) {
  stop("No matched-patient MRDetect rows passed the export filters.", call. = FALSE)
}

duplicate_keys <- patient_features %>%
  count(
    baseline_source,
    Patient,
    sample_id_tested,
    mutation_source_sample_id,
    mutation_source_vcf,
    name = "n"
  ) %>%
  filter(n > 1)

if (nrow(duplicate_keys) > 0) {
  stop(
    "Duplicate MRDetect rows detected after filtering; inspect before export. ",
    "Number of duplicate keys: ", nrow(duplicate_keys),
    call. = FALSE
  )
}

write_csv(patient_features, output_csv, na = "")
saveRDS(patient_features, output_rds)

message("Wrote: ", output_csv)
message("Wrote: ", output_rds)
message("Rows: ", nrow(patient_features))
message("Patients: ", n_distinct(patient_features$Patient))
message("Rows by baseline source:")
print(patient_features %>% count(baseline_source), n = Inf)
