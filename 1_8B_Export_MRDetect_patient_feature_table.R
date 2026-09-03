# =============================================================================
# 1_8B_Export_MRDetect_patient_feature_table.R
# Project: cfWGS MRDetect
#
# Purpose:
#   Create a non-destructive, analysis-ready table of matched-patient MRDetect
#   values using both tumor/BM-derived mutations and cfDNA-derived mutations as
#   the patient-specific baseline/diagnosis mutation source.
#
# Unit of analysis:
#   One row is one queried cfDNA sample evaluated against one matched-patient
#   baseline/diagnosis mutation-list VCF and one mutation-source compartment.
#
# Interpretation boundary:
#   This is a retrospective feature export. Mutation-source and queried-sample
#   dates are retained so temporal ordering can be audited; it is not itself a
#   prospective-validation table.
#
# Input:
#   MRDetect_output_winter_2025/Processed_R_outputs/
#     cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.rds
#   Output_tables_2025/clinical_support/cfwgs_sample_identity_map.csv
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
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

input_file <- file.path(
  "MRDetect_output_winter_2025",
  "Processed_R_outputs",
  "cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.rds"
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

sample_identity_path <- file.path(
  "Output_tables_2025", "clinical_support", "cfwgs_sample_identity_map.csv"
)
if (!file.exists(sample_identity_path)) {
  stop("cfWGS sample identity map not found: ", sample_identity_path, call. = FALSE)
}

clinical_sample_dates <- readr::read_csv(
  sample_identity_path,
  show_col_types = FALSE
) %>%
  transmute(
    Patient = as.character(Patient),
    sample_id_tested = as.character(Sample_ID),
    tested_sample_date = as.Date(Date_of_sample_collection)
  ) %>%
  filter(
    !is.na(Patient),
    !is.na(sample_id_tested),
    nzchar(sample_id_tested),
    !is.na(tested_sample_date)
  ) %>%
  distinct()

duplicate_clinical_sample_dates <- clinical_sample_dates %>%
  count(Patient, sample_id_tested, name = "n") %>%
  filter(n > 1)
if (nrow(duplicate_clinical_sample_dates) > 0L) {
  stop(
    "Clinical metadata has duplicate Patient/Sample_ID date mappings; ",
    "cannot enforce mutation-source temporal ordering.",
    call. = FALSE
  )
}

required_cols <- c(
  "Study", "Patient", "Patient_Bam", "Sample_ID", "Sample_ID_Bam", "BAM", "VCF",
  "VCF_clean", "Sample_type", "Sample_type_Bam", "timepoint_info",
  "timepoint_info_Bam", "Mut_source", "Filter_source", "plotting_type",
  "Date_of_sample_collection",
  "sites_checked", "reads_checked", "sites_detected", "reads_detected",
  "total_reads", "detection_rate",
  "detection_rate_as_reads_detected_over_reads_checked",
  "detection_rate_as_reads_detected_over_total_reads",
  "sites_detection_rate", "detection_rate_zscore_charm",
  "detection_rate_zscore_reads_checked_charm",
  "detection_rate_zscore_total_reads_charm", "sites_rate_zscore_charm",
  "mean_det_charm", "sd_det_charm", "mean_det_checked_charm",
  "sd_det_checked_charm", "mean_det_total_charm", "sd_det_total_charm",
  "mean_sites_charm", "sd_sites_charm", "Relapsed",
  "Num_days_to_closest_relapse"
)

missing_cols <- setdiff(required_cols, names(mrdetect))
if (length(missing_cols) > 0) {
  stop(
    "Input MRDetect table is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

patient_feature_candidates <- mrdetect %>%
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
    mutation_source_date = as.Date(Date_of_sample_collection),
    sample_id_tested = as.character(Sample_ID_Bam)
  ) %>%
  left_join(
    clinical_sample_dates,
    by = c("Patient", "sample_id_tested")
  )

# The table included with the manuscript source data is limited to mutations
# measured at baseline/diagnosis. Preserve excluded later-source rows in an audit rather than silently
# including maintenance, relapse, or progression mutation lists.
nonbaseline_mutation_source_rows <- patient_feature_candidates %>%
  filter(is.na(.data$timepoint_info) |
           !.data$timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  transmute(
    Study,
    Patient,
    mutation_source_sample_id = Sample_ID,
    mutation_source_timepoint = timepoint_info,
    mutation_source_date,
    sample_id_tested,
    timepoint_tested = timepoint_info_Bam,
    tested_sample_date,
    Mut_source,
    VCF,
    exclusion_reason = "mutation source is not labelled Baseline or Diagnosis"
  ) %>%
  distinct() %>%
  arrange(Patient, mutation_source_date, tested_sample_date)

write_csv(
  nonbaseline_mutation_source_rows,
  file.path(output_dir, "MRDetect_excluded_nonbaseline_mutation_sources.csv"),
  na = ""
)

patient_feature_candidates <- patient_feature_candidates %>%
  filter(.data$timepoint_info %in% c("Baseline", "Diagnosis"))

temporally_invalid_pairs <- patient_feature_candidates %>%
  filter(
    Patient == "SPORE_0012",
    Sample_ID == "SPORE_0012_T4_BM_cells",
    !is.na(mutation_source_date),
    !is.na(tested_sample_date),
    tested_sample_date < mutation_source_date
  ) %>%
  transmute(
    Patient,
    mutation_source_sample_id = Sample_ID,
    mutation_source_date,
    sample_id_tested,
    tested_sample_date,
    Mut_source,
    exclusion_reason = "tested plasma predates mutation-source sample"
  ) %>%
  distinct() %>%
  arrange(Patient, mutation_source_date, tested_sample_date)

write_csv(
  temporally_invalid_pairs,
  file.path(output_dir, "MRDetect_temporally_invalid_source_query_pairs.csv"),
  na = ""
)

patient_features <- patient_feature_candidates %>%
  filter(
    !(
      Patient == "SPORE_0012" &
        Sample_ID == "SPORE_0012_T4_BM_cells" &
        !is.na(mutation_source_date) &
        !is.na(tested_sample_date) &
        tested_sample_date < mutation_source_date
    )
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
    tested_sample_date,
    timepoint_tested = timepoint_info_Bam,
    mutation_source_sample_id = Sample_ID,
    mutation_source_type = Sample_type,
    mutation_source_date,
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
