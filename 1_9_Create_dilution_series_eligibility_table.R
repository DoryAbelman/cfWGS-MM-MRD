# =============================================================================
# Create_dilution_series_eligibility_table.R
# Project: cfWGS MRDetect (Winter 2025)
# Author: Dory Abelman
# Date: February 2026
#
# Purpose 
#   This script finds patients who have:
#     • One “tumor-high” timepoint (strong signal) AND
#     • One “tumor-low” timepoint (near-background signal)
#   so we can use those paired samples to build a dilution series.
#
# Key definitions 
#   • detection_rate_as_reads_detected_over_reads_checked is a decimal.
#       - 0.01  = 1%
#       - 1e-4 = 0.01%
#   • “tumor-high”  timepoint: detection_rate >= 0.05  (>=0.5%)
#   • “tumor-low”   timepoint: detection_rate <= 5e-4 (<=0.05%)
#
# Outputs:
#   1) eligible_dilution_pairs_Feb2026.csv
#        A table of patients with one high and one low timepoint.
#   2) dilution_plan_Feb2026.csv
#        For each eligible pair, a plan for 10^-1 to 10^-6 target fractions
#        with 3 technical replicates per level.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ──────────────────────────────────────────────────────────────────────────────
# Inputs / outputs
#   input_file = the BM MRDetect table that already includes healthy control means
#   output_dir = folder where we write eligibility + dilution plan tables
# ──────────────────────────────────────────────────────────────────────────────
input_file <- "MRDetect_output_winter_2025/Processed_R_outputs/BM_muts_plots_baseline/cfWGS_MRDetect_BM_data_updated_Feb2026.csv"
output_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/Dilution_series/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# Load and standardize columns
#   We compute the difference to healthy controls for the detection_rate metric.
#   We also flag if a sample is “significantly above healthy controls” using
#   z-score or mean+2*sd if z-score is missing.
# ──────────────────────────────────────────────────────────────────────────────
raw_df <- read_csv(input_file, show_col_types = FALSE)

# Map possible column names for healthy control mean/sd
mean_hc <- raw_df$mean_detection_rate_reads_checked_charm

sd_hc <- raw_df$sd_detection_rate_reads_checked_charm


rate <- raw_df$detection_rate_as_reads_detected_over_reads_checked
zscore <- raw_df$detection_rate_zscore_reads_checked_charm

working_df <- raw_df %>%
  mutate(
    mean_detection_rate_reads_checked_charm = mean_hc,
    sd_detection_rate_reads_checked_charm = sd_hc,
    detection_rate_diff_vs_hc = rate - mean_hc,
    not_sig_above_hc = dplyr::case_when(
      # Z-score >= 2 means significantly above healthy controls
      !is.na(zscore) ~ zscore < 2,
      !is.na(mean_hc) & !is.na(sd_hc) ~ rate <= (mean_hc + 2 * sd_hc),
      TRUE ~ NA
    )
  )

# ──────────────────────────────────────────────────────────────────────────────
# Identify tumor-high and tumor-low timepoints per patient
#   “High” = potential tumor-rich sample to spike in
#   “Low”  = background-like sample to dilute into
#   NOTE: Criteria loosened to capture more eligible cases.
# ──────────────────────────────────────────────────────────────────────────────
high_threshold <- 0.005  # 0.5% in decimal (loosened from 1%)
low_threshold  <- 5e-4   # 0.05% in decimal (loosened from 0.01%)

join_vars <- intersect(c("Patient", "Mut_source", "Filter_source"), names(working_df))

high_df <- working_df %>%
  filter(detection_rate_as_reads_detected_over_reads_checked >= high_threshold) %>%
  transmute(
    across(all_of(join_vars)),
    high_sample_id = Sample_ID_Bam,
    high_timepoint = timepoint_info_Bam,
    high_detection_rate = detection_rate_as_reads_detected_over_reads_checked,
    high_mean_hc = mean_detection_rate_reads_checked_charm,
    high_diff_vs_hc = detection_rate_diff_vs_hc,
    high_zscore = detection_rate_zscore_reads_checked_charm,
    high_not_sig_above_hc = not_sig_above_hc
  )

low_df <- working_df %>%
  filter(
    detection_rate_as_reads_detected_over_reads_checked <= low_threshold,
    not_sig_above_hc == TRUE
  ) %>%
  transmute(
    across(all_of(join_vars)),
    low_sample_id = Sample_ID_Bam,
    low_timepoint = timepoint_info_Bam,
    low_detection_rate = detection_rate_as_reads_detected_over_reads_checked,
    low_mean_hc = mean_detection_rate_reads_checked_charm,
    low_diff_vs_hc = detection_rate_diff_vs_hc,
    low_zscore = detection_rate_zscore_reads_checked_charm,
    low_not_sig_above_hc = not_sig_above_hc
  )

# Keep only one tumor-high and one tumor-low sample per patient (per join_vars)
#   - tumor-high: highest detection rate
#   - tumor-low: lowest detection rate
high_df_top <- high_df %>%
  group_by(across(all_of(join_vars))) %>%
  slice_max(order_by = high_detection_rate, n = 1, with_ties = FALSE) %>%
  ungroup()

low_df_top <- low_df %>%
  group_by(across(all_of(join_vars))) %>%
  slice_min(order_by = low_detection_rate, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join high and low samples by patient (and other shared metadata)
eligible_pairs <- high_df_top %>%
  inner_join(low_df_top, by = join_vars) %>%
  distinct()

# ──────────────────────────────────────────────────────────────────────────────
# Create dilution plan for each pair
#   For each eligible pair, we compute the fraction of the “tumor-high” sample
#   required to reach target tumor fractions of:
#   10%, 1%, 0.1%, 0.01%, 0.001%, 0.0001% (10^-1 to 10^-6)
#   with 3 technical replicates each.
#
#   IMPORTANT: We must account for the tumor-low background.
#     Mixed detection rate = (high_rate * f) + (low_rate * (1 - f))
#     Solve for f (fraction of tumor-high sample):
#       f = (target - low_rate) / (high_rate - low_rate)
#     If target <= low_rate, then the target is NOT achievable for that pair.
#
#   Example:
#     high_rate = 0.05 (5%), low_rate = 0.0004 (0.04%), target = 0.001 (0.1%)
#     f = (0.001 - 0.0004) / (0.05 - 0.0004) ≈ 0.0121 (1.21%)
# ──────────────────────────────────────────────────────────────────────────────

target_fractions <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)

if (nrow(eligible_pairs) > 0) {
  dilution_plan <- eligible_pairs %>%
    crossing(
      target_fraction = target_fractions
    ) %>%
    mutate(
      required_tumor_fraction = (target_fraction - low_detection_rate) / (high_detection_rate - low_detection_rate),
      required_normal_fraction = 1 - required_tumor_fraction,
      feasible = !is.na(required_tumor_fraction) &
        high_detection_rate > low_detection_rate &
        target_fraction > low_detection_rate &
        required_tumor_fraction > 0 & required_tumor_fraction < 1,
      mix_ratio = if_else(
        feasible,
        paste0(round(required_tumor_fraction, 6), ":", round(required_normal_fraction, 6)),
        NA_character_
      )
    )
} else {
  dilution_plan <- eligible_pairs
}

# ──────────────────────────────────────────────────────────────────────────────
# ALTERNATIVE ANALYSIS: Dilution plan using difference-to-healthy-controls
#   Instead of raw detection rates, use the background-corrected values
#   (difference vs healthy controls) to compute the required fractions.
#   This removes the contribution of the baseline healthy control signal.
# ──────────────────────────────────────────────────────────────────────────────

if (nrow(eligible_pairs) > 0) {
  dilution_plan_diff_vs_hc <- eligible_pairs %>%
    crossing(
      target_fraction = target_fractions
    ) %>%
    mutate(
      required_tumor_fraction = (target_fraction - low_diff_vs_hc) / (high_diff_vs_hc - low_diff_vs_hc),
      required_normal_fraction = 1 - required_tumor_fraction,
      feasible = !is.na(required_tumor_fraction) &
        high_diff_vs_hc > low_diff_vs_hc &
        target_fraction > low_diff_vs_hc &
        required_tumor_fraction > 0 & required_tumor_fraction < 1,
      mix_ratio = if_else(
        feasible,
        paste0(round(required_tumor_fraction, 6), ":", round(required_normal_fraction, 6)),
        NA_character_
      )
    )
} else {
  dilution_plan_diff_vs_hc <- eligible_pairs
}

# Filter to only patients/samples that can reach 1e-6 using difference-to-HC approach
patients_reaching_1e6 <- dilution_plan_diff_vs_hc %>%
  filter(
    target_fraction == 1e-6,
    feasible == TRUE
  ) %>%
  select(
    -target_fraction, -required_tumor_fraction,
    -required_normal_fraction, -feasible, -mix_ratio
  ) %>%
  distinct()

# Extract just the sample IDs for patients reaching 1e-6
patients_reaching_1e6_samples <- patients_reaching_1e6 %>%
  select(
    all_of(join_vars),
    high_sample_id,
    high_timepoint,
    low_sample_id,
    low_timepoint
  ) %>%
  distinct()

# Filter dilution plan to only 1e-6 rows where feasible
dilution_plan_diff_vs_hc_1e6_only <- dilution_plan_diff_vs_hc %>%
  filter(
    target_fraction == 1e-6,
    feasible == TRUE
  )

# Keep all rows for patients where 1e-6 is feasible
patients_with_feasible_1e6 <- dilution_plan_diff_vs_hc %>%
  filter(
    target_fraction == 1e-6,
    feasible == TRUE
  ) %>%
  select(all_of(join_vars)) %>%
  distinct()

dilution_plan_diff_vs_hc_for_1e6_patients <- dilution_plan_diff_vs_hc %>%
  semi_join(patients_with_feasible_1e6, by = join_vars)

# ──────────────────────────────────────────────────────────────────────────────
# Write outputs
#   1) eligible_dilution_pairs_Feb2026.csv
#        All eligible patient/sample pairs
#   2) dilution_plan_raw_Feb2026.csv
#        Dilution plan using raw detection rates (one row per target fraction)
#   3) dilution_plan_diff_vs_hc_Feb2026.csv
#        Dilution plan using difference-to-HC (one row per target fraction)
#   4) patients_reaching_1e6_diff_vs_hc_Feb2026.csv
#        Only patients/samples that can feasibly reach 1e-6 with diff-to-HC approach
# ──────────────────────────────────────────────────────────────────────────────
write_csv(eligible_pairs, file.path(output_dir, "eligible_dilution_pairs_Feb2026.csv"))
write_csv(dilution_plan, file.path(output_dir, "dilution_plan_raw_Feb2026.csv"))
write_csv(dilution_plan_diff_vs_hc, file.path(output_dir, "dilution_plan_diff_vs_hc_Feb2026.csv"))
write_csv(dilution_plan_diff_vs_hc_for_1e6_patients, file.path(output_dir, "dilution_plan_diff_vs_hc_1e6_capable_patients_Feb2026.csv"))
write_csv(patients_reaching_1e6, file.path(output_dir, "patients_reaching_1e6_diff_vs_hc_Feb2026.csv"))
write_csv(patients_reaching_1e6_samples, file.path(output_dir, "patients_reaching_1e6_sample_ids_Feb2026.csv"))

message("Eligible pairs written to: ", file.path(output_dir, "eligible_dilution_pairs_Feb2026.csv"))
message("Dilution plan (raw rates) written to: ", file.path(output_dir, "dilution_plan_raw_Feb2026.csv"))
message("Dilution plan (diff vs HC) written to: ", file.path(output_dir, "dilution_plan_diff_vs_hc_Feb2026.csv"))
message("Dilution plan (all rows for 1e-6 capable patients) written to: ", file.path(output_dir, "dilution_plan_diff_vs_hc_1e6_capable_patients_Feb2026.csv"))
message("Patients reaching 1e-6 (diff vs HC) written to: ", file.path(output_dir, "patients_reaching_1e6_diff_vs_hc_Feb2026.csv"))
message("Sample IDs for patients reaching 1e-6 written to: ", file.path(output_dir, "patients_reaching_1e6_sample_ids_Feb2026.csv"))
