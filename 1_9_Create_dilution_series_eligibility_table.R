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
# STEP 1: Load data and compute healthy-control-relative metrics
#
#   The core metric is detection_rate_as_reads_detected_over_reads_checked:
#   the fraction of cfDNA reads (at patient-specific mutation sites) that
#   show the somatic mutation. This is a proxy for circulating tumor fraction.
#
#   We add two derived columns:
#     • detection_rate_diff_vs_hc: the patient's detection rate minus the mean
#       healthy-control detection rate at those same sites. This removes the
#       background "noise floor" that exists even in true negatives.
#
#     • not_sig_above_hc: a Boolean flag. TRUE means this sample's signal is
#       NOT significantly elevated above healthy controls - i.e., it looks like
#       background. We use z-score < 2 as the criterion when available;
#       otherwise fall back to raw_rate <= mean_hc + 2*SD.
#       This flag is key for identifying genuine MRD-negative timepoints.
# ──────────────────────────────────────────────────────────────────────────────
raw_df <- read_csv(input_file, show_col_types = FALSE)

# Extract HC mean/SD and per-sample rate/z-score as vectors for clarity
mean_hc <- raw_df$mean_detection_rate_reads_checked_charm
sd_hc   <- raw_df$sd_detection_rate_reads_checked_charm
rate    <- raw_df$detection_rate_as_reads_detected_over_reads_checked
zscore  <- raw_df$detection_rate_zscore_reads_checked_charm

working_df <- raw_df %>%
  mutate(
    mean_detection_rate_reads_checked_charm = mean_hc,
    sd_detection_rate_reads_checked_charm   = sd_hc,
    # Background-corrected signal: how far above HC noise floor is this sample?
    detection_rate_diff_vs_hc = rate - mean_hc,
    # Is this sample statistically indistinguishable from healthy controls?
    # Preferred: use pre-computed z-score (z < 2 = not significantly elevated).
    # Fallback: raw rate within mean + 2*SD of HC distribution.
    not_sig_above_hc = dplyr::case_when(
      !is.na(zscore) ~ zscore < 2,
      !is.na(mean_hc) & !is.na(sd_hc) ~ rate <= (mean_hc + 2 * sd_hc),
      TRUE ~ NA
    )
  )

# ──────────────────────────────────────────────────────────────────────────────
# STEP 2: Identify tumor-high and tumor-low candidate timepoints per patient
#
#   The goal is to find one sample to use as the cfDNA SPIKE-IN SOURCE
#   (tumor-high) and one to use as the BACKGROUND MATRIX (tumor-low) to
#   dilute into. Together they form one patient's dilution pair.
#
#   -- Tumor-HIGH candidates: the spike-in source --
#   Criterion: detection_rate >= 0.5% (loosened from 1% to capture more pairs)
#   These are timepoints where the patient had measurable circulating tumor
#   cfDNA -- typically at diagnosis, progression, or relapse.
#
#   -- Tumor-LOW candidates: the background matrix --
#   This is NOT simply "lowest detection rate". There are TWO gates:
#     Gate 1: detection_rate <= 0.05%  (raw rate floor)
#     Gate 2: not_sig_above_hc == TRUE  (z-score < 2 vs healthy controls)
#   Gate 2 is critical: a sample can have a numerically low detection rate
#   but still be statistically elevated above the HC noise floor. Both gates
#   must pass so the background matrix is genuinely MRD-negative.
#
#   NOTE: Thresholds were loosened from the original (1% / 0.01%) to maximize
#   the number of eligible patient pairs.
# ──────────────────────────────────────────────────────────────────────────────
high_threshold <- 0.005  # 0.5% in decimal (loosened from 1%)
low_threshold  <- 5e-4   # 0.05% in decimal (loosened from 0.01%)

# join_vars defines the grouping level for "per patient" selection.
# Includes Mut_source and Filter_source so that pairs are matched on the
# same mutation set and filtering strategy (apples-to-apples comparison).
join_vars <- intersect(c("Patient", "Mut_source", "Filter_source"), names(working_df))

# Tumor-high candidates: all timepoints above the high threshold
high_df <- working_df %>%
  filter(detection_rate_as_reads_detected_over_reads_checked >= high_threshold) %>%
  transmute(
    across(all_of(join_vars)),
    high_sample_id   = Sample_ID_Bam,
    high_timepoint   = timepoint_info_Bam,
    high_detection_rate = detection_rate_as_reads_detected_over_reads_checked,
    high_mean_hc     = mean_detection_rate_reads_checked_charm,
    high_diff_vs_hc  = detection_rate_diff_vs_hc,
    high_zscore      = detection_rate_zscore_reads_checked_charm,
    high_not_sig_above_hc = not_sig_above_hc
  )

# Tumor-low candidates: must pass BOTH gates (rate floor + HC z-score gate).
# A patient who achieved deep remission will have a timepoint where their
# cfDNA detection rate is low AND their z-score vs HC is < 2 (not elevated).
low_df <- working_df %>%
  filter(
    detection_rate_as_reads_detected_over_reads_checked <= low_threshold,
    not_sig_above_hc == TRUE  # z-score < 2: statistically at HC background level
  ) %>%
  transmute(
    across(all_of(join_vars)),
    low_sample_id   = Sample_ID_Bam,
    low_timepoint   = timepoint_info_Bam,
    low_detection_rate = detection_rate_as_reads_detected_over_reads_checked,
    low_mean_hc     = mean_detection_rate_reads_checked_charm,
    low_diff_vs_hc  = detection_rate_diff_vs_hc,
    low_zscore      = detection_rate_zscore_reads_checked_charm,
    low_not_sig_above_hc = not_sig_above_hc
  )

# STEP 3: Select the single best high and low timepoint per patient
#
#   If a patient has multiple qualifying timepoints for either category, we
#   pick the single most extreme one. This maximises the signal-to-background
#   ratio, giving the widest possible dynamic range:
#     • tumor-high: HIGHEST detection rate (strongest signal = best spike-in)
#     • tumor-low:  LOWEST detection rate (cleanest background = best matrix)
high_df_top <- high_df %>%
  group_by(across(all_of(join_vars))) %>%
  slice_max(order_by = high_detection_rate, n = 1, with_ties = FALSE) %>%
  ungroup()

low_df_top <- low_df %>%
  group_by(across(all_of(join_vars))) %>%
  slice_min(order_by = low_detection_rate, n = 1, with_ties = FALSE) %>%
  ungroup()

# STEP 4: Eligibility gate -- patients must have BOTH a high AND a low sample
#
#   inner_join keeps only patients present in both tables. Patients who were
#   always MRD-positive (no qualifying low sample) or who have no high-signal
#   timepoints (no qualifying high sample) are excluded.
eligible_pairs <- high_df_top %>%
  inner_join(low_df_top, by = join_vars) %>%
  distinct()

# ──────────────────────────────────────────────────────────────────────────────
# STEP 5: Compute the dilution plan (using raw detection rates)
#
#   For each eligible pair, we calculate the physical mixing fractions needed
#   to achieve a series of TARGET tumor fractions spanning 6 orders of magnitude:
#   10^-1 (10%) down to 10^-6 (0.0001%).
#
#   KEY INSIGHT: the tumor-low sample is not perfectly zero - it has a small
#   residual detection rate (the background). If we ignore this and treat the
#   tumor-low sample as a pure zero, our target fractions will be off.
#   Instead we model the mixture explicitly:
#
#     Mixed detection rate = (high_rate × f) + (low_rate × (1 - f))
#
#   Where f = fraction of the tumor-high sample in the final mix.
#   Solving for f to hit a given target:
#
#     f = (target - low_rate) / (high_rate - low_rate)
#
#   This means:
#     • If target ≤ low_rate  → not achievable (you can't go below background)
#     • If f ≤ 0 or f ≥ 1    → not feasible (outside 0–100% mixing range)
#     • feasible == TRUE      → this dilution level can actually be made
#
#   Example with real numbers:
#     high_rate = 0.05 (5%), low_rate = 0.0004 (0.04%), target = 0.001 (0.1%)
#     f = (0.001 - 0.0004) / (0.05 - 0.0004) = 0.0006 / 0.0496 ≈ 1.21% tumor-high
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
# STEP 6: Alternative dilution plan using HC-corrected (diff vs HC) rates
#
#   The raw detection rate includes a universal background noise floor that
#   exists even in healthy controls (technical noise from sequencing errors,
#   mapping artifacts, etc.). Using the raw rate to compute mixing fractions
#   means we may be over-estimating how much of the "low" sample's signal
#   is actually tumor-derived.
#
#   This alternative uses detection_rate_diff_vs_hc = (sample_rate - mean_HC_rate)
#   for both the high and low samples. After subtraction, the shared HC noise
#   floor cancels out, and the values represent only the tumor-specific signal.
#
#   The mixing formula is the same:
#     f = (target - low_diff_vs_hc) / (high_diff_vs_hc - low_diff_vs_hc)
#
#   This is generally the preferred approach for determining feasibility at
#   very low target fractions (e.g., 10^-5 or 10^-6) where the noise floor
#   would otherwise dominate the raw rate.
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

# ──────────────────────────────────────────────────────────────────────────────
# STEP 7: Identify the "gold standard" subset - patients whose pairs span
#         far enough to feasibly reach the 10^-6 (0.0001%) target level
#
#   Reaching 10^-6 requires both:
#     • A very high tumor-high detection rate (lots of signal to dilute down)
#     • A very low tumor-low detection rate (minimal background floor)
#   Not all eligible pairs can achieve this. We filter to only those where
#   the 1e-6 level is feasible under the HC-corrected model, then extract
#   all dilution levels for those patients (not just the 1e-6 row).
# ──────────────────────────────────────────────────────────────────────────────

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
