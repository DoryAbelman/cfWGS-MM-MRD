# =============================================================================
# MRDetect_processing.R
# Project:  cfWGS MRDetect (Winter 2025)
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_8_Process_Cumulative_VAFs_MRDetect.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script processes MRDetect outputs for
#     patient samples that feed model training, clinical concordance, and
#     survival analyses.
#
# Author:   Dory Abelman
# Date:     January 2025
# Last Updated: May 2025
#
# Purpose:
#   1. Read all MRDetect CSV outputs.
#   2. Annotate each record with source file, sample metadata, and z-scores
#      based on CHARM_healthy controls.
#   3. Filter to cfDNA timepoints and export
#      both raw and z-scored tables for downstream plotting.
#
# Dependencies:
#   • readr, data.table, tidyverse (dplyr, tidyr, stringr), openxlsx
#   • ggplot2, ggbreak, patchwork, scales, conflicted
#
# Input Files:
#   • ../MRDetect_output_winter_2025/MRDetect_outputs/*.csv
#   • ../combined_clinical_data_updated_April2025.csv
#
# Output Directory (created if necessary):
#   • MRDetect_output_winter_2025/Processed_R_outputs/
#   • Writes:
#       - cfWGS_Winter2025All_MRDetect_May2025.txt
#       - cfWGS_Winter2025All_MRDetect_May2025.rds
#       - cfWGS_Winter2025All_MRDetect_with_Zscore_May2025.txt
#       - cfWGS_Winter2025All_MRDetect_with_Zscore_May2025.rds
#       - cfWGS MRDetect BM data updated May.csv
#       - cfWGS MRDetect Blood data updated May.csv
#
# Usage:
#   Rscript MRDetect_processing.R
# =============================================================================
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#

# ──────────────────────────────────────────────────────────────────────────────
# 1) Load libraries
# ──────────────────────────────────────────────────────────────────────────────
library(readr)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(ggplot2)
library(ggbreak)
library(patchwork)
library(scales)
#library(conflicted)
library(lubridate)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

restrict_to_single_baseline_mutation_source <- function(dat, audit_dir, audit_prefix) {
  baseline_labels <- c("Baseline", "Diagnosis")

  baseline_sources <- dat %>%
    filter(timepoint_info %in% baseline_labels) %>%
    mutate(
      source_key = coalesce(
        na_if(as.character(VCF_factor), ""),
        na_if(as.character(VCF_clean), ""),
        na_if(as.character(VCF), ""),
        na_if(as.character(Sample_ID), ""),
        paste(as.character(Patient), as.character(Mut_source), as.character(Timepoint), as.character(timepoint_info), sep = "|")
      ),
      source_date = suppressWarnings(as.Date(Date_of_sample_collection_Sample_ID)),
      timepoint_rank = case_when(
        str_to_upper(Timepoint) %in% c("T0", "TP0", "D0", "0", "01", "1") ~ 1,
        str_to_upper(Timepoint) %in% c("T1", "TP1") ~ 2,
        str_detect(Timepoint, "^[Tt]?[0-9]+$") ~ as.numeric(str_remove(str_to_upper(Timepoint), "^T")) + 10,
        TRUE ~ 999
      ),
      label_rank = case_when(
        timepoint_info == "Diagnosis" ~ 1,
        timepoint_info == "Baseline" ~ 2,
        TRUE ~ 9
      )
    ) %>%
    distinct(
      Patient, Mut_source, source_key, Sample_ID, Timepoint, timepoint_info,
      source_date, timepoint_rank, label_rank
    )

  duplicate_groups <- baseline_sources %>%
    count(Patient, Mut_source, name = "n_baseline_diagnosis_mutation_sources") %>%
    filter(n_baseline_diagnosis_mutation_sources > 1)

  if (nrow(duplicate_groups) == 0L) {
    return(dat)
  }

  selected_sources <- baseline_sources %>%
    semi_join(duplicate_groups, by = c("Patient", "Mut_source")) %>%
    arrange(
      Patient, Mut_source,
      coalesce(source_date, as.Date("9999-12-31")),
      timepoint_rank,
      label_rank,
      source_key
    ) %>%
    group_by(Patient, Mut_source) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      Patient, Mut_source, selected_source_key = source_key,
      selected_source = paste(Sample_ID, Timepoint, timepoint_info, source_date, sep = "|")
    )

  baseline_decisions <- baseline_sources %>%
    semi_join(duplicate_groups, by = c("Patient", "Mut_source")) %>%
    left_join(selected_sources, by = c("Patient", "Mut_source")) %>%
    mutate(
      retained = source_key == selected_source_key,
      mutation_source = paste(Sample_ID, Timepoint, timepoint_info, source_date, sep = "|")
    ) %>%
    arrange(Patient, Mut_source, desc(retained), timepoint_rank, source_key)

  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    baseline_decisions,
    file.path(audit_dir, paste0(audit_prefix, "_baseline_mutation_source_deduplication_audit.csv"))
  )

  dat %>%
    mutate(
      source_key = coalesce(
        na_if(as.character(VCF_factor), ""),
        na_if(as.character(VCF_clean), ""),
        na_if(as.character(VCF), ""),
        na_if(as.character(Sample_ID), ""),
        paste(as.character(Patient), as.character(Mut_source), as.character(Timepoint), as.character(timepoint_info), sep = "|")
      )
    ) %>%
    left_join(selected_sources, by = c("Patient", "Mut_source")) %>%
    filter(is.na(selected_source_key) | source_key == selected_source_key) %>%
    select(-source_key, -selected_source_key, -selected_source)
}

# resolve common conflicts
#conflicted::conflicts_prefer("dplyr::mutate")
#conflicted::conflicts_prefer("dplyr::filter")
#conflicted::conflicts_prefer("dplyr::select")
#conflicted::conflicts_prefer("dplyr::summarize")


# ──────────────────────────────────────────────────────────────────────────────
# 2) Set up paths and create output directory
# ──────────────────────────────────────────────────────────────────────────────
input_root <- "MRDetect_output_winter_2025/MRDetect_outputs/"
outdir     <- "MRDetect_output_winter_2025/Processed_R_outputs/"
project    <- "cfWGS_Winter2025"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}


# ──────────────────────────────────────────────────────────────────────────────
# 3) Read and combine all CSVs, label with ‘source_file’, ‘Mut_source’ & ‘Filter_source’
# ──────────────────────────────────────────────────────────────────────────────
csv_files <- list.files(input_root, pattern = "\\.csv$", full.names = TRUE)
spring2026_mrdetect_files <- spring2026_revision_files(
  "MRDetect_outputs",
  "^MRDetect_all_RESULTS_combined_with_source[.]csv$"
)
csv_files <- unique(c(csv_files, spring2026_mrdetect_files))
if (!length(csv_files)) {
  stop("No MRDetect CSV files found in historical or Spring 2026 revision inputs.", call. = FALSE)
}

read_and_label <- function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  if (!"filename" %in% names(df)) {
    if ("source_file" %in% names(df)) {
      df$filename <- df$source_file
    } else {
      df$filename <- basename(file)
    }
  }
  if (!"source_file" %in% names(df)) {
    df$source_file <- basename(file)
  }
  df$input_source_file <- basename(file)
  return(df)
}
all_files <- bind_rows(lapply(csv_files, read_and_label))

all_files <- all_files %>%
  filter(
    !str_detect(source_file, "^M4CHIP_"),
    !str_detect(filename, "^M4CHIP_"),
    !str_detect(input_source_file, "dilution_series")
  )

all_files <- all_files %>%
  mutate(
    Mut_source = case_when(
      str_detect(input_source_file, "BM_muts") | str_detect(source_file, "BM_muts") ~ "BM_cells",
      str_detect(input_source_file, "Blood_muts") | str_detect(source_file, "Blood_muts") ~ "Blood",
      TRUE                                   ~ NA_character_
    ),
    Filter_source = case_when(
      str_detect(source_file, "encode_only") ~ "Encode_only",
      str_detect(source_file, "STR_encode") |
        str_detect(source_file, "encode_STR_removed") ~ "STR_encode",
      TRUE                                   ~ NA_character_
    )
  ) %>%
  # drop 'VCF' if it already exists
  select(-any_of("VCF")) %>%
  mutate(
    # 1) remove leading "./" from filename
    filename = str_remove(filename, "^\\./"),
    # 2) rebuild VCF from whatever follows "_VS_" in the (cleaned) filename
    VCF = sub(".*_VS_", "\\1", filename),
    # 3) strip “.mutect2…” and “.fil…” from the VCF string itself
    VCF_clean = VCF %>%
      str_remove("\\.mutect2.*") %>%
      str_remove("\\.fil.*") %>%
      str_remove("\\.somatic.*")
  )

# fix one truncated VCF_clean
all_files <- all_files %>%
  mutate(
    VCF_clean = ifelse(
      VCF_clean == "TFRIM4_0189_Bm_P_WG_ZC-02",
      "TFRIM4_0189_Bm_P_WG_ZC-02-01-O-DNA",
      VCF_clean
    )
  )


# ──────────────────────────────────────────────────────────────────────────────
# 4) Load clinical metadata and unify ‘VCF_clean’ naming
# ──────────────────────────────────────────────────────────────────────────────
cfWGS_metadata <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv",
  include_revision_extra = TRUE
) %>%
  mutate(
    VCF_clean_merge = str_remove(Bam, "\\.filter.*"),
    # correct internal PG→WG for those five patient‐specific IDs
    VCF_clean_merge = if_else(
      VCF_clean_merge %in% c(
        "TFRIM4_0031_Bm_P_PG_M4-CA-02-01-O-DNA",
        "TFRIM4_0032_Bm_P_PG_M4-HP-01-01-O-DNA",
        "TFRIM4_0033_Bm_P_PG_M4-MJ-06-01-O-DNA",
        "TFRIM4_0034_Bm_P_PG_M4-VA-02-01-O-DNA",
        "TFRIM4_0035_Bm_P_PG_M4-VA-06-01-O-DNA"
      ),
      str_replace(VCF_clean_merge, "PG", "WG"),
      VCF_clean_merge
    )
  )

spring2026_panel_metadata <- cfWGS_metadata %>%
  filter(!is.na(mutect2_pair_id), nzchar(mutect2_pair_id)) %>%
  transmute(
    VCF_clean = mutect2_pair_id,
    VCF_panel_patient = Patient,
    VCF_panel_sample_id = Sample_ID,
    VCF_panel_sample_type = Sample_type,
    VCF_panel_timepoint_info = timepoint_info
  ) %>%
  distinct(VCF_clean, .keep_all = TRUE)

if (nrow(spring2026_panel_metadata) > 0) {
  all_files <- all_files %>%
    left_join(spring2026_panel_metadata, by = "VCF_clean") %>%
    mutate(
      Mut_source = case_when(
        VCF_panel_sample_type == "BM_cells" ~ "BM_cells",
        VCF_panel_sample_type == "Blood_plasma_cfDNA" ~ "Blood",
        TRUE ~ Mut_source
      ),
      VCF_panel_is_allowed_panel_timepoint = case_when(
        is.na(VCF_panel_timepoint_info) ~ TRUE,
        VCF_panel_timepoint_info %in% c("Baseline", "Diagnosis") ~ TRUE,
        TRUE ~ FALSE
      )
    )

  excluded_disallowed_panel <- all_files %>%
    filter(
      !is.na(VCF_panel_timepoint_info),
      !VCF_panel_is_allowed_panel_timepoint
    )
  if (nrow(excluded_disallowed_panel) > 0) {
    readr::write_csv(
      excluded_disallowed_panel,
      file.path(outdir, "spring2026_mrdetect_excluded_disallowed_vcf_panels.csv")
    )
  }

  all_files <- all_files %>%
    filter(VCF_panel_is_allowed_panel_timepoint) %>%
    select(-VCF_panel_is_allowed_panel_timepoint)
}

# ──────────────────────────────────────────────────────────────────────────────
# 5) Standardize column names: replace spaces & dots with underscores
# ──────────────────────────────────────────────────────────────────────────────
# 1) remove trailing spaces
colnames(all_files) <- gsub("\\s+$", "", colnames(all_files))
# 2) replace any runs of spaces with "_"
colnames(all_files) <- gsub("\\s+", "_", colnames(all_files))
# 3) replace "." with "_"
colnames(all_files) <- gsub("\\.", "_", colnames(all_files))


# -----------------------------------------------------------------------------
# 6) MRDetect parser compatibility correction
# -----------------------------------------------------------------------------
#
# The Winter 2025 MRDetect CSV files used for the manuscript were generated by a
# parser version that assigned three count columns in the wrong order. The
# manuscript analysis has always corrected that ordering before calculating
# detection-rate features. Keep this TRUE for the preserved manuscript input set.
# If future MRDetect CSVs are produced by a corrected parser, set this to FALSE
# only after confirming the input column definitions against parser documentation.
#
# Key definitions after the correction:
#   sites_detected  = number of patient-specific mutant sites with at least one read
#   reads_detected  = number of individual reads carrying a somatic mutation
#   total_reads     = total reads inspected at those sites (denominator for detection_rate)
apply_mrdetect_parser_column_correction <- TRUE
required_mrdetect_count_cols <- c("sites_detected", "reads_detected", "total_reads")
missing_mrdetect_count_cols <- setdiff(required_mrdetect_count_cols, names(all_files))
if (length(missing_mrdetect_count_cols) > 0) {
  stop(
    "MRDetect count columns are missing after column-name standardization: ",
    paste(missing_mrdetect_count_cols, collapse = ", "),
    call. = FALSE
  )
}

if (isTRUE(apply_mrdetect_parser_column_correction)) {
  temp_sites <- all_files$sites_detected
  temp_reads <- all_files$reads_detected
  temp_total <- all_files$total_reads

  all_files <- all_files %>%
    mutate(
      sites_detected = temp_reads,
      reads_detected = temp_total,
      total_reads    = temp_sites
    )

  rm(temp_sites, temp_reads, temp_total)
}



# ──────────────────────────────────────────────────────────────────────────────
# 7) Annotate each row with ‘Study’ based on BAM prefix
# ──────────────────────────────────────────────────────────────────────────────# BAM file naming conventions used to identify cohort:
#   TFRI...  = M4 cohor)
#   MY...    = MyP cohort 
#   EGA...   = Landau et al. external validation cohort (from EGA repository)
#   HCCCFD.. = HCC patients used as non-MM cancer controls
#   HCC-...  = HCC healthy donors used as cancer-free controls
#   TGL...   = CHARM healthy donors (TGL49 panel; PRIMARY HEALTHY CONTROL POPULATION
#              used for z-score normalization throughout this analysis)
#   SPORE... = SPORE cohort (secondary MM cohort)
all_files <- all_files %>%
  mutate(
    Study = case_when(
      is_xplus_charm_healthy_bam(BAM) ~ "XPLUS_CHARM_healthy",
      str_detect(BAM, "^TFRI")       ~ "M4",
      str_detect(BAM, "^MY")         ~ "MyP",
      str_detect(BAM, "^EGA")        ~ "Landau",
      str_detect(BAM, "^HCCCFD")     ~ "HCC",
      str_detect(BAM, "^HCC\\-")      ~ "HCC_healthy",
      str_detect(BAM, "^TGL")        ~ "CHARM_healthy",
      str_detect(BAM, "^SPORE")      ~ "SPORE",
      TRUE                            ~ NA_character_
    )
  )


# ──────────────────────────────────────────────────────────────────────────────
# 8) Merge MRDetect outputs with clinical metadata
# ──────────────────────────────────────────────────────────────────────────────
tmp_meta <- cfWGS_metadata %>%
  select(-Bam) %>%
  rename(Study_VCF = Study)

Merged_MRDetect <- all_files %>%
  left_join(tmp_meta, by = c("VCF_clean" = "VCF_clean_merge"))

if ("mutect2_pair_id" %in% names(tmp_meta)) {
  tmp_meta_pair <- tmp_meta %>%
    filter(!is.na(mutect2_pair_id), nzchar(mutect2_pair_id)) %>%
    distinct(mutect2_pair_id, .keep_all = TRUE)

  pair_joined <- all_files %>%
    left_join(tmp_meta_pair, by = c("VCF_clean" = "mutect2_pair_id"))

  fill_cols <- intersect(names(Merged_MRDetect), names(pair_joined))
  metadata_cols <- setdiff(fill_cols, names(all_files))
  for (col in metadata_cols) {
    Merged_MRDetect[[col]] <- dplyr::coalesce(Merged_MRDetect[[col]], pair_joined[[col]])
  }
}

## Fix issue with ones that don't copy over

# ──────────────────────────────────────────────────────────────────────────────
# 9) Add BAM‐level sample info (Sample_ID, Patient, Sample_type, timepoint_info)
# ──────────────────────────────────────────────────────────────────────────────
# Make a helper to normalise BAM strings
normalize_bam <- function(x) {
  x %>%
    # remove internal _PG_ or _WG_ tokens
    str_replace_all("_[PW]G_", "_") %>%
    str_replace("^([PW]G)_", "") %>%
    str_replace("_([PW]G)$", "") %>%
    # remove the .filter. segment if present
    str_replace("\\.filter\\.", ".")
}

# Build bam_info with normalised key
bam_info <- cfWGS_metadata %>%
  mutate(Bam_norm = normalize_bam(Bam)) %>%
  select(Bam, Bam_norm, Sample_ID, Patient, Sample_type, timepoint_info) %>%
  rename_with(~ paste0(.x, "_Bam"),
              .cols = c(Sample_ID, Patient, Sample_type, timepoint_info))

# Add the same normalized key to Merged_MRDetect
Merged_MRDetect <- Merged_MRDetect %>%
  mutate(BAM_norm = str_replace_all(BAM, "_[PW]G_", "_") |> 
           str_replace("_([PW]G)$", "") |> 
           str_replace("\\.filter\\.", ".") |> 
           str_replace("^([PW]G)_", ""))

# Join using the normalized key
Merged_MRDetect <- Merged_MRDetect %>%
  left_join(bam_info, by = c("BAM_norm" = "Bam_norm"))

# optional: drop the helper columns if you don’t need them in the final table
Merged_MRDetect <- Merged_MRDetect %>% select(-BAM_norm)

## Old way
# bam_info <- cfWGS_metadata %>%
#   select(Bam, Sample_ID, Patient, Sample_type, timepoint_info) %>%
#   rename_with(~ paste0(.x, "_Bam"), .cols = -Bam)
# 
# Merged_MRDetect <- Merged_MRDetect %>%
#   left_join(bam_info, by = c("BAM" = "Bam"))
# 

Merged_MRDetect %>%
  filter(is.na(Patient_Bam)) %>%
  pull(BAM) %>% 
  unique() # seems good
  
# ──────────────────────────────────────────────────────────────────────────────
# 10) Classify matched vs. unmatched plasma; force CHARM_healthy→cfDNA
# ──────────────────────────────────────────────────────────────────────────────
# MRDetect runs each cfDNA sample (BAM) against the VCF from one patient's
# bone-marrow mutations. That VCF is the "patient" the sample is being tested against.
#   Matched_plasma   = BAM and VCF come from the SAME patient (the true signal)
#   Unmatched_plasma = BAM and VCF come from DIFFERENT patients
#                      (used as cross-patient negative controls)
# All CHARM_healthy donors are always run against a patient's VCF, so their
# Sample_type_Bam field is set to Blood_plasma_cfDNA here.
Merged_MRDetect <- Merged_MRDetect %>%
  mutate(
    plotting_type       = if_else(Patient_Bam == Patient, "Matched_plasma", "Unmatched_plasma"),
    Sample_type_Bam     = if_else(Study == "CHARM_healthy", "Blood_plasma_cfDNA", Sample_type_Bam)
  )


# ──────────────────────────────────────────────────────────────────────────────
# 11) Subset to cfDNA timepoints only.
# The queried BAM may be a longitudinal cfDNA sample, but the tumor-mutation VCF
# panel above is restricted to Diagnosis/Baseline labels. Relapse-context samples
# that are the baseline for a follow-up series should be labelled Baseline in the
# metadata before entering this analysis.
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect_cfDNA <- Merged_MRDetect %>%
  filter(Sample_type_Bam == "Blood_plasma_cfDNA",
         timepoint_info_Bam %in% c("Baseline","Diagnosis","Progression"))


# ──────────────────────────────────────────────────────────────────────────────
# 12) Compute CHARM_healthy means & SDs → join back for z-scores
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect$VCF_factor <- factor(Merged_MRDetect$VCF, levels = unique(Merged_MRDetect$VCF))

# Z-SCORE NORMALIZATION STRATEGY:
# For each patient VCF (VCF_factor), we compute a healthy-control reference
# distribution using CHARM_healthy donors run against that same VCF.
# This is very important: each patient has a different mutation set (different VCF),
# so the background noise level differs per patient. By grouping on VCF_factor,
# we get a VCF-specific mean and SD for healthy controls, ensuring that the
# z-score for a patient sample is relative to healthy donors tested against
# that exact same set of patient mutations.
#   z = (patient_rate - mean_HC_rate) / sd_HC_rate
# High z-score (>>2) → detection rate significantly above healthy-control noise,
# suggesting true circulating tumor DNA.
healthy_reference_candidates <- Merged_MRDetect %>%
  filter(Study %in% c("CHARM_healthy", "XPLUS_CHARM_healthy", "HCC_healthy")) %>%
  mutate(
    healthy_reference_tier = case_when(
      Study == "CHARM_healthy" ~ "CHARM_healthy",
      TRUE ~ "non_mm_healthy_fallback"
    )
  )

healthy_reference_choice <- healthy_reference_candidates %>%
  count(VCF_factor, Mut_source, Filter_source, healthy_reference_tier, name = "n_healthy_reference_rows") %>%
  mutate(healthy_reference_priority = case_when(
    healthy_reference_tier == "CHARM_healthy" ~ 1L,
    healthy_reference_tier == "non_mm_healthy_fallback" ~ 2L,
    TRUE ~ 99L
  )) %>%
  group_by(VCF_factor, Mut_source, Filter_source) %>%
  slice_min(order_by = healthy_reference_priority, n = 1, with_ties = FALSE) %>%
  ungroup()

write_csv(
  healthy_reference_choice,
  file.path(outdir, "spring2026_mrdetect_healthy_reference_choice_audit.csv")
)

zscore_lookup <- healthy_reference_candidates %>%
  inner_join(
    healthy_reference_choice %>%
      select(VCF_factor, Mut_source, Filter_source, healthy_reference_tier),
    by = c("VCF_factor", "Mut_source", "Filter_source", "healthy_reference_tier")
  ) %>%
  select(VCF_factor, Mut_source, Filter_source,
         detection_rate,
         detection_rate_as_reads_detected_over_reads_checked,
         detection_rate_as_reads_detected_over_total_reads,
         sites_detection_rate) %>%
  group_by(VCF_factor, Mut_source, Filter_source) %>%
  summarize(
    mean_det_charm                     = mean(detection_rate, na.rm = TRUE),
    sd_det_charm                       = sd(detection_rate, na.rm = TRUE),
    mean_det_checked_charm             = mean(detection_rate_as_reads_detected_over_reads_checked, na.rm = TRUE),
    sd_det_checked_charm               = sd(detection_rate_as_reads_detected_over_reads_checked, na.rm = TRUE),
    mean_det_total_charm               = mean(detection_rate_as_reads_detected_over_total_reads, na.rm = TRUE),
    sd_det_total_charm                 = sd(detection_rate_as_reads_detected_over_total_reads, na.rm = TRUE),
    mean_sites_charm                   = mean(sites_detection_rate, na.rm = TRUE),
    sd_sites_charm                     = sd(sites_detection_rate, na.rm = TRUE)
  ) %>%
  ungroup()

Merged_MRDetect_zscore <- Merged_MRDetect %>%
  left_join(zscore_lookup, by = c("VCF_factor","Mut_source","Filter_source")) %>%
  mutate(
    detection_rate_zscore_charm                   = (detection_rate - mean_det_charm) / sd_det_charm,
    detection_rate_zscore_reads_checked_charm     = (detection_rate_as_reads_detected_over_reads_checked - mean_det_checked_charm) / sd_det_checked_charm,
    detection_rate_zscore_total_reads_charm       = (detection_rate_as_reads_detected_over_total_reads - mean_det_total_charm) / sd_det_total_charm,
    sites_rate_zscore_charm                       = (sites_detection_rate - mean_sites_charm) / sd_sites_charm,
    # Difference vs healthy controls (absolute values)
    detection_rate_diff_vs_hc                     = detection_rate - mean_det_charm,
    detection_rate_checked_diff_vs_hc             = detection_rate_as_reads_detected_over_reads_checked - mean_det_checked_charm,
    detection_rate_total_diff_vs_hc               = detection_rate_as_reads_detected_over_total_reads - mean_det_total_charm,
    sites_rate_diff_vs_hc                         = sites_detection_rate - mean_sites_charm,
    # Difference vs healthy controls (as percentages)
    detection_rate_checked_diff_vs_hc_pct         = (detection_rate_checked_diff_vs_hc) * 100
  )


# ──────────────────────────────────────────────────────────────────────────────
# 13) Identify and drop any duplicate runs (keep max total_reads)
# ──────────────────────────────────────────────────────────────────────────────
dup_keys <- c("BAM","Mut_source","Filter_source","VCF","VCF_clean","Study",
              "Sample_ID_Bam","Patient_Bam","VCF_factor")
Merged_MRDetect_zscore <- Merged_MRDetect_zscore %>%
  group_by(across(all_of(dup_keys))) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup()


# ──────────────────────────────────────────────────────────────────────────────
# 14) Save both “raw” and “z-scored” MRDetect tables
# ──────────────────────────────────────────────────────────────────────────────
base_name_raw    <- paste0(project, "All_MRDetect_Sep2025")
base_name_zscore <- paste0(project, "All_MRDetect_with_Zscore_Sep2025")

write.table(Merged_MRDetect,
            file = file.path(outdir, paste0(base_name_raw, ".txt")),
            sep = "\t", row.names = FALSE)

saveRDS(Merged_MRDetect,
        file = file.path(outdir, paste0(base_name_raw, ".rds")))

write.table(Merged_MRDetect_zscore,
            file = file.path(outdir, paste0(base_name_zscore, ".txt")),
            sep = "\t", row.names = FALSE)

saveRDS(Merged_MRDetect_zscore,
        file = file.path(outdir, paste0(base_name_zscore, ".rds")))

message("→ MRDetect tables (raw & z-scored) written to: ", outdir)


# ──────────────────────────────────────────────────────────────────────────────
# 14a) Export individual detection rates for baseline/diagnosis + healthy controls
# ──────────────────────────────────────────────────────────────────────────────
# Extract all baseline/diagnosis samples from MM cohort + all CHARM_healthy samples
detection_rates_all_samples <- Merged_MRDetect_zscore %>%
  filter(
    # Include MM baseline/diagnosis samples
    (Timepoint %in% c("01", "T0") & Study == "CHARM_healthy")
  ) %>%
  select(
    Study, Patient, Patient_Bam, Sample_ID_Bam, BAM, VCF, VCF_clean,
    Sample_type_Bam, timepoint_info_Bam, Mut_source, Filter_source,
    sites_checked, reads_checked, sites_detected, reads_detected, total_reads,
    detection_rate, 
    detection_rate_as_reads_detected_over_reads_checked,
    detection_rate_as_reads_detected_over_total_reads,
    sites_detection_rate,
    detection_rate_zscore_charm,
    detection_rate_zscore_reads_checked_charm,
    detection_rate_zscore_total_reads_charm,
    sites_rate_zscore_charm
  ) %>%
  arrange(Study, Patient_Bam, Mut_source, Filter_source)

# Export as CSV and RDS
write_csv(detection_rates_all_samples,
          file = file.path(outdir, "All_detection_rates_baseline_and_controls_Feb2026.csv"))
saveRDS(detection_rates_all_samples,
        file = file.path(outdir, "All_detection_rates_baseline_and_controls_Feb2026.rds"))

message("→ Individual detection rates (baseline + healthy controls) written to: ", outdir)
message("  Total samples exported: ", nrow(detection_rates_all_samples))


## See what is NA
# show BAM file paths where Patient_Bam is NA
Merged_MRDetect_zscore %>%
  filter(is.na(Patient_Bam)) %>%
  pull(BAM) %>% unique()

# ──────────────────────────────────────────────────────────────────────────────
# 15) Additional Processing - BM muts
# ──────────────────────────────────────────────────────────────────────────────
### Now Process Further
### First BM Muts

# 1) Start from your z-scored MRDetect table
df <- Merged_MRDetect_zscore %>%
  filter(plotting_type == "Matched_plasma",
         Mut_source      == "BM_cells",
         Filter_source   == "STR_encode")

# Rename df's columns to the names num_cols vector expects
rename_map <- c(
  mean_detection_rate_charm                       = "mean_det_charm",
  sd_detection_rate_charm                         = "sd_det_charm",
  mean_detection_rate_reads_checked_charm         = "mean_det_checked_charm",
  sd_detection_rate_reads_checked_charm           = "sd_det_checked_charm",
  mean_detection_rate_total_reads_charm           = "mean_det_total_charm",
  sd_detection_rate_total_reads_charm             = "sd_det_total_charm",
  mean_sites_rate_charm                           = "mean_sites_charm",
  sd_sites_rate_charm                             = "sd_sites_charm"
)

df <- df %>% rename(!!!rename_map)

# 2) Collapse duplicate BAM runs by averaging all numeric QC/MRD metrics
num_cols <- c(
  "sites_checked", "reads_checked", "sites_detected", "reads_detected",
  "total_reads", "detection_rate",
  "detection_rate_as_reads_detected_over_reads_checked",
  "detection_rate_as_reads_detected_over_total_reads",
  "sites_detection_rate", "mean_detection_rate_charm", "sd_detection_rate_charm",
  "mean_detection_rate_reads_checked_charm", "sd_detection_rate_reads_checked_charm",
  "mean_detection_rate_total_reads_charm", "sd_detection_rate_total_reads_charm",
  "mean_sites_rate_charm", "sd_sites_rate_charm",
  "detection_rate_zscore_charm","sites_rate_zscore_charm",
  "detection_rate_zscore_reads_checked_charm",
  "detection_rate_zscore_total_reads_charm"
)

df <- df %>%
  group_by(
    Sample_ID, Sample_ID_Bam, VCF, VCF_clean, Study, Patient,
    Date_of_sample_collection, Sample_type, Timepoint, Study_VCF,
    timepoint_info, Sample_type_Bam, timepoint_info_Bam,
    Mut_source, Filter_source, plotting_type, VCF_factor
  ) %>%
  summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 3) Recompute MRD status & cumulative VAF
# Mrd_by_WGS threshold of 4.5 z-scores: sites_rate_zscore_charm measures how many
# standard deviations above the CHARM healthy-control mean this sample's
# site-detection rate falls. A z-score > 4.5 was chosen as a conservative positive
# threshold to maximize specificity at this screening step. The final model-based
# probability thresholds in script 3_1 are more precisely optimized to the
# frontline training data and should be used for clinical calls.
#
# Cumulative_VAF: the mutation detection rate (reads_detected / reads_checked)
# for MRD-positive samples only; set to 0 for negatives. This is the primary
# quantitative cfDNA feature passed to downstream analysis scripts.
df <- df %>%
  mutate(
    Mrd_by_WGS    = if_else(sites_rate_zscore_charm > 4.5, "Positive", "Negative"),
    Cumulative_VAF = if_else(Mrd_by_WGS == "Positive",
                             detection_rate_as_reads_detected_over_reads_checked,
                             0)
  )

## Check for dups 
# Diagnostics: where are the dupes?
cfWGS_metadata %>%
  count(Sample_ID) %>% filter(n > 1) %>% arrange(desc(n)) -> meta_dupes

df %>%
  count(Sample_ID, Sample_ID_Bam) %>% filter(n > 1) %>% arrange(desc(n)) -> df_pair_dupes

message("# meta dupes (by Sample_ID): ", nrow(meta_dupes))
message("# df pair dupes (by Sample_ID, Sample_ID_Bam): ", nrow(df_pair_dupes))

# Make metadata unique per Sample_ID (break ties deterministically)
cc <- cfWGS_metadata %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection)) %>%
  arrange(Sample_ID, desc(Date_of_sample_collection)) %>%
  group_by(Sample_ID) %>%
  # keep only ONE row per Sample_ID (latest date; if ties, take the first after arrange)
  slice_head(n = 1) %>%
  ungroup() %>%
  select(Sample_ID, Date_of_sample_collection, Timepoint)

# sanity check
stopifnot(nrow(cc) == dplyr::n_distinct(cc$Sample_ID))

# Also ensure df has unique pairs before joining dates
df_pairs <- df %>%
  select(Sample_ID, Sample_ID_Bam, Patient) %>%
  distinct()

# 4) Join dates for both IDs
start_dates <- df_pairs %>%
  left_join(cc, by = "Sample_ID") %>%
  rename(Date_of_sample_collection_Sample_ID = Date_of_sample_collection,
         Timepoint_Sample_ID                  = Timepoint) %>%
  left_join(cc, by = c("Sample_ID_Bam" = "Sample_ID")) %>%
  rename(Date_of_sample_collection_Sample_ID_Bam = Date_of_sample_collection,
         Timepoint_Sample_ID_Bam                  = Timepoint) %>%
  mutate(
    num_days  = as.numeric(
      difftime(Date_of_sample_collection_Sample_ID_Bam,
               Date_of_sample_collection_Sample_ID,
               units = "days")
    ),
    num_weeks = num_days / 7
  )

# Final guard: confirm uniqueness using the same key used for the downstream join.
start_date_dupes <- start_dates %>%
  count(Sample_ID, Sample_ID_Bam, Patient) %>%
  filter(n > 1)
if (nrow(start_date_dupes) > 0) {
  readr::write_csv(
    start_date_dupes,
    file.path(outdir, "start_date_duplicate_key_audit.csv")
  )
  stop(
    "Duplicate start-date rows remain for the Sample_ID/Sample_ID_Bam/Patient join key. ",
    "Audit written to: ", file.path(outdir, "start_date_duplicate_key_audit.csv"),
    call. = FALSE
  )
}

# 5) Attach start_dates back onto df
combined_data_plot <- df %>%
  left_join(start_dates, by = c("Sample_ID", "Sample_ID_Bam", "Patient"))

# 6) Flag “Good_baseline_marrow” from your All_feature_data
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_August2025.rds")
cohort_df <- readRDS("cohort_assignment_table_updated.rds")
good_pts <- All_feature_data %>%
  filter(Sample_type == "BM_cells",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>% unique()

combined_data_plot <- combined_data_plot %>%
  mutate(
    Good_baseline_marrow = if_else(Patient %in% good_pts, "Yes", "No")
  )

# 7) Calculate START_DATE and percent_change from baseline
combined_data_plot <- combined_data_plot %>%
  rename(START_DATE = num_days) %>%
  filter(Good_baseline_marrow == "Yes",
         Sample_type_Bam == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  mutate(
    baseline_date = min(START_DATE[timepoint_info %in% c("Baseline","Diagnosis")], na.rm = TRUE),
    baseline_rate = detection_rate_as_reads_detected_over_reads_checked[START_DATE == baseline_date][1],
    absolute_change = detection_rate_as_reads_detected_over_reads_checked - baseline_rate,
    percent_change  = if_else(!is.na(baseline_rate),
                              (absolute_change / baseline_rate) * 100,
                              NA_real_),
    Weeks_since_baseline = pmax(0, START_DATE) / 7
  ) %>%
  ungroup() 

## Get the change since the first treatment timepoint 
# Define the treatment sample categories
treatment_samples <- c("Post_induction", "Post_transplant", "Maintenance", "1.5yr maintenance")

combined_data_plot <- combined_data_plot %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(START_DATE) %>%  # Ensure entries are sorted by days since baseline within each patient
  dplyr::mutate(
    # Find the first treatment sample date if available
    first_treatment_date = ifelse(any(timepoint_info_Bam %in% treatment_samples),
                                  START_DATE[timepoint_info_Bam %in% treatment_samples][1],
                                  NA_real_),
    # Calculate the detection rate for the first treatment timepoint
    treatment_detection_rate = ifelse(!is.na(first_treatment_date),
                                      detection_rate_as_reads_detected_over_reads_checked[START_DATE == first_treatment_date],
                                      NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    # Calculate the absolute and percent change in detection rate since the first treatment timepoint
    absolute_change_detection_rate_treatment = detection_rate_as_reads_detected_over_reads_checked - treatment_detection_rate,
    percent_change_detection_rate_treatment = (absolute_change_detection_rate_treatment / treatment_detection_rate) * 100,
    # Calculate weeks since the first treatment date
    Weeks_since_first_treatment = (START_DATE - first_treatment_date) / 7
  )

# For consistency 
combined_data_plot$Weeks_since_second_start <- combined_data_plot$Weeks_since_first_treatment
combined_data_plot$percent_change_detection_rate_second_timepoint <- combined_data_plot$percent_change_detection_rate_treatment


# 8) Write out 
export_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/BM_muts_plots_baseline/"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

combined_data_plot <- restrict_to_single_baseline_mutation_source(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "bm"
)

bm_source_duplicate_audit <- combined_data_plot %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
  count(Patient, Mut_source, name = "n_baseline_diagnosis_mutation_sources") %>%
  filter(n_baseline_diagnosis_mutation_sources > 1) %>%
  left_join(
    combined_data_plot %>%
      filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
      distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
      group_by(Patient, Mut_source) %>%
      summarise(
        mutation_sources = paste(
          sort(unique(paste(Sample_ID, Timepoint, timepoint_info, sep = "|"))),
          collapse = "; "
        ),
        .groups = "drop"
      ),
    by = c("Patient", "Mut_source")
  )

write_csv(
  bm_source_duplicate_audit,
  file.path(export_dir, "spring2026_bm_mrdetect_duplicate_baseline_source_audit.csv")
)

if (nrow(bm_source_duplicate_audit) > 0L) {
  warning(
    "Patients with more than one baseline/diagnosis BM-mutation source were found. ",
    "Audit written to: ",
    file.path(export_dir, "spring2026_bm_mrdetect_duplicate_baseline_source_audit.csv")
  )
}

readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS_MRDetect_BM_data_updated_Sep.csv")
)



# ──────────────────────────────────────────────────────────────────────────────
# 15) Additional Processing - Blood muts
# ──────────────────────────────────────────────────────────────────────────────
### Next for blood-derived muts
#Merged_MRDetect_zscore <- readRDS(file = "MRDetect_output_winter_2025/Processed_R_outputs/cfWGS_Winter2025All_MRDetect_with_Zscore_May2025.rds")
# 1) Start from your z-scored MRDetect table
df <- Merged_MRDetect_zscore %>%
  filter(plotting_type == "Matched_plasma",
         Mut_source      == "Blood",
         Filter_source   == "STR_encode")

# Rename df's columns to the names num_cols vector expects
rename_map <- c(
  mean_detection_rate_charm                       = "mean_det_charm",
  sd_detection_rate_charm                         = "sd_det_charm",
  mean_detection_rate_reads_checked_charm         = "mean_det_checked_charm",
  sd_detection_rate_reads_checked_charm           = "sd_det_checked_charm",
  mean_detection_rate_total_reads_charm           = "mean_det_total_charm",
  sd_detection_rate_total_reads_charm             = "sd_det_total_charm",
  mean_sites_rate_charm                           = "mean_sites_charm",
  sd_sites_rate_charm                             = "sd_sites_charm"
)

df <- df %>% rename(!!!rename_map)

# 2) Collapse duplicate BAM runs by averaging all numeric QC/MRD metrics
num_cols <- c(
  "sites_checked", "reads_checked", "sites_detected", "reads_detected",
  "total_reads", "detection_rate",
  "detection_rate_as_reads_detected_over_reads_checked",
  "detection_rate_as_reads_detected_over_total_reads",
  "sites_detection_rate", "mean_detection_rate_charm", "sd_detection_rate_charm",
  "mean_detection_rate_reads_checked_charm", "sd_detection_rate_reads_checked_charm",
  "mean_detection_rate_total_reads_charm", "sd_detection_rate_total_reads_charm",
  "mean_sites_rate_charm", "sd_sites_rate_charm",
  "detection_rate_zscore_charm","sites_rate_zscore_charm",
  "detection_rate_zscore_reads_checked_charm",
  "detection_rate_zscore_total_reads_charm"
)

df <- df %>%
  group_by(
    Sample_ID, Sample_ID_Bam, VCF, VCF_clean, Study, Patient,
    Date_of_sample_collection, Sample_type, Timepoint, Study_VCF,
    timepoint_info, Sample_type_Bam, timepoint_info_Bam,
    Mut_source, Filter_source, plotting_type, VCF_factor
  ) %>%
  summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 3) Recompute MRD status & cumulative VAF
df <- df %>%
  mutate(
    Mrd_by_WGS    = if_else(sites_rate_zscore_charm > 4.5, "Positive", "Negative"),
    Cumulative_VAF = if_else(Mrd_by_WGS == "Positive",
                             detection_rate_as_reads_detected_over_reads_checked,
                             0)
  ) ## Old version

# 4) Pull in clinical dates & timepoints for both Sample_ID and Sample_ID_Bam
cc <- cfWGS_metadata %>%
  group_by(Sample_ID) %>%
  slice_max(Date_of_sample_collection, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Sample_ID, Date_of_sample_collection, Timepoint)

# Build start_dates _including_ Patient
start_dates <- df %>%
  select(Sample_ID, Sample_ID_Bam, Patient) %>% 
  left_join(cc,                    by = "Sample_ID")           %>% 
  rename(Date_of_sample_collection_Sample_ID = Date_of_sample_collection,
         Timepoint_Sample_ID            = Timepoint)           %>% 
  left_join(cc, by = c("Sample_ID_Bam" = "Sample_ID"))         %>% 
  rename(Date_of_sample_collection_Sample_ID_Bam = Date_of_sample_collection,
         Timepoint_Sample_ID_Bam            = Timepoint)       %>% 
  mutate(
    num_days  = as.numeric(difftime(
      Date_of_sample_collection_Sample_ID_Bam,
      Date_of_sample_collection_Sample_ID,
      units = "days")),
    num_weeks = num_days / 7
  ) %>% 
  distinct()


# 5) Attach start_dates back onto df
combined_data_plot <- df %>%
  left_join(start_dates, by = c("Sample_ID", "Sample_ID_Bam", "Patient"))

# 6) Flag “Good_baseline_sample” from your All_feature_data
good_pts <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>% unique()

combined_data_plot <- combined_data_plot %>%
  mutate(
    Good_baseline_sample = if_else(Patient %in% good_pts, "Yes", "No")
  )

filtered_good_pts <- tibble(Patient = good_pts) %>%
  inner_join(cohort_df, by = "Patient")

filtered_good_pts %>%
  pull(Patient) %>% unique()

# 7) Calculate START_DATE and percent_change from baseline
combined_data_plot <- combined_data_plot %>%
  rename(START_DATE = num_days) %>%
  filter(Good_baseline_sample == "Yes",
         Sample_type_Bam == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  mutate(
    baseline_date = min(START_DATE[timepoint_info %in% c("Baseline","Diagnosis")], na.rm = TRUE),
    baseline_rate = detection_rate_as_reads_detected_over_reads_checked[START_DATE == baseline_date][1],
    absolute_change = detection_rate_as_reads_detected_over_reads_checked - baseline_rate,
    percent_change  = if_else(!is.na(baseline_rate),
                              (absolute_change / baseline_rate) * 100,
                              NA_real_),
    Weeks_since_baseline = pmax(0, START_DATE) / 7
  ) %>%
  ungroup() 


## Now edit this to show the time difference since the second timepoint 
## First calculate the percent change in detection rate and put dates

## Get the change since the first treatment timepoint 
# Define the treatment sample categories
treatment_samples <- c("Post_induction", "Post_transplant", "Maintenance", "1.5yr maintenance")

combined_data_plot <- combined_data_plot %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(START_DATE) %>%  # Ensure entries are sorted by days since baseline within each patient
  dplyr::mutate(
    # Find the first treatment sample date if available
    first_treatment_date = ifelse(any(timepoint_info_Bam %in% treatment_samples),
                                  START_DATE[timepoint_info_Bam %in% treatment_samples][1],
                                  NA_real_),
    # Calculate the detection rate for the first treatment timepoint
    treatment_detection_rate = ifelse(!is.na(first_treatment_date),
                                      detection_rate_as_reads_detected_over_reads_checked[START_DATE == first_treatment_date],
                                      NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    # Calculate the absolute and percent change in detection rate since the first treatment timepoint
    absolute_change_detection_rate_treatment = detection_rate_as_reads_detected_over_reads_checked - treatment_detection_rate,
    percent_change_detection_rate_treatment = (absolute_change_detection_rate_treatment / treatment_detection_rate) * 100,
    # Calculate weeks since the first treatment date
    Weeks_since_first_treatment = (START_DATE - first_treatment_date) / 7
  )

# For consistency 
combined_data_plot$Weeks_since_second_start <- combined_data_plot$Weeks_since_first_treatment
combined_data_plot$percent_change_detection_rate_second_timepoint <- combined_data_plot$percent_change_detection_rate_treatment

## Get rid of marrows that are not at baseline 
combined_data_plot <- combined_data_plot %>%
  filter(Timepoint %in% c("01", "T0", "T1", "1"))
combined_data_plot <- combined_data_plot %>% 
  filter(Sample_type == "Blood_plasma_cfDNA")


# 8) Write out 
export_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

combined_data_plot <- restrict_to_single_baseline_mutation_source(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_good_baseline"
)

readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS MRDetect Blood data updated Sep.csv")
)


## Redo for all bloods available
# 5) Attach start_dates back onto df
combined_data_plot <- df %>%
  left_join(start_dates, by = c("Sample_ID", "Sample_ID_Bam", "Patient"))

# 6) Flag “Good_baseline_sample” from your All_feature_data
good_pts <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>% unique()

combined_data_plot <- combined_data_plot %>%
  mutate(
    Good_baseline_sample = if_else(Patient %in% good_pts, "Yes", "No")
  )

# 7) Calculate START_DATE and percent_change from baseline
combined_data_plot <- combined_data_plot %>%
  rename(START_DATE = num_days) %>%
  filter(#Good_baseline_sample == "Yes",
         Sample_type_Bam == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  mutate(
    baseline_date = min(START_DATE[timepoint_info %in% c("Baseline","Diagnosis")], na.rm = TRUE),
    baseline_rate = detection_rate_as_reads_detected_over_reads_checked[START_DATE == baseline_date][1],
    absolute_change = detection_rate_as_reads_detected_over_reads_checked - baseline_rate,
    percent_change  = if_else(!is.na(baseline_rate),
                              (absolute_change / baseline_rate) * 100,
                              NA_real_),
    Weeks_since_baseline = pmax(0, START_DATE) / 7
  ) %>%
  ungroup() 


## Now edit this to show the time difference since the second timepoint 
## First calculate the percent change in detection rate and put dates

## Get the change since the first treatment timepoint 
# Define the treatment sample categories
treatment_samples <- c("Post_induction", "Post_transplant", "Maintenance", "1.5yr maintenance")

combined_data_plot <- combined_data_plot %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(START_DATE) %>%  # Ensure entries are sorted by days since baseline within each patient
  dplyr::mutate(
    # Find the first treatment sample date if available
    first_treatment_date = ifelse(any(timepoint_info_Bam %in% treatment_samples),
                                  START_DATE[timepoint_info_Bam %in% treatment_samples][1],
                                  NA_real_),
    # Calculate the detection rate for the first treatment timepoint
    treatment_detection_rate = ifelse(!is.na(first_treatment_date),
                                      detection_rate_as_reads_detected_over_reads_checked[START_DATE == first_treatment_date],
                                      NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    # Calculate the absolute and percent change in detection rate since the first treatment timepoint
    absolute_change_detection_rate_treatment = detection_rate_as_reads_detected_over_reads_checked - treatment_detection_rate,
    percent_change_detection_rate_treatment = (absolute_change_detection_rate_treatment / treatment_detection_rate) * 100,
    # Calculate weeks since the first treatment date
    Weeks_since_first_treatment = (START_DATE - first_treatment_date) / 7
  )

# For consistency 
combined_data_plot$Weeks_since_second_start <- combined_data_plot$Weeks_since_first_treatment
combined_data_plot$percent_change_detection_rate_second_timepoint <- combined_data_plot$percent_change_detection_rate_treatment

## Get rid of marrows that are not at baseline 
combined_data_plot <- combined_data_plot %>%
  filter(Timepoint %in% c("01", "T0", "T1", "1"))
combined_data_plot <- combined_data_plot %>% 
  filter(Sample_type == "Blood_plasma_cfDNA")


# 8) Write out 
export_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

combined_data_plot <- restrict_to_single_baseline_mutation_source(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_all_patients"
)

blood_source_duplicate_audit <- combined_data_plot %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
  count(Patient, Mut_source, name = "n_baseline_diagnosis_mutation_sources") %>%
  filter(n_baseline_diagnosis_mutation_sources > 1) %>%
  left_join(
    combined_data_plot %>%
      filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
      distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
      group_by(Patient, Mut_source) %>%
      summarise(
        mutation_sources = paste(
          sort(unique(paste(Sample_ID, Timepoint, timepoint_info, sep = "|"))),
          collapse = "; "
        ),
        .groups = "drop"
      ),
    by = c("Patient", "Mut_source")
  )

write_csv(
  blood_source_duplicate_audit,
  file.path(export_dir, "spring2026_blood_mrdetect_duplicate_baseline_source_audit.csv")
)

if (nrow(blood_source_duplicate_audit) > 0L) {
  warning(
    "Patients with more than one baseline/diagnosis blood-mutation source were found. ",
    "Audit written to: ",
    file.path(export_dir, "spring2026_blood_mrdetect_duplicate_baseline_source_audit.csv")
  )
}

readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS MRDetect Blood data updated Sep with all patients.csv")
)


## Script complete.
