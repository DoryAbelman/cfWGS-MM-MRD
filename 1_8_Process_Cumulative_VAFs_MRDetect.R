# =============================================================================
# MRDetect_processing.R
# Project:  cfWGS MRDetect (Winter 2025)
# Author:   Dory Abelman
# Date:     January 2025
# Last Updated: May 2025
#
# Purpose:
#   1. Read all MRDetect CSV outputs.
#   2. Annotate each record with source file, sample metadata, and z-scores
#      based on CHARM_healthy controls.
#   3. Filter to cfDNA timepoints (Diagnosis/Baseline/Progression) and export
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
#
# Usage:
#   Rscript MRDetect_processing.R
# =============================================================================

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
library(conflicted)

# resolve common conflicts
conflicted::conflicts_prefer("dplyr::mutate")
conflicted::conflicts_prefer("dplyr::filter")
conflicted::conflicts_prefer("dplyr::select")
conflicted::conflicts_prefer("dplyr::summarize")


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

read_and_label <- function(file) {
  df <- read_csv(file)
  df$source_file <- basename(file)
  return(df)
}
all_files <- bind_rows(lapply(csv_files, read_and_label))

all_files <- all_files %>%
  mutate(
    Mut_source = case_when(
      str_detect(source_file, "BM_muts")    ~ "BM_cells",
      str_detect(source_file, "Blood_muts") ~ "Blood",
      TRUE                                   ~ NA_character_
    ),
    Filter_source = case_when(
      str_detect(source_file, "encode_only") ~ "Encode_only",
      str_detect(source_file, "STR_encode")  ~ "STR_encode",
      TRUE                                   ~ NA_character_
    )
  ) %>%
  # drop ‘VCF’ column if present
  select(-any_of("VCF")) %>%
  # clean up ‘filename’ field
  mutate(
    filename   = str_remove(filename, "^\\./"),
    VCF        = sub(".*_VS_", "\\1", filename),
    VCF_clean  = filename %>%
      str_remove("\\.mutect2.*") %>%
      str_remove("\\.fil.*")
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
cfWGS_metadata <- read_csv("../combined_clinical_data_updated_April2025.csv") %>%
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


# ──────────────────────────────────────────────────────────────────────────────
# 5) OPTIONAL: Swap mis‐labeled columns (sites_detected ↔ reads_detected ↔ total_reads)
#    (only needed if you used an older parser version that mis-assigned these)
# ──────────────────────────────────────────────────────────────────────────────

# Temporarily store each column
temp_sites <- all_files$sites_detected
temp_reads <- all_files$reads_detected
temp_total <- all_files$total_reads

# Reassign them in the correct order
all_files <- all_files %>%
  mutate(
    sites_detected = temp_reads,
    reads_detected = temp_total,
    total_reads    = temp_sites
  )

# Remove the temporary objects
rm(temp_sites, temp_reads, temp_total)



# ──────────────────────────────────────────────────────────────────────────────
# 6) Annotate each row with ‘Study’ based on BAM prefix
# ──────────────────────────────────────────────────────────────────────────────
all_files <- all_files %>%
  mutate(
    Study = case_when(
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
# 7) Standardize column names: replace spaces & dots with underscores
# ──────────────────────────────────────────────────────────────────────────────
colnames(all_files) <- colnames(all_files) %>%
  str_trim(end = "right") %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("\\.", "_")


# ──────────────────────────────────────────────────────────────────────────────
# 8) Merge MRDetect outputs with clinical metadata
# ──────────────────────────────────────────────────────────────────────────────
tmp_meta <- cfWGS_metadata %>%
  select(-Bam) %>%
  rename(Study_VCF = Study)

Merged_MRDetect <- all_files %>%
  left_join(tmp_meta, by = c("VCF_clean" = "VCF_clean_merge"))


# ──────────────────────────────────────────────────────────────────────────────
# 9) Add BAM‐level sample info (Sample_ID, Patient, Sample_type, timepoint_info)
# ──────────────────────────────────────────────────────────────────────────────
bam_info <- cfWGS_metadata %>%
  select(Bam, Sample_ID, Patient, Sample_type, timepoint_info) %>%
  rename_with(~ paste0(.x, "_Bam"), .cols = -Bam)

Merged_MRDetect <- Merged_MRDetect %>%
  left_join(bam_info, by = c("BAM" = "Bam"))


# ──────────────────────────────────────────────────────────────────────────────
# 10) Classify matched vs. unmatched plasma; force CHARM_healthy→cfDNA
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect <- Merged_MRDetect %>%
  mutate(
    plotting_type       = if_else(Patient_Bam == Patient, "Matched_plasma", "Unmatched_plasma"),
    Sample_type_Bam     = if_else(Study == "CHARM_healthy", "Blood_plasma_cfDNA", Sample_type_Bam)
  )


# ──────────────────────────────────────────────────────────────────────────────
# 11) Subset to cfDNA timepoints only (Baseline/Diagnosis/Progression)
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect_cfDNA <- Merged_MRDetect %>%
  filter(Sample_type_Bam == "Blood_plasma_cfDNA",
         timepoint_info_Bam %in% c("Baseline","Diagnosis","Progression"))


# ──────────────────────────────────────────────────────────────────────────────
# 12) Compute CHARM_healthy means & SDs → join back for z-scores
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect$VCF_factor <- factor(Merged_MRDetect$VCF, levels = unique(Merged_MRDetect$VCF))

zscore_lookup <- Merged_MRDetect %>%
  filter(Study == "CHARM_healthy") %>%
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
    sites_rate_zscore_charm                       = (sites_detection_rate - mean_sites_charm) / sd_sites_charm
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
base_name_raw    <- paste0(project, "All_MRDetect_May2025")
base_name_zscore <- paste0(project, "All_MRDetect_with_Zscore_May2025")

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
