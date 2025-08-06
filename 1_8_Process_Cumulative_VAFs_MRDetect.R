# =============================================================================
# 1_8_Process_Cumulative_VAFs_MRDetect.R
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
#       - cfWGS MRDetect BM data updated May.csv
#       - cfWGS MRDetect Blood data updated May.csv
#
# Usage:
#   Rscript 1_8_Process_Cumulative_VAFs_MRDetect.R
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
#library(conflicted)
library(lubridate)

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
cfWGS_metadata <- read_csv("combined_clinical_data_updated_April2025.csv") %>%
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
# 5) Standardize column names: replace spaces & dots with underscores
# ──────────────────────────────────────────────────────────────────────────────
# 1) remove trailing spaces
colnames(all_files) <- gsub("\\s+$", "", colnames(all_files))
# 2) replace any runs of spaces with "_"
colnames(all_files) <- gsub("\\s+", "_", colnames(all_files))
# 3) replace "." with "_"
colnames(all_files) <- gsub("\\.", "_", colnames(all_files))


# ──────────────────────────────────────────────────────────────────────────────
# 6) OPTIONAL: Swap mis‐labeled columns (sites_detected ↔ reads_detected ↔ total_reads)
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
# 7) Annotate each row with ‘Study’ based on BAM prefix
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
  summarise(across(all_of(num_cols), mean, na.rm = TRUE), .groups = "drop")

# 3) Recompute MRD status & cumulative VAF
df <- df %>%
  mutate(
    Mrd_by_WGS    = if_else(sites_rate_zscore_charm > 4.5, "Positive", "Negative"),
    Cumulative_VAF = if_else(Mrd_by_WGS == "Positive",
                             detection_rate_as_reads_detected_over_reads_checked,
                             0)
  )

# 4) Pull in clinical dates & timepoints for both Sample_ID and Sample_ID_Bam
cc <- cfWGS_metadata %>%
  group_by(Sample_ID) %>%
  slice_max(Date_of_sample_collection, n = 1) %>%
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

# 6) Flag “Good_baseline_marrow” from your All_feature_data
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_May2025.rds")
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

readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS_MRDetect_BM_data_updated_May.csv")
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
  summarise(across(all_of(num_cols), mean, na.rm = TRUE), .groups = "drop")

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
  slice_max(Date_of_sample_collection, n = 1) %>%
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
  filter(Timepoint %in% c("01", "T0", "1"))
combined_data_plot <- combined_data_plot %>% 
  filter(Sample_type == "Blood_plasma_cfDNA")


# 8) Write out 
export_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/"
readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS MRDetect Blood data updated June.csv")
)


## Script complete.
