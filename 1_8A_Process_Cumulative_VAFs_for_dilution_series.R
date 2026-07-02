# =============================================================================
# 1_8A_Process_Cumulative_VAFs_for_dilution_series.R
# Project:  cfWGS MRDetect (Winter 2025)
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_8A_Process_Cumulative_VAFs_for_dilution_series.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script processes MRDetect outputs for
#     dilution-series samples used by 3_1_part2 limit-of-detection analyses.
#
# Author:   Dory Abelman
# Date:     January 2025
# Last Updated: May 2025
#
# Purpose:
#   Identical MRDetect processing pipeline to 1_8_Process_Cumulative_VAFs_MRDetect.R
#   but applied exclusively to the EXPERIMENTAL DILUTION SERIES samples rather than
#   the main patient cohort. Reads MRDetect CSV outputs from the dilution series
#   input directory, annotates records with z-scores relative to CHARM healthy
#   controls, and exports processed tables for use in the LOD (limit-of-detection)
#   analysis (script 3_1_part2).
#
# Dependencies:
#   • readr, data.table, tidyverse (dplyr, tidyr, stringr), openxlsx
#   • ggplot2, ggbreak, patchwork, scales, conflicted
#
# Input Files:
#   • MRDetect_output_winter_2025/MRDetect_outputs/Dilution_series/*.csv
#   • combined_clinical_data_updated_April2025.csv
#
# Output Directory (created if necessary):
#   • MRDetect_output_winter_2025/Processed_R_outputs/
#   • Writes:
#       - cfWGS_Winter2025Dilution_series_May2025.rds
#       - cfWGS_Winter2025Dilution_series_May2025_with_zscore.rds
#
# Usage:
#   Rscript 1_8A_Process_Cumulative_VAFs_for_dilution_series.R
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
library(conflicted)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

# resolve common conflicts
conflicted::conflicts_prefer(
  dplyr::mutate,
  dplyr::filter,
  dplyr::select,
  dplyr::summarise,
  dplyr::summarize,
  .quiet = TRUE
)


# ──────────────────────────────────────────────────────────────────────────────
# 2) Set up paths and create output directory
# ──────────────────────────────────────────────────────────────────────────────
input_root <- "MRDetect_output_winter_2025/MRDetect_outputs/Dilution_series/"
outdir     <- "MRDetect_output_winter_2025/Processed_R_outputs/"
project    <- "cfWGS_Winter2025"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}


# ──────────────────────────────────────────────────────────────────────────────
# 3) Read and combine all CSVs, label with ‘source_file’, ‘Mut_source’ & ‘Filter_source’
# ──────────────────────────────────────────────────────────────────────────────
csv_files <- list.files(input_root, pattern = "\\.csv$", full.names = TRUE)
spring2026_dilution_mrdetect_files <- spring2026_revision_files(
  "MRDetect_outputs",
  "^(MRDetect_all_RESULTS_combined_with_source_dilution_series|MRDetect_dilution_series_updated)[.]csv$"
)
spring2026_combined_mrdetect_files <- spring2026_revision_files(
  "MRDetect_outputs",
  "^MRDetect_all_RESULTS_combined_with_source[.]csv$"
)
csv_files <- unique(c(
  csv_files,
  spring2026_dilution_mrdetect_files,
  spring2026_combined_mrdetect_files
))
if (!length(csv_files)) {
  stop("No dilution-series MRDetect CSV files found in historical or Spring 2026 revision inputs.", call. = FALSE)
}

read_mrdetect_csv <- function(file) {
  header_line <- readLines(file, n = 1, warn = FALSE)
  header_cols <- trimws(strsplit(header_line, ",", fixed = TRUE)[[1]])
  probe_lines <- readLines(file, n = 25, warn = FALSE)
  data_probe <- probe_lines[-1]
  nonempty_data_probe <- data_probe[nzchar(data_probe)]
  field_counts <- count.fields(file, sep = ",", quote = "", blank.lines.skip = TRUE)
  header_field_count <- field_counts[1]
  max_data_field_count <- max(field_counts[-1], na.rm = TRUE)
  has_one_extra_data_field <- is.finite(max_data_field_count) &&
    max_data_field_count == header_field_count + 1
  has_terminal_empty_data_field <- length(nonempty_data_probe) > 0 &&
    !endsWith(header_line, ",") &&
    any(endsWith(nonempty_data_probe, ","))
  has_legacy_empty_field_before_filename <- has_one_extra_data_field &&
    tail(header_cols, 1) == "filename" &&
    !has_terminal_empty_data_field

  if (!has_terminal_empty_data_field && !has_legacy_empty_field_before_filename) {
    return(read_csv(file, show_col_types = FALSE))
  }

  col_names <- if (has_legacy_empty_field_before_filename) {
    c(head(header_cols, -1), "legacy_empty_field", tail(header_cols, 1))
  } else {
    c(header_cols, "terminal_empty_field")
  }
  df <- read_csv(
    file,
    col_names = col_names,
    skip = 1,
    show_col_types = FALSE
  )
  empty_field <- if (has_legacy_empty_field_before_filename) {
    "legacy_empty_field"
  } else {
    "terminal_empty_field"
  }
  empty_values <- as.character(df[[empty_field]])
  empty_values[is.na(empty_values)] <- ""
  if (any(nzchar(empty_values))) {
    stop(
      "MRDetect CSV has an unexpected non-empty extra field: ",
      file,
      call. = FALSE
    )
  }
  df %>% select(-all_of(empty_field))
}

read_and_label <- function(file) {
  df <- read_mrdetect_csv(file)
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

pwgval_dilution_metadata_for_mrdetect <- load_spring2026_pwgval_dilution_metadata(required = FALSE)
pwgval_dilution_bams <- if (is.null(pwgval_dilution_metadata_for_mrdetect)) {
  character()
} else {
  unique(pwgval_dilution_metadata_for_mrdetect$BAM)
}

if (length(spring2026_combined_mrdetect_files) > 0 && length(pwgval_dilution_bams) > 0) {
  spring2026_combined_basenames <- basename(spring2026_combined_mrdetect_files)
  all_files <- all_files %>%
    filter(
      !.data$input_source_file %in% spring2026_combined_basenames |
        .data$BAM %in% pwgval_dilution_bams |
        str_detect(.data$BAM, "^M4CHIP_") |
        str_detect(.data$source_file, "^M4CHIP_") |
        str_detect(.data$filename, "^M4CHIP_")
    )
}

all_files <- all_files %>%
  mutate(
    Mut_source = case_when(
      str_detect(input_source_file, "BM_muts") | str_detect(source_file, "BM_muts") ~ "BM_cells",
      str_detect(input_source_file, "Blood_muts") | str_detect(source_file, "Blood_muts") ~ "Blood",
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
      str_remove("\\.fil.*") %>%
      str_remove("\\.somatic.*")
  )

extract_mrdetect_mutation_list_patient <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "VA-[0-9]+") ~ stringr::str_extract(x, "VA-[0-9]+"),
    stringr::str_detect(x, "IMG-[0-9]+") ~ stringr::str_extract(x, "IMG-[0-9]+"),
    TRUE ~ NA_character_
  )
}

extract_mrdetect_mutation_list_timepoint <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "VA-[0-9]+-[0-9]+-[A-Z]-DNA") ~
      stringr::str_match(x, "VA-[0-9]+-([0-9]+)-[A-Z]-DNA")[, 2],
    stringr::str_detect(x, "VA-[0-9]+-R-[A-Z]-DNA") ~ "R",
    stringr::str_detect(x, "IMG-[0-9]+-T[0-9]+-[A-Z]") ~
      stringr::str_match(x, "IMG-[0-9]+-(T[0-9]+)-[A-Z]")[, 2],
    TRUE ~ NA_character_
  )
}

is_baseline_or_diagnosis_mutation_list <- function(timepoint) {
  timepoint %in% c("01", "T0", "T1")
}

if (length(pwgval_dilution_bams) > 0) {
  pwgval_bam_patient_lookup <- pwgval_dilution_metadata_for_mrdetect %>%
    transmute(
      BAM,
      PWGVAL_dilution_patient = .data$Patient,
      PWGVAL_dilution_sample_id = .data$Sample_ID,
      PWGVAL_dilution_LOD = .data$LOD
    ) %>%
    distinct(BAM, .keep_all = TRUE)

  all_files <- all_files %>%
    mutate(
      is_pwgval_m4chip_query_bam = .data$BAM %in% pwgval_dilution_bams | str_detect(.data$BAM, "^M4CHIP_"),
      VCF_mutation_list_patient_for_dilution = extract_mrdetect_mutation_list_patient(.data$VCF),
      VCF_mutation_list_timepoint_for_dilution = extract_mrdetect_mutation_list_timepoint(.data$VCF),
      VCF_mutation_list_is_baseline_or_diagnosis =
        is_baseline_or_diagnosis_mutation_list(.data$VCF_mutation_list_timepoint_for_dilution),
      VCF_mutation_list_matches_PWGVAL_patient = NA,
      VCF_mutation_list_usable_for_PWGVAL = NA
    ) %>%
    left_join(pwgval_bam_patient_lookup, by = "BAM") %>%
    mutate(
      VCF_mutation_list_matches_PWGVAL_patient = if_else(
        .data$is_pwgval_m4chip_query_bam,
        !is.na(.data$PWGVAL_dilution_patient) &
          !is.na(.data$VCF_mutation_list_patient_for_dilution) &
          .data$VCF_mutation_list_patient_for_dilution == .data$PWGVAL_dilution_patient,
        NA
      ),
      VCF_mutation_list_usable_for_PWGVAL = if_else(
        .data$is_pwgval_m4chip_query_bam,
        .data$VCF_mutation_list_matches_PWGVAL_patient %in% TRUE &
          .data$VCF_mutation_list_is_baseline_or_diagnosis %in% TRUE,
        NA
      ),
      PWGVAL_mutation_list_match_audit_note = case_when(
        !.data$is_pwgval_m4chip_query_bam ~ NA_character_,
        is.na(.data$PWGVAL_dilution_patient) ~
          "M4CHIP/PWGVAL-like queried BAM is absent from PWGVAL dilution metadata",
        is.na(.data$VCF_mutation_list_patient_for_dilution) ~
          "Could not infer personalized mutation-list patient from MRDetect source filename",
        !.data$VCF_mutation_list_matches_PWGVAL_patient ~
          "Excluded from PWGVAL dilution scoring because queried BAM patient and personalized mutation-list patient differ",
        !.data$VCF_mutation_list_is_baseline_or_diagnosis ~
          "Excluded from PWGVAL dilution scoring because the personalized mutation list is not baseline/diagnosis",
        .data$VCF_mutation_list_usable_for_PWGVAL ~
          "Kept: PWGVAL queried BAM matched to same-patient baseline/diagnosis personalized mutation list",
        TRUE ~ "Excluded from PWGVAL dilution scoring"
      )
    )

  pwgval_mutation_list_match_audit <- all_files %>%
    filter(.data$is_pwgval_m4chip_query_bam) %>%
    select(
      input_source_file,
      source_file,
      BAM,
      VCF,
      VCF_clean,
      PWGVAL_dilution_patient,
      PWGVAL_dilution_sample_id,
      PWGVAL_dilution_LOD,
      VCF_mutation_list_patient_for_dilution,
      VCF_mutation_list_timepoint_for_dilution,
      VCF_mutation_list_is_baseline_or_diagnosis,
      VCF_mutation_list_matches_PWGVAL_patient,
      VCF_mutation_list_usable_for_PWGVAL,
      PWGVAL_mutation_list_match_audit_note
    ) %>%
    arrange(BAM, VCF)

  if (nrow(pwgval_mutation_list_match_audit) > 0) {
    readr::write_csv(
      pwgval_mutation_list_match_audit,
      file.path(outdir, "spring2026_pwgval_mrdetect_mutation_list_match_audit.csv")
    )
    readr::write_csv(
      pwgval_mutation_list_match_audit %>%
        filter(.data$VCF_mutation_list_usable_for_PWGVAL %in% TRUE),
      file.path(outdir, "spring2026_pwgval_mrdetect_kept_baseline_diagnosis_mutation_lists.csv")
    )
  }

  all_files <- all_files %>%
    filter(
      !.data$is_pwgval_m4chip_query_bam |
        .data$VCF_mutation_list_usable_for_PWGVAL %in% TRUE
    )
}

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
    left_join(spring2026_panel_metadata, by = "VCF_clean")
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
# The dilution-series MRDetect CSV files used for the manuscript were generated
# by the same parser version as the main-sample files. That parser assigned
# three count columns in the wrong order, so the manuscript analysis corrects
# that ordering before calculating detection-rate features. Keep this TRUE for
# the preserved manuscript input set. If future MRDetect CSVs are produced by a
# corrected parser, set this to FALSE only after confirming the input column
# definitions against parser documentation.
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
# ──────────────────────────────────────────────────────────────────────────────
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
      str_detect(BAM, "^M4CHIP_")    ~ "M4",
      TRUE                            ~ NA_character_
    )
  )

## Add the values if missing from the fields 
# If Mut_source is empty, assign based on VCF_clean patterns;
# then set Filter_source to "STR_encode" for all rows

all_files <- all_files %>%
  mutate(
    Mut_source = ifelse(
      is.na(Mut_source) | Mut_source == "",
      case_when(
        grepl("-O-DNA|_Bm_", VCF_clean) ~ "BM_cells",
        grepl("-P-DNA|_Pl_", VCF_clean) ~ "Blood",
        TRUE                        ~ Mut_source
      ),
      Mut_source
    ),
    Filter_source = "STR_encode" # This is what was run prior to loading
  )


# ──────────────────────────────────────────────────────────────────────────────
# 8) Merge MRDetect outputs with clinical metadata
# ──────────────────────────────────────────────────────────────────────────────
tmp_meta <- cfWGS_metadata %>%
  select(-Bam) %>%
  rename(Study_VCF = Study)

Merged_MRDetect_dilution <- all_files %>%
  left_join(tmp_meta, by = c("VCF_clean" = "VCF_clean_merge"))

if ("mutect2_pair_id" %in% names(tmp_meta)) {
  tmp_meta_pair <- tmp_meta %>%
    filter(!is.na(mutect2_pair_id), nzchar(mutect2_pair_id)) %>%
    distinct(mutect2_pair_id, .keep_all = TRUE)

  pair_joined <- all_files %>%
    left_join(tmp_meta_pair, by = c("VCF_clean" = "mutect2_pair_id"))

  fill_cols <- intersect(names(Merged_MRDetect_dilution), names(pair_joined))
  metadata_cols <- setdiff(fill_cols, names(all_files))
  for (col in metadata_cols) {
    Merged_MRDetect_dilution[[col]] <- dplyr::coalesce(
      Merged_MRDetect_dilution[[col]],
      pair_joined[[col]]
    )
  }
}

if ("VCF_panel_sample_type" %in% names(Merged_MRDetect_dilution)) {
  Merged_MRDetect_dilution <- Merged_MRDetect_dilution %>%
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

  excluded_disallowed_dilution_panel <- Merged_MRDetect_dilution %>%
    filter(
      !is.na(VCF_panel_timepoint_info),
      !VCF_panel_is_allowed_panel_timepoint
    )
  if (nrow(excluded_disallowed_dilution_panel) > 0) {
    readr::write_csv(
      excluded_disallowed_dilution_panel,
      file.path(outdir, "spring2026_pwgval_mrdetect_excluded_disallowed_vcf_mutation_lists.csv")
    )
  }

  Merged_MRDetect_dilution <- Merged_MRDetect_dilution %>%
    filter(VCF_panel_is_allowed_panel_timepoint) %>%
    select(-VCF_panel_is_allowed_panel_timepoint)
}

# ──────────────────────────────────────────────────────────────────────────────
# 9) Add BAM‐level sample info (Sample_ID, Patient, Sample_type, timepoint_info)
# ──────────────────────────────────────────────────────────────────────────────
bam_info <- read_dilution_metadata_with_spring2026("Metadata_dilution_series.csv")

# 3. Join back into the main Merged_MRDetect_dilution table
Merged_MRDetect_dilution <- Merged_MRDetect_dilution %>%
  left_join(bam_info, by = "BAM")

if ("Patient.x" %in% names(Merged_MRDetect_dilution)) {
  Merged_MRDetect_dilution$Patient <- Merged_MRDetect_dilution$Patient.x
  Merged_MRDetect_dilution$Patient.x <- NULL
}
if ("Patient.y" %in% names(Merged_MRDetect_dilution)) {
  Merged_MRDetect_dilution$Patient_Bam <- dplyr::coalesce(
    Merged_MRDetect_dilution$Patient_Bam,
    Merged_MRDetect_dilution$Patient.y
  )
  Merged_MRDetect_dilution$Patient.y <- NULL
}

pwgval_mrdetect_availability <- bam_info %>%
  filter(.data$dilution_series == "PWGVAL_M4CHIP") %>%
  distinct(BAM, Patient_Bam, LOD, Sample_ID, Merge) %>%
  mutate(
    present_in_mrdetect_results = .data$BAM %in% unique(all_files$BAM),
    audit_note = if_else(
      .data$present_in_mrdetect_results,
      "PWGVAL dilution BAM present in loaded MRDetect outputs",
      "PWGVAL dilution BAM absent from loaded MRDetect outputs"
    )
  )
if (nrow(pwgval_mrdetect_availability) > 0) {
  readr::write_csv(
    pwgval_mrdetect_availability,
    file.path(outdir, "spring2026_pwgval_mrdetect_dilution_availability_audit.csv")
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# 10) Classify matched vs. unmatched plasma; force CHARM_healthy→cfDNA
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect_dilution <- Merged_MRDetect_dilution %>%
  mutate(
    plotting_type       = if_else(Patient_Bam == Patient, "Matched_plasma", "Unmatched_plasma"),
    Sample_type_Bam     = if_else(Study == "CHARM_healthy", "Blood_plasma_cfDNA", Sample_type_Bam)
  )


# ──────────────────────────────────────────────────────────────────────────────
# 11) Laod in the matched healthy control data and join
# ──────────────────────────────────────────────────────────────────────────────
Healthy_reference <- readRDS("MRDetect_output_winter_2025/Processed_R_outputs/cfWGS_Winter2025All_MRDetect_May2025.rds")
needed_healthy_reference_vcfs <- Merged_MRDetect_dilution %>%
  filter(!is.na(.data$VCF_clean), nzchar(.data$VCF_clean)) %>%
  distinct(VCF_clean)

Healthy_reference <- Healthy_reference %>%
  filter(
    .data$Study == "CHARM_healthy",
    .data$VCF_clean %in% needed_healthy_reference_vcfs$VCF_clean
  )

healthy_reference_audit <- needed_healthy_reference_vcfs %>%
  left_join(
    Healthy_reference %>%
      count(VCF_clean, Mut_source, Filter_source, name = "healthy_reference_rows"),
    by = "VCF_clean"
  ) %>%
  mutate(
    healthy_reference_rows = coalesce(.data$healthy_reference_rows, 0L),
    healthy_reference_status = if_else(
      .data$healthy_reference_rows > 0,
      "available",
      "missing"
    )
  ) %>%
  arrange(.data$VCF_clean, .data$Mut_source, .data$Filter_source)
readr::write_csv(
  healthy_reference_audit,
  file.path(outdir, "spring2026_dilution_mrdetect_healthy_reference_audit.csv")
)

## Ensure matching column types
Merged_MRDetect_dilution <- Merged_MRDetect_dilution %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

Healthy_reference <- Healthy_reference %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

Merged_MRDetect_dilution_joined <- bind_rows(Merged_MRDetect_dilution, Healthy_reference)
# ──────────────────────────────────────────────────────────────────────────────
# 12) Compute CHARM_healthy means & SDs → join back for z-scores
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect_dilution_joined$VCF_factor <- factor(Merged_MRDetect_dilution_joined$VCF, levels = unique(Merged_MRDetect_dilution_joined$VCF))

zscore_lookup <- Merged_MRDetect_dilution_joined %>%
  filter(Study == "CHARM_healthy") %>%
  dplyr::select(VCF_factor, Mut_source, Filter_source,
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

Merged_MRDetect_dilution_zscore <- Merged_MRDetect_dilution_joined %>%
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
Merged_MRDetect_dilution_zscore <- Merged_MRDetect_dilution_zscore %>%
  group_by(across(all_of(dup_keys))) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup()


# ──────────────────────────────────────────────────────────────────────────────
# 14) Save both “raw” and “z-scored” MRDetect tables
# ──────────────────────────────────────────────────────────────────────────────
base_name_raw    <- paste0(project, "Dilution_series_May2025")
base_name_zscore <- paste0(project, "Dilution_series_May2025_with_zscore")

write.table(Merged_MRDetect_dilution,
            file = file.path(outdir, paste0(base_name_raw, ".txt")),
            sep = "\t", row.names = FALSE)

saveRDS(Merged_MRDetect_dilution,
        file = file.path(outdir, paste0(base_name_raw, ".rds")))

write.table(Merged_MRDetect_dilution_zscore,
            file = file.path(outdir, paste0(base_name_zscore, ".txt")),
            sep = "\t", row.names = FALSE)

saveRDS(Merged_MRDetect_dilution_zscore,
        file = file.path(outdir, paste0(base_name_zscore, ".rds")))

message("→ MRDetect tables (raw & z-scored) written to: ", outdir)
