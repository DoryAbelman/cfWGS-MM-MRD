source("setup_packages.R")
source("config.R")
source("helpers.R")

# =============================================================================
# sample_availability_summary.R
# Project:  cfWGS MRD Detection (M4 / IMMAGINE / SPORE)
# Author:   Dory Abelman
# Date:     January 2025
# Last Updated: May 2025
#
# Purpose:
#   1. Load raw BM sample list and processing log; clean and compute DNA availability per sample.
#   2. Merge with clinical metadata and feature‐level “Evidence_of_Disease” to flag
#      baseline (Diagnosis/Baseline) and progression (Progression/Relapse) BM or cfDNA samples.
#   3. Build a patient‐level summary of sample availability, quality flags, study membership, and eligibility
#      for MRDetect (requires ≥2 cfDNA timepoints or ≥1 BM timepoint).
#   4. Export:
#       • Filtered_TFRIM4_Processing_Log.csv
#       • summary_table_of_samples_and_patient_availability_cfWGS.txt
#       • high_quality_patients_list_for_baseline_mut_calling.{csv,rds}
#       • patient_cohort_assignment.{csv,rds}
#   5. Print console summaries of overall counts and per‐study breakdowns.
#
# Dependencies:
#   • readxl, dplyr, readr, stringr, glue
#
# Input Files:
#   • Clinical data/M4/M4 V1 BM processed at baseline.xlsx
#   • TFRIM4_Processing Log_Nov2024.xlsx   (sheet 6)
#   • combined_clinical_data_updated_April2025.csv
#   • All_feature_data_Feb2025.rds
#   • summary_table_of_samples_and_patient_availability_cfWGS - for making the flow chart of samples.xlsx
#
# Output Files (written to working directory or “Output_tables_2025”):
#   • Filtered_TFRIM4_Processing_Log.csv
#   • summary_table_of_samples_and_patient_availability_cfWGS.txt
#   • high_quality_patients_list_for_baseline_mut_calling.csv
#   • high_quality_patients_list_for_baseline_mut_calling.rds
#   • patient_cohort_assignment.csv
#   • patient_cohort_assignment.rds
#
# Usage:
#   Rscript sample_availability_summary.R
# =============================================================================

# =============================================================================
# sample_availability_summary.R
# Project: cfWGS MRD detection (M4 / IMMAGINE / SPORE)
# Author:  Dory Abelman
# Date:    January 2025 (last update May 2025)
# =============================================================================


# ──────────────────────────────────────────────────────────────────────────────
# 1. FILE PATHS & SETUP
# ──────────────────────────────────────────────────────────────────────────────
bm_list_path        <- "Clinical data/M4/M4 V1 BM processed at baseline.xlsx"
processing_log_path <- "TFRIM4_Processing Log_Nov2024.xlsx"
clinical_csv_path   <- "combined_clinical_data_updated_April2025.csv"
features_rds_path   <- "All_feature_data_Feb2025.rds"
failed_info_path    <- "summary_table_of_samples_and_patient_availability_cfWGS - for making the flow chart of samples.xlsx"
export_dir          <- "Output_tables_2025"

if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 2. CLEAN PROCESSING LOG
# ──────────────────────────────────────────────────────────────────────────────

# 2.1 Read raw files
bm_data  <- read_excel(bm_list_path)
log_data <- read_excel(processing_log_path, sheet = 6)

# 2.2 Tag “sent to OICR” if any column mentions “OICR|30X|40X”
log_data <- log_data %>%
  mutate(
    sent_to_OICR = if_else(
      apply(log_data, 1, function(row) {
        any(str_detect(row, regex("OICR|30X|40X", ignore_case = TRUE)), na.rm = TRUE)
      }),
      1, 0
    )
  )

# 2.3 Keep only relevant columns
keep_columns <- c(
  "Specimen ID (red= 35 patient list)",
  "DNA Qubit (ng/ul)",
  "DNA total (ng)",
  "Pugh Lab Re-Qubit",
  "Pugh Lab vol Check",
  "Pugh Lab total ng available",
  "Notes",
  "Library Prep Date",
  "Library Prep Technician",
  "Sample Vol used for Library Prep (uL)",
  "Amount of DNA used for Library Prep (ng)",
  "Amount of DNA remaining post-library prep (ng)",
  "Library ID",
  "Total Library Yield (ng)"
)

filtered_log <- log_data %>%
  select(all_of(keep_columns), sent_to_OICR) %>%
  mutate(
    across(
      c(
        `Pugh Lab total ng available`,
        `Amount of DNA used for Library Prep (ng)`,
        `Amount of DNA remaining post-library prep (ng)`,
        `Total Library Yield (ng)`
      ),
      as.numeric
    )
  ) %>%
  mutate(
    DNA_available = if_else(
      !is.na(`Pugh Lab total ng available`),
      `Pugh Lab total ng available`,
      rowSums(
        select(
          ., 
          `Amount of DNA used for Library Prep (ng)`, 
          `Amount of DNA remaining post-library prep (ng)`
        ),
        na.rm = TRUE
      )
    )
  )

write_csv(filtered_log, "Filtered_TFRIM4_Processing_Log.csv")
message("Saved → Filtered_TFRIM4_Processing_Log.csv")


# ──────────────────────────────────────────────────────────────────────────────
# 3. HELPER TO FETCH PATIENT LISTS
# ──────────────────────────────────────────────────────────────────────────────

get_patient_list <- function(df, sample_type, timepoints, evidence_col = NULL, evidence_val = NULL) {
  df2 <- df %>%
    filter(
      Sample_type %in% sample_type,
      timepoint_info %in% timepoints
    )
  if (!is.null(evidence_col) && !is.null(evidence_val)) {
    df2 <- df2 %>% filter(!!sym(evidence_col) == evidence_val)
  }
  df2 %>% distinct(Patient) %>% pull(Patient)
}


# ──────────────────────────────────────────────────────────────────────────────
# 4. LOAD CLINICAL & FEATURE DATA
# ──────────────────────────────────────────────────────────────────────────────

combined_clinical <- read_csv(clinical_csv_path)
All_feature_data   <- readRDS(features_rds_path)

all_patients <- combined_clinical %>% distinct(Patient) %>% pull(Patient)


# ──────────────────────────────────────────────────────────────────────────────
# 5. PATIENT‐LEVEL AVAILABILITY & QUALITY
# ──────────────────────────────────────────────────────────────────────────────

bm_baseline   <- get_patient_list(combined_clinical, "BM_cells",           c("Diagnosis","Baseline"))
cfDNA_baseline<- get_patient_list(combined_clinical, "Blood_plasma_cfDNA", c("Diagnosis","Baseline"))

qual_bm       <- get_patient_list(All_feature_data,   "BM_cells",           c("Diagnosis","Baseline"), "Evidence_of_Disease", 1)
qual_cfDNA    <- get_patient_list(All_feature_data,   "Blood_plasma_cfDNA", c("Diagnosis","Baseline"), "Evidence_of_Disease", 1)

cfDNA_counts  <- combined_clinical %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  summarise(total_cfDNA_samples = n(), .groups = "drop")

# Build base summary table
summary_table <- tibble(Patient = all_patients) %>%
  mutate(
    BM_diagnosis_baseline    = if_else(Patient %in% bm_baseline,    1, 0),
    cfDNA_diagnosis_baseline = if_else(Patient %in% cfDNA_baseline, 1, 0),
    Qualifying_BM            = case_when(
      !Patient %in% bm_baseline       ~ NA_real_,
      Patient %in% qual_bm            ~ 1,
      TRUE                             ~ 0
    ),
    Qualifying_cfDNA         = case_when(
      !Patient %in% cfDNA_baseline    ~ NA_real_,
      Patient %in% qual_cfDNA         ~ 1,
      TRUE                             ~ 0
    )
  ) %>%
  left_join(cfDNA_counts, by = "Patient")


# Earliest BM Evidence_of_Disease
bm_first_EoD <- All_feature_data %>%
  filter(Sample_type == "BM_cells") %>%
  group_by(Patient) %>%
  slice_min(order_by = Date_of_sample_collection, with_ties = FALSE) %>%
  summarise(BM_first_sample_EoD = first(Evidence_of_Disease), .groups = "drop")

# Earliest cfDNA Evidence_of_Disease
cfDNA_first_EoD <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  slice_min(order_by = Date_of_sample_collection, with_ties = FALSE) %>%
  summarise(cfDNA_first_sample_EoD = first(Evidence_of_Disease), .groups = "drop")

summary_table <- summary_table %>%
  left_join(bm_first_EoD, by = "Patient") %>%
  left_join(cfDNA_first_EoD, by = "Patient")


# High-quality progression: BM & cfDNA
bm_hq_prog <- All_feature_data %>%
  filter(Sample_type == "BM_cells", timepoint_info %in% c("Progression","Relapse")) %>%
  group_by(Patient) %>%
  summarise(high_quality_progression_sample_BM = max(Evidence_of_Disease, na.rm = TRUE), .groups = "drop") %>%
  mutate(high_quality_progression_sample_BM = na_if(high_quality_progression_sample_BM, Inf))

bm_prog_available <- All_feature_data %>%
  filter(Sample_type == "BM_cells", timepoint_info %in% c("Progression","Relapse")) %>%
  group_by(Patient) %>%
  summarise(BM_progression_available = 1, .groups = "drop")

cfDNA_hq_prog <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% c("Progression","Relapse")) %>%
  group_by(Patient) %>%
  summarise(high_quality_progression_sample_cfDNA = max(Evidence_of_Disease, na.rm = TRUE), .groups = "drop") %>%
  mutate(high_quality_progression_sample_cfDNA = na_if(high_quality_progression_sample_cfDNA, Inf))

cfDNA_prog_available <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% c("Progression","Relapse")) %>%
  group_by(Patient) %>%
  summarise(cfDNA_progression_available = 1, .groups = "drop")

summary_table <- summary_table %>%
  left_join(bm_hq_prog,      by = "Patient") %>%
  left_join(bm_prog_available, by = "Patient") %>%
  left_join(cfDNA_hq_prog,   by = "Patient") %>%
  left_join(cfDNA_prog_available, by = "Patient") %>%
  mutate(
    BM_progression_available    = replace_na(BM_progression_available, 0),
    cfDNA_progression_available = replace_na(cfDNA_progression_available, 0)
  )


# Eligibility for MRDetect (baseline ∪ progression) & high-quality baseline
summary_table <- summary_table %>%
  mutate(
    eligible_in_MRDetect_study = if_else(
      (
        (BM_diagnosis_baseline == 1 | cfDNA_diagnosis_baseline == 1) &
          total_cfDNA_samples >= 2
      ) |
        (
          (BM_progression_available == 1 | cfDNA_progression_available == 1) &
            total_cfDNA_samples >= 2
        ),
      1, 0
    ),
    High_quality_baseline = if_else(
      (Qualifying_BM == 1 | Qualifying_cfDNA == 1) & total_cfDNA_samples >= 2, 1, 0
    )
  )

high_quality_patients <- summary_table %>%
  filter(High_quality_baseline == 1) %>%
  pull(Patient)

write_csv(high_quality_patients, "high_quality_patients_list_for_baseline_mut_calling.csv")
saveRDS(high_quality_patients,   "high_quality_patients_list_for_baseline_mut_calling.rds")
message("Saved → high_quality_patients_list_for_baseline_mut_calling.{csv,rds}")


# ──────────────────────────────────────────────────────────────────────────────
# 6. REPORT & SUMMARIES
# ──────────────────────────────────────────────────────────────────────────────

# 6.1 Study assignment
summary_table <- summary_table %>%
  mutate(
    Study = case_when(
      str_starts(Patient, "IMG")   ~ "IMMAGINE",
      str_starts(Patient, "SPORE") ~ "SPORE",
      TRUE                           ~ "M4"
    )
  )

study_counts <- summary_table %>%
  filter(eligible_in_MRDetect_study == 1) %>%
  count(Study, name = "Total_Eligible_Patients")

cat("Eligible patients per Study:\n")
print(study_counts); cat("\n")


# 6.2 Combined high-quality summary (baseline ∪ progression ∪ failure)
summary_filtered <- summary_table %>%
  filter(
    BM_diagnosis_baseline == 1 |
      BM_progression_available == 1 |
      BM_failed == 1 |
      cfDNA_diagnosis_baseline == 1 |
      cfDNA_progression_available == 1 |
      cfDNA_failed == 1
  ) %>%
  drop_na(total_cfDNA_samples)

final_summary <- summary_filtered %>%
  mutate(
    any_BM        = BM_diagnosis_baseline == 1 | BM_progression_available == 1 | BM_failed == 1,
    any_cfDNA     = cfDNA_diagnosis_baseline == 1 | cfDNA_progression_available == 1 | cfDNA_failed == 1,
    HQ_BM_bl      = Qualifying_BM == 1,
    HQ_cfDNA_bl   = Qualifying_cfDNA == 1,
    HQ_BM_prog    = high_quality_progression_sample_BM == 1,
    HQ_cfDNA_prog = high_quality_progression_sample_cfDNA == 1
  )

patient_totals <- final_summary %>%
  summarise(
    Total_Patients                = n_distinct(Patient),
    Total_BM                      = sum(any_BM, na.rm = TRUE),
    Total_cfDNA                   = sum(any_cfDNA, na.rm = TRUE),
    Total_Baseline_Samples        = sum(BM_diagnosis_baseline == 1, na.rm = TRUE) +
      sum(cfDNA_diagnosis_baseline == 1, na.rm = TRUE) +
      sum(BM_failed == 1, na.rm = TRUE) +
      sum(cfDNA_failed == 1, na.rm = TRUE),
    High_Quality_Baseline_BM      = sum(HQ_BM_bl, na.rm = TRUE),
    High_Quality_Baseline_cfDNA   = sum(HQ_cfDNA_bl, na.rm = TRUE),
    Total_Progression_Samples     = sum(BM_progression_available == 1, na.rm = TRUE) +
      sum(cfDNA_progression_available == 1, na.rm = TRUE),
    High_Quality_Progression_BM   = sum(HQ_BM_prog, na.rm = TRUE),
    High_Quality_Progression_cfDNA= sum(HQ_cfDNA_prog, na.rm = TRUE),
    Num_Failed_BM                 = sum(BM_failed == 1, na.rm = TRUE),
    Num_Failed_cfDNA              = sum(cfDNA_failed == 1, na.rm = TRUE)
  )

cat("Overall patient/sample summary:\n")
print(patient_totals); cat("\n")

study_totals <- final_summary %>%
  filter(any_BM | any_cfDNA) %>%
  group_by(Study) %>%
  summarise(
    Patients_in_Study     = n_distinct(Patient),
    BM_Samples_in_Study   = sum(any_BM, na.rm = TRUE),
    cfDNA_Samples_in_Study= sum(any_cfDNA, na.rm = TRUE),
    .groups = "drop"
  )

cat("Patients & samples per Study (filtered):\n")
print(study_totals); cat("\n")


# ──────────────────────────────────────────────────────────────────────────────
# 7. COHORT ASSIGNMENT FOR BASELINE PASSING SAMPLES
# ──────────────────────────────────────────────────────────────────────────────

failed_info <- read_excel(failed_info_path) %>%
  rename(
    BM_status    = `BM Status`,
    cfDNA_status = `cfDNA Status`
  ) %>%
  filter(!Patient %in% c("IMG-146", "IMG-163")) # exclude missing‐date cases

stats <- failed_info %>%
  summarise(
    bm_collected = sum(BM_status %in% c("Sequenced_pass","Sequenced_fail_QC","Failed_extraction","Failed_insufficient","Missed_release","Depleted_prestudy")),
    bm_seq_pass  = sum(BM_status == "Sequenced_pass"),
    bm_failures  = sum(BM_status %in% c("Sequenced_fail_QC","Failed_extraction","Failed_insufficient")),
    cf_collected = sum(cfDNA_status %in% c("Sequenced_pass","Sequenced_fail_QC","Failed_extraction","Failed_insufficient","Missed_release","Depleted_prestudy")),
    cf_seq_pass  = sum(cfDNA_status == "Sequenced_pass"),
    cf_failures  = sum(cfDNA_status %in% c("Sequenced_fail_QC","Failed_extraction","Failed_insufficient"))
  ) %>%
  mutate(
    bm_attempts     = bm_collected,
    bm_pass_rate    = round(bm_seq_pass / bm_attempts * 100, 1),
    bm_failure_rate = round(bm_failures  / bm_attempts * 100, 1),
    cf_attempts     = cf_collected,
    cf_pass_rate    = round(cf_seq_pass / cf_attempts * 100, 1),
    cf_failure_rate = round(cf_failures  / cf_attempts * 100, 1)
  )

cat(
  glue(
    "Of {bm_attempts} BM specimens attempted:\n",
    "- {bm_seq_pass} ({bm_pass_rate}%) passed QC; {bm_failures} ({bm_failure_rate}%) failed.\n\n",
    "Of {cf_attempts} cfDNA libraries attempted:\n",
    "- {cf_seq_pass} ({cf_pass_rate}%) passed QC; {cf_failures} ({cf_failure_rate}%) failed.\n"
  ),
  "\n"
)

patient_cohort <- failed_info %>%
  filter(BM_status == "Sequenced_pass" | cfDNA_status == "Sequenced_pass") %>%
  distinct(Patient, Study) %>%
  mutate(
    Cohort = case_when(
      Study == "SPORE"                    ~ "Non-frontline",
      Patient %in% c("IMG-098","IMG-159") ~ "Non-frontline",
      TRUE                                 ~ "Frontline"
    )
  ) %>%
  select(Patient, Cohort)

write_csv(patient_cohort, file.path(export_dir, "patient_cohort_assignment.csv"))
saveRDS(patient_cohort,   file.path(export_dir, "patient_cohort_assignment.rds"))
message("Saved → patient_cohort_assignment.{csv,rds} in ", export_dir)

# ──────────────────────────────────────────────────────────────────────────────
# End of file
# =============================================================================
