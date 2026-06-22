# =============================================================================
# 1_6_Identify_High_Quality_Patient_Pairs.R
# Project:  cfWGS MRD Detection (M4 / IMMAGINE / SPORE)
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_6_Identify_High_Quality_Patient_Pairs.R
#
# Manuscript outputs created/updated:
#   - Figure 1B source table: patient/sample availability and quality-control
#     counts used to assemble the final study-flow diagram.
#
# Pipeline role:
#   This script defines which patients and sample pairs are eligible for the
#   downstream baseline analyses. It does not draw the final flowchart panel
#   itself; instead, it exports the auditable source table copied into
#   final_manuscript_objects/Figure_1B.
#
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
#   4. Export the Figure 1B source table to `Final Tables and Figures/`, and
#      copy the same source table into the manuscript-labeled
#      `final_manuscript_objects/Figure_1B/` folder.
#   5. Export reusable pipeline intermediates for downstream mutation calling
#      and cohort assignment under `Output_tables_2025/`.
#   6. Export support-only sample-processing QC under
#      `Output_tables_2025/sample_qc_support/`.
#   7. Print console summaries of overall counts and per-study breakdowns.
#
# Dependencies:
#   • readxl, dplyr, readr, stringr, glue
#
# Input Files:
#   • Clinical data/M4/M4 V1 BM processed at baseline.xlsx
#   • TFRIM4_Processing Log_Nov2024.xlsx   (sheet 6)
#   • combined_clinical_data_updated_April2025.csv
#   • Jan2025_exported_data/All_feature_data_August2025.rds
#   • summary_table_of_samples_and_patient_availability_cfWGS - for making the flow chart of samples.xlsx
#
# Output Files:
#   Manuscript-facing:
#   - Final Tables and Figures/Table for creating sample flowchart updated3.csv
#     (Figure 1B source table copied into final_manuscript_objects/)
#
#   Pipeline intermediates:
#   - Output_tables_2025/high_quality_patients_list_for_baseline_mut_calling2.csv
#   - Output_tables_2025/high_quality_patients_list_for_baseline_mut_calling2.rds
#   - Output_tables_2025/patient_cohort_assignment.csv
#   - Output_tables_2025/patient_cohort_assignment.rds
#
#   Support-only QC:
#   - Output_tables_2025/sample_qc_support/Filtered_TFRIM4_Processing_Log.csv
#
# Usage:
#   Rscript 1_6_Identify_High_Quality_Patient_Pairs.R
# =============================================================================
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

library(readxl)
library(dplyr)
library(readr)
library(stringr)
library(glue)
library(tidyr)

# Shared manuscript-output helpers.
# Figure 1B itself is assembled outside R, but this script creates the exact
# source table used for that assembly and copies it into final_manuscript_objects/.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)

# ──────────────────────────────────────────────────────────────────────────────
# 1. FILE PATHS & SETUP
# ──────────────────────────────────────────────────────────────────────────────
bm_list_path        <- "Clinical data/M4/M4 V1 BM processed at baseline.xlsx"
processing_log_path <- "TFRIM4_Processing Log_Nov2024.xlsx"
clinical_csv_path   <- "combined_clinical_data_updated_April2025.csv"
features_rds_path   <- "Jan2025_exported_data/All_feature_data_August2025.rds"
failed_info_path    <- "summary_table_of_samples_and_patient_availability_cfWGS - for making the flow chart of samples.xlsx"
export_dir          <- "Output_tables_2025"
final_tables_dir    <- "Final Tables and Figures"
support_dir         <- file.path(export_dir, "sample_qc_support")

if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
if (!dir.exists(final_tables_dir)) dir.create(final_tables_dir, recursive = TRUE)
if (!dir.exists(support_dir)) dir.create(support_dir, recursive = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 2. CLEAN PROCESSING LOG
# ──────────────────────────────────────────────────────────────────────────────

# 2.1 Read raw files
bm_data  <- read_excel(bm_list_path, .name_repair = "unique_quiet")
log_data <- read_excel(processing_log_path, sheet = 6, .name_repair = "unique_quiet")

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
      ~ suppressWarnings(as.numeric(.x))
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

filtered_log_path <- file.path(support_dir, "Filtered_TFRIM4_Processing_Log.csv")
write_csv(filtered_log, filtered_log_path)
message("Saved support QC table -> ", filtered_log_path)


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

max_flag_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  out <- max(x, na.rm = TRUE)
  if (is.infinite(out)) NA_real_ else out
}


# ──────────────────────────────────────────────────────────────────────────────
# 4. LOAD CLINICAL & FEATURE DATA
# ──────────────────────────────────────────────────────────────────────────────

combined_clinical <- read_csv(clinical_csv_path, show_col_types = FALSE)
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
  summarise(total_cfDNA_samples = dplyr::n(), .groups = "drop")

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
  summarise(high_quality_progression_sample_BM = max_flag_or_na(Evidence_of_Disease), .groups = "drop")

bm_prog_available <- All_feature_data %>%
  filter(Sample_type == "BM_cells", timepoint_info %in% c("Progression","Relapse")) %>%
  group_by(Patient) %>%
  summarise(BM_progression_available = 1, .groups = "drop")

cfDNA_hq_prog <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% c("Progression","Relapse")) %>%
  group_by(Patient) %>%
  summarise(high_quality_progression_sample_cfDNA = max_flag_or_na(Evidence_of_Disease), .groups = "drop")

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

high_quality_csv <- file.path(export_dir, "high_quality_patients_list_for_baseline_mut_calling2.csv")
high_quality_rds <- file.path(export_dir, "high_quality_patients_list_for_baseline_mut_calling2.rds")
write_csv(tibble(Patient = high_quality_patients), high_quality_csv)
saveRDS(high_quality_patients, high_quality_rds)
message("Saved high-quality mutation patient list -> ", high_quality_csv, " and ", high_quality_rds)


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
failed_flags <- read_excel(failed_info_path, .name_repair = "unique_quiet") %>%
  rename(
    BM_status    = `BM Status`,
    cfDNA_status = `cfDNA Status`
  ) %>%
  filter(!Patient %in% c("IMG-146", "IMG-163")) %>%
  group_by(Patient) %>%
  summarise(
    BM_failed = as.integer(any(BM_status %in% c("Sequenced_fail_QC", "Failed_extraction", "Failed_insufficient"), na.rm = TRUE)),
    cfDNA_failed = as.integer(any(cfDNA_status %in% c("Sequenced_fail_QC", "Failed_extraction", "Failed_insufficient"), na.rm = TRUE)),
    .groups = "drop"
  )

summary_table <- summary_table %>%
  left_join(failed_flags, by = "Patient") %>%
  mutate(
    BM_failed = replace_na(BM_failed, 0L),
    cfDNA_failed = replace_na(cfDNA_failed, 0L)
  )

summary_filtered <- summary_table %>%
  filter(
    BM_diagnosis_baseline == 1 |
      BM_progression_available == 1 |
#      BM_failed == 1 |
      cfDNA_diagnosis_baseline == 1 |
      cfDNA_progression_available == 1 
 #     cfDNA_failed == 1
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
 #   Num_Failed_BM                 = sum(BM_failed == 1, na.rm = TRUE),
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

failed_info <- read_excel(failed_info_path, .name_repair = "unique_quiet") %>%
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
  glue::glue_data(
    stats,
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
message("Saved patient cohort assignment -> ", file.path(export_dir, "patient_cohort_assignment.{csv,rds}"))



# Analyst note:
#   `cohort_assignment_table_updated.rds` is the reviewed cohort assignment used
#   by the submitted analysis. The earlier `patient_cohort_assignment.rds` export
#   above is retained as a reproducible pipeline intermediate, while this curated
#   assignment is used to annotate the manuscript Figure 1B source table.
cohort_df <- readRDS("cohort_assignment_table_updated.rds")

failed_info <- failed_info %>% 
  left_join(cohort_df)

# Identify the cohort patients used as denominators for the disease-evidence
# coverage summaries below.
cohort_patients <- cohort_df %>% distinct(Patient)

# Which cohort patients have a baseline BM sample with evidence of disease?
bm_good_patients <- All_feature_data %>%
  filter(
    Sample_type == "BM_cells",
    timepoint_info %in% c("Diagnosis","Baseline"),
    Evidence_of_Disease == 1
  ) %>%
  distinct(Patient)

# Which cohort patients have a baseline cfDNA sample with evidence of disease?
cfDNA_good_patients <- All_feature_data %>%
  filter(
    Sample_type == "Blood_plasma_cfDNA",
    timepoint_info %in% c("Diagnosis","Baseline"),
    Evidence_of_Disease == 1
  ) %>%
  distinct(Patient)

# Add sample-quality flags to the Figure 1B source table.
failed_info <- failed_info %>%
  mutate(
    high_quality_BM    = as.integer(Patient %in% bm_good_patients$Patient),
    high_quality_cfDNA = as.integer(Patient %in% cfDNA_good_patients$Patient)
  )

failed_info <- failed_info %>%
  mutate(
    Cohort = case_when(
      Study == "M4" & is.na(Cohort)        ~ "Frontline_omitted",
      TRUE                                  ~ Cohort
    )
  )

fig1b_source_table <- file.path(final_tables_dir, "Table for creating sample flowchart updated3.csv")
write.csv(failed_info, file = fig1b_source_table)

# MANUSCRIPT OUTPUT: Figure 1B source table
# The final Figure 1B visual panel is created from this table in the manuscript
# figure file. This is the command-line-regenerated source table that documents
# which samples/patients pass sequencing, quality, cohort, and disease-evidence
# filters for the sample-flow diagram.
ms_copy_artifact(
  source_path = fig1b_source_table,
  artifact_id = "FIG1B",
  role = "figure_source_csv",
  description = "Figure 1B source table: patient/sample flowchart counts and cohort/QC annotations.",
  script_name = "1_6_Identify_High_Quality_Patient_Pairs.R"
)

# Support-only QA summary:
#   These proportions are printed to the console so the analyst can confirm the
#   fraction of cohort patients with disease-evaluable baseline BM/cfDNA samples.
#   They are not copied to `final_manuscript_objects/` because no final
#   manuscript figure or table is mapped to these summaries.
prop_bm <- cohort_patients %>%
  mutate(has_good_BM = Patient %in% bm_good_patients$Patient) %>%
  summarise(
    n_total = dplyr::n(),
    n_good  = sum(has_good_BM),
    prop_good_BM = n_good / n_total
  )

prop_cfDNA <- cohort_patients %>%
  mutate(has_good_cfDNA = Patient %in% cfDNA_good_patients$Patient) %>%
  summarise(
    n_total   = dplyr::n(),
    n_good    = sum(has_good_cfDNA),
    prop_good_cfDNA = n_good / n_total
  )

# Print the overall disease-evaluable proportions.
prop_bm
prop_cfDNA

# BM proportion by cohort and overall.
prop_bm_by_cohort <- failed_info %>%
  filter(BM_status == "Sequenced_pass") %>%
  group_by(Cohort) %>%
  summarise(
    n_total = dplyr::n(),
    n_good  = sum(high_quality_BM, na.rm = TRUE),
    prop_good_BM = n_good / n_total,
    .groups = "drop"
  ) %>%
  bind_rows(
    failed_info %>%
      filter(BM_status == "Sequenced_pass") %>%
      summarise(
        Cohort = "All",
        n_total = dplyr::n(),
        n_good  = sum(high_quality_BM, na.rm = TRUE),
        prop_good_BM = n_good / n_total
      )
  )

# cfDNA proportion by cohort and overall.
prop_cfDNA_by_cohort <- failed_info %>%
  filter(cfDNA_status == "Sequenced_pass") %>%
  group_by(Cohort) %>%
  summarise(
    n_total = dplyr::n(),
    n_good  = sum(high_quality_cfDNA, na.rm = TRUE),
    prop_good_cfDNA = n_good / n_total,
    .groups = "drop"
  ) %>%
  bind_rows(
    failed_info %>%
      filter(cfDNA_status == "Sequenced_pass") %>%
      summarise(
        Cohort = "All",
        n_total = dplyr::n(),
        n_good  = sum(high_quality_cfDNA, na.rm = TRUE),
        prop_good_cfDNA = n_good / n_total
      )
  )

# Print the cohort-stratified disease-evaluable proportions.
prop_bm_by_cohort
prop_cfDNA_by_cohort


# ──────────────────────────────────────────────────────────────────────────────
# End of file
# =============================================================================
