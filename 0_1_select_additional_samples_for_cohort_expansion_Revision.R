# Eligible sample selection and classification script

# This R script loads a multi‐sheet Excel workbook containing
# patient diagnosis dates, sample collection information, MRD test dates
# and treatment line data.  It identifies patients with a baseline
# sample (bone marrow or peripheral blood) collected within ±30 days of
# their documented date of diagnosis.  For the patients meeting this
# baseline criterion, it then classifies all samples into one of three
# major categories (Baseline, MRD or Treatment) based on the
# relationship between the sample collection date and key clinical
# events.  The script also flags whether each patient ultimately
# progressed, and it records the number of days between baseline and
# MRD samples relative to the associated diagnosis or MRD test dates.
#
# Author: Dory Abelman
# Date: 2026-02-05

## Load required libraries
suppressPackageStartupMessages({
  library(readxl)    # for reading Excel files
  library(dplyr)     # for data manipulation
  library(lubridate) # for handling dates
  library(stringr)   # for string operations
  library(tidyr)     # for data tidying
})

#-------------------------------------------------------------------
# 0. Load and process patient ID mapping and complete sample inventory
#-------------------------------------------------------------------

cat("\n========== LOADING PATIENT ID MAPPING AND SAMPLE INVENTORY ==========\n\n")

# Load patient ID mapping file (links PATIENT_ID to STUDY_NAME)
id_mapping_file <- "Clinical data/IMMAGINE/IMMAGINE  LIBERATE_ID_noMRNorDoB_6Feb2026.xlsx"

if (!file.exists(id_mapping_file)) {
  stop("ERROR: ID mapping file not found at: ", id_mapping_file)
}

id_mapping <- read_excel(id_mapping_file) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), str_trim))

cat("✓ Loaded ID mapping file with", nrow(id_mapping), "rows\n")
cat("  Columns:", paste(names(id_mapping), collapse = ", "), "\n\n")

# Load complete sample inventory (all samples from both studies)
sample_inventory_file <- "Clinical data/IMMAGINE/2025-05-01_LIB-IMMAGINE_sample_report.xlsx"

if (!file.exists(sample_inventory_file)) {
  stop("ERROR: Sample inventory file not found at: ", sample_inventory_file)
}

# Get sheet names and load both sheets
sheet_names <- excel_sheets(sample_inventory_file)
cat("✓ Found", length(sheet_names), "sheets in sample inventory:", paste(sheet_names, collapse = ", "), "\n")

# Load both sheets
sheet1 <- read_excel(sample_inventory_file, sheet = sheet_names[1]) %>%
  mutate(source_sheet = sheet_names[1])

sheet2 <- read_excel(sample_inventory_file, sheet = sheet_names[2]) %>%
  mutate(source_sheet = sheet_names[2])

# Merge sheets together
sample_inventory <- bind_rows(sheet1, sheet2)

cat("  Sheet 1:", nrow(sheet1), "samples\n")
cat("  Sheet 2:", nrow(sheet2), "samples\n")
cat("  Combined:", nrow(sample_inventory), "samples\n\n")

# Match Event Comment to STUDY_NAME in ID mapping to add PATIENT_ID
# Clean the Event Comment field for matching
sample_inventory_with_ids <- sample_inventory %>%
  mutate(
    Event_Comment_clean = str_trim(`Event Comment`),
    STUDY_NAME = Event_Comment_clean
  ) %>%
  left_join(
    id_mapping %>% select(PATIENT_ID, STUDY_NAME),
    by = "STUDY_NAME"
  )

# Report matching statistics
n_matched <- sum(!is.na(sample_inventory_with_ids$PATIENT_ID))
n_unmatched <- sum(is.na(sample_inventory_with_ids$PATIENT_ID))

cat("========== SAMPLE-TO-PATIENT MATCHING RESULTS ==========\n")
cat("  Matched to patient ID:", n_matched, "samples\n")
cat("  Unmatched (no patient ID):", n_unmatched, "samples\n")
cat("  Match rate:", sprintf("%.1f%%", 100 * n_matched / nrow(sample_inventory_with_ids)), "\n\n")

# Show unique Event Comments that didn't match
if (n_unmatched > 0) {
  unmatched_events <- sample_inventory_with_ids %>%
    filter(is.na(PATIENT_ID)) %>%
    count(Event_Comment_clean) %>%
    arrange(desc(n))

  cat("Top unmatched Event Comments:\n")
  print(head(unmatched_events, 10))
  cat("\n")

  # Export unmatched Event Comment values as a CSV and an RDS list
  unmatched_event_list <- unmatched_events$Event_Comment_clean

  # CSV with counts
  write.csv(as.data.frame(unmatched_events), file = "unmatched_event_comments.csv", row.names = FALSE)

  # RDS containing character vector of unmatched Event Comment strings
  saveRDS(unmatched_event_list, file = "unmatched_event_comments.rds")

  cat("✓ Exported unmatched Event Comments to 'unmatched_event_comments.csv' and 'unmatched_event_comments.rds'\n\n")
}

# Export the complete sample inventory with patient IDs
write.csv(
  sample_inventory_with_ids,
  file = "complete_sample_inventory_with_patient_IDs.csv",
  row.names = FALSE
)
cat("✓ Complete sample inventory written to 'complete_sample_inventory_with_patient_IDs.csv'\n\n")

# Create summary table by patient
patient_sample_summary <- sample_inventory_with_ids %>%
  filter(!is.na(PATIENT_ID)) %>%
  group_by(PATIENT_ID) %>%
  summarise(
    n_samples = n(),
    sample_types = paste(unique(`Specimen Type`), collapse = "; "),
    dates = paste(unique(`Specimen collection date`), collapse = "; "),
    .groups = 'drop'
  )

write.csv(
  patient_sample_summary,
  file = "patient_sample_summary.csv",
  row.names = FALSE
)
cat("✓ Patient sample summary written to 'patient_sample_summary.csv'\n")
cat("  Total patients with samples:", nrow(patient_sample_summary), "\n\n")

cat("========== ID MAPPING AND INVENTORY COMPLETE ==========\n\n")

## Now work on getting the info of available samples for the comparison
# Define the path to the Excel workbook.  The workbook must reside in
# the working directory or you must adjust the path accordingly.
file_path <- "Clinical data/IMMAGINE/IMMAGINE  LIBERATE_MRD_withoutMRN_or_DOB_18Dec2025.xlsx"

# Check if file exists, provide helpful error message if not
if (!file.exists(file_path)) {
  stop(
    "ERROR: Excel file not found at: ", file_path, "\n",
    "Current working directory: ", getwd(), "\n",
    "Please ensure the file exists at the specified path or update the file_path variable."
  )
}

#-------------------------------------------------------------------
# 1. Import data from each worksheet
#-------------------------------------------------------------------

# Read diagnosis dates
date_dx <- read_excel(
  file_path, sheet = "IMM_MRD+_DateDx",
  col_types = c(
    PATIENT_ID = "text",            # preserve patient identifiers as text
    STUDY_NAME = "text",
    GENDER = "text",
    DATE_DIAGNOSIS = "date",
    PRIMARY_PATHOLOGY = "text",
    VITAL_STATUS = "text",
    CAUSE_OF_DEATH = "text",
    DATE_OF_DEATH = "date",
    DATE_OF_LAST_FOLLOWUP = "date"
  )
)

# Read sample collection records
samples_raw <- read_excel(
  file_path, sheet = "IMMAGINE & LIBERATE SAMPLES",
  col_types = c(
    STUDY_NAME = "text",
    PATIENT_ID = "text",
    GENDER = "text",
    COLLECTION_DATE = "date",
    PROCEDURE_TYPE_TEXT = "text",
    SAMPLE_COLLECTION_DATE = "date",
    TISSUE_REPOSITORY_TEXT = "text"
  )
)

# Read MRD test dates
mrd_raw <- read_excel(
  file_path, sheet = "MRD",
  col_types = c(
    STUDY_NAME = "text",
    PATIENT_ID = "text",
    GENDER = "text",
    BM_COLLCTN_DATE = "date",
    MALIGNANT_PC = "text",
    `RESIDUAL_PC%` = "numeric"
  )
)

# Read treatment line data
line_tx_raw <- read_excel(
  file_path, sheet = "IMM_MRD+_Line_Tx",
  col_types = c(
    PATIENT_ID = "text",
    STUDY_NAME = "text",
    REGIMEN_NAME = "text",
    LINE_OF_TREATMENT = "numeric",
    START_DATE = "date",
    END_DATE = "date",
    STUDY_DRUG = "text",
    INTENT = "text",
    BEST_RESPONSE_DATE = "date",
    BEST_RESPONSE = "text",
    PROGRESSION = "text",
    PROGRESSION_DATE = "date",
    PROGRESSION_DETAILS = "text",
    ADDITIONAL_DETAILS = "text"
  )
)

#-------------------------------------------------------------------
# 2. Data preparation and cleaning
#-------------------------------------------------------------------

# Prepare diagnosis data
date_dx_clean <- date_dx %>%
  mutate(
    PATIENT_ID = str_trim(PATIENT_ID),
    DATE_DIAGNOSIS = as.Date(DATE_DIAGNOSIS)
  ) %>%
  select(PATIENT_ID, DATE_DIAGNOSIS)

# Prepare sample data: unify collection date and assign sample ID
samples_clean <- samples_raw %>%
  mutate(
    PATIENT_ID = str_trim(PATIENT_ID),
    # choose SAMPLE_COLLECTION_DATE when available; otherwise fall back to COLLECTION_DATE
    sample_date = coalesce(SAMPLE_COLLECTION_DATE, COLLECTION_DATE),
    sample_date = as.Date(sample_date),
    sample_type = if_else(is.na(PROCEDURE_TYPE_TEXT), "Unknown", PROCEDURE_TYPE_TEXT),
    # Standardize sample types: treat Plasma as Peripheral blood sample
    sample_type = if_else(sample_type == "Plasma", "Peripheral blood sample", sample_type),
    sample_id = row_number(),
    data_source = "clinical_excel"
  ) %>%
  select(sample_id, PATIENT_ID, STUDY_NAME, GENDER, sample_type, sample_date, data_source)

# Transform sample inventory to match samples_clean structure
cat("\n========== INTEGRATING SAMPLE INVENTORY WITH CLINICAL SAMPLES ==========\n")

inventory_for_merge <- sample_inventory_with_ids %>%
  filter(!is.na(PATIENT_ID)) %>%  # Only keep samples with matched patient IDs
  mutate(
    PATIENT_ID = str_trim(PATIENT_ID),
    sample_date = as.Date(`Specimen collection date`),
    sample_type = if_else(is.na(`Specimen Type`), "Unknown", `Specimen Type`),
    # Standardize sample types: treat Plasma as Peripheral blood sample
    sample_type = if_else(sample_type == "Plasma", "Peripheral blood sample", sample_type),
    GENDER = NA_character_,  # Not available in inventory
    data_source = "sample_inventory"
  ) %>%
  select(PATIENT_ID, STUDY_NAME, GENDER, sample_type, sample_date, data_source)

# Assign sample_id after combining
inventory_for_merge <- inventory_for_merge %>%
  mutate(sample_id = row_number() + max(samples_clean$sample_id))

cat("  Clinical Excel samples:", nrow(samples_clean), "\n")
cat("  Sample Inventory (matched):", nrow(inventory_for_merge), "\n")

# Combine the datasets
samples_combined <- bind_rows(
  samples_clean,
  inventory_for_merge
)

cat("  Combined total:", nrow(samples_combined), "\n\n")

# Check for potential duplicates (same patient, date, and type)
duplicates <- samples_combined %>%
  group_by(PATIENT_ID, sample_date, sample_type) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(duplicates) > 0) {
  cat("⚠️  Found", nrow(duplicates), "duplicate samples (same patient, date, type)\n")
  cat("   Removing duplicates and keeping first occurrence. Review 'removed_duplicate_samples.csv' for details.\n\n")
  
  write.csv(duplicates, file = "removed_duplicate_samples.csv", row.names = FALSE)
  
  # Remove duplicates by keeping only first occurrence
  samples_combined <- samples_combined %>%
    distinct(PATIENT_ID, sample_date, sample_type, .keep_all = TRUE)
  
  cat("✓ After deduplication:", nrow(samples_combined), "samples\n")
  cat("  Removed:", nrow(duplicates) - n_distinct(duplicates$PATIENT_ID, duplicates$sample_date, duplicates$sample_type), "duplicate rows\n\n")
} else {
  cat("✓ No duplicate samples detected\n\n")
}

# Replace samples_clean with the combined (and deduplicated) dataset
samples_clean <- samples_combined

cat("✓ Using combined dataset with", nrow(samples_clean), "samples from", 
    n_distinct(samples_clean$PATIENT_ID), "patients\n")
cat("  - From clinical Excel:", sum(samples_clean$data_source == "clinical_excel"), "\n")
cat("  - From sample inventory:", sum(samples_clean$data_source == "sample_inventory"), "\n\n")

cat("========== INTEGRATION COMPLETE ==========\n\n")

# Prepare MRD data
mrd_clean <- mrd_raw %>%
  mutate(
    PATIENT_ID = str_trim(PATIENT_ID),
    mrd_date = as.Date(BM_COLLCTN_DATE)
  ) %>%
  select(PATIENT_ID, mrd_date)

# Prepare treatment line data
line_tx_clean <- line_tx_raw %>%
  mutate(
    PATIENT_ID = str_trim(PATIENT_ID),
    start_date = as.Date(START_DATE),
    end_date = as.Date(END_DATE)
  ) %>%
  select(
    PATIENT_ID, STUDY_NAME, REGIMEN_NAME, LINE_OF_TREATMENT,
    start_date, end_date, PROGRESSION
  )

#-------------------------------------------------------------------
# 3. Summarize progression at the patient level
#-------------------------------------------------------------------

# Define a helper that interprets progression strings
interpret_progression <- function(progression) {
  ifelse(
    is.na(progression), NA,
    case_when(
      str_detect(str_to_lower(progression), "progressive") ~ TRUE,
      str_detect(str_to_lower(progression), "no evidence") ~ FALSE,
      TRUE ~ NA
    )
  )
}

progression_summary <- line_tx_clean %>%
  mutate(progressed_flag = interpret_progression(PROGRESSION)) %>%
  group_by(PATIENT_ID) %>%
  summarise(
    progressed = case_when(
      any(progressed_flag == TRUE, na.rm = TRUE) ~ "Yes",
      any(progressed_flag == FALSE, na.rm = TRUE) & !any(progressed_flag == TRUE, na.rm = TRUE) ~ "No",
      TRUE ~ "Unknown"
    )
  ) %>%
  ungroup()

#-------------------------------------------------------------------
# 4. Identify baseline samples (within ±30 days of diagnosis)
#-------------------------------------------------------------------

samples_with_dx <- samples_clean %>%
  left_join(date_dx_clean, by = "PATIENT_ID") %>%
  mutate(
    days_from_diagnosis = as.numeric(sample_date - DATE_DIAGNOSIS),
    baseline_flag = if_else(!is.na(days_from_diagnosis) & abs(days_from_diagnosis) <= 30, TRUE, FALSE)
  )

#-------------------------------------------------------------------
# 5. Identify MRD timepoint samples (within ±14 days of any MRD test)
#-------------------------------------------------------------------

# IMPORTANT: The MRD sheet contains only bone marrow (BM_COLLCTN_DATE).
# Blood samples at MRD timepoints are identified by proximity to these BM dates.
# 
# Potential limitations:
# - If blood samples are drawn at different times than BM during MRD testing,
#   they may be misclassified or missed
# - The ±14 day window may be too strict/loose for your use case
# - Consider validating that matched samples truly represent the same clinical event

# For each sample, find the nearest MRD test date and compute the days difference.
# We perform a left join to replicate samples for all MRD dates of the same patient,
# compute absolute differences, and then summarise to keep the nearest event.

if (nrow(mrd_clean) > 0) {
  mrd_diffs <- samples_with_dx %>%
    left_join(mrd_clean, by = "PATIENT_ID", relationship = "many-to-many") %>%
    mutate(
      diff_to_mrd = as.numeric(sample_date - mrd_date),
      abs_diff_to_mrd = abs(diff_to_mrd)
    ) %>%
    group_by(sample_id) %>%
    summarise(
      nearest_abs_diff_mrd = if_else(any(!is.na(abs_diff_to_mrd)), min(abs_diff_to_mrd, na.rm = TRUE), Inf),
      nearest_mrd_diff = {
        min_idx <- which.min(abs_diff_to_mrd)
        if_else(length(min_idx) == 0 || all(is.na(mrd_date)), NA_real_, diff_to_mrd[min_idx[1]])
      },
      nearest_mrd_date = {
        min_idx <- which.min(abs_diff_to_mrd)
        if_else(length(min_idx) == 0 || all(is.na(mrd_date)), as.Date(NA), mrd_date[min_idx[1]])
      },
      .groups = 'drop'
    ) %>%
    mutate(
      mrd_flag = if_else(nearest_abs_diff_mrd <= 14, TRUE, FALSE)
    )
} else {
  mrd_diffs <- samples_with_dx %>%
    select(sample_id) %>%
    mutate(
      nearest_abs_diff_mrd = Inf,
      nearest_mrd_diff = NA_real_,
      nearest_mrd_date = as.Date(NA),
      mrd_flag = FALSE
    )
}

# Merge MRD summary back to samples
samples_with_dx_mrd <- samples_with_dx %>%
  left_join(mrd_diffs, by = "sample_id")

#-------------------------------------------------------------------
# 6. Identify treatment samples (samples collected during any treatment line)
#-------------------------------------------------------------------

if (nrow(line_tx_clean) > 0) {
  treatment_flags <- samples_with_dx_mrd %>%
    select(sample_id, PATIENT_ID, sample_date) %>%
    left_join(line_tx_clean %>% select(PATIENT_ID, start_date, end_date, LINE_OF_TREATMENT, REGIMEN_NAME), by = "PATIENT_ID") %>%
    mutate(
      in_treatment = !is.na(start_date) & !is.na(sample_date) & sample_date >= start_date & (is.na(end_date) | sample_date <= end_date)
    ) %>%
    group_by(sample_id) %>%
    summarise(
      in_treatment = any(in_treatment, na.rm = TRUE),
      tx_lines = paste(unique(LINE_OF_TREATMENT[in_treatment]), collapse = ";"),
      tx_regimens = paste(unique(REGIMEN_NAME[in_treatment]), collapse = ";")
    ) %>%
    ungroup()
} else {
  treatment_flags <- samples_with_dx_mrd %>%
    select(sample_id) %>%
    mutate(
      in_treatment = FALSE,
      tx_lines = "",
      tx_regimens = ""
    )
}

# Merge treatment flags back to samples and add progression info
samples_final <- samples_with_dx_mrd %>%
  left_join(treatment_flags, by = "sample_id") %>%
  left_join(progression_summary, by = "PATIENT_ID") %>%
  mutate(
    sample_class = case_when(
      baseline_flag ~ "Baseline",
      !baseline_flag & mrd_flag ~ "MRD",
      !baseline_flag & !mrd_flag & in_treatment ~ "Treatment",
      TRUE ~ "Other"
    ),
    days_from_dx = if_else(sample_class == "Baseline", days_from_diagnosis, NA_real_),
    days_from_mrd = if_else(sample_class == "MRD", nearest_mrd_diff, NA_real_)
  )

#-------------------------------------------------------------------
# 7. Define two cohorts based on MRD sample type requirements
#-------------------------------------------------------------------

# Cohort 1: Baseline (BM or blood) + ANY MRD sample (BM or blood)
eligible_patients_any_mrd <- samples_final %>%
  group_by(PATIENT_ID) %>%
  summarise(
    has_baseline_bm_or_blood = any(
      baseline_flag & sample_type %in% c("Bone marrow sample", "Peripheral blood sample"),
      na.rm = TRUE
    ),
    # ANY MRD sample (BM or blood)
    has_any_mrd = any(
      sample_class == "MRD",
      na.rm = TRUE
    )
  ) %>%
  filter(has_baseline_bm_or_blood & has_any_mrd) %>%
  pull(PATIENT_ID)

# Cohort 2: Baseline (BM or blood) + AT LEAST ONE blood sample at MRD timepoint
eligible_patients_mrd_blood_only <- samples_final %>%
  group_by(PATIENT_ID) %>%
  summarise(
    has_baseline_bm_or_blood = any(
      baseline_flag & sample_type %in% c("Bone marrow sample", "Peripheral blood sample"),
      na.rm = TRUE
    ),
    # At least one MRD blood sample
    has_mrd_blood = any(
      sample_class == "MRD" & sample_type == "Peripheral blood sample",
      na.rm = TRUE
    )
  ) %>%
  filter(has_baseline_bm_or_blood & has_mrd_blood) %>%
  pull(PATIENT_ID)

eligible_samples_any_mrd <- samples_final %>%
  filter(PATIENT_ID %in% eligible_patients_any_mrd) %>%
  arrange(PATIENT_ID, sample_date, sample_class)

eligible_samples_mrd_blood_only <- samples_final %>%
  filter(PATIENT_ID %in% eligible_patients_mrd_blood_only) %>%
  arrange(PATIENT_ID, sample_date, sample_class)

#-------------------------------------------------------------------
# 8. Diagnostic output - analyze blood MRD eligibility failures
#-------------------------------------------------------------------

cat("\n===== DIAGNOSTIC REPORT: BLOOD MRD ELIGIBILITY ANALYSIS =====\n\n")

# Identify patients with blood MRD samples
patients_with_blood_mrd <- samples_final %>%
  filter(sample_class == "MRD" & sample_type == "Peripheral blood sample") %>%
  pull(PATIENT_ID) %>%
  unique()

cat("Total patients with blood MRD samples:", length(patients_with_blood_mrd), "\n\n")

# For each patient, analyze why they might fail eligibility
blood_mrd_analysis <- samples_final %>%
  filter(PATIENT_ID %in% patients_with_blood_mrd) %>%
  group_by(PATIENT_ID) %>%
  summarise(
    # Check eligibility criteria
    has_baseline_bm_or_blood = any(
      baseline_flag & sample_type %in% c("Bone marrow sample", "Peripheral blood sample"),
      na.rm = TRUE
    ),
    has_blood_mrd = any(
      sample_class == "MRD" & sample_type == "Peripheral blood sample",
      na.rm = TRUE
    ),
    # Additional diagnostic info
    total_samples = n(),
    n_timepoints = n_distinct(sample_class),
    baseline_samples = sum(baseline_flag, na.rm = TRUE),
    mrd_samples = sum(sample_class == "MRD", na.rm = TRUE),
    blood_mrd_samples = sum(sample_class == "MRD" & sample_type == "Peripheral blood sample", na.rm = TRUE),
    bm_mrd_samples = sum(sample_class == "MRD" & sample_type == "Bone marrow sample", na.rm = TRUE),
    baseline_blood = sum(baseline_flag & sample_type == "Peripheral blood sample", na.rm = TRUE),
    baseline_bm = sum(baseline_flag & sample_type == "Bone marrow sample", na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    is_eligible = has_baseline_bm_or_blood & has_blood_mrd,
    failure_reason = case_when(
      is_eligible ~ "ELIGIBLE",
      !has_baseline_bm_or_blood & has_blood_mrd ~ "No baseline BM or blood sample",
      has_baseline_bm_or_blood & !has_blood_mrd ~ "No blood MRD sample (should not appear here)",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(failure_reason, PATIENT_ID)

# Split into eligible and ineligible
blood_mrd_eligible <- blood_mrd_analysis %>% filter(is_eligible)
blood_mrd_ineligible <- blood_mrd_analysis %>% filter(!is_eligible)

# Output both to CSV
write.csv(blood_mrd_analysis, file = "blood_MRD_patient_analysis.csv", row.names = FALSE)
cat("✓ Full blood MRD analysis written to 'blood_MRD_patient_analysis.csv'\n\n")

cat("========== BLOOD MRD PATIENTS: ELIGIBLE ==========\n")
cat("Count:", nrow(blood_mrd_eligible), "\n\n")
if (nrow(blood_mrd_eligible) > 0) {
  print(blood_mrd_eligible %>% select(PATIENT_ID, total_samples, n_timepoints, baseline_blood, baseline_bm, blood_mrd_samples, bm_mrd_samples))
}

cat("\n========== BLOOD MRD PATIENTS: INELIGIBLE ==========\n")
cat("Count:", nrow(blood_mrd_ineligible), "\n\n")
if (nrow(blood_mrd_ineligible) > 0) {
  print(blood_mrd_ineligible %>% select(PATIENT_ID, failure_reason, total_samples, n_timepoints, baseline_blood, baseline_bm, blood_mrd_samples, bm_mrd_samples))
  
  cat("\n----- Failure Breakdown -----\n")
  failure_summary <- blood_mrd_ineligible %>%
    group_by(failure_reason) %>%
    summarise(count = n(), .groups = 'drop')
  print(failure_summary)
  
  cat("\n----- Patient IDs for Manual Review -----\n")
  for (reason in unique(blood_mrd_ineligible$failure_reason)) {
    patients <- blood_mrd_ineligible %>%
      filter(failure_reason == reason) %>%
      pull(PATIENT_ID) %>%
      paste(collapse = ", ")
    cat(sprintf("%s (%d patients):\n  %s\n\n", 
                reason, 
                length(blood_mrd_ineligible$PATIENT_ID[blood_mrd_ineligible$failure_reason == reason]),
                patients))
  }
}

# Blood samples by classification (original diagnostic)
blood_samples_summary <- samples_final %>%
  filter(sample_type == "Peripheral blood sample") %>%
  group_by(sample_class) %>%
  summarise(count = n(), .groups = 'drop')

cat("\n===== BLOOD SAMPLES CLASSIFICATION SUMMARY =====\n")
print(blood_samples_summary)

# Blood samples that are marked "Other" - potentially missed
blood_other <- samples_final %>%
  filter(sample_type == "Peripheral blood sample" & sample_class == "Other") %>%
  select(PATIENT_ID, sample_date, sample_type, days_from_diagnosis, nearest_abs_diff_mrd)

if (nrow(blood_other) > 0) {
  cat("\n⚠️  WARNING: Peripheral blood samples classified as 'Other':\n")
  print(blood_other)
  cat("\nThese may be MRD blood samples if ±14 day threshold is too strict.\n")
} else {
  cat("\n✓ No peripheral blood samples classified as 'Other'.\n")
}

# MRD samples by type
mrd_by_type <- samples_final %>%
  filter(sample_class == "MRD") %>%
  group_by(sample_type) %>%
  summarise(count = n(), .groups = 'drop')

cat("\nMRD samples by type:\n")
print(mrd_by_type)

cat("\n===== END DIAGNOSTIC REPORT =====\n\n")

#-------------------------------------------------------------------
# 9. Output results for both cohorts
#-------------------------------------------------------------------

# Cohort 1: Any MRD (BM or blood)
final_output_any_mrd <- eligible_samples_any_mrd %>%
  select(
    sample_id,
    PATIENT_ID,
    STUDY_NAME,
    sample_date,
    sample_type,
    sample_class,
    days_from_dx,
    days_from_mrd,
    tx_lines,
    tx_regimens,
    progressed
  )

write.csv(final_output_any_mrd, file = "eligible_samples_ANY_MRD.csv", row.names = FALSE)
cat("✓ Results written to 'eligible_samples_ANY_MRD.csv'\n")

# Cohort 2: Blood MRD only
final_output_mrd_blood <- eligible_samples_mrd_blood_only %>%
  select(
    sample_id,
    PATIENT_ID,
    STUDY_NAME,
    sample_date,
    sample_type,
    sample_class,
    days_from_dx,
    days_from_mrd,
    tx_lines,
    tx_regimens,
    progressed
  )

write.csv(final_output_mrd_blood, file = "eligible_samples_BLOOD_MRD_ONLY.csv", row.names = FALSE)
cat("✓ Results written to 'eligible_samples_BLOOD_MRD_ONLY.csv'\n\n")

#-------------------------------------------------------------------
# Summary statistics for both cohorts
#-------------------------------------------------------------------

cat("========== COHORT 1: ANY MRD (BM or Blood) ==========\n\n")
cat("Total eligible patients:", n_distinct(final_output_any_mrd$PATIENT_ID), "\n")
cat("Total eligible samples:", nrow(final_output_any_mrd), "\n\n")

class_counts_any <- final_output_any_mrd %>% count(sample_class)
cat("Sample classification breakdown:\n")
print(class_counts_any)

sample_type_counts_any <- final_output_any_mrd %>% count(sample_type)
cat("\nSample type breakdown:\n")
print(sample_type_counts_any)

# Identify ineligible blood MRD patients with multiple timepoints
multi_timepoint_patients <- blood_mrd_ineligible %>%
  filter(n_timepoints >= 2) %>%
  pull(PATIENT_ID)

cat("Patients with blood MRD but no baseline AND multiple timepoints:", length(multi_timepoint_patients), "\n\n")


cat("\n========== COHORT 2: BLOOD MRD ONLY ==========\n\n")
cat("Total eligible patients:", n_distinct(final_output_mrd_blood$PATIENT_ID), "\n")
cat("Total eligible samples:", nrow(final_output_mrd_blood), "\n\n")

class_counts_blood <- final_output_mrd_blood %>% count(sample_class)
cat("Sample classification breakdown:\n")
print(class_counts_blood)

sample_type_counts_blood <- final_output_mrd_blood %>% count(sample_type)
cat("\nSample type breakdown:\n")
print(sample_type_counts_blood)

cat("\n========== SUMMARY ==========\n")
cat("Patients with baseline + ANY MRD:", n_distinct(final_output_any_mrd$PATIENT_ID), "\n")
cat("Patients with baseline + BLOOD MRD:", n_distinct(final_output_mrd_blood$PATIENT_ID), "\n")
cat("Patients needing blood MRD collection:", 
    n_distinct(final_output_any_mrd$PATIENT_ID) - n_distinct(final_output_mrd_blood$PATIENT_ID), "\n\n")

cat("========== BLOOD MRD EXPANSION OPPORTUNITY ==========\n")
cat("Total patients with blood MRD samples:", length(patients_with_blood_mrd), "\n")
cat("  - Already eligible (have baseline):", nrow(blood_mrd_eligible), "\n")
cat("  - Need baseline collection:", nrow(blood_mrd_ineligible), "\n\n")

if (nrow(blood_mrd_ineligible) > 0) {
  multi_tp_count <- length(multi_timepoint_patients)
  single_tp_count <- nrow(blood_mrd_ineligible) - multi_tp_count
  
  cat("========== SUMMARY REPORT FOR COHORT EXPANSION ==========\n\n")
  
  cat("COHORT STATUS:\n")
  cat(sprintf("- Primary cohort (baseline + blood MRD): %d patients, %d samples\n",
              n_distinct(final_output_mrd_blood$PATIENT_ID),
              nrow(final_output_mrd_blood)))
  
  cat(sprintf("- Extended cohort (baseline + any MRD): %d patients, %d samples\n",
              n_distinct(final_output_any_mrd$PATIENT_ID),
              nrow(final_output_any_mrd)))
  
  cat("\nBLOOD MRD EXPANSION TARGETS:\n")
  cat(sprintf("- Total patients with blood MRD: %d\n", length(patients_with_blood_mrd)))
  cat(sprintf("  ├─ Already eligible: %d patients (have baseline)\n", nrow(blood_mrd_eligible)))
  cat(sprintf("  └─ Missing baseline: %d patients\n", nrow(blood_mrd_ineligible)))
  
  if (multi_tp_count > 0) {
    cat(sprintf("\nMULTI-TIMEPOINT OPPORTUNITIES:\n"))
    cat(sprintf("- %d patients with blood MRD AND multiple sample timepoints (but no baseline)\n", multi_tp_count))
    cat("  These patients have additional samples at other timepoints beyond baseline/MRD,\n")
    cat("  suggesting rich longitudinal data that could be valuable if baseline samples are obtained.\n")
    cat("  See 'blood_MRD_multi_timepoint_expansion_targets.csv' for details.\n\n")
  }
  
  cat("RECOMMENDATIONS:\n")
  cat("1. Prioritize obtaining baseline samples for the ", nrow(blood_mrd_ineligible), " patients with blood MRD\n")
  if (multi_tp_count > 0) {
    cat("2. Give priority to the ", multi_tp_count, " multi-timepoint patients (IDs: ",
        paste(multi_timepoint_patients, collapse = ", "), ")\n")
  }
  cat("3. Review blood_MRD_patient_analysis.csv for complete diagnostic information\n\n")
}


if (length(multi_timepoint_patients) > 0) {
  # Get detailed sample info for these patients
  multi_timepoint_details <- samples_final %>%
    filter(PATIENT_ID %in% multi_timepoint_patients) %>%
    select(
      PATIENT_ID, STUDY_NAME, sample_type, sample_date, sample_class,
      days_from_diagnosis, baseline_flag, nearest_abs_diff_mrd, nearest_mrd_diff,
      mrd_flag, in_treatment, progressed
    ) %>%
    arrange(PATIENT_ID, sample_date)
  
  # Export to CSV
  write.csv(multi_timepoint_details, file = "blood_MRD_multi_timepoint_expansion_targets.csv", row.names = FALSE)
  cat("✓ Detailed sample info written to 'blood_MRD_multi_timepoint_expansion_targets.csv'\n\n")
  
  # Display in console
  cat("Sample details by patient:\n\n")
  for (pid in multi_timepoint_patients) {
    patient_samples <- multi_timepoint_details %>% filter(PATIENT_ID == pid)
    cat(sprintf("--- Patient %s ---\n", pid))
    cat(sprintf("Total samples: %d | Timepoints: %d\n", nrow(patient_samples), n_distinct(patient_samples$sample_class)))
    print(patient_samples %>% select(sample_type, sample_date, sample_class, days_from_diagnosis, mrd_flag, in_treatment))
    cat("\n")
  }
  
  # Summary of what timepoints they have
  timepoint_summary <- multi_timepoint_details %>%
    group_by(PATIENT_ID) %>%
    summarise(
      total_samples = n(),
      timepoints = paste(unique(sample_class), collapse = ", "),
      has_mrd = any(sample_class == "MRD"),
      has_treatment = any(sample_class == "Treatment"),
      has_other = any(sample_class == "Other"),
      .groups = 'drop'
    )
  
  cat("\n----- Timepoint Breakdown -----\n")
  print(timepoint_summary)
} else {
  cat("No patients with multiple timepoints found.\n")
}

cat("\n===== END DETAILED ANALYSIS =====\n\n")

cat("\n✓ Script completed successfully!\n")

