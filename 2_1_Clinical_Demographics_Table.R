# =============================================================================
# Script: 2_1_Clinical_Demographics_Table.R
#
# Description:
#   This script processes the final aggregated cfWGS feature table (with clinical
#   and demographic annotations), selects each patient’s earliest Baseline/Diagnosis
#   sample, defines clinical cohorts (Newly diagnosed vs Pre-treated), and builds
#   a “Table 1” of baseline characteristics.  It also exports a cohort assignment
#   file for downstream merging.
#
# Steps:
#   1. Load the aggregated RDS of cfWGS features + clinical/demographics
#   2. Identify patients with both BM and blood data
#   3. Subset to Diagnosis/Baseline, resolve duplicates (CA-02, SPORE_0009, etc.)
#   4. Define cohorts by patient ID patterns (Newly diagnosed vs Pre-treated)
#   5. Clean and recode categorical and continuous variables
#   6. Build categorical Table 1 with gtsummary + overall column
#   7. (Optional) Build continuous Table 1
#   8. Export tables to Word and save cohort assignment as TXT/RDS
#
# Inputs:
#   • RDS: Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.rds
#   • CSV: Output_tables_2025/patient_cohort_assignment.csv
#   • RDS: cohort_assignment_table_updated.rds
#
# Outputs:
#   • Word: table1_clinical_categorical_updated_final_3.docx
#   • Optional Word helper: baseline_characteristics_updated.docx
#   • Optional exploratory cohort/fragmentomics audit output printed to console
#
# Dependencies:
#   library(tidyverse)
#   library(gtsummary)
#   library(officer)
#   library(flextable)
#
# Usage:
#   Rscript Scripts_2025/Final_Scripts/2_1_Clinical_Demographics_Table.R
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/2_1_Clinical_Demographics_Table.R
#
# Manuscript outputs created/updated:
#   - Table 1: baseline clinical/demographic characteristics, exported as the
#     generated Word source and, when present, the manually edited final DOCX/PDF
#     used for the manuscript.
#
# Pipeline role:
#   Table 1 is built from the baseline/diagnosis sample per patient after
#   resolving known duplicate baseline records. The cohort assignment exported
#   here is reused by later scripts so clinical grouping is consistent across
#   figures, supplementary tables, and survival analyses.
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

# -----------------------------------------------------------
# 0.  PACKAGES  (install once, then keep only library() lines)
# -----------------------------------------------------------
# install.packages(c("tidyverse", "gtsummary", "gt", "officer"))   # ← run once
library(tidyverse)
library(gtsummary)
library(officer)
library(flextable)

# Shared manuscript-output helpers.
# These copy the exact final Table 1 artifacts into final_manuscript_objects/
# while retaining this script as the place where the Table 1 source is built.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)


# -----------------------------------------------------------
# 1. Load baseline table inputs
# -----------------------------------------------------------
# The Table 1 source table is intentionally generated from the same historical
# aggregate RDS used by the original script. See the input-version note above
# before replacing this with a newer aggregate table.
dat <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.rds")


### Load the curated patient cohort list
# `patient_cohort_assignment.csv` defines which patients enter Table 1 and the
# cohort grouping shown as columns. The RDS cohort file is retained for the
# exploratory audit block at the end of the script.
# Define export directory
export_dir <- "Output_tables_2025"

# Load from CSV (if you want to view/edit easily)
patient_cohort_tbl_csv <- read.csv(file.path(export_dir, "patient_cohort_assignment.csv"))

cohort_df <- readRDS("cohort_assignment_table_updated.rds")


### Select baseline/diagnosis rows
# Table 1 is a patient-level baseline table. This section keeps only baseline
# or diagnosis rows and then resolves known duplicate baseline records.
# 1. Subset to only Diagnosis or Baseline timepoints
dat_tb <- dat %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline"))


# 2. Check for duplicate rows per patient in this subset
dup_patients <- dat_tb %>%
  count(Patient) %>%
  filter(n > 1) %>%
  pull(Patient)

# ensure Date is Date class
dat_tb <- dat_tb %>%
  mutate(Date = as.Date(Date))

# 1) Remove CA-03 timepoint 02
dat_tb2 <- dat_tb %>%
  filter(!(Patient == "CA-03" & timepoint_info == "02"))

# 2) Consolidate the two CA-02 rows
resp_CA02 <- dat_tb2 %>%
  filter(Patient == "CA-02") %>%
  # order so timepoint “01” comes before “02”
  arrange(factor(timepoint_info, levels = c("01","02"))) %>%
  summarise(across(everything(), ~{
    vals <- .
    # first non-NA in order
    first_val <- vals[which(!is.na(vals))[1]]
    # any other non-NA ≠ first_val
    other    <- vals[!is.na(vals) & vals != first_val][1]
    # if first_val is Unknown/Other but other is a “real” value, use other
    if (!is.na(first_val) &&
        first_val %in% c("Unknown","Other") &&
        !is.na(other) &&
        !other %in% c("Unknown","Other")) {
      other
    } else {
      first_val
    }
  }))

# 3) Drop all CA-02 originals
dat_tb3 <- dat_tb2 %>% filter(Patient != "CA-02")

# 4) For SPORE_0009, keep only the 2016-08-24 Baseline row
dat_tb4 <- dat_tb3 %>%
  filter(!(Patient == "SPORE_0009" &
             !(Date == as.Date("2016-08-24") & timepoint_info == "Baseline")))

# 5) Re-bind the collapsed CA-02 row
dat_tb_final <- bind_rows(dat_tb4, resp_CA02) %>%
  arrange(Patient)

dat_base <- dat_tb_final %>%
  filter(Patient %in% patient_cohort_tbl_csv$Patient)

## Remove dup
dat_base <- dat_base %>%
  filter(!(Patient == "CA-03" & Timepoint == "02"))

# 6) Quick duplicate check
dups <- dat_base %>%
  count(Patient, timepoint_info) %>%
  filter(n > 1)

if (nrow(dups)) {
  warning("Still multiple Diagnosis/Baseline rows for: ",
          paste(unique(dups$Patient), collapse = ", "))
} else {
  message("All patients now have at most one Diagnosis/Baseline row.")
}


# -----------------------------------------------------------
# 2. Define cohorts for Table 1
# -----------------------------------------------------------
dat_base <- dat_base %>%
  left_join(patient_cohort_tbl_csv, by = "Patient")


dat_base <- dat_base %>%
  mutate(cohort = case_when(
    Cohort == "Frontline"     ~ "Frontline induction-transplant",
    TRUE                      ~ Cohort
  ))


# -----------------------------------------------------------
# 3. Recode variables shown in Table 1
# -----------------------------------------------------------
# ------------------------- 2.1  Convert continuous vars to numeric
# 1. Collapse Subtype
dat_base <- dat_base %>%
  mutate(Subtype = fct_collapse(
    Subtype,
    IgG              = c("IgG Kappa", "IgG Lambda", "IgG (Subtype Unknown)"),
    IgA              = c("IgA Kappa", "IgA Lambda"),
    `Light Chain Only` = "Light Chain Only (Kappa)",
    Other            = c("Other", NA),
    .default         = "Other"
  ))

# Fix issue 
dat_base <- dat_base %>%
  mutate(Subtype = na_if(as.character(Subtype), ".default")) %>%
  mutate(Subtype = factor(Subtype))  # re-factor after conversion


# 2. Define your final var lists
vars_cat <- c("Gender", "AGE_GROUP", "ISS_STAGE",
              "Cytogenetic_Risk", "Subtype", "ECOG_SCORE")

### Most things NA, so trying best with what isn't NA for everything at least across both groups 
## Plasma cll percent NA for most

# 2. Recode “Unknown” and “Unknown/Missing” to NA and drop them
dat_base <- dat_base %>%
  # 1) Turn all your cats into character
  mutate(across(all_of(vars_cat), as.character)) %>%
  
  # 2) Replace "Unknown" or "Unknown/Missing" with NA
  mutate(across(
    all_of(vars_cat),
    ~ if_else(.x %in% c("Unknown", "Unknown/Missing"), NA_character_, .x)
  )) %>%
  
  # 3) Back to factor and drop any empty levels
  mutate(across(all_of(vars_cat), as.factor)) %>%
  mutate(across(all_of(vars_cat), fct_drop))

### Fill missing cytogenetic risk where FISH markers are sufficient
# This does not overwrite an existing risk label. It only derives risk for rows
# where Cytogenetic_Risk is missing and the high-risk FISH markers are present.
dat_base <- dat_base %>%
  mutate(
    # 1. Calculate risk from FISH results
    Cytogenetic_Risk_calc = case_when(
      # Any high-risk marker positive → High Risk
      T_4_14  == "Positive" |
        T_14_16 == "Positive" |
        DEL_17P == "Positive"    ~ "High Risk",
      
      # All three markers explicitly NEGATIVE → Standard Risk
      !is.na(T_4_14)  & T_4_14  == "Negative" &
        !is.na(T_14_16) & T_14_16 == "Negative" &
        !is.na(DEL_17P) & DEL_17P == "Negative" ~ "Standard Risk",
      
      # Otherwise (any missing or contradictory data) → remain NA
      TRUE                      ~ NA_character_
    ),
    
    # 2. Fill only the NA entries in your original Cytogenetic_Risk
    Cytogenetic_Risk = coalesce(Cytogenetic_Risk, Cytogenetic_Risk_calc),
    
    # 3. (re)factor with the two levels
    Cytogenetic_Risk = factor(
      Cytogenetic_Risk,
      levels = c("Standard Risk", "High Risk")
    )
  ) %>%
  select(-Cytogenetic_Risk_calc)   # drop the helper column if you like

### Export deterministic Table 1 source-data companions
# These CSV/TSV files are not the formatted manuscript table. They are a
# reviewer/developer-friendly audit trail for the categorical counts and
# denominators that feed the gtsummary Table 1 object below. Keeping them in the
# original Table 1 script makes the final_manuscript_objects source-data folder
# reproducible without relying on the separate reproducible_workflow generator.
format_table1_count_percent <- function(n, denominator) {
  if (is.na(n) || is.na(denominator) || denominator == 0) {
    return("0 (0%)")
  }
  paste0(n, " (", round(100 * n / denominator), "%)")
}

table1_labels <- c(
  Gender = "Gender",
  AGE_GROUP = "Age Group",
  ISS_STAGE = "ISS Stage",
  Cytogenetic_Risk = "Cytogenetic Risk",
  Subtype = "Myeloma Ig Subtype",
  ECOG_SCORE = "ECOG Performance Status"
)
table1_cohort_levels <- c("Frontline induction-transplant", "Non-frontline")

table1_source_rows <- list()
for (variable in vars_cat) {
  values <- as.character(dat_base[[variable]])
  display_values <- ifelse(is.na(values), "(Missing)", values)
  levels_to_show <- unique(display_values)
  levels_to_show <- levels_to_show[order(levels_to_show == "(Missing)", levels_to_show)]

  for (level in levels_to_show) {
    source_row <- data.frame(
      variable = unname(table1_labels[[variable]]),
      level = level,
      stringsAsFactors = FALSE
    )
    for (cohort_level in table1_cohort_levels) {
      in_cohort <- dat_base$cohort == cohort_level
      n_level <- sum(in_cohort & display_values == level, na.rm = TRUE)
      denominator <- sum(in_cohort, na.rm = TRUE)
      source_row[[cohort_level]] <- format_table1_count_percent(n_level, denominator)
    }
    source_row[["Total"]] <- format_table1_count_percent(
      sum(display_values == level, na.rm = TRUE),
      nrow(dat_base)
    )
    table1_source_rows[[length(table1_source_rows) + 1]] <- source_row
  }
}

table1_source_counts <- bind_rows(table1_source_rows)
dir.create("Output_tables_2025", recursive = TRUE, showWarnings = FALSE)
table1_source_counts_path <- file.path(
  "Output_tables_2025",
  "Table_1_clinical_demographics_computed_source_counts.csv"
)
write.csv(table1_source_counts, table1_source_counts_path, row.names = FALSE, quote = TRUE)

table1_source_qc <- data.frame(
  input_rds = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.rds",
  cohort_csv = "Output_tables_2025/patient_cohort_assignment.csv",
  n_baseline_table_patients = nrow(dat_base),
  n_frontline = sum(dat_base$cohort == "Frontline induction-transplant", na.rm = TRUE),
  n_non_frontline = sum(dat_base$cohort == "Non-frontline", na.rm = TRUE),
  output_csv = table1_source_counts_path,
  stringsAsFactors = FALSE
)
table1_source_qc_path <- file.path(
  "Output_tables_2025",
  "Table_1_clinical_demographics_computed_source_qc.tsv"
)
write.table(
  table1_source_qc,
  file = table1_source_qc_path,
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

ms_copy_artifact(
  source_path = table1_source_counts_path,
  artifact_id = "TABLE1_DOCX",
  role = "main_table_source_counts_csv",
  description = "Table 1 computed categorical count/percent source data used to build the gtsummary table.",
  script_name = "2_1_Clinical_Demographics_Table.R"
)

ms_copy_artifact(
  source_path = table1_source_qc_path,
  artifact_id = "TABLE1_DOCX",
  role = "main_table_source_qc_tsv",
  description = "Table 1 source-data QC summary with input paths and cohort denominators.",
  script_name = "2_1_Clinical_Demographics_Table.R"
)

### Build and export categorical manuscript Table 1 source
tbl1_cat <- dat_base %>%
  select(all_of(c(vars_cat, "cohort"))) %>%
  tbl_summary(
    by           = cohort,
    missing      = "ifany",                # show NA counts
    missing_text = "(Missing)",
    statistic    = all_categorical() ~ "{n} ({p}%)",
    label        = list(
      Gender            ~ "Gender",
      AGE_GROUP         ~ "Age Group",
      ISS_STAGE         ~ "ISS Stage",
      Cytogenetic_Risk  ~ "Cytogenetic Risk",
      Subtype           ~ "Myeloma Ig Subtype",
      ECOG_SCORE        ~ "ECOG Performance Status"
    )
  ) %>%
  add_overall(last = TRUE, col_label = "Total") %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

# Export the scripted DOCX source. The final manuscript DOCX/PDF in the output
# map appear to be edited/exported versions of this source table.
tbl1_flex <- as_flex_table(tbl1_cat)
doc <- read_docx() %>%
  body_add_flextable(tbl1_flex) %>%
  body_end_section_portrait()
table1_generated_docx <- "table1_clinical_categorical_updated_final_3.docx"
print(doc, target = table1_generated_docx)

# MANUSCRIPT OUTPUT: Table 1
# This script regenerates the categorical Table 1 source DOCX above. The mapped
# final manuscript artifacts are the manually DA-edited DOCX and exported PDF in
# Final Tables and Figures/. Those frozen files are copied when present because
# they are the exact files used in the current manuscript. If the edited DOCX is
# absent, we still export the regenerated source DOCX so a command-line run has
# a usable Table 1 artifact, but the manifest will make the source filename clear.
table1_final_docx <- file.path(
  "Final Tables and Figures",
  "table1_clinical_categorical_updated_final_3_DA_edited.docx"
)
table1_final_pdf <- file.path(
  "Final Tables and Figures",
  "table1_clinical_categorical_updated_final_3_DA_edited.pdf"
)

table1_docx_to_copy <- if (file.exists(table1_final_docx)) {
  table1_final_docx
} else {
  warning(
    "Final edited Table 1 DOCX was not found; copying regenerated source DOCX instead: ",
    table1_generated_docx
  )
  table1_generated_docx
}

ms_copy_artifact(
  source_path = table1_docx_to_copy,
  artifact_id = "TABLE1_DOCX",
  role = "main_table_docx",
  description = "Table 1: clinical and demographic baseline characteristics. Edited manuscript DOCX is preferred when present; otherwise the regenerated source DOCX is copied.",
  script_name = "2_1_Clinical_Demographics_Table.R"
)

if (file.exists(table1_final_pdf)) {
  ms_copy_artifact(
    source_path = table1_final_pdf,
    artifact_id = "TABLE1_PDF",
    role = "main_table_pdf",
    description = "Table 1: edited/exported manuscript PDF corresponding to the final DOCX table.",
    script_name = "2_1_Clinical_Demographics_Table.R"
  )
} else {
  warning(
    "Final edited Table 1 PDF was not found. The PDF is a manuscript/export artifact and is not regenerated by this R script: ",
    table1_final_pdf
  )
}

### Optional continuous Table 1 companion
# This block is retained from the original script for traceability, but it is
# not the mapped final manuscript Table 1. Set the flag below to TRUE only when
# you intentionally want to regenerate the historical continuous-variable helper
# table (`baseline_characteristics_updated.docx`).
export_exploratory_continuous_table1 <- FALSE

if (isTRUE(export_exploratory_continuous_table1)) {
# ------------------------- 2.1  Convert continuous vars to numeric
vars_cont <- c(
  "AGE", "dFLC", "Albumin", "B2_micro", "Calcium",
  "Creatinine", "Hemoglobin", "LDH", "Plasma_pct"
)

dat_base <- dat_base %>%
  mutate(across(all_of(vars_cont), as.numeric))

# ------------------------- 2.2  Make categorical vars factors
vars_cat <- c(
  "Gender", "AGE_GROUP", "ISS_STAGE", "R_ISS_STAGE",
  "Cytogenetic_Risk", "Subtype", "ECOG_SCORE", "KPS_SCORE"
)

dat_base <- dat_base %>%
  mutate(across(all_of(vars_cat), ~ factor(.x)))


# -----------------------------------------------------------
# 4.  BUILD TABLE 1  ---------------------------------------
# -----------------------------------------------------------
tbl1 <- dat_base %>%
  select(all_of(c(vars_cat, vars_cont, "cohort"))) %>%
  tbl_summary(
    by        = cohort,
    type      = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",
    digits    = all_continuous() ~ 1,
    missing   = "always",                 # <‑‑ show missing counts explicitly
    label     = list(
      AGE        ~ "Age (years)",
      dFLC       ~ "dFLC (mg/L)",
      B2_micro   ~ "β₂‑microglobulin (mg/L)",
      Plasma_pct ~ "BM plasma cells (%)"
    )
  ) %>%
  add_overall(last = TRUE, col_label = "Total") %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

# -----------------------------------------------------------
# 5.  EXPORT TO WORD  --------------------------------------
# -----------------------------------------------------------
# Convert to flextable
tbl1_flex <- as_flex_table(tbl1)

# Export using flextable functions
doc <- read_docx() %>%
  body_add_flextable(tbl1_flex) %>%
  body_end_section_portrait()

print(doc, target = "baseline_characteristics_updated.docx")
}



### Exploratory fragmentomics/clinical MRD audit
# This optional audit does not create a manuscript figure or table. It lists
# possible additional non-IMG/non-SPORE patients with fragmentomics plus
# clinical MRD data who are not in the current cohort assignment file. Keep
# FALSE for the manuscript pipeline; set TRUE during cohort-expansion review.
run_fragmentomics_clinical_mrd_audit <- FALSE

if (isTRUE(run_fragmentomics_clinical_mrd_audit)) {

# Define which columns count as “clinical MRD”
clinical_mrd_cols <- c(
  "MRD_by_clinical_testing",
  "MRD_by_clinical_testing_stringent",
  "MRD_Clinical_Binary",
  "MRD_Clinical_Stringent_Binary"
)

patients_with_FS_and_clinical_MRD <- dat %>%
  # 1) require FS not missing
  filter(!is.na(FS)) %>%
  # 2) require at least one clinical‐MRD column not missing
  filter(if_any(all_of(clinical_mrd_cols), ~ !is.na(.))) %>%
  # 3) Check not at baseline 
  filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  # 4) pull unique patient IDs
  distinct(Patient) %>%
  pull(Patient)

# Print the result
cat("Patients with FS and clinical MRD data:\n")
print(patients_with_FS_and_clinical_MRD)
cat("Total:", length(patients_with_FS_and_clinical_MRD), "patients\n")

# Get patients that do not start with IMG or SPORE
non_img_spore_patients <- patients_with_FS_and_clinical_MRD[
  !grepl("^(IMG|SPORE)", patients_with_FS_and_clinical_MRD)
]

# Which of these are NOT in the cohort_df Patient list?
new_possible_patients <- setdiff(non_img_spore_patients, cohort_df$Patient)

# Show them
cat("Patients with FS+clinical MRD, not in cohort_df, and not IMG/SPORE:\n")
print(new_possible_patients)
cat("Total:", length(new_possible_patients), "patients\n")
}


### Archived idea: excluded from the manuscript workflow
# The original working script had considered adding more complexity after this
# point, but that path is excluded from all final manuscript figure/table outputs.
