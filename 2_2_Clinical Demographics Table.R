# =============================================================================
# Script: generate_baseline_characteristics_table.R
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
#   • RDS: Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.rds
#
# Outputs:
#   • Word: table1_categorical_updated_final.docx
#   • Word: baseline_characteristics_updated.docx
#   • TXT/RDS: cohort_assignment_table.{txt,rds}
#
# Dependencies:
#   library(tidyverse)
#   library(gtsummary)
#   library(officer)
#   library(flextable)
#
# Usage:
#   source("generate_baseline_characteristics_table.R")
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================

# -----------------------------------------------------------
# 0.  PACKAGES  (install once, then keep only library() lines)
# -----------------------------------------------------------
# install.packages(c("tidyverse", "gtsummary", "gt", "officer"))   # ← run once
library(tidyverse)
library(gtsummary)
library(officer)
library(flextable)


# -----------------------------------------------------------
# 1.  DATA  
# -----------------------------------------------------------
# keep ONLY the earliest sample per patient
dat <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.rds")


### Get the list of qualifying samples 
# Identify patients with ≥2 rows where either zscore_bm or zscore_blood is not NA
patients_have_data_on <- dat %>%
  filter(!is.na(zscore_BM) | !is.na(zscore_blood)) %>%  # keep rows with data in at least one of the two columns
  group_by(Patient) %>%
  filter(n() >= 2) %>%                                    # require at least 2 such rows per patient
  distinct(Patient) %>%                                   # one entry per patient
  pull(Patient)                                           # extract as a simple vector


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
  filter(Patient %in% patients_have_data_on)

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
# 2.  DEFINE COHORTS  ------------
# -----------------------------------------------------------
# Define the two groups
new_dx_ids <- patients_have_data_on[
  grepl("^[A-Z]{2}-\\d{2}$", patients_have_data_on) |
    patients_have_data_on %in% c("IMG-060", "IMG-127", "IMG-181")
]

pretreated_ids <- patients_have_data_on[
  grepl("SPORE", patients_have_data_on) |
    patients_have_data_on %in% c("IMG-098", "IMG-159")
]

dat_base <- dat_base %>%
  mutate(cohort = case_when(
    Patient %in% new_dx_ids     ~ "Newly diagnosed",
    Patient %in% pretreated_ids ~ "Pre‑treated",
    TRUE                        ~ NA_character_
  ))

# -----------------------------------------------------------
# 3.  VARIABLES TO SHOW IN TABLE 1 --------------------------
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

### Add the risk if it is NA 
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


# 2. Build the categorical-only Table 1
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

# 3. (Optional) Export to Word
tbl1_flex <- as_flex_table(tbl1_cat)
doc <- read_docx() %>%
  body_add_flextable(tbl1_flex) %>%
  body_end_section_portrait()
print(doc, target = "table1_categorical_updated_final.docx")

# -----------------------------------------------------------
# 6.  DONE!  -----------------------------------------------
# -----------------------------------------------------------
# The file 'baseline_characteristics.docx' is now in your working directory.
# Open it in Word to fine‑tune column widths, font, or add footnotes.



### Export important things 
# Create a cohort assignment dataframe
cohort_df <- data.frame(
  Patient = c(new_dx_ids, pretreated_ids),
  NewDx = c(rep(1, length(new_dx_ids)), rep(0, length(pretreated_ids)))
)

# Optional: sort for clarity
cohort_df <- cohort_df[order(cohort_df$Patient), ]

# Export to file (choose your desired path)
write.table(cohort_df, file = "cohort_assignment_table.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Also optionally save as RDS
saveRDS(cohort_df, file = "cohort_assignment_table.rds")







##### Below here is testing 
### Removed most of these since didn't have info for them
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
  mutate(across(all_of(vars_cat), ~ factor(.x, exclude = NULL)))


# -----------------------------------------------------------
# 4.  BUILD TABLE 1  ---------------------------------------
# -----------------------------------------------------------
tbl1 <- dat_base %>%
  select(all_of(c(vars_cat, "cohort"))) %>%
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

