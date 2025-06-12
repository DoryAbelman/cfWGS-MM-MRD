# ──────────────────────────────────────────────────────────────────────────────
#  process_clinical_metadata.R
#
#  Purpose:
#    • Ingest and harmonize multiple clinical/extract metadata sources (SPORE, M4, IMMAGINE, MyC).
#    • Clean and normalize sample identifiers, sample types, and collection dates.
#    • Annotate each sample with a unified timepoint_info (Diagnosis, Baseline, Post-induction,
#      Post-transplant, Maintenance, Relapse/Progression, etc.).
#    • Integrate FISH flags, lab timepoints, relapse/progression dates, and OS follow-up data.
#    • Compute time-to-relapse (absolute and non-absolute) and maximum follow-up durations for PFS analyses.
#    • Generate summary tables (counts by study, sample type, timepoint) and export QC plots.
#    • Export intermediary and final “combined_clinical_data_updated” master table along with
#      derived tables: baseline dates, PFS_days, latest_dates_per_patient, patient_counts, BAM lists, etc.
#
#  Author:        Dory Abelman
#  Last updated:  2025-06-06
#
#  Notes:
#    • All file paths are assumed relative to your project root.
#    • Before running, ensure these inputs exist under “Clinical data” and “M4_CMRG_Data” directories:
#         – cfWGS_MS_Integrated_Metadata.xlsx
#         – SPORE: Spore_bams_with_collected_dates.txt,
#                  Extracted_clinical_data_cfDNA_project.xlsx,
#                  Spore_timepoint_info_updated.xlsx,
#                  SPORE_OS_info.xlsx
#         – M4: M4_COHORT_STAGING.xlsx,
#               M4_COHORT_STEM_CELL_TRANSPLANT.xlsx,
#               M4_COHORT_LABS.xlsx,
#               Relapse_dates.xlsx,
#               M4 Status and best response update Dec 2 24 just dates.xlsx,
#               M4_COHORT_DEMO.xlsx,
#               M4_COHORT_CHEMOTHERAPY.xlsx,
#               TFRIM4_Processing Log_Nov2024.xlsx
#         – IMMAGINE: IMG_request_20241009 (2).xlsx,
#                     Cleaned transplant dates just dates.xlsx,
#                     Extracted_clinical_MRD_data.xlsx,
#                     Cleaned_Patient_Follow-Up_Table_IMMAGINE.csv,
#                     NAs in timepoint info - IMMAGINE additional timepoint info.csv
#
#    • Outputs:
#         – “combined_clinical_data_updated_April2025.csv” (master metadata table)
#         – Various exported tables under “Exported_data_tables_clinical” and “Oct2024_exported_data”
#         – Summary CSVs: patient_counts_summary.csv, patient_counts_detailed.csv, BAM lists, etc.
#         – QC/barplots under “plots_Oct”#
#
#   *This script is intended to be a one‐shot cleanup for downstream analysis.*
# ──────────────────────────────────────────────────────────────────────────────



# ─── 1.  Load libraries ─────────────────────────────────────────────────────────
library(readxl)       # reading Excel files
library(data.table)   # fast fread/fwrite if needed
library(tidyverse)    # dplyr, tidyr, ggplot2, etc.
library(lubridate)    # date parsing


# ─── 2.  Configuration & helper functions (optional) ───────────────────────────
# [ Single place to change paths ]
#
# CFG <- list(
#   root        = here::here(),
#   out_tables  = "out/tables",
#   out_plots   = "out/plots",
#   match_window = 60   # use when matching SPORE sample dates to clinical dates
# )
# 
# dir.create(CFG$out_tables, showWarnings = FALSE, recursive = TRUE)
# dir.create(CFG$out_plots,  showWarnings = FALSE, recursive = TRUE)
# 
# pth <- function(...) fs::path(CFG$root, ...)
# write_out <- function(df, name) {
#   write_csv(df, pth(CFG$out_tables, paste0(name, ".csv")))
# }


# ─── 3.  Load raw clinical files ────────────────────────────────────────────────
##  3.1  cf/WGS integrated metadata (M4, MyC etc.)
clinical_data <- read_excel("cfWGS_MS_Integrated_Metadata.xlsx", sheet = 1)

##  3.2  SPORE “bams with collected dates” spreadsheet
spore_data <- read_tsv("Clinical data/SPORE/Spore_bams_with_collected_dates.txt")

##  3.3  M4 staging (diagnosis dates)
M4_diagnosis_date      <- read_excel("Clinical data/M4/M4_COHORT_STAGING.xlsx")

##  3.4  SPORE relapse/extracted info (long list of diagnosis/relapse dates)
spore_data_extracted   <- read_excel("Clinical data/SPORE/Extracted_clinical_data_cfDNA_project.xlsx")

##  3.5  SPORE timepoint info (dates of interest)
spore_timepoint_info   <- read_excel("Clinical data/SPORE/Spore_timepoint_info_updated.xlsx")

##  3.6  SPORE overall survival info (for later OS annotations)
spore_OS_info          <- read_excel("Clinical data/SPORE/SPORE_OS_info.xlsx")


# ─── 4.  Tidy SPORE “bams” table ────────────────────────────────────────────────
##  4.1  Create Sample_type from “Bams have” patterns
spore_data <- spore_data %>%
  mutate(
    Sample_type = case_when(
      grepl("Bm_T", `Bams have`) ~ "BM_cells",
      grepl("Pl_T", `Bams have`) ~ "Blood_plasma_cfDNA",
      grepl("Pb_R", `Bams have`) ~ "Blood_Buffy_coat",
      TRUE ~ NA_character_
    )
  )

##  4.2  Sort by Patient & date to compute “Timepoint” (≥14 days → new timepoint)
spore_data <- spore_data %>%
  arrange(Patient, Date_of_sample_collection) %>%
  group_by(Patient) %>%
  mutate(
    Timepoint = {
      date_diff      <- c(0, as.numeric(diff(Date_of_sample_collection)))
      timepoint_grp  <- cumsum(date_diff > 14) + 1
      as.character(timepoint_grp)
    }
  ) %>%
  ungroup()

##  4.3  Add “Study = SPORE” + unified “Sample_ID” and rename “Bams have” → “Bam”
spore_data <- spore_data %>%
  mutate(
    Study     = "SPORE",
    Sample_ID = paste0(Patient, "_T", Timepoint, "_", Sample_type)
  ) %>%
  rename(Bam = `Bams have`)


# ─── 5.  Tidy cf/WGS “clinical_data” (M4 + MyC) ─────────────────────────────────
clinical_data <- clinical_data %>%
  mutate(
    Sample_type = case_when(
      Sample_type == "P" ~ "Blood_plasma_cfDNA",
      Sample_type == "O" ~ "BM_cells",
      Sample_type == "B" ~ "Blood_Buffy_coat",
      TRUE ~ Sample_type
    )
  )

# (Note: at this point clinical_data already has its own Study column 
#  or you’ll add one manually, e.g., mutate(Study="M4") or mutate(Study="MyC").)


# ─── 6.  Combine SPORE + clinical_data into one “combined_clinical_data” ───────
combined_clinical_data <- bind_rows(spore_data, clinical_data) %>%
  filter(!is.na(Timepoint))   # drop any rows that lack a Timepoint


# ─── 7.  Annotate M4 timepoint_info (diagnosis, post_induction, etc.) ──────────
combined_clinical_data <- combined_clinical_data %>%
  mutate(
    timepoint_info = case_when(
      Study == "M4" & Timepoint == "01" ~ "Diagnosis",
      Study == "M4" & Timepoint == "03" ~ "Post_induction",
      Study == "M4" & Timepoint == "05" ~ "Post_transplant",
      Study == "M4" & Timepoint == "07" ~ "Maintenance",
      Timepoint == "R-"                 ~ "Relapse",
      TRUE                               ~ NA_character_
    )
  )


# ─── 8.  Build SPORE timepoints long → integrate with combined_clinical_data ─────
##  8.1  Ensure SPORE timepoint info dates are Date class
spore_timepoint_info <- spore_timepoint_info %>%
  mutate(`Timepoint of interest` = as.Date(`Timepoint of interest`))

##  8.2  Transform “spore_data_extracted” to long: Actual_Diagnosis → Diagnosis,
##       Relapse_date* → Progression + add Status_at_timepoint = NA
spore_data_long <- spore_data_extracted %>%
  select(Patient, Actual_Diagnosis_Date, Relapse_date, Relapse_date_2, Relapse_date_3,
         Relapse_date_4, Relapse_date_5, Relapse_date_6) %>%
  pivot_longer(
    cols      = c(Actual_Diagnosis_Date, Relapse_date, Relapse_date_2, Relapse_date_3,
                  Relapse_date_4, Relapse_date_5, Relapse_date_6),
    names_to  = "timepoint_info",
    values_to = "Timepoint of interest"
  ) %>%
  filter(!is.na(`Timepoint of interest`)) %>%
  mutate(
    timepoint_info = case_when(
      timepoint_info == "Actual_Diagnosis_Date" ~ "Diagnosis",
      grepl("Relapse_date", timepoint_info)    ~ "Progression",
      TRUE                                      ~ timepoint_info
    ),
    Status_at_timepoint = NA_character_
  )

##  8.3  Bind the two SPORE sources & remove duplicates
combined_spore_timepoint_info <- bind_rows(
  spore_timepoint_info,
  spore_data_long
) %>%
  select(Patient, `Timepoint of interest`, timepoint_info) %>%
  mutate(timepoint_info = gsub("Relapse", "Progression", timepoint_info)) %>%
  distinct()


# ─── 9.  Join combined_clinical_data (only SPORE rows) → timepoint_info by Patient ─
##      and pick the row with smallest date_diff ≤ 60 days; else “Awaiting_integration”
combined_clinical_data <- combined_clinical_data %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

temp <- combined_clinical_data %>%
  filter(Study == "SPORE") %>%
  select(-timepoint_info) %>%
  left_join(combined_spore_timepoint_info, by = "Patient") %>%
  mutate(
    date_diff = abs(as.numeric(
      difftime(Date_of_sample_collection, `Timepoint of interest`, units = "days")
    ))
  ) %>%
  mutate(
    timepoint_info = ifelse(!is.na(date_diff) & date_diff <= 60, timepoint_info, "Awaiting_integration")
  ) %>%
  filter(timepoint_info != "Awaiting_integration") %>%
  group_by(Bam) %>%
  slice_min(date_diff, with_ties = FALSE) %>%
  select(-`Timepoint of interest`, -date_diff) %>%
  distinct() %>%
  ungroup()

##  9.1  All SPORE rows without a match get “Awaiting_integration”
combined_spore_data <- combined_clinical_data %>%
  filter(Study == "SPORE") %>%
  select(-timepoint_info) %>%
  left_join(temp, by = c("Bam","Patient","Sample_type","Timepoint")) %>%
  mutate(timepoint_info = coalesce(timepoint_info, "Awaiting_integration"))

##  9.2  Re‐bind SPORE + non‐SPORE (to update combined_clinical_data_updated)
combined_clinical_data_updated <- combined_clinical_data %>%
  filter(Study != "SPORE") %>%
  bind_rows(combined_spore_data) %>%
  mutate(timepoint_info = if_else(is.na(timepoint_info), "Awaiting_integration", timepoint_info))


# ─── 10.  Fix special M4 cases (CA-08, IMG-060, IMG-098, etc.) ───────────────────
##   (e.g. “CA-08” mislabeled “09” → “08”, IMG-060 T0 → Diagnosis, etc.)
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    timepoint_info = case_when(
      Study == "M4" & Timepoint == "08" ~ "1.5yr maintenance",
      Study == "M4" & Timepoint == "09" & Patient == "CA-08" ~ "1.5yr maintenance",
      Study == "M4" & Timepoint == "09"              ~ "2yr maintenance",
      Study == "M4" & Timepoint == "10"              ~ "2.5yr maintenance",
      Patient == "IMG-060" & Timepoint == "T0"       ~ "Diagnosis",
      Patient == "IMG-098" & (Timepoint == "T0" | Timepoint == "T1") ~ "Diagnosis",
      grepl("IMG", Patient) & Timepoint == "T0"      ~ "Diagnosis",
      TRUE                                           ~ timepoint_info
    )
  )

## Fix IMG cases
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    Date_of_sample_collection = ifelse(Patient == "IMG-181" & Timepoint == "T0", as.Date("2021-04-27"), Date_of_sample_collection),
    timepoint_info = case_when(
      Patient == "IMG-181" & Timepoint == "T0" ~ "Diagnosis",
      Patient == "IMG-181" & Timepoint %in% c("T6", "T12", "T15") ~ "Maintenance",
      TRUE ~ timepoint_info
    ),
    Timepoint = ifelse(Patient == "IMG-181" & Timepoint == "T6", "07", Timepoint)
  )

# ─── 11.  Write out the first “combined_clinical_data” intermiate file ─────────
write.csv(
  combined_clinical_data_updated,
  "combined_clinical_data_updated_tmp.csv",
  row.names = FALSE
)


# ─── 12.  QC: get counts & write summary tables ────────────────────────────────
##  12.1  Number of unique patients
num_patients <- combined_clinical_data_updated %>%
  summarise(num_patients = n_distinct(Patient))

##  12.2  Samples by Sample_type
samples_by_type      <- combined_clinical_data_updated %>%
  group_by(Sample_type) %>%
  summarise(num_samples = n())

##  12.3  Samples by timepoint_info
samples_by_timepoint <- combined_clinical_data_updated %>%
  group_by(timepoint_info) %>%
  summarise(num_samples = n())

##  12.4  Samples by Study
samples_by_study     <- combined_clinical_data_updated %>%
  group_by(Study) %>%
  summarise(num_samples = n())

##  12.5  Print to screen & write CSVs
print(num_patients)
print(samples_by_type)
print(samples_by_timepoint)
print(samples_by_study)

# Export the number of unique patients
write.csv(num_patients, "num_patients.csv", row.names = FALSE)

# Export the number of samples by sample type
write.csv(samples_by_type, "samples_by_type.csv", row.names = FALSE)

# Export the number of samples by timepoint
write.csv(samples_by_timepoint, "samples_by_timepoint.csv", row.names = FALSE)

# Export the number of samples by study
write.csv(samples_by_study, "samples_by_study.csv", row.names = FALSE)


# ─── 13.  Make simple barplots (Sample_type, timepoint_info, Study) ─────────────
output_dir <- "plots_Oct"
dir.create(output_dir, showWarnings = FALSE)

##  13.1  Samples by Sample_type
plot_samples_by_type <- ggplot(samples_by_type,
                               aes(x = Sample_type, y = num_samples, fill = Sample_type)) +
  geom_col(stat = "identity") +
  labs(title = "Number of Samples by Sample Type",
       x = "Sample Type", y = "Number of Samples") +
  theme_minimal()
ggsave(file.path(output_dir, "samples_by_type.png"), width = 6, height = 4)

##  13.2  Samples by Timepoint (factor levels in order)
samples_by_timepoint$timepoint_info <- factor(
  samples_by_timepoint$timepoint_info,
  levels = c("Diagnosis","Baseline","Post_induction","Post_transplant",
             "Maintenance","1.5yr maintenance","2yr maintenance",
             "Progression","Relapse","Awaiting_integration","Unknown")
)
plot_samples_by_timepoint <- ggplot(samples_by_timepoint,
                                    aes(x = timepoint_info, y = num_samples, fill = timepoint_info)) +
  geom_col(stat = "identity") +
  geom_text(aes(label = paste0("n=", num_samples)), vjust = -0.5, size = 3.5) +
  labs(title = "Number of Samples by Timepoint",
       x = "Timepoint", y = "Number of Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "samples_by_timepoint.png"), width = 6, height = 4)

##  13.3  Stacked bar: Samples by Study & Sample_type
samples_by_study_and_type <- combined_clinical_data %>%
  group_by(Study, Sample_type) %>%
  summarise(num_samples = n(), .groups = "drop")

plot_samples_by_study_and_type <- ggplot(samples_by_study_and_type,
                                         aes(x = Study, y = num_samples, fill = Sample_type)) +
  geom_col(stat = "identity", position = "stack") +
  labs(title = "Samples by Study and Sample Type",
       x = "Study", y = "Number of Samples") +
  theme_minimal()
ggsave(file.path(output_dir, "samples_by_study_and_type.png"), width = 6, height = 4)


##  13.4  Number of Patients by Study
patients_by_study <- combined_clinical_data %>%
  group_by(Study) %>%
  summarise(num_patients = n_distinct(Patient), .groups = "drop")

plot_patients_by_study <- ggplot(patients_by_study,
                                 aes(x = Study, y = num_patients, fill = Study)) +
  geom_col(stat = "identity") +
  geom_text(aes(label = paste0("n=", num_patients)), vjust = -0.5, size = 5) +
  labs(title = "Number of Patients by Study",
       x = "Study", y = "Number of Patients") +
  theme_minimal()
ggsave(file.path(output_dir, "patients_by_study.png"), width = 6, height = 4)


# ─── 14.  Tidy M4 labs & extract timepoints → combine back into main table ──────
##  14.1  Read M4 lab workbook
M4_Labs <- read_excel("M4_CMRG_Data/M4_COHORT_LABS.xlsx")

##  14.2  Select columns
M4_Labs <- M4_Labs %>%
  select(M4_id, study_patient_id, LAB_TYPE, LAB_DATE, PURPOSE, LAB_NAME, LAB_VALUE)

##  14.3  Define extract_timepoint() (purpose → “01”, “03”, … “R”)
# Define the function to extract the timepoint
extract_timepoint <- function(purpose) {
  if (is.na(purpose) || purpose == "") return(NA)
  
  # Convert to lower case for case-insensitive matching
  purpose_lower <- tolower(purpose)
  
  # Remove leading and trailing whitespace
  purpose_lower <- trimws(purpose_lower)
  
  # Initialize timepoint as NA
  timepoint <- NA
  
  # Check for 'Diagnosis' related terms
  if (grepl("diagnosis|at diagnosis|baseline|at first consult|baseline line|dara baseline", purpose_lower)) {
    timepoint <- "01"
  } 
  # Check for 'Relapse'
  else if (grepl("relapse", purpose_lower)) {
    timepoint <- "R"
  } 
  else {
    # Try to extract visit number patterns
    match <- regexpr("v(isit)?[ .-]*#?[ ]*(\\d+|r)", purpose_lower, perl=TRUE)
    if (match[1] != -1) {
      visit_str <- regmatches(purpose_lower, match)
      visit_num <- sub("v(isit)?[ .-]*#?[ ]*", "", visit_str, perl=TRUE)
      visit_num <- toupper(visit_num)
      if (visit_num == "R") {
        timepoint <- "R"
      } else if (grepl("^\\d+$", visit_num)) {
        timepoint <- sprintf("%02d", as.integer(visit_num))
      }
    }
  }
  return(timepoint)
}

##  14.4  Apply extract_timepoint() → TIMEPOINT + map to descriptions
description_mapping <- c(
  "01" = "Diagnosis", "03" = "Post-induction", "05" = "Post-transplant",
  "07" = "1yr maintenance", "08" = "1.5yr maintenance", "09" = "2yr maintenance",
  "10" = "2.5yr maintenance", "11" = "3yr maintenance", "12" = "3.5yr maintenance",
  "13" = "4yr maintenance", "14" = "4.5yr maintenance", "15" = "5yr maintenance",
  "R"  = "Relapse"
)

M4_Labs <- M4_Labs %>%
  mutate(
    TIMEPOINT       = sapply(PURPOSE, extract_timepoint),
    TIMEPOINT_DESC  = description_mapping[TIMEPOINT]
  )

write.csv(M4_labs_cleaned, file = "M4_labs_cleaned.csv")

##  14.5  Make M4_dates_cmrg with unique M4_id, LAB_DATE, TIMEPOINT, TIMEPOINT_DESC
M4_dates_cmrg <- M4_Labs %>%
  select(M4_id, LAB_DATE, PURPOSE, TIMEPOINT, TIMEPOINT_DESC) %>%
  distinct()


# ─── 15.  (Legacy) “M4 Status & best response” Excel → long ─────────────────────
data <- read_excel("Clinical data/M4/M4 Status and best response update Dec 2 24 just dates.xlsx")
names(data) <- make.names(names(data))
data <- data %>% select(-Best.Response.CR.or.VGPR, -Date)

colnames(data) <- c(
  "SITE", "PATIENT", "Current.Status", "CONSENT",
  "Timepoint_01_Prior_TX", "Timepoint_03_Post_Chemo",
  "Timepoint_05_Post_ASCT_100_or_1-3_DYS", "Timepoint_07_Post_Maintenance_Therapy_12_MONTHS",
  "Timepoint_08_6_Months_Post_07", "Timepoint_09_Every_6_Months",
  "Timepoint_10_Every_6_Months", "Timepoint_11_Every_6_Months",
  "Timepoint_12_Every_6_Months", "Timepoint_13_Every_6_Months",
  "Timepoint_R_Relapse"
)

# Transform data to long format and clean it to standardize timepoint names
data_long <- data %>%
  pivot_longer(
    cols      = -c(SITE, PATIENT, Current.Status),
    names_to  = "Timepoint",
    values_to = "Date"
  ) %>%
  mutate(
    Timepoint = str_trim(Timepoint),
    Timepoint = gsub("\\s+", " ", Timepoint),
    Date      = sapply(Date, function(x) {
      parsed <- parse_date_time(as.character(x), orders = c("dmy","mdy","ymd"), quiet = TRUE)
      as_date(parsed)
    }),
    Timepoint_Code = case_when(
      Timepoint == "Timepoint_01_Prior_TX"                  ~ "01",
      Timepoint == "Timepoint_03_Post_Chemo"                ~ "03",
      Timepoint == "Timepoint_05_Post_ASCT_100_or_1-3_DYS"   ~ "05",
      Timepoint == "Timepoint_07_Post_Maintenance_Therapy_12_MONTHS" ~ "07",
      Timepoint == "Timepoint_08_6_Months_Post_07"           ~ "08",
      Timepoint == "Timepoint_09_Every_6_Months"             ~ "09",
      Timepoint == "Timepoint_10_Every_6_Months"             ~ "10",
      Timepoint == "Timepoint_11_Every_6_Months"             ~ "11",
      Timepoint == "Timepoint_12_Every_6_Months"             ~ "12",
      Timepoint == "Timepoint_13_Every_6_Months"             ~ "13",
      Timepoint == "Timepoint_R_Relapse"                     ~ "R",
      Timepoint == "CONSENT"                                 ~ "CONSENT",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(PATIENT = str_pad(as.character(PATIENT), width = 2, pad = "0")) %>%
  filter(!is.na(Date), Timepoint_Code != "CONSENT") %>%
  rename(Patient_Code = PATIENT) %>%
  mutate(Patient = paste(SITE, Patient_Code, sep = "-")) %>%
  distinct()

M4_dates <- data_long
rm(data_long)

##  15.1  Convert “R” → “R-” in M4_dates
M4_dates <- M4_dates %>%
  mutate(Timepoint_Code = ifelse(Timepoint_Code == "R", "R-", Timepoint_Code))


# ─── 16.  Merge M4_dates into combined_clinical_data_updated ────────────────────
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(
    M4_dates %>% select(Patient, Timepoint_Code, Date),
    by = c("Patient", "Timepoint" = "Timepoint_Code")
  ) %>%
  mutate(
    Date_of_sample_collection = coalesce(Date, Date_of_sample_collection)
  ) %>%
  select(-Date)

##  16.1  Write intermediate file + M4_dates to disk)
write.csv(M4_dates,
          "Exported_data_tables_clinical/M4_dates_Nov2024.csv",
          row.names = FALSE)


# ─── 17.  Use M4_Labs_earliest to fill any remaining missing M4 dates ───────────
M4_Labs_filtered <- M4_Labs %>%
  select(M4_id, PURPOSE, LAB_DATE, TIMEPOINT, TIMEPOINT_DESC) %>%
  filter(!grepl("Baseline|baseline", PURPOSE)) %>%
  distinct()

M4_Labs_filtered <- M4_Labs_filtered %>%
  mutate(Timepoint_Code = ifelse(TIMEPOINT == "R", "R-", TIMEPOINT)) %>%
  filter(!is.na(TIMEPOINT))

M4_Labs_earliest <- M4_Labs_filtered %>%
  group_by(M4_id, TIMEPOINT) %>%
  filter(LAB_DATE == min(LAB_DATE, na.rm = TRUE)) %>%
  ungroup()

##  17.1  Fix CA-08 mislabeled “08”
M4_Labs_earliest <- M4_Labs_earliest %>%
  mutate(Timepoint_Code = if_else(M4_id == "CA-08" & Timepoint_Code == "08", "09", Timepoint_Code))

##  17.2  Read additional relapse dates from Excel → Relapse_dates_M4_clean
Relapse_dates_M4 <- read_excel("M4_CMRG_Data/Relapse_dates.xlsx")

Relapse_dates_M4_clean <- Relapse_dates_M4 %>%
  select(SITE, PATIENT, `Current Status`, `CMRG Relapse Date`, `M4 Relapse Visit`) %>%
  rename(
    Site            = SITE,
    Patient         = PATIENT,
    Status          = `Current Status`,
    CMRG_Relapse_Date = `CMRG Relapse Date`,
    M4_Relapse_Visit  = `M4 Relapse Visit`
  ) %>%
  mutate(
    CMRG_Relapse_Date = as.numeric(CMRG_Relapse_Date),
    M4_Relapse_Visit  = as.numeric(M4_Relapse_Visit),
    CMRG_Relapse_Date = as.Date(CMRG_Relapse_Date, origin = "1899-12-30"),
    M4_Relapse_Visit  = as.Date(M4_Relapse_Visit,  origin = "1899-12-30")
  )

##  17.3  Correct error in sample IDs (diagnosis_bam_values, unknown_timepoint_bam)
diagnosis_bam_values   <- c(
  "TFRIM4_0166_Bm_P_WG_CA-06.filter.deduped.recalibrated.bam",
  "TFRIM4_0168_Bm_P_WG_CA-15.filter.deduped.recalibrated.bam",
  "TFRIM4_0169_Bm_P_WG_CA-16.filter.deduped.recalibrated.bam",
  "TFRIM4_0175_Bm_P_WG_EK-09.filter.deduped.recalibrated.bam"
)
unknown_timepoint_bam  <- "TFRIM4_0178_Bm_P_WG_FZ-07.filter.deduped.recalibrated.bam"

# Update Timepoint, timepoint_info, and Date_of_sample_collection accordingly
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    Timepoint = case_when(
      Bam %in% diagnosis_bam_values ~ "01",
      Bam == unknown_timepoint_bam ~ NA_character_,
      TRUE ~ Timepoint
    ),
    timepoint_info = ifelse(Bam %in% diagnosis_bam_values, "Diagnosis", timepoint_info),
    Date_of_sample_collection = ifelse(Bam %in% c(diagnosis_bam_values, unknown_timepoint_bam), NA, Date_of_sample_collection)
  )

## Edit the IDs to match
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(Sample_ID = ifelse(Bam %in% diagnosis_bam_values, gsub("-R-O", "-01-O", Sample_ID), Sample_ID))


## Now add to previous dataframe 
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(M4_dates %>% select(Patient, Timepoint_Code, Date), 
            by = c("Patient", "Timepoint" = "Timepoint_Code")) %>%
  mutate(Date_of_sample_collection = coalesce(Date, Date_of_sample_collection)) %>%
  select(-Date)  # Remove the Date column after updating

## If still NA, go off of labs
## Merge 
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(M4_Labs_earliest %>% dplyr::select(M4_id, Timepoint_Code, LAB_DATE),
            by = c("Patient" = "M4_id", "Timepoint" = "Timepoint_Code")) %>%
  # Update Date_of_sample_collection only if it is NA, using LAB_DATE where available
  mutate(Date_of_sample_collection = if_else(is.na(Date_of_sample_collection), LAB_DATE, Date_of_sample_collection)) %>%
  # Drop the temporary LAB_DATE column added from the join
  select(-LAB_DATE) %>% 
  unique()

##  17.4  Merge M4_dates back (again) + then M4_Labs_earliest to fill any remaining NAs
## Get all unique dates for PFS
# Step 1: Separate SITE and Patient_Code in M4_Labs_earliest
tmp <- M4_Labs_earliest %>%
  separate(M4_id, into = c("SITE", "Patient_Code"), sep = "-")

# Step 2: Identify new combinations in M4_Labs_earliest not present in M4_dates
new_combinations <- tmp %>%
  rename(Date = LAB_DATE) %>%
  select(SITE, Patient_Code, Timepoint_Code, Date) %>%
  anti_join(M4_dates %>% select(SITE, Patient_Code, Timepoint_Code), by = c("SITE", "Patient_Code", "Timepoint_Code"))

# Step 3: Combine new combinations with M4_dates, prioritizing dates in M4_dates
combined_dates <- M4_dates %>%
  select(SITE, Patient_Code, Timepoint_Code, Date) %>%
  bind_rows(new_combinations)

# Step 4: Calculate the date range for each patient
patient_date_range <- combined_dates %>%
  group_by(SITE, Patient_Code) %>%
  summarise(
    Earliest_Date = min(Date, na.rm = TRUE),
    Latest_Date = max(Date, na.rm = TRUE),
    Days_Between = as.numeric(Latest_Date - Earliest_Date)
  ) %>%
  ungroup()

# Create the Patient column by concatenating SITE and Patient_Code
patient_date_range <- patient_date_range %>%
  mutate(Patient = paste(SITE, Patient_Code, sep = "-")) %>%
  select(Patient, everything())  # Optional: Move Patient column to the front

# Export the dataframe to a CSV file
write.csv(patient_date_range, "Oct2024_exported_data/patient_date_range_M4_PFS.csv", row.names = FALSE)


##  17.5 Update the VA-05 with updated info 
# Define the specific Bam value to update
target_bam <- "TFRIM4_0181_Cf_P_PG_VA-07-05-P-DNA.filter.deduped.recalibrated.bam"

# Update Date_of_sample_collection, Timepoint, timepoint_info, and modify Sample_ID
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection, origin = "1970-01-01"))

combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    Date_of_sample_collection = ifelse(Bam == target_bam, as.Date("2021-01-01", format = "%Y-%m-%d"), Date_of_sample_collection),
    Timepoint = ifelse(Bam == target_bam, "R-", Timepoint),
    timepoint_info = ifelse(Bam == target_bam, "Relapse", timepoint_info),
    Sample_ID = ifelse(Bam == target_bam, gsub("05-P", "R-P", Sample_ID), Sample_ID)
  )

combined_clinical_data_updated$Date_of_sample_collection <- as.Date(combined_clinical_data_updated$Date_of_sample_collection)


## 17.6 Correct IMG cases
# Set Plot_timepoint_numeric to 6 for the specific BAM entry
combined_clinical_data_updated$Timepoint[combined_clinical_data_updated$Bam == "MYC_0181_Pl_T_PG_MyP-181-T2-P.filter.deduped.recalibrated.bam"] <- "T6"

## 17.7 Fix FZ-07 BM sample in combined_clinical_data_updated
source_row <- combined_clinical_data_updated %>%
  filter(Bam == "TFRIM4_0178_Cf_P_PG_FZ-07-01-P-DNA.filter.deduped.recalibrated.bam")

new_row <- source_row
new_row$Bam         <- "TFRIM4_0178_Bm_P_WG_FZ-07.filter.deduped.recalibrated.bam"
new_row$Sample_ID   <- "FZ-07-01-O"
new_row$Sample_type <- "BM_cells"

combined_clinical_data_updated <- combined_clinical_data_updated %>%
  filter(Bam != "TFRIM4_0178_Bm_P_WG_FZ-07.filter.deduped.recalibrated.bam") %>%
  bind_rows(new_row)




# ─── 18.  Get progression info  ─────────────────

# 18.1 Calculate the progression info 
# Add Sample_ID column with Site and zero-padded Patient
Relapse_dates_M4_clean <- Relapse_dates_M4_clean %>%
  mutate(Sample_ID = paste0(Site, "-", sprintf("%02d", Patient)))

### Add in missing dates of progressive disease from CMRG
CMRG_progressive_disease <- read_excel("M4_CMRG_Data/M4_COHORT_CHEMOTHERAPY.xlsx")

# Perform a left join to bring PROGRESSION_DATE into Relapse_dates_M4_clean based on Sample_ID and M4_id
Relapse_dates_M4_clean <- Relapse_dates_M4_clean %>%
  left_join(
    CMRG_progressive_disease %>% filter(!is.na(PROGRESSION_DATE)) %>% unique() %>% dplyr::select(M4_id, PROGRESSION_DATE),
    by = c("Sample_ID" = "M4_id")
  )

### Add another column called integrated progression date which puts the value of M4_Relapse Visit if Progression is NA
# Create the Integrated_Progression_Date column
Relapse_dates_M4_clean <- Relapse_dates_M4_clean %>%
  mutate(
    Integrated_Progression_Date = as.Date(ifelse(
      is.na(PROGRESSION_DATE),
      as.Date(M4_Relapse_Visit),
      as.Date(PROGRESSION_DATE))
    )
  )

## Add everything together 
# Step 1: Filter combined_spore_timepoint_info for rows with timepoint_info == "Progression"
spore_progression <- combined_spore_timepoint_info %>%
  filter(timepoint_info == "Progression") %>%
  rename(Progression_date = `Timepoint of interest`)

## Add additional dates 
# Step 1: Modify spore_OS_info
tmp <- spore_OS_info %>%
  mutate(
    timepoint_info = ifelse(`Current status` == "Deceased", "Progression", NA),
    Progression_date = `Status last follow up`
  ) %>%
  dplyr::select(Patient, Progression_date, timepoint_info) %>%
  filter(!is.na(timepoint_info))

# Step 2: Bind rows with spore_progression
spore_progression <- bind_rows(spore_progression, tmp)


# Step 2: Select Patient and PROGRESSION_DATE columns from Relapse_dates_M4_clean
m4_progression <- Relapse_dates_M4_clean %>%
  select(Sample_ID, Integrated_Progression_Date) %>%
  rename(Patient = Sample_ID, Progression_date = Integrated_Progression_Date)

# Step 3: Combine the two dataframes into Relapse_dates_full
Relapse_dates_full <- bind_rows(spore_progression, m4_progression) %>%
  filter(!is.na(Progression_date)) %>% select(-timepoint_info)

### Add the IMMAGINE dates here 
IMMAGINE_progression <- read_excel("Clinical data/IMMAGINE/Extracted_clinical_MRD_data.xlsx", sheet = 1)

# Fix dates
IMMAGINE_progression <- IMMAGINE_progression %>%
  # Split the 'Relapse Dates' column into separate rows
  separate_rows(`Relapse Dates`, sep = ", ") %>%
  # Rename columns to the desired names
  rename(Diagnosis_date = `Diagnosis Date`, Relapse_date = `Relapse Dates`) %>%
  # Convert 'Diagnosis_date' and 'Relapse_date' columns to Date format
  mutate(
    Diagnosis_date = as.Date(Diagnosis_date),
    Relapse_date = as.Date(Relapse_date, format = "%Y-%m-%d")
  )

## Join to relapse_dates 
temp <- IMMAGINE_progression %>% 
  select(ID, Relapse_date) %>%
  rename(Patient = ID, Progression_date = Relapse_date) %>% 
  filter(!is.na(Progression_date)) %>% 
  unique()

Relapse_dates_full <- bind_rows(Relapse_dates_full, temp)
Relapse_dates_full <- bind_rows(Relapse_dates_full, tmp %>% select(Patient, Progression_date))
## Now calculate the time to closest progression for each sample in the combined_clinical_info

# Ensure Date_of_sample_collection and Progression_date are in Date format
combined_clinical_data_updated$Date_of_sample_collection <- as.Date(combined_clinical_data_updated$Date_of_sample_collection)
Relapse_dates_full$Progression_date <- as.Date(Relapse_dates_full$Progression_date)

# Step 1: Add 'Relapsed' column
## Need to edit for those that relapsed before sample collected
test <- combined_clinical_data_updated %>%
  mutate(
    Relapsed = ifelse(Patient %in% Relapse_dates_full$Patient, "Y", "N")
  )

## Add the other dates got from David 
# Create a data frame with the new entries
new_dates <- data.frame(
  Patient = c("IMG-060", "IMG-098"),
  Progression_date = as.Date(c("2022-08-29", "2024-06-03")) # Check these
)

# Add the new dates to the existing data frame
Relapse_dates_full <- bind_rows(Relapse_dates_full, new_dates) %>% unique()

## Export relapse dates 
write.csv(Relapse_dates_full, file = "Relapse dates cfWGS updated.csv", row.names = F)
write.csv(Relapse_dates_M4_clean, "Relapse_dates_M4_clean.csv", row.names = FALSE)


# 18.2  Find the maximum number of days available per patient 
# Join on Patient and calculate the days difference
# Only relapse dates AFTER the sample was collected, unless within 35 days
## Then count that as a relapse timepoint
Time_to_relapse2 <- test %>%
  left_join(
    Relapse_dates_full,
    by = c("Patient")
  ) %>%
  group_by(Patient, Date_of_sample_collection) %>%
  mutate(
    # Calculate non-absolute days to relapse only for valid dates
    days_to_relapse_non_absolute = case_when(
      Progression_date > Date_of_sample_collection ~ as.numeric(Date_of_sample_collection - Progression_date), 
      Progression_date <= Date_of_sample_collection & Progression_date >= Date_of_sample_collection - 35 ~ 0,
      TRUE ~ NA_real_  # Set to NA for dates that do not meet conditions
    )
  ) %>%
  filter(!is.na(days_to_relapse_non_absolute)) %>%
  summarise(
    Relapsed = first(Relapsed),
    
    # Calculate the absolute days to relapse for valid dates
    Num_days_to_closest_relapse_absolute = min(abs(days_to_relapse_non_absolute), na.rm = TRUE),
    
    # Get the closest relapse date that came strictly after the sample collection
    Num_days_to_closest_relapse_non_absolute = days_to_relapse_non_absolute[which.min(abs(days_to_relapse_non_absolute))]
  ) %>%
  ungroup()

## Now add this back 
tmp <- combined_clinical_data_updated
combined_clinical_data_updated <- combined_clinical_data_updated %>% 
  dplyr::select(-Num_days_to_closest_relapse, -Relapsed, -Num_days_to_closest_relapse_absolute, -Num_days_to_closest_relapse_non_absolute) %>%
  left_join(Time_to_relapse2)

## Add other timepoint_info
tmp <- read.csv("Clinical data/IMMAGINE/NAs in timepoint info - IMMAGINE additional timepoint info.csv")
tmp <- tmp %>% 
  rename(timepoint_info_updated = timepoint_info) %>% 
  select(-Date_of_sample_collection)

# Merge combined_clinical_data_updated with specimen_collection, replacing `timepoint_info` with `timepoint_info_updated` wherever there is a match
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(tmp, by = c("Patient", "Sample_ID")) %>%
  dplyr::mutate(
    # Replace `timepoint_info` with `timepoint_info_updated` wherever there is a match
    timepoint_info = if_else(!is.na(timepoint_info_updated), timepoint_info_updated, timepoint_info)
  ) %>%
  dplyr::select(-timepoint_info_updated)  # Drop `timepoint_info_updated` after merging

# Replace "MRD" with "Maintenance" in the `timepoint_info` column
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    timepoint_info = if_else(timepoint_info == "MRD", "Maintenance", timepoint_info)
  )


### Take the date of the relapse from the baseline marrow 
# Filter baseline bone marrow samples with "Diagnosis" or "baseline" in timepoint_info
## Edit the relapse dates of additional SPORE patients 
baseline_BM_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "BM_cells", timepoint_info %in% c("Diagnosis", "Baseline"))

## Export 
write.csv(baseline_with_relapse, file = "Oct2024_exported_data/Days to relapse baseline BM cell samples updated.txt", row.names = F)

combined_clinical_data_updated$Num_days_to_closest_relapse <- combined_clinical_data_updated$Num_days_to_closest_relapse_absolute

# Update Relapsed to "No" only if it is NA
combined_clinical_data_updated$Relapsed[is.na(combined_clinical_data_updated$Relapsed)] <- "N"

# View the result to confirm changes
table(combined_clinical_data_updated$Relapsed)


### Calculate longest follwowup time have on patients 
### Need additional data on immagine for this to see how long we have
# Step 1: Identify baseline marrow sample date for each patient
baseline_dates <- combined_clinical_data_updated %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  group_by(Patient) %>%
  summarise(Baseline_Date = min(Date_of_sample_collection, na.rm = TRUE)) %>%
  ungroup()

# Step 2: Join baseline dates with main data
combined_with_baseline <- combined_clinical_data_updated %>% filter(Sample_type == "Blood_plasma_cfDNA") %>%
  left_join(baseline_dates, by = "Patient") %>%
  # Step 3: Calculate follow-up days only for dates after baseline
  mutate(
    Followup_Days = ifelse(Date_of_sample_collection > Baseline_Date,
                           as.numeric(Date_of_sample_collection - Baseline_Date),
                           NA)
  )

# Step 4: Summarize maximum follow-up days for each patient
patient_followup <- combined_with_baseline %>%
  group_by(Patient) %>%
  summarise(Max_Followup_Days = max(Followup_Days, na.rm = TRUE)) %>%
  ungroup()

patient_followup <- patient_followup %>% filter(Max_Followup_Days > 0)

## Now join with relapse info to get the PFS
Time_to_relapse2 <- combined_clinical_data_updated %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  left_join(
    Relapse_dates_full,
    by = c("Patient")
  ) %>%
  group_by(Patient, Date_of_sample_collection) %>%
  mutate(
    # Calculate non-absolute days to relapse only for valid dates
    days_to_relapse_non_absolute = case_when(
      Progression_date > Date_of_sample_collection ~ as.numeric(Date_of_sample_collection - Progression_date), 
      Progression_date <= Date_of_sample_collection & Progression_date >= Date_of_sample_collection - 35 ~ 0,
      TRUE ~ NA_real_  # Set to NA for dates that do not meet conditions
    )
  ) %>%
  filter(!is.na(days_to_relapse_non_absolute)) %>%
  summarise(
    Relapsed = first(Relapsed),
    
    # Calculate the absolute days to relapse for valid dates
    Num_days_to_closest_relapse_absolute = min(abs(days_to_relapse_non_absolute), na.rm = TRUE),
    
    # Get the closest relapse date that came strictly after the sample collection
    Num_days_to_closest_relapse_non_absolute = days_to_relapse_non_absolute[which.min(abs(days_to_relapse_non_absolute))]
  ) %>%
  ungroup()

# Step to keep the highest Num_days_to_closest_relapse_absolute for each patient
Time_to_relapse2 <- Time_to_relapse2 %>% 
  group_by(Patient) %>%
  slice_max(Num_days_to_closest_relapse_absolute, with_ties = FALSE) %>%
  ungroup()

## Now add this back 
patient_followup <- left_join(patient_followup, Time_to_relapse2)

### Have to think about multiple relapse cases, what to do for that for the kaplan myer curves 
# Keep only the max num_days for each patient
max_num_days_per_patient <- patient_followup %>%
  group_by(Patient) %>%
  summarise(
    max_num_days = max(Max_Followup_Days, na.rm = TRUE),
    max_num_weeks = max_num_days / 7
  ) %>%
  ungroup()


## Export
write.table(max_num_days_per_patient, file = "Oct2024_exported_data/max_num_days_per_patient_Jan28_updated.txt", sep = "\t", quote = F, row.names = F)



### 18.3 Clean processing log to get additional dates and export 
file_path <- "TFRIM4_Processing Log_Nov2024.xlsx"  # Replace with your file path
processing_log <- read_excel(file_path)

# Step 3: Select only the first three columns
processing_log_cleaned <- processing_log %>%
  select(
    `Plasma Pugh Lab Identifier`,
    `Red is on 35-Patient list`,
    `Blood Drawn Date [YYYY/MM/DD]`
  )

# Ensure naming is consistent
processing_log_cleaned <- processing_log_cleaned %>%
  mutate(`Plasma Pugh Lab Identifier` = str_replace_all(`Plasma Pugh Lab Identifier`, "_", "-"))


# Step 4: Rename and process Sample_Code
processing_log_cleaned <- processing_log_cleaned %>%
  rename(Sample_Code = `Plasma Pugh Lab Identifier`) %>%  # Rename the column
  mutate(
    Sample_Code = str_remove(Sample_Code, "M4-"),  # Remove "M4-" prefix
    Sample_Code = str_remove_all(Sample_Code, "(_[0-9]|-[0-9])$"),  # Remove endings like -1, -2, _1, _2
    Patient = str_sub(Sample_Code, 1, 5),  # Extract the first 5 characters as Patient
    Timepoint = str_extract(Sample_Code, "[^-]+$")  # Extract everything after the last "-"
  )

## Make date consistent
# Define a function to process mixed date formats
process_dates <- function(date_column) {
  date_column <- gsub("\\?", NA, date_column) # Replace question marks with NA
  date_column <- trimws(date_column)  # Remove leading/trailing whitespace
  
  processed_dates <- sapply(date_column, function(date) {
    if (is.na(date) || date == "") {
      return(NA)
    } else if (grepl("^\\d{5}$", date)) {
      # Convert Excel serial numbers (Julian dates)
      as.character(as.Date(as.numeric(date), origin = "1899-12-30"))
    } else if (grepl("^\\d{4}-\\d{2}-\\d{2}$", date)) {
      # ISO 8601 format: YYYY-MM-DD
      as.character(as.Date(date))
    } else if (grepl("^\\d{2}/\\d{1,2}/\\d{4}$", date)) {
      # Format: DD/MM/YYYY or D/M/YYYY
      as.character(dmy(date))
    } else if (grepl("^\\d{1,2}-\\d{1,2}-\\d{4}$", date)) {
      # Format: D-M-YYYY or DD-MM-YYYY
      as.character(dmy(date))
    } else if (grepl("^\\d{4}/\\d{1,2}/\\d{1,2}$", date)) {
      # Format: YYYY/MM/DD or YYYY/M/D
      as.character(ymd(date))
    } else {
      # Invalid or unknown format
      return(NA)
    }
  })
  
  return(processed_dates)
}

# Apply the function to clean the dates
processing_log_cleaned <- processing_log_cleaned %>%
  mutate(
    Cleaned_Date = process_dates(`Blood Drawn Date [YYYY/MM/DD]`),
    Cleaned_Date = as.Date(Cleaned_Date)  # Convert to Date type
  )

## Fix +1 Notation 
processing_log_cleaned <- processing_log_cleaned %>%
  mutate(Timepoint = str_replace_all(Timepoint, "\\+.*", ""))

## Extract relevant columns 
# Process the data
processing_log_cleaned <- processing_log_cleaned %>%
  select(Patient, Timepoint, Cleaned_Date) %>%  # Keep only the specified columns
  distinct() %>%  # Remove duplicate rows
  mutate(Sample_type = "Blood_plasma_cfDNA") %>% # Add the new column with the specified value 
  filter(!is.na(Cleaned_Date))

### Export
# Ensure the folder exists
dir.create("Exported_data_tables_clinical", showWarnings = FALSE)

# Export the processed data as a CSV file
write.csv(processing_log_cleaned, file = "Exported_data_tables_clinical/M4_processing_log_dates.csv", row.names = FALSE)


### 18.4 Get fill dates table for PFS

## Make M4_dates_full to use for the PFS table 
### Join all of the ones available 
# 1) Start with the values in M4_dates

# Replace 'R-' with 'R' in the Timepoint column
M4_dates <- M4_dates %>%
  mutate(Timepoint_Code = gsub("^R-", "R", Timepoint_Code))

df1 <- M4_dates %>%
  select(
    Patient = Patient, 
    Timepoint = Timepoint_Code, 
    Date = Date
  ) %>%
  arrange(Patient, Timepoint, Date) %>%
  group_by(Patient, Timepoint) %>%
  slice(1) %>%  # Take the first date if multiple exist
  ungroup()

# Ensure Date is of Date class
df1$Date <- as.Date(df1$Date)

# 2) Add additional timepoints or patients from M4_dates_cmrg not included in M4_dates
df2 <- M4_dates_cmrg %>%
  select(
    Patient = M4_id, 
    Timepoint = TIMEPOINT, 
    Date = LAB_DATE
  ) %>%
  mutate(
    Date = as.Date(Date), 
    Timepoint = as.character(Timepoint),
    Timepoint = unname(Timepoint)  # Remove names attribute from Timepoint
  ) %>%
  arrange(Patient, Timepoint, Date) %>%
  group_by(Patient, Timepoint) %>%
  slice(1) %>%  # Take the first date if multiple exist
  ungroup()

# Find entries in df2 not already in df1
df2_new <- anti_join(df2, df1, by = c("Patient", "Timepoint"))

# Combine df1 and df2_new
df_combined <- bind_rows(df1, df2_new)

# 3) Add remaining entries from processing_log_cleaned not included in previous data
df3 <- processing_log_cleaned %>%
  select(
    Patient = Patient, 
    Timepoint = Timepoint, 
    Date = Cleaned_Date
  ) %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Patient, Timepoint, Date) %>%
  group_by(Patient, Timepoint) %>%
  slice(1) %>%  # Take the first date if multiple exist
  ungroup()

# Find entries in df3 not already in df_combined
df3_new <- anti_join(df3, df_combined, by = c("Patient", "Timepoint"))

# Combine all data into the super table
df_combined <- bind_rows(df_combined, df3_new)

df_combined <- df_combined %>% filter(!is.na(Timepoint))

## Now make a table with the longest dates have per pt and then the time to relapse for PGS

## Use this for the PFS curves
# Calculate the number of days from the baseline (Timepoint == "01") for each sample
df_combined <- df_combined %>%
  group_by(Patient) %>%
  mutate(Days_from_baseline = as.numeric(Date - min(Date[Timepoint == "01"], na.rm = TRUE))) %>%
  ungroup()

# Determine the maximum number of days to a sample for each patient
max_days <- df_combined %>%
  group_by(Patient) %>%
  summarize(Max_Days = max(Days_from_baseline, na.rm = TRUE)) %>%
  ungroup()

# Identify if the patient relapsed and the number of days to the relapse sample (Timepoint == "R")
relapse_info <- df_combined %>%
  filter(Timepoint == "R", !is.na(Days_from_baseline)) %>%
  group_by(Patient) %>%
  summarize(Relapse_Days = first(Days_from_baseline), Relapsed = 1) %>%
  ungroup()


# Merge the maximum days and relapse information
PFS_days <- max_days %>%
  left_join(relapse_info, by = "Patient") %>%
  mutate(
    Relapse_Days = ifelse(is.na(Relapse_Days), 0, Relapse_Days),
    Relapsed = ifelse(is.na(Relapsed), 0, Relapsed),
    Final_Days = ifelse(Relapsed == 1, Relapse_Days, Max_Days)  # Add a new column for the final days
  )

# Export the PFS_dates table to a CSV file
write.csv(PFS_days, file = "Exported_data_tables_clinical/PFS_days.csv", row.names = FALSE)

### Get the latest date per patient 
latest_dates <- df_combined %>%
  group_by(Patient) %>%
  summarise(latest_date = max(Date, na.rm = TRUE))

#### Add SPORE to this and IMG 
# Extract the Patient and 'Status last follow up' from spore_OS_info,
# renaming it to latest_date
tmp <- spore_OS_info %>%
  select(Patient, latest_date = `Status last follow up`)

# Combine the two tables
latest_dates <- bind_rows(latest_dates, tmp)

## Export this 
write.csv(latest_dates, file = "Exported_data_tables_clinical/latest_dates_per_patient.csv", row.names = FALSE)


# ─── 19.  Filter/compare patient sets for cfDNA vs. BM longitudinal → export counts ─

df <- combined_clinical_data_updated

# 19.1  Identify who has ≥2 distinct cfDNA timepoints (“longitudinal cfDNA”)
longitudinal_cfDNA <- df %>%
  filter(Sample_type=="Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  summarise(n_timepoints=n_distinct(Timepoint)) %>%
  filter(n_timepoints>=2) %>%
  pull(Patient)

# 19.2  Baseline/Diagnosis BM & cfDNA sets, Progression/Relapse BM & cfDNA sets
baseline_BM       <- df %>% filter(Sample_type=="BM_cells", timepoint_info=="Diagnosis")         %>% distinct(Patient) %>% pull(Patient)
baseline_cfDNA    <- df %>% filter(Sample_type=="Blood_plasma_cfDNA", timepoint_info=="Diagnosis") %>% distinct(Patient) %>% pull(Patient)
progression_BM    <- df %>% filter(Sample_type=="BM_cells", timepoint_info %in% c("Progression","Relapse")) %>% distinct(Patient) %>% pull(Patient)
progression_cfDNA <- df %>% filter(Sample_type=="Blood_plasma_cfDNA", timepoint_info %in% c("Progression","Relapse")) %>% distinct(Patient) %>% pull(Patient)

# 19.3  Helper to count intersections
count_intersect <- function(A,B) length(intersect(A,B))

n_a <- count_intersect(baseline_BM, longitudinal_cfDNA)
n_b <- count_intersect(baseline_cfDNA, longitudinal_cfDNA)
n_c <- count_intersect(progression_BM, longitudinal_cfDNA)
n_d <- count_intersect(progression_cfDNA, longitudinal_cfDNA)

no_baseline_BM    <- setdiff(progression_BM, baseline_BM)
n_e <- count_intersect(no_baseline_BM, longitudinal_cfDNA)

no_baseline_cfDNA <- setdiff(progression_cfDNA, baseline_cfDNA)
n_f <- count_intersect(no_baseline_cfDNA, longitudinal_cfDNA)

no_baseline_BM_cfDNA <- setdiff(baseline_cfDNA, baseline_BM)
n_g <- count_intersect(no_baseline_BM_cfDNA, longitudinal_cfDNA)

# 19.4  Build results + detailed lists → write CSVs
results <- tibble(
  description = c(
    "Baseline/diagnosis BM + longitudinal cfDNA",
    "Baseline/diagnosis cfDNA + longitudinal cfDNA",
    "Progression/relapse BM + longitudinal cfDNA",
    "Progression/relapse cfDNA + longitudinal cfDNA",
    "Prog/relapse BM BUT NOT baseline BM + longitudinal cfDNA",
    "Prog/relapse cfDNA BUT NOT baseline cfDNA + longitudinal cfDNA",
    "Baseline cfDNA BUT NOT baseline BM + longitudinal cfDNA"
  ),
  count = c(n_a, n_b, n_c, n_d, n_e, n_f, n_g)
)

find_intersect <- function(A,B) intersect(A,B)

patients_a <- find_intersect(baseline_BM, longitudinal_cfDNA)
patients_b <- find_intersect(baseline_cfDNA, longitudinal_cfDNA)
patients_c <- find_intersect(progression_BM, longitudinal_cfDNA)
patients_d <- find_intersect(progression_cfDNA, longitudinal_cfDNA)
patients_e <- find_intersect(setdiff(progression_BM, baseline_BM), longitudinal_cfDNA)
patients_f <- find_intersect(setdiff(progression_cfDNA, baseline_cfDNA), longitudinal_cfDNA)
patients_g <- find_intersect(no_baseline_BM_cfDNA, longitudinal_cfDNA)

detailed_results <- tibble(
  description = rep(results$description, times = c(length(patients_a), length(patients_b),
                                                   length(patients_c), length(patients_d),
                                                   length(patients_e), length(patients_f),
                                                   length(patients_g))),
  patient     = c(patients_a, patients_b, patients_c, patients_d, patients_e, patients_f, patients_g)
)

write.csv(results,          "patient_counts_summary.csv",  row.names = FALSE)
write.csv(detailed_results, "patient_counts_detailed.csv",  row.names = FALSE)


# ─── 20.  Filter BAM lists (IMG only, SPORE only) → write out CSVs ─────────────
bams_filtered_IMG <- combined_clinical_data_updated %>%
  filter(grepl("IMG", Patient, ignore.case = TRUE)) %>%
  select(Bam, Patient, Sample_ID)
write.csv(bams_filtered_IMG, "bams_with_IMG_patients.csv", row.names = FALSE)

bams_filtered_SPORE <- combined_clinical_data_updated %>%
  filter(grepl("SPORE", Patient, ignore.case = TRUE)) %>%
  select(Bam, Patient, Sample_ID)
write.csv(bams_filtered_SPORE, "SPORE_bams.csv", row.names = FALSE)

cat("Filtered BAM list exported to bams_with_IMG_patients.csv and SPORE_bams.csv\n")


# ─── 21.  Final tweaks & write “combined_clinical_data_updated_April2025.csv” ──
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(timepoint_info = if_else(
    Patient == "SPORE_0007" & Date_of_sample_collection == as.Date("2019-08-21"),
    "Post_transplant", timepoint_info
  ))

write.csv(combined_clinical_data_updated,
          "combined_clinical_data_updated_April2025.csv",
          row.names = FALSE)


# ─── 22.  Build “baseline_dates_all” (SPORE, M4, IMMAGINE) ───────────────────────
###  22.1  SPORE: pick the first BM “Baseline” for SPORE_0009, else first “Diagnosis”
baseline_0009 <- combined_spore_timepoint_info %>%
  filter(Patient=="SPORE_0009", timepoint_info=="Baseline") %>%
  slice_min(`Timepoint of interest`, with_ties=FALSE)

diagnosis_others <- combined_spore_timepoint_info %>%
  filter(Patient!="SPORE_0009", timepoint_info=="Diagnosis") %>%
  group_by(Patient) %>%
  slice_min(`Timepoint of interest`, with_ties=FALSE) %>%
  ungroup()

spore_baseline_dates <- bind_rows(baseline_0009, diagnosis_others)

###  22.2  M4: from M4_diagnosis_date (rename M4_id→Patient)
diagnosis_dates_m4 <- M4_diagnosis_date %>%
  select(Patient = M4_id, Diagnosis_date = DIAGNOSIS_DATE) %>%
  filter(!is.na(Diagnosis_date)) %>%
  distinct()

###  22.3  IMMAGINE: from IMMAGINE_progression → Diagnosis_date
diagnosis_dates_immagine <- IMMAGINE_progression %>%
  select(Patient = ID, Diagnosis_date) %>%
  filter(!is.na(Diagnosis_date)) %>%
  distinct()

###  22.4  Bind all three sources, drop NA
spore_baseline_clean <- spore_baseline_dates %>%
  select(Patient, Diagnosis_date = `Timepoint of interest`) %>%
  distinct()

baseline_dates_all <- bind_rows(
  spore_baseline_clean,
  diagnosis_dates_m4,
  diagnosis_dates_immagine
) %>%
  distinct() %>%
  filter(!is.na(Diagnosis_date))

n_distinct(baseline_dates_all$Patient)

dir.create("exported_clinical_data_April2025", showWarnings = FALSE)
write.csv(baseline_dates_all,
          file.path("exported_clinical_data_April2025", "all_baseline_dates.csv"),
          row.names = FALSE)



# ─── 23.  Get the patient info, number of patients with each feature ───────────────────────

# Define baseline and monitoring groups
baseline_timepoints <- c("Diagnosis", "Baseline")
monitoring_timepoints <- c("Post_transplant", "Post_induction", "Maintenance", "1.5yr maintenance", "2yr maintenance")

# Filter for baseline cfDNA samples
baseline_cfDNA <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% baseline_timepoints) %>%
  dplyr::select(Patient) %>%
  distinct()

# Count unique patients with baseline cfDNA samples
num_baseline_patients <- n_distinct(baseline_cfDNA$Patient)

# Filter for monitoring cfDNA samples
monitoring_cfDNA <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% monitoring_timepoints) %>%
  dplyr::select(Patient) %>%
  distinct()

# Count unique patients with monitoring cfDNA samples
num_monitoring_patients <- n_distinct(monitoring_cfDNA$Patient)

# Find patients with both baseline and monitoring samples
patients_with_both <- intersect(baseline_cfDNA$Patient, monitoring_cfDNA$Patient)
num_patients_with_both <- length(patients_with_both)

# Print results
cat("Number of patients with baseline cfDNA samples:", num_baseline_patients, "\n")
cat("Number of patients with monitoring cfDNA samples:", num_monitoring_patients, "\n")
cat("Number of patients with both baseline and monitoring cfDNA samples:", num_patients_with_both, "\n")

# Should I also keep only ones that were high in baseline to begin with? 
# I think best to keep ones with both 

# Number of cfDNA samples for patients with both baseline and monitoring cfDNA types
num_samples_with_both <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", Patient %in% patients_with_both) %>%
  nrow()

cat("Number of cfDNA samples for patients with both baseline and monitoring cfDNA types:", num_samples_with_both, "\n")

write.csv(patients_with_both, file = "patients_with_cfDNA_at_baseline_and_monitoring.csv")

samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", Patient %in% patients_with_both)

# Total unique patients with cfDNA samples
total_unique_patients <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  dplyr::select(Patient) %>%
  distinct() %>%
  nrow()

# Total number of cfDNA samples
total_cfDNA_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  nrow()

# Print results
cat("Total number of unique patients with cfDNA samples:", total_unique_patients, "\n")
cat("Total number of cfDNA samples:", total_cfDNA_samples, "\n")

# Number of baseline samples and patients
baseline_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% baseline_timepoints)

num_baseline_samples <- nrow(baseline_samples)
num_baseline_patients <- n_distinct(baseline_samples$Patient)

# Number of monitoring samples and patients
monitoring_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% monitoring_timepoints)

num_monitoring_samples <- nrow(monitoring_samples)
num_monitoring_patients <- n_distinct(monitoring_samples$Patient)

# Print results
cat("Number of baseline cfDNA samples:", num_baseline_samples, "\n")
cat("Number of unique patients with baseline cfDNA samples:", num_baseline_patients, "\n")
cat("Number of monitoring cfDNA samples:", num_monitoring_samples, "\n")
cat("Number of unique patients with monitoring cfDNA samples:", num_monitoring_patients, "\n")

# Combine patients from baseline and monitoring groups
patients_baseline <- baseline_samples$Patient
patients_monitoring <- monitoring_samples$Patient

# Union of baseline and monitoring patients
patients_with_either <- union(patients_baseline, patients_monitoring)

# Number of unique patients with either baseline or monitoring samples
num_patients_with_either <- length(patients_with_either)

# Print the result
cat("Number of unique patients with either baseline or monitoring cfDNA samples:", num_patients_with_either, "\n")

# Extract site identifiers from patients_with_either
site_identifiers <- sapply(patients_with_either, function(patient) {
  strsplit(patient, "[-_]")[[1]][1]  # Split by hyphen or underscore and take the first part
})

# Exclude IMG from the count
filtered_sites <- unique(site_identifiers[!grepl("^IMG", site_identifiers)])

# Number of unique sites
num_sites <- length(filtered_sites)

# Print the result
cat("Number of unique sites (excluding IMG):", num_sites, "\n")

write.csv(patients_with_either, file = "patients_with_either_cfDNA_at_baseline_and_monitoring.csv", row.names = FALSE)
