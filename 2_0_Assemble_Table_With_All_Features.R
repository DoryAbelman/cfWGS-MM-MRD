# ==============================================================================
# 2_Assemble_Table_With_All_Features
#
# Description:
#   Loads MRD assay results (flow‐cytometry, adaptive counts, Rapid Novor, PET),
#   clinical metadata (M4, SPORE, IMMAGINE), relapse dates, fragmentomics metrics,
#   WGS‐derived tumor fractions & structural features, and mutation count tables.
#   Cleans and harmonizes all inputs by Patient + Timepoint (deduplicating,
#   coalescing NAs, fixing sample codes), computes derived MRD flags, and merges
#   everything into one comprehensive table.  Finally, writes out CSV and RDS
#   versions for downstream analyses and figure generation.
#
# Author:    Dory Abelman
# Last update: May 2025
# ==============================================================================



### MRD comparison figures updated Winter 2025
### Author: Dory Abelman 

library(tidyverse)
library(readxl)
library(gridExtra)


### Load in all data:
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_May2025.rds") # MRD sample purity
M4_MRD_filtered <- read_tsv("M4_MRD_filtered.txt") # The MRD data from labs

path <- file.path("MRDetect_output_winter_2025", "Processed_R_outputs", "BM_muts_data")
#MRD_cfWGS_backup <- MRD_cfWGS
MRD_cfWGS_BM <- read.csv(file.path("MRDetect_output_winter_2025/Processed_R_outputs/BM_muts_plots_baseline/cfWGS MRDetect BM data updated May.csv")) # MRD data from MRDetect
MRD_cfWGS_blood <- read.csv(file.path("MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/cfWGS MRDetect Blood data updated June.csv")) # MRD data from MRDetect


Relapse_dates_M4_clean <- read_csv("Relapse_dates_M4_clean.csv") # Relapse dates
combined_clinical_data_updated <- read_csv("combined_clinical_data_updated_April2025.csv") ## Aggregated clinical info
SPORE_MRD <- read_excel("Clinical data/SPORE/SPORE_pct_flow_extracted.xlsx") %>% filter(!is.na(Patient))
IMMAGINE_MRD <- read_excel("Clinical data/IMMAGINE/Extracted_clinical_MRD_data.xlsx", sheet = 3)


### Make a full MRD table 
### Issue is M4 is by timepoint and others by date. Can do the others by date
## and then try to merge all together with the combined clinical data 
## or join in merged_table 
## need to clearly indicate where everything is - ie all key clinical files being used
# Step 1: Clean and transform `% Cells` in table2
table2 <- IMMAGINE_MRD %>%
  mutate(
    # Convert `% Cells` to numeric where "Not detected" is 0
    `% Cells` = ifelse(`% Cells` == "Not detected", 0, as.numeric(`% Cells`))
  )

table2 <- table2 %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))

# Step 2: Rename columns for consistency
table1 <- SPORE_MRD %>%
  rename(ID = Patient, Date = `Timepoint of interest`, Flow_Pct_Cells = `Edited_Pct_Flow`)

table2 <- table2 %>%
  rename(Date = Date, Flow_Pct_Cells = `% Cells`)

# Step 3: Merge the tables
# Join by 'ID' and 'Date'
MRD_SPORE_IMMAGINE <- bind_rows(table1, table2) %>% select(ID, Date, Flow_Pct_Cells, `MRD Status`)
MRD_SPORE_IMMAGINE <- MRD_SPORE_IMMAGINE %>%
  rename(MRD_Status = `MRD Status`) %>%
  mutate(
    MRD_Status = ifelse(is.na(Flow_Pct_Cells) | Flow_Pct_Cells > 0, "Positive", "Negative")
  )

M4_MRD_clinical <- M4_MRD_filtered %>%
  rename(MRD_Results_FLOW = MRD_Results)

## Now add to the other table 
MRD_SPORE_IMMAGINE <- MRD_SPORE_IMMAGINE %>% 
  rename(Patient = ID, Flow_pct_cells = Flow_Pct_Cells, MRD_Results_FLOW = MRD_Status)

## Multiple by 100 to be consistent with other table for FLOW 
MRD_SPORE_IMMAGINE$Flow_pct_cells <- MRD_SPORE_IMMAGINE$Flow_pct_cells*100

M4_MRD_clinical <- M4_MRD_clinical %>% 
  rename(Patient = M4_id)

cfWGS_Clinical_MRD <- bind_rows(M4_MRD_clinical, MRD_SPORE_IMMAGINE)

## Calculate MRD rates 
# Create MRD_by_clinical_testing column based on the specified conditions
cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(
    MRD_by_clinical_testing = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Mrd1E5 is Positive
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 1
      Rapid_Novor > 1 ~ "Positive",
      
      # If none of the conditions are met, set as Negative
      TRUE ~ "Negative"
    )
  )

# Create MRD_by_clinical_testing_stringent column based on the specified conditions
cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(
    MRD_by_clinical_testing_stringent = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Mrd1E5 is Positive
      Mrd1E6 == "Positive" ~ "Positive",
      
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 0.1
      Rapid_Novor > 0.1 ~ "Positive",
      
      # If none of the conditions are met, set as Negative
      TRUE ~ "Negative"
    )
  )


### Add dates to M4 for use in the patient tracker 
M4_dates <- read_csv("M4_dates_Oct2024.csv") 
tmp <- M4_dates %>% 
  select(Patient, Timepoint_Code, Date) %>% 
  unique() %>% 
  rename(Visit_Number = Timepoint_Code) %>% 
  mutate(Visit_Number = gsub("R-", "R", Visit_Number))

cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>% 
  left_join(tmp %>% select(Patient, Visit_Number, Date), by = c("Patient", "Visit_Number")) %>%
  mutate(Date = coalesce(Date.x, Date.y)) %>%
  select(-Date.x, -Date.y)


## Add the IDs if there is a matched sample 
tmp2 <- combined_clinical_data_updated %>% 
  #  filter(Study != "M4") %>%
  select(Patient, Date_of_sample_collection, Timepoint, Sample_ID, timepoint_info) %>% 
  unique() %>%
  mutate(
    Sample_Code = gsub("(-P|-B|-O|-OZ-DNA|-OZ|_BM_cells|_Blood_plasma_cfDNA|_Blood_Buffy_coat)$", "", Sample_ID),  # Remove specified suffixes
    Sample_Code = gsub("MyP|MyC|MyR", "IMG", Sample_Code)  # Remove specified suffixes
  ) %>%
  select(-Sample_ID) %>%  # Remove the original Sample_ID column
  unique()  # Ensure rows are unique

cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>% select(-Patient_Code, -Site_Code) %>% 
  rename(Timepoint = Visit_Number)

## Merge 
# Ensure dates are in Date format
cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(Date = as.Date(Date))

tmp2 <- tmp2 %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

# Join cfWGS_Clinical_MRD with tmp2 by Patient, matching dates within one month
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD %>% filter(is.na(Timepoint)) %>%
  left_join(tmp2, by = "Patient") %>%
  # Calculate the absolute difference in days
  mutate(date_diff = abs(difftime(Date, Date_of_sample_collection, units = "days"))) %>%
  # Filter to keep matches within 2 weeks
  filter(is.na(Timepoint.x) & date_diff <= 14 | !is.na(Timepoint.x)) %>%
  # For each row where Timepoint.x or Sample_Code.x is NA, keep the closest date match
  group_by(Patient, Date) %>%
  mutate(min_diff = min(date_diff, na.rm = TRUE)) %>%
  filter(date_diff == min_diff | !is.na(Timepoint.x)) %>%
  ungroup() %>%
  # Fill in missing Sample_Code and Timepoint from tmp2 where needed
  mutate(
    Sample_Code = ifelse(is.na(Sample_Code.x), Sample_Code.y, Sample_Code.x),
    Timepoint = ifelse(is.na(Timepoint.x), Timepoint.y, Timepoint.x)
  ) %>%
  # Select final columns and remove temporary columns
  select(-Sample_Code.x, -Sample_Code.y, -Timepoint.x, -Timepoint.y, -Date_of_sample_collection, -date_diff, -min_diff) %>% 
  unique()

## Bind rows to the other 

# Make sure both data frames have Sample_Code (and Timepoint) as character
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    Sample_Code = as.character(Sample_Code),
    Timepoint   = as.character(Timepoint)
  )

cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(
    Sample_Code = as.character(Sample_Code),
    Timepoint   = as.character(Timepoint)
  )
cfWGS_Clinical_MRD_filled <- bind_rows(cfWGS_Clinical_MRD_filled, cfWGS_Clinical_MRD %>% filter(!is.na(Timepoint)))

## Add timepoint info if NA
# 1) make a clean lookup table
tmp_info <- tmp2 %>%
  ungroup() %>% 
  select(Sample_Code, timepoint_info) %>%
  distinct()

# 2) join + coalesce
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  left_join(tmp_info, by = "Sample_Code", suffix = c("", ".new")) %>%
  mutate(
    timepoint_info = coalesce(timepoint_info, timepoint_info.new)
  ) %>%
  select(-timepoint_info.new)

#cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% select(-Sample_Code)

## Add timepoint values 
# Create the mapping of timepoints to descriptions
description_mapping <- c(
  "01" = "Diagnosis",
  "03" = "Post_induction",
  "05" = "Post_transplant",
  "07" = "1yr maintenance",
  "08" = "1.5yr maintenance",
  "09" = "2yr maintenance",
  "10" = "2.5yr maintenance",
  "11" = "3yr maintenance",
  "12" = "3.5yr maintenance",
  "13" = "4yr maintenance",
  "14" = "4.5yr maintenance",
  "15" = "5yr maintenance",
  "R"  = "Relapse"
)

# Add timepoint_info based on the mapping
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    timepoint_info = if_else(
      !grepl("^(SPORE|IMG)", Patient),  # Check if Patient does not start with SPORE or IMG
      description_mapping[Timepoint],  # Apply the mapping
      timepoint_info  # Keep the existing value otherwise
    )
  )


### Add the new Rapid Novor data got from Aimee later 
Rapid_Novor_update <- read_excel("Clinical data/M4/Updated Rapid Novor Data.xlsx")

# Restructure the data
# Step 1: Standardize and Clean Non-Numeric Values
Rapid_Novor_cleaned <- Rapid_Novor_update %>%
  mutate(across(
    -`UHN Serial`, 
    ~ case_when(
      . %in% c("Negative", "Not Quantifiable", "NQ") ~ 0,  # Set these to 0
      TRUE ~ suppressWarnings(as.numeric(.))  # Convert others to numeric
    )
  )) %>%
  rename(M4_id = `UHN Serial`) # Rename UHN Serial to M4_id

# Step 2: Restructure the Data
Rapid_Novor_restructured <- Rapid_Novor_cleaned %>%
  pivot_longer(
    cols = -M4_id, # All columns except M4_id
    names_to = "Timepoint", # New column for Timepoint
    values_to = "Value" # New column for the values
  ) 

# Step 1: Remove rows where Value is NA
Rapid_Novor_restructured <- Rapid_Novor_restructured %>%
  filter(!is.na(Value))

# Step 2: Consolidate 9 and 9P into 9, and R and RP into R
Rapid_Novor_consolidated <- Rapid_Novor_restructured %>%
  mutate(Timepoint = case_when(
    Timepoint %in% c("9", "9P") ~ "9",  # Consolidate 9 and 9P into 9
    Timepoint %in% c("R", "RP") ~ "R",  # Consolidate R and RP into R
    TRUE ~ Timepoint  # Keep other timepoints as they are
  )) %>%
  group_by(M4_id, Timepoint) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")  # Average if duplicates exist

## Prep for join 
Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    Patient = M4_id,  # Rename M4_id to Patient
  ) %>%
  select(-M4_id)  # Optionally remove the original M4_id column if not needed


Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    Patient = case_when(
      Patient %in% c("VA-018", "VA-017", "VA-016") ~ sub("^VA-0", "VA-", Patient),  # Remove leading 0
      TRUE ~ Patient  # Keep all other values unchanged
    )
  )

Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    Timepoint = if_else(Timepoint == "R", Timepoint, str_pad(Timepoint, width = 2, pad = "0")),  # Add leading 0 only if not "R"
    Sample_Code = paste0(Patient, "-", Timepoint)  # Create Sample_Code as Patient-Timepoint
  )

## Correct errors
# Create the mapping of timepoints to descriptions
description_mapping <- c(
  "01" = "Diagnosis",
  "03" = "Post_induction",
  "05" = "Post_transplant",
  "07" = "1yr maintenance",
  "08" = "1.5yr maintenance",
  "09" = "2yr maintenance",
  "10" = "2.5yr maintenance",
  "11" = "3yr maintenance",
  "12" = "3.5yr maintenance",
  "13" = "4yr maintenance",
  "14" = "4.5yr maintenance",
  "15" = "5yr maintenance",
  "R"  = "Relapse"
)

# Add timepoint_info based on the mapping
Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    timepoint_info = description_mapping[Timepoint]  # Map Timepoint to description
  ) 


# Check for unmapped timepoints
unmapped_timepoints <- Rapid_Novor_consolidated %>%
  filter(is.na(timepoint_info)) %>%
  pull(Timepoint) %>%
  unique()

if (length(unmapped_timepoints) > 0) {
  warning("Unmapped timepoints: ", paste(unmapped_timepoints, collapse = ", "))
}

Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  rename(Rapid_Novor = Value)

Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>% unique()

# Perform an inner join to compare Rapid_Novor values for matching Patient and Timepoint
# Set a tolerance threshold (e.g., 0.01 for small differences)
tolerance <- 1

# Perform an inner join to compare Rapid_Novor values for matching Patient and Timepoint
# Some samples are really different
discrepancies <- Rapid_Novor_consolidated %>%
  inner_join(
    cfWGS_Clinical_MRD_filled %>% select(Patient, Timepoint, Rapid_Novor),
    by = c("Patient", "Timepoint"),
    suffix = c("_consolidated", "_clinical")
  ) %>%
  filter(!is.na(Rapid_Novor_consolidated), !is.na(Rapid_Novor_clinical)) %>% # Ensure both values are non-missing
  mutate(
    difference = abs(Rapid_Novor_consolidated - Rapid_Novor_clinical), # Calculate absolute difference
    within_tolerance = difference <= tolerance # Check if difference is within tolerance
  )


## Now add them to the original table and overwrite original since I think was errors 
# Step 1: Remove the Rapid_Novor and Rapid_Novor_Binary columns from cfWGS_Clinical_MRD_filled
tmp <- cfWGS_Clinical_MRD_filled %>%
  select(-Rapid_Novor)

# Step 2: Check for duplicates in tmp
dups_by_pt <- tmp %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()


# Step 2: Full join with Rapid_Novor_consolidated to ensure all rows are included
cfWGS_Clinical_MRD_filled <- tmp %>%
  full_join(
    Rapid_Novor_consolidated %>% 
      select(Patient, Timepoint, Rapid_Novor, Sample_Code, timepoint_info)
  ) %>% 
  unique()

# Identify duplicate rows based on Sample_Code
duplicate_rows <- cfWGS_Clinical_MRD_filled %>%
  group_by(Sample_Code) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()



### Add in the PET results
PET_M4 <- read_excel("Clinical data/M4/PET_dataframe.xlsx")

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% 
  full_join(PET_M4 %>% select(PET, Patient_Timepoint), by = c("Sample_Code" = "Patient_Timepoint"))


## Clean up 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    MRD_Results_FLOW = case_when(
      Flow_pct_cells == 0 ~ "Negative",
      Flow_pct_cells > 0 & MRD_Results_FLOW == "-" ~ "Positive",
      TRUE ~ MRD_Results_FLOW  # Keep the original value if neither condition applies
    )
  )

## Recalculate rates
# Create MRD_by_clinical_testing column based on the specified conditions
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    MRD_by_clinical_testing = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Mrd1E5 is Positive
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 1
      Rapid_Novor > 1 ~ "Positive",
      
      # If none of the conditions are met, set as Negative
      TRUE ~ "Negative"
    )
  )

# Create MRD_by_clinical_testing_stringent column based on the specified conditions
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    MRD_by_clinical_testing_stringent = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Mrd1E5 is Positive
      Mrd1E6 == "Positive" ~ "Positive",
      
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 0.1
      Rapid_Novor > 0.05 ~ "Positive",
      
      # If none of the conditions are met, set as Negative
      TRUE ~ "Negative"
    )
  )

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()
### Now assess correlates

# Prepare binary columns for concordance analysis with is.na checks
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    # Binary indicator for Flow: 1 if Positive, 0 if Negative, NA if missing
    Flow_Binary = case_when(
      !is.na(MRD_Results_FLOW) & MRD_Results_FLOW == "Positive" ~ 1,
      !is.na(MRD_Results_FLOW) & MRD_Results_FLOW == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for Adaptive: 1 if Mrd1E6 or Mrd1E5 is Positive, 0 if Negative, NA if missing
    Adaptive_Binary = case_when(
      (!is.na(Mrd1E6) & Mrd1E6 == "Positive") | (!is.na(Mrd1E5) & Mrd1E5 == "Positive") ~ 1,
      (!is.na(Mrd1E6) & Mrd1E6 == "Indeterminate") & (!is.na(Mrd1E5) & Mrd1E5 == "Negative") ~ 0,
      (!is.na(Mrd1E6) & Mrd1E6 == "Negative") | (!is.na(Mrd1E5) & Mrd1E5 == "Negative") ~ 0,
      (!is.na(Mrd1E6) & Mrd1E6 == "Positive") ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for Rapid Novor if positive (> 0), 0 if not, NA if missing
    Rapid_Novor_Binary = case_when(
      !is.na(Rapid_Novor) & Rapid_Novor > 0 ~ 1,
      !is.na(Rapid_Novor) & Rapid_Novor <= 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for PET: 1 if Positive, 0 if Negative, NA if missing
    PET_Binary = case_when(
      !is.na(PET) & PET == "Positive" ~ 1,
      !is.na(PET) & PET == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for MRD by clinical testing (standard): 1 if Positive, 0 if Negative, NA if missing
    MRD_Clinical_Binary = case_when(
      !is.na(MRD_by_clinical_testing) & MRD_by_clinical_testing == "Positive" ~ 1,
      !is.na(MRD_by_clinical_testing) & MRD_by_clinical_testing == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for MRD by clinical testing (stringent): 1 if Positive, 0 if Negative, NA if missing
    MRD_Clinical_Stringent_Binary = case_when(
      !is.na(MRD_by_clinical_testing_stringent) & MRD_by_clinical_testing_stringent == "Positive" ~ 1,
      !is.na(MRD_by_clinical_testing_stringent) & MRD_by_clinical_testing_stringent == "Negative" ~ 0,
      TRUE ~ NA_real_
    )
  )

### Add cfWGS
# --- 1. Make a “BM” lookup table -----------------------------

tmp_bm <- MRD_cfWGS_BM %>%
  # (you already filtered to baseline timepoints and selected only
  # the columns you care about)
  select(
    Patient,
    Date  = Date_of_sample_collection_Sample_ID_Bam,
    Timepoint = Timepoint_Sample_ID_Bam,
    Mrd_by_WGS_BM    = Mrd_by_WGS,
    Cumulative_VAF_BM        = Cumulative_VAF,
    detect_rate_BM   = detection_rate_as_reads_detected_over_reads_checked,
    zscore_BM        = sites_rate_zscore_charm,
    z_score_detection_rate_BM = detection_rate_zscore_reads_checked_charm, 
    PercentChangeFromBaseline_BM       = percent_change_detection_rate, #percent_change
    PercentChangeAtSecondTimepoint_BM       = percent_change_detection_rate_second_timepoint
  ) %>%
  distinct()                # just in case you have duplicates

# --- 2. Make a “Blood” lookup table -------------------------

tmp_blood <- MRD_cfWGS_blood %>%
  # keep only the same columns you want, but rename for clarity
  select(
    Patient,
    Date  = Date_of_sample_collection_Sample_ID_Bam,
    Timepoint = Timepoint_Sample_ID_Bam,
    Mrd_by_WGS_blood = Mrd_by_WGS,
    Cumulative_VAF_blood     = Cumulative_VAF,
    detect_rate_blood= detection_rate_as_reads_detected_over_reads_checked,
    zscore_blood     = sites_rate_zscore_charm,
    z_score_detection_rate_blood = detection_rate_zscore_reads_checked_charm, 
    PercentChangeFromBaseline_blood    = percent_change,
    PercentChangeAtSecondTimepoint_blood    = percent_change_detection_rate_second_timepoint
  ) %>%
  distinct()

#  Parse the Date columns in your look-up tables
tmp_bm <- tmp_bm %>%
  mutate(Date = as.Date(Date)) 

tmp_blood <- tmp_blood %>%
  mutate(Date = as.Date(Date))

### check for duplicates 
# 1. Find any Patient/Date/Timepoint combos that appear >1× in tmp_bm or tmp_blood
bm_dups <- tmp_bm %>% 
  count(Patient, Date, Timepoint) %>% 
  filter(n > 1)

blood_dups <- tmp_blood %>% 
  count(Patient, Date, Timepoint) %>% 
  filter(n > 1)

bm_dups
blood_dups


# --- 3. Join both into your clinical table ------------------
# Ensure all Date columns are of type Date
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(Date = as.Date(Date))

tmp_bm <- tmp_bm %>%
  mutate(Date = as.Date(Date))

tmp_blood <- tmp_blood %>%
  mutate(Date = as.Date(Date))
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  # first bring in the BM calls
  full_join(tmp_bm,    by = c("Patient","Date","Timepoint")) %>%
  # then bring in the blood calls
  full_join(tmp_blood, by = c("Patient","Date","Timepoint"))

# --- 4. Create two binary MRD flags -------------------------

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    cfWGS_BM_Binary = case_when(
      Mrd_by_WGS_BM    == "Positive" ~ 1,
      Mrd_by_WGS_BM    == "Negative" ~ 0,
      TRUE                          ~ NA_real_
    ),
    cfWGS_blood_Binary = case_when(
      Mrd_by_WGS_blood == "Positive" ~ 1,
      Mrd_by_WGS_blood == "Negative" ~ 0,
      TRUE                          ~ NA_real_
    )
  )

# --- 5. Quick sanity-check -------------------------

# How many rows got both BM & blood info?
cfWGS_Clinical_MRD_filled %>% 
  summarize(
    total = dplyr::n(),
    with_BM    = sum(!is.na(Mrd_by_WGS_BM)),
    with_blood = sum(!is.na(Mrd_by_WGS_blood))
  )


#### Reorganize
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  select(
    # 1) Core identifiers
    Patient,
    Date,
    timepoint_info,
    Sample_Code,
    Timepoint,
    
    # 2) Clinical MRD assays
    Flow_pct_cells,
    Flow_Binary,
    MRD_Results_FLOW,
    
    PerMillionCount_Adaptive,
    Adaptive_Frequency,
    Adaptive_Binary,
    
    Rapid_Novor,
    Rapid_Novor_Binary,
    
    PET,
    PET_Binary,
    
    MRD_by_clinical_testing,
    MRD_by_clinical_testing_stringent,
    MRD_Clinical_Binary,
    MRD_Clinical_Stringent_Binary,
    
    # 3) WGS‐based MRD, BM mutations
    Mrd_by_WGS_BM,
    cfWGS_BM_Binary,
    Cumulative_VAF_BM,
    detect_rate_BM,
    zscore_BM,
    z_score_detection_rate_BM,
    PercentChangeFromBaseline_BM,
    PercentChangeAtSecondTimepoint_BM,
    
    # 4) WGS‐based MRD, blood mutations
    Mrd_by_WGS_blood,
    cfWGS_blood_Binary,
    Cumulative_VAF_blood,
    detect_rate_blood,
    zscore_blood,
    z_score_detection_rate_blood,
    PercentChangeFromBaseline_blood,
    PercentChangeAtSecondTimepoint_blood,
    
    # 5) Anything else you might have missed
    everything()
  )


cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()


## Get duplicates
# 1) Pull out all the rows where (Patient, Date) is duplicated
dups <- cfWGS_Clinical_MRD_filled %>%
  group_by(Patient, Date, Sample_Code) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()

# Inspect
dups %>% 
  select(Patient, Date, Sample_Code) %>%
  arrange(Patient, Date)

## Remove and consolidate 
#cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
#  filter(
#    is.na(Sample_ID) | !str_detect(Sample_ID, "(_Blood_Buffy_coat|_BM_cells|-B$|-O$)")
#  )

### Add in additional dates for missing cases 
M4_processing_log_dates <- read.csv("Exported_data_tables_clinical/M4_processing_log_dates.csv")

M4_processing_log_dates <- M4_processing_log_dates %>%
  mutate(Sample_Code = paste0(Patient, "-", Timepoint))

# Compare dates for matching Sample_Codes and calculate days difference
date_comparison <- cfWGS_Clinical_MRD_filled %>%
  # make sure Sample_Code is character, not logical
  mutate(Sample_Code = as.character(Sample_Code)) %>%
  select(Sample_Code, Date) %>%
  inner_join(
    M4_processing_log_dates %>%
      # likewise ensure it’s character
      mutate(Sample_Code = as.character(Sample_Code)) %>%
      select(Sample_Code, Cleaned_Date),
    by = "Sample_Code"
  ) %>%
  mutate(
    # coerce your cleaned date to Date
    Cleaned_Date   = as.Date(Cleaned_Date),
    # calculate days difference
    Days_Difference = as.numeric(Date - Cleaned_Date)
  ) %>%
  # drop any zero or NA differences
  filter(!is.na(Days_Difference) & Days_Difference != 0)

# Find duplicate date entries
duplicate_dates <- M4_processing_log_dates %>%
  group_by(Sample_Code) %>%  # Group by the sample code column
  filter(dplyr::n() > 1)  # Keep only rows where the Cleaned_Date appears more than once

# Keep only those closest to lab 
# Merge the two datasets to compare dates
M4_processing_log_dates <- M4_processing_log_dates %>%
  mutate(Cleaned_Date = as.Date(Cleaned_Date))

## From previous script
M4_Labs <- read.csv(file = "M4_labs_cleaned.csv") %>% # generated in script 1_0
  mutate(LAB_DATE = as.Date(LAB_DATE))

# Merge the two datasets with a left join to retain all processing log dates
merged_data <- M4_processing_log_dates %>%
  left_join(
    M4_Labs %>%
      select(Patient, TIMEPOINT, LAB_DATE), 
    by = c("Patient" = "Patient", "Timepoint" = "TIMEPOINT")
  )

# Calculate the absolute difference in days (assign NA for rows without LAB_DATE)
merged_data <- merged_data %>%
  mutate(Date_Difference = ifelse(is.na(LAB_DATE), NA, abs(as.numeric(Cleaned_Date - LAB_DATE))))

# For each Sample_Code, keep the row with the smallest Date_Difference (or retain rows without LAB_DATE)
closest_date_entries <- merged_data %>%
  group_by(Sample_Code) %>%
  filter(is.na(Date_Difference) | Date_Difference == min(Date_Difference, na.rm = TRUE)) %>%
  ungroup() %>%
  unique()

# Keep only relevant columns
final_processing_log <- closest_date_entries %>%
  select(Patient, Timepoint, Sample_Code, Cleaned_Date) %>%
  unique()

# Find duplicate date entries
duplicate_dates <- final_processing_log %>%
  group_by(Sample_Code) %>%  # Group by the sample code column
  filter(dplyr::n() > 1)  # Keep only rows where the Cleaned_Date appears more than once


# Define the rows to exclude
duplicates_to_exclude <- data.frame(
  Sample_Code = c("MJ-03-07", "VA-05-05", "ZC-06-01"),
  Cleaned_Date = as.Date(c("2020-11-04", "2019-12-05", "2021-11-12"))
)

# Filter out the duplicates from the main dataframe
M4_processing_log_dates <- M4_processing_log_dates %>%
  anti_join(duplicates_to_exclude, by = c("Sample_Code", "Cleaned_Date"))

## Recode VA-07 since 05 was actually relapse timepoint 
M4_processing_log_dates <- M4_processing_log_dates %>%
  mutate(
    # define your condition once
    is_va07_05 = Patient == "VA-07" & Timepoint == "05",
    
    # change Timepoint where needed
    Timepoint  = if_else(is_va07_05, "R-", Timepoint),
    
    # change Sample_Code where needed
    Sample_Code = if_else(is_va07_05, "VA-07-R", Sample_Code)
  ) %>%
  select(-is_va07_05)

## Take the Date value for labs since may have been typos in processing log 
# Update NA dates or sample codes in cfWGS_Clinical_MRD_filled with values from M4_processing_log_dates
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  # 1) fix types in your main table
  mutate(
    Patient     = as.character(Patient),
    Timepoint   = as.character(Timepoint),
    Sample_Code = as.character(Sample_Code)
  ) %>%
  left_join(
    M4_processing_log_dates %>%
      # 2) also fix types in your lookup
      mutate(
        Patient     = as.character(Patient),
        Timepoint   = as.character(Timepoint),
        Sample_Code = as.character(Sample_Code)
      ) %>%
      # if there really are duplicate log rows, pick one
      distinct(Patient, Timepoint, .keep_all = TRUE) %>%
      select(Patient, Timepoint, Cleaned_Date, Sample_Code),
    by     = c("Patient", "Timepoint"),
    suffix = c("", ".log")   # keeps your original Sample_Code, brings in Sample_Code.log
  ) %>%
  mutate(
    # 3a) fill in Date from Cleaned_Date only when Date is NA
    Date = coalesce(Date, as.Date(Cleaned_Date)),
    
    # 3b) fill in Sample_Code from the log when the original is missing/empty
    Sample_Code = if_else(
      is.na(Sample_Code) | Sample_Code == "",
      Sample_Code.log,
      Sample_Code
    )
  ) %>%
  # 4) drop the helper columns
  select(-Cleaned_Date, -Sample_Code.log)



# Convert numeric dates back to proper Date format
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(Date = as.Date(Date, origin = "1970-01-01"))


### Add sample codes where missing 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% filter(!is.na(Patient))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    Sample_Code = if_else(
      # only rewrite when it's missing or blank
      is.na(Sample_Code) | Sample_Code == "",
      case_when(
        # SPORE → “PATIENT_T<Timepoint>”
        str_detect(Patient, "^SPORE") ~ paste0(Patient, "_T", Timepoint),
        # IMG → “PATIENT-<Timepoint>”
        str_detect(Patient, "^IMG")   ~ paste0(Patient, "-", Timepoint),
        # all others → default “PATIENT-<Timepoint>”
        TRUE                          ~ paste0(Patient, "-", Timepoint)
      ),
      Sample_Code
    )
  )

## Make corrections 

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    # 1) fix any Timepoint that is exactly "R-"
    Timepoint = str_replace(Timepoint, "^R-$", "R"),
    
    # 2) fix any Patient ending in "R-"
    Sample_Code   = str_replace(Sample_Code,   "R-$",    "R")
  )

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()

## Add missing timepoint_info 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  # 1) compute the “new” labels based purely on Timepoint (and that special CA-08 tweak)
  mutate(
    new_timepoint_info = case_when(
      Timepoint == "01"                       ~ "Diagnosis",
      Timepoint == "03"                       ~ "Post_induction",
      Timepoint == "05"                       ~ "Post_transplant",
      Timepoint == "07"                       ~ "Maintenance",
      Timepoint == "08"                       ~ "1.5yr maintenance",
      Timepoint == "09" & Patient == "CA-08"  ~ "1.5yr maintenance",  # CA-08 special
      Timepoint == "09"                       ~ "2yr maintenance",
      Timepoint == "10"                       ~ "2.5yr maintenance",
      Timepoint == "11"                       ~ "3yr maintenance",
      Timepoint == "12"                       ~ "3.5yr maintenance",
      Timepoint == "13"                       ~ "4yr maintenance",
      Timepoint == "14"                       ~ "4.5yr maintenance",
      Timepoint == "15"                       ~ "5yr maintenance",
      Timepoint == "R"                        ~ "Relapse",
      TRUE                                    ~ NA_character_
    )
  ) %>%
  # 2) only overwrite the old timepoint_info when it was NA
  mutate(
    timepoint_info = coalesce(timepoint_info, new_timepoint_info)
  ) %>%
  # 3) drop the helper column
  select(-new_timepoint_info)



# Find duplicate date entries
duplicate_dates <- cfWGS_Clinical_MRD_filled %>%
  group_by(Sample_Code) %>%  # Group by the sample code column
  filter(dplyr::n() > 1)  # Keep only rows where the Cleaned_Date appears more than once

# Find sample codes with NA in the Date column
missing_date_samples <- cfWGS_Clinical_MRD_filled %>%
  dplyr::filter(is.na(Date)) %>%
  distinct(Sample_Code)

# Export the missing date sample codes to a CSV
write.csv(
  missing_date_samples,
  file = "Exported_data_tables_clinical/missing_date_samples_May2025_2.csv",
  row.names = FALSE
)



### Need to add relapse info  
### After redo confirmation for IMG-098

## Export relapse dates 

Relapse_dates_full <- read_csv("Exported_data_tables_clinical/Relapse dates cfWGS updated.csv") 

# Step 2: Calculate Num_days_to_closest_relapse
# Join on Patient and calculate the days difference
Time_to_relapse <- cfWGS_Clinical_MRD_filled %>% select(Patient, Date) %>% unique() %>%
  left_join(
    Relapse_dates_full,
    by = c("Patient")
  ) %>%
  group_by(Patient, Date) %>%
  mutate(
    # Calculate non-absolute days to relapse only for valid dates
    days_to_relapse_non_absolute = case_when(
      Progression_date > Date ~ as.numeric(Date - Progression_date), 
      Progression_date <= Date & Progression_date >= Date - 35 ~ 0,
      TRUE ~ NA_real_  # Set to NA for dates that do not meet conditions
    )
  ) %>%
  filter(!is.na(days_to_relapse_non_absolute)) %>%
  summarise(
    # Calculate the absolute days to relapse for valid dates
    Num_days_to_closest_relapse_absolute = min(abs(days_to_relapse_non_absolute), na.rm = TRUE),
    Num_days_to_closest_relapse = Num_days_to_closest_relapse_absolute,
    # Get the closest relapse date that came strictly after the sample collection
    Num_days_to_closest_relapse_non_absolute = days_to_relapse_non_absolute[which.min(abs(days_to_relapse_non_absolute))]
  ) %>%
  ungroup() %>%
  unique()

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  left_join(Time_to_relapse %>% select(Patient, Date, Num_days_to_closest_relapse) %>% unique())

## Add a column for Relapsed 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% 
  mutate(Relapsed = ifelse(!is.na(Num_days_to_closest_relapse), "Yes", "No"))

## Some editing 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(timepoint_info = ifelse(Sample_Code == "IMG-060-T36", "Maintenance", timepoint_info))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(timepoint_info = ifelse(Timepoint == "R", "Relapse", timepoint_info))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(timepoint_info = ifelse(Sample_Code == "IMG-159-T11", "Maintenance", timepoint_info))


### Now get the relapse date as the time since the baseline marrow and use this for the PFS plot


### Export 
write.csv(cfWGS_Clinical_MRD_filled, file = "cfWGS clinical MRD values with timepoint and dates updated May 2025.csv", row.names = F)




#### Updated May 2025 - now add lab values like free kappa, lambda, M-protein, ext, and use in the final classifier
#### Also add the fragmentomicss info 

#### First add in fragmentomics 

# 1) Read in your fragmentomics + cfWGS data
##### Ensure is correct with everything
frag <- read_csv("Results_Fragmentomics/Key_fragmentomics_data_updated2.csv")
## Filter the frag df 

# 1) Load fragmentomics + cfWGS features, keeping only the new columns you need
# 2) Keep only the fragmentomics‐specific metrics plus join keys
frag_small <- frag %>%
  select(
    Patient,
    Date_of_sample_collection,
    FS,
    Proportion.Short,
    Site,
    Mean.Coverage,
    Midpoint.Coverage,
    Midpoint.normalized,
    Amplitude,
    Zscore.Coverage,
    Zscore.Midpoint,
    Zscore.Amplitude,
    Threshold.Coverage,
    Threshold.Midpoint,
    Threshold.Amplitude
  ) %>%
  rename(
    Date        = Date_of_sample_collection
  )

frag_small <- frag_small %>% filter(!is.na(Patient))
# 2) If you only want one row per Patient/Date/Sample_Code—drop any duplicates
# 2) Collapse duplicates by averaging all numeric columns
frag_small_unique <- frag_small %>%
  group_by(Patient, Date) %>%
  summarise(
    across(
      c(
        FS,
        Proportion.Short,
        Mean.Coverage,
        Midpoint.Coverage,
        Midpoint.normalized,
        Amplitude,
        Zscore.Coverage,
        Zscore.Midpoint,
        Zscore.Amplitude,
        Threshold.Coverage,
        Threshold.Midpoint,
        Threshold.Amplitude
      ),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )



# 3) Join onto your clinical MRD table
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  left_join(
    frag_small_unique
  )





##### Next add other info like clinical data
#### Clinical data 
clinical_data_integrate <- read.csv("Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv") ### Generated in 1_1B

# make sure both date columns are Date class
clinical_data_integrate <- clinical_data_integrate %>%
  mutate(Date = as.Date(Date))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(Date = as.Date(Date))

## See duplicates 
dups <- clinical_data_integrate %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()

## Add the sample codes 
clinical_data_integrate <- clinical_data_integrate %>%
  mutate(
    Sample_Code =
      case_when(
        # SPORE → “PATIENT_T<Timepoint>”
        str_detect(Patient, "^SPORE") ~ paste0(Patient, "_T", Timepoint),
        # IMG → “PATIENT-<Timepoint>”
        str_detect(Patient, "^IMG")   ~ paste0(Patient, "-", Timepoint),
        # all others → default “PATIENT-<Timepoint>”
        TRUE                          ~ paste0(Patient, "-", Timepoint)
      ),
    Sample_Code
  )


# rather than fuzzy join, join by study timepoint 
joined <- clinical_data_integrate %>%
  full_join(
    cfWGS_Clinical_MRD_filled,
    by = c("Patient" = "Patient", "Sample_Code" = "Sample_Code"))

# now coalesce the two Timepoint / timepoint_info pairs,
# and drop the .y (cfWGS) helper columns
joined <- joined %>%
  # 1) Coalesce every .x/.y pair into a single column
  {
    # find all the bases that have both .x and .y
    xy_bases <- intersect(
      sub("\\.x$", "", grep("\\.x$", names(.), value = TRUE)),
      sub("\\.y$", "", grep("\\.y$", names(.), value = TRUE))
    )
    
    # for each base, create a single coalesced column
    coalesced <- map_dfc(xy_bases, function(base) {
      xcol <- paste0(base, ".x")
      ycol <- paste0(base, ".y")
      tibble(!!base := coalesce(.[[ycol]], .[[xcol]]))
    })
    
    # drop all the old .x/.y and bind the new coalesced ones
    select(., -matches("\\.(x|y)$")) %>%
      bind_cols(coalesced)
  }

## Fill demo info 
# pull out one row per patient with their non-NA demographics
demo_map <- joined %>%
  select(Patient, AGE, Gender, AGE_GROUP) %>%
  distinct(Patient, .keep_all = TRUE)

# 2) drop the old (possibly‐sparse) demo columns, and re-attach the full map
cfWGS_Clinical_MRD_filled_final <- joined %>%
  select(-AGE, -Gender, -AGE_GROUP) %>%
  left_join(demo_map, by = "Patient")


### The dates don't match exactly - strange...
cfWGS_Clinical_MRD_filled_final <- cfWGS_Clinical_MRD_filled_final %>% unique() 
dups <- cfWGS_Clinical_MRD_filled_final %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()

write.csv(dups, file = "duplicate_check.csv")
### Now add in the tumor fraction, other things identified and CNA from WGS

### Dedup the table to have one row per patient-timepoint combo
## helper: coalesce down an arbitrary vector (character OR numeric)
first_non_na <- function(x) {
  # treat "" and "-" as NA too
  x[x %in% c("", "-")] <- NA
  # return the first non‑NA (or NA if none)
  x[which(!is.na(x))[1]]
}

cfWGS_dedup <- cfWGS_Clinical_MRD_filled_final %>%
  # ------- 1. add priority flags --------------------------------------------
mutate(
  has_z     = !is.na(zscore_BM) | !is.na(zscore_blood),
  relapse_y = str_to_upper(Relapsed) == "YES"
) %>%
  
  # ------- 2. sort so best rows come first ----------------------------------
arrange(Patient, Timepoint,
        desc(has_z),          # 1st priority: has a z‑score, date used for rest of analysis 
        desc(relapse_y)) %>%  # 2nd priority: relapsed == "Yes", since this was not held in some older tables
  
  # ------- 3. collapse each Patient‑Timepoint group -------------------------
group_by(Patient, Timepoint) %>%
  summarise(
    across(
      # drop helper columns before summarising
      -c(has_z, relapse_y),
      first_non_na,
      .names = "{.col}"
    ),
    .groups = "drop"
  )


dups <- cfWGS_dedup %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()



#### Now add info from all feature data 
# 1. From your WGS feature table, select the key and all of the WGS columns you care about
wgs_subset <- All_feature_data %>%
  select(
    Patient,
    Date_of_sample_collection, 
    Sample_type,
    Timepoint,
    Tumor_Fraction,
    # CNAs
    del1p, amp1q, del13q, del17p, hyperdiploid,
    # IGH translocations
    IGH_MAF, IGH_CCND1, IGH_MYC, IGH_FGFR3,
    # mutation summary
    Mut_identified, Mut_genes, Mut_highest_VAF, Mut_type, Evidence_of_Disease
  )

# 2. Rename them with a “WGS_” prefix so they don’t overwrite your clinical calls
wgs_prefixed <- wgs_subset %>%
  rename_with(
    ~ paste0("WGS_", .),
    -c(Patient, Date_of_sample_collection, Sample_type, Timepoint)
  )

## Replace NA to 0
wgs_prefixed <- wgs_prefixed %>%
  mutate(across(
    all_of(c(
      "WGS_del1p", "WGS_amp1q", "WGS_del13q", "WGS_del17p",
      "WGS_IGH_MAF", "WGS_IGH_CCND1", "WGS_IGH_MYC", "WGS_IGH_FGFR3"
    )),
    ~ replace_na(.x, 0)
  ))

## Make wide 

## See what is dup
wgs_prefixed %>% 
  filter(Sample_type %in% c("BM_cells","Blood_plasma_cfDNA")) %>%
  group_by(Patient, Timepoint, Date_of_sample_collection, Sample_type) %>%
  filter(dplyr::n()>1) %>%
  ungroup() -> dups

# how many duplicate groups?
dups %>% 
  summarise(count = dplyr::n(), .by = c(Patient, Timepoint, Date_of_sample_collection, Sample_type)) %>% 
  arrange(desc(count))

## Correct 
wgs_prefixed <- wgs_prefixed %>%
  filter(Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA")) %>%
  group_by(Patient, Timepoint, Date_of_sample_collection, Sample_type) %>%
  summarise(
    # 1) mean of Tumor_Fraction
    WGS_Tumor_Fraction = mean(WGS_Tumor_Fraction, na.rm = TRUE),
    # 2) first() of all the other WGS_ columns
    across(
      .cols = starts_with("WGS_") & 
        !all_of("WGS_Tumor_Fraction"),
      .fns  = first
    ),
    .groups = "drop"
  )

# Normalize timepoints
wgs_prefixed <- wgs_prefixed %>%
  mutate(
    Timepoint = str_replace(Timepoint, "^R-$", "R")
  )

wgs_wide <- wgs_prefixed %>%
  # 1. keep only BM_cells and Blood_plasma_cfDNA
  filter(Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA")) %>%
  
  # 2. pivot so each WGS_… column splits into two: one per sample type
  pivot_wider(
    id_cols        = c(Patient, Timepoint, Date_of_sample_collection),
    names_from     = Sample_type,
    values_from    = starts_with("WGS_"),
    names_sep      = "_"
  )


# 3. Join onto  clinical MRD table by matching modifiers 
# make sure both date columns are Date objects
# Normalize the timepoints 
cfWGS_Clinical_MRD_filled_final <- cfWGS_dedup %>%
  mutate(
    Date      = as.Date(Date),
    Timepoint = str_replace(Timepoint, "^R-$", "R"),
    Sample_Code = str_replace(Sample_Code, "R-$", "R")
  )


wgs_wide <- wgs_wide %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

# fuzzy‐join on Patient, Sample_type, Timepoint (exact) 
# and Date within ±14 days
joined <- cfWGS_Clinical_MRD_filled_final %>%
  fuzzy_full_join(
    wgs_wide,
    by = c(
      "Patient" = "Patient",
      "Timepoint" = "Timepoint",
      "Date" = "Date_of_sample_collection"
    ),
    match_fun = list(
      `==`,            # Patient must match exactly
      `==`,            # Timepoint must match exactly
      function(d1, d2) abs(d1 - d2) <= 7  # dates within 7 days
    )
  )  


## reorganize
# 1. coalesce + drop the .x/.y duplicates
joined2 <- joined %>%
  mutate(
    Patient   = coalesce(Patient.x, Patient.y),
    Timepoint = coalesce(Timepoint.x, Timepoint.y)
  ) %>%
  select(-Patient.x, -Patient.y, -Timepoint.x, -Timepoint.y)

# 2. define your blocks
id_cols      <- c("Patient", "Timepoint", "Sample_Code", "timepoint_info",
                  "Date", "Date_of_sample_collection", "AGE", "Gender", "AGE_GROUP")

ext_cols     <- c("DIAGNOSIS_DATE","ECOG_SCORE","ISS_STAGE","KPS_SCORE",
                  "DETAILS","ISS_stage_notes","R_ISS_STAGE",
                  "Cytogenetic_Risk","Abnormality_Count","Abnormality_Status")

lab_cols     <- c("Study","dFLC","Albumin","B2_micro","Calcium","Creatinine",
                  "Hemoglobin","LDH","Plasma_pct","Plasma_count",
                  "IgA","IgG","IgM","Kappa","Kappa_Lambda_Ratio",
                  "Lambda","M_Protein","Subtype")

fish_cols    <- c("T_4_14","T_11_14","T_14_16",
                  "DEL_17P","DEL_1P","AMP_1Q","Fish_abnormalities",
                  "DEL_13","DEL_14",
                  "trisomy_3","trisomy_7","trisomy_9","trisomy_15",
                  "hyperdiploid","T_14_20","CKS1B_DUP")

flow_cols    <- c("Flow_pct_cells","Flow_Binary","MRD_Results_FLOW")
adaptive_cols<- c("PerMillionCount_Adaptive","Adaptive_Frequency","Adaptive_Binary")
rapid_cols   <- c("Rapid_Novor","Rapid_Novor_Binary")
pet_cols     <- c("PET","PET_Binary")

clinical_mrd_cols <- c(
  "MRD_by_clinical_testing","MRD_by_clinical_testing_stringent",
  "MRD_Clinical_Binary","MRD_Clinical_Stringent_Binary"
)

wgs_mrd_cols <- c(
  "Mrd_by_WGS_BM","cfWGS_BM_Binary","Cumulative_VAF_BM",
  "detect_rate_BM","zscore_BM", "z_score_detection_rate_BM",
  "PercentChangeFromBaseline_BM","PercentChangeAtSecondTimepoint_BM",
  "Mrd_by_WGS_blood","cfWGS_blood_Binary","Cumulative_VAF_blood",
  "detect_rate_blood","zscore_blood", "z_score_detection_rate_blood",
  "PercentChangeFromBaseline_blood","PercentChangeAtSecondTimepoint_blood",
  "Mrd1E5","Mrd1E6"
)

coverage_cols <- c(
  "FS", "Proportion.Short",                    
  "Mean.Coverage","Midpoint.Coverage","Midpoint.normalized","Amplitude",
  "Zscore.Coverage","Zscore.Midpoint","Zscore.Amplitude",
  "Threshold.Coverage","Threshold.Midpoint","Threshold.Amplitude"
)

# grab the remaining WGS_… columns
wgs_cols <- grep("^WGS_", colnames(joined2), value = TRUE)

# 3. put it all together
joined_clean <- joined2 %>%
  select(
    # identifiers & dates
    all_of(id_cols),
    # extrinsic
    all_of(ext_cols),
    # labs
    all_of(lab_cols),
    # FISH / cytogenetics
    all_of(fish_cols),
    # flow / adaptive / rapid / PET
    all_of(flow_cols),
    all_of(adaptive_cols),
    all_of(rapid_cols),
    all_of(pet_cols),
    # clinical MRD
    all_of(clinical_mrd_cols),
    # WGS MRD
    all_of(wgs_mrd_cols),
    # coverage (incl FS)
    all_of(coverage_cols),
    # all the rest of the WGS_… columns
    all_of(wgs_cols),
    # anything else not yet selected
    everything()
  )


joined_clean <- joined_clean %>% select(-Date_of_sample_collection)

## get duplicates 
dups <- joined_clean %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()


### There are still dups, remove 

# helper to grab the first non‐NA (or NA if none)
first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) x[[1]] else NA
}

joined_clean <- joined_clean %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    across(everything(), first_non_na),
    .groups = "drop"
  )

## get duplicates 
dups <- joined_clean %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()



#### Now add in the mutation counts and response info 

## Load in the metadata and try to match to it 
cfWGS_metadata <- read.csv("combined_clinical_data_updated_April2025.csv")

# Create a VCF_clean_merge column from the 'Bam' column in cfWGS_metadata
cfWGS_metadata <- cfWGS_metadata %>%
  mutate(VCF_clean_merge = gsub("\\.filter.*", "", Bam))  # Remove '.filter' and anything following it

# Modify select non-overlapping cases - OICR put PG but should be WG to match (just internal naming conventions)
cfWGS_metadata <- cfWGS_metadata %>%
  mutate(VCF_clean_merge = ifelse(
    VCF_clean_merge %in% c(
      "TFRIM4_0031_Bm_P_PG_M4-CA-02-01-O-DNA",
      "TFRIM4_0032_Bm_P_PG_M4-HP-01-01-O-DNA",
      "TFRIM4_0033_Bm_P_PG_M4-MJ-06-01-O-DNA",
      "TFRIM4_0034_Bm_P_PG_M4-VA-02-01-O-DNA",
      "TFRIM4_0035_Bm_P_PG_M4-VA-06-01-O-DNA"
    ), 
    gsub("PG", "WG", VCF_clean_merge),  # Replace WG with PG for these specific values
    VCF_clean_merge  # Keep other values unchanged
  ))


### Now load in the counts files 
# 1. Load the counts files
Blood_muts <- read.delim("Mutation_counts/mutation_counts_table_Blood_no_RSID.txt",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BM_muts    <- read.delim("Mutation_counts/mutation_counts_table_BM_no_RSID.txt",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Clean the 'File' column: remove '.filter' and anything after
Blood_muts$File <- sub("\\.fil.*$", "", Blood_muts$File)
Blood_muts <- Blood_muts[!grepl("\\.PASS", Blood_muts$File), ]
Blood_muts <- Blood_muts[!grepl("\\Bm_P", Blood_muts$File), ]


BM_muts$File    <- sub("\\.fil.*$", "", BM_muts$File)
BM_muts$File <- sub("\\.mutect2.*$", "", BM_muts$File)


# Update specific value in VCF_clean that was truncated
BM_muts <- BM_muts %>%
  mutate(File = ifelse(
    File == "TFRIM4_0189_Bm_P_WG_ZC-02", 
    "TFRIM4_0189_Bm_P_WG_ZC-02-01-O-DNA", 
    File  # Keep other values unchanged
  ))


# join BM counts
cfWGS_with_counts <- cfWGS_metadata %>%
  left_join(
    BM_muts %>% 
      rename(BM_Mutation_Count = Mutation_Count),
    by = c("VCF_clean_merge" = "File")
  ) %>%
  # join Blood counts
  left_join(
    Blood_muts %>% 
      rename(Blood_Mutation_Count = Mutation_Count),
    by = c("VCF_clean_merge" = "File")
  )


## Add this to the main table 
# 1) make sure both date columns are Date
joined_clean <- joined_clean %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

counts_df <- cfWGS_with_counts %>%
  mutate(collection_date = as.Date(Date_of_sample_collection, "%Y-%m-%d")) %>%
  filter(Sample_type != "Blood_Buffy_coat") %>%
  select(Patient, Timepoint, BM_Mutation_Count, Blood_Mutation_Count) %>% unique()

## Check for dups
counts_df %>%
  count(Patient, Timepoint) %>%
  filter(n > 1)

# Consolidate
counts_consolidated <- counts_df %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    BM_Mutation_Count = {
      bm <- BM_Mutation_Count[!is.na(BM_Mutation_Count)]
      if (length(bm)) bm[1] else NA_integer_
    },
    Blood_Mutation_Count = {
      bl <- Blood_Mutation_Count[!is.na(Blood_Mutation_Count)]
      if (length(bl)) bl[1] else NA_integer_
    },
    .groups = "drop"
  )

counts_consolidated <- counts_consolidated %>%
  mutate(Timepoint = if_else(Timepoint == "R-", "R", Timepoint))


# 2) fuzzy‐join within a 14-day window
joined_with_counts <- joined_clean %>%
  left_join(
    counts_consolidated,
    by = c(
      "Patient" = "Patient",
      "Timepoint"        = "Timepoint"
    )
  )  


### Now correct the duplicate rows in joined with counts 
first_non_na <- function(x) {
  x[x %in% c("", "-")] <- NA      # treat blanks and dashes as missing
  vec <- x[!is.na(x)]             # drop all NA
  if (length(vec)) vec[[1]] else NA
}

joined_consolidated <- joined_with_counts %>%
  # normalize your “R-” codes up front
  mutate(
    Timepoint   = str_replace(Timepoint,   "^R-$", "R"),
    Sample_Code = str_replace(Sample_Code, "^R-$", "R")
  ) %>%
  group_by(Patient, Timepoint, Date) %>%   # or just (Patient, Timepoint)
  reframe(
    across(everything(), first_non_na)
  )


## Now have everything need to make plots
# Write to CSV (for Excel/sharing)
write.csv(joined_consolidated, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated.csv", row.names = FALSE)

# Write to RDS (for loading back into R with full structure)
saveRDS(joined_consolidated, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated.rds")

### Do a bit more cleaning for duplicate or missing timepoints

# 1) two lookup tables from your clinical master (de-duplicated if needed)
cc_by_info <- combined_clinical_data_updated %>% 
  select(Patient, timepoint_info, 
         TP_info = Timepoint, 
         Date_info = Date_of_sample_collection)

cc_by_date <- combined_clinical_data_updated %>% 
  select(Patient, Date_of_sample_collection, 
         TP_date = Timepoint, 
         Info_date = timepoint_info)

cc_by_date <- cc_by_date %>%
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection, origin = "1970-01-01")
  )

# 2) left-join both onto joined_consolidated
joined_filled2 <- joined_consolidated %>% 
  # join by info
  left_join(cc_by_info, 
            by = c("Patient","timepoint_info")) %>% 
  # join by date
  left_join(cc_by_date, 
            by = c("Patient","Date" = "Date_of_sample_collection")) %>% 
  
  mutate(
    # only replace missing Timepoint
    Timepoint = if_else(
      is.na(Timepoint),
      coalesce(TP_info, TP_date),
      Timepoint
    ),
    
    # only replace missing timepoint_info
    timepoint_info = if_else(
      is.na(timepoint_info),
      coalesce(Info_date, timepoint_info),
      timepoint_info
    ),
    
    # only replace missing Date
    Date = if_else(
      is.na(Date),
      coalesce(Date_info, Date),
      Date
    )
  ) %>%
  # drop the helper columns
  select(-TP_info, -Date_info, -TP_date, -Info_date)

# 3) if you want to toss any that still lack a Timepoint
joined_filled2 <- joined_filled2 %>% 
  filter(!is.na(Timepoint))

# 4) finally de-duplicate exactly as before
joined_clean2 <- joined_filled2 %>%
  group_by(Patient, Timepoint, Date) %>%
  reframe(
    across(everything(), first_non_na)
  )


## Now have everything need to make plots
# Write to CSV (for Excel/sharing)
write.csv(joined_clean2, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated2.csv", row.names = FALSE)

# Write to RDS (for loading back into R with full structure)
saveRDS(joined_clean2, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated2.rds")


## Get unique pts 
patients_with_zscores <- joined_clean2 %>%
  filter(!is.na(zscore_BM) | !is.na(zscore_blood)) %>%
  distinct(Patient)


### Add missing tumor fraction to table 
missing_tumor_fraction <- read.csv("Missing_tumor_fraction.csv")

# 3) join back onto dat and fill your two WGS‐tumor‐fraction columns
joined_clean2 <- joined_clean2 %>%
  # join in the new BM / blood TF values
  left_join(
    missing_tumor_fraction %>%
      select(Patient, Timepoint,
             tf_BM_new = WGS_Tumor_Fraction_BM_cells,
             tf_blood_new = WGS_Tumor_Fraction_Blood_plasma_cfDNA),
    by = c("Patient", "Timepoint")
  ) %>%
  # overwrite only where the originals were NA
  mutate(
    WGS_Tumor_Fraction_BM_cells = coalesce(WGS_Tumor_Fraction_BM_cells, tf_BM_new),
    WGS_Tumor_Fraction_Blood_plasma_cfDNA = coalesce(WGS_Tumor_Fraction_Blood_plasma_cfDNA, tf_blood_new)
  ) %>%
  # drop the helper cols
  select(-tf_BM_new, -tf_blood_new)


### One more duplicate check 
# 1) Identify which (Patient, timepoint_info) groups have a “-NA” sample_code
dup_keys <- joined_clean2 %>%
  filter(str_detect(Sample_Code, "-NA$")) %>%
  distinct(Patient, timepoint_info)

# 2) For those groups, coalesce the two rows into one—
#    preferring any value coming from the non-“-NA” Sample_Code row
collapsed <- joined_clean2 %>%
  semi_join(dup_keys, by = c("Patient","timepoint_info")) %>%
  group_by(Patient, timepoint_info) %>%
  summarise(across(everything(), ~ {
    col <- cur_column()
    vals <- .x
    
    if (col == "Sample_Code") {
      # drop the “-NA” code, then take the first remaining
      clean <- vals[!str_detect(vals, "-NA$")]
      if (length(clean)>0) clean[1] else NA_character_
    } else {
      # take the first non-NA value
      non_na <- vals[!is.na(vals)]
      if (length(non_na)>0) non_na[1] else NA
    }
  }), .groups="drop")

# 3) Everything *not* in those duplicate groups just stays as-is
unchanged <- joined_clean2 %>%
  anti_join(dup_keys, by = c("Patient","timepoint_info"))

# 4) Re-combine
joined_clean2_fixed <- bind_rows(unchanged, collapsed) %>%
  arrange(Patient, Date)   # or however you like to sort

joined_clean2 <- joined_clean2_fixed %>%
  mutate(
    Sample_Code = if_else(
      is.na(Sample_Code),
      paste0(Patient, "-", Timepoint),
      Sample_Code
    )
  )


## Another check
# helper to treat "" or "-" as NA, then pick first non-NA
first_non_na <- function(x) {
  x[x %in% c("", "-")] <- NA
  vals <- x[!is.na(x)]
  if (length(vals)) vals[1] else NA
}

joined_clean3 <- joined_clean2 %>%
  # 1) fix any Sample_Code ending in "R-"
  mutate(
    Sample_Code = str_replace(Sample_Code, "R-$", "R")
  ) %>%
  
  # 2) group on Patient, Date, Timepoint and collapse
  group_by(Patient, Date, Timepoint) %>%
  reframe(
    across(everything(), first_non_na)
  )

## Clean dates 
df <- joined_clean3

# 1) For each Patient/Timepoint/Sample_Code, compute exactly one fill_date:
fill_dates <- df %>%
  group_by(Patient, Timepoint, Sample_Code) %>%
  summarise(
    fill_date = {
      ud <- unique(Date[!is.na(Date)])
      if (length(ud) == 1) {
        # exactly one date → keep it
        ud
      } else {
        # zero dates or multiple conflicting → NA
        as.Date(NA)
      }
    },
    .groups = "drop"
  )

# 2) Join back and coalesce only the missing Dates
filled_df <- df %>%
  left_join(fill_dates, by = c("Patient","Timepoint","Sample_Code")) %>%
  mutate(
    Date = coalesce(Date, fill_date)
  ) %>%
  select(-fill_date)

# 3) (Optional) sanity‐check that you didn’t miss any fillable gaps
residual <- filled_df %>%
  group_by(Patient, Timepoint, Sample_Code) %>%
  filter(any(is.na(Date)) & any(!is.na(Date))) %>%
  distinct(Patient, Timepoint, Sample_Code)

if (nrow(residual)) {
  warning("These groups still have mixed NA/non-NA dates:\n",
          paste0(sprintf("%s / %s / %s", residual$Patient, residual$Timepoint, residual$Sample_Code),
                 collapse = "\n"))
} else {
  message("All missing dates successfully filled where possible.")
}

filled_df <- filled_df %>% unique()

filled_df %>%
  # for safety, make sure Date is Date class
  mutate(Date = as.Date(Date)) %>%
  # group by Patient + Timepoint
  group_by(Patient, Timepoint) %>%
  # collect unique dates and count them
  summarise(
    n_dates = n_distinct(Date),
    dates   = paste(sort(unique(as.character(Date))), collapse = ", "),
    .groups = "drop"
  ) %>%
  # keep only those with more than one distinct date
  filter(n_dates > 1)

## One last pass to coalesce
filled_df <- filled_df %>%
  #group on Patient, Date and collapse
  group_by(Patient, Date) %>%
  reframe(
    across(everything(), first_non_na)
  )

## Now have everything need to make plots

# Write to CSV (for Excel/sharing)
write.csv(filled_df, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.csv", row.names = FALSE)

# Write to RDS (for loading back into R with full structure)
saveRDS(filled_df, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.rds")


## Quick check 
filled_df %>%
  filter(is.na(FS) & !is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA))

filled_df %>%
  filter(!is.na(FS) & is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA))



### Check labs are ok 

# 1) Define your “lab” columns
lab_vars <- c(
  "Albumin", "B2_micro", "Calcium", "Creatinine", "Hemoglobin", "LDH",
  "IgA", "IgG", "IgM", "Kappa", "Kappa_Lambda_Ratio", "Lambda", "M_Protein"
)

# Keep only those that actually exist
lab_vars <- intersect(lab_vars, names(filled_df))

# 2) Compute min / Q1 / median / mean / Q3 / max for each
lab_summary <- filled_df %>%
  summarise(across(
    all_of(lab_vars),
    list(
      min    = ~ min(.x, na.rm=TRUE),
      Q1     = ~ quantile(.x, 0.25, na.rm=TRUE),
      median = ~ median(.x, na.rm=TRUE),
      mean   = ~ mean(.x, na.rm=TRUE),
      Q3     = ~ quantile(.x, 0.75, na.rm=TRUE),
      max    = ~ max(.x, na.rm=TRUE)
    ),
    .names = "{col}_{fn}"
  ))

print(lab_summary, width = Inf)

# 3) (Optional) reshape for easier inspection
tmp <- lab_summary %>%
  pivot_longer(everything(),
               names_to = c("Lab", "Stat"),
               names_sep = "_") %>%
  pivot_wider(names_from = Stat, values_from = value) %>%
  arrange(Lab) %>%
  print(n = Inf)

# 4) Quick pairwise scatter/diag plot to catch weird units
#    This will be heavy if you have many samples—consider sampling 100 pts if slow.
library(GGally)
filled_df %>%
  select(all_of(lab_vars)) %>%
  # drop raws with any NAs across labs
  drop_na() %>%
  # pick a sample if too big
  slice_sample(n = min(200, dplyr::n())) %>%
  ggpairs(title = "Lab Variables Pairwise Check")


filled_df %>% 
  summarise(
    Ca_min  = min(Calcium,    na.rm=TRUE),
    Ca_max  = max(Calcium,    na.rm=TRUE),
    B2_min  = min(B2_micro,   na.rm=TRUE),
    B2_max  = max(B2_micro,   na.rm=TRUE),
    LDH_min = min(LDH,        na.rm=TRUE),
    LDH_max = max(LDH,        na.rm=TRUE)
  )
