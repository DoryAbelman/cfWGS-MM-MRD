# =============================================================================
# Script: 2_1_Part2_Cohort_Swim_Plot.R
# Purpose: Build cohort swim plot of treatment timelines
# =============================================================================
#### Now make swim plot 

library(tidyverse)
library(readxl)
library(lubridate)
library(patchwork)

## --------------------------
## 1. Load all three sources
## --------------------------

# 1a. “tidy_treatments.csv” you already have
tt_raw <- read_csv("Clinical data/SPORE/tidy_treatments.csv",
                   col_types = cols(.default = "c"))   # keep all as character; we’ll parse later



# 1b. M4 cohort chemotherapy Excel
m4_raw <- read_excel("M4_CMRG_Data/M4_COHORT_CHEMOTHERAPY.xlsx",
                     col_types = "text")                # read everything as character

# 1c. IMMAGINE chemotherapy CSV
imm_raw <- read_csv("Clinical data/IMMAGINE/Cleaned_IMMAGINE_chemotherapy.csv",
                    col_types = cols(.default = "c"))

## -------------------------------------------------
## 2. Tidy each data set into (patient, event, date)
## -------------------------------------------------

## First SPORE 

# ─────────────────────────────────────────────────────────────────────────────
# A. SPORE “tidy_treatments.csv”
# ─────────────────────────────────────────────────────────────────────────────
## 2a. tidy_treatments  (already long, just rename / parse), includes stransplant
tt <- tt_raw %>%
  transmute(
    patient = patient,                        # your CSV’s patient column
    event   = str_to_sentence(line),          # “line” → e.g. “Diagnosis”, “Transplant”
    start   = ymd(start_date),                # parse yyyy-mm-dd
    end     = NA,                          # single-day events, set as same
    details = regimen                         # whatever you’d like to show here
  )


# ─────────────────────────────────────────────────────────────────────────────
# B. M4 cohort 
# ─────────────────────────────────────────────────────────────────────────────
### Next M4
## 2b. M4 chemotherapy rows  → “Chemotherapy” events
m4_chemo <- m4_raw %>%
  mutate(
    # convert Excel serial → R Date
    start = as.Date(as.numeric(START_DATE), origin = "1899-12-30"),
    end   = as.Date(as.numeric(END_DATE),   origin = "1899-12-30"),
    details = REGIMEN_NAME
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Chemotherapy",
    start,
    end,
    details
  )

m4_chemo <- m4_chemo %>% filter(!is.na(start))
write.csv(m4_chemo %>% filter(patient %in% cohort_df$Patient) %>% filter(is.na(end)), "m4_chemo.csv")

### Add transplant info for M4 
# Read the raw sheet
m4_trans_raw <- read_excel(
  "M4_CMRG_Data/M4_COHORT_STEM_CELL_TRANSPLANT.xlsx",
  col_types = "text"
)

# Robust date parser with NA‐guards
parse_m4_date <- function(x) {
  n <- length(x)
  out <- rep(as.Date(NA), n)
  
  # only non-NA entries
  not_na <- !is.na(x)
  
  # 2a. Excel serials: purely digits
  is_serial <- not_na & str_detect(x, "^[0-9]+$")
  out[is_serial] <- as.Date(as.numeric(x[is_serial]), origin = "1899-12-30")
  
  # 2b. Remaining non-blank, non-serial text
  is_text <- not_na & !is_serial & x != ""
  out[is_text] <- ymd(x[is_text])
  
  out
}

# 2c. Transplant events
m4_transplant_events <- m4_trans_raw %>%
  mutate(
    tx_date = parse_m4_date(TRANSPLANT_DATE),
    details = str_c(
      PROCEDURE_TYPE_TEXT,
      TRANSPLANT_TYPE,
      if_else(is.na(INJECTED_CD34) | INJECTED_CD34 == "",
              "",
              str_c("CD34:", INJECTED_CD34)),
      sep = " | "
    )
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Transplant",
    start   = tx_date,
    end     = tx_date,
    details
  )

# 2d. Check for any remaining failures in parsing
bad_dates <- m4_trans_raw %>%
  filter(!is.na(TRANSPLANT_DATE) & TRANSPLANT_DATE != "" &
           is.na(parse_m4_date(TRANSPLANT_DATE))) %>%
  pull(TRANSPLANT_DATE) %>%
  unique()

if (length(bad_dates)) {
  message("These raw TRANSPLANT_DATE values still failed to parse:\n",
          paste(bad_dates, collapse = ", "))
}



## 2e. Immune‐response event
m4_immun_resp <- m4_trans_raw %>%
  mutate(IMMUN_RESPONSE_DATE = as.numeric(IMMUN_RESPONSE_DATE)) %>%
  # keep only rows with a non‐blank IMMUN_RESPONSE_DATE
  filter(!is.na(IMMUN_RESPONSE_DATE) & IMMUN_RESPONSE_DATE != "") %>%
  mutate(
    # parse with the same helper
    resp_date = parse_m4_date(IMMUN_RESPONSE_DATE)
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Immune response",
    start   = resp_date,
    end     = resp_date,
    details = IMMUN_RESPONSE_TYPE
  )

## 2f. Best‐response event
m4_best_resp <- m4_trans_raw %>%
  mutate(BEST_RESPONSE_DATE = as.numeric(BEST_RESPONSE_DATE)) %>%
  filter(!is.na(BEST_RESPONSE_DATE) & BEST_RESPONSE_DATE != "") %>%
  mutate(
    best_date = parse_m4_date(BEST_RESPONSE_DATE)
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Best response",
    start   = best_date,
    end     = best_date,
    details = BEST_RESPONSE
  )




# ─────────────────────────────────────────────────────────────────────────────
#  C. IMMAGINE 
# ─────────────────────────────────────────────────────────────────────────────
### Lastly IMMAGINE
## 2d. IMMAGINE rows  (already long but rename fields)
imm_clean <- imm_raw %>%
  rename(all = `Patient,Event,Date,Details`) %>%
  separate(
    col  = all,
    into = c("patient", "event", "date", "details"),
    sep  = ",",
    fill = "right"      # in case some rows end with a comma
  )

## Fix dates 
imm_clean <- imm_clean %>%
  # Capitalize event names if you like
  mutate(event = str_to_sentence(event)) %>%
  
  # Count dashes, pad to first-of-month or first-of-year
  mutate(
    n_dash   = str_count(date, fixed("-")),
    date_full = case_when(
      n_dash == 2 ~ date,                   # “YYYY-MM-DD”
      n_dash == 1 ~ paste0(date, "-01"),    # “YYYY-MM”    → “YYYY-MM-01”
      n_dash == 0 ~ paste0(date, "-01-01"), # “YYYY”       → “YYYY-01-01”
      TRUE        ~ date
    )
  ) %>%
  
  # Parse into real Date, set end = start
  mutate(
    start = ymd(date_full),
    end   = NA
  ) %>%
  
  # Drop helper columns
  select(patient, event, start, end, details)




# ─────────────────────────────────────────────────────────────────────────────
# D. Relapse, Baseline & Last follow‐up events
# ─────────────────────────────────────────────────────────────────────────────
#### Add diagnosis dates to this and the censor dates
# a) relapse dates (one row per patient, many relapse‐cols)
relapse_dates_full <- read_csv(
  "Relapse dates cfWGS updated.csv",
  col_types = cols(.default = "c")
)

# b) sample‐collection / censor dates per patient
censor_tbl <- readRDS("Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds")
censor_tbl$Baseline_Date <- censor_tbl$baseline_date # for consistency

# 2) Build a “Relapse” events table ----------------------------------------
relapse_events <- relapse_dates_full %>%
  # if your patient ID column is called something else, rename it:
  rename(patient = Patient) %>%  
  # pivot all the relapse columns (e.g. Relapse1, Relapse2, …) into long form:
  pivot_longer(
    cols      = -patient,
    names_to  = "which_relapse",
    values_to = "date"
  ) %>%
  filter(!is.na(date) & date != "") %>%            # drop missing
  mutate(
    start   = ymd(date),                           # parse to Date
    end     = start,
    event   = "Relapse",
    details = which_relapse                        # e.g. “Relapse1”
  ) %>%
  select(patient, event, start, end, details)


# 3) Build a baseline and date of last followup events table -----------------------------
# 1. Baseline events from the Baseline_Date column
baseline_events <- censor_tbl %>%
  transmute(
    patient = Patient,
    event   = "Baseline",
    # Baseline_Date is POSIXct; convert to Date
    start   = as_date(Baseline_Date),
    end     = start,
    details = NA_character_
  )

# 2. Last follow-up events from the Censor_date column
followup_events <- censor_tbl %>%
  transmute(
    patient = Patient,
    event   = "Last follow-up",
    start   = censor_date,
    end     = start,
    details = paste0("Relapsed: ", relapsed)
  )




# ─────────────────────────────────────────────────────────────────────────────
# E. Sample collections (BM vs cfDNA)
# ─────────────────────────────────────────────────────────────────────────────
### Now add all the sample collection dates 
combined_clinical_data_updated <- read.csv("combined_clinical_data_updated_April2025.csv")

# 1. BM sample collection events
bm_events <- combined_clinical_data_updated %>%
  filter(Sample_type == "BM_cells") %>%
  transmute(
    patient = Patient,
    event   = "BM sample collection",
    start   = Date_of_sample_collection,
    end     = start,
    details = Sample_ID
  )

# 2. cfDNA (plasma) sample collection events
cfDNA_events <- combined_clinical_data_updated %>%
  # assuming plasma‐derived cfDNA are labeled "Blood_plasma"
  filter(str_detect(Sample_type, regex("plasma", ignore_case = TRUE))) %>%
  transmute(
    patient = Patient,
    event   = "cfDNA sample collection",
    start   = Date_of_sample_collection,
    end     = start,
    details = Sample_ID
  )


# ─────────────────────────────────────────────────────────────────────────────
# F. MRD test dates
# ─────────────────────────────────────────────────────────────────────────────
dat <- read.csv("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.csv")


mfc_events <- dat %>%
  filter(!is.na(Flow_Binary)) %>%
  mutate(
    start   = as_date(Date),
    end     = start,
    details = paste0("Flow result: ", Flow_Binary)
  ) %>%
  transmute(
    patient = Patient,
    event   = "MRD (MFC)",
    start,
    end,
    details
  )

# 2. clonoSEQ (adaptive) MRD test events
clonoseq_events <- dat %>%
  filter(!is.na(Adaptive_Binary)) %>%
  mutate(
    start   = as_date(Date),
    end     = start,
    details = paste0("Adaptive result: ", Adaptive_Binary)
  ) %>%
  transmute(
    patient = Patient,
    event   = "MRD (clonoSEQ)",
    start,
    end,
    details
  )

# Get everythong as date 
# e.g. if bm_events has character dates in “YYYY‑MM‑DD” form:
bm_events <- bm_events %>% 
  mutate(
    start = as_date(start),      # parse “2020‑07‑23” into Date
    end   = as_date(end)         # if end is also character
  )

# repeat for any other dfs with char dates...
cfDNA_events   <- cfDNA_events   %>% mutate(start = as_date(start), end = as_date(end))
followup_events<- followup_events%>% mutate(start = as_date(start), end = as_date(end))


# ─────────────────────────────────────────────────────────────────────────────
# 3. Combine Everything & Check
# ─────────────────────────────────────────────────────────────────────────────
all_events <- bind_rows(
  tt,
  m4_chemo,
  m4_transplant_events,
  m4_immun_resp,
  m4_best_resp,
  m4_progression,
  imm_clean,
  relapse_events,
  baseline_events,
  followup_events,
  bm_events,
  cfDNA_events,
  clonoseq_events, 
  mfc_events
) %>%
  arrange(patient, start)

all_events <- all_events %>% filter(!is.na(start))

# Quick parse check
bad <- all_events %>%
  filter(is.na(start)) %>%
  distinct(patient, event, start) 

if (nrow(bad)) {
  message("Warning: the following rows have NA start dates:\n")
  print(bad)
} else {
  message("All events parsed successfully!")
}


### Now edit 
all_events <- all_events %>%
  mutate(
    event = case_when(
      # 1) For any SPORE patient whose event name has a digit → Chemotherapy
      str_detect(patient, "^SPORE") & str_detect(event, "\\d+") ~ "Chemotherapy",
      # 2) For any SPORE patient whose details mention ASCT → Transplant
      str_detect(patient, "^SPORE") & 
        str_detect(details, regex("ASCT", ignore_case = TRUE)) ~ "Transplant",
      # 3) Globally, any “Induction” → Chemotherapy
      event == "Induction" ~ "Chemotherapy",
      # 4) Otherwise, keep whatever was there
      TRUE ~ event
    )
  )

all_events <- all_events %>%
  mutate(
    # Collapse both “Best response” and “Immune response” into “Response”
    event = case_when(
      event %in% c("Best response", "Immune response") ~ "Response",
      TRUE                                             ~ event
    )
  ) %>%
  unique() %>%
  # Drop any stray Excel‐origin rows that parsed to 1899-12-31
  filter(start != as.Date("1899-12-31"))

cohort_df <- readRDS("cohort_assignment_table_updated.rds")

## Filter to in cohort df 
all_events <- all_events %>% filter(patient %in% cohort_df$Patient)

## add end to be start if unsure 
all_events <- all_events %>%
  mutate(
    end = if_else(is.na(end), start, end)
  )

all_events <- all_events %>% select(-Patient, -Progression_date)


### Add new info from Sarah and Esther 
M4_new <- read_excel("Clinical data/M4/Updated dates from Sarah for swim plot - DA edited.xlsx")
IMG_new <- read_excel("Clinical data/IMMAGINE/Updated_data_esther.xlsx")

  # 1) Parse the M4_new dates into Date class
  M4_new2 <- M4_new %>%
  mutate(
    start = dmy(start),
    end   = dmy(end)
  )

# 2) Drop from all_events any rows that share patient+event+start with M4_new2
all_events_clean <- all_events %>%
  anti_join(
    M4_new2 %>% select(patient, event, start),
    by = c("patient", "event", "start")
  )

# 3) Bind them back together and (optionally) re‑order
all_events_updated <- bind_rows(all_events_clean, M4_new2) %>%
  arrange(patient, start)


# 1) Convert IMG_new’s POSIXct columns to Date
IMG_new2 <- IMG_new %>%
  mutate(
    start = as_date(start),
    end   = as_date(end)
  )

# 2) Drop all “Transplant” rows for IMG‑181 and IMG‑098 from your master table
all_events_clean <- all_events_updated %>%
  filter(!(patient %in% c("IMG-181","IMG-098") & event == "Transplant"))

# 3) (Optional) also guard against exact dupes on patient+event+start
all_events_clean <- all_events_clean %>%
  anti_join(
    IMG_new2 %>% filter(event=="Transplant") %>% select(patient, event, start),
    by = c("patient","event","start")
  )

# 4) Finally bind your cleaned master with the new rows
all_events <- bind_rows(all_events_clean, IMG_new2) %>%
  arrange(patient, start)

## Ensure ASCT treated as transplant 
all_events <- all_events %>%
  mutate(
    event = if_else(
      details == "ASCT",    # when TRUE…
      "Transplant",         # …use this
      event,                # when FALSE…
      missing = event       # …and when NA, also use this
    )
  )

## Now if ongoing put to the latest timepoint have on patient for chemotherapy or progression 
# define which events to fall back on if there is no next Chemo/Transplant
fallback_events <- c("Relapse", "Last follow-up")

# 2) For each patient, stretch only Chemo rows whose end == start
all_events_updated <- all_events %>%
  group_by(patient) %>%
  arrange(start) %>%
  group_modify(~ {
    df <- .
    
    for (i in seq_len(nrow(df))) {
      if (
        !is.na(df$event[i]) &&
        df$event[i] == "Chemotherapy" &&
        !is.na(df$start[i]) &&
        !is.na(df$end[i]) &&
        df$end[i] == df$start[i]
      ) {
        # 1) Next Chemo or Transplant
        next_chemo_tx <- df$start[
          df$event %in% c("Chemotherapy", "Transplant") &
            df$start > df$start[i]
        ]
        
        if (length(next_chemo_tx) > 0) {
          df$end[i] <- min(next_chemo_tx)
        } else {
          # 2) Fallback: earliest Relapse/Last followup AFTER this chemo start
          future_fallback <- df$start[
            df$event %in% fallback_events &
              df$start > df$start[i]
          ]
          if (length(future_fallback) > 0) {
            df$end[i] <- min(future_fallback)
          }
          # else leave df$end[i] as-is (or set to NA if you prefer)
        }
      }
    }
    
    df
  }) %>%
  ungroup()



all_events <- all_events_updated
# Export all_events to CSV
write_csv(all_events %>% select(-details_2), "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated.csv")

# Export all_events to RDS
saveRDS(all_events, "Final Tables and Figures/all_events_for_swim_plot_combined_updated2.rds")

#all_events <- readRDS("Final Tables and Figures/all_events_for_swim_plot_combined_updated2.rds")
### Export clean version with correct id 

### this is a supplementary table used in manuscript 

# 1. Load the ID map (must have columns Patient, New_ID)
id_map <- readRDS("id_map.rds") %>% distinct(Patient, New_ID)

# 2. Merge all_events with id_map and replace Patient with New_ID
all_events_updated <- all_events %>%
  rename(Patient = patient) %>%
  left_join(id_map, by = "Patient") %>%
  mutate(Patient = coalesce(New_ID, Patient),  # use New_ID when available
  details = if_else(
    event %in% c("BM sample collection", "cfDNA sample collection"),
    NA_character_,
    details
  )) %>%
  select(-New_ID, -details_2)                      # drop temp + details_2 col

# 3. Write to CSV
write_csv(all_events_updated,
          "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated_with_new_patient_ID.csv")


## Export for Esteban to get more info 
spore_chemo_events <- all_events %>%
  filter(event == "Chemotherapy", grepl("^SPORE", patient)) %>% 
  mutate(end = NA) %>% 
  filter(details != "ASCT")

# Export to CSV
write_csv(spore_chemo_events, "spore_chemo_events.csv")

## Export for Esteban to get more info 
spore_all_events <- all_events %>%
  filter(grepl("^SPORE", patient)) %>% 
  mutate(end = NA) %>% 
  filter(details != "ASCT")

# Export to CSV
write_csv(spore_all_events, "spore_all_events.csv")


 #### Now assemble plot 
# ──────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ──────────────────────────────────────────────────────────────────────
events   <- read_csv("Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated.csv",
                     col_types = cols(
                       patient = col_character(),
                       event   = col_character(),
                       start   = col_date(),
                       end     = col_date(),
                       details = col_character()
                     ))

events <- events %>%
  # recode any “Relapse” rows into “Progression”
  mutate(
    event = if_else(event == "Relapse", "Progression", event)
  )

cohort_df <- readRDS("cohort_assignment_table_updated.rds")

# ──────────────────────────────────────────────────────────────────────
# 2. MERGE COHORT INFO  &  KEEP PATIENTS WITH A BASELINE
# ──────────────────────────────────────────────────────────────────────
events <- events %>%
  left_join(cohort_df,  by = c("patient" = "Patient")) %>%
  mutate(cohort = if_else(Cohort == "Frontline", "Front-line cohort",
                          "Non-front-line cohort"))


# (a) Front‐line patients use the "Baseline" event date
baseline_front <- events %>%
  filter(cohort == "Front-line cohort", event == "Baseline") %>%
  select(patient, baseline_date = start)

# (b) Non‐front‐line patients use the first BM or blood draw
#     where Timepoint is “Diagnosis” or “Baseline” in your sample table
non_ids <- cohort_df %>%
  filter(Cohort != "Frontline") %>%
  pull(Patient)

# then compute baseline_non only for those non-frontline patients
baseline_non <- combined_clinical_data_updated %>%
  filter(
    Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA"),
    timepoint_info %in% c("Diagnosis", "Baseline")
  ) %>%
  group_by(Patient) %>%
  summarize(
    baseline_date = min(Date_of_sample_collection, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(patient = Patient) %>%
  filter(patient %in% non_ids)

# bind the two sets of baseline dates
# parse baseline_date in baseline_non to Date
baseline_non <- baseline_non %>%
  mutate(baseline_date = ymd(baseline_date))

# (Optional) if baseline_front ended up as character too, parse it likewise:
baseline_front <- baseline_front %>%
  mutate(baseline_date = as_date(baseline_date))

# now bind
baseline_tbl <- bind_rows(baseline_front, baseline_non)

# re‐join and drop any patients w/o a computed baseline
events <- events %>%
  left_join(baseline_tbl, by = "patient") %>%
  filter(!is.na(baseline_date))

# ──────────────────────────────────────────────────────────────────────
# 3. DAYS-FROM-BASELINE & INTERVAL FLAG
# ──────────────────────────────────────────────────────────────────────
events <- events %>%
  mutate(
    start_day   = as.numeric(start   - baseline_date),
    end_day     = as.numeric(end     - baseline_date),
    is_interval = start_day != end_day
  )

## edit typo
events <- events %>%
  # if the details column is exactly "ASCT" (or contains ASCT), recode event
  mutate(
    event = if_else(
      str_detect(details, regex("^ASCT$", ignore_case = TRUE)),
      "Transplant",
      event
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4. EXTEND ONE-DAY CHEMO BARS → 30 d BEFORE NEXT CHEMO START
# ─────────────────────────────────────────────────────────────────────────────
#chemo_adj <- events %>%
#  filter(event == "Chemotherapy") %>%
#  arrange(patient, start_day) %>%
#  group_by(patient) %>%
#  mutate(next_start = lead(start_day)) %>%
#  mutate(end_day = if_else(end_day == start_day & !is.na(next_start) &
#                             (next_start - 30 > start_day),
#                           next_start - 30, end_day)) %>%
#  select(patient, start_day, end_day)

#events <- events %>%
#  left_join(chemo_adj, by = c("patient", "start_day"),
#            suffix = c("", "_new")) %>%
#  mutate(
#    end_day = coalesce(end_day_new, end_day),
#    end     = baseline_date + days(end_day),
#    is_interval = if_else(event == "Chemotherapy" & end_day > start_day,
#                          TRUE, is_interval)
#  ) %>%
#  select(-ends_with("_new"))

# ─────────────────────────────────────────────────────────────────────────────
# 5. REGIMEN → COLOUR GROUP
# ─────────────────────────────────────────────────────────────────────────────

## See events
all_chemo <- events %>%
  filter(event == "Chemotherapy") %>%
  distinct(details) %>%
  arrange(details)

events <- events %>%
  mutate(
    chemo_group = case_when(
      event != "Chemotherapy"                                 ~ NA_character_,
      str_detect(details, regex("CY[Bb]OR.?[DP]", ignore_case = TRUE))  
      ~ "CyBorD", 
      str_detect(details, regex("\\bV?RD\\b|Lenalidomide|Rev", TRUE))
      ~ "R/VRD ± Len",
      str_detect(details, regex("Dara", TRUE))                ~ "Dara-based",
      str_detect(details, regex("Carfilzomib|\\bKD\\b|CAR", TRUE))
      ~ "Carfilzomib-based",
      str_detect(details, regex("Ixazomib|\\bIxa\\b", TRUE))  ~ "Ixazomib-based",
      str_detect(details, regex("Elranatamab", TRUE))         ~ "Elranatamab",
      str_detect(details, regex("Iberdomide", TRUE))          ~ "Iberdomide",
      TRUE                                                    ~ "Other"
    )
  )

## Updated 
events <- events %>%
  mutate(
    chemo_group = case_when(
      event != "Chemotherapy"                                            ~ NA_character_,
      # CyBorD (any D or P variant)
      str_detect(details, regex("CY[Bb]OR.?[DP]", ignore_case=TRUE))    ~ "CyBorD",
      # VRD / RVD / Lenalidomide / Revlimid variants
      str_detect(details, regex("\\bV?RD\\b|Lenalidomide|Rev", TRUE))   ~ "Lenalidomide-based",
      # Daratumumab‐based combos (even DaraPom, DaraRVD, etc.)
      str_detect(details, regex("Dara|Daratumumab", TRUE))              ~ "Dara-based",
      # Pomalidomide‐based (PomDex, KPomD, CAR POM, ELO POM, DAR POM, PomDex)
      str_detect(details, regex("Pom(Dex|D)|POM|KPomD|CAR POM|ELO POM|DAR POM|PomDex", TRUE))
      ~ "Pomalidomide-based",
      # Carfilzomib combinations (Carfil-, KD-, but avoid “DAR CAR” which is caught above)
      str_detect(details, regex("Carfilozomib|\\bKD\\b|CAR(?! *REV)", TRUE))
      ~ "Carfilzomib-based",
      # Ixazomib combos
      str_detect(details, regex("Ixa|Ixazomib", TRUE))                  ~ "Ixazomib-based",
      # Elranatamab (including the /daratumumab trial)
      str_detect(details, regex("Elranatamab", TRUE))                   ~ "Elranatamab",
      # Iberdomide combos
      str_detect(details, regex("Iberdomide", TRUE))                    ~ "Iberdomide",
      # Cyclophosphamide + dexamethasone (“CYCLO DEX” or “CYCLONE”)
      str_detect(details, regex("CYCLO DEX|CYCLONE", TRUE))             ~ "Cyclophosphamide-based",
      # Early-phase/trial drugs
 #     str_detect(details, regex("MEDI|TAK-|VSV", TRUE))                 ~ "Clinical trial",
      # Catch-all for everything else
      TRUE                                                              ~ "Other"
    )
  )

chemo_cols <- c(
  "CyBorD"              = "#4477AA",
  "R/VRD ± Len"         = "#CC6677",
  "Dara-based"          = "#228833",
  "Carfilzomib-based"   = "#AA3377",
  "Ixazomib-based"      = "#66C2A5",
  "Elranatamab"         = "#EE7733",
  "Iberdomide"          = "#994F00",
  "Other"               = "#999999"
)


### Updated
chemo_cols <- c(
  "CyBorD"                    = "#0072B2",  # blue
  "Lenalidomide-based"        = "#CC6677",  # pink 
  "Dara-based"                = "#009E73",  # green
  "Pomalidomide-based"        = "#D55E00",  # vermillon
  "Carfilzomib-based"         = "#E69F00",  # orange
  "Ixazomib-based"            = "#56B4E9",  # sky blue
  "Elranatamab"               = "#F0E442",  # yellow
  "Iberdomide"                = "#800080",  # purple
  "Cyclophosphamide-based"    = "#994F00",  # brown-orange
#  "Clinical trial"            = "#F7C6C7",  # light pink
  "Other"                     = "#999999"   # grey
)

# ─────────────────────────────────────────────────────────────────────────────
# 6. SHAPE MAP FOR ALL POINT-EVENTS
# ─────────────────────────────────────────────────────────────────────────────
shape_map <- c(
  "BM sample collection"   = 3,   # +
  "cfDNA sample collection"= 4,   # x
  "MRD (MFC)"              = 1,   # circle
  "MRD (clonoSEQ)"         = 18,  # diamond
  "Transplant"             = 7,  # square with x
  "Progression"            = 16   # filled circle
)

# Treat transplants & progression as points
events <- events %>%
  mutate(is_interval = if_else(event %in% c("Transplant", "Progression"),
                               FALSE, is_interval))

# ─────────────────────────────────────────────────────────────────────────────
# 7. PATIENT ORDER (earliest baseline → top)
# ─────────────────────────────────────────────────────────────────────────────
## Change to based on tumor fraction change from ichorCNA
dat_tf <- dat %>%
  mutate(Date = as_date(Date)) %>%
  transmute(
    Patient,
    Date,
    timepoint_info,
    TF = suppressWarnings(as.numeric(WGS_Tumor_Fraction_Blood_plasma_cfDNA))
  )

# 1) First non-NA baseline TF where timepoint_info is Baseline/Diagnosis
baseline_df <- dat_tf %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis"),
         !is.na(TF)) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            baseline_TF   = TF,
            baseline_date = Date)

# 2) First non-NA TF strictly AFTER that baseline date
followup_df <- dat_tf %>%
  inner_join(baseline_df, by = "Patient") %>%
  filter(!is.na(TF),
         Date > baseline_date) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            first_post_TF   = TF,
            first_post_date = Date)

# 3) Combine + percent change (baseline → first post)
tf_change <- baseline_df %>%
  left_join(followup_df, by = "Patient") %>%
  mutate(
    pct_change = if_else(!is.na(baseline_TF) & !is.na(first_post_TF) & baseline_TF > 0,
                         100 * (first_post_TF - baseline_TF) / baseline_TF,
                         NA_real_)
  )

## Add cohort 
tf_change <- tf_change %>%
  left_join(events %>% distinct(patient, cohort) %>% rename(Patient = patient),
            by = "Patient")

## Adjust for NAs 
tf_change <- tf_change %>%
  mutate(
    pct_for_plot = case_when(
      !is.na(pct_change) ~ pct_change,
      is.na(pct_change) & baseline_TF == 0 & first_post_TF > 0 ~ 200,  # big value for sorting
      TRUE ~ NA_real_
    ),
    pct_label = case_when(
      !is.na(pct_change) ~ sprintf("%.0f%%", pct_change),
      is.na(pct_change) & baseline_TF == 0 & first_post_TF > 0 ~ ">100%",
      TRUE ~ NA_character_
    )
  )


### Now try with sites zscore 
## Add more info to patients 
Additional_info <- read_csv("MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/cfWGS MRDetect Blood data updated Sep with all patients.csv")

add_lookup <- Additional_info %>%
  filter(Sample_type == "Blood_plasma_cfDNA",
         !is.na(sites_rate_zscore_charm)) %>%
  transmute(
    Patient,
    add_Date  = as_date(Date_of_sample_collection_Sample_ID_Bam),
    zscore_blood_from_additional = as.numeric(sites_rate_zscore_charm)
  )

# keep a row id so no row is ever lost
dat_clean <- dat %>%
  mutate(Date = as_date(Date),
         .row_id = row_number())

dat2 <- dat_clean %>%
  # many-to-many join on Patient
  left_join(add_lookup, by = "Patient", relationship = "many-to-many") %>%
  mutate(
    add_Date = as_date(add_Date),
    # if there's no candidate date, set diff = Inf so the original row survives
    day_diff = ifelse(is.na(add_Date), Inf, abs(as.integer(add_Date - Date)))
  ) %>%
  # pick the closest Additional_info date per ORIGINAL ROW (not per Patient+Date)
  group_by(.row_id) %>%
  slice_min(order_by = day_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # fill only if original was NA and the nearest candidate is within 7 days
  mutate(
    zscore_blood = if_else(
      is.na(zscore_blood) &
        is.finite(day_diff) & day_diff <= 7 &
        !is.na(zscore_blood_from_additional),
      zscore_blood_from_additional,
      zscore_blood
    )
  ) %>%
  select(-add_Date, -day_diff, -zscore_blood_from_additional, -.row_id)

# 3) Quick report of how many were filled
filled_n <- sum(is.na(dat$zscore_blood) & !is.na(dat2$zscore_blood))
total_na_before <- sum(is.na(dat$zscore_blood))
total_na_after  <- sum(is.na(dat2$zscore_blood))

message("zscore_blood backfilled: ", filled_n,
        " (NA before: ", total_na_before, ", NA after: ", total_na_after, ")")

# 1) First non-NA baseline TF where timepoint_info is Baseline/Diagnosis
baseline_df <- dat2 %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis"),
         !is.na(zscore_BM)) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            baseline_TF   = zscore_BM,
            baseline_date = Date)

# 2) First non-NA TF strictly AFTER that baseline date
followup_df <- dat2 %>%
  inner_join(baseline_df, by = "Patient") %>%
  filter(!is.na(zscore_BM),
         Date > baseline_date) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            first_post_TF   = zscore_BM,
            first_post_date = Date)

# 3) Combine + percent change (baseline → first post)
BM_change <- baseline_df %>%
  left_join(followup_df, by = "Patient") %>%
  mutate(
    pct_change = if_else(!is.na(baseline_TF) & !is.na(first_post_TF) & baseline_TF > 0,
                         100 * (first_post_TF - baseline_TF) / baseline_TF,
                         NA_real_)
  )

## Add cohort 
BM_change <- BM_change %>%
  left_join(events %>% distinct(patient, cohort) %>% rename(Patient = patient),
            by = "Patient")

### Now for blood
# 1) First non-NA baseline TF where timepoint_info is Baseline/Diagnosis
baseline_df <- dat2 %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis"),
         !is.na(zscore_blood)) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            baseline_TF   = zscore_blood,
            baseline_date = Date)

# 2) First non-NA TF strictly AFTER that baseline date
followup_df <- dat2 %>%
  inner_join(baseline_df, by = "Patient") %>%
  filter(!is.na(zscore_blood),
         Date > baseline_date) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            first_post_TF   = zscore_blood,
            first_post_date = Date)

# 3) Combine + percent change (baseline → first post)
Blood_change <- baseline_df %>%
  left_join(followup_df, by = "Patient") %>%
  mutate(
    pct_change = if_else(!is.na(baseline_TF) & !is.na(first_post_TF) & baseline_TF > 0,
                         100 * (first_post_TF - baseline_TF) / baseline_TF,
                         NA_real_)
  )

## Add cohort 
Blood_change <- Blood_change %>%
  left_join(events %>% distinct(patient, cohort) %>% rename(Patient = patient),
            by = "Patient")

#  Add suffixes to every column EXCEPT Patient and cohort
blood_renamed <- Blood_change %>%
  rename_with(~ paste0(.x, "_blood"),
              .cols = -c(Patient, cohort))

bm_renamed <- BM_change %>%
  rename_with(~ paste0(.x, "_bm"),
              .cols = -c(Patient, cohort))

# Join side-by-side on Patient + cohort
change_combined <- full_join(
  blood_renamed,
  bm_renamed,
  by = c("Patient", "cohort")
)

### Use the tumor fraction instead since too many with NAs for the zscore due to poor quality lists

## Original method, by length of tracking 
# patient_order <- events %>%
#   group_by(cohort, patient) %>%
#   summarize(first_day = min(start_day), .groups = "drop") %>%
#   arrange(cohort, first_day) %>%
#   mutate(y = rev(row_number()))
# 
# events <- events %>% left_join(patient_order, by = c("cohort", "patient"))


## Updated, by tumor fraciton 
patient_order <- tf_change %>%
  mutate(
    cohort = factor(cohort, levels = c("Front-line cohort", "Non-front-line cohort"))
  ) %>%
  arrange(cohort, is.na(pct_for_plot), pct_for_plot) %>%     # NAs last, then ascending %Δ
  group_by(cohort) %>%
  mutate(y = row_number()) %>%                               # 1,2,3... within cohort
  ungroup() %>%
  transmute(
    cohort,
    patient = Patient,                                       # match events$patient
    y,
    pct_for_plot,
    pct_label
  )

patient_order_tf <- patient_order

# 2) Join onto eventpatient_order# 2) Join onto events for plotting
events <- events %>%
  left_join(patient_order, by = c("cohort", "patient"))


# ─────────────────────────────────────────────────────────────────────────────
# 8. FRONT-LINE COHORT PLOT
# ─────────────────────────────────────────────────────────────────────────────
front_data <- events %>% filter(cohort == "Front-line cohort")

p_front <- ggplot() +
  geom_segment(
    data = front_data %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y,
        colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  geom_segment(
    data = front_data %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y),
    colour = "black", size = 5, lineend = "round"
  ) +
  geom_point(
    data = front_data %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day, y = y, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
 # scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen") +
 # scale_shape_manual(values = shape_map, name = "Point events") +
  scale_colour_manual(
    values = chemo_cols,
    name   = "Chemotherapy regimen",
    drop   = FALSE            # <-- keep all colours in legend
  ) +
  scale_shape_manual(
    values = shape_map,
    name   = "Point events",
    drop   = FALSE            # <-- keep all shapes in legend
  ) +
  scale_x_continuous("Days from baseline", expand = expansion(mult = 0.02)) +
  scale_y_continuous(
    NULL,
    breaks = patient_order$y[patient_order$cohort == "Front-line cohort"],
    labels = patient_order$patient[patient_order$cohort == "Front-line cohort"]
  ) +
  ggtitle("Front-line cohort") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 6),
    legend.position    = "bottom",
    legend.box         = "vertical",
    plot.title         = element_text(face = "bold", hjust = 0)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 9. NON-FRONT-LINE COHORT PLOT
# ─────────────────────────────────────────────────────────────────────────────
non_data <- events %>% filter(cohort == "Non-front-line cohort")

p_non <- ggplot() +
  geom_segment(
    data = non_data %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y,
        colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  geom_segment(
    data = non_data %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y),
    colour = "black", size = 5, lineend = "round"
  ) +
  geom_point(
    data = non_data %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day, y = y, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
  scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen", guide = "none") +
  scale_shape_manual(values = shape_map, name = "Point events" , guide = "none")  +
  scale_x_continuous("Days from baseline", expand = expansion(mult = 0.02))+  scale_y_continuous(
    NULL,
    breaks = patient_order$y[patient_order$cohort == "Non-front-line cohort"],
    labels = patient_order$patient[patient_order$cohort == "Non-front-line cohort"]
  ) +
  ggtitle("Non-front-line cohort") +
  guides(colour = FALSE, shape = FALSE) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 6),
    legend.position    = "none",
  #  legend.box         = "vertical",
    plot.title         = element_text(face = "bold", hjust = 0)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 10. COMBINE & EXPORT
# ─────────────────────────────────────────────────────────────────────────────
combined_plot <- p_front / p_non +
  plot_layout(guides = "collect", heights = c(3, 1)) &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.spacing.x = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) &
  guides(
    colour = guide_legend(
      title   = "Chemotherapy regimen",
      nrow    = 2,
      byrow   = TRUE
    ),
    shape = guide_legend(
      title   = "Point events",
      nrow    = 2,
      byrow   = TRUE
    )
  )

ggsave("swimplot_by_cohort_v4.pdf", combined_plot,
       width = 12, height = 10, device = cairo_pdf)

ggsave("swimplot_by_cohort_v4.png", combined_plot,
       width = 12, height = 10, dpi = 300)



## Seperately 
# Export the front-line cohort plot
# Export the non-front-line cohort plot
ggsave(
  filename = "swimplot_nonfrontline_cohort_updated.png",
  plot     = p_non,
  width    = 10,
  height   = 2,
  dpi      = 500
)

ggsave(
  filename = "swimplot_frontline_cohort_updated.png",
  plot     = p_front,
  width    = 10,
  height   = 10,
  dpi      = 500
)




# ─────────────────────────────────────────────────────────────────────────────
# 11. SINGLE PLOT WITH BOTH COHORTS, FRONT‐LINE FIRST
# ─────────────────────────────────────────────────────────────────────────────
# A) Build a combined patient ordering: front-line patients first, then non-front-line,
#    each ordered by their first start_day.
# patient_order_combined <- events %>%
#   group_by(patient, cohort) %>%
#   summarize(first_day = min(start_day), .groups = "drop") %>%
#   mutate(
#     cohort = factor(cohort, levels = c("Front-line cohort", "Non-front-line cohort"))
#   ) %>%
#   arrange(cohort, first_day) %>%
#   mutate(y = row_number())

## To use tumor fraction 
patient_order_combined <- patient_order

# B) Join the new y positions back into `events`
events_combined <- events %>%
  left_join(patient_order_combined, by = c("patient","cohort")) %>%
  mutate(y = y.y) %>%      # pull in the new y
  select(-y.x, -y.y)       # drop the old and the suffixed copy

# C) Single swim‐plot
p_combined <- ggplot() +
  # 1) Chemo duration bars
  geom_segment(
    data = events_combined %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y, colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  # 2) Other intervals (e.g. maintenance) in black
  geom_segment(
    data = events_combined %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y),
    colour = "black", size = 5, lineend = "round"
  ) +
  # 3) Point events (samples, MRD, transplant, progression)
  geom_point(
    data = events_combined %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day, y = y, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
  # 4) Scales: chemo colours & shapes, keep all keys
  scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen", drop = FALSE) +
  scale_shape_manual(values = shape_map,   name = "Point events",          drop = FALSE) +
  # 5) X-axis: days from baseline
  scale_x_continuous("Days from baseline", expand = expansion(mult = c(0.02, 0.02))) +
  # 6) Y-axis: patient names in the combined order, top-down
  scale_y_reverse(
    breaks = patient_order_combined$y,
    labels = patient_order_combined$patient,
    name   = "Patient"
  ) +
  # 7) Theme & legend at bottom
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 10),
    axis.title.x        = element_text(size = 12, face = "bold", vjust = -0.5),
    axis.title.y        = element_text(size = 12, face = "bold", vjust =  1),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.spacing.x   = unit(0.3, "cm"),
    plot.margin        = margin(5, 10, 5, 10)
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 2, byrow = TRUE)
  )

p_combined <- p_combined +
  labs(
    title = "Swim Plot of Treatment Timelines") +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# D) Display / Save
print(p_combined)
ggsave("swimplot_combined_cohorts_v3.png", p_combined,
       width = 10, height = 12, dpi = 500)
ggsave("swimplot_combined_cohorts_wide_v3.png", p_combined,
       width = 16, height = 12, dpi = 500)



## Above is figure 1A swim plot 




#### Remove EK-09 pre-treatment and instead add star
## Add other features to make nice 
## Need to specify which patients have baseline BM and not 
### ──────────────────────────────────────────────────────────────
### 1.  Prepare the tumour-fraction table  (ord_df)
### ──────────────────────────────────────────────────────────────
# A) Standardise patient names (strip the "_Baseline")
ord_df <- read.csv("ordering_df_for_Figure_1.csv")
ord_df <- ord_df %>%
  mutate(
    patient = sub("_Baseline$", "", Sample),
    
    # map your 'Cohort' to the swim-plot cohort names
    cohort  = recode(Cohort,
                     "Train" = "Front-line cohort",
                     "Test"    = "Non-front-line cohort"),
    
    # sample-type flag for the y-axis label
    sample_type = case_when(
      Paired                     ~ "Paired",
      !Paired & is.na(TumourFraction) ~ "BM only",
      TRUE                       ~ "Blood only"
    )
  )

### ──────────────────────────────────────────────────────────────
### 2.  Re-order patients by cohort → descending tumour-fraction
### ──────────────────────────────────────────────────────────────
patient_order_combined <- ord_df %>%
  mutate(
    cohort = factor(Cohort,
                    levels = c("Train", "Test"))
  ) %>%
  group_by(cohort) %>%
  mutate(y = row_number()) %>%
  ungroup() %>%
  select(patient, cohort, y, Paired, sample_type, TumourFraction)


# tack the new 'y' onto events
events_combined <- events %>%
  left_join(patient_order_combined, by = c("patient", "cohort")) %>%
  mutate(y = y.y) %>%      # pull in the new y
  select(-y.x, -y.y)       # drop the old and the suffixed copy

patient_levels <- patient_order_combined$patient
patient_order_combined <- patient_order_combined %>%
  mutate(patient = factor(patient, levels = patient_levels))

events_combined <- events_combined %>%
  mutate(patient = factor(patient, levels = patient_levels))

# ─────────────────────────────────────────────────────────────────────────────
# 2) Tumour‐fraction strip
# ─────────────────────────────────────────────────────────────────────────────
ann_tf <- ggplot(patient_order_combined,
                 aes(x = TumourFraction, y = patient, group = 1)) +
  geom_path(colour = "grey70", size = 0.4) +
  geom_point(colour = "grey20", size = 2) +
  scale_x_continuous(
    name   = "cfDNA\ntumour fraction",
    limits = c(0, max(patient_order_combined$TumourFraction, na.rm = TRUE) * 1.05),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    limits = rev(patient_levels),
    expand = c(0, 0)    # shrink vertical padding between rows
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_text(size = 10, hjust = 1, lineheight = 0.8),
    plot.margin     = margin(0, 1, 0, 1)
  )




# ─────────────────────────────────────────────────────────────────────────────
# 3) Cohort colour bar
# ─────────────────────────────────────────────────────────────────────────────
# 1) Recode the cohort labels
#patient_order_combined <- patient_order_combined %>%
#  mutate(cohort = recode(cohort,
#                         "Front-line cohort"     = "Train",
#                         "Non-front-line cohort" = "Test"))

# 2) Define your new colour mapping
cohort_cols <- c(
  "Train" = "#1f77b4",
  "Test"  = "#e6550d"
)


ann_cohort <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = cohort)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = cohort_cols, name = "Cohort") +
  scale_x_continuous(name   = "Cohort", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(expand = c(0, 0),   limits = rev(patient_levels)) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    plot.margin     = margin(0, 0, 0, 0)
  )




# ─────────────────────────────────────────────────────────────────────────────
# 4) Paired‐status bar
# ─────────────────────────────────────────────────────────────────────────────
paired_cols <- c(
  "Paired"     = "#9467bd",  # purple
  "BM only"    = "#1f77b4",  # blue (same as Front‑line cohort)
  "Blood only" = "#CC6677"   # red (same as LEN‑based chemo)
)


ann_paired <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = factor(sample_type))) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = paired_cols, name = "Baseline samples available") +
  scale_x_continuous(name   = "Samples\navailable", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(
    limits = rev(patient_levels),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.margin     = margin(0, 0, 0, 0)
  )



# ─────────────────────────────────────────────────────────────────────────────
# 5) Swim‐plot (rebuilt on events_combined)
# ─────────────────────────────────────────────────────────────────────────────
events_combined <- events_combined %>%
  mutate(
    start_day_plot = pmax(start_day, -500),
    end_day_plot   = pmax(end_day,   -500)
  )

xmax <- max(events_combined$end_day_plot, na.rm=TRUE) * 1.05

p_swim <- ggplot() +
  geom_segment(
    data = events_combined %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day_plot, xend = end_day_plot, y = patient, yend = patient, colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  geom_segment(
    data = events_combined %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day_plot, xend = end_day_plot, y = patient, yend = patient),
    colour = "black", size = 5, lineend = "round"
  ) +
  geom_point(
    data = events_combined %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day_plot, y = patient, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
  scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen", drop = FALSE) +
  scale_shape_manual(values = shape_map, name = "Point events", drop = FALSE) +
  scale_x_continuous(
    "Days from baseline",
    limits = c(-500, xmax),
    # pick whatever breaks you like; here's an example every 500 days
    breaks = c(-500, seq(0, ceiling(xmax/500)*500, by = 500)),
    # label -750 as "<-750", everything else as its value
    labels = function(x) ifelse(x == -500, "<-500", as.character(x)),
    expand = expansion(mult = c(0.02,0.02))
  ) +
  scale_y_discrete(labels = NULL, limits = rev(levels(events_combined$patient))) +  # reverse so first is on top
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_blank(),  # we show them in ann_paired
    axis.text.x        = element_text(size = 10),             # ← bigger tick labels
    axis.title.x       = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12, face = "bold"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    legend.spacing.x   = unit(0.3, "cm"),
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(title = "Swim Plot of Treatment Timelines", y = NULL)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Patchwork: three left panels + swim plot
# ─────────────────────────────────────────────────────────────────────────────
final_plot <- ann_tf + ann_cohort + ann_paired + p_swim +
  plot_layout(widths = c(0.05, 0.035, 0.035, 0.88),
              guides = "collect") &
  theme(
    # position & spacing
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    # unify title/text/key sizes
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.8, "lines")
  ) &
  guides(
    fill   = guide_legend(
      nrow          = 1,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    colour = guide_legend(
      nrow          = 4,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    shape  = guide_legend(
      nrow          = 3,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8),
      override.aes  = list(size = 2.6)  # match your point size
    )
  )

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_updated.png",
       final_plot,
       width  = 15,
       height = 10,
       dpi    = 500)

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_wide_updated.png",
       final_plot,
       width  = 18,
       height = 10,
       dpi    = 500)



### Edit to be simpler - less shapes, clearer
# 1) Collapse to two point types

lastfu_tbl <- followup_events %>%
  filter(
    event   == "Last follow-up",
    details == "Relapsed: 0"
  ) %>%
  group_by(patient) %>%
  summarise(
    last_followup = max(start, na.rm = TRUE),
    .groups = "drop"
  )

events_combined2 <- events_combined %>%
  left_join(lastfu_tbl, by = "patient") %>%
  mutate(
    is_ongoing = is_interval &
      event == "Chemotherapy" &
      ( is.na(end) | end == last_followup )
  )

# Turn off is_ongoing for VA-02 since treatment stopped
events_combined2 <- events_combined2 %>%
  left_join(lastfu_tbl, by = "patient") %>%
  mutate(
      is_ongoing = if_else(patient == "VA-02", FALSE, is_ongoing)
  )

events_combined2 <- events_combined2 %>%
  mutate(
    point_type = case_when(
      event %in% c("BM sample collection", "cfDNA sample collection") ~ "Sample collection",
      event %in% c("MRD (MFC)", "MRD (clonoSEQ)")                       ~ "Clinical MRD assay",
      TRUE                                                             ~ NA_character_
    ),
    is_progression = event == "Progression",
    is_transplant  = event == "Transplant"
    )  

events_combined2 <- events_combined2 %>%
  mutate(
    # if event is NA, replace with "Baseline"
    event = if_else(is.na(event), "Baseline", event)
  )

events_combined2 <- events_combined2 %>%
  mutate(
    across(
      starts_with("is_"),       # pick every column whose name starts with "is_"
      ~ replace_na(.x, FALSE)   # turn any NA into FALSE
    )
  )

# 2) Shapes for the two point types
shape_map_simple <- c(
  "Sample collection" = 16,  # circle
  "Clinical MRD assay"         = 17   #  triangle
)

events_combined2 <- events_combined2 %>%
  mutate(
    chemo_group_simple = case_when(
      event != "Chemotherapy" ~ NA_character_,
      # IMiDs
      str_detect(details, regex("Lenalidomide|Rev|Pom(Dex|D)|Iberdomide", TRUE)) ~ "IMiD-based",
      # Proteasome inhibitors
      str_detect(details, regex("CY[Bb]OR|Bort|Carfil|KD|Ixa|Ixazomib", TRUE))   ~ "PI-based",
      # Antibodies
      str_detect(details, regex("Dara|Daratumumab|Elranatamab", TRUE))           ~ "Antibody-based",
      # fallback
      TRUE                                                                      ~ "Other"
    )
  )

# Previously treated 
pretx_list <- c("EK-09", "IMG-098", "SPORE_0009")
events_combined2 <- events_combined2 %>%
  mutate(
    previously_treated = patient %in% pretx_list
  )

# Limit early chemo bars
events_combined2 <- events_combined2 %>%
  mutate(
    start_day_plot = if_else(
      event == "Chemotherapy",
      pmax(start_day, 0),  # cap at 0 only for chemo
      start_day_plot               # otherwise leave as‐is
    ),
    end_day_plot = if_else(
      event == "Chemotherapy",
      pmax(end_day, 0),    # cap at 0 only for chemo
      start_day_plot                 # otherwise leave as‐is
    )
  )

## Add ongoing to patient who stopped maintenance early 
#events_combined2 <- events_combined2 %>%
#  mutate(
#    is_ongoing = if_else(
#      patient == "VA-02" & event == "Last follow-up",
#      TRUE,
#      is_ongoing
#    )
#  )

chemo_cols_simple <- c(
  "IMiD-based"     = "#009E73",  # the same green you use for MRD+/maintenance curves
  "PI-based"       = "#E69F00",  # that classic orange from your cohort bars
  "Antibody-based" = "#0072B2",  # the deep blue used on your panel A ROC lines
  "Other"          = "#999999"   # neutral grey as a fallback
)

## Other 
chemo_cols_simple <- c(
  "IMiD-based"     = "#35B779FF",  # the same green you use for MRD+/maintenance curves
  "PI-based"       = "#E69F00FF",  # that classic orange from your cohort bars
  "Antibody-based" = "#9467bd",  # the deep blue used on your panel A ROC lines
  "Other"          = "#999999"   # neutral grey as a fallback
)

# Original
chemo_cols_simple <- c(
  "IMiD-based"         = "#1b9e77",
  "PI-based"           = "#d95f02",
  "Antibody-based"     = "#7570b3",
  "Other"              = "#999999"
)

## Change to weeks 
events_combined2 <- events_combined2 %>%
  mutate(
    start_week_plot = start_day_plot/7,
    end_week_plot = end_day_plot/7
    )

## Months
events_combined2 <- events_combined2 %>%
  mutate(
    start_month_plot = start_day_plot/30.44,
    end_month_plot = end_day_plot/30.44
  )

xmax <- max(events_combined2$end_week_plot, na.rm=TRUE) * 1.05
max_months <- ceiling(xmax * 7 / 30.44)  # if xmax is in weeks, convert to months
month_breaks <- seq(0, max_months, by = 3)

## Consolidate Points
events_combined2 <- events_combined2 %>%
  mutate(
    event_type = case_when(
      is_transplant      ~ "Transplant",
      is_progression     ~ "Progression",
      is_ongoing         ~ "Ongoing",
      point_type == "Sample collection" ~ "Sample collection",
      point_type == "MRD assay"         ~ "Clinical MRD assay",
      TRUE               ~ NA_character_
    )
  )


## Reorder 
# 1) Compute total follow‑up per patient
patient_order_tbl <- events_combined2 %>% filter(!(patient == "VA-02" & event == "Last follow-up")) %>% # remove this since not plotted
  group_by(Cohort, patient) %>%
  summarise(
    # use end_day if available, otherwise start_day
    followup = max(coalesce(end_day, start_day), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # within each cohort, sort descending so longest follow‑up is first
  arrange(Cohort, desc(followup))

# 2) Pull out the ordered patient vector
patient_order <- patient_order_tbl %>% pull(patient)

## Change back to tumor fraction 
patient_order <- patient_order_tf

# 3) Re‑factor your patient column
events_combined2 <- events_combined2 %>%
  mutate(
    patient = factor(patient, levels = patient_order)
  )  



### Edit the patient labels for de-identification in figure
df <- cohort_df
make_patient_id_map <- function(df, patient_col = "Patient") {
  p <- sym(patient_col)
  
  df %>%
    distinct(!!p, .keep_all = FALSE) %>%         # keep first occurrence, preserve row order
    mutate(
      .prefix = case_when(
        str_starts(!!p, "IMG")   ~ "IMG",
        str_starts(!!p, "SPORE") ~ "SPORE",
        TRUE                     ~ "M4"
      )
    ) %>%
    group_by(.prefix) %>%
    mutate(
      .idx   = row_number(),
      # width: 2 digits unless group size >= 100, then compute digits of n()
      .width = if_else(dplyr::n() >= 100L,
                       as.integer(floor(log10(dplyr::n())) + 1L),
                       2L)
    ) %>%
    ungroup() %>%
    mutate(New_ID = paste0(.prefix, "-", str_pad(.idx, width = .width, pad = "0"))) %>%
    select(!!p, New_ID) %>%
    rename(Patient = !!p)
}


id_map <- make_patient_id_map(df)

# sanity checks
stopifnot(anyDuplicated(id_map$Patient) == 0)
stopifnot(anyDuplicated(id_map$New_ID)  == 0)

cohort_df_anon <- cohort_df %>% left_join(id_map, by = "Patient")

## Re-export supp table with new IDs 
all_events_tmp <- all_events %>%
  left_join(id_map, by = c("patient" = "Patient")) %>%
  mutate(patient = New_ID) %>%
  select(names(all_events))  # keep original column order

write_csv(all_events_tmp %>% select(-details_2), "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated_with_new_IDs.csv")

# Save as CSV
write.csv(id_map,
          file = "id_map.csv",
          row.names = FALSE)

# Save as RDS
saveRDS(id_map,
        file = "id_map.rds")

## Get labels vector
# build label lookup (id_map: Patient, New_ID)
lab_vec <- id_map$New_ID; names(lab_vec) <- id_map$Patient

lab_fun <- function(y) {
  yy  <- as.character(y)
  out <- unname(lab_vec[yy])             # remove names → pure character
  nas <- is.na(out); out[nas] <- yy[nas] # fallback to original
  out
}


# 3) Plot
p_swim <- ggplot() +
  # Chemo intervals (coloured)
  geom_segment(
    data = events_combined2 %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_month_plot, xend = end_month_plot, y = patient, yend = patient, colour = chemo_group_simple),
    linewidth = 5, lineend = "round"
  ) +
  # Other intervals (black)
  geom_segment(
    data = events_combined2 %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_month_plot, xend = end_month_plot, y = patient, yend = patient),
    colour = "black", linewidth = 5, lineend = "round"
  ) +
  # Point events (ONLY two shapes)
  geom_point(
    data = events_combined2 %>% filter(!is_interval, !is.na(point_type)),
    aes(x = start_month_plot, y = patient, shape = point_type),
    size = 2.6, stroke = 0.35, colour = "black", fill = "white"
  ) +
  # Progression tick (vertical line across the bar)
  geom_text(
    data   = events_combined2 %>% filter(is_progression),
    aes(x = start_month_plot, y = patient),
    label    = "|",          # single vertical bar
    fontface = "bold",       # bold weight
    size     = 5,            # tweak to taste
    colour   = "red"
  ) +
  # Last follow-up arrow
  geom_segment(
    data = events_combined2 %>% filter(is_ongoing),
    aes(
      x    = end_month_plot,
      xend = end_month_plot + 50/30.44,      # arrow extends 50 days to the right
      y    = patient,
      yend = patient
    ),
    colour = "black",
    linewidth = 0.6,
    lineend   = "butt",                # so the bar doesn’t cap‑over the arrow
    arrow     = arrow(
      length = unit(2, "mm"),          # size of the arrow head
      ends   = "last",                 # only draw the head at the end
      type   = "closed"                # filled triangle
    )
  ) +
  # Transplant "T"
  geom_text(
    data = events_combined2 %>% filter(is_transplant),
    aes(x = start_month_plot, y = patient, label = "T"),
    fontface = "bold", size = 3
  ) +
  # Previously treated 
  geom_point(
    data = filter(events_combined2, previously_treated),
    aes(x = -40/30.44, y = patient),
    shape  = 8,        # asterisk/star glyph
    size   = 2,        # adjust as needed
    colour = "black"
  ) +
  # Scales
  scale_colour_manual(values = chemo_cols_simple, name = "Chemotherapy regimen", drop = FALSE) +
  scale_shape_manual(values = shape_map_simple, name = "Point events", drop = FALSE) +
  scale_x_continuous(
    "\nMonths since baseline",
    limits = c(-1.4, max(month_breaks)),
    breaks = month_breaks,
    labels = as.character,
    minor_breaks = minor_breaks, # minor ticks every 3 months
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(labels = NULL, limits = rev(levels(events_combined2$patient))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),  # ← turns off vertical grid lines
    panel.grid.minor.x = element_blank(),  # ← turns off minor vertical grid lines
    axis.text.y        = element_blank(),
    axis.text.x        = element_text(size = 12),
    axis.title.x       = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12, face = "bold"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    legend.spacing.x   = unit(0.3, "cm"),
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(title = "Swim Plot of Treatment Timelines", y = NULL)

p_swim <- p_swim +
  scale_y_discrete(
    # by default it will use the patient factor levels
    limits = rev(levels(events_combined2$patient))
  ) 

## Change order 
p_swim <- p_swim + 
  scale_y_discrete(
    limits = rev(levels(events_combined2$patient)),
    breaks = NULL    # no ticks or labels
  )

ggsave("Final Tables and Figures/Test_swim2.png",
       p_swim,
       width  = 15,
       height = 10,
       dpi    = 500)

### Add back cohort and paired status 
# 2) Define your new colour mapping
cohort_cols <- c(
  "Train" = "#1f77b4",
  "Test"  = "#e6550d"
)

ann_cohort <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = cohort)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = cohort_cols, name = "Cohort") +
  scale_x_continuous(name = "Cohort", limits = c(0.5, 1.5), expand = c(0,0), breaks = NULL) +
  scale_y_discrete(
    limits = rev(patient_order),  # character levels, not numeric
    labels = lab_fun,
    expand = c(0,0)
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_text(size = 10, hjust = 1, lineheight = 0.8),
    plot.margin  = margin(0,0,0,0)
  )

# ann_cohort <- ggplot(patient_order_combined,
#                      aes(x = 1, y = patient, fill = cohort)) +
#   geom_tile(width = 0.9, height = 0.9) +
#   scale_fill_manual(values = cohort_cols, name = "Cohort") +
#   scale_x_continuous(name   = "Cohort", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
#  # scale_y_discrete(expand = c(0, 0),   limits = rev(patient_order)) +
#   scale_y_discrete(
#     limits = rev(patient_order),
#     labels = function(y) {  # anonymized labels
#       out <- lab_vec[y]
#       out[is.na(out)] <- y[is.na(out)]   # fallback: original if any ID missing
#       out
#     },
#     expand = c(0, 0)
#   ) +
#   theme_minimal(base_size = 8) +
#   theme(
#     panel.grid      = element_blank(),
#     axis.title.y    = element_blank(),
#     axis.title.x    = element_text(size = 10, hjust = 0.5, lineheight = 0.8),
#     axis.ticks.x    = element_blank(),
#     axis.ticks.y    = element_blank(),
#     axis.text.y     = element_text(size = 10, hjust = 1, lineheight = 0.8), # For patients
#     plot.margin     = margin(0, 0, 0, 0)
#   )

paired_cols <- c(
  "Paired"     = "#9467bd",  # purple
  "BM only"    = "#1f77b4",  # blue (same as Front‑line cohort)
  "Blood only" = "#CC6677"   # red (same as LEN‑based chemo)
)


ann_paired <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = factor(sample_type))) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = paired_cols, name = "Baseline samples available") +
  scale_x_continuous(name   = "Samples\navailable", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(
    limits = rev(patient_order),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.margin     = margin(0, 0, 0, 0)
  )



## Assemble final plot
final_plot <- ann_cohort + p_swim +
  plot_layout(widths = c(0.05, 0.95),
              guides = "collect") &
  theme(
    # position & spacing
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    # unify title/text/key sizes
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.8, "lines")
  ) &
  guides(
    fill   = guide_legend(
      nrow          = 1,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    colour = guide_legend(
      nrow          = 4,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    shape  = guide_legend(
      nrow          = 3,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8),
      override.aes  = list(size = 2.6)  # match your point size
    )
  )

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_updated3.png",
       final_plot,
       width  = 15,
       height = 10,
       dpi    = 500)

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_wide_updated3.png",
       final_plot,
       width  = 16.5,
       height = 10,
       dpi    = 500)



## Add legend 

# 1) Build a dummy data.frame with one row per symbol
legend_df <- data.frame(
  y     = c(3, 2, 1),
  label = c("Transplant", "Progression", "Ongoing")
)

# 2) Make the legend plot
p_symbols <- ggplot(legend_df, aes(y = y)) +
  # Title
  ggtitle("Symbols") +
  # 2a) Transplant: a bold "T"
  geom_text(
    data = subset(legend_df, label=="Transplant"),
    aes(x = 0, label = "T"),
    fontface = "bold",
    size     = 6
  ) +
  # 2b) Progression: red vertical line
  geom_text(
    data = subset(legend_df, label=="Progression"),
    aes(x = 0, y = y),
    label    = "|",          # single vertical bar
    fontface = "bold",       # bold weight
    size     = 6,            # match your other legend sizes
    colour   = "red"
  ) +
  # 2c) Ongoing: black right arrow
  geom_segment(
    data = subset(legend_df, label=="Ongoing"),
    aes(x = 0, xend = 0.6, yend = y),
    arrow   = arrow(length = unit(4, "mm"), ends="last", type="closed"),
    colour  = "black",
    size    = 1
  ) +
  # text labels
  geom_text(
    aes(x = 1.2, label = label),
    hjust = 0,
    size  = 5
  ) +
  # clean up
  scale_y_continuous(limits = c(0.5, 3.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(-0.2, 2), expand = c(0,0)) +
  theme_void() +
  theme(
    plot.title      = element_text(size=14, face="bold", hjust=0),
    plot.margin     = margin(5,5,5,5)
  )

# 3) (Optional) view it
print(p_symbols)











#### Below here is testing 
## Just legend - not used in final code
# 1. Dummy data for chemo colours
df_chemo <- data.frame(
  chemo_group = factor(names(chemo_cols), levels = names(chemo_cols)),
  x = seq_along(names(chemo_cols)),
  y = 1
)

# 2. Dummy data for point-event shapes
df_shape <- data.frame(
  event = factor(names(shape_map), levels = names(shape_map)),
  x = seq_along(names(shape_map)),
  y = 2
)

# 3. Build dummy plot to generate full legend
legend_plot <- ggplot() +
  geom_point(
    data = df_chemo,
    aes(x = x, y = y, colour = chemo_group),
    size = 5
  ) +
  geom_point(
    data = df_shape,
    aes(x = x, y = y, shape = event),
    size = 5
  ) +
  scale_colour_manual(
    name   = "Chemotherapy regimen",
    values = chemo_cols,
    drop   = FALSE
  ) +
  scale_shape_manual(
    name   = "Point events",
    values = shape_map,
    drop   = FALSE
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.spacing.x = unit(0.3, "cm"),
    legend.spacing.y = unit(0.3, "cm")
  )

# 4. Extract the legend grob
full_legend <- get_legend(legend_plot)

# 5. (Optional) Draw to the current device
grid.newpage()
grid.draw(full_legend)

# 6. Save it as its own file
ggsave("swimplot_full_legend.png", full_legend,
       width = 14, height = 2, dpi = 500)

# 3. Draw the legend on a blank page (for interactive viewing)
grid::grid.newpage()
grid::grid.draw(full_legend)

# 4. Save the legend alone to file
ggsave("full_swimplot_legend.pdf", full_legend,
       width = 8, height = 2, device = cairo_pdf)
ggsave("full_swimplot_legend.png", full_legend,
       width = 8, height = 2, dpi = 300)


