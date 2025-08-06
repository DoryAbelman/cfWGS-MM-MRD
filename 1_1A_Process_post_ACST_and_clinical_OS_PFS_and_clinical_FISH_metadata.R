# ──────────────────────────────────────────────────────────────────────────────
#  1_1A_Process_post_ACST_and_clinical_OS_PFS_and_clinical_FISH_metadata.R
#
#  Purpose:
#    • Ingest and clean clinical/transplant/progression data for SPORE, M4, IMMAGINE cohorts
#    • Normalize column names, extract timepoints, treatments, and FISH flags
#    • Merge ASCT (autologous stem‐cell transplant) dates with relapse/progression dates
#    • Compute time‐to‐relapse metrics (days, flags within 12/18/24 months)
#    • Integrate overall survival (OS) follow‐up and vital‐status information
#    • Export tidy “tidy_timepoints.csv”, “tidy_treatments.csv”, “tidy_progression.csv”,
#      “tidy_fish.csv”, “ASCT_relapse_summary.csv”, and a combined OS/PFS master table
#
#  Usage:
#    1. Place this script in your project root.
#    2. Make sure the following source files exist (paths relative to project root):
#         • Clinical data/SPORE/cfDNA project and clinical data just clinical dates treatment and progression.xlsx
#         • M4_CMRG_Data/June update/M4_COHORT_STEM_CELL_TRANSPLANT.xlsx
#         • Clinical data/IMMAGINE/IMG_request_20241009 (2).xlsx
#         • Clinical data/IMMAGINE/Cleaned transplant dates just dates.xlsx
#         • Clinical data/IMMAGINE/Extracted clinical MRD data.xlsx (sheet 1)
#         • Clinical data/IMMAGINE/Cleaned_Patient_Follow-Up_Table_IMMAGINE.csv
#         • Clinical data/SPORE/SPORE_OS_info.xlsx
#         • M4_CMRG_Data/M4_COHORT_DEMO.xlsx
#
#    3. From R (in project root), run:
#         source("process_post_ACST_and_clinical_OS_PFS.R")
#
#  *This script is intended as a one-shot cleaning script to be run after the previous 1_0A script* 
#
#  Author:        Dory Abelman
#  Last updated:  2025-05-13
#
#  Required packages:
#    tidyverse, readxl, lubridate, janitor, purrr, glue
# ──────────────────────────────────────────────────────────────────────────────


# ─── 0.  Library loading & file checks ────────────────────────────────────────
library(readxl)      # for reading Excel files
library(dplyr)
library(tidyr)
library(lubridate)   # date parsing
library(janitor)     # cleaning column names, removing empty rows/cols
library(stringr)
library(purrr)       # for functional programming (map, etc.)
library(glue)        # for building file‐paths and messages

# Helper to verify that required input files exist
ensure_exists <- function(path) {
  if (!file.exists(path)) {
    stop(glue("❌ Required file not found: {path}"), call. = FALSE)
  }
  invisible(path)
}

# List all input files needed by this script:
input_files <- c(
  "Clinical data/SPORE/cfDNA project and clinical data just clinical dates treatment and progression.xlsx",
  "M4_CMRG_Data/June update/M4_COHORT_STEM_CELL_TRANSPLANT.xlsx",
  "Clinical data/IMMAGINE/IMG_request_20241009 (2).xlsx",
  "Clinical data/IMMAGINE/Cleaned transplant dates just dates.xlsx",
  "Clinical data/IMMAGINE/Extracted_clinical_MRD_data.xlsx",
  "Clinical data/IMMAGINE/Cleaned_Patient_Follow-Up_Table_IMMAGINE.csv",
  "Clinical data/SPORE/SPORE_OS_info.xlsx",
  "M4_CMRG_Data/M4_COHORT_DEMO.xlsx"
)

# Verify existence of each input
walk(input_files, ensure_exists)

message("✅ All required input files are present.")





#### Now get the post-ACST dates for each project  

## For SPORE
# ──────────────────────────────────────────────────────────────────────────────
#  1.  Ingest the raw sheet ----------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
raw   <- read_excel("Clinical data/SPORE/cfDNA project and clinical data just clinical dates treatment and progression.xlsx",
                    sheet = 1,                          # adjust if needed
                    na    = c("NA","N/A","", " ")) |>   # treat blank‑likes as NA
  janitor::remove_empty("rows")                 # drop fully‑empty rows

# readxl will automatically deduplicate identical headers by appending “…1”, “…2”.
# Keep a copy for inspection if you like:
write_lines(names(raw), "original_headers.txt")

# ──────────────────────────────────────────────────────────────────────────────
#  2.  Normalise column names --------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
# make syntactically valid & lower‑snake‑case names *once*
names(raw) <- names(raw) |> 
  janitor::make_clean_names()              # "treatment_1_prior_timepoint"

# collapse verbose headers to simple patterns we can pivot on
names(raw) <- names(raw) |>
  str_replace("^treatment_(\\d+)_prior_timepoint$",   "treatment_\\1") |>
  str_replace("^start_date_(\\d+)$",                  "treatment_start_\\1") |>
  str_replace("^relapse_progression_date_(\\d+)$",    "relapse_\\1")

# quick sanity check for accidental duplicates created by the collapse:
stopifnot(!anyDuplicated(names(raw)))

# ──────────────────────────────────────────────────────────────────────────────
#  3.  Core identifiers --------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
key_cols <- c(
  "patient",
  "fish",
  "fish_date",
  "timepoint_of_interest",
  "status_at_timepoint",
  "flow_monotypic_pc_per_total_events_percent"  # ← updated
)

core <- raw |> 
  select(all_of(key_cols))

# ──────────────────────────────────────────────────────────────────────────────
#  4.  Time‑point table --------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
timepoints <- core |>
  select(patient,
         timepoint_date  = timepoint_of_interest,
         status          = status_at_timepoint,
         flow_pc_total   = flow_monotypic_pc_per_total_events_percent) |>
  mutate(across(contains("date"), ~ suppressWarnings(as_date(.x)))) |>
  drop_na(timepoint_date) |>
  distinct()

# ──────────────────────────────────────────────────────────────────────────────
#  5.  Treatment table (wide → long) ------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
treat_long <- raw |>
  select(patient,
         matches("^treatment_\\d+$|^treatment_start_\\d+$")) |>
  pivot_longer(
    -patient,
    names_to  = c(".value", "line"),                   # .value keeps column group
    names_pattern = "^(treatment|treatment_start)_(\\d+)$"
  ) |>
  rename(regimen = treatment, start_date = treatment_start) |>
  mutate(
    line       = parse_integer(line),
    start_date = suppressWarnings(as_date(start_date))
  ) |>
  drop_na(regimen, start_date) |>
  distinct()

# ──────────────────────────────────────────────────────────────────────────────
#  6.  Progression / relapse table --------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
# 1. identify your relapse columns
relapse_cols <- names(raw)[ grepl("^relapse_\\d+$", names(raw)) ]

# 2. coerce them all to character (so pivot_longer can combine them)
raw <- raw %>%
  mutate(across(all_of(relapse_cols), ~ as.character(.x)))

# 3. now pivot and then parse back into a Date
### This doesn't take everything but have already from above so ok
progression <- raw %>%
  select(patient, all_of(relapse_cols)) %>%
  pivot_longer(
    -patient,
    names_to      = "event",
    values_to     = "relapse_date",
    names_pattern = "^relapse_(\\d+)$"
  ) %>%
  mutate(
    event        = as.integer(event),
    relapse_date = ymd(relapse_date)      # lubridate::ymd handles "YYYY-MM-DD" strings
  ) %>%
  drop_na(relapse_date) %>%
  distinct()

# ──────────────────────────────────────────────────────────────────────────────
#  7.  FISH table (plus quick binary flags) -----------------------------------
# ──────────────────────────────────────────────────────────────────────────────
fish <- core |>
  select(patient,
         fish_abnormalities = fish,
         fish_date) |>
  mutate(
    fish_date     = suppressWarnings(as_date(fish_date)),
    # 1p, 13, 14, 17p deletions / monosomies
    del_1p        = str_detect(fish_abnormalities, regex("1p deletion|monosomy 1p", TRUE)),
    del_13        = str_detect(fish_abnormalities, regex("13q deletion|monosomy 13", TRUE)),
    del_14        = str_detect(fish_abnormalities, regex("monosomy 14|14q deletion", TRUE)),
    del_17p       = str_detect(fish_abnormalities, regex("17p deletion|monosomy 17|tp53", TRUE)),
    
    # 1q gains
    gain_1q       = str_detect(fish_abnormalities, regex("1q gain|1q duplication", TRUE)),
    
    # common trisomies
    trisomy_3     = str_detect(fish_abnormalities, regex("trisomy 3", TRUE)),
    trisomy_7     = str_detect(fish_abnormalities, regex("trisomy 7", TRUE)),
    trisomy_9     = str_detect(fish_abnormalities, regex("trisomy 9", TRUE)),
    trisomy_15    = str_detect(fish_abnormalities, regex("trisomy 15", TRUE)),
    
    # hyperdiploid / tetraploid signature
    hyperdiploid  = str_detect(fish_abnormalities, regex("hyperdiploid|tetraploid", TRUE)),
    
    # canonical translocations
    t_11_14       = str_detect(fish_abnormalities, regex("t\\s*\\(11;14", TRUE)),
    t_4_14        = str_detect(fish_abnormalities, regex("t\\s*\\(4;14",  TRUE)),
    t_14_16       = str_detect(fish_abnormalities, regex("t\\s*\\(14;16", TRUE)),
    t_14_20       = str_detect(fish_abnormalities, regex("t\\s*\\(14;20", TRUE)),
    
    # CKS1B region duplication
    cks1b_dup     = str_detect(fish_abnormalities, regex("cks1b", TRUE))
  ) |>
  distinct()

fish <- fish %>% filter(!is.na(fish_abnormalities))

# ──────────────────────────────────────────────────────────────────────────────
#  8.  Tidy outputs ------------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
list(
  timepoints  = timepoints,
  treatments  = treat_long,
  progression = progression,
  fish        = fish
) |> purrr::iwalk(\(df, nm) {
  write_csv(df, glue::glue("tidy_{nm}.csv"))
})








#### Now get all the ACST dates for M4
# ──────────────────────────────────────────────────────────────────────────────
#  1.  Read the raw table ------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
df_raw <- read_excel(
  "M4_CMRG_Data/June update/M4_COHORT_STEM_CELL_TRANSPLANT.xlsx",    # ← your file here
  sheet = 1,                          # adjust if needed
  na    = c("NA","N/A","", " ")) |>   # treat blank‑likes as NA
  janitor::remove_empty("rows")                 # drop fully‑empty rows

# ──────────────────────────────────────────────────────────────────────────────
#  2.  Clean names & drop empty columns ---------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
df <- df_raw |>
  clean_names() |>
  remove_empty("cols")


df <- df %>%
  filter(
    # keep rows where either
    #  • m4_id is ≤ 5 chars long, OR
    #  • it matches two letters, a dash, two digits
    nchar(m4_id) <= 5 |
      str_detect(m4_id, "^[A-Z]{2}-\\d{2}$")
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3.  Parse dates & numbers ---------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
date_cols <- c(
  "transplant_date",
  "admission_date",
  "discharge_date",
  "immun_response_date",
  "best_response_date",
  "progression_date"
)

num_cols <- c(
  "plasma_cell_percent",
  "plasma_cell_col_res_percent",
  "injected_cd34"
)

df_clean <- df %>%
  mutate(
    across(all_of(date_cols), ~ {
      x <- .x
      # if it’s pure numeric, treat it as Excel serial
      if (is.numeric(x)) {
        # Excel’s “1900” origin (adjust for the 1900‑leap‑year bug if needed):
        excel_numeric_to_date(x)
      } else {
        # otherwise parse character dates flexibly
        parse_date_time(as.character(x),
                        orders = c("dmy", "ymd"),
                        quiet  = TRUE)
      }
    }),
    # now fix injected_cd34: extract the coefficient then multiply by 10^6
    injected_cd34 = parse_number(as.character(injected_cd34),
                                 locale = locale(decimal_mark = ".")) * 1e6,
    
    # re‑parse your other numeric %% columns as before
    plasma_cell_percent          = parse_number(as.character(plasma_cell_percent)),
    plasma_cell_col_res_percent  = parse_number(as.character(plasma_cell_col_res_percent))
  )


# ──────────────────────────────────────────────────────────────────────────────
# 4.  Rename & factorize key columns (optional) --------------------------------
# ──────────────────────────────────────────────────────────────────────────────
df_clean <- df_clean |>
  rename_with(~ str_replace_all(.x, "_date$", "_dt")) |>
  mutate(
    procedure_type  = as_factor(procedure_type_text),
    transplant_type = as_factor(transplant_type),
    line_of_treatment = as.integer(line_of_treatment)
  )







#### Do for IMMAGINE 
# ──────────────────────────────────────────────────────────────────────────────
#  1.  read the worksheet ------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
raw <- read_excel(
  "Clinical data/IMMAGINE/IMG_request_20241009 (2).xlsx",   # ← your file
  sheet = 1,
  na    = c("", "NA", "N/A", " ", "Pending")
) %>% 
  remove_empty("cols")              %>%  # drop blank columns
  clean_names()                     %>%  # nicer snake‑case names
  mutate(across(everything(), as.character))  # keep everything character

# ──────────────────────────────────────────────────────────────────────────────
#  2.  helper: extract & parse *all* dates in a string -------------------------
# ──────────────────────────────────────────────────────────────────────────────
extract_dates <- function(x) {
  # 1) coerce to character
  x_chr <- as.character(x)
  
  # 2) pull Excel‑serial candidates (five digits) and textual dates
  serials   <- str_extract_all(x_chr, "\\b\\d{5}\\b", simplify = FALSE) |> unlist()
  text_dates <- str_extract_all(
    x_chr,
    "(\\d{1,2}[/-]\\d{1,2}[/-]\\d{2,4})|([A-Za-z]{3,9}[ -]\\d{2,4})",
    simplify = FALSE
  ) |> unlist()
  
  # 3) combine, unique & drop empty
  candidates <- unique(c(serials, text_dates))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) == 0) return(as_date(NA))
  
  # 4) parse each candidate
  parsed <- map(candidates, function(d) {
    if (str_detect(d, "^\\d{5}$")) {
      excel_numeric_to_date(as.numeric(d))
    } else {
      parse_date_time(d,
                      orders = c("dmy", "dmY", "dmyy", "ymd", "my", "BdY"),
                      quiet  = TRUE
      ) %>% as_date()
    }
  })
  
  # 5) drop NAs, unlist, sort & unique
  parsed <- purrr::compact(parsed)    # <-- explicitly purrr::compact
  parsed <- unlist(parsed)
  sort(unique(parsed))
}

# ──────────────────────────────────────────────────────────────────────────────
#  3.  grab the first date in each time‑point column ---------------------------
# ──────────────────────────────────────────────────────────────────────────────
# 2) read in the Excel file (treat em‑dash or blank as NA)
df <- read_excel(
  path  = "Clinical data/IMMAGINE/Cleaned transplant dates just dates.xlsx",
  sheet = 1, 
  na    = c("", "—", "-", "NA")
)

# 3) pivot to long, clean names, split multi‑dates, parse into Date class
df_long <- df %>%
  pivot_longer(
    cols      = -Patient,
    names_to  = "Treatment",
    values_to = "Date"
  ) %>%
  mutate(Treatment = recode(Treatment,
                            "Diagnosis date(s)"                 = "Diagnosis",
                            "Induction date(s)"                 = "Induction",
                            "Transplant date(s)"                = "Transplant",
                            "Maintenance date(s)"               = "Maintenance",
                            "MRD tests (date)"                  = "MRD test",
                            "Relapse / Progression date(s)"     = "Relapse/Progression"
  )) %>%
  separate_rows(Date, sep = ",\\s*") %>%
  filter(!is.na(Date), Date != "") %>%
  mutate(Date = parse_date_time(Date,
                                orders = c("Y-m-d", "Y-m", "Y")))




# ──────────────────────────────────────────────────────────────────────────────
#  4.  FISH table  -------------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
fish <- raw %>% 
  select(id, fish) %>% 
  rename(fish_text = fish)

fish_flags <- fish %>% 
  mutate(
    del_13q     = str_detect(fish_text, regex("13q|monosomy 13", TRUE)),
    del_17p     = str_detect(fish_text, regex("17p|tp53", TRUE)),
    loss_1p     = str_detect(fish_text, regex("1p", TRUE)),
    gain_1q     = str_detect(fish_text, regex("1q( gain|dup)", TRUE)),
    trisomy_3   = str_detect(fish_text, regex("trisomy 3", TRUE)),
    trisomy_7   = str_detect(fish_text, regex("trisomy 7", TRUE)),
    trisomy_9   = str_detect(fish_text, regex("trisomy 9", TRUE)),
    trisomy_15  = str_detect(fish_text, regex("trisomy 15", TRUE)),
    t_11_14     = str_detect(fish_text, regex("t\\s*\\(11;14", TRUE)),
    t_4_14      = str_detect(fish_text, regex("t\\s*\\(4;14",  TRUE)),
    t_14_16     = str_detect(fish_text, regex("t\\s*\\(14;16", TRUE)),
    t_14_20     = str_detect(fish_text, regex("t\\s*\\(14;20", TRUE)),
    cks1b_dup   = str_detect(fish_text, regex("cks1b", TRUE))
  )

fish_flags <- fish_flags %>% filter(!is.na(fish_text))





#### Now join into one big table for all the ASCT dates
# 1) From df_long: pick the “Transplant” Treatment
df_long_xplant <- df_long %>%
  filter(Treatment == "Transplant") %>%
  select(Patient, Date) %>%
  mutate(
    Source = "df_long",
    Date   = as_date(Date)
  )

# 2) From df_clean: grab transplant_dt
df_clean_xplant <- df_clean %>%
  rename(Patient = m4_id) %>%
  select(Patient, Date = transplant_dt) %>%
  filter(!is.na(Date)) %>%
  mutate(
    Source = "df_clean",
    Date   = as_date(Date)
  )

# 3) From treat_long: pick ASCT or ASTC #2 regimens
treat_long_xplant <- treat_long %>%
  filter(regimen %in% c("ASCT", "ASTC #2")) %>%
  rename(Patient = patient,
         Date    = start_date) %>%
  filter(!is.na(Date)) %>%
  mutate(
    Source = "treat_long",
    Date   = as_date(Date)
  )

# 4) Combine all three and sort
all_xplants <- bind_rows(df_long_xplant,
                         df_clean_xplant,
                         treat_long_xplant) %>%
  # optional: drop exact duplicates 
  distinct(Patient, Date, Source, .keep_all = TRUE) %>%
  arrange(Patient, Date)


### Now calculate time to relapse 

# 1) Keep only the transplant_date
acst_tbl <- all_xplants %>%
  select(Patient, transplant_date = Date)

# 2) Join to all progression dates on or after each transplant
acst_joined <- acst_tbl %>%
  left_join(Relapse_dates_full, by = "Patient") %>%
  filter(Progression_date >= transplant_date)

# 3) For each transplant, pick the *earliest* progression_date
next_relapse <- acst_joined %>%
  group_by(Patient, transplant_date) %>%
  slice_min(Progression_date, with_ties = FALSE) %>%
  ungroup() %>%
  rename(nearest_relapse_date = Progression_date) %>%
  mutate(
    days_to_relapse = as.integer(nearest_relapse_date - transplant_date)
  )

# 4) Patients/transplants with *no* later relapse
never_relapsed <- acst_tbl %>%
  anti_join(next_relapse, by = c("Patient", "transplant_date")) %>%
  mutate(
    nearest_relapse_date = as_date(NA),
    days_to_relapse       = NA_integer_
  )

# 5) Combine and add flags
final_tbl <- bind_rows(next_relapse, never_relapsed) %>%
  arrange(Patient, transplant_date) %>%
  mutate(
    relapsed      = if_else(!is.na(nearest_relapse_date), 1, 0),
    within_12mo   = if_else(relapsed == 1 & days_to_relapse <= 365, 1, 0),
    within_18mo   = if_else(relapsed == 1 & days_to_relapse <= 548, 1, 0),
    within_24mo   = if_else(relapsed == 1 & days_to_relapse <= 730, 1, 0)
  )

# Inspect
print(final_tbl)

# Export as CSV
write.csv(final_tbl, "relapse_flags_table.csv", row.names = FALSE)

# Export as RDS
saveRDS(final_tbl, file = "relapse_flags_table.rds")


### Re-correct FISH flags 
## --- 1) Regex patterns -----------------------------------------------------
fish_patterns <- list(
  del_13q   = c("monosomy 13",       "del\\s*13"),
  del_17p   = c("17p",               "del\\s*17"),
  loss_1p   = c("1p32\\.3 loss",     "loss 1p", "1p loss",  "del\\s*1p"),
  gain_1q   = c("1q21\\.3 gain",     "1q gain",   "gain 1q", "cks1b gain"),
  trisomy_3 = "trisomy\\s*3",
  trisomy_7 = "trisomy\\s*7",
  trisomy_9 = "trisomy\\s*9",
  trisomy_15= "trisomy\\s*15",
  t_11_14   = "11;\\s*14",
  t_4_14 = c(
    "4;\\s*14",
    "igh\\s*,\\s*fgfr3\\s*translocation"   # can stay mixed‐case
  ),  
  t_14_16   = "14;\\s*16",
  t_14_20   = "14;\\s*20",
  cks1b_dup = c("cks1b dup", "cks1b gain")
)

## --- 2) Helper that returns a logical column -------------------------------
flag_it <- function(text, pat) {
  out <- str_detect(
    string      = str_to_lower(text %||% ""),      # handle NA safely
    pattern     = str_c("(", str_c(pat, collapse = "|"), ")")
  )
  replace_na(out, FALSE)
}

## --- 3) Re‑score every column in one mutate --------------------------------
fish_flags_corrected <- fish_flags %>%          # keeps original columns
  mutate(across(
    .cols = names(fish_patterns),
    .fns  = ~ flag_it(fish_text, fish_patterns[[cur_column()]]),
    .names = "{.col}"                                # overwrite / create
  )) %>%
  # re‑order to your preferred column order -------------
relocate(id, fish_text, everything())


#  ──────────────────────────────────────────────────────────────────────────────
#   Export key tables for each cohort                                         

export_dir <- "Clinical data/Exported clinical data April 2025"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# SPORE
write_csv(timepoints, file.path(export_dir, "SPORE_timepoints.csv"))
write_csv(treat_long,  file.path(export_dir, "SPORE_treatments.csv"))
write_csv(progression,  file.path(export_dir, "SPORE_progression.csv"))
write_csv(fish,         file.path(export_dir, "SPORE_fish_flags.csv"))

# M4
write_csv(df_clean,     file.path(export_dir, "M4_cohort_df_clean.csv"))
# if you have a distilled ASCT table for M4 only:
write_csv(df_clean %>% 
            select(m4_id, transplant_dt),
          file.path(export_dir, "M4_ASCT_dates.csv"))

# IMMAGINE
write_csv(df_long,      file.path(export_dir, "IMMAGINE_dates_long.csv"))
write_csv(fish_flags_corrected,   file.path(export_dir, "IMMAGINE_fish_flags.csv")) ### This is wrong, edit

# Combined ASCT/relapse summary
write_csv(all_xplants,  file.path(export_dir, "all_cohorts_ASCT_dates.csv"))
write_csv(final_tbl,     file.path(export_dir, "ASCT_relapse_summary.csv"))


#  ──────────────────────────────────────────────────────────────────────────────
#  Optional plots                                                         
# A. Histogram of days-to-relapse across all cohorts
p1 <- final_tbl %>% 
  filter(!is.na(days_to_relapse)) %>%
  ggplot(aes(x = days_to_relapse)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribution of Days to Relapse Post‑Transplant",
    x     = "Days to Relapse",
    y     = "Count"
  )
ggsave("Exported clinical data April 2025/days_to_relapse_histogram.png", p1, width = 6, height = 4)

# 10b. Proportion relapsed within each window
p2 <- final_tbl %>%
  summarise(
    `≤ 12 mo` = mean(within_12mo),
    `≤ 18 mo` = mean(within_18mo),
    `≤ 24 mo` = mean(within_24mo)
  ) %>%
  pivot_longer(everything(), names_to = "Window", values_to = "Proportion") %>%
  ggplot(aes(x = Window, y = Proportion)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Transplants Followed by Relapse Within Time Windows",
    x     = "Time Window After ACST",
    y     = "Proportion Relapsed"
  )
ggsave("Exported clinical data April 2025/relapse_window_proportions.png", p2, width = 6, height = 4)



# Filter and select columns for SPORE cases
spore_clinical_subset <- combined_clinical_data_updated %>%
  filter(Study == "SPORE") %>%
  select(Patient, Date_of_sample_collection, timepoint_info)

# Export to CSV
write_csv(spore_clinical_subset, "SPORE_Patient_Timeline_ClinicalInfo.csv")




### To the final table, add the date of last followup and update relapse 

## Load for IMMAGINE 
IMMAGINE_OS <- read.csv("Clinical data/IMMAGINE/Cleaned_Patient_Follow-Up_Table_IMMAGINE.csv")

# 1) Clean IMMAGINE_OS
# 1) Compute each patient’s first sample date
first_sample_dates <- combined_clinical_data_updated %>%
  group_by(Patient) %>%
  summarise(
    first_sample_date = min(Date_of_sample_collection, na.rm = TRUE),
    .groups = "drop"
  )

# 2) Clean your IMMAGINE_OS table (empty → NA, parse dates)
immagine_clean <- IMMAGINE_OS %>%
  mutate(across(c(Maintenance_Start_Date, Relapse1_Date:Relapse3_Date, Last_Followup_Date),
                ~ na_if(., "") %>% as.Date())) %>%
  rename(Patient = Patient_ID)

# 3) Pivot out only the relapse‐date columns, join the first‐sample date, 
#    and filter to keep only relapse_dates *after* that first_sample_date
immagine_relapse <- immagine_clean %>%
  pivot_longer(
    cols      = c(Relapse1_Date, Relapse2_Date, Relapse3_Date),
    names_to  = "which_relapse",
    values_to = "relapse_date"
  ) %>%
  filter(!is.na(relapse_date)) %>%
  left_join(first_sample_dates, by = "Patient") %>%
  filter(relapse_date > first_sample_date) %>%
  group_by(Patient) %>%
  summarise(
    nearest_relapse_date_img = min(relapse_date),
    .groups = "drop"
  )

## Correct the "-"
# Fix the Patient codes in your main table:
final_tbl <- final_tbl %>%
  mutate(
    Patient = str_replace_all(Patient, "[\u2010\u2011\u2012\u2013]", "-")
  )

# And also in the IMMAGINE data:
immagine_clean <- immagine_clean %>%
  mutate(
    Patient = str_replace_all(Patient, "[\u2010\u2011\u2012\u2013]", "-")
  )

combined_tbl <- final_tbl %>%
  full_join(immagine_relapse, by = "Patient") %>%
  mutate(
    # overwrite only for those IMG‐ patients
    nearest_relapse_date = coalesce(nearest_relapse_date_img, nearest_relapse_date),
    days_to_relapse      = as.integer(nearest_relapse_date - transplant_date),
    relapsed             = if_else(!is.na(nearest_relapse_date), 1, 0),
    within_12mo          = if_else(relapsed == 1 & days_to_relapse <= 365, 1, 0),
    within_18mo          = if_else(relapsed == 1 & days_to_relapse <= 548, 1, 0),
    within_24mo          = if_else(relapsed == 1 & days_to_relapse <= 730, 1, 0)
  ) %>%
  select(-nearest_relapse_date_img)


## Now add the date of last followup 
# 1) Make sure immagine_clean$Last_Followup_Date is Date
immagine_clean <- immagine_clean %>%
  mutate(Last_Followup_Date = as_date(Last_Followup_Date))

# 2) Left-join those two fields onto combined_tbl
combined_tbl <- combined_tbl %>%
  left_join(
    immagine_clean %>% 
      select(Patient, Last_Followup_Date, Status_Last_FU),
    by = "Patient"
  ) 


### Pull fo M4 
M4_OS_info <- read_excel("M4_CMRG_Data/M4_COHORT_DEMO.xlsx")

### See if any later than latest dates table from above

# 1) Make sure both date‐columns are real Dates
latest_dates2 <- latest_dates %>%
  mutate(latest_date = as_date(latest_date))

M4_OS2 <- M4_OS_info %>%
  mutate(DATE_OF_LAST_FOLLOWUP = as_date(DATE_OF_LAST_FOLLOWUP))

# 2) Join them by patient ID
#    Option A: if your ID is stored in M4_OS_info$M4_id
joined_A <- latest_dates2 %>%
  inner_join(M4_OS2, by = c("Patient" = "M4_id"))

# 3) Flag and filter
out_of_date_A <- joined_A %>%
  mutate(later_than_fu = latest_date > DATE_OF_LAST_FOLLOWUP) %>%
  filter(later_than_fu)

out_of_date_B <- joined_A %>%
  mutate(later_than_fu = latest_date < DATE_OF_LAST_FOLLOWUP) %>%
  filter(later_than_fu)

out_of_date_A %>%
  select(Patient, latest_date, DATE_OF_LAST_FOLLOWUP)

out_of_date_B %>%
  select(Patient, latest_date, DATE_OF_LAST_FOLLOWUP)

## Reconcile to always have the latest date 
# 1. Prepare latest_dates: rename Patient → M4_id, ensure it’s a Date
latest_dates2 <- latest_dates %>%
  rename(M4_id = Patient) %>%
  mutate(latest_date = as_date(latest_date))

# 2. Update M4_OS_info
M4_OS_updated <- M4_OS_info %>%
  # make sure your OS dates are Date class
  mutate(DATE_OF_LAST_FOLLOWUP = as_date(DATE_OF_LAST_FOLLOWUP)) %>%
  
  # join in the new latest_date
  left_join(latest_dates2, by = "M4_id") %>%
  
  # take the later of the two dates whenever latest_date is non-NA
  mutate(
    DATE_OF_LAST_FOLLOWUP = pmax(DATE_OF_LAST_FOLLOWUP, latest_date, na.rm = TRUE)
  ) %>%
  
  # clean up
  select(-latest_date)

M4_OS_updated <- M4_OS_updated %>% select(-DOB, -CONSENT_OBTAINED, -CONSENT_DATE, -study_patient_id) ## Remove unnecessary columns

### Now integrate this with the relapse info in final table 
# 1) Prep the OS table: rename and ensure dates are Date class
M4_OS2 <- M4_OS_updated %>%
  rename(Patient = M4_id) %>%
  mutate(
    DATE_OF_LAST_FOLLOWUP = as_date(DATE_OF_LAST_FOLLOWUP),
    DATE_OF_DEATH          = as_date(DATE_OF_DEATH)
  ) %>%
  select(
    Patient,
    VITAL_STATUS,
    CAUSE_OF_DEATH,
    DATE_OF_DEATH,
    DATE_OF_LAST_FOLLOWUP
  )

# 2) Join onto your final_tbl
combined_tbl <- combined_tbl %>%
  left_join(M4_OS2, by = "Patient")


# 2) Left-join those two fields onto combined_tbl
combined_tbl <- combined_tbl %>%
  # 3) Coalesce so M4_OS_info values stay when present
  mutate(
    DATE_OF_LAST_FOLLOWUP = coalesce(DATE_OF_LAST_FOLLOWUP, Last_Followup_Date),
    VITAL_STATUS          = coalesce(VITAL_STATUS, Status_Last_FU)
  ) %>%
  # 4) Drop the helper cols
  select(-Last_Followup_Date, -Status_Last_FU)

### Now add for SPORE 
## Spore already have in the SPORE OS info

# 1) Clean and rename the SPORE OS info so its columns match combined_tbl
spore_OS2 <- spore_OS_info %>%
  rename(
    VITAL_STATUS_SPORE       = `Current status`,
    DATE_OF_LAST_FOLLOWUP_SPORE = `Status last follow up`
  ) %>%
  mutate(
    DATE_OF_LAST_FOLLOWUP_SPORE = as_date(DATE_OF_LAST_FOLLOWUP_SPORE)
  )

# 2) Left‐join onto your existing combined_tbl
combined_tbl <- combined_tbl %>%
  left_join(spore_OS2, by = "Patient") %>%
  
  # 3) Coalesce so you keep existing M4/IMG values and fill in SPORE where missing
  mutate(
    DATE_OF_LAST_FOLLOWUP = coalesce(DATE_OF_LAST_FOLLOWUP, DATE_OF_LAST_FOLLOWUP_SPORE),
    VITAL_STATUS          = coalesce(VITAL_STATUS, VITAL_STATUS_SPORE)
  ) %>%
  
  # 4) Drop the helper columns
  select(-DATE_OF_LAST_FOLLOWUP_SPORE, -VITAL_STATUS_SPORE)

# Check date of death 
combined_tbl <- combined_tbl %>%
  # ensure your death column is Date
  mutate(DATE_OF_DEATH = as_date(DATE_OF_DEATH)) %>%
  
  # for anyone flagged as Deceased but missing DATE_OF_DEATH,
  # set DATE_OF_DEATH = DATE_OF_LAST_FOLLOWUP
  mutate(
    DATE_OF_DEATH = coalesce(
      DATE_OF_DEATH,
      if_else(VITAL_STATUS == "Deceased", DATE_OF_LAST_FOLLOWUP, as.Date(NA))
    )
  )


## Now export - have everything needed for OS and PFS curves now 
# Export to CSV
write.csv(combined_tbl, file = "combined_clinical_MRD_OS_table_updated_May2025.csv", row.names = FALSE)

# Export to RDS
saveRDS(combined_tbl, file = "combined_clinical_MRD_OS_table_updated_May2025.rds")
