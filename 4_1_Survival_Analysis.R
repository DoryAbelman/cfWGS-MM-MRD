################################################################################
##  Survival Analysis and Relapse Detection Sensitivity
##  
##  Purpose: 
##    Generate Kaplan-Meier survival curves stratified by MRD status at key
##    clinical timepoints, and calculate sensitivity of each assay for detecting
##    future relapse among frontline-treated multiple myeloma patients.
##    Includes head-to-head comparison of cfWGS, clinical assays (MFC, clonoSEQ),
##    and EasyM proteomic MRD.
##  
##  Main Analyses:
##    1. Progression-free survival (PFS) curves stratified by MRD status
##       at landmark timepoints (Post-ASCT, 1yr Maintenance, etc.)
##    2. Sensitivity calculation: % of patients who relapsed that were MRD+
##       at each timepoint, for each assay independently
##    3. Comparative sensitivity barplots (BM vs Blood-derived cfWGS models)
##    4. Optional: Time-window prediction analysis in non-frontline cohort
##  
##  Input Data Sources:
##    - all_patients_with_BM_and_blood_calls_updated5.rds 
##      (from 3_1_Optimize_cfWGS_thresholds.R)
##    - EasyM_all_samples_with_optimized_calls.csv 
##      (from 3_1_A_Process_and_optimize_EasyM.R)
##    - Censor_dates_per_patient_for_PFS_updated.rds 
##      (clinical outcome tracking)
##  
##  Output Locations:
##    - Output_tables_2025/detection_progression_updated6/
##      └─ Kaplan-Meier curves (PNG, organized by timepoint)
##      └─ Sensitivity tables (CSV)
##      └─ Comparative barplots (Supp_6A, Supp_8A)
##  
##  Scripts this Script Depends On:
##    1. 3_1_Optimize_cfWGS_thresholds.R - cfWGS model optimization/thresholds
##    2. 3_1_A_Process_and_optimize_EasyM.R - EasyM model development
##    3. 2_0_Assemble_Table_With_All_Features.R - Feature integration
##  
##  Author: Dory Abelman
##  Updated: February 2026
##  
################################################################################

## ── 0. SETUP: Load Packages and Configure Paths ──────────────────────────────
##
##  This section:
##    - Loads all required R packages for survival analysis & plotting
##    - Defines input/output file paths
##    - Creates output directories if they don't exist
##
## ────────────────────────────────────────────────────────────────────────────

# Required packages for survival analysis, visualization, and data wrangling
library(tidyverse)       # dplyr, ggplot2, tidyr
library(lubridate)       # Date/time operations
library(survival)        # Survival objects, survfit()
library(survminer)       # ggsurvplot() for KM curves
library(broom)           # Tidy model outputs
library(patchwork)       # Combine multiple plots
library(tableone)        # Create summary tables (optional)
library(timeROC)         # Time-dependent ROC analysis (optional)

## ─────────────────────────────────────────────────────────────────────────────
## INPUT FILES: Clinical outcomes and cfWGS results
## ─────────────────────────────────────────────────────────────────────────────

# File 1: Clinical metadata - censor dates, relapse status per patient
#         Used to compute time-to-event and determine event status for survival analysis
final_tbl_rds <- "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"

# File 2: Main data table with all cfWGS model calls and clinical MRD results
#         Each row = one sample (patient + timepoint)
#         Contains: MFC calls, clonoSEQ calls, cfWGS calls (multiple models)
dat_rds       <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated5.rds"

## ─────────────────────────────────────────────────────────────────────────────
## OUTPUT DIRECTORY: All results saved here
## ─────────────────────────────────────────────────────────────────────────────

outdir <- "Output_tables_2025/detection_progression_updated6"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# SOURCE DATA DIRECTORY: Stores figure source data for manuscript submission
outdir_source_data <- "Output_tables_2025/Source_data"
dir.create(outdir_source_data, showWarnings = FALSE, recursive = TRUE)

# ═════════════════════════════════════════════════════════════════════════════
# FILE VERSIONING: Prevents overwriting previous results
# ═════════════════════════════════════════════════════════════════════════════
# All output files include today's date in their names:
#   KM_assay_timepoint_updated_no_CI_2026-02-24.png
#   frontline_postASCT_sensitivity_2026-02-24.csv
# This ensures each run creates new files without overwriting previous versions.
# To compare runs across different dates, simply check the date suffix.
# ═════════════════════════════════════════════════════════════════════════════

date_tag <- format(Sys.Date(), "%Y-%m-%d")

# ═════════════════════════════════════════════════════════════════════════════
# SOURCE DATA EXPORTS: Figure source data for manuscript
# ═════════════════════════════════════════════════════════════════════════════
# All figure source data saved to: Output_tables_2025/Source_data/
#
# KEY EXPORTS (exported immediately after data frame creation):
#   • Supp_6A_BM_sensitivity_barplot_source_data_YYYY-MM-DD.csv
#     └─ Source: sens_df_bm (BM subset sensitivity by assay × timepoint)
#
#   • Supp_8A_blood_sensitivity_barplot_source_data_YYYY-MM-DD.csv
#     └─ Source: sens_df_blood (Blood subset sensitivity by assay × timepoint)
#
#   • SuppFig8B_blood_HR_plot_source_data_YYYY-MM-DD.csv
#     └─ Source: hr_plot_df_blood (Blood-derived HR by landmark × assay)
#
#   • Supp_Figure_6B_BM_HR_plot_source_data_YYYY-MM-DD.csv
#     └─ Source: hr_plot_df_bm (BM-derived HR by landmark × assay)
#
# ADDITIONAL EXPORTS (analytical summaries):
#   • frontline_followup_summary_YYYY-MM-DD.csv (see line ~1046)
#   • frontline_postASCT_sensitivity_YYYY-MM-DD.csv (see line ~1175)
#   • frontline_1yr_sensitivity_YYYY-MM-DD.csv
#   • All progression metrics CSVs (see lines ~1900, 2330, 2763)
#
# NOTE: KM curve source data = filtered survival_df for each assay/timepoint
#       (generated dynamically in loop; not exported separately)
# ═════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("═", 80), "\n")
cat("SURVIVAL ANALYSIS: Frontline Cohort\n")
cat("Output directory:", outdir, "\n")
cat("Date tag for versioning:", date_tag, "\n")
cat(strrep("═", 80), "\n\n")

## ── 1. LOAD AND TIDY CORE TABLES ─────────────────────────────────────────────
##
##  This section:
##    - Loads clinical outcome metadata (follow-up dates, relapse status)
##    - Loads main cfWGS data table with all model predictions
##    - Filters to frontline cohort only (primary analysis population)
##    - Creates additional computed columns as needed
##
## ────────────────────────────────────────────────────────────────────────────

cat("1. Loading and preparing clinical data tables...\n")

# CLINICAL OUTCOMES TABLE
# Standardizes column names to lowercase and selects key columns:
#   - Patient: unique patient identifier
#   - baseline_date: treatment start date
#   - censor_date: last known follow-up date (event or censoring)
#   - relapsed: binary indicator (1=relapsed, 0=censored without relapse)
final_tbl <- readRDS(final_tbl_rds) %>%
  rename_with(
    tolower,
    any_of(c("Baseline_Date", "Censor_date", "Relapsed"))
  ) %>%
  transmute(
    Patient       = as.character(Patient),
    baseline_date = as.Date(baseline_date),
    censor_date   = as.Date(censor_date),
    relapsed      = as.integer(relapsed)
  )

cat(sprintf("  ✓ Clinical outcomes: %d patients loaded\n", n_distinct(final_tbl$Patient)))

# MAIN cfWGS + CLINICAL ASSAYS DATA TABLE
# Each row represents one sample (patient + timepoint)
# Contains predictions from all MRD models and timepoint classifications
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)  # standardize timepoint names for consistent grouping
  )

cat(sprintf("  ✓ Sample data: %d samples from %d patients loaded\n", 
            nrow(dat), n_distinct(dat$Patient)))

# FILTER TO FRONTLINE COHORT
# Primary analysis uses only patients in frontline/induction-to-transplant cohort
dat <- dat %>% filter(Cohort == "Frontline")

cat(sprintf("  ✓ Filtered to Frontline cohort: %d samples from %d patients\n", 
            nrow(dat), n_distinct(dat$Patient)))

# CREATE ADDITIONAL BM cfWGS MODELS (alternative decision rules)
# Screen call: different threshold optimized for sensitivity (detect 95% of relapsers)
#   Uses probability >= 0.350 instead of default optimized threshold
dat <- dat %>%
  mutate(
    BM_zscore_only_detection_rate_screen_call = as.integer(BM_zscore_only_detection_rate_prob >= 0.350),
  )

cat("  ✓ Data preparation complete\n\n")

## ── 1A. LOAD EasyM PROTEOMIC MRD DATA ─────────────────────────────────────────
##
##  This section:
##    - Loads EasyM M-protein measurements and optimized binary calls
##    - Merges EasyM data by patient and timepoint
##    - Gracefully handles missing EasyM data with placeholder columns
##
##  EasyM Data Sources (from script 3_1_A):
##    - EasyM_value: continuous M-protein measure (%)
##    - EasyM_optimized_binary: binary call (1=positive, 0=negative) using
##                            optimized thresholds from script 3_1_A
##    - EasyM_optimized_call: character label of call for reporting
##
## ────────────────────────────────────────────────────────────────────────────

cat("1A. Loading EasyM proteomic MRD data...\n")

OUTPUT_DIR_EASYM <- "Output_EasyM_MRD_analysis_2025"
EasyM_file <- file.path(OUTPUT_DIR_EASYM, "EasyM_all_samples_with_optimized_calls.csv")

if (file.exists(EasyM_file)) {
  # Load EasyM predictions from script 3_1_A output
  EasyM_data <- readr::read_csv(EasyM_file, show_col_types = FALSE)
  cat(sprintf("  ✓ Loaded EasyM data: %d samples\n", nrow(EasyM_data)))
  
  # Merge EasyM data into main dataset
  # Join key: Patient (character) + Timepoint (character)
  # relationship = "many-to-one": multiple samples per patient, but one EasyM value per timepoint
  dat <- dat %>%
    mutate(Patient = as.character(Patient), Timepoint = as.character(Timepoint)) %>%
    left_join(
      EasyM_data %>%
        select(Patient, Timepoint, EasyM_value, EasyM_optimized_binary, EasyM_optimized_call) %>%
        mutate(Patient = as.character(Patient), Timepoint = as.character(Timepoint)),
      by = c("Patient", "Timepoint"),
      relationship = "many-to-one"
    )
  
  n_easym_matched <- sum(!is.na(dat$EasyM_optimized_binary))
  cat(sprintf("  ✓ Merged EasyM data: %d samples with EasyM calls\n", n_easym_matched))
  
} else {
  # If EasyM file not found, create placeholder columns (survival analysis will skip EasyM with NA filter)
  cat(sprintf("  ⚠ EasyM file not found at: %s\n", EasyM_file))
  cat("    Creating placeholder columns (EasyM will be skipped in downstream analyses)\n")
  
  dat <- dat %>%
    mutate(
      EasyM_value = NA_real_,
      EasyM_optimized_binary = NA_integer_,
      EasyM_optimized_call = NA_character_
    )
}

cat("  ✓ EasyM data preparation complete\n\n")

## ── 2. BUILD SURVIVAL DATA TABLE (sample-level) ───────────────────────────────
##
##  This section:
##    - Merges clinical outcomes (relapse, follow-up dates) with sample-level data
##    - Computes time-to-event: days from sample draw to relapse or censoring
##    - Creates binary event indicator (1=relapsed, 0=censored)
##    - Selects only columns needed for downstream survival analysis
##
##  Key computations:
##    - Time_to_event = censor_date - sample_date 
##      (how long from THIS sample until event/censoring)
##    - Relapsed_Binary = 1 if patient relapsed (regardless of when), 0 if censored
##
##  Output: survival_df
##    - One row per sample (patient × timepoint)
##    - Ready for stratified KM analysis by MRD status
##
## ────────────────────────────────────────────────────────────────────────────

cat("2. Building survival analysis table...\n")

survival_df <- dat %>%
  mutate(
    Patient     = as.character(Patient),
    sample_date = as.Date(Date)  # standardize date column name and type
  ) %>%
  # Merge in patient-level outcomes (censor_date, relapsed status)
  left_join(
    final_tbl %>% 
      select(Patient, censor_date, relapsed),
    by = "Patient"
  ) %>%
  # Calculate time from sample draw to event/censoring
  mutate(
    Time_to_event   = as.numeric(censor_date - sample_date),  # days
    Relapsed_Binary = as.integer(relapsed)                     # 1=event, 0=censored
  ) %>%
  # Select only columns needed for survival analyses
  # (This reduces memory footprint and makes code more readable)
  select(
    Patient, Timepoint, sample_date, censor_date, timepoint_info,
    Time_to_event, Relapsed_Binary,
    # Clinical MRD assays
    Flow_Binary, Adaptive_Binary, Rapid_Novor_Binary,
    Flow_pct_cells, Adaptive_Frequency,
    # PET imaging (if needed)
    PET_Binary,
    # cfWGS BM-derived models
    BM_zscore_only_detection_rate_call, BM_zscore_only_detection_rate_prob, 
    BM_zscore_only_detection_rate_screen_call,
    # cfWGS blood-derived models (multiple variants)
    Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob, 
    Blood_base_prob, Blood_base_call,
    Blood_plus_fragment_prob, Blood_plus_fragment_call,
    Blood_plus_fragment_min_prob, Blood_plus_fragment_min_call, 
    # Fragmentomics models
    Fragmentomics_mean_coverage_only_prob, Fragmentomics_mean_coverage_only_call,
    # EasyM proteomic MRD
    EasyM_optimized_binary, EasyM_value
  )

cat(sprintf("  ✓ Survival table created: %d samples from %d patients\n", 
            nrow(survival_df), n_distinct(survival_df$Patient)))

# QUICK DATA VALIDATION
# Verify the survival data looks reasonable
cat("\n  Data validation checks:\n")
cat(sprintf("    - Follow-up time: %.1f to %.1f days (median: %.1f days)\n",
            min(survival_df$Time_to_event, na.rm=TRUE),
            max(survival_df$Time_to_event, na.rm=TRUE),
            median(survival_df$Time_to_event, na.rm=TRUE)))

relapse_counts <- table(survival_df$Relapsed_Binary, useNA="ifany")
cat(sprintf("    - Relapse events: %d / %d (%.1f%%)\n",
            relapse_counts["1"], sum(!is.na(survival_df$Relapsed_Binary)),
            relapse_counts["1"] / sum(!is.na(survival_df$Relapsed_Binary)) * 100))

cat(sprintf("    - Timepoints represented: %s\n",
            paste(unique(survival_df$timepoint_info), collapse=", ")))

cat("  ✓ Survival table validation complete\n\n")

## ── 3. KAPLAN-MEIER SURVIVAL ANALYSIS CONFIGURATION ──────────────────────────
##
##  This section configures settings for survival curve generation:
##    - Define which assays/models to include in KM analyses
##    - Specify timepoints to analyze
##    - Set visualization parameters (colors, DPI, minimum group size)
##    - Map human-readable labels for each timepoint
##
## ────────────────────────────────────────────────────────────────────────────

cat("3. Configuring Kaplan-Meier survival analyses...\n")

# LIST OF MRD TECHNOLOGIES/MODELS TO ANALYZE
# Each entry: column_name = "Display Label for Plots"
# Includes: clinical assays (MFC, clonoSEQ), cfWGS models, and EasyM
techs <- c(
  # Clinical MRD assays
  Flow_Binary        = "MFC",
  Adaptive_Binary    = "clonoSEQ",
  EasyM_optimized_binary = "EasyM (Proteomic MRD)",
  # Bone Marrow-derived WGS mutations
  BM_zscore_only_detection_rate_call    = "cfWGS of BM-Derived Mutations (cVAF Model)", 
  BM_zscore_only_detection_rate_screen_call = "cfWGS of BM-derived mutations (High Sensitivity)", 
  # Blood Plasma-derived WGS mutations
  Blood_zscore_only_sites_call = "cfWGS of cfDNA-Derived Mutations (Sites Model)",
  Blood_plus_fragment_call = "cfWGS of cfDNA-Derived Mutations (Combined Model)"
  # NOTE: Fragmentomics models could be added here if desired
)

cat(sprintf("  - Analyzing %d MRD technologies/models\n", length(techs)))
cat(sprintf("    %s\n", paste("   ", names(techs), sep=" ✓ ", collapse="\n    "))[1:3])
cat("    ...\n")

# CLINICAL TIMEPOINTS TO ANALYZE
# Extract all unique timepoints from the data
# These represent key clinical decision points (diagnosis, post-transplant, maintenance, etc.)
tps <- unique(survival_df$timepoint_info)
cat(sprintf("  - Timepoints to analyze: %s\n", paste(tps, collapse=", ")))

# VISUALIZATION PARAMETERS
dpi_target <- 500  # Resolution for PNG output of KM curves

# GROUP SIZE THRESHOLDS
# Skip KM curve if < this many patients in either MRD+/MRD- group
# (Small sample sizes lead to unreliable survival estimates)
min_n <- 5

# COLOR PALETTE FOR KM CURVES
# MRD negative = black, MRD positive = red
pal_2 <- c("black", "red")

# HUMAN-READABLE LABELS FOR TIMEPOINTS
# Maps internal names (as stored in timepoint_info column) to plot labels
tp_labels <- c(
  `diagnosis`          = "Diagnosis",
  `post_transplant`    = "Post‑ASCT",
  `1yr maintenance`    = "One-Year Maintenance", 
  `post_induction`     = "Post‑Induction",
  `post_asct`          = "Post‑ASCT",
  `maintenance`        = "Maintenance",
  `1yr maint`          = "One-Year Maintenance"
)

# PREPARE BASELINE DIAGNOSIS DATES
# Compute earliest sample date per patient (used for relative time calculations if needed)
dx_tbl <- survival_df %>%
  group_by(Patient) %>%
  summarise(
    diagnosis_date = suppressWarnings(min(sample_date[timepoint_info == "diagnosis"], na.rm = TRUE)),
    .groups = "drop"
  )

cat("  ✓ Configuration complete\n\n")


## ── 4. GENERATE KAPLAN-MEIER CURVES: Timepoint × Technology Analysis ────────
##
##  This section contains nested loops that:
##    1. Iterate through each clinical timepoint
##    2. For each timepoint, iterate through each MRD technology
##    3. For each combination, generate a stratified KM curve:
##       - Subjects split into MRD+ vs MRD- groups based on assay result
##       - Curve shows PFS probability over follow-up time
##       - Risk table shows number at risk per group at each time
##    4. Saves PNG at 500 dpi for manuscript inclusion
##
##  Nested Loop Structure:
##    for (timepoint in all timepoints)
##      for (assay in all technologies)
##        Generate KM curve for that timepoint + assay combo
##        Save to: outdir/timepoint/KM_assay_timepoint.png
##
##  Output Directory Organization:
##    detection_progression_updated6/
##    ├── diagnosis/
##    │   ├── KM_MFC_Diagnosis_updated_no_CI.png
##    │   ├── KM_clonoSEQ_Diagnosis_updated_no_CI.png
##    │   ├── KM_EasyM_Diagnosis_updated_no_CI.png
##    │   └── ...
##    ├── post_transplant/
##    │   ├── KM_MFC_Post‑ASCT_updated_no_CI.png
##    │   ├── KM_clonoSEQ_Post‑ASCT_updated_no_CI.png
##    │   ├── KM_EasyM_Post‑ASCT_updated_no_CI.png
##    │   ├── KM_cfWGS of BM-Derived Mutations (cVAF Model)_Post‑ASCT_updated_no_CI.png
##    │   └── ...
##    └── 1yr_maintenance/
##        └── ...
##
##  Notes on Curve Generation:
##    - Skips assays with < min_n (5) patients in either group
##    - Uses Kaplan-Meier non-parametric estimator
##    - Includes log-rank p-value testing MRD+ vs MRD- groups
##    - Risk table shows N at risk below x-axis at selected timepoints
##    - X-axis = months from MRD assessment, Y-axis = PFS probability
##
## ────────────────────────────────────────────────────────────────────────────

cat("4. Generating Kaplan-Meier survival curves...\n")
cat("   (This may take a minute - creating curves for all timepoint×assay combinations)\n\n")

# Loop 1: Iterate through each timepoint
for(tp in tps) {
  # Get nice label for this timepoint
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  # Create subdirectory for this timepoint's curves
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Loop 2: Iterate through each MRD technology/model
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_no_CI_", date_tag, ".png"))
    
    # ─────────────────────────────────────────────────────────────────────────
    # PREPARE DATA FOR THIS TIMEPOINT × ASSAY
    # Filter criteria:
    #   - timepoint_info == tp: Only samples from this timepoint
    #   - !is.na(Time_to_event): Remove if follow-up time missing
    #   - !is.na(Relapsed_Binary): Remove if relapse status unknown
    #   - !is.na(.data[[var]]): Remove if MRD assay result missing
    # arrange() + group_by() + slice(1) ensures we keep only the FIRST sample per patient
    # at that timepoint (in case multiple draws on same date)
    # ─────────────────────────────────────────────────────────────────────────
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%           # keep just the first draw per patient if multiple at that timepoint
      ungroup() %>%
      # Create stratification group: MRD+ (assay value = 1) vs MRD- (assay value = 0)
      # This is the key predictor variable for the KM curves
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    # Convert time-to-event from days to months (multiply by 12, divide by 365)
    # Using 30.44 = average days per month (365.25 / 12) for accuracy
    df_sub <- df_sub %>% 
      mutate(Time_to_event = Time_to_event/30.44) # divide days by 30.44 to get months
    
    # ─────────────────────────────────────────────────────────────────────────
    # VALIDATION: Skip assays with insufficient data
    # ─────────────────────────────────────────────────────────────────────────
    
    # Skip if total sample size is too small (min_n = 5)
    if(nrow(df_sub) < min_n) next
    
    # Skip if we don't have both MRD+ and MRD- groups (need both for comparison)
    if(n_distinct(df_sub$Group) < 2) next
    
    # ─────────────────────────────────────────────────────────────────────────
    # FIT KAPLAN-MEIER MODEL AND GENERATE SURVIVAL PLOT
    # ─────────────────────────────────────────────────────────────────────────
    
    # Surv() creates survival object with time and event indicator
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    
    # survfit() fits KM curves stratified by MRD status (one curve per group)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    # ggsurvplot generates publication-quality KM plot with:
    #   - Log-rank p-value comparing groups
    #   - Risk table showing number of subjects at risk over time
    #   - Customized colors, labels, and formatting
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = TRUE,
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = FALSE,
      risk.table      = TRUE,
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),  # ← left‑align
      palette         = pal_2,
      # legend.title    = paste0(assay_lab, " MRD"),
      legend.title    = "MRD status",
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Time since MRD assessment (months)",
      ylab            = "Progression-free survival",
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position = "top",
          axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),          # no grid
          panel.grid.minor = element_blank(),
          #  Make the tick‑labels (the numbers) larger:
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y      = element_text(size = 15),
          axis.title.x      = element_text(size = 14)
        ),
    )
    
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title      = element_text(hjust = 0, face = "plain"),
      )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}


## Redo with CI 
# 5) Loop
for(tp in tps) {
  #  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_with_CI.png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%           # keep just the first draw per patient if multiple at that timepoint
      ungroup() %>%
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    df_sub <- df_sub %>% 
      mutate(Time_to_event = Time_to_event/30.44) # divide days by 30.44 to get months
    
    # skip if too few pts
    if(nrow(df_sub) < min_n) next
    
    # skip if only one group present
    if(n_distinct(df_sub$Group) < 2) next
    
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = TRUE,
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = TRUE,
      risk.table      = TRUE,
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),  # ← left‑align
      palette         = pal_2,
      # legend.title    = paste0(assay_lab, " MRD"),
      legend.title    = "MRD status",
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Time since MRD assessment (months)",
      ylab            = "Progression-free survival",
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position = "top",
          axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),          # no grid
          panel.grid.minor = element_blank(),
          #  Make the tick‑labels (the numbers) larger:
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y      = element_text(size = 15),
          axis.title.x      = element_text(size = 14)
        ),
    )
    
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title      = element_text(hjust = 0, face = "plain"),
      )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}


### Set the CI instead to be to 90% 
for(tp in tps) {
  #  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_with_CI_90.png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%           # keep just the first draw per patient if multiple at that timepoint
      ungroup() %>%
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    df_sub <- df_sub %>% 
      mutate(Time_to_event = Time_to_event/30.44) # divide days by 30.44 to get months
    
    # skip if too few pts
    if(nrow(df_sub) < min_n) next
    
    # skip if only one group present
    if(n_distinct(df_sub$Group) < 2) next
    
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    # 90% CI from survfit; "log" (Greenwood on log scale) is common, "log-log" is also fine
    fit <- survfit(
      surv_obj ~ Group,
      data      = df_sub,
      conf.int  = 0.90     # ← 90% CI
    )
    
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = TRUE,
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = TRUE,
      conf.int.alpha  = 0.1,    
      risk.table      = TRUE,
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),  # ← left‑align
      palette         = pal_2,
      # legend.title    = paste0(assay_lab, " MRD"),
      legend.title    = "MRD status",
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Time since MRD assessment (months)",
      ylab            = "Progression-free survival",
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position = "top",
          axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),          # no grid
          panel.grid.minor = element_blank(),
          #  Make the tick‑labels (the numbers) larger:
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y      = element_text(size = 15),
          axis.title.x      = element_text(size = 14)
        ),
    )
    
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title      = element_text(hjust = 0, face = "plain"),
      )
    
    # # Convert survfit object to a data.frame
    # fit_df <- broom::tidy(fit)  # gives time, n.risk, n.event, surv, std.err, conf.low, conf.high, strata
    # 
    # # Overlay ribbons for CIs
    # km$plot <- km$plot +
    #   geom_ribbon(
    #     data = fit_df,
    #     aes(x = time, ymin = conf.low, ymax = conf.high, fill = strata, group = strata),
    #     inherit.aes = FALSE, alpha = 0.2
    #   )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}

 



#### Redo but show the plot starting from diagnosis
### To complex since the default risk table doesn't support delayed entry - skipped
### delayed entry (a.k.a. left truncation) to avoid immortal-bias
# pretty p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-4) return("p < 1e-4")
  if (p < 1e-3) return("p < 0.001")
  if (p < 1e-2) return("p < 0.01")
  sprintf("p = %.2f", p)
}

# risk table for delayed-entry (start–stop) data
make_risktable <- function(df, breaks) {
  tmp <- df %>%
    mutate(Group_label = dplyr::recode(as.character(Group),
                                       "Negative" = "MRD–",
                                       "Positive" = "MRD+"))
  purrr::map_dfr(breaks, function(ti) {
    tmp %>%
      dplyr::group_by(Group_label) %>%
      dplyr::summarise(n = sum(entry_m <= ti & exit_m > ti), .groups = "drop") %>%
      dplyr::mutate(time = ti)
  }) %>%
    tidyr::pivot_wider(names_from = time, values_from = n) %>%
    dplyr::arrange(factor(Group_label, levels = c("MRD–","MRD+"))) %>%
    dplyr::rename(`MRD status` = Group_label)
}


for (tp in tps) {
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_from_diagnosis.png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%
      ungroup() %>%
      left_join(dx_tbl, by = "Patient") %>%
      mutate(
        Group   = factor(ifelse(.data[[var]] == 1, "Positive", "Negative"),
                         levels = c("Negative","Positive")),
        # entry = months from diagnosis to MRD test
        entry_m = as.numeric(sample_date - diagnosis_date) / 30.44,
        # exit = entry + observed time after MRD test
        exit_m  = entry_m + (Time_to_event / 30.44)
      ) %>%
      filter(!is.na(entry_m), !is.na(exit_m), exit_m >= entry_m)
    
    if (nrow(df_sub) < min_n) next
    if (n_distinct(df_sub$Group) < 2) next
    
    # Cox with delayed entry
    surv_obj <- Surv(time = df_sub$entry_m, time2 = df_sub$exit_m, event = df_sub$Relapsed_Binary)
    fit     <- survfit(surv_obj ~ Group, data = df_sub)
    
    cox_fit <- coxph(surv_obj ~ Group, data = df_sub, ties = "breslow")
    s <- summary(cox_fit)
    
    # Pull p-values (log-rank / score / wald) safely
    p_lrt   <- suppressWarnings(as.numeric(s$logtest[3]))  # likelihood-ratio p
    p_score <- suppressWarnings(as.numeric(s$sctest[3]))   # score test p
    p_wald  <- suppressWarnings(as.numeric(s$waldtest[3])) # Wald p
    
    # Fallback if all above are NA (e.g., separation): grab coefficient p if present
    coef_p <- NA_real_
    cs <- try(coef(summary(cox_fit)), silent = TRUE)
    if (!inherits(cs, "try-error")) {
      # find the row for Group (works even if level names change)
      r <- grep("^Group", rownames(cs), value = FALSE)
      if (length(r) >= 1) coef_p <- cs[r[1], "Pr(>|z|)"]
    }
    
    # Choose the best available p
    p_any <- if (!is.na(p_lrt)) p_lrt else if (!is.na(p_score)) p_score else if (!is.na(p_wald)) p_wald else coef_p
    
    # Format: show thresholds instead of tiny decimals
    fmt_p <- function(p) {
      if (is.na(p)) return("p = NA")
      if (p < 1e-4) return("p < 1e-4")
      if (p < 1e-3) return("p < 0.001")
      if (p < 1e-2) return("p < 0.01")
      sprintf("p = %.2f", p)
    }
    pval_str <- fmt_p(p_any)
    
    km <- ggsurvplot(
      fit, data = df_sub,
      pval              = pval_str,
      break.time.by     = 12,
      conf.int          = TRUE,
      risk.table        = TRUE,
      risk.table.title  = "Number at risk",                    # ← add
      risk.table.title.theme = element_text(hjust = 0),        # ← add (left-align)
      palette           = pal_2,
      legend.title      = "MRD status",
      legend.labs       = c("MRD-","MRD+"),  # use ASCII hyphen to avoid file/device issues
      xlab              = "Months since diagnosis",
      ylab              = "Progression-free survival",
      title             = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25,
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position  = "top",
          axis.line        = element_line(colour = "black"),   # ← add
          panel.grid.major = element_blank(),                  # ← add
          panel.grid.minor = element_blank(),                  # ← add
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y     = element_text(size = 15),
          axis.title.x     = element_text(size = 14)
        )
    )
    
    # match the old post-processing of table/plot
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title   = element_text(hjust = 0, face = "plain")
      )
    
    km$plot <- km$plot +
      theme(axis.title.x = element_blank())
    
    combined <- ggarrange(km$plot, km$table, ncol = 1, heights = c(3,1))
    ggsave(fname, plot = combined, width = 7, height = 7, dpi = dpi_target)
  }
}




## ── 5. SENSITIVITY ANALYSIS: Frontline Cohort MRD Performance ───────────────
##
##  This section calculates sensitivity and specificity metrics for each MRD assay
##  in the Frontline cohort at two critical clinical timepoints:
##    1. Post-ASCT (post-autologous stem cell transplant)
##    2. Maintenance-1yr (1-year maintenance therapy)
##
##  For each assay and timepoint:
##    - Filter to Frontline cohort only
##    - Stratify by MRD status (positive/negative)
##    - Calculate sensitivity = % of relapsed with MRD+
##    - Calculate specificity = % of non-relapsed with MRD-
##    - Calculate PPV/NPV = predictive values
##
##  Output: Summary table with MRD performance metrics
##
## ────────────────────────────────────────────────────────────────────────────

cat("5. Calculating MRD sensitivity/specificity metrics (Frontline cohort)...\n\n")

#### Now get stats for results 
### First on frontline 
# 1) Define frontline cohort for analysis
front_patients <- dat %>%
  filter(Cohort == "Frontline") %>%
  distinct(Patient)

# ─────────────────────────────────────────────────────────────────────────────
# 2) COMPUTE FOLLOW-UP TIME AND RELAPSE STATUS FOR FRONTLINE COHORT
# Time = number of days from baseline date (diagnosis) to censor/relapse date
# This is the time-to-event used for KM curves and outcome calculations
# ─────────────────────────────────────────────────────────────────────────────

pfs_front <- final_tbl %>%
  filter(Patient %in% front_patients$Patient) %>%
  # compute time from baseline to censor/relapse
  mutate(time_days = as.numeric(censor_date - baseline_date)) 

## check one row per patient (sanity check to ensure data integrity)
pfs_front %>% 
  count(Patient) %>% 
  filter(n > 1) -> dups
if(nrow(dups)) stop("Duplicate patients found: ", paste(dups$Patient, collapse = ", "))

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY STATISTICS: Frontline cohort follow-up and relapse
# ─────────────────────────────────────────────────────────────────────────────

median_fu_mo <- median(pfs_front$time_days / 30.44, na.rm = TRUE)
n_front      <- nrow(pfs_front)
n_rel        <- sum(pfs_front$relapsed)
pct_rel      <- n_rel / n_front * 100

cat(sprintf(
  "Frontline cohort summary:\n  Median follow-up: %.1f months\n  Relapses: %d/%d (%.0f%%) patients\n\n",
  median_fu_mo, n_rel, n_front, pct_rel
))

# ─────────────────────────────────────────────────────────────────────────────
# 3) DETAILED FOLLOW-UP STATISTICS
# Comprehensive summary of follow-up times in days and months
# Includes: min/Q1/median/Q3/max/mean/sd for all subjects
# ─────────────────────────────────────────────────────────────────────────────

followup_stats <- pfs_front %>%
  summarise(
    N_patients      = dplyr::n(),
    N_relapses      = sum(relapsed),
    Relapse_rate    = N_relapses / N_patients * 100,
    
    # Follow-up time in days
    min_days        = min(time_days, na.rm = TRUE),
    q1_days         = quantile(time_days, 0.25, na.rm = TRUE),
    median_days     = median(time_days, na.rm = TRUE),
    q3_days         = quantile(time_days, 0.75, na.rm = TRUE),
    max_days        = max(time_days, na.rm = TRUE),
    mean_days       = mean(time_days, na.rm = TRUE),
    sd_days         = sd(time_days, na.rm = TRUE),
    
    # Follow-up time in months (for clinical interpretation)
    min_months      = min(time_days, na.rm = TRUE) / 30.44,
    q1_months       = quantile(time_days / 30.44, 0.25, na.rm = TRUE),
    median_months   = median(time_days / 30.44, na.rm = TRUE),
    q3_months       = quantile(time_days / 30.44, 0.75, na.rm = TRUE),
    max_months      = max(time_days, na.rm = TRUE) / 30.44,
    mean_months     = mean(time_days, na.rm = TRUE) / 30.44,
    sd_months       = sd(time_days, na.rm = TRUE) / 30.44
  )

cat("Follow-up time summary (Frontline cohort):\n")
print(followup_stats)
cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4) OPTIONAL: Visualization of follow-up time distribution
# Histogram showing how follow-up times are distributed across cohort
# ─────────────────────────────────────────────────────────────────────────────

followup_hist <- pfs_front %>%
  mutate(followup_months = time_days / 30.44) %>%
  ggplot(aes(x = followup_months)) +
  geom_histogram(binwidth = 3, boundary = 0) +
  labs(
    x = "Follow-up time (months)",
    y = "Number of patients",
    title = "Distribution of follow-up times in frontline cohort"
  )

# If you want to export the stats to CSV
write_csv(followup_stats, file.path(outdir, paste0("frontline_followup_summary_", date_tag, ".csv")))
## Can use 1A for this as well 

# 4.  Assays & timepoint definitions ------------------------------------------
assays <- c(
  EasyM    = "EasyM_optimized_binary",
  clonoSEQ = "Adaptive_Binary",
  Flow     = "Flow_Binary",
  cfWGS_BM    = "BM_zscore_only_detection_rate_call",
  cfWGS_BM_screen    = "BM_zscore_only_detection_rate_screen_call",
  cfWGS_Blood_Sites    = "Blood_zscore_only_sites_call",
  cfWGS_Blood_Combined = "Blood_plus_fragment_call"
)

post_labels   <- c("post_transplant")
one_year_labels <- c("1yr maintenance")

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTION: compute_sens()
# Calculates sensitivity for a given MRD assay
#
# Sensitivity = % of RELAPSED patients who were MRD+ at baseline timepoint
# Formula: N_positive / N_tested (where N_tested = samples with assay result)
#
# Input:
#   - df: data frame with relapsed patients and assay columns
#   - col: column name of assay (e.g., "Adaptive_Binary", "EasyM_optimized_binary")
# Output:
#   - tibble with N_tested, N_positive, and computed Sensitivity
# ─────────────────────────────────────────────────────────────────────────────

compute_sens <- function(df, col) {
  df2      <- df %>% filter(!is.na(.data[[col]]))  # Remove rows with missing assay result
  n_tested <- nrow(df2)                            # Total relapsed patients with assay result
  n_pos    <- sum(df2[[col]] == 1, na.rm = TRUE)   # Number of those who were MRD+
  tibble(
    N_tested      = n_tested,
    N_positive    = n_pos,
    Sensitivity   = n_pos / n_tested              # Sensitivity = % MRD+ among relapsed
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. POST-ASCT SENSITIVITY CALCULATIONS
# Sensitivity = % of relapsed patients who had detectable MRD at post-ASCT
# This answers: "Which assays best predicted relapse at post-ASCT?"
# ─────────────────────────────────────────────────────────────────────────────

# Filter to:
#   - Frontline cohort only
#   - timepoint_info contains "post_transplant" (post-ASCT samples)
#   - Take first (earliest) sample per patient
# Then restrict to relapsed patients only for sensitivity calculation

post_df <- dat %>%
  filter(
    Patient        %in% front_patients$Patient,
    str_detect(timepoint_info, paste(post_labels, collapse = "|"))
  ) %>%
  arrange(Patient, sample_date) %>%
  group_by(Patient) %>%
  slice(1) %>%   # earliest post-ASCT sample per patient
  ungroup() %>%
  select(Patient, one_of(assays)) %>%
  left_join(final_tbl %>% select(Patient, relapsed), by = "Patient") %>%
  filter(relapsed == 1)  # Keep only relapsed patients for sensitivity calc

# Calculate sensitivity for each assay
post_stats <- map_dfr(names(assays), ~ {
  col <- assays[.x]
  stats <- compute_sens(post_df, col)
  stats %>% mutate(Assay = .x)
}, .id = NULL) %>%
  select(Assay, everything())

# ─────────────────────────────────────────────────────────────────────────────
# 6. ONE-YEAR MAINTENANCE SENSITIVITY CALCULATIONS
# Sensitivity = % of relapsed patients who had detectable MRD at 1-year maintenance
# This answers: "Which assays best predicted relapse at 1-year maintenance?"
# ─────────────────────────────────────────────────────────────────────────────

# Filter to:
#   - Frontline cohort only
#   - timepoint_info contains "1yr maintenance" (1-year maintenance samples)
#   - Take first (earliest) sample per patient at that timepoint
# Then restrict to relapsed patients only for sensitivity calculation

year_df <- dat %>%
  filter(
    Patient        %in% front_patients$Patient,
    str_detect(timepoint_info, paste(one_year_labels, collapse = "|"))
  ) %>%
  arrange(Patient, sample_date) %>%
  group_by(Patient) %>%
  slice(1) %>%   # earliest 1-year maintenance sample per patient
  ungroup() %>%
  select(Patient, one_of(assays)) %>%
  left_join(final_tbl %>% select(Patient, relapsed), by = "Patient") %>%
  filter(relapsed == 1)  # Keep only relapsed patients for sensitivity calc

# Calculate sensitivity for each assay at 1-year maintenance
year_stats <- map_dfr(names(assays), ~ {
  col <- assays[.x]
  stats <- compute_sens(year_df, col)
  stats %>% mutate(Assay = .x)
}, .id = NULL) %>%
  select(Assay, everything())

# ─────────────────────────────────────────────────────────────────────────────
# 7. PRINT SENSITIVITY RESULTS TO CONSOLE
# ─────────────────────────────────────────────────────────────────────────────

cat("\n")
cat("═════════════════════════════════════════════════════════════════════════\n")
cat("POST-ASCT SENSITIVITY: % of relapsed patients with MRD+ at post-ASCT\n")
cat("═════════════════════════════════════════════════════════════════════════\n")
print(post_stats)

cat("\n")
cat("═════════════════════════════════════════════════════════════════════════\n")
cat("1-YEAR MAINTENANCE SENSITIVITY: % of relapsed patients with MRD+ at 1yr\n")
cat("═════════════════════════════════════════════════════════════════════════\n")
print(year_stats)

# ─────────────────────────────────────────────────────────────────────────────
# 8. SAVE SENSITIVITY RESULTS TO CSV
# ─────────────────────────────────────────────────────────────────────────────

cat("\nSaving sensitivity results...\n")
write_csv(
  post_stats,
  file.path(outdir, paste0("frontline_postASCT_sensitivity_", date_tag, ".csv"))
)
write_csv(
  year_stats,
  file.path(outdir, paste0("frontline_1yr_sensitivity_", date_tag, ".csv"))
)


#### Now do seperately only amongst patients who got a cfWGS test - seperate for BM and blood 
# -- after you’ve built `post_df` and `year_df` as before… --------------

# 5a.  Head-to-head in the BM-cfWGS subset (only patients with BM Z-score) ---

# define the BM-cfWGS column name
bm_col <- assays["cfWGS_BM"]

post_df_BM <- post_df %>% filter(!is.na(.data[[bm_col]]))
year_df_BM <- year_df %>% filter(!is.na(.data[[bm_col]]))

post_stats_BM <- map_dfr(names(assays), function(a) {
  compute_sens(post_df_BM, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

year_stats_BM <- map_dfr(names(assays), function(a) {
  compute_sens(year_df_BM, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

message("Post-ASCT sensitivities among those with BM-cfWGS:")
print(post_stats_BM)

message("1-yr sensitivities among those with BM-cfWGS:")
print(year_stats_BM)

# 5b.  Head-to-head in the blood-cfWGS subset (only patients with blood Z-score) -

# Define the blood-derived cfWGS column name
blood_col <- assays["cfWGS_Blood_Sites"]

# Filter to subjects with blood-cfWGS data
post_df_blood <- post_df %>% filter(!is.na(.data[[blood_col]]))
year_df_blood <- year_df %>% filter(!is.na(.data[[blood_col]]))

# Calculate sensitivity for each assay within blood subset
post_stats_blood <- map_dfr(names(assays), function(a) {
  compute_sens(post_df_blood, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

year_stats_blood <- map_dfr(names(assays), function(a) {
  compute_sens(year_df_blood, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

cat("BLOOD-SUBSET: Post-ASCT sensitivities (among those with blood-cfWGS):\n")
print(post_stats_blood)

cat("\nBLOOD-SUBSET: 1-yr sensitivities (among those with blood-cfWGS):\n")
print(year_stats_blood)

# ─────────────────────────────────────────────────────────────────────────────
# 11. SAVE SUBSET SENSITIVITY RESULTS TO CSV
# Saves BM and blood subsets separately for comparison
# ─────────────────────────────────────────────────────────────────────────────

write_csv(post_stats_BM,    file.path(outdir, paste0("frontline_postASCT_sens_BMcfWGS_", date_tag, ".csv")))
write_csv(year_stats_BM,    file.path(outdir, paste0("frontline_1yr_sens_BMcfWGS_", date_tag, ".csv")))
write_csv(post_stats_blood, file.path(outdir, paste0("frontline_postASCT_sens_bloodcfWGS_", date_tag, ".csv")))
write_csv(year_stats_blood, file.path(outdir, paste0("frontline_1yr_sens_bloodcfWGS_", date_tag, ".csv")))



# ─────────────────────────────────────────────────────────────────────────────
# 12. SENSITIVITY BARPLOT FOR MANUSCRIPT SUPPLEMENT
# Comparative visualization of sensitivity across assays and timepoints
# ─────────────────────────────────────────────────────────────────────────────

# Combine BM-subset results (post-ASCT and maintenance) into long format
# Filter to exclude blood-only cfWGS columns
# Recodes assay names to short manuscript labels
# Enforces desired display order: clinical assays → cfWGS

sens_df_bm <- bind_rows(
  post_stats_BM  %>% mutate(Timepoint = "Post-ASCT"),
  year_stats_BM  %>% mutate(Timepoint = "Maintenance-1yr")
) %>%
  # Exclude blood-derived assays and screening variant from BM-subset comparison
  filter(Assay != "cfWGS_Blood_Sites") %>%
  filter(Assay != "cfWGS_Blood_Combined") %>%
  filter(Assay != "cfWGS_BM_screen") %>%
  mutate(
    # Convert to percentage for labeling
    Sens_pct   = Sensitivity * 100,
    # Manuscript-friendly assay names
    Assay      = recode(Assay,
                        EasyM          = "EasyM (Opt)",
                        clonoSEQ       = "clonoSEQ",
                        Flow           = "MFC",
                        cfWGS_BM       = "cfWGS"),
    # Enforce display order: cfWGS, clonoSEQ, MFC, then EasyM (Opt) rightmost
    Assay = factor(Assay, levels = c("cfWGS", "clonoSEQ", "MFC", "EasyM (Opt)")),
    # Enforce timepoint order for legend and grouping
    Timepoint = factor(Timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: Supp_6A (BM sensitivity barplot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  sens_df_bm,
  file.path(outdir_source_data, paste0("Supp_6A_BM_sensitivity_barplot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: Supp_6A (BM sensitivity)\n")

# ─────────────────────────────────────────────────────────────────────────────
# 13. CONFIGURE BARPLOT COLORS AND THEME
# ─────────────────────────────────────────────────────────────────────────────
# Color scheme:
#   - Post-ASCT: Deep teal (#31688E) - early therapeutic assessment
#   - Maintenance-1yr: Bright green (#35B779) - longer-term surveillance
# Theme minimizes clutter while maintaining publication quality

custom_cols <- c(
  "Post-ASCT"       = "#31688E",  # deep teal
  "Maintenance-1yr" = "#35B779"   # bright green
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title      = element_text(size = 11),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top",
    plot.margin     = margin(10, 10, 30, 10)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 14. BUILD GROUPED BARPLOT: Sensitivity by Assay × Timepoint
# ─────────────────────────────────────────────────────────────────────────────
# Grouped bars allow visual comparison of:
#   - Sensitivity across assays (x-axis: clonoSEQ, MFC, EasyM, cfWGS)
#   - Timepoint effect (bars grouped by color: Post-ASCT vs Maintenance-1yr)
# Bar height = sensitivity %; text labels show exact percentages

p_sens <- ggplot(sens_df_bm,
                 aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  # position_dodge separates bars for same assay; width controls bar width
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",        # black outline for clarity
           size     = 0.3) +
  # Add percentage labels on top of bars
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,           # position above bar
            size     = 3.5) +
  # Use custom color mapping with explicit order
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols[c("Post-ASCT", "Maintenance-1yr")],  # enforce mapping
    limits = c("Post-ASCT", "Maintenance-1yr")                # enforce order
  ) +
  # Y-axis: 0-100% scale
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of cfDNA and Clinical MRD Assays in Relapsing Patients",
    x     = "Technology",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

print(p_sens)

# SAVE BM-SUBSET SENSITIVITY BARPLOT (500 DPI for manuscripts)
# Saves to: Final Tables and Figures/Supp_6A_Fig_sensitivity_by_tech_training3_{date}.png

ggsave(
  filename = paste0("Final Tables and Figures/Supp_6A_Fig_sensitivity_by_tech_training3_", date_tag, ".png"),
  plot     = p_sens,
  width    = 6,
  height   = 4,
  dpi      = 500
)

cat("  ✓ Saved BM-subset sensitivity barplot\n\n")

# BLOOD-SUBSET SENSITIVITY BARPLOT - Same structure as BM analysis
# Restricted to blood-derived cfWGS samples for head-to-head comparison


## Now for blood
sens_df_blood <- bind_rows(
  post_stats_blood  %>% mutate(Timepoint = "Post-ASCT"),
  year_stats_blood  %>% mutate(Timepoint = "Maintenance-1yr")
) %>%
  filter(Assay != "cfWGS_BM") %>%
  filter(Assay != "cfWGS_BM_screen") %>%
  mutate(
    # Convert to percentage for labeling
    Sens_pct   = Sensitivity * 100,
    # Manuscript-friendly assay names (with line breaks for blood cfWGS variants)
    Assay      = recode(Assay,
                        EasyM                 = "EasyM (Opt)",
                        clonoSEQ              = "clonoSEQ",
                        Flow                  = "MFC",
                        cfWGS_Blood_Sites     = "cfWGS\n(Sites Model)",
                        cfWGS_Blood_Combined  = "cfWGS\n(Combined Model)"),
    # Enforce timepoint order: Post-ASCT → Maintenance-1yr
    Timepoint = factor(Timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
  ) 

# Enforce assay ordering: cfWGS variants (blood) first, then clonoSEQ, MFC, EasyM (Opt) rightmost
# Note: Blood models use "\n" for line break in x-axis labels; EasyM (Opt) is rightmost for consistency
sens_df_blood <- sens_df_blood %>%
  mutate(Assay = factor(Assay,
                        levels = c("cfWGS\n(Sites Model)",
                                   "cfWGS\n(Combined Model)", 
                                   "clonoSEQ", "MFC", "EasyM (Opt)")))

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: Supp_8A (Blood sensitivity barplot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  sens_df_blood,
  file.path(outdir_source_data, paste0("Supp_8A_blood_sensitivity_barplot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: Supp_8A (Blood sensitivity)\n")

# Build blood-subset barplot (same structure as BM version)
p_sens_blood <- ggplot(sens_df_blood,
                       aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols[c("Post-ASCT", "Maintenance-1yr")],  # enforce mapping
    limits = c("Post-ASCT", "Maintenance-1yr")                # enforce order
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of cfDNA and Clinical MRD Assays in Relapsing Patients",
    x     = "Technology",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

# Display the blood-subset barplot
p_sens_blood

# SAVE BLOOD-SUBSET SENSITIVITY BARPLOT (500 DPI for manuscripts)
# Saves to: Final Tables and Figures/Supp_8A_Fig_sensitivity_by_tech_training_blood2_{date}.png

ggsave(
  filename = paste0("Final Tables and Figures/Supp_8A_Fig_sensitivity_by_tech_training_blood2_", date_tag, ".png"),
  plot     = p_sens_blood,
  width    = 6,
  height   = 4,
  dpi      = 500
)

cat("  ✓ Saved blood-subset sensitivity barplot\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 17. NON-FRONTLINE COHORT: Time-Window Prediction Analysis
# (Non-annotated cohort with time-window predictions from 3_1_A)
# ─────────────────────────────────────────────────────────────────────────────

## Time-window analysis: Prediction of relapse within specific time windows
## Using models trained on frontline cohort, evaluate on non-frontline subjects
## Subset to post-transplant and BM-cfWGS tested

# Filter to 1-year maintenance timepoint with BM-cfWGS result
df_km <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(BM_zscore_only_detection_rate_call)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 18. CALCULATE 24-MONTH RELAPSE-FREE SURVIVAL (RFS) BY CFWGS BM
# ─────────────────────────────────────────────────────────────────────────────
# RFS at 24 months = probability of not relapsing within 24 months
# Stratified by BM-cfWGS MRD status (positive vs negative)

fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
  data = df_km
)

# Extract RFS probability at 24 months (convert months to days: 24 * 30.44)
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# Extract RFS % for MRD- (stratum 1) and MRD+ (stratum 2)
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

# Calculate hazard ratio and 95% CI for cfWGS BM positivity
cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

# ─────────────────────────────────────────────────────────────────────────────
# 19. CALCULATE MEDIAN RELAPSE-FREE SURVIVAL BY CFWGS BM
# Median RFS = time at which 50% of patients have relapsed
# ─────────────────────────────────────────────────────────────────────────────

med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44   # Convert days to months
med_pos_cf <- med_cf[2] / 30.44   # Convert days to months

# ─────────────────────────────────────────────────────────────────────────────
# 20. CALCULATE 24-MONTH RFS BY FLOW CYTOMETRY
# Comparable time-window analysis using MFC/flow cytometry data
# ─────────────────────────────────────────────────────────────────────────────

# Filter to subjects with both 1-year maintenance and Flow_Binary data
df_km_fl <- df_km %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Flow_Binary)
  )

# Fit KM curve stratified by Flow MRD status
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)

# Extract RFS at 24 months for flow cytometry
sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

# ─────────────────────────────────────────────────────────────────────────────
# 21. CALCULATE MEDIAN RFS BY FLOW CYTOMETRY
# ─────────────────────────────────────────────────────────────────────────────

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km_fl),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high


## check clonoSEQ
# fit the Kaplan–Meier curve
df_km_clonoSEQ <- df_km %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km_clonoSEQ),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high


## check EasyM (Proteomic MRD)
# Filter to subjects with EasyM data at 1-year maintenance
df_km_easym <- df_km %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(EasyM_optimized_binary)
  )

# fit the Kaplan–Meier curve
fit_easym <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
  data = df_km_easym
)

sum_easym24   <- summary(fit_easym, times = t24)
rfs_neg_easym <- sum_easym24$surv[1] * 100
rfs_pos_easym <- sum_easym24$surv[2] * 100

## Median RFS by EasyM ----
fit_easym <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
  data = df_km_easym
)
med_easym <- surv_median(fit_easym)$median      # vector of two values
med_neg_easym <- med_easym[1] / 30.44           # convert days→months
med_pos_easym <- med_easym[2] / 30.44

cox_easym <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
        data = df_km_easym),
  exponentiate = TRUE, conf.int = TRUE
)
hr_easym      <- cox_easym$estimate
ci_lo_easym   <- cox_easym$conf.low
ci_hi_easym   <- cox_easym$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "BM_zscore_only_detection_rate_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$BM_zscore_only_detection_rate_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## --- Ns (head-to-head subset: requires cfWGS present) -------------------
n_cfWGS <- df_km %>%
  filter(!is.na(BM_zscore_only_detection_rate_call)) %>%
  distinct(Patient) %>%
  nrow()

n_MFC <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>%
  nrow()

n_clonoSEQ <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>%
  nrow()

n_easym <- df_km %>%
  filter(!is.na(EasyM_optimized_binary)) %>%
  distinct(Patient) %>%
  nrow()

## 6. Draft paragraph ----
paragraph <- glue(
  "After one year of maintenance therapy, among patients with cfWGS available (n={n_cfWGS}), ",
  "BM-cfWGS MRD-negative patients had {round(rfs_neg_cf)}% relapse-free survival at 24 months ",
  "versus {round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {ifelse(is.na(med_neg_cf),'NR',round(med_neg_cf,1))} vs ",
  "{ifelse(is.na(med_pos_cf),'NR',round(med_pos_cf,1))} months. ",
  "For MFC (n={n_MFC}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months ",
  "versus {round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {ifelse(is.na(med_neg_fl),'NR',round(med_neg_fl,1))} vs ",
  "{ifelse(is.na(med_pos_fl),'NR',round(med_pos_fl,1))} months. ",
  "clonoSEQ (n={n_clonoSEQ}) showed a similar direction of effect ",
  "(HR = {round(hr_clonoSEQ,2)}; 95% CI [{round(ci_lo_clonoSEQ,2)}–{round(ci_hi_clonoSEQ,2)}]). ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), comparable to flow ",
  "(ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients in the head-to-head subset, the minimum ",
  "detectable HR for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)


writeLines(paragraph)


## Get tibble 
# 7) Assemble a single summary table of all key metrics ---------------------

# gather into a single-row tibble
metrics_1yr <- tibble(
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  RFS24_em_neg   = rfs_neg_easym,
  RFS24_em_pos   = rfs_pos_easym,
  MedRFS_em_neg  = med_neg_easym,
  MedRFS_em_pos  = med_pos_easym,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
  HR_em          = hr_easym,
  CI_low_em      = ci_lo_easym,
  CI_high_em     = ci_hi_easym,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS,
  N_MFC          = n_MFC,
  N_clonoSEQ     = n_clonoSEQ,
  N_easym        = n_easym
)



### Redo for post-transplant 
## 1. Subset to post-transplant & BM-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(BM_zscore_only_detection_rate_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# we already have `fit_fl` from before
df_km_fl <- df_km %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Flow_Binary)
  )

# fit the Kaplan–Meier curve
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)

sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km_fl),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

## Check clonoSEQ 
# fit the Kaplan–Meier curve
df_km_clonoSEQ <- df_km %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km_clonoSEQ),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high


## check EasyM (Proteomic MRD)
# Filter to subjects with EasyM data at post-transplant
df_km_easym <- df_km %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(EasyM_optimized_binary)
  )

# fit the Kaplan–Meier curve
fit_easym <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
  data = df_km_easym
)

sum_easym24   <- summary(fit_easym, times = t24)
rfs_neg_easym <- sum_easym24$surv[1] * 100
rfs_pos_easym <- sum_easym24$surv[2] * 100

## Median RFS by EasyM ----
fit_easym <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
  data = df_km_easym
)
med_easym <- surv_median(fit_easym)$median      # vector of two values
med_neg_easym <- med_easym[1] / 30.44           # convert days→months
med_pos_easym <- med_easym[2] / 30.44

cox_easym <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
        data = df_km_easym),
  exponentiate = TRUE, conf.int = TRUE
)
hr_easym      <- cox_easym$estimate
ci_lo_easym   <- cox_easym$conf.low
ci_hi_easym   <- cox_easym$conf.high



## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "BM_zscore_only_detection_rate_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$BM_zscore_only_detection_rate_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## Get N
n_cfWGS <- df_km %>%
  filter(!is.na(BM_zscore_only_detection_rate_call)) %>%
  distinct(Patient) %>%
  nrow()

n_MFC <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>%
  nrow()

n_clonoSEQ <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>%
  nrow()

n_easym <- df_km %>%
  filter(!is.na(EasyM_optimized_binary)) %>%
  distinct(Patient) %>%
  nrow()

## 6. Draft paragraph ----
paragraph <- glue(
  "At post-transplant, BM-cfWGS (n={n_cfWGS}) MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC (n={n_MFC}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "clonoSEQ (n={n_clonoSEQ}) showed a similar direction of effect, ",
  "though with greater uncertainty due to fewer patients. ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), ",
  "comparable to flow cytometry (ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients, the minimum detectible HR ",
  "for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph)

## Compile 
metrics_post_transplant <- tibble(
  Landmark        = "post_transplant",
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  RFS24_em_neg   = rfs_neg_easym,
  RFS24_em_pos   = rfs_pos_easym,
  MedRFS_em_neg  = med_neg_easym,
  MedRFS_em_pos  = med_pos_easym,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
  HR_em          = hr_easym,
  CI_low_em      = ci_lo_easym,
  CI_high_em     = ci_hi_easym,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS,
  N_MFC          = n_MFC,
  N_clonoSEQ     = n_clonoSEQ,
  N_easym        = n_easym
)

metrics_1yr <- metrics_1yr %>%
  mutate(Landmark = "1yr_maintenance") %>%
  select(Landmark, everything())  # move Landmark to front

progression_metrics <- bind_rows(
  metrics_post_transplant,
  metrics_1yr
)

# 8) Export summary table --------------------------------------------
write_csv(
  progression_metrics,
  file.path(outdir, paste0("cfWGS_vs_flow_progression_summary_", date_tag, ".csv"))
)

# (Optional) also save as RDS for later use
saveRDS(
  progression_metrics,
  file.path(outdir, paste0("cfWGS_vs_flow_progression_summary_updated_", date_tag, ".rds"))
)



##### Now get the stats on cfWGS from blood derived muts as well 
#### Now get other results and do power analysis 

## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Blood_zscore_only_sites_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# fit the Kaplan–Meier curve
df_km_fl <- df_km %>%
  filter(
    !is.na(Flow_Binary)
  )

# fit the Kaplan–Meier curve
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)

sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

## 3b. Assess EasyM (Proteomic MRD) ----
# Filter patients with EasyM data available
df_km_easym <- df_km %>%
  filter(!is.na(EasyM_optimized_binary))

# Kaplan-Meier curve for EasyM
fit_easym <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
  data = df_km_easym
)

# 24-month RFS by EasyM
sum_easym24 <- summary(fit_easym, times = t24)
rfs_neg_easym <- sum_easym24$surv[1] * 100
rfs_pos_easym <- sum_easym24$surv[2] * 100

# Median RFS by EasyM (days → months)
med_easym <- surv_median(fit_easym)$median
med_neg_easym <- med_easym[1] / 30.44
med_pos_easym <- med_easym[2] / 30.44

# Cox regression for EasyM hazard ratio
cox_easym <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_easym     <- cox_easym$estimate
ci_lo_easym  <- cox_easym$conf.low
ci_hi_easym  <- cox_easym$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "Blood_zscore_only_sites_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_zscore_only_sites_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_zscore_only_sites_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

n_easym_blood <- df_km %>%
  filter(!is.na(EasyM_optimized_binary)) %>%
  distinct(Patient) %>% nrow()

## Blood paragraph with Ns (and NR-safe medians)
paragraph_blood <- glue(
  "After one year of maintenance therapy, among patients with Blood-cfWGS available (n={n_cfWGS_blood}), ",
  "Blood-cfWGS MRD-negative patients had {round(rfs_neg_cf)}% relapse-free survival at 24 months ",
  "versus {round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by Blood-cfWGS was {ifelse(is.na(med_neg_cf),'NR',round(med_neg_cf,1))} vs ",
  "{ifelse(is.na(med_pos_cf),'NR',round(med_pos_cf,1))} months. ",
  "For MFC (n={n_MFC_blood}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months ",
  "versus {round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {ifelse(is.na(med_neg_fl),'NR',round(med_neg_fl,1))} vs ",
  "{ifelse(is.na(med_pos_fl),'NR',round(med_pos_fl,1))} months. ",
  "clonoSEQ (n={n_clonoSEQ_blood}) showed a similar direction of effect ",
  "(HR = {round(hr_clonoSEQ,2)}; 95% CI [{round(ci_lo_clonoSEQ,2)}–{round(ci_hi_clonoSEQ,2)}]). ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), comparable to flow ",
  "(ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients in the head-to-head subset, the minimum ",
  "detectable HR for 80% power is {round(hr80,1)} (and power to detect HR = {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph_blood)


## Get tibble 
# 7) Assemble a single summary table of all key metrics ---------------------

# gather into a single-row tibble
metrics_1yr <- tibble(
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  RFS24_em_neg   = rfs_neg_easym,
  RFS24_em_pos   = rfs_pos_easym,
  MedRFS_em_neg  = med_neg_easym,
  MedRFS_em_pos  = med_pos_easym,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  HR_em          = hr_easym,
  CI_low_em      = ci_lo_easym,
  CI_high_em     = ci_hi_easym,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood,
  N_easym        = n_easym_blood
)



### Redo for post-transplant 
## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Blood_zscore_only_sites_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# fit the Kaplan–Meier curve
df_km_fl <- df_km %>%
  filter(
    !is.na(Flow_Binary)
  )

# fit the Kaplan–Meier curve
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)

sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

## 3b. Assess EasyM (Proteomic MRD) ----
# Filter patients with EasyM data available
df_km_easym <- df_km %>%
  filter(!is.na(EasyM_optimized_binary))

# Kaplan-Meier curve for EasyM
fit_easym <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
  data = df_km_easym
)

# 24-month RFS by EasyM
sum_easym24 <- summary(fit_easym, times = t24)
rfs_neg_easym <- sum_easym24$surv[1] * 100
rfs_pos_easym <- sum_easym24$surv[2] * 100

# Median RFS by EasyM (days → months)
med_easym <- surv_median(fit_easym)$median
med_neg_easym <- med_easym[1] / 30.44
med_pos_easym <- med_easym[2] / 30.44

# Cox regression for EasyM hazard ratio
cox_easym <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ EasyM_optimized_binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_easym     <- cox_easym$estimate
ci_lo_easym  <- cox_easym$conf.low
ci_hi_easym  <- cox_easym$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "Blood_zscore_only_sites_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_zscore_only_sites_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_zscore_only_sites_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

n_easym_blood <- df_km %>%
  filter(!is.na(EasyM_optimized_binary)) %>%
  distinct(Patient) %>% nrow()

## 6. Draft paragraph ----
paragraph <- glue(
  "At post-transplant, Blood-cfWGS MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC, MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), ",
  "comparable to flow cytometry (ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients, the minimum detectible HR ",
  "for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph)

## Compile 
metrics_post_transplant <- tibble(
  Landmark        = "post_transplant",
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  RFS24_em_neg   = rfs_neg_easym,
  RFS24_em_pos   = rfs_pos_easym,
  MedRFS_em_neg  = med_neg_easym,
  MedRFS_em_pos  = med_pos_easym,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  HR_em          = hr_easym,
  CI_low_em      = ci_lo_easym,
  CI_high_em     = ci_hi_easym,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood,
  N_easym        = n_easym_blood
)

metrics_1yr <- metrics_1yr %>%
  mutate(Landmark = "1yr_maintenance") %>%
  select(Landmark, everything())  # move Landmark to front

progression_metrics_blood <- bind_rows(
  metrics_post_transplant,
  metrics_1yr
)

# 8) Export summary table --------------------------------------------
write_csv(
  progression_metrics_blood,
  file.path(outdir, paste0("cfWGS_vs_flow_progression_summary_blood_muts_", date_tag, ".csv"))
)

# (Optional) also save as RDS for later use
saveRDS(
  progression_metrics_blood,
  file.path(outdir, paste0("cfWGS_vs_flow_progression_summaryy_blood_muts_updated_", date_tag, ".rds"))
)




#### Now redo with the combined model for blood 

## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Blood_plus_fragment_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# fit the Kaplan–Meier curve
df_km_fl <- df_km %>%
  filter(
    !is.na(Flow_Binary)
  )

# fit the Kaplan–Meier curve
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)

sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "Blood_plus_fragment_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_plus_fragment_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_plus_fragment_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

## Blood paragraph with Ns (and NR-safe medians)
paragraph_blood <- glue(
  "After one year of maintenance therapy, among patients with Blood-cfWGS available (n={n_cfWGS_blood}), ",
  "Blood-cfWGS MRD-negative patients had {round(rfs_neg_cf)}% relapse-free survival at 24 months ",
  "versus {round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by Blood-cfWGS was {ifelse(is.na(med_neg_cf),'NR',round(med_neg_cf,1))} vs ",
  "{ifelse(is.na(med_pos_cf),'NR',round(med_pos_cf,1))} months. ",
  "For MFC (n={n_MFC_blood}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months ",
  "versus {round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {ifelse(is.na(med_neg_fl),'NR',round(med_neg_fl,1))} vs ",
  "{ifelse(is.na(med_pos_fl),'NR',round(med_pos_fl,1))} months. ",
  "clonoSEQ (n={n_clonoSEQ_blood}) showed a similar direction of effect ",
  "(HR = {round(hr_clonoSEQ,2)}; 95% CI [{round(ci_lo_clonoSEQ,2)}–{round(ci_hi_clonoSEQ,2)}]). ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), comparable to flow ",
  "(ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients in the head-to-head subset, the minimum ",
  "detectable HR for 80% power is {round(hr80,1)} (and power to detect HR = {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph_blood)


## Get tibble 
# 7) Assemble a single summary table of all key metrics ---------------------

# gather into a single-row tibble
metrics_1yr <- tibble(
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood
)



### Redo for post-transplant 
## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Blood_plus_fragment_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# fit the Kaplan–Meier curve
df_km_fl <- df_km %>%
  filter(
    !is.na(Flow_Binary)
  )

# fit the Kaplan–Meier curve
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)

sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "Blood_plus_fragment_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_plus_fragment_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_plus_fragment_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

## 6. Draft paragraph ----
paragraph <- glue(
  "At post-transplant, Blood-cfWGS MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC, MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), ",
  "comparable to flow cytometry (ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients, the minimum detectible HR ",
  "for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph)

## Compile 
metrics_post_transplant <- tibble(
  Landmark        = "post_transplant",
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood
)

metrics_1yr <- metrics_1yr %>%
  mutate(Landmark = "1yr_maintenance") %>%
  select(Landmark, everything())  # move Landmark to front

progression_metrics_blood_combined <- bind_rows(
  metrics_post_transplant,
  metrics_1yr
)

# 8) Export summary table --------------------------------------------
write_csv(
  progression_metrics_blood_combined,
  file.path(outdir, paste0("cfWGS_vs_flow_progression_summary_blood_muts_combined_model_", date_tag, ".csv"))
)

# (Optional) also save as RDS for later use
saveRDS(
  progression_metrics_blood_combined,
  file.path(outdir, paste0("cfWGS_vs_flow_progression_summary_blood_muts_updated_combined_model_", date_tag, ".rds"))
)










### Make HR figure 
# 1. reshape into long format
hr_plot_df_blood <- progression_metrics_blood %>%
  select(Landmark,
         HR_cf,   CI_low_cf,   CI_high_cf,
         HR_fl,   CI_low_fl,   CI_high_fl,
         HR_seq,   CI_low_seq,   CI_high_seq,
         HR_em,   CI_low_em,   CI_high_em) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf|fl|seq|em)"
  ) %>%
  mutate(
    Assay = recode(Assay,
                   cf = "cfWGS (Sites Model)",
                   fl = "MFC",
                   seq = "clonoSEQ",
                   em = "EasyM (Opt)"),
    Assay = factor(Assay,
                   levels = c("cfWGS (Sites Model)", "clonoSEQ", "MFC", "EasyM (Opt)")),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post‑ASCT", "Maintenance-1yr"))
  )

# reshape combined model (cfWGS only)
hr_plot_df_combined <- progression_metrics_blood_combined %>%
  select(Landmark,
         HR_cf, CI_low_cf, CI_high_cf) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf)"
  ) %>%
  mutate(
    Assay = "cfWGS (Combined Model)",
    Assay = factor(Assay,
                   levels = c("cfWGS (Sites Model)", "cfWGS (Combined Model)", "clonoSEQ", "MFC", "EasyM (Opt)")),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post‑ASCT", "Maintenance-1yr"))
  )


# bind together Sites + Combined models
hr_plot_df_blood <- bind_rows(hr_plot_df_blood, hr_plot_df_combined) %>%
  mutate(Assay = factor(Assay,
                        levels = c("cfWGS (Sites Model)", "cfWGS (Combined Model)", "clonoSEQ", "MFC", "EasyM (Opt)")))

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: SuppFig8B (Blood HR plot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  hr_plot_df_blood,
  file.path(outdir_source_data, paste0("SuppFig8B_blood_HR_plot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: SuppFig8B (Blood HR)\n")

p_hr <- ggplot(hr_plot_df_blood,
               aes(x = HR, y = fct_rev(Landmark), colour = Assay)) +
  # reference line
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # 1) horizontal CIs
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    position = position_dodge(width = 0.6),
    size     = 0.5
  ) +
  
  # 2) dots at the HR
  geom_point(
    position = position_dodge(width = 0.6),
    size     = 3
  ) +
  
  # log scale axis
  scale_x_continuous(
    "Hazard ratio (log scale)",
    trans        = "log10",
    limits       = c(0.05, 50),
    breaks       = c(0.1, 0.5, 1, 2, 5, 10, 20, 50),
    minor_breaks = c(
      0.05, 0.06, 0.08,           # between 0.05 & 0.1
      0.15, 0.2, 0.3, 0.4,        # between 0.1 & 0.5
      0.6, 0.8,                   # between 0.5 & 1
      1.5, 3,                      # between 1 & 5
      6, 8,                        # between 5 & 10
      15, 30                       # between 10 & 50
    ),
    labels = label_number(accuracy = .1)
  )+
  annotation_logticks(
    sides  = "b",
    short  = unit(2, "pt"),
    mid    = unit(4, "pt"),
    long   = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    name   = NULL,
    values = c("cfWGS (Sites Model)" = "#35608DFF",
               "cfWGS (Combined Model)" = "#440154FF",
               "MFC"   = "#43BF71FF",
               "clonoSEQ"= "#E69F00FF",   # orange for clonoSEQ
               "EasyM (Opt)" = "#D81B60FF")  # magenta for EasyM
  ) +
  
  labs(
    y        = NULL,
    title    = "Relapse hazard ratios stratified by MRD assay\nand landmark timepoint",
    #   subtitle = "cfWGS vs. MFC (95% CI)"
  ) +
  
  # classic theme with no gridlines
  theme_classic(base_size = 11) +
  theme(
    panel.grid         = element_blank(),   # no grid at all
    plot.title         = element_text(face = "bold",
                                      hjust = 0.5),  # bold + centered
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8)
  )

ggsave(paste0("Final Tables and Figures/SuppFig8B_cfWGS_blood_HR_updated3_", date_tag, ".png"),
       p_hr, width = 6, height = 4, dpi = 600)



### Now for BM derived muts
hr_plot_df_bm <- progression_metrics %>%
  select(Landmark,
         HR_cf,   CI_low_cf,   CI_high_cf,
         HR_fl,   CI_low_fl,   CI_high_fl,
         HR_seq,   CI_low_seq,   CI_high_seq,
         HR_em,   CI_low_em,   CI_high_em) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf|fl|seq|em)"
  ) %>%
  mutate(
    Assay = recode(Assay,
                   cf = "cfWGS",
                   fl = "MFC",
                   seq = "clonoSEQ",
                   em = "EasyM (Opt)"),
    Assay = factor(Assay,
                   levels = c("cfWGS", "clonoSEQ", "MFC", "EasyM (Opt)")),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post‑ASCT", "Maintenance-1yr"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: Supp_Figure_6B (BM HR plot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  hr_plot_df_bm,
  file.path(outdir_source_data, paste0("Supp_Figure_6B_BM_HR_plot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: Supp_Figure_6B (BM HR)\n")

p_hr_bm <- ggplot(hr_plot_df_bm,
                  aes(x = HR, y = fct_rev(Landmark), colour = Assay)) +
  # reference line
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # 1) horizontal CIs
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    position = position_dodge(width = 0.6),
    size     = 0.5
  ) +
  
  # 2) dots at the HR
  geom_point(
    position = position_dodge(width = 0.6),
    size     = 3
  ) +
  
  # log scale axis
  scale_x_log10(
    "Hazard ratio (log scale)",
    limits       = c(0.19, 250),           # extended to 250 to accommodate ~220 CI
    breaks       = c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 250),
    minor_breaks = c(
      0.25, 0.3, 0.4, 0.6, 0.8,    # between 0.2 & 1
      1.5,           # between 1 & 2
      3, 4,          # between 2 & 5
      6, 8,          # between 5 & 10
      15, 30,        # between 10 & 50
      40, 60, 80,    # between 20 & 100
      150            # between 100 & 200
    ),
    labels = function(x) {
      sapply(x, function(xx) {
        if (xx > 1) {
          sprintf("%.0f", xx)
        } else {
          sprintf("%.2f", xx)
        }
      })
    }
  ) +
  annotation_logticks(
    sides = "b",
    short = unit(2, "pt"),
    mid   = unit(4, "pt"),
    long  = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    name   = NULL,
    values = c("cfWGS" = "#35608DFF",
               "MFC"   = "#43BF71FF",
               "clonoSEQ"= "#E69F00FF",   # orange for clonoSEQ
               "EasyM (Opt)" = "#D81B60FF")    # magenta/pink for EasyM
  ) +
  
  labs(
    y        = NULL,
    title    = "Relapse hazard ratios stratified by MRD assay\nand landmark timepoint",
    #  subtitle = "(95% CI)"
  ) +
  
  # classic theme with no gridlines
  theme_classic(base_size = 11) +
  theme(
    panel.grid         = element_blank(),   # no grid at all
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold",
                                      hjust = 0.5),  # bold + centered
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8)
  )

p_hr_bm

ggsave(paste0("Final Tables and Figures/Supp_Figure_6B_cfWGS_BM_HR_updated3_", date_tag, ".png"),
       p_hr_bm, width = 6, height = 4, dpi = 600)






### Now make time to relapse figure 
df <- survival_df %>%                           # <- your tibble
  # keep samples beyond baseline / diagnosis
  filter(!str_detect(timepoint_info, regex("Diagnosis|Baseline", TRUE))) %>%
  
  # drop rows with missing probability or time
  filter(!is.na(BM_zscore_only_detection_rate_prob),
         !is.na(Time_to_event)) %>%
  
  # enforce non-negative time to event for relapse event visits a few days off from CMRG date
  mutate(
    days_before_event = pmax(Time_to_event, 0), # set negative values to 0 
    mrd_status      = factor(
      BM_zscore_only_detection_rate_call,
      levels = c(0, 1),
      labels = c("MRD-", "MRD+")
    ),
    progress_status = factor(
      Relapsed_Binary,
      levels = c(0, 1),
      labels = c("No relapse", "Relapse")
    ),
    
    # time *before* the anchor (positive value) → plot reversed
    days_before_event = Time_to_event,    # keep positive for clarity
    months_before_event = days_before_event/30.44
  )

## Only multiple points 
df_slim <- df %>%
  group_by(Patient) %>% 
  filter(dplyr::n() > 1) %>%   # keep only patients with >1 row
  ungroup()

# ────────────────────────────────────────────────────────────────
# 2.  Plot  ──────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
youden_thresh <- 0.4215524
#youden_thresh2 <- 0.35

max_mo <- max(df$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df, aes(months_before_event, BM_zscore_only_detection_rate_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
   # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status),
  shape  = 21,
  colour = "black",
  size   = 2
  ) +
  
  # 4) event line
  #geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("cVAF Model Probability",
                     limits = c(0,1),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  # scale_colour_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  # scale_fill_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  
  # # 7) stroke scale for MRD call
  # scale_discrete_manual(
  #   aesthetics = "stroke",
  #   values     = c("MRD-" = 0, "MRD+" = 1),
  #   guide      = guide_legend(
  #     title = "MRD call",
  #     override.aes = list(
  #       shape  = 21,
  #       fill   = "white",   # white interior in legend, makes stroke obvious
  #       size   = 4,
  #       colour = "black",
  #       stroke = c(0, 1)
  #     )
  #   )
  # ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "Longitudinal cfWGS MRD Probability by Patient Outcome\nUsing BM-Derived Mutation Lists"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "bottom", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 10)    # even smaller
  )

print(p_prob)

## Add the samples 
# Make sure the outcome labels match your scales
df <- df %>% mutate(progress_status = factor(progress_status,
                                             levels = c("No relapse","Relapse")))

n_patients_by <- df %>%
  distinct(Patient, progress_status) %>%
  count(progress_status, name = "n_patients")

n_timepoints_by <- df %>%
  count(progress_status, name = "n_timepoints")

# pull counts (0 if a group is absent)
get_n <- function(tbl, lvl, col) {
  val <- tbl %>% filter(progress_status == lvl) %>% pull({{col}})
  if (length(val) == 0) 0 else val
}

n_pat_nr  <- get_n(n_patients_by,  "No relapse", n_patients)
n_time_nr <- get_n(n_timepoints_by, "No relapse", n_timepoints)
n_pat_rl  <- get_n(n_patients_by,  "Relapse",    n_patients)
n_time_rl <- get_n(n_timepoints_by, "Relapse",   n_timepoints)

# text to print
lab_nr <- paste0("No relapse: n=", n_pat_nr, " patients; ", n_time_nr, " samples")
lab_rl <- paste0("Relapse: n=", n_pat_rl, " patients; ", n_time_rl, " samples")

# place labels at bottom-left (remember: x is reversed, so 'left' == large x)
x_left <- max_mo - 0.02 * max_mo  # a small inset from the left border

p_prob2 <- p_prob +
  scale_colour_manual(
    name = "Patient outcome",
#    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
    values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  ) +
  scale_fill_manual(
    name = "Patient outcome",
 #   values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
   values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  )

print(p_prob2)

# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated5.png",
       p_prob, width = 6, height = 4.5, dpi = 600)

ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated5_label.png",
       p_prob2, width = 6, height = 4.5, dpi = 600)



### Instead change the scale to match what Trevor was thinking 
df_plot <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%   # months → days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                       # keep all progressors
      row_number() == which.min(days_before_event)       # keep *latest* censor
  ) %>%
  ungroup()

max_days <- ceiling(max(df_plot$days_before_event, na.rm = TRUE) / 180) * 180

df_plot <- df_plot %>%
  mutate(Time_to_event = if_else(Time_to_event >= -30 & Time_to_event < 0, 0, Time_to_event))

## ─────────────────────────────────────────────────────────────
## 1)  Build the scatter plot                                  
## ─────────────────────────────────────────────────────────────
p_time <- ggplot(df_plot,
                 aes(x = BM_zscore_only_detection_rate_prob,
                     y = days_before_event)) +
  
  # ① dashed horizontal line at “event” (0 days)
  #  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  
  # ② points – colour = outcome, stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape  = 21,
             size   = 3,
             colour = "black") +
  
  # ③ axes
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_y_reverse(
    "Days until relapse (or censor)",
    limits = c(max_days, 0),
    breaks = seq(0, max_days, by = 180),      # every ~6 months
    minor_breaks = seq(0, max_days, by = 90)  # every 3 months
  ) +
  
  # ④ colours for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # ⑤ stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1.1),   # ring only if MRD+
    guide = guide_legend(
      title          = "MRD call",
      override.aes   = list(shape = 21,
                            size  = 4,
                            colour = "black",
                            fill   = "white",
                            stroke = c(0, 1.1))
    )
  ) +
  
  # ⑥ theme / labels
  labs(
    title = "Time to relapse vs. cfWGS MRD probability"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  )

print(p_time)

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob2.png"),
       p_time, width = 6, height = 4, dpi = 600)


### Show non-relapsers as infinity
# 1) compute days & a “days_plot” that sends non-relapsers to ∞
plot_df2 <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%        # months→days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                                 # keep all relapsers
      row_number() == which.min(days_before_event)                   # for non-relapsers, keep their last sample
  ) %>%
  ungroup()

# define axis maximum and “infinity” sentinel
max_days <- ceiling(max(plot_df2$days_before_event, na.rm = TRUE) / 180) * 180
overflow <- max_days + 180   # a little beyond the longest follow‑up

plot_df2 <- plot_df2 %>%
  mutate(
    days_plot = if_else(progress_status == "No relapse", overflow, days_before_event)
  )


## Check corrs 
# 1) subset to relapsers
rel_df <- plot_df2 %>%
  filter(progress_status == "Relapse")

# 2) run cor.test for each metric
spearman_prob <- cor.test(
  rel_df$BM_zscore_only_detection_rate_prob,
  rel_df$days_before_event,
  method = "spearman"
)

spearman_zscore <- cor.test(
  rel_df$zscore_BM,
  rel_df$days_before_event,
  method = "spearman"
)

spearman_detect_rate <- cor.test(
  rel_df$detect_rate_BM,
  rel_df$days_before_event,
  method = "spearman"
)

# 3) print them
cat("\n--- cfWGS model probability vs days ---\n")
print(spearman_prob)

cat("\n--- zscore_BM vs days ---\n")
print(spearman_zscore)

cat("\n--- detect_rate_BM vs days ---\n")
print(spearman_detect_rate)


# Edit days
plot_df2 <- plot_df2 %>%
  mutate(days_plot = if_else(days_plot >= -30 & days_plot <= 0, 0, days_plot))

# 2) compute Spearman rho on the relapsers
spearman_res <- with(
  filter(plot_df2, progress_status == "Relapse"),
  cor.test(BM_zscore_only_detection_rate_prob,
           days_before_event,
           method = "spearman")
)

## Check other metrics
# 2A) pull out estimate + p‑value
rho_all  <- spearman_res$estimate
pval_all <- spearman_res$p.value

# B) excluding relapse samples (days_plot > 0 only)
spearman_pre <- with(
  filter(plot_df2, progress_status == "Relapse", days_plot > 0),
  cor.test(BM_zscore_only_detection_rate_prob,
           days_before_event,
           method = "spearman")
)

rho_pre  <- spearman_pre$estimate
pval_pre <- spearman_pre$p.value

pval_all_str <- ifelse(pval_all < 0.001, "<0.001", sprintf("%.3f", pval_all))
pval_pre_str <- ifelse(pval_pre < 0.001, "<0.001", sprintf("%.3f", pval_pre))

annot_text <- sprintf("All relapse samples:\nρ=%.2f, p=%s\nPre-relapse only:\nρ=%.2f, p=%s",
                      rho_all, pval_all_str,
                      rho_pre, pval_pre_str)

# 3) make the scatter
p_time_inf <- ggplot(plot_df2,
                     aes(x = BM_zscore_only_detection_rate_prob,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = 0.35, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0,1),
    breaks = seq(0,1,by=0.2),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # annotate Spearman
  # annotate("text",
  #          x = 0.01,      # left margin
  #          y = 1650,
  #          label = sprintf("ρ = %.2f\np = %s", rho, pval_str),
  #          hjust = 0,
  #          size = 3.5) +
  annotate("text",
           x = 0.01,
           y = 1450,
           label = annot_text,
           hjust = 0,
           size = 3.5) +

  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using BM-Derived Mutations")

print(p_time_inf)

# 4) render / save
ggsave(file.path("Final Tables and Figures/Fig_4D_time_to_relapse_infinity_no_reverse2_BM_muts_updated3.png"),
       p_time_inf, width = 5.5, height = 4.5, dpi = 600)


## Try different legend layout 
library(cowplot)    # get_legend()
library(patchwork)  # easy assembly

# --- 1) build the main plot *without* a legend (legend handled below) ---
p_main <- ggplot(plot_df2,
                 aes(x = BM_zscore_only_detection_rate_prob, y = days_plot)) +
  geom_vline(xintercept = 0.35, linetype = "dotted", colour = "grey40") +
  geom_point(aes(colour = progress_status, fill = progress_status),
             shape = 21, size = 3, colour = "black") +
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0.1,1),
    breaks = seq(0.1,1,by=0.2),
    labels = scales::percent_format(accuracy=1)
  )  +# colours
scale_colour_manual(
  "Patient outcome",
  values = c("No relapse" = "black", "Relapse" = "red")
) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  guides(fill = "none") +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none"        # <- hide here; we'll place it below
  ) +
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using BM-Derived Mutations")

# --- 2) LEGEND-ONLY PLOT (do NOT inherit guides(fill = 'none')) ---
# -----build a dummy legend that always shows both levels -----
p_legend_only <- ggplot(legend_df, aes(x, y, fill = progress_status)) +
  # make the plotting layer invisible *in the panel* …
  geom_point(shape = 21, size = 0, colour = "black", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = pal_vals,
    breaks = names(pal_vals)
  ) +
  guides(fill = guide_legend(
    ncol = 1,                 # <- stack items vertically
    byrow = TRUE,
    title.position = "top",   # title on its own line
    label.hjust = 0,          # left-align labels
    override.aes = list(shape = 21, size = 3, alpha = 1, colour = "black")
  )) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "vertical",         # <- vertical legend
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(6, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  )

# --- 3) make two small “text boxes” for the right columns ---
txt_all <- sprintf("All relapse samples:\nρ=%.2f, p=%s", rho_all, pval_all_str)
txt_pre <- sprintf("Pre-relapse only:\nρ=%.2f, p=%s", rho_pre, pval_pre_str)

mini_box <- function(s) {
  ggplot() +
    annotate("label", x = 0, y = 1, label = s,
             hjust = 0, vjust = 1, size = 3.5,
             label.size = 0, fill = scales::alpha("white", 0.7)) +
    xlim(0,1) + ylim(0,1) +
    theme_void()
}

col2 <- mini_box(txt_all)
col3 <- mini_box(txt_pre)

# --- 4) assemble: plot on top; 3 columns underneath ---
bottom_row <- p_legend_only | col2 | col3
final_plot <- p_main / bottom_row + plot_layout(heights = c(1, 0.22))

# show and save
print(final_plot)
ggsave("Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_BM_muts.png",
       final_plot, width = 5.5, height = 5.5, dpi = 600)



plot_df2 %>%
  filter(
    is.na(BM_zscore_only_detection_rate_prob) | is.na(days_plot) |
      BM_zscore_only_detection_rate_prob < 0 | BM_zscore_only_detection_rate_prob > 1 |
      days_plot < 0 | days_plot > overflow
  )

time_to_relapse_BM <- plot_df2



### Redo now for blood derived muts 
### Now redo for blood-derived muts
df <- survival_df %>%                           # <- your tibble
  # keep samples beyond baseline / diagnosis
  filter(!str_detect(timepoint_info, regex("Diagnosis|Baseline", TRUE))) %>%
  
  # drop rows with missing probability or time
  filter(!is.na(Blood_zscore_only_sites_prob),
         !is.na(Time_to_event)) %>%
  
  # enforce non-negative time to event for relapse event visits a few days off from CMRG date
  mutate(
    days_before_event = pmax(Time_to_event, 0), # set negative values to 0 
    mrd_status      = factor(
      Blood_zscore_only_sites_call,
      levels = c(0, 1),
      labels = c("MRD-", "MRD+")
    ),
    progress_status = factor(
      Relapsed_Binary,
      levels = c(0, 1),
      labels = c("No relapse", "Relapse")
    ),
    
    # time *before* the anchor (positive value) → plot reversed
    days_before_event = Time_to_event,    # keep positive for clarity
    months_before_event = days_before_event/30.44
  )

## Only multiple points 
df_slim <- df %>%
  group_by(Patient) %>% 
  filter(dplyr::n() > 1) %>%   # keep only patients with >1 row
  ungroup()

# ────────────────────────────────────────────────────────────────
# 2.  Plot  ──────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
youden_thresh <- 0.5166693
youden_thr <- youden_thresh # for consistency
max_mo <- max(df_slim$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df, aes(months_before_event, Blood_zscore_only_sites_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
  
  # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status),
    shape  = 21,
    colour = "black",
    size   = 2
  ) +
  
  # 4) event line
  #geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("Sites Model Probability",
                     limits = c(0.3,1),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  # scale_colour_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  # scale_fill_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  # # 7) stroke scale for MRD call
  # scale_discrete_manual(
  #   aesthetics = "stroke",
  #   values     = c("MRD-" = 0, "MRD+" = 1),
  #   guide      = guide_legend(
  #     title = "MRD call",
  #     override.aes = list(
  #       shape  = 21,
  #       fill   = "white",   # white interior in legend, makes stroke obvious
  #       size   = 4,
  #       colour = "black",
  #       stroke = c(0, 1)
  #     )
  #   )
  # ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "Longitudinal cfWGS MRD Probability by Patient Outcome\nUsing cfDNA-Derived Mutation Lists"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "bottom", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 10)    # even smaller
  )

print(p_prob)

### Add label 
# Make sure the outcome labels match your scales
df <- df %>% mutate(progress_status = factor(progress_status,
                                             levels = c("No relapse","Relapse")))

n_patients_by <- df %>%
  distinct(Patient, progress_status) %>%
  count(progress_status, name = "n_patients")

n_timepoints_by <- df %>%
  count(progress_status, name = "n_timepoints")

# pull counts (0 if a group is absent)
get_n <- function(tbl, lvl, col) {
  val <- tbl %>% filter(progress_status == lvl) %>% pull({{col}})
  if (length(val) == 0) 0 else val
}

n_pat_nr  <- get_n(n_patients_by,  "No relapse", n_patients)
n_time_nr <- get_n(n_timepoints_by, "No relapse", n_timepoints)
n_pat_rl  <- get_n(n_patients_by,  "Relapse",    n_patients)
n_time_rl <- get_n(n_timepoints_by, "Relapse",   n_timepoints)

# text to print
lab_nr <- paste0("No relapse: n=", n_pat_nr, " patients; ", n_time_nr, " samples")
lab_rl <- paste0("Relapse: n=", n_pat_rl, " patients; ", n_time_rl, " samples")

# place labels at bottom-left (remember: x is reversed, so 'left' == large x)
x_left <- max_mo - 0.02 * max_mo  # a small inset from the left border

p_prob2 <- p_prob +
  scale_colour_manual(
    name = "Patient outcome",
#    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
    values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  ) +
  scale_fill_manual(
    name = "Patient outcome",
   # values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
   values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  )



print(p_prob2)
# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated3_blood4.png",
       p_prob, width = 6, height = 4.5, dpi = 600)

ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated3_blood2_labelled3.png",
       p_prob2, width = 6, height = 4.5, dpi = 600)


## Way Trevor was thinking 
### Instead change the scale to match what Trevor was thinking 
df_plot <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%   # months → days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                       # keep all progressors
      row_number() == which.min(days_before_event)       # keep *latest* censor
  ) %>%
  ungroup()

max_days <- ceiling(max(df_plot$days_before_event, na.rm = TRUE) / 180) * 180

df_plot <- df_plot %>%
  mutate(Time_to_event = if_else(Time_to_event >= -30 & Time_to_event < 0, 0, Time_to_event))

## ─────────────────────────────────────────────────────────────
## 1)  Build the scatter plot                                  
## ─────────────────────────────────────────────────────────────
p_time <- ggplot(df_plot,
                 aes(x = Blood_zscore_only_sites_prob,
                     y = days_before_event)) +
  
  # ① dashed horizontal line at “event” (0 days)
  #  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  
  # ② points – colour = outcome, stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape  = 21,
             size   = 3,
             colour = "black") +
  
  # ③ axes
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_y_reverse(
    "Days until relapse (or censor)",
    limits = c(max_days, 0),
    breaks = seq(0, max_days, by = 180),      # every ~6 months
    minor_breaks = seq(0, max_days, by = 90)  # every 3 months
  ) +
  
  # ④ colours for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # ⑤ stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1.1),   # ring only if MRD+
    guide = guide_legend(
      title          = "MRD call",
      override.aes   = list(shape = 21,
                            size  = 4,
                            colour = "black",
                            fill   = "white",
                            stroke = c(0, 1.1))
    )
  ) +
  
  # ⑥ theme / labels
  labs(
    title = "Time to relapse vs. cfWGS MRD probability"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  )

print(p_time)

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob2.png"),
       p_time, width = 6, height = 4, dpi = 600)


### Show non-relapsers as infinity
# 1) compute days & a “days_plot” that sends non-relapsers to ∞
plot_df2 <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%        # months→days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                                 # keep all relapsers
      row_number() == which.min(days_before_event)                   # for non-relapsers, keep their last sample
  ) %>%
  ungroup()

# define axis maximum and “infinity” sentinel
max_days <- ceiling(max(plot_df2$days_before_event, na.rm = TRUE) / 180) * 180
overflow <- max_days + 180   # a little beyond the longest follow‑up

plot_df2 <- plot_df2 %>%
  mutate(
    days_plot = if_else(progress_status == "No relapse", overflow, days_before_event)
  )


## Check corrs 
# 1) subset to relapsers
rel_df <- plot_df2 %>%
  filter(progress_status == "Relapse")

# 2) run cor.test for each metric
spearman_prob <- cor.test(
  rel_df$Blood_zscore_only_sites_prob,
  rel_df$days_before_event,
  method = "spearman"
)

# 3) print them
cat("\n--- cfWGS model probability vs days ---\n")
print(spearman_prob)


# Edit days
plot_df2 <- plot_df2 %>%
  mutate(days_plot = if_else(days_plot >= -35 & days_plot <= 0, 0, days_plot))

# 2) compute Spearman rho on the relapsers
spearman_res <- with(
  filter(plot_df2, progress_status == "Relapse"),
  cor.test(Blood_zscore_only_sites_prob,
           days_before_event,
           method = "spearman")
)

## Check other metrics
# 2A) pull out estimate + p‑value
# A) including relapse samples (your current spearman_res)
rho_all  <- spearman_res$estimate
pval_all <- spearman_res$p.value

# B) excluding relapse samples (days_plot > 0 only)
spearman_pre <- with(
  filter(plot_df2, progress_status == "Relapse", days_plot > 0),
  cor.test(Blood_zscore_only_sites_prob,
           days_before_event,
           method = "spearman")
)

rho_pre  <- spearman_pre$estimate
pval_pre <- spearman_pre$p.value

pval_all_str <- ifelse(pval_all < 0.001, "<0.001", sprintf("%.3f", pval_all))
pval_pre_str <- ifelse(pval_pre < 0.001, "<0.001", sprintf("%.3f", pval_pre))

annot_text <- sprintf("All relapse samples:\nρ=%.2f, p=%s\nPre-relapse only:\nρ=%.2f, p=%s",
                      rho_all, pval_all_str,
                      rho_pre, pval_pre_str)

## See range
plot_df2 %>%
  summarise(
    min   = min(Blood_zscore_only_sites_prob, na.rm = TRUE),
    max   = max(Blood_zscore_only_sites_prob, na.rm = TRUE),
    range = max(Blood_zscore_only_sites_prob, na.rm = TRUE) -
      min(Blood_zscore_only_sites_prob, na.rm = TRUE)
  )

# 3) make the scatter
p_time_inf <- ggplot(plot_df2,
                     aes(x = Blood_zscore_only_sites_prob,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = youden_thr, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability (%)",
    limits = c(0.3,1),
    breaks = seq(0,1,by=0.1),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # annotate Spearman
  # annotate("text",
  #          x = 0.90,      # left margin
  #          y = 1650,
  #          label = sprintf("ρ = %.2f\np = %s", rho, pval_fmt),
  #          hjust = 0,
  #          size = 3.5) +
  annotate("text",
           x = 0.8,
           y = 1450,
           label = annot_text,
           hjust = 0,
           size = 3.5) +
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using cfDNA-Derived Mutations")

print(p_time_inf)

# 4) render / save
ggsave(file.path("Final Tables and Figures/Fig_4_D_time_to_relapse_infinity_no_reverse2_blood_muts4_small2.png"),
       p_time_inf, width = 5.5, height = 4.5, dpi = 600)


time_to_relapse_blood <- plot_df2

## Check what is not included
plot_df2 %>%
  filter(
    is.na(Blood_zscore_only_sites_prob) | is.na(days_plot) |
      Blood_zscore_only_sites_prob < 0.3 | Blood_zscore_only_sites_prob > 1 |
      days_plot < 0 | days_plot > overflow
  )



### Do different legend layout 

## Try different legend layout 

# --- 1) build the main plot *without* a legend (legend handled below) ---
p_main <-  ggplot(plot_df2,
                  aes(x = Blood_zscore_only_sites_prob,
                      y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = youden_thr, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability (%)",
    limits = c(0.3,1),
    breaks = seq(0,1,by=0.1),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using cfDNA-Derived Mutations")

# --- 2) LEGEND-ONLY PLOT (do NOT inherit guides(fill = 'none')) ---
# -----build a dummy legend that always shows both levels -----
p_legend_only <- ggplot(legend_df, aes(x, y, fill = progress_status)) +
  # make the plotting layer invisible *in the panel* …
  geom_point(shape = 21, size = 0, colour = "black", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = pal_vals,
    breaks = names(pal_vals)
  ) +
  guides(fill = guide_legend(
    ncol = 1,                 # <- stack items vertically
    byrow = TRUE,
    title.position = "top",   # title on its own line
    label.hjust = 0,          # left-align labels
    override.aes = list(shape = 21, size = 3, alpha = 1, colour = "black")
  )) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "vertical",         # <- vertical legend
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(6, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  )

# --- 3) make two small “text boxes” for the right columns ---
txt_all <- sprintf("All relapse samples:\nρ=%.2f, p=%s", rho_all, pval_all_str)
txt_pre <- sprintf("Pre-relapse only:\nρ=%.2f, p=%s", rho_pre, pval_pre_str)

mini_box <- function(s) {
  ggplot() +
    annotate("label", x = 0, y = 1, label = s,
             hjust = 0, vjust = 1, size = 3.5,
             label.size = 0, fill = scales::alpha("white", 0.7)) +
    xlim(0,1) + ylim(0,1) +
    theme_void()
}

col2 <- mini_box(txt_all)
col3 <- mini_box(txt_pre)

# --- 4) assemble: plot on top; 3 columns underneath ---
bottom_row <- p_legend_only | col2 | col3
final_plot <- p_main / bottom_row + plot_layout(heights = c(1, 0.22))

# show and save
print(final_plot)
ggsave("Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_blood_muts.png",
       final_plot, width = 5.5, height = 5.5, dpi = 600)













### Now evaluate on non-frontline too 
## Re-introdue dat since previously limited 
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )


## Do rescored 
dat <- dat %>%
  ## Add the screen column 
  mutate(
    BM_zscore_only_detection_rate_screen_call  = as.integer(BM_zscore_only_detection_rate_prob >= 0.350),
  )

################################################################################
##  Time-window prediction performance in Non-frontline cohort
################################################################################

# 1) Your assays vector
assays <- c(
  EasyM        = "EasyM_optimized_binary",
  Flow         = "Flow_Binary",
  cfWGS_BM     = "BM_zscore_only_detection_rate_call",
  cfWGS_Blood  = "Blood_zscore_only_sites_call", 
  cfWGS_Blood_Combined = "Blood_plus_fragment_call"
)

## Get additional dates
Relapse_dates_full <- read_csv(
  "Exported_data_tables_clinical/Relapse dates cfWGS updated2.csv",
  col_types = cols(
    Patient          = col_character(),
    Progression_date = col_date(format = "%Y-%m-%d")
  )
)

# 2) Build df_sf: non-frontline samples + assays + relapse info
df_sf <- dat %>%
  filter(tolower(Cohort) == "non-frontline") %>%
  transmute(
    Patient,
    sample_date = as.Date(Date),
    Flow_Binary,
    BM_zscore_only_detection_rate_call,
    Blood_zscore_only_sites_call,
    Blood_plus_fragment_call
  ) %>%
  # keep rows with at least one assay result
  filter(if_any(all_of(assays), ~ !is.na(.x))) %>%
  # join in per‐patient relapse_date + relapsed flag
  left_join(
    final_tbl %>% 
      transmute(
        Patient,
        relapse_date = as.Date(censor_date),
        relapsed     = as.integer(relapsed)
      ),
    by = "Patient"
  )

## Edit if progression date earlier than sample collected, since some patients had multiple
df_sf2 <- df_sf %>%
  select(-relapse_date, -relapsed) %>%           # ① drop the old pair
  left_join(Relapse_dates_full, by = "Patient") %>%
  filter(Progression_date >= sample_date) %>%
  group_by(Patient, sample_date,
           Flow_Binary,
           BM_zscore_only_detection_rate_call,
           Blood_zscore_only_sites_call) %>%
  slice_min(Progression_date, with_ties = FALSE) %>%
  ungroup() %>%
  rename(relapse_date = Progression_date) %>%    # ② now safe to rename
  mutate(relapsed = as.integer(!is.na(relapse_date)))

# 3) Define windows (days) to evaluate
windows <- c(90, 180, 365, 730)

# 4) Helper to compute metrics for one assay + one window
calc_metrics <- function(df, assay_label, col_name, win_d) {
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[col_name]])) %>%
    dplyr::mutate(
      test_pos        = (.data[[col_name]] == 1),
      event_in_window = relapsed == 1 &
        !is.na(relapse_date) &
        relapse_date <= sample_date + lubridate::days(win_d)   # << here
    )
  
  tab <- table(
    factor(df2$test_pos,        levels = c(FALSE, TRUE)),
    factor(df2$event_in_window, levels = c(FALSE, TRUE))
  )
  
  tp <- tab["TRUE","TRUE"]; fn <- tab["FALSE","TRUE"]
  fp <- tab["TRUE","FALSE"]; tn <- tab["FALSE","FALSE"]
  
  tibble::tibble(
    Window_days = win_d,
    Assay       = assay_label,
    N_samples   = nrow(df2),
    N_patients  = dplyr::n_distinct(df2$Patient),
    TP = tp, FN = fn, FP = fp, TN = tn,
    Sensitivity = if((tp+fn)>0) tp/(tp+fn) else NA_real_,
    Specificity = if((tn+fp)>0) tn/(tn+fp) else NA_real_,
    PPV         = if((tp+fp)>0) tp/(tp+fp) else NA_real_,
    NPV         = if((tn+fn)>0) tn/(tn+fn) else NA_real_
  )
}


# 6) Inspect results
results <- imap_dfr(assays, 
                    # .x = column name, .y = assay label
                    .f = function(col_name, assay_label) {
                      map_dfr(windows, function(win_d) {
                        calc_metrics(df_sf2, assay_label, col_name, win_d)
                      })
                    }
)

results %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity)) %>%
  print(n = Inf)

#### Now do only for those who got the cfWGS test done 
# --- 0) use the refined table everywhere ---------------------------------
bm_col    <- assays["cfWGS_BM"]      # "BM_zscore_only_detection_rate_call"
blood_col <- assays["cfWGS_Blood"]   # "Blood_zscore_only_sites_call"

# Rebuild subsets FROM df_sf2 (not df_sf)
df_sf_BM    <- df_sf2 %>% dplyr::filter(!is.na(.data[[bm_col]]))
df_sf_blood <- df_sf2 %>% dplyr::filter(!is.na(.data[[blood_col]]))

# Recompute results on the consistent base
results_BM <- purrr::map_dfr(windows, function(w) {
  purrr::map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_BM, a, assays[[a]], w)
  })
}) %>% dplyr::arrange(Window_days, dplyr::desc(Sensitivity), dplyr::desc(Specificity))

results_blood <- purrr::map_dfr(windows, function(w) {
  purrr::map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_blood, a, assays[[a]], w)
  })
}) %>% dplyr::arrange(Window_days, dplyr::desc(Sensitivity), dplyr::desc(Specificity))


# --- 1) narrative counts that match the EXACT subsets above ---------------
wins <- c(180, 365)

# Helper: per-window event counts on a given DF (sample-level)
count_relapses_by_window <- function(df, wins = c(180, 365)) {
  purrr::map_dfr(wins, function(w) {
    df %>%
      dplyr::mutate(event_in_window = relapsed == 1 &
                      !is.na(relapse_date) &
                      relapse_date <= sample_date + lubridate::days(w)) %>%
      dplyr::summarise(
        Window_days       = w,
        Samples_relapsed  = sum(event_in_window, na.rm = TRUE),
        Patients_relapsed = dplyr::n_distinct(Patient[event_in_window])
      )
  })
}

# A) BM-only wording (matches results_BM denominators for the cfWGS_BM rows)
summ_bm <- count_relapses_by_window(df_sf_BM, wins)
txt_bm <- glue::glue(
  "Among bone-marrow cfWGS samples, {summ_bm$Samples_relapsed[summ_bm$Window_days==180]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_bm$Samples_relapsed[summ_bm$Window_days==365]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_bm, "\n")

# B) Blood-only wording (matches results_blood denominators for the cfWGS_Blood rows)
summ_blood <- count_relapses_by_window(df_sf_blood, wins)
txt_blood <- glue::glue(
  "Among blood cfWGS samples, {summ_blood$Samples_relapsed[summ_blood$Window_days==180]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_blood$Samples_relapsed[summ_blood$Window_days==365]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_blood, "\n")

# C) “Any cfWGS” wording = union of rows that have BM OR Blood assays
cfwgs_any <- df_sf2 %>% dplyr::filter(!is.na(.data[[bm_col]]) | !is.na(.data[[blood_col]]))
summ_any  <- count_relapses_by_window(cfwgs_any, wins)
txt_any <- glue::glue(
  "In the test cohort, sampling times were heterogeneous, so we assessed each assay’s ",
  "ability to predict progression within fixed time windows (180 and 365 days). ",
  "Among cfWGS samples, {summ_any$Samples_relapsed[summ_any$Window_days==180]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_any$Samples_relapsed[summ_any$Window_days==365]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_any, "\n")



# 5a) Restrict to BM-cfWGS subset, then loop ------------------------------
bm_col     <- assays["cfWGS_BM"]    # "BM_zscore_only_detection_rate_call"
df_sf_BM   <- df_sf2 %>% filter(!is.na(.data[[bm_col]]))

results_BM <- map_dfr(windows, function(w) {
  map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_BM, a, assays[[a]], w)
  })
}) %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity))

# 5b) Restrict to blood-cfWGS subset, then loop ---------------------------
blood_col   <- assays["cfWGS_Blood"]  # "Blood_zscore_only_sites_call"
df_sf_blood <- df_sf2 %>% filter(!is.na(.data[[blood_col]]))

results_blood <- map_dfr(windows, function(w) {
  map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_blood, a, assays[[a]], w)
  })
}) %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity))

# 6) Inspect both ------------------------------------------
message("=== BM-cfWGS subset ===")
print(results_BM)

message("=== Blood-cfWGS subset ===")
print(results_blood)


## Get info
library(scales)   # for percent()

results_BM %>%
  # pick the assays & windows you want to narrate
  filter(Assay %in% c("cfWGS_BM", "Flow"), Window_days %in% c(180, 365)) %>%
  # build a sentence for each row
  rowwise() %>%
  mutate(
    sentence = glue(
      "{Assay} at {Window_days}-day window detected {TP}/{TP + FN} progressors ",
      "(sensitivity {percent(Sensitivity)}, specificity {percent(Specificity)})."
    )
  ) %>%
  ungroup() %>%
  # print them to the console
  pull(sentence) %>%
  cat(sep = "\n")

## For blood
results_blood %>%
  # pick the assays & windows you want to narrate
  filter(Assay %in% c("cfWGS_Blood", "Flow", "cfWGS_Blood_Combined"), Window_days %in% c(180, 365)) %>%
  # build a sentence for each row
  rowwise() %>%
  mutate(
    sentence = glue(
      "{Assay} at {Window_days}-day window detected {TP}/{TP + FN} progressors ",
      "(sensitivity {percent(Sensitivity)}, specificity {percent(Specificity)})."
    )
  ) %>%
  ungroup() %>%
  # print them to the console
  pull(sentence) %>%
  cat(sep = "\n")

## Get event counts 
# define your windows of interest
windows <- c(90, 180, 365, 730)

# this will count, for each window:
# - how many distinct patients relapsed within that window
# - how many samples fall into that window

event_counts_BM <- map_dfr(windows, function(w) {
  df_w <- df_sf_BM %>%
    filter(
      relapsed == 1,
      !is.na(relapse_date),
      relapse_date <= sample_date + days(w)
    )
  
  tibble(
    Window_days = w,
    N_patients  = n_distinct(df_w$Patient),
    N_samples   = nrow(df_w)
  )
})

# view the result
print(event_counts_BM)

## Get number of progressors
wins <- c(180, 365)

# Helper that, given a data frame of samples, returns counts per window
count_relapses_by_window <- function(df, wins = c(180, 365)) {
  purrr::map_dfr(wins, function(w) {
    df %>%
      dplyr::mutate(
        event_in_window = relapsed == 1 &
          !is.na(relapse_date) &
          relapse_date <= sample_date + lubridate::days(w)
      ) %>%
      dplyr::summarise(
        Window_days        = w,
        Samples_relapsed   = sum(event_in_window, na.rm = TRUE),
        Patients_relapsed  = dplyr::n_distinct(Patient[event_in_window])
      )
  })
}

# A) "cfWGS samples" = any sample with BM or Blood cfWGS available
cfwgs_any <- df_sf2 %>%
  dplyr::filter(!is.na(.data[[bm_col]]) | !is.na(.data[[blood_col]]))
summ_any <- count_relapses_by_window(cfwgs_any, wins)

txt_any <- glue::glue(
  "In the test cohort, sampling times were heterogeneous, so we assessed each assay’s ",
  "ability to predict progression within fixed time windows (180 and 365 days). ",
  "Among cfWGS samples, {summ_any$Samples_relapsed[summ_any$Window_days==180]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_any$Samples_relapsed[summ_any$Window_days==365]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_any, "\n")

# (Optional) BM-only wording
cfwgs_bm   <- df_sf2 %>% dplyr::filter(!is.na(.data[[bm_col]]))
summ_bm    <- count_relapses_by_window(cfwgs_bm, wins)
txt_bm <- glue::glue(
  "Among bone-marrow cfWGS samples, {summ_bm$Samples_relapsed[summ_bm$Window_days==180]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_bm$Samples_relapsed[summ_bm$Window_days==365]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_bm, "\n")

# (Optional) Blood-only wording
cfwgs_blood <- df_sf2 %>% dplyr::filter(!is.na(.data[[blood_col]]))
summ_blood  <- count_relapses_by_window(cfwgs_blood, wins)
txt_blood <- glue::glue(
  "Among blood cfWGS samples, {summ_blood$Samples_relapsed[summ_blood$Window_days==180]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_blood$Samples_relapsed[summ_blood$Window_days==365]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==365]} patients relapsed within 365 days."
)
 cat(txt_blood, "\n")



## Make figure 
# ────────────────────────────────────────────────────────────────────────────
# 1) Prepare the data
# ────────────────────────────────────────────────────────────────────────────
sens_BM_df <- results_BM %>%
  # if you want to drop the blood‑only assay, uncomment:
  # filter(Assay != "cfWGS_Blood") %>%
  
  # turn Window_days into a nice factor
  mutate(
    Timepoint = factor(
      Window_days,
      levels = c(90, 180, 365, 730),
      labels = c("90 days", "180 days", "365 days", "730 days")
    ),
    Sens_pct = Sensitivity * 100,
    Assay = recode(
      Assay,
      Flow        = "MFC",
      cfWGS_BM    = "cfWGS",
    )
  )

sens_BM_df <- sens_BM_df %>% filter(Assay != "cfWGS_Blood")
# ────────────────────────────────────────────────────────────────────────────
# 2) Colours & theme (match your existing style)
# ────────────────────────────────────────────────────────────────────────────
custom_cols <- c(
  "90 days"  = "#440154FF",
  "180 days" = "#31688EFF",
  "365 days" = "#35B779FF",
  "730 days" = "#E69F00FF"
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title      = element_text(size = 11),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top",
    plot.margin     = margin(10, 10, 30, 10)
  )

# ────────────────────────────────────────────────────────────────────────────
# 3) Build the grouped bar‑plot
# ────────────────────────────────────────────────────────────────────────────
p_sens_bm <- ggplot(sens_BM_df,
                    aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Window",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays over\nfollow‑up windows (Test Cohort)",
    x     = "Assay",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 14)
  )

# ────────────────────────────────────────────────────────────────────────────
# 4) (Optional) Save
# ────────────────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/Supp_Fig_6_Fig_sensitivity_windows_BM_test_cohort_updated2.png",
       plot = p_sens_bm,
       width = 4.75, height = 6.25, dpi = 500)


### Now remake for blood muts
sens_blood_df <- results_blood %>%
  # if you want to drop the blood‑only assay, uncomment:
  # filter(Assay != "cfWGS_Blood") %>%
  
  # turn Window_days into a nice factor
  mutate(
    Timepoint = factor(
      Window_days,
      levels = c(90, 180, 365, 730),
      labels = c("90 days", "180 days", "365 days", "730 days")
    ),
    Sens_pct = Sensitivity * 100,
    Assay = recode(
      Assay,
      Flow        = "MFC",
      cfWGS_Blood    = "cfWGS",
    )
  )

sens_blood_df <- sens_blood_df %>% filter(Assay != "cfWGS_BM")

p_sens_blood <- ggplot(sens_blood_df,
                       aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Window",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays over\nfollow‑up windows (Test Cohort)",
    x     = "Assay",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ────────────────────────────────────────────────────────────────────────────
# 4) (Optional) Save
# ────────────────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/Supp_Fig_8_Fig_sensitivity_windows_blood_test_cohort3.png",
       plot = p_sens_blood,
       width = 5, height = 5, dpi = 500)





### Export this
# full results (all patients with any assay)
#write_csv(
#  results,
#  file.path(outdir, "all_assays_timewindow_results.csv")
#)

# BM‐cfWGS subset
write_csv(
  results_BM,
  file.path(outdir, "BM_cfWGS_timewindow_results2.csv")
)

# blood‐cfWGS subset
write_csv(
  results_blood,
  file.path(outdir, "blood_cfWGS_timewindow_results2.csv")
)

## export 
# --- Clean BM results ---
results_BM_clean <- results_BM %>%
  # remove blood assays
  filter(!Assay %in% c("cfWGS_Blood", "cfWGS_Blood_Combined")) %>%
  # rename Flow → MFC
  mutate(
    Assay = case_when(
      Assay == "Flow" ~ "MFC",
      TRUE ~ Assay
    )
  )

# --- Clean Blood results ---
results_blood_clean <- results_blood %>%
  # remove BM assays
  filter(!Assay %in% c("cfWGS_BM")) %>%
  # rename assays
  mutate(
    Assay = case_when(
      Assay == "Flow" ~ "MFC",
      Assay == "cfWGS_Blood" ~ "cfWGS_Blood (Sites Model)",
      Assay == "cfWGS_Blood_Combined" ~ "cfWGS_Blood (Combined Model)",
      TRUE ~ Assay
    )
  )

# --- Export to Excel ---
export_list <- list(
  "BM_models"    = results_BM_clean,
  "Blood_models" = results_blood_clean
)

# Write to Excel
write_xlsx(
  export_list,
  path = file.path("Final Tables and Figures/Supplementary_Table_9_timewindow_results_test_cohort.xlsx")
)

# event counts for BM‐cfWGS
write_csv(
  event_counts_BM,
  file.path(outdir, "BM_cfWGS_event_counts2.csv")
)

# For blood
event_counts_blood <- map_dfr(windows, function(w) {
  df_w <- df_sf_blood %>%
    filter(
      relapsed == 1,
      !is.na(relapse_date),
      relapse_date <= sample_date + days(w)
    )
  tibble(
    Window_days = w,
    N_patients  = n_distinct(df_w$Patient),
    N_samples   = nrow(df_w)
  )
})

write_csv(
  event_counts_blood,
  file.path(outdir, "blood_cfWGS_event_counts2.csv")
)




### See at what pont an increase occured 
d2m <- function(days, digits = 1) round(as.numeric(days) / 30.44, digits)

# 1) choose the probability column you want to analyze:
assay_prob <- "BM_zscore_only_detection_rate_prob"  # or "BM_zscore_only_detection_rate_prob"
# Optional: restrict to surveillance timepoints only (post-induction / ASCT / maintenance)
surveillance_only <- FALSE
surv_regex <- "(post[_ -]?induction|pre[_ -]?(asct|transplant)|post[_ -]?(asct|transplant)|maintenance)"

# ==== BASE (relapse patients only, consistent with figure) ====
df0 <- time_to_relapse_BM %>%
  mutate(timepoint_info = tolower(timepoint_info)) %>%
  filter(progress_status == "Relapse",
         !is.na(.data[[assay_prob]]))

if (surveillance_only) {
  df0 <- df0 %>% filter(grepl(surv_regex, timepoint_info))
}

# ==== A) STRICT pre-progression monitoring set ====
# Use days_before_event > 0 to enforce strictly pre-progression (same logic as sample_date < censor_date)
df_pre <- df0 %>%
  filter(days_before_event > 0) %>%
  arrange(Patient, sample_date)

n_patients <- n_distinct(df_pre$Patient)
n_samples  <- nrow(df_pre)

sample_timing_stats <- df_pre %>%
  summarise(
    median_days_before = median(days_before_event, na.rm = TRUE),
    iqr_days_before    = IQR(days_before_event,    na.rm = TRUE),
    min_days_before    = min(days_before_event,    na.rm = TRUE),
    max_days_before    = max(days_before_event,    na.rm = TRUE),
    .groups = "drop"
  )

sample_timing_stats <- sample_timing_stats %>%
  mutate(
    median_days_before_mo = d2m(median_days_before),
    iqr_days_before_mo    = d2m(iqr_days_before),
    min_days_before_mo    = d2m(min_days_before),
    max_days_before_mo    = d2m(max_days_before)
  )

# ==== B) Per-patient nadir & first increase (pre-progression only) ====
advance_df <- df0 %>%
  # Keep rows up to relapse (+30d tolerance) so patients with slight date mismatches are still present
  dplyr::filter(sample_date <= (censor_date + 30)) %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(sample_date, .by_group = TRUE) %>%
  dplyr::group_modify(~{
    df <- .
    prog_date <- df$censor_date[1]
    
    # STRICT pre-relapse rows for nadir / first-increase logic
    df_pre <- dplyr::filter(df, sample_date < prog_date)
    
    if (nrow(df_pre) == 0L) {
      return(tibble::tibble(
        nadir_date = as.Date(NA),  nadir_prob = NA_real_,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    # Nadir = minimum probability before relapse
    nadir_row <- df_pre %>%
      dplyr::slice_min(.data[[assay_prob]], with_ties = FALSE)
    
    nadir_date <- nadir_row$sample_date
    nadir_prob <- nadir_row[[assay_prob]]
    
    # First increase after nadir (strictly later in time AND strictly higher prob)
    post_nadir <- df_pre %>%
      dplyr::filter(sample_date > nadir_date,
                    .data[[assay_prob]] > nadir_prob) %>%
      dplyr::slice_head(n = 1)
    
    if (nrow(post_nadir) == 0L) {
      return(tibble::tibble(
        nadir_date = nadir_date,  nadir_prob = nadir_prob,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase  = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    first_inc_date <- post_nadir$sample_date
    first_inc_prob <- post_nadir[[assay_prob]]
    
    tibble::tibble(
      nadir_date = nadir_date,
      nadir_prob = nadir_prob,
      first_inc_date = first_inc_date,
      first_inc_prob = first_inc_prob,
      days_to_first_increase  = as.numeric(first_inc_date - nadir_date),
      days_before_progression = as.numeric(prog_date - first_inc_date),
      prog_date = prog_date
    )
  }) %>%
  dplyr::ungroup()


# Nadir timing relative to progression
nadir_timing_stats <- advance_df %>%
  mutate(days_nadir_before = as.numeric(prog_date - nadir_date)) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before,    na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before,    na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before,    na.rm = TRUE),
    .groups = "drop"
  )

nadir_timing_stats <- nadir_timing_stats %>%
  mutate(
    median_nadir_mo = d2m(median_nadir_days),
    iqr_nadir_mo    = d2m(iqr_nadir_days),
    min_nadir_mo    = d2m(min_nadir_days),
    max_nadir_mo    = d2m(max_nadir_days)
  )

# First increase timing (after nadir, and its lead time to progression)
increase_stats2 <- advance_df %>%
  summarise(
    median_after_nadir = median(days_to_first_increase,   na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,      na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,      na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,      na.rm = TRUE),
    median_before_prog = median(days_before_progression,  na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,     na.rm = TRUE),
    min_before_prog    = min(days_before_progression,     na.rm = TRUE),
    max_before_prog    = max(days_before_progression,     na.rm = TRUE),
    .groups            = "drop"
  )

increase_stats2 <- increase_stats2 %>%
  mutate(
    median_after_nadir_mo = d2m(median_after_nadir),
    iqr_after_nadir_mo    = d2m(iqr_after_nadir),
    min_after_nadir_mo    = d2m(min_after_nadir),
    max_after_nadir_mo    = d2m(max_after_nadir),
    median_before_prog_mo = d2m(median_before_prog),
    iqr_before_prog_mo    = d2m(iqr_before_prog),
    min_before_prog_mo    = d2m(min_before_prog),
    max_before_prog_mo    = d2m(max_before_prog)
  )

# ==== C) One-paragraph sentence  ====
cat(glue(
  "We analyzed {n_patients} patients (total {n_samples} samples) collected a median ",
  "{sample_timing_stats$median_days_before} days (~{sample_timing_stats$median_days_before_mo} mo) before progression ",
  "(IQR {sample_timing_stats$iqr_days_before} days, ~{sample_timing_stats$iqr_days_before_mo} mo; ",
  "range {sample_timing_stats$min_days_before}–{sample_timing_stats$max_days_before} days, ",
  "~{sample_timing_stats$min_days_before_mo}–{sample_timing_stats$max_days_before_mo} mo). ",
  "The nadir detection probability occurred a median ",
  "{nadir_timing_stats$median_nadir_days} days (~{nadir_timing_stats$median_nadir_mo} mo) before progression ",
  "(IQR {nadir_timing_stats$iqr_nadir_days} days, ~{nadir_timing_stats$iqr_nadir_mo} mo; ",
  "range {nadir_timing_stats$min_nadir_days}–{nadir_timing_stats$max_nadir_days} days, ",
  "~{nadir_timing_stats$min_nadir_mo}–{nadir_timing_stats$max_nadir_mo} mo). ",
  "From that nadir, the first increase occurred a median ",
  "{increase_stats2$median_after_nadir} days (~{increase_stats2$median_after_nadir_mo} mo) later ",
  "(IQR {increase_stats2$iqr_after_nadir} days, ~{increase_stats2$iqr_after_nadir_mo} mo; ",
  "range {increase_stats2$min_after_nadir}–{increase_stats2$max_after_nadir} days, ",
  "~{increase_stats2$min_after_nadir_mo}–{increase_stats2$max_after_nadir_mo} mo), ",
  "which was a median {increase_stats2$median_before_prog} days ",
  "(~{increase_stats2$median_before_prog_mo} mo) before IMWG-defined clinical progression ",
  "(IQR {increase_stats2$iqr_before_prog} days, ~{increase_stats2$iqr_before_prog_mo} mo; ",
  "range {increase_stats2$min_before_prog}–{increase_stats2$max_before_prog} days, ",
  "~{increase_stats2$min_before_prog_mo}–{increase_stats2$max_before_prog_mo} mo).\n"
))



##### Now redo for blood muts 
# 1) choose the probability column you want to analyze:
assay_prob <- "Blood_zscore_only_sites_call"  # or "BM_zscore_only_detection_rate_prob"
# Optional: restrict to surveillance timepoints only (post-induction / ASCT / maintenance)
surveillance_only <- FALSE
surv_regex <- "(post[_ -]?induction|pre[_ -]?(asct|transplant)|post[_ -]?(asct|transplant)|maintenance)"

# ==== BASE (relapse patients only, consistent with figure) ====
df0 <- time_to_relapse_blood %>%
  mutate(timepoint_info = tolower(timepoint_info)) %>%
  filter(progress_status == "Relapse",
         !is.na(.data[[assay_prob]]))

if (surveillance_only) {
  df0 <- df0 %>% filter(grepl(surv_regex, timepoint_info))
}

# ==== A) STRICT pre-progression monitoring set ====
# Use days_before_event > 0 to enforce strictly pre-progression (same logic as sample_date < censor_date)
df_pre <- df0 %>%
  filter(days_before_event > 0) %>%
  arrange(Patient, sample_date)

n_patients <- n_distinct(df_pre$Patient)
n_samples  <- nrow(df_pre)

sample_timing_stats <- df_pre %>%
  summarise(
    median_days_before = median(days_before_event, na.rm = TRUE),
    iqr_days_before    = IQR(days_before_event,    na.rm = TRUE),
    min_days_before    = min(days_before_event,    na.rm = TRUE),
    max_days_before    = max(days_before_event,    na.rm = TRUE),
    .groups = "drop"
  )

sample_timing_stats <- sample_timing_stats %>%
  mutate(
    median_days_before_mo = d2m(median_days_before),
    iqr_days_before_mo    = d2m(iqr_days_before),
    min_days_before_mo    = d2m(min_days_before),
    max_days_before_mo    = d2m(max_days_before)
  )

# ==== B) Per-patient nadir & first increase (pre-progression only) ====
advance_df <- df0 %>%
  # Keep rows up to relapse (+30d tolerance) so patients with slight date mismatches are still present
  dplyr::filter(sample_date <= (censor_date + 30)) %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(sample_date, .by_group = TRUE) %>%
  dplyr::group_modify(~{
    df <- .
    prog_date <- df$censor_date[1]
    
    # STRICT pre-relapse rows for nadir / first-increase logic
    df_pre <- dplyr::filter(df, sample_date < prog_date)
    
    if (nrow(df_pre) == 0L) {
      return(tibble::tibble(
        nadir_date = as.Date(NA),  nadir_prob = NA_real_,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    # Nadir = minimum probability before relapse
    nadir_row <- df_pre %>%
      dplyr::slice_min(.data[[assay_prob]], with_ties = FALSE)
    
    nadir_date <- nadir_row$sample_date
    nadir_prob <- nadir_row[[assay_prob]]
    
    # First increase after nadir (strictly later in time AND strictly higher prob)
    post_nadir <- df_pre %>%
      dplyr::filter(sample_date > nadir_date,
                    .data[[assay_prob]] > nadir_prob) %>%
      dplyr::slice_head(n = 1)
    
    if (nrow(post_nadir) == 0L) {
      return(tibble::tibble(
        nadir_date = nadir_date,  nadir_prob = nadir_prob,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase  = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    first_inc_date <- post_nadir$sample_date
    first_inc_prob <- post_nadir[[assay_prob]]
    
    tibble::tibble(
      nadir_date = nadir_date,
      nadir_prob = nadir_prob,
      first_inc_date = first_inc_date,
      first_inc_prob = first_inc_prob,
      days_to_first_increase  = as.numeric(first_inc_date - nadir_date),
      days_before_progression = as.numeric(prog_date - first_inc_date),
      prog_date = prog_date
    )
  }) %>%
  dplyr::ungroup()


# Nadir timing relative to progression
nadir_timing_stats <- advance_df %>%
  mutate(days_nadir_before = as.numeric(prog_date - nadir_date)) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before,    na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before,    na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before,    na.rm = TRUE),
    .groups = "drop"
  )

nadir_timing_stats <- nadir_timing_stats %>%
  mutate(
    median_nadir_mo = d2m(median_nadir_days),
    iqr_nadir_mo    = d2m(iqr_nadir_days),
    min_nadir_mo    = d2m(min_nadir_days),
    max_nadir_mo    = d2m(max_nadir_days)
  )

# First increase timing (after nadir, and its lead time to progression)
increase_stats2 <- advance_df %>%
  summarise(
    median_after_nadir = median(days_to_first_increase,   na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,      na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,      na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,      na.rm = TRUE),
    median_before_prog = median(days_before_progression,  na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,     na.rm = TRUE),
    min_before_prog    = min(days_before_progression,     na.rm = TRUE),
    max_before_prog    = max(days_before_progression,     na.rm = TRUE),
    .groups            = "drop"
  )

increase_stats2 <- increase_stats2 %>%
  mutate(
    median_after_nadir_mo = d2m(median_after_nadir),
    iqr_after_nadir_mo    = d2m(iqr_after_nadir),
    min_after_nadir_mo    = d2m(min_after_nadir),
    max_after_nadir_mo    = d2m(max_after_nadir),
    median_before_prog_mo = d2m(median_before_prog),
    iqr_before_prog_mo    = d2m(iqr_before_prog),
    min_before_prog_mo    = d2m(min_before_prog),
    max_before_prog_mo    = d2m(max_before_prog)
  )

# ==== C) One-paragraph sentence  ====
cat(glue(
  "We analyzed {n_patients} patients (total {n_samples} samples) collected a median ",
  "{sample_timing_stats$median_days_before} days (~{sample_timing_stats$median_days_before_mo} mo) before progression ",
  "(IQR {sample_timing_stats$iqr_days_before} days, ~{sample_timing_stats$iqr_days_before_mo} mo; ",
  "range {sample_timing_stats$min_days_before}–{sample_timing_stats$max_days_before} days, ",
  "~{sample_timing_stats$min_days_before_mo}–{sample_timing_stats$max_days_before_mo} mo). ",
  "The nadir detection probability occurred a median ",
  "{nadir_timing_stats$median_nadir_days} days (~{nadir_timing_stats$median_nadir_mo} mo) before progression ",
  "(IQR {nadir_timing_stats$iqr_nadir_days} days, ~{nadir_timing_stats$iqr_nadir_mo} mo; ",
  "range {nadir_timing_stats$min_nadir_days}–{nadir_timing_stats$max_nadir_days} days, ",
  "~{nadir_timing_stats$min_nadir_mo}–{nadir_timing_stats$max_nadir_mo} mo). ",
  "From that nadir, the first increase occurred a median ",
  "{increase_stats2$median_after_nadir} days (~{increase_stats2$median_after_nadir_mo} mo) later ",
  "(IQR {increase_stats2$iqr_after_nadir} days, ~{increase_stats2$iqr_after_nadir_mo} mo; ",
  "range {increase_stats2$min_after_nadir}–{increase_stats2$max_after_nadir} days, ",
  "~{increase_stats2$min_after_nadir_mo}–{increase_stats2$max_after_nadir_mo} mo), ",
  "which was a median {increase_stats2$median_before_prog} days ",
  "(~{increase_stats2$median_before_prog_mo} mo) before IMWG-defined clinical progression ",
  "(IQR {increase_stats2$iqr_before_prog} days, ~{increase_stats2$iqr_before_prog_mo} mo; ",
  "range {increase_stats2$min_before_prog}–{increase_stats2$max_before_prog} days, ",
  "~{increase_stats2$min_before_prog_mo}–{increase_stats2$max_before_prog_mo} mo).\n"
))










### Maybe follow up with Esteban to see if clinical or biochemical progression since don't know for non-frontline
### Leave as blank for now
# Filter out patients whose ID starts with "IMG-" or "SPORE-"

# ═════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY: Script Output and Results
# ═════════════════════════════════════════════════════════════════════════════
#
# This comprehensive survival analysis script has produced the following:
#
# ─────────────────────────────────────────────────────────────────────────────
# PRIMARY OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────
#
# 1. KAPLAN-MEIER SURVIVAL CURVES (PNG files @ 500 DPI)
#    Location: detection_progression_updated6/[timepoint]/
#    Files: KM_[Assay]_[Timepoint]_updated_no_CI.png
#    - Diagnosis: 4 assays (EasyM, clonoSEQ, MFC, cfWGS)
#    - Post-ASCT: 7 assays (+ BM/blood variants)
#    - Maintenance-1yr: 7 assays (+ BM/blood variants)
#    Each plot shows: PFS curves, risk table, log-rank p-value
#
# 2. SENSITIVITY TABLES (CSV files)
#    Location: Output directory root
#    - frontline_postASCT_sensitivity.csv: % MRD+ among relapsed at post-ASCT
#    - frontline_1yr_sensitivity.csv: % MRD+ among relapsed at 1-year
#    - frontline_postASCT_sens_BMcfWGS.csv: Sensitivity in BM-cfWGS subset
#    - frontline_1yr_sens_BMcfWGS.csv: Sensitivity in BM-cfWGS subset
#    - frontline_postASCT_sens_bloodcfWGS.csv: Sensitivity in blood-cfWGS subset
#    - frontline_1yr_sens_bloodcfWGS.csv: Sensitivity in blood-cfWGS subset
#
# 3. SENSITIVITY BARPLOTS (PNG files @ 500 DPI)
#    Location: Final Tables and Figures/
#    - Supp_6A_Fig_sensitivity_by_tech_training3.png (BM-subset analysis)
#    - Supp_8A_Fig_sensitivity_by_tech_training_blood2.png (blood-subset analysis)
#    Each plot shows: Grouped bars comparing assay sensitivities at post-ASCT vs 1-year
#
# ─────────────────────────────────────────────────────────────────────────────
# SECONDARY ANALYSES
# ─────────────────────────────────────────────────────────────────────────────
#
# 4. NON-FRONTLINE COHORT TIME-WINDOW ANALYSES
#    - 24-month RFS by cfWGS BM at 1-year maintenance
#    - Median RFS and hazard ratios by assay
#    - Comparable metrics for Flow, clonoSEQ, EasyM
#
# 5. SPEARMAN CORRELATION ANALYSIS
#    - Pairwise correlations between assays
#    - Both full cohort and subset analyses
#
# 6. NARRATIVE SUMMARY STATISTICS
#    - Patient/sample counts and demographics
#    - Timing analysis: Days from baseline to collection/relapse
#    - Nadir timing and relapse progression timing
#    - Formatted paragraphs suitable for methods/results sections
#
# ─────────────────────────────────────────────────────────────────────────────
# KEY STATISTICAL FINDINGS EXPORTED FOR VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────
#
# - PFS probability at each timepoint (KM curves)
# - Log-rank p-values testing MRD+/- difference
# - Hazard ratios with 95% confidence intervals
# - Sensitivity metrics (% detection among relapsers)
# - Time-to-event distributions (days and months)
#
# ─────────────────────────────────────────────────────────────────────────────
# NEXT STEPS
# ─────────────────────────────────────────────────────────────────────────────
#
# After running this script:
# 1. Review KM curves in output folders for publication-ready figures
# 2. Check sensitivity tables for clinical assay performance comparison
# 3. Examine barplots for supplementary figures
# 4. Use narrative statistics for manuscript methods/results sections
# 5. Consult statistical summaries for hypothesis testing and effect sizes
#
# For questions or modifications, reference upstream scripts:
#   - 3_1_Optimize_cfWGS_thresholds.R (cfWGS optimization)
#   - 3_1_A_Process_and_optimize_EasyM.R (EasyM model)
#   - 2_0_Assemble_Table_With_All_Features.R (feature integration)
#
# ═════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("═" %*% 80, "\n")
cat("SURVIVAL ANALYSIS COMPLETE\n")
cat("Scripts completed: All KM curves, sensitivity analyses, and barplots generated\n")
cat("═" %*% 80, "\n\n")


