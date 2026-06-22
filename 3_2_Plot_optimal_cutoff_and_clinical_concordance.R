# =============================================================================
# Script: 3_2_Plot_optimal_cutoff_and_clinical_concordance.R
# Project:  cfWGS MRD detection (M4 / SPORE / IMMAGINE)
# How to run:
#   Rscript Scripts_2025/Final_Scripts/3_2_Plot_optimal_cutoff_and_clinical_concordance.R
#
# Manuscript outputs created/updated:
#   - Figure 3D-E: BM cfWGS clinical-concordance and probability panels.
#   - Figure 4C-D: blood cfWGS clinical-concordance and probability panels.
#   - Extended Data Figure 5E-G: BM supplementary concordance/diagnostic panels.
#   - Extended Data Figure 7F-H: blood supplementary concordance/diagnostic panels.
#   - Supplementary Table 8: BM model concordance/source-data workbook.
#   - Supplementary Table 10: blood model concordance/source-data workbook.
#
# Pipeline role:
#   This script starts from the preserved model calls/probabilities created in
#   3_1_Optimize_cfWGS_thresholds.R and compares them with clinical MRD assays
#   and EasyM. It is the main script for interpreting model calls against
#   independent clinical measurements, not for changing model weights.
#
# Author:   Dory Abelman
# Date:     May 28, 2025
# Updated Feb 2026
#
# Script purpose:
#   1. Read the processed cfWGS call/probability table produced by
#      3_1_Optimize_cfWGS_thresholds.R.
#   2. Join EasyM MRD measurements where available.
#   3. Compare BM-informed and blood/cfDNA-informed cfWGS calls with clinical
#      MRD comparators (MFC, clonoSEQ, and EasyM).
#   4. Generate the manuscript-facing concordance plots, confusion-matrix
#      panels, and supplementary concordance workbooks listed above.
#
# Execution note:
#   The historical working filenames in this script do not always match final
#   manuscript numbering. The authoritative mapping is the `ms_copy_artifact()`
#   block beside each final output. Those blocks copy the final figures/tables
#   into final_manuscript_objects/ using names such as Figure_3D,
#   Figure_4D, Extended_Data_Figure_7F, and Supplementary_Table_10.
#
# Inputs:
#   • Output_tables_2025/all_patients_with_BM_and_blood_calls_updated5.rds
#     from 3_1_Optimize_cfWGS_thresholds.R.
#   • Output_EasyM_MRD_analysis_2025/EasyM_all_samples_with_optimized_calls.csv
#     from 3_1_A_Process_and_optimize_EasyM.R, when EasyM comparison panels are
#     available.
#   • Optional preserved model/threshold RDS files in Output_tables_2025/ for
#     model-definition labels and benchmarking; this script does not refit
#     models.
#
# Outputs:
#   Historical working outputs are written to Output_figures_2025/,
#   Output_tables_2025/, and Final Tables and Figures/. Final manuscript
#   outputs are copied into final_manuscript_objects/ by artifact ID.
#
# =============================================================================
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

# ===========================================================================
# SECTION 0: SETUP
# ===========================================================================
# Load required packages for data wrangling, visualization, and analysis

library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping (pivot_*, separate_*, etc.)
library(ggplot2)     # Publication-quality graphics
library(pROC)        # ROC curves and AUC calculation
library(patchwork)   # Multi-panel figure assembly
library(janitor)     # Quick crosstabs (tabyl) and data cleaning
library(gt)          # Beautiful HTML/image tables for reports
library(glue)        # String interpolation for informative messages
library(rmda)        # Decision curve analysis (DCA) for clinical utility
library(lubridate)   # Date/time manipulation
library(scales)      # Formatting helpers (percent_format,, etc.)
library(readr)       # Fast CSV import/export
library(stringr)     # String normalization for timepoint labels and sample IDs
library(tibble)      # column_to_rownames for confusion matrix helpers
library(broom)       # Tidy summaries for logistic model outputs
library(purrr)       # Functional iteration over model/list columns

# Shared helper for final manuscript-organized outputs.
# The script keeps its historical output filenames while copying final
# manuscript-facing figure/table components into
# Scripts_2025/Final_Scripts/final_manuscript_objects with final labels such as
# Figure_3D, Figure_4D, Extended_Data_Figure_5E, and Supplementary_Table_10.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)


# ===========================================================================
# SECTION 1: DATA INPUT AND IDENTIFICATION OF SOURCE FILES
# ===========================================================================
# This script integrates cfWGS results with EasyM MRD measurements
# Output directory structure from prior analyses

outdir <- "Output_tables_2025"
OUTPUT_DIR_FIGURES <- "Output_figures_2025"
OUTPUT_DIR_EASYM <- "Output_EasyM_MRD_analysis_2025"  # Output from script 3_1_A

# ─────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT SETUP
# ─────────────────────────────────────────────────────────────────────────
# Create subdirectory for figure source data (filtered dataframes & metadata)
outdir_source_data <- "Output_tables_2025/Source_data"
if (!dir.exists(outdir_source_data)) {
  dir.create(outdir_source_data, showWarnings = FALSE)
  cat(sprintf("✓ Created source data directory: %s\n", outdir_source_data))
}

# Version stamp for consistent file naming across tables and figures
version_date <- as.character(Sys.Date())

# ===========================================================================
# SECTION 1A: LOAD EXISTING cfWGS DATA
# ===========================================================================
# Read previously processed cfWGS dataset with binary calls and probabilities

cat("\nLoading cfWGS data from prior analysis...\n")

dat <- readRDS(file.path(outdir, "all_patients_with_BM_and_blood_calls_updated5.rds"))
cat(sprintf("✓ Loaded %d rows from cfWGS dataset\n", nrow(dat)))

# Optional model and threshold definitions used only for labels/benchmarking.
# Keep these paths project-relative for Code Ocean portability. The scored calls
# and probabilities are read from all_patients_with_BM_and_blood_calls_updated5.rds;
# this script does not retrain the models.
PATH_MODEL_LIST <- file.path(outdir, "selected_combo_models_2025-09-17.rds")
PATH_THRESHOLD_LIST <- file.path(outdir, "selected_combo_thresholds_2025-09-17.rds")

if (file.exists(PATH_MODEL_LIST)) {
  selected_models <- readRDS(PATH_MODEL_LIST)
  cat("✓ Loaded cfWGS model definitions\n")
} else {
  selected_models <- NULL
  cat("⚠ cfWGS model file not found, proceeding without model details\n")
}

if (file.exists(PATH_THRESHOLD_LIST)) {
  selected_thr <- readRDS(PATH_THRESHOLD_LIST)
  cat("✓ Loaded cfWGS threshold definitions\n")
} else {
  selected_thr <- NULL
  cat("⚠ cfWGS threshold file not found, proceeding without threshold details\n")
}

# ===========================================================================
# SECTION 1B: LOAD EasyM DATA FROM SCRIPT 3_1_A
# ===========================================================================
# 3_1_A produces the following key files:
#   1. EasyM_all_samples_with_optimized_calls.csv
#      Contains: Patient, Timepoint, EasyM_value, EasyM_clinician_binary,
#                EasyM_optimized_binary, EasyM_optimized_call, threshold_method
#   2. EasyM_threshold_values_by_timepoint.csv
#      Contains: Timepoint and corresponding threshold cutoff (%) used

cat("\nLoading EasyM data from script 3_1_A...\n")

EasyM_file <- file.path(OUTPUT_DIR_EASYM, "EasyM_all_samples_with_optimized_calls.csv")
EasyM_thresholds <- NULL

if (file.exists(EasyM_file)) {
  EasyM_data <- readr::read_csv(EasyM_file, show_col_types = FALSE)
  cat(sprintf("✓ Loaded EasyM data: %d samples across %d timepoints\n", 
              n_distinct(EasyM_data$Patient), n_distinct(EasyM_data$Timepoint)))
  
  # Merge EasyM data into cfWGS dataset by matching Patient and Timepoint
  dat <- dat %>%
    mutate(Patient = as.character(Patient), Timepoint = as.character(Timepoint)) %>%
    left_join(
      EasyM_data %>%
        select(
          Patient, Timepoint, 
          EasyM_value,
          EasyM_clinician_binary,
          EasyM_optimized_binary,
          EasyM_optimized_call,
          threshold_method
        ) %>%
        mutate(Patient = as.character(Patient), Timepoint = as.character(Timepoint)),
      by = c("Patient", "Timepoint"),
      relationship = "many-to-one"
    )
  
  n_easy_m_matched <- sum(!is.na(dat$EasyM_value))
  cat(sprintf("✓ Merged EasyM data: %d samples matched with cfWGS\n", n_easy_m_matched))
  cat(sprintf("  - EasyM_clinician_binary: %d samples\n", sum(!is.na(dat$EasyM_clinician_binary))))
  cat(sprintf("  - EasyM_optimized_binary: %d samples\n", sum(!is.na(dat$EasyM_optimized_binary))))
  cat(sprintf("  - Using optimized thresholds: %d samples\n", 
              sum(dat$threshold_method == "optimized", na.rm = TRUE)))
  cat(sprintf("  - Using clinician calls: %d samples\n", 
              sum(dat$threshold_method == "clinician", na.rm = TRUE)))
  
  # Load threshold reference table
  threshold_ref_file <- file.path(OUTPUT_DIR_EASYM, "tables", "EasyM_threshold_values_by_timepoint.csv")
  if (file.exists(threshold_ref_file)) {
    EasyM_thresholds <- readr::read_csv(threshold_ref_file, show_col_types = FALSE)
    cat(sprintf("✓ Loaded threshold reference: %d timepoints\n", nrow(EasyM_thresholds)))
  }
  
} else {
  cat(sprintf("⚠ EasyM file not found at: %s\n", EasyM_file))
  cat("   Creating placeholder columns to prevent downstream errors...\n")
  
  # Create placeholder columns so downstream code doesn't break
  dat <- dat %>%
    mutate(
      EasyM_value = NA_real_,
      EasyM_clinician_binary = NA_integer_,
      EasyM_optimized_binary = NA_integer_,
      EasyM_optimized_call = NA_character_,
      threshold_method = NA_character_
    )
  
  EasyM_thresholds <- NULL
}

cat("\n")

# ===========================================================================
# SECTION 2: DATA PREPARATION AND HARMONIZATION
# ===========================================================================
# Prepare combined cfWGS + EasyM dataset for visualization and analysis

# Helper function: Calculate positivity rates
# Used to summarize MRD+ frequencies within subgroups
summarise_pos <- function(df) {
  df %>%
    summarise(
      n_total  = n(),
      n_pos    = sum(cfWGS_BM_Binary == 1, na.rm = TRUE),
      pos_rate = n_pos / n_total
    )
}

# Harmonize timepoint labels for consistent visualization
# (standardize naming across different assays)
dat <- dat %>%
  mutate(
    landmark_tp = case_when(
      timepoint_info == "Post_transplant"  ~ "Post-ASCT",
      timepoint_info == "1yr maintenance"  ~ "Maintenance-1yr",
      Timepoint == "07"                    ~ "Maintenance-1yr",
      TRUE                                 ~ NA_character_
    )
  )

# ---------------------------------------------------------------------------
#  3.  FRONTLINE cohort: positivity by landmark - BM-derived mutations -------
# Figure 3D input table: positivity rates for cfWGS, clinical MRD assays, and
# EasyM at the frontline landmark timepoints.
front_tbl <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    !is.na(BM_zscore_only_detection_rate_call)
  ) %>%
  ## Compare both EasyM methods (any-detect and optimized) with other modalities
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_zscore_only_detection_rate_call, 
                  EasyM_clinician_binary, EasyM_optimized_binary),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  dplyr::summarise(
    n_total  = dplyr::n(),
    n_pos    = sum(Result == 1L, na.rm = TRUE),
    pos_rate = n_pos / n_total,
    .groups  = "drop"
  ) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary         = "MFC",
      Adaptive_Binary     = "clonoSEQ",
      BM_zscore_only_detection_rate_call = "cfWGS",
      EasyM_clinician_binary = "EasyM (any)",
      EasyM_optimized_binary = "EasyM (opt)"
    )
  )

# Analyst note: quantify patients who are negative at both frontline landmark
# timepoints by clinical assays and by BM-informed cfWGS. These values support
# manuscript text but are not separately staged as a figure/table artifact.
front_long <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    landmark_tp %in% c("Post-ASCT", "Maintenance-1yr")
  ) %>%
  filter(!is.na(BM_zscore_only_detection_rate_call)) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_zscore_only_detection_rate_call, 
                  EasyM_clinician_binary, EasyM_optimized_binary),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary                     = "MFC",
      Adaptive_Binary                 = "clonoSEQ",
      BM_zscore_only_detection_rate_call = "cfWGS",
      EasyM_clinician_binary          = "EasyM (any)",
      EasyM_optimized_binary          = "EasyM (opt)"
    )
  )

# Find patients negative at both timepoints, per method
neg_both_tbl <- front_long %>%
  group_by(Technology, Patient) %>%
  summarise(
    n_tp = n_distinct(landmark_tp),
    all_neg = (n_tp == 2 & all(Result == 0L)),
    .groups = "drop"
  ) %>%
  filter(all_neg)

# Pull method-specific patient sets for overlap summaries.
neg_sets <- split(neg_both_tbl$Patient, neg_both_tbl$Technology)

mfc_neg   <- neg_sets[["MFC"]]
seq_neg   <- neg_sets[["clonoSEQ"]]
cfwgs_neg <- neg_sets[["cfWGS"]]

# Overlaps between clinical-assay-negative and cfWGS-negative sets.
mfc_also_cfwgs <- intersect(mfc_neg, cfwgs_neg)
seq_also_cfwgs <- intersect(seq_neg, cfwgs_neg)

# Build analyst-facing prose snippets for manuscript text checks.
mfc_sentence <- sprintf(
  "Of %d patients that were negative at both timepoints by MFC, %d (%.1f%%) were also negative by cfWGS.",
  length(mfc_neg),
  length(mfc_also_cfwgs),
  100 * length(mfc_also_cfwgs) / max(1, length(mfc_neg))
)

seq_sentence <- sprintf(
  "Of %d patients that were negative at both timepoints by clonoSEQ, %d (%.1f%%) were also negative by cfWGS.",
  length(seq_neg),
  length(seq_also_cfwgs),
  100 * length(seq_also_cfwgs) / max(1, length(seq_neg))
)

mfc_sentence
seq_sentence

# Analyst note: paired-assay availability counts used to verify the sample
# denominators that appear in concordance text.
tp_labels <- c(
  "Post-ASCT" = "post-ASCT",
  "Maintenance-1yr" = "1-year maintenance"
)

front_counts <- dat %>%
  filter(Cohort == "Frontline", !is.na(landmark_tp)) %>%
  mutate(
    has_screen   = !is.na(BM_zscore_only_detection_rate_call),
    has_clinical = !is.na(Flow_Binary) | !is.na(Adaptive_Binary)
  ) %>%
  group_by(landmark_tp) %>%
  summarise(
    n_total_cfWGS_screen = sum(has_screen, na.rm = TRUE),
    n_with_clinical      = sum(has_screen & has_clinical, na.rm = TRUE),
    frac_string          = paste0(n_with_clinical, "/", n_total_cfWGS_screen),
    .groups = "drop"
  )

concord_sentence <- front_counts %>%
  mutate(tp_label = recode(landmark_tp, !!!tp_labels)) %>%
  arrange(match(tp_label, c("post-ASCT", "1-year maintenance"))) %>%
  summarise(
    sentence = glue(
      "We assessed concordance between cfDNA-based MRD and clinical assays at {glue_collapse(glue('{tp_label} in samples with at least one clinical MRD result (n = {n_with_clinical}/{n_total_cfWGS_screen})'), sep = ' and ')}."
    )
  ) %>%
  pull(sentence)

concord_sentence

# ---------------------------------------------------------------------------
#  4.  NON-FRONTLINE cohort: pooled positivity -------------------------------
# Figure 3D test-cohort input table: positivity rates pooled across eligible
# non-frontline follow-up samples.
non_tbl <- dat %>%
  mutate(landmark_tp = "All timepoints") %>%
  filter(!timepoint_info %in% c("Baseline", "Diagnosis")) %>% 
  filter(
    Cohort == "Non-frontline",
    !is.na(BM_zscore_only_detection_rate_call),
    !is.na(MRD_truth) # restrict to only ones with MRD for fair comparison
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_zscore_only_detection_rate_call),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  dplyr::summarise(
    n_total  = dplyr::n(),
    n_pos    = sum(Result == 1L, na.rm = TRUE),
    pos_rate = n_pos / n_total,
    .groups  = "drop"
  ) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary         = "MFC",
      Adaptive_Binary     = "clonoSEQ",
      BM_zscore_only_detection_rate_call = "cfWGS"
  )
  )


# Export source tables used by the BM-informed positivity panels.
readr::write_csv(
  front_tbl,
  file.path(outdir, "Positivity_by_Landmark_TimePoint_BoneMarrow_Frontline_updated5.csv")
)

readr::write_csv(
  non_tbl,
  file.path(outdir, "Positivity_All_TimePoints_BoneMarrow_NonFrontline_updated5.csv")
)


# Sanity check: count complete versus partially paired assay availability.
# The figure keeps samples with cfWGS plus at least one clinical MRD comparator,
# matching the manuscript denominator rather than requiring all three assays.
dat %>%
  filter(Cohort == "Frontline", !is.na(landmark_tp)) %>%
  group_by(landmark_tp) %>%
  summarise(
    total            = n_distinct(Sample_Code),
    at_least_one     = sum(!is.na(BM_zscore_only_detection_rate_call) | !is.na(Flow_Binary) | !is.na(Adaptive_Binary)),
    all_three        = sum(!is.na(BM_zscore_only_detection_rate_call) & !is.na(Flow_Binary) & !is.na(Adaptive_Binary)),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#  0.  Aesthetics ------------------------------------------------------------
custom_cols <- c("Post-ASCT"       = "#E41A1C",
                 "Maintenance-1yr" = "#377EB8",   
                 "All timepoints" = "#999999")   # grey (non‑frontline single bar)

plot_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title.x  = element_text(size = 11),
    axis.title.y  = element_text(size = 11),
    plot.title    = element_text(hjust = .5, face = "bold", size = 12),
    axis.line     = element_line(colour = "black"),
    panel.grid    = element_blank(),      # remove gridlines
    panel.background = element_blank(),
    legend.position  = "top",
    plot.margin      = margin(10, 10, 30, 10)      # t, r, b, l  (pt)
  )

# ===========================================================================
# SECTION 3: FIGURE GENERATION - POSITIVITY RATE COMPARISONS BY TIMEPOINT
# ===========================================================================
# Compare MRD+ detection rates across cfWGS, EasyM clinician, and EasyM optimized
# Stratified by treatment phase (post-ASCT vs maintenance)

# ---------------------------------------------------------------------------
#  SUBSECTION 3.1: Frontline/Post-ASCT Positivity Rates
# ---------------------------------------------------------------------------
# Compare MRD+ detection at post-ASCT landmark between methods

# Clean and re-factor timepoint column
front_tbl <- front_tbl %>%
  # Normalize unicode hyphen to ASCII hyphen-minus for consistent rendering
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) %>%
  # Convert to factor with standardized order
  mutate(
    landmark_tp = factor(landmark_tp,
                         levels = c("Post-ASCT", "Maintenance-1yr"))
  ) %>%
  # Reorder Technology for consistent display order across all plots
  mutate(
    Technology = factor(Technology,
                        levels = c("cfWGS", "clonoSEQ", "MFC", "EasyM (opt)", "EasyM (any)"))
  )

# Build grouped barplot comparing MRD+ rates across technologies and timepoints
p_front_grouped <- ggplot(front_tbl, 
                          aes(x    = Technology,
                              y    = pos_rate * 100,
                              fill = landmark_tp)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black") +
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)    # forces legend order
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.08)),
    limits = c(0, 100)
  ) +
  labs(
    title = "cfWGS Positivity: Frontline cohort, BM-derived muts",
    x     = NULL,
    y     = "Positivity rate"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(aes(label = sprintf("%.0f%%", pos_rate * 100)),
            position = position_dodge(width = 0.8),
            vjust    = -0.4,
            size     = 3.5)

# Save the first standalone frontline BM positivity plot for provenance.
ggsave(
    filename = file.path(OUTPUT_DIR_FIGURES, "Fig_BM_positivity_by_tech_updated5.png"),
  plot     = p_front_grouped,
  width    = 6,
  height   = 4,
  dpi      = 500
)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: BM positivity by technology (frontline)
# ══════════════════════════════════════════════════════════════════════════
readr::write_csv(
  front_tbl %>% mutate(Figure = "Fig_BM_positivity_by_tech_updated5"),
  file.path(outdir_source_data, "Fig_BM_positivity_by_tech_updated5_source_data.csv")
)


# Intermediate BM positivity plot retained for historical/provenance output.
# The final manuscript-staged BM positivity component is the faceted Figure 3D
# panel built below.
custom_cols <- c(
  "Post-ASCT"       = "#31688E",  # deep teal-blue leaning toward purple
  "Maintenance-1yr" = "#35B779",  # bright forest green
  "All timepoints"  = "#999999"   # neutral grey for the single-bar group
)

p_front_grouped <- ggplot(front_tbl, 
                          aes(x    = Technology,
                              y    = pos_rate * 100,
                              fill = landmark_tp)) +
  # ① solid bars with black outline
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  
  # Manual fill to match the manuscript color palette.
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  
  # ③ y axis as percent, 0–100, with adequate spacing for 100% label
  scale_y_continuous(
    labels = percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.08)),
    limits = c(0, 100)
  ) +
  
  # ④ labels on top of bars
  geom_text(aes(label = sprintf("%d%%", round(pos_rate * 100))),
            position = position_dodge(width = 0.8),
            vjust    = -0.4,
            size     = 3.5,
            family   = "sans") +
  
  # ⑤ titles
  labs(
    title = "MRD Positivity by Technology (Training Cohort)",
    x     = "Technology",
    y     = "Positivity rate"
  ) +
  plot_theme

ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_updated5.png"),
  plot     = p_front_grouped,
  width    = 6,
  height   = 4,
  dpi      = 500
)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Training cohort BM positivity (updated5, same as Fig_updated5)
# ══════════════════════════════════════════════════════════════════════════
# (Already exported above; duplicates omitted)


# ---------------------------------------------------------------------------
#  2. Build non-frontline/test-cohort BM positivity plot ----------------------

# Clean and re-factor the timepoint column for consistent plotting.
non_tbl <- non_tbl %>%
  # 1) normalize that unicode hyphen to the ASCII hyphen-minus
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) %>%
  # Reorder Technology for consistent display order across all plots
  mutate(
    Technology = factor(Technology,
                        levels = c("cfWGS", "clonoSEQ", "MFC", "EasyM (opt)", "EasyM (any)"))
  )

# Build a standalone non-frontline/test-cohort plot for provenance.
p_non_grouped <- ggplot(non_tbl, 
                          aes(x    = Technology,
                              y    = pos_rate * 100,
                              fill = landmark_tp)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black") +
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)    # forces legend order
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.08)),
    limits = c(0, 100)
  ) +
  labs(
    title = "cfWGS Positivity: Later-line cohort, BM-derived muts",
    x     = NULL,
    y     = "Positivity rate"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(aes(label = sprintf("%.0f%%", pos_rate * 100)),
            position = position_dodge(width = 0.8),
            vjust    = -0.4,
            size     = 3.5)

# Save the first standalone non-frontline/test-cohort BM positivity plot for
# provenance.
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_BM_positivity_by_tech_later_line5.png"),
  plot     = p_non_grouped,
  width    = 6,
  height   = 4,
  dpi      = 500
)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Non-frontline BM positivity
# ══════════════════════════════════════════════════════════════════════════
readr::write_csv(
  non_tbl %>% mutate(Figure = "Fig_BM_positivity_by_tech_later_line5"),
  file.path(outdir_source_data, "Fig_BM_positivity_by_tech_later_line5_source_data.csv")
)


# Intermediate themed non-frontline/test-cohort plot retained for provenance.
p_non_grouped <- ggplot(non_tbl, 
                          aes(x    = Technology,
                              y    = pos_rate * 100,
                              fill = landmark_tp)) +
  # ① solid bars with black outline
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  
  # Manual fill to match the manuscript color palette.
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  
  # ③ y axis as percent, 0–100, with adequate spacing for 100% label
  scale_y_continuous(
    labels = percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.08)),
    limits = c(0, 100)
  ) +
  
  # ④ labels on top of bars
  geom_text(aes(label = sprintf("%d%%", round(pos_rate * 100))),
            position = position_dodge(width = 0.8),
            vjust    = -0.4,
            size     = 3.5,
            family   = "sans") +
  
  # ⑤ titles
  labs(
    title = "MRD Positivity by Technology (Test Cohort)",
    x     = "Technology",
    y     = "Positivity rate"
  ) +
  plot_theme

ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_updated_non_frontline6.png"),
  plot     = p_non_grouped,
  width    = 3.5,
  height   = 4,
  dpi      = 500
)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Non-frontline BM positivity (test cohort version)
# ══════════════════════════════════════════════════════════════════════════
# (Already exported above; same data as Fig_BM_positivity_by_tech_later_line5)


# Main Figure 3D: combine training/frontline and test/non-frontline tables into
# one faceted BM-informed positivity panel.
front_tbl2 <- front_tbl  %>% mutate(Cohort = "Training Cohort")
non_tbl2   <- non_tbl    %>% mutate(Cohort = "Test Cohort")
combo_tbl  <- bind_rows(front_tbl2, non_tbl2)

# Set display order for cohorts, timepoints, and technologies.
combo_tbl <- combo_tbl %>%
  mutate(
    Cohort      = factor(Cohort,    levels = c("Training Cohort","Test Cohort")),
    landmark_tp = factor(landmark_tp,
                         levels = c("Post-ASCT",
                                    "Maintenance-1yr",
                                    "All timepoints")),
    Technology = factor(Technology,
                        levels = c("cfWGS", "clonoSEQ", "MFC", "EasyM (opt)", "EasyM (any)"))
  )

# Single faceted plot used as the final Figure 3D component.
p_pos_by_tech <- ggplot(combo_tbl, 
                        aes(x    = Technology,
                            y    = pos_rate * 100,
                            fill = landmark_tp)) +
  # ① solid bars with black outline
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  
  # ② labels on top of bars
  geom_text(aes(label = sprintf("%d%%", round(pos_rate * 100))),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 2.5,
            family   = "sans") +
  
  # Manual fill to match the manuscript color palette.
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  
  # ④ free‐x per facet so “clonoSEQ” drops out of Test, and percent y‐axis
  facet_wrap(~ Cohort, nrow = 1, scales = "free_x") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, 100)
  ) +
  
  # ⑤ titles
  labs(
    title = "MRD Positivity by Technology (BM-Derived mutation lists)",
    x     = "Technology",
    y     = "Positivity Rate"
  ) +
  
  # ⑥ clean up the theme to match p_perf
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y     = element_text(size = 9),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "top"
  )

# Save the final Figure 3D working PNG before staging it to the manuscript
# output directory.
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_facet7.png"),
  plot     = p_pos_by_tech,
  width    = 7.5,    # wider to accommodate two facets
  height   = 4.5,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 3D
#
# What this is:
#   BM-informed cfWGS positivity by clinical-assay technology and cohort.
#
# Why it is here:
#   This PNG is the final Main Figure 3D component. The source-data filename
#   below still says facet6, while the final plotted component is facet7.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_facet7.png"),
  artifact_id = "FIG3D",
  role = "figure_panel_png",
  description = "BM-informed cfWGS positivity by technology panel used as Main Figure 3D.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Faceted BM positivity (training + test cohorts)
# ══════════════════════════════════════════════════════════════════════════
readr::write_csv(
  combo_tbl %>% mutate(Figure = "Fig_4I_BM_positivity_by_tech_facet6"),
  file.path(outdir_source_data, "Fig_4I_BM_positivity_by_tech_facet6_source_data.csv")
)








# Main Figure 4C: repeat the positivity-rate analysis using blood/cfDNA-derived
# mutation lists and the blood cfWGS model calls.
#
# Analyst note: count paired assay availability at each frontline landmark
# timepoint before building the plot.
front_counts <- dat %>%
  filter(Cohort == "Frontline", !is.na(landmark_tp)) %>%
  mutate(
    has_screen   = !is.na(Blood_zscore_only_sites_call),
    has_clinical = !is.na(Flow_Binary) | !is.na(Adaptive_Binary)
  ) %>%
  group_by(landmark_tp) %>%
  summarise(
    n_total_cfWGS_screen = sum(has_screen, na.rm = TRUE),
    n_with_clinical      = sum(has_screen & has_clinical, na.rm = TRUE),
    frac_string          = paste0(n_with_clinical, "/", n_total_cfWGS_screen),
    .groups = "drop"
  )

concord_sentence <- front_counts %>%
  mutate(tp_label = recode(landmark_tp, !!!tp_labels)) %>%
  arrange(match(tp_label, c("post-ASCT", "1-year maintenance"))) %>%
  summarise(
    sentence = glue(
      "We assessed concordance between cfDNA-based MRD and clinical assays at {glue_collapse(glue('{tp_label} in samples with at least one clinical MRD result (n = {n_with_clinical}/{n_total_cfWGS_screen})'), sep = ' and ')}."
    )
  ) %>%
  pull(sentence)

concord_sentence

front_tbl <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    !is.na(Blood_zscore_only_sites_call),
    !is.na(Flow_Binary) | !is.na(Adaptive_Binary)
  ) %>%
  # Add the blood screening-threshold call used in Figure 4C.
  mutate(
    Blood_zscore_screen_call  = as.integer(Blood_zscore_only_sites_prob >= 0.380)
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call, Blood_zscore_screen_call, 
                  EasyM_clinician_binary, EasyM_optimized_binary),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  dplyr::summarise(
    n_total  = dplyr::n(),
    n_pos    = sum(Result == 1L, na.rm = TRUE),
    pos_rate = n_pos / n_total,
    .groups  = "drop"
  ) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary         = "MFC",
      Adaptive_Binary     = "clonoSEQ",
      Blood_zscore_only_sites_call = "cfWGS (confirm)",
      Blood_zscore_screen_call = "cfWGS (screen)",
      EasyM_clinician_binary = "EasyM (any)",
      EasyM_optimized_binary = "EasyM (opt)"
    )
  )
#      Blood_zscore_only_sites_call = "cfWGS (confirm)",
#      Blood_zscore_screen_call = "cfWGS (screen)"
#    )
#  )

# Analyst note: quantify patients negative at both landmark timepoints by
# clinical assays and by blood/cfDNA-informed cfWGS. These values support
# manuscript text but are not separately staged as a figure/table artifact.
front_long <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    landmark_tp %in% c("Post-ASCT", "Maintenance-1yr")
  ) %>%
  filter(!is.na(Blood_zscore_only_sites_call)) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call, 
                  EasyM_clinician_binary, EasyM_optimized_binary),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary                     = "MFC",
      Adaptive_Binary                 = "clonoSEQ",
      Blood_zscore_only_sites_call = "cfWGS",
      EasyM_clinician_binary          = "EasyM (any)",
      EasyM_optimized_binary          = "EasyM (opt)"
    )
  )

# Find patients negative at both timepoints, per method
neg_both_tbl <- front_long %>%
  group_by(Technology, Patient) %>%
  summarise(
    n_tp = n_distinct(landmark_tp),
    all_neg = (n_tp == 2 & all(Result == 0L)),
    .groups = "drop"
  ) %>%
  filter(all_neg)

# Pull sets of patients
neg_sets <- split(neg_both_tbl$Patient, neg_both_tbl$Technology)

mfc_neg   <- neg_sets[["MFC"]]
seq_neg   <- neg_sets[["clonoSEQ"]]
cfwgs_neg <- neg_sets[["cfWGS"]]

# Overlaps
mfc_also_cfwgs <- intersect(mfc_neg, cfwgs_neg)
seq_also_cfwgs <- intersect(seq_neg, cfwgs_neg)

# Build sentences
mfc_sentence <- sprintf(
  "Of %d patients that were negative at both timepoints by MFC, %d (%.1f%%) were also negative by cfWGS.",
  length(mfc_neg),
  length(mfc_also_cfwgs),
  100 * length(mfc_also_cfwgs) / max(1, length(mfc_neg))
)

seq_sentence <- sprintf(
  "Of %d patients that were negative at both timepoints by clonoSEQ, %d (%.1f%%) were also negative by cfWGS.",
  length(seq_neg),
  length(seq_also_cfwgs),
  100 * length(seq_also_cfwgs) / max(1, length(seq_neg))
)

mfc_sentence
seq_sentence


# ---------------------------------------------------------------------------
#  4.  NON‑FRONTLINE cohort: pooled positivity -------------------------------
non_tbl <- dat %>%
  mutate(landmark_tp = "All timepoints") %>%
  filter(!timepoint_info %in% c("Baseline", "Diagnosis")) %>% 
  filter(
    Cohort == "Non-frontline",
    !is.na(Blood_zscore_only_detection_rate_call),
    !is.na(MRD_truth) # restrict to only ones with MRD for fair comparison
  ) %>%
  # Add the blood screening-threshold call used in Figure 4C.
  mutate(
    Blood_zscore_screen_call  = as.integer(Blood_zscore_only_sites_prob >= 0.380)
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call, Blood_zscore_screen_call, 
                  EasyM_clinician_binary, EasyM_optimized_binary),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  dplyr::summarise(
    n_total  = dplyr::n(),
    n_pos    = sum(Result == 1L, na.rm = TRUE),
    pos_rate = n_pos / n_total,
    .groups  = "drop"
  ) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary         = "MFC",
      Adaptive_Binary     = "clonoSEQ",
      Blood_zscore_only_sites_call = "cfWGS (confirm)",
      Blood_zscore_screen_call = "cfWGS (screen)",
      EasyM_clinician_binary = "EasyM (any)",
      EasyM_optimized_binary = "EasyM (opt)"
    )
  )


# Clean and re-factor the frontline timepoint column for plotting.
front_tbl <- front_tbl %>%
  # Normalize unicode hyphen to ASCII hyphen-minus for consistent rendering.
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) %>%
  # Convert to a factor with Post-ASCT first.
  mutate(
    landmark_tp = factor(landmark_tp,
                         levels = c("Post-ASCT", "Maintenance-1yr"))
  )

# Clean and re-factor the non-frontline/test-cohort timepoint column.
non_tbl <- non_tbl %>%
  # Normalize unicode hyphen to ASCII hyphen-minus for consistent rendering.
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) 



# Combine training/frontline and test/non-frontline tables into one faceted
# blood/cfDNA-informed positivity panel.
front_tbl_blood <- front_tbl  %>% mutate(Cohort = "Training Cohort")
non_tbl_blood   <- non_tbl    %>% mutate(Cohort = "Test Cohort")
combo_tbl  <- bind_rows(front_tbl_blood, non_tbl_blood)

combo_tbl <- combo_tbl %>%
  mutate(
    Cohort      = factor(Cohort,    levels = c("Training Cohort","Test Cohort")),
    landmark_tp = factor(landmark_tp,
                         levels = c("Post-ASCT",
                                    "Maintenance-1yr",
                                    "All timepoints")),
    Technology  = factor(Technology,
                         levels = c(
                           "cfWGS (screen)",
                           "cfWGS (confirm)",
                           "clonoSEQ",
                           "MFC",
                           "EasyM (opt)",
                           "EasyM (any)"
                         ))
    
  )

# Analyst note: build prose snippets for checking the Figure 4C percentages
# against the plotted denominators.
fig_ref   <- "Figure 4C"
cohort_in <- "Training Cohort"                    # change if you want the Test Cohort

# Helper: fetch row for a (timepoint, technology, cohort) and format "XX% (a/b)"
pull_fmt <- function(df, tp, tech, cohort = cohort_in, digits = 0) {
  row <- df %>%
    filter(landmark_tp == tp, Technology == tech, Cohort == cohort) %>%
    slice(1)
  
  if (nrow(row) == 0 || is.na(row$n_total) || row$n_total == 0) return("NA")
  
  pct <- round(100 * row$pos_rate, digits)
  sprintf("%d%% (%d/%d)", pct, row$n_pos, row$n_total)
}

# Helper: build one sentence for a given timepoint
build_sentence <- function(df, tp, fig = fig_ref) {
  scr  <- pull_fmt(df, tp, "cfWGS (screen)")
  conf <- pull_fmt(df, tp, "cfWGS (confirm)")
  seqv <- pull_fmt(df, tp, "clonoSEQ")
  mfc  <- pull_fmt(df, tp, "MFC")
  
  sprintf(
    "At %s, the cfWGS screening threshold identified %s as MRD-positive, whereas the confirmatory threshold was positive in %s, with clinical assays showing %s by clonoSEQ and %s by MFC (%s).",
    tp, scr, conf, seqv, mfc, fig
  )
}

tp1 <- "Post-ASCT"
tp2 <- "Maintenance-1yr"

sentence_postASCT   <- build_sentence(combo_tbl, tp1)
sentence_maint1year <- build_sentence(combo_tbl, tp2)

sentence_postASCT
sentence_maint1year


p_pos_by_tech <- ggplot(combo_tbl, 
                        aes(x    = Technology,
                            y    = pos_rate * 100,
                            fill = landmark_tp)) +
  # ① solid bars with black outline
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  
  # ② labels on top of bars
  geom_text(aes(label = sprintf("%d%%", round(pos_rate * 100))),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 2.5,
            family   = "sans") +
  
  # Manual fill to match the manuscript color palette.
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  
  # ④ free‐x per facet so “clonoSEQ” drops out of Test, and percent y‐axis
  facet_wrap(~ Cohort, nrow = 1, scales = "free_x") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, 100)
  ) +
  
  # ⑤ titles
  labs(
    title = "MRD Positivity by Technology (cfDNA-derived mutation lists)",
    x     = "Technology",
    y     = "Positivity Rate"
  ) +
  
  # ⑥ clean up the theme to match p_perf
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y     = element_text(size = 9),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "top"
  )

# Save the final Figure 4C working PNG before staging it to the manuscript
# output directory.
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_5I_Blood_positivity_by_tech_facet_updated6.png"),
  plot     = p_pos_by_tech,
  width    = 7.5,    # wider to accommodate two facets
  height   = 4.5,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 4C
#
# What this is:
#   Blood/cfDNA-informed cfWGS positivity by clinical-assay technology and
#   cohort.
#
# Why it is here:
#   This PNG is the final Main Figure 4C component.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = file.path(OUTPUT_DIR_FIGURES, "Fig_5I_Blood_positivity_by_tech_facet_updated6.png"),
  artifact_id = "FIG4C",
  role = "figure_panel_png",
  description = "Blood/cfDNA-informed cfWGS positivity by technology panel used as Main Figure 4C.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)






# Concordance and PPV/NPV summaries for BM-informed calls. These tables support
# manuscript text and downstream supplementary-table construction.
# ---------------------------------------------------------------------------
# 1.  Helper to compute pairwise concordance ---------------------------------
pair_concord <- function(df, test_a, test_b) {
  tmp <- df %>%
    select(Sample_Code, !!test_a, !!test_b) %>%
    drop_na() %>%
    mutate(Concordant = (.data[[test_a]] == .data[[test_b]]))
  
  tibble(
    test_a      = test_a,
    test_b      = test_b,
    n_total     = nrow(tmp),
    n_conc      = sum(tmp$Concordant),
    conc_rate   = n_conc / n_total,
    # discordant categories
    a_pos_b_neg = sum(tmp[[test_a]] == 1 & tmp[[test_b]] == 0),
    a_neg_b_pos = sum(tmp[[test_a]] == 0 & tmp[[test_b]] == 1)
  )
}

# ---------------------------------------------------------------------------
# 2.  Landmark labels --------------------------------------------------------
dat <- dat %>%
  mutate(landmark = recode(timepoint_info,
                           "Post_transplant"  = "Post_ASCT",
                           "1yr maintenance"  = "Maintenance",
                           .default           = NA_character_))

# ---------------------------------------------------------------------------
# 3.  FRONTLINE: concordance & positivity  -----------------------------------
front <- dat %>% filter(Cohort == "Frontline", !is.na(landmark)) %>% filter(!is.na(BM_base_zscore_call))

# --- 3a. Pairwise concordance at Post‑ASCT ----------------------------------
pa   <- front %>% filter(landmark == "Post_ASCT")
post_conc <- bind_rows(
  pair_concord(pa, "BM_zscore_only_detection_rate_call", "Adaptive_Binary"),
  pair_concord(pa, "BM_zscore_only_detection_rate_call", "Flow_Binary"),
  pair_concord(pa, "Adaptive_Binary",  "Flow_Binary")
)

# --- 3b. Pairwise concordance at Maintenance --------------------------------
ma   <- front %>% filter(landmark == "Maintenance")
maint_conc <- bind_rows(
  pair_concord(ma, "BM_zscore_only_detection_rate_call", "Adaptive_Binary"),
  pair_concord(ma, "BM_zscore_only_detection_rate_call", "Flow_Binary")
)

# --- 3c. Positivity counts ---------------------------------------------------
pos_tbl <- front %>%
  filter(landmark %in% c("Post_ASCT", "Maintenance")) %>%
  pivot_longer(
    cols      = c(BM_zscore_only_detection_rate_call, Adaptive_Binary, Flow_Binary),
    names_to  = "Test",
    values_to = "Result"
  ) %>%
  drop_na(Result) %>%
  group_by(landmark, Test) %>%
  dplyr::summarise(
    pos = sum(Result == 1L, na.rm = TRUE),
    tot = dplyr::n(),
    .groups = "drop"
  )

# Non-frontline/test-cohort BM concordance summary, aggregated across eligible
# follow-up samples.
non <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(BM_zscore_only_detection_rate_call)) %>% filter(timepoint_info != "Baseline") %>% 
  filter(timepoint_info != "Diagnosis")

non_conc <- pair_concord(non, "BM_zscore_only_detection_rate_call", "Flow_Binary")




# ---------------------------------------------------------------------------
# 4.  FRONTLINE PPV / NPV for BM_zscore_only_detection_rate_call --------------------------------
ppv_npv <- function(df, pred_col = "BM_zscore_only_detection_rate_call") {
  # Build a 2×2 table of prediction vs truth
  tbl <- table(
    Pred  = df[[pred_col]],
    Truth = df$MRD_truth,
    useNA = "no"
  )
  
  # Ensure the table has both rows 0/1 and columns 0/1
  all_lv <- c("0", "1")
  tbl     <- tbl[all_lv, all_lv, drop = FALSE]  # missing rows/cols become NA
  tbl[is.na(tbl)] <- 0                           # convert those NA counts back to 0
  
  # Extract TP, FP, TN, FN as scalars
  TP <- tbl["1", "1"]
  FP <- tbl["1", "0"]
  TN <- tbl["0", "0"]
  FN <- tbl["0", "1"]
  
  tibble(
    PPV = TP / (TP + FP),
    NPV = TN / (TN + FN)
  )
}

# PPV/NPV versus the combined clinical MRD truth for the two frontline
# landmark timepoints.
ppv_post  <- ppv_npv(pa %>% filter(!is.na(MRD_truth)))
ppv_maint <- ppv_npv(ma %>% filter(!is.na(MRD_truth)))

print(ppv_post)
print(ppv_maint)


# ---------------------------------------------------------------------------
# 5.  NON-FRONTLINE cohort ----------------------------------------------------
non <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(BM_zscore_only_detection_rate_call)) %>% filter(timepoint_info != "Diagnosis")  %>% filter(timepoint_info != "Baseline")

non_cm <- non %>%
  filter(Cohort != "Frontline",
         !is.na(BM_zscore_only_detection_rate_call),
         !is.na(MRD_truth)) %>%
  tabyl(BM_zscore_only_detection_rate_call, MRD_truth) %>%
  # ensure integer rows 0 and 1 exist
  complete(
    BM_zscore_only_detection_rate_call = c(0L, 1L), 
    fill = list(`0` = 0, `1` = 0)
  ) %>%
  column_to_rownames("BM_zscore_only_detection_rate_call")


# extract counts
TPn <- non_cm["1", "1"]
FPn <- non_cm["1", "0"]
TNn <- non_cm["0", "0"]
FNn <- non_cm["0", "1"]

sens_non <- TPn / (TPn + FNn)
spec_non <- TNn / (TNn + FPn)

sens_non
spec_non

# Overall positivity non‑frontline
non_pos <- non %>%
  pivot_longer(
    cols      = c(BM_zscore_only_detection_rate_call, Flow_Binary),
    names_to  = "Test",
    values_to = "Result"
  ) %>%
  drop_na(Result) %>%
  group_by(Test) %>%
  dplyr::summarise(
    pos = sum(Result == 1L, na.rm = TRUE),
    tot = dplyr::n(),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------
# 6.  Collect outputs --------------------------------------------------------
stats_out <- list(
  post_concordance   = post_conc,
  maint_concordance  = maint_conc,
  frontline_pos      = pos_tbl,
  PPV_post_ASCT      = ppv_post,
  PPV_maintenance    = ppv_maint,
  nonfront_sens      = sens_non,
  nonfront_spec      = spec_non,
  nonfront_pos       = non_pos
)

print(stats_out)

# ---------------------------------------------------------------------------
# 7.  Manuscript-text denominator check --------------------------------------
# Print a paragraph-style summary so analysts can cross-check values in the
# manuscript text. This does not create a final manuscript artifact.
write_para <- TRUE

if (write_para) {
  # helpers to pull numbers
  g <- function(a,b,df) df %>% filter(test_a==a, test_b==b) %>% pull(n_conc)
  n <- function(a,b,df) df %>% filter(test_a==a, test_b==b) %>% pull(n_total)
  r <- function(a,b,df) df %>% filter(test_a==a, test_b==b) %>% pull(conc_rate)
  
  # post‑ASCT numbers
  X  <- g("BM_zscore_only_detection_rate_call","Adaptive_Binary", post_conc)
  Y  <- n("BM_zscore_only_detection_rate_call","Adaptive_Binary", post_conc)
  XX <- sprintf("%.0f", 100*r("BM_zscore_only_detection_rate_call","Adaptive_Binary", post_conc))
  Xp <- g("BM_zscore_only_detection_rate_call","Flow_Binary", post_conc)
  Yp <- n("BM_zscore_only_detection_rate_call","Flow_Binary", post_conc)
  XXp<- sprintf("%.0f", 100*r("BM_zscore_only_detection_rate_call","Flow_Binary", post_conc))
  Z  <- g("Adaptive_Binary","Flow_Binary", post_conc)
  W  <- n("Adaptive_Binary","Flow_Binary", post_conc)
  YY <- sprintf("%.0f",100*r("Adaptive_Binary","Flow_Binary", post_conc))
  
  # discordant counts
  n_cf_pos_cl_neg <- post_conc %>%
    filter(test_a=="BM_zscore_only_detection_rate_call", test_b=="Adaptive_Binary") %>%
    pull(a_pos_b_neg)
  m_cf_neg_cl_pos <- post_conc %>%
    filter(test_a=="BM_zscore_only_detection_rate_call", test_b=="Adaptive_Binary") %>%
    pull(a_neg_b_pos)
  
  # maintenance
  A  <- g("BM_zscore_only_detection_rate_call","Adaptive_Binary", maint_conc)
  B  <- n("BM_zscore_only_detection_rate_call","Adaptive_Binary", maint_conc)
  AA <- sprintf("%.0f",100*r("BM_zscore_only_detection_rate_call","Adaptive_Binary", maint_conc))
  C  <- g("BM_zscore_only_detection_rate_call","Flow_Binary", maint_conc)
  D  <- n("BM_zscore_only_detection_rate_call","Flow_Binary", maint_conc)
  BB <- sprintf("%.0f",100*r("BM_zscore_only_detection_rate_call","Flow_Binary", maint_conc))
  
  p <- ppv_post$PPV; q <- ppv_post$NPV
  p2<- ppv_maint$PPV; q2<- ppv_maint$NPV
  
  para <- glue("
    At post-ASCT, cfWGS agreed with clonoSEQ in {X}/{Y} ({XX}%) samples and with MFC in {Xp}/{Yp} ({XXp}%). 
    clonoSEQ vs. MFC were concordant in {Z}/{W} ({YY}%) paired samples. 
    Of the discordant post-ASCT samples, cfWGS was positive/ clonoSEQ negative in {n_cf_pos_cl_neg} cases and negative/ clonoSEQ positive in {m_cf_neg_cl_pos}. 
    At the 1-year maintenance timepoint, cfWGS agreed with clonoSEQ in {A}/{B} ({AA}%) samples and with MFC in {C}/{D} ({BB}%). 
    The PPV and NPV of cfWGS were {sprintf('%.0f',p*100)}% and {sprintf('%.0f',q*100)}% at post-ASCT, and {sprintf('%.0f',p2*100)}% and {sprintf('%.0f',q2*100)}% at maintenance. 
    In the non-frontline cohort, sensitivity and specificity of cfWGS were {sprintf('%.0f',stats_out$nonfront_sens*100)}% and {sprintf('%.0f',stats_out$nonfront_spec*100)}%, with an overall positivity rate of {stats_out$nonfront_pos %>% filter(Test=='BM_zscore_only_detection_rate_call') %>% summarise(sprintf('%.0f%%', 100*pos/tot)) %>% pull()}.
  ")
  
  cat(para)
}


# Assay-specific PPV/NPV summaries.
# These compare cfWGS to each clinical assay separately, rather than to the
# combined MRD truth variable.
# 1.  General PPV/NPV helper that takes any truth column  -------------------
ppv_npv_any <- function(df, pred_col = "BM_zscore_only_detection_rate_call", truth_col) {
  tbl <- table(
    Pred  = df[[pred_col]],
    Truth = df[[truth_col]],
    useNA = "no"
  )
  # ensure both levels exist
  all_lv <- c("0", "1")
  tbl     <- tbl[all_lv, all_lv, drop = FALSE]
  tbl[is.na(tbl)] <- 0
  
  TP <- tbl["1","1"]
  FP <- tbl["1","0"]
  TN <- tbl["0","0"]
  FN <- tbl["0","1"]
  
  tibble(
    truth         = truth_col,
    PPV           = TP / (TP + FP),
    NPV           = TN / (TN + FN),
    n_pos_calls   = TP + FP,
    n_neg_calls   = TN + FN
  )
}

# 2.  Filter post-ASCT frontline samples -------------------------------
pa <- dat %>%
  filter(Cohort == "Frontline", landmark == "Post_ASCT")

# 3.  Compute PPV/NPV versus clonoSEQ  ----------------------------------
ppv_clono <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Adaptive_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Adaptive_Binary"
)

# 4.  Compute PPV/NPV versus MFC  ---------------------------------------
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Flow_Binary"
)

# 5.  Bind together and print -------------------------------------------
bind_rows(ppv_clono, ppv_mfc)


# Repeat assay-specific PPV/NPV at the frontline maintenance landmark.
pa <- dat %>%
  filter(Cohort == "Frontline", landmark == "Maintenance")
ppv_clono <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Adaptive_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Adaptive_Binary"
)
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_clono, ppv_mfc)

# Repeat assay-specific PPV/NPV for the non-frontline/test cohort.
pa <- dat %>%
  filter(Cohort == "Non-frontline") %>% filter(timepoint_info != "Diagnosis")

ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_mfc)


# Manuscript-text denominator check for assay-specific PPV/NPV summaries.
# ---- compact wrappers to fetch %PPV/%NPV for a timepoint vs a given truth assay ----
get_ppvnpv <- function(df, landmark_value, truth_col) {
  tmp <- df %>%
    filter(Cohort == "Frontline", landmark == landmark_value, !is.na(.data[[truth_col]]))
  if (nrow(tmp) == 0) return(c(PPV = NA_real_, NPV = NA_real_))
  out <- ppv_npv_any(tmp, pred_col = "BM_zscore_only_detection_rate_call", truth_col = truth_col)
  c(PPV = out$PPV, NPV = out$NPV)
}

fmt_pct0 <- function(x) ifelse(is.na(x), "NA", sprintf("%.0f%%", 100*x))

# ---- pull numbers ----
# Post-ASCT
ppvnpv_post_seq <- get_ppvnpv(dat, "Post_ASCT", "Adaptive_Binary")  # clonoSEQ
ppvnpv_post_mfc <- get_ppvnpv(dat, "Post_ASCT", "Flow_Binary")      # MFC

# Maintenance-1yr
ppvnpv_maint_seq <- get_ppvnpv(dat, "Maintenance", "Adaptive_Binary")
ppvnpv_maint_mfc <- get_ppvnpv(dat, "Maintenance", "Flow_Binary")

para2 <- glue("
To further assess assay performance, we examined concordance between cfWGS and clinical MRD tests at each time point. 
At post-ASCT, cfWGS demonstrated high concordance with both clinical assays, agreeing with clonoSEQ in {X}/{Y} ({XX}%) samples (PPV = {fmt_pct0(ppvnpv_post_seq['PPV'])}, NPV = {fmt_pct0(ppvnpv_post_seq['NPV'])}) 
and with MFC in {Xp}/{Yp} ({XXp}%) samples (PPV = {fmt_pct0(ppvnpv_post_mfc['PPV'])}, NPV = {fmt_pct0(ppvnpv_post_mfc['NPV'])}). 
At the 1-year maintenance timepoint, cfWGS agreed with clonoSEQ in {A}/{B} ({AA}%) samples (PPV = {fmt_pct0(ppvnpv_maint_seq['PPV'])}, NPV = {fmt_pct0(ppvnpv_maint_seq['NPV'])}) 
and with MFC in {C}/{D} ({BB}%) samples (PPV = {fmt_pct0(ppvnpv_maint_mfc['PPV'])}, NPV = {fmt_pct0(ppvnpv_maint_mfc['NPV'])}).
")

cat(para2)


# Export BM-informed concordance helper tables used for review and source-data
# tracing. Final manuscript-staged workbooks are generated later in this script.
# 1. Export post-ASCT pairwise concordance (frontline BM)
readr::write_csv(
  post_conc,
  file.path(outdir, "Frontline_BoneMarrow_PostASCT_Pairwise_Concordance3.csv")
)

# 2. Export maintenance-timepoint pairwise concordance (frontline BM)
readr::write_csv(
  maint_conc,
  file.path(outdir, "Frontline_BoneMarrow_Maintenance_Pairwise_Concordance3.csv")
)

# 3. Export frontline positivity counts by test & landmark (Post_ASCT + Maintenance)
readr::write_csv(
  pos_tbl,
  file.path(outdir, "Frontline_BoneMarrow_Positivity_PostASCT_and_Maintenance3.csv")
)

# 4. Export PPV/NPV at Post-ASCT for BM_zscore_only_detection_rate_call
readr::write_csv(
  ppv_post,
  file.path(outdir, "Frontline_BoneMarrow_PostASCT_PPV_NPV3.csv")
)

# 5. Export PPV/NPV at Maintenance for BM_zscore_only_detection_rate_call
readr::write_csv(
  ppv_maint,
  file.path(outdir, "Frontline_BoneMarrow_Maintenance_PPV_NPV3.csv")
)




# Extended Data Figure 5E-G: BM-informed confusion-matrix inputs against the
# combined clinical MRD truth.
# helper: build a tidy 2 × 2 contingency table ----------------------------
make_ct <- function(df,
                    pred  = "BM_zscore_only_detection_rate_call",   # cfWGS
                    truth = "MRD_truth") {                          # reference
  out <- df %>%
    select(all_of(c(pred, truth))) %>%
    drop_na() %>%                                   # keep only paired calls
    mutate(across(everything(), as.integer)) %>%    # ensure 0/1 integers
    count(!!sym(pred), !!sym(truth), name = "n") %>%
    complete(!!sym(pred) := 0:1,
             !!sym(truth):= 0:1,
             fill = list(n = 0)) %>%                # add missing cells
    mutate(row = ifelse(!!sym(pred)==1, "Pred_Pos", "Pred_Neg"),
           col = ifelse(!!sym(truth)==1,"Truth_Pos","Truth_Neg")) %>%
    select(row, col, n)
  
  # add summary columns
  tp <- out %>% filter(row=="Pred_Pos",  col=="Truth_Pos")  %>% pull(n)
  fp <- out %>% filter(row=="Pred_Pos",  col=="Truth_Neg")  %>% pull(n)
  fn <- out %>% filter(row=="Pred_Neg",  col=="Truth_Pos")  %>% pull(n)
  tn <- out %>% filter(row=="Pred_Neg",  col=="Truth_Neg")  %>% pull(n)
  
  tibble(
    TP = tp, FP = fp, FN = fn, TN = tn,
    Sensitivity = TP/(TP+FN),
    Specificity = TN/(TN+FP),
    PPV = TP/(TP+FP),
    NPV = TN/(TN+FN)
  )
}

# 1.  FRONTLINE – post-ASCT -----------------------------------------------
ct_post_ASCT <- front %>%
  filter(landmark == "Post_ASCT") %>%
  make_ct()

# 2.  FRONTLINE – 1-year maintenance --------------------------------------
ct_maint <- front %>%
  filter(landmark == "Maintenance") %>%
  make_ct()

# 3.  NON-FRONTLINE – all baseline / follow-up samples --------------------
ct_nonfront <- non %>% make_ct()

# Export compact contingency-table workbook for traceability.
writexl::write_xlsx(
  list(
    "Contingency_Post_ASCT"     = ct_post_ASCT,
    "Contingency_Maintenance"   = ct_maint,
    "Contingency_NonFrontline"  = ct_nonfront
  ),
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables3.xlsx")
)

# Build comparator-specific BM contingency tables against MFC and clonoSEQ.
ct_post_Flow   <- front  %>% filter(landmark=="Post_ASCT") %>% 
  make_ct(truth = "Flow_Binary")
ct_post_Clono  <- front  %>% filter(landmark=="Post_ASCT") %>% 
  make_ct(truth = "Adaptive_Binary")

ct_maint_Flow  <- front  %>% filter(landmark=="Maintenance") %>% 
  make_ct(truth = "Flow_Binary")
ct_maint_Clono <- front  %>% filter(landmark=="Maintenance") %>% 
  make_ct(truth = "Adaptive_Binary")

ct_non_Flow    <- non %>% make_ct(truth = "Flow_Binary")

writexl::write_xlsx(
  list(
    "Post_ASCT_vs_Flow"        = ct_post_Flow,
    "Post_ASCT_vs_clonoSEQ"    = ct_post_Clono,
    "Maintenance_vs_Flow"      = ct_maint_Flow,
    "Maintenance_vs_clonoSEQ"  = ct_maint_Clono,
    "NonFront_vs_Flow"         = ct_non_Flow
  ),
  path = file.path(outdir, "cfWGS_contingency_vs_Flow_clonoSEQ3.xlsx")
)


# Plotting helper for Extended Data Figure 5E-G: turn a one-row TP/FP/FN/TN
# table into the long format expected by ggplot.
ct_to_long <- function(ct_row, label){
  with(ct_row, tibble(
    Obs  = rep(c("neg","pos"), each = 2),          # rows
    Pred = rep(c("neg","pos"), times = 2),         # cols
    Count = c(TN, FP, FN, TP),
    model = label,
    PPV   = PPV,
    NPV   = NPV
  ))
}

## ───────────────────────────────────────────────────────────────
## B.  Build the three data frames to plot
## ───────────────────────────────────────────────────────────────
# 1) Post-ASCT (Frontline)
cm_post <- bind_rows(
  ct_to_long(ct_post_Flow ,  "Flow (MFC)"),
  ct_to_long(ct_post_Clono,  "clonoSEQ")
) |> mutate(model = factor(model, levels = c("Flow (MFC)","clonoSEQ")))

# 2) Maintenance (Frontline)
cm_maint <- bind_rows(
  ct_to_long(ct_maint_Flow , "Flow (MFC)"),
  ct_to_long(ct_maint_Clono, "clonoSEQ")
) |> mutate(model = factor(model, levels = c("Flow (MFC)","clonoSEQ")))

# 3) Non-frontline (one comparator only)
cm_non <- ct_to_long(ct_non_Flow, "Flow (MFC)")

## ───────────────────────────────────────────────────────────────
## C.  Confusion-matrix plotting helper
## ───────────────────────────────────────────────────────────────
# Active plotting helper used for Extended Data Figure 5E-G and Extended Data
# Figure 7F-H. Earlier draft versions of this helper were removed from the
# active script to keep the manuscript path unambiguous.
plot_cm <- function(df, main_title,
                    col_low = "#f2f2f2", col_high = "#4a4a4a") {
  
  df <- df %>%
    mutate(
      # keep corners consistent
      Obs  = factor(Obs,  levels = c("neg","pos")),
      Pred = factor(Pred, levels = c("pos","neg"))
    ) %>%
    group_by(model) %>%                             # choose white text on dark cells (per facet)
    mutate(text_col = if_else(Count >= 0.7 * max(Count), "white", "black")) %>%
    ungroup()
  
  ggplot(df, aes(x = Obs, y = Pred, fill = Count)) +
    geom_tile(colour = "white", linewidth = 0) +
    geom_text(aes(label = Count, color = text_col), size = 4, show.legend = FALSE) +
    facet_wrap(~ model, nrow = 1) +
    scale_fill_gradient(
      low    = col_low,
      high   = col_high,
      limits = c(0, max(df$Count)),
      oob    = scales::squish,
      guide  = "none"
    ) +
    scale_color_identity() +
    scale_x_discrete(position = "top") +
    labs(
      x = "Clinical MRD",
      y = "cfWGS MRD",
      title = main_title
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text      = element_text(face = "bold"),
      axis.text.y     = element_text(size = 9),
      axis.text.x     = element_text(size = 9, vjust = 0),
      axis.title      = element_text(size = 10),
      panel.grid      = element_blank(),
      legend.position = "none",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
}



## ───────────────────────────────────────────────────────────────
## D.  Draw & save the three panels
## ───────────────────────────────────────────────────────────────
p_post   <- plot_cm(cm_post ,  "Confusion Matrix at Post-ASCT (Training Cohort)")
p_maint  <- plot_cm(cm_maint,  "Confusion Matrix at 1‑Year Maintenance (Training Cohort)")
p_non    <- plot_cm(cm_non ,   "Confusion Matrix of Test Cohort")

ggsave("Final Tables and Figures/Fig4_confmat_post_ASCT_updated5.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_maintenance5.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_nonfront5.png",
       p_non,   width = 3, height = 2.75, dpi = 600)   # single facet – narrower

# -------------------------------------------------------------------------
# Manuscript outputs: Extended Data Figure 5E-G
#
# What these are:
#   BM-informed cfWGS confusion-matrix panels against clinical comparators at
#   post-ASCT, 1-year maintenance, and non-frontline/test-cohort settings.
#
# Why they are here:
#   The retained final components use older updated4 filenames, while this
#   script currently writes updated5 versions from the same plotted section.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig4_confmat_post_ASCT_updated5.png",
  artifact_id = "EDFIG5E",
  role = "figure_panel_png",
  description = "BM-informed post-ASCT clinical-comparator confusion matrix used as Extended Data Figure 5E.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig4_confmat_maintenance5.png",
  artifact_id = "EDFIG5F",
  role = "figure_panel_png",
  description = "BM-informed one-year maintenance clinical-comparator confusion matrix used as Extended Data Figure 5F.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig4_confmat_nonfront5.png",
  artifact_id = "EDFIG5G",
  role = "figure_panel_png",
  description = "BM-informed non-frontline clinical-comparator confusion matrix used as Extended Data Figure 5G.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Confusion matrices (BM) - individual panels
# ══════════════════════════════════════════════════════════════════════════
readr::write_csv(
  cm_post %>% mutate(Figure = "Fig4_confmat_post_ASCT_updated5"),
  file.path(outdir_source_data, "Fig4_confmat_post_ASCT_updated5_source_data.csv")
)
readr::write_csv(
  cm_maint %>% mutate(Figure = "Fig4_confmat_maintenance5"),
  file.path(outdir_source_data, "Fig4_confmat_maintenance5_source_data.csv")
)
readr::write_csv(
  cm_non %>% mutate(Figure = "Fig4_confmat_nonfront5"),
  file.path(outdir_source_data, "Fig4_confmat_nonfront5_source_data.csv")
)

# Auxiliary combined BM layout retained for provenance. This is not staged as a
# final manuscript artifact; the final ED5E-G components are copied separately
# by the `ms_copy_artifact()` calls above.
combined_cm <- (p_post  +
                  theme(
                    panel.spacing = unit(1, "lines"),
                    plot.margin   = margin(5,5,5,5),
                    # allow annotations to overflow the panel
                    panel.clip    = "off"
                  )
) /
  (p_maint +
     theme(
       panel.spacing = unit(1, "lines"),
       plot.margin   = margin(5,5,5,5),
       panel.clip    = "off"
     )
  ) /
  (p_non   +
     theme(
       panel.spacing = unit(1, "lines"),
       plot.margin   = margin(5,5,5,5),
       panel.clip    = "off"
     )
  ) + 
  plot_layout(ncol = 1, heights = c(1,1,1)) & 
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# optionally give it an overall title
#combined_cm <- combined_cm + plot_annotation(
#  title = "Confusion tables at Post-ASCT, Maintenance & Non-frontline"
#)


# Save auxiliary stacked BM layout.
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three_4.png",
       combined_cm,
       width  = 4,
       height = 7,      # three panels tall
       dpi    = 600)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Combined confusion matrices (stacked)
# ══════════════════════════════════════════════════════════════════════════
# (Source data already exported individually above)

# Auxiliary side-by-side BM layout retained for provenance. This is not staged
# as a final manuscript artifact.
combined_cm <- (p_post +
                  theme(
                    panel.spacing = unit(1, "lines"),
                    plot.margin   = margin(5, 5, 5, 5),
                    panel.clip    = "off"
                  )
) |  # <-- use | instead of / to put them side by side
  (p_maint +
     theme(
       panel.spacing = unit(1, "lines"),
       plot.margin   = margin(5, 5, 5, 5),
       panel.clip    = "off"
     )
  ) |
  (p_non +
     theme(
       panel.spacing = unit(1, "lines"),
       plot.margin   = margin(5, 5, 5, 5),
       panel.clip    = "off"
     )
  ) +
  plot_layout(ncol = 3, widths = c(1, 1, 0.5)) &  # three columns
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# Save auxiliary side-by-side BM layout.
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three_side_by_side4.png",
       combined_cm,
       width  = 15,
       height = 3,      # three panels tall
       dpi    = 600)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Combined confusion matrices (side-by-side)
# ══════════════════════════════════════════════════════════════════════════
# (Source data already exported individually above)






# Extended Data Figure 7F-H: blood/cfDNA-informed confusion-matrix inputs.
front_blood <- dat %>% filter(Cohort == "Frontline", !is.na(landmark)) %>% filter(!is.na(Blood_zscore_only_sites_call))
non_blood <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(Blood_zscore_only_sites_call)) %>% filter(!timepoint_info %in% c("Diagnosis", "Baseline"))


make_ct <- function(df,
                    pred  = "Blood_zscore_only_sites_call",   # cfWGS
                    truth = "MRD_truth") {                          # reference
  out <- df %>%
    select(all_of(c(pred, truth))) %>%
    drop_na() %>%                                   # keep only paired calls
    mutate(across(everything(), as.integer)) %>%    # ensure 0/1 integers
    count(!!sym(pred), !!sym(truth), name = "n") %>%
    complete(!!sym(pred) := 0:1,
             !!sym(truth):= 0:1,
             fill = list(n = 0)) %>%                # add missing cells
    mutate(row = ifelse(!!sym(pred)==1, "Pred_Pos", "Pred_Neg"),
           col = ifelse(!!sym(truth)==1,"Truth_Pos","Truth_Neg")) %>%
    select(row, col, n)
  
  # add summary columns
  tp <- out %>% filter(row=="Pred_Pos",  col=="Truth_Pos")  %>% pull(n)
  fp <- out %>% filter(row=="Pred_Pos",  col=="Truth_Neg")  %>% pull(n)
  fn <- out %>% filter(row=="Pred_Neg",  col=="Truth_Pos")  %>% pull(n)
  tn <- out %>% filter(row=="Pred_Neg",  col=="Truth_Neg")  %>% pull(n)
  
  N      <- tp + fp + fn + tn
  Agree  <- tp + tn
  tibble(
    TP = tp, FP = fp, FN = fn, TN = tn,
    N = N, Agree = Agree,
    Concordance = Agree / N,
    Sensitivity = tp/(tp+fn),
    Specificity = tn/(tn+fp),
    PPV = tp/(tp+fp),
    NPV = tn/(tn+fn)
  )
}

# 1.  FRONTLINE – post-ASCT -----------------------------------------------
ct_post_ASCT <- front_blood %>%
  filter(landmark == "Post_ASCT") %>%
  make_ct()

# 2.  FRONTLINE – 1-year maintenance --------------------------------------
ct_maint <- front_blood %>%
  filter(landmark == "Maintenance") %>%
  make_ct()

# 3.  NON-FRONTLINE – all baseline / follow-up samples --------------------
ct_nonfront <- non_blood %>% make_ct()

# Export compact blood contingency-table workbook for traceability.
writexl::write_xlsx(
  list(
    "Contingency_Post_ASCT"     = ct_post_ASCT,
    "Contingency_Maintenance"   = ct_maint,
    "Contingency_NonFrontline"  = ct_nonfront
  ),
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables_blood3.xlsx")
)

# Build comparator-specific blood/cfDNA contingency tables against MFC and
# clonoSEQ.
ct_post_Flow   <- front_blood  %>% filter(landmark=="Post_ASCT") %>% 
  make_ct(truth = "Flow_Binary")
ct_post_Clono  <- front_blood  %>% filter(landmark=="Post_ASCT") %>% 
  make_ct(truth = "Adaptive_Binary")

ct_maint_Flow  <- front_blood  %>% filter(landmark=="Maintenance") %>% 
  make_ct(truth = "Flow_Binary")
ct_maint_Clono <- front_blood  %>% filter(landmark=="Maintenance") %>% 
  make_ct(truth = "Adaptive_Binary")

ct_non_Flow    <- non_blood %>% make_ct(truth = "Flow_Binary")

writexl::write_xlsx(
  list(
    "Post_ASCT_vs_Flow"        = ct_post_Flow,
    "Post_ASCT_vs_clonoSEQ"    = ct_post_Clono,
    "Maintenance_vs_Flow"      = ct_maint_Flow,
    "Maintenance_vs_clonoSEQ"  = ct_maint_Clono,
    "NonFront_vs_Flow"         = ct_non_Flow
  ),
  path = file.path(outdir, "cfWGS_contingency_vs_Flow_clonoSEQ_blood_calls_3.xlsx") ### these are the metrics
)


# Format an analyst-facing prose check using rounded manuscript-style
# percentages.
fmt_pct <- function(x) sprintf("%.0f%%", 100*x)

post_sentence <- glue(
  "At post-ASCT, confirmatory cfDNA-based MRD demonstrated strong concordance with clonoSEQ ",
  "({ct_post_Clono$Agree}/{ct_post_Clono$N}, {fmt_pct(ct_post_Clono$Concordance)}; ",
  "PPV {fmt_pct(ct_post_Clono$PPV)}, NPV {fmt_pct(ct_post_Clono$NPV)}) ",
  "and moderate concordance with MFC ",
  "({ct_post_Flow$Agree}/{ct_post_Flow$N}, {fmt_pct(ct_post_Flow$Concordance)}; ",
  "PPV {fmt_pct(ct_post_Flow$PPV)}, NPV {fmt_pct(ct_post_Flow$NPV)})."
)

post_sentence

maint_sentence <- glue(
  "At maintenance, confirmatory cfDNA-based MRD demonstrated strong concordance with clonoSEQ ",
  "({ct_maint_Clono$Agree}/{ct_maint_Clono$N}, {fmt_pct(ct_maint_Clono$Concordance)}; ",
  "PPV {fmt_pct(ct_maint_Clono$PPV)}, NPV {fmt_pct(ct_maint_Clono$NPV)}) ",
  "and moderate concordance with MFC ",
  "({ct_maint_Flow$Agree}/{ct_maint_Flow$N}, {fmt_pct(ct_maint_Flow$Concordance)}; ",
  "PPV {fmt_pct(ct_maint_Flow$PPV)}, NPV {fmt_pct(ct_maint_Flow$NPV)})."
)

maint_sentence

# Rebuild long-format confusion-matrix rows for the blood/cfDNA panels.
ct_to_long <- function(ct_row, label){
  with(ct_row, tibble(
    Obs  = rep(c("neg","pos"), each = 2),          # rows
    Pred = rep(c("neg","pos"), times = 2),         # cols
    Count = c(TN, FP, FN, TP),
    model = label,
    PPV   = PPV,
    NPV   = NPV
  ))
}

## ───────────────────────────────────────────────────────────────
## B.  Build the three data frames to plot
## ───────────────────────────────────────────────────────────────
# 1) Post-ASCT (Frontline)
cm_post <- bind_rows(
  ct_to_long(ct_post_Flow ,  "Flow (MFC)"),
  ct_to_long(ct_post_Clono,  "clonoSEQ")
) |> mutate(model = factor(model, levels = c("Flow (MFC)","clonoSEQ")))

# 2) Maintenance (Frontline)
cm_maint <- bind_rows(
  ct_to_long(ct_maint_Flow , "Flow (MFC)"),
  ct_to_long(ct_maint_Clono, "clonoSEQ")
) |> mutate(model = factor(model, levels = c("Flow (MFC)","clonoSEQ")))

# 3) Non-frontline (one comparator only)
cm_non <- ct_to_long(ct_non_Flow, "Flow (MFC)")

## ───────────────────────────────────────────────────────────────
##  Draw & save the three panels
## ───────────────────────────────────────────────────────────────
p_post   <- plot_cm(cm_post ,  "Confusion Matrix at Post-ASCT (Training Cohort)")
p_maint  <- plot_cm(cm_maint,  "Confusion Matrix at 1‑Year Maintenance (Training Cohort)")
p_non    <- plot_cm(cm_non ,   "Confusion Matrix of Test Cohort")

ggsave("Final Tables and Figures/Fig5_confmat_post_ASCT_blood_updated6.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_maintenance_blood_updated6.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_nonfront_blood_updated6.png",
       p_non,   width = 3, height = 2.75, dpi = 600)   # single facet – narrower

# -------------------------------------------------------------------------
# Manuscript outputs: Extended Data Figure 7F-H
#
# What these are:
#   Blood/cfDNA-informed cfWGS confusion-matrix panels against clinical
#   comparators at post-ASCT, 1-year maintenance, and non-frontline/test-cohort
#   settings.
#
# Why they are here:
#   The retained final components use older updated5 filenames, while this
#   script currently writes updated6 versions from the same plotted section.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig5_confmat_post_ASCT_blood_updated6.png",
  artifact_id = "EDFIG7F",
  role = "figure_panel_png",
  description = "Blood/cfDNA-informed post-ASCT clinical-comparator confusion matrix used as Extended Data Figure 7F.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig5_confmat_maintenance_blood_updated6.png",
  artifact_id = "EDFIG7G",
  role = "figure_panel_png",
  description = "Blood/cfDNA-informed one-year maintenance clinical-comparator confusion matrix used as Extended Data Figure 7G.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig5_confmat_nonfront_blood_updated6.png",
  artifact_id = "EDFIG7H",
  role = "figure_panel_png",
  description = "Blood/cfDNA-informed non-frontline clinical-comparator confusion matrix used as Extended Data Figure 7H.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: Blood confusion matrices (cfWGS vs clinical - individual panels)
# ══════════════════════════════════════════════════════════════════════════
readr::write_csv(
  cm_post %>% mutate(Figure = "Fig5_confmat_post_ASCT_blood_updated6"),
  file.path(outdir_source_data, "Fig5_confmat_post_ASCT_blood_updated6_source_data.csv")
)
readr::write_csv(
  cm_maint %>% mutate(Figure = "Fig5_confmat_maintenance_blood_updated6"),
  file.path(outdir_source_data, "Fig5_confmat_maintenance_blood_updated6_source_data.csv")
)
readr::write_csv(
  cm_non %>% mutate(Figure = "Fig5_confmat_nonfront_blood_updated6"),
  file.path(outdir_source_data, "Fig5_confmat_nonfront_blood_updated6_source_data.csv")
)

# Auxiliary combined blood layout retained for provenance. This is not staged
# as a final manuscript artifact; the final ED7F-H components are copied
# separately by the `ms_copy_artifact()` calls above.
combined_cm <- (p_post  +
                  theme(
                    panel.spacing = unit(1, "lines"),
                    plot.margin   = margin(5,5,5,5),
                    # allow annotations to overflow the panel
                    panel.clip    = "off"
                  )
) /
  (p_maint +
     theme(
       panel.spacing = unit(1, "lines"),
       plot.margin   = margin(5,5,5,5),
       panel.clip    = "off"
     )
  ) /
  (p_non   +
     theme(
       panel.spacing = unit(1, "lines"),
       plot.margin   = margin(5,5,5,5),
       panel.clip    = "off"
     )
  ) + 
  plot_layout(ncol = 1, heights = c(1,1,1)) & 
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# optionally give it an overall title
#combined_cm <- combined_cm + plot_annotation(
#  title = "Confusion tables at Post-ASCT, Maintenance & Non-frontline"
#)


# Save auxiliary stacked blood layout.
ggsave("Final Tables and Figures/Fig5J_confusion_matrices_all_three_blood.png",
       combined_cm,
       width  = 4,
       height = 7,      # three panels tall
       dpi    = 600)



# Supplementary analysis: repeat the blood confusion-matrix summaries using the
# blood plus fragmentomics combined model. These outputs are not currently
# staged as final manuscript figure panels, but the metrics support exploratory
# checks and Supplementary Table 10 construction.
front_blood <- dat %>% filter(Cohort == "Frontline", !is.na(landmark)) %>% filter(!is.na(Blood_plus_fragment_call))
non_blood <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(Blood_plus_fragment_call)) %>% filter(!timepoint_info %in% c("Diagnosis", "Baseline"))


make_ct <- function(df,
                    pred  = "Blood_plus_fragment_call",   # cfWGS
                    truth = "MRD_truth") {                          # reference
  out <- df %>%
    select(all_of(c(pred, truth))) %>%
    drop_na() %>%                                   # keep only paired calls
    mutate(across(everything(), as.integer)) %>%    # ensure 0/1 integers
    count(!!sym(pred), !!sym(truth), name = "n") %>%
    complete(!!sym(pred) := 0:1,
             !!sym(truth):= 0:1,
             fill = list(n = 0)) %>%                # add missing cells
    mutate(row = ifelse(!!sym(pred)==1, "Pred_Pos", "Pred_Neg"),
           col = ifelse(!!sym(truth)==1,"Truth_Pos","Truth_Neg")) %>%
    select(row, col, n)
  
  # add summary columns
  tp <- out %>% filter(row=="Pred_Pos",  col=="Truth_Pos")  %>% pull(n)
  fp <- out %>% filter(row=="Pred_Pos",  col=="Truth_Neg")  %>% pull(n)
  fn <- out %>% filter(row=="Pred_Neg",  col=="Truth_Pos")  %>% pull(n)
  tn <- out %>% filter(row=="Pred_Neg",  col=="Truth_Neg")  %>% pull(n)
  
  N      <- tp + fp + fn + tn
  Agree  <- tp + tn
  tibble(
    TP = tp, FP = fp, FN = fn, TN = tn,
    N = N, Agree = Agree,
    Concordance = Agree / N,
    Sensitivity = tp/(tp+fn),
    Specificity = tn/(tn+fp),
    PPV = tp/(tp+fp),
    NPV = tn/(tn+fn)
  )
}

# 1.  FRONTLINE – post-ASCT -----------------------------------------------
ct_post_ASCT <- front_blood %>%
  filter(landmark == "Post_ASCT") %>%
  make_ct()

# 2.  FRONTLINE – 1-year maintenance --------------------------------------
ct_maint <- front_blood %>%
  filter(landmark == "Maintenance") %>%
  make_ct()

# 3.  NON-FRONTLINE – all baseline / follow-up samples --------------------
ct_nonfront <- non_blood %>% make_ct()

# Export compact combined-model contingency-table workbook for traceability.
writexl::write_xlsx(
  list(
    "Contingency_Post_ASCT"     = ct_post_ASCT,
    "Contingency_Maintenance"   = ct_maint,
    "Contingency_NonFrontline"  = ct_nonfront
  ),
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables_blood2_fragment_combined.xlsx")
)

# Build comparator-specific combined-model contingency tables.
ct_post_Flow   <- front_blood  %>% filter(landmark=="Post_ASCT") %>% 
  make_ct(truth = "Flow_Binary")
ct_post_Clono  <- front_blood  %>% filter(landmark=="Post_ASCT") %>% 
  make_ct(truth = "Adaptive_Binary")

ct_maint_Flow  <- front_blood  %>% filter(landmark=="Maintenance") %>% 
  make_ct(truth = "Flow_Binary")
ct_maint_Clono <- front_blood  %>% filter(landmark=="Maintenance") %>% 
  make_ct(truth = "Adaptive_Binary")

ct_non_Flow    <- non_blood %>% make_ct(truth = "Flow_Binary")

writexl::write_xlsx(
  list(
    "Post_ASCT_vs_Flow"        = ct_post_Flow,
    "Post_ASCT_vs_clonoSEQ"    = ct_post_Clono,
    "Maintenance_vs_Flow"      = ct_maint_Flow,
    "Maintenance_vs_clonoSEQ"  = ct_maint_Clono,
    "NonFront_vs_Flow"         = ct_non_Flow
  ),
  path = file.path(outdir, "cfWGS_contingency_vs_Flow_clonoSEQ_blood_calls_3_fragment_combined_model.xlsx") ### these are the metrics
)


# Format an analyst-facing prose check using rounded manuscript-style
# percentages.
fmt_pct <- function(x) sprintf("%.0f%%", 100*x)

post_sentence <- glue(
  "At post-ASCT, confirmatory cfDNA-based MRD demonstrated strong concordance with clonoSEQ ",
  "({ct_post_Clono$Agree}/{ct_post_Clono$N}, {fmt_pct(ct_post_Clono$Concordance)}; ",
  "PPV {fmt_pct(ct_post_Clono$PPV)}, NPV {fmt_pct(ct_post_Clono$NPV)}) ",
  "and moderate concordance with MFC ",
  "({ct_post_Flow$Agree}/{ct_post_Flow$N}, {fmt_pct(ct_post_Flow$Concordance)}; ",
  "PPV {fmt_pct(ct_post_Flow$PPV)}, NPV {fmt_pct(ct_post_Flow$NPV)})."
)

post_sentence

maint_sentence <- glue(
  "At maintenance, confirmatory cfDNA-based MRD demonstrated strong concordance with clonoSEQ ",
  "({ct_maint_Clono$Agree}/{ct_maint_Clono$N}, {fmt_pct(ct_maint_Clono$Concordance)}; ",
  "PPV {fmt_pct(ct_maint_Clono$PPV)}, NPV {fmt_pct(ct_maint_Clono$NPV)}) ",
  "and moderate concordance with MFC ",
  "({ct_maint_Flow$Agree}/{ct_maint_Flow$N}, {fmt_pct(ct_maint_Flow$Concordance)}; ",
  "PPV {fmt_pct(ct_maint_Flow$PPV)}, NPV {fmt_pct(ct_maint_Flow$NPV)})."
)

maint_sentence





# Supplementary Table 10: all cfWGS/fragmentomics call metrics against clinical
# comparators across frontline landmarks and non-frontline/test-cohort samples.
## ── Config ─────────────────────────────────────────────────────────────
TP_FRONTLINE <- c("Post_ASCT", "Maintenance")
COMPARATORS  <- c(MFC = "Flow_Binary", clonoSEQ = "Adaptive_Binary")
pct0 <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.0f%%", 100*x))
pretty_pred <- function(col) {
  if (col %in% names(techs)) {
    techs[[col]]
  } else {
    col
  }
}

# Human-readable labels for the model/call columns exported in Supplementary
# Table 10.
techs <- c(
  BM_zscore_only_sites_call            = "BM Sites Z-score",
  BM_zscore_only_detection_rate_call   = "BM cVAF Z-score",
  BM_rate_only_call                    = "BM cVAF",
  BM_base_call                         = "BM All Mut Features",
  BM_base_zscore_call                  = "BM Sites + cVAF Z-score",
  BM_plus_fragment_call                = "BM + Fragmentomics",
  BM_plus_fragment_min_call            = "BM + Fragments (min)",
  BM_base_zscore_screen_call           = "BM Sites + cVAF Z-score (screening)",
  
  Blood_zscore_only_sites_call         = "Blood Sites Z-score",
  Blood_zscore_only_detection_rate_call= "Blood cVAF Z-score",
  Blood_rate_only_call                 = "Blood cVAF",
  Blood_base_call                      = "Blood All Mut Features",
  Blood_base_zscore_call               = "Blood Sites + cVAF Z-score",
  Blood_plus_fragment_call             = "Blood + Fragmentomics",
  Blood_plus_fragment_min_call         = "Blood + Fragments (min)",
  
  Fragmentomics_full_call              = "Fragmentomics: FS + MeanCov + TF + PropShort",
  Fragmentomics_min_call               = "Fragmentomics: FS + MeanCov",
  Fragmentomics_FS_only_call           = "Fragmentomics: Fragment Size only",
  Fragmentomics_mean_coverage_only_call= "Fragmentomics: Mean coverage only",
  Fragmentomics_prop_short_only_call   = "Fragmentomics: Prop. short fragments only",
  Fragmentomics_tumor_fraction_only_call = "Fragmentomics: Tumor fraction only"
)


## ── Core contingency (one row) ─────────────────────────────────────────
compute_ct <- function(df, pred_col, truth_col) {
  if (!all(c(pred_col, truth_col) %in% names(df))) return(NULL)

  # Build an explicit 2 x 2 table on numeric 0/1 calls.  The original
  # table()/complete() implementation was sensitive to retained column names
  # and factor levels when called from Rscript; this version makes the
  # prediction/truth roles unambiguous and preserves zero-count cells.
  dd <- df %>%
    transmute(
      Pred = as.integer(.data[[pred_col]]),
      Truth = as.integer(.data[[truth_col]])
    ) %>%
    filter(!is.na(Pred), !is.na(Truth))

  if (nrow(dd) == 0) return(NULL)

  tbl <- dd %>%
    count(Pred, Truth, name = "Freq") %>%
    complete(Pred = 0:1, Truth = 0:1, fill = list(Freq = 0)) %>%
    arrange(Pred, Truth)
  
  tp <- tbl$Freq[tbl$Pred == 1 & tbl$Truth == 1]
  fp <- tbl$Freq[tbl$Pred == 1 & tbl$Truth == 0]
  fn <- tbl$Freq[tbl$Pred == 0 & tbl$Truth == 1]
  tn <- tbl$Freq[tbl$Pred == 0 & tbl$Truth == 0]
  
  N <- tp + fp + fn + tn; Agree <- tp + tn
  tibble(
    N = N, Agree = Agree,
    Concordance = ifelse(N > 0, Agree/N, NA_real_),
    TP = tp, FP = fp, FN = fn, TN = tn,
    Sensitivity = ifelse((tp + fn) > 0, tp/(tp+fn), NA_real_),
    Specificity = ifelse((tn + fp) > 0, tn/(tn+fp), NA_real_),
    PPV = ifelse((tp + fp) > 0, tp/(tp+fp), NA_real_),
    NPV = ifelse((tn + fn) > 0, tn/(tn+fn), NA_real_)
  )
}

## ── Build metrics for ALL *_call columns ───────────────────────────────
build_metrics_frontline_vs_nonfront <- function(dat, pred_regex = "_call$") {
  
  pred_cols <- names(dat)[grepl(pred_regex, names(dat))]
  if (length(pred_cols) == 0) stop("No *_call columns found.")
  
  out <- list()
  
  ## A) FRONTLINE: split by timepoint
  df_fl <- dat %>% filter(Cohort == "Frontline")
  for (tp in TP_FRONTLINE) {
    df_tp <- df_fl %>% filter(landmark == tp)
    if (nrow(df_tp) == 0) next
    
    for (pred in pred_cols) {
      base_df <- df_tp %>% filter(!is.na(.data[[pred]]))
      if (nrow(base_df) == 0) next
      
      for (cmp_name in names(COMPARATORS)) {
        truth_col <- COMPARATORS[[cmp_name]]
        df_pair   <- base_df %>% filter(!is.na(.data[[truth_col]]))
        if (nrow(df_pair) == 0) next
        
        ct <- compute_ct(df_pair, pred_col = pred, truth_col = truth_col)
        if (is.null(ct)) next
        
        out[[length(out) + 1]] <-
          ct %>%
          mutate(
            Cohort      = "Frontline",
            Timepoint   = tp,
            Pred_Column = pred,
            Pred_Label  = pretty_pred(pred),
            Comparator  = cmp_name,
            Agree_str   = sprintf("%d/%d", Agree, N),
            Conc_pct    = pct0(Concordance),
            PPV_pct     = pct0(PPV),
            NPV_pct     = pct0(NPV)
          ) %>%
          select(Cohort, Timepoint, Pred_Column, Pred_Label, Comparator,
                 Agree_str, Conc_pct, PPV_pct, NPV_pct,
                 N, Agree, TP, FP, FN, TN, Concordance, PPV, NPV, Sensitivity, Specificity)
      }
    }
  }
  
  ## B) NON-FRONTLINE: aggregate across timepoints (exclude Dx/Baseline)
  df_nf <- dat %>%
    filter(Cohort == "Non-frontline",
           !timepoint_info %in% c("Diagnosis","Baseline"))
  
  if (nrow(df_nf) > 0) {
    for (pred in pred_cols) {
      base_df <- df_nf %>% filter(!is.na(.data[[pred]]))
      if (nrow(base_df) == 0) next
      
      for (cmp_name in names(COMPARATORS)) {
        truth_col <- COMPARATORS[[cmp_name]]
        df_pair   <- base_df %>% filter(!is.na(.data[[truth_col]]))
        if (nrow(df_pair) == 0) next
        
        ct <- compute_ct(df_pair, pred_col = pred, truth_col = truth_col)
        if (is.null(ct)) next
        
        out[[length(out) + 1]] <-
          ct %>%
          mutate(
            Cohort      = "Non-frontline",
            Timepoint   = "All (excl. Dx/Baseline)",
            Pred_Column = pred,
            Pred_Label  = pretty_pred(pred),
            Comparator  = cmp_name,
            Agree_str   = sprintf("%d/%d", Agree, N),
            Conc_pct    = pct0(Concordance),
            PPV_pct     = pct0(PPV),
            NPV_pct     = pct0(NPV)
          ) %>%
          select(Cohort, Timepoint, Pred_Column, Pred_Label, Comparator,
                 Agree_str, Conc_pct, PPV_pct, NPV_pct,
                 N, Agree, TP, FP, FN, TN, Concordance, PPV, NPV, Sensitivity, Specificity)
      }
    }
  }
  
  if (length(out) == 0) return(tibble())
  bind_rows(out) %>% arrange(Cohort, Timepoint, Pred_Label, Comparator)
}

## ── Sentence helper (handles both cohorts) ─────────────────────────────
row_to_sentence <- function(row) {
  tp_phrase <- if (row$Cohort == "Frontline") {
    if (row$Timepoint == "Maintenance") "At maintenance"
    else if (row$Timepoint == "Post_ASCT") "At post-ASCT"
    else glue("At {row$Timepoint}")
  } else {
    "Overall in the non-frontline cohort"
  }
  glue(
    "{tp_phrase}, {row$Pred_Label} showed concordance with {row$Comparator} ",
    "({row$Agree_str}, {row$Conc_pct}; PPV {row$PPV_pct}, NPV {row$NPV_pct})."
  )
}

# Build the final Supplementary Table 10 metrics table.
metrics_tbl <- build_metrics_frontline_vs_nonfront(dat)

metrics_tbl <- metrics_tbl %>%
  mutate(Cohort = recode(Cohort,
                         "Frontline" = "Train",
                         "Non-frontline" = "Test"))

# Export workbook using the historical filename. The final manuscript copy is
# staged below as Supplementary Table 10.
writexl::write_xlsx(list("All_Call_Metrics" = metrics_tbl),
                    path = file.path(outdir, "Supplementary_Table_9_All_call_metrics_against_clinical_metrics.xlsx"))

# -------------------------------------------------------------------------
# Manuscript output: Supplementary Table 10
#
# What this is:
#   All call metrics against clinical comparators across cohort/timepoint
#   strata.
#
# Why it is here:
#   The historical script filename says Supplementary Table 9, but the audited
#   manuscript map identifies this workbook as final Supplementary Table 10.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = file.path(outdir, "Supplementary_Table_9_All_call_metrics_against_clinical_metrics.xlsx"),
  artifact_id = "STABLE10",
  role = "workbook_xlsx",
  description = "All call metrics against clinical metrics workbook used as Supplementary Table 10.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)



# ------------------------------------------------------------------------------
# SECTION: SUPPLEMENTARY TABLE 8 - DISCORDANCE AND CONCORDANCE SOURCE WORKBOOK
# ------------------------------------------------------------------------------
# This section builds the row-level BM and blood/cfDNA comparison tables used in
# Supplementary Table 8. Rows include both concordant and discordant comparisons,
# plus baseline mutation counts and selected fragmentomics quality/context fields
# that help interpret discordant cases.

# Add baseline mutation counts.
# For each patient, pull the BM and blood mutation counts from their first
# diagnosis/baseline sample, then join those baseline counts back onto every
# row for that patient.
baseline_counts <- dat %>%
  # Keep only diagnosis/baseline visits
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  # Use the earliest diagnosis/baseline date when multiple baseline-like rows
  # exist for the same patient.
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_min(order_by = Date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Select the two mutation‐count columns and rename them for clarity
  select(
    Patient,
    BM_MutCount_Baseline    = BM_Mutation_Count,
    Blood_MutCount_Baseline = Blood_Mutation_Count
  )

# Merge those baseline counts back onto the full dataset.
dat <- dat %>%
  left_join(baseline_counts, by = "Patient")

# Analyst-facing preview of the new columns. This is printed during command-line
# execution for QC but is not exported as a manuscript artifact.
dat %>%
  select(Patient, timepoint_info, BM_Mutation_Count, BM_MutCount_Baseline,
         Blood_Mutation_Count, Blood_MutCount_Baseline) %>%
  head(10)






# Add fragmentomics outlier flags using healthy-control reference ranges.
# These flags help annotate whether discordant samples also have unusual
# fragmentomics features.
fs_cutoffs_tbl <- read_csv(
  file.path(outdir, "FS_cutoffs_table.csv"),
  col_types = cols(
    Method       = col_character(),
    Lower_Cutoff = col_double(),
    Upper_Cutoff = col_double()
  )
)

mm_ranges <- read_csv(
  file.path(outdir, "HC_Ranges_Selected_MM_DARs_Metrics.csv"),
  col_types = cols(
    metric           = col_character(),
    mean_value       = col_double(),
    sd_value         = col_double(),
    lower_gaussian   = col_double(),
    upper_gaussian   = col_double(),
    lower_empirical  = col_double(),
    upper_empirical  = col_double()
  )
)

# Turn cutoffs into named vectors for easy lookup during mutation.
lower_fs <- fs_cutoffs_tbl$Lower_Cutoff
upper_fs <- fs_cutoffs_tbl$Upper_Cutoff

# Use the Gaussian healthy-control bounds from the MM-DARs reference table.
emp_ranges <- mm_ranges %>%
  select(metric, lower = lower_gaussian, upper = upper_gaussian)

lower_mm <- setNames(emp_ranges$lower, emp_ranges$metric)
upper_mm <- setNames(emp_ranges$upper, emp_ranges$metric)

# Flag outliers in the working data table.
dat <- dat %>%
  # 3a. FS outlier
  mutate(
    FS_outlier = case_when(
      !is.na(FS) & (FS < lower_fs | FS > upper_fs) ~ 1L,
      TRUE                                         ~ 0L
    )
  ) %>%
  # 3b. MM-DARs outliers for each metric
  mutate(
    across(
      .cols = c("Mean.Coverage", "Midpoint.Coverage",
                "Midpoint.normalized", "Amplitude"),
      .fns  = ~ as.integer(.x < lower_mm[cur_column()] |
                             .x > upper_mm[cur_column()]),
      .names = "{.col}_outlier"
    )
  )

# Save the augmented data table for traceability and downstream auditing.
readr::write_csv(
  dat,
  file.path(outdir, "dat_with_fragment_and_DARs_outlier_flags.csv")
)


# Columns retained in the discordance workbook for sample identity and
# interpretation.
id_cols   <- c("Patient", "Date", "Sample_Code", "Timepoint", "timepoint_info")

# Columns that help explain why cfWGS and clinical calls differ.
aux_cols  <- c("Adaptive_Frequency",              # clonoSEQ cumulative VAF
               "Flow_pct_cells", grep("^BM.*(_prob|_call)$", names(dat), value = TRUE),
               "FS", "Mean.Coverage", "detect_rate_BM", "zscore_BM", 
               "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
               "BM_MutCount_Baseline", "Blood_MutCount_Baseline", "FS_outlier", "Mean.Coverage_outlier")

# Warn if any expected explanatory columns are absent.
missing <- setdiff(aux_cols, names(dat))
if (length(missing)) warning("These columns are missing: ", paste(missing, collapse = ", "))
# BM-informed frontline comparison table for Supplementary Table 8.
combined_discord_tbl <- dat %>%
  filter(
    !is.na(BM_zscore_only_detection_rate_call),
    !is.na(landmark),                             # only landmark timepoints
    (!is.na(Adaptive_Binary) | !is.na(Flow_Binary))
  ) %>%
  pivot_longer(
    cols      = c(Adaptive_Binary, Flow_Binary),
    names_to  = "Comparator",
    values_to = "Reference"
  ) %>%
  filter(!is.na(Reference)) %>%
  mutate(
    # Label the comparator
    Comparator = recode(
      Comparator,
      Adaptive_Binary = "clonoSEQ",
      Flow_Binary     = "MFC"
    ),
    # Define discordance/concordance category
    category = case_when(
      BM_zscore_only_detection_rate_call == 1L & Reference == 0L ~
        paste0("cfWGS_pos / ", Comparator, "_neg"),
      BM_zscore_only_detection_rate_call == 0L & Reference == 1L ~
        paste0("cfWGS_neg / ", Comparator, "_pos"),
      TRUE ~ "concordant"
    )
  ) %>%
  select(
    all_of(id_cols),
    Cohort,
    landmark,
    Comparator,
    category,
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)
  ) %>%
  arrange(landmark, Patient, Comparator)

combined_discord_tbl_slim <- combined_discord_tbl %>% filter(category != "concordant")

# Export a CSV companion for audit/review.
write.csv(
  combined_discord_tbl,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_BM2.csv"),
  row.names = FALSE
)



# BM-informed non-frontline/test-cohort comparison table for Supplementary
# Table 8, aggregated across eligible timepoints.
combined_discord_tbl_non_frontline <- dat %>%
  filter(
    !is.na(BM_zscore_only_detection_rate_call),
    (!is.na(Adaptive_Binary) | !is.na(Flow_Binary))
  ) %>%
  filter(Cohort == "Non-frontline") %>%
  pivot_longer(
    cols      = c(Adaptive_Binary, Flow_Binary),
    names_to  = "Comparator",
    values_to = "Reference"
  ) %>%
  filter(!is.na(Reference)) %>%
  mutate(
    # Label the comparator
    Comparator = recode(
      Comparator,
      Adaptive_Binary = "clonoSEQ",
      Flow_Binary     = "MFC"
    ),
    # Define discordance/concordance category
    category = case_when(
      BM_zscore_only_detection_rate_call == 1L & Reference == 0L ~
        paste0("cfWGS_pos / ", Comparator, "_neg"),
      BM_zscore_only_detection_rate_call == 0L & Reference == 1L ~
        paste0("cfWGS_neg / ", Comparator, "_pos"),
      TRUE ~ "concordant"
    )
  ) %>%
  select(
    all_of(id_cols),
    Cohort,
    Comparator,
    category,
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)
  ) %>%
  arrange(Patient, Comparator)


# Export a CSV companion for audit/review.
write.csv(
  combined_discord_tbl_non_frontline,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_BM_non_frontline2.csv"),
  row.names = FALSE
)





# Supplementary analysis: exploratory association checks for BM-informed
# discordance direction. These regression summaries are printed for analyst QC
# and are not staged as final manuscript tables.

# 1)  Flag discordance (1 = discordant, 0 = concordant) ----------------------
tbl <- combined_discord_tbl 

# Numeric predictors assessed in the exploratory discordance models.
num_vars <- c(
  "Adaptive_Frequency",
  "Flow_pct_cells",
  "detect_rate_BM",
  "zscore_BM",
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
  "FS",
  "Mean.Coverage",
  "BM_MutCount_Baseline",
  "Blood_MutCount_Baseline"
)

# 2) Tag each row as 'missed', 'captured', or 'concordant'
tbl <- tbl %>%
  mutate(
    direction = case_when(
      str_detect(category, "^cfWGS_neg") ~ "missed",
      str_detect(category, "^cfWGS_pos") ~ "captured",
      TRUE                                ~ "concordant"
    ),
    # for binary modeling
    is_missed   = as.integer(direction == "missed"),
    is_captured = as.integer(direction == "captured")
  )


# Descriptive summary by discordance direction.
tbl %>%
  group_by(direction) %>%
  summarise(
    across(all_of(num_vars),
           list(median = ~ median(.x, na.rm=TRUE),
                IQR    = ~ IQR(.x,   na.rm=TRUE)),
           .names = "{.col}_{.fn}")
  ) %>%
  print(n = Inf)

# 4) Helper to fit a logistic model for a given binary outcome
fit_disc_model <- function(outcome) {
  formula <- as.formula(paste(outcome, "~", paste(num_vars, collapse = " + ")))
  glm(formula, data = tbl %>% drop_na(all_of(num_vars)), family = binomial)
}

# Sparse discordance strata can produce complete/quasi-complete separation.
# Profile-likelihood intervals are preferred when available, but they can fail
# when there are too few non-NA support points.  In that case we keep the model
# summary command-line-runnable by reporting Wald intervals and recording the
# interval method explicitly.
safe_tidy_glm <- function(model, conf.level = 0.95, exponentiate = TRUE) {
  profiled <- tryCatch(
    broom::tidy(model, conf.int = TRUE, conf.level = conf.level, exponentiate = exponentiate),
    error = function(e) NULL
  )

  if (!is.null(profiled)) {
    return(profiled %>% mutate(ci_method = "profile"))
  }

  z <- qnorm((1 + conf.level) / 2)
  wald <- broom::tidy(model, conf.int = FALSE, exponentiate = FALSE) %>%
    mutate(
      conf.low = estimate - z * std.error,
      conf.high = estimate + z * std.error,
      ci_method = "wald_fallback"
    )

  if (isTRUE(exponentiate)) {
    wald <- wald %>%
      mutate(
        estimate = exp(estimate),
        conf.low = exp(conf.low),
        conf.high = exp(conf.high)
      )
  }

  wald
}

# Fit pooled logistic models. Interpret cautiously because discordance strata
# are small and may show complete/quasi-complete separation.
glm_missed   <- fit_disc_model("is_missed")   # 1 = missed, 0 = else
glm_captured <- fit_disc_model("is_captured") # 1 = captured, 0 = else

# 6) Tidy results with odds‐ratios and 95% CIs
tidy_missed   <- safe_tidy_glm(glm_missed)
tidy_captured <- safe_tidy_glm(glm_captured)

print(tidy_missed)
print(tidy_captured)

# Optional comparator-stratified exploratory models.
by_comp <- tbl %>%
  group_by(Comparator) %>%
  nest() %>%
  mutate(
    missed_mod   = map(data, ~ glm(is_missed   ~ ., data = select(.x, all_of(num_vars), is_missed),   family=binomial)),
    captured_mod = map(data, ~ glm(is_captured ~ ., data = select(.x, all_of(num_vars), is_captured), family=binomial)),
    missed_tidy   = map(missed_mod, safe_tidy_glm),
    captured_tidy = map(captured_mod, safe_tidy_glm)
  ) %>%
  select(Comparator, missed_tidy, captured_tidy) %>%
  unnest(c(missed_tidy, captured_tidy), names_sep = "_")

print(by_comp)




# Blood/cfDNA-informed frontline comparison table for Supplementary Table 8.
aux_cols  <- c("Adaptive_Frequency",              # clonoSEQ cumulative VAF
               "Flow_pct_cells",                  # MFC percent tumor cells
               "FS", "Mean.Coverage",  grep("^Blood.*(_prob|_call)$", names(dat), value = TRUE),
               "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
               "BM_MutCount_Baseline", "Blood_MutCount_Baseline", "FS_outlier", "Mean.Coverage_outlier")

combined_discord_tbl2 <- dat %>%
  filter(
    !is.na(Blood_zscore_only_sites_call),
    !is.na(landmark),                             # only landmark timepoints
    (!is.na(Adaptive_Binary) | !is.na(Flow_Binary))
  ) %>%
  pivot_longer(
    cols      = c(Adaptive_Binary, Flow_Binary),
    names_to  = "Comparator",
    values_to = "Reference"
  ) %>%
  filter(!is.na(Reference)) %>%
  mutate(
    # Label the comparator
    Comparator = recode(
      Comparator,
      Adaptive_Binary = "clonoSEQ",
      Flow_Binary     = "MFC"
    ),
    # Define discordance/concordance category
    category = case_when(
      Blood_zscore_only_sites_call == 1L & Reference == 0L ~
        paste0("cfWGS_pos / ", Comparator, "_neg"),
      Blood_zscore_only_sites_call == 0L & Reference == 1L ~
        paste0("cfWGS_neg / ", Comparator, "_pos"),
      TRUE ~ "concordant"
    )
  ) %>%
  select(
    all_of(id_cols),
    Cohort,
    landmark,
    Comparator,
    category,
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)
  ) %>%
  arrange(landmark, Patient, Comparator)

combined_discord_tbl_slim <- combined_discord_tbl2 %>% filter(category != "concordant")

# Export a CSV companion for audit/review.
write.csv(
  combined_discord_tbl2,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_Blood2.csv"),
  row.names = FALSE
)



# Blood/cfDNA-informed non-frontline/test-cohort comparison table for
# Supplementary Table 8, aggregated across eligible timepoints.
combined_discord_tbl_non_frontline2 <- dat %>%
  filter(
    !is.na(Blood_zscore_only_sites_call),
    (!is.na(Adaptive_Binary) | !is.na(Flow_Binary))
  ) %>%
  filter(Cohort == "Non-frontline") %>%
  pivot_longer(
    cols      = c(Adaptive_Binary, Flow_Binary),
    names_to  = "Comparator",
    values_to = "Reference"
  ) %>%
  filter(!is.na(Reference)) %>%
  mutate(
    # Label the comparator
    Comparator = recode(
      Comparator,
      Adaptive_Binary = "clonoSEQ",
      Flow_Binary     = "MFC"
    ),
    # Define discordance/concordance category
    category = case_when(
      Blood_zscore_only_sites_call == 1L & Reference == 0L ~
        paste0("cfWGS_pos / ", Comparator, "_neg"),
      Blood_zscore_only_sites_call == 0L & Reference == 1L ~
        paste0("cfWGS_neg / ", Comparator, "_pos"),
      TRUE ~ "concordant"
    )
  ) %>%
  select(
    all_of(id_cols),
    Cohort,
    Comparator,
    category,
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)
  ) %>%
  arrange(Patient, Comparator)


# Export a CSV companion for audit/review.
write.csv(
  combined_discord_tbl_non_frontline2,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_blood_non_frontline2.csv"),
  row.names = FALSE
)

# Create the final multi-sheet Supplementary Table 8 workbook.
library(openxlsx)

# Load de-identified patient ID map and baseline dates for exported workbook
# columns.
id_map <- readRDS("id_map.rds") %>% distinct(Patient, New_ID)

Baseline_dates <- read_csv("Final Tables and Figures/Baseline dates for samples.csv")
# Make a joinable version using the same Patient key as the exported workbook
# when a de-identified New_ID is available.
baseline_join <- Baseline_dates %>%
  # assume columns: patient (original ID) and start (baseline date); adapt if needed
  mutate(start = as_date(start)) %>%
  left_join(id_map, by = c("patient" = "Patient")) %>%
  mutate(Patient = coalesce(New_ID, patient)) %>%
  select(Patient, start)

# Prepare each sheet: rename clinical assay columns, remap patient IDs, add
# days/months from diagnosis, and remove raw date/sample-code fields before
# export.
prepare_tbl <- function(df, id_map, baseline_join){
  if (is.null(df)) stop("Input table is NULL.")
  
  df %>%
    tibble::as_tibble() %>%
    ungroup() %>%
    rename(
      clonoSEQ_Tumor_Ig_Frequency = any_of("Adaptive_Frequency"),
      MFC_Pct_Tumor_Cells         = any_of("Flow_pct_cells")
    ) %>%
    # Map Patient -> New_ID (like before)
    { if ("Patient" %in% names(.)) {
      left_join(., id_map, by = "Patient") %>%
        mutate(Patient = coalesce(New_ID, Patient)) %>%
        select(-any_of("New_ID"))
    } else . } %>%
    # Add baseline start date and compute days_from_dx when Date exists
    left_join(baseline_join, by = "Patient") %>%
    mutate(
      Date = if ("Date" %in% names(.)) as_date(Date) else Date,
      days_from_dx = if ("Date" %in% names(.))
        as.integer(difftime(Date, start, units = "days"))
      else NA_integer_,
      months_from_dx = if ("Date" %in% names(.))
        as.numeric(interval(start, Date) / months(1))
      else NA_real_
    ) %>%
    relocate(any_of(c("days_from_dx","months_from_dx")), .after = "Patient") %>%
    # flatten any list-cols for Excel safety
    mutate(across(where(is.list), ~ map_chr(., ~ paste0(as.character(.x), collapse = "; ")))) %>%
    # remove identifiers you don’t want to export
    select(-any_of(c("Sample_Code", "Date", "start"))) %>%
    as.data.frame(check.names = FALSE)
}

# Build the four workbook sheets: BM/blood by training/test cohort.
BM_Train    <- prepare_tbl(combined_discord_tbl,                id_map, baseline_join)
BM_Test     <- prepare_tbl(combined_discord_tbl_non_frontline,  id_map, baseline_join)
Blood_Train <- prepare_tbl(combined_discord_tbl2,               id_map, baseline_join)
Blood_Test  <- prepare_tbl(combined_discord_tbl_non_frontline2, id_map, baseline_join)

# Write the final workbook with filters and readable column widths.
add_sheet_with_style <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  if (is.null(ncol(data)) || is.na(ncol(data)) || ncol(data) == 0) {
    writeData(wb, sheet_name, data.frame(`(empty)` = character(0)))
    return(invisible(NULL))
  }
  writeData(wb, sheet_name, data, headerStyle = createStyle(textDecoration = "bold"))
  addFilter(wb, sheet = sheet_name, rows = 1, cols = 1:ncol(data))
  setColWidths(wb, sheet = sheet_name, cols = 1:ncol(data), widths = "auto")
}

wb <- createWorkbook()
add_sheet_with_style(wb, "BM_Train",    BM_Train)
add_sheet_with_style(wb, "BM_Test",     BM_Test)
add_sheet_with_style(wb, "Blood_Train", Blood_Train)
add_sheet_with_style(wb, "Blood_Test",  Blood_Test)

saveWorkbook(wb, "Final Tables and Figures/Supplementary_Table_8_model_comparisons_to_clinical_metrics3.xlsx",
             overwrite = TRUE)

# -------------------------------------------------------------------------
# Manuscript output: Supplementary Table 8
#
# What this is:
#   Multi-sheet workbook comparing cfWGS model outputs to clinical metrics,
#   including discordance/concordance summaries used with Main Figure 3E and
#   Main Figure 4D.
#
# Why it is here:
#   This is the code-generated workbook mapped to final Supplementary Table 8.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supplementary_Table_8_model_comparisons_to_clinical_metrics3.xlsx",
  artifact_id = "STABLE8",
  role = "workbook_xlsx",
  description = "Model comparisons to clinical metrics workbook used as Supplementary Table 8.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)

# ------------------------------------------------------------------------------
# SECTION: MAIN FIGURE 3E AND MAIN FIGURE 4D - PROBABILITY VS CLINICAL ASSAYS
# ------------------------------------------------------------------------------
# This section compares cfWGS model probabilities with quantitative clinical
# MRD assay values. The historical working filenames include "Fig4K/Fig5K",
# but the final manuscript-staged outputs are:
#   - Figure 3E: BM-informed cfWGS probability vs MFC, clonoSEQ, and EasyM
#   - Figure 4D: blood/cfDNA-informed cfWGS probability vs MFC, clonoSEQ, and EasyM
#
# LOD/reference constants used for plotting.
lod_cfWGS   <- 0.00011   # 0.011 % 
lod_cfWGS_blood   <- 0.00061   # 0.061 %

lod_clonoMF <- 1e-5      # 10-5  (for both MFC & clonoSEQ)

shape_pal <- c(
  detected     = 16,    # filled circle
  `not detected` = 4    # open cross
)

# Figure 3E input: build BM-informed plotting data for frontline landmark
# samples with at least one quantitative clinical comparator.
plot_df <- dat %>%
  mutate(Flow_pct_cells = Flow_pct_cells/100) %>% # to be consistent
  filter(
    Cohort == "Frontline",
    !is.na(BM_zscore_only_detection_rate_call),
    !is.na(landmark_tp),
    !is.na(Adaptive_Frequency) | !is.na(Flow_pct_cells)
  ) %>%
  pivot_longer(
    cols      = c(Adaptive_Frequency, Flow_pct_cells),
    names_to  = "Comparator",
    values_to = "x_val"
  ) %>%
  drop_na(x_val) %>%
  rowwise() %>%
  mutate(
    # Recode comparator names to match plot labels and facet levels.
    Comparator = recode(Comparator,
                        Adaptive_Frequency = "clonoSEQ",
                        Flow_pct_cells     = "MFC"),
    
    # Build reference binary calls from the original assay-specific binary
    # columns, while treating clonoSEQ zero frequency as not detected.
    ref_binary = case_when(
      Comparator == "clonoSEQ" & x_val == 0 ~ 0L,
      Comparator == "clonoSEQ" ~ Adaptive_Binary,
      Comparator == "MFC"      ~ Flow_Binary,
      TRUE                     ~ NA_integer_
    ),
    
    cfwgs_bin = BM_zscore_only_detection_rate_call,
    concord   = (cfwgs_bin == ref_binary),
    detected  = if_else(BM_zscore_only_detection_rate_call == 1,
                        "detected", "not detected")
  ) %>%
  ungroup() %>%
  mutate(
    x_plot   = if_else(x_val <= 1e-6, 1e-6, x_val),
    y_plot   = if_else(BM_zscore_only_detection_rate_prob <= 1e-5, 1e-5, BM_zscore_only_detection_rate_prob),
    category = if_else(concord, "concordant", "discordant")
  ) %>% 
  select(
    Patient, Sample_Code, landmark_tp, Comparator, x_val, x_plot,
    BM_zscore_only_detection_rate_call, y_plot,
    ref_binary, cfwgs_bin, concord, detected, category,
    Num_days_to_closest_relapse, Flow_Binary, Adaptive_Binary
  )

# Add detection-pattern and relapse-within-1-year labels for plotting.
plot_df2 <- plot_df %>%
  mutate(
    shape_cat = case_when(
      cfwgs_bin == 1 & ref_binary == 0 ~ "cfWGS only",
      cfwgs_bin == 0 & ref_binary == 1 ~ "Comparator only",
      cfwgs_bin == 1 & ref_binary == 1 ~ "Both",
      cfwgs_bin == 0 & ref_binary == 0 ~ "Neither"
    ),   # color by relapse status
    relapse_cat = if_else(
      Num_days_to_closest_relapse <= 365,
      "Relapsed ≤365 d",
      "No relapse ≤365 d",
      missing = "No relapse ≤365 d"
    )
  )

# Normalize labels before plotting.
hyphen_rx <- "[\u2010\u2011\u2012\u2013\u2014\u2212]"  # all the usual dash culprits

plot_df2 <- plot_df2 %>%
  mutate(
    # Normalize dashes and whitespace.
    landmark_tp = str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
    
    # Factor with ASCII hyphens so facet order is deterministic.
    landmark_timepoint = factor(
      landmark_tp,
      levels = c("Post-ASCT", "Maintenance-1yr")
    ),
    
    detect_cat = factor(shape_cat,
                        levels = c("cfWGS only", "Comparator only", "Both", "Neither")),
    relapse_cat = factor(relapse_cat,
                         levels = c("Relapsed ≤365 d", "No relapse ≤365 d")),
    relapse_flag = if_else(relapse_cat == "Relapsed ≤365 d", "Relapse", "No relapse")
  )


# QC check: any rows with missing detection-pattern labels should be reviewed.
plot_df2 %>% filter(is.na(shape_cat)) %>% select(Patient, Comparator, cfwgs_bin, ref_binary)

# Historical detection-pattern palette retained for provenance; the active plot
# below uses relapse status for point fill.
detect_cols <- c(
  `cfWGS only`      = "#1b9e77",
  `Comparator only` = "#d95f02",
  `Both`            = "#7570b3",
  `Neither`         = "#999999"
)

# Viridis palette retained from the original plotting workflow.
library(viridisLite)

pal5 <- viridis(5, option = "D") 

# Assign colors to detection-pattern categories.
detect_cols <- c(
  `cfWGS only`      = pal5[1],   "#440154FF",
  `Comparator only` = pal5[2],   "#31688EFF",
  `Both`            = pal5[4],   "#73D055FF",
  `Neither`         = "#999999"   # "#FDE725FF"
)


# Correlations annotated in the BM-informed clinical-comparator panel.
corr_df <- plot_df2 %>%
  group_by(landmark_timepoint, Comparator) %>%
  summarize(
    rho = cor(x_plot, y_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(x_plot, y_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,
    y = 0.99
  )


# Build the historical two-comparator BM panel. The final manuscript panel adds
# EasyM below and is staged as Figure 3E.
p_scatter_simple <- ggplot(plot_df2,
                           aes(x = x_plot, y = y_plot,
                               fill   = relapse_cat)) +  
  # LOD reference lines
  geom_hline(yintercept = 0.4215524,   linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = lod_clonoMF, linetype = "dashed", colour = "grey80") +
  
  # points
  geom_point(shape = 21, size = 2, alpha = 0.9,
             colour = "black") +   # outline colour (same for all)
  
  # colour legend for detection pattern
  scale_fill_manual(name = "Detection pattern",
                    values = detect_cols) +
  
  # log axes
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0.12, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator)) +
  
  labs(title = "cfWGS of BM-Derived Mutations MRD\nProbability vs. Clinical Assays",
       x = "Comparator MRD level",
       y = "cVAF Model Probability") +
  # start from a white‐background theme with borders
  theme_bw(base_size = 11) +    
  
  # Color points by relapse within one year.
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black",
      "Unknown"           = "#bbbbbb"
    )
  ) +
  
  # Add Spearman correlation labels.
  geom_text(
    data = corr_df,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +
  # Final theme adjustments.
  theme(
    # draw a thin black border around each facet
    panel.border      = element_rect(colour = "black", fill = NA, size = 0.5),
    
    # kill all the internal grid lines
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    
    # if you want a white strip header with black outline:
    strip.background  = element_rect(fill = "white", colour = "black"),
    strip.text        = element_text(face = "bold"),
    
    # axis titles bold
    axis.title        = element_text(size = 11),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  ) 

p_scatter_simple

# Save the historical two-comparator BM panel for provenance.
ggsave("Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ_clean_BM_muts_updated4.png",
       p_scatter_simple,
       width  = 6.5, height = 5, dpi = 600)


# ===========================================================================
# SECTION 3B: SCATTER PLOTS WITH EASYМ ADDED AS THIRD COLUMN (BM-DERIVED)
# ===========================================================================
# Final Figure 3E version: add EasyM as a third comparator column.

# Load EasyM data and merge into the working cfWGS table.
cat("\n=== Loading EasyM data for scatter plots ===\n")
EasyM_data_file <- file.path(OUTPUT_DIR_EASYM, "EasyM_all_samples_with_optimized_calls.csv")

if (file.exists(EasyM_data_file)) {
  EasyM_dat <- readr::read_csv(EasyM_data_file, show_col_types = FALSE)
  
  # Remove any existing EasyM_value columns before joining to avoid .x/.y
  # suffixes in command-line runs.
  dat_prep <- dat %>%
    select(-starts_with("EasyM_value"))
  
  dat_with_easym <- dat_prep %>%
    left_join(
      EasyM_dat %>%
        select(Patient, Timepoint, EasyM_value) %>%
        distinct(Patient, Timepoint, .keep_all = TRUE),
      by = c("Patient" = "Patient", "Timepoint" = "Timepoint"),
      relationship = "many-to-one"
    )
  
  cat(sprintf("✓ EasyM data merged: %d rows with EasyM_value\n", 
              sum(!is.na(dat_with_easym$EasyM_value))))
} else {
  warning("EasyM data file not found - EasyM columns will be empty")
  # If EasyM data are absent, create an all-NA column so downstream code fails
  # gracefully by producing empty EasyM facets rather than stopping.
  dat_with_easym <- dat %>%
    select(-starts_with("EasyM_value")) %>%
    mutate(EasyM_value = NA_real_)
}

# Load EasyM thresholds by timepoint (for reference lines)
easyM_thresholds <- readr::read_csv(
  file.path(OUTPUT_DIR_EASYM, "EasyM_threshold_values_by_timepoint.csv"),
  show_col_types = FALSE
)

# Build Figure 3E plotting data with EasyM as an additional comparator.
plot_df_with_easym <- bind_rows(
  # Original clonoSEQ and MFC rows.
  plot_df2 %>%
    mutate(Comparator = as.character(Comparator)) %>%
    select(Patient, Sample_Code, landmark_timepoint, Comparator, 
           x_plot, y_plot, relapse_cat, Num_days_to_closest_relapse),
  # EasyM comparator rows.
  dat_with_easym %>%
    mutate(Flow_pct_cells = Flow_pct_cells/100) %>%
    filter(
      Cohort == "Frontline",
      !is.na(BM_zscore_only_detection_rate_call),
      !is.na(landmark_tp),
      !is.na(EasyM_value)
    ) %>%
    mutate(
      Comparator = "EasyM",
      x_val = EasyM_value,
      # Cap EasyM at 100% and set a minimum floor for log-scale plotting.
      x_plot = pmin(pmax(EasyM_value, 1e-6), 1.0),
      y_plot = if_else(BM_zscore_only_detection_rate_prob <= 1e-5, 
                       1e-5, BM_zscore_only_detection_rate_prob),
      landmark_tp = str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
      landmark_timepoint = factor(landmark_tp, levels = c("Post-ASCT", "Maintenance-1yr")),
      relapse_cat = if_else(
        Num_days_to_closest_relapse <= 365,
        "Relapsed ≤365 d",
        "No relapse ≤365 d",
        missing = "No relapse ≤365 d"
      ),
      relapse_cat = factor(relapse_cat, levels = c("Relapsed ≤365 d", "No relapse ≤365 d"))
    ) %>%
    select(Patient, Sample_Code, landmark_timepoint, Comparator, 
           x_plot, y_plot, relapse_cat, Num_days_to_closest_relapse)
)

# Set comparator and timepoint factor order so the panel layout is stable.
plot_df_with_easym <- plot_df_with_easym %>%
  mutate(
    # Trim whitespace before applying explicit factor levels.
    Comparator = str_trim(as.character(Comparator)),
    Comparator = factor(Comparator, levels = c("clonoSEQ", "MFC", "EasyM"), ordered = FALSE),
    landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
  )

# Calculate Spearman correlations for all three comparators.
corr_df_with_easym <- plot_df_with_easym %>%
  group_by(landmark_timepoint, Comparator) %>%
  summarize(
    rho = cor(x_plot, y_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(x_plot, y_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,
    y = 0.99,
    # Match comparator factor levels to the plot data.
    Comparator = factor(Comparator, levels = c("clonoSEQ", "MFC", "EasyM"))
  )

# Map EasyM thresholds by timepoint for reference lines.
# Timepoint codes: 05 = Post-ASCT (Post-Transplant), 07 = Maintenance-1yr (1-Year Maintenance)
easyM_threshold_lines <- easyM_thresholds %>%
  filter(Timepoint %in% c("05", "07")) %>%
  mutate(
    landmark_timepoint = case_when(
      Timepoint == "05" ~ "Post-ASCT",
      Timepoint == "07" ~ "Maintenance-1yr",
      TRUE ~ NA_character_
    ),
    # Convert percentage (0-100 scale) to proportion (0-1 scale).
    xintercept = Threshold_raw_percent / 100,
    Comparator = "EasyM",
    # Match comparator factor levels to the plot data.
    Comparator = factor(Comparator, levels = c("clonoSEQ", "MFC", "EasyM"))
  ) %>%
  select(landmark_timepoint, Comparator, xintercept)

# Re-apply factor ordering directly before plotting to avoid ordering drift from
# row binding or joins.
plot_df_with_easym <- plot_df_with_easym %>%
  mutate(landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr")))

corr_df_with_easym <- corr_df_with_easym %>%
  mutate(landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr")))

easyM_threshold_lines <- easyM_threshold_lines %>%
  mutate(landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr")))

# Build the final Figure 3E panel with EasyM added.
p_scatter_with_easym_bm <- ggplot(plot_df_with_easym,
                                   aes(x = x_plot, y = y_plot,
                                       fill = relapse_cat)) +
  # Reference lines for cfWGS and comparator thresholds.
  geom_hline(yintercept = 0.4215524, linetype = "dashed", colour = "grey80") +
  # Clinical-comparator LOD line for clonoSEQ and MFC facets.
  geom_vline(aes(xintercept = xintercept), linetype = "dashed", colour = "grey80",
             data = data.frame(
               landmark_timepoint = factor(rep(c("Post-ASCT", "Maintenance-1yr"), 2),
                                          levels = c("Post-ASCT", "Maintenance-1yr")),
               Comparator = factor(c("clonoSEQ", "clonoSEQ", "MFC", "MFC"), 
                                   levels = c("clonoSEQ", "MFC", "EasyM")),
               xintercept = lod_clonoMF
             )) +
  # EasyM thresholds by landmark timepoint.
  geom_vline(aes(xintercept = xintercept), linetype = "dashed", colour = "grey80",
             data = easyM_threshold_lines) +
  
  # Points
  geom_point(shape = 21, size = 2, alpha = 0.9, colour = "black") +
  
  # Color points by relapse within one year.
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black"
    )
  ) +
  
  # Log-scale x-axis with fixed labels across comparators.
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%", "100%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0.12, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  # Facet rows by timepoint and columns by comparator.
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator),
             drop = FALSE,
             as.table = TRUE) +
  
  # Annotations
  geom_text(
    data = corr_df_with_easym,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +
  
  # Labels
  labs(title = "cfWGS of BM-Derived Mutations MRD\nProbability vs. Clinical Assays and EasyM",
       x = "Comparator MRD level",
       y = "cfWGS Model Probability") +
  
  # Theme
  theme_bw(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8)
  )

p_scatter_with_easym_bm

# Save the final Figure 3E working PNG before staging it to the manuscript
# output directory.
ggsave("Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ_EasyM_BM_muts_updated5.png",
       p_scatter_with_easym_bm,
       width = 8.5, height = 5, dpi = 600)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 3E
#
# What this is:
#   BM-informed cfWGS comparison to MFC, clonoSEQ, and EasyM.
#
# Why it is here:
#   This PNG is the final Main Figure 3E component and also supports
#   Supplementary Table 8.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ_EasyM_BM_muts_updated5.png",
  artifact_id = "FIG3E",
  role = "figure_panel_png",
  description = "BM-informed cfWGS comparison to MFC, clonoSEQ, and EasyM used as Main Figure 3E.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)

# Export source data
readr::write_csv(
  plot_df_with_easym %>% mutate(Figure = "Fig4K_cfWGS_vs_clinical_assays_EasyM_BM"),
  file.path(outdir_source_data, "Fig4K_cfWGS_vs_clinical_assays_EasyM_BM_source_data.csv")
)


# ===========================================================================
# SECTION 4: AUXILIARY CROSS-PLATFORM EASYM VISUALIZATIONS
# ===========================================================================
# Generate auxiliary scatter plots comparing EasyM residual M-protein directly
# to cfWGS probabilities. These are source/QC outputs and are not currently
# staged as final manuscript figure panels.

# ---------------------------------------------------------------------------
# SUBSECTION 4.1: EasyM vs Bone Marrow cfWGS Model
# ---------------------------------------------------------------------------
# Visualize relationship between EasyM M-protein % and cfWGS BM model probability
# Points colored by relapse status, faceted by treatment timepoint

# Prepare data for EasyM vs BM cfWGS comparison
plot_df_easym_bm <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(BM_zscore_only_detection_rate_call),
    !is.na(landmark_tp),
    !is.na(EasyM_value)
  ) %>%
  mutate(
    # Set minimum detection floor for log-scale visualization
    EasyM_plot = if_else(EasyM_value <= 1e-6, 1e-6, EasyM_value),
    BM_prob_plot = if_else(BM_zscore_only_detection_rate_prob <= 1e-5, 
                           1e-5, BM_zscore_only_detection_rate_prob),
    
    # Stratify by clinical outcome: relapse within 1 year vs not
    relapse_cat = if_else(
      Num_days_to_closest_relapse <= 365,
      "Relapsed ≤365 d",
      "No relapse ≤365 d",
      missing = "No relapse ≤365 d"
    ),
    
    # Standardize landmark timepoint labels
    landmark_timepoint = factor(
      str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
      levels = c("Post-ASCT", "Maintenance-1yr")
    )
  )

# Calculate correlations
corr_easym_bm <- plot_df_easym_bm %>%
  group_by(landmark_timepoint) %>%
  summarize(
    rho = cor(EasyM_plot, BM_prob_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(EasyM_plot, BM_prob_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,
    y = 0.99
  )

# Create scatter plot
p_easym_bm <- ggplot(plot_df_easym_bm,
                     aes(x = EasyM_plot, y = BM_prob_plot, fill = relapse_cat)) +
  # Reference lines
  geom_hline(yintercept = 0.4215524, linetype = "dashed", colour = "grey80") +
  
  # Points
  geom_point(shape = 21, size = 2, alpha = 0.9, colour = "black") +
  
  # Scales
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0.12, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_wrap(~ landmark_timepoint, nrow = 2) +
  
  # Labels and theme
  labs(title = "cfWGS of BM-Derived Mutations MRD\nProbability vs. EasyM",
       x = "EasyM Proteomic MRD level",
       y = "cVAF Model Probability") +
  theme_bw(base_size = 11) +
  
  # Colors
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black",
      "Unknown"           = "#bbbbbb"
    )
  ) +
  
  # Add correlation text
  geom_text(
    data = corr_easym_bm,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +
  
  # Theme adjustments
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "white", colour = "black"),
    strip.text        = element_text(face = "bold"),
    axis.title        = element_text(size = 11),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8)
  )

p_easym_bm

# Save the auxiliary EasyM-vs-BM probability plot for provenance.
ggsave("Final Tables and Figures/FigS_EasyM_vs_BM_cfWGS_prob.png",
       p_easym_bm,
       width = 5, height = 6, dpi = 600)


# Figure 4D input: repeat the probability-vs-comparator workflow using
# blood/cfDNA-derived mutation lists and blood cfWGS model probabilities.
plot_df <- dat %>%
  mutate(Flow_pct_cells = Flow_pct_cells/100) %>% # to be consistent
  filter(
    Cohort == "Frontline",
    !is.na(Blood_zscore_only_sites_call),
    !is.na(landmark_tp),
    !is.na(Adaptive_Frequency) | !is.na(Flow_pct_cells)
  ) %>%
  pivot_longer(
    cols      = c(Adaptive_Frequency, Flow_pct_cells),
    names_to  = "Comparator",
    values_to = "x_val"
  ) %>%
  drop_na(x_val) %>%
  rowwise() %>%
  mutate(
    # Recode comparator names to match plot labels and facet levels.
    Comparator = recode(Comparator,
                        Adaptive_Frequency = "clonoSEQ",
                        Flow_pct_cells     = "MFC"),
    
    # Build reference binary calls from assay-specific binary columns, treating
    # clonoSEQ zero frequency as not detected.
    ref_binary = case_when(
      Comparator == "clonoSEQ" & x_val == 0 ~ 0L,
      Comparator == "clonoSEQ" ~ Adaptive_Binary,
      Comparator == "MFC"      ~ Flow_Binary,
      TRUE                     ~ NA_integer_
    ),
    
    cfwgs_bin = Blood_zscore_only_sites_call,
    concord   = (cfwgs_bin == ref_binary),
    detected  = if_else(Blood_zscore_only_sites_call == 1,
                        "detected", "not detected")
  ) %>%
  ungroup() %>%
  mutate(
    x_plot   = if_else(x_val <= 1e-6, 1e-6, x_val),
    y_plot   = if_else(Blood_zscore_only_sites_prob <= 1e-5, 1e-5, Blood_zscore_only_sites_prob),
    category = if_else(concord, "concordant", "discordant")
  ) %>% 
  select(
    Patient, Sample_Code, landmark_tp, Comparator, x_val, x_plot,
    Blood_zscore_only_sites_call, y_plot,
    ref_binary, cfwgs_bin, concord, detected, category,
    Num_days_to_closest_relapse, Flow_Binary, Adaptive_Binary
  )

# Add detection-pattern and relapse-within-1-year labels for plotting.
plot_df2 <- plot_df %>%
  mutate(
    shape_cat = case_when(
      cfwgs_bin == 1 & ref_binary == 0 ~ "cfWGS only",
      cfwgs_bin == 0 & ref_binary == 1 ~ "Comparator only",
      cfwgs_bin == 1 & ref_binary == 1 ~ "Both",
      cfwgs_bin == 0 & ref_binary == 0 ~ "Neither"
    ),   # color by relapse status
    relapse_cat = if_else(
      Num_days_to_closest_relapse <= 365,
      "Relapsed ≤365 d",
      "No relapse ≤365 d",
      missing = "No relapse ≤365 d"
    )
  )

# Normalize labels before plotting.
hyphen_rx <- "[\u2010\u2011\u2012\u2013\u2014\u2212]"  # all the usual dash culprits

plot_df2 <- plot_df2 %>%
  mutate(
    # Normalize dashes and whitespace.
    landmark_tp = str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
    
    # Factor with ASCII hyphens so facet order is deterministic.
    landmark_timepoint = factor(
      landmark_tp,
      levels = c("Post-ASCT", "Maintenance-1yr")
    ),
    
    detect_cat = factor(shape_cat,
                        levels = c("cfWGS only", "Comparator only", "Both", "Neither")),
    relapse_cat = factor(relapse_cat,
                         levels = c("Relapsed ≤365 d", "No relapse ≤365 d")),
    relapse_flag = if_else(relapse_cat == "Relapsed ≤365 d", "Relapse", "No relapse")
  )


# QC check: any rows with missing detection-pattern labels should be reviewed.
plot_df2 %>% filter(is.na(shape_cat)) %>% select(Patient, Comparator, cfwgs_bin, ref_binary)

# ------------------------------------------------------------
# 2.  Palette for detection pattern.
# ------------------------------------------------------------
detect_cols <- c(
  `cfWGS only`      = "#1b9e77",
  `Comparator only` = "#d95f02",
  `Both`            = "#7570b3",
  `Neither`         = "#999999"
)

# Viridis palette retained from the original plotting workflow.
library(viridisLite)

pal5 <- viridis(5, option = "D") 

# Assign colors to detection-pattern categories.
detect_cols <- c(
  `cfWGS only`      = pal5[1],   "#440154FF",
  `Comparator only` = pal5[2],   "#31688EFF",
  `Both`            = pal5[4],   "#73D055FF",
  `Neither`         = "#999999"   # "#FDE725FF"
)


# Correlations annotated in the blood/cfDNA-informed clinical-comparator panel.
corr_df <- plot_df2 %>%
  group_by(landmark_timepoint, Comparator) %>%
  summarize(
    rho = cor(x_plot, y_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(x_plot, y_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,
    y = 0.99
  )


# Build the historical two-comparator blood/cfDNA panel. The final manuscript
# panel adds EasyM below and is staged as Figure 4D.
p_scatter_simple_blood <- ggplot(plot_df2,
                                 aes(x = x_plot, y = y_plot,
                                     fill   = relapse_cat)) +  
  # LOD reference lines
  # 0.5166693 is the Youden's J optimal probability threshold for the blood combo
  # model (maximises sensitivity + specificity simultaneously). Derived from the
  # ROC analysis in 3_1_Optimize_cfWGS_thresholds.R and stored in selected_thr.
  geom_hline(yintercept = 0.5166693,   linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = lod_clonoMF, linetype = "dashed", colour = "grey80") +
  
  # points
  geom_point(shape = 21, size = 2, alpha = 0.9,
             colour = "black") +   # outline colour (same for all)
  
  # colour legend for detection pattern
  scale_fill_manual(name = "Detection pattern",
                    values = detect_cols) +
  
  # log axes
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0.34, 1),
    breaks = seq(0.4, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator)) +
  
  labs(title = "cfWGS of cfDNA-Derived Mutations MRD\nProbability vs. Clinical Assays",
       x = "Comparator MRD level",
       y = "Sites Model Probability") +
  # start from a white‐background theme with borders
  theme_bw(base_size = 11) +    
  
  # Color points by relapse within one year.
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black",
      "Unknown"           = "#bbbbbb"
    )
  ) +
  
  # Add Spearman correlation labels.
  geom_text(
    data = corr_df,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +
  # Final theme adjustments.
  theme(
    # draw a thin black border around each facet
    panel.border      = element_rect(colour = "black", fill = NA, size = 0.5),
    
    # kill all the internal grid lines
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    
    # if you want a white strip header with black outline:
    strip.background  = element_rect(fill = "white", colour = "black"),
    strip.text        = element_text(face = "bold"),
    
    # axis titles bold
    axis.title        = element_text(size = 11),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  ) 

p_scatter_simple_blood

# Save the historical two-comparator blood/cfDNA panel for provenance.
ggsave("Final Tables and Figures/Fig5K_cfWGS_vs_MFC_clonoSEQ_clean_Blood_muts_updated3.png",
       p_scatter_simple_blood,
       width  = 6.5, height = 5, dpi = 600)


# ===========================================================================
# SECTION 3C: SCATTER PLOTS WITH EASYМ ADDED AS THIRD COLUMN (BLOOD-DERIVED)
# ===========================================================================
# Final Figure 4D version: add EasyM as a third comparator column.

# Ensure EasyM data are available, reusing the merged table from the BM section
# when possible.
if (!exists("dat_with_easym")) {
  cat("\n=== Loading EasyM data for blood scatter plots ===\n")
  EasyM_data_file <- file.path(OUTPUT_DIR_EASYM, "EasyM_all_samples_with_optimized_calls.csv")
  
  if (file.exists(EasyM_data_file)) {
    EasyM_dat <- readr::read_csv(EasyM_data_file, show_col_types = FALSE)
    
    # Remove any existing EasyM_value columns before joining to avoid .x/.y
    # suffixes in command-line runs.
    dat_prep <- dat %>%
      select(-starts_with("EasyM_value"))
    
    dat_with_easym <- dat_prep %>%
      left_join(
        EasyM_dat %>%
          select(Patient, Timepoint, EasyM_value) %>%
          distinct(Patient, Timepoint, .keep_all = TRUE),
        by = c("Patient" = "Patient", "Timepoint" = "Timepoint"),
        relationship = "many-to-one"
      )
    
    cat(sprintf("✓ EasyM data merged: %d rows with EasyM_value\n", 
                sum(!is.na(dat_with_easym$EasyM_value))))
  } else {
    warning("EasyM data file not found - EasyM columns will be empty")
    # If EasyM data are absent, create an all-NA column so downstream code fails
    # gracefully by producing empty EasyM facets rather than stopping.
    dat_with_easym <- dat %>%
      select(-starts_with("EasyM_value")) %>%
      mutate(EasyM_value = NA_real_)
  }
}

# Blood-specific comparator LOD floor used for plotting.
lod_blood <- 1e-5

# Map EasyM thresholds by timepoint for reference lines.
easyM_threshold_lines_blood <- easyM_thresholds %>%
  filter(Timepoint %in% c("05", "07")) %>%
  mutate(
    landmark_timepoint = case_when(
      Timepoint == "05" ~ "Post-ASCT",
      Timepoint == "07" ~ "Maintenance-1yr",
      TRUE ~ NA_character_
    ),
    # Convert percentage (0-100 scale) to proportion (0-1 scale).
    xintercept = Threshold_raw_percent / 100,
    Comparator = "EasyM",
    # Match comparator factor levels to the plot data.
    Comparator = factor(Comparator, levels = c("clonoSEQ", "MFC", "EasyM"))
  ) %>%
  select(landmark_timepoint, Comparator, xintercept)

plot_df_blood_with_easym <- bind_rows(
  # Original clonoSEQ and MFC rows.
  plot_df2 %>%
    mutate(Comparator = as.character(Comparator)) %>%
    select(Patient, Sample_Code, landmark_timepoint, Comparator,
           x_plot, y_plot, relapse_cat, Num_days_to_closest_relapse),
  # EasyM comparator rows.
  dat_with_easym %>%
    mutate(Flow_pct_cells = Flow_pct_cells/100) %>%
    filter(
      Cohort == "Frontline",
      !is.na(Blood_zscore_only_sites_call),
      !is.na(landmark_tp),
      !is.na(EasyM_value)
    ) %>%
    mutate(
      Comparator = "EasyM",
      x_val = EasyM_value,
      # Cap EasyM at 100% and set a minimum floor for log-scale plotting.
      x_plot = pmin(pmax(EasyM_value, 1e-6), 1.0),
      y_plot = if_else(Blood_zscore_only_sites_prob <= 1e-5,
                       1e-5, Blood_zscore_only_sites_prob),
      landmark_tp = str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
      landmark_timepoint = factor(landmark_tp, levels = c("Post-ASCT", "Maintenance-1yr")),
      relapse_cat = if_else(
        Num_days_to_closest_relapse <= 365,
        "Relapsed ≤365 d",
        "No relapse ≤365 d",
        missing = "No relapse ≤365 d"
      ),
      relapse_cat = factor(relapse_cat, levels = c("Relapsed ≤365 d", "No relapse ≤365 d"))
    ) %>%
    select(Patient, Sample_Code, landmark_timepoint, Comparator,
           x_plot, y_plot, relapse_cat, Num_days_to_closest_relapse)
)

# Set comparator and timepoint factor order so the panel layout is stable.
plot_df_blood_with_easym <- plot_df_blood_with_easym %>%
  mutate(
    # Trim whitespace before applying explicit factor levels.
    Comparator = str_trim(as.character(Comparator)),
    Comparator = factor(Comparator, levels = c("clonoSEQ", "MFC", "EasyM"), ordered = FALSE),
    landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
  )

# Calculate Spearman correlations for all three comparators.
corr_df_blood_with_easym <- plot_df_blood_with_easym %>%
  group_by(landmark_timepoint, Comparator) %>%
  summarize(
    rho = cor(x_plot, y_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(x_plot, y_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,
    y = 0.99,
    # Match comparator factor levels to the plot data.
    Comparator = factor(Comparator, levels = c("clonoSEQ", "MFC", "EasyM"))
  )

# Re-apply factor ordering directly before plotting to avoid ordering drift from
# row binding or joins.
plot_df_blood_with_easym <- plot_df_blood_with_easym %>%
  mutate(landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr")))

corr_df_blood_with_easym <- corr_df_blood_with_easym %>%
  mutate(landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr")))

easyM_threshold_lines_blood <- easyM_threshold_lines_blood %>%
  mutate(landmark_timepoint = factor(landmark_timepoint, levels = c("Post-ASCT", "Maintenance-1yr")))

# Build the final Figure 4D panel with EasyM added.
p_scatter_with_easym_blood <- ggplot(plot_df_blood_with_easym,
                                      aes(x = x_plot, y = y_plot,
                                          fill = relapse_cat)) +
  # Reference lines for cfWGS and comparator thresholds.
  geom_hline(yintercept = 0.4215524, linetype = "dashed", colour = "grey80") +
  # Clinical-comparator LOD line for clonoSEQ and MFC facets.
  geom_vline(
    data = data.frame(
      landmark_timepoint = factor(rep(c("Post-ASCT", "Maintenance-1yr"), 2),
                                  levels = c("Post-ASCT", "Maintenance-1yr")),
      Comparator = factor(c("clonoSEQ", "clonoSEQ", "MFC", "MFC"),
                          levels = c("clonoSEQ", "MFC", "EasyM")),
      xintercept = lod_blood
    ),
    aes(xintercept = xintercept),
    linetype = "dashed", colour = "grey80"
  ) +
  # EasyM thresholds by landmark timepoint.
  geom_vline(
    data = easyM_threshold_lines_blood,
    aes(xintercept = xintercept),
    linetype = "dashed", colour = "grey80"
  ) +

  # Points
  geom_point(shape = 21, size = 2, alpha = 0.9, colour = "black") +

  # Color points by relapse within one year.
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black"
    )
  ) +

  # Log-scale x-axis with fixed labels across comparators.
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%", "100%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0.12, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +

  # Facet rows by timepoint and columns by comparator.
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator),
             as.table = TRUE) +

  # Annotations
  geom_text(
    data = corr_df_blood_with_easym,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +

  # Labels
  labs(title = "cfWGS of cfDNA-Derived Mutations MRD\nProbability vs. Clinical Assays and EasyM",
       x = "Comparator MRD level",
       y = "cfWGS Model Probability") +

  # Theme
  theme_bw(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8)
  )

p_scatter_with_easym_blood

# Save the final Figure 4D working PNG before staging it to the manuscript
# output directory.
ggsave("Final Tables and Figures/Fig5K_cfWGS_vs_MFC_clonoSEQ_EasyM_Blood_muts_updated5.png",
       p_scatter_with_easym_blood,
       width = 8.5, height = 5, dpi = 600)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 4D
#
# What this is:
#   Blood/cfDNA-informed cfWGS comparison to MFC, clonoSEQ, and EasyM.
#
# Why it is here:
#   This PNG is the final Main Figure 4D component and also supports
#   Supplementary Table 8.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig5K_cfWGS_vs_MFC_clonoSEQ_EasyM_Blood_muts_updated5.png",
  artifact_id = "FIG4D",
  role = "figure_panel_png",
  description = "Blood/cfDNA-informed cfWGS comparison to MFC, clonoSEQ, and EasyM used as Main Figure 4D.",
  script_name = "3_2_Plot_optimal_cutoff_and_clinical_concordance.R"
)

# Export source data
readr::write_csv(
  plot_df_blood_with_easym %>% mutate(Figure = "Fig5K_cfWGS_vs_clinical_assays_EasyM_blood"),
  file.path(outdir_source_data, "Fig5K_cfWGS_vs_clinical_assays_EasyM_blood_source_data.csv")
)


#### EasyM vs Blood model probability  ────────────────────────────────────────
# Build plotting data for EasyM comparison with Blood model
plot_df_easym_blood <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(Blood_zscore_only_sites_call),
    !is.na(landmark_tp),
    !is.na(EasyM_value)
  ) %>%
  mutate(
    # Set minimum detection floor for plotting
    EasyM_plot = if_else(EasyM_value <= 1e-6, 1e-6, EasyM_value),
    Blood_prob_plot = if_else(Blood_zscore_only_sites_prob <= 1e-5, 
                              1e-5, Blood_zscore_only_sites_prob),
    
    # Relapse category
    relapse_cat = if_else(
      Num_days_to_closest_relapse <= 365,
      "Relapsed ≤365 d",
      "No relapse ≤365 d",
      missing = "No relapse ≤365 d"
    ),
    
    # Factor landmark timepoint
    landmark_timepoint = factor(
      str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
      levels = c("Post-ASCT", "Maintenance-1yr")
    )
  )

# Calculate correlations
corr_easym_blood <- plot_df_easym_blood %>%
  group_by(landmark_timepoint) %>%
  summarize(
    rho = cor(EasyM_plot, Blood_prob_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(EasyM_plot, Blood_prob_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,
    y = 0.99
  )

# Create scatter plot
p_easym_blood <- ggplot(plot_df_easym_blood,
                        aes(x = EasyM_plot, y = Blood_prob_plot, fill = relapse_cat)) +
  # Reference lines
  geom_hline(yintercept = 0.380, linetype = "dashed", colour = "grey80") +
  
  # Points
  geom_point(shape = 21, size = 2, alpha = 0.9, colour = "black") +
  
  # Scales
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0.34, 1),
    breaks = seq(0.4, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_wrap(~ landmark_timepoint, nrow = 2) +
  
  # Labels and theme
  labs(title = "cfWGS of cfDNA-Derived Mutations MRD\nProbability vs. EasyM",
       x = "EasyM Proteomic MRD level",
       y = "Sites Model Probability") +
  theme_bw(base_size = 11) +
  
  # Colors
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black",
      "Unknown"           = "#bbbbbb"
    )
  ) +
  
  # Add correlation text
  geom_text(
    data = corr_easym_blood,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +
  
  # Theme adjustments
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "white", colour = "black"),
    strip.text        = element_text(face = "bold"),
    axis.title        = element_text(size = 11),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8)
  )

p_easym_blood

# Save figure to publication-quality resolution
ggsave("Final Tables and Figures/FigS_EasyM_vs_Blood_cfWGS_prob_updated6.png",
       p_easym_blood,
       width = 5, height = 6, dpi = 600)

# ══════════════════════════════════════════════════════════════════════════
# SOURCE DATA: EasyM vs blood-derived cfWGS scatter plot
# ══════════════════════════════════════════════════════════════════════════
readr::write_csv(
  plot_df_easym_blood %>% 
    select(Patient, Timepoint, EasyM_plot, Blood_prob_plot, relapse_cat) %>%
    mutate(Figure = "FigS_EasyM_vs_Blood_cfWGS_prob_updated6"),
  file.path(outdir_source_data, "FigS_EasyM_vs_Blood_cfWGS_prob_updated6_source_data.csv")
)


# ===========================================================================
# COMPLETION SUMMARY AND OUTPUT DOCUMENTATION
# ===========================================================================
# Command-line completion summary. The historical working files remain in their
# original output folders for traceability. Final manuscript-facing copies are
# staged by `ms_copy_artifact()` into final_manuscript_objects/.

cat("\n", strrep("=", 80), "\n")
cat("SCRIPT 3_2 COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 80), "\n\n")

cat("FINAL MANUSCRIPT ARTIFACTS STAGED BY THIS SCRIPT\n")
cat("----------------------------------------------------------------------------\n")
cat("  Figure 3D: BM-informed positivity by technology and cohort\n")
cat("  Figure 3E: BM-informed cfWGS probability vs MFC, clonoSEQ, and EasyM\n")
cat("  Figure 4C: blood/cfDNA-informed positivity by technology and cohort\n")
cat("  Figure 4D: blood/cfDNA-informed cfWGS probability vs MFC, clonoSEQ, and EasyM\n")
cat("  Extended Data Figure 5E-G: BM-informed clinical-comparator confusion matrices\n")
cat("  Extended Data Figure 7F-H: blood/cfDNA-informed clinical-comparator confusion matrices\n")
cat("  Supplementary Table 8: row-level model comparisons to clinical metrics\n")
cat("  Supplementary Table 10: all call metrics against clinical comparators\n\n")

cat("DATA INTEGRATION CHECKS\n")
cat("----------------------------------------------------------------------------\n")
cat(sprintf("  cfWGS rows loaded: %d\n", nrow(dat)))
cat(sprintf(
  "  EasyM rows merged: %d (%.1f%% of cfWGS rows)\n",
  sum(!is.na(dat$EasyM_value)),
  100 * sum(!is.na(dat$EasyM_value)) / nrow(dat)
))
cat(sprintf(
  "  EasyM optimized-threshold rows: %d\n",
  sum(dat$threshold_method == "optimized", na.rm = TRUE)
))
cat(sprintf(
  "  EasyM clinician-call rows: %d\n",
  sum(dat$threshold_method == "clinician", na.rm = TRUE)
))
cat("\n")

cat("EASYM THRESHOLD REFERENCE\n")
cat("----------------------------------------------------------------------------\n")
threshold_summary_tbl <- NULL
if (exists("EasyM_thresholds", inherits = FALSE) && !is.null(EasyM_thresholds)) {
  threshold_summary_tbl <- EasyM_thresholds
} else if (exists("easyM_thresholds", inherits = FALSE) && !is.null(easyM_thresholds)) {
  threshold_summary_tbl <- easyM_thresholds
}

if (!is.null(threshold_summary_tbl) && nrow(threshold_summary_tbl) > 0) {
  threshold_col <- dplyr::case_when(
    "opt_cut_raw" %in% names(threshold_summary_tbl) ~ "opt_cut_raw",
    "Threshold_raw_percent" %in% names(threshold_summary_tbl) ~ "Threshold_raw_percent",
    TRUE ~ NA_character_
  )

  for (i in seq_len(nrow(threshold_summary_tbl))) {
    row <- threshold_summary_tbl[i, ]
    tp <- row$Timepoint
    cut <- if (!is.na(threshold_col)) row[[threshold_col]] else NA_real_
    if (is.na(cut)) {
      cat(sprintf("  TP%s: threshold value column not found\n", tp))
    } else {
      cat(sprintf("  TP%s: threshold = %.4f%%\n", tp, cut))
    }
  }
} else {
  cat("  (Threshold table not available)\n")
}
cat("\n")
