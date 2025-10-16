# =============================================================================
# Script: 3_2_Plot_optimal_cutoff_and_clinical_concordance.R
# Project:  cfWGS MRD detection (M4 / SPORE / IMMAGINE)
# Author:   Dory Abelman
# Date:     May 28, 2025
#
# Purpose:
#   1. Read the processed dataset produced by `optimize_cfWGS_thresholds.R`
#      (contains *_prob and *_call columns plus selected threshold table).
#   2. Generate all publication‑quality visualisations:
#        • Density plots and smoothed ROC curves for BM and blood models
#        • Dual‑threshold blood panels
#        • Calibration curves
#        • Decision‑curve analysis
#        • Water‑fall and longitudinal (“spaghetti”) plots
#   3. Build and save contingency tables (gt images) for BM and blood rules.
#   4. Save all figures (PNG, 300–500 dpi) to the `Output_tables_2025/` folder.
#
# Inputs  (produced by the analysis script):
#   • Output_tables_2025/all_patients_with_BM_and_blood_calls.rds
#   • Output_tables_2025/STable_selected_thresholds.csv
#
# Outputs:
#   • FigX_cfWGS_MRD_Panel.png
#   • FigX_cfWGS_MRD_BloodDualThresholds.png
#   • FigX_Calibration*.png, FigX_DecisionCurve*.png, …
#   • tbl_*_contingency.png (gt images)
#
# =============================================================================

# -------- 0.  Load packages --------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)        # ROC + AUC
library(patchwork)   # multi‑panel figures
library(janitor)     # tabyl + adorn_totals
library(gt)          # pretty contingency tables
library(glue)
library(rmda)        # decision‑curve analysis
library(lubridate)
library(scales)   # for percent_format()


# -------- 1.  Read processed data & thresholds ------------------------------
outdir   <- "Output_tables_2025"
OUTPUT_DIR_FIGURES    <- "Output_figures_2025"
dat      <- readRDS(file.path(outdir, "all_patients_with_BM_and_blood_calls_updated5.rds"))
PATH_MODEL_LIST       <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_models_2025-09-17.rds"
PATH_THRESHOLD_LIST   <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_thresholds_2025-09-17.rds"

selected_models <- readRDS(PATH_MODEL_LIST)
selected_thr    <- readRDS(PATH_THRESHOLD_LIST)

### Now get positivity rates compared to other tech
# ---------------------------------------------------------------------------
#  1.  Helper: summarise positivity rate -------------------------------------
summarise_pos <- function(df) {
  df %>%
    summarise(
      n_total  = n(),
      n_pos    = sum(cfWGS_BM_Binary == 1, na.rm = TRUE),
      pos_rate = n_pos / n_total
    )
}

# ---------------------------------------------------------------------------
#  2.  Recoding of landmark labels --------------------------------------
# If your `timepoint_info` already uses exactly "Post‑ASCT" and "Maintenance‑1yr"
# you can drop this mutate() block.
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
#  3.  FRONTLINE cohort: positivity by landmark ------------------------------
front_tbl <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    !is.na(BM_zscore_only_detection_rate_call)
  ) %>%
  ## Add the screen column 
#  mutate(
#    BM_base_zscore_screen_call  = as.integer(BM_base_zscore_call >= 0.350),
#  ) %>%
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
      BM_zscore_only_detection_rate_call = "cfWGS",
   #   BM_base_zscore_screen_call = "cfWGS (screen)"
      
    )
  )

## Get double negatives 
# Reshape to long format for all three methods
front_long <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    landmark_tp %in% c("Post-ASCT", "Maintenance-1yr")
  ) %>%
  filter(!is.na(BM_zscore_only_detection_rate_call)) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_zscore_only_detection_rate_call),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary                     = "MFC",
      Adaptive_Binary                 = "clonoSEQ",
      BM_zscore_only_detection_rate_call = "cfWGS"
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

## Counts 
# 1) Counts per landmark timepoint
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
#  4.  NON‑FRONTLINE cohort: pooled positivity -------------------------------
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


## Export
readr::write_csv(
  front_tbl,
  file.path(outdir, "Positivity_by_Landmark_TimePoint_BoneMarrow_Frontline_updated4.csv")
)

readr::write_csv(
  non_tbl,
  file.path(outdir, "Positivity_All_TimePoints_BoneMarrow_NonFrontline_updated4.csv")
)


### Decide if to restrict for 3 or keep as just 2
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

# ---------------------------------------------------------------------------
#  1.  Frontline plot  ------------------------------------------------------------

# 1) Clean and re-factor your time-point column
front_tbl <- front_tbl %>%
  # 1) normalize that unicode hyphen to the ASCII hyphen-minus
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) %>%
  # 2) now convert to a factor with Post-ASCT first
  mutate(
    landmark_tp = factor(landmark_tp,
                         levels = c("Post-ASCT", "Maintenance-1yr"))
  )

# 2) Rebuild the grouped barplot for frontline
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
    expand = expansion(mult = c(0, 0.05)),
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
            vjust    = -0.3,
            size     = 3.5)

# 4) save
ggsave(
    filename = file.path(OUTPUT_DIR_FIGURES, "Fig_BM_positivity_by_tech_updated.png"),
  plot     = p_front_grouped,
  width    = 6,
  height   = 4,
  dpi      = 500
)


## Edit theme - figure 4I
# Updated custom palette to match the two-tone “teal → green” plus grey
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
  
  # ② manual fill so it matches your other panels
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  
  # ③ y axis as percent, 0–100, small padding
  scale_y_continuous(
    labels = percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.02)),
    limits = c(0, 100)
  ) +
  
  # ④ labels on top of bars
  geom_text(aes(label = sprintf("%d%%", round(pos_rate * 100))),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_updated2.png"),
  plot     = p_front_grouped,
  width    = 6,
  height   = 4,
  dpi      = 500
)


# ---------------------------------------------------------------------------
#  2. Now create NON‑FRONTLINE plot ----------------------------------------------------

# 1) Clean and re-factor your time-point column
non_tbl <- non_tbl %>%
  # 1) normalize that unicode hyphen to the ASCII hyphen-minus
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) 

# 2) Rebuild the grouped barplot for frontline
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
    expand = expansion(mult = c(0, 0.05)),
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
            vjust    = -0.3,
            size     = 3.5)

# 4) save
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_BM_positivity_by_tech_later_line2.png"),
  plot     = p_front_grouped,
  width    = 6,
  height   = 4,
  dpi      = 500
)


## Updated theme 
p_non_grouped <- ggplot(non_tbl, 
                          aes(x    = Technology,
                              y    = pos_rate * 100,
                              fill = landmark_tp)) +
  # ① solid bars with black outline
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  
  # ② manual fill so it matches your other panels
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  
  # ③ y axis as percent, 0–100, small padding
  scale_y_continuous(
    labels = percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.02)),
    limits = c(0, 100)
  ) +
  
  # ④ labels on top of bars
  geom_text(aes(label = sprintf("%d%%", round(pos_rate * 100))),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_updated_non_frontline5.png"),
  plot     = p_non_grouped,
  width    = 3.5,
  height   = 4,
  dpi      = 500
)


### Do a facetted figure with both 
# 1) Combine the two tables into one and add a Cohort column
front_tbl2 <- front_tbl  %>% mutate(Cohort = "Training Cohort")
non_tbl2   <- non_tbl    %>% mutate(Cohort = "Test Cohort")
combo_tbl  <- bind_rows(front_tbl2, non_tbl2)

# Set order
combo_tbl <- combo_tbl %>%
  mutate(
    Cohort      = factor(Cohort,    levels = c("Training Cohort","Test Cohort")),
    landmark_tp = factor(landmark_tp,
                         levels = c("Post-ASCT",
                                    "Maintenance-1yr",
                                    "All timepoints"))
  )

# 3) Single facetted plot
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
            size     = 3.5,
            family   = "sans") +
  
  # ③ manual fill so it matches your other panels
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  
  # ④ free‐x per facet so “clonoSEQ” drops out of Test, and percent y‐axis
  facet_wrap(~ Cohort, nrow = 1, scales = "free_x") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.02)),
    limits = c(0, 100)
  ) +
  
  # ⑤ titles
  labs(
    title = "MRD Positivity by Technology (BM-informed mutation lists)",
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

# 4) Save
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_facet5.png"),
  plot     = p_pos_by_tech,
  width    = 6.5,    # wider to accommodate two facets
  height   = 3.85,
  dpi      = 500
)








### Now redo using blood derived muts
## See number available 
# 1) Counts per landmark timepoint
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
  ## Add the screen column 
  mutate(
    Blood_zscore_screen_call  = as.integer(Blood_zscore_only_sites_prob >= 0.380),
  ) %>%
  
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call, Blood_zscore_screen_call),
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
      Blood_zscore_screen_call = "cfWGS (screen)"
    )
  )

## Get double negatives 
# Reshape to long format for all three methods
front_long <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    landmark_tp %in% c("Post-ASCT", "Maintenance-1yr")
  ) %>%
  filter(!is.na(Blood_zscore_only_sites_call)) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  mutate(
    Technology = recode(
      Technology,
      Flow_Binary                     = "MFC",
      Adaptive_Binary                 = "clonoSEQ",
      Blood_zscore_only_sites_call = "cfWGS"
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
  ## Add the screen column 
  mutate(
    Blood_zscore_screen_call  = as.integer(Blood_zscore_only_sites_prob >= 0.380),
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call, Blood_zscore_screen_call),
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
      Blood_zscore_screen_call = "cfWGS (screen)"
    )
  )


# 1) Clean and re-factor your time-point column
front_tbl <- front_tbl %>%
  # 1) normalize that unicode hyphen to the ASCII hyphen-minus
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) %>%
  # 2) now convert to a factor with Post-ASCT first
  mutate(
    landmark_tp = factor(landmark_tp,
                         levels = c("Post-ASCT", "Maintenance-1yr"))
  )

# 1) Clean and re-factor your time-point column
non_tbl <- non_tbl %>%
  # 1) normalize that unicode hyphen to the ASCII hyphen-minus
  mutate(
    landmark_tp = str_replace_all(landmark_tp, "\u2011", "-")
  ) 



# 3) Single facetted plot
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
                           "MFC"
                         ))
    
  )

## Make sentence
# --- CONFIG ---
fig_ref   <- "Figure 3D"
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

# --- RUN ---
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
            size     = 3,
            family   = "sans") +
  
  # ③ manual fill so it matches your other panels
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

# 4) Save
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_5I_Blood_positivity_by_tech_facet_updated6.png"),
  plot     = p_pos_by_tech,
  width    = 6.5,    # wider to accommodate two facets
  height   = 4,
  dpi      = 500
)






### Now get the important metrics 
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

## for non-frontline 
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

## To get the PPV and NPV 
# Now can run:
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
# 7.  (OPTIONAL)  Auto‑generate paragraph -----------------------------------
write_para <- TRUE   # set TRUE if you want it printed

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


### Get PPV and NPV seperately across technologies rather than on MRD truth
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

# 2.  Filter post-ASCT FRONTLINE samples -------------------------------
pa <- dat %>%
  filter(Cohort == "Frontline", landmark == "Post_ASCT")

# 3.  Compute PPV/NPV vs. clonoSEQ  --------------------------------------
ppv_clono <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Adaptive_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Adaptive_Binary"
)

# 4.  Compute PPV/NPV vs. MFC  ------------------------------------------
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Flow_Binary"
)

# 5.  Bind together and print -------------------------------------------
bind_rows(ppv_clono, ppv_mfc)


## Now for maintenance 
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

## Now for non-frontline 
pa <- dat %>%
  filter(Cohort == "Non-frontline") %>% filter(timepoint_info != "Diagnosis")

ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_mfc)


### Make big paragraph
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


## Export 
# 1. Export post-ASCT pairwise concordance (frontline BM)
readr::write_csv(
  post_conc,
  file.path(outdir, "Frontline_BoneMarrow_PostASCT_Pairwise_Concordance2.csv")
)

# 2. Export maintenance-timepoint pairwise concordance (frontline BM)
readr::write_csv(
  maint_conc,
  file.path(outdir, "Frontline_BoneMarrow_Maintenance_Pairwise_Concordance2.csv")
)

# 3. Export frontline positivity counts by test & landmark (Post_ASCT + Maintenance)
readr::write_csv(
  pos_tbl,
  file.path(outdir, "Frontline_BoneMarrow_Positivity_PostASCT_and_Maintenance2.csv")
)

# 4. Export PPV/NPV at Post-ASCT for BM_zscore_only_detection_rate_call
readr::write_csv(
  ppv_post,
  file.path(outdir, "Frontline_BoneMarrow_PostASCT_PPV_NPV2.csv")
)

# 5. Export PPV/NPV at Maintenance for BM_zscore_only_detection_rate_call
readr::write_csv(
  ppv_maint,
  file.path(outdir, "Frontline_BoneMarrow_Maintenance_PPV_NPV2.csv")
)




#### Now make contingency table 
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

# 4.  Export  ----------------------------------------------------
writexl::write_xlsx(
  list(
    "Contingency_Post_ASCT"     = ct_post_ASCT,
    "Contingency_Maintenance"   = ct_maint,
    "Contingency_NonFrontline"  = ct_nonfront
  ),
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables3.xlsx")
)

### Now do to MFC and clonoSEQ seperately 
# build the six tables -----------------------------------------------------
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


## Now plot the contingency tables 
## A.  Helper: turn a 1-row TP/FP/FN/TN tibble into a long tibble
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
## C.  A small plotting helper (avoids repeated code)
## ───────────────────────────────────────────────────────────────
plot_cm <- function(df, main_title){
  df <- df %>%
    mutate(
      # now x = Obs (clinical), y = Pred (cfWGS)
      Obs  = factor(Obs,  levels = c("neg","pos")),
      Pred = factor(Pred, levels = c("pos","neg"))  # note: flipped so “neg” sits in the same corner
    )
  
  ggplot(df, aes(x = Obs, y = Pred, fill = Count)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = Count), size = 4) +
    facet_wrap(~ model, nrow = 1) +
    scale_fill_viridis_c(
      option = "D", begin = 0.3, end = 0.9, guide = "none"
    ) +
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

ggsave("Final Tables and Figures/Fig4_confmat_post_ASCT_updated3.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_maintenance3.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_nonfront3.png",
       p_non,   width = 3, height = 2.75, dpi = 600)   # single facet – narrower

## As one 
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


# save
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three_3.png",
       combined_cm,
       width  = 4,
       height = 7,      # three panels tall
       dpi    = 600)

## Side by side 
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


# save
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three_side_by_side3.png",
       combined_cm,
       width  = 15,
       height = 3,      # three panels tall
       dpi    = 600)






### Now do for the blood samples 
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

# 4.  Export  ----------------------------------------------------
writexl::write_xlsx(
  list(
    "Contingency_Post_ASCT"     = ct_post_ASCT,
    "Contingency_Maintenance"   = ct_maint,
    "Contingency_NonFrontline"  = ct_nonfront
  ),
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables_blood2.xlsx")
)

### Now do to MFC and clonoSEQ seperately 
# build the six tables -----------------------------------------------------
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


# format the sentence you quoted (rounded like your example)
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

## Now plot the contingency tables 
## A.  Helper: turn a 1-row TP/FP/FN/TN tibble into a long tibble
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
## C.  A small plotting helper (avoids repeated code)
## ───────────────────────────────────────────────────────────────
plot_cm <- function(df, main_title){
  df <- df %>%
    mutate(
      # now x = Obs (clinical), y = Pred (cfWGS)
      Obs  = factor(Obs,  levels = c("neg","pos")),
      Pred = factor(Pred, levels = c("pos","neg"))  # note: flipped so “neg” sits in the same corner
    )
  
  ggplot(df, aes(x = Obs, y = Pred, fill = Count)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = Count), size = 4) +
    facet_wrap(~ model, nrow = 1) +
    scale_fill_viridis_c(
      option = "D", begin = 0.3, end = 0.9, guide = "none"
    ) +
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

ggsave("Final Tables and Figures/Fig5_confmat_post_ASCT_blood_updated4.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_maintenance_blood_updated4.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_nonfront_blood_updated4.png",
       p_non,   width = 3, height = 2.75, dpi = 600)   # single facet – narrower

## As one 
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


# save
ggsave("Final Tables and Figures/Fig5J_confusion_matrices_all_three_blood.png",
       combined_cm,
       width  = 4,
       height = 7,      # three panels tall
       dpi    = 600)



### Redo blood with the fragment combined model 

### Now do for the blood samples 
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

# 4.  Export  ----------------------------------------------------
writexl::write_xlsx(
  list(
    "Contingency_Post_ASCT"     = ct_post_ASCT,
    "Contingency_Maintenance"   = ct_maint,
    "Contingency_NonFrontline"  = ct_nonfront
  ),
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables_blood2_fragment_combined.xlsx")
)

### Now do to MFC and clonoSEQ seperately 
# build the six tables -----------------------------------------------------
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


# format the sentence you quoted (rounded like your example)
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





### Just make for everything - PPV/NPV to clinical
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

# your mapping vector
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
  dd <- df %>% select(all_of(c(pred_col, truth_col))) %>% drop_na() %>% mutate(across(everything(), as.integer))
  if (nrow(dd) == 0) return(NULL)
  
  tbl <- as.data.frame(table(Pred = dd[[pred_col]], Truth = dd[[truth_col]]))
  tbl <- complete(tbl,
                  Pred  = factor(c(0,1), levels = c(0,1)),
                  Truth = factor(c(0,1), levels = c(0,1)),
                  fill  = list(Freq = 0)) %>% arrange(Pred, Truth)
  
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

## ── Run ─────────────────────────────────────────────────────────────────
metrics_tbl <- build_metrics_frontline_vs_nonfront(dat)

metrics_tbl <- metrics_tbl %>%
  mutate(Cohort = recode(Cohort,
                         "Frontline" = "Train",
                         "Non-frontline" = "Test"))

# Example: sentences for Frontline–Maintenance
# maint_sents <- metrics_tbl %>%
#   filter(Cohort == "Frontline", Timepoint == "Maintenance") %>%
#   rowwise() %>% mutate(sentence = row_to_sentence(cur_data())) %>%
#   ungroup() %>% select(Cohort, Timepoint, Pred_Label, Comparator, sentence)

# Example: sentences for Non-frontline (aggregated)
# nf_sents <- metrics_tbl %>%
#   filter(Cohort == "Non-frontline") %>%
#   rowwise() %>% mutate(sentence = row_to_sentence(cur_data())) %>%
#   ungroup() %>% select(Cohort, Timepoint, Pred_Label, Comparator, sentence)

# Export
writexl::write_xlsx(list("All_Call_Metrics" = metrics_tbl),
                    path = file.path(outdir, "Supplementary_Table_9_All_call_metrics_against_clinical_metrics.xlsx"))


# Example: filter to the Combined Model column if it’s named "Blood_plus_fragment_call"
# metrics_tbl %>% filter(Pred_Column == "Blood_plus_fragment_call")

# Example: export the full table
# writexl::write_xlsx(list("All_Call_Metrics" = metrics_tbl),
#   path = file.path(outdir, "All_call_metrics_by_cohort_timepoint.xlsx"))

# Example: make sentences for all rows at Maintenance for Frontline
# maint_sentences <- metrics_tbl %>%
#   filter(Cohort == "Frontline", Timepoint == "Maintenance") %>%
#   rowwise() %>%
#   mutate(sentence = row_to_sentence(cur_data())) %>%
#   ungroup() %>%
#   select(Cohort, Timepoint, Pred_Label, Comparator, sentence)






# ------------------------------------------------------------------------------
# SECTION: CHARACTERISTICS OF MISCLASSIFIED SAMPLES
# ------------------------------------------------------------------------------

### Look at characteristics of the samples that were misclassified

### Add baseline mutation count 
# For each patient, pull the BM and blood mutation counts from their first “Diagnosis” or “Baseline” sample,
# then join those baseline counts back onto every row for that patient.
baseline_counts <- dat %>%
  # Keep only diagnosis/baseline visits
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  # Ensure we pick the earliest by date (in case a patient has both “Diagnosis” and “Baseline” or multiple baseline rows)
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

# Now merge those baseline counts back onto the full dataset
dat <- dat %>%
  left_join(baseline_counts, by = "Patient")

# Inspect the new columns
dat %>%
  select(Patient, timepoint_info, BM_Mutation_Count, BM_MutCount_Baseline,
         Blood_Mutation_Count, Blood_MutCount_Baseline) %>%
  head(10)






### Add in the expected positions of fragmentomics data in healthy controls 
# 1) Load the cut-offs tables
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

# 2) Turn those into named vectors for easy lookup
lower_fs <- fs_cutoffs_tbl$Lower_Cutoff
upper_fs <- fs_cutoffs_tbl$Upper_Cutoff

# pick the guassian bounds from mm_ranges
emp_ranges <- mm_ranges %>%
  select(metric, lower = lower_gaussian, upper = upper_gaussian)

lower_mm <- setNames(emp_ranges$lower, emp_ranges$metric)
upper_mm <- setNames(emp_ranges$upper, emp_ranges$metric)

# ────────────────────────────────────────────────────────────────────────────
# 3) Flag outliers in dat
# ────────────────────────────────────────────────────────────────────────────
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

# 4) Save the augmented dat
readr::write_csv(
  dat,
  file.path(outdir, "dat_with_fragment_and_DARs_outlier_flags.csv")
)


#### Now report on 
# Columns you definitely need for inspection
id_cols   <- c("Patient", "Date", "Sample_Code", "Timepoint", "timepoint_info")

# Columns that explain why calls differ
aux_cols  <- c("Adaptive_Frequency",              # clonoSEQ cumulative VAF (rename to your actual column name)
               "Flow_pct_cells", grep("^BM.*(_prob|_call)$", names(dat), value = TRUE),
               "FS", "Mean.Coverage", "detect_rate_BM", "zscore_BM", 
               "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
               "BM_MutCount_Baseline", "Blood_MutCount_Baseline", "FS_outlier", "Mean.Coverage_outlier")

# Make sure they exist
missing <- setdiff(aux_cols, names(dat))
if (length(missing)) warning("These columns are missing: ", paste(missing, collapse = ", "))



### Now do one big table with everything discordant

# 2.  Build one combined table ---------------------------------------------
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
    all_of(id_cols),     # Patient, Sample_Code, Timepoint, timepoint_info
    Cohort,              # frontline vs non-frontline
    landmark,            # e.g. "Post_ASCT" or "Maintenance"
    Comparator,          # "clonoSEQ" or "MFC"
    category,            # your discordance/concordance label
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)     # all the cfWGS, clonoSEQ & MFC metrics you specified
  ) %>%
  arrange(landmark, Patient, Comparator)

combined_discord_tbl_slim <- combined_discord_tbl %>% filter(category != "concordant")

# 3.  (Optional) write out to CSV ------------------------------------------
write.csv(
  combined_discord_tbl,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_BM2.csv"),
  row.names = FALSE
)



### Do for non-frontline now at all timepoints
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
    all_of(id_cols),     # Patient, Sample_Code, Timepoint, timepoint_info
    Cohort,              # frontline vs non-frontline
    Comparator,          # "clonoSEQ" or "MFC"
    category,            # your discordance/concordance label
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)     # all the cfWGS, clonoSEQ & MFC metrics you specified
  ) %>%
  arrange(Patient, Comparator)


## Export 
write.csv(
  combined_discord_tbl_non_frontline,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_BM_non_frontline2.csv"),
  row.names = FALSE
)





### Now check what is associated with discordance 

# 1)  Flag discordance (1 = discordant, 0 = concordant) ----------------------
tbl <- combined_discord_tbl 

# Select numeric predictors you want to test -----------------------------
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


# 3) Quick descriptive summary by direction
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

# 5) Fit pooled logistic models
glm_missed   <- fit_disc_model("is_missed")   # 1 = missed, 0 = else
glm_captured <- fit_disc_model("is_captured") # 1 = captured, 0 = else

# 6) Tidy results with odds‐ratios and 95% CIs
tidy_missed   <- tidy(glm_missed,   conf.int = TRUE, exponentiate = TRUE)
tidy_captured <- tidy(glm_captured, conf.int = TRUE, exponentiate = TRUE)

print(tidy_missed)
print(tidy_captured)

# 7) (Optional) By comparator breakdown -----------------------
by_comp <- tbl %>%
  group_by(Comparator) %>%
  nest() %>%
  mutate(
    missed_mod   = map(data, ~ glm(is_missed   ~ ., data = select(.x, all_of(num_vars), is_missed),   family=binomial)),
    captured_mod = map(data, ~ glm(is_captured ~ ., data = select(.x, all_of(num_vars), is_captured), family=binomial)),
    missed_tidy   = map(missed_mod,   ~ tidy(.x, conf.int=TRUE, exponentiate=TRUE)),
    captured_tidy = map(captured_mod, ~ tidy(.x, conf.int=TRUE, exponentiate=TRUE))
  ) %>%
  select(Comparator, missed_tidy, captured_tidy) %>%
  unnest(c(missed_tidy, captured_tidy), names_sep = "_")

print(by_comp)




### Show discordance of blood-derived muts 
### Now do one big table with everything discordant
aux_cols  <- c("Adaptive_Frequency",              # clonoSEQ cumulative VAF (rename to your actual column name)
               "Flow_pct_cells",                     # MFC % cells; rename if needed
               "FS", "Mean.Coverage",  grep("^Blood.*(_prob|_call)$", names(dat), value = TRUE),
               "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
               "BM_MutCount_Baseline", "Blood_MutCount_Baseline", "FS_outlier", "Mean.Coverage_outlier")

# 2.  Build one combined table ---------------------------------------------
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
    all_of(id_cols),     # Patient, Sample_Code, Timepoint, timepoint_info
    Cohort,              # frontline vs non-frontline
    landmark,            # e.g. "Post_ASCT" or "Maintenance"
    Comparator,          # "clonoSEQ" or "MFC"
    category,            # your discordance/concordance label
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)     # all the cfWGS, clonoSEQ & MFC metrics you specified
  ) %>%
  arrange(landmark, Patient, Comparator)

combined_discord_tbl_slim <- combined_discord_tbl2 %>% filter(category != "concordant")

# 3.  (Optional) write out to CSV ------------------------------------------
write.csv(
  combined_discord_tbl2,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_Blood2.csv"),
  row.names = FALSE
)



### Do for non-frontline now at all timepoints
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
    all_of(id_cols),     # Patient, Sample_Code, Timepoint, timepoint_info
    Cohort,              # frontline vs non-frontline
    Comparator,          # "clonoSEQ" or "MFC"
    category,            # your discordance/concordance label
    Relapsed,
    Num_days_to_closest_relapse,
    all_of(aux_cols)     # all the cfWGS, clonoSEQ & MFC metrics you specified
  ) %>%
  arrange(Patient, Comparator)


## Export 
write.csv(
  combined_discord_tbl_non_frontline2,
  file.path(outdir, "Supplementary_Table_combined_discordance_table_blood_non_frontline2.csv"),
  row.names = FALSE
)



### Create full supplementary table 
library(openxlsx)

# --- 0) Load ID map
id_map <- readRDS("id_map.rds") %>% distinct(Patient, New_ID)

Baseline_dates <- read_csv("Final Tables and Figures/Baseline dates for samples.csv")
# Make a joinable version that uses the same Patient key as your output (New_ID if available)
baseline_join <- Baseline_dates %>%
  # assume columns: patient (original ID) and start (baseline date); adapt if needed
  mutate(start = as_date(start)) %>%
  left_join(id_map, by = c("patient" = "Patient")) %>%
  mutate(Patient = coalesce(New_ID, patient)) %>%
  select(Patient, start)

# --- 1) Prep helper: rename, remap Patient, add days_from_dx, drop raw dates --
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

# --- 2) Build the four tables ----------------------------------------------
BM_Train    <- prepare_tbl(combined_discord_tbl,                id_map, baseline_join)
BM_Test     <- prepare_tbl(combined_discord_tbl_non_frontline,  id_map, baseline_join)
Blood_Train <- prepare_tbl(combined_discord_tbl2,               id_map, baseline_join)
Blood_Test  <- prepare_tbl(combined_discord_tbl_non_frontline2, id_map, baseline_join)

# --- 3) Write workbook ------------------------------------------------------
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



### Now make figure showing LOD of each

## 4L: LOD of cfWGS vs MFC/clonoSEQ discordances

# 1.  Build plotting data  (front-line only, both comparators)  ───────────────
lod_cfWGS   <- 0.00011   # 0.011 % 
lod_cfWGS_blood   <- 0.00061   # 0.061 %

lod_clonoMF <- 1e-5      # 10-5  (for both MFC & clonoSEQ)

shape_pal <- c(
  detected     = 16,    # filled circle
  `not detected` = 4    # open cross
)

# 2.  Build the plotting data  ────────────────────────────────────────────────
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
    # First, recode your comparator names so you can match them later:
    Comparator = recode(Comparator,
                        Adaptive_Frequency = "clonoSEQ",
                        Flow_pct_cells     = "MFC"),
    
    # Now build ref_binary by looking at x_val for zeros:
    ref_binary = case_when(
      # If this is a clonoSEQ row with x_val == 0, it must be “not detected”
      Comparator == "clonoSEQ" & x_val == 0 ~ 0L,
      # Otherwise use the original binary calls
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

# 2) Add shape & relapse flags
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


#### Clean version 
# ------------------------------------------------------------
# 1.  Tidy / recode (same as you already have)
# ------------------------------------------------------------
hyphen_rx <- "[\u2010\u2011\u2012\u2013\u2014\u2212]"  # all the usual dash culprits

plot_df2 <- plot_df2 %>%
  mutate(
    # normalize dashes and whitespace
    landmark_tp = str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
    
    # now factor with ASCII hyphens
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


## See what is missing 
plot_df2 %>% filter(is.na(shape_cat)) %>% select(Patient, Comparator, cfwgs_bin, ref_binary)

# ------------------------------------------------------------
# 2.  Palette for detection pattern (your choice of colours)
# ------------------------------------------------------------
detect_cols <- c(
  `cfWGS only`      = "#1b9e77",
  `Comparator only` = "#d95f02",
  `Both`            = "#7570b3",
  `Neither`         = "#999999"
)

# Updated colors
library(viridisLite)

# grab a 5‑colour Viridis ramp
pal5 <- viridis(5, option = "D") 

# pick:
#  • darkest purple    = cfWGS only
#  • deep blue         = Comparator only
#  • mid‑green (mix)   = Both
#  • bright yellow     = Neither or grey 
detect_cols <- c(
  `cfWGS only`      = pal5[1],   "#440154FF",
  `Comparator only` = pal5[2],   "#31688EFF",
  `Both`            = pal5[4],   "#73D055FF",
  `Neither`         = "#999999"   # "#FDE725FF"
)


## Get correlations
corr_df <- plot_df2 %>%
  group_by(landmark_timepoint, Comparator) %>%
  summarize(
    rho = cor(x_plot, y_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(x_plot, y_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,   # choose x/y annotation positions for your data
    y = 0.99
  )


# ------------------------------------------------------------
# 3.  Build the plot (shape fixed to 21: filled circle with outline)
# ------------------------------------------------------------
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
  
  ## add color 
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black",
      "Unknown"           = "#bbbbbb"
    )
  ) +
  
  ## add text
  geom_text(
    data = corr_df,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +
  # now tweak:
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
    
    # any other theme tweaks you like…
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  ) 

p_scatter_simple

# save
ggsave("Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ_clean_BM_muts_updated4.png",
       p_scatter_simple,
       width  = 6.5, height = 5, dpi = 600)


#### Now redo for blood-derived muts
# 2.  Build the plotting data  ────────────────────────────────────────────────
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
    # First, recode your comparator names so you can match them later:
    Comparator = recode(Comparator,
                        Adaptive_Frequency = "clonoSEQ",
                        Flow_pct_cells     = "MFC"),
    
    # Now build ref_binary by looking at x_val for zeros:
    ref_binary = case_when(
      # If this is a clonoSEQ row with x_val == 0, it must be “not detected”
      Comparator == "clonoSEQ" & x_val == 0 ~ 0L,
      # Otherwise use the original binary calls
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

# 2) Add shape & relapse flags
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


#### Clean version 
# ------------------------------------------------------------
# 1.  Tidy / recode 
# ------------------------------------------------------------
hyphen_rx <- "[\u2010\u2011\u2012\u2013\u2014\u2212]"  # all the usual dash culprits

plot_df2 <- plot_df2 %>%
  mutate(
    # normalize dashes and whitespace
    landmark_tp = str_replace_all(landmark_tp, hyphen_rx, "-") |> trimws(),
    
    # now factor with ASCII hyphens
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


## See what is missing 
plot_df2 %>% filter(is.na(shape_cat)) %>% select(Patient, Comparator, cfwgs_bin, ref_binary)

# ------------------------------------------------------------
# 2.  Palette for detection pattern (your choice of colours)
# ------------------------------------------------------------
detect_cols <- c(
  `cfWGS only`      = "#1b9e77",
  `Comparator only` = "#d95f02",
  `Both`            = "#7570b3",
  `Neither`         = "#999999"
)

# Updated colors
library(viridisLite)

# grab a 5‑colour Viridis ramp
pal5 <- viridis(5, option = "D") 

# pick:
#  • darkest purple    = cfWGS only
#  • deep blue         = Comparator only
#  • mid‑green (mix)   = Both
#  • bright yellow     = Neither or grey 
detect_cols <- c(
  `cfWGS only`      = pal5[1],   "#440154FF",
  `Comparator only` = pal5[2],   "#31688EFF",
  `Both`            = pal5[4],   "#73D055FF",
  `Neither`         = "#999999"   # "#FDE725FF"
)


## Get correlations
corr_df <- plot_df2 %>%
  group_by(landmark_timepoint, Comparator) %>%
  summarize(
    rho = cor(x_plot, y_plot, method = "spearman", use = "complete.obs"),
    p   = cor.test(x_plot, y_plot, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("ρ = %.2f\np = %.2f", rho, p),
    x = 0.035,   # choose x/y annotation positions for your data
    y = 0.99
  )


# ------------------------------------------------------------
# 3.  Build the plot (shape fixed to 21: filled circle with outline)
# ------------------------------------------------------------
p_scatter_simple_blood <- ggplot(plot_df2,
                                 aes(x = x_plot, y = y_plot,
                                     fill   = relapse_cat)) +  
  # LOD reference lines
  geom_hline(yintercept = 0.5166693,   linetype = "dashed", colour = "grey80") + # youden threshold of the used model
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
  
  ## add color 
  scale_fill_manual(
    name = "Relapse ≤1 year",
    values = c(
      "Relapsed ≤365 d"   = "red",
      "No relapse ≤365 d" = "black",
      "Unknown"           = "#bbbbbb"
    )
  ) +
  
  ## add text
  geom_text(
    data = corr_df,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 2.5,
    inherit.aes = FALSE
  ) +
  # now tweak:
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
    
    # any other theme tweaks you like…
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  ) 

p_scatter_simple_blood

# save
ggsave("Final Tables and Figures/Fig5K_cfWGS_vs_MFC_clonoSEQ_clean_Blood_muts_updated3.png",
       p_scatter_simple_blood,
       width  = 6.5, height = 5, dpi = 600)

