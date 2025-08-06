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
dat      <- readRDS(file.path(outdir, "all_patients_with_BM_and_blood_calls_updated3.rds"))
PATH_MODEL_LIST       <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_models_2025-07-25.rds"
PATH_THRESHOLD_LIST   <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_thresholds_2025-07-25.rds"

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
    landmark_tp = recode(timepoint_info,
                         "Post_transplant"  = "Post‑ASCT",
                         "1yr maintenance"  = "Maintenance‑1yr",
                         .default           = NA_character_)
  )

# ---------------------------------------------------------------------------
#  3.  FRONTLINE cohort: positivity by landmark ------------------------------
front_tbl <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    !is.na(BM_base_zscore_call)
  ) %>%
  ## Add the screen column 
  mutate(
    BM_base_zscore_screen_call  = as.integer(BM_base_zscore_call >= 0.350),
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_base_zscore_call),
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
      BM_base_zscore_call = "cfWGS",
   #   BM_base_zscore_screen_call = "cfWGS (screen)"
      
    )
  )


# ---------------------------------------------------------------------------
#  4.  NON‑FRONTLINE cohort: pooled positivity -------------------------------
non_tbl <- dat %>%
  mutate(landmark_tp = "All timepoints") %>%
  filter(!timepoint_info %in% c("Baseline", "Diagnosis")) %>% 
  filter(
    Cohort == "Non-frontline",
    !is.na(BM_base_zscore_call),
    !is.na(MRD_truth) # restrict to only ones with MRD for fair comparison
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_base_zscore_call),
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
      BM_base_zscore_call = "cfWGS" 
  )
  )


## Export
readr::write_csv(
  front_tbl,
  file.path(outdir, "Positivity_by_Landmark_TimePoint_BoneMarrow_Frontline_updated2.csv")
)

readr::write_csv(
  non_tbl,
  file.path(outdir, "Positivity_All_TimePoints_BoneMarrow_NonFrontline_updated2.csv")
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_updated_non_frontline4.png"),
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
    title = "MRD Positivity by Technology",
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_facet2.png"),
  plot     = p_pos_by_tech,
  width    = 6,    # wider to accommodate two facets
  height   = 4,
  dpi      = 500
)








### Now redo using blood derived muts
front_tbl <- dat %>%
  filter(
    Cohort == "Frontline",
    !is.na(landmark_tp),
    !is.na(Blood_zscore_only_sites_call)
  ) %>%
  ## Add the screen column 
  mutate(
    Blood_zscore_screen_call  = as.integer(Blood_zscore_only_sites_prob >= 0.457),
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
    Blood_zscore_screen_call  = as.integer(Blood_zscore_only_sites_prob >= 0.457),
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
    title = "MRD Positivity by Technology",
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_5I_Blood_positivity_by_tech_facet_updated3.png"),
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
  pair_concord(pa, "BM_base_zscore_call", "Adaptive_Binary"),
  pair_concord(pa, "BM_base_zscore_call", "Flow_Binary"),
  pair_concord(pa, "Adaptive_Binary",  "Flow_Binary")
)

# --- 3b. Pairwise concordance at Maintenance --------------------------------
ma   <- front %>% filter(landmark == "Maintenance")
maint_conc <- bind_rows(
  pair_concord(ma, "BM_base_zscore_call", "Adaptive_Binary"),
  pair_concord(ma, "BM_base_zscore_call", "Flow_Binary")
)

# --- 3c. Positivity counts ---------------------------------------------------
pos_tbl <- front %>%
  filter(landmark %in% c("Post_ASCT", "Maintenance")) %>%
  pivot_longer(
    cols      = c(BM_base_zscore_call, Adaptive_Binary, Flow_Binary),
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
non <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(BM_base_zscore_call)) %>% filter(timepoint_info != "Baseline") %>% 
  filter(timepoint_info != "Diagnosis")

non_conc <- pair_concord(non, "BM_base_zscore_call", "Flow_Binary")




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
non <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(BM_base_zscore_call)) %>% filter(timepoint_info != "Diagnosis")  %>% filter(timepoint_info != "Baseline")

non_cm <- non %>%
  filter(Cohort != "Frontline",
         !is.na(BM_base_zscore_call),
         !is.na(MRD_truth)) %>%
  tabyl(BM_base_zscore_call, MRD_truth) %>%
  # ensure integer rows 0 and 1 exist
  complete(
    BM_base_zscore_call = c(0L, 1L), 
    fill = list(`0` = 0, `1` = 0)
  ) %>%
  column_to_rownames("BM_base_zscore_call")


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
    cols      = c(BM_base_zscore_call, Flow_Binary),
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
  X  <- g("BM_base_zscore_call","Adaptive_Binary", post_conc)
  Y  <- n("BM_base_zscore_call","Adaptive_Binary", post_conc)
  XX <- sprintf("%.0f", 100*r("BM_base_zscore_call","Adaptive_Binary", post_conc))
  Xp <- g("BM_base_zscore_call","Flow_Binary", post_conc)
  Yp <- n("BM_base_zscore_call","Flow_Binary", post_conc)
  XXp<- sprintf("%.0f", 100*r("BM_base_zscore_call","Flow_Binary", post_conc))
  Z  <- g("Adaptive_Binary","Flow_Binary", post_conc)
  W  <- n("Adaptive_Binary","Flow_Binary", post_conc)
  YY <- sprintf("%.0f",100*r("Adaptive_Binary","Flow_Binary", post_conc))
  
  # discordant counts
  n_cf_pos_cl_neg <- post_conc %>%
    filter(test_a=="BM_base_zscore_call", test_b=="Adaptive_Binary") %>%
    pull(a_pos_b_neg)
  m_cf_neg_cl_pos <- post_conc %>%
    filter(test_a=="BM_base_zscore_call", test_b=="Adaptive_Binary") %>%
    pull(a_neg_b_pos)
  
  # maintenance
  A  <- g("BM_base_zscore_call","Adaptive_Binary", maint_conc)
  B  <- n("BM_base_zscore_call","Adaptive_Binary", maint_conc)
  AA <- sprintf("%.0f",100*r("BM_base_zscore_call","Adaptive_Binary", maint_conc))
  C  <- g("BM_base_zscore_call","Flow_Binary", maint_conc)
  D  <- n("BM_base_zscore_call","Flow_Binary", maint_conc)
  BB <- sprintf("%.0f",100*r("BM_base_zscore_call","Flow_Binary", maint_conc))
  
  p <- ppv_post$PPV; q <- ppv_post$NPV
  p2<- ppv_maint$PPV; q2<- ppv_maint$NPV
  
  para <- glue("
    At post-ASCT, cfWGS agreed with clonoSEQ in {X}/{Y} ({XX}%) samples and with MFC in {Xp}/{Yp} ({XXp}%). 
    clonoSEQ vs. MFC were concordant in {Z}/{W} ({YY}%) paired samples. 
    Of the discordant post-ASCT samples, cfWGS was positive/ clonoSEQ negative in {n_cf_pos_cl_neg} cases and negative/ clonoSEQ positive in {m_cf_neg_cl_pos}. 
    At the 1-year maintenance timepoint, cfWGS agreed with clonoSEQ in {A}/{B} ({AA}%) samples and with MFC in {C}/{D} ({BB}%). 
    The PPV and NPV of cfWGS were {sprintf('%.0f',p*100)}% and {sprintf('%.0f',q*100)}% at post-ASCT, and {sprintf('%.0f',p2*100)}% and {sprintf('%.0f',q2*100)}% at maintenance. 
    In the non-frontline cohort, sensitivity and specificity of cfWGS were {sprintf('%.0f',stats_out$nonfront_sens*100)}% and {sprintf('%.0f',stats_out$nonfront_spec*100)}%, with an overall positivity rate of {stats_out$nonfront_pos %>% filter(Test=='BM_base_zscore_call') %>% summarise(sprintf('%.0f%%', 100*pos/tot)) %>% pull()}.
  ")
  
  cat(para)
}


### Get PPV and NPV seperately across technologies rather than on MRD truth
# 1.  General PPV/NPV helper that takes any truth column  -------------------
ppv_npv_any <- function(df, pred_col = "BM_base_zscore_call", truth_col) {
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
  pred_col  = "BM_base_zscore_call",
  truth_col = "Adaptive_Binary"
)

# 4.  Compute PPV/NPV vs. MFC  ------------------------------------------
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_base_zscore_call",
  truth_col = "Flow_Binary"
)

# 5.  Bind together and print -------------------------------------------
bind_rows(ppv_clono, ppv_mfc)


## Now for maintenance 
pa <- dat %>%
  filter(Cohort == "Frontline", landmark == "Maintenance")
ppv_clono <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Adaptive_Binary)),
  pred_col  = "BM_base_zscore_call",
  truth_col = "Adaptive_Binary"
)
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_base_zscore_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_clono, ppv_mfc)

## Now for non-frontline 
pa <- dat %>%
  filter(Cohort == "Non-frontline") %>% filter(timepoint_info != "Diagnosis")

ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_base_zscore_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_mfc)


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

# 4. Export PPV/NPV at Post-ASCT for BM_base_zscore_call
readr::write_csv(
  ppv_post,
  file.path(outdir, "Frontline_BoneMarrow_PostASCT_PPV_NPV2.csv")
)

# 5. Export PPV/NPV at Maintenance for BM_base_zscore_call
readr::write_csv(
  ppv_maint,
  file.path(outdir, "Frontline_BoneMarrow_Maintenance_PPV_NPV2.csv")
)




#### Now make contingency table 
# helper: build a tidy 2 × 2 contingency table ----------------------------
make_ct <- function(df,
                    pred  = "BM_base_zscore_call",   # cfWGS
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
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables2.xlsx")
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
  path = file.path(outdir, "cfWGS_contingency_vs_Flow_clonoSEQ2.xlsx")
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

ggsave("Final Tables and Figures/Fig4_confmat_post_ASCT_updated2.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_maintenance2.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_nonfront2.png",
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
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three_2.png",
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
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three_side_by_side2.png",
       combined_cm,
       width  = 15,
       height = 3,      # three panels tall
       dpi    = 600)






### Now do for the blood samples 
front_blood <- dat %>% filter(Cohort == "Frontline", !is.na(landmark)) %>% filter(!is.na(Blood_zscore_only_sites_call))
non_blood <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(Blood_zscore_only_sites_call))

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
  
  tibble(
    TP = tp, FP = fp, FN = fn, TN = tn,
    Sensitivity = TP/(TP+FN),
    Specificity = TN/(TN+FP),
    PPV = TP/(TP+FP),
    NPV = TN/(TN+FN)
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
  path = file.path(outdir, "cfWGS_contingency_vs_Flow_clonoSEQ_blood_calls_2.xlsx")
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

ggsave("Final Tables and Figures/Fig5_confmat_post_ASCT_blood_updated2.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_maintenance_blood_updated2.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_nonfront_blood_updated2.png",
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
id_cols   <- c("Patient", "Sample_Code", "Timepoint", "timepoint_info")

# Columns that explain why calls differ
aux_cols  <- c("Adaptive_Frequency",              # clonoSEQ cumulative VAF (rename to your actual column name)
               "Flow_pct_cells",                     # MFC % cells; rename if needed
               "BM_base_zscore_call",  "BM_zscore_only_sites_call",  "BM_base_zscore_call", "BM_base_zscore_prob",           # cfWGS probability (before threshold)
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
    !is.na(BM_base_zscore_call),
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
      BM_base_zscore_call == 1L & Reference == 0L ~
        paste0("cfWGS_pos / ", Comparator, "_neg"),
      BM_base_zscore_call == 0L & Reference == 1L ~
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
  file.path(outdir, "Supplementary_Table_combined_discordance_table_BM.csv"),
  row.names = FALSE
)



### Do for non-frontline now at all timepoints
combined_discord_tbl_non_frontline <- dat %>%
  filter(
    !is.na(BM_base_zscore_call),
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
      BM_base_zscore_call == 1L & Reference == 0L ~
        paste0("cfWGS_pos / ", Comparator, "_neg"),
      BM_base_zscore_call == 0L & Reference == 1L ~
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
  file.path(outdir, "Supplementary_Table_combined_discordance_table_BM_non_frontline.csv"),
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
    !is.na(BM_base_zscore_prob),
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
    
    cfwgs_bin = BM_base_zscore_call,
    concord   = (cfwgs_bin == ref_binary),
    detected  = if_else(BM_base_zscore_call == 1,
                        "detected", "not detected")
  ) %>%
  ungroup() %>%
  mutate(
    x_plot   = if_else(x_val <= 1e-6, 1e-6, x_val),
    y_plot   = if_else(BM_base_zscore_prob <= 1e-5, 1e-5, BM_base_zscore_prob),
    category = if_else(concord, "concordant", "discordant")
  ) %>% 
  select(
    Patient, Sample_Code, landmark_tp, Comparator, x_val, x_plot,
    BM_base_zscore_prob, y_plot,
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
plot_df2 <- plot_df2 %>%                       # << only this line changes
  mutate(
    landmark_timepoint = factor(landmark_tp,
                                levels = c("Post‑ASCT", "Maintenance‑1yr")),
    detect_cat = factor(shape_cat,
                        levels = c("cfWGS only", "Comparator only",
                                   "Both", "Neither")),
    relapse_cat = factor(relapse_cat,
                         levels = c("Relapsed ≤365 d",
                                    "No relapse ≤365 d")),
    relapse_flag = if_else(relapse_cat == "Relapsed ≤365 d",
                           "Relapse", "No relapse")   
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
    x = 0.035,   # choose your x/y annotation positions for your data
    y = 0.99
  )


# ------------------------------------------------------------
# 3.  Build the plot (shape fixed to 21: filled circle with outline)
# ------------------------------------------------------------
p_scatter_simple <- ggplot(plot_df2,
                           aes(x = x_plot, y = y_plot,
                               fill   = relapse_cat)) +  
  # LOD reference lines
  geom_hline(yintercept = 0.410,   linetype = "dashed", colour = "grey80") +
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
    limits = c(0.18, 1),
    breaks = seq(0.2, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator)) +
  
  labs(title = "Comparison of cfWGS MRD Probability\nand Clinical Assays Across Landmark Timepoints",
       x = "Comparator MRD level",
       y = "Combined Model Probability") +
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

# save
ggsave("Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ_clean_BM_muts_updated.png",
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
# 1.  Tidy / recode (same as you already have)
# ------------------------------------------------------------
plot_df2 <- plot_df2 %>%                       # << only this line changes
  mutate(
    landmark_timepoint = factor(landmark_tp,
                                levels = c("Post‑ASCT", "Maintenance‑1yr")),
    detect_cat = factor(shape_cat,
                        levels = c("cfWGS only", "Comparator only",
                                   "Both", "Neither")),
    relapse_cat = factor(relapse_cat,
                         levels = c("Relapsed ≤365 d",
                                    "No relapse ≤365 d")),
    relapse_flag = if_else(relapse_cat == "Relapsed ≤365 d",
                           "Relapse", "No relapse")   
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
    x = 0.0019,   # choose your x/y annotation positions for your data
    y = 0.99
  )


# ------------------------------------------------------------
# 3.  Build the plot (shape fixed to 21: filled circle with outline)
# ------------------------------------------------------------
p_scatter_simple_blood <- ggplot(plot_df2,
                                 aes(x = x_plot, y = y_plot,
                                     fill   = relapse_cat)) +  
  # LOD reference lines
  geom_hline(yintercept = 0.5226406,   linetype = "dashed", colour = "grey80") + # youden threshold
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
    limits = c(0.38, 1),
    breaks = seq(0.4, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator)) +
  
  labs(title = "Comparison of cfWGS MRD Probability\nand Clinical Assays Across Landmark Timepoints",
       x = "Comparator MRD level",
       y = "Combined Model Probability") +
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
ggsave("Final Tables and Figures/Fig5K_cfWGS_vs_MFC_clonoSEQ_clean_Bloof_muts_updated.png",
       p_scatter_simple_blood,
       width  = 6.5, height = 5, dpi = 600)
