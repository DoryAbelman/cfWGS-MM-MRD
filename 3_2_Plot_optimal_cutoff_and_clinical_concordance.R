# =============================================================================
# Script:   cfWGS_MRD_make_figures.R
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
dat      <- readRDS(file.path(outdir, "all_patients_with_BM_and_blood_calls_updated2.rds"))
PATH_MODEL_LIST       <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_models_2025-06-17.rds"
PATH_THRESHOLD_LIST   <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_thresholds_2025-06-17.rds"

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
    !is.na(BM_zscore_only_detection_rate_call)
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
  file.path(outdir, "Positivity_by_Landmark_TimePoint_BoneMarrow_Frontline_updated.csv")
)

readr::write_csv(
  non_tbl,
  file.path(outdir, "Positivity_All_TimePoints_BoneMarrow_NonFrontline_updated.csv")
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_updated.png"),
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_BM_positivity_by_tech_later_line.png"),
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_updated_non_frontline.png"),
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
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y     = element_text(size = 9),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "top"
  )

# 4) Save
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_4I_BM_positivity_by_tech_facet.png"),
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
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y     = element_text(size = 9),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "top"
  )

# 4) Save
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_5I_Blood_positivity_by_tech_facet_updated.png"),
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
front <- dat %>% filter(Cohort == "Frontline", !is.na(landmark)) %>% filter(!is.na(BM_zscore_only_detection_rate_call))

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
non <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(BM_zscore_only_detection_rate_call))

non_cm <- dat %>%
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
  filter(Cohort == "Non-frontline")
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "BM_zscore_only_detection_rate_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_mfc)


## Export 
# 1. Export post-ASCT pairwise concordance (frontline BM)
readr::write_csv(
  post_conc,
  file.path(outdir, "Frontline_BoneMarrow_PostASCT_Pairwise_Concordance.csv")
)

# 2. Export maintenance-timepoint pairwise concordance (frontline BM)
readr::write_csv(
  maint_conc,
  file.path(outdir, "Frontline_BoneMarrow_Maintenance_Pairwise_Concordance.csv")
)

# 3. Export frontline positivity counts by test & landmark (Post_ASCT + Maintenance)
readr::write_csv(
  pos_tbl,
  file.path(outdir, "Frontline_BoneMarrow_Positivity_PostASCT_and_Maintenance.csv")
)

# 4. Export PPV/NPV at Post-ASCT for BM_zscore_only_detection_rate_call
readr::write_csv(
  ppv_post,
  file.path(outdir, "Frontline_BoneMarrow_PostASCT_PPV_NPV.csv")
)

# 5. Export PPV/NPV at Maintenance for BM_zscore_only_detection_rate_call
readr::write_csv(
  ppv_maint,
  file.path(outdir, "Frontline_BoneMarrow_Maintenance_PPV_NPV.csv")
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
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables.xlsx")
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
  path = file.path(outdir, "cfWGS_contingency_vs_Flow_clonoSEQ.xlsx")
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

ggsave("Final Tables and Figures/Fig4_confmat_post_ASCT_updated.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_maintenance.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig4_confmat_nonfront.png",
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
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three.png",
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
ggsave("Final Tables and Figures/Fig4J_confusion_matrices_all_three_side_by_side.png",
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
  path = file.path(outdir, "cfWGS_vs_MRD_truth_contingency_tables_blood.xlsx")
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
  path = file.path(outdir, "cfWGS_contingency_vs_Flow_clonoSEQ_blood_calls.xlsx")
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

ggsave("Final Tables and Figures/Fig5_confmat_post_ASCT_blood_updated.png",
       p_post,  width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_maintenance_blood_updated.png",
       p_maint, width = 5, height = 2.75, dpi = 600)
ggsave("Final Tables and Figures/Fig5_confmat_nonfront_blood_updated.png",
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
               "BM_zscore_only_detection_rate_call",  "BM_zscore_only_sites_call",              # cfWGS probability (before threshold)
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
  file.path(outdir_discordances, "combined_discordance_table.csv"),
  row.names = FALSE
)



### Do for non-frontline now at all timepoints
combined_discord_tbl_non_frontline <- dat %>%
  filter(
    !is.na(BM_zscore_only_detection_rate_call),
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
  file.path(outdir_discordances, "combined_discordance_table_all_timepoints.csv"),
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

# ──────────────────────────────────────────────────────────────────────────────
# 1.  Build plotting data  (front-line only, both comparators)  ───────────────
# ──────────────────────────────────────────────────────────────────────────────
lod_cfWGS   <- 0.0001   # 0.011 % 
lod_cfWGS_blood   <- 0.00061   # 0.061 %

lod_clonoMF <- 1e-5      # 10-5  (for both MFC & clonoSEQ)

shape_pal <- c(
  detected     = 16,    # filled circle
  `not detected` = 4    # open cross
)

# ──────────────────────────────────────────────────────────────────────────────
# 2.  Build the plotting data  ────────────────────────────────────────────────
plot_df <- dat %>%
  mutate(Flow_pct_cells = Flow_pct_cells/100) %>% # to be consistent
  filter(
    Cohort == "Frontline",
    !is.na(z_score_detection_rate_BM),
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
    detected  = if_else(x_val > 0 & detect_rate_BM > 0,
                        "detected", "not detected")
  ) %>%
  ungroup() %>%
  mutate(
    x_plot   = if_else(x_val <= 1e-6, 1e-6, x_val),
    y_plot   = if_else(detect_rate_BM <= 1e-5, 1e-5, detect_rate_BM),
    category = if_else(concord, "concordant", "discordant")
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3.  Draw  ───────────────────────────────────────────────────────────────────
# relapse palette
col_relapse <- c(
  `Relapsed ≤365 d`     = "red",
  `No relapse ≤365 d`   = "black"
)

# ──────────────────────────────────────────────────────────────────────────────
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
    
# ──────────────────────────────────────────────────────────────────────────────
# 3) Plot
p_scatter2 <- ggplot(plot_df2,
                     aes(x = x_plot, y = y_plot,
                         shape  = shape_cat,
                         colour = relapse_cat)) +
  
  # LOD lines
  geom_hline(yintercept = lod_cfWGS,   linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = lod_clonoMF, linetype = "dashed", colour = "grey50") +
  
  # points
  geom_point(size = 3, alpha = 0.85) +
  
  # Shape legend
  scale_shape_manual(
    name = "Detection pattern",
    values = c(
      "cfWGS only"       = 24,  # ▲
      "Comparator only"  = 25,  # ▼
      "Both detected"    = 16,  # ●
      "Neither detected" = 15   # ■
    )
  ) +
  
  # relapse colours
  scale_colour_manual(values = col_relapse, name = "Relapse") +
  
  # log–log scales with “Not detected” at half‐LOD
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%")
  ) +
  scale_y_log10(
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2),
    labels = c("Not detected", "0.01%", "0.1%", "1%")
  ) +
  
  # facets
  facet_wrap(~ Comparator, nrow = 1) +
  
  # labels & theme
  labs(
    title = "cfWGS vs clinical MRD assays (Frontline cohort)",
    x     = "Comparator MRD level",
    y     = "cfWGS cVAF"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ.png",
       p_scatter2,
       width  = 7,
       height = 3.5,
       dpi    = 600)
       

### Clean up this plot to look prettier and match other styles more
# Make sure landmark_timepoint is a factor in the order you want:
plot_df2 <- plot_df2 %>%
  mutate(
    landmark_timepoint = factor(landmark_tp,
                                levels = c("Post‑ASCT", "Maintenance‑1yr")),
    shape_cat = factor(shape_cat, 
                       levels = c("cfWGS only", "Comparator only", "Both", "Neither")),
    relapse_cat = factor(relapse_cat, levels = c("Relapsed ≤365 d", "No relapse ≤365 d"))
  )

p_scatter2 <- ggplot(plot_df2,
                     aes(x      = x_plot,
                         y      = y_plot,
                         shape  = shape_cat,
                         colour = relapse_cat)) +
  # LOD lines
  geom_hline(yintercept = lod_cfWGS,   linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = lod_clonoMF, linetype = "dashed", colour = "grey80") +
  
  # points
  geom_point(size = 3, alpha = 0.9) +
  
  # shape legend
  scale_shape_manual(
    name   = "Detection pattern",
    values = c(
      "cfWGS only"       = 24,  # ▲
      "Comparator only"  = 25,  # ▼
      "Both"    = 16,  # ●
      "Neither" = 15   # ■
    ),
    guide = guide_legend(override.aes = list(size = 4, alpha = 1))
  ) +
  
  # relapse colours
  scale_colour_manual(
    name   = "Relapse",
    values = col_relapse,
    guide  = guide_legend(order = 2, override.aes = list(shape = 16, size = 4))
  ) +
  
  # log–log axes with a little breathing room
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_log10(
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2),
    labels = c("Not detected", "0.01%", "0.1%", "1%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  # two-dimensional faceting
  facet_grid(
    rows = vars(landmark_timepoint),
    cols = vars(Comparator),
    labeller = label_wrap_gen(10)    # wrap long facet labels if needed
  ) +
  
  # titles & theme
  labs(
    title = "cfWGS vs clinical MRD assays (Training cohort)",
    x     = "Comparator MRD level",
    y     = "cfWGS cVAF"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title        = element_text(face = "bold", hjust = 0.5),
    strip.background  = element_rect(fill = "grey90", colour = NA),
    strip.text        = element_text(face = "bold"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    panel.grid.minor  = element_blank()
  )

p_scatter2 <- p_scatter2 +
  # start from a white‐background theme with borders
  theme_bw(base_size = 11) +    
  
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
    axis.title        = element_text(face = "bold", size = 12),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    
    # any other theme tweaks you like…
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  )
# save
ggsave("Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ_updated.png",
       p_scatter2,
       width  = 8.25,
       height = 4.5,
       dpi    = 600)



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

# ------------------------------------------------------------
# 3.  Build the plot (shape fixed to 21: filled circle with outline)
# ------------------------------------------------------------
p_scatter_simple <- ggplot(plot_df2,
                           aes(x = x_plot, y = y_plot,
                               fill   = detect_cat,   # interior colour = detection pattern
                               stroke = relapse_flag)) +  # outline width = relapse
  # LOD reference lines
  geom_hline(yintercept = lod_cfWGS,   linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = lod_clonoMF, linetype = "dashed", colour = "grey80") +
  
  # points
  geom_point(shape = 21, size = 3, alpha = 0.9,
             colour = "black") +   # outline colour (same for all)
  
  # colour legend for detection pattern
  scale_fill_manual(name = "Detection pattern",
                    values = detect_cols) +
  
  # fix stroke scale so legend shows outline vs. none
  scale_discrete_manual(                     # works now
    aesthetics = "stroke",
    values = c("Relapse" = 1.3, "No relapse" = 0),
    guide  = guide_legend(
      title = "Relapse",
      override.aes = list(fill   = detect_cols["Both"],
                          shape  = 21,
                          size   = 4,
                          colour = "black")
    )
  ) +
  
  # log‑log axes
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_log10(
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2),
    labels = c("Not detected", "0.01%", "0.1%", "1%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator)) +
  
  labs(title = "cfWGS vs clinical MRD assays (Training cohort)",
       x = "Comparator MRD level",
       y = "cfWGS cVAF") +
  # start from a white‐background theme with borders
  theme_bw(base_size = 11) +    
  
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
    axis.title        = element_text(size = 12),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    
    # any other theme tweaks you like…
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  )

# save
ggsave("Final Tables and Figures/Fig4K_cfWGS_vs_MFC_clonoSEQ_clean_BM_muts.png",
       p_scatter_simple,
       width  = 8.25, height = 4.5, dpi = 600)














#### Now redo for blood-derived muts
# 2.  Build the plotting data  ────────────────────────────────────────────────
plot_df <- dat %>%
  mutate(Flow_pct_cells = Flow_pct_cells/100) %>% # to be consistent
  filter(
    Cohort == "Frontline",
    !is.na(z_score_detection_rate_blood),
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
    detected  = if_else(x_val > 0 & detect_rate_BM > 0,
                        "detected", "not detected")
  ) %>%
  ungroup() %>%
  mutate(
    x_plot   = if_else(x_val <= 1e-6, 1e-6, x_val),
    y_plot   = if_else(detect_rate_blood <= 1e-5, 1e-5, detect_rate_blood),
    category = if_else(concord, "concordant", "discordant")
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3.  Draw  ───────────────────────────────────────────────────────────────────
# relapse palette
col_relapse <- c(
  `Relapsed ≤365 d`     = "red",
  `No relapse ≤365 d`   = "black"
)

# ──────────────────────────────────────────────────────────────────────────────
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

## Clean up
plot_df2 <- plot_df2 %>%
  mutate(
    landmark_timepoint = factor(landmark_tp,
                                levels = c("Post‑ASCT", "Maintenance‑1yr")),
    shape_cat = factor(shape_cat, 
                       levels = c("cfWGS only", "Comparator only", "Both", "Neither")),
    relapse_cat = factor(relapse_cat, levels = c("Relapsed ≤365 d", "No relapse ≤365 d"))
  )

p_scatter2 <- ggplot(plot_df2,
                     aes(x      = x_plot,
                         y      = y_plot,
                         shape  = shape_cat,
                         colour = relapse_cat)) +
  # LOD lines
  geom_hline(yintercept = lod_cfWGS_blood,   linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = lod_clonoMF, linetype = "dashed", colour = "grey80") +
  
  # points
  geom_point(size = 3, alpha = 0.9) +
  
  # shape legend
  scale_shape_manual(
    name   = "Detection pattern",
    values = c(
      "cfWGS only"       = 24,  # ▲
      "Comparator only"  = 25,  # ▼
      "Both"    = 16,  # ●
      "Neither" = 15   # ■
    ),
    guide = guide_legend(override.aes = list(size = 4, alpha = 1))
  ) +
  
  # relapse colours
  scale_colour_manual(
    name   = "Relapse",
    values = col_relapse,
    guide  = guide_legend(order = 2, override.aes = list(shape = 16, size = 4))
  ) +
  
  # log–log axes with a little breathing room
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%"),
    expand = expansion(mult = c(0.05, 0))
  ) +
  scale_y_log10(
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  # two-dimensional faceting
  facet_grid(
    rows = vars(landmark_timepoint),
    cols = vars(Comparator),
    labeller = label_wrap_gen(10)    # wrap long facet labels if needed
  ) +
  
  # titles & theme
  labs(
    title = "cfWGS vs clinical MRD assays (Training cohort)",
    x     = "Comparator MRD level",
    y     = "cfWGS cVAF"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title        = element_text(face = "bold", hjust = 0.5),
    strip.background  = element_rect(fill = "grey90", colour = NA),
    strip.text        = element_text(face = "bold"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    panel.grid.minor  = element_blank()
  )

p_scatter2 <- p_scatter2 +
  # start from a white‐background theme with borders
  theme_bw(base_size = 11) +    
  
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
    axis.title        = element_text(face = "bold", size = 12),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    
    # any other theme tweaks you like…
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  )
# save
ggsave("Final Tables and Figures/Fig5K_cfWGS_vs_MFC_clonoSEQ_blood_updated.png",
       p_scatter2,
       width  = 8.25,
       height = 4.5,
       dpi    = 600)




#### Clean updated version 
plot_df2 <- plot_df2 %>%                      
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
# 2.  Build the plot (shape fixed to 21: filled circle with outline)
# ------------------------------------------------------------
p_scatter_simple_blood <- ggplot(plot_df2,
                           aes(x = x_plot, y = y_plot,
                               fill   = detect_cat,   # interior colour = detection pattern
                               stroke = relapse_flag)) +  # outline width = relapse
  # LOD reference lines
#  geom_hline(yintercept = lod_cfWGS_blood,   linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = lod_clonoMF, linetype = "dashed", colour = "grey80") +
  
  # points
  geom_point(shape = 21, size = 3, alpha = 0.9,
             colour = "black") +   # outline colour (same for all)
  
  # colour legend for detection pattern
  scale_fill_manual(name = "Detection pattern",
                    values = detect_cols) +
  
  # fix stroke scale so legend shows outline vs. none
  scale_discrete_manual(                     # works now
    aesthetics = "stroke",
    values = c("Relapse" = 1.3, "No relapse" = 0),
    guide  = guide_legend(
      title = "Relapse",
      override.aes = list(fill   = detect_cols["Both"],
                          shape  = 21,
                          size   = 4,
                          colour = "black")
    )
  ) +
  
  # log‑log axes
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_log10(
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("Not detected", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(mult = c(0.5, 0.1))
  ) +
  
  facet_grid(rows = vars(landmark_timepoint),
             cols = vars(Comparator)) +
  
  labs(title = "cfWGS vs clinical MRD assays (Training cohort)",
       x = "Comparator MRD level",
       y = "cfWGS cVAF") +
  # start from a white‐background theme with borders
  theme_bw(base_size = 11) +    
  
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
    axis.title        = element_text(size = 12),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    
    # any other theme tweaks you like…
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),    # ← smaller legend labels
    
  )

# save
ggsave("Final Tables and Figures/Fig5K_cfWGS_vs_MFC_clonoSEQ_clean_Blood_muts.png",
       p_scatter_simple_blood,
       width  = 8.25, height = 4.5, dpi = 600)









### Below here is supplement

# ────────────────────────────────────────────────────────────────────────────
# 1.  Prep helper ------------------------------------------------------------
density_panel <- function(df, prob_col, thr, title_txt) {
  ggplot(df, aes(x = .data[[prob_col]],
                 fill = MRD_label)) +
    geom_density(alpha = 0.4, adjust = 0.8) +
    geom_vline(xintercept = thr, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = title_txt,
         x     = "Predicted probability",
         y     = "Density",
         fill  = "MRD truth") +
    theme_classic()
}

roc_panel <- function(df, prob_col, title_txt, col = "#0072B2") {
  roc_obj <- roc(df$MRD_label, df[[prob_col]], levels = c("Negative","Positive"))
  roc_sm  <- smooth(roc_obj, method = "binormal")
  auc_val <- as.numeric(auc(roc_obj))
  
  df_roc <- tibble(
    fpr = 1 - roc_sm$specificities,
    tpr = roc_sm$sensitivities
  )
  
  ggplot(df_roc, aes(fpr, tpr)) +
    geom_line(size = 1, colour = col) +
    geom_abline(linetype = "dashed") +
    annotate("text", x = 0.6, y = 0.2,
             label = paste0("AUC = ", round(auc_val, 3))) +
    labs(title = title_txt,
         x = "1 – Specificity", y = "Sensitivity") +
    theme_classic()
}

make_gt_cm <- function(tabyl_tbl, title_text){
  # drop the Totals row/col for clarity in the body
  df <- tabyl_tbl %>%
    filter(truth != "Total") %>%
    select(-Total)
  
  df %>%
    gt(rowname_col = "truth") %>%
    tab_header(
      title    = md(glue::glue("**{title_text}**")),
      subtitle = "Contingency table vs. MRD_truth"
    ) %>%
    fmt_number(
      columns  = c("neg", "pos"),
      decimals = 0
    ) %>%
    cols_label(
      neg = "Predicted neg",
      pos = "Predicted pos"
    ) %>%
    tab_spanner(
      label   = "Prediction",
      columns = c("neg", "pos")
    ) %>%
    opt_align_table_header(align = "center")
}

# ────────────────────────────────────────────────────────────────────────────
# 2.  Data frames ------------------------------------------------------------
dat_fig <- dat %>%
  filter(!is.na(MRD_truth)) %>%
  mutate(MRD_label = factor(MRD_truth,
                            levels = c(0,1),
                            labels = c("Negative","Positive")))

dat_bm    <- dat_fig %>% filter(!is.na(BM_zscore_only_prob))
dat_bl_s  <- dat_fig %>% filter(!is.na(BloodSensPriority_prob))
dat_bl_a  <- dat_fig %>% filter(!is.na(BloodAccPriority_prob))

# thresholds (pulled from your earlier objects)
thr_bm   <- selected_rows %>% filter(combo == "BM_zscore_only") %>% pull(threshold)
thr_bl_s <- selected_rows %>% filter(combo == "Blood_all_extras",
                                     trial == "sens>=0.9") %>%
  slice(1) %>% pull(threshold)   # sens‑priority row
thr_bl_a <- selected_rows %>% filter(combo == "Blood_all_extras",
                                     trial == "accuracy") %>%
  slice(1) %>% pull(threshold)   # acc‑priority row

# ────────────────────────────────────────────────────────────────────────────
# 3.  Panels -----------------------------------------------------------------
pA <- density_panel(dat_bm,   "BM_zscore_only_prob",
                    thr_bm,   "A) BM: sites score probability")
pB <- density_panel(dat_bl_s, "BloodSensPriority_prob",
                    thr_bl_s, "B) Blood (sens priority) probability")
pC <- roc_panel(dat_bm,   "BM_zscore_only_prob",   "C) ROC: BM sites score")
pD <- roc_panel(dat_bl_s, "BloodSensPriority_prob","D) ROC: Blood sens priority",
                col = "#D55E00")

fig_panel <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "cfWGS MRD scoring: probability distributions & ROC",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig_panel)

ggsave("FigX_cfWGS_MRD_Panel.png", fig_panel,
       width = 12, height = 8, dpi = 500)

## Show for blood seperately with senseitivity or accuracy 
# ────────────────────────────────────────────────────────────────────────────
# 3b.  Blood density + both thresholds ---------------------------------------
pB2 <- ggplot(
  dat_fig %>% filter(!is.na(BloodSensPriority_prob)),
  aes(x = BloodSensPriority_prob, fill = MRD_label)
) +
  geom_density(alpha = 0.4, adjust = 0.8) +
  # sens-priority threshold
  geom_vline(xintercept = thr_bl_s,
             linetype = "dashed", color = "#D55E00", size = 0.8) +
  # acc-priority threshold
  geom_vline(xintercept = thr_bl_a,
             linetype = "dashed", color = "#0072B2", size = 0.8) +
  coord_cartesian(xlim = c(0,1)) +
  labs(
    title = "B) Blood: predicted MRD probability",
    x = "Predicted probability", y = "Density", fill = "MRD truth"
  ) +
  theme_classic()

# ────────────────────────────────────────────────────────────────────────────
# 3d.  Combined ROC for both blood rules --------------------------------------
# compute raw ROCs
roc_s <- pROC::roc(dat_fig$MRD_label, dat_fig$BloodSensPriority_prob,
                   levels = c("Negative","Positive"), direction = "<")
roc_a <- pROC::roc(dat_fig$MRD_label, dat_fig$BloodAccPriority_prob,
                   levels = c("Negative","Positive"), direction = "<")

# smooth them
roc_s_s <- pROC::smooth(roc_s, method = "binormal")
roc_a_s <- pROC::smooth(roc_a, method = "binormal")

# build a combined data frame
df_roc <- bind_rows(
  tibble(fpr = 1 - roc_s_s$specificities,
         tpr =    roc_s_s$sensitivities,
         method = "sens>=0.9"),
  tibble(fpr = 1 - roc_a_s$specificities,
         tpr =    roc_a_s$sensitivities,
         method = "accuracy")
)

# extract AUCs
auc_s <- as.numeric(pROC::auc(roc_s))
auc_a <- as.numeric(pROC::auc(roc_a))

pD2 <- ggplot(df_roc, aes(x = fpr, y = tpr, color = method)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.2,
           label = paste0("sens>=0.9 AUC=", round(auc_s, 3)),
           color = "#D55E00") +
  annotate("text", x = 0.6, y = 0.1,
           label = paste0("accuracy AUC=", round(auc_a, 3)),
           color = "#0072B2") +
  scale_color_manual(values = c("sens>=0.9" = "#D55E00",
                                "accuracy"  = "#0072B2")) +
  labs(
    title = "D) ROC: Blood MRD call",
    x = "1 – Specificity", y = "Sensitivity", color = "Rule"
  ) +
  theme_classic()

# ────────────────────────────────────────────────────────────────────────────
# 4.  Render & save updated panel --------------------------------------------
fig_panel2 <- (pA | pB2) / (pC | pD2) +
  plot_annotation(
    title = "cfWGS MRD scoring — blood & BM combined",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig_panel2)

ggsave("FigX_cfWGS_MRD_BloodDualThresholds.png",
       fig_panel2, width = 12, height = 8, dpi = 500)


# ────────────────────────────────────────────────────────────────────────────
# 4.  Contingency tables -----------------------------------------------------
tbl_bm <- dat %>%
  filter(!is.na(MRD_truth), !is.na(BM_zscore_only_detection_rate_call)) %>%
  mutate(truth = factor(MRD_truth, levels = c(0,1), labels = c("neg","pos")),
         call  = factor(BM_zscore_only_detection_rate_call, levels = c(0,1), labels = c("neg","pos"))) %>%
  tabyl(truth, call) %>%
  adorn_totals(where = "both")

# Build the sensitivity-priority blood contingency table --------------------
tbl_bl <- dat %>%
  filter(!is.na(MRD_truth), !is.na(BloodSensPriority_call)) %>%
  mutate(truth = factor(MRD_truth, levels = c(0,1), labels = c("neg","pos")),
         call  = factor(BloodSensPriority_call, levels = c(0,1), labels = c("neg","pos"))) %>%
  tabyl(truth, call) %>%
  adorn_totals(where = "both")

# Build the accuracy-priority blood contingency table --------------------
tbl_bl_acc <- dat %>%
  filter(!is.na(MRD_truth), !is.na(BloodAccPriority_call)) %>%
  mutate(
    truth   = factor(MRD_truth, levels = c(0,1), labels = c("neg","pos")),
    call    = factor(BloodAccPriority_call, levels = c(0,1), labels = c("neg","pos"))
  ) %>%
  tabyl(truth, call) %>%
  adorn_totals(where = "both")


# 2. Generate GT tables -----------------------------------------------------
gt_bm    <- make_gt_cm(tbl_bm,    "BM z-score-only call")
gt_bl_s  <- make_gt_cm(tbl_bl,    "Blood sens-priority call")
gt_bl_a  <- make_gt_cm(tbl_bl_acc,"Blood acc-priority call")


# ────────────────────────────────────────────────────────────────────────────
# 2.  Save each table on its own ---------------------------------------------
gtsave(
  data = gt_bm,
  filename = file.path(outdir, "tbl_BM_zscore_only_contingency.png"),
  zoom = 8
)

gtsave(
  data = gt_bl_s,
  filename = file.path(outdir, "tbl_Blood_sens_priority_contingency.png"),
  zoom = 8
)

gtsave(
  data = gt_bl_a,
  filename = file.path(outdir, "tbl_Blood_acc_priority_contingency.png"),
  zoom = 8
)


#  A) CALIBRATION CURVE  -----------------------------------------------------
cal_df <- dat %>%
  filter(!is.na(MRD_truth), !is.na(BloodSensPriority_prob)) %>%
  mutate(bin = ntile(BloodSensPriority_prob, 10)) %>%            # deciles
  group_by(bin) %>%
  summarise(
    mean_pred = mean(BloodSensPriority_prob),
    obs_rate  = mean(MRD_truth)
  )

p_cal <- ggplot(cal_df, aes(mean_pred, obs_rate)) +
  geom_line() +
  geom_point(size = 2) +
  geom_abline(linetype = "dashed") +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "Calibration: BloodSensPriority",
       x = "Mean predicted probability (per decile)",
       y = "Observed MRD positivity rate") +
  theme_classic()

# ---------------------------------------------------------------------------
#  B) DECISION‑CURVE ANALYSIS  ----------------------------------------------
dca_obj <- decision_curve(
  MRD_truth ~ BloodSensPriority_prob,
  data       = dat,
  family     = binomial(link = "logit"),
  thresholds = seq(0.05, 0.95, by = 0.05),
  confidence.intervals = 0.95 # 95% CI 
)

p_dca <- plot_decision_curve(
  list(dca_obj),
  curve.names = "BloodSensPriority",
  standardize = TRUE
) +
  ggtitle("Decision‑curve analysis (standardised net benefit)")

# ---------------------------------------------------------------------------
#  C) WATER‑FALL  (baseline samples only) -----------------------------------
dat_base <- dat %>%
  filter(timepoint_info %in% c("Baseline","Diagnosis"),
         !is.na(BloodSensPriority_prob)) %>%
  arrange(desc(BloodSensPriority_prob)) %>%
  mutate(sample_id = row_number(),
         truth_f   = factor(MRD_truth, labels = c("Negative","Positive")))

p_water <- ggplot(dat_base,
                  aes(sample_id, BloodSensPriority_prob,
                      fill = truth_f)) +
  geom_col() +
  scale_fill_manual(values = c("#999999","#D55E00")) +
  labs(title = "Baseline blood MRD probabilities (ordered)",
       x = "Sample (ordered)",
       y = "Predicted probability",
       fill = "MRD truth") +
  theme_classic()

# ---------------------------------------------------------------------------
#  D) LONGITUDINAL SPAGHETTI  -----------------------------------------------
dat_long <- dat %>%
  filter(!is.na(BloodSensPriority_prob),
         !is.na(MRD_truth)) %>%
  group_by(Patient) %>%
  filter(n() >= 2) %>%            # need at least 2 time‑points
  ungroup() %>%
  mutate(MRD_label = factor(MRD_truth,
                            labels = c("Neg","Pos")))

p_long <- ggplot(dat_long,
                 aes(Date, BloodSensPriority_prob,
                     colour = MRD_label, group = Patient)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c("#0072B2","#D55E00")) +
  labs(title = "Longitudinal blood MRD probability per patient",
       x = "Collection date",
       y = "Predicted probability",
       colour = "MRD truth") +
  theme_classic()

# ---------------------------------------------------------------------------
#  SAVE ALL  -----------------------------------------------------------------
ggsave("FigX_Calibration.png",          p_cal,   width = 5, height = 5, dpi = 500)
ggsave("FigX_DecisionCurve.png",        p_dca,   width = 6, height = 4, dpi = 500)
ggsave("FigX_Waterfall.png",            p_water, width = 7, height = 4, dpi = 500)
ggsave("FigX_LongitudinalSpaghetti.png",p_long,  width = 7, height = 4, dpi = 500)



### Do same for BM 
# ---------------------------------------------------------------------------
# 1.  Pull BM threshold ------------------------------------------------------
thr_bm <- selected_rows %>%
  filter(combo == "BM_zscore_only") %>%
  pull(threshold)

# ---------------------------------------------------------------------------
# 2.  Prepare BM‐specific data frames ----------------------------------------
dat_bm_cal  <- dat %>%
  filter(!is.na(MRD_truth), !is.na(BM_zscore_only_prob)) %>%
  mutate(bin = ntile(BM_zscore_only_prob, 10),
         MRD_label = factor(MRD_truth,
                            levels = c(0,1),
                            labels = c("Negative","Positive")))

dat_bm_base <- dat %>%
  filter(timepoint_info %in% c("Baseline","Diagnosis"),
         !is.na(BM_zscore_only_prob)) %>%
  arrange(desc(BM_zscore_only_prob)) %>%
  mutate(sample_id = row_number(),
         truth_f   = factor(MRD_truth,
                            levels = c(0,1),
                            labels = c("Negative","Positive")))

dat_bm_long <- dat %>%
  filter(!is.na(BM_zscore_only_prob),
         !is.na(Date), !is.na(MRD_truth)) %>%
  group_by(Patient) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  mutate(MRD_label = factor(MRD_truth,
                            levels = c(0,1),
                            labels = c("Neg","Pos")))

# ---------------------------------------------------------------------------
# 3.  Calibration curve for BM_zscore_only ----------------------------------
cal_bm <- dat_bm_cal %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(BM_zscore_only_prob),
    obs_rate  = mean(MRD_truth)
  )

p_cal_bm <- ggplot(cal_bm, aes(mean_pred, obs_rate)) +
  geom_line() +
  geom_point(size = 2) +
  geom_abline(linetype = "dashed") +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "Calibration: BM_zscore_only",
       x = "Mean predicted probability (deciles)",
       y = "Observed MRD rate") +
  theme_classic()

# ---------------------------------------------------------------------------
# 4.  Decision-curve analysis for BM_zscore_only ----------------------------
dca_bm <- decision_curve(
  MRD_truth ~ BM_zscore_only_prob,
  data       = dat,
  family     = binomial(link = "logit"),
  thresholds = seq(0.05, 0.95, by = 0.05),
  confidence.intervals = 0L
)

p_dca_bm <- plot_decision_curve(
  list(dca_bm),
  curve.names  = "BM_zscore_only",
  standardize  = TRUE
) +
  ggtitle("Decision-curve: BM_zscore_only")

# ---------------------------------------------------------------------------
# 5.  Water-fall plot (BM baseline) -----------------------------------------
p_water_bm <- ggplot(dat_bm_base,
                     aes(sample_id, BM_zscore_only_prob,
                         fill = truth_f)) +
  geom_col() +
  scale_fill_manual(values = c("#999999","#0072B2")) +
  labs(title = "Baseline BM probabilities (BM_zscore_only)",
       x = "Sample (ordered)",
       y = "Predicted probability",
       fill = "MRD truth") +
  theme_classic()

# ---------------------------------------------------------------------------
# 6.  Longitudinal spaghetti (BM) -------------------------------------------
p_long_bm <- ggplot(dat_bm_long,
                    aes(Date, BM_zscore_only_prob,
                        colour = MRD_label, group = Patient)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c("#0072B2","#D55E00")) +
  labs(title = "Longitudinal BM probability per patient",
       x = "Collection date",
       y = "Predicted probability",
       colour = "MRD truth") +
  theme_classic()

# ---------------------------------------------------------------------------
# 7.  Save all BM plots ------------------------------------------------------
ggsave("FigX_Calibration_BM.png",          p_cal_bm,   width = 5, height = 5, dpi = 500)
ggsave("FigX_DecisionCurve_BM.png",        p_dca_bm,   width = 6, height = 4, dpi = 500)
ggsave("FigX_Waterfall_BM.png",            p_water_bm, width = 7, height = 4, dpi = 500)
ggsave("FigX_LongitudinalSpaghetti_BM.png",p_long_bm,  width = 7, height = 4, dpi = 500)



