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


#### Rescore based on new established limit in the dilution serires 
dat <- dat %>%
  mutate(
    Blood_zscore_only_sites_call_rescored = if_else(
      Blood_zscore_only_sites_prob >= 0.457,
      1L,    # “positive” call
      0L     # “negative” call
    )
  )

### Now re-get sensetivity, specificity and accuracy 
metrics_frontline <- dat %>%
  # 1) restrict to frontline patients
  filter(Cohort == "Frontline") %>%
  # 2) compute the four cells of the confusion matrix
  summarise(
    TP = sum(Blood_zscore_only_sites_call_rescored == 1 & MRD_truth == 1, na.rm = TRUE),
    TN = sum(Blood_zscore_only_sites_call_rescored == 0 & MRD_truth == 0, na.rm = TRUE),
    FP = sum(Blood_zscore_only_sites_call_rescored == 1 & MRD_truth == 0, na.rm = TRUE),
    FN = sum(Blood_zscore_only_sites_call_rescored == 0 & MRD_truth == 1, na.rm = TRUE)
  ) %>%
  # 3) derive the performance metrics
  mutate(
    Sensitivity = TP / (TP + FN),
    Specificity = TN / (TN + FP),
    Accuracy    = (TP + TN) / (TP + TN + FP + FN)
  )

print(metrics_frontline)



### Next get the high specificity 'confirm' threshold 
# 1) Compute ROC on frontline cohort
roc_front <- dat %>%
  filter(Cohort == "Frontline") %>%
  with( roc(response  = MRD_truth,
            predictor = Blood_zscore_only_sites_prob,
            levels    = c(0,1),
            direction = "<") )

# 2) Get all thresholds with sens, spec, PPV, NPV
roc_df <- coords(
  roc_front,
  x         = "all",
  ret       = c("threshold","sensitivity","specificity","ppv","npv"),
  transpose = FALSE
) %>%
  as_tibble()

# 3) Compute prevalence in frontline cohort
prev <- dat %>%
  filter(Cohort == "Frontline") %>%
  summarise(p = mean(MRD_truth == 1, na.rm = TRUE)) %>%
  pull(p)

# 4) Add accuracy, balanced accuracy, F1
roc_metrics <- roc_df %>%
  mutate(
    accuracy      = sensitivity * prev + specificity * (1 - prev),
    bal_accuracy  = (sensitivity + specificity) / 2,
    f1            = 2 * sensitivity * ppv / (sensitivity + ppv)
  ) %>%
  arrange(desc(specificity), desc(accuracy))  # or sort by your preferred metric


# 3) Find thresholds with specificity ≥ 0.80, pick the one with highest sensitivity
best <- roc_df %>%
  filter(specificity >= 0.80) %>%
  arrange(desc(sensitivity)) %>%
  slice(1)

confirm_thr <- best$threshold
message("Chosen confirm threshold (≥80% spec): ", round(confirm_thr, 3))
message("At this cutoff: spec = ", round(best$specificity,2),
        ", sens = ", round(best$sensitivity,2))

# 4) Create confirm call column
dat <- dat %>%
  mutate(
    Blood_zscore_only_sites_call_confirm =
      if_else(Blood_zscore_only_sites_prob >= confirm_thr, 1L, 0L)
  )

# 5) Recompute metrics on frontline
metrics_confirm <- dat %>%
  filter(Cohort == "Frontline") %>%
  summarise(
    TP = sum(Blood_zscore_only_sites_call_confirm == 1 & MRD_truth == 1, na.rm = TRUE),
    TN = sum(Blood_zscore_only_sites_call_confirm == 0 & MRD_truth == 0, na.rm = TRUE),
    FP = sum(Blood_zscore_only_sites_call_confirm == 1 & MRD_truth == 0, na.rm = TRUE),
    FN = sum(Blood_zscore_only_sites_call_confirm == 0 & MRD_truth == 1, na.rm = TRUE)
  ) %>%
  mutate(
    Sensitivity = TP / (TP + FN),
    Specificity = TN / (TN + FP),
    Accuracy    = (TP + TN) / (TP + TN + FP + FN)
  )

print(metrics_confirm)



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
    !is.na(Blood_zscore_only_sites_call)
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call, Blood_zscore_only_sites_call_rescored),
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
      Blood_zscore_only_sites_call = "cfWGS_blood_confirm",
      Blood_zscore_only_sites_call_rescored = "cfWGS_blood_screen"
    )
  )

# ---------------------------------------------------------------------------
#  4.  NON‑FRONTLINE cohort: pooled positivity -------------------------------
non_tbl <- dat %>%
  mutate(landmark_tp = "All timepoints") %>%
  filter(!timepoint_info %in% c("Baseline", "Diagnosis")) %>% 
  filter(
    Cohort == "Non-frontline",
    !is.na(Blood_zscore_only_sites_call)
  ) %>%
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, Blood_zscore_only_sites_call),
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
      Blood_zscore_only_sites_call = "cfWGS_BM" 
  )
  )


## Export
readr::write_csv(
  front_tbl,
  file.path(outdir, "Positivity_by_Landmark_TimePoint_PB_cfDNA_Frontline_updated.csv")
)

readr::write_csv(
  non_tbl,
  file.path(outdir, "Positivity_All_TimePoints_PB_cfDNA_NonFrontline_updated.csv")
)


### Decide if to restrict for 3 or keep as just 2
dat %>%
  filter(Cohort == "Frontline", !is.na(landmark_tp)) %>%
  group_by(landmark_tp) %>%
  summarise(
    total            = n_distinct(Sample_Code),
    at_least_one     = sum(!is.na(Blood_zscore_only_sites_call) | !is.na(Flow_Binary) | !is.na(Adaptive_Binary)),
    all_three        = sum(!is.na(Blood_zscore_only_sites_call) & !is.na(Flow_Binary) & !is.na(Adaptive_Binary)),
    .groups = "drop"
  )



### see how many have no BM 
frontline <- dat %>%
  filter(Cohort == "Frontline")

# 1) Samples & patients with a blood‐based sites z‐score call but missing the BM call
no_bm <- frontline %>%
  filter(
    is.na(BM_zscore_only_sites_call),
    !is.na(Blood_zscore_only_sites_call),
    ! timepoint_info %in% c("Baseline", "Diagnosis", "Relapse")
  )

n_samples_no_bm  <- nrow(no_bm)
n_patients_no_bm <- no_bm %>% distinct(Patient) %>% nrow()

total <- frontline %>%
  filter(
    !is.na(Blood_zscore_only_sites_call),
    ! timepoint_info %in% c("Baseline", "Diagnosis", "Relapse")
  )

n_samples_total  <- nrow(total)
n_patients_total <- total %>% distinct(Patient) %>% nrow()

# 2) Of those, how many are MRD‐truth positive?
truth_samples_no_bm  <- no_bm %>% filter(!is.na(MRD_truth)) %>% nrow()
truth_patients_no_bm <- no_bm %>% filter(!is.na(MRD_truth)) %>% distinct(Patient) %>% nrow()

tibble(
  n_samples_total = n_samples_total, 
  n_patients_total = n_patients_total,
  samples_without_bm   = n_samples_no_bm,
  patients_without_bm  = n_patients_no_bm,
  truth_pos_samples    = truth_samples_no_bm,
  truth_pos_patients   = truth_patients_no_bm
)




## Plot 
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#  0.  Aesthetics ------------------------------------------------------------
custom_cols <- c("Post-ASCT"       = "#E41A1C",
                 "Maintenance-1yr" = "#377EB8",   
                 "All timepoints" = "#999999")   # grey (non‑frontline single bar)

plot_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    axis.title.x  = element_text(size = 12),
    axis.title.y  = element_text(size = 12),
    plot.title    = element_text(hjust = .5, face = "bold", size = 14),
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
    title = "cfWGS Positivity: Frontline cohort, PB_cfDNA-derived muts",
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
    filename = file.path(OUTPUT_DIR_FIGURES, "Fig_PB_cfDNA_positivity_by_tech_updated.png"),
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
    title = "cfWGS Positivity: Later-line cohort, PB_cfDNA-derived muts",
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_PB_cfDNA_positivity_by_tech_later_line.png"),
  plot     = p_front_grouped,
  width    = 6,
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
front <- dat %>% filter(Cohort == "Frontline", !is.na(landmark)) %>% filter(!is.na(Blood_zscore_only_sites_call))

# --- 3a. Pairwise concordance at Post‑ASCT ----------------------------------
pa   <- front %>% filter(landmark == "Post_ASCT")
post_conc <- bind_rows(
  pair_concord(pa, "Blood_zscore_only_sites_call", "Adaptive_Binary"),
  pair_concord(pa, "Blood_zscore_only_sites_call", "Flow_Binary"),
  pair_concord(pa, "Adaptive_Binary",  "Flow_Binary")
)

# --- 3b. Pairwise concordance at Maintenance --------------------------------
ma   <- front %>% filter(landmark == "Maintenance")
maint_conc <- bind_rows(
  pair_concord(ma, "Blood_zscore_only_sites_call", "Adaptive_Binary"),
  pair_concord(ma, "Blood_zscore_only_sites_call", "Flow_Binary")
)

# --- 3c. Positivity counts ---------------------------------------------------
pos_tbl <- front %>%
  filter(landmark %in% c("Post_ASCT", "Maintenance")) %>%
  pivot_longer(
    cols      = c(Blood_zscore_only_sites_call, Adaptive_Binary, Flow_Binary),
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
# 4.  FRONTLINE PPV / NPV for Blood_zscore_only_sites_call --------------------------------
ppv_npv <- function(df, pred_col = "Blood_zscore_only_sites_call") {
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
non <- dat %>% filter(Cohort == "Non-frontline") %>% filter(!is.na(Blood_zscore_only_sites_call))

non_cm <- dat %>%
  filter(Cohort != "Frontline",
         !is.na(Blood_zscore_only_sites_call),
         !is.na(MRD_truth)) %>%
  tabyl(Blood_zscore_only_sites_call, MRD_truth) %>%
  # ensure integer rows 0 and 1 exist
  complete(
    Blood_zscore_only_sites_call = c(0L, 1L), 
    fill = list(`0` = 0, `1` = 0)
  ) %>%
  column_to_rownames("Blood_zscore_only_sites_call")


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
    cols      = c(Blood_zscore_only_sites_call, Flow_Binary),
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
  X  <- g("Blood_zscore_only_sites_call","Adaptive_Binary", post_conc)
  Y  <- n("Blood_zscore_only_sites_call","Adaptive_Binary", post_conc)
  XX <- sprintf("%.0f", 100*r("Blood_zscore_only_sites_call","Adaptive_Binary", post_conc))
  Xp <- g("Blood_zscore_only_sites_call","Flow_Binary", post_conc)
  Yp <- n("Blood_zscore_only_sites_call","Flow_Binary", post_conc)
  XXp<- sprintf("%.0f", 100*r("Blood_zscore_only_sites_call","Flow_Binary", post_conc))
  Z  <- g("Adaptive_Binary","Flow_Binary", post_conc)
  W  <- n("Adaptive_Binary","Flow_Binary", post_conc)
  YY <- sprintf("%.0f",100*r("Adaptive_Binary","Flow_Binary", post_conc))
  
  # discordant counts
  n_cf_pos_cl_neg <- post_conc %>%
    filter(test_a=="Blood_zscore_only_sites_call", test_b=="Adaptive_Binary") %>%
    pull(a_pos_b_neg)
  m_cf_neg_cl_pos <- post_conc %>%
    filter(test_a=="Blood_zscore_only_sites_call", test_b=="Adaptive_Binary") %>%
    pull(a_neg_b_pos)
  
  # maintenance
  A  <- g("Blood_zscore_only_sites_call","Adaptive_Binary", maint_conc)
  B  <- n("Blood_zscore_only_sites_call","Adaptive_Binary", maint_conc)
  AA <- sprintf("%.0f",100*r("Blood_zscore_only_sites_call","Adaptive_Binary", maint_conc))
  C  <- g("Blood_zscore_only_sites_call","Flow_Binary", maint_conc)
  D  <- n("Blood_zscore_only_sites_call","Flow_Binary", maint_conc)
  BB <- sprintf("%.0f",100*r("Blood_zscore_only_sites_call","Flow_Binary", maint_conc))
  
  p <- ppv_post$PPV; q <- ppv_post$NPV
  p2<- ppv_maint$PPV; q2<- ppv_maint$NPV
  
  para <- glue("
    At post-ASCT, cfWGS agreed with clonoSEQ in {X}/{Y} ({XX}%) samples and with MFC in {Xp}/{Yp} ({XXp}%). 
    clonoSEQ vs. MFC were concordant in {Z}/{W} ({YY}%) paired samples. 
    Of the discordant post-ASCT samples, cfWGS was positive/ clonoSEQ negative in {n_cf_pos_cl_neg} cases and negative/ clonoSEQ positive in {m_cf_neg_cl_pos}. 
    At the 1-year maintenance timepoint, cfWGS agreed with clonoSEQ in {A}/{B} ({AA}%) samples and with MFC in {C}/{D} ({BB}%). 
    The PPV and NPV of cfWGS were {sprintf('%.0f',p*100)}% and {sprintf('%.0f',q*100)}% at post-ASCT, and {sprintf('%.0f',p2*100)}% and {sprintf('%.0f',q2*100)}% at maintenance. 
    In the non-frontline cohort, sensitivity and specificity of cfWGS were {sprintf('%.0f',stats_out$nonfront_sens*100)}% and {sprintf('%.0f',stats_out$nonfront_spec*100)}%, with an overall positivity rate of {stats_out$nonfront_pos %>% filter(Test=='Blood_zscore_only_sites_call') %>% summarise(sprintf('%.0f%%', 100*pos/tot)) %>% pull()}.
  ")
  
  cat(para)
}


### Get PPV and NPV seperately across technologies rather than on MRD truth
# 1.  General PPV/NPV helper that takes any truth column  -------------------
ppv_npv_any <- function(df, pred_col = "Blood_zscore_only_sites_call", truth_col) {
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
  pred_col  = "Blood_zscore_only_sites_call",
  truth_col = "Adaptive_Binary"
)

# 4.  Compute PPV/NPV vs. MFC  ------------------------------------------
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "Blood_zscore_only_sites_call",
  truth_col = "Flow_Binary"
)

# 5.  Bind together and print -------------------------------------------
bind_rows(ppv_clono, ppv_mfc)


## Now for maintenance 
pa <- dat %>%
  filter(Cohort == "Frontline", landmark == "Maintenance")
ppv_clono <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Adaptive_Binary)),
  pred_col  = "Blood_zscore_only_sites_call",
  truth_col = "Adaptive_Binary"
)
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "Blood_zscore_only_sites_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_clono, ppv_mfc)

## Now for non-frontline 
pa <- dat %>%
  filter(Cohort == "Non-frontline")
ppv_mfc <- ppv_npv_any(
  df        = pa %>% filter(!is.na(Flow_Binary)),
  pred_col  = "Blood_zscore_only_sites_call",
  truth_col = "Flow_Binary"
)

bind_rows(ppv_mfc)


## Export 
# 1. Export post-ASCT pairwise concordance (frontline PB_cfDNA)
readr::write_csv(
  post_conc,
  file.path(outdir, "Frontline_PB_cfDNA_PostASCT_Pairwise_Concordance.csv")
)

# 2. Export maintenance-timepoint pairwise concordance (frontline PB_cfDNA)
readr::write_csv(
  maint_conc,
  file.path(outdir, "Frontline_PB_cfDNA_Maintenance_Pairwise_Concordance.csv")
)

# 3. Export frontline positivity counts by test & landmark (Post_ASCT + Maintenance)
readr::write_csv(
  pos_tbl,
  file.path(outdir, "Frontline_PB_cfDNA_Positivity_PostASCT_and_Maintenance.csv")
)

# 4. Export PPV/NPV at Post-ASCT for Blood_zscore_only_sites_call
readr::write_csv(
  ppv_post,
  file.path(outdir, "Frontline_PB_cfDNA_PostASCT_PPV_NPV.csv")
)

# 5. Export PPV/NPV at Maintenance for Blood_zscore_only_sites_call
readr::write_csv(
  ppv_maint,
  file.path(outdir, "Frontline_PB_cfDNA_Maintenance_PPV_NPV.csv")
)








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
  file.path(outdir, "dat_with_fragment_and_DARs_outlier_flags_scored.csv")
)


#### Now report on 
# Columns you definitely need for inspection
id_cols   <- c("Patient", "Sample_Code", "Timepoint", "timepoint_info")

# Columns that explain why calls differ
aux_cols  <- c("Adaptive_Frequency",              # clonoSEQ cumulative VAF (rename to your actual column name)
               "Flow_pct_cells",                     # MFC % cells; rename if needed
               "Blood_zscore_only_sites_prob",  "Blood_zscore_only_sites_call", "Blood_zscore_only_sites_call_rescored",             # cfWGS probability (before threshold)
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

combined_discord_tbl_slim <- combined_discord_tbl %>% filter(category != "concordant")

# 3.  (Optional) write out to CSV ------------------------------------------
write.csv(
  combined_discord_tbl,
  file.path(outdir_discordances, "combined_discordance_table_blood_calls.csv"),
  row.names = FALSE
)



### Do for non-frontline now at all timepoints
combined_discord_tbl_non_frontline <- dat %>%
  filter(
    !is.na(Blood_zscore_only_sites_call),
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









### Now plot




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
  filter(!is.na(MRD_truth), !is.na(Blood_zscore_only_sites_call)) %>%
  mutate(truth = factor(MRD_truth, levels = c(0,1), labels = c("neg","pos")),
         call  = factor(Blood_zscore_only_sites_call, levels = c(0,1), labels = c("neg","pos"))) %>%
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



