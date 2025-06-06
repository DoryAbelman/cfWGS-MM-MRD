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
dat      <- readRDS(file.path(outdir, "all_patients_with_BM_and_blood_calls.rds"))
selected <- read.csv(file.path(outdir, "STable_selected_thresholds.csv"))

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
#  2.  Recoding of your landmark labels --------------------------------------
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
  filter(Cohort == "Frontline",
         !is.na(landmark_tp),
         !is.na(BM_zscore_only_call)) %>%
  # pivot technologies into long form
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_zscore_only_call),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  summarise(
    n_total  = n(),
    n_pos    = sum(Result == 1, na.rm = TRUE),
    pos_rate = n_pos / n_total,
    .groups  = "drop"
  ) %>%
  # rename for display
  mutate(
    Technology = recode(Technology,
                        "Flow_Binary"       = "MFC",
                        "Adaptive_Binary"   = "clonoSEQ",
                        "BM_zscore_only_call"   = "cfWGS_BM")
  )
# ---------------------------------------------------------------------------
#  4.  NON‑FRONTLINE cohort: pooled positivity -------------------------------
non_tbl <- dat %>%
  mutate(landmark_tp = "All time‑points") %>%
  filter(Cohort == "Non-frontline",
         !is.na(BM_zscore_only_call)) %>%
  # pivot technologies into long form
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BM_zscore_only_call),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  summarise(
    n_total  = n(),
    n_pos    = sum(Result == 1, na.rm = TRUE),
    pos_rate = n_pos / n_total,
    .groups  = "drop"
  ) %>%
  # rename for display
  mutate(
    Technology = recode(Technology,
                        "Flow_Binary"       = "MFC",
                        "Adaptive_Binary"   = "clonoSEQ",
                        "BM_zscore_only_call"   = "cfWGS_BM")) 

## Export 
readr::write_csv(
  front_tbl,
  "Positivity_by_Landmark_TimePoint_BoneMarrow_Frontline.csv"
)
readr::write_csv(
  non_tbl,
  "Positivity_All_TimePoints_BoneMarrow_NonFrontline.csv"
)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#  0.  Aesthetics ------------------------------------------------------------
custom_cols <- c("Post‑ASCT"       = "#E41A1C",   # red  (same as poster)
                 "Maintenance‑1yr" = "#377EB8",   # blue (same as poster)
                 "All time‑points" = "#999999")   # grey (non‑frontline single bar)

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

# helper for % above bar + n under bar --------------------------------------
add_bar_labels <- function() {
  list(
    geom_text(                                 # % label above bar
      aes(label = sprintf("%.0f%%", pos_rate*100)),
      vjust = -0.45, size = 3.5
    ),
    geom_text(                                 # “n =” label under bar
      aes(label = sprintf("n = %d", n_total), y = 0.01),
      vjust = 1.8, size = 3.5
    )
  )
}

# ---------------------------------------------------------------------------
#  1.  FRONTLINE plot --------------------------------------------------------
p_front <- ggplot(front_tbl,
                  aes(x = landmark_tp,
                      y = pos_rate,
                      fill = landmark_tp)) +
  geom_bar(stat = "identity", width = 0.7, colour = "black") +
  scale_fill_manual(values = custom_cols, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     limits = c(0, 1.05)) +
  labs(title = "cfWGS BM positivity: Frontline cohort",
       x = NULL, y = "Positivity rate") +
  plot_theme +
  add_bar_labels()

ggsave("Figures_May2025/MRD comparisons/BM_pos_rate_Frontline.png",
       p_front, width = 5, height = 4, dpi = 500)

# ---------------------------------------------------------------------------
#  2.  NON‑FRONTLINE plot ----------------------------------------------------
p_non <- ggplot(non_tbl,
                aes(x = landmark_tp,
                    y = pos_rate,
                    fill = landmark_tp)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black") +
  scale_fill_manual(values = custom_cols, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     limits = c(0, 1.05)) +
  labs(title = "cfWGS BM positivity — Non‑frontline cohort",
       x = NULL, y = "Positivity rate") +
  plot_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  add_bar_labels()

ggsave("Figures_May2025/MRD comparisons/BM_pos_rate_NonFrontline.png",
       p_non, width = 4, height = 4, dpi = 500)





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
front <- dat %>% filter(Cohort == "Frontline", !is.na(landmark))

# --- 3a. Pairwise concordance at Post‑ASCT ----------------------------------
pa   <- front %>% filter(landmark == "Post_ASCT")
post_conc <- bind_rows(
  pair_concord(pa, "BM_zscore_only_call", "Adaptive_Binary"),
  pair_concord(pa, "BM_zscore_only_call", "Flow_Binary"),
  pair_concord(pa, "Adaptive_Binary",  "Flow_Binary")
)

# --- 3b. Pairwise concordance at Maintenance --------------------------------
ma   <- front %>% filter(landmark == "Maintenance")
maint_conc <- bind_rows(
  pair_concord(ma, "BM_zscore_only_call", "Adaptive_Binary"),
  pair_concord(ma, "BM_zscore_only_call", "Flow_Binary")
)

# --- 3c. Positivity counts ---------------------------------------------------
pos_tbl <- front %>%
  filter(landmark %in% c("Post_ASCT","Maintenance")) %>%
  pivot_longer(cols = c(BM_zscore_only_call, Adaptive_Binary, Flow_Binary),
               names_to = "Test", values_to = "Result") %>%
  drop_na(Result) %>%
  group_by(landmark, Test) %>%
  summarise(
    pos = sum(Result == 1),
    tot = n(),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------
# 4.  FRONTLINE PPV / NPV for BM_zscore_only_call --------------------------------
ppv_npv <- function(df, pred_col = "BM_zscore_only_call") {
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

# Now you can run:
ppv_post  <- ppv_npv(pa %>% filter(!is.na(MRD_truth)))
ppv_maint <- ppv_npv(ma %>% filter(!is.na(MRD_truth)))

print(ppv_post)
print(ppv_maint)


# ---------------------------------------------------------------------------
# 5.  NON-FRONTLINE cohort ----------------------------------------------------
non <- dat %>% filter(Cohort != "Frontline")

# Build the 2×2 confusion table
non_cm <- non %>%
  filter(!is.na(BM_zscore_only_call), !is.na(MRD_truth)) %>%
  tabyl(BM_zscore_only_call, MRD_truth)

non_cm <- non_cm %>%
  mutate(BM_zscore_only_call = as.character(BM_zscore_only_call))

# Ensure both levels 0 and 1 appear as rows and columns:
if (!all(c("0","1") %in% rownames(non_cm))) {
  non_cm <- complete(non_cm, BM_zscore_only_call = c("0","1"), fill = list(`0`=0, `1`=0))
}
if (!all(c("0","1") %in% colnames(non_cm))) {
  non_cm <- cbind(non_cm, `0` = 0, `1` = 0)[, c("0","1")]
}

# Extract counts as scalars
TPn <- non_cm["1", "1", drop = TRUE]
FPn <- non_cm["1", "0", drop = TRUE]
TNn <- non_cm["0", "0", drop = TRUE]
FNn <- non_cm["0", "1", drop = TRUE]

# Now compute sensitivity and specificity
sens_non <- TPn / (TPn + FNn)
spec_non <- TNn / (TNn + FPn)

sens_non
spec_non


# Overall positivity non‑frontline
non_pos <- non %>%
  pivot_longer(cols = c(BM_zscore_only_call, Flow_Binary),
               names_to = "Test", values_to = "Result") %>%
  drop_na(Result) %>%
  group_by(Test) %>%
  summarise(pos = sum(Result==1), tot=n(), .groups="drop")

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
  X  <- g("BM_zscore_only_call","Adaptive_Binary", post_conc)
  Y  <- n("BM_zscore_only_call","Adaptive_Binary", post_conc)
  XX <- sprintf("%.0f", 100*r("BM_zscore_only_call","Adaptive_Binary", post_conc))
  Xp <- g("BM_zscore_only_call","Flow_Binary", post_conc)
  Yp <- n("BM_zscore_only_call","Flow_Binary", post_conc)
  XXp<- sprintf("%.0f", 100*r("BM_zscore_only_call","Flow_Binary", post_conc))
  Z  <- g("Adaptive_Binary","Flow_Binary", post_conc)
  W  <- n("Adaptive_Binary","Flow_Binary", post_conc)
  YY <- sprintf("%.0f",100*r("Adaptive_Binary","Flow_Binary", post_conc))
  
  # discordant counts
  n_cf_pos_cl_neg <- post_conc %>%
    filter(test_a=="BM_zscore_only_call", test_b=="Adaptive_Binary") %>%
    pull(a_pos_b_neg)
  m_cf_neg_cl_pos <- post_conc %>%
    filter(test_a=="BM_zscore_only_call", test_b=="Adaptive_Binary") %>%
    pull(a_neg_b_pos)
  
  # maintenance
  A  <- g("BM_zscore_only_call","Adaptive_Binary", maint_conc)
  B  <- n("BM_zscore_only_call","Adaptive_Binary", maint_conc)
  AA <- sprintf("%.0f",100*r("BM_zscore_only_call","Adaptive_Binary", maint_conc))
  C  <- g("BM_zscore_only_call","Flow_Binary", maint_conc)
  D  <- n("BM_zscore_only_call","Flow_Binary", maint_conc)
  BB <- sprintf("%.0f",100*r("BM_zscore_only_call","Flow_Binary", maint_conc))
  
  p <- ppv_post$PPV; q <- ppv_post$NPV
  p2<- ppv_maint$PPV; q2<- ppv_maint$NPV
  
  para <- glue("
    At post-ASCT, cfWGS agreed with clonoSEQ in {X}/{Y} ({XX}%) samples and with MFC in {Xp}/{Yp} ({XXp}%). 
    clonoSEQ vs. MFC were concordant in {Z}/{W} ({YY}%) paired samples. 
    Of the discordant post-ASCT samples, cfWGS was positive/ clonoSEQ negative in {n_cf_pos_cl_neg} cases and negative/ clonoSEQ positive in {m_cf_neg_cl_pos}. 
    At the 1-year maintenance timepoint, cfWGS agreed with clonoSEQ in {A}/{B} ({AA}%) samples and with MFC in {C}/{D} ({BB}%). 
    The PPV and NPV of cfWGS were {sprintf('%.0f',p*100)}% and {sprintf('%.0f',q*100)}% at post-ASCT, and {sprintf('%.0f',p2*100)}% and {sprintf('%.0f',q2*100)}% at maintenance. 
    In the non-frontline cohort, sensitivity and specificity of cfWGS were {sprintf('%.0f',stats_out$nonfront_sens*100)}% and {sprintf('%.0f',stats_out$nonfront_spec*100)}%, with an overall positivity rate of {stats_out$nonfront_pos %>% filter(Test=='BM_zscore_only_call') %>% summarise(sprintf('%.0f%%', 100*pos/tot)) %>% pull()}.
  ")
  
  cat(para)
}

## Export 
# 1. Export post‐ASCT pairwise concordance for bone marrow (frontline cohort)
readr::write_csv(
  post_conc,
  "Frontline_BoneMarrow_PostASCT_Pairwise_Concordance.csv"
)

# 2. Export maintenance‐timepoint pairwise concordance for bone marrow (frontline cohort)
readr::write_csv(
  maint_conc,
  "Frontline_BoneMarrow_Maintenance_Pairwise_Concordance.csv"
)

# 3. Export frontline positivity counts by test and landmark for bone marrow
#    (this includes Post_ASCT and Maintenance combined in one table)
readr::write_csv(
  pos_tbl,
  "Frontline_BoneMarrow_Positivity_PostASCT_and_Maintenance.csv"
)

# 4. Export PPV/NPV at post‐ASCT for BM_zscore_only_call (frontline bone marrow)
readr::write_csv(
  ppv_post,
  "Frontline_BoneMarrow_PostASCT_PPV_NPV.csv"
)

# 5. Export PPV/NPV at maintenance for BM_zscore_only_call (frontline bone marrow)
readr::write_csv(
  ppv_maint,
  "Frontline_BoneMarrow_Maintenance_PPV_NPV.csv"
)



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


# Columns you definitely need for inspection
id_cols   <- c("Patient", "Sample_Code", "Timepoint", "timepoint_info")

# Columns that explain why calls differ
aux_cols  <- c("Adaptive_Frequency",              # clonoSEQ cumulative VAF (rename to your actual column name)
               "Flow_pct_cells",                     # MFC % cells; rename if needed
               "BM_zscore_only_prob",                # cfWGS probability (before threshold)
               "FS", "Mean.Coverage", "detect_rate_BM", "zscore_BM", 
               "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
               "BM_MutCount_Baseline", "Blood_MutCount_Baseline")

# Make sure they exist
missing <- setdiff(aux_cols, names(dat))
if (length(missing)) warning("These columns are missing: ", paste(missing, collapse = ", "))

## Slice discordances
discord_tbl <- dat %>%
  filter(!is.na(BM_zscore_only_call),           # cfWGS call present
         !is.na(Adaptive_Binary)) %>%           # clonoSEQ call present
  mutate(category = case_when(
    BM_zscore_only_call == 1 & Adaptive_Binary == 0 ~ "cfWGS_pos / clonoSEQ_neg",
    BM_zscore_only_call == 0 & Adaptive_Binary == 1 ~ "cfWGS_neg / clonoSEQ_pos",
    TRUE                                            ~ "concordant"
  ))

disc_cf_pos <- discord_tbl %>% filter(category == "cfWGS_pos / clonoSEQ_neg")
disc_cf_neg <- discord_tbl %>% filter(category == "cfWGS_neg / clonoSEQ_pos")


outdir <- "discordant_tables"
dir.create(outdir, showWarnings = FALSE)

disc_cf_pos %>%
  select(all_of(id_cols), BM_zscore_only_call, Adaptive_Binary,
         # Flow_Binary if you want it too
         all_of(aux_cols)) %>%
  arrange(Adaptive_Frequency) %>%                        # lowest VAF at top
  write.csv(file.path(outdir, "cfWGSpos_clSEQneg.csv"), row.names = FALSE)

disc_cf_neg %>%
  select(all_of(id_cols), BM_zscore_only_call, Adaptive_Binary,
         all_of(aux_cols)) %>%
  arrange(Adaptive_Frequency) %>%
  write.csv(file.path(outdir, "cfWGSneg_clSEQpos.csv"), row.names = FALSE)


## Now for MFC 

# -----------------------------------------------------------------------------
# 2.  Build a “discordance” table for cfWGS (BM_zscore_only_call) vs MFC  -----
discord_flow_vs_bm <- dat %>%
  filter(!is.na(BM_zscore_only_call),    # cfWGS call present
         !is.na(Flow_Binary))  %>%          # MFC call present
mutate(category = case_when(
  BM_zscore_only_call == 1 & Flow_Binary == 0 ~ "cfWGS_pos / MFC_neg",
  BM_zscore_only_call == 0 & Flow_Binary == 1 ~ "cfWGS_neg / MFC_pos",
  TRUE                                          ~ "concordant"
))

disc_cf_pos_flow <- discord_flow_vs_bm %>% 
  filter(category == "cfWGS_pos / MFC_neg")

disc_cf_neg_flow <- discord_flow_vs_bm %>% 
  filter(category == "cfWGS_neg / MFC_pos")

# -----------------------------------------------------------------------------
# 3.  Export CSVs for inspection ---------------------------------------------
# Disconnect where cfWGS positive but MFC negative
disc_cf_pos_flow %>%
  select(
    all_of(id_cols),
    BM_zscore_only_call, Flow_Binary,
    all_of(aux_cols)
  ) %>%
  arrange(Flow_pct_cells) %>%   # sort by MFC percentage (lowest first)
  write.csv(
    file.path(outdir, "cfWGSpos_MFCneg.csv"),
    row.names = FALSE
  )

# Disconnect where cfWGS negative but MFC positive
disc_cf_neg_flow %>%
  select(
    all_of(id_cols),
    BM_zscore_only_call, Flow_Binary,
    all_of(aux_cols)
  ) %>%
  arrange(Flow_pct_cells) %>%   # MFC percentage (lowest first)
  write.csv(
    file.path(outdir, "cfWGSneg_MFCpos.csv"),
    row.names = FALSE
  )



### Run stats on discordant cases 
library(exact2x2)      # for exact McNemar 
# library(logistf)     # if you decide on Firth logistic

## 3.1  Build flags -----------------------------------------------------------
df_disc <- dat %>%
  filter(!is.na(BM_zscore_only_call)) %>%  # only require cfWGS to be non-NA
  mutate(
    disc_cf_vs_clono = case_when(
      !is.na(Adaptive_Binary) & BM_zscore_only_call != Adaptive_Binary ~ 1L,
      !is.na(Adaptive_Binary)                                          ~ 0L,
      TRUE                                                             ~ NA_integer_
    ),
    disc_cf_vs_flow = case_when(
      !is.na(Flow_Binary) & BM_zscore_only_call != Flow_Binary ~ 1L,
      !is.na(Flow_Binary)                                      ~ 0L,
      TRUE                                                     ~ NA_integer_
    )
  )


## 3.2  Exact McNemar for cfWGS vs clonoSEQ at post‑ASCT ----------------------
post <- df_disc %>% filter(timepoint_info == "Post_transplant")

tbl_mc <- table(
  cfWGS = post$BM_zscore_only_call,
  clono = post$Adaptive_Binary
)

# Assuming tbl_mc is your 2x2 matrix of paired results
exact2x2::exact2x2(tbl_mc, paired = TRUE)

## 3.3  Example: Wilcoxon of clonoSEQ VAF between concordant vs discordant ----
post %>%
  mutate(concord = ifelse(disc_cf_vs_clono == 1, "discordant", "concordant")) %>%
  wilcox.test(Adaptive_Frequency ~ concord, exact = TRUE)

## 3.4  Write a compact table of all discordant samples -----------------------
post %>%
  filter(disc_cf_vs_clono == 1 | disc_cf_vs_flow == 1) %>%
  select(Patient, Sample_Code, BM_zscore_only_call, Adaptive_Binary,
         Flow_Binary, Adaptive_Frequency, Flow_pct_cells,
         BM_zscore_only_prob, FS, WGS_Tumor_Fraction_Blood_plasma_cfDNA,
         Mean.Coverage, BM_MutCount_Baseline, Blood_MutCount_Baseline) %>%
  arrange(desc(Adaptive_Frequency)) %>%
  write.csv("discordant_postASCT_summary.csv", row.names = FALSE)



### Now do for maintenance 
# -----------------------------------------------------------------------------
# 3.2  Exact McNemar + Wilcoxon + export for MAINTENANCE ----------------------
maint <- df_disc %>% filter(timepoint_info == "1yr maintenance")

# 3.2a  McNemar’s test (cfWGS vs. clonoSEQ) at maintenance
tbl_mc_maint <- table(
  cfWGS = maint$BM_zscore_only_call,
  clono = maint$Adaptive_Binary
)
# Exact McNemar (install exact2x2 if not already)
# install.packages("exact2x2")
exact2x2::exact2x2(tbl_mc_maint, paired = TRUE)

# 3.2b  Wilcoxon: clonoSEQ VAF in concordant vs. discordant pairs
maint %>%
  mutate(concord = ifelse(disc_cf_vs_clono == 1, "discordant", "concordant")) %>%
  wilcox.test(Adaptive_Frequency ~ concord, exact = TRUE)

# 3.2c  Export a CSV of all discordant (cfWGS vs. clonoSEQ or cfWGS vs. MFC)
maint %>%
  filter(disc_cf_vs_clono == 1 | disc_cf_vs_flow == 1) %>%
  select(
    Patient, Sample_Code, timepoint_info,
    BM_zscore_only_call, Adaptive_Binary, Flow_Binary,
    Adaptive_Frequency, Flow_pct_cells,
    BM_zscore_only_prob, FS, WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    Mean.Coverage, BM_MutCount_Baseline, Blood_MutCount_Baseline
  ) %>%
  arrange(desc(Adaptive_Frequency)) %>%
  write.csv("discordant_maintenance_summary.csv", row.names = FALSE)








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
  filter(!is.na(MRD_truth), !is.na(BM_zscore_only_call)) %>%
  mutate(truth = factor(MRD_truth, levels = c(0,1), labels = c("neg","pos")),
         call  = factor(BM_zscore_only_call, levels = c(0,1), labels = c("neg","pos"))) %>%
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


















# --- Prepare data -------------------------------------------------------
# — Prep data ----------------------------------------------------------
#— 2.  Prepare separate data frames ---------------------------------------
dat_bm <- dat %>%
  filter(!is.na(MRD_truth), !is.na(zscore_BM)) %>%
  mutate(
    MRD_label    = factor(MRD_truth, levels = c(0,1), labels = c("Negative","Positive")),
    zscore_capped = pmin(pmax(zscore_BM, -5), 20)  # cap at [-5, 20]
  )

dat_blood <- dat %>%
  filter(!is.na(MRD_truth), !is.na(combo_blood_prob)) %>%
  mutate(
    MRD_label = factor(MRD_truth,
                       levels = c(0,1),
                       labels = c("Negative","Positive"))
  )

# --- 2.  Panel A: BM z-score density with caps --------------------------
#— Cap the threshold too -----------------------------------------------
thr_zbm_capped <- if (thr_zbm < -5) {
  -5
} else if (thr_zbm > 20) {
  20
} else thr_zbm

#— Panel A: BM z-score density with custom x-axis -----------------------
p1 <- ggplot(dat_bm, aes(x = zscore_capped, fill = MRD_label)) +
  geom_density(alpha = 0.4, adjust = 1) +
  geom_vline(xintercept = thr_zbm_capped,
             linetype = "dashed", size = 0.8) +
  scale_x_continuous(
    limits = c(-5, 20),
    breaks = c(-5, 0, 5, 10, 15, 20),
    labels = c("<-5", "0", "5", "10", "15", ">20")
  ) +
  labs(
    title = "A) BM z-score distribution",
    x     = "BM z-score",
    y     = "Density",
    fill  = "MRD truth"
  ) +
  theme_classic()


#— 4.  Panel B: Blood combo density --------------------------------------
p2 <- ggplot(dat_blood, aes(x = combo_blood_prob, fill = MRD_label)) +
  geom_density(alpha = 0.4, adjust = 0.5) +
  geom_vline(xintercept = thr_cbp, linetype = "dashed", size = 0.8) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "B) Blood combo probability",
    x     = "Blood combo probability",
    y     = "Density",
    fill  = "MRD truth"
  ) +
  theme_classic()

#— 5.  Panel C: Smoothed ROC for BM z-score -------------------------------
roc_bm_raw <- pROC::roc(
  dat_bm$MRD_label,
  dat_bm$zscore_BM,
  levels    = c("Negative", "Positive"),
  direction = "<"
)

roc_bm_s   <- pROC::smooth(roc_bm_raw, method = "binormal")
roc_bm_df  <- tibble(
  fpr = 1 - roc_bm_s$specificities,
  tpr =    roc_bm_s$sensitivities
)
auc_bm     <- as.numeric(pROC::auc(roc_bm_raw))

p3 <- ggplot(roc_bm_df, aes(x = fpr, y = tpr)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.2,
           label = paste0("AUC = ", round(auc_bm, 3))) +
  labs(
    title = "C) ROC: BM z-score",
    x     = "1 – Specificity",
    y     = "Sensitivity"
  ) +
  theme_classic()

#— 6.  Panel D: Smoothed ROC for Blood combo -----------------------------
roc_cb_raw <- pROC::roc(
  dat_blood$MRD_label,
  dat_blood$combo_blood_prob,
  levels    = c("Negative", "Positive"),
  direction = "<"
)

roc_cb_s   <- pROC::smooth(roc_cb_raw, method = "binormal")
roc_cb_df  <- tibble(
  fpr = 1 - roc_cb_s$specificities,
  tpr =    roc_cb_s$sensitivities
)
auc_cb     <- as.numeric(pROC::auc(roc_cb_raw))

p4 <- ggplot(roc_cb_df, aes(x = fpr, y = tpr)) +
  geom_line(size = 1, color = "#D55E00") +
  geom_abline(linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.2,
           label = paste0("AUC = ", round(auc_cb, 3))) +
  labs(
    title = "D) ROC: Blood combo probability",
    x     = "1 – Specificity",
    y     = "Sensitivity"
  ) +
  theme_classic()

#— 7.  Combine & save -----------------------------------------------------
fig_panel <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "cfWGS MRD calling: distributions & smoothed ROC",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  )

print(fig_panel)

ggsave(
  "FigX_cfWGS_MRD_Panel_smoothed.png",
  fig_panel,
  width  = 12,
  height = 8,
  dpi    = 500
)



### Get the contingency tables
# 1) Contingency tables ---------------------------------------------------

# a) BM z-score call
tbl_bm <- dat %>%
  filter(!is.na(MRD_truth), !is.na(call_zscore_BM)) %>%
  mutate(
    truth    = factor(MRD_truth, levels=c(0,1), labels=c("neg","pos")),
    call_bm  = factor(call_zscore_BM, levels=c(0,1), labels=c("neg","pos"))
  ) %>%
  tabyl(truth, call_bm) %>%
  adorn_totals(where="both")
print(tbl_bm)

# b) Blood combo call
tbl_bl <- dat %>%
  filter(!is.na(MRD_truth), !is.na(call_combo_blood)) %>%
  mutate(
    truth    = factor(MRD_truth, levels=c(0,1), labels=c("neg","pos")),
    call_bl  = factor(call_combo_blood, levels=c(0,1), labels=c("neg","pos"))
  ) %>%
  tabyl(truth, call_bl) %>%
  adorn_totals(where="both")
print(tbl_bl)

library(janitor)
library(gt)

make_gt_cm <- function(tabyl_tbl, title_text){
  # drop the Totals row/col for clarity in the body
  df <- tabyl_tbl %>%
    filter(truth != "Total") %>%
    select(-Total)
  
  df %>%
    # convert truth into rownames
    gt(rowname_col = "truth") %>%
    tab_header(
      title    = md(glue::glue("**{title_text}**")),
      subtitle = "Contingency table vs. MRD_truth"
    ) %>%
    fmt_number(
      columns  = vars(neg, pos),
      decimals = 0
    ) %>%
    cols_label(
      neg = "Predicted neg",
      pos = "Predicted pos"
    ) %>%
    tab_spanner(
      label   = "Prediction",
      columns = vars(neg, pos)
    ) %>%
    opt_align_table_header(align = "center")
}

# BM z-score table
gt_bm <- make_gt_cm(tbl_bm,    "BM z-score MRD call")
# Blood combo table
gt_bl <- make_gt_cm(tbl_bl,    "Blood combo MRD call")




#### Now get the negativity rate and contingency tables 

##For consistency 
dat$cfWGS_BM_Binary <- dat$call_zscore_BM