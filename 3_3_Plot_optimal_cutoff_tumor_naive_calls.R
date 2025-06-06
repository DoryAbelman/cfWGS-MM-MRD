# =============================================================================
# Script:   cfWGS_MRD_make_figures_cfDNA_calls.R
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
      n_total  = dplyr::n(),
      n_pos    = sum(BloodSensPriority_call == 1, na.rm = TRUE),
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
         !is.na(BloodSensPriority_call)) %>%
  # pivot technologies into long form
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BloodSensPriority_call),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  summarise(
    n_total  = dplyr::n(),                # force dplyr’s dplyr::n()
    n_pos    = sum(Result == 1, na.rm = TRUE),
    pos_rate = n_pos / dplyr::n(),        # or refer back to n_total
    .groups  = "drop"
  ) %>%
  # rename for display
  mutate(
    Technology = recode(Technology,
                        "Flow_Binary"       = "MFC",
                        "Adaptive_Binary"   = "clonoSEQ",
                        "BloodSensPriority_call"   = "cfWGS_PB_cfDNA")
  )
# ---------------------------------------------------------------------------
#  4.  NON‑FRONTLINE cohort: pooled positivity -------------------------------
non_tbl <- dat %>%
  mutate(landmark_tp = "All time‑points") %>%
  filter(Cohort == "Non-frontline",
         !is.na(BloodSensPriority_call)) %>%
  # pivot technologies into long form
  pivot_longer(
    cols      = c(Flow_Binary, Adaptive_Binary, BloodSensPriority_call),
    names_to  = "Technology",
    values_to = "Result"
  ) %>%
  filter(!is.na(Result)) %>%
  group_by(landmark_tp, Technology) %>%
  summarise(
    n_total  = dplyr::n(),
    n_pos    = sum(Result == 1, na.rm = TRUE),
    pos_rate = n_pos / n_total,
    .groups  = "drop"
  ) %>%
  # rename for display
  mutate(
    Technology = recode(Technology,
                        "Flow_Binary"       = "MFC",
                        "Adaptive_Binary"   = "clonoSEQ",
                        "BloodSensPriority_call"   = "cfWGS_PB_cfDNA")) 

## Export 
readr::write_csv(
  front_tbl,
  "Positivity_by_Landmark_TimePoint_PB_cfDNA_Frontline.csv"
)
readr::write_csv(
  non_tbl,
  "Positivity_All_TimePoints_PB_cfDNA_NonFrontline.csv"
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
#  Helper to compute pairwise concordance ---------------------------------
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

# ---------------------------------------------------------------------------
# 3.  FRONTLINE: concordance & positivity  -----------------------------------
front <- dat %>% filter(Cohort == "Frontline", !is.na(landmark))

# --- 3a. Pairwise concordance at Post‑ASCT ----------------------------------
pa   <- front %>% filter(landmark == "Post_ASCT")
post_conc <- bind_rows(
  pair_concord(pa, "BloodSensPriority_call", "Adaptive_Binary"),
  pair_concord(pa, "BloodSensPriority_call", "Flow_Binary"),
  pair_concord(pa, "Adaptive_Binary",  "Flow_Binary")
)

# --- 3b. Pairwise concordance at Maintenance --------------------------------
ma   <- front %>% filter(landmark == "Maintenance")
maint_conc <- bind_rows(
  pair_concord(ma, "BloodSensPriority_call", "Adaptive_Binary"),
  pair_concord(ma, "BloodSensPriority_call", "Flow_Binary")
)

# --- 3c. Positivity counts ---------------------------------------------------
pos_tbl <- front %>%
  filter(landmark %in% c("Post_ASCT","Maintenance")) %>%
  pivot_longer(cols = c(BloodSensPriority_call, Adaptive_Binary, Flow_Binary),
               names_to = "Test", values_to = "Result") %>%
  drop_na(Result) %>%
  group_by(landmark, Test) %>%
  summarise(
    pos = sum(Result == 1),
    tot = dplyr::n(),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------
# 4.  FRONTLINE PPV / NPV for BloodSensPriority_call --------------------------------
ppv_npv <- function(df, pred_col = "BloodSensPriority_call") {
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
  filter(!is.na(BloodSensPriority_call), !is.na(MRD_truth)) %>%
  tabyl(BloodSensPriority_call, MRD_truth)

non_cm <- non_cm %>%
  mutate(BloodSensPriority_call = as.character(BloodSensPriority_call))

# Ensure both levels 0 and 1 appear as rows and columns:
if (!all(c("0","1") %in% rownames(non_cm))) {
  non_cm <- complete(non_cm, BloodSensPriority_call = c("0","1"), fill = list(`0`=0, `1`=0))
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
  pivot_longer(cols = c(BloodSensPriority_call, Flow_Binary),
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
  X  <- g("BloodSensPriority_call","Adaptive_Binary", post_conc)
  Y  <- n("BloodSensPriority_call","Adaptive_Binary", post_conc)
  XX <- sprintf("%.0f", 100*r("BloodSensPriority_call","Adaptive_Binary", post_conc))
  Xp <- g("BloodSensPriority_call","Flow_Binary", post_conc)
  Yp <- n("BloodSensPriority_call","Flow_Binary", post_conc)
  XXp<- sprintf("%.0f", 100*r("BloodSensPriority_call","Flow_Binary", post_conc))
  Z  <- g("Adaptive_Binary","Flow_Binary", post_conc)
  W  <- n("Adaptive_Binary","Flow_Binary", post_conc)
  YY <- sprintf("%.0f",100*r("Adaptive_Binary","Flow_Binary", post_conc))
  
  # discordant counts
  n_cf_pos_cl_neg <- post_conc %>%
    filter(test_a=="BloodSensPriority_call", test_b=="Adaptive_Binary") %>%
    pull(a_pos_b_neg)
  m_cf_neg_cl_pos <- post_conc %>%
    filter(test_a=="BloodSensPriority_call", test_b=="Adaptive_Binary") %>%
    pull(a_neg_b_pos)
  
  # maintenance
  A  <- g("BloodSensPriority_call","Adaptive_Binary", maint_conc)
  B  <- n("BloodSensPriority_call","Adaptive_Binary", maint_conc)
  AA <- sprintf("%.0f",100*r("BloodSensPriority_call","Adaptive_Binary", maint_conc))
  C  <- g("BloodSensPriority_call","Flow_Binary", maint_conc)
  D  <- n("BloodSensPriority_call","Flow_Binary", maint_conc)
  BB <- sprintf("%.0f",100*r("BloodSensPriority_call","Flow_Binary", maint_conc))
  
  p <- ppv_post$PPV; q <- ppv_post$NPV
  p2<- ppv_maint$PPV; q2<- ppv_maint$NPV
  
  para <- glue("
    At post-ASCT, cfWGS agreed with clonoSEQ in {X}/{Y} ({XX}%) samples and with MFC in {Xp}/{Yp} ({XXp}%). 
    clonoSEQ vs. MFC were concordant in {Z}/{W} ({YY}%) paired samples. 
    Of the discordant post-ASCT samples, cfWGS was positive/ clonoSEQ negative in {n_cf_pos_cl_neg} cases and negative/ clonoSEQ positive in {m_cf_neg_cl_pos}. 
    At the 1-year maintenance timepoint, cfWGS agreed with clonoSEQ in {A}/{B} ({AA}%) samples and with MFC in {C}/{D} ({BB}%). 
    The PPV and NPV of cfWGS were {sprintf('%.0f',p*100)}% and {sprintf('%.0f',q*100)}% at post-ASCT, and {sprintf('%.0f',p2*100)}% and {sprintf('%.0f',q2*100)}% at maintenance. 
    In the non-frontline cohort, sensitivity and specificity of cfWGS were {sprintf('%.0f',stats_out$nonfront_sens*100)}% and {sprintf('%.0f',stats_out$nonfront_spec*100)}%, with an overall positivity rate of {stats_out$nonfront_pos %>% filter(Test=='BloodSensPriority_call') %>% summarise(sprintf('%.0f%%', 100*pos/tot)) %>% pull()}.
  ")
  
  cat(para)
}

## Export 
# 1. Export post‐ASCT pairwise concordance for bone marrow (frontline cohort)
readr::write_csv(
  post_conc,
  "Frontline_PB_cfDNA_PostASCT_Pairwise_Concordance.csv"
)

# 2. Export maintenance‐timepoint pairwise concordance for bone marrow (frontline cohort)
readr::write_csv(
  maint_conc,
  "Frontline_PB_cfDNA_Maintenance_Pairwise_Concordance.csv"
)

# 3. Export frontline positivity counts by test and landmark for bone marrow
#    (this includes Post_ASCT and Maintenance combined in one table)
readr::write_csv(
  pos_tbl,
  "Frontline_PB_cfDNA_Positivity_PostASCT_and_Maintenance.csv"
)

# 4. Export PPV/NPV at post‐ASCT for BM_zscore_only_call (frontline bone marrow)
readr::write_csv(
  ppv_post,
  "Frontline_PB_cfDNA_PostASCT_PPV_NPV.csv"
)

# 5. Export PPV/NPV at maintenance for BM_zscore_only_call (frontline bone marrow)
readr::write_csv(
  ppv_maint,
  "Frontline_PB_cfDNA_Maintenance_PPV_NPV.csv"
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
  filter(!is.na(BloodSensPriority_call),           # cfWGS call present
         !is.na(Adaptive_Binary)) %>%           # clonoSEQ call present
  mutate(category = case_when(
    BloodSensPriority_call == 1 & Adaptive_Binary == 0 ~ "cfWGS_pos / clonoSEQ_neg",
    BloodSensPriority_call == 0 & Adaptive_Binary == 1 ~ "cfWGS_neg / clonoSEQ_pos",
    TRUE                                            ~ "concordant"
  ))

disc_cf_pos <- discord_tbl %>% filter(category == "cfWGS_pos / clonoSEQ_neg")
disc_cf_neg <- discord_tbl %>% filter(category == "cfWGS_neg / clonoSEQ_pos")


outdir <- "discordant_tables"
dir.create(outdir, showWarnings = FALSE)

disc_cf_pos %>%
  select(all_of(id_cols), BloodSensPriority_call, Adaptive_Binary,
         # Flow_Binary if you want it too
         all_of(aux_cols)) %>%
  arrange(Adaptive_Frequency) %>%                        # lowest VAF at top
  write.csv(file.path(outdir, "cfWGSpos_clSEQneg_PB_cfDNA.csv"), row.names = FALSE)

disc_cf_neg %>%
  select(all_of(id_cols), BloodSensPriority_call, Adaptive_Binary,
         all_of(aux_cols)) %>%
  arrange(Adaptive_Frequency) %>%
  write.csv(file.path(outdir, "cfWGSneg_clSEQpos_PB_cfDNA.csv"), row.names = FALSE)


## Now for MFC 

# -----------------------------------------------------------------------------
# 2.  Build a “discordance” table for cfWGS (BloodSensPriority_call) vs MFC  -----
discord_flow_vs_bm <- dat %>%
  filter(!is.na(BloodSensPriority_call),    # cfWGS call present
         !is.na(Flow_Binary))  %>%          # MFC call present
  mutate(category = case_when(
    BloodSensPriority_call == 1 & Flow_Binary == 0 ~ "cfWGS_pos / MFC_neg",
    BloodSensPriority_call == 0 & Flow_Binary == 1 ~ "cfWGS_neg / MFC_pos",
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
    BloodSensPriority_call, Flow_Binary,
    all_of(aux_cols)
  ) %>%
  arrange(Flow_pct_cells) %>%   # sort by MFC percentage (lowest first)
  write.csv(
    file.path(outdir, "cfWGSpos_MFCneg_PB_cfDNA.csv"),
    row.names = FALSE
  )

# Disconnect where cfWGS negative but MFC positive
disc_cf_neg_flow %>%
  select(
    all_of(id_cols),
    BloodSensPriority_call, Flow_Binary,
    all_of(aux_cols)
  ) %>%
  arrange(Flow_pct_cells) %>%   # MFC percentage (lowest first)
  write.csv(
    file.path(outdir, "cfWGSneg_MFCpos_PB_cfDNA.csv"),
    row.names = FALSE
  )



### Run stats on discordant cases 
library(exact2x2)      # for exact McNemar 
# library(logistf)     # if you decide on Firth logistic

## 3.1  Build flags -----------------------------------------------------------
df_disc <- dat %>%
  filter(!is.na(BloodSensPriority_call)) %>%  # only require cfWGS to be non-NA
  mutate(
    disc_cf_vs_clono = case_when(
      !is.na(Adaptive_Binary) & BloodSensPriority_call != Adaptive_Binary ~ 1L,
      !is.na(Adaptive_Binary)                                          ~ 0L,
      TRUE                                                             ~ NA_integer_
    ),
    disc_cf_vs_flow = case_when(
      !is.na(Flow_Binary) & BloodSensPriority_call != Flow_Binary ~ 1L,
      !is.na(Flow_Binary)                                      ~ 0L,
      TRUE                                                     ~ NA_integer_
    )
  )


## 3.2  Exact McNemar for cfWGS vs clonoSEQ at post‑ASCT ----------------------
post <- df_disc %>% filter(timepoint_info == "Post_transplant")

tbl_mc <- table(
  cfWGS = post$BloodSensPriority_call,
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
  select(Patient, Sample_Code, BloodSensPriority_call, Adaptive_Binary,
         Flow_Binary, Adaptive_Frequency, Flow_pct_cells,
         BM_zscore_only_prob, FS, WGS_Tumor_Fraction_Blood_plasma_cfDNA,
         Mean.Coverage, BM_MutCount_Baseline, Blood_MutCount_Baseline) %>%
  arrange(desc(Adaptive_Frequency)) %>%
  write.csv("discordant_postASCT_summary_PB_cfDNA.csv", row.names = FALSE)



### Now do for maintenance 
# -----------------------------------------------------------------------------
# 3.2  Exact McNemar + Wilcoxon + export for MAINTENANCE ----------------------
maint <- df_disc %>% filter(timepoint_info == "1yr maintenance")

# 3.2a  McNemar’s test (cfWGS vs. clonoSEQ) at maintenance
tbl_mc_maint <- table(
  cfWGS = maint$BloodSensPriority_call,
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
    BloodSensPriority_call, Adaptive_Binary, Flow_Binary,
    Adaptive_Frequency, Flow_pct_cells,
    BM_zscore_only_prob, FS, WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    Mean.Coverage, BM_MutCount_Baseline, Blood_MutCount_Baseline
  ) %>%
  arrange(desc(Adaptive_Frequency)) %>%
  write.csv("discordant_maintenance_summary_PB_cfDNA.csv", row.names = FALSE)
