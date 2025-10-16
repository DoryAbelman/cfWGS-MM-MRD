# ===============================================================================
# Script: 2_4_Longitudinal_features_analysis.R
# Author: Dory Abelman
# Date: 2025-06
# Purpose:
#   - Compute longitudinal MRD feature statistics and visualizations
#   - Handles NA values and dynamic column presence
#   - Identifies baseline and first follow-up samples using date logic
#   - Summarizes cohort-level counts and paired feature changes
#   - Generates publication-quality plots (violin/box, spaghetti, correlation)
# Dependencies:
#   tidyverse, lubridate, rstatix, patchwork, ggpubr, glue
# Inputs:
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.rds
#   - cohort_assignment_table_updated.rds
# Outputs:
#   - CSV tables in 'Output_tables_2025'
#   - PDF figures: longitudinal pairs, spaghetti, correlation plots
# ===============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(rstatix)
  library(patchwork)
  library(ggpubr)
  library(glue)
  library(GGally)
  library(scales)
  library(viridis)
})

## ───── 1. Load data ──────────────────────────────────────────────────────────
### Set paths 
outdir   <- "Output_tables_2025"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

### Load data 
file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds")
cohort_df <- readRDS("cohort_assignment_table_updated.rds")

dat <- file 

## Load good patients 
Good_pts <- read.csv("baseline_high_quality_patients_updated.csv",
                    stringsAsFactors = FALSE) ## From 2_0 script

# 1.  Join cohort_df and keep frontline only -------------------------------------
dat <- dat %>%                # <‑‑ your master data
  left_join(cohort_df, by = "Patient") 

dat <- dat %>% mutate(Date = as.Date(Date))   # safety

## Filter to baseline for simplicity and ones have data on  
dat <- dat %>% filter(Cohort == "Frontline")

dat <- dat %>%
  filter(!is.na(zscore_BM) | !is.na(zscore_blood) | !is.na(FS))

## Set na the ones without evidence of disease? 
# Load the RDS file and compute the good-patient sets
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds")

BM_good_pts <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_BM_cells == 1) %>%
  pull(Patient) %>%
  unique()

cfDNA_good_pts <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
  pull(Patient) %>%
  unique()


# Filter dat to rows with evidence of disease
bm_feats    <- c("zscore_BM", "z_score_detection_rate_BM", "detect_rate_BM", "sites_rate_BM")
blood_feats <- c("zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood", "sites_rate_blood")

dat <- dat %>%
  # keep rows with at least one key measurement
  filter(!is.na(zscore_BM) | !is.na(zscore_blood) | !is.na(FS)) %>%
  mutate(
    across(all_of(bm_feats),
           ~ if_else(Patient %in% BM_good_pts,    .x, NA_real_)),
    across(all_of(blood_feats),
           ~ if_else(Patient %in% cfDNA_good_pts, .x, NA_real_))
  )

# Remove IMG-127-T15 ichorCNA since not a complete case? But have BM for it so keep in count? 
dat2 <- dat %>%
  mutate(
    WGS_Tumor_Fraction_Blood_plasma_cfDNA =
      if_else(Patient == "IMG-127" & Date == as.Date("2022-08-15"),
              NA_real_,
              WGS_Tumor_Fraction_Blood_plasma_cfDNA)
  )

## ───── 2. Identify BASELINE and FIRST-FOLLOW rows (start-date logic) ────────
baseline_tbl <- dat2 %>%
  group_by(Patient) %>%
  filter(
    timepoint_info %in% c("Diagnosis", "Baseline") |
      Date == min(Date, na.rm = TRUE)
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(Patient, Date,
         everything())                  # keep all cols - easier later

follow1_tbl <- dat2 %>%
  anti_join(baseline_tbl, by = c("Patient", "Date")) %>%
  filter(timepoint_info != "Relapse") %>% # omit these 
  group_by(Patient) %>%
  filter(Date == min(Date, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  select(Patient, Date, everything())

## prefix helpers
rename_pref <- function(df, suff) {
  rename_with(df, ~ paste0(., "_", suff), -c(Patient, Date)) %>%
    rename(!!paste0("Date_", suff) := Date)
}

paired_master <- rename_pref(baseline_tbl, "baseline") %>%
  inner_join(rename_pref(follow1_tbl, "follow1"), by = "Patient")

## ───── 3. Cohort-level summary counts  ────────────────
n_patients     <- n_distinct(paired_master$Patient)
n_cfDNA_draws  <- dat %>% filter(!is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA)) %>% nrow()
n_cfDNA_draws_matched_baseline <- dat %>% filter(!is.na(zscore_BM) | !is.na(zscore_blood)) %>% nrow()

median_follow <- dat %>%
  filter(!is.na(zscore_BM) | !is.na(zscore_blood)) %>%
  group_by(Patient) %>%
  summarise(months = interval(min(Date, na.rm = TRUE),
                              max(Date, na.rm = TRUE)) %/% months(1)) %>%
  summarise(med = median(months, na.rm = TRUE),
            lo  = min(months,   na.rm = TRUE),
            hi  = max(months,   na.rm = TRUE))

sentence <- glue(
  "Across the longitudinal series ({n_cfDNA_draws_matched_baseline} cfDNA draws with matched baseline samples ",
  "from {n_patients} patients, ",
  "median follow-up = {median_follow$med} months, range ",
  "{median_follow$lo}–{median_follow$hi}),"
)

cat(sentence, "\n")

### Do seperately by feature 
# ----- define the feature set you consider "cfDNA features"
features_any <- c(
  "zscore_BM", "z_score_detection_rate_BM", "detect_rate_BM",
  "zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood",
  "FS", "Mean.Coverage", "Proportion.Short",
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA"
)
features_any <- intersect(features_any, names(dat))  # only keep columns that exist

# =============  A) ANY FEATURE (>=1 present on a draw)  =============
dat_any <- dat %>%
  filter(if_any(all_of(features_any), ~ !is.na(.x)))

n_patients_any <- n_distinct(dat_any$Patient)
n_draws_any    <- nrow(dat_any)

median_follow_any <- dat_any %>%
  group_by(Patient) %>%
  summarise(months = lubridate::interval(min(Date, na.rm = TRUE),
                                         max(Date, na.rm = TRUE)) %/% months(1),
            .groups = "drop") %>%
  summarise(med = median(months, na.rm = TRUE),
            lo  = min(months,   na.rm = TRUE),
            hi  = max(months,   na.rm = TRUE))

sentence_any <- glue::glue(
  "Across the longitudinal series ({n_draws_any} cfDNA draws with ≥1 feature ",
  "from {n_patients_any} patients, median follow-up = {median_follow_any$med} months, ",
  "range {median_follow_any$lo}–{median_follow_any$hi})."
)

# =============  B) zscore_BM universe  =============
dat_bm <- dat %>% filter(!is.na(zscore_BM))

n_patients_bm <- n_distinct(dat_bm$Patient)
n_draws_bm    <- nrow(dat_bm)

median_follow_bm <- dat_bm %>%
  group_by(Patient) %>%
  summarise(months = lubridate::interval(min(Date, na.rm = TRUE),
                                         max(Date, na.rm = TRUE)) %/% months(1),
            .groups = "drop") %>%
  summarise(med = median(months, na.rm = TRUE),
            lo  = min(months,   na.rm = TRUE),
            hi  = max(months,   na.rm = TRUE))

sentence_bm <- glue::glue(
  "z-score (BM lists): {n_draws_bm} draws from {n_patients_bm} patients; ",
  "follow-up median {median_follow_bm$med} months (range {median_follow_bm$lo}–{median_follow_bm$hi})."
)

# =============  C) zscore_blood universe  =============
dat_blood <- dat %>% filter(!is.na(zscore_blood))

n_patients_blood <- n_distinct(dat_blood$Patient)
n_draws_blood    <- nrow(dat_blood)

median_follow_blood <- dat_blood %>%
  group_by(Patient) %>%
  summarise(months = lubridate::interval(min(Date, na.rm = TRUE),
                                         max(Date, na.rm = TRUE)) %/% months(1),
            .groups = "drop") %>%
  summarise(med = median(months, na.rm = TRUE),
            lo  = min(months,   na.rm = TRUE),
            hi  = max(months,   na.rm = TRUE))

sentence_blood <- glue::glue(
  "z-score (blood lists): {n_draws_blood} draws from {n_patients_blood} patients; ",
  "follow-up median {median_follow_blood$med} months (range {median_follow_blood$lo}–{median_follow_blood$hi})."
)

# Print the three sentences for quick copy/paste into the manuscript
cat(sentence_any, "\n")
cat(sentence_bm,  "\n")
cat(sentence_blood, "\n")




## ───── 4. Feature list  ────────────────────────────────────────
features <- c("zscore_BM", "z_score_detection_rate_BM", "detect_rate_BM",
              "zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood",
              "FS", "Mean.Coverage", "Proportion.Short",
              "WGS_Tumor_Fraction_Blood_plasma_cfDNA")

## convenience: only keep those that actually exist
features <- features[features %in% names(dat)]

## ───── 5.Paired stats per feature ────────────────────────────────
pairwise_stats <- map_dfr(features, function(f) {
  # 1. Pull out only the non-missing draws for this feature, per patient
  df <- dat %>%
    select(Patient, Date, value = !!sym(f)) %>%
    filter(!is.na(value)) %>%
    group_by(Patient) %>%
    mutate(cnt = dplyr::n()) %>%       # count per patient
    filter(cnt >= 2) %>%        # need at least 2 draws
    arrange(Date) %>%
    slice_head(n = 2) %>%       # first two timepoints
    ungroup()
  
  if (nrow(df) == 0) return(NULL)
  
  # 2. Spread into baseline / follow1
  df2 <- df %>%
    group_by(Patient) %>%
    arrange(Date) %>%
    mutate(Time = c("baseline","follow1")) %>%
    ungroup() %>%
    select(Patient, Time, value) %>%
    pivot_wider(
      id_cols     = Patient,
      names_from  = Time,
      values_from = value,
      names_glue  = "{.value}_{Time}"
    ) %>%
    drop_na()  
  
  # 3. Extract vectors, counts, medians
  x_bl   <- df2$value_baseline
  x_f1   <- df2$value_follow1
  n_pairs <- length(x_bl)
  med_bl  <- median(x_bl, na.rm = TRUE)
  med_f1  <- median(x_f1, na.rm = TRUE)
  
  # 4. Wilcoxon only if at least 2 pairs
  pval <- if (n_pairs >= 2) {
    wilcox.test(df2$value_baseline, df2$value_follow1,
                paired = TRUE, exact = FALSE)$p.value
  } else {
    NA_real_
  }
  
  tibble::tibble(
    Feature         = f,
    N_pairs         = n_pairs,
    Mean_baseline    = mean(x_bl, na.rm = TRUE),
    SD_baseline      = sd(x_bl,   na.rm = TRUE),
    Median_baseline = med_bl,
    Mean_follow1     = mean(x_f1, na.rm = TRUE),
    SD_follow1       = sd(x_f1,   na.rm = TRUE),
    Median_follow1  = med_f1,
    Percent_change   = (median(x_f1) - median(x_bl)) / median(x_bl) * 100,
    Fold_change   = (median(x_f1)) / median(x_bl),
    p_Wilcoxon      = pval
  )
}) %>%
  mutate(q_BH = p.adjust(p_Wilcoxon, method = "BH"))

## Write updated stats to outdir
write_csv(pairwise_stats, file.path(outdir, "placeholder_stats_v4.csv"))

### See which patients are included 
# Return the patients contributing baseline→first-post pairs for a given column
pairwise_included <- purrr::map(features, function(f) {
  df <- dat %>%
    select(Patient, Date, value = !!sym(f)) %>%
    filter(!is.na(value)) %>%
    group_by(Patient) %>%
    mutate(cnt = dplyr::n()) %>%
    filter(cnt >= 2) %>%          # need at least 2 draws
    arrange(Date) %>%
    slice_head(n = 2) %>%         # first two timepoints (as in your stats)
    ungroup()
  
  if (nrow(df) == 0) {
    return(tibble(Feature = f, N_pairs = 0L, Included_Patients = list(character(0))))
  }
  
  # Spread to baseline/follow1 as in your stats, then keep only complete pairs
  df2 <- df %>%
    group_by(Patient) %>%
    arrange(Date) %>%
    mutate(Time = c("baseline","follow1")) %>%
    ungroup() %>%
    select(Patient, Time, value) %>%
    pivot_wider(
      id_cols     = Patient,
      names_from  = Time,
      values_from = value,
      names_glue  = "{.value}_{Time}"
    ) %>%
    drop_na()
  
  tibble(
    Feature          = f,
    N_pairs          = nrow(df2),
    Included_Patients = list(sort(unique(df2$Patient)))
  )
}) %>%
  bind_rows()

## (Optional) join onto your pairwise_stats to see N agree side-by-side
pairwise_stats_with_ids <- pairwise_stats %>%
  left_join(pairwise_included %>% select(Feature, Included_Patients),
            by = "Feature")

for (i in seq_len(nrow(pairwise_included))) {
  f   <- pairwise_included$Feature[i]
  n   <- pairwise_included$N_pairs[i]
  ids <- pairwise_included$Included_Patients[[i]]
  cat("\n=== Feature:", f, "===\nN_pairs =", n, "\nPatients:\n")
  if (length(ids)) {
    cat(paste(ids, collapse = ", "), "\n")
  } else {
    cat("(none)\n")
  }
}

## Check labs 
### Assemble paragraph 
# helper to format p-values
format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# 5. Extract values for paragraph
format_p <- function(p) {
  if (is.na(p)) return("NA"); if (p<0.001) return("<0.001"); sprintf("%.3f",p)
}
bm_site <- filter(pairwise_stats,Feature=="zscore_BM")
bm_vafz <- filter(pairwise_stats,Feature=="z_score_detection_rate_BM")
bm_rate <- filter(pairwise_stats,Feature=="detect_rate_BM")
cf_site <- filter(pairwise_stats,Feature=="zscore_blood")
cf_vafz <- filter(pairwise_stats,Feature=="z_score_detection_rate_blood")
cf_rate <- filter(pairwise_stats,Feature=="detect_rate_blood")
fs_feat <- filter(pairwise_stats,Feature=="FS")
mc_feat <- filter(pairwise_stats,Feature=="Mean.Coverage")
ps_feat <- filter(pairwise_stats,Feature=="Proportion.Short")
tf_feat <- filter(pairwise_stats,Feature=="WGS_Tumor_Fraction_Blood_plasma_cfDNA")

# 6. Key Spearman correlations
clin_cor <- cor.test(dat$detect_rate_blood, dat$MRD_Clinical_Binary,
                     method="spearman",use="pairwise.complete.obs")
rho_clin <- round(clin_cor$estimate,2); p_clin <- format_p(clin_cor$p.value)

tf_bm_cor <- cor.test(dat$WGS_Tumor_Fraction_Blood_plasma_cfDNA, dat$Blood_Mutation_Count,
                      method="spearman",use="pairwise.complete.obs")
rho_tf_bm <- round(tf_bm_cor$estimate,2); p_tf_bm <- format_p(tf_bm_cor$p.value)

tf_fs_cor <- cor.test(dat$WGS_Tumor_Fraction_Blood_plasma_cfDNA, dat$FS,
                      method="spearman",use="pairwise.complete.obs")
rho_tf_fs <- round(tf_fs_cor$estimate,2); p_tf_fs <- format_p(tf_fs_cor$p.value)

# 7. Assemble paragraph
paragraph <- glue(
  "In BM, the site-detection z-score declined from {round(bm_site$Median_baseline,2)} to {round(bm_site$Median_follow1,2)} ",
  "(fold = {round(bm_site$Fold_change,2)}×; p = {format_p(bm_site$p_Wilcoxon)}), and the cumulative VAF z-score from ",
  "{round(bm_vafz$Median_baseline,2)} to {round(bm_vafz$Median_follow1,2)} ",
  "(fold = {round(bm_vafz$Fold_change,2)}×; p = {format_p(bm_vafz$p_Wilcoxon)}). The raw detection rate also decreased significantly, ",
  "dropping by {round((bm_rate$Median_follow1 - bm_rate$Median_baseline)*100,2)} percentage points ",
  "({round(bm_rate$Median_baseline*100,2)}%→{round(bm_rate$Median_follow1*100,2)}%; ",
  "fold = {round(bm_rate$Fold_change,2)}×; p = {format_p(bm_rate$p_Wilcoxon)}).",
  
  "\n\nIn cfDNA, the blood site-detection z-score dropped from {round(cf_site$Median_baseline,2)} to {round(cf_site$Median_follow1,2)} ",
  "(fold = {round(cf_site$Fold_change,2)}×; p = {format_p(cf_site$p_Wilcoxon)}), the cumulative VAF z-score from ",
  "{round(cf_vafz$Median_baseline,2)} to {round(cf_vafz$Median_follow1,2)} ",
  "(fold = {round(cf_vafz$Fold_change,2)}×; p = {format_p(cf_vafz$p_Wilcoxon)}), and the raw blood detection rate from ",
  "{round(cf_rate$Median_baseline*100,2)}% to {round(cf_rate$Median_follow1*100,2)}% ",
  "(fold = {round(cf_rate$Fold_change,2)}×; p = {format_p(cf_rate$p_Wilcoxon)}). ",
  "These tracked closely with clinical response: Spearman ρ = {rho_clin} (p = {p_clin}).",
  
  "\n\nFragmentomic signals moved in tandem: the fragment-size score contracted by ",
  "{round((fs_feat$Median_follow1 - fs_feat$Median_baseline)/abs(fs_feat$Median_baseline)*100,2)}% ",
  "(fold = {round(fs_feat$Fold_change,2)}×; p = {format_p(fs_feat$p_Wilcoxon)}), the mean coverage over MM enhancers fell from ",
  "{round(mc_feat$Median_baseline,2)}× to {round(mc_feat$Median_follow1,2)}× ",
  "(fold = {round(mc_feat$Fold_change,2)}×; p = {format_p(mc_feat$p_Wilcoxon)}), and the proportion of short fragments ",
  "decreased from {round(ps_feat$Median_baseline*100,2)}% to {round(ps_feat$Median_follow1*100,2)}% ",
  "(fold = {round(ps_feat$Fold_change,2)}×; p = {format_p(ps_feat$p_Wilcoxon)}). ",
  "Meanwhile, ichorCNA tumour fraction dropped from {round(tf_feat$Median_baseline*100,2)}% to ",
  "{round(tf_feat$Median_follow1*100,2)}% (fold = {round(tf_feat$Fold_change,2)}×; p = {format_p(tf_feat$p_Wilcoxon)}), ",
  "remaining strongly correlated with blood mutation burden (ρ = {rho_tf_bm}; p = {p_tf_bm}) and FS (ρ = {rho_tf_fs}; p = {p_tf_fs})."
)

# print it
cat(paragraph, "\n")



## ───── 6. More correlations and with labs (pairwise complete obs already) ────
# ── 0. Clean the calcium outlier ──────────────────────────────────────────────
dat_clean <- dat %>%
  mutate(Calcium = if_else(Calcium > 10, NA_real_, Calcium))

# ── 1. Optionally coerce stage columns to numeric levels ----------------------
dat_clean <- dat_clean %>%
  mutate(
    ISS_STAGE_num  = as.numeric(factor(ISS_STAGE,  levels = c("Stage I","Stage II","Stage III"))),
    R_ISS_STAGE_num = as.numeric(factor(R_ISS_STAGE, ordered = TRUE))
  )


extras <- c("AGE", "ISS_STAGE_num", "R_ISS_STAGE_num", "ECOG_SCORE", "KPS_SCORE",
            "Albumin", "B2_micro", "Calcium", "Creatinine", "Hemoglobin",
            "LDH", "dFLC", "M_Protein", "Kappa_Lambda_Ratio",
            "Flow_pct_cells", "Adaptive_Frequency")

core_feats <- c("zscore_blood", "detect_rate_blood", "z_score_detection_rate_blood",
                "zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM",
                "WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Blood_Mutation_Count",
                "FS", "Mean.Coverage", "Proportion.Short")

vars_to_corr <- intersect(c(core_feats, extras), names(dat_clean))

# ── 3. Keep ONLY numeric columns ---------------------------------------------
numeric_vars <- vars_to_corr[sapply(dat_clean[vars_to_corr], is.numeric)]

if (length(numeric_vars) < 2) stop("Need ≥2 numeric variables to correlate.")

# ── 4. All-against-all correlations ------------------------------------------
var_pairs <- combn(numeric_vars, 2, simplify = FALSE)

# 4) Compute correlations + normality checks
corr_table <- map_dfr(var_pairs, function(v) {
  x_all <- dat_clean[[v[1]]]
  y_all <- dat_clean[[v[2]]]
  ok    <- complete.cases(x_all,y_all)
  n     <- sum(ok)
  if(n < 3) return(NULL)
  
  x <- x_all[ok]; y <- y_all[ok]
  
  # 4a) Normality tests if n <= 5000
  normal_x <- if(n <= 5000) shapiro.test(x)$p.value > 0.05 else FALSE
  normal_y <- if(n <= 5000) shapiro.test(y)$p.value > 0.05 else FALSE
  
  # 4b) Spearman and Pearson estimates + p-values
  rho_sp    <- suppressWarnings(cor(x,y,method="spearman"))
  p_sp      <- suppressWarnings(cor.test(x,y,method="spearman")$p.value)
  rho_pe    <- suppressWarnings(cor(x,y,method="pearson"))
  test_pe   <- suppressWarnings(cor.test(x,y,method="pearson"))
  p_pe      <- test_pe$p.value
  
  # 4c) Decide which method to use
  method_use <- if(normal_x && normal_y) "Pearson" else "Spearman"
  
  tibble(
    var1          = v[1],
    var2          = v[2],
    n_pairs       = n,
    normal_x      = normal_x,
    normal_y      = normal_y,
    rho_spearman  = round(rho_sp, 3),
    p_spearman    = signif(p_sp, 3),
    rho_pearson   = round(rho_pe, 3),
    p_pearson     = signif(p_pe, 3),
    method_use    = method_use
  )
}) %>%
  ungroup() %>%
  # 5) Optionally FDR-correct separately for each method
  mutate(
    q_BH_spearman = p.adjust(p_spearman, method="BH"),
    q_BH_pearson = p.adjust(p_pearson,  method="BH")
  )

### Get significant 
# 1) Identify MRD metrics vs lab correlations with q_BH < 0.05
core_feats <- c("zscore_BM","z_score_detection_rate_BM","detect_rate_BM",
                "zscore_blood","z_score_detection_rate_blood","detect_rate_blood",
                "FS","Mean.Coverage","Proportion.Short",
                "WGS_Tumor_Fraction_Blood_plasma_cfDNA")

lab_feats  <- c("Albumin","B2_micro","Calcium","Creatinine","Hemoglobin",
                "LDH","dFLC","M_Protein","Kappa_Lambda_Ratio",
                "Flow_pct_cells","Adaptive_Frequency")

sig_corrs <- corr_table %>%
  filter(
    q_BH_spearman < 0.05,
    ((var1 %in% core_feats) & (var2 %in% lab_feats)) |
      ((var2 %in% core_feats) & (var1 %in% lab_feats))
  ) %>%
  # normalize so var1 always = MRD metric
  mutate(
    m   = if_else(var1 %in% core_feats, var1, var2),
    lab = if_else(var1 %in% lab_feats, var1, var2)
  ) %>%
  select(m, lab, n_pairs, rho_spearman, rho_pearson, p_spearman, p_pearson, q_BH_spearman)

print(sig_corrs, n = Inf)




### Export what have so far
write_csv(pairwise_stats, file.path(outdir, "paired_stats_summary2.csv"))
write_csv(corr_table, file.path(outdir, "mrd_metric_lab_all_correlations2.csv"))
write_csv(sig_corrs, file.path(outdir, "mrd_metric_lab_sig_correlations2.csv"))







## ───── 7. Figures – all filtered on !is.na for the plotted feature ──────────
## Helper for paired violin/box with spaghetti
## 9.1  Utilities ------------------------------------------------------------
baseline_dates <- readRDS("Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds")
baseline_dates$Baseline_Date <- baseline_dates$baseline_date # for consistency

  dat <- dat_clean %>%
  left_join(baseline_dates, by = "Patient") %>%
  mutate(
    Weeks_Since_Baseline = as.numeric(difftime(Date, Baseline_Date, units = "weeks")),
    Weeks_Since_Baseline = case_when(
      Weeks_Since_Baseline >= -2 & Weeks_Since_Baseline < 0 ~ 0, # minor fluctuations
      TRUE                                                ~ Weeks_Since_Baseline
    )
  )

# 1. Pivot to long form
plot_df <- dat %>%
  select(
    Patient,
    Weeks_Since_Baseline,
    Num_days_to_closest_relapse,
    cVAF       = detect_rate_BM,
    cVAF_z     = z_score_detection_rate_BM,
    sites      = sites_rate_BM,
    sites_z    = zscore_BM
  ) %>%
  mutate(
    relapse_within_180 = if_else(
      Num_days_to_closest_relapse <= 180,
      TRUE,
      FALSE,
      missing = FALSE    # recode any NA here to FALSE since didn't relapse
    )) %>%
  pivot_longer(
    cols = c(cVAF, cVAF_z, sites, sites_z),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  drop_na(Value)

# 2. Plot
# 1. Build the plot and save to a variable
custom_labels <- c(
  cVAF    = "Cumulative VAF",
  cVAF_z  = "Cumulative VAF Z-score",
  sites   = "Proportion of Sites Detected",
  sites_z = "Proportion of Sites Detected Z-score"
)

p_traj <- ggplot(plot_df, aes(x = Weeks_Since_Baseline, y = Value, group = Patient)) +
  geom_line(colour = "grey70", alpha = 0.4) +
  geom_point(aes(color = relapse_within_180), size = 1.5, alpha = 0.8) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol  = 4,
    labeller = labeller(Metric = custom_labels)
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180 d"),
    name   = "Relapse status",
    na.value = "black",         # <- this makes NA points black
    na.translate = FALSE      # ← drop NA from the legend since these patients didn't progress
    
  ) +
  labs(
    x     = "Weeks Since Baseline",
    y     = "Value",
    title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

# 2. Save as PNG in outdir
ggsave(
  filename = file.path(outdir, "Fig3A_metrics_trajectories_BM.png"),
  plot     = p_traj,
  device   = "png",
  width    = 12,    # inches
  height   = 4,     # inches
  dpi      = 600
)


### Edit style 
# ─────────────────────────────────────────────────────────────
# Prep: derive a patient-level relapse flag (any TP ≤180d)
# ─────────────────────────────────────────────────────────────
plot_df_pat <- plot_df %>%
  group_by(Patient) %>%
  mutate(
    patient_relapse180 = any(replace_na(relapse_within_180, FALSE))
  ) %>%
  ungroup() %>%
  mutate(
    # make it an ordered factor for consistent legend / facet order
    patient_relapse180 = factor(patient_relapse180,
                                levels = c(FALSE, TRUE))
  )

# nice labels for the patient-level factor
patient_relapse_labs <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d (last TP)"
)
patient_relapse_labs_short <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d"
)



# ===================================================================
# VERSION 1
# Lines coloured by *patient-level* relapse (any TP ≤180d)
# (All points inherit the same colour as the line.)
# ===================================================================
p_traj_patientline <- ggplot(
  plot_df_pat,
  aes(x = Weeks_Since_Baseline, y = Value, group = Patient,
      colour = patient_relapse180)
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(
    ~ Metric,
    scales   = "free_y",
    ncol     = 4,
    labeller = labeller(Metric = custom_labels)
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = patient_relapse_labs,
    name   = "Patient relapse status"
  ) +
  labs(
    x     = "Weeks Since Baseline",
    y     = "Value",
    title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

ggsave(
  filename = file.path(outdir, "Fig3A_metrics_trajectories_BM_byPatientLine_2.png"),
  plot     = p_traj_patientline,
  device   = "png",
  width    = 12,
  height   = 4,
  dpi      = 600
)



# ===================================================================
# VERSION 2
# Facet by patient-level relapse status *and* by Metric.
# ===================================================================
patient_relapse_labs_short <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d"
)

# 1) Prep: patient-level flag
plot_df_pat <- plot_df %>%
  group_by(Patient) %>%
  mutate(patient_relapse180 = any(replace_na(relapse_within_180, FALSE))) %>%
  ungroup() %>%
  mutate(
    patient_relapse180 = factor(patient_relapse180, levels = c(FALSE, TRUE)),
    relapse_within_180 = replace_na(relapse_within_180, FALSE)
  )

# 2) Build the two mini‑plots

## A) Cumulative VAF
df_cvaf <- filter(plot_df_pat, Metric == "cVAF")

## If wanted segmented line 
df_cvaf_seg <- df_cvaf %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_cvaf <- ggplot(df_cvaf, aes(Weeks_Since_Baseline, Value, group = Patient)) +
  #geom_line(color = "black", alpha = 0.6) +
  # coloured segments
  geom_segment(
    data = df_cvaf_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["cVAF"],
    x     = "Weeks Since Baseline",
    y     = "Cumulative VAF"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

## B) Proportion of Sites Detected
df_sites <- filter(plot_df_pat, Metric == "sites")
df_sites_seg <- df_sites %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_sites <- ggplot(df_sites, aes(Weeks_Since_Baseline, Value, group = Patient)) +
 # geom_line(color = "black", alpha = 0.6) +
  # coloured segments
  geom_segment(
    data = df_sites_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["sites"],
    x     = "Weeks Since Baseline",
    y     = "Prop. Mutant Sites Detected"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

# 3) Combine them side by side
p_combined <- p_cvaf + p_sites + 
  plot_annotation(
    title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  ) & 
  theme(legend.position = "bottom")

# 4) Save
ggsave(
  filename = file.path(outdir, "Fig3A_sideBySide_lockedY_segmented2.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)


# 2) Build the two mini‑plots for other features

## A) Cumulative VAF
df_cvaf <- filter(plot_df_pat, Metric == "cVAF_z")
df_cvaf_seg <- df_cvaf %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_cvaf <- ggplot(df_cvaf, aes(Weeks_Since_Baseline, Value, group = Patient)) +
 # geom_line(color = "black", alpha = 0.6) +
  # coloured segments
  geom_segment(
    data = df_cvaf_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["cVAF_z"],
    x     = "Weeks Since Baseline",
    y     = "Cumulative VAF (Z)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

# cap y values at 2500 for plotting only
## Flag capped points and keep their true values for labeling
df_cvaf_plot <- df_cvaf %>%
  mutate(
    overcap    = Value > 2500,                # TRUE if > 2500
    Value_plot = pmin(Value, 2500),           # cap for plotting
    label_val  = ifelse(overcap, round(Value), NA_real_)  # show true value
  ) %>%
  filter(!is.na(Weeks_Since_Baseline), !is.na(Value_plot))

## Segments as before, but capped for plotting
df_cvaf_seg_plot <- df_cvaf %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  ungroup() %>%
  mutate(
    y_plot    = pmin(pmax(y,    0), 2500),
    yend_plot = pmin(pmax(yend, 0), 2500)
  ) %>%
  filter(!is.na(x), !is.na(xend), !is.na(y_plot), !is.na(yend_plot))

p_cvaf_modified <- ggplot(df_cvaf_plot, aes(Weeks_Since_Baseline, Value_plot, group = Patient)) +
  geom_segment(
    data = df_cvaf_seg_plot,
    aes(x = x, y = y_plot, xend = xend, yend = yend_plot, colour = seg_relapse),
    size = 0.4, alpha = 0.6
  ) +
  # points: circle for ≤2500, triangle for >2500
  geom_point(aes(color = relapse_within_180, shape = overcap),
             size = 1.8, alpha = 0.8) +
  # labels for capped points: print the true value
  geom_text(
    data = ~ dplyr::filter(., overcap),
    aes(label = label_val),
    vjust = -0.75, size = 2
  ) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  # shape 16 = filled circle; 17 = filled triangle
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 17),
    labels = c(`FALSE` = "≤2500", `TRUE` = ">2500 (capped)"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["cVAF_z"],
    x     = "Weeks Since Baseline",
    y     = "Cumulative VAF (Z)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA),
        plot.margin = margin(5.5, 5.5, 12, 5.5)) +
  # give a little headroom so labels above 2500 aren't clipped
  coord_cartesian(ylim = c(0, 2600), clip = "on") +
  scale_y_continuous(
    breaks = c(0, 500, 1000, 1500, 2000, 2500),
    labels = function(x) ifelse(x == 2500, ">2500", x)
  )

## B) Proportion of Sites Detected
df_sites <- filter(plot_df_pat, Metric == "sites_z")
df_sites_seg <- df_sites %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_sites <- ggplot(df_sites, aes(Weeks_Since_Baseline, Value, group = Patient)) +
 # geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_sites_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["sites_z"],
    x     = "Weeks Since Baseline",
    y     = "Prop. Mutant Sites Detected (Z)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

# 3) Combine them side by side
p_combined <- p_cvaf_modified + p_sites + 
     plot_annotation(
      title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
    ) & 
      theme(legend.position = "bottom")
    
    
# 4) Save
ggsave(
  filename = file.path(outdir, "Fig3A_sideBySide_lockedY_zscore_segmented4.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)


### Redo as function but color based on change 

# ---- knobs to adjust ----
q_flat <- 0.90     # 90th percentile of non-progressor follow-up values
tol    <- 0        # require strict increase (set >0 to ignore tiny jitters)
cap    <- 2500      # no capping for Z-scores; set to 2500 for capped metrics
week_cutoff <- 250  # define where we look for the patient’s trough

# flag the *actual* clinical progression timepoint (0 days from progression)
plot_df_pat <- plot_df_pat %>%
  mutate(prog_tp = !is.na(Num_days_to_closest_relapse) &
           Num_days_to_closest_relapse == 0)

# learn thresholds from NON-progressors (patients with no relapse ≤180d)
thr_tbl <- plot_df_pat %>%
  filter(as.character(patient_relapse180) == "FALSE",
         Weeks_Since_Baseline > 52) %>%               # follow-up only
  group_by(Metric) %>%
  summarise(thr = quantile(Value, q_flat, na.rm = TRUE), .groups = "drop")

thr_named <- setNames(thr_tbl$thr, thr_tbl$Metric)

## Function to plot
make_panel_minlogic <- function(metric, data, ylab,
                                col_after = "red", col_before = "black",
                                week_cutoff = 75, tol = 0, cap = Inf,
                                pad_before_weeks = 0, plot_title = NULL) {   # <-- added here
  
  df0 <- data %>% dplyr::filter(Metric == metric)
  
  # threshold (fallback if not in thr_named)
  thr_val <- suppressWarnings(as.numeric(thr_named[[metric]]))
  if (is.na(thr_val)) {
    thr_val <- stats::quantile(df0$Value[df0$Weeks_Since_Baseline > 0],
                               q_flat, na.rm = TRUE)
  }
  
  
  # ---- mark what to show as percent ----
  percent_metrics <- c("TF_ichorCNA", "Proportion.Short", "cVAF_blood", "cVAF", "sites", "sites_blood")
  is_percent <- metric %in% percent_metrics
  
  ## --- low-tail handling ----------------------------------------------------
  low_tail_metrics <- c("Mean.Coverage")   # add more if needed
  is_low <- metric %in% low_tail_metrics
  
  # anchor (prefer ≤ week_cutoff if available): trough for high-tail, peak for low-tail
  anchor_tbl <- df0 %>%
    dplyr::group_by(Patient) %>%
    dplyr::group_modify(~{
      d     <- .x
      d_pre <- d %>% dplyr::filter(Weeks_Since_Baseline <= week_cutoff)
      pick  <- function(dd) {
        if (is_low) dd[which.max(dd$Value), , drop = FALSE]  # peak for low-tail
        else        dd[which.min(dd$Value), , drop = FALSE]  # trough for high-tail
      }
      r <- if (nrow(d_pre) > 0) pick(d_pre) else pick(d)
      tibble::tibble(anchor_time = r$Weeks_Since_Baseline[1],
                     anchor_val  = r$Value[1])
    }) %>%
    dplyr::ungroup()
  
  # crossing time:
  #  - high-tail: require any *rise* after trough; first time >= threshold
  #  - low-tail : require any *fall* after peak;   first time <= threshold
  eps       <- 1e-8
  step_size <- 0.01
  
  cross_tbl <- df0 %>%
    dplyr::left_join(anchor_tbl, by = "Patient") %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(
      anchor_time = dplyr::first(anchor_time),
      anchor_val  = dplyr::first(anchor_val),
      
      # did the trajectory move in the expected direction after the anchor?
      move_ok = if (is_low) {
        any(Weeks_Since_Baseline > anchor_time &
              Value < (anchor_val - step_size), na.rm = TRUE)
      } else {
        any(Weeks_Since_Baseline > anchor_time &
              Value > (anchor_val + step_size), na.rm = TRUE)
      },
      
      # first time it crosses the appropriate side of the threshold
      cross_time = {
        if (is.na(anchor_time) || !move_ok) NA_real_ else {
          t_after <- if (is_low) {
            Weeks_Since_Baseline[Weeks_Since_Baseline > anchor_time &
                                   Value <= (thr_val + eps)]
          } else {
            Weeks_Since_Baseline[Weeks_Since_Baseline > anchor_time &
                                   Value >= (thr_val - eps)]
          }
          if (length(t_after)) min(t_after) else NA_real_
        }
      },
      .groups = "drop"
    )
  
  # ---- PLOT-FRIENDLY COPIES (cap applied here) ----
  df_plot <- df0 %>%
    dplyr::left_join(cross_tbl, by = "Patient") %>%
    dplyr::mutate(
      after_cross_pt = !is.na(cross_time) & Weeks_Since_Baseline >= cross_time,
      Value_plot     = if (is.finite(cap)) pmin(pmax(Value, 0), cap) else Value,
      overcap        = is.finite(cap) & Value > cap,
      # label: use % for percent metrics, otherwise numeric
      label_val = dplyr::case_when(
        is_percent & Value < 1 ~ sprintf("%.1f%%", Value * 100),
        is_percent & Value >= 1 ~ sprintf("%.0f%%", Value * 100),
        !is_percent & Value < 1 ~ sprintf("%.2f", Value),
        TRUE ~ as.character(round(Value))
      ),
      tp_shape = dplyr::case_when(
        prog_tp ~ "Progression",
        TRUE    ~ "Routine follow-up"
      )
    ) %>%
    dplyr::filter(!is.na(Weeks_Since_Baseline), !is.na(Value_plot))
  
  seg_plot <- df0 %>%
    dplyr::left_join(cross_tbl, by = "Patient") %>%
    dplyr::arrange(Patient, Weeks_Since_Baseline) %>%
    dplyr::group_by(Patient) %>%
    dplyr::mutate(
      x    = Weeks_Since_Baseline,
      y    = Value,
      xend = dplyr::lead(Weeks_Since_Baseline),
      yend = dplyr::lead(Value),
      y_plot    = if (is.finite(cap)) pmin(pmax(y,    0), cap) else y,
      yend_plot = if (is.finite(cap)) pmin(pmax(yend, 0), cap) else yend,
      after_cross_seg = !is.na(cross_time) &
        (xend >= cross_time - pad_before_weeks)  # color crossing seg too
    ) %>%
    dplyr::filter(!is.na(x), !is.na(xend), !is.na(y_plot), !is.na(yend_plot)) %>%
    dplyr::ungroup()
  
  mk_y_breaks <- function(cap) {
    if (!is.finite(cap)) return(waiver())
    if (cap >= 1000) {
      c(0, 500, 1000, 1500, 2000, cap)
    } else if (cap >= 10) {
      unique(c(0, pretty(c(0, cap), n = 6), cap))
    } else if (cap >= 1) {
      unique(c(0, pretty(c(0, cap), n = 6), cap))
    } else {
      # very small caps (e.g., 0.25): a few nicely spaced decimals
      br <- pretty(c(0, cap), n = 5)
      unique(c(0, br[br > 0 & br < cap], cap))
    }
  }
  
  mk_y_labels <- function(cap) {
    if (!is.finite(cap)) return(waiver())
    function(x) {
      if (cap < 1) {
        labs    <- sprintf("%.2f", x)   # two decimals for tiny ranges
        cap_lab <- sprintf("%.2f", cap)
      } else if (cap < 10) {
        labs    <- sprintf("%.1f", x)   # one decimal for mid ranges
        cap_lab <- sprintf("%.1f", cap)
      } else {
        # integers for large caps (no thousands separator, no .0)
        labs    <- scales::number_format(accuracy = 1, big.mark = "", decimal.mark = ".")(x)
        cap_lab <- formatC(cap, format = "f", digits = 0)  # -> "2500"
      }
      ix <- which(abs(x - cap) < 1e-8)
      if (length(ix)) labs[ix] <- paste0(">", cap_lab)
      labs
    }
  }
  
  top_expand <- function(cap) {
    if (!is.finite(cap)) return(0.05)
    if (cap < 1) 0.05 else 0.05   # a bit more headroom for tiny caps
  }
  
  # ---- PLOT ----
  ggplot(df_plot, aes(Weeks_Since_Baseline, Value_plot, group = Patient)) +
    geom_segment(
      data = seg_plot,
      aes(x = x, y = y_plot, xend = xend, yend = yend_plot, colour = after_cross_seg),
      linewidth = 0.55, alpha = 0.85
    ) +
    geom_point(
      aes(colour = after_cross_pt, shape = tp_shape, fill = after_cross_pt),
      size = 2, alpha = 0.9, stroke = 0.8
    ) +
    # labels for capped points (print true value just above the triangle)
    geom_text(
      data = dplyr::filter(df_plot, overcap),
      aes(label = label_val, colour = after_cross_pt),
      nudge_y = if (is.finite(cap)) cap * 0.03 else 0.03,   # ~3% headroom
      size = 2.4, show.legend = FALSE
    ) +
    facet_wrap(~ patient_relapse180, nrow = 1,
               labeller = labeller(patient_relapse180 = c(`FALSE`="No relapse ≤180d",
                                                          `TRUE` ="Relapse ≤180d"))) +
    labs(
      title = if (!is.null(plot_title)) plot_title
      else paste0(metric, " threshold from non-progressors: ", signif(thr_val, 3)),
      x = "Weeks Since Baseline", y = ylab
    ) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA),
          plot.margin = margin(t = 10, r = 8, b = 8, l = 8)) +
    # color = state relative to threshold
    scale_color_manual(
      values = c(`FALSE` = col_before, `TRUE` = col_after),
      labels = c(`FALSE` = "Below cfDNA threshold", `TRUE` = "Elevated cfDNA signal"),
      name   = "State"
    ) +
    # shape = timepoint type (triangle shows capped)
    scale_shape_manual(
      values = c("Routine follow-up" = 16, "Progression" = 23, "Capped (>cap)" = 17),
      name   = "Timepoint"
    ) +
    scale_fill_manual(values = c(`FALSE` = col_before, `TRUE` = col_after), guide = "none") +
    # if Mean.Coverage, flip y-axis and set custom limits
    # add y-scale/coord conditionally as a list of layers
    # ---- y-axis scaling ----
  (
    if (metric == "Mean.Coverage") {
      list(
        scale_y_reverse()
      )
    } else if (is_percent) {
      list(
        scale_y_continuous(
          labels = scales::percent_format(accuracy = 1),
          expand = expansion(mult = c(0.04, top_expand(cap)))
        ),
        coord_cartesian(ylim = if (is.finite(cap)) c(0, cap) else NULL, clip = "off"),
        labs(y = paste0(ylab, " (%)"))
      )
    } else {
      list(
        scale_y_continuous(
          breaks = mk_y_breaks(cap),
          labels = mk_y_labels(cap),
          expand = expansion(mult = c(0.04, top_expand(cap)))
        ),
        coord_cartesian(ylim = if (is.finite(cap)) c(0, cap) else NULL, clip = "off")
      )
    }
  ) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(shape = 16)),
      shape  = guide_legend(order = 2)
    )
  
}



# Example calls (Z-score panels: no cap)
p_sites  <- make_panel_minlogic("sites_z", "Proportion of Mutant Sites Detected (Z)", data = plot_df_pat, plot_title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles")
p_cvaf  <- make_panel_minlogic("cVAF_z",  "Cumulative VAF (Z)", cap = 2500,  data = plot_df_pat, plot_title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles")    # cap on
p_sites
p_cvaf


# Combine side-by-side with a single legend at the bottom
p_cvaf <- p_cvaf + ggtitle(NULL)
p_sites <- p_sites + ggtitle(NULL)

p_combined <- (p_cvaf + p_sites) +
  plot_annotation(
    title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles"
  ) +
  plot_layout(guides = "collect") &                           # <-- collect legends
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) &
  guides(                                                     # optional: one-row legends
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  )

p_combined

# Save
ggsave(
  filename = file.path(outdir, "BM_Zscore_Longitudinal_Monitoring4.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)

## Redo for non z-score BM features
# Example calls (Z-score panels: no cap)
p_sites  <- make_panel_minlogic("sites", "Proportion of Mutant Sites Detected", data = plot_df_pat, plot_title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles")
p_cvaf  <- make_panel_minlogic("cVAF",  "Cumulative VAF", cap = 0.25, data = plot_df_pat, plot_title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles")    # cap on
p_sites
p_cvaf


# Combine side-by-side with a single legend at the bottom
p_cvaf <- p_cvaf + ggtitle(NULL)
p_sites <- p_sites + ggtitle(NULL)

p_combined <- (p_cvaf + p_sites) +
  plot_annotation(
    title = "Longitudinal trajectories of MRD metrics from baseline BM mutation profiles"
  ) +
  plot_layout(guides = "collect") &                           # <-- collect legends
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) &
  guides(                                                     # optional: one-row legends
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  )

p_combined

ggsave(
  filename = file.path(outdir, "BM_raw_features_Longitudinal_Monitoring3.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)


## Do for blood and fragmentomic features
## Create unified plotting and threshold df to use with function 
## Ensure to flip the regulatory coverage fragment feature
## Just turn everything below 1 into a percent potentially 

# ------------------------------
# Helper: sample-level relapse flag
# ------------------------------
relapse_flag <- function(days) {
  dplyr::if_else(days <= 180, TRUE, FALSE, missing = FALSE)
}

# ------------------------------
# 1) BM mutation metrics (keep original names)
# ------------------------------
plot_df_bm <- dat %>%
  transmute(
    Patient,
    Weeks_Since_Baseline,
    Num_days_to_closest_relapse,
    cVAF   = detect_rate_BM,
    cVAF_z = z_score_detection_rate_BM,
    sites  = sites_rate_BM,
    sites_z = zscore_BM,
    relapse_within_180 = relapse_flag(Num_days_to_closest_relapse),
    Panel = "BM"
  ) %>%
  pivot_longer(
    cols = c(cVAF, cVAF_z, sites, sites_z),
    names_to = "Metric", values_to = "Value"
  ) %>%
  drop_na(Value)

# ------------------------------
# 2) Blood mutation metrics (rename so they’re distinct)
# ------------------------------
plot_df_blood <- dat %>%
  transmute(
    Patient,
    Weeks_Since_Baseline,
    Num_days_to_closest_relapse,
    cVAF_blood   = detect_rate_blood,
    cVAF_z_blood = z_score_detection_rate_blood,
    sites_blood  = sites_rate_blood,
    sites_z_blood = zscore_blood,
    relapse_within_180 = relapse_flag(Num_days_to_closest_relapse),
    Panel = "Blood"
  ) %>%
  pivot_longer(
    cols = c(cVAF_blood, cVAF_z_blood, sites_blood, sites_z_blood),
    names_to = "Metric", values_to = "Value"
  ) %>%
  drop_na(Value)

# ------------------------------
# 3) Fragmentomics metrics
#     - flip Mean.Coverage (inverse relationship)
#     - rename ichorCNA TF for a cleaner Metric name
# ------------------------------
plot_df_fragmentomics <- dat %>%
  transmute(
    Patient,
    Weeks_Since_Baseline,
    Num_days_to_closest_relapse,
    FS,
    Mean.Coverage,
    Proportion.Short,
    TF_ichorCNA = WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    relapse_within_180 = relapse_flag(Num_days_to_closest_relapse),
    Panel = "Fragmentomics"
  ) %>%
  pivot_longer(
    cols = c(FS, Mean.Coverage, Proportion.Short, TF_ichorCNA),
    names_to = "Metric", values_to = "Value"
  ) %>%
  drop_na(Value)

# ------------------------------
# 4) Combine all panels
# ------------------------------
plot_df_all <- bind_rows(plot_df_bm, plot_df_blood, plot_df_fragmentomics) %>%
  # optional: add patient-level relapse factor for faceting
  group_by(Patient) %>%
  mutate(patient_relapse180 = any(relapse_within_180, na.rm = TRUE)) %>%
  ungroup()

plot_df_all <- plot_df_all %>%
  mutate(prog_tp = !is.na(Num_days_to_closest_relapse) &
           Num_days_to_closest_relapse == 0)

## Now get thresholds 
# which metrics should use the *lower* tail?
low_tail_metrics <- c("Mean.Coverage")   # add others here if needed

# learn thresholds from NON-progressors, follow-up only
thr_tbl <- plot_df_all %>%
  dplyr::filter(as.character(patient_relapse180) == "FALSE",
                Weeks_Since_Baseline > 52) %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarise(
    thr = quantile(
      Value,
      # 0.10 for low-tail metrics; 0.90 for the rest
      ifelse(dplyr::first(Metric) %in% low_tail_metrics, 1 - q_flat, q_flat),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

thr_named <- setNames(thr_tbl$thr, thr_tbl$Metric)


# ------------------------------
# 5) Labels for plotting
# ------------------------------
custom_labels <- c(
  # BM mutation metrics
  cVAF    = "Cumulative VAF",
  cVAF_z  = "Cumulative VAF Z-score",
  sites   = "Proportion of Sites Detected",
  sites_z = "Proportion of Sites Detected Z-score",
  
  # Blood mutation metrics (renamed)
  cVAF_blood    = "Cumulative VAF",
  cVAF_z_blood  = "Cumulative VAF Z-score",
  sites_blood   = "Proportion of Sites Detected",
  sites_z_blood = "Proportion of Sites Detected Z-score",
  
  # Fragmentomics
  FS                        = "Fragment-size score",
  Mean.Coverage.Flipped     = "cfDNA coverage at MM active regulatory sites",
  Proportion.Short          = "Short-fragment proportion",
  TF_ichorCNA               = "cfDNA tumor fraction (ichorCNA)"
)

# optional: lock Metric factor order to labels order
#plot_df_all <- plot_df_all %>%
#  mutate(Metric = factor(Metric, levels = names(custom_labels)))


## Now plot 
## Blood features
# Raw calls 
p_sites  <- make_panel_minlogic("sites_blood", "Proportion of Mutant Sites Detected", data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")
p_cvaf  <- make_panel_minlogic("cVAF_blood",  "Cumulative VAF", cap = 0.22,  data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")    # cap on
p_sites
p_cvaf


# Combine side-by-side with a single legend at the bottom
p_cvaf <- p_cvaf + ggtitle(NULL)
p_sites <- p_sites + ggtitle(NULL)

p_combined <- (p_cvaf + p_sites) +
  plot_annotation(
    title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles"
  ) +
  plot_layout(guides = "collect") &                           # <-- collect legends
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) &
  guides(                                                     # optional: one-row legends
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  )

p_combined

ggsave(
  filename = file.path(outdir, "Blood_raw_features_Longitudinal_Monitoring3.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)


# Example calls (Z-score panels: no cap)
p_sites  <- make_panel_minlogic("sites_z_blood", "Proportion of Mutant Sites Detected (Z)", data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")
p_cvaf  <- make_panel_minlogic("cVAF_z_blood",  "Cumulative VAF (Z)", cap = 160,  data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")    # cap on
p_sites
p_cvaf


# Combine side-by-side with a single legend at the bottom
p_cvaf <- p_cvaf + ggtitle(NULL)
p_sites <- p_sites + ggtitle(NULL)

p_combined <- (p_cvaf + p_sites) +
  plot_annotation(
    title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles"
  ) +
  plot_layout(guides = "collect") &                           # <-- collect legends
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) &
  guides(                                                     # optional: one-row legends
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  )

p_combined

ggsave(
  filename = file.path(outdir, "Blood_zscore_features_Longitudinal_Monitoring3.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)


## Fragmentomic features
# Example calls (Z-score panels: no cap)
p_frag2 <- make_panel_minlogic("TF_ichorCNA",  "cfDNA Tumor Fraction", cap = NA,  data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")    # cap on
p_frag1  <- make_panel_minlogic("Proportion.Short", "Short cfDNA Fragments", data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")
p_frag1
p_frag2

# Combine side-by-side with a single legend at the bottom
p_frag1 <- p_frag1 + ggtitle(NULL)
p_frag2 <- p_frag2 + ggtitle(NULL)

p_combined <- (p_frag1 + p_frag2) +
  plot_annotation(
    title = "Longitudinal trajectories of fragmentomic features"
  ) +
  plot_layout(guides = "collect") &                           # <-- collect legends
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) &
  guides(                                                     # optional: one-row legends
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  )

p_combined

ggsave(
  filename = file.path(outdir, "Fragmentation_Features_1_Longitudinal_Monitoring2.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)



## Other features
p_frag1 <- make_panel_minlogic("FS",  "Fragment-size score", cap = NA,  data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")    # cap on
p_frag2  <- make_panel_minlogic("Mean.Coverage", "Mean cfDNA coverage (MM regs)", cap = NA, data = plot_df_all, plot_title = "Longitudinal trajectories of MRD metrics from baseline cfDNA mutation profiles")
p_frag1
p_frag2

# Combine side-by-side with a single legend at the bottom
p_frag1 <- p_frag1 + ggtitle(NULL)
p_frag2 <- p_frag2 + ggtitle(NULL)

p_combined <- (p_frag1 + p_frag2) +
  plot_annotation(
    title = "Longitudinal trajectories of fragmentomic features"
  ) +
  plot_layout(guides = "collect") &                           # <-- collect legends
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) &
  guides(                                                     # optional: one-row legends
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  )

p_combined

ggsave(
  filename = file.path(outdir, "Fragmentation_Features_2_Longitudinal_Monitoring2.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)






############################################################
# Section: Progressor vs Non-Progressor Example
# Purpose: Add illustrative comparison of representative cases
# Notes: Highlight differences in cfDNA or MRD trajectories
############################################################


#### Add specific example of progressor vs non-progressor 
### Need the full healthy control info for this 
Merged_MRDetect_zscore <- readRDS("MRDetect_output_winter_2025/Processed_R_outputs/cfWGS_Winter2025All_MRDetect_with_Zscore_May2025.rds")

## -----------------------------
## USER PARAMS
## -----------------------------
pid        <- "CA-08"
to_percent <- TRUE   # change to FALSE if Cumulative_VAF_BM already in % units
hc_timept  <- c("01")  # acceptable baseline codes in healthy file


## -----------------------------
## 1.  Patient CA-08 longitudinal cVAF (from `dat`)
## -----------------------------
ca08_df <- dat %>%
  filter(Patient == pid) %>%
  select(Patient, Date, timepoint_info, Cumulative_VAF_BM) %>%
  filter(!is.na(timepoint_info),
         !is.na(Cumulative_VAF_BM)) %>%
  arrange(Date) %>%
  # Ensure unique timepoint_info order (chronological)
  mutate(timepoint_info = factor(timepoint_info, levels = unique(timepoint_info))) %>%
  # convert to % if required
  mutate(
    cVAF = if (to_percent) Cumulative_VAF_BM * 100 else Cumulative_VAF_BM
  ) %>%
  mutate(plot_group = "Patient CA-08")


## capture ordered levels for later re-use (healthy facet)
tp_levels <- levels(ca08_df$timepoint_info)


## -----------------------------
## 2.  Healthy-control background values
##     (detection_rate_as_reads_detected_over_reads_checked)
## -----------------------------
hc_vals <- Merged_MRDetect_zscore %>%
  filter(
    Mut_source    == "BM_cells",
    Filter_source == "STR_encode",
    Study         == "CHARM_healthy",
    Patient       == pid,
    Sample_type   == "BM_cells",
    Timepoint %in% hc_timept
  ) %>%
  pull(detection_rate_as_reads_detected_over_reads_checked)

if (length(hc_vals) == 0) {
  warning("No healthy-control rows found for CA-08 under supplied filters.")
}

# Convert to % if the patient values were scaled
if (to_percent) hc_vals <- hc_vals * 100

# Repeat the healthy-control distribution at *each* CA-08 timepoint
# so the two facets have aligned x-axes & equal panel widths.
hc_df <- tibble(
  timepoint_info = factor(rep(tp_levels, each = length(hc_vals)), levels = tp_levels),
  cVAF           = rep(hc_vals, times = length(tp_levels)),
  plot_group     = "Healthy controls"
)


## -----------------------------
## 3.  Bind and plot
## -----------------------------
# 1) Compute healthy‑control median (in same units as ca08_df$cVAF)
median_hc <- median(hc_df$cVAF, na.rm = TRUE)

# 2) Define custom, prettier x‑axis labels
x_labels <- c(
  "Diagnosis"         = "Diagnosis",
  "Post_induction"    = "Post‑induction",
  "1yr maintenance"   = "1‑year maintenance",
  "1.5yr maintenance" = "1.5‑year maintenance",
  "Relapse"           = "Relapse",
  "Healthy_controls"  = "Healthy controls"
)

# 3) Re‑factor your patient df to include the Healthy_controls level
#    (so that tick appears, but no red points/lines show there)
all_levels <- c(levels(ca08_df$timepoint_info), "Healthy_controls")
ca08_df2 <- ca08_df %>%
  mutate(timepoint_info = factor(timepoint_info, levels = all_levels))

# 4) Plot
p_single3 <- ggplot(ca08_df, aes(x = timepoint_info, y = cVAF, group = 1)) +
  # patient trajectory
  geom_line(color = "red", size = 0.6) +
  geom_point(color = "red", size = 2) +
  # healthy‐control median as a dotted line
  geom_hline(yintercept = median_hc,
             linetype     = "dotted",
             colour       = "black",
             size         = 0.5) +
  # healthy controls at the new "Healthy_controls" tick
  geom_point(data = hc_df,
             aes(x = timepoint_info, y = cVAF),
             colour = "black", size = 1.2, alpha = 0.6,
             position = position_jitter(width = 0.1, height = 0)
  ) +
  # nicer x labels
  scale_x_discrete(labels = x_labels) +
  # y‐axis percent (or raw) as before
  scale_y_continuous(
    name   = if (to_percent) "Tumour‑informed cumulative VAF (%)"
    else "Tumour‑informed cumulative VAF",
    labels = if (to_percent) scales::label_number(accuracy = 0.1)
    else scales::waiver()
  ) +
  labs(
    x     = NULL,
    title = paste0("Tumour‑informed longitudinal tracking: Progressor")
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.x   = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(outdir, "Fig2A_CA08_cVAF_vsHealthy_median_dotted.png"),
  p_single3, width = 5, height = 3, dpi = 600
)




## Add new patient ID 
# Load id_map first
id_map <- readRDS("id_map.rds")  


### As function
## From this 2A and B is used
library(rlang)
make_cVAF_hc_plot <- function(pid,
                              dat,
                              merged,
                              hc_timept    = c("01"),
                              to_percent   = TRUE,
                              tp_order     = NULL,
                              show_hc_points = TRUE,
                              hc_filters   = list(Mut_source="BM_cells",
                                                  Filter_source="STR_encode",
                                                  Study="CHARM_healthy",
                                                  Sample_type="BM_cells"),
                              return_data  = FALSE,
                              progressor_color    = "red",
                              nonprogressor_color = "black") {
  
  # ------------------ Patient longitudinal data ------------------
  pat_df <- dat %>%
    filter(.data$Patient == pid) %>%
    select(Patient, Date, timepoint_info, Cumulative_VAF_BM,
           Num_days_to_closest_relapse = any_of("Num_days_to_closest_relapse")) %>%
    filter(!is.na(timepoint_info), !is.na(Cumulative_VAF_BM)) %>%
    arrange(Date)
  
  if (nrow(pat_df) == 0) {
    abort(paste0("No data found in `dat` for patient ", pid, "."))
  }
  
  # Determine relapse / progressor status (≤180 days at ANY timepoint)
  if (!"Num_days_to_closest_relapse" %in% names(pat_df)) {
    pat_is_prog <- NA
  } else {
    pat_is_prog <- any(pat_df$Num_days_to_closest_relapse <= 180, na.rm = TRUE)
  }
  
  # Choose colour & title suffix
  pat_col <- if (isTRUE(pat_is_prog)) progressor_color else nonprogressor_color
  title_suffix <- if (isTRUE(pat_is_prog)) "Progressor" else "Non‑progressor"
  
  # x-axis ordering
  if (is.null(tp_order)) {
    tp_levels <- unique(pat_df$timepoint_info)  # chronological because pat_df arranged by Date
  } else {
    # use supplied order but keep only levels that appear
    tp_levels <- intersect(tp_order, unique(pat_df$timepoint_info))
    # ensure we didn't lose any because of naming mismatches
    missing_in_order <- setdiff(unique(pat_df$timepoint_info), tp_levels)
    if (length(missing_in_order)) {
      tp_levels <- c(tp_levels, missing_in_order)
    }
  }
  
  # add Healthy_controls at the end
  tp_levels_all <- c(tp_levels, "Healthy_controls")
  
  pat_df <- pat_df %>%
    mutate(
      cVAF = if (to_percent) Cumulative_VAF_BM * 100 else Cumulative_VAF_BM,
      timepoint_info = factor(timepoint_info, levels = tp_levels_all)
    )
  
  # ------------------ Healthy-control values ------------------
  # build filter expression from hc_filters
  mf <- merged
  if ("Mut_source"    %in% names(hc_filters)) mf <- mf %>% filter(.data$Mut_source    == hc_filters$Mut_source)
  if ("Filter_source" %in% names(hc_filters)) mf <- mf %>% filter(.data$Filter_source == hc_filters$Filter_source)
  if ("Study"         %in% names(hc_filters)) mf <- mf %>% filter(.data$Study         == hc_filters$Study)
  if ("Sample_type"   %in% names(hc_filters)) mf <- mf %>% filter(.data$Sample_type   == hc_filters$Sample_type)
  
  # patient match + baseline timepoint(s)
  mf <- mf %>%
    filter(.data$Patient == pid,
           .data$timepoint_info %in% hc_timept)
  
  hc_vals <- mf$detection_rate_as_reads_detected_over_reads_checked
  
  if (length(hc_vals) == 0) {
    warning("No healthy-control rows found for patient ", pid, 
            " under current filters; using NA median.")
  }
  
  if (to_percent) hc_vals <- hc_vals * 100
  
  median_hc <- median(hc_vals, na.rm = TRUE)
  
  # Data frame for plotting raw HC points (if requested)
  hc_df <- tibble(
    timepoint_info = factor("Healthy_controls", levels = tp_levels_all),
    cVAF           = hc_vals
  )
  
  # ------------------ label mapping ------------------
  # Use global object tp_label_map if exists; else auto-camel-case
  if (exists("tp_label_map", inherits = TRUE)) {
    lbl_map <- get("tp_label_map", inherits = TRUE)
  } else {
    lbl_map <- setNames(gsub("_", " ", tp_levels_all), tp_levels_all)
  }
  # ensure we have label for Healthy_controls
  if (!"Healthy_controls" %in% names(lbl_map)) {
    lbl_map["Healthy_controls"] <- "Healthy controls"
  }
  
  # ------------------ Plot ------------------
  new_pid <- id_map$New_ID[match(pid, id_map$Patient)]
  if (is.na(new_pid)) new_pid <- pid  # fallback if not found
  
  p <- ggplot(pat_df, aes(x = timepoint_info, y = cVAF, group = 1)) +
    geom_line(color = pat_col, linewidth = 0.6) +
    geom_point(color = pat_col, size = 2) +
    # dotted median (only if median not NA)
    { if (!is.na(median_hc)) 
      geom_hline(yintercept = median_hc,
                 linetype = "dotted", colour = "black", linewidth = 0.5) 
      else NULL } +
    # HC jittered points at right tick (optional)
    { if (show_hc_points && length(hc_vals)) 
      geom_point(data = hc_df,
                 aes(x = timepoint_info, y = cVAF),
                 inherit.aes = FALSE,
                 colour = "black", size = 1.2, alpha = 0.6,
                 position = position_jitter(width = 0.1, height = 0)) 
      else NULL } +
    scale_x_discrete(labels = lbl_map) +
    scale_y_continuous(
      name   = if (to_percent) "Tumour‑informed cumulative VAF (%)" else "Tumour‑informed cumulative VAF",
      labels = if (to_percent) label_number(accuracy = 0.1) else waiver()
    ) +
    labs(
      x     = NULL,
      title = paste0("Tumour‑informed longitudinal tracking: ", new_pid, " (", title_suffix, ")")
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title  = element_text(size = 10.5, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  
  if (return_data) {
    out <- list(plot = p,
                patient_df = pat_df,
                hc_df = hc_df,
                median_hc = median_hc,
                is_progressor = pat_is_prog)
    return(out)
  } else {
    return(p)
  }
}

 
all_pids <- cohort_df %>%
  filter(Cohort == "Frontline") %>%
  pull(Patient) %>%
  unique() %>%
  sort()


# Optionally save each:
for(pid in all_pids){
  # Attempt to make the plot; if it errors, skip
  p <- tryCatch(
    make_cVAF_hc_plot(pid, dat, Merged_MRDetect_zscore,
                      hc_timept = "Diagnosis", to_percent = TRUE),
    error = function(e){
      message("Skipping ", pid, ": ", e$message)
      return(NULL)
    }
  )
  
  # If the plot was successfully created, save it
  if(!is.null(p)){
    ggsave(
      filename = file.path(outdir, paste0(pid, "_cVAF_vsHealthy_updated.png")),
      plot     = p,
      width    = 5,
      height   = 3,
      dpi      = 600
    )
  }
}
















#### Below here is testing and old code no longer needed 

## Now work on blood muts
# 1. Pivot to long form
plot_df_blood <- dat %>%
  select(
    Patient,
    Weeks_Since_Baseline,
    Num_days_to_closest_relapse,
    cVAF       = detect_rate_blood,
    cVAF_z     = z_score_detection_rate_blood,
    sites      = sites_rate_blood,
    sites_z    = zscore_blood
  ) %>%
  mutate(
    relapse_within_180 = if_else(
      Num_days_to_closest_relapse <= 180,
      TRUE,
      FALSE,
      missing = FALSE    # recode any NA here to FALSE since didn't relapse
    )) %>%
  pivot_longer(
    cols = c(cVAF, cVAF_z, sites, sites_z),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  drop_na(Value)

# 2. Plot
p_traj <- ggplot(plot_df_blood, aes(x = Weeks_Since_Baseline, y = Value, group = Patient)) +
  geom_line(colour = "grey70", alpha = 0.4) +
  geom_point(aes(color = relapse_within_180), size = 1.5, alpha = 0.8) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol  = 4,
    labeller = labeller(Metric = custom_labels)
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180 d"),
    name   = "Relapse status",
    na.value = "black",         # <- this makes NA points black
    na.translate = FALSE      # ← drop NA from the legend since these patients didn't progress
    
  ) +
  labs(
    x     = "Weeks Since Baseline",
    y     = "Value",
    title = "Longitudinal trajectories of MRD metrics from baseline PB cfDNA mutation profiles"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

# 2. Save as PNG in outdir
ggsave(
  filename = file.path(outdir, "Fig3A_metrics_trajectories_Blood.png"),
  plot     = p_traj,
  device   = "png",
  width    = 12,    # inches
  height   = 4,     # inches
  dpi      = 600
)


### Add other versions 
### Edit style 
# ─────────────────────────────────────────────────────────────
# Prep: derive a patient-level relapse flag (any TP ≤180d)
# ─────────────────────────────────────────────────────────────
plot_df_pat <- plot_df_blood %>%
  group_by(Patient) %>%
  mutate(
    patient_relapse180 = any(replace_na(relapse_within_180, FALSE))
  ) %>%
  ungroup() %>%
  mutate(
    # make it an ordered factor for consistent legend / facet order
    patient_relapse180 = factor(patient_relapse180,
                                levels = c(FALSE, TRUE))
  )

# nice labels for the patient-level factor
patient_relapse_labs <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d (last TP)"
)
patient_relapse_labs_short <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d"
)



# ===================================================================
# VERSION 1
# Lines coloured by *patient-level* relapse (any TP ≤180d)
# (All points inherit the same colour as the line.)
# ===================================================================
p_traj_patientline <- ggplot(
  plot_df_pat,
  aes(x = Weeks_Since_Baseline, y = Value, group = Patient,
      colour = patient_relapse180)
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(
    ~ Metric,
    scales   = "free_y",
    ncol     = 4,
    labeller = labeller(Metric = custom_labels)
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = patient_relapse_labs,
    name   = "Patient relapse status"
  ) +
  labs(
    x     = "Weeks Since Baseline",
    y     = "Value",
    title = "Longitudinal trajectories of MRD metrics from baseline PB cfDNA mutation profiles"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

ggsave(
  filename = file.path(outdir, "Fig3B_metrics_trajectories_blood_byPatientLine.png"),
  plot     = p_traj_patientline,
  device   = "png",
  width    = 12,
  height   = 4,
  dpi      = 600
)



# ===================================================================
# VERSION 2
# Facet by patient-level relapse status *and* by Metric.
# ===================================================================
patient_relapse_labs_short <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d"
)

# 1) Prep: patient-level flag
plot_df_pat <- plot_df_blood %>%
  group_by(Patient) %>%
  mutate(patient_relapse180 = any(replace_na(relapse_within_180, FALSE))) %>%
  ungroup() %>%
  mutate(
    patient_relapse180 = factor(patient_relapse180, levels = c(FALSE, TRUE)),
    relapse_within_180 = replace_na(relapse_within_180, FALSE)
  )

# 2) Build the two mini‑plots

## A) Cumulative VAF
df_cvaf <- filter(plot_df_pat, Metric == "cVAF")
df_cvaf_seg <- df_cvaf %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_cvaf <- ggplot(df_cvaf, aes(Weeks_Since_Baseline, Value, group = Patient)) +
#  geom_line(color = "black", alpha = 0.6) +
  # coloured segments
  geom_segment(
    data = df_cvaf_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["cVAF"],
    x     = "Weeks Since Baseline",
    y     = "Cumulative VAF"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

## B) Proportion of Sites Detected
df_sites <- filter(plot_df_pat, Metric == "sites")
df_sites_seg <- df_sites %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()


p_sites <- ggplot(df_sites, aes(Weeks_Since_Baseline, Value, group = Patient)) +
  #geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_sites_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["sites"],
    x     = "Weeks Since Baseline",
    y     = "Prop. Mutant Sites Detected"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

# 3) Combine them side by side
p_combined <- p_cvaf + p_sites + 
  plot_annotation(
    title = "Longitudinal trajectories of MRD metrics from baseline PB cfDNA mutation profiles",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  ) & 
  theme(legend.position = "bottom")

# 4) Save
ggsave(
  filename = file.path(outdir, "Fig3B_sideBySide_lockedY_blood_segmented2.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)


# 2) Build the two mini‑plots

## A) Cumulative VAF
df_cvaf <- filter(plot_df_pat, Metric == "cVAF_z")
df_cvaf_seg <- df_cvaf %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_cvaf <- ggplot(df_cvaf, aes(Weeks_Since_Baseline, Value, group = Patient)) +
#  geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_cvaf_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["cVAF_z"],
    x     = "Weeks Since Baseline",
    y     = "Cumulative VAF (Z)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

## B) Proportion of Sites Detected
df_sites <- filter(plot_df_pat, Metric == "sites_z")
df_sites_seg <- df_sites %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()


p_sites <- ggplot(df_sites, aes(Weeks_Since_Baseline, Value, group = Patient)) +
#  geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_sites_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["sites_z"],
    x     = "Weeks Since Baseline",
    y     = "Prop. Mutant Sites Detected (Z)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

# 3) Combine them side by side
p_combined <- p_cvaf + p_sites + 
  plot_annotation(
    title = "Longitudinal trajectories of MRD metrics from baseline PB cfDNA mutation profiles",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  ) & 
  theme(legend.position = "bottom")


# 4) Save
ggsave(
  filename = file.path(outdir, "Fig3B_sideBySide_lockedY_zscore_blood_segmented2.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)






## Now for fragmentomics 
# 1. Pivot to long form
plot_df_fragmentomics <- dat %>%
  select(
    Patient,
    Weeks_Since_Baseline,
    Num_days_to_closest_relapse,
    FS,
    Mean.Coverage,
    Proportion.Short,
    WGS_Tumor_Fraction_Blood_plasma_cfDNA
  ) %>%
  mutate(
    relapse_within_180 = if_else(
      Num_days_to_closest_relapse <= 180,
      TRUE,
      FALSE,
      missing = FALSE    # recode any NA here to FALSE since didn't relapse
    )) %>%
  pivot_longer(
    cols = c(FS, Mean.Coverage, Proportion.Short, WGS_Tumor_Fraction_Blood_plasma_cfDNA),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  drop_na(Value)

### Flip the mean.coverage since has an inverse relationship 
# 1) after your pivot_longer/drop_na, add a new y‐column:
plot_df_fragmentomics2 <- plot_df_fragmentomics %>%
  mutate(
    # flip only Mean.Coverage by negating it
    Value = if_else(Metric == "Mean.Coverage",
                     -Value,
                     Value))


# 2. Plot
# 1. Build the plot and save to a variable
custom_labels <- c(
  FS                             = "Fragment-size score",
  Mean.Coverage                  = "cfDNA coverage at MM active regulatory sites",
  Mean.Coverage.Flipped                  = "cfDNA coverage at MM active regulatory sites",
  Proportion.Short               = "Short-fragment proportion",
  WGS_Tumor_Fraction_Blood_plasma_cfDNA = "cfDNA tumor fraction (ichorCNA)"
)

p_traj <- ggplot(plot_df_fragmentomics, aes(x = Weeks_Since_Baseline, y = Value, group = Patient)) +
  geom_line(colour = "grey70", alpha = 0.4) +
  geom_point(aes(color = relapse_within_180), size = 1.5, alpha = 0.8) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol  = 4,
    labeller = labeller(Metric = custom_labels)
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180 d"),
    name   = "Relapse status",
    na.value = "black",         # <- this makes NA points black
    na.translate = FALSE      # ← drop NA from the legend since these patients didn't progress
    
  ) +
  labs(
    x     = "Weeks Since Baseline",
    y     = "Value",
    title = "Longitudinal trajectories of fragmentomic features"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

# 2. Save as PNG in outdir
ggsave(
  filename = file.path(outdir, "Fig3A_metrics_trajectories_fragmentomics.png"),
  plot     = p_traj,
  device   = "png",
  width    = 12,    # inches
  height   = 4,     # inches
  dpi      = 600
)

# Flipped version
p_traj <- ggplot(plot_df_fragmentomics2, aes(x = Weeks_Since_Baseline, y = Value, group = Patient)) +
  geom_line(colour = "grey70", alpha = 0.4) +
  geom_point(aes(color = relapse_within_180), size = 1.5, alpha = 0.8) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol  = 4,
    labeller = labeller(Metric = custom_labels)
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180 d"),
    name   = "Relapse status",
    na.value = "black",         # <- this makes NA points black
    na.translate = FALSE      # ← drop NA from the legend since these patients didn't progress
    
  ) +
  labs(
    x     = "Weeks Since Baseline",
    y     = "Value",
    title = "Longitudinal trajectories of fragmentomic features"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

# 2. Save as PNG in outdir
ggsave(
  filename = file.path(outdir, "Fig3C_metrics_trajectories_fragmentomics_flipped.png"),
  plot     = p_traj,
  device   = "png",
  width    = 12,    # inches
  height   = 4,     # inches
  dpi      = 600
)


### Now other versions 
### Edit style 
# ─────────────────────────────────────────────────────────────
# Prep: derive a patient-level relapse flag (any TP ≤180d)
# ─────────────────────────────────────────────────────────────
plot_df_pat <- plot_df_fragmentomics2 %>%
  group_by(Patient) %>%
  mutate(
    patient_relapse180 = any(replace_na(relapse_within_180, FALSE))
  ) %>%
  ungroup() %>%
  mutate(
    # make it an ordered factor for consistent legend / facet order
    patient_relapse180 = factor(patient_relapse180,
                                levels = c(FALSE, TRUE))
  )

# nice labels for the patient-level factor
patient_relapse_labs <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d (last TP)"
)
patient_relapse_labs_short <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d"
)



# ===================================================================
# VERSION 1
# Lines coloured by *patient-level* relapse (any TP ≤180d)
# (All points inherit the same colour as the line.)
# ===================================================================
p_traj_patientline <- ggplot(
  plot_df_pat,
  aes(x = Weeks_Since_Baseline, y = Value, group = Patient,
      colour = patient_relapse180)
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(
    ~ Metric,
    scales   = "free_y",
    ncol     = 4,
    labeller = labeller(Metric = custom_labels)
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = patient_relapse_labs,
    name   = "Patient relapse status"
  ) +
  labs(
    x     = "Weeks Since Baseline",
    y     = "Value",
    title = "Longitudinal trajectories of fragmentomic features"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

ggsave(
  filename = file.path(outdir, "Fig3C_metrics_trajectories_blood_byPatientLine.png"),
  plot     = p_traj_patientline,
  device   = "png",
  width    = 12,
  height   = 4,
  dpi      = 600
)



# ===================================================================
# VERSION 2
# Facet by patient-level relapse status *and* by Metric.
# ===================================================================
patient_relapse_labs_short <- c(
  `FALSE` = "No relapse ≤180d",
  `TRUE`  = "Relapse ≤180d"
)

# 1) Prep: patient-level flag
plot_df_pat <- plot_df_fragmentomics2 %>%
  group_by(Patient) %>%
  mutate(patient_relapse180 = any(replace_na(relapse_within_180, FALSE))) %>%
  ungroup() %>%
  mutate(
    patient_relapse180 = factor(patient_relapse180, levels = c(FALSE, TRUE)),
    relapse_within_180 = replace_na(relapse_within_180, FALSE)
  )

# 2) Build the two mini‑plots

## A) Cumulative VAF
df_cvaf <- filter(plot_df_pat, Metric == "FS")
df_cvaf_seg <- df_cvaf %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_cvaf <- ggplot(df_cvaf, aes(Weeks_Since_Baseline, Value, group = Patient)) +
  # geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_cvaf_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["FS"],
    x     = "Weeks Since Baseline",
    y     = "Fragment-size score"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

## B) Mean coverage (flipped)
df_sites <- filter(plot_df_pat, Metric == "Mean.Coverage")
df_sites_seg <- df_sites %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_sites <- ggplot(df_sites, aes(Weeks_Since_Baseline, -Value, group = Patient)) +
 # geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_sites_seg,
    aes(
      x      = x,
      y      = -y,
      xend   = xend,
      yend   = -yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  scale_y_reverse() +    
  labs(
    title = custom_labels["Mean.Coverage"],
    x     = "Weeks Since Baseline",
    y     = "Mean cfDNA coverage (MM regs)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

# 3) Combine them side by side
p_combined <- p_cvaf + p_sites + 
  plot_annotation(
    title = "Longitudinal trajectories of fragmentomic features",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  ) & 
  theme(legend.position = "bottom")

# 4) Save
ggsave(
  filename = file.path(outdir, "Fig3C_sideBySide_lockedY_fragmentomics_segmented2.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)


# 2) Now for next

## A) Cumulative VAF
df_cvaf <- filter(plot_df_pat, Metric == "Proportion.Short")
df_cvaf_seg <- df_cvaf %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p_cvaf <- ggplot(df_cvaf, aes(Weeks_Since_Baseline, Value, group = Patient)) +
  # geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_cvaf_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["Proportion.Short"],
    x     = "Weeks Since Baseline",
    y     = "Short cfDNA Fragments (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

## B) Proportion of Sites Detected
df_sites <- filter(plot_df_pat, Metric == "WGS_Tumor_Fraction_Blood_plasma_cfDNA")
df_sites_seg <- df_sites %>%
  arrange(Patient, Weeks_Since_Baseline) %>%
  group_by(Patient) %>%
  mutate(
    x    = Weeks_Since_Baseline,
    y    = Value,
    xend = lead(Weeks_Since_Baseline),
    yend = lead(Value),
    # segment should turn red if its *end* point is a relapse
    seg_relapse = lead(relapse_within_180, default = FALSE)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()


p_sites <- ggplot(df_sites, aes(Weeks_Since_Baseline, Value, group = Patient)) +
 #  geom_line(color = "black", alpha = 0.6) +
  geom_segment(
    data = df_sites_seg,
    aes(
      x      = x,   y    = y,
      xend   = xend, yend = yend,
      colour = seg_relapse
    ),
    size  = 0.4,
    alpha = 0.6
  ) +
  geom_point(aes(color = relapse_within_180),
             size = 1.8, alpha = 0.8) +
  facet_wrap(~ patient_relapse180, nrow = 1,
             labeller = labeller(patient_relapse180 = patient_relapse_labs_short)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    labels = c(`FALSE` = "No relapse ≤180d", `TRUE` = "Relapse ≤180d"),
    name   = NULL
  ) +
  labs(
    title = custom_labels["WGS_Tumor_Fraction_Blood_plasma_cfDNA"],
    x     = "Weeks Since Baseline",
    y     = "cfDNA Tumor Fraction"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))

# 3) Combine them side by side
p_combined <- p_cvaf + p_sites + 
  plot_annotation(
    title = "Longitudinal trajectories of fragmentomic features",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  ) & 
  theme(legend.position = "bottom")


# 4) Save
ggsave(
  filename = file.path(outdir, "Fig3B_sideBySide_lockedY_2_fragmentomics_segmented2.png"),
  plot     = p_combined,
  width    = 12,   # 4 panels across
  height   = 4,
  dpi      = 600
)










### Now correlation with lab values 
# 1. Build a data frame of all relevant columns
# Make slightly smaller to avoid being overcrowded
pair_df <- dat_clean %>%
  select(
    Patient,
    zscore_blood,                   # sites z-score (blood)
    z_score_detection_rate_blood,   # cVAF z-score (blood)
    detect_rate_blood,              # raw cVAF (%) (blood)
    zscore_BM,                      # sites z-score (BM)
    z_score_detection_rate_BM,      # cVAF z-score (BM)
    detect_rate_BM,                 # raw cVAF (%) (BM)
    FS,                             # fragment-size score
    Mean.Coverage,                  # mean regulatory coverage (×)
    Proportion.Short,               # short-fragment proportion (%)
    WGS_Tumor_Fraction_Blood_plasma_cfDNA,  # ichorCNA tumor fraction (%)
    M_Protein,                      # M-protein level
    Calcium,                        # serum calcium
    B2_micro,                       # β2-microglobulin
    Kappa_Lambda_Ratio, 
    Flow_pct_cells,                 # % aberrant plasma cells by MFC
    Adaptive_Frequency              # clonoSEQ dominant rearrangement frequency
  ) 

# 2. Rename for Plot-Friendly Labels
pair_df <- pair_df %>%
  rename(
    `Blood cVAF z-score`          = z_score_detection_rate_blood,
    `Blood sites z-score`            = zscore_blood,
    `Blood cVAF (%)`              = detect_rate_blood,
    `BM sites z-score`            = zscore_BM,
    `BM cVAF z-score`             = z_score_detection_rate_BM,
    `BM cVAF (%)`                 = detect_rate_BM,
    `Fragment-size score`         = FS,
    `Regulatory coverage (×)`     = Mean.Coverage,
    `Short-fragment proportion` = Proportion.Short,
    `Tumor fraction (%)`          = WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    `M-protein (g/L)`             = M_Protein,
    `Calcium (mmol/L)`            = Calcium,
    `β2-microglobulin (mg/L)`     = B2_micro,
    `Kappa:Lambda ratio`          = Kappa_Lambda_Ratio,
    `% aberrant plasma cells (MFC)`     = Flow_pct_cells,
    `clonoSEQ frequency`          = Adaptive_Frequency
  )

# 3. Generate the GGally pairs plot
wrap_labels <- function(labels, width = 20) {
  sapply(labels, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

wrapped_labels <- wrap_labels(names(pair_df), width = 15)


pair_df <- pair_df %>% select(-Patient)
p_pairs <- ggpairs(
  pair_df,
  columns    = 1:ncol(pair_df),
  upper      = list(continuous = wrap("cor", size = 2.5, method = "spearman", use = "pairwise.complete.obs")),
  lower      = list(continuous = wrap("points", alpha = 0.4, size = 0.8)),
  diag       = list(continuous = wrap("densityDiag")),
  columnLabels = wrapped_labels
) +
  theme_classic(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text  = element_text(size = 9)
  )

# 4. Save as high-res PNG
ggsave(
  filename = file.path(outdir, "Supplementary_Figure_pairs_metrics_clinical_3B_2.png"),
  plot     = p_pairs,
  width    = 18,    # inches
  height   = 18,    # inches
  dpi      = 600
)





#### Now make heatmap 
# 1. Compute Spearman correlation matrix
corr_mat <- cor(pair_df, method = "spearman", use = "pairwise.complete.obs")

# 2. Turn into long format for ggplot
corr_df_updated <- as.data.frame(as.table(corr_mat)) %>%
  set_names(c("Metric1", "Metric2", "rho")) %>%
  # ensure the ordering matches the matrix
  mutate(
    Metric1 = factor(Metric1, levels = colnames(corr_mat)),
    Metric2 = factor(Metric2, levels = colnames(corr_mat))
  )

# 3. (Optional) order by hierarchical clustering for prettier blocks
ord <- hclust(dist(corr_mat))$order
levs <- colnames(corr_mat)[ord]
corr_df_updated <- corr_df_updated %>%
  mutate(
    Metric1 = factor(Metric1, levels = levs),
    Metric2 = factor(Metric2, levels = levs)
  )

# 4. Plot heatmap
p_heatmap <- ggplot(corr_df_updated, aes(x = Metric1, y = Metric2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  scale_fill_viridis_c(option = "D", name = expression(rho~"(Spearman)")) +
  coord_equal() +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.grid      = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Spearman correlation heatmap of cfDNA features and biomarkers"
  )

# 5. Save it
ggsave(
  file.path(outdir, "Fig_heatmap_spearman_updated.png"),
  plot   = p_heatmap,
  width  = 6,
  height = 5,
  dpi    = 600
)

# 6. Keep only the upper‐triangle (including diagonal)
corr_df_tri <- corr_df_updated %>%
  mutate(
    r = as.integer(Metric1),
    c = as.integer(Metric2)
  ) %>%
  filter(r < c) %>%     # strictly upper, drops diagonal
  select(-r, -c)

# 7. Plot only the upper triangle
p_heatmap_tri <- ggplot(corr_df_tri, aes(x = Metric1, y = Metric2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2) +
  scale_fill_viridis_c(option = "D", name = expression(rho~"(Spearman)")) +
  coord_equal() +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.grid      = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Spearman correlation heatmap"
  )

p_heatmap_tri <- p_heatmap_tri +
  theme(
    # place legend at x=0.85 (85% from left), y=0.15 (15% from bottom)
    legend.position     = c(0.85, 0.15),
    # anchor the legend box’s top-right corner at that point
    legend.justification = c(1, 0),
    # give it a semi-opaque background so tiles underneath don’t show through
    legend.background   = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size     = unit(0.8, "lines"),
    legend.text         = element_text(size = 6),
    legend.title        = element_text(size = 7)
  )

p_heatmap_tri <- ggplot(corr_df_tri, aes(x = Metric1, y = Metric2, fill = rho)) +
  geom_tile(color = "white") +
  
  # draw text in white when rho < -0.3, else black
  geom_text(aes(label = sprintf("%.2f", rho), color = rho < -0.3),
            size = 2) +
  scale_color_manual(
    values = c(`TRUE` = "white", `FALSE` = "black"),
    guide  = FALSE
  ) +
  
  # force legend from -1 to +1, with breaks including 1
  scale_fill_viridis_c(
    option = "D",
    name   = expression(rho~"(Spearman)"),
    limits = c(-0.6, 1),
    breaks = c(-0.5, 0, 0.5, 1),
    labels = c("-0.5", "0.0", "0.5", "1.0")
  ) +
  
  coord_equal() +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y     = element_text(size = 6),
    panel.grid      = element_blank(),
    plot.title      = element_text(size = 8, face = "bold"),
    legend.position = c(0.9, 0.2),
    legend.justification = c(1, 0),
    legend.background   = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size     = unit(0.6, "lines"),
    legend.text         = element_text(size = 6),
    legend.title        = element_text(size = 7)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Spearman correlation heatmap"
  )



# 8. Save the triangular heatmap
### Figure 3B
ggsave(
  filename = file.path(outdir, "Fig_heatmap_spearman_upper_triangle2.png"),
  plot     = p_heatmap_tri,
  width    = 6,
  height   = 5,
  dpi      = 600
)








################################################################################
