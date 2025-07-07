# ===============================================================================
# Script: longitudinal_MRD_feature_amalusis.R
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
  library(viridis)
})

## ───── 1. Load data ──────────────────────────────────────────────────────────
### Set paths 
outdir   <- "Output_tables_2025"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

### Load data 
file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.rds")
cohort_df <- readRDS("cohort_assignment_table_updated.rds")

dat <- file 

# 1.  Join cohort_df and keep frontline only -------------------------------------
dat <- dat %>%                # <‑‑ your master data
  left_join(cohort_df, by = "Patient") 

dat <- dat %>% mutate(Date = as.Date(Date))   # safety

## Filter to baseline for simplicity and ones have data on  
dat <- dat %>% filter(Cohort == "Frontline")

dat <- dat %>%
  filter(!is.na(zscore_BM) | !is.na(zscore_blood) | !is.na(FS))

## ───── 2. Identify BASELINE and FIRST-FOLLOW rows (start-date logic) ────────
baseline_tbl <- dat %>%
  group_by(Patient) %>%
  filter(
    timepoint_info %in% c("Diagnosis", "Baseline") |
      Date == min(Date, na.rm = TRUE)
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(Patient, Date,
         everything())                  # keep all cols - easier later

follow1_tbl <- dat %>%
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
n_cfDNA_draws  <- dat %>% filter(!is.na(detect_rate_blood)) %>% nrow()
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
write_csv(pairwise_stats, file.path(outdir, "placeholder_stats_v3.csv"))


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

# ── 2. Define core and extra variable names -----------------------------------
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
write_csv(pairwise_stats, file.path(outdir, "paired_stats_summary.csv"))
write_csv(corr_table, file.path(outdir, "mrd_metric_lab_all_correlations.csv"))
write_csv(sig_corrs, file.path(outdir, "mrd_metric_lab_sig_correlations.csv"))







## ───── 7. Figures – all filtered on !is.na for the plotted feature ──────────
## Helper for paired violin/box with spaghetti
## 9.1  Utilities ------------------------------------------------------------
baseline_dates <- readRDS("Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS.rds")

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
    cVAF       = Cumulative_VAF_BM,
    cVAF_z     = z_score_detection_rate_BM,
    sites      = detect_rate_BM,
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


## Now for blood 
# 1. Pivot to long form
plot_df_blood <- dat %>%
  select(
    Patient,
    Weeks_Since_Baseline,
    Num_days_to_closest_relapse,
    cVAF       = Cumulative_VAF_blood,
    cVAF_z     = z_score_detection_rate_blood,
    sites      = detect_rate_blood,
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
# 1. Build the plot and save to a variable
custom_labels <- c(
  cVAF    = "Cumulative VAF",
  cVAF_z  = "Cumulative VAF Z-score",
  sites   = "Proportion of Sites Detected",
  sites_z = "Proportion of Sites Detected Z-score"
)

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

# 2. Plot
# 1. Build the plot and save to a variable
custom_labels <- c(
  FS                             = "Fragment-size score",
  Mean.Coverage                  = "cfDNA coverage at MM active regulatory sites",
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
  filename = file.path(outdir, "Supplementary_Figure_pairs_metrics_clinical_3B.png"),
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
  file.path(outdir, "Fig_heatmap_spearman.png"),
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
  filename = file.path(outdir, "Fig_heatmap_spearman_upper_triangle.png"),
  plot     = p_heatmap_tri,
  width    = 6,
  height   = 5,
  dpi      = 600
)








################################################################################
