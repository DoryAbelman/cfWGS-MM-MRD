################################################################################
## longitudinal_MRD_features_v2.R
##  – NA-robust, column-robust, start-date aware                        (2025-06)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(rstatix)
  library(patchwork)
  library(ggpubr)
  library(glue)
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
write_csv(pairwise_stats, file.path(outdir, "placeholder_stats_v2.csv"))


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
paired_plot <- function(var, ylab) {
  
  bl <- paste0(var, "_baseline")
  f1 <- paste0(var, "_follow1")
  if (!all(c(bl,f1) %in% names(paired_master))) return(NULL)
  
  paired_master %>%
    select(Patient, all_of(c(bl,f1))) %>%
    pivot_longer(-Patient, names_to="Time", values_to="value") %>%
    drop_na() %>%
    mutate(Time = factor(Time,
                         levels=c(bl,f1),
                         labels=c("Baseline","Post-Tx#1"))) %>%
    ggplot(aes(Time, value, group=Patient))+
    geom_line(alpha=.4)+
    geom_violin(aes(fill=Time), width=.9, alpha=.4, colour=NA)+
    geom_boxplot(width=.25, outlier.shape=NA)+
    labs(y=ylab, x=NULL)+
    theme_bw(base_size=11)+
    theme(legend.position="none")
}

p1 <- paired_plot("zscore_blood", "Blood mutation z-score")
p2 <- paired_plot("detect_rate_blood", "Blood detection rate")
p3 <- paired_plot("FS", "Fragment-size score")

plot_combined <- (p1 | p2 | p3) + plot_annotation(tag_levels = "A")

ggsave(file.path(outdir, "Fig_longitudinal_pairs_v2.pdf"),
       plot = plot_combined,
       width = 10, height = 4, dpi = 600, useDingbats = FALSE)

## Longitudinal spaghetti – keep rows with at least one var not-NA
traj_vars <- intersect(c("zscore_blood","z_score_detection_rate_blood",
                         "detect_rate_blood","FS",
                         "Mean.Coverage","Proportion.Short",
                         "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
                       names(dat))

traj_long <- dat %>%
  select(Patient, Date, all_of(traj_vars)) %>%
  pivot_longer(-c(Patient,Date), names_to="Feature", values_to="value") %>%
  drop_na()

ggplot(traj_long, aes(Date,value,group=Patient,colour=Patient))+
  geom_line(alpha=.6)+
  facet_wrap(~Feature, scales="free_y", ncol=3)+
  theme_bw(base_size=10)+
  theme(legend.position="none")+
  labs(x=NULL,y=NULL) %>%
  ggsave("Fig_longitudinal_spaghetti_v2.pdf",
         width=8, height=10, dpi=600, useDingbats=FALSE)

## Correlation figure (only when both columns exist)
corr_plot <- function(xvar,yvar,xlab,ylab){
  if(!(xvar %in% names(dat) && yvar %in% names(dat))) return(NULL)
  dat %>% select(all_of(c(xvar,yvar))) %>% drop_na() %>%
    ggplot(aes(.data[[xvar]], .data[[yvar]]))+
    geom_point(alpha=.7, size=1.5)+
    geom_smooth(method="lm",se=FALSE, linetype="dashed")+
    stat_cor(method="spearman", label.x.npc="left", label.y.npc=.9, size=3)+
    theme_bw(base_size=11)+
    labs(x=xlab,y=ylab)
}

c1 <- corr_plot("zscore_blood","z_score_detection_rate_blood",
                "Blood mutation z-score","Blood detection-rate z-score")
c2 <- corr_plot("WGS_Tumor_Fraction_Blood_plasma_cfDNA","Blood_Mutation_Count",
                "Tumour fraction","Blood mutation count")
c3 <- corr_plot("WGS_Tumor_Fraction_Blood_plasma_cfDNA","FS",
                "Tumour fraction","Fragment-size score")

((c1|c2|c3)+plot_annotation(tag_levels="A")) %>%
  ggsave("Fig_correlations_v2.pdf",
         width=10, height=3.5, dpi=600, useDingbats=FALSE)

## ───── 8. Console read-out for quick copy-paste ─────────────────────────────
cat("\n──── Study-level counts ────\n")
cat(glue::glue("cfDNA draws = {n_cfDNA_draws}\n",
               "BM aspirates = {n_BM_aspirates}\n",
               "Patients = {n_patients}\n",
               "Median follow-up = {median_follow$med} months ",
               "(range {median_follow$lo}–{median_follow$hi})\n\n"))

print(pairwise_stats, width=Inf)
print(corr_pairs, width=Inf)

sessionInfo()
################################################################################
