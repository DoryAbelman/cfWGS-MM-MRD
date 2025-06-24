################################################################################
##  Detection-Rate & Progression Analysis
##  cfWGS MRD manuscript – Dory A.
################################################################################

## ── 0.  SETUP ────────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(tableone)      # optional – baseline table
library(timeROC)       # optional – AUC vs time

## INPUT (already saved from your pipeline) ------------------------------
final_tbl_rds <- "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS.rds"
dat_rds       <- "output_tables_2025/all_patients_with_BM_and_blood_calls_updated2.rds"

## OUTPUT ----------------------------------------------------------------------
outdir <- "Output_tables_2025/detection_progression"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ── 1.  LOAD & TIDY CORE TABLES ──────────────────────────────────────────────
final_tbl <- readRDS(final_tbl_rds) %>%
  transmute(
    Patient       = as.character(Patient),
    baseline_date = as.Date(Baseline_Date),
    censor_date   = as.Date(Censor_date),
    relapsed      = as.integer(Relapsed)
  )

dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )

## Limit to frontline 
dat <- dat %>% filter(Cohort == "Frontline")

## ── 2.  BUILD PFS TABLE (patient-level) ──────────────────────────────────────
survival_df <- dat %>%
  mutate(
    Patient     = as.character(Patient),
    sample_date = as.Date(Date)                     # the draw date for each row
  ) %>%
  # bring in censor_date + relapsed from final_tbl
  left_join(
    final_tbl %>% 
      select(Patient, censor_date, relapsed),
    by = "Patient"
  ) %>%
  # compute days from the sample to the event/censor
  mutate(
    Time_to_event   = as.numeric(censor_date - sample_date),
    Relapsed_Binary = as.integer(relapsed)
  ) %>%
  # keep only the columns your KM‐loop needs
  select(
    Patient, Timepoint, sample_date, censor_date, timepoint_info,
    Time_to_event, Relapsed_Binary,
    Flow_Binary, Adaptive_Binary, Rapid_Novor_Binary,
    Flow_pct_cells, Adaptive_Frequency,
    PET_Binary,
    BM_zscore_only_detection_rate_call, BM_zscore_only_detection_rate_prob,
    Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob
  )

# sanity checks
table(survival_df$Relapsed_Binary, useNA="ifany")
summary(survival_df$Time_to_event)
table(survival_df$timepoint_info)


#### get function ready 
# 2) Friendly tech names
techs <- c(
  Flow_Binary        = "MFC",
  Adaptive_Binary    = "clonoSEQ",
  BM_zscore_only_detection_rate_call    = "cfWGS_BM", 
  Blood_zscore_only_sites_call = "cfWGS_PB_cfDNA"
)

# 3) All timepoints to cover
tps <- unique(survival_df$timepoint_info)
dpi_target <- 500

# 4) Minimum n per group to plot
min_n <- 5

# 5) Loop
for(tp in tps) {
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, ".png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%           # keep just the first draw per patient if multiple at that timepoint
      ungroup() %>%
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    df_sub <- df_sub %>% 
      mutate(Time_to_event = Time_to_event/30.44) # divide days by 30.44 to get months
    
    # skip if too few pts
    if(nrow(df_sub) < min_n) next
    
    # skip if only one group present
    if(n_distinct(df_sub$Group) < 2) next
    
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = TRUE,
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = TRUE,
      risk.table      = TRUE,
      palette         = c("#E7B800","#2E9FDF"),
      legend.title    = paste0(assay_lab, " MRD"),
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Months since sample",
      ylab            = "PFS",
      title           = paste0(assay_lab, " MRD at ", tp),
      risk.table.height = 0.25
    )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3.5,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 8, 
      height   = 9,
      dpi      = dpi_target
    )
  }
}


### To do 
### Check that days since sample ok 
## See why goes up to 1 at start 
## Add other analyses







#### Now get stats for results 
### First on frontline 
# 2.  Define frontline cohort -------------------------------------------------
front_patients <- dat %>%
  filter(Cohort == "Frontline") %>%
  distinct(Patient)

# 3.  Median follow-up & relapse rate -----------------------------------------
pfs_front <- final_tbl %>%
  filter(Patient %in% front_patients$Patient) %>%
  # compute time from baseline to censor/relapse
  mutate(time_days = as.numeric(censor_date - baseline_date)) 

## check one row per patient 
pfs_front %>% 
  count(Patient) %>% 
  filter(n > 1) -> dups
if(nrow(dups)) stop("Duplicate patients found: ", paste(dups$Patient, collapse = ", "))


median_fu_mo <- median(pfs_front$time_days / 30.44, na.rm = TRUE)
n_front      <- nrow(pfs_front)
n_rel        <- sum(pfs_front$relapsed)
pct_rel      <- n_rel / n_front * 100

message(glue::glue(
  "Median follow-up: {round(median_fu_mo,1)} months\n",
  "Relapses: {n_rel}/{n_front} ({round(pct_rel)}%) front-line patients"
))

## More info 
# 1) Basic summary statistics for follow-up time (in days and months)
followup_stats <- pfs_front %>%
  summarise(
    N_patients      = dplyr::n(),
    N_relapses      = sum(relapsed),
    Relapse_rate    = N_relapses / N_patients * 100,
    
    min_days        = min(time_days, na.rm = TRUE),
    q1_days         = quantile(time_days, 0.25, na.rm = TRUE),
    median_days     = median(time_days, na.rm = TRUE),
    q3_days         = quantile(time_days, 0.75, na.rm = TRUE),
    max_days        = max(time_days, na.rm = TRUE),
    mean_days       = mean(time_days, na.rm = TRUE),
    sd_days         = sd(time_days, na.rm = TRUE),
    
    min_months      = min(time_days, na.rm = TRUE) / 30.44,
    q1_months       = quantile(time_days / 30.44, 0.25, na.rm = TRUE),
    median_months   = median(time_days / 30.44, na.rm = TRUE),
    q3_months       = quantile(time_days / 30.44, 0.75, na.rm = TRUE),
    max_months      = max(time_days, na.rm = TRUE) / 30.44,
    mean_months     = mean(time_days, na.rm = TRUE) / 30.44,
    sd_months       = sd(time_days, na.rm = TRUE) / 30.44
  )

print(followup_stats)

# 2) Optional: distribution of follow-up times
followup_hist <- pfs_front %>%
  mutate(followup_months = time_days / 30.44) %>%
  ggplot(aes(x = followup_months)) +
  geom_histogram(binwidth = 3, boundary = 0) +
  labs(
    x = "Follow-up time (months)",
    y = "Number of patients",
    title = "Distribution of follow-up times in frontline cohort"
  )

# If you want to export the stats to CSV
write_csv(followup_stats, file.path(outdir, "frontline_followup_summary.csv"))

# 4.  Assays & timepoint definitions ------------------------------------------
assays <- c(
#  EasyM  = "Rapid_Novor_Binary",
  clonoSEQ = "Adaptive_Binary",
  Flow     = "Flow_Binary",
  cfWGS_BM    = "BM_zscore_only_detection_rate_call",
  cfWGS_Blood    = "Blood_zscore_only_sites_call"
)

post_labels   <- c("post_transplant")
one_year_labels <- c("1yr maintenance")

compute_sens <- function(df, col) {
  df2      <- df %>% filter(!is.na(.data[[col]]))
  n_tested <- nrow(df2)
  n_pos    <- sum(df2[[col]] == 1, na.rm = TRUE)
  tibble(
    N_tested      = n_tested,
    N_positive    = n_pos,
    Sensitivity   = n_pos / n_tested
  )
}

# 5.  Post-ASCT sensitivities -----------------------------------------------
post_df <- dat %>%
  filter(
    Patient        %in% front_patients$Patient,
    str_detect(timepoint_info, paste(post_labels, collapse = "|"))
  ) %>%
  arrange(Patient, sample_date) %>%
  group_by(Patient) %>%
  slice(1) %>%   # earliest post-ASCT sample
  ungroup() %>%
  select(Patient, one_of(assays)) %>%
  left_join(final_tbl %>% select(Patient, relapsed), by = "Patient") %>%
  filter(relapsed == 1)

post_stats <- map_dfr(names(assays), ~ {
  col <- assays[.x]
  stats <- compute_sens(post_df, col)
  stats %>% mutate(Assay = .x)
}, .id = NULL) %>%
  select(Assay, everything())

# 6.  One-year sensitivities -----------------------------------------------
year_df <- dat %>%
  filter(
    Patient        %in% front_patients$Patient,
    str_detect(timepoint_info, paste(one_year_labels, collapse = "|"))
  ) %>%
  arrange(Patient, sample_date) %>%
  group_by(Patient) %>%
  slice(1) %>%   # earliest 1-yr maintenance sample
  ungroup() %>%
  select(Patient, one_of(assays)) %>%
  left_join(final_tbl %>% select(Patient, relapsed), by = "Patient") %>%
  filter(relapsed == 1)

year_stats <- map_dfr(names(assays), ~ {
  col <- assays[.x]
  stats <- compute_sens(year_df, col)
  stats %>% mutate(Assay = .x)
}, .id = NULL) %>%
  select(Assay, everything())

# 7.  Print tables ------------------------------------------------------------
message("Post-ASCT sensitivity among relapsers:")
print(post_stats)

message("1-year sensitivity among relapsers:")
print(year_stats)

# 8.  (Optional) write out results -------------------------------------------
write_csv(
  post_stats,
  file.path(outdir, "frontline_postASCT_sensitivity.csv")
)
write_csv(
  year_stats,
  file.path(outdir, "frontline_1yr_sensitivity.csv")
)


#### Now do seperately only amongst patients who got a cfWGS test - seperate for BM and blood 
# —— after you’ve built `post_df` and `year_df` as before… ——————————————

# 5a.  Head-to-head in the BM-cfWGS subset (only patients with BM Z-score) ———

# define the BM-cfWGS column name
bm_col <- assays["cfWGS_BM"]

post_df_BM <- post_df %>% filter(!is.na(.data[[bm_col]]))
year_df_BM <- year_df %>% filter(!is.na(.data[[bm_col]]))

post_stats_BM <- map_dfr(names(assays), function(a) {
  compute_sens(post_df_BM, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

year_stats_BM <- map_dfr(names(assays), function(a) {
  compute_sens(year_df_BM, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

message("Post-ASCT sensitivities among those with BM-cfWGS:")
print(post_stats_BM)

message("1-yr sensitivities among those with BM-cfWGS:")
print(year_stats_BM)

# 5b.  Head-to-head in the blood-cfWGS subset (only patients with blood Z-score) —

blood_col <- assays["cfWGS_Blood"]

post_df_blood <- post_df %>% filter(!is.na(.data[[blood_col]]))
year_df_blood <- year_df %>% filter(!is.na(.data[[blood_col]]))

post_stats_blood <- map_dfr(names(assays), function(a) {
  compute_sens(post_df_blood, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

year_stats_blood <- map_dfr(names(assays), function(a) {
  compute_sens(year_df_blood, assays[a]) %>% 
    mutate(Assay = a)
}) %>% select(Assay, everything())

message("Post-ASCT sensitivities among those with blood-cfWGS:")
print(post_stats_blood)

message("1-yr sensitivities among those with blood-cfWGS:")
print(year_stats_blood)

# 6.  (Optional) write them out ———————————————————————————————

write_csv(post_stats_BM,    file.path(outdir, "frontline_postASCT_sens_BMcfWGS.csv"))
write_csv(year_stats_BM,    file.path(outdir, "frontline_1yr_sens_BMcfWGS.csv"))
write_csv(post_stats_blood, file.path(outdir, "frontline_postASCT_sens_bloodcfWGS.csv"))
write_csv(year_stats_blood, file.path(outdir, "frontline_1yr_sens_bloodcfWGS.csv"))





#### Now get other results and do power analysis 

## 1. Subset to post-transplant & BM-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(BM_zscore_only_detection_rate_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# we already have `fit_fl` from before
sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "BM_zscore_only_detection_rate_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$BM_zscore_only_detection_rate_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## 6. Draft paragraph ----
paragraph <- glue(
  "After one year of maintenance therapy, BM-cfWGS MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC, MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), ",
  "comparable to flow cytometry (ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients, the minimum detectible HR ",
  "for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph)


## Get tibble 
# 7) Assemble a single summary table of all key metrics ---------------------

# gather into a single-row tibble
metrics_1yr <- tibble(
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100
)



### Redo for post-transplant 
## 1. Subset to post-transplant & BM-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(BM_zscore_only_detection_rate_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ BM_zscore_only_detection_rate_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# we already have `fit_fl` from before
sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "BM_zscore_only_detection_rate_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$BM_zscore_only_detection_rate_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## 6. Draft paragraph ----
paragraph <- glue(
  "At post-transplant, BM-cfWGS MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC, MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), ",
  "comparable to flow cytometry (ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients, the minimum detectible HR ",
  "for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph)

## Compile 
metrics_post_transplant <- tibble(
  Landmark        = "post_transplant",
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100
)

metrics_1yr <- metrics_1yr %>%
  mutate(Landmark = "1yr_maintenance") %>%
  select(Landmark, everything())  # move Landmark to front

progression_metrics <- bind_rows(
  metrics_post_transplant,
  metrics_1yr
)

# 8) Export summary table --------------------------------------------
write_csv(
  progression_metrics,
  file.path(outdir, "cfWGS_vs_flow_progression_summary.csv")
)

# (Optional) also save as RDS for later use
saveRDS(
  progression_metrics,
  file.path(outdir, "cfWGS_vs_flow_progression_summary.rds")
)





##### Now get the stats on cfWGS from blood derived muts as well 
#### Now get other results and do power analysis 

## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Blood_zscore_only_sites_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS BM (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# we already have `fit_fl` from before
sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "Blood_zscore_only_sites_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_zscore_only_sites_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## 6. Draft paragraph ----
paragraph <- glue(
  "After one year of maintenance therapy, Blood-cfWGS MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by Blood-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC, MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), ",
  "comparable to flow cytometry (ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients, the minimum detectible HR ",
  "for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph)


## Get tibble 
# 7) Assemble a single summary table of all key metrics ---------------------

# gather into a single-row tibble
metrics_1yr <- tibble(
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100
)



### Redo for post-transplant 
## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Blood_zscore_only_sites_call)
  )

## 2. 24-month RFS by cfWGS Blood ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_zscore_only_sites_call,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_cf      <- cox_cf$estimate
ci_lo_cf   <- cox_cf$conf.low
ci_hi_cf   <- cox_cf$conf.high

## 2a) Median RFS by cfWGS Blood (days → months) ----
med_cf <- surv_median(fit_cf)$median
med_neg_cf <- med_cf[1] / 30.44
med_pos_cf <- med_cf[2] / 30.44

## 2b) 24-month RFS by flow cytometry ----
# we already have `fit_fl` from before
sum_fl24   <- summary(fit_fl, times = t24)
rfs_neg_fl <- sum_fl24$surv[1] * 100
rfs_pos_fl <- sum_fl24$surv[2] * 100

## 3. Median RFS by flow ----
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "Blood_zscore_only_sites_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_zscore_only_sites_call==1)
zα     <- qnorm(1-0.05/2); zβ <- qnorm(0.80)
hr80   <- exp(2*(zα+zβ) / sqrt(d * prop_p * (1-prop_p)))

# power to detect HR = 2.0
## --- power to detect HR = 2 with Schoenfeld formula --------------------
hr_target <- 2
ln_hr     <- log(hr_target)              # ln HR
p         <- prop_p                      # proportion in MRD-positive group
z_alpha   <- qnorm(1 - 0.05/2)           # 1.96 for two-sided α = .05

z_beta    <- sqrt(d * p * (1 - p)) * ln_hr - z_alpha
pw2       <- pnorm(z_beta)               # ≈ power for HR = 2

## 6. Draft paragraph ----
paragraph <- glue(
  "At post-transplant, Blood-cfWGS MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC, MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), ",
  "comparable to flow cytometry (ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients, the minimum detectible HR ",
  "for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph)

## Compile 
metrics_post_transplant <- tibble(
  Landmark        = "post_transplant",
  RFS24_cf_neg   = rfs_neg_cf,
  RFS24_cf_pos   = rfs_pos_cf,
  MedRFS_cf_neg  = med_neg_cf,
  MedRFS_cf_pos  = med_pos_cf,
  RFS24_fl_neg   = rfs_neg_fl,
  RFS24_fl_pos   = rfs_pos_fl,
  MedRFS_fl_neg  = med_neg_fl,
  MedRFS_fl_pos  = med_pos_fl,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  Spearman_prob  = rho1,
  Spearman_flow  = rho2,
  Events         = d,
  Patients       = nrow(df_km),
  HR_80pct       = hr80,
  Power_HR2_pct  = pw2 * 100
)

metrics_1yr <- metrics_1yr %>%
  mutate(Landmark = "1yr_maintenance") %>%
  select(Landmark, everything())  # move Landmark to front

progression_metrics_blood <- bind_rows(
  metrics_post_transplant,
  metrics_1yr
)

# 8) Export summary table --------------------------------------------
write_csv(
  progression_metrics_blood,
  file.path(outdir, "cfWGS_vs_flow_progression_summary_blood_muts.csv")
)

# (Optional) also save as RDS for later use
saveRDS(
  progression_metrics_blood,
  file.path(outdir, "cfWGS_vs_flow_progression_summaryy_blood_muts.rds")
)









### Now evaluate on non-frontline too 
## Re-introdue dat since previously limited 
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )

################################################################################
##  Time-window prediction performance in Non-frontline cohort
################################################################################

# 1) Your assays vector
assays <- c(
  Flow         = "Flow_Binary",
  cfWGS_BM     = "BM_zscore_only_detection_rate_call",
  cfWGS_Blood  = "Blood_zscore_only_sites_call"
)

# 2) Build df_sf: non-frontline samples + assays + relapse info
df_sf <- dat %>%
  filter(tolower(Cohort) == "non-frontline") %>%
  transmute(
    Patient,
    sample_date = as.Date(Date),
    Adaptive_Binary,
    Flow_Binary,
    BM_zscore_only_detection_rate_call,
    Blood_zscore_only_sites_call
  ) %>%
  # keep rows with at least one assay result
  filter(if_any(all_of(assays), ~ !is.na(.x))) %>%
  # join in per‐patient relapse_date + relapsed flag
  left_join(
    final_tbl %>% 
      transmute(
        Patient,
        relapse_date = as.Date(censor_date),
        relapsed     = as.integer(relapsed)
      ),
    by = "Patient"
  )

# 2) Build the sample‐level df for non‐frontline -----------------------------
df_sf <- dat %>%
  filter(tolower(Cohort) == "non-frontline") %>%
  transmute(
    Patient,
    sample_date = as.Date(Date),
    Adaptive_Binary,
    Flow_Binary,
    BM_zscore_only_detection_rate_call,
    Blood_zscore_only_sites_call
  ) %>%
  # keep rows with at least one assay result
  filter(if_any(all_of(assays), ~ !is.na(.x))) %>%
  # join in per‐patient relapse_date + relapsed flag
  left_join(
    final_tbl %>% 
      transmute(
        Patient,
        relapse_date = as.Date(censor_date),
        relapsed     = as.integer(relapsed)
      ),
    by = "Patient"
  )

# 3) Define windows (days) to evaluate
windows <- c(90, 180, 365, 730)

# 4) Helper to compute metrics for one assay + one window
calc_metrics <- function(df, assay_label, col_name, win_d) {
  df2 <- df %>%
    filter(!is.na(.data[[col_name]])) %>%
    mutate(
      test_pos        = (.data[[col_name]] == 1),
      event_in_window = relapsed == 1 &
        !is.na(relapse_date) &
        relapse_date <= sample_date + days(win_d)
    )
  
  # force both TRUE/FALSE levels even if absent
  tab <- table(
    factor(df2$test_pos,        levels = c(FALSE, TRUE)),
    factor(df2$event_in_window, levels = c(FALSE, TRUE))
  )
  
  tp <- tab["TRUE","TRUE"]
  fn <- tab["FALSE","TRUE"]
  fp <- tab["TRUE","FALSE"]
  tn <- tab["FALSE","FALSE"]
  
  tibble(
    Window_days = win_d,
    Assay       = assay_label,
    N_samples   = nrow(df2),
    N_patients   = n_distinct(df2$Patient), 
    TP = tp, FN = fn, FP = fp, TN = tn,
    Sensitivity = if((tp+fn)>0) tp/(tp+fn) else NA_real_,
    Specificity = if((tn+fp)>0) tn/(tn+fp) else NA_real_,
    PPV         = if((tp+fp)>0) tp/(tp+fp) else NA_real_,
    NPV         = if((tn+fn)>0) tn/(tn+fn) else NA_real_
  )
}


# 6) Inspect results
results %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity)) %>%
  print(n = Inf)



#### Now do only for those who got the cfWGS test done 

# 5a) Restrict to BM-cfWGS subset, then loop ------------------------------
bm_col     <- assays["cfWGS_BM"]    # "BM_zscore_only_detection_rate_call"
df_sf_BM   <- df_sf %>% filter(!is.na(.data[[bm_col]]))

results_BM <- map_dfr(windows, function(w) {
  map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_BM, a, assays[[a]], w)
  })
}) %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity))

# 5b) Restrict to blood-cfWGS subset, then loop ---------------------------
blood_col   <- assays["cfWGS_Blood"]  # "Blood_zscore_only_sites_call"
df_sf_blood <- df_sf %>% filter(!is.na(.data[[blood_col]]))

results_blood <- map_dfr(windows, function(w) {
  map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_blood, a, assays[[a]], w)
  })
}) %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity))

# 6) Inspect both ------------------------------------------
message("=== BM-cfWGS subset ===")
print(results_BM)

message("=== Blood-cfWGS subset ===")
print(results_blood)

## Get event counts 
# define your windows of interest
windows <- c(90, 180, 365, 730)

# assuming you’ve already built df_sf (with Patient, sample_date, relapsed, relapse_date)
# this will count, for each window:
# - how many distinct patients relapsed within that window
# - how many samples fall into that window

event_counts_BM <- map_dfr(windows, function(w) {
  df_w <- df_sf_BM %>%
    filter(
      relapsed == 1,
      !is.na(relapse_date),
      relapse_date <= sample_date + days(w)
    )
  
  tibble(
    Window_days = w,
    N_patients  = n_distinct(df_w$Patient),
    N_samples   = nrow(df_w)
  )
})

# view the result
print(event_counts)


### Export this
# full results (all patients with any assay)
write_csv(
  results,
  file.path(outdir, "all_assays_timewindow_results.csv")
)

# BM‐cfWGS subset
write_csv(
  results_BM,
  file.path(outdir, "BM_cfWGS_timewindow_results.csv")
)

# blood‐cfWGS subset
write_csv(
  results_blood,
  file.path(outdir, "blood_cfWGS_timewindow_results.csv")
)

# event counts for BM‐cfWGS
write_csv(
  event_counts_BM,
  file.path(outdir, "BM_cfWGS_event_counts.csv")
)

# (if you also computed event_counts for the blood subset:)
event_counts_blood <- map_dfr(windows, function(w) {
  df_w <- df_sf_blood %>%
    filter(
      relapsed == 1,
      !is.na(relapse_date),
      relapse_date <= sample_date + days(w)
    )
  tibble(
    Window_days = w,
    N_patients  = n_distinct(df_w$Patient),
    N_samples   = nrow(df_w)
  )
})

write_csv(
  event_counts_blood,
  file.path(outdir, "blood_cfWGS_event_counts.csv")
)



### Maybe follow up with Esteban to see if clinical or biochemical progression since don't know for non-frontline
### Leave as blank for now
# Filter out patients whose ID starts with "IMG-" or "SPORE-"
relapse_filtered <- Relapse_dates_full %>%
  filter(!grepl("^IMG|^SPORE", Patient))

# Export to RDS
saveRDS(relapse_filtered, file = "Relapse_dates_full_filtered.rds")

# Export to CSV (optional)
write.csv(relapse_filtered, file = "Relapse_dates_full_filtered.csv", row.names = FALSE)








### Below here is testing 

## ── 5.  FUNCTIONS TO ANALYSE & PLOT ONE ASSAY ───────────────────────────────
plot_km <- function(df, assay_col,
                    palette = c("#1b9e77", "#d95f02")) {
  # 1) filter & create a simple marker factor --------------------------------
  df2 <- df %>%
    filter(!is.na(.data[[assay_col]])) %>%
    mutate(marker = case_when(
      .data[[assay_col]] %in% c("pos","positive",1,"1",TRUE)  ~ "MRD-positive",
      .data[[assay_col]] %in% c("neg","negative",0,"0",FALSE) ~ "MRD-negative",
      TRUE                                                      ~ NA_character_
    )) %>%
    filter(!is.na(marker)) %>%
    mutate(marker = factor(marker, levels = c("MRD-negative","MRD-positive")))
  
  if (nrow(df2) == 0) {
    stop("No data for assay: ", assay_col)
  }
  
  # 2) fit Kaplan–Meier ------------------------------------------------------
  fit <- survfit(Surv(time_mo, event) ~ marker, data = df2)
  
  # 3) extract log-rank p-value yourself ------------------------------------
  pval <- surv_pvalue(fit, data = df2)$pval
  
  # 4) plot (pval passed as label string) -----------------------------------
  ggs <- ggsurvplot(
    fit, data       = df2,
    risk.table      = TRUE,
    pval            = paste0("p=", signif(pval, 2)),
    conf.int        = FALSE,
    palette         = palette,
    xlab            = "Months from transplant",
    ylab            = "Relapse-free survival",
    legend.title    = NULL,
    legend.labs     = levels(df2$marker),
    title           = paste(assay_col, "(post-ASCT)")
  )
  
  # 5) univariable Cox (HR + 95% CI) -----------------------------------------
  cox_mod <- coxph(Surv(time_mo, event) ~ marker, data = df2)
  cox_tab <- broom::tidy(
    cox_mod,
    exponentiate = TRUE,
    conf.int     = TRUE
  ) %>%
    select(HR = estimate, CI_low = conf.low, CI_high = conf.high, p.value)
  
  list(plot = ggs, cox = cox_tab)
}



## list of assays to iterate over
assays <- c("BM_zscore_only_detection_rate_call",
            "Blood_zscore_only_sites_call",
            "Flow_Binary",
            "Adaptive_Binary")

# 1. Are the columns really in surv_df?
names(surv_df)
# 2. How many non-NA values per assay?
sapply(surv_df[assays], function(x) sum(!is.na(x)))


# pick only assays with data
non_na <- sapply(assays, function(a) sum(!is.na(surv_df[[a]])))
valid <- assays[non_na > 0]
if (length(valid) < length(assays)) {
  message("Dropping assays with no data: ", paste(setdiff(assays, valid), collapse = ", "))
}

# run
km_res <- map(valid, ~plot_km(surv_df, .x))

## export individual PDFs
walk2(km_res, valid, ~ggsave(
  file.path(outdir, paste0(.y, "_KM.pdf")),
  plot = .x$plot$plot, width = 5, height = 4
))


## combined figure (2×2 grid)
combined <- wrap_plots(map(km_res, "plot"), ncol = 2)
ggsave(file.path(outdir, "Figure_KM_4panel.pdf"),
       combined, width = 8, height = 6)

## HR table
hr_tbl <- map2_dfr(km_res, assays,
                   \(x, nm) mutate(x$cox, Assay = nm)) %>%
  select(Assay, HR, CI_low, CI_high, p.value)

write_csv(hr_tbl, file.path(outdir, "Cox_HR_table.csv"))

## ── 6.  DETECTION-RATE & LEAD-TIME STATISTICS ──────────────────────────────
## helper to compute sensitivity etc. at landmark
detection_metrics <- function(df, assay_col) {
  df %>%
    filter(!is.na(.data[[assay_col]])) %>%
    summarise(
      Assay                = assay_col,
      N_tested             = n(),
      N_positive           = sum(.data[[assay_col]] == "pos"),
      Prevalence_pos       = N_positive / N_tested,
      Events_tested        = sum(event),
      Sensitivity          = sum(event & (.data[[assay_col]] == "pos")) / Events_tested,
      Specificity          = sum(!event & (.data[[assay_col]] == "neg")) /
        sum(!event),
      Relapse_rate_neg     = sum(event & (.data[[assay_col]] == "neg")) /
        sum(.data[[assay_col]] == "neg"),
      Relapse_rate_pos     = sum(event & (.data[[assay_col]] == "pos")) /
        N_positive
    )
}

det_tbl <- map_dfr(assays, ~detection_metrics(surv_df, .x))
write_csv(det_tbl, file.path(outdir, "Detection_rate_summary.csv"))

## ── 7.  LEAD TIME TO RELAPSE (earliest positive sample vs relapse) ─────────
## For each assay, find earliest positive date AFTER transplant,
## then compute (relapse_date − first_pos_date)
lead_tbl <- map_dfr(assays, function(a) {
  
  ## pull all positive samples for that assay
  pos_tbl <- dat %>%
    filter(.data[[a]] %in% c("pos", "positive", 1, "1", TRUE)) %>%
    select(Patient, first_pos_date = sample_date) %>%
    group_by(Patient) %>%
    summarise(first_pos_date = min(first_pos_date), .groups = "drop")
  
  ## join relapse date
  df <- final_tbl %>%
    select(Patient, relapse_date = censor_date, relapsed) %>%  # censor_date = relapse if relapsed==1
    filter(relapsed == 1) %>%                                  # only real relapses
    left_join(pos_tbl, by = "Patient") %>%
    mutate(
      lead_days = as.numeric(relapse_date - first_pos_date)
    )
  
  summarise(df,
            Assay         = a,
            N_relp_w_pos  = sum(!is.na(lead_days)),
            Median_lead_d = median(lead_days, na.rm = TRUE),
            IQR_lead_d    = IQR(lead_days,  na.rm = TRUE))
})

write_csv(lead_tbl, file.path(outdir, "Lead_time_summary.csv"))

## ── 8.  (OPTIONAL) TIME-DEPENDENT AUC AT 24 & 36 MONTHS ────────────────────
auc_tbl <- map_dfr(assays, function(a) {
  marker <- as.numeric(surv_df[[a]] == "pos")
  roc24  <- timeROC(T = surv_df$time_mo,
                    delta = surv_df$event,
                    marker = marker,
                    cause = 1,
                    times = c(24, 36),
                    iid = TRUE)
  tibble(
    Assay = a,
    AUC24 = roc24$AUC[1],
    AUC36 = roc24$AUC[2]
  )
})

write_csv(auc_tbl, file.path(outdir, "AUC_timeROC_summary.csv"))

## ── 9.  BASELINE CHARACTERISTICS TABLE (optional) ──────────────────────────
vars <- c("AGE", "Gender", "ISS_STAGE", "Cytogenetic_Risk")
baseline_tbl <- dat %>%
  filter(timepoint_info == "diagnosis") %>%              # baseline sample
  select(Patient, all_of(vars)) %>%
  unique()

CreateTableOne(vars = vars, data = baseline_tbl) %>%
  print(varLabels = TRUE, showAllLevels = TRUE)

################################################################################
##  END
################################################################################
