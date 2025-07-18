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
final_tbl_rds <- "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"
dat_rds       <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated2.rds"

## OUTPUT ----------------------------------------------------------------------
outdir <- "Output_tables_2025/detection_progression"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ── 1.  LOAD & TIDY CORE TABLES ──────────────────────────────────────────────
final_tbl <- readRDS(final_tbl_rds) %>%
  rename_with(
    tolower,
    any_of(c("Baseline_Date", "Censor_date", "Relapsed"))
  ) %>%
  transmute(
    Patient       = as.character(Patient),
    baseline_date = as.Date(baseline_date),
    censor_date   = as.Date(censor_date),
    relapsed      = as.integer(relapsed)
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
  BM_zscore_only_detection_rate_call    = "cfWGS of BM-derived mutations", 
  Blood_zscore_only_sites_call = "cfWGS of PB‑cfDNA‑derived mutations"
)

# 3) All timepoints to cover
tps <- unique(survival_df$timepoint_info)
dpi_target <- 500

# 4) Minimum n per group to plot
min_n <- 5

pal_2 <- viridis(2, option = "D", begin = 0.3, end = 0.7)   # nice mid‑range hues
# old colors: palette = c("#E7B800","#2E9FDF")

tp_labels <- c(
  `diagnosis` = "diagnosis",
  `post_transplant`    = "post‑ASCT",
  `1yr maintenance` = "one year maintenance", 
  `post_induction`    = "post‑induction"
  )

# 5) Loop
for(tp in tps) {
  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated.png"))
    
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
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),  # ← left‑align
      palette         = pal_2,
     # legend.title    = paste0(assay_lab, " MRD"),
     legend.title    = "MRD status",
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Time since MRD assessment (months)",
      ylab            = "Progression-free survival",
      title = str_wrap(paste0("PFS stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 16),
          legend.position = "top",
          axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),          # no grid
          panel.grid.minor = element_blank(),
          #  Make the tick‑labels (the numbers) larger:
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y      = element_text(size = 15),
          axis.title.x      = element_text(size = 14)
        ),
    )
    
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title      = element_text(hjust = 0, face = "plain"),
      )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3.5,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 8,
      dpi      = dpi_target
    )
  }
}





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
# fit the Kaplan–Meier curve
df_km_fl <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Flow_Binary)
  )

# fit the Kaplan–Meier curve
fit_fl <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
  data = df_km_fl
)

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

## Get additional dates
Relapse_dates_full <- read_csv(
  "Relapse dates cfWGS updated.csv",
  col_types = cols(
    Patient          = col_character(),
    Progression_date = col_date(format = "%Y-%m-%d")
  )
)

# 2) Build df_sf: non-frontline samples + assays + relapse info
df_sf <- dat %>%
  filter(tolower(Cohort) == "non-frontline") %>%
  transmute(
    Patient,
    sample_date = as.Date(Date),
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

## Edit if progression date earlier than sample collected, since some patients had multiple
df_sf2 <- df_sf %>%
  select(-relapse_date, -relapsed) %>%           # ① drop the old pair
  left_join(Relapse_dates_full, by = "Patient") %>%
  filter(Progression_date >= sample_date) %>%
  group_by(Patient, sample_date,
           Flow_Binary,
           BM_zscore_only_detection_rate_call,
           Blood_zscore_only_sites_call) %>%
  slice_min(Progression_date, with_ties = FALSE) %>%
  ungroup() %>%
  rename(relapse_date = Progression_date) %>%    # ② now safe to rename
  mutate(relapsed = as.integer(!is.na(relapse_date)))


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
results <- imap_dfr(assays, 
                    # .x = column name, .y = assay label
                    .f = function(col_name, assay_label) {
                      map_dfr(windows, function(win_d) {
                        calc_metrics(df_sf2, assay_label, col_name, win_d)
                      })
                    }
)

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


## Get info
library(scales)   # for percent()

results_BM %>%
  # pick the assays & windows you want to narrate
  filter(Assay %in% c("cfWGS_BM", "Flow"), Window_days %in% c(180, 365)) %>%
  # build a sentence for each row
  rowwise() %>%
  mutate(
    sentence = glue(
      "{Assay} at {Window_days}-day window detected {TP}/{TP + FN} progressors ",
      "(sensitivity {percent(Sensitivity)}, specificity {percent(Specificity)})."
    )
  ) %>%
  ungroup() %>%
  # print them to the console
  pull(sentence) %>%
  cat(sep = "\n")

## For blood
results_blood %>%
  # pick the assays & windows you want to narrate
  filter(Assay %in% c("cfWGS_Blood", "Flow"), Window_days %in% c(180, 365)) %>%
  # build a sentence for each row
  rowwise() %>%
  mutate(
    sentence = glue(
      "{Assay} at {Window_days}-day window detected {TP}/{TP + FN} progressors ",
      "(sensitivity {percent(Sensitivity)}, specificity {percent(Specificity)})."
    )
  ) %>%
  ungroup() %>%
  # print them to the console
  pull(sentence) %>%
  cat(sep = "\n")

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
print(event_counts_BM)


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




### See at what pont an increase occured 
# 1) choose the probability column you want to analyze:
assay_prob <- "BM_zscore_only_detection_rate_prob"  # or "BM_zscore_only_detection_rate_prob"

## Reload dat 
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )


### Now see how early it was to progression
survival_df <- dat %>%
  filter(Cohort == "Frontline") %>%
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
    Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob, Cumulative_VAF_BM, Cumulative_VAF_blood
  )


# 1) build advance_df: per‐patient nadir & first increase before progression
advance_df <- survival_df %>%
  filter(
    Relapsed_Binary == 1,
    !is.na(.data[[assay_prob]]),
    # include samples up to and including the progression date
    sample_date <= censor_date
  ) %>%
  group_by(Patient) %>%
  arrange(sample_date) %>%
  group_modify(~ {
    df <- .
    # find the nadir (min prob)
    nadir_idx   <- which.min(df[[assay_prob]])
    nadir_date  <- df$sample_date[nadir_idx]
    nadir_prob  <- df[[assay_prob]][nadir_idx]
    # progression date
    prog_date   <- df$censor_date[1]
    # find first sample after nadir with prob > nadir
    inc_idx     <- which(df$sample_date > nadir_date &
                           df[[assay_prob]] > nadir_prob)
    if (length(inc_idx)) {
      first_inc_date <- df$sample_date[inc_idx[1]]
      first_inc_prob <- df[[assay_prob]][inc_idx[1]]
    } else {
      first_inc_date <- as.Date(NA)
      first_inc_prob <- NA_real_
    }
    tibble(
      nadir_date             = nadir_date,
      nadir_prob             = nadir_prob,
      first_inc_date         = first_inc_date,
      first_inc_prob         = first_inc_prob,
      days_to_first_increase = as.numeric(first_inc_date - nadir_date),
      days_before_progression= as.numeric(prog_date - first_inc_date),
      prog_date              = prog_date
    )
  }) %>%
  ungroup()

# A) progressing‐patient dataset with non‐missing assay values
df_prog <- survival_df %>%
  filter(
    Relapsed_Binary == 1,
    !is.na(.data[[assay_prob]]),
    sample_date <= censor_date    # ← restrict to pre-progression here
  ) %>%
  arrange(Patient, sample_date) %>%
  # compute days before progression for each sample
  mutate(days_before_prog = as.numeric(censor_date - sample_date))

# B) total samples & timing relative to progression
n_patients <- df_prog %>% pull(Patient) %>% unique() %>% length()
n_samples  <- nrow(df_prog)
sample_timing_stats <- df_prog %>%
  summarise(
    median_days_before = median(days_before_prog),
    iqr_days_before    = IQR(days_before_prog),
    min_days_before    = min(days_before_prog),
    max_days_before    = max(days_before_prog),
    .groups = "drop"
  )

# C) nadir timing stats (days before progression)
nadir_timing_stats <- advance_df %>%
  mutate(
    days_nadir_before = as.numeric(prog_date - nadir_date)
  ) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before, na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before, na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before, na.rm = TRUE),
    .groups = "drop"
  )

# D) first‐increase timing stats (already in advance_df)
increase_stats <- advance_df %>%
  summarise(
    median_to_increase = median(days_to_first_increase, na.rm = TRUE),
    iqr_to_increase    = IQR(days_to_first_increase, na.rm = TRUE),
    min_inc_days       = min(days_to_first_increase, na.rm = TRUE),
    max_inc_days       = max(days_to_first_increase, na.rm = TRUE),
    .groups = "drop"
  )

# E) print results‐section sentences
cat(glue(
  "We analyzed {n_patients} patients (total {n_samples} samples) collected a median ",
  "{sample_timing_stats$median_days_before} days before progression (IQR ",
  "{sample_timing_stats$iqr_days_before} days; range ",
  "{sample_timing_stats$min_days_before}–{sample_timing_stats$max_days_before} days). ",
  "The nadir detection probability occurred a median ",
  "{nadir_timing_stats$median_nadir_days} days before progression ",
  "(IQR {nadir_timing_stats$iqr_nadir_days} days; range ",
  "{nadir_timing_stats$min_nadir_days}–{nadir_timing_stats$max_nadir_days} days). ",
  "From that nadir, the first increase occurred a median ",
  "{increase_stats$median_to_increase} days later ",
  "(IQR {increase_stats$iqr_to_increase} days; range ",
  "{increase_stats$min_inc_days}–{increase_stats$max_inc_days} days) before progression.\n"
))

# summarise both “after nadir” and “before progression”
increase_stats2 <- advance_df %>%
  summarise(
    # after nadir → first increase
    median_after_nadir = median(days_to_first_increase, na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,   na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,   na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,   na.rm = TRUE),
    # first increase → progression
    median_before_prog = median(days_before_progression, na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,   na.rm = TRUE),
    min_before_prog    = min(days_before_progression,   na.rm = TRUE),
    max_before_prog    = max(days_before_progression,   na.rm = TRUE),
    .groups            = "drop"
  )

# then one sentence with both
cat(glue(
  "From the nadir, the first increase occurred a median ",
  "{increase_stats2$median_after_nadir} days later ",
  "(IQR {increase_stats2$iqr_after_nadir} days; range ",
  "{increase_stats2$min_after_nadir}–{increase_stats2$max_after_nadir} days), ",
  "which was a median {increase_stats2$median_before_prog} days before progression ",
  "(IQR {increase_stats2$iqr_before_prog} days; range ",
  "{increase_stats2$min_before_prog}–{increase_stats2$max_before_prog} days).\n"
))





### Now for blood muts
# 1) choose the probability column you want to analyze:
assay_prob <- "Blood_zscore_only_sites_prob"  

## Reload dat 
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )


### Now see how early it was to progression
survival_df <- dat %>%
  filter(Cohort == "Frontline") %>%
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
    Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob, Cumulative_VAF_BM, Cumulative_VAF_blood
  )


# 1) build advance_df: per‐patient nadir & first increase before progression
advance_df <- survival_df %>%
  filter(
    Relapsed_Binary == 1,
    !is.na(.data[[assay_prob]]),
    # include samples up to and including the progression date
    sample_date <= censor_date
  ) %>%
  group_by(Patient) %>%
  arrange(sample_date) %>%
  group_modify(~ {
    df <- .
    # find the nadir (min prob)
    nadir_idx   <- which.min(df[[assay_prob]])
    nadir_date  <- df$sample_date[nadir_idx]
    nadir_prob  <- df[[assay_prob]][nadir_idx]
    # progression date
    prog_date   <- df$censor_date[1]
    # find first sample after nadir with prob > nadir
    inc_idx     <- which(df$sample_date > nadir_date &
                           df[[assay_prob]] > nadir_prob)
    if (length(inc_idx)) {
      first_inc_date <- df$sample_date[inc_idx[1]]
      first_inc_prob <- df[[assay_prob]][inc_idx[1]]
    } else {
      first_inc_date <- as.Date(NA)
      first_inc_prob <- NA_real_
    }
    tibble(
      nadir_date             = nadir_date,
      nadir_prob             = nadir_prob,
      first_inc_date         = first_inc_date,
      first_inc_prob         = first_inc_prob,
      days_to_first_increase = as.numeric(first_inc_date - nadir_date),
      days_before_progression= as.numeric(prog_date - first_inc_date),
      prog_date              = prog_date
    )
  }) %>%
  ungroup()

# A) progressing‐patient dataset with non‐missing assay values
df_prog <- survival_df %>%
  filter(
    Relapsed_Binary == 1,
    !is.na(.data[[assay_prob]]),
    sample_date <= censor_date    # ← restrict to pre-progression here
  ) %>%
  arrange(Patient, sample_date) %>%
  # compute days before progression for each sample
  mutate(days_before_prog = as.numeric(censor_date - sample_date))

# B) total samples & timing relative to progression
n_patients <- df_prog %>% pull(Patient) %>% unique() %>% length()
n_samples  <- nrow(df_prog)
sample_timing_stats <- df_prog %>%
  summarise(
    median_days_before = median(days_before_prog),
    iqr_days_before    = IQR(days_before_prog),
    min_days_before    = min(days_before_prog),
    max_days_before    = max(days_before_prog),
    .groups = "drop"
  )

# C) nadir timing stats (days before progression)
nadir_timing_stats <- advance_df %>%
  mutate(
    days_nadir_before = as.numeric(prog_date - nadir_date)
  ) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before, na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before, na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before, na.rm = TRUE),
    .groups = "drop"
  )

# D) first‐increase timing stats (already in advance_df)
increase_stats <- advance_df %>%
  summarise(
    median_to_increase = median(days_to_first_increase, na.rm = TRUE),
    iqr_to_increase    = IQR(days_to_first_increase, na.rm = TRUE),
    min_inc_days       = min(days_to_first_increase, na.rm = TRUE),
    max_inc_days       = max(days_to_first_increase, na.rm = TRUE),
    .groups = "drop"
  )

# E) print results‐section sentences
cat(glue(
  "We analyzed {n_patients} patients (total {n_samples} samples) collected a median ",
  "{sample_timing_stats$median_days_before} days before progression (IQR ",
  "{sample_timing_stats$iqr_days_before} days; range ",
  "{sample_timing_stats$min_days_before}–{sample_timing_stats$max_days_before} days). ",
  "The nadir detection probability occurred a median ",
  "{nadir_timing_stats$median_nadir_days} days before progression ",
  "(IQR {nadir_timing_stats$iqr_nadir_days} days; range ",
  "{nadir_timing_stats$min_nadir_days}–{nadir_timing_stats$max_nadir_days} days). ",
  "From that nadir, the first increase occurred a median ",
  "{increase_stats$median_to_increase} days later ",
  "(IQR {increase_stats$iqr_to_increase} days; range ",
  "{increase_stats$min_inc_days}–{increase_stats$max_inc_days} days) before progression.\n"
))


# summarise both “after nadir” and “before progression”
increase_stats2 <- advance_df %>%
  summarise(
    # after nadir → first increase
    median_after_nadir = median(days_to_first_increase, na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,   na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,   na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,   na.rm = TRUE),
    # first increase → progression
    median_before_prog = median(days_before_progression, na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,   na.rm = TRUE),
    min_before_prog    = min(days_before_progression,   na.rm = TRUE),
    max_before_prog    = max(days_before_progression,   na.rm = TRUE),
    .groups            = "drop"
  )

# then one sentence with both
cat(glue(
  "From the nadir, the first increase occurred a median ",
  "{increase_stats2$median_after_nadir} days later ",
  "(IQR {increase_stats2$iqr_after_nadir} days; range ",
  "{increase_stats2$min_after_nadir}–{increase_stats2$max_after_nadir} days), ",
  "which was a median {increase_stats2$median_before_prog} days before progression ",
  "(IQR {increase_stats2$iqr_before_prog} days; range ",
  "{increase_stats2$min_before_prog}–{increase_stats2$max_before_prog} days).\n"
))






### Maybe follow up with Esteban to see if clinical or biochemical progression since don't know for non-frontline
### Leave as blank for now
# Filter out patients whose ID starts with "IMG-" or "SPORE-"
relapse_filtered <- Relapse_dates_full %>%
  filter(!grepl("^IMG|^SPORE", Patient))

# Export to RDS
saveRDS(relapse_filtered, file = "Relapse_dates_full_filtered.rds")

# Export to CSV (optional)
write.csv(relapse_filtered, file = "Relapse_dates_full_filtered.csv", row.names = FALSE)

