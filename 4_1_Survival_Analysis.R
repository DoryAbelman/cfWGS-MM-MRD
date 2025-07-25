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
dat_rds       <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated3.rds"

## OUTPUT ----------------------------------------------------------------------
outdir <- "Output_tables_2025/detection_progression_updated"
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

## Do rescored 
dat <- dat %>%
  ## Add the screen column 
  mutate(
    BM_base_zscore_screen_call  = as.integer(BM_base_zscore_prob >= 0.350),
  )

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
    BM_base_zscore_call, BM_base_zscore_prob, BM_base_zscore_screen_call,
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
  BM_base_zscore_call    = "cfWGS of BM-derived mutations", 
  BM_base_zscore_screen_call    = "cfWGS of BM-derived mutations (high sensetivity)", 
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
## Can use 1A for this as well 

# 4.  Assays & timepoint definitions ------------------------------------------
assays <- c(
  #  EasyM  = "Rapid_Novor_Binary",
  clonoSEQ = "Adaptive_Binary",
  Flow     = "Flow_Binary",
  cfWGS_BM    = "BM_base_zscore_call",
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



### Make barplot for supplement 
# ─────────────────────────────────────────────────────────────
# 0.  Combine the two tables  → long format
# ─────────────────────────────────────────────────────────────
sens_df <- bind_rows(
  post_stats_BM  %>% mutate(Timepoint = "Post-ASCT"),
  year_stats_BM  %>% mutate(Timepoint = "Maintenance-1yr")
) %>%
  filter(Assay != "cfWGS_Blood") %>%
  mutate(
    # percentages for labelling
    Sens_pct   = Sensitivity * 100,
    # nicer assay labels for the x‑axis
    Assay      = recode(Assay,
                        clonoSEQ       = "clonoSEQ",
                        Flow           = "MFC",
                        cfWGS_BM       = "cfWGS")
  )

# ─────────────────────────────────────────────────────────────
# 1.  Custom palette & base theme
# ─────────────────────────────────────────────────────────────
custom_cols <- c(
  "Post-ASCT"       = "#31688E",  # deep teal
  "Maintenance-1yr" = "#35B779"   # bright green
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title      = element_text(size = 11),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top",
    plot.margin     = margin(10, 10, 30, 10)
  )

# ─────────────────────────────────────────────────────────────
# 2.  Build the grouped bar‑plot
# ─────────────────────────────────────────────────────────────
p_sens <- ggplot(sens_df,
                 aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays among relapsing patients",
    x     = "Technology",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

print(p_sens)

# ─────────────────────────────────────────────────────────────
# 3.  Save
# ─────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Supp_6A_Fig_sensitivity_by_tech_training.png",
  plot     = p_sens,
  width    = 6,
  height   = 4,
  dpi      = 500
)


## Now for blood
sens_df <- bind_rows(
  post_stats_blood  %>% mutate(Timepoint = "Post-ASCT"),
  year_stats_blood  %>% mutate(Timepoint = "Maintenance-1yr")
) %>%
  filter(Assay != "cfWGS_BM") %>%
  mutate(
    # percentages for labelling
    Sens_pct   = Sensitivity * 100,
    # nicer assay labels for the x‑axis
    Assay      = recode(Assay,
                        clonoSEQ       = "clonoSEQ",
                        Flow           = "MFC",
                        cfWGS_Blood       = "cfWGS")
  )

p_sens_blood <- ggplot(sens_df,
                       aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays among relapsing patients",
    x     = "Technology",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

# ─────────────────────────────────────────────────────────────
# 3.  Save
# ─────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Supp_8A_Fig_sensitivity_by_tech_training_blood.png",
  plot     = p_sens,
  width    = 6,
  height   = 4,
  dpi      = 500
)




#### Now get other results and do power analysis 

## 1. Subset to post-transplant & BM-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(BM_base_zscore_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ BM_base_zscore_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ BM_base_zscore_call,
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


## check clonoSEQ
# fit the Kaplan–Meier curve
df_km_clonoSEQ <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "BM_base_zscore_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$BM_base_zscore_call==1)
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
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
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
    !is.na(BM_base_zscore_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ BM_base_zscore_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ BM_base_zscore_call,
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

## Check clonoSEQ 
# fit the Kaplan–Meier curve
df_km_clonoSEQ <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high



## 4. Spearman correlations ----
# replace with your actual probability column:
prob_var <- "BM_base_zscore_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$BM_base_zscore_call==1)
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
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  HR_cf          = hr_cf,
  CI_low_cf      = ci_lo_cf,
  CI_high_cf     = ci_hi_cf,
  HR_fl          = hr_fl,
  CI_low_fl      = ci_lo_fl,
  CI_high_fl     = ci_hi_fl,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
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
  file.path(outdir, "cfWGS_vs_flow_progression_summary_updated.rds")
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

# Now clonoSEQ
df_km_clonoSEQ <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

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
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
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

df_km_clonoSEQ <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Adaptive_Binary)
  )

# fit the Kaplan–Meier curve
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km_clonoSEQ
)

sum_clonoSEQ24   <- summary(fit_clonoSEQ, times = t24)
rfs_neg_clonoSEQ <- sum_clonoSEQ24$surv[1] * 100
rfs_pos_clonoSEQ <- sum_clonoSEQ24$surv[2] * 100

## 3. Median RFS by clonoSEQ ----
fit_clonoSEQ <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
  data = df_km
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

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
  RFS24_seq_neg  = rfs_neg_clonoSEQ,
  RFS24_seq_pos  = rfs_pos_clonoSEQ,
  MedRFS_seq_neg = med_neg_clonoSEQ,
  MedRFS_seq_pos = med_pos_clonoSEQ,
  HR_seq         = hr_clonoSEQ,
  CI_low_seq     = ci_lo_clonoSEQ,
  CI_high_seq    = ci_hi_clonoSEQ,
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
  file.path(outdir, "cfWGS_vs_flow_progression_summaryy_blood_muts_updated.rds")
)



### Make HR figure 
# 1. reshape into long format
hr_plot_df <- progression_metrics_blood %>%
  select(Landmark,
         HR_cf,   CI_low_cf,   CI_high_cf,
         HR_fl,   CI_low_fl,   CI_high_fl,
         HR_seq,   CI_low_seq,   CI_high_seq) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf|fl|seq)"
  ) %>%
  mutate(
    Assay = recode(Assay,
                   cf = "cfWGS",
                   fl = "MFC",
                   seq = "clonoSEQ"),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post‑ASCT", "Maintenance-1yr"))
  )

p_hr <- ggplot(hr_plot_df,
               aes(x = HR, y = fct_rev(Landmark), colour = Assay)) +
  # reference line
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # 1) horizontal CIs
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    position = position_dodge(width = 0.6),
    size     = 0.5
  ) +
  
  # 2) dots at the HR
  geom_point(
    position = position_dodge(width = 0.6),
    size     = 3
  ) +
  
  # log scale axis
  scale_x_continuous(
    "Hazard ratio (log scale)",
    trans        = "log10",
    limits       = c(0.05, 16),           # start at 0.05
    breaks       = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16),
    minor_breaks = c(
      0.06, 0.08,      # between 0.05 & 0.1
      0.15, 0.2,       # between 0.1 & 0.25
      0.3, 0.4,        # between 0.25 & 0.5
      0.6, 0.8,        # between 0.5 & 1
      1.5,             # between 1 & 2
      3, 6,            # between 2 & 4 & 8
      12               # between 8 & 16
    ),
    labels = label_number(accuracy = .01)
  ) +
  annotation_logticks(
    sides  = "b",
    short  = unit(2, "pt"),
    mid    = unit(4, "pt"),
    long   = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    name   = NULL,
    values = c("cfWGS" = "#35608DFF",
               "MFC"   = "#43BF71FF",
               "clonoSEQ"= "#E69F00FF"   # orange for clonoSEQ
    )
  ) +
  
  labs(
    y        = NULL,
    title    = "Relapse hazard ratios stratified by MRD assay\nand landmark timepoint",
    #   subtitle = "cfWGS vs. MFC (95% CI)"
  ) +
  
  # classic theme with no gridlines
  theme_classic(base_size = 11) +
  theme(
    panel.grid         = element_blank(),   # no grid at all
    plot.title         = element_text(face = "bold",
                                      hjust = 0.5),  # bold + centered
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8)
  )

ggsave("Final Tables and Figures/F5X_cfWGS_blood_HR_updated.png",
       p_hr, width = 6, height = 4, dpi = 600)



### Now for BM derived muts
hr_plot_df <- progression_metrics %>%
  select(Landmark,
         HR_cf,   CI_low_cf,   CI_high_cf,
         HR_fl,   CI_low_fl,   CI_high_fl,
         HR_seq,   CI_low_seq,   CI_high_seq) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf|fl|seq)"
  ) %>%
  mutate(
    Assay = recode(Assay,
                   cf = "cfWGS",
                   fl = "MFC",
                   seq = "clonoSEQ"),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post‑ASCT", "Maintenance-1yr"))
  )

p_hr_bm <- ggplot(hr_plot_df,
                  aes(x = HR, y = fct_rev(Landmark), colour = Assay)) +
  # reference line
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # 1) horizontal CIs
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    position = position_dodge(width = 0.6),
    size     = 0.5
  ) +
  
  # 2) dots at the HR
  geom_point(
    position = position_dodge(width = 0.6),
    size     = 3
  ) +
  
  # log scale axis
  scale_x_log10(
    "Hazard ratio (log scale)",
    limits       = c(0.19, 200),           # now starts at 0.2
    breaks       = c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200),
    minor_breaks = c(
      0.25, 0.3, 0.4, 0.6, 0.8,    # between 0.2 & 1
      1.5,           # between 1 & 2
      3, 4,          # between 2 & 5
      6, 8,          # between 5 & 10
      15, 30,        # between 10 & 50
      40, 60, 80,    # between 20 & 100
      150            # between 100 & 200
    ),
    labels = function(x) {
      sapply(x, function(xx) {
        if (xx > 1) {
          sprintf("%.0f", xx)
        } else {
          sprintf("%.2f", xx)
        }
      })
    }
  ) +
  annotation_logticks(
    sides = "b",
    short = unit(2, "pt"),
    mid   = unit(4, "pt"),
    long  = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    name   = NULL,
    values = c("cfWGS" = "#35608DFF",
               "MFC"   = "#43BF71FF",
               "clonoSEQ"= "#E69F00FF")   # orange for clonoSEQ
  ) +
  
  labs(
    y        = NULL,
    title    = "Relapse hazard ratios stratified by MRD assay\nand landmark timepoint",
    #  subtitle = "(95% CI)"
  ) +
  
  # classic theme with no gridlines
  theme_classic(base_size = 11) +
  theme(
    panel.grid         = element_blank(),   # no grid at all
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold",
                                      hjust = 0.5),  # bold + centered
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8)
  )

ggsave("Final Tables and Figures/Figure_4X_cfWGS_BM_HR_updated.png",
       p_hr_bm, width = 6, height = 4, dpi = 600)










### Now make time to relapse figure 
df <- survival_df %>%                           # <- your tibble
  # keep samples beyond baseline / diagnosis
  filter(!str_detect(timepoint_info, regex("Diagnosis|Baseline", TRUE))) %>%
  
  # drop rows with missing probability or time
  filter(!is.na(BM_base_zscore_prob),
         !is.na(Time_to_event)) %>%
  
  # enforce non-negative time to event for relapse event visits a few days off from CMRG date
  mutate(
    days_before_event = pmax(Time_to_event, 0), # set negative values to 0 
    mrd_status      = factor(
      BM_base_zscore_call,
      levels = c(0, 1),
      labels = c("MRD-", "MRD+")
    ),
    progress_status = factor(
      Relapsed_Binary,
      levels = c(0, 1),
      labels = c("No relapse", "Relapse")
    ),
    
    # time *before* the anchor (positive value) → plot reversed
    days_before_event = Time_to_event,    # keep positive for clarity
    months_before_event = days_before_event/30.44
  )

## Only multiple points 
df_slim <- df %>%
  group_by(Patient) %>% 
  filter(dplyr::n() > 1) %>%   # keep only patients with >1 row
  ungroup()

# ────────────────────────────────────────────────────────────────
# 2.  Plot  ──────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
youden_thresh <- 0.436
max_mo <- max(df_slim$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df_slim, aes(months_before_event, BM_base_zscore_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
  
  # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status,
    stroke = mrd_status
  ),
  shape  = 21,
  colour = "black",
  size   = 2
  ) +
  
  # 4) event line
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("cfWGS MRD probability",
                     limits = c(0,1),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # 7) stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0,    # no ring
                   "MRD+" = 1), # visible ring
    guide      = guide_legend(
      title = "MRD call",
      override.aes = list(
        fill   = detect_cols["Both"] %||% c("MRD-"="#35608DFF","MRD+"="#43BF71FF"),
        shape  = 21,
        fill   = "white",    # <–– white interior, not green
        size   = 4,
        colour = "black",
        stroke = c(0,1)
      )
    )
  ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "cfWGS probability trajectories vs. time to relapse"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "right", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 7)    # even smaller
  )

print(p_prob)


# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F4I_cfWGS_prob_vs_time_updated.png",
       p_prob, width = 8, height = 4.5, dpi = 600)



### Instead change the scale to match what Trevor was thinking 
df_plot <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%   # months → days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                       # keep all progressors
      row_number() == which.min(days_before_event)       # keep *latest* censor
  ) %>%
  ungroup()

max_days <- ceiling(max(df_plot$days_before_event, na.rm = TRUE) / 180) * 180

## ─────────────────────────────────────────────────────────────
## 1)  Build the scatter plot                                  
## ─────────────────────────────────────────────────────────────
p_time <- ggplot(df_plot,
                 aes(x = BM_base_zscore_prob,
                     y = days_before_event)) +
  
  # ① dashed horizontal line at “event” (0 days)
  #  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  
  # ② points – colour = outcome, stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape  = 21,
             size   = 3,
             colour = "black") +
  
  # ③ axes
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_y_reverse(
    "Days until relapse (or censor)",
    limits = c(max_days, 0),
    breaks = seq(0, max_days, by = 180),      # every ~6 months
    minor_breaks = seq(0, max_days, by = 90)  # every 3 months
  ) +
  
  # ④ colours for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # ⑤ stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1.1),   # ring only if MRD+
    guide = guide_legend(
      title          = "MRD call",
      override.aes   = list(shape = 21,
                            size  = 4,
                            colour = "black",
                            fill   = "white",
                            stroke = c(0, 1.1))
    )
  ) +
  
  # ⑥ theme / labels
  labs(
    title = "Time to relapse vs. cfWGS MRD probability"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  )

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob.png"),
       p_time, width = 6, height = 4, dpi = 600)


### Show non-relapsers as infinity
# 1) compute days & a “days_plot” that sends non-relapsers to ∞
plot_df2 <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%        # months→days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                                 # keep all relapsers
      row_number() == which.min(days_before_event)                   # for non-relapsers, keep their last sample
  ) %>%
  ungroup()

# define axis maximum and “infinity” sentinel
max_days <- ceiling(max(plot_df2$days_before_event, na.rm = TRUE) / 180) * 180
overflow <- max_days + 30   # a little beyond the longest follow‑up

plot_df2 <- plot_df2 %>%
  mutate(
    days_plot = if_else(progress_status == "No relapse", overflow, days_before_event)
  )


## Check corrs 
# 1) subset to relapsers
rel_df <- plot_df2 %>%
  filter(progress_status == "Relapse")

# 2) run cor.test for each metric
spearman_prob <- cor.test(
  rel_df$BM_base_zscore_prob,
  rel_df$days_before_event,
  method = "spearman"
)

spearman_zscore <- cor.test(
  rel_df$zscore_BM,
  rel_df$days_before_event,
  method = "spearman"
)

spearman_detect_rate <- cor.test(
  rel_df$detect_rate_BM,
  rel_df$days_before_event,
  method = "spearman"
)

# 3) print them
cat("\n--- cfWGS model probability vs days ---\n")
print(spearman_prob)

cat("\n--- zscore_BM vs days ---\n")
print(spearman_zscore)

cat("\n--- detect_rate_BM vs days ---\n")
print(spearman_detect_rate)


# 2) compute Spearman rho on the relapsers
spearman_res <- with(
  filter(plot_df2, progress_status == "Relapse"),
  cor.test(BM_base_zscore_prob,
           days_before_event,
           method = "spearman")
)

## Check other metrics
# 2A) pull out estimate + p‑value
rho  <- spearman_res$estimate
pval <- spearman_res$p.value

# 3) make the scatter
p_time_inf <- ggplot(plot_df2,
                     aes(x = BM_base_zscore_prob,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_hline(yintercept = overflow, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),  # now increasing from 0 to overflow
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(as.character(seq(0, 1620, by = 180)), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0,1),
    breaks = seq(0,1,by=0.1),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF")
  ) +
  
  # stroke legend for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1),
    guide      = guide_legend(
      title = "MRD call",
      override.aes = list(
        shape  = 21,
        size   = 4,
        fill   = "white",
        colour = "black",
        stroke = c(0,1)
      )
    )
  ) +
  
  # annotate Spearman
  annotate("text",
           x = 0.02,      # left margin
           y = 1400,
           label = sprintf("Spearman ρ = %.2f,\np = %.3f", rho, pval),
           hjust = 0,
           size = 3) +
  
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  ) +
  
  labs(title = "Time to relapse vs. cfWGS MRD probability")

# 4) render / save
ggsave(file.path(outdir, "Fig_time_to_relapse_infinity_no_reverse.png"),
       p_time_inf, width = 6, height = 4, dpi = 600)



### Other metrics 
# 3) make the scatter
rho  <- spearman_detect_rate$estimate
pval <- spearman_detect_rate$p.value


p_time_inf <- ggplot(plot_df2,
                     aes(x = detect_rate_BM,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_hline(yintercept = overflow, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),  # now increasing from 0 to overflow
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(as.character(seq(0, 1620, by = 180)), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_log10(
    "cfWGS detection rate (log scale)",
    limits      = c(1e-4, 0.07),                   # can’t include zero
    breaks      = c(1e-4, 1e-3, 1e-2, 1e-1),       # your major decades
    minor_breaks= c(
      2e-4, 5e-4,   # between 1e-4 & 1e-3
      2e-3, 5e-3,   # between 1e-3 & 1e-2
      2e-2, 5e-2    # between 1e-2 & 1e-1
    ),
    labels      = scales::percent_format(accuracy = 0.1)
  ) +
  annotation_logticks(
    sides = "b",     # draw ticks on the bottom axis
    short = unit(2, "pt"),
    mid   = unit(4, "pt"),
    long  = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF")
  ) +
  
  # stroke legend for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1),
    guide      = guide_legend(
      title = "MRD call",
      override.aes = list(
        shape  = 21,
        size   = 4,
        fill   = "white",
        colour = "black",
        stroke = c(0,1)
      )
    )
  ) +
  
  # annotate Spearman
  annotate("text",
           x = 0.01,      # left margin
           y = 1500,
           label = sprintf("Spearman ρ = %.2f,\np = %.3f", rho, pval),
           hjust = 0,
           size = 3) +
  
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  ) +
  
  labs(title = "Time to relapse vs. detection rate")

# 4) render / save
ggsave(file.path(outdir, "Fig_time_to_relapse_infinity_no_reverse_detection_rate.png"),
       p_time_inf, width = 6, height = 4, dpi = 600)




### Now z-score of sites
# 3) make the scatter
rho  <- spearman_zscore$estimate
pval <- spearman_zscore$p.value


p_time_inf <- ggplot(plot_df2,
                     aes(x = zscore_BM,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_hline(yintercept = overflow, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),  # now increasing from 0 to overflow
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(as.character(seq(0, 1620, by = 180)), "∞")
  ) +
  
  scale_x_continuous(
    "cfWGS z-score of proportion of sites detected",
  ) +
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF")
  ) +
  
  # stroke legend for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1),
    guide      = guide_legend(
      title = "MRD call",
      override.aes = list(
        shape  = 21,
        size   = 4,
        fill   = "white",
        colour = "black",
        stroke = c(0,1)
      )
    )
  ) +
  
  # annotate Spearman
  annotate("text",
           x = 230,      # left margin
           y = 1500,
           label = sprintf("Spearman ρ = %.2f,\np = %.3f", rho, pval),
           hjust = 0,
           size = 3) +
  
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  ) +
  
  labs(title = "Time to relapse vs. z-score of proportion of sites detected")

# 4) render / save
ggsave(file.path(outdir, "Fig_time_to_relapse_infinity_no_reverse_sites.png"),
       p_time_inf, width = 6, height = 4, dpi = 600)






### Redo other way for blood derived muts 
### Now redo for blood-derived muts
df <- survival_df %>%                           # <- your tibble
  # keep samples beyond baseline / diagnosis
  filter(!str_detect(timepoint_info, regex("Diagnosis|Baseline", TRUE))) %>%
  
  # drop rows with missing probability or time
  filter(!is.na(Blood_zscore_only_sites_prob),
         !is.na(Time_to_event)) %>%
  
  # enforce non-negative time to event for relapse event visits a few days off from CMRG date
  mutate(
    days_before_event = pmax(Time_to_event, 0), # set negative values to 0 
    mrd_status      = factor(
      Blood_zscore_only_sites_call,
      levels = c(0, 1),
      labels = c("MRD-", "MRD+")
    ),
    progress_status = factor(
      Relapsed_Binary,
      levels = c(0, 1),
      labels = c("No relapse", "Relapse")
    ),
    
    # time *before* the anchor (positive value) → plot reversed
    days_before_event = Time_to_event,    # keep positive for clarity
    months_before_event = days_before_event/30.44
  )

## Only multiple points 
df_slim <- df %>%
  group_by(Patient) %>% 
  filter(dplyr::n() > 1) %>%   # keep only patients with >1 row
  ungroup()

# ────────────────────────────────────────────────────────────────
# 2.  Plot  ──────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
youden_thresh <- 0.523
max_mo <- max(df_slim$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df_slim, aes(months_before_event, Blood_zscore_only_sites_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
  
  # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status,
    stroke = mrd_status
  ),
  shape  = 21,
  colour = "black",
  size   = 2
  ) +
  
  # 4) event line
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("cfWGS MRD probability",
                     limits = c(0.4,0.9),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # 7) stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0,    # no ring
                   "MRD+" = 1), # visible ring
    guide      = guide_legend(
      title = "MRD call",
      override.aes = list(
        fill   = detect_cols["Both"] %||% c("MRD-"="#35608DFF","MRD+"="#43BF71FF"),
        shape  = 21,
        fill   = "white",    # <–– white interior, not green
        size   = 4,
        colour = "black",
        stroke = c(0,1)
      )
    )
  ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "cfWGS probability trajectories vs. time to relapse"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "right", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 7)    # even smaller
  )

print(p_prob)


# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F5I_cfWGS_prob_vs_time_updated_blood_muts.png",
       p_prob, width = 8, height = 4.5, dpi = 600)




## Way Trevor was thinking 
### Instead change the scale to match what Trevor was thinking 
df_plot <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%   # months → days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                       # keep all progressors
      row_number() == which.min(days_before_event)       # keep *latest* censor
  ) %>%
  ungroup()

max_days <- ceiling(max(df_plot$days_before_event, na.rm = TRUE) / 180) * 180


## 1)  Build the scatter plot                                  
p_time <- ggplot(df_plot,
                 aes(x = Blood_zscore_only_sites_call,
                     y = days_before_event)) +
  
  # ① dashed horizontal line at “event” (0 days)
  #  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  
  # ② points – colour = outcome, stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape  = 21,
             size   = 3,
             colour = "black") +
  
  # ③ axes
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_y_reverse(
    "Days until relapse (or censor)",
    limits = c(max_days, 0),
    breaks = seq(0, max_days, by = 180),      # every ~6 months
    minor_breaks = seq(0, max_days, by = 90)  # every 3 months
  ) +
  
  # ④ colours for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # ⑤ stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1.1),   # ring only if MRD+
    guide = guide_legend(
      title          = "MRD call",
      override.aes   = list(shape = 21,
                            size  = 4,
                            colour = "black",
                            fill   = "white",
                            stroke = c(0, 1.1))
    )
  ) +
  
  # ⑥ theme / labels
  labs(
    title = "Time to relapse vs. cfWGS MRD probability for blood derived muts"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  )

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob_blood.png"),
       p_time, width = 6, height = 4, dpi = 600)


### Show non-relapsers as infinity
# 1) compute days & a “days_plot” that sends non-relapsers to ∞
plot_df2 <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%        # months→days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                                 # keep all relapsers
      row_number() == which.min(days_before_event)                   # for non-relapsers, keep their last sample
  ) %>%
  ungroup()

# define axis maximum and “infinity” sentinel
max_days <- ceiling(max(plot_df2$days_before_event, na.rm = TRUE) / 180) * 180
overflow <- max_days + 30   # a little beyond the longest follow‑up

plot_df2 <- plot_df2 %>%
  mutate(
    days_plot = if_else(progress_status == "No relapse", overflow, days_before_event)
  )


## Check corrs 
# 1) subset to relapsers
rel_df <- plot_df2 %>%
  filter(progress_status == "Relapse")

# 2) run cor.test for each metric
spearman_prob <- cor.test(
  rel_df$BM_base_zscore_prob,
  rel_df$days_before_event,
  method = "spearman"
)

spearman_zscore <- cor.test(
  rel_df$zscore_BM,
  rel_df$days_before_event,
  method = "spearman"
)

spearman_detect_rate <- cor.test(
  rel_df$detect_rate_BM,
  rel_df$days_before_event,
  method = "spearman"
)

# 3) print them
cat("\n--- cfWGS model probability vs days ---\n")
print(spearman_prob)

cat("\n--- zscore_BM vs days ---\n")
print(spearman_zscore)

cat("\n--- detect_rate_BM vs days ---\n")
print(spearman_detect_rate)













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
  cfWGS_BM     = "BM_base_zscore_call",
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
    BM_base_zscore_call,
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
           BM_base_zscore_call,
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
bm_col     <- assays["cfWGS_BM"]    # "BM_base_zscore_call"
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


## Make figure 
# ────────────────────────────────────────────────────────────────────────────
# 1) Prepare the data
# ────────────────────────────────────────────────────────────────────────────
sens_BM_df <- results_BM %>%
  # if you want to drop the blood‑only assay, uncomment:
  # filter(Assay != "cfWGS_Blood") %>%
  
  # turn Window_days into a nice factor
  mutate(
    Timepoint = factor(
      Window_days,
      levels = c(90, 180, 365, 730),
      labels = c("90 days", "180 days", "365 days", "730 days")
    ),
    Sens_pct = Sensitivity * 100,
    Assay = recode(
      Assay,
      Flow        = "MFC",
      cfWGS_BM    = "cfWGS",
    )
  )

sens_BM_df <- sens_BM_df %>% filter(Assay != "cfWGS_Blood")
# ────────────────────────────────────────────────────────────────────────────
# 2) Colours & theme (match your existing style)
# ────────────────────────────────────────────────────────────────────────────
custom_cols <- c(
  "90 days"  = "#440154FF",
  "180 days" = "#31688EFF",
  "365 days" = "#35B779FF",
  "730 days" = "#E69F00FF"
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title      = element_text(size = 11),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top",
    plot.margin     = margin(10, 10, 30, 10)
  )

# ────────────────────────────────────────────────────────────────────────────
# 3) Build the grouped bar‑plot
# ────────────────────────────────────────────────────────────────────────────
p_sens_bm <- ggplot(sens_BM_df,
                    aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Window",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays over\nfollow‑up windows (Test Cohort)",
    x     = "Assay",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 14)
  )

# ────────────────────────────────────────────────────────────────────────────
# 4) (Optional) Save
# ────────────────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/Supp_Fig_6_Fig_sensitivity_windows_BM_test_cohort.png",
       plot = p_sens_bm,
       width = 4.75, height = 6, dpi = 500)


### Now remake for blood muts
sens_blood_df <- results_blood %>%
  # if you want to drop the blood‑only assay, uncomment:
  # filter(Assay != "cfWGS_Blood") %>%
  
  # turn Window_days into a nice factor
  mutate(
    Timepoint = factor(
      Window_days,
      levels = c(90, 180, 365, 730),
      labels = c("90 days", "180 days", "365 days", "730 days")
    ),
    Sens_pct = Sensitivity * 100,
    Assay = recode(
      Assay,
      Flow        = "MFC",
      cfWGS_Blood    = "cfWGS",
    )
  )

sens_blood_df <- sens_blood_df %>% filter(Assay != "cfWGS_BM")

p_sens_blood <- ggplot(sens_blood_df,
                       aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Window",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays over\nfollow‑up windows (Test Cohort)",
    x     = "Assay",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ────────────────────────────────────────────────────────────────────────────
# 4) (Optional) Save
# ────────────────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/Supp_Fig_8_Fig_sensitivity_windows_blood_test_cohort.png",
       plot = p_sens_blood,
       width = 5, height = 5, dpi = 500)





### Export this
# full results (all patients with any assay)
#write_csv(
#  results,
#  file.path(outdir, "all_assays_timewindow_results.csv")
#)

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
assay_prob <- "BM_base_zscore_prob"  # or "BM_base_zscore_prob"

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
    BM_base_zscore_call, BM_base_zscore_prob,
    Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob, Cumulative_VAF_BM, Cumulative_VAF_blood,
    zscore_BM, zscore_blood, z_score_detection_rate_BM, z_score_detection_rate_blood, detect_rate_BM, detect_rate_blood
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
    BM_base_zscore_call, BM_base_zscore_prob,
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

