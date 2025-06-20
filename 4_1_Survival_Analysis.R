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
    PET_Binary,
    BM_zscore_only_detection_rate_call,
    Blood_zscore_only_sites_call
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
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    # skip if too few pts
    if(nrow(df_sub) < min_n) next
    
    # skip if only one group present
    if(n_distinct(df_sub$Group) < 2) next
    
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = TRUE,
      conf.int        = TRUE,
      risk.table      = TRUE,
      palette         = c("#E7B800","#2E9FDF"),
      legend.title    = paste0(assay_lab, " MRD"),
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Days since sample",
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
