### Script for assembling what is needed for survival analysis 
## author: Dory 
## Last updated: June 18, 2025



# ── libraries ──────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)   # pretty KM curves + log-rank p
library(broom)       # clean model outputs
library(timeROC)     # time dependent ROC

# ── file paths (EDIT for your project layout) ──────────────────────────────
## raw objects you already have on disk or in memory
# 1) sample-level clinical + MRD + dates info
#    (this is the 'combined_clinical_data_updated' you created earlier)
clinical_rds      <- "combined_clinical_data_updated_April2025.csv"

# 2) a table with *one* row per patient giving the first documented relapse/
#    progression date (can be NA) – you called it 'Relapse_dates_full'
relapse_csv       <- "Relapse_dates_full.csv"

# 3) table of “last-seen” dates for patients that have **not** yet relapsed
#    (you generated this as 'latest_dates_per_patient.csv')
lastfu_csv        <- "latest_dates_per_patient.csv"

# 4) long-format MRD results – one row per sample per assay
#    MUST contain: Patient, Sample_ID, timepoint_info, assay, status
#    If you stored MRD results in several objects, rbind them first.
mrd_long_rds      <- "data/MRD_results_long.rds"

# ── output folder ──────────────────────────────────────────────────────────
outdir <- "output_tables_2025"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


## 1A  clinical sample table  ───────────────────────────────────────────────
clin <- readRDS(clinical_rds) %>%
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection),   # ensure Date
    Patient                   = as.character(Patient)
  )

## 1B  relapse dates  ───────────────────────────────────────────────────────
relapse_dates <- read_csv(relapse_csv, show_col_types = FALSE) %>%
  transmute(
    Patient         = as.character(Patient),
    relapse_date    = as.Date(Progression_date)         # rename for clarity
  )

## 1C  last follow-up (censoring) dates  ────────────────────────────────────
last_fu <- read_csv(lastfu_csv, show_col_types = FALSE) %>%
  transmute(
    Patient      = as.character(Patient),
    last_fu_date = as.Date(latest_date)
  )


## Now build PFS table 

## 2A  baseline date per patient -------------------------------------------
baseline_tbl <- clin %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  group_by(Patient) %>%
  summarise(baseline_date = min(Date_of_sample_collection, na.rm = TRUE),
            .groups = "drop")

## 2B  merge baseline, relapse and last-FU dates ----------------------------
pfs_tbl <- baseline_tbl %>%
  left_join(relapse_dates, by = "Patient") %>%
  left_join(last_fu,      by = "Patient") %>%
  mutate(
    # choose relapse date if it exists, otherwise last FU
    final_date = coalesce(relapse_date, last_fu_date),
    event      = if_else(!is.na(relapse_date), 1L, 0L),
    time_days  = as.numeric(final_date - baseline_date)
  ) %>%
  select(Patient, time_days, event)

## 2C  sanity check ---------------------------------------------------------
stopifnot(all(pfs_tbl$time_days >= 0))
pfs_tbl %>% summarise(median_follow_up_mo = median(time_days[event == 0] / 30.44))

## 2D  export ---------------------------------------------------------------
write_csv(pfs_tbl, file.path(outdir, "PFS_days.csv"))


### Now get MRD landmark from combined clinical data updated






#### Now plot 

## Build dataframe
## --- join PFS table with MRD landmark data -----------------------------
surv_df <- pfs_tbl %>%                       # <-- your PFS_days renamed
  select(Patient, time_days, event) %>%
  left_join(mrd_landmark, by = "Patient") %>%
  mutate(time_months = time_days / 30.44)    # convert to months for plotting


### Function to get KM curves
plot_km <- function(df, assay_col, landmark = "Post-ASCT") {
  form <- as.formula(paste0("Surv(time_months, event) ~ ", assay_col))
  fit  <- survfit(form, data = df)
  
  ggsurv <- ggsurvplot(
    fit,
    data          = df,
    risk.table    = TRUE,
    pval          = TRUE,          # log-rank p
    conf.int      = FALSE,
    palette       = c("#1b9e77", "#d95f02"),  # green / orange
    legend.title  = assay_col,
    legend.labs   = c("MRD- negative", "MRD- positive"),
    xlab          = "Months from transplant",
    title         = paste(assay_col, "at", landmark, "— RFS")
  )
  
  ## HR from univariable Cox
  cox <- coxph(form, data = df)
  tidy_cox <- broom::tidy(cox, exp = TRUE) %>%        # exp() gives HR
    select(HR = estimate, CI_low = conf.low, CI_high = conf.high, p.value)
  
  list(plot = ggsurv, hr = tidy_cox)
}


## Run every assay
assays <- c("EasyM", "cfWGS", "clonoSEQ", "MFC")   # update list as needed

km_results <- map(assays, ~plot_km(surv_df, .x))

## --- Export plots -------------------------------------------------------
dir.create("figures/KM", showWarnings = FALSE)
walk2(km_results, assays, ~ggsave(
  filename = file.path("figures/KM", paste0(.y, "_km.pdf")),
  plot = .x$plot$plot,
  width = 5, height = 4))

## --- Compile HR table for manuscript -----------------------------------
hr_tbl <- map2_dfr(km_results, assays, ~mutate(.x$hr, Assay = .y)) %>%
  select(Assay, HR, CI_low, CI_high, p.value)

write_csv(hr_tbl, "figures/KM/HR_table.csv")




### Get other stats 
## Concordance (Harrell's C) for each assay
c_stats <- map2_dfr(assays, assays, function(a, name) {
  coxph(as.formula(paste0("Surv(time_months, event) ~ ", a)), data = surv_df) |>
    summary() |>
    {\(x) tibble(Assay = name, C_index = x$concordance[1])}()
})

## Time-dependent ROC at 24 months
roc24 <- map_dfr(assays, function(a) {
  roc <- timeROC(
    T       = surv_df$time_months,
    delta   = surv_df$event,
    marker  = as.numeric(surv_df[[a]] == "pos"),
    cause   = 1,
    times   = 24,
    iid     = TRUE
  )
  tibble(Assay = a, AUC_24mo = roc$AUC)
})
