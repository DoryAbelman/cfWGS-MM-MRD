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
dat_rds       <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated5.rds"

## OUTPUT ----------------------------------------------------------------------
outdir <- "Output_tables_2025/detection_progression_updated6"
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
    BM_zscore_only_detection_rate_screen_call  = as.integer(BM_zscore_only_detection_rate_prob >= 0.350),
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
    BM_zscore_only_detection_rate_call, BM_zscore_only_detection_rate_prob, BM_zscore_only_detection_rate_screen_call,
    Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob, Blood_base_prob, Blood_base_call,
    Blood_plus_fragment_prob, Blood_plus_fragment_call,
    Blood_plus_fragment_min_prob, Blood_plus_fragment_min_call, Fragmentomics_mean_coverage_only_prob, Fragmentomics_mean_coverage_only_call
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
  BM_zscore_only_detection_rate_call    = "cfWGS of BM-Derived Mutations (cVAF Model)", 
  BM_zscore_only_detection_rate_screen_call    = "cfWGS of BM-derived mutations (high sensetivity)", 
  Blood_zscore_only_sites_call = "cfWGS of cfDNA-Derived Mutations (Sites Model)",
  Blood_plus_fragment_min_call = "cfWGS of cfDNA-Derived Mutations (Combined Model)",
  Fragmentomics_mean_coverage_only_call = "Fragmentomics model"
)

# techs <- c(
#   Blood_base_call = "cfWGS of PB‑cfDNA‑derived mutations (combined model)"
# )

# 3) All timepoints to cover
tps <- unique(survival_df$timepoint_info)
dpi_target <- 500

# 4) Minimum n per group to plot
min_n <- 5

#pal_2 <- viridis(2, option = "D", begin = 0.3, end = 0.7)   # nice mid‑range hues
# old colors: palette = c("#E7B800","#2E9FDF")
pal_2 <- c("black", "red") # updated 

tp_labels <- c(
  `diagnosis` = "Diagnosis",
  `post_transplant`    = "Post‑ASCT",
  `1yr maintenance` = "One-Year Maintenance", 
  `post_induction`    = "Post‑Induction"
)

# once before the loops
dx_tbl <- survival_df %>%
  group_by(Patient) %>%
  summarise(
    diagnosis_date = suppressWarnings(min(sample_date[timepoint_info == "diagnosis"], na.rm = TRUE)),
    .groups = "drop"
  )


# 5) Loop
for(tp in tps) {
#  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_no_CI.png"))
    
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
      conf.int        = FALSE,
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
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
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
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}


## Redo with CI 
# 5) Loop
for(tp in tps) {
  #  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_with_CI.png"))
    
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
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
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
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}


### Set the CI instead to be to 90% 
for(tp in tps) {
  #  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_with_CI_90.png"))
    
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
    
    # 90% CI from survfit; "log" (Greenwood on log scale) is common, "log-log" is also fine
    fit <- survfit(
      surv_obj ~ Group,
      data      = df_sub,
      conf.int  = 0.90     # ← 90% CI
    )
    
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = TRUE,
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = TRUE,
      conf.int.alpha  = 0.1,    
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
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
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
    
    # # Convert survfit object to a data.frame
    # fit_df <- broom::tidy(fit)  # gives time, n.risk, n.event, surv, std.err, conf.low, conf.high, strata
    # 
    # # Overlay ribbons for CIs
    # km$plot <- km$plot +
    #   geom_ribbon(
    #     data = fit_df,
    #     aes(x = time, ymin = conf.low, ymax = conf.high, fill = strata, group = strata),
    #     inherit.aes = FALSE, alpha = 0.2
    #   )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}

 



#### Redo but show the plot starting from diagnosis
### To complex since the default risk table doesn't support delayed entry - skipped
### delayed entry (a.k.a. left truncation) to avoid immortal-bias
# pretty p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-4) return("p < 1e-4")
  if (p < 1e-3) return("p < 0.001")
  if (p < 1e-2) return("p < 0.01")
  sprintf("p = %.2f", p)
}

# risk table for delayed-entry (start–stop) data
make_risktable <- function(df, breaks) {
  tmp <- df %>%
    mutate(Group_label = dplyr::recode(as.character(Group),
                                       "Negative" = "MRD–",
                                       "Positive" = "MRD+"))
  purrr::map_dfr(breaks, function(ti) {
    tmp %>%
      dplyr::group_by(Group_label) %>%
      dplyr::summarise(n = sum(entry_m <= ti & exit_m > ti), .groups = "drop") %>%
      dplyr::mutate(time = ti)
  }) %>%
    tidyr::pivot_wider(names_from = time, values_from = n) %>%
    dplyr::arrange(factor(Group_label, levels = c("MRD–","MRD+"))) %>%
    dplyr::rename(`MRD status` = Group_label)
}


for (tp in tps) {
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_from_diagnosis.png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%
      ungroup() %>%
      left_join(dx_tbl, by = "Patient") %>%
      mutate(
        Group   = factor(ifelse(.data[[var]] == 1, "Positive", "Negative"),
                         levels = c("Negative","Positive")),
        # entry = months from diagnosis to MRD test
        entry_m = as.numeric(sample_date - diagnosis_date) / 30.44,
        # exit = entry + observed time after MRD test
        exit_m  = entry_m + (Time_to_event / 30.44)
      ) %>%
      filter(!is.na(entry_m), !is.na(exit_m), exit_m >= entry_m)
    
    if (nrow(df_sub) < min_n) next
    if (n_distinct(df_sub$Group) < 2) next
    
    # Cox with delayed entry
    surv_obj <- Surv(time = df_sub$entry_m, time2 = df_sub$exit_m, event = df_sub$Relapsed_Binary)
    fit     <- survfit(surv_obj ~ Group, data = df_sub)
    
    cox_fit <- coxph(surv_obj ~ Group, data = df_sub, ties = "breslow")
    s <- summary(cox_fit)
    
    # Pull p-values (log-rank / score / wald) safely
    p_lrt   <- suppressWarnings(as.numeric(s$logtest[3]))  # likelihood-ratio p
    p_score <- suppressWarnings(as.numeric(s$sctest[3]))   # score test p
    p_wald  <- suppressWarnings(as.numeric(s$waldtest[3])) # Wald p
    
    # Fallback if all above are NA (e.g., separation): grab coefficient p if present
    coef_p <- NA_real_
    cs <- try(coef(summary(cox_fit)), silent = TRUE)
    if (!inherits(cs, "try-error")) {
      # find the row for Group (works even if level names change)
      r <- grep("^Group", rownames(cs), value = FALSE)
      if (length(r) >= 1) coef_p <- cs[r[1], "Pr(>|z|)"]
    }
    
    # Choose the best available p
    p_any <- if (!is.na(p_lrt)) p_lrt else if (!is.na(p_score)) p_score else if (!is.na(p_wald)) p_wald else coef_p
    
    # Format: show thresholds instead of tiny decimals
    fmt_p <- function(p) {
      if (is.na(p)) return("p = NA")
      if (p < 1e-4) return("p < 1e-4")
      if (p < 1e-3) return("p < 0.001")
      if (p < 1e-2) return("p < 0.01")
      sprintf("p = %.2f", p)
    }
    pval_str <- fmt_p(p_any)
    
    km <- ggsurvplot(
      fit, data = df_sub,
      pval              = pval_str,
      break.time.by     = 12,
      conf.int          = TRUE,
      risk.table        = TRUE,
      risk.table.title  = "Number at risk",                    # ← add
      risk.table.title.theme = element_text(hjust = 0),        # ← add (left-align)
      palette           = pal_2,
      legend.title      = "MRD status",
      legend.labs       = c("MRD-","MRD+"),  # use ASCII hyphen to avoid file/device issues
      xlab              = "Months since diagnosis",
      ylab              = "Progression-free survival",
      title             = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25,
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position  = "top",
          axis.line        = element_line(colour = "black"),   # ← add
          panel.grid.major = element_blank(),                  # ← add
          panel.grid.minor = element_blank(),                  # ← add
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y     = element_text(size = 15),
          axis.title.x     = element_text(size = 14)
        )
    )
    
    # match the old post-processing of table/plot
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title   = element_text(hjust = 0, face = "plain")
      )
    
    km$plot <- km$plot +
      theme(axis.title.x = element_blank())
    
    combined <- ggarrange(km$plot, km$table, ncol = 1, heights = c(3,1))
    ggsave(fname, plot = combined, width = 7, height = 7, dpi = dpi_target)
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
  cfWGS_BM    = "BM_zscore_only_detection_rate_call",
  cfWGS_BM_screen    = "BM_zscore_only_detection_rate_screen_call",
  cfWGS_Blood_Sites    = "Blood_zscore_only_sites_call",
  cfWGS_Blood_Combined = "Blood_plus_fragment_call"
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

blood_col <- assays["cfWGS_Blood_Sites"]

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
  filter(Assay != "cfWGS_Blood_Sites") %>%
  filter(Assay != "cfWGS_Blood_Combined") %>%
  mutate(
    # percentages for labelling
    Sens_pct   = Sensitivity * 100,
    # nicer assay labels for the x‑axis
    Assay      = recode(Assay,
                        clonoSEQ       = "clonoSEQ",
                        Flow           = "MFC",
                        cfWGS_BM       = "cfWGS"),
    # enforce desired order
    Timepoint = factor(Timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
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
p_sens <- ggplot(sens_df %>% filter(Assay != "cfWGS_BM_screen"),
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
    values = custom_cols[c("Post-ASCT", "Maintenance-1yr")],  # enforce mapping
    limits = c("Post-ASCT", "Maintenance-1yr")                # enforce order
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of cfDNA and Clinical MRD Assays in Relapsing Patients",
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
  filename = "Final Tables and Figures/Supp_6A_Fig_sensitivity_by_tech_training3.png",
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
  filter(Assay != "cfWGS_BM_screen") %>%
  mutate(
    # percentages for labelling
    Sens_pct   = Sensitivity * 100,
    # nicer assay labels for the x-axis
    Assay      = recode(Assay,
                        clonoSEQ              = "clonoSEQ",
                        Flow                  = "MFC",
                        cfWGS_Blood_Sites     = "cfWGS\n(Sites Model)",
                        cfWGS_Blood_Combined  = "cfWGS\n(Combined Model)"),
    # enforce desired order
    Timepoint = factor(Timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
  ) 

sens_df <- sens_df %>%
  mutate(Assay = factor(Assay,
                        levels = c("cfWGS\n(Sites Model)",
                                   "cfWGS\n(Combined Model)", 
                                   "clonoSEQ", "MFC")))

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
    values = custom_cols[c("Post-ASCT", "Maintenance-1yr")],  # enforce mapping
    limits = c("Post-ASCT", "Maintenance-1yr")                # enforce order
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of cfDNA and Clinical MRD Assays in Relapsing Patients",
    x     = "Technology",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_sens_blood
# ─────────────────────────────────────────────────────────────
# 3.  Save
# ─────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Supp_8A_Fig_sensitivity_by_tech_training_blood2.png",
  plot     = p_sens_blood,
  width    = 6,
  height   = 4,
  dpi      = 500
)




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
df_km_fl <- df_km %>%
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
  data = df_km_fl
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km_fl),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high


## check clonoSEQ
# fit the Kaplan–Meier curve
df_km_clonoSEQ <- df_km %>%
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
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km_clonoSEQ),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high

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

## --- Ns (head-to-head subset: requires cfWGS present) -------------------
n_cfWGS <- df_km %>%
  filter(!is.na(BM_zscore_only_detection_rate_call)) %>%
  distinct(Patient) %>%
  nrow()

n_MFC <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>%
  nrow()

n_clonoSEQ <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>%
  nrow()

## 6. Draft paragraph ----
paragraph <- glue(
  "After one year of maintenance therapy, among patients with cfWGS available (n={n_cfWGS}), ",
  "BM-cfWGS MRD-negative patients had {round(rfs_neg_cf)}% relapse-free survival at 24 months ",
  "versus {round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {ifelse(is.na(med_neg_cf),'NR',round(med_neg_cf,1))} vs ",
  "{ifelse(is.na(med_pos_cf),'NR',round(med_pos_cf,1))} months. ",
  "For MFC (n={n_MFC}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months ",
  "versus {round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {ifelse(is.na(med_neg_fl),'NR',round(med_neg_fl,1))} vs ",
  "{ifelse(is.na(med_pos_fl),'NR',round(med_pos_fl,1))} months. ",
  "clonoSEQ (n={n_clonoSEQ}) showed a similar direction of effect ",
  "(HR = {round(hr_clonoSEQ,2)}; 95% CI [{round(ci_lo_clonoSEQ,2)}–{round(ci_hi_clonoSEQ,2)}]). ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), comparable to flow ",
  "(ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients in the head-to-head subset, the minimum ",
  "detectable HR for 80% power is {round(hr80,1)} (and power to detect HR = 2.0 is {round(pw2*100,1)}%). ",
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
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS,
  N_MFC          = n_MFC,
  N_clonoSEQ     = n_clonoSEQ
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
df_km_fl <- df_km %>%
  filter(
    timepoint_info == "post_transplant",
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
  data = df_km_fl
)
med_fl <- surv_median(fit_fl)$median      # vector of two values
med_neg_fl <- med_fl[1] / 30.44           # convert days→months
med_pos_fl <- med_fl[2] / 30.44

cox_fl <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Flow_Binary,
        data = df_km_fl),
  exponentiate = TRUE, conf.int = TRUE
)
hr_fl      <- cox_fl$estimate
ci_lo_fl   <- cox_fl$conf.low
ci_hi_fl   <- cox_fl$conf.high

## Check clonoSEQ 
# fit the Kaplan–Meier curve
df_km_clonoSEQ <- df_km %>%
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
  data = df_km_clonoSEQ
)
med_clonoSEQ <- surv_median(fit_clonoSEQ)$median      # vector of two values
med_neg_clonoSEQ <- med_clonoSEQ[1] / 30.44           # convert days→months
med_pos_clonoSEQ <- med_clonoSEQ[2] / 30.44

cox_clonoSEQ <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Adaptive_Binary,
        data = df_km_clonoSEQ),
  exponentiate = TRUE, conf.int = TRUE
)
hr_clonoSEQ      <- cox_clonoSEQ$estimate
ci_lo_clonoSEQ   <- cox_clonoSEQ$conf.low
ci_hi_clonoSEQ   <- cox_clonoSEQ$conf.high



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

## Get N
n_cfWGS <- df_km %>%
  filter(!is.na(BM_zscore_only_detection_rate_call)) %>%
  distinct(Patient) %>%
  nrow()

n_MFC <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>%
  nrow()

n_clonoSEQ <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>%
  nrow()

## 6. Draft paragraph ----
paragraph <- glue(
  "At post-transplant, BM-cfWGS (n={n_cfWGS}) MRD-negative patients had ",
  "{round(rfs_neg_cf)}% relapse-free survival at 24 months versus ",
  "{round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by BM-cfWGS was {round(med_neg_cf,1)} vs {round(med_pos_cf,1)} months. ",
  "For MFC (n={n_MFC}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months versus ",
  "{round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {round(med_neg_fl,1)} vs {round(med_pos_fl,1)} months. ",
  "clonoSEQ (n={n_clonoSEQ}) showed a similar direction of effect, ",
  "though with greater uncertainty due to fewer patients. ",
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
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS,
  N_MFC          = n_MFC,
  N_clonoSEQ     = n_clonoSEQ
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
# fit the Kaplan–Meier curve
df_km_fl <- df_km %>%
  filter(
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

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
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
  data = df_km_clonoSEQ
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

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_zscore_only_sites_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

## Blood paragraph with Ns (and NR-safe medians)
paragraph_blood <- glue(
  "After one year of maintenance therapy, among patients with Blood-cfWGS available (n={n_cfWGS_blood}), ",
  "Blood-cfWGS MRD-negative patients had {round(rfs_neg_cf)}% relapse-free survival at 24 months ",
  "versus {round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by Blood-cfWGS was {ifelse(is.na(med_neg_cf),'NR',round(med_neg_cf,1))} vs ",
  "{ifelse(is.na(med_pos_cf),'NR',round(med_pos_cf,1))} months. ",
  "For MFC (n={n_MFC_blood}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months ",
  "versus {round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {ifelse(is.na(med_neg_fl),'NR',round(med_neg_fl,1))} vs ",
  "{ifelse(is.na(med_pos_fl),'NR',round(med_pos_fl,1))} months. ",
  "clonoSEQ (n={n_clonoSEQ_blood}) showed a similar direction of effect ",
  "(HR = {round(hr_clonoSEQ,2)}; 95% CI [{round(ci_lo_clonoSEQ,2)}–{round(ci_hi_clonoSEQ,2)}]). ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), comparable to flow ",
  "(ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients in the head-to-head subset, the minimum ",
  "detectable HR for 80% power is {round(hr80,1)} (and power to detect HR = {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph_blood)


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
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood
)



### Redo for post-transplant 
## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
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
# fit the Kaplan–Meier curve
df_km_fl <- df_km %>%
  filter(
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

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
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
  data = df_km_clonoSEQ
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

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_zscore_only_sites_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

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
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood
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




#### Now redo with the combined model for blood 

## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "1yr maintenance",
    !is.na(Blood_plus_fragment_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
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
df_km_fl <- df_km %>%
  filter(
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

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
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
  data = df_km_clonoSEQ
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
prob_var <- "Blood_plus_fragment_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_plus_fragment_call==1)
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

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_plus_fragment_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

## Blood paragraph with Ns (and NR-safe medians)
paragraph_blood <- glue(
  "After one year of maintenance therapy, among patients with Blood-cfWGS available (n={n_cfWGS_blood}), ",
  "Blood-cfWGS MRD-negative patients had {round(rfs_neg_cf)}% relapse-free survival at 24 months ",
  "versus {round(rfs_pos_cf)}% for MRD-positive patients ",
  "(HR = {round(hr_cf,2)}; 95% CI [{round(ci_lo_cf,2)}–{round(ci_hi_cf,2)}]). ",
  "Median RFS by Blood-cfWGS was {ifelse(is.na(med_neg_cf),'NR',round(med_neg_cf,1))} vs ",
  "{ifelse(is.na(med_pos_cf),'NR',round(med_pos_cf,1))} months. ",
  "For MFC (n={n_MFC_blood}), MRD-negative patients had {round(rfs_neg_fl)}% RFS at 24 months ",
  "versus {round(rfs_pos_fl)}% for MRD-positive patients ",
  "(HR = {round(hr_fl,2)}; 95% CI [{round(ci_lo_fl,2)}–{round(ci_hi_fl,2)}]), ",
  "with median RFS of {ifelse(is.na(med_neg_fl),'NR',round(med_neg_fl,1))} vs ",
  "{ifelse(is.na(med_pos_fl),'NR',round(med_pos_fl,1))} months. ",
  "clonoSEQ (n={n_clonoSEQ_blood}) showed a similar direction of effect ",
  "(HR = {round(hr_clonoSEQ,2)}; 95% CI [{round(ci_lo_clonoSEQ,2)}–{round(ci_hi_clonoSEQ,2)}]). ",
  "We also examined continuous MRD levels (model probability) against time-to-relapse ",
  "and found Spearman’s ρ = {round(rho1,2)} (p = {signif(p1,2)}), comparable to flow ",
  "(ρ = {round(rho2,2)}; p = {signif(p2,2)}). ",
  "With only {d} events among {nrow(df_km)} patients in the head-to-head subset, the minimum ",
  "detectable HR for 80% power is {round(hr80,1)} (and power to detect HR = {round(pw2*100,1)}%). ",
  "Accordingly, these analyses are presented as descriptive, hypothesis-generating results."
)

writeLines(paragraph_blood)


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
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood
)



### Redo for post-transplant 
## 1. Subset to post-transplant & Blood-cfWGS tested ----
df_km <- survival_df %>%
  filter(
    timepoint_info == "post_transplant",
    !is.na(Blood_plus_fragment_call)
  )

## 2. 24-month RFS by cfWGS BM ----
fit_cf <- survfit(
  Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
  data = df_km
)
# survival probabilities at 24 months:
t24    <- 24 * 30.44
sum_cf <- summary(fit_cf, times = t24)

# extract: strata 1 = negative, 2 = positive
rfs_neg_cf <- sum_cf$surv[1] * 100
rfs_pos_cf <- sum_cf$surv[2] * 100

cox_cf <- tidy(
  coxph(Surv(Time_to_event, Relapsed_Binary) ~ Blood_plus_fragment_call,
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
df_km_fl <- df_km %>%
  filter(
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

# Now clonoSEQ
df_km_clonoSEQ <- df_km %>%
  filter(
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
  data = df_km_clonoSEQ
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
prob_var <- "Blood_plus_fragment_prob"  
ct1 <- cor.test(df_km[[prob_var]], df_km$Time_to_event, method = "spearman")
rho1 <- ct1$estimate; p1 <- ct1$p.value

ct2 <- cor.test(df_km$Flow_pct_cells, df_km$Time_to_event, method = "spearman")
rho2 <- ct2$estimate; p2 <- ct2$p.value

## 5. Power diagnostics ----
d      <- sum(df_km$Relapsed_Binary)
prop_p <- mean(df_km$Blood_plus_fragment_call==1)
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

## Ns for BLOOD head-to-head subset (requires Blood cfWGS present)
n_cfWGS_blood <- df_km %>%
  filter(!is.na(Blood_plus_fragment_call)) %>%
  distinct(Patient) %>% nrow()

n_MFC_blood <- df_km %>%
  filter(!is.na(Flow_Binary)) %>%
  distinct(Patient) %>% nrow()

n_clonoSEQ_blood <- df_km %>%
  filter(!is.na(Adaptive_Binary)) %>%
  distinct(Patient) %>% nrow()

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
  Power_HR2_pct  = pw2 * 100,
  N_cfWGS        = n_cfWGS_blood,
  N_MFC          = n_MFC_blood,
  N_clonoSEQ     = n_clonoSEQ_blood
)

metrics_1yr <- metrics_1yr %>%
  mutate(Landmark = "1yr_maintenance") %>%
  select(Landmark, everything())  # move Landmark to front

progression_metrics_blood_combined <- bind_rows(
  metrics_post_transplant,
  metrics_1yr
)

# 8) Export summary table --------------------------------------------
write_csv(
  progression_metrics_blood_combined,
  file.path(outdir, "cfWGS_vs_flow_progression_summary_blood_muts_combined_model.csv")
)

# (Optional) also save as RDS for later use
saveRDS(
  progression_metrics_blood_combined,
  file.path(outdir, "cfWGS_vs_flow_progression_summary_blood_muts_updated_combined_model.rds")
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
                   cf = "cfWGS (Sites Model)",
                   fl = "MFC",
                   seq = "clonoSEQ"),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post‑ASCT", "Maintenance-1yr"))
  )

# reshape combined model (cfWGS only)
hr_plot_df_combined <- progression_metrics_blood_combined %>%
  select(Landmark,
         HR_cf, CI_low_cf, CI_high_cf) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf)"
  ) %>%
  mutate(
    Assay = "cfWGS (Combined Model)",
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post‑ASCT", "Maintenance-1yr"))
  )


# bind together
hr_plot_df <- bind_rows(hr_plot_df, hr_plot_df_combined)

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
    limits       = c(0.05, 34),   # extend upper limit
    breaks       = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32),
    minor_breaks = c(
      0.06, 0.08,      # between 0.05 & 0.1
      0.15, 0.2,       # between 0.1 & 0.25
      0.3, 0.4,        # between 0.25 & 0.5
      0.6, 0.8,        # between 0.5 & 1
      1.5,             # between 1 & 2
      3, 6,            # between 2 & 4 & 8
      12, 24           # new minor breaks between 8–16 and 16–32
    ),
    labels = label_number(accuracy = .01)
  )+
  annotation_logticks(
    sides  = "b",
    short  = unit(2, "pt"),
    mid    = unit(4, "pt"),
    long   = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    name   = NULL,
    values = c("cfWGS (Sites Model)" = "#35608DFF",
               "cfWGS (Combined Model)" = "#440154FF",
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

ggsave("Final Tables and Figures/SuppFig8B_cfWGS_blood_HR_updated3.png",
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

p_hr_bm

ggsave("Final Tables and Figures/Supp_Figure_6B_cfWGS_BM_HR_updated3.png",
       p_hr_bm, width = 6, height = 4, dpi = 600)






### Now make time to relapse figure 
df <- survival_df %>%                           # <- your tibble
  # keep samples beyond baseline / diagnosis
  filter(!str_detect(timepoint_info, regex("Diagnosis|Baseline", TRUE))) %>%
  
  # drop rows with missing probability or time
  filter(!is.na(BM_zscore_only_detection_rate_prob),
         !is.na(Time_to_event)) %>%
  
  # enforce non-negative time to event for relapse event visits a few days off from CMRG date
  mutate(
    days_before_event = pmax(Time_to_event, 0), # set negative values to 0 
    mrd_status      = factor(
      BM_zscore_only_detection_rate_call,
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
youden_thresh <- 0.4215524
#youden_thresh2 <- 0.35

max_mo <- max(df$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df, aes(months_before_event, BM_zscore_only_detection_rate_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
   # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status),
  shape  = 21,
  colour = "black",
  size   = 2
  ) +
  
  # 4) event line
  #geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("cVAF Model Probability",
                     limits = c(0,1),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  # scale_colour_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  # scale_fill_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  
  # # 7) stroke scale for MRD call
  # scale_discrete_manual(
  #   aesthetics = "stroke",
  #   values     = c("MRD-" = 0, "MRD+" = 1),
  #   guide      = guide_legend(
  #     title = "MRD call",
  #     override.aes = list(
  #       shape  = 21,
  #       fill   = "white",   # white interior in legend, makes stroke obvious
  #       size   = 4,
  #       colour = "black",
  #       stroke = c(0, 1)
  #     )
  #   )
  # ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "Longitudinal cfWGS MRD Probability by Patient Outcome\nUsing BM-Derived Mutation Lists"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "bottom", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 10)    # even smaller
  )

print(p_prob)

## Add the samples 
# Make sure the outcome labels match your scales
df <- df %>% mutate(progress_status = factor(progress_status,
                                             levels = c("No relapse","Relapse")))

n_patients_by <- df %>%
  distinct(Patient, progress_status) %>%
  count(progress_status, name = "n_patients")

n_timepoints_by <- df %>%
  count(progress_status, name = "n_timepoints")

# pull counts (0 if a group is absent)
get_n <- function(tbl, lvl, col) {
  val <- tbl %>% filter(progress_status == lvl) %>% pull({{col}})
  if (length(val) == 0) 0 else val
}

n_pat_nr  <- get_n(n_patients_by,  "No relapse", n_patients)
n_time_nr <- get_n(n_timepoints_by, "No relapse", n_timepoints)
n_pat_rl  <- get_n(n_patients_by,  "Relapse",    n_patients)
n_time_rl <- get_n(n_timepoints_by, "Relapse",   n_timepoints)

# text to print
lab_nr <- paste0("No relapse: n=", n_pat_nr, " patients; ", n_time_nr, " samples")
lab_rl <- paste0("Relapse: n=", n_pat_rl, " patients; ", n_time_rl, " samples")

# place labels at bottom-left (remember: x is reversed, so 'left' == large x)
x_left <- max_mo - 0.02 * max_mo  # a small inset from the left border

p_prob2 <- p_prob +
  scale_colour_manual(
    name = "Patient outcome",
#    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
    values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  ) +
  scale_fill_manual(
    name = "Patient outcome",
 #   values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
   values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  )

print(p_prob2)

# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated5.png",
       p_prob, width = 6, height = 4.5, dpi = 600)

ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated5_label.png",
       p_prob2, width = 6, height = 4.5, dpi = 600)



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

df_plot <- df_plot %>%
  mutate(Time_to_event = if_else(Time_to_event >= -30 & Time_to_event < 0, 0, Time_to_event))

## ─────────────────────────────────────────────────────────────
## 1)  Build the scatter plot                                  
## ─────────────────────────────────────────────────────────────
p_time <- ggplot(df_plot,
                 aes(x = BM_zscore_only_detection_rate_prob,
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

print(p_time)

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob2.png"),
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
overflow <- max_days + 180   # a little beyond the longest follow‑up

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
  rel_df$BM_zscore_only_detection_rate_prob,
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


# Edit days
plot_df2 <- plot_df2 %>%
  mutate(days_plot = if_else(days_plot >= -30 & days_plot <= 0, 0, days_plot))

# 2) compute Spearman rho on the relapsers
spearman_res <- with(
  filter(plot_df2, progress_status == "Relapse"),
  cor.test(BM_zscore_only_detection_rate_prob,
           days_before_event,
           method = "spearman")
)

## Check other metrics
# 2A) pull out estimate + p‑value
rho_all  <- spearman_res$estimate
pval_all <- spearman_res$p.value

# B) excluding relapse samples (days_plot > 0 only)
spearman_pre <- with(
  filter(plot_df2, progress_status == "Relapse", days_plot > 0),
  cor.test(BM_zscore_only_detection_rate_prob,
           days_before_event,
           method = "spearman")
)

rho_pre  <- spearman_pre$estimate
pval_pre <- spearman_pre$p.value

pval_all_str <- ifelse(pval_all < 0.001, "<0.001", sprintf("%.3f", pval_all))
pval_pre_str <- ifelse(pval_pre < 0.001, "<0.001", sprintf("%.3f", pval_pre))

annot_text <- sprintf("All relapse samples:\nρ=%.2f, p=%s\nPre-relapse only:\nρ=%.2f, p=%s",
                      rho_all, pval_all_str,
                      rho_pre, pval_pre_str)

# 3) make the scatter
p_time_inf <- ggplot(plot_df2,
                     aes(x = BM_zscore_only_detection_rate_prob,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = 0.35, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0,1),
    breaks = seq(0,1,by=0.2),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # annotate Spearman
  # annotate("text",
  #          x = 0.01,      # left margin
  #          y = 1650,
  #          label = sprintf("ρ = %.2f\np = %s", rho, pval_str),
  #          hjust = 0,
  #          size = 3.5) +
  annotate("text",
           x = 0.01,
           y = 1450,
           label = annot_text,
           hjust = 0,
           size = 3.5) +

  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using BM-Derived Mutations")

print(p_time_inf)

# 4) render / save
ggsave(file.path("Final Tables and Figures/Fig_4D_time_to_relapse_infinity_no_reverse2_BM_muts_updated3.png"),
       p_time_inf, width = 5.5, height = 4.5, dpi = 600)


## Try different legend layout 
library(cowplot)    # get_legend()
library(patchwork)  # easy assembly

# --- 1) build the main plot *without* a legend (legend handled below) ---
p_main <- ggplot(plot_df2,
                 aes(x = BM_zscore_only_detection_rate_prob, y = days_plot)) +
  geom_vline(xintercept = 0.35, linetype = "dotted", colour = "grey40") +
  geom_point(aes(colour = progress_status, fill = progress_status),
             shape = 21, size = 3, colour = "black") +
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0.1,1),
    breaks = seq(0.1,1,by=0.2),
    labels = scales::percent_format(accuracy=1)
  )  +# colours
scale_colour_manual(
  "Patient outcome",
  values = c("No relapse" = "black", "Relapse" = "red")
) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  guides(fill = "none") +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none"        # <- hide here; we'll place it below
  ) +
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using BM-Derived Mutations")

# --- 2) LEGEND-ONLY PLOT (do NOT inherit guides(fill = 'none')) ---
# -----build a dummy legend that always shows both levels -----
p_legend_only <- ggplot(legend_df, aes(x, y, fill = progress_status)) +
  # make the plotting layer invisible *in the panel* …
  geom_point(shape = 21, size = 0, colour = "black", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = pal_vals,
    breaks = names(pal_vals)
  ) +
  guides(fill = guide_legend(
    ncol = 1,                 # <- stack items vertically
    byrow = TRUE,
    title.position = "top",   # title on its own line
    label.hjust = 0,          # left-align labels
    override.aes = list(shape = 21, size = 3, alpha = 1, colour = "black")
  )) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "vertical",         # <- vertical legend
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(6, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  )

# --- 3) make two small “text boxes” for the right columns ---
txt_all <- sprintf("All relapse samples:\nρ=%.2f, p=%s", rho_all, pval_all_str)
txt_pre <- sprintf("Pre-relapse only:\nρ=%.2f, p=%s", rho_pre, pval_pre_str)

mini_box <- function(s) {
  ggplot() +
    annotate("label", x = 0, y = 1, label = s,
             hjust = 0, vjust = 1, size = 3.5,
             label.size = 0, fill = scales::alpha("white", 0.7)) +
    xlim(0,1) + ylim(0,1) +
    theme_void()
}

col2 <- mini_box(txt_all)
col3 <- mini_box(txt_pre)

# --- 4) assemble: plot on top; 3 columns underneath ---
bottom_row <- p_legend_only | col2 | col3
final_plot <- p_main / bottom_row + plot_layout(heights = c(1, 0.22))

# show and save
print(final_plot)
ggsave("Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_BM_muts.png",
       final_plot, width = 5.5, height = 5.5, dpi = 600)



plot_df2 %>%
  filter(
    is.na(BM_zscore_only_detection_rate_prob) | is.na(days_plot) |
      BM_zscore_only_detection_rate_prob < 0 | BM_zscore_only_detection_rate_prob > 1 |
      days_plot < 0 | days_plot > overflow
  )

time_to_relapse_BM <- plot_df2



### Redo now for blood derived muts 
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
youden_thresh <- 0.5166693
youden_thr <- youden_thresh # for consistency
max_mo <- max(df_slim$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df, aes(months_before_event, Blood_zscore_only_sites_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
  
  # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status),
    shape  = 21,
    colour = "black",
    size   = 2
  ) +
  
  # 4) event line
  #geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("Sites Model Probability",
                     limits = c(0.3,1),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  # scale_colour_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  # scale_fill_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  # # 7) stroke scale for MRD call
  # scale_discrete_manual(
  #   aesthetics = "stroke",
  #   values     = c("MRD-" = 0, "MRD+" = 1),
  #   guide      = guide_legend(
  #     title = "MRD call",
  #     override.aes = list(
  #       shape  = 21,
  #       fill   = "white",   # white interior in legend, makes stroke obvious
  #       size   = 4,
  #       colour = "black",
  #       stroke = c(0, 1)
  #     )
  #   )
  # ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "Longitudinal cfWGS MRD Probability by Patient Outcome\nUsing cfDNA-Derived Mutation Lists"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "bottom", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 10)    # even smaller
  )

print(p_prob)

### Add label 
# Make sure the outcome labels match your scales
df <- df %>% mutate(progress_status = factor(progress_status,
                                             levels = c("No relapse","Relapse")))

n_patients_by <- df %>%
  distinct(Patient, progress_status) %>%
  count(progress_status, name = "n_patients")

n_timepoints_by <- df %>%
  count(progress_status, name = "n_timepoints")

# pull counts (0 if a group is absent)
get_n <- function(tbl, lvl, col) {
  val <- tbl %>% filter(progress_status == lvl) %>% pull({{col}})
  if (length(val) == 0) 0 else val
}

n_pat_nr  <- get_n(n_patients_by,  "No relapse", n_patients)
n_time_nr <- get_n(n_timepoints_by, "No relapse", n_timepoints)
n_pat_rl  <- get_n(n_patients_by,  "Relapse",    n_patients)
n_time_rl <- get_n(n_timepoints_by, "Relapse",   n_timepoints)

# text to print
lab_nr <- paste0("No relapse: n=", n_pat_nr, " patients; ", n_time_nr, " samples")
lab_rl <- paste0("Relapse: n=", n_pat_rl, " patients; ", n_time_rl, " samples")

# place labels at bottom-left (remember: x is reversed, so 'left' == large x)
x_left <- max_mo - 0.02 * max_mo  # a small inset from the left border

p_prob2 <- p_prob +
  scale_colour_manual(
    name = "Patient outcome",
#    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
    values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  ) +
  scale_fill_manual(
    name = "Patient outcome",
   # values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
   values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  )



print(p_prob2)
# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated3_blood4.png",
       p_prob, width = 6, height = 4.5, dpi = 600)

ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated3_blood2_labelled3.png",
       p_prob2, width = 6, height = 4.5, dpi = 600)


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

df_plot <- df_plot %>%
  mutate(Time_to_event = if_else(Time_to_event >= -30 & Time_to_event < 0, 0, Time_to_event))

## ─────────────────────────────────────────────────────────────
## 1)  Build the scatter plot                                  
## ─────────────────────────────────────────────────────────────
p_time <- ggplot(df_plot,
                 aes(x = Blood_zscore_only_sites_prob,
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

print(p_time)

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob2.png"),
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
overflow <- max_days + 180   # a little beyond the longest follow‑up

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
  rel_df$Blood_zscore_only_sites_prob,
  rel_df$days_before_event,
  method = "spearman"
)

# 3) print them
cat("\n--- cfWGS model probability vs days ---\n")
print(spearman_prob)


# Edit days
plot_df2 <- plot_df2 %>%
  mutate(days_plot = if_else(days_plot >= -35 & days_plot <= 0, 0, days_plot))

# 2) compute Spearman rho on the relapsers
spearman_res <- with(
  filter(plot_df2, progress_status == "Relapse"),
  cor.test(Blood_zscore_only_sites_prob,
           days_before_event,
           method = "spearman")
)

## Check other metrics
# 2A) pull out estimate + p‑value
# A) including relapse samples (your current spearman_res)
rho_all  <- spearman_res$estimate
pval_all <- spearman_res$p.value

# B) excluding relapse samples (days_plot > 0 only)
spearman_pre <- with(
  filter(plot_df2, progress_status == "Relapse", days_plot > 0),
  cor.test(Blood_zscore_only_sites_prob,
           days_before_event,
           method = "spearman")
)

rho_pre  <- spearman_pre$estimate
pval_pre <- spearman_pre$p.value

pval_all_str <- ifelse(pval_all < 0.001, "<0.001", sprintf("%.3f", pval_all))
pval_pre_str <- ifelse(pval_pre < 0.001, "<0.001", sprintf("%.3f", pval_pre))

annot_text <- sprintf("All relapse samples:\nρ=%.2f, p=%s\nPre-relapse only:\nρ=%.2f, p=%s",
                      rho_all, pval_all_str,
                      rho_pre, pval_pre_str)

## See range
plot_df2 %>%
  summarise(
    min   = min(Blood_zscore_only_sites_prob, na.rm = TRUE),
    max   = max(Blood_zscore_only_sites_prob, na.rm = TRUE),
    range = max(Blood_zscore_only_sites_prob, na.rm = TRUE) -
      min(Blood_zscore_only_sites_prob, na.rm = TRUE)
  )

# 3) make the scatter
p_time_inf <- ggplot(plot_df2,
                     aes(x = Blood_zscore_only_sites_prob,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = youden_thr, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability (%)",
    limits = c(0.3,1),
    breaks = seq(0,1,by=0.1),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # annotate Spearman
  # annotate("text",
  #          x = 0.90,      # left margin
  #          y = 1650,
  #          label = sprintf("ρ = %.2f\np = %s", rho, pval_fmt),
  #          hjust = 0,
  #          size = 3.5) +
  annotate("text",
           x = 0.8,
           y = 1450,
           label = annot_text,
           hjust = 0,
           size = 3.5) +
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using cfDNA-Derived Mutations")

print(p_time_inf)

# 4) render / save
ggsave(file.path("Final Tables and Figures/Fig_4_D_time_to_relapse_infinity_no_reverse2_blood_muts4_small2.png"),
       p_time_inf, width = 5.5, height = 4.5, dpi = 600)


time_to_relapse_blood <- plot_df2

## Check what is not included
plot_df2 %>%
  filter(
    is.na(Blood_zscore_only_sites_prob) | is.na(days_plot) |
      Blood_zscore_only_sites_prob < 0.3 | Blood_zscore_only_sites_prob > 1 |
      days_plot < 0 | days_plot > overflow
  )



### Do different legend layout 

## Try different legend layout 

# --- 1) build the main plot *without* a legend (legend handled below) ---
p_main <-  ggplot(plot_df2,
                  aes(x = Blood_zscore_only_sites_prob,
                      y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = youden_thr, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # ∞‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or ∞ for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "∞")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability (%)",
    limits = c(0.3,1),
    breaks = seq(0,1,by=0.1),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using cfDNA-Derived Mutations")

# --- 2) LEGEND-ONLY PLOT (do NOT inherit guides(fill = 'none')) ---
# -----build a dummy legend that always shows both levels -----
p_legend_only <- ggplot(legend_df, aes(x, y, fill = progress_status)) +
  # make the plotting layer invisible *in the panel* …
  geom_point(shape = 21, size = 0, colour = "black", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = pal_vals,
    breaks = names(pal_vals)
  ) +
  guides(fill = guide_legend(
    ncol = 1,                 # <- stack items vertically
    byrow = TRUE,
    title.position = "top",   # title on its own line
    label.hjust = 0,          # left-align labels
    override.aes = list(shape = 21, size = 3, alpha = 1, colour = "black")
  )) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "vertical",         # <- vertical legend
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(6, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  )

# --- 3) make two small “text boxes” for the right columns ---
txt_all <- sprintf("All relapse samples:\nρ=%.2f, p=%s", rho_all, pval_all_str)
txt_pre <- sprintf("Pre-relapse only:\nρ=%.2f, p=%s", rho_pre, pval_pre_str)

mini_box <- function(s) {
  ggplot() +
    annotate("label", x = 0, y = 1, label = s,
             hjust = 0, vjust = 1, size = 3.5,
             label.size = 0, fill = scales::alpha("white", 0.7)) +
    xlim(0,1) + ylim(0,1) +
    theme_void()
}

col2 <- mini_box(txt_all)
col3 <- mini_box(txt_pre)

# --- 4) assemble: plot on top; 3 columns underneath ---
bottom_row <- p_legend_only | col2 | col3
final_plot <- p_main / bottom_row + plot_layout(heights = c(1, 0.22))

# show and save
print(final_plot)
ggsave("Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_blood_muts.png",
       final_plot, width = 5.5, height = 5.5, dpi = 600)













### Now evaluate on non-frontline too 
## Re-introdue dat since previously limited 
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )


## Do rescored 
dat <- dat %>%
  ## Add the screen column 
  mutate(
    BM_zscore_only_detection_rate_screen_call  = as.integer(BM_zscore_only_detection_rate_prob >= 0.350),
  )

################################################################################
##  Time-window prediction performance in Non-frontline cohort
################################################################################

# 1) Your assays vector
assays <- c(
  Flow         = "Flow_Binary",
  cfWGS_BM     = "BM_zscore_only_detection_rate_call",
  cfWGS_Blood  = "Blood_zscore_only_sites_call", 
  cfWGS_Blood_Combined = "Blood_plus_fragment_call"
)

## Get additional dates
Relapse_dates_full <- read_csv(
  "Exported_data_tables_clinical/Relapse dates cfWGS updated2.csv",
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
    Blood_zscore_only_sites_call,
    Blood_plus_fragment_call
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
    dplyr::filter(!is.na(.data[[col_name]])) %>%
    dplyr::mutate(
      test_pos        = (.data[[col_name]] == 1),
      event_in_window = relapsed == 1 &
        !is.na(relapse_date) &
        relapse_date <= sample_date + lubridate::days(win_d)   # << here
    )
  
  tab <- table(
    factor(df2$test_pos,        levels = c(FALSE, TRUE)),
    factor(df2$event_in_window, levels = c(FALSE, TRUE))
  )
  
  tp <- tab["TRUE","TRUE"]; fn <- tab["FALSE","TRUE"]
  fp <- tab["TRUE","FALSE"]; tn <- tab["FALSE","FALSE"]
  
  tibble::tibble(
    Window_days = win_d,
    Assay       = assay_label,
    N_samples   = nrow(df2),
    N_patients  = dplyr::n_distinct(df2$Patient),
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
# --- 0) use the refined table everywhere ---------------------------------
bm_col    <- assays["cfWGS_BM"]      # "BM_zscore_only_detection_rate_call"
blood_col <- assays["cfWGS_Blood"]   # "Blood_zscore_only_sites_call"

# Rebuild subsets FROM df_sf2 (not df_sf)
df_sf_BM    <- df_sf2 %>% dplyr::filter(!is.na(.data[[bm_col]]))
df_sf_blood <- df_sf2 %>% dplyr::filter(!is.na(.data[[blood_col]]))

# Recompute results on the consistent base
results_BM <- purrr::map_dfr(windows, function(w) {
  purrr::map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_BM, a, assays[[a]], w)
  })
}) %>% dplyr::arrange(Window_days, dplyr::desc(Sensitivity), dplyr::desc(Specificity))

results_blood <- purrr::map_dfr(windows, function(w) {
  purrr::map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_blood, a, assays[[a]], w)
  })
}) %>% dplyr::arrange(Window_days, dplyr::desc(Sensitivity), dplyr::desc(Specificity))


# --- 1) narrative counts that match the EXACT subsets above ---------------
wins <- c(180, 365)

# Helper: per-window event counts on a given DF (sample-level)
count_relapses_by_window <- function(df, wins = c(180, 365)) {
  purrr::map_dfr(wins, function(w) {
    df %>%
      dplyr::mutate(event_in_window = relapsed == 1 &
                      !is.na(relapse_date) &
                      relapse_date <= sample_date + lubridate::days(w)) %>%
      dplyr::summarise(
        Window_days       = w,
        Samples_relapsed  = sum(event_in_window, na.rm = TRUE),
        Patients_relapsed = dplyr::n_distinct(Patient[event_in_window])
      )
  })
}

# A) BM-only wording (matches results_BM denominators for the cfWGS_BM rows)
summ_bm <- count_relapses_by_window(df_sf_BM, wins)
txt_bm <- glue::glue(
  "Among bone-marrow cfWGS samples, {summ_bm$Samples_relapsed[summ_bm$Window_days==180]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_bm$Samples_relapsed[summ_bm$Window_days==365]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_bm, "\n")

# B) Blood-only wording (matches results_blood denominators for the cfWGS_Blood rows)
summ_blood <- count_relapses_by_window(df_sf_blood, wins)
txt_blood <- glue::glue(
  "Among blood cfWGS samples, {summ_blood$Samples_relapsed[summ_blood$Window_days==180]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_blood$Samples_relapsed[summ_blood$Window_days==365]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_blood, "\n")

# C) “Any cfWGS” wording = union of rows that have BM OR Blood assays
cfwgs_any <- df_sf2 %>% dplyr::filter(!is.na(.data[[bm_col]]) | !is.na(.data[[blood_col]]))
summ_any  <- count_relapses_by_window(cfwgs_any, wins)
txt_any <- glue::glue(
  "In the test cohort, sampling times were heterogeneous, so we assessed each assay’s ",
  "ability to predict progression within fixed time windows (180 and 365 days). ",
  "Among cfWGS samples, {summ_any$Samples_relapsed[summ_any$Window_days==180]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_any$Samples_relapsed[summ_any$Window_days==365]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_any, "\n")



# 5a) Restrict to BM-cfWGS subset, then loop ------------------------------
bm_col     <- assays["cfWGS_BM"]    # "BM_zscore_only_detection_rate_call"
df_sf_BM   <- df_sf2 %>% filter(!is.na(.data[[bm_col]]))

results_BM <- map_dfr(windows, function(w) {
  map_dfr(names(assays), function(a) {
    calc_metrics(df_sf_BM, a, assays[[a]], w)
  })
}) %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity))

# 5b) Restrict to blood-cfWGS subset, then loop ---------------------------
blood_col   <- assays["cfWGS_Blood"]  # "Blood_zscore_only_sites_call"
df_sf_blood <- df_sf2 %>% filter(!is.na(.data[[blood_col]]))

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
  filter(Assay %in% c("cfWGS_Blood", "Flow", "cfWGS_Blood_Combined"), Window_days %in% c(180, 365)) %>%
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

## Get number of progressors
wins <- c(180, 365)

# Helper that, given a data frame of samples, returns counts per window
count_relapses_by_window <- function(df, wins = c(180, 365)) {
  purrr::map_dfr(wins, function(w) {
    df %>%
      dplyr::mutate(
        event_in_window = relapsed == 1 &
          !is.na(relapse_date) &
          relapse_date <= sample_date + lubridate::days(w)
      ) %>%
      dplyr::summarise(
        Window_days        = w,
        Samples_relapsed   = sum(event_in_window, na.rm = TRUE),
        Patients_relapsed  = dplyr::n_distinct(Patient[event_in_window])
      )
  })
}

# A) "cfWGS samples" = any sample with BM or Blood cfWGS available
cfwgs_any <- df_sf2 %>%
  dplyr::filter(!is.na(.data[[bm_col]]) | !is.na(.data[[blood_col]]))
summ_any <- count_relapses_by_window(cfwgs_any, wins)

txt_any <- glue::glue(
  "In the test cohort, sampling times were heterogeneous, so we assessed each assay’s ",
  "ability to predict progression within fixed time windows (180 and 365 days). ",
  "Among cfWGS samples, {summ_any$Samples_relapsed[summ_any$Window_days==180]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_any$Samples_relapsed[summ_any$Window_days==365]} ",
  "from {summ_any$Patients_relapsed[summ_any$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_any, "\n")

# (Optional) BM-only wording
cfwgs_bm   <- df_sf2 %>% dplyr::filter(!is.na(.data[[bm_col]]))
summ_bm    <- count_relapses_by_window(cfwgs_bm, wins)
txt_bm <- glue::glue(
  "Among bone-marrow cfWGS samples, {summ_bm$Samples_relapsed[summ_bm$Window_days==180]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_bm$Samples_relapsed[summ_bm$Window_days==365]} ",
  "from {summ_bm$Patients_relapsed[summ_bm$Window_days==365]} patients relapsed within 365 days."
)
cat(txt_bm, "\n")

# (Optional) Blood-only wording
cfwgs_blood <- df_sf2 %>% dplyr::filter(!is.na(.data[[blood_col]]))
summ_blood  <- count_relapses_by_window(cfwgs_blood, wins)
txt_blood <- glue::glue(
  "Among blood cfWGS samples, {summ_blood$Samples_relapsed[summ_blood$Window_days==180]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==180]} patients relapsed within 180 days ",
  "and {summ_blood$Samples_relapsed[summ_blood$Window_days==365]} ",
  "from {summ_blood$Patients_relapsed[summ_blood$Window_days==365]} patients relapsed within 365 days."
)
 cat(txt_blood, "\n")



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
    limits = c(0, 105),
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
ggsave("Final Tables and Figures/Supp_Fig_6_Fig_sensitivity_windows_BM_test_cohort_updated2.png",
       plot = p_sens_bm,
       width = 4.75, height = 6.25, dpi = 500)


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
ggsave("Final Tables and Figures/Supp_Fig_8_Fig_sensitivity_windows_blood_test_cohort3.png",
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
  file.path(outdir, "BM_cfWGS_timewindow_results2.csv")
)

# blood‐cfWGS subset
write_csv(
  results_blood,
  file.path(outdir, "blood_cfWGS_timewindow_results2.csv")
)

## export 
# --- Clean BM results ---
results_BM_clean <- results_BM %>%
  # remove blood assays
  filter(!Assay %in% c("cfWGS_Blood", "cfWGS_Blood_Combined")) %>%
  # rename Flow → MFC
  mutate(
    Assay = case_when(
      Assay == "Flow" ~ "MFC",
      TRUE ~ Assay
    )
  )

# --- Clean Blood results ---
results_blood_clean <- results_blood %>%
  # remove BM assays
  filter(!Assay %in% c("cfWGS_BM")) %>%
  # rename assays
  mutate(
    Assay = case_when(
      Assay == "Flow" ~ "MFC",
      Assay == "cfWGS_Blood" ~ "cfWGS_Blood (Sites Model)",
      Assay == "cfWGS_Blood_Combined" ~ "cfWGS_Blood (Combined Model)",
      TRUE ~ Assay
    )
  )

# --- Export to Excel ---
export_list <- list(
  "BM_models"    = results_BM_clean,
  "Blood_models" = results_blood_clean
)

# Write to Excel
write_xlsx(
  export_list,
  path = file.path("Final Tables and Figures/Supplementary_Table_9_timewindow_results_test_cohort.xlsx")
)

# event counts for BM‐cfWGS
write_csv(
  event_counts_BM,
  file.path(outdir, "BM_cfWGS_event_counts2.csv")
)

# For blood
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
  file.path(outdir, "blood_cfWGS_event_counts2.csv")
)




### See at what pont an increase occured 
d2m <- function(days, digits = 1) round(as.numeric(days) / 30.44, digits)

# 1) choose the probability column you want to analyze:
assay_prob <- "BM_zscore_only_detection_rate_prob"  # or "BM_zscore_only_detection_rate_prob"
# Optional: restrict to surveillance timepoints only (post-induction / ASCT / maintenance)
surveillance_only <- FALSE
surv_regex <- "(post[_ -]?induction|pre[_ -]?(asct|transplant)|post[_ -]?(asct|transplant)|maintenance)"

# ==== BASE (relapse patients only, consistent with figure) ====
df0 <- time_to_relapse_BM %>%
  mutate(timepoint_info = tolower(timepoint_info)) %>%
  filter(progress_status == "Relapse",
         !is.na(.data[[assay_prob]]))

if (surveillance_only) {
  df0 <- df0 %>% filter(grepl(surv_regex, timepoint_info))
}

# ==== A) STRICT pre-progression monitoring set ====
# Use days_before_event > 0 to enforce strictly pre-progression (same logic as sample_date < censor_date)
df_pre <- df0 %>%
  filter(days_before_event > 0) %>%
  arrange(Patient, sample_date)

n_patients <- n_distinct(df_pre$Patient)
n_samples  <- nrow(df_pre)

sample_timing_stats <- df_pre %>%
  summarise(
    median_days_before = median(days_before_event, na.rm = TRUE),
    iqr_days_before    = IQR(days_before_event,    na.rm = TRUE),
    min_days_before    = min(days_before_event,    na.rm = TRUE),
    max_days_before    = max(days_before_event,    na.rm = TRUE),
    .groups = "drop"
  )

sample_timing_stats <- sample_timing_stats %>%
  mutate(
    median_days_before_mo = d2m(median_days_before),
    iqr_days_before_mo    = d2m(iqr_days_before),
    min_days_before_mo    = d2m(min_days_before),
    max_days_before_mo    = d2m(max_days_before)
  )

# ==== B) Per-patient nadir & first increase (pre-progression only) ====
advance_df <- df0 %>%
  # Keep rows up to relapse (+30d tolerance) so patients with slight date mismatches are still present
  dplyr::filter(sample_date <= (censor_date + 30)) %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(sample_date, .by_group = TRUE) %>%
  dplyr::group_modify(~{
    df <- .
    prog_date <- df$censor_date[1]
    
    # STRICT pre-relapse rows for nadir / first-increase logic
    df_pre <- dplyr::filter(df, sample_date < prog_date)
    
    if (nrow(df_pre) == 0L) {
      return(tibble::tibble(
        nadir_date = as.Date(NA),  nadir_prob = NA_real_,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    # Nadir = minimum probability before relapse
    nadir_row <- df_pre %>%
      dplyr::slice_min(.data[[assay_prob]], with_ties = FALSE)
    
    nadir_date <- nadir_row$sample_date
    nadir_prob <- nadir_row[[assay_prob]]
    
    # First increase after nadir (strictly later in time AND strictly higher prob)
    post_nadir <- df_pre %>%
      dplyr::filter(sample_date > nadir_date,
                    .data[[assay_prob]] > nadir_prob) %>%
      dplyr::slice_head(n = 1)
    
    if (nrow(post_nadir) == 0L) {
      return(tibble::tibble(
        nadir_date = nadir_date,  nadir_prob = nadir_prob,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase  = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    first_inc_date <- post_nadir$sample_date
    first_inc_prob <- post_nadir[[assay_prob]]
    
    tibble::tibble(
      nadir_date = nadir_date,
      nadir_prob = nadir_prob,
      first_inc_date = first_inc_date,
      first_inc_prob = first_inc_prob,
      days_to_first_increase  = as.numeric(first_inc_date - nadir_date),
      days_before_progression = as.numeric(prog_date - first_inc_date),
      prog_date = prog_date
    )
  }) %>%
  dplyr::ungroup()


# Nadir timing relative to progression
nadir_timing_stats <- advance_df %>%
  mutate(days_nadir_before = as.numeric(prog_date - nadir_date)) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before,    na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before,    na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before,    na.rm = TRUE),
    .groups = "drop"
  )

nadir_timing_stats <- nadir_timing_stats %>%
  mutate(
    median_nadir_mo = d2m(median_nadir_days),
    iqr_nadir_mo    = d2m(iqr_nadir_days),
    min_nadir_mo    = d2m(min_nadir_days),
    max_nadir_mo    = d2m(max_nadir_days)
  )

# First increase timing (after nadir, and its lead time to progression)
increase_stats2 <- advance_df %>%
  summarise(
    median_after_nadir = median(days_to_first_increase,   na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,      na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,      na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,      na.rm = TRUE),
    median_before_prog = median(days_before_progression,  na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,     na.rm = TRUE),
    min_before_prog    = min(days_before_progression,     na.rm = TRUE),
    max_before_prog    = max(days_before_progression,     na.rm = TRUE),
    .groups            = "drop"
  )

increase_stats2 <- increase_stats2 %>%
  mutate(
    median_after_nadir_mo = d2m(median_after_nadir),
    iqr_after_nadir_mo    = d2m(iqr_after_nadir),
    min_after_nadir_mo    = d2m(min_after_nadir),
    max_after_nadir_mo    = d2m(max_after_nadir),
    median_before_prog_mo = d2m(median_before_prog),
    iqr_before_prog_mo    = d2m(iqr_before_prog),
    min_before_prog_mo    = d2m(min_before_prog),
    max_before_prog_mo    = d2m(max_before_prog)
  )

# ==== C) One-paragraph sentence  ====
cat(glue(
  "We analyzed {n_patients} patients (total {n_samples} samples) collected a median ",
  "{sample_timing_stats$median_days_before} days (~{sample_timing_stats$median_days_before_mo} mo) before progression ",
  "(IQR {sample_timing_stats$iqr_days_before} days, ~{sample_timing_stats$iqr_days_before_mo} mo; ",
  "range {sample_timing_stats$min_days_before}–{sample_timing_stats$max_days_before} days, ",
  "~{sample_timing_stats$min_days_before_mo}–{sample_timing_stats$max_days_before_mo} mo). ",
  "The nadir detection probability occurred a median ",
  "{nadir_timing_stats$median_nadir_days} days (~{nadir_timing_stats$median_nadir_mo} mo) before progression ",
  "(IQR {nadir_timing_stats$iqr_nadir_days} days, ~{nadir_timing_stats$iqr_nadir_mo} mo; ",
  "range {nadir_timing_stats$min_nadir_days}–{nadir_timing_stats$max_nadir_days} days, ",
  "~{nadir_timing_stats$min_nadir_mo}–{nadir_timing_stats$max_nadir_mo} mo). ",
  "From that nadir, the first increase occurred a median ",
  "{increase_stats2$median_after_nadir} days (~{increase_stats2$median_after_nadir_mo} mo) later ",
  "(IQR {increase_stats2$iqr_after_nadir} days, ~{increase_stats2$iqr_after_nadir_mo} mo; ",
  "range {increase_stats2$min_after_nadir}–{increase_stats2$max_after_nadir} days, ",
  "~{increase_stats2$min_after_nadir_mo}–{increase_stats2$max_after_nadir_mo} mo), ",
  "which was a median {increase_stats2$median_before_prog} days ",
  "(~{increase_stats2$median_before_prog_mo} mo) before IMWG-defined clinical progression ",
  "(IQR {increase_stats2$iqr_before_prog} days, ~{increase_stats2$iqr_before_prog_mo} mo; ",
  "range {increase_stats2$min_before_prog}–{increase_stats2$max_before_prog} days, ",
  "~{increase_stats2$min_before_prog_mo}–{increase_stats2$max_before_prog_mo} mo).\n"
))



##### Now redo for blood muts 
# 1) choose the probability column you want to analyze:
assay_prob <- "Blood_zscore_only_sites_call"  # or "BM_zscore_only_detection_rate_prob"
# Optional: restrict to surveillance timepoints only (post-induction / ASCT / maintenance)
surveillance_only <- FALSE
surv_regex <- "(post[_ -]?induction|pre[_ -]?(asct|transplant)|post[_ -]?(asct|transplant)|maintenance)"

# ==== BASE (relapse patients only, consistent with figure) ====
df0 <- time_to_relapse_blood %>%
  mutate(timepoint_info = tolower(timepoint_info)) %>%
  filter(progress_status == "Relapse",
         !is.na(.data[[assay_prob]]))

if (surveillance_only) {
  df0 <- df0 %>% filter(grepl(surv_regex, timepoint_info))
}

# ==== A) STRICT pre-progression monitoring set ====
# Use days_before_event > 0 to enforce strictly pre-progression (same logic as sample_date < censor_date)
df_pre <- df0 %>%
  filter(days_before_event > 0) %>%
  arrange(Patient, sample_date)

n_patients <- n_distinct(df_pre$Patient)
n_samples  <- nrow(df_pre)

sample_timing_stats <- df_pre %>%
  summarise(
    median_days_before = median(days_before_event, na.rm = TRUE),
    iqr_days_before    = IQR(days_before_event,    na.rm = TRUE),
    min_days_before    = min(days_before_event,    na.rm = TRUE),
    max_days_before    = max(days_before_event,    na.rm = TRUE),
    .groups = "drop"
  )

sample_timing_stats <- sample_timing_stats %>%
  mutate(
    median_days_before_mo = d2m(median_days_before),
    iqr_days_before_mo    = d2m(iqr_days_before),
    min_days_before_mo    = d2m(min_days_before),
    max_days_before_mo    = d2m(max_days_before)
  )

# ==== B) Per-patient nadir & first increase (pre-progression only) ====
advance_df <- df0 %>%
  # Keep rows up to relapse (+30d tolerance) so patients with slight date mismatches are still present
  dplyr::filter(sample_date <= (censor_date + 30)) %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(sample_date, .by_group = TRUE) %>%
  dplyr::group_modify(~{
    df <- .
    prog_date <- df$censor_date[1]
    
    # STRICT pre-relapse rows for nadir / first-increase logic
    df_pre <- dplyr::filter(df, sample_date < prog_date)
    
    if (nrow(df_pre) == 0L) {
      return(tibble::tibble(
        nadir_date = as.Date(NA),  nadir_prob = NA_real_,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    # Nadir = minimum probability before relapse
    nadir_row <- df_pre %>%
      dplyr::slice_min(.data[[assay_prob]], with_ties = FALSE)
    
    nadir_date <- nadir_row$sample_date
    nadir_prob <- nadir_row[[assay_prob]]
    
    # First increase after nadir (strictly later in time AND strictly higher prob)
    post_nadir <- df_pre %>%
      dplyr::filter(sample_date > nadir_date,
                    .data[[assay_prob]] > nadir_prob) %>%
      dplyr::slice_head(n = 1)
    
    if (nrow(post_nadir) == 0L) {
      return(tibble::tibble(
        nadir_date = nadir_date,  nadir_prob = nadir_prob,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase  = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    first_inc_date <- post_nadir$sample_date
    first_inc_prob <- post_nadir[[assay_prob]]
    
    tibble::tibble(
      nadir_date = nadir_date,
      nadir_prob = nadir_prob,
      first_inc_date = first_inc_date,
      first_inc_prob = first_inc_prob,
      days_to_first_increase  = as.numeric(first_inc_date - nadir_date),
      days_before_progression = as.numeric(prog_date - first_inc_date),
      prog_date = prog_date
    )
  }) %>%
  dplyr::ungroup()


# Nadir timing relative to progression
nadir_timing_stats <- advance_df %>%
  mutate(days_nadir_before = as.numeric(prog_date - nadir_date)) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before,    na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before,    na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before,    na.rm = TRUE),
    .groups = "drop"
  )

nadir_timing_stats <- nadir_timing_stats %>%
  mutate(
    median_nadir_mo = d2m(median_nadir_days),
    iqr_nadir_mo    = d2m(iqr_nadir_days),
    min_nadir_mo    = d2m(min_nadir_days),
    max_nadir_mo    = d2m(max_nadir_days)
  )

# First increase timing (after nadir, and its lead time to progression)
increase_stats2 <- advance_df %>%
  summarise(
    median_after_nadir = median(days_to_first_increase,   na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,      na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,      na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,      na.rm = TRUE),
    median_before_prog = median(days_before_progression,  na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,     na.rm = TRUE),
    min_before_prog    = min(days_before_progression,     na.rm = TRUE),
    max_before_prog    = max(days_before_progression,     na.rm = TRUE),
    .groups            = "drop"
  )

increase_stats2 <- increase_stats2 %>%
  mutate(
    median_after_nadir_mo = d2m(median_after_nadir),
    iqr_after_nadir_mo    = d2m(iqr_after_nadir),
    min_after_nadir_mo    = d2m(min_after_nadir),
    max_after_nadir_mo    = d2m(max_after_nadir),
    median_before_prog_mo = d2m(median_before_prog),
    iqr_before_prog_mo    = d2m(iqr_before_prog),
    min_before_prog_mo    = d2m(min_before_prog),
    max_before_prog_mo    = d2m(max_before_prog)
  )

# ==== C) One-paragraph sentence  ====
cat(glue(
  "We analyzed {n_patients} patients (total {n_samples} samples) collected a median ",
  "{sample_timing_stats$median_days_before} days (~{sample_timing_stats$median_days_before_mo} mo) before progression ",
  "(IQR {sample_timing_stats$iqr_days_before} days, ~{sample_timing_stats$iqr_days_before_mo} mo; ",
  "range {sample_timing_stats$min_days_before}–{sample_timing_stats$max_days_before} days, ",
  "~{sample_timing_stats$min_days_before_mo}–{sample_timing_stats$max_days_before_mo} mo). ",
  "The nadir detection probability occurred a median ",
  "{nadir_timing_stats$median_nadir_days} days (~{nadir_timing_stats$median_nadir_mo} mo) before progression ",
  "(IQR {nadir_timing_stats$iqr_nadir_days} days, ~{nadir_timing_stats$iqr_nadir_mo} mo; ",
  "range {nadir_timing_stats$min_nadir_days}–{nadir_timing_stats$max_nadir_days} days, ",
  "~{nadir_timing_stats$min_nadir_mo}–{nadir_timing_stats$max_nadir_mo} mo). ",
  "From that nadir, the first increase occurred a median ",
  "{increase_stats2$median_after_nadir} days (~{increase_stats2$median_after_nadir_mo} mo) later ",
  "(IQR {increase_stats2$iqr_after_nadir} days, ~{increase_stats2$iqr_after_nadir_mo} mo; ",
  "range {increase_stats2$min_after_nadir}–{increase_stats2$max_after_nadir} days, ",
  "~{increase_stats2$min_after_nadir_mo}–{increase_stats2$max_after_nadir_mo} mo), ",
  "which was a median {increase_stats2$median_before_prog} days ",
  "(~{increase_stats2$median_before_prog_mo} mo) before IMWG-defined clinical progression ",
  "(IQR {increase_stats2$iqr_before_prog} days, ~{increase_stats2$iqr_before_prog_mo} mo; ",
  "range {increase_stats2$min_before_prog}–{increase_stats2$max_before_prog} days, ",
  "~{increase_stats2$min_before_prog_mo}–{increase_stats2$max_before_prog_mo} mo).\n"
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

