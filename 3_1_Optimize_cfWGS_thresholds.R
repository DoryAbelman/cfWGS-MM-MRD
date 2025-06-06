# =============================================================================
# Script:   optimize_cfWGS_thresholds.R
# Project:  cfWGS MRD detection (M4 / SPORE / IMMAGINE)
# Author:   Dory Abelman
# Date:     May 28, 2025
#
# Purpose:
#   1. Restrict analysis to frontline induction–transplant cohort.
#   2. Derive a ground‐truth MRD label from clinical assays (MFC or clonoSEQ).
#   3. Optimize thresholds for selected cfWGS metrics (z-score, detect_rate, etc.)
#      via ROC analysis against MRD_truth.
#   4. Evaluate performance (sensitivity, specificity, accuracy, precision, NPV, AUC).
#   5. Fit ridge-penalized logistic combination predictors and compare via
#      accuracy and sens≥0.9 criteria.
#   6. Save optimized thresholds, performance tables, and updated cfWGS calls.
#
# Outputs:
#   • combo_results.csv / combo_results.rds
#   • all_patients_with_BM_and_blood_calls.csv / .rds
#   • STable_* performance tables (*.csv)
#
# =============================================================================

# -------- 0.  Load packages --------------------------------------------------
library(dplyr)
library(tidyr)
library(pROC)       # ROC curves, coords(), auc()
library(yardstick)  # confusion-matrix metrics (if still used; otherwise drop)
library(purrr)      # map_dfr()
library(glmnet)     # ridge-penalized CV
library(Matrix)     # sparse model matrices (dgCMatrix)
library(janitor)    # tabyl() + adorn_totals()
library(glue)       # string glue for messages




### Set paths 
outdir   <- "Output_tables_2025"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

### Load data 
file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.rds")
cohort_df <- readRDS("cohort_assignment_table_updated.rds")

dat <- file 



# ────────────────────────────────────────────────────────────────────────────
# 1.  Join cohort_df and keep frontline only -------------------------------------
dat <- dat %>%                # <‑‑ your master data
  left_join(cohort_df, by = "Patient") 


# ────────────────────────────────────────────────────────────────────────────
# 2.  Create MRD_truth reference label --------------------------------------

# Choose which external test(s) you trust most:
# "Flow_Binary"    -> MFC
# "Adaptive_Binary"-> clonoSEQ
# Go with adaptive first, and then if NA take clonoSEQ
#truth_choice <- "either"

dat <- dat %>% 
  mutate(
    MRD_truth = case_when(
      !is.na(Adaptive_Binary)            ~ Adaptive_Binary,             # use clonoSEQ if available
      is.na(Adaptive_Binary)
      & !is.na(Flow_Binary)              ~ Flow_Binary,                 # else use MFC 
      TRUE                               ~ NA_real_                     # missing if neither assay run
    )
  )

## Keep only MRD timepoints to optimize timepoint 
dat_mrd <- dat %>%
  filter(!timepoint_info %in% c("Diagnosis", "Baseline"))


# create train and test df - train on frontline, apply on hold out df
# frontline = training set, everything else = hold‑out
train_df <- dat_mrd %>% filter(Cohort == "Frontline")
hold_df  <- dat_mrd %>% filter(Cohort != "Frontline")

# ────────────────────────────────────────────────────────────────────────────
# 2.  Helper: build sparse X & y ---------------------------------------------
make_sparse_xy <- function(df, predictors) {
  mm <- model.matrix(
    reformulate(predictors, response = "MRD_truth"),
    data = df
  )
  list(
    x = as(mm[, -1], "dgCMatrix"),  # drop intercept
    y = df$MRD_truth
  )
}

# 3.  Train ridge‐penalized combos ------------------------------------------
# 3a) Classic BM combo
bm_complete    <- train_df %>% drop_na(MRD_truth, zscore_BM, detect_rate_BM)
bm_xy          <- make_sparse_xy(bm_complete, c("zscore_BM", "detect_rate_BM"))
cv_bm          <- cv.glmnet(bm_xy$x, bm_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 3b) Extended BM combo (+ FS, coverage, BM_tf)
bm_ext_complete <- train_df %>% 
  drop_na(MRD_truth, zscore_BM, detect_rate_BM, FS, Mean.Coverage, WGS_Tumor_Fraction_Blood_plasma_cfDNA)
bm_ext_xy       <- make_sparse_xy(
  bm_ext_complete,
  c("zscore_BM", "detect_rate_BM", "FS", "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA")
)
cv_bm_ext      <- cv.glmnet(bm_ext_xy$x, bm_ext_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 3c) Classic blood combo
bl_complete    <- train_df %>% drop_na(MRD_truth, zscore_blood, detect_rate_blood)
bl_xy          <- make_sparse_xy(bl_complete, c("zscore_blood", "detect_rate_blood"))
cv_bl          <- cv.glmnet(bl_xy$x, bl_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 3d) Extended blood combo (+ FS, coverage, blood_tf)
bl_ext_complete <- train_df %>%
  drop_na(MRD_truth, zscore_blood, detect_rate_blood, FS, Mean.Coverage, WGS_Tumor_Fraction_Blood_plasma_cfDNA)
bl_ext_xy       <- make_sparse_xy(
  bl_ext_complete,
  c("zscore_blood", "detect_rate_blood", "FS",
    "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA")
)
cv_bl_ext      <- cv.glmnet(bl_ext_xy$x, bl_ext_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 4.  Add combo probabilities to all patients -------------------------------
add_prob <- function(df, predictors, fit) {
  # logical mask: rows with ALL predictors present
  keep <- stats::complete.cases(df[, predictors])
  
  # initialise output vector
  prob <- rep(NA_real_, nrow(df))
  
  if (any(keep)) {
    mm   <- model.matrix(reformulate(predictors), data = df[keep, ])
    newx <- as(mm[, -1], "dgCMatrix")
    prob[keep] <- as.numeric(
      predict(fit, newx = newx, s = "lambda.min", type = "response")
    )
  }
  prob
}


## Add it to the dataframe
dat <- dat %>%
  mutate(
    combo_BM_prob        = add_prob(., c("zscore_BM", "detect_rate_BM"),          cv_bm),
    combo_BM_ext_prob    = add_prob(., c("zscore_BM", "detect_rate_BM", "FS",
                                         "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"), cv_bm_ext),
    combo_blood_prob     = add_prob(., c("zscore_blood", "detect_rate_blood"),    cv_bl),
    combo_blood_ext_prob = add_prob(., c("zscore_blood", "detect_rate_blood", "FS",
                                         "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"), cv_bl_ext)
  )

# 5.  Define metric sets ----------------------------------------------------
base_metrics     <- c("zscore_BM","detect_rate_BM",
                      "zscore_blood","detect_rate_blood",
                      "combo_BM_prob","combo_blood_prob")
extended_metrics <- c(base_metrics,
                      "FS","Mean.Coverage",
                      "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
                      "combo_BM_ext_prob","combo_blood_ext_prob")

# 6.  Optimiser function ----------------------------------------------------
optimise_metric <- function(df, metric, truth_col="MRD_truth",
                            trial = c("accuracy","sens>=0.9")) {
  trial <- match.arg(trial)
  tmp   <- df %>% select(all_of(c(metric, truth_col))) %>% drop_na()
  roc   <- pROC::roc(tmp[[truth_col]], tmp[[metric]])
  best  <- if (trial=="accuracy") {
    pROC::coords(roc, "best", best.method="youden",
                 ret=c("threshold","sensitivity","specificity","accuracy"),
                 transpose=FALSE)
  } else {
    pROC::coords(roc, "all",
                 ret=c("threshold","sensitivity","specificity","accuracy"),
                 transpose=FALSE) %>%
      filter(sensitivity>=0.9) %>%
      slice_max(specificity, n=1)
  }
  tibble(metric   = metric,
         trial    = trial,
         threshold= best[["threshold"]],
         sens     = best[["sensitivity"]],
         spec     = best[["specificity"]],
         accuracy = best[["accuracy"]],
         auc      = as.numeric(pROC::auc(roc)))
}

# -------------- Re-slice out train & hold‑out ----------------------------
train_df <- dat %>%
  filter(Cohort == "Frontline",
         !timepoint_info %in% c("Diagnosis", "Baseline"))   # keep MRD visits only

hold_df  <- dat %>%
  filter(Cohort != "Frontline",
         !timepoint_info %in% c("Diagnosis", "Baseline"))

run_trial <- function(metrics, trial_tag, df_train) {
  map_dfr(metrics, optimise_metric,
          df        = df_train,
          truth_col = "MRD_truth",
          trial     = trial_tag)
}

# 7.  Run trials on FRONTLINE training set ---------------------------------
results <- bind_rows(
  run_trial(base_metrics,     "accuracy",  train_df) %>% mutate(set="base"),
  run_trial(base_metrics,     "sens>=0.9", train_df) %>% mutate(set="base"),
  run_trial(extended_metrics, "accuracy",  train_df) %>% mutate(set="extended"),
  run_trial(extended_metrics, "sens>=0.9", train_df) %>% mutate(set="extended")
)
print(results)




#### See if different combination of features would be better
# 1.  Define candidate combos ------------------------------------------------
combos <- list(
  # — BM combos with zscore —
  BM_zscore_only         = c("zscore_BM"),
  BM_base             = c("zscore_BM",    "detect_rate_BM"),
  BM_plus_FS          = c("zscore_BM",    "detect_rate_BM", "FS"),
  BM_plus_MeanCov     = c("zscore_BM",    "detect_rate_BM", "Mean.Coverage"),
  BM_plus_TF          = c("zscore_BM",    "detect_rate_BM", "WGS_Tumor_Fraction_BM_cells"),
  BM_all_extras       = c("zscore_BM",    "detect_rate_BM", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_BM_cells"),
  # — BM combos without zscore —
  BM_rate_base        = c("detect_rate_BM"),
  BM_rate_plus_FS     = c("detect_rate_BM", "FS"),
  BM_rate_plus_MeanCov= c("detect_rate_BM", "Mean.Coverage"),
  BM_rate_plus_TF     = c("detect_rate_BM", "WGS_Tumor_Fraction_BM_cells"),
  BM_rate_all_extras  = c("detect_rate_BM", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_BM_cells"),

  # — Blood combos with zscore —
  Blood_zscore_only      = c("zscore_blood"),
  Blood_base          = c("zscore_blood", "detect_rate_blood"),
  Blood_plus_FS       = c("zscore_blood", "detect_rate_blood", "FS"),
  Blood_plus_MeanCov  = c("zscore_blood", "detect_rate_blood", "Mean.Coverage"),
  Blood_plus_TF       = c("zscore_blood", "detect_rate_blood",
                          "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_all_extras    = c("zscore_blood", "detect_rate_blood", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  # — Blood combos without zscore —
  Blood_rate_base     = c("detect_rate_blood"),
  Blood_rate_plus_FS  = c("detect_rate_blood", "FS"),
  Blood_rate_plus_MeanCov = c("detect_rate_blood", "Mean.Coverage"),
  Blood_rate_plus_TF  = c("detect_rate_blood",
                          "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_rate_all_extras   = c("detect_rate_blood", "FS",
                              "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  
  ### BM and blood together 
  BM_blood_ultimate = c("zscore_BM", "zscore_blood", "detect_rate_blood", "FS",
                        "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA")
)


# 2.  Utility: evaluate one combo ------------------------------------------
eval_combo <- function(preds, label) {
  # complete‑case slice
  df_cc <- train_df %>% drop_na(MRD_truth, all_of(preds))
  
  # skip if too small or only one class
  if (nrow(df_cc) < 20 || n_distinct(df_cc$MRD_truth) < 2)
    return(tibble())
  
  # --------------------------------------------------------------------
  # ➊  FIT MODEL  (glmnet if ≥2 preds, else glm) ------------------------
  # --------------------------------------------------------------------
  if (length(preds) >= 2) {
    # ridge CV
    mm  <- model.matrix(reformulate(preds, response = "MRD_truth"), df_cc)
    X   <- as(mm[, -1, drop = FALSE], "dgCMatrix")
    y   <- df_cc$MRD_truth
    
    fit <- cv.glmnet(X, y,
                     family  = "binomial",
                     alpha   = 0,
                     nfolds  = 5,
                     nlambda = 30)
    
    prob <- as.numeric(
      predict(fit, newx = X, s = "lambda.min", type = "response")
    )
  } else {
    # single predictor → plain logistic regression
    form <- reformulate(preds, response = "MRD_truth")
    fit  <- glm(form, data = df_cc, family = binomial)
    prob <- predict(fit, type = "response")
  }
  
  # --------------------------------------------------------------------
  # ➋  ROC & performance metrics ---------------------------------------
  # --------------------------------------------------------------------
  roc_obj <- pROC::roc(df_cc$MRD_truth, prob)
  
  # Accuracy‑optimal threshold
  acc_best <- pROC::coords(
    roc_obj, "best", best.method = "youden",
    ret = c("threshold","sensitivity","specificity","accuracy"),
    transpose = FALSE
  )
  # Highest specificity with sens ≥ 0.9
  sens90 <- pROC::coords(
    roc_obj, "all",
    ret = c("threshold","sensitivity","specificity","accuracy"),
    transpose = FALSE
  ) %>%
    filter(sensitivity >= 0.9) %>%
    slice_max(specificity, n = 1)
  
  bind_rows(
    tibble(combo = label, trial = "accuracy",  !!!acc_best),
    tibble(combo = label, trial = "sens>=0.9", !!!sens90)
  ) %>%
    mutate(auc = as.numeric(pROC::auc(roc_obj)))
}


# 3.  Run across all combos --------------------------------------------------
combo_results <- imap_dfr(combos, eval_combo)

print(combo_results)

# Save as an RDS:
saveRDS(combo_results,
        file = file.path(outdir, "combo_results.rds"))

# Export as CSV (no row names):
write.csv(combo_results,
          file = file.path(outdir, "combo_results.csv"),
          row.names = FALSE)







### Now select best one and apply to everything 
# ────────────────────────────────────────────────────────────────────────────
# 4.  SELECT THE THREE RULES YOU WANT ----------------------------------------
# helper to pull top‑n accuracy rows for a given combo
pick_rows <- function(cmb, n = 1) {
  combo_results %>%
    filter(combo == cmb) %>%
    slice_max(accuracy, n = n, with_ties = FALSE)
}

bm_row        <- pick_rows("BM_zscore_only", n = 1)
blood_rows    <- pick_rows("Blood_all_extras", n = 2)
blood_row_add <- pick_rows("Blood_base", n = 1)
selected_rows <- bind_rows(bm_row, blood_rows, blood_row_add)   # 4 rows in total
print(selected_rows)

# ────────────────────────────────────────────────────────────────────────────
# 5.  Function to refit & apply a single rule --------------------------------
# ────────────────────────────────────────────────────────────────────────────
fit_and_apply <- function(row, ix, dat, train_df, combos) {
  cmb   <- row$combo
  thr   <- row$threshold
  preds <- combos[[cmb]]
  
  # ---- (re)fit -------------------------------------------------------------
  train_cc <- train_df %>% drop_na(MRD_truth, all_of(preds))
  
  if (length(preds) >= 2) {
    mm <- model.matrix(reformulate(preds, response = "MRD_truth"), train_cc)
    X  <- as(mm[, -1, drop = FALSE], "dgCMatrix")
    y  <- train_cc$MRD_truth
    
    fit <- cv.glmnet(X, y, family = "binomial",
                     alpha = 0, nfolds = 5, nlambda = 30)
    
    predict_fun <- function(df_sub) {
      mm  <- model.matrix(reformulate(preds), data = df_sub)
      newx <- as(mm[, -1, drop = FALSE], "dgCMatrix")
      as.numeric(predict(fit, newx = newx, s = "lambda.min",
                         type = "response"))
    }
  } else {
    fit <- glm(reformulate(preds, response = "MRD_truth"),
               data = train_cc, family = binomial)
    
    predict_fun <- function(df_sub) {
      predict(fit, newdata = df_sub, type = "response")
    }
  }
  
  # ---- apply to all rows ---------------------------------------------------
  keep <- complete.cases(dat[, preds])
  prob <- rep(NA_real_, nrow(dat))
  if (any(keep)) prob[keep] <- predict_fun(dat[keep, ])
  
  # suffix: only for the two Blood_all_extras rows
  suffix <- if (cmb == "Blood_all_extras") paste0("_acc", ix) else ""
  
  prob_col  <- paste0(cmb, suffix, "_prob")
  call_col  <- paste0(cmb, suffix, "_call")
  
  dat[[prob_col]]  <- prob
  dat[[call_col]]  <- if_else(prob >= thr, 1L, 0L, NA_integer_)
  dat
}

# ────────────────────────────────────────────────────────────────────────────
# 6.  LOOP OVER THE 3 SELECTED ROWS -----------------------------------------
for (i in seq_len(nrow(selected_rows))) {
  dat <- fit_and_apply(selected_rows[i, ], i, dat, train_df, combos)
}



# ────────────────────────────────────────────────────────────────────────────
# 7.  Quick checks -----------------------------------------------------------
message("BM_zscore_only positives: ",
        sum(dat$BM_zscore_only_call == 1, na.rm = TRUE))
message("Blood_all_extras_acc2 positives: ",
        sum(dat$Blood_all_extras_acc2_call == 1, na.rm = TRUE))
message("Blood_all_extras_acc3 positives: ",
        sum(dat$Blood_all_extras_acc3_call == 1, na.rm = TRUE))


### Rename for clarity 
dat <- dat %>%
  rename(
    BloodSensPriority_prob = Blood_all_extras_acc2_prob,
    BloodSensPriority_call = Blood_all_extras_acc2_call,
    BloodAccPriority_prob   = Blood_all_extras_acc3_prob,
    BloodAccPriority_call   = Blood_all_extras_acc3_call
  )

message("BloodSensPriority positives: ",
        sum(dat$BloodSensPriority_call == 1, na.rm = TRUE))
message("BloodAccPriority positives: ",
        sum(dat$BloodAccPriority_call == 1, na.rm = TRUE))
message("BloodBase positives: ",
        sum(dat$Blood_base_call == 1, na.rm = TRUE))



# ────────────────────────────────────────────────────────────────────────────
# 8.  Export -----------------------------------------------------------------
write.csv(dat,
          file = file.path(outdir,
                           "all_patients_with_BM_and_blood_calls.csv"),
          row.names = FALSE)
saveRDS(dat,
        file = file.path(outdir,
                         "all_patients_with_BM_and_blood_calls.rds"))





### Now get some proportions for the manuscript
# ────────────────────────────────────────────────────────────────────────────

## First get the number of patients and samples available as well as their metrics 

# define  call‐columns and a helper to summarise one
rules <- c("BM_zscore_only_call", 
           "BloodSensPriority_call", 
           "Blood_base_call")

summarise_rule <- function(df, call_col) {
  df %>%
    # keep only rows where we have both a truth and a call
    filter(!is.na(MRD_truth), !is.na(.data[[call_col]])) %>%
    filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>% # remove baseline
    group_by(Cohort) %>%
    summarise(
      samples      = n(),                             # total samples
      patients     = n_distinct(Patient),             # unique patients
      positives    = sum(.data[[call_col]] == 1),     # call==1
      negatives    = sum(.data[[call_col]] == 0),     # call==0
      sens         = round(
        sum(.data[[call_col]] == 1 & MRD_truth == 1) /
          sum(MRD_truth == 1),
        3),
      spec         = round(
        sum(.data[[call_col]] == 0 & MRD_truth == 0) /
          sum(MRD_truth == 0),
        3),
      accuracy     = round((positives + negatives) / samples, 3)
    ) %>%
    mutate(rule = call_col) %>%
    select(rule, everything())
}

# 1) Overall by cohort
overall_tbl <- bind_rows(
  lapply(rules, function(cl) summarise_rule(dat, cl))
)
print(overall_tbl)

# 2) By cohort × timepoint_info
by_tp_tbl <- bind_rows(
  lapply(rules, function(call_col) {
    dat %>%
      filter(!is.na(MRD_truth), !is.na(.data[[call_col]])) %>%
      group_by(Cohort, timepoint_info) %>%
      summarise(
        samples   = n(),
        patients  = n_distinct(Patient),
        positives = sum(.data[[call_col]] == 1),
        negatives = sum(.data[[call_col]] == 0),
        sens      = round(sum(.data[[call_col]] == 1 & MRD_truth == 1) / sum(MRD_truth == 1), 3),
        spec      = round(sum(.data[[call_col]] == 0 & MRD_truth == 0) / sum(MRD_truth == 0), 3),
        accuracy  = round((positives + negatives) / samples, 3)
      ) %>%
      mutate(rule = call_col) %>%
      select(rule, everything())
  })
)
print(by_tp_tbl)

# 3) Export for manuscript
write.csv(overall_tbl, "MRD_performance_by_cohort.csv", row.names = FALSE)
write.csv(by_tp_tbl,  "MRD_performance_by_cohort_and_timepoint.csv", row.names = FALSE)



### See how many additional patients had blood samples but not BM 
dat %>%
  filter(!is.na(zscore_blood), is.na(zscore_BM)) %>%
  filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  distinct(Patient, Cohort) %>%
  count(Cohort, name = "n_patients_with_blood_only")

dat %>%
  filter(!is.na(zscore_blood), is.na(zscore_BM)) %>%
  filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  filter(!is.na(MRD_truth)) %>%
  distinct(Patient, Cohort) %>%
  count(Cohort, name = "n_patients_with_blood_only")


# ────────────────────────────────────────────────────────────────────────────
# A. 1) OVERALL performance by Cohort --------------------------------------
# 1) Define your three “rules” (call columns)
rules  <- c("BM_zscore_only_call",
            "BloodSensPriority_call",
            "BloodAccPriority_call",
            "Blood_base_call")



cohorts <- unique(dat$Cohort)
tps     <- unique(dat$timepoint_info)

# — 1. overall by cohort ----------------------------------------------------
overall_list <- list()
ctr <- 1
for(cl in rules){
  for(co in cohorts){
    dfsub <- dat %>%
      filter(Cohort == co,
             !is.na(MRD_truth),
             !is.na(.data[[cl]]))
    n   <- nrow(dfsub)
    tp  <- sum(dfsub[[cl]] == 1 & dfsub$MRD_truth == 1)
    tn  <- sum(dfsub[[cl]] == 0 & dfsub$MRD_truth == 0)
    fp  <- sum(dfsub[[cl]] == 1 & dfsub$MRD_truth == 0)
    fn  <- sum(dfsub[[cl]] == 0 & dfsub$MRD_truth == 1)
    sens <- if((tp+fn)>0) tp/(tp+fn) else NA
    spec <- if((tn+fp)>0) tn/(tn+fp) else NA
    acc  <- if(n>0) (tp+tn)/n else NA
    ppv  <- if((tp+fp)>0) tp/(tp+fp) else NA
    npv  <- if((tn+fn)>0) tn/(tn+fn) else NA
    
    overall_list[[ctr]] <- tibble(
      rule        = cl,
      cohort      = co,
      n           = n,
      tp, tn, fp, fn,
      sensitivity = round(sens,3),
      specificity = round(spec,3),
      accuracy    = round(acc,3),
      ppv         = round(ppv,3),
      npv         = round(npv,3)
    )
    ctr <- ctr + 1
  }
}
tbl_overall <- bind_rows(overall_list)

write.csv(
  tbl_overall,
  file = file.path(outdir, "STable_overall_perf_by_cohort.csv"),
  row.names = FALSE
)

# — 2. performance by cohort × timepoint_info ------------------------------
by_tp_list <- list()
ctr <- 1
for(cl in rules){
  for(co in cohorts){
    for(tp in tps){
      dfsub <- dat %>%
        filter(Cohort == co,
               timepoint_info == tp,
               !is.na(MRD_truth),
               !is.na(.data[[cl]]))
      n   <- nrow(dfsub)
      tp_  <- sum(dfsub[[cl]] == 1 & dfsub$MRD_truth == 1)
      tn_  <- sum(dfsub[[cl]] == 0 & dfsub$MRD_truth == 0)
      fp_  <- sum(dfsub[[cl]] == 1 & dfsub$MRD_truth == 0)
      fn_  <- sum(dfsub[[cl]] == 0 & dfsub$MRD_truth == 1)
      sens <- if((tp_+fn_)>0) tp_/(tp_+fn_) else NA
      spec <- if((tn_+fp_)>0) tn_/(tn_+fp_) else NA
      acc  <- if(n>0) (tp_+tn_)/n else NA
      
      by_tp_list[[ctr]] <- tibble(
        rule           = cl,
        cohort         = co,
        timepoint_info = tp,
        n              = n,
        tp   = tp_, tn = tn_, fp = fp_, fn = fn_,
        sensitivity = round(sens,3),
        specificity = round(spec,3),
        accuracy    = round(acc,3)
      )
      ctr <- ctr + 1
    }
  }
}
tbl_by_timepoint <- bind_rows(by_tp_list)

tbl_by_timepoint <- tbl_by_timepoint %>% filter(!is.na(accuracy))
write.csv(
  tbl_by_timepoint,
  file = file.path(outdir, "STable_perf_by_cohort_and_timepoint.csv"),
  row.names = FALSE
)

# — 3. positivity / negativity rates by cohort × timepoint_info ------------
rates_list <- list()
ctr <- 1
for(cl in rules){
  for(co in cohorts){
    for(tp in tps){
      dfsub <- dat %>%
        filter(Cohort == co,
               timepoint_info == tp,
               !is.na(.data[[cl]]))
      n_obs    <- nrow(dfsub)
      n_pos    <- sum(dfsub[[cl]] == 1)
      n_neg    <- sum(dfsub[[cl]] == 0)
      pos_rate <- if(n_obs>0) n_pos / n_obs else NA
      neg_rate <- if(n_obs>0) n_neg / n_obs else NA
      
      rates_list[[ctr]] <- tibble(
        rule           = cl,
        cohort         = co,
        timepoint_info = tp,
        n_obs, n_pos, n_neg,
        pos_rate = round(pos_rate,3),
        neg_rate = round(neg_rate,3)
      )
      ctr <- ctr + 1
    }
  }
}
tbl_rates_tp <- bind_rows(rates_list)

write.csv(
  tbl_rates_tp,
  file = file.path(outdir, "STable_pos_neg_rates_by_cohort_and_timepoint.csv"),
  row.names = FALSE
)


# ────────────────────────────────────────────────────────────────────────────
# C.  2) POS / NEG counts by time‑point & cohort ----------------------------
tbl_by_tp <- dat %>%
  mutate(Call = factor(BloodSensPriority_call, levels = c(0,1),
                       labels = c("Negative","Positive")),
         Cohort_group = ifelse(Cohort == "Frontline","Frontline","Non‑frontline")) %>%
  group_by(Cohort_group, timepoint_info, Call) %>%
  summarise(N = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Call, values_from = N, values_fill = 0) %>%
  mutate(total = Positive + Negative,
         pos_rate = round(Positive / total, 3))

write.csv(tbl_by_tp,
          file = file.path(outdir,
                           "STable_counts_by_timepoint.csv"),
          row.names = FALSE)

# ────────────────────────────────────────────────────────────────────────────
# D.  3) Landmark MRD‑negativity rates --------------------------------------
landmark_tp <- c("Post‑ASCT","Maintenance‑1yr")

tbl_landmark <- dat %>%
  filter(timepoint_info %in% landmark_tp) %>%
  mutate(cfWGS = BloodSensPriority_call,
         clonoSEQ = Adaptive_Binary,
         MFC      = Flow_Binary) %>%
  tidyr::pivot_longer(cfWGS:MFC, names_to = "Assay", values_to = "Result") %>%
  mutate(Neg = Result == 0) %>%
  group_by(Cohort, timepoint_info, Assay) %>%
  summarise(N = n(), Neg = sum(Neg, na.rm = TRUE),
            Neg_rate = round(Neg / N, 3),
            .groups = "drop")

write.csv(tbl_landmark,
          file = file.path(outdir,
                           "STable_landmark_negativity.csv"),
          row.names = FALSE)

# ────────────────────────────────────────────────────────────────────────────
# E.  4) save the threshold table you selected ------------------------------
write.csv(selected_rows,
          file = file.path(outdir,
                           "STable_selected_thresholds.csv"),
          row.names = FALSE)

# ────────────────────────────────────────────────────────────────────────────
# F.  5) simple cohort totals for PRISMA / flow chart -----------------------
tbl_cohort_totals <- dat %>%
  group_by(Cohort) %>%
  summarise(
    N_rows            = n(),
    N_BM              = sum(!is.na(BM_zscore_only_call)),
    N_Blood           = sum(!is.na(BloodSensPriority_call)),
    BM_pos            = sum(BM_zscore_only_call == 1, na.rm = TRUE),
    Blood_pos         = sum(BloodSensPriority_call == 1, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(tbl_cohort_totals,
          file = file.path(outdir,
                           "STable_cohort_flowchart_counts.csv"),
          row.names = FALSE)

message("Exported 5 summary tables to ", outdir)
