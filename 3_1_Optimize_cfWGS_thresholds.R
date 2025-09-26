# =============================================================================
# Script:   optimize_cfWGS_thresholds.R
# Project:  cfWGS MRD detection (M4 / SPORE / IMMAGINE)
# Author:   Dory Abelman
# Date:     May 28, 2025
#
# Purpose:
#   • Load & preprocess clinical and cfWGS feature data
#   • Define MRD ground truth from clonoSEQ/MFC assays
#   • Optimize univariate thresholds for cfWGS metrics
#   • Train & nested-CV elastic-net classifiers on BM, blood, and fragmentomics
#   • Select top models and thresholds; save them for scoring
#   • Score full cohort and dilution series; export performance tables
#
# Outputs:
#   ├ combo_results_*.{rds,csv}
#   ├ models_cfWGS_final.rds
#   ├ selected_combo_{models,thresholds}.rds
#   ├ all_patients_with_*_calls.{rds,csv}
#   └ STable_*.csv
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup & Package Loading
# -----------------------------------------------------------------------------

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
library(Matrix)
library(caret)
library(viridis) 
library(patchwork)
library(DescTools)   # for CalibrationPlot()
library(PRROC)      # pr.curve()
library(scales)     # percent_format()
library(stringr)
library(forcats)

# -----------------------------------------------------------------------------
# 2. Data Import & Cohort Definition
# -----------------------------------------------------------------------------

### Set paths 
outdir   <- "Output_tables_2025"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

### Load data 
file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated7.rds")
cohort_df <- readRDS("cohort_assignment_table_updated.rds")

dat <- file 

# 1.  Join cohort_df and keep frontline only -------------------------------------
dat <- dat %>%                # <‑‑ your master data
  left_join(cohort_df, by = "Patient") 


# -----------------------------------------------------------------------------
# 3. Ground-Truth MRD Label Construction
# -----------------------------------------------------------------------------
#    Create MRD_truth reference label --------------------------------------

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
hold_df  <- dat_mrd %>% filter(Cohort == "Non-frontline")


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

# -----------------------------------------------------------------------------
# 4. Univariate Threshold Optimization - not used in final manuscript
# -----------------------------------------------------------------------------

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
  filter(Cohort == "Non-frontline",
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



# -----------------------------------------------------------------------------
# 5. Single 5-fold CV - not used in manuscript but for testing only 
# -----------------------------------------------------------------------------

#### See if different combination of features would be better
# 1.  Define candidate combos ------------------------------------------------
combos <- list(
  # — BM combos with zscore —
  BM_zscore_only         = c("zscore_BM"),
  BM_base             = c("zscore_BM",    "detect_rate_BM"),
  BM_plus_FS          = c("zscore_BM",    "detect_rate_BM", "FS", "Proportion.Short"),
  BM_plus_MeanCov     = c("zscore_BM",    "detect_rate_BM", "Mean.Coverage"),
  BM_plus_TF          = c("zscore_BM",    "detect_rate_BM", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  BM_all_extras       = c("zscore_BM",    "detect_rate_BM", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  # — BM combos without zscore —
  BM_rate_base        = c("detect_rate_BM"),
  BM_rate_plus_FS     = c("detect_rate_BM", "FS", "Proportion.Short"),
  BM_rate_plus_MeanCov= c("detect_rate_BM", "Mean.Coverage"),
  BM_rate_plus_TF     = c("detect_rate_BM", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  BM_rate_all_extras  = c("detect_rate_BM", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),

  # — Blood combos with zscore —
  Blood_zscore_only      = c("zscore_blood"),
  Blood_base          = c("zscore_blood", "detect_rate_blood"),
  Blood_plus_FS       = c("zscore_blood", "detect_rate_blood", "FS", "Proportion.Short"),
  Blood_plus_MeanCov  = c("zscore_blood", "detect_rate_blood", "Mean.Coverage"),
  Blood_plus_TF       = c("zscore_blood", "detect_rate_blood",
                          "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_all_extras    = c("zscore_blood", "detect_rate_blood", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  # — Blood combos without zscore —
  Blood_rate_base     = c("detect_rate_blood"),
  Blood_rate_plus_FS  = c("detect_rate_blood", "FS"),
  Blood_rate_plus_MeanCov = c("detect_rate_blood", "Mean.Coverage"),
  Blood_rate_plus_TF  = c("detect_rate_blood",
                          "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_rate_all_extras   = c("detect_rate_blood", "FS",
                              "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  
  ### BM and blood together 
  BM_blood_ultimate = c("zscore_BM", "zscore_blood", "detect_rate_blood", "FS",
                        "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA"), 
  
  ### Just with fragmentomics 
  Just_fragmentomics_full = c("FS", "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  Just_fragmentomics_small = c("FS", "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Just_fragmentomics_small2 = c("FS", "Mean.Coverage")
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
set.seed(2025)
combo_results <- imap_dfr(combos, eval_combo)

print(combo_results)

# Save as an RDS:
saveRDS(combo_results,
        file = file.path(outdir, "combo_results_updated.rds"))

# Export as CSV (no row names):
write.csv(combo_results,
          file = file.path(outdir, "combo_results_updated.csv"),
          row.names = FALSE)



#  Save fitted combo models for later scoring ─────────────────────────

# Collect the exact fits you just trained
model_list <- list(
  BM_zscore_only   = cv_bm,       # classic BM combo
  BM_ext           = cv_bm_ext,   # extended BM combo
  Blood_base       = cv_bl,       # classic blood combo
  Blood_ext        = cv_bl_ext    # extended blood combo
)

# Write them to disk so your dilution-series script can just load & predict
saveRDS(model_list,
        file = file.path(outdir, "models_cfWGS_final.rds"))





# -----------------------------------------------------------------------------
# 6. Repeated CV - not used in manuscript but for testing only 
# -----------------------------------------------------------------------------

###### Now do a repeated 5x5 cross validation to improve generalizability 
###### Final numbers reported in the manuscript

# parallel backend (optional, speeds up train())
library(doParallel)
registerDoParallel()

# your data splits
train_df <- dat_mrd %>% filter(Cohort == "Frontline") 
hold_df  <- dat_mrd %>% filter(Cohort == "Non-frontline")

train_df <- train_df %>% filter(!is.na(MRD_truth))
hold_df <- hold_df %>% filter(!is.na(MRD_truth))

# convert binary to factor with levels "neg"/"pos"
train_df <- train_df %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))
hold_df  <- hold_df  %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))

positive_class <- "pos"
# ──────────────────────────────────────────────────────────────────────────────
# 1.  caret::train control (5×5 repeated CV) ------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
ctrl <- trainControl(
  method            = "repeatedcv",
  number            = 5,
  repeats           = 5,
  classProbs        = TRUE,
  summaryFunction   = twoClassSummary,
  savePredictions   = "final",   # store OOF preds
  allowParallel     = TRUE
)

# ──────────────────────────────────────────────────────────────────────────────
# 2.  Fit each combo, get CV + hold‐out metrics -------------------------------
# ──────────────────────────────────────────────────────────────────────────────
## Helper function 
youden_metrics <- function(response, prob, thr) {
  roc_obj <- roc(response, prob, quiet = TRUE, levels = c("neg", "pos"))
  tibble(
    auc = as.numeric(auc(roc_obj)),
    sens = coords(roc_obj, thr, ret = "sensitivity"),
    spec = coords(roc_obj, thr, ret = "specificity"),
    acc  = mean((prob >= thr) == (response == positive_class))
  )
}

## Run on cohort 
results <- imap(combos, function(preds, label) {
  
  # complete cases in train
  df_cc <- train_df %>% drop_na(all_of(preds), MRD_truth)
  if(nrow(df_cc) < 20 || n_distinct(df_cc$MRD_truth) < 2) {
    message("Skipping ", label, " (n=", nrow(df_cc), ")")
    return(NULL)
  }
  
  # fit the model (glm or glmnet)
  form <- reformulate(preds, "MRD_truth")
  fit  <- suppressWarnings(
    if (length(preds) == 1) {
      train(form, df_cc, method="glm", family=binomial,
            metric="ROC", trControl=ctrl, preProc=c("center","scale"))
    } else {
      train(form, df_cc, method="glmnet", metric="ROC",
            tuneLength=10, trControl=ctrl, preProc=c("center","scale"))
    }
  )
  
  # pull out OOF preds
  oof <- if (length(preds)==1) fit$pred else {
    fit$pred %>%
      filter(alpha==fit$bestTune$alpha,
             lambda==fit$bestTune$lambda)
  }
  roc_oof <- roc(oof$obs, oof[[positive_class]], quiet=TRUE,
                 levels=c("neg","pos"))
  
  # helper to pick threshold by criterion
  pick_thr <- function(mode){
    df <- coords(
      roc_oof, "all",
      ret = c("threshold","sensitivity","specificity","accuracy"),
      transpose = FALSE
    ) %>% as_tibble()
    if(mode=="youden"){
      coords(
        roc_oof, "best", best.method="youden",
        ret="threshold", transpose=FALSE
      )[[1]]
    } else if(mode=="sens90"){
      df %>%
        filter(sensitivity >= .9) %>%
        slice_max(specificity, n = 1, with_ties = FALSE) %>%
        pull(threshold)
    } else {
      df %>%
        filter(specificity >= .9) %>%
        slice_max(sensitivity, n = 1, with_ties = FALSE) %>%
        pull(threshold)
    }
  }
  
  # now build a numeric vector of exactly three thresholds:
  thr_list <- map_dbl(c("youden","sens90","spec90"), pick_thr)
  
  
  # build a little function to compute metrics at a given thr
  eval_at <- function(roc_obj, response, prob, thr){
    coords(roc_obj, x=thr,
           ret=c("sensitivity","specificity","accuracy"),
           transpose=FALSE) %>%
      as_tibble() %>%
      mutate(
        auc = as.numeric(auc(roc_obj)),
        thr = thr
      )
  }
  
  # 1) CV metrics for each mode
  cv_metrics <- map2_dfr(
    thr_list,
    c("youden","sens>=0.9","spec>=0.9"),
    ~ eval_at(roc_oof, oof$obs, oof[[positive_class]], .x) %>%
      mutate(
        combo = label,    # <— the feature-combo name
        mode  = .y,       # <— which thresholding mode
        .before = 1
      )
  )
  
  # now do hold-out
  hold_cc    <- hold_df %>% drop_na(all_of(preds), MRD_truth)
  probs_hold <- predict(fit, newdata=hold_cc, type="prob")[[positive_class]]
  roc_hold   <- roc(hold_cc$MRD_truth, probs_hold, quiet=TRUE,
                    levels=c("neg","pos"))
  
  hold_metrics <- map2_dfr(
    thr_list,
    c("youden","sens>=0.9","spec>=0.9"),
    ~ eval_at(roc_hold, hold_cc$MRD_truth, probs_hold, .x) %>%
      mutate(
        combo = label,    # <— again, the feature-combo name
        mode  = .y,
        .before = 1
      )
  )
  
  # return
  list(cv   = cv_metrics,
       hold = hold_metrics,
       fit  = fit                   # <<< return the caret model object
  )
})


# bind into tables
nested_metrics <- map_dfr(results, "cv")
hold_metrics   <- map_dfr(results, "hold")

# Save for manuscript ------------------------------------------------------
write_csv(nested_metrics, "nested_metrics.csv")
write_csv(hold_metrics,   "hold_metrics.csv")


### Now select best one and apply to everything 
# ────────────────────────────────────────────────────────────────────────────
# 4.  SELECT THE THREE RULES YOU WANT ----------------------------------------
# helper to pull top‑n accuracy rows for a given combo
pick_rows <- function(cmb, n = 1) {
  nested_metrics %>%
    filter(combo == cmb) %>%
    slice_max(accuracy, n = n, with_ties = FALSE)
}

bm_row        <- pick_rows("BM_zscore_only", n = 1)
blood_rows    <- pick_rows("Blood_rate_all_extras", n = 1)
blood_row_add <- pick_rows("Blood_base", n = 1)
selected_rows <- bind_rows(bm_row, blood_rows, blood_row_add)   # 4 rows in total
print(selected_rows)


# ────────────────────────────────────────────────────────────────────────────
# 4b.  Fit & save only the selected combo models ----------------------------

# Extract the combo names you picked
selected_combos <- selected_rows$combo

# Get the models from before
selected_models <- map(results[selected_combos], "fit")

# Save both the models and the thresholds
saveRDS(selected_models,
        file = file.path(outdir, "selected_combo_models_updated.rds"))
saveRDS(selected_rows,
        file = file.path(outdir, "selected_combo_thresholds_updated.rds"))



# ────────────────────────────────────────────────────────────────────────────
# 5. NOW APPLY TO REST OF COHORT  --------------------------------------------

# Define function 
apply_selected <- function(row, idx, dat, combos, models, positive_class){
  cmb   <- row$combo
  thr   <- row$threshold      # make sure your selected_rows has exactly this column!
  preds <- combos[[cmb]]
  fit   <- models[[cmb]]
  
  # build a prediction function
  if (inherits(fit, "cv.glmnet")) {
    predict_fun <- function(df_sub) {
      mm   <- model.matrix(reformulate(preds), data = df_sub)
      Xnew <- as(mm[,-1, drop=FALSE], "dgCMatrix")
      as.numeric(predict(fit, newx = Xnew, s="lambda.min", type="response"))
    }
  } else if (inherits(fit, "glm")) {
    predict_fun <- function(df_sub) {
      predict(fit, newdata = df_sub, type="response")
    }
  } else if (inherits(fit, "train")) {
    predict_fun <- function(df_sub) {
      prob_df <- predict(fit, newdata = df_sub, type="prob")
      prob_df[[ positive_class ]]
    }
  } else {
    stop("Unknown model type for combo: ", cmb)
  }
  
  # apply to all rows with non-missing preds
  keep <- complete.cases(dat[, preds, drop=FALSE])
  prob <- rep(NA_real_, nrow(dat))
  prob[keep] <- predict_fun(dat[keep, ])
  
  # optional suffix if you have multiples of same combo
  suffix <- if(cmb=="Blood_all_extras") paste0("_acc", idx) else ""
  pcol   <- paste0(cmb, suffix, "_prob")
  ccol   <- paste0(cmb, suffix, "_call")
  
  dat[[pcol]] <- prob
  dat[[ccol]] <- if_else(prob >= thr, 1L, 0L, NA_integer_)
  
  dat
}

selected_rows <- selected_rows %>%
  mutate(threshold = thr)    # or rename(threshold = threshold_cv)

# Loop through all selected rows
for(i in seq_len(nrow(selected_rows))) {
  dat <- apply_selected(
    selected_rows[i, ], 
    i, 
    dat, 
    combos, 
    selected_models, 
    positive_class
  )
}





### This is main model used in manuscript
# -----------------------------------------------------------------------------
# 7. Elastic-Net Classifier Training & Nested CV - main model
# -----------------------------------------------------------------------------

### Next do a full nested cross validation 5x5 with hyperparameter tuning 
### Even more robust with less data leakage
positive_class <- "pos"

### Since doing hyperparameter tuning, can use less features 
combos_small <- list(
  # Bone marrow signal comparisons
  BM_zscore_only_sites     = c("zscore_BM"),
  BM_zscore_only_detection_rate     = c("z_score_detection_rate_BM"),
  BM_rate_only = c("detect_rate_BM"),
  BM_base            = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM"),
  BM_base_zscore            = c("zscore_BM", "z_score_detection_rate_BM"),
  BM_plus_fragment   = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM", "FS", "Mean.Coverage", "Proportion.Short", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  BM_plus_fragment_min   = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM", "FS", "Mean.Coverage"),
  
  # Blood signal comparisons
  Blood_zscore_only_sites  = c("zscore_blood"),
  Blood_zscore_only_detection_rate  = c("z_score_detection_rate_blood"),
  Blood_rate_only = c("detect_rate_blood"),
  Blood_base         = c("zscore_blood", "detect_rate_blood", "z_score_detection_rate_blood"),
  Blood_base_zscore         = c("zscore_blood", "z_score_detection_rate_blood"),
  Blood_plus_fragment = c("zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood", "FS", "Mean.Coverage", "Proportion.Short", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_plus_fragment_min = c("zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood", "FS", "Mean.Coverage"),
  
  # Combined signals
  BM_blood_combo     = c("zscore_BM", "zscore_blood", "detect_rate_blood", "FS", "Mean.Coverage"),
  
  # Fragmentomics-only
  Fragmentomics_full = c("FS", "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  Fragmentomics_min  = c("FS", "Mean.Coverage"),
  Fragmentomics_FS_only = c("FS"), 
  Fragmentomics_mean_coverage_only = c("Mean.Coverage"), 
  Fragmentomics_prop_short_only = c("Proportion.Short"), 
  Fragmentomics_tumor_fraction_only = c("WGS_Tumor_Fraction_Blood_plasma_cfDNA")
)

# 1) Identify which predictors belong to BM vs Blood combos
# BM combos: anything starting “BM_” but drop “BM_blood_combo”
bm_names <- names(combos_small)[
  startsWith(names(combos_small), "BM_") &
    names(combos_small) != "BM_blood_combo"
]
combos_bm    <- combos_small[bm_names]

# Blood combos: anything starting “Blood_” OR containing “Just_fragmentomics”
blood_names <- names(combos_small)[
  startsWith(names(combos_small), "Blood_") |
    grepl("Fragmentomics_", names(combos_small))
]
combos_blood <- combos_small[blood_names]

# Combos just fragmentomics - train seperately to get bigger cohort 
fragmentomics_names <- names(combos_small)[
    grepl("Fragmentomics_", names(combos_small))
]
combos_fragmentomics <- combos_small[fragmentomics_names]


# 2) Build train_bm / train_blood by dropping any row missing *any* BM or Blood predictors
bm_preds    <- unique(unlist(combos_bm))
blood_preds <- unique(unlist(combos_blood))

## Add filter that evidence of disease at baseline 
# 2) Which cohort patients have a “good” baseline BM sample?
Good_pts <- read.csv("baseline_high_quality_patients_updated.csv",
                     stringsAsFactors = FALSE) ## From 2_0 script

bm_good_patients <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_BM_cells == 1) %>%
  select(Patient) %>%         # keep as a tibble with a Patient column
  distinct()                  # remove duplicates

cfDNA_good_patients <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
  select(Patient) %>%         # keep as a tibble with a Patient column
  distinct()                  # remove duplicates

# train_bm_original <- train_bm
# train_blood_original <- train_blood
# hold_bm_original <- hold_bm
# hold_blood_original <- hold_blood
# your data splits
train_df <- dat_mrd %>% filter(Cohort == "Frontline") 
hold_df  <- dat_mrd %>% filter(Cohort == "Non-frontline")

train_df <- train_df %>% filter(!is.na(MRD_truth))
hold_df <- hold_df %>% filter(!is.na(MRD_truth))

# convert binary to factor with levels "neg"/"pos"
train_df <- train_df %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))
hold_df  <- hold_df  %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))

positive_class <- "pos"
train_bm    <- train_df %>% drop_na(all_of(c("MRD_truth", bm_preds))) %>% filter(Patient %in% bm_good_patients$Patient)
train_blood <- train_df %>% drop_na(all_of(c("MRD_truth", blood_preds))) %>% filter(Patient %in% cfDNA_good_patients$Patient)

hold_bm    <- hold_df %>% drop_na(all_of(c("MRD_truth", bm_preds)))  %>% filter(Patient %in% bm_good_patients$Patient)
hold_blood <- hold_df %>% drop_na(all_of(c("MRD_truth", blood_preds))) %>% filter(Patient %in% cfDNA_good_patients$Patient)


### see what has blood but not BM 
# Restrict to frontline samples in each training set
fb <- train_blood %>% filter(Cohort == "Frontline")
fbm <- train_bm    %>% filter(Cohort == "Frontline")

# 1) Identify blood samples with no BM call
no_bm <- fb %>%
  # drop any that have a matching (Patient, Timepoint) in the BM table
  anti_join(fbm, by = c("Patient", "Timepoint")) %>%
  # keep only non‐baseline/diagnosis timepoints
  filter(! timepoint_info %in% c("Baseline", "Diagnosis")) 


n_samples_no_bm  <- nrow(no_bm)
n_patients_no_bm <- no_bm %>% distinct(Patient) %>% nrow()

# Summarize
tibble(
  samples_without_bm   = n_samples_no_bm,
  patients_without_bm  = n_patients_no_bm
)

## Get reason 
failures <- read_csv("Table for creating sample flowchart updated3.csv")
no_bm <- left_join(no_bm)

patients_no_bm <- no_bm %>%
  select(Patient) %>%
  distinct()

failures_joined <- patients_no_bm %>%
  left_join(failures, by = "Patient")

failures_joined

## Get total number of samples
dat_mrd %>%
  filter(Cohort == "Frontline") %>%
  drop_na(all_of(c(blood_preds))) %>%
  filter(Patient %in% cfDNA_good_patients$Patient) %>%
  summarise(
    total_patients = n_distinct(Patient),
    total_samples  = dplyr::n()
  )


# 3) Build train/hold sets for fragmentomics-only models
frag_preds       <- unique(unlist(combos_fragmentomics))

train_fragmentomics <- train_df %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

hold_fragmentomics  <- hold_df %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

## RUN the nested cross-validation
# Outer 5-fold nested CV
run_nested_with_validation <- function(train_data,
                                       valid_data,           # ◀ NEW: pass your hold-out cohort here
                                       combo_list,
                                       positive_class = "pos") {
  # 0) Make sure your `ctrl` (trainControl) is defined in the parent environment
  #    exactly as before (5×5 repeated CV, classProbs=TRUE, summaryFunction=twoClassSummary, etc.)
  
  # 1) Outer 5‐fold indices (stratified on MRD_truth)
  outer_folds <- createFolds(train_data$MRD_truth,
                             k = 5,
                             returnTrain = TRUE)
  
  outer_preds <- list()    # one element per combo
  
  # 2) Loop over each feature-combo
  results <- purrr::imap(combo_list, function(preds, label) {
    
    combo_fold_preds <- vector("list", length(outer_folds))  # collect per fold
    
    # 2a) Collect per‐fold metrics
    combo_metrics <- purrr::map_dfr(seq_along(outer_folds), function(i) {
      # split out vs. in
      train_out <- train_data[ outer_folds[[i]], ]
      test_out  <- train_data[-outer_folds[[i]], ]
      
      # complete‐case
      train_cc <- train_out %>% drop_na(MRD_truth, all_of(preds))
      test_cc  <- test_out  %>% drop_na(MRD_truth, all_of(preds))
      
      # Export diagnostics
      message(sprintf("[%s | fold %d] n_train=%d  n_test=%d  classes(train)=%d  classes(test)=%d",
                      label, i,
                      nrow(train_cc), nrow(test_cc),
                      n_distinct(train_cc$MRD_truth),
                      n_distinct(test_cc$MRD_truth)))
      
      # SKIP if either side is too small or single‐class
    #  if ( # nrow(train_cc) < 20 || removed this criteria
    #      n_distinct(train_cc$MRD_truth) < 2 ||
    #      n_distinct(test_cc$MRD_truth)  < 2) {
    #    return(NULL)
    #  }
      
      if (n_distinct(train_cc$MRD_truth) < 2 ||
          n_distinct(test_cc$MRD_truth)  < 2) {
        message("‑‑ fold dropped: missing class after drop_na()")
        return(NULL)
      }
      
      message("‑‑ fold kept: fitting model now")
      
      # build formula
      f <- reformulate(preds, response = "MRD_truth")
      
      # 2b) inner CV via caret::train
      fit <- if (length(preds) == 1) {
        train(f, data    = train_cc,
              method  = "glm",
              family  = binomial,
              metric  = "ROC",
              trControl = ctrl,
              preProc = c("center","scale"))
      } else {
        train(f, data    = train_cc,
              method    = "glmnet",
              metric    = "ROC",
              # ◀ EXPANDED GRID over alpha ∈ {0,0.25,0.5,0.75,1} and 50 log‐lambda values
              tuneGrid  = expand.grid(
                alpha  = seq(0,1, by=0.25),
                lambda = 10^seq(-3, 1, length.out = 50)
              ),
              trControl = ctrl,
              preProc   = c("center","scale"))
      }
      
      # 2c) predict on outer test
      probs    <- predict(fit, newdata = test_cc, type = "prob")[[positive_class]]
      message("    --> length(probs) = ", length(probs))
      
      roc_o    <- pROC::roc(test_cc$MRD_truth, probs,
                            quiet  = TRUE,
                            levels = c("neg","pos"))
      youden_t <- unlist(
        pROC::coords(roc_o, "best",
                     best.method = "youden",
                     ret         = "threshold")
      )
      youden_t <- as.numeric(youden_t[1])
      
      ## --- store the raw predictions for later plotting ----
      combo_fold_preds[[i]] <<- tibble(
        combo = label,
        fold  = i,
        truth = test_cc$MRD_truth,
        prob  = probs
      )
      
      # pull out single‐value sens/spec
      sens_o <- pROC::coords(roc_o, x = youden_t,
                             ret       = "sensitivity",
                             transpose = FALSE)[[1]]
      spec_o <- pROC::coords(roc_o, x = youden_t,
                             ret       = "specificity",
                             transpose = FALSE)[[1]]
      
      # sensitivity at 95% specificity
      sens95_o <- pROC::coords(roc_o,
                               x      = 0.95,
                               input  = "specificity",
                               ret    = "sensitivity",
                               transpose = FALSE)[[1]]
      
      # specificity at 95% sensitivity
      spec95_o <- pROC::coords(roc_o,
                               x      = 0.95,
                               input  = "sensitivity",
                               ret    = "specificity",
                               transpose = FALSE)[[1]]
      
      # PPV / NPV at Youden threshold
      ppv_o    <- pROC::coords(roc_o,
                               x      = youden_t,
                               input  = "threshold",
                               ret    = "ppv",
                               transpose = FALSE)[[1]]
      npv_o    <- pROC::coords(roc_o,
                               x      = youden_t,
                               input  = "threshold",
                               ret    = "npv",
                               transpose = FALSE)[[1]]
      
      # F1 and balanced accuracy
      f1_o     <- 2 * ppv_o * sens_o / (ppv_o + sens_o)
      balacc_o <- (sens_o + spec_o) / 2
      
      # Brier score
      brier_o  <- mean((probs - as.numeric(test_cc$MRD_truth == positive_class))^2)
      
      # partial AUC over specificity ≥ 0.9
      pauc90_o <- as.numeric(
        pROC::auc(roc_o,
                  partial.auc = c(0.9, 1),
                  partial.auc.correct = TRUE)
      )
      
      # fold sizes
      n_train_i <- nrow(train_cc)
      n_test_i  <- nrow(test_cc)
      
      tibble(
        combo       = label,
        auc         = as.numeric( pROC::auc(roc_o) ),
        sensitivity = sens_o,
        specificity = spec_o,
        sens_at_95_spec       = sens95_o,
        spec_at_95_sen       = spec95_o,
        bal_accuracy  = balacc_o,
        ppv           = ppv_o,
        npv           = npv_o,
        f1            = f1_o,
        brier         = brier_o,
        pAUC90        = pauc90_o,
        n_train       = n_train_i,
        n_test        = n_test_i,
        accuracy    = mean((probs >= youden_t) ==
                             (test_cc$MRD_truth == positive_class))
      )
    })
    
    # --- NEW: bind the fold‑level predictions for this combo
    preds_df <- bind_rows(combo_fold_preds)
    message("   --> preds_df rows for ", label, " = ", nrow(preds_df))
    
    # 2d) summarize across *kept* folds
    nested_summary <- combo_metrics %>%
      summarise(
        combo            = first(combo),
        n_folds          = dplyr::n(),
        auc_mean         = mean(auc,           na.rm = TRUE),
        auc_sd           = sd(auc,             na.rm = TRUE),
        sens_mean        = mean(sensitivity,   na.rm = TRUE),
        sens_sd          = sd(sensitivity,     na.rm = TRUE),
        spec_mean        = mean(specificity,   na.rm = TRUE),
        spec_sd          = sd(specificity,     na.rm = TRUE),
        balacc_mean      = mean(bal_accuracy,   na.rm = TRUE),
        balacc_sd        = sd(bal_accuracy,     na.rm = TRUE),
        ppv_mean         = mean(ppv,            na.rm = TRUE),
        npv_mean         = mean(npv,            na.rm = TRUE),
        f1_mean          = mean(f1,             na.rm = TRUE),
        brier_mean       = mean(brier,          na.rm = TRUE),
        pAUC90_mean      = mean(pAUC90,         na.rm = TRUE),
        sens95_mean      = mean(sens_at_95_spec,         na.rm = TRUE),
        spec95_mean      = mean(spec_at_95_sen,         na.rm = TRUE),
        mean_n_train     = mean(n_train,        na.rm = TRUE),
        mean_n_test      = mean(n_test,         na.rm = TRUE),
        acc_mean         = mean(accuracy,       na.rm = TRUE),
        .groups          = "drop"
      )
    
    # 3) Re‐fit on the *entire* train_data
    full_cc   <- train_data %>% drop_na(MRD_truth, all_of(preds))
    f_full    <- reformulate(preds, "MRD_truth")
    final_fit <- if (length(preds) == 1) {
      train(f_full, data    = full_cc,
            method  = "glm", family = binomial,
            metric  = "ROC", trControl = ctrl,
            preProc = c("center","scale"))
    } else {
      train(f_full, data    = full_cc,
            method    = "glmnet",
            metric    = "ROC",
            tuneGrid  = expand.grid(
              alpha  = seq(0,1, by=0.25),
              lambda = 10^seq(-3, 1, length.out = 10)
            ),
            trControl = ctrl,
            preProc   = c("center","scale"))
    }
    
    # 4) Compute Youden threshold on full‐train ROC
    full_probs <- predict(final_fit, newdata = full_cc, type="prob")[[positive_class]]
    roc_full   <- pROC::roc(full_cc$MRD_truth, full_probs,
                            quiet  = TRUE,
                            levels = c("neg","pos"))
    youden_full <- unlist(
      pROC::coords(roc_full, "best",
                   best.method = "youden",
                   ret         = "threshold")
    )
    youden_full <- as.numeric(youden_full[1])
    
    # keep that threshold so we can reuse it later 
    threshold_full <- youden_full
    
    # 5) EVALUATE on the *external* validation set
    valid_cc   <- valid_data %>% drop_na(MRD_truth, all_of(preds))
    valid_probs<- predict(final_fit, newdata = valid_cc, type="prob")[[positive_class]]
    roc_val    <- pROC::roc(valid_cc$MRD_truth, valid_probs,
                            quiet  = TRUE,
                            levels = c("neg","pos"))
    sens_val   <- pROC::coords(roc_val, x = youden_full,
                               ret = "sensitivity")[[1]]
    spec_val   <- pROC::coords(roc_val, x = youden_full,
                               ret = "specificity")[[1]]
    acc_val    <- mean((valid_probs >= youden_full) ==
                         (valid_cc$MRD_truth == positive_class))
    
    # sensitivity/spec at 95% on the hold‐out
    sens95_val <- pROC::coords(roc_val,
                               x      = 0.95,
                               input  = "specificity",
                               ret    = "sensitivity",
                               transpose = FALSE)[1]
    
    spec95_val <- pROC::coords(roc_val,
                               x      = 0.95,
                               input  = "sensitivity",
                               ret    = "specificity",
                               transpose = FALSE)[1]
    
    # PPV / NPV at the full-train Youden cutoff
    ppv_val    <- pROC::coords(roc_val,
                               x      = youden_full,
                               input  = "threshold",
                               ret    = "ppv",
                               transpose = FALSE)[1]
    npv_val    <- pROC::coords(roc_val,
                               x      = youden_full,
                               input  = "threshold",
                               ret    = "npv",
                               transpose = FALSE)[1]
    
    # F1 and balanced accuracy
    f1_val     <- 2 * ppv_val * sens_val / (ppv_val + sens_val)
    balacc_val <- (sens_val + spec_val) / 2
    
    # Brier score
    brier_val  <- mean((valid_probs - as.numeric(valid_cc$MRD_truth == positive_class))^2)
    
    # Partial AUC over specificity ≥ 0.9
    pauc90_val <- as.numeric(
      pROC::auc(roc_val,
                partial.auc           = c(0.9, 1),
                partial.auc.correct   = TRUE)
    )
    
    
    validation_summary <- tibble(
      combo       = label,
      auc_valid   = as.numeric( pROC::auc(roc_val) ),
      sens_valid  = sens_val,
      spec_valid  = spec_val,
      sens_at_95_spec_valid  = sens95_val,
      spec_at_95_sens_valid  = spec95_val,
      bal_accuracy   = balacc_val,
      ppv            = ppv_val,
      npv            = npv_val,
      f1             = f1_val,
      brier          = brier_val,
      pAUC90         = pauc90_val,
      accuracy   = acc_val
    )
    
    # return all pieces
    list(
      nested     = nested_summary,
      final_fit  = final_fit,
      bestTune   = final_fit$bestTune,
      validation = validation_summary,
      threshold  = threshold_full,
      preds_df   = preds_df           # <‑‑  NEW: stash the preds in the result
    )
  })
  
  # pull them back out
  nested_metrics     <- map_dfr(results, "nested")
  models             <- map(results, "final_fit")
  best_tunes         <- map(results, "bestTune")
  validation_metrics <- map_dfr(results, "validation")
  thresholds         <- map_dbl(results, "threshold") # collect thresholds
  outer_predictions  <- map_dfr(results, "preds_df") # ---  NEW: bind all combos’ outer‑fold predictions
  
  message("Rows in outer_predictions = ", nrow(outer_predictions))
  
  list(
    nested_metrics     = nested_metrics,
    models             = models,
    best_tunes         = best_tunes,
    validation_metrics = validation_metrics,
    thresholds         = thresholds,
    outer_predictions  = outer_predictions   ### NEW  –– outer predictions are here
  )
}


# 4) Run nested CV separately
set.seed(2025)
nested_bm_validation_updated3 <- run_nested_with_validation(
  train_data = train_bm,
  valid_data = hold_bm,
  combo_list = combos_bm,
  positive_class = "pos"
)

set.seed(2025)
nested_blood_validation_updated3 <- run_nested_with_validation(
  train_data = train_blood,
  valid_data = hold_blood,
  combo_list = combos_blood,
  positive_class = "pos"
)

set.seed(2025)
nested_fragmentomics_validation_updated3 <- run_nested_with_validation(
  train_data    = train_fragmentomics,
  valid_data    = hold_fragmentomics,
  combo_list    = combos_fragmentomics,
  positive_class = "pos"
)


bm_preds <- nested_bm_validation_updated3$outer_predictions


## For testing - just when need quick iteration on script 
# 1) Minimal combos list with only BM_base
combos_small_test <- list(
  BM_base = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM")
)

# 2) For clarity we’ll treat all entries in combos_small_test as BM combos
combos_bm_test <- combos_small_test

# 3) Pull out the predictor names and build train / hold sets
bm_preds_test   <- unique(unlist(combos_bm_test))
train_bm_test   <- train_df %>% drop_na(all_of(c("MRD_truth", bm_preds_test)))
hold_bm_test    <- hold_df  %>% drop_na(all_of(c("MRD_truth", bm_preds_test)))

# 4) Run nested CV just on BM_base
test_out <- run_nested_with_validation(
  train_data     = train_bm_test,
  valid_data     = hold_bm_test,
  combo_list     = combos_bm_test,
  positive_class = positive_class
)

# 5) Export
## Export the models that are updated
# saveRDS(nested_bm_validation_updated2, file = "nested_bm_validation_updated2.rds")
# saveRDS(nested_blood_validation_updated2, file = "nested_blood_validation_updated2.rds")
# saveRDS(nested_fragmentomics_validation_updated2, file = "nested_fragmentomics_validation_updated2.rds")

## Export the models that are updated
saveRDS(nested_bm_validation_updated3, file = "nested_bm_validation_updated3.rds")
saveRDS(nested_blood_validation_updated3, file = "nested_blood_validation_updated3_run2.rds")
saveRDS(nested_fragmentomics_validation_updated3, file = "nested_fragmentomics_validation_updated3.rds")

## Load back in (optional)
# nested_bm_validation_updated3 <- readRDS("nested_bm_validation_updated3.rds")
# nested_blood_validation_updated3 <- readRDS("nested_blood_validation_updated3.rds")
# nested_fragmentomics_validation_updated3 <- readRDS("nested_fragmentomics_validation_updated3.rds")


## For consistency 
nested_bm_validation_updated2 <- nested_bm_validation_updated3
nested_blood_validation_updated2 <- nested_blood_validation_updated3
nested_fragmentomics_validation_updated2 <- nested_fragmentomics_validation_updated3

# -----------------------------------------------------------------------------
# 8.Get Model Metrics Together and Exported
# -----------------------------------------------------------------------------

# All validation metrics in one tibble
all_val <- bind_rows(
  nested_bm_validation_updated2$validation_metrics,
  nested_blood_validation_updated2$validation_metrics
)

## All primary
all_primary <- bind_rows(
  nested_bm_validation_updated2$nested_metrics,
  nested_blood_validation_updated2$nested_metrics
)

# One named list of all models
all_models <- c(
  nested_bm_validation_updated2$models,
  nested_blood_validation_updated2$models
)

# One named vector of the corresponding thresholds
all_thresholds <- c(
  nested_bm_validation_updated2$thresholds,
  nested_blood_validation_updated2$thresholds
)


### Remove and re-add fragmentomics from full training 
# the fragmentomics combo names
frag_names <- names(combos_fragmentomics)

# 1) Filter out fragmentomics rows from the blood‐based tibbles
all_val_clean <- all_val %>%
  filter(!combo %in% frag_names)

all_primary_clean <- all_primary %>%
  filter(!combo %in% frag_names)

# 2) Bind in the fragmentomics validation and nested metrics
all_val <- bind_rows(
  all_val_clean,
  nested_fragmentomics_validation_updated2$validation_metrics
)

all_primary <- bind_rows(
  all_primary_clean,
  nested_fragmentomics_validation_updated2$nested_metrics
)

# 3) Rebuild the models list: drop fragmentomics from the old and add the new
all_models_clean <- all_models[ ! names(all_models) %in% frag_names ]

all_models <- c(
  all_models_clean,
  nested_fragmentomics_validation_updated2$models
)

## Add additional metrics

# ---- 1. gather all outer-fold predictions -----------------------------------
# 1) pool everything *except* fragmentomics combos from the BM / blood sets
outer_preds <- bind_rows(
  # keep all fragmentomics
  nested_fragmentomics_validation_updated2$outer_predictions,
  # drop any frag_names from BM before binding
  nested_bm_validation_updated2$outer_predictions %>%
    filter(! combo %in% frag_names),
  # drop any frag_names from blood before binding
  nested_blood_validation_updated2$outer_predictions %>%
    filter(! combo %in% frag_names)
) %>%
  mutate(truth = factor(truth, levels = c("neg","pos"))) # make sure the response is a factor with neg first, pos second

# ---- 2. pooled ROC, AUC and 95 % CI per combo ------------------------------
pooled_auc_tbl <- outer_preds %>%
  group_by(combo) %>%
  summarise(
    roc_obj     = list(roc(truth, prob,
                           levels    = c("neg","pos"),
                           direction = "<",
                           quiet     = TRUE)),
    auc_pooled  = as.numeric(auc(roc_obj[[1]])),
    ci_low      = ci.auc(roc_obj[[1]], method = "delong")[1],
    ci_high     = ci.auc(roc_obj[[1]], method = "delong")[3],
    .groups     = "drop"
  ) %>%
  mutate(
    auc_pooled_ci = sprintf("%.2f (%.2f–%.2f)", auc_pooled, ci_low, ci_high)
  ) %>%
  select(-roc_obj)   # drop the heavy objects

# ---- 3. add to your primary metrics table ----------------------------------
all_primary <- all_primary %>%
  left_join(pooled_auc_tbl, by = "combo")


## Pick winners after reviewin the tables
wanted <- c(
  # BM-only models
  "BM_zscore_only_sites",              # best BM by AUC
  "BM_zscore_only_detection_rate",     # best BM by sensitivity
  "BM_base_zscore",     # best BM by sensitivity
  
  
  # Blood-only models
  "Blood_base",   #  best blood by AUC
  "Blood_rate_only",                   # second best blood by AUC
  "Blood_plus_fragment_min",           # best blood by accuracy
  "Blood_zscore_only_detection_rate",  # best blood by sensitivity
  "Blood_zscore_only_sites",  # best blood in validation cohort
  
  # Fragmentomics-only models
  "Fragmentomics_tumor_fraction_only",             # best fragmentomics in validation cohort and good sens
  "Fragmentomics_mean_coverage_only"   # best auc in primary and accuracy
)

# Then pull just those rows
selected_primary <- all_primary # %>%
#  filter(combo %in% wanted)

## Save
selected_combos   <- selected_primary$combo
selected_models   <- all_models[selected_combos]
selected_thr      <- all_thresholds[selected_combos]

# saveRDS(selected_models,
#         file = file.path(outdir, "selected_combo_models_2025-07-25.rds"))
# saveRDS(selected_thr,
#         file = file.path(outdir, "selected_combo_thresholds_2025-07-25.rds"))
# 
# saveRDS(selected_models,
#         file = file.path(outdir, "selected_combo_models_2025-07-31.rds"))
# saveRDS(selected_thr,
#         file = file.path(outdir, "selected_combo_thresholds_2025-07-31.rds"))


saveRDS(selected_models,
        file = file.path(outdir, "selected_combo_models_2025-08-15.rds"))
saveRDS(selected_thr,
        file = file.path(outdir, "selected_combo_thresholds_2025-08-15.rds"))


# -----------------------------------------------------------------------------
# 8b. Export Combined Metrics & Models
# -----------------------------------------------------------------------------
# 1) Validation‐set metrics (hold‐out cohort)
saveRDS(all_val,
        file = file.path(outdir, "all_validation_metrics.rds"))
write.csv(all_val,
          file = file.path(outdir, "all_validation_metrics.csv"),
          row.names = FALSE)

# 2) Nested‐CV metrics (primary/frontline cohort)
saveRDS(all_primary,
        file = file.path(outdir, "all_nested_cv_metrics.rds"))
write.csv(all_primary,
          file = file.path(outdir, "all_nested_cv_metrics.csv"),
          row.names = FALSE)

# 3) Full list of fitted models
saveRDS(all_models,
        file = file.path(outdir, "all_cfWGS_models_list.rds"))

# (Optional) If you also want to export the thresholds vector:
saveRDS(all_thresholds,
        file = file.path(outdir, "all_model_thresholds.rds"))
write.csv(
  tibble(
    combo     = names(all_thresholds),
    threshold = unname(all_thresholds)
  ),
  file = file.path(outdir, "all_model_thresholds.csv"),
  row.names = FALSE
)

message("Exported validation metrics, nested‐CV metrics, models list, and thresholds to ", outdir)


### Export cleaned version for supplementary tables
# 1) pull out only the non‑all‑NA columns from your training results
primary_clean <- all_primary %>%
  select(where(~ !all(is.na(.)))) %>%   # drop any column that is entirely NA
  mutate(cohort = "Training")           # tag it

# 2) do the same for the validation set, and rename its cols to match the *_mean convention
val_clean <- all_val %>%
  select(where(~ !all(is.na(.)))) %>%
  rename(
    auc_mean       = auc_valid,
    sens_mean      = sens_valid,
    spec_mean      = spec_valid,
    balacc_mean    = bal_accuracy,
    ppv_mean       = ppv,
    npv_mean       = npv,
    f1_mean        = f1,
    brier_mean     = brier,
    pAUC90_mean    = pAUC90,
    acc_mean       = accuracy,
    sens95_mean    = sens_at_95_spec_valid,
    spec95_mean    = spec_at_95_sens_valid
  ) %>%
  mutate(cohort = "Testing")



# 1. Helper to pull the first (and only) column out of any data.frame columns
flatten_df_cols <- function(df) {
  is_df <- vapply(df, is.data.frame, logical(1))
  df[is_df] <- lapply(df[is_df], function(x) x[[1]])
  df
}

# 2. Helper to clean up the sens95/spec95 names
fix_names <- function(df) {
  names(df) <- sub("sens95_mean\\$sensitivity", "sens95_mean", names(df))
  names(df) <- sub("spec95_mean\\$specificity", "spec95_mean", names(df))
  df
}

# Apply to both data sets
primary_clean <- primary_clean %>%
  flatten_df_cols() %>%
  fix_names()

val_clean <- val_clean %>%
  flatten_df_cols() %>%
  fix_names()

# 3) stack them into one combined table
combined_metrics <- bind_rows(primary_clean, val_clean) %>%
  select(
    -sens95_mean,
    -spec95_mean,
    -ppv_mean,
    -npv_mean,
    -f1_mean
  )

# 4) Export only the Fragmentomics models
frag_tbl <- combined_metrics %>%
  filter(startsWith(combo, "Fragmentomics_"))

write_csv(
  frag_tbl,
  "Output_tables_2025/Supplementary_Table_4_Fragmentomics_Performance3.csv"
)

# 2) Export all models
write_csv(
  primary_clean %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))),
  "Final Tables and Figures/Supplementary_Table_3_All_Model_performance_nested_CV_updated3.csv"
)

write_csv(
  val_clean %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))),
  "Final Tables and Figures/Supplementary_Table_5_All_Model_performance_nested_CV_updated_on_testing_cohort3_used_inMS.csv"
)




##### Now see sensetivity at 95% specificity in the all models file 
# all_models is your named list of caret::train objects
sens_at_95spec <- map_dbl(all_models, function(model) {
  # 1) Grab the best‐tuning parameters
  best <- model$bestTune
  
  # 2) Subset the out‐of‐fold predictions to only those rows matching bestTune
  preds_best <- inner_join(
    model$pred,
    best,
    by = names(best)
  )
  
  # 3) Build ROC on those held‐out probabilities
  roc_obj <- roc(
    response  = preds_best$obs,    # observed neg/pos
    predictor = preds_best$pos,    # predicted probability of “pos”
    levels    = c("neg", "pos"),
    direction = "<"
  )
  
  # 4) Extract the sensitivity at 95% specificity
  sens <- coords(
    roc_obj,
    x         = 0.95,
    input     = "specificity",
    ret       = "sensitivity",
    transpose = FALSE
  )
  
  as.numeric(sens)
})

# Combine into a dataframe for reporting
tibble(
  model            = names(sens_at_95spec),
  sens_at_95_spec  = sens_at_95spec
)


# Compute specificity at 95% sensitivity for each model
spec_at_95sen <- map_dbl(all_models, function(model) {
  # 1) Optimal tuning parameters
  best <- model$bestTune
  
  # 2) Subset to out-of-fold predictions at bestTune
  preds_best <- inner_join(
    model$pred,
    best,
    by = names(best)
  )
  
  # 3) Build ROC on those held-out predictions
  roc_obj <- roc(
    response  = preds_best$obs,   # true class: “neg”/“pos”
    predictor = preds_best$pos,   # predicted probability for “pos”
    levels    = c("neg","pos"),
    direction = "<"
  )
  
  # 4) Extract specificity when sensitivity = 0.95
  spec <- coords(
    roc_obj,
    x         = 0.95,
    input     = "sensitivity",
    ret       = "specificity",
    transpose = FALSE
  )
  
  as.numeric(spec)
})

# Create a summary table
table <- tibble(
  model                  = names(spec_at_95sen),
  spec_at_95_sensitivity = spec_at_95sen
)



# -----------------------------------------------------------------------------
# 9A. Full-Cohort & Prep for Dilution-Series Scoring
# -----------------------------------------------------------------------------
### Apply on entire dataframe 
## Define function
apply_selected <- function(dat, models, thresholds, positive_class = "pos") {
  out <- dat            # start with your full data.frame
  n   <- nrow(dat)
  
  # Check that models and thresholds align
  if (!all(names(models) %in% names(thresholds))) {
    stop("Model names and threshold names do not match!")
  }
  common <- intersect(names(models), names(thresholds))
  
  for (cmb in common) {
    fit <- models[[cmb]]
    thr <- thresholds[[cmb]]
    
    # which predictors the model actually needs
    preds <- setdiff(names(fit$trainingData), ".outcome")
    
    # figure out which rows we can actually score
    can_score <- complete.cases(dat[preds])
    
    # pre‐allocate
    probs <- rep(NA_real_,    n)
    calls <- rep(NA_integer_, n)
    
    # only predict where allowed
    if (any(can_score)) {
      probs[can_score] <- predict(fit,
                                  newdata = dat[can_score, preds, drop = FALSE],
                                  type    = "prob")[[positive_class]]
      calls[can_score] <- as.integer(probs[can_score] >= thr)
    }
    
    # stick them back on
    out[[paste0(cmb, "_prob")]] <- probs
    out[[paste0(cmb, "_call")]] <- calls
  }
  
  out
}

## Run function
### Primary cohort
data_scored <- apply_selected(
  dat           = dat,            # your big master data frame
  models        = selected_models,
  thresholds    = selected_thr,
  positive_class= "pos"
)


### Now restrict to elidgible patients 
## Patient eligibility vectors from earlier 
good_blood_patients <- cfDNA_good_patients$Patient
good_bm_patients    <- bm_good_patients$Patient

## All probability and score columns produced by apply_selected
score_cols <- grep("_(prob|call)$", names(data_scored), value = TRUE, ignore.case = TRUE)

## Map prob column -> model base name (strip the trailing "_prob")
model_base <- tibble(
  model_col  = score_cols,
  model_base = str_remove(score_cols, "_(?i:prob|call)$")   # strip _prob or _call
)

## Heuristic family mapping based on the names you provided
model_map <- model_base %>%
  mutate(family = case_when(
    str_detect(model_base, regex("^BM_", ignore_case = TRUE))            ~ "BM",
    str_detect(model_base, regex("^Blood_", ignore_case = TRUE))         ~ "BLOOD",
    str_detect(model_base, regex("^Fragmentomics_", ignore_case = TRUE)) ~ "FRAG",
    TRUE                                                                 ~ "OTHER"
  ))

# --- 4) Mask disallowed patients for BOTH prob+call columns --------
data_scored_masked <- data_scored

typed_na <- function(x) {
  if (is.factor(x))         return(NA)                 # factor NA
  if (is.logical(x))        return(NA)                 # logical NA
  if (is.integer(x))        return(NA_integer_)        # integer NA
  if (is.double(x))         return(NA_real_)           # numeric NA
  if (is.numeric(x))        return(NA_real_)
  if (is.character(x))      return(NA_character_)      # character NA
  if (inherits(x, "Date"))  return(as.Date(NA))        # Date NA
  if (inherits(x, "POSIXt"))return(as.POSIXct(NA))     # POSIXct NA
  NA
}

walk2(model_map$model_col, model_map$family, function(col, fam) {
  allowed <- switch(
    fam,
    "BM"    = good_bm_patients,
    "BLOOD" = good_blood_patients,
    "FRAG"  = unique(data_scored$Patient),
    "OTHER" = unique(data_scored$Patient)
  )
  
  x <- data_scored[[col]]
  mask <- data_scored$Patient %in% allowed
  x[!mask] <- typed_na(x)        # preserve the column's original type
  data_scored_masked[[col]] <<- x
})

# --- 5) Simple diagnostics: effective N before/after --------------
message("=== Effective N BEFORE masking ===")
for (m in score_cols) {
  n_eff <- data_scored %>% tidyr::drop_na(all_of(m)) %>% nrow()
  message(sprintf("%-40s N = %d", m, n_eff))
}

message("=== Effective N AFTER masking ===")
for (m in score_cols) {
  n_eff <- data_scored_masked %>% tidyr::drop_na(all_of(m)) %>% nrow()
  message(sprintf("%-40s N = %d", m, n_eff))
}


### Save this 
# saveRDS(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated3.rds"))
# write_csv(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated3.csv"))

## Another export
saveRDS(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated4_full.rds"))
write_csv(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated4_full.csv"))

saveRDS(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated4.rds"))
write_csv(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated4.csv"))





# -----------------------------------------------------------------------------
# 9B. Now re-get model metrics on the whole cohort 
# -----------------------------------------------------------------------------

## See what the metrics look like at 95% specificity 

# 1) Filter to frontline cohort and compute prevalence
frontline <- data_scored_masked %>% 
  filter(Cohort == "Frontline") %>% 
  filter(! timepoint_info %in% c("Baseline", "Diagnosis")) 

test_cohort <- data_scored_masked %>% 
  filter(Cohort == "Non-frontline") %>%
  filter(! timepoint_info %in% c("Baseline", "Diagnosis")) 

prev <- mean(frontline$MRD_truth == 1, na.rm = TRUE)

# 2) Identify all “_prob” columns (i.e. your different models)
prob_cols <- grep("_prob$", colnames(data_scored), value = TRUE)

### 3A) First get metrics for entire cohort 
# 3) Loop over models, build ROC & compute metrics at every threshold
all_metrics_rescored_primary <- map_dfr(prob_cols, function(prob_col) {
  # build ROC object
  roc_obj <- roc(
    response  = frontline$MRD_truth,
    predictor = frontline[[prob_col]],
    levels    = c(0, 1),
    direction = "<"
  )
  
  # extract thresholds with sens, spec, ppv, npv
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  ) %>% 
    as_tibble()
  
  # compute additional metrics
  roc_df %>% 
    mutate(
      model        = prob_col,
      accuracy     = sensitivity * prev + specificity * (1 - prev),
      bal_accuracy = (sensitivity + specificity) / 2,
      f1           = 2 * sensitivity * ppv / (sensitivity + ppv)
    ) %>%
    select(model, threshold, sensitivity, specificity, ppv, npv, accuracy, bal_accuracy, f1)
})

# Get Youden cutoff
threshold_df <- tibble(
  model = paste0(names(all_thresholds), "_prob"),
  youden = unname(all_thresholds)
)

# Filter your all_metrics_rescored_primary to only the Youden‐index row for each model
metrics_youden <- all_metrics_rescored_primary %>%
  inner_join(threshold_df, by = "model") %>%
  group_by(model) %>%
  # pick the row whose threshold minimizes |threshold - youden|
  slice_min(order_by = abs(threshold - youden), n = 1) %>%
  ungroup() %>%
  select(-youden)

print(metrics_youden)

### 3B) Now for each model, find the threshold with ≥95% specificity that maximizes sensitivity,
#    then compute sensitivity, specificity and accuracy at that threshold
metrics_at_95spec <- map_dfr(prob_cols, function(prob_col) {
  # build ROC
  roc_obj <- roc(
    response  = frontline$MRD_truth,
    predictor = frontline[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # get all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x        = "all",
    ret      = c("threshold","sensitivity","specificity"),
    transpose= FALSE
  ) %>% as_tibble()
  
  # filter to spec ≥ 0.95
  cand <- roc_df %>%
    filter(specificity >= 0.95)
  
  # if no threshold meets the bar, return NA
  if (nrow(cand) == 0) {
    return(tibble(
      model        = prob_col,
      threshold    = NA_real_,
      sensitivity  = NA_real_,
      specificity  = NA_real_,
      accuracy     = NA_real_
    ))
  }
  
  # pick the one with highest sensitivity
  best <- cand %>%
    arrange(desc(sensitivity)) %>%
    slice(1)
  
  # compute accuracy = sens*prev + spec*(1-prev)
  best <- best %>%
    mutate(
      accuracy = sensitivity * prev + specificity * (1 - prev),
      model    = prob_col
    ) %>%
    select(model, threshold, sensitivity, specificity, accuracy)
  
  return(best)
})

print(metrics_at_95spec)

## Now for sens


### Get more metrics 
# Prevalence: proportion of positives in data
prev <- mean(frontline$MRD_truth == 1, na.rm = TRUE)

metrics_at_95sens <- map_dfr(prob_cols, function(prob_col) {
  # 1) Build ROC
  roc_obj <- roc(
    response  = frontline$MRD_truth,
    predictor = frontline[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # 2) Extract all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold","sensitivity","specificity"),
    transpose = FALSE
  ) %>% as_tibble()
  
  # 3) Filter to sens ≥ 0.95
  cand <- roc_df %>%
    filter(sensitivity >= 0.94)
  
  # 4) If none meet the bar, fall back to highest observed sens
  if (nrow(cand) == 0) {
    best <- roc_df %>% slice_max(sensitivity, n = 1)
  } else {
    # 5) Otherwise pick the one with highest specificity
    best <- cand %>% arrange(desc(specificity)) %>% slice(1)
  }
  
  # 6) Compute precision, F1, and balanced accuracy at that threshold
  best <- best %>%
    mutate(
      model      = prob_col,
      accuracy   = sensitivity * prev + specificity * (1 - prev),
      precision  = ifelse(
        (sensitivity * prev + (1 - specificity) * (1 - prev)) == 0,
        NA,
        (sensitivity * prev) / (sensitivity * prev + (1 - specificity) * (1 - prev))
      ),
      f1         = ifelse(
        (precision + sensitivity) == 0,
        NA,
        2 * (precision * sensitivity) / (precision + sensitivity)
      ),
      bal_accuracy = (sensitivity + specificity) / 2
    )
  
  best %>% select(model, threshold, sensitivity, specificity, accuracy, precision, f1, bal_accuracy)
})

print(metrics_at_95sens)


### Save this 
# Save metrics evaluated at Youden threshold
write_rds(metrics_youden, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold.rds"))
write_csv(metrics_youden, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold.csv"))

# Save metrics evaluated at fixed 95% specificity (if applicable)
write_rds(metrics_at_95spec, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec.rds"))
write_csv(metrics_at_95spec, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec.csv"))

write_rds(metrics_at_95sens, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens.rds"))
write_csv(metrics_at_95sens, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens_updated.csv"))


## Tidy into one thing 
# 1. Tidy up and rename columns on each
youden_tbl <- metrics_youden %>%
  rename(
    Model                 = model,
    Threshold_Youden      = threshold,
    Sensitivity_Youden    = sensitivity,
    Specificity_Youden    = specificity,
    PPV_Youden            = ppv,
    NPV_Youden            = npv,
    Accuracy_Youden       = accuracy,
    BalancedAcc_Youden    = bal_accuracy,
    F1_Youden             = f1
  )

spec95_tbl <- metrics_at_95spec %>%
  rename(
    Model                   = model,
    Threshold_95Spec        = threshold,
    Sensitivity_95Spec      = sensitivity,
    Specificity_95Spec      = specificity,
    Accuracy_95Spec         = accuracy
  )

sens95_tbl <- metrics_at_95sens %>%
  rename(
    Model                   = model,
    Threshold_95Sens        = threshold,
    Sensitivity_95Sens      = sensitivity,
    Specificity_95Sens      = specificity,
    Accuracy_95Sens         = accuracy,
    Precision_95Sens        = precision,
    F1_95Sens               = f1,
    BalancedAccuracy_95Sens = bal_accuracy
  )


# 2. Join them all together
all_perf_metrics <- youden_tbl %>%
  left_join(spec95_tbl, by = "Model") %>%
  left_join(sens95_tbl, by = "Model")

all_perf_metrics <- all_perf_metrics %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# 3. Export as a single supplementary Excel (or CSV)
library(writexl)    # or use write.csv if you prefer CSV
write_xlsx(
  list(`All Performance Metrics` = all_perf_metrics),
  path = "Final Tables and Figures/Supplementary_Table_4_All_Model_Metrics_Refit4.xlsx"
)

# — or, if prefer CSV:
write.csv(all_perf_metrics, 
          "Final Tables and Figures/Supplementary_Table_4_All_Model_Metrics_Refit4.csv", 
          row.names = FALSE)




##### Now redo for non-frontline 
### 3A) First get metrics for entire cohort 
# 3) Loop over models, build ROC & compute metrics at every threshold
all_metrics_rescored_testing <- map_dfr(prob_cols, function(prob_col) {
  # build ROC object
  roc_obj <- roc(
    response  = test_cohort$MRD_truth,
    predictor = test_cohort[[prob_col]],
    levels    = c(0, 1),
    direction = "<"
  )
  
  # extract thresholds with sens, spec, ppv, npv
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  ) %>% 
    as_tibble()
  
  # compute additional metrics
  roc_df %>% 
    mutate(
      model        = prob_col,
      accuracy     = sensitivity * prev + specificity * (1 - prev),
      bal_accuracy = (sensitivity + specificity) / 2,
      f1           = 2 * sensitivity * ppv / (sensitivity + ppv)
    ) %>%
    select(model, threshold, sensitivity, specificity, ppv, npv, accuracy, bal_accuracy, f1)
})

# Filter your all_metrics_rescored_primary to only the Youden‐index row for each model
metrics_youden_testing <- all_metrics_rescored_testing %>%
  inner_join(threshold_df, by = "model") %>%
  group_by(model) %>%
  # pick the row whose threshold minimizes |threshold - youden|
  slice_min(order_by = abs(threshold - youden), n = 1) %>%
  ungroup() %>%
  select(-youden)


### 3B) Now for each model, find the threshold with ≥95% specificity that maximizes sensitivity,
#    then compute sensitivity, specificity and accuracy at that threshold
metrics_at_95spec_test <- map_dfr(prob_cols, function(prob_col) {
  # build ROC
  roc_obj <- roc(
    response  = test_cohort$MRD_truth,
    predictor = test_cohort[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # get all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x        = "all",
    ret      = c("threshold","sensitivity","specificity"),
    transpose= FALSE
  ) %>% as_tibble()
  
  # filter to spec ≥ 0.95
  cand <- roc_df %>%
    filter(specificity >= 0.95)
  
  # if no threshold meets the bar, return NA
  if (nrow(cand) == 0) {
    return(tibble(
      model        = prob_col,
      threshold    = NA_real_,
      sensitivity  = NA_real_,
      specificity  = NA_real_,
      accuracy     = NA_real_
    ))
  }
  
  # pick the one with highest sensitivity
  best <- cand %>%
    arrange(desc(sensitivity)) %>%
    slice(1)
  
  # compute accuracy = sens*prev + spec*(1-prev)
  best <- best %>%
    mutate(
      accuracy = sensitivity * prev + specificity * (1 - prev),
      model    = prob_col
    ) %>%
    select(model, threshold, sensitivity, specificity, accuracy)
  
  return(best)
})


## Now for sens

### Get more metrics 
# Prevalence: proportion of positives in data
prev <- mean(test_cohort$MRD_truth == 1, na.rm = TRUE)

metrics_at_95sens_test <- map_dfr(prob_cols, function(prob_col) {
  # 1) Build ROC
  roc_obj <- roc(
    response  = test_cohort$MRD_truth,
    predictor = test_cohort[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # 2) Extract all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold","sensitivity","specificity"),
    transpose = FALSE
  ) %>% as_tibble()
  
  # 3) Filter to sens ≥ 0.95
  cand <- roc_df %>%
    filter(sensitivity >= 0.94)
  
  # 4) If none meet the bar, fall back to highest observed sens
  if (nrow(cand) == 0) {
    best <- roc_df %>% slice_max(sensitivity, n = 1)
  } else {
    # 5) Otherwise pick the one with highest specificity
    best <- cand %>% arrange(desc(specificity)) %>% slice(1)
  }
  
  # 6) Compute precision, F1, and balanced accuracy at that threshold
  best <- best %>%
    mutate(
      model      = prob_col,
      accuracy   = sensitivity * prev + specificity * (1 - prev),
      precision  = ifelse(
        (sensitivity * prev + (1 - specificity) * (1 - prev)) == 0,
        NA,
        (sensitivity * prev) / (sensitivity * prev + (1 - specificity) * (1 - prev))
      ),
      f1         = ifelse(
        (precision + sensitivity) == 0,
        NA,
        2 * (precision * sensitivity) / (precision + sensitivity)
      ),
      bal_accuracy = (sensitivity + specificity) / 2
    )
  
  best %>% select(model, threshold, sensitivity, specificity, accuracy, precision, f1, bal_accuracy)
})

print(metrics_at_95sens_test)


### Save this 
# Save metrics evaluated at Youden threshold
write_rds(metrics_youden_testing, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold_test_cohort.rds"))
write_csv(metrics_youden_testing, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold_test_cohort.csv"))

# Save metrics evaluated at fixed 95% specificity (if applicable)
write_rds(metrics_at_95spec_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec_test_cohort.rds"))
write_csv(metrics_at_95spec_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec_test_cohort.csv"))

write_rds(metrics_at_95sens_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens_test_cohort.rds"))
write_csv(metrics_at_95sens_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens_updated_test_cohort.csv"))


## Tidy into one thing 
# 1. Tidy up and rename columns on each
youden_tbl <- metrics_youden_testing %>%
  rename(
    Model                 = model,
    Threshold_Youden      = threshold,
    Sensitivity_Youden    = sensitivity,
    Specificity_Youden    = specificity,
    PPV_Youden            = ppv,
    NPV_Youden            = npv,
    Accuracy_Youden       = accuracy,
    BalancedAcc_Youden    = bal_accuracy,
    F1_Youden             = f1
  )

spec95_tbl <- metrics_at_95spec_test %>%
  rename(
    Model                   = model,
    Threshold_95Spec        = threshold,
    Sensitivity_95Spec      = sensitivity,
    Specificity_95Spec      = specificity,
    Accuracy_95Spec         = accuracy
  )

sens95_tbl <- metrics_at_95sens_test %>%
  rename(
    Model                   = model,
    Threshold_95Sens        = threshold,
    Sensitivity_95Sens      = sensitivity,
    Specificity_95Sens      = specificity,
    Accuracy_95Sens         = accuracy,
    Precision_95Sens        = precision,
    F1_95Sens               = f1,
    BalancedAccuracy_95Sens = bal_accuracy
  )


# 2. Join them all together
all_perf_metrics <- youden_tbl %>%
  left_join(spec95_tbl, by = "Model") %>%
  left_join(sens95_tbl, by = "Model")

# 3. Export as a single supplementary Excel (or CSV)
library(writexl)    # or use write.csv if you prefer CSV

all_perf_metrics <- all_perf_metrics %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write_xlsx(
  list(`All Performance Metrics` = all_perf_metrics),
  path = "Final Tables and Figures/Supplementary_Table_5_All_Model_Metrics_Refit_Test_Cohort4.xlsx"
)

# — or, if prefer CSV:
write.csv(all_perf_metrics, 
          "Final Tables and Figures/Supplementary_Table_5_All_Model_Metrics_Refit_Test_Cohort3.csv", 
          row.names = FALSE)








### Now get whole cohort performance on rescored dataframe - skipped
# Define the specific model and a second fixed threshold
# model_name <- "BM_zscore_only_sites_prob"
# fixed_thr  <- 0.2211
# 
# # Pull out the Youden threshold for this model
# youden_thr <- threshold_df %>%
#   filter(model == model_name) %>%
#   pull(youden)
# 
# # Subset your full metrics table to just this model
# metrics_model <- all_metrics_rescored_primary %>%
#   filter(model == model_name)
# 
# # Find the row nearest the Youden cutoff
# metrics_youden <- metrics_model %>%
#   slice_min(order_by = abs(threshold - youden_thr), n = 1) %>%
#   mutate(type = "youden")
# 
# # Find the row nearest your fixed cutoff
# metrics_fixed <- metrics_model %>%
#   slice_min(order_by = abs(threshold - fixed_thr), n = 1) %>%
#   mutate(type = "fixed")
# 
# # Combine and print
# metrics_two <- bind_rows(metrics_youden, metrics_fixed) %>%
#   select(model, type, threshold, sensitivity, specificity, ppv, npv, accuracy, bal_accuracy, f1)
# 
# print(metrics_two)
# 
# write_rds(metrics_two, file = file.path(outdir, "cfWGS_model_metrics_fixed_z_score_model_deeper_LOD.rds"))
# write_csv(metrics_two, file = file.path(outdir, "cfWGS_model_metrics_fixed_z_score_model_deeper_LOD.csv"))
# write.csv(all_metrics_rescored_primary, file = "Final Tables and Figures/Suplementary_Table_3_All_BM_Zscore_Model_Performance_Cutoffs.csv")


# -----------------------------------------------------------------------------
# 10A. Performance Summaries & Exports - BM 
# -----------------------------------------------------------------------------
#### Plot to see what the model looks like
# For convenience
bm_obj      <- nested_bm_validation_updated2
models_list <- bm_obj$models            # named list of caret models
valid_df    <- train_bm                  # must contain MRD_truth + all predictors
bm_preds <- nested_bm_validation_updated2$outer_predictions



# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })

# 1) Compute one ROC curve per combo from your outer‑fold preds to get the nested CV version 
roc_dfs <- bm_preds %>% 
  group_by(combo) %>% 
  group_map(~{
    roc_obj <- roc(.x$truth, .x$prob,
                   levels    = c("neg","pos"),
                   direction = "<",
                   quiet     = TRUE)
    tibble(
      combo = .y$combo,
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities,
      auc   = as.numeric(auc(roc_obj))
    )
  }) %>% 
  bind_rows()


roc_df  <- bind_rows(roc_dfs)
auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 2) nice human-readable labels
pretty_combo_names <- c(
  BM_base                          = "All Mut Features",
  BM_base_zscore                   = "BM Sites + cVAF Z-score",
  BM_plus_fragment                 = "BM + Fragmentomics",
  BM_plus_fragment_min             = "BM + Fragments (min)",
  BM_rate_only                     = "BM cVAF",
  BM_zscore_only_detection_rate    = "BM cVAF Z-score",
  BM_zscore_only_sites             = "BM Sites Z-score"
)

## Get colors 
okabe_ito8 <- c(
  "#000000",  # black
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7"   # reddish purple
)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[levels(roc_df$combo)], 
    ", AUC = ", 
    formatC(auc_tbl$auc, digits = 2, format = "f")
  ),
  levels(roc_df$combo)
)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Pooled Outer‑Fold ROC Curve\n(nested 5×5 CV)"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(roc_df$combo))],
    labels = legend_labels
  )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- bm_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # force the same ordering of combos
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height = 0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width = 0.015) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x     = "Mean sensitivity (CV)",
    y     = "Mean specificity (CV)",
    title = "Fold‑Wise Sensitivity & Specificity\n(mean ± SD)"
  ) +
  # reuse the Okabe-Ito palette
  scale_colour_manual(
    values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 11)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_bm_updated.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

## Above is figure 4A and B 

### Now additional stats 

# ──────────────────────────────────────────────────────────────────────────────
# 1. Whole cohort summary: pick just the two BM z-score models, relabel, and rename cols
# supp figure 5A
cv_tbl <- metrics_youden %>%
  filter(model %in% c("BM_zscore_only_sites_prob", 
                      "BM_zscore_only_detection_rate_prob", 
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob          = "Sites model",
                   BM_zscore_only_detection_rate_prob = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
  #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens %>%
  filter(model %in% c("BM_zscore_only_sites_prob",
                      "BM_zscore_only_detection_rate_prob",
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob           = "Sites model",
                   BM_zscore_only_detection_rate_prob  = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
   # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Train Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)

  )
# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp5A_classifier_performance_bar_updated2.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)


 ## do for validation as well 
# ──────────────────────────────────────────────────────────────────────────────
cv_tbl <- metrics_youden_testing %>%
  filter(model %in% c("BM_zscore_only_sites_prob", 
                      "BM_zscore_only_detection_rate_prob", 
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob          = "Sites model",
                   BM_zscore_only_detection_rate_prob = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens_test %>%
  filter(model %in% c("BM_zscore_only_sites_prob",
                      "BM_zscore_only_detection_rate_prob",
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob           = "Sites model",
                   BM_zscore_only_detection_rate_prob  = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Test Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )

# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp5A_classifier_performance_bar_test_cohort_updated2.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)


### Add points on ROC where models selected - skipped
# # ────────────────────────────────────────────────────────────────
# # A.  Build a tibble of (model, threshold) pairs to mark
# # ────────────────────────────────────────────────────────────────
# mark_tbl <- tribble(
#   ~combo,                          ~threshold,
#   "BM_zscore_only_sites",           0.501,
#   "BM_zscore_only_sites",           0.221,
#   "BM_zscore_only_detection_rate",  0.436
# )
# 
# # ────────────────────────────────────────────────────────────────
# # B.  For each row in mark_tbl, compute FPR & TPR
# # ────────────────────────────────────────────────────────────────
# marker_pts <- pmap_dfr(mark_tbl, function(combo, threshold) {
#   
#   # retrieve fitted caret model & its predicted probs
#   prob <- predict(models_list[[combo]], newdata = valid_df, type = "prob")[ , positive_class]
#   
#   # truth vector (pos/neg as in earlier code)
#   truth <- valid_df$MRD_truth
#   
#   # binary calls at this threshold
#   pred_pos <- prob >= threshold
#   truth_pos <- truth == "pos"
#   
#   # confusion parts
#   tp <- sum(pred_pos & truth_pos)
#   fn <- sum(!pred_pos & truth_pos)
#   fp <- sum(pred_pos & !truth_pos)
#   tn <- sum(!pred_pos & !truth_pos)
#   
#   tibble(
#     combo      = combo,
#     threshold  = threshold,
#     fpr        = fp / (fp + tn),
#     tpr        = tp / (tp + fn)
#   )
# })
# 
# # ensure same factor order as curves
# marker_pts$combo <- factor(marker_pts$combo, levels = levels(roc_df$combo))
# 
# # ────────────────────────────────────────────────────────────────
# # C.  Add to roc_plot
# # ────────────────────────────────────────────────────────────────
# selected_models <- c("BM_zscore_only_sites", "BM_zscore_only_detection_rate")
# 
# roc_plot_tmp <- roc_df %>%
#   filter(combo %in% selected_models) %>%               # keep only those two
#   ggplot(aes(x = fpr, y = tpr, colour = combo)) +
#   geom_line(size = 1) +
#   geom_abline(lty = 2, colour = "grey60") +
#   labs(
#     x      = "False-positive rate (1 − specificity)",
#     y      = "True-positive rate (sensitivity)",
#     colour = NULL,
#     title  = "Cross-validated ROC curves (nested 5×5 folds)"
#   ) +
#   theme_bw(12) +
#   theme(
#     legend.position    = c(0.63, 0.15),
#     legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
#     legend.key.size    = unit(0.8, "lines"),
#     legend.text        = element_text(size = 10),
#     panel.grid         = element_blank()
#   ) +
#   scale_colour_manual(
#     # pick the same palette entries but only for your two models
#     values = okabe_ito8[ match(selected_models, levels(roc_df$combo)) ],
#     labels = legend_labels[selected_models]
#   )
# roc_plot_final <- roc_plot_tmp +
#   geom_point(
#     data    = marker_pts,
#     aes(x = fpr, y = tpr, shape = combo),
#     fill    = "white",
#     colour  = "black",
#     stroke  = 0.9,
#     size    = 3,
#     inherit.aes = FALSE
#   ) +
#   scale_shape_manual(
#     name   = "Marked thresholds",
#     values = c(
#       BM_zscore_only_sites           = 21,  # circle
#       BM_zscore_only_detection_rate  = 24   # triangle
#     ),
#     labels = c(
#       BM_zscore_only_sites          = "Sites model",
#       BM_zscore_only_detection_rate = "cVAF model"
#     )
#   ) +
#   guides(
#     # a) the shape legend for marked thresholds
#     shape = guide_legend(
#       title       = "Marked thresholds",
#       order       = 1,       # draw this first
#       ncol        = 1,
#       byrow       = FALSE,
#       keywidth    = unit(1.2, "lines"),
#       keyheight   = unit(1.2, "lines"),
#       override.aes = list(
#         fill   = "white",
#         colour = "black",
#         size   = 4,
#         stroke = 0.9
#       )
#     ),
#     # b) the colour legend for ROC curves
#     colour = guide_legend(
#       title       = "Model (AUC)",
#       order       = 2,       # draw this second
#       ncol        = 1,
#       byrow       = FALSE,
#       keywidth    = unit(2,   "lines"),
#       keyheight   = unit(0.5, "lines"),
#       override.aes = list(
#         size = 1.5           # slightly thicker lines in the key
#       )
#     )
#   ) +
#   theme(
#     # remove legend background box
#     legend.background = element_blank(),
#     legend.position       = c(0.72, 0.24),   # x=0.85 (near right), y=0.25 (higher)
#     # shrink the space between items
#     legend.spacing.y  = unit(0.2, "cm"),
#     legend.key        = element_blank(),
#     # style legend titles & text
#     legend.title      = element_text(face = "bold", size = 10),
#     legend.text       = element_text(size = 9), 
#     # title bold 
#     plot.title   = element_text(face = "bold", size = 13)
#   )
# 
# 
# ggsave(
#   file.path("Final Tables and Figures/Supp5A_ROC_performance.png"),
#   plot   = roc_plot_final,
#   width  = 5,
#   height = 4.25,
#   dpi    = 600
# )
# 
# 





## add contingency table 
# ──────────────────────────────────────────────────────────────────────────────
# 2) pull thresholds and models
bm_obj <- nested_bm_validation_updated2
thresh <- bm_obj$thresholds
mods   <- bm_obj$models[c("BM_zscore_only_sites",
                          "BM_base_zscore")]

# nice labels
model_labs <- c(
  "Sites model" = "BM_zscore_only_sites",
  "Combined model"  = "BM_base_zscore"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = train_bm, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, train_bm$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 5) +
  facet_wrap(~ model) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Count",
    begin  = 0.3,      # shift palette toward its lighter end
    end    = 0.9       # avoid the very darkest purples
  ) +  
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = c("pos", "neg")) +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Training Cohort"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig4C_confusion_tables_primary_updated3.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)


## Now validation
# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = hold_bm, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, hold_bm$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 5) +
  facet_wrap(~ model) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Count",
    begin  = 0.3,      # shift palette toward its lighter end
    end    = 0.9       # avoid the very darkest purples
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Test Cohort"
  ) +
  scale_y_discrete(limits = c("pos", "neg")) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig4C_confusion_tables_test3.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)




#### Now do the performance plot on validation cohort
# For convenience
valid_df    <- hold_bm                  # must contain MRD_truth + all predictors

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

roc_df  <- bind_rows(roc_dfs)
auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Test Samples ROC Curves"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(roc_df$combo))],
    labels = legend_labels
  )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- bm_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # force the same ordering of combos
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# Make a named vector of labels
labels_perf <- perf_df %>%
  transmute(
    combo,
    pretty = pretty_combo_names[combo],
    lbl = paste0(
      pretty,
      " (", percent(sens_mean, 1),
      " sens, ", percent(spec_mean, 1),
      " spec)"
    )
  ) %>%
  { setNames(.$lbl, .$combo) }
    
# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(sens_mean, spec_mean, colour = combo)) +
  geom_point(size=3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height=0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width=0.015) +
  geom_vline(xintercept=0.5, lty=2, colour="grey80") +
  geom_hline(yintercept=0.5, lty=2, colour="grey80") +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(perf_df$combo))],
    labels = labels_perf
  ) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    x     = "Mean sensitivity",
    y     = "Mean specificity",
    title = "Sensitivity vs. Specificity of cfWGS\nModels in the Test Cohort",
    colour = NULL
  ) +
  theme_bw(12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0,0),
    legend.background = element_rect(fill = alpha("white",0.7), colour=NA),
    plot.title      = element_text(hjust=0.5, face = "bold", size = 16)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_bm_validation_updated2.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/4E_performance_nested_folds_bm_validation_updated2.png",
  plot     = perf_plot,
  width    = 5,
  height   = 4,
  dpi      = 500
)






#### Now re-do for blood derived muts 
## add contingency table 
# ──────────────────────────────────────────────────────────────────────────────
# 2) pull thresholds and models
blood_obj      <- nested_blood_validation_updated2

thresh <- blood_obj$thresholds
mods   <- blood_obj$models[c("Blood_plus_fragment_min",
                          "Blood_zscore_only_sites")]

# nice labels
model_labs <- c(
  "Combined model" = "Blood_plus_fragment_min",
  "Sites model"  = "Blood_zscore_only_sites"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = train_blood, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, train_blood$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 5) +
  facet_wrap(~ model) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Count",
    begin  = 0.3,      # shift palette toward its lighter end
    end    = 0.9       # avoid the very darkest purples
  ) +  
  scale_x_discrete(position = "top") +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Training Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig5C_confusion_tables_primary_blood4.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)


## Now validation
# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = hold_blood, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, hold_blood$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 5) +
  facet_wrap(~ model) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Count",
    begin  = 0.3,      # shift palette toward its lighter end
    end    = 0.9       # avoid the very darkest purples
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Test Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig5C_confusion_tables_test_blood4.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)



### Now do the sensetivity vs specificity barplots for blood muts
cv_tbl <- metrics_youden %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_min_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_min_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_min_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_min_prob = "Combined model")  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Train Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )
# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp7A_classifier_performance_bar_updated_blood_muts.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)


## do for validation as well 
# ──────────────────────────────────────────────────────────────────────────────
cv_tbl <- metrics_youden_testing %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_min_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_min_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens_test %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_min_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_min_prob = "Combined model")
    ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc")) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Test Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )

# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp_7B_classifier_performance_bar_test_cohort_updated2.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)




### Make barplots for fragmentomic features
cv_tbl <- metrics_youden %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_prob", 
                      "Fragmentomics_prop_short_only_prob", 
                      "Fragmentomics_full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_prob = "Prop. short model",
                   Fragmentomics_full_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_prob", 
                      "Fragmentomics_prop_short_only_prob", 
                      "Fragmentomics_full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_prob = "Prop. short model",
                   Fragmentomics_full_prob = "Combined model")
  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Train Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )
# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp9D_classifier_performance_bar_updated_frag.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)


## do for validation as well 
# ──────────────────────────────────────────────────────────────────────────────
cv_tbl <- metrics_youden_testing %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_prob", 
                      "Fragmentomics_prop_short_only_prob", 
                      "Fragmentomics_full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_prob = "Prop. short model",
                   Fragmentomics_full_prob = "Combined model")
    ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens_test %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_prob", 
                      "Fragmentomics_prop_short_only_prob", 
                      "Fragmentomics_full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_prob = "Prop. short model",
                   Fragmentomics_full_prob = "Combined model")
    ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Test Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )

# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp_Fig9F_classifier_performance_bar_test_cohort_updated2_frag.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)




##### Now make the contingency tables
## add contingency table 
# ──────────────────────────────────────────────────────────────────────────────
# 2) pull thresholds and models
fragmentomics_obj      <- nested_fragmentomics_validation_updated3
models_list <- fragmentomics_obj$models            # named list of caret models
valid_df    <- train_fragmentomics                  # must contain MRD_truth + all predictors
fragmentomics_preds <- nested_fragmentomics_validation_updated3$outer_predictions

mods   <- fragmentomics_obj$models[c("Fragmentomics_mean_coverage_only",                             "Fragmentomics_min")]

# nice labels
model_labs <- c(
  "MM-sites coverage" = "Fragmentomics_mean_coverage_only",
  "Fragment score + MM-sites"  = "Fragmentomics_min"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = train_fragmentomics, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, train_fragmentomics$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 5) +
  facet_wrap(~ model) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Count",
    begin  = 0.3,      # shift palette toward its lighter end
    end    = 0.9       # avoid the very darkest purples
  ) +  
  scale_x_discrete(position = "top") +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Training Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Supp_Fig9E_confusion_tables_primary_fragmentomics.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)


## Now validation
# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = hold_fragmentomics, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, hold_fragmentomics$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 5) +
  facet_wrap(~ model) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Count",
    begin  = 0.3,      # shift palette toward its lighter end
    end    = 0.9       # avoid the very darkest purples
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Test Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Supp_Fig9F_confusion_tables_val_fragmentomics.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)









## DO ROC curve on primary blood muts
# -----------------------------------------------------------------------------
# 10B. Performance Summaries & Exports - Blood 
# -----------------------------------------------------------------------------
#### Plot to see what the model looks like
# For convenience
models_list <- blood_obj$models            # named list of caret models
valid_df    <- train_blood                  # must contain MRD_truth + all predictors
blood_preds <- blood_obj$outer_predictions

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })

roc_dfs <- blood_preds %>% 
  group_by(combo) %>% 
  group_map(~{
    roc_obj <- roc(.x$truth, .x$prob,
                   levels    = c("neg","pos"),
                   direction = "<",
                   quiet     = TRUE)
    tibble(
      combo = .y$combo,
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities,
      auc   = as.numeric(auc(roc_obj))
    )
  }) %>% 
  bind_rows()

# Rebuild roc_df without any Fragmentomics_ combos as re-trained later
roc_df <- bind_rows(roc_dfs) %>%
  filter(!grepl("^Fragmentomics_", combo))

# Then recompute AUC table
auc_tbl <- roc_df %>%
  distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 2) nice human-readable labels
pretty_combo_names <- c(
  Blood_base                          = "All Mut Features",
  Blood_base_zscore                   = "Blood Sites + cVAF Z-score",
  Blood_plus_fragment                 = "Blood + Fragmentomics",
  Blood_plus_fragment_min             = "Blood + Fragments (min)",
  Blood_rate_only                     = "Blood cVAF",
  Blood_zscore_only_detection_rate    = "Blood cVAF Z-score",
  Blood_zscore_only_sites             = "Blood Sites Z-score"
)

## Get colors 
okabe_ito8 <- c(
  "#000000",  # black
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7"   # reddish purple
)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[levels(roc_df$combo)], 
    ", AUC = ", 
    formatC(auc_tbl$auc, digits = 2, format = "f")
  ),
  levels(roc_df$combo)
)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Pooled Outer‑Fold ROC Curve\n(nested 5×5 CV)"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  ) +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(roc_df$combo))],
    labels = legend_labels
  )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- blood_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # keep only combos that actually made it into your roc_df
  filter(combo %in% levels(roc_df$combo)) %>%
  # now drop any leftover unused levels
  mutate(combo = factor(combo, levels = levels(roc_df$combo))) 

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height = 0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width = 0.015) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x     = "Mean sensitivity (CV)",
    y     = "Mean specificity (CV)",
    title = "Fold‑Wise Sensitivity & Specificity\n(mean ± SD)"
  ) +
  # reuse the Okabe-Ito palette
  scale_colour_manual(
    values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 11)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Fig_5A_updated_combined_ROC_and_performance_nested_folds_blood_features_updated.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)




### Add the metrics from dilution series
### Add points on ROC where models selected 
analysis_df <- data_scored_masked %>% filter(Sample_Code %in% train_blood$Sample_Code) %>%
  filter(!timepoint_info %in% c("Baseline", "Diagnosis"))

# Identify probability columns present in analysis_df
prob_cols <- grep("_prob$", names(analysis_df), value = TRUE)

# Build long-form predictions: one row per sample × model
# Map "Blood_zscore_only_sites_prob" -> combo "Blood_zscore_only_sites", etc.
pred_long <- analysis_df %>%
  select(MRD_truth, all_of(prob_cols)) %>%
  pivot_longer(cols = all_of(prob_cols),
               names_to = "model_prob",
               values_to = "prob") %>%
  mutate(
    combo = sub("_prob$", "", model_prob),
    # robust truth: accept 1/"pos" as positive
    truth_pos = MRD_truth %in% c(1, "pos"),
    truth_neg = MRD_truth %in% c(0, "neg")
  ) %>%
  select(combo, prob, truth_pos, truth_neg)

# Do ROC df on whole cohort refit
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

# Rebuild roc_df without any Fragmentomics_ combos as re-trained later
roc_df <- bind_rows(roc_dfs) %>%
  filter(!grepl("^Fragmentomics_", combo))

roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# Keep factor levels consistent with your ROC plot (if roc_df exists already)
if (exists("roc_df")) {
  pred_long$combo <- factor(pred_long$combo, levels = levels(roc_df$combo))
}

# ------------------------------------------------------------------
# Thresholds to mark (same names as in your ROC 'combo' column)
# ------------------------------------------------------------------
mark_tbl <- tribble(
  ~combo,                      ~threshold,
  "Blood_zscore_only_sites",      0.432,
  "Blood_zscore_only_sites",      0.380,
  "Blood_plus_fragment_min",      0.435
)

# ------------------------------------------------------------------
# Compute sensitivity/specificity at each (combo, threshold)
# ------------------------------------------------------------------
marker_pts <- pmap_dfr(mark_tbl, function(combo, threshold) {
  df <- pred_long %>% filter(combo == !!combo)
  
  # Confusion counts using >= threshold rule
  pred_pos <- df$prob >= threshold
  tp <- sum(pred_pos & df$truth_pos, na.rm = TRUE)
  fn <- sum(!pred_pos & df$truth_pos, na.rm = TRUE)
  fp <- sum(pred_pos & df$truth_neg, na.rm = TRUE)
  tn <- sum(!pred_pos & df$truth_neg, na.rm = TRUE)
  
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  
  tibble(
    combo       = combo,
    threshold   = threshold,
    sensitivity = sens,
    specificity = spec,
    tpr         = sens,
    fpr         = 1 - spec
  )
})

# Align factor levels for plotting aesthetics
if (exists("roc_df")) {
  marker_pts$combo <- factor(marker_pts$combo, levels = levels(roc_df$combo))
}

# Optional: round for reporting
marker_pts <- marker_pts %>%
  mutate(across(c(sensitivity, specificity, tpr, fpr), ~ round(.x, 3)))

print(marker_pts)

selected_models <- c("Blood_zscore_only_sites", "Blood_plus_fragment_min")

roc_plot_tmp <- roc_df %>%
  filter(combo %in% selected_models) %>%               # keep only those two
  ggplot(aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,
    title  = "Full-Cohort ROC (refit on all training samples)"
  ) +
  theme_bw(12) +
  theme(
    legend.position    = c(0.63, 0.15),
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid         = element_blank()
  ) +
  scale_colour_manual(
    # pick the same palette entries but only for your two models
    values = okabe_ito8[ match(selected_models, levels(roc_df$combo)) ],
    labels = c(
      Blood_zscore_only_sites          = "Sites model",
      Blood_plus_fragment_min = "Combined model"
    )
  )
roc_plot_final <- roc_plot_tmp +
  geom_point(
    data    = marker_pts,
    aes(x = fpr, y = tpr),
    fill    = "white",
    colour  = "darkgrey",
    stroke  = 0.9,
    size    = 2,
    inherit.aes = FALSE
  ) +
  guides(
    # a) the shape legend for marked thresholds
    shape = guide_legend(
      title       = "Marked thresholds",
      order       = 1,       # draw this first
      ncol        = 1,
      byrow       = FALSE,
      keywidth    = unit(1.2, "lines"),
      keyheight   = unit(1.2, "lines"),
      override.aes = list(
        fill   = "white",
        colour = "black",
        size   = 4,
        stroke = 0.9
      )
    ),
    # b) the colour legend for ROC curves
    colour = guide_legend(
      title       = "Model",
      order       = 2,       # draw this second
      ncol        = 1,
      byrow       = FALSE,
      keywidth    = unit(2,   "lines"),
      keyheight   = unit(0.5, "lines"),
      override.aes = list(
        size = 1.5           # slightly thicker lines in the key
      )
    )
  ) +
  theme(
    # remove legend background box
    legend.background = element_blank(),
    #  legend.position       = c(0.72, 0.24),   # x=0.85 (near right), y=0.25 (higher)
    legend.position       = "right",
    # shrink the space between items
    #    legend.spacing.y  = unit(0.2, "cm"),
    legend.key        = element_blank(),
    # style legend titles & text
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9), 
    # title bold 
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.title.position = "plot"
  )


ggsave(
  file.path("Final Tables and Figures/Supp7D_ROC_performance_blood_updated3.png"),
  plot   = roc_plot_final,
  width  = 5.5,
  height = 4,
  dpi    = 600
)


### Do ROC curve on validation for blood derived muts
#### Now do the performance plot on validation cohort
# For convenience
valid_df    <- hold_blood                # must contain MRD_truth + all predictors

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

roc_df  <- bind_rows(roc_dfs)

# Rebuild roc_df without any Fragmentomics_ combos as re-trained later
roc_df <- bind_rows(roc_dfs) %>%
  filter(!grepl("^Fragmentomics_", combo))

auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Test Samples ROC Curves"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(roc_df$combo))],
    labels = legend_labels
  )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- blood_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # keep only combos that actually made it into your roc_df
  filter(combo %in% levels(roc_df$combo)) %>%
  # now drop any leftover unused levels
  mutate(combo = factor(combo, levels = levels(roc_df$combo))) 

# Make a named vector of labels
labels_perf <- perf_df %>%
  transmute(
    combo,
    pretty = pretty_combo_names[combo],
    lbl = paste0(
      pretty,
      " (", percent(sens_mean, 1),
      " sens, ", percent(spec_mean, 1),
      " spec)"
    )
  ) %>%
  { setNames(.$lbl, .$combo) }

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(sens_mean, spec_mean, colour = combo)) +
  geom_point(size=3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height=0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width=0.015) +
  geom_vline(xintercept=0.5, lty=2, colour="grey80") +
  geom_hline(yintercept=0.5, lty=2, colour="grey80") +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(perf_df$combo))],
    labels = labels_perf
  ) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    x     = "Mean sensitivity",
    y     = "Mean specificity",
    title = "Sensitivity vs. Specificity of cfWGS\nModels in the Test Cohort",
    colour = NULL
  ) +
  theme_bw(12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0,0),
    legend.background = element_rect(fill = alpha("white",0.7), colour=NA),
    plot.title      = element_text(hjust=0.5, face = "bold", size = 16)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))


# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Supp_Fig_7C_combined_ROC_and_performance_nested_folds_bm_validation.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/Supp_Fig_7C_performance_nested_folds_blood_validation_updated.png",
  plot     = perf_plot,
  width    = 5,
  height   = 4,
  dpi      = 500
)


















#### Now do ROC curve on validation cohort for BM muts - old 
# For convenience
valid_df    <- hold_bm                  # must contain MRD_truth + all predictors

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

roc_df  <- bind_rows(roc_dfs)
auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Test Samples ROC Curves"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(roc_df$combo))],
    labels = legend_labels
  )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- bm_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # force the same ordering of combos
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# Make a named vector of labels
labels_perf <- perf_df %>%
  transmute(
    combo,
    pretty = pretty_combo_names[combo],
    lbl = paste0(
      pretty,
      " (", percent(sens_mean, 1),
      " sens, ", percent(spec_mean, 1),
      " spec)"
    )
  ) %>%
  { setNames(.$lbl, .$combo) }

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(sens_mean, spec_mean, colour = combo)) +
  geom_point(size=3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height=0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width=0.015) +
  geom_vline(xintercept=0.5, lty=2, colour="grey80") +
  geom_hline(yintercept=0.5, lty=2, colour="grey80") +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(perf_df$combo))],
    labels = labels_perf
  ) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    x     = "Mean sensitivity",
    y     = "Mean specificity",
    title = "Sensitivity vs. specificity of cfWGS models\nin the test cohort",
    colour = NULL
  ) +
  theme_bw(12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = c(0.05, 0.05),
    legend.justification = c(0,0),
    legend.background = element_rect(fill = alpha("white",0.7), colour=NA),
    plot.title      = element_text(hjust=0.5, face = "bold", size = 14)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_bm_validation.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/4E_performance_nested_folds_bm_validation_updated.png",
  plot     = perf_plot,
  width    = 6,
  height   = 4,
  dpi      = 500
)














## Add the other metrics
#### Now do on validation for blood
# # For convenience
# valid_df    <- hold_blood                  # must contain MRD_truth + all predictors
# 
# # ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })
# 
# # Rebuild roc_df without any Fragmentomics_ combos as re-trained later
# roc_df <- bind_rows(roc_dfs) %>%
#   filter(!grepl("^Fragmentomics_", combo))
# 
# auc_tbl <- roc_df %>% distinct(combo, auc)
# 
# # 1) reorder combos by AUC
# auc_tbl <- auc_tbl %>% arrange(desc(auc))
# roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)
# 
# # 4) make the plot, using a qualitative brewer palette
# roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
#   geom_line(size = 1) +
#   geom_abline(lty = 2, colour = "grey60") +
#   labs(
#     x      = "False-positive rate (1 − specificity)",
#     y      = "True-positive rate (sensitivity)",
#     colour = NULL,              # no title above the legend
#     title  = "Hold-out samples ROC curves"
#   ) +
#   theme_bw(14) +
#   theme(
#     legend.position    = c(0.63, 0.15),  # inside, bottom-right
#     legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
#     legend.key.size    = unit(0.8, "lines"),
#     legend.text        = element_text(size = 10),
#     panel.grid = element_blank()
#   ) +
#   scale_colour_manual(
#     values = okabe_ito8[1:length(levels(roc_df$combo))],
#     labels = legend_labels
#   )
# 
# 
# # ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
# perf_df <- bm_obj$nested_metrics %>%
#   select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
#   # force the same ordering of combos
#   mutate(combo = factor(combo, levels = levels(roc_df$combo)))
# 
# # 1) Build performance plot
# perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(
#     xmin = pmax(0, sens_mean - sens_sd),
#     xmax = pmin(1, sens_mean + sens_sd)
#   ), height = 0.015) +
#   geom_errorbar(aes(
#     ymin = pmax(0, spec_mean - spec_sd),
#     ymax = pmin(1, spec_mean + spec_sd)
#   ), width = 0.015) +
#   geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
#   geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
#   scale_x_continuous(limits = c(0,1)) +
#   scale_y_continuous(limits = c(0,1)) +
#   labs(
#     x     = "Mean sensitivity",
#     y     = "Mean specificity",
#     title = "Hold-out samples performance"
#   ) +
#   # reuse the Okabe-Ito palette
#   scale_colour_manual(
#     values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
#   ) +
#   theme_bw(base_size = 14) +
#   theme(
#     panel.grid      = element_blank(),
#     legend.position = "none",
#     plot.title      = element_text(hjust = 0.5),
#     axis.title      = element_text(size = 13),
#     axis.text       = element_text(size = 11)
#   )
# 
# # ── 3) Combine with roc_plot ────────────────────────────────────────────────
# combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))
# 
# # ── 4) Export ───────────────────────────────────────────────────────────────
# ggsave(
#   filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_blood_features_validation.png",
#   plot     = combined_plot,
#   width    = 12,
#   height   = 6,
#   dpi      = 500
# )






# -----------------------------------------------------------------------------
# 10C. Performance Summaries & Exports - Fragmentomics
# -----------------------------------------------------------------------------
#### Plot to see what the model looks like
# For convenience
fragmentomics_obj      <- nested_fragmentomics_validation_updated3
models_list <- fragmentomics_obj$models            # named list of caret models
valid_df    <- train_fragmentomics                  # must contain MRD_truth + all predictors
fragmentomics_preds <- nested_fragmentomics_validation_updated3$outer_predictions

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })

## On the pooled outer fold predictions 
# 1) Compute one ROC curve per combo from your outer‑fold preds
# 1) Compute one ROC curve per combo from your outer‑fold preds to get the nested CV version 
roc_dfs <- fragmentomics_preds %>% 
  group_by(combo) %>% 
  group_map(~{
    roc_obj <- roc(.x$truth, .x$prob,
                   levels    = c("neg","pos"),
                   direction = "<",
                   quiet     = TRUE)
    tibble(
      combo = .y$combo,
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities,
      auc   = as.numeric(auc(roc_obj))
    )
  }) %>% 
  bind_rows()


roc_df  <- bind_rows(roc_dfs)
auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 2) nice human-readable labels
pretty_combo_names <- c(
  Fragmentomics_full                  = "FS + MeanCov + TF + PropShort",
  Fragmentomics_min                   = "FS + MeanCov",
  Fragmentomics_FS_only               = "Fragment Size Only",
  Fragmentomics_mean_coverage_only    = "Mean Coverage Only",
  Fragmentomics_prop_short_only       = "Prop. Short Fragments Only",
  Fragmentomics_tumor_fraction_only   = "Tumor Fraction Only"
)

## Get colors 
okabe_ito8 <- c(
  "#000000",  # black
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7"   # reddish purple
)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[levels(roc_df$combo)], 
    ", AUC = ", 
    formatC(auc_tbl$auc, digits = 2, format = "f")
  ),
  levels(roc_df$combo)
)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Pooled Outer‑Fold ROC Curve\n(nested 5×5 CV)"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_colour_manual(
    values = okabe_ito8[1:length(levels(roc_df$combo))],
    labels = legend_labels
  )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- fragmentomics_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # force the same ordering of combos
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height = 0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width = 0.015) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x     = "Mean sensitivity (CV)",
    y     = "Mean specificity (CV)",
    title = "Fold‑Wise Sensitivity & Specificity\n(mean ± SD)"
  ) +
  # reuse the Okabe-Ito palette
  scale_colour_manual(
    values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 11)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_fragmentomics_updated3.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)






# -----------------------------------------------------------------------------
# 11. Additional Diagnostic Plots
# -----------------------------------------------------------------------------

fig_dir <- file.path(outdir, "Diagnostic_Figures_2025")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Helper: pretty label for each combo
pretty_combo_names <- c(
  BM_zscore_only_sites            = "BM Sites Z-score",
  BM_zscore_only_detection_rate   = "BM cVAF Z-score",
  Blood_rate_only                 = "Blood cVAF",
  Blood_plus_fragment_min         = "Blood + Fragments (min)",
  Blood_zscore_only_detection_rate= "Blood cVAF Z-score",
  Fragmentomics_min               = "Fragmentomics (FS + Cov.)",
  Fragmentomics_mean_coverage_only= "Fragmentomics Coverage"
)

# choose the same frontline MRD set you trained on
valid_df <- train_df          # contains MRD_truth + ALL predictors
positive_class <- "pos"

# loop over every wanted combo
for (cmb in wanted) {
  
  fit <- selected_models[[cmb]]
  thr <- selected_thr[[cmb]]
  
  # figure out which predictors this model requires
  preds <- setdiff(colnames(fit$trainingData), ".outcome")
  
  # restrict to rows that have all predictors
  keep <- complete.cases(valid_df[, preds])
  if (!any(keep)) {
    message("Skipping ", cmb, " (no rows with complete predictors).")
    next
  }
  
  probs <- predict(fit,
                   newdata = valid_df[keep, preds, drop = FALSE],
                   type = "prob")[[positive_class]]
  
  truth_fac <- valid_df$MRD_truth[keep]          # factor "neg"/"pos"
  truth     <- ifelse(truth_fac == "pos", 1, 0)  # numeric 0/1
  
  # ── 11.1 Calibration plot (deciles) ──────────────────────────────────────
  cal_df <- tibble(truth = truth, pred = probs) %>%
    mutate(dec = ntile(pred, 10)) %>%
    group_by(dec) %>%
    summarise(mean_pred = mean(pred), obs_rate = mean(truth), .groups = "drop")
  
  cal_plot <- ggplot(cal_df, aes(mean_pred, obs_rate)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_abline(lty = 2, col = "grey60") +
    scale_x_continuous(labels = percent_format(1)) +
    scale_y_continuous(labels = percent_format(1)) +
    labs(
      x = "Mean predicted probability",
      y = "Observed MRD-positive rate",
      title = paste0("Calibration: ", pretty_combo_names[cmb])
    ) +
    theme_bw(14) + theme(panel.grid = element_blank())
  
  ggsave(
    filename = file.path(fig_dir,
                         paste0("calibration_", cmb, ".png")),
    plot = cal_plot, width = 6, height = 5, dpi = 300
  )
  
  # ── 11.2 Precision–Recall curve ─────────────────────────────────────────
  pr <- pr.curve(scores.class0 = probs[truth == 1],
                 scores.class1 = probs[truth == 0],
                 curve = TRUE)
  pr_df <- data.frame(recall = pr$curve[, 1],
                      precision = pr$curve[, 2])
  
  pr_plot <- ggplot(pr_df, aes(recall, precision)) +
    geom_line(size = 1, col = "#E69F00") +                 # orange
    labs(
      x = "Recall (sensitivity)",
      y = "Precision (PPV)",
      title = paste0("PR curve: ", pretty_combo_names[cmb]),
      subtitle = sprintf("AUPRC = %.2f", pr$auc.integral)
    ) +
    theme_bw(14) + theme(panel.grid = element_blank())
  
  ggsave(
    filename = file.path(fig_dir,
                         paste0("pr_curve_", cmb, ".png")),
    plot = pr_plot, width = 6, height = 5, dpi = 300
  )
  
  # ── 11.3 Sensitivity & specificity vs cutoff ────────────────────────────
  roc_obj    <- roc(response = truth_fac, predictor = probs,
                    quiet = TRUE, levels = c("neg", "pos"))
  thresholds <- seq(0, 1, by = 0.01)
  
  # coords() returns a matrix with columns "sensitivity" & "specificity"
  metric_df  <- as_tibble(
    coords(
      roc_obj,
      x         = thresholds,
      input     = "threshold",
      ret       = c("sensitivity", "specificity"),
      transpose = FALSE
    )
  )
  
  # bind the numeric threshold vector back on
  coords_df <- tibble(threshold = thresholds) %>%
    bind_cols(metric_df)
  
  thr_plot <- ggplot(coords_df, aes(x = threshold)) +
    geom_line(aes(y = sensitivity,   col = "Sensitivity"),   size = 1) +
    geom_line(aes(y = specificity,   col = "Specificity"),   size = 1) +
    scale_color_manual(
      NULL,
      values = c("Sensitivity" = "#009E73", "Specificity" = "#0072B2")
    ) +
    scale_x_continuous(labels = percent_format(0.1)) +
    scale_y_continuous(labels = percent_format(1)) +
    labs(
      x     = "Probability cutoff",
      y     = "Value",
      title = paste0("Sens/Spec vs cutoff: ", pretty_combo_names[cmb])
    ) +
    theme_bw(14) +
    theme(panel.grid = element_blank())
  
  ggsave(
    filename = file.path(fig_dir,
                         paste0("sens_spec_vs_cutoff_", cmb, ".png")),
    plot   = thr_plot,
    width  = 6,
    height = 5,
    dpi    = 300
  )
  
  message("Exported diagnostics for ", cmb)
}














#### Below here is testing - not used in main analysis 

### Now try smoothed 
# ── 1. Compute binormal-smoothed ROC curves ────────────────────────────────
smoothed_rocs <- imap_dfr(models_list, function(fit, label) {
  # 1a) get predicted probabilities on your primary cohort
  prob <- predict(fit, newdata = valid_df, type = "prob")[, positive_class]
  
  # 1b) raw ROC object
  roc_obj <- roc(
    response  = valid_df$MRD_truth,
    predictor = prob,
    levels    = c("neg","pos"),
    direction = "<",
    quiet     = TRUE
  )
  
  # 1c) binormal smoothing
  roc_s <- smooth(roc_obj, method = "binormal")
  
  # 1d) pack into a tibble
  tibble(
    combo = label,
    fpr   = 1 - roc_s$specificities,
    tpr   =     roc_s$sensitivities,
    auc   = as.numeric(auc(roc_obj))
  )
})

# ── 2. Extract one AUC per combo for legend labels ─────────────────────────
auc_tbl <- smoothed_rocs %>%
  distinct(combo, auc) %>%
  arrange(combo)

# ── 3. Plot ────────────────────────────────────────────────────────────────
ggplot(smoothed_rocs, aes(x = fpr, y = tpr, color = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, color = "grey60") +
  scale_color_manual(
    values = scales::hue_pal()(nrow(auc_tbl)),
    labels = setNames(
      sprintf("%s | AUC = %.2f", auc_tbl$combo, auc_tbl$auc),
      auc_tbl$combo
    )
  ) +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = "Combo",
    title  = "ROC curves on primary cohort\n(binormal smoothing)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )





### Get more info on total number of samples
### Now check sample counts 

# Define export directory
export_dir <- "Output_tables_2025"

## Edit to ones we have cfDNA analyzed on in primary cohorts
dat_valid <- dat %>%
  filter(Patient %in% cohort_df$Patient) %>% filter(!is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA)) %>% filter(!is.na(FS))

## 2. Per‑patient summaries
patient_summ <- dat_valid %>%
  group_by(Patient) %>%
  summarise(
    n_samples   = dplyr::n(),                                   # how many cfDNA samples
    first_date  = min(Date, na.rm = TRUE),
    last_date   = max(Date, na.rm = TRUE),
    follow_days = as.numeric(difftime(last_date, first_date, units = "days")),
    .groups = "drop"
  )

## 4. Study‑level numbers
n_patients        <- n_distinct(dat_valid$Patient)        # should be 51 (check)
total_samples     <- nrow(dat_valid)
median_per_pt     <- median(patient_summ$n_samples)
range_per_pt      <- range(patient_summ$n_samples)
median_follow     <- median(patient_summ$follow_days)
range_follow      <- range(patient_summ$follow_days)

## 5. Auto‑generate the sentence
glue(
  "This study comprises cfWGS analysis of {total_samples} blood cfDNA samples ",
  "from {n_patients} MM patients (median {median_per_pt} samples per patient, ",
  "range {range_per_pt[1]}–{range_per_pt[2]}), collected over a median follow‑up ",
  "of {median_follow} days (range {range_follow[1]}–{range_follow[2]})."
)


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
  distinct(Sample_Code, Cohort) %>%
  count(Cohort, name = "n_samples_with_blood_only")

