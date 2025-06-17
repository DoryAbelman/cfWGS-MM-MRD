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
library(Matrix)
library(caret)



### Set paths 
outdir   <- "Output_tables_2025"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

### Load data 
file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated4.rds")
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







### Next do a full nested cross validation 5x5 with hyperparameter tuning 
### Even more robust with less data leakage
positive_class <- "pos"

# 1) Identify which predictors belong to BM vs Blood combos
# BM combos: anything starting “BM_” but drop “BM_blood_ultimate”
bm_names <- names(combos)[
  startsWith(names(combos), "BM_") &
    names(combos) != "BM_blood_ultimate"
]
combos_bm    <- combos[bm_names]

# Blood combos: anything starting “Blood_” OR containing “Just_fragmentomics”
blood_names <- names(combos)[
  startsWith(names(combos), "Blood_") |
    grepl("Just_fragmentomics", names(combos))
]
combos_blood <- combos[blood_names]

# 2) Build train_bm / train_blood by dropping any row missing *any* BM or Blood predictors
bm_preds    <- unique(unlist(combos_bm))
blood_preds <- unique(unlist(combos_blood))

train_bm    <- train_df %>% drop_na(all_of(c("MRD_truth", bm_preds)))
train_blood <- train_df %>% drop_na(all_of(c("MRD_truth", blood_preds)))


train_bm    <- train_df %>% drop_na(all_of(c("MRD_truth", bm_preds)))
train_blood <- train_df %>% drop_na(all_of(c("MRD_truth", blood_preds)))

hold_bm    <- hold_df %>% drop_na(all_of(c("MRD_truth", bm_preds)))
hold_blood <- hold_df %>% drop_na(all_of(c("MRD_truth", blood_preds)))

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
  
  # 2) Loop over each feature-combo
  results <- purrr::imap(combo_list, function(preds, label) {
    
    # 2a) Collect per‐fold metrics
    combo_metrics <- purrr::map_dfr(seq_along(outer_folds), function(i) {
      # split out vs. in
      train_out <- train_data[ outer_folds[[i]], ]
      test_out  <- train_data[-outer_folds[[i]], ]
      
      # complete‐case
      train_cc <- train_out %>% drop_na(MRD_truth, all_of(preds))
      test_cc  <- test_out  %>% drop_na(MRD_truth, all_of(preds))
      
      # SKIP if either side is too small or single‐class
      if (nrow(train_cc) < 20 ||
          n_distinct(train_cc$MRD_truth) < 2 ||
          n_distinct(test_cc$MRD_truth)  < 2) {
        return(NULL)
      }
      
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
      roc_o    <- pROC::roc(test_cc$MRD_truth, probs,
                            quiet  = TRUE,
                            levels = c("neg","pos"))
      youden_t <- pROC::coords(roc_o, "best",
                               best.method = "youden",
                               ret         = "threshold")[[1]]
      
      # pull out single‐value sens/spec
      sens_o <- pROC::coords(roc_o, x = youden_t,
                             ret       = "sensitivity",
                             transpose = FALSE)[[1]]
      spec_o <- pROC::coords(roc_o, x = youden_t,
                             ret       = "specificity",
                             transpose = FALSE)[[1]]
      
      tibble(
        combo       = label,
        auc         = as.numeric( pROC::auc(roc_o) ),
        sensitivity = sens_o,
        specificity = spec_o,
        accuracy    = mean((probs >= youden_t) ==
                             (test_cc$MRD_truth == positive_class))
      )
    })
    
    # 2d) summarize across *kept* folds
    nested_summary <- combo_metrics %>%
      summarise(
        combo     = first(combo),
        auc_mean  = mean(auc,           na.rm = TRUE),
        sens_mean = mean(sensitivity,   na.rm = TRUE),
        spec_mean = mean(specificity,   na.rm = TRUE),
        acc_mean  = mean(accuracy,      na.rm = TRUE),
        .groups   = "drop"
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
    youden_full <- pROC::coords(roc_full, "best",
                                best.method = "youden",
                                ret         = "threshold")[[1]]
    
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
    
    validation_summary <- tibble(
      combo       = label,
      auc_valid   = as.numeric( pROC::auc(roc_val) ),
      sens_valid  = sens_val,
      spec_valid  = spec_val,
      acc_valid   = acc_val
    )
    
    # return all pieces
    list(
      nested     = nested_summary,
      final_fit  = final_fit,
      bestTune   = final_fit$bestTune,
      validation = validation_summary
    )
  })
  
  # pull them back out
  nested_metrics     <- map_dfr(results, "nested")
  models             <- map(results, "final_fit")
  best_tunes         <- map(results, "bestTune")
  validation_metrics <- map_dfr(results, "validation")
  
  list(
    nested_metrics     = nested_metrics,
    models             = models,
    best_tunes         = best_tunes,
    validation_metrics = validation_metrics
  )
}


# 4) Run nested CV separately
nested_bm_validation <- run_nested_with_validation(
  train_data = train_bm,
  valid_data = hold_bm,
  combo_list = combos_bm,
  positive_class = "pos"
)

nested_blood_validation <- run_nested_with_validation(
  train_data = train_blood,
  valid_data = hold_blood,
  combo_list = combos_blood,
  positive_class = "pos"
)


# 5) Export
saveRDS(nested_bm_validation, file = "nested_bm_validation.rds")
saveRDS(nested_blood_validation, file = "nested_blood_validation.rds")


#### Other things to add: Bootstrapping for confidence intervals and permutation testing to show AUC drop  



# ──────────────────────────────────────────────────────────────────────────────
# 4.  Save for manuscript ------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
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




# ────────────────────────────────────────────────────────────────────────────
# 7.  Quick checks -----------------------------------------------------------
message("BM_zscore_only positives: ",
        sum(dat$BM_zscore_only_call == 1, na.rm = TRUE))
message("Blood_base_call positives: ",
        sum(dat$Blood_base_call == 1, na.rm = TRUE))
message("Blood_rate_all_extras_call positives: ",
        sum(dat$Blood_rate_all_extras_call == 1, na.rm = TRUE))


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





#### Get more info on total number of samples
### Now check sample counts 

# Define export directory
export_dir <- "Output_tables_2025"

# Load from CSV (if you want to view/edit easily)
patient_cohort_tbl_csv <- read.csv(file.path(export_dir, "patient_cohort_assignment.csv"))

## Edit to ones we have cfDNA analyzed on in primary cohorts
dat_valid <- dat %>%
  filter(Patient %in% patient_cohort_tbl_csv$Patient) %>% filter(!is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA))

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










#### Below here is testing - not used in main analysis 
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
