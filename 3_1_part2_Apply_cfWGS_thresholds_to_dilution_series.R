# ---------------------------------------------------------------------------
# Script:   score_dilution_series_cfWGS.R
# Purpose:  Apply the previously optimised cfWGS MRD classifier to the dilution series
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(Matrix)
library(readr)     # read_rds / write_csv
library(glmnet)
library(pROC)
library(stringr)


# ── 1. FILE PATHS ───────────────────────────────────────────────────────────
PATH_MODEL_LIST       <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_models.rds"
PATH_THRESHOLD_LIST   <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_thresholds.rds"
PATH_DILUTION_FRAGMENTOMICS         <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Results_Fragmentomics/Dilution_series/key_fragmentomics_info_dilution_series.rds"
PATH_DILUTION_PROCESSED_MRDetect    <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/MRDetect_output_winter_2025/Processed_R_outputs/cfWGS_Winter2025Dilution_series_May2025_with_zscore.rds"
PATH_DILUTION_CLINICAL  <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series/Metadata_dilution_series.csv"
PATH_TUMOR_FRACTION <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series/tumor_fraction_dilution_series.txt" # To add

OUTPUT_DIR            <- "Dilution_Series_Scoring_2025"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ── 2. LOAD SAVED MODELS & THRESHOLDS ───────────────────────────────────────
model_list    <- readRDS(PATH_MODEL_LIST)
selected_rows <- readRDS(PATH_THRESHOLD_LIST)

# ── 3. REDEFINE YOUR COMBOS ─────────────────────────────────────────────────
combos <- list(
  BM_zscore_only   = c("zscore_BM"),
  Blood_all_extras = c("zscore_blood","detect_rate_blood","FS",
                       "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_base       = c("zscore_blood","detect_rate_blood")
)

# ── 4. LOAD AND ASSEMBLE DILUTION SERIES DATA ────────────────────────────────
frag_df     <- read_rds(PATH_DILUTION_FRAGMENTOMICS)
mrdetect_df <- read_rds(PATH_DILUTION_PROCESSED_MRDetect)
clinical_df <- read_csv(PATH_DILUTION_CLINICAL)
tumor_fraction_dilution <- read_tsv(PATH_TUMOR_FRACTION)

# A) Combine features into one table 

# Merge based on BAM (make sure both are character columns)
mrdetect_df <- mrdetect_df %>%
  left_join(
    clinical_df %>% dplyr::select(BAM, Merge),
    by = "BAM"
  )

# ── Filter MRDetect df to only timepoint 0 (diagnosis) and select/rename ─────────
mrdetect_wide <- mrdetect_df %>%
  filter(Timepoint %in% c("01")) %>%
  filter(!is.na(Patient_Bam)) %>% # omit healthy controls since not plotted 
  dplyr::select(
    Patient,
    Date_of_sample_collection,
    Mut_source,
    Merge,
    detection_rate,
    sites_rate_zscore_charm
  ) %>%
  mutate(
    arm = case_when(
      Mut_source == "BM_cells"               ~ "BM",
      Mut_source == "Blood"      ~ "blood",
      TRUE                                             ~ NA_character_
    )
  ) %>%
  dplyr::select(-Mut_source) %>%
  pivot_wider(
    names_from  = arm,
    values_from = c(detection_rate, sites_rate_zscore_charm),
    names_sep   = "_"
  ) %>%
  rename(
    detect_rate_BM            = detection_rate_BM,
    detect_rate_blood         = detection_rate_blood,
    zscore_BM                 = sites_rate_zscore_charm_BM,
    zscore_blood              = sites_rate_zscore_charm_blood
  )


# Ensure tumor_fraction_dilution$Bam has .bam to match clinical_df$Bam
tumor_fraction_dilution <- tumor_fraction_dilution %>%
  mutate(BAM_full = paste0(Bam, ".bam")) %>%
  left_join(clinical_df, by = c("BAM_full" = "Bam"))

# ── Join everything together ─────────────────────────────────────────────
dilution_df <- mrdetect_wide %>%
  left_join(
    frag_df,
    by = c("Merge" = "Sample", "Patient" = "Patient")
  )

dilution_df <- dilution_df %>%
  rename(Sample = Merge)

tumor_fraction_dilution <- tumor_fraction_dilution %>% dplyr::select(-Bam) %>%
  rename(Bam = BAM)

dilution_df <- dilution_df %>% left_join(tumor_fraction_dilution %>% dplyr::select(Tumor_fraction, Bam)) 

dilution_df <- dilution_df %>%
  rename(WGS_Tumor_Fraction_Blood_plasma_cfDNA = Tumor_fraction) # for consistency

# ── 5. PREDICTION FUNCTION based on pre-trained model ────────────────────
predict_combo <- function(df, cmb, model) {
  preds <- combos[[cmb]]
  keep  <- complete.cases(df[, preds])
  prob  <- rep(NA_real_, nrow(df))
  
  if (any(keep)) {
    if (inherits(model, "cv.glmnet")) {
      mm <- model.matrix(reformulate(preds), df[keep, ])
      prob[keep] <- as.numeric(
        predict(model,
                newx = as(mm[, -1], "dgCMatrix"),
                s    = "lambda.min",
                type = "response")
      )
    } else {
      prob[keep] <- predict(model, newdata = df[keep, ], type = "response")
    }
  }
  prob
}

# ── 6. APPLY EACH SELECTED RULE TO DILUTION DATA ────────────────────────────
for (i in seq_len(nrow(selected_rows))) {
  cmb      <- selected_rows$combo[i]
  thr      <- selected_rows$threshold[i]
  fit_obj  <- model_list[[cmb]]
  suffix   <- if (cmb == "Blood_all_extras") paste0("_acc", i) else ""
  
  prob_col <- paste0(cmb, suffix, "_prob")
  call_col <- paste0(cmb, suffix, "_call")
  
  dilution_df[[prob_col]] <- predict_combo(dilution_df, cmb, fit_obj)
  dilution_df[[call_col]] <- if_else(
    dilution_df[[prob_col]] >= thr, 1L, 0L, NA_integer_
  )
}

# ── 7. SAVE SCORED DILUTION SERIES ──────────────────────────────────────────
write_rds(dilution_df, file.path(OUTPUT_DIR, "dilution_series_scored.rds"))
write_csv(dilution_df, file.path(OUTPUT_DIR, "dilution_series_scored.csv"))

message("Finished: results written to ", OUTPUT_DIR)
