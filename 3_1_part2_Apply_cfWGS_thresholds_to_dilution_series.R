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
library(broom)       # tidying model outputs
library(patchwork)   # optional – for combining plots, if needed


# ── 1. FILE PATHS ───────────────────────────────────────────────────────────
PATH_MODEL_LIST       <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_models.rds"
PATH_THRESHOLD_LIST   <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_thresholds.rds"
PATH_DILUTION_FRAGMENTOMICS         <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Results_Fragmentomics/Dilution_series/key_fragmentomics_info_dilution_series.rds"
PATH_DILUTION_PROCESSED_MRDetect    <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/MRDetect_output_winter_2025/Processed_R_outputs/cfWGS_Winter2025Dilution_series_May2025_with_zscore.rds"
PATH_DILUTION_CLINICAL  <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series/Metadata_dilution_series.csv"
PATH_TUMOR_FRACTION <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series/tumor_fraction_dilution_series.txt" # To add

OUTPUT_DIR            <- "Dilution_Series_Scoring_2025"
OUTPUT_DIR_TABLES     <- "Output_tables_2025"
OUTPUT_DIR_FIGURES    <- "Output_figures_2025"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_DIR_FIGURES)) dir.create(OUTPUT_DIR_FIGURES, recursive = TRUE)


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

### Rename for clarity 
dilution_df <- dilution_df %>%
  rename(
    BloodSensPriority_prob = Blood_all_extras_acc2_prob,
    BloodSensPriority_call = Blood_all_extras_acc2_call,
    BloodAccPriority_prob   = Blood_all_extras_acc3_prob,
    BloodAccPriority_call   = Blood_all_extras_acc3_call
  )

# ── 7. SAVE SCORED DILUTION SERIES ──────────────────────────────────────────
write_rds(dilution_df, file.path(OUTPUT_DIR, "dilution_series_scored.rds"))
write_csv(dilution_df, file.path(OUTPUT_DIR, "dilution_series_scored.csv"))
write_csv(dilution_df, file.path(OUTPUT_DIR_TABLES, "dilution_series_scored.csv"))

message("Finished: results written to ", OUTPUT_DIR)



#### Now make the plots 
## ----------------------------------------------------------------------
## 1.  Candidate feature list  ------------------------------------------
##     (add / drop columns as you wish)
## ----------------------------------------------------------------------
select <- dplyr::select ## for conisistency

feature_cols <- c(
  "detect_rate_blood", "detect_rate_BM",
  "zscore_blood", "zscore_BM",
  "Mean.Coverage", "Midpoint.Coverage", "Midpoint.normalized",
  "Amplitude", "Zscore.Coverage", "Zscore.Midpoint", "Zscore.Amplitude",
  "Proportion.Short", "FS",
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
  "BM_zscore_only_prob",
  "BloodSensPriority_prob", "BloodAccPriority_prob", "Blood_base_prob"
)

## ----------------------------------------------------------------------
## 2.  Long‑format data frame & correlation screen  ---------------------
## ----------------------------------------------------------------------
plot_df <- dilution_df %>%
  select(LOD, all_of(feature_cols)) %>%
  pivot_longer(-LOD, names_to = "feature", values_to = "value") %>%
  filter(!is.na(value))

corr_tbl <- plot_df %>%
  group_by(feature) %>%
  summarise(
    r  = cor(value, LOD, method = "spearman"),
    p  = cor.test(value, LOD, method = "spearman")$p.value,
    r2 = r^2,
    .groups = "drop"
  ) %>%
  arrange(desc(r2))

## save the full table for the supplement
write_csv(corr_tbl, file.path(OUTPUT_DIR_TABLES, "Supp_Table_LOD_feature_correlations.csv"))

## ----------------------------------------------------------------------
## 3.  Pick “interesting” features  -------------------------------------
##    (|r| > 0.80 & p < 0.05 by default; tweak as needed)
## ----------------------------------------------------------------------
sig_features <- corr_tbl %>%
  filter(abs(r) > 0.80, p < 0.05) %>%
  pull(feature)

## ----------------------------------------------------------------------
## 4.  Build the figure  -------------------------------------------------
## ----------------------------------------------------------------------
sig_plot_df <- plot_df %>% filter(feature %in% sig_features)

# annotate each panel with the correlation coefficient
annot_df <- corr_tbl %>%
  filter(feature %in% sig_features) %>%
  mutate(label = sprintf("r = %.2f", r))

p <- ggplot(sig_plot_df, aes(x = LOD, y = value)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_wrap(~ feature, scales = "free_y") +
  geom_text(
    data = annot_df, aes(x = Inf, y = Inf, label = label),
    hjust = 1.05, vjust = 1.3, size = 3.8
  ) +
  labs(
    x = "Dilution (estimated tumour fraction, %)",
    y = "Feature value"
  ) +
  theme_bw(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))

ggsave("Fig_LOD_feature_correlations.svg", p,
       width = 8, height = 6, dpi = 500)

## Make nicer figure for publication
# 1) define your BM features of interest
bm_features <- c(
  "detect_rate_BM",           # BM mutation detection rate
  "zscore_BM",                # BM z-score
  "BM_zscore_only_prob"       # BM MRD probability
)

# 2) prepare long data
plot_df_bm <- dilution_df %>%
  select(LOD, all_of(bm_features)) %>%
  pivot_longer(-LOD, names_to = "feature", values_to = "value")

# 3) get correlations & annotate
corr_bm <- plot_df_bm %>%
  group_by(feature) %>%
  summarise(
    r   = cor(value, LOD, method = "spearman"),
    r2  = r^2,
    p   = cor.test(value, LOD, method = "spearman")$p.value,
    .groups = "drop"
  )

annot_df <- corr_bm %>%
  mutate(label = sprintf("r² = %.2f", r2))

# 4) nicer facet labels
facet_labels <- c(
  detect_rate_BM         = "BM mutation\ndetection rate",
  zscore_BM              = "BM sites z-score",
  BM_zscore_only_prob    = "BM MRD probability"
)

# 5) clean ggplot - pearson 
p_bm <- ggplot(plot_df_bm, aes(x = LOD, y = value)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, size = 0.7) +
  facet_wrap(~ feature,
             scales = "free_y",
             labeller = as_labeller(facet_labels)) +
  geom_text(
    data = annot_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.2,
    size = 3.5
  ) +
  labs(
    x = "Dilution tumour fraction (%)",
    y = "Values"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text        = element_text(face = "bold"),
    panel.border      = element_rect(color = "black", fill = NA),
    axis.title.x      = element_text(margin = margin(t = 8)),
    axis.ticks        = element_line(color = "black")
  )

## For Spearman
# 1) Compute Spearman correlation (on raw data, using ranks internally)
corr_bm_spearman <- plot_df_bm %>%
  group_by(feature) %>%
  summarise(
    rho = cor(value, LOD, method = "spearman"),
    .groups = "drop"
  )

# 2) Prepare annotation
annot_spear <- corr_bm_spearman %>%
  mutate(label = sprintf("ρ = %.2f", rho))

# 3) Create plot with actual values and Spearman rho
p_bm_spearman_actual <- ggplot(plot_df_bm, aes(x = LOD, y = value)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, size = 0.7) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels)) +
  geom_text(
    data    = annot_spear,
    aes(x = Inf, y = Inf, label = label),
    hjust   = 1.1, vjust = 1.2,
    size    = 3.5
  ) +
  labs(
    x = "Tumour fraction (%)",
    y = "Feature value"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 4) Save to output directory
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_BM_metrics_spearman_actual.png"),
  plot     = p_bm_spearman_actual,
  width    = 10,
  height   = 5,
  dpi      = 500
)

