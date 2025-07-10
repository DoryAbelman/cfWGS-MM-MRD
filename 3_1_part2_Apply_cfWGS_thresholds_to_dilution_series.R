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
library(virdis)

# ── 1. FILE PATHS ───────────────────────────────────────────────────────────
PATH_MODEL_LIST       <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_models_2025-06-19.rds"
PATH_THRESHOLD_LIST   <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Output_tables_2025/selected_combo_thresholds_2025-06-19.rds"
PATH_DILUTION_FRAGMENTOMICS         <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Results_Fragmentomics/Dilution_series/key_fragmentomics_info_dilution_series.rds"
PATH_DILUTION_PROCESSED_MRDetect    <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/MRDetect_output_winter_2025/Processed_R_outputs/cfWGS_Winter2025Dilution_series_May2025_with_zscore.rds"
PATH_DILUTION_CLINICAL  <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series/Metadata_dilution_series.csv"
PATH_TUMOR_FRACTION <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series/tumor_fraction_dilution_series.txt" # To add

OUTPUT_DIR            <- "Dilution_Series_Scoring_2025"
OUTPUT_DIR_TABLES     <- "Output_tables_2025"
OUTPUT_DIR_FIGURES    <- "Final Tables and Figures"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_DIR_FIGURES)) dir.create(OUTPUT_DIR_FIGURES, recursive = TRUE)


# ── 2. LOAD SAVED MODELS & THRESHOLDS ───────────────────────────────────────
selected_models <- readRDS(PATH_MODEL_LIST)
selected_thr    <- readRDS(PATH_THRESHOLD_LIST)

# ── 3. LOAD AND ASSEMBLE DILUTION SERIES DATA ────────────────────────────────
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
    sites_rate_zscore_charm,
    detection_rate_zscore_reads_checked_charm
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
    values_from = c(detection_rate, sites_rate_zscore_charm, detection_rate_zscore_reads_checked_charm),
    names_sep   = "_"
  ) %>%
  rename(
    detect_rate_BM            = detection_rate_BM,
    detect_rate_blood         = detection_rate_blood,
    z_score_detection_rate_BM = detection_rate_zscore_reads_checked_charm_BM, 
    zscore_BM                 = sites_rate_zscore_charm_BM,
    zscore_blood              = sites_rate_zscore_charm_blood,
    z_score_detection_rate_blood = detection_rate_zscore_reads_checked_charm_blood
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

# ── 6. APPLY EACH SELECTED RULE TO DILUTION DATA ────────────────────────────
dilution_df <- apply_selected(
  dat           = dilution_df,
  models        = selected_models,
  thresholds    = selected_thr,
  positive_class= "pos"
)

### Now edit the dilution LOD to the correct one 
# these are your true measured VAFs (%)
baseline_vaf <- 0.688    # baseline‐only library
neg_vaf      <- 0.0087   # MRD‐negative library
orig_full    <- 1.300    # the original “100%” point

dilution_df <- dilution_df %>%
  mutate(
    LOD_updated = (LOD / orig_full) * baseline_vaf +
      (1 - LOD / orig_full) * neg_vaf
  )

## For consistency
dilution_df$LOD_original <- dilution_df$LOD
dilution_df$LOD <- dilution_df$LOD_updated

# ── 7. SAVE SCORED DILUTION SERIES ──────────────────────────────────────────
write_rds(dilution_df, file.path(OUTPUT_DIR, "dilution_series_scored_updated2.rds"))
write_csv(dilution_df, file.path(OUTPUT_DIR, "dilution_series_scored_updated2.csv"))
write_csv(dilution_df, file.path(OUTPUT_DIR_TABLES, "dilution_series_scored_updated2.csv"))

message("Finished: results written to ", OUTPUT_DIR)



#### Now make the plots 
## ----------------------------------------------------------------------
## 1.  Candidate feature list  ------------------------------------------
##     (add / drop columns as you wish)
## ----------------------------------------------------------------------
select <- dplyr::select ## for conisistency

feature_cols <- c(
  "detect_rate_blood", "detect_rate_BM",
  "zscore_blood", "zscore_BM", "z_score_detection_rate_BM", "z_score_detection_rate_blood",
  "Mean.Coverage", "Midpoint.Coverage", "Midpoint.normalized",
  "Amplitude", "Zscore.Coverage", "Zscore.Midpoint", "Zscore.Amplitude",
  "Proportion.Short", "FS",
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
  "BM_zscore_only_sites_prob",
  "BM_zscore_only_detection_rate_prob",
  "Blood_zscore_only_detection_rate_prob",
  "Blood_rate_only_prob",
  "Blood_zscore_only_sites_prob",
  "Blood_plus_fragment_min_prob",
  "Fragmentomics_min_prob",
  "Fragmentomics_mean_coverage_only_prob"
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
write_csv(corr_tbl, file.path(OUTPUT_DIR_TABLES, "Supp_Table_LOD_feature_correlations_updated2.csv"))

## ----------------------------------------------------------------------
## 3.  Pick “interesting” features  -------------------------------------
##    (|r| > 0.80 & p < 0.05 by default; tweak as needed)
## ----------------------------------------------------------------------
sig_features <- corr_tbl %>%
  filter(abs(r) > 0.80, p < 0.05) %>%
  pull(feature)

sig_features_df <- corr_tbl %>%
  filter(abs(r) > 0.80, p < 0.05) 
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

ggsave("Final Tables and Figures/Fig_LOD_feature_correlations_updated.svg", p,
       width = 8, height = 6, dpi = 500)

## Make nicer figure for publication

# ───────────────────────────────────────────
# 1. BM Features
# ───────────────────────────────────────────
# 1) define your BM features of interest
bm_features <- c(
  "detect_rate_BM",           # BM mutation detection rate
  "zscore_BM",                # BM z-score
  "z_score_detection_rate_BM"           # BM mutation detection rate
#  "BM_zscore_only_detection_rate_prob",       # BM MRD probability
#  "BM_zscore_only_sites_prob"       # BM MRD probability
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
  BM_zscore_only_sites_prob                 = "BM sites model prob.",
  BM_zscore_only_detection_rate_prob             = "BM cVAF model prob.",
  detect_rate_BM                     = "BM cVAF",
  z_score_detection_rate_BM    = "BM cVAF Z-score",
  zscore_BM             = "BM sites Z-score"
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
    p   = cor.test(value, LOD, method = "spearman")$p.value,
    .groups = "drop"
  )

# 2) Prepare annotation
annot_spear <- corr_bm_spearman %>%
  mutate(
    p_text = if_else(p < 0.01, "p < 0.01", paste0("p = ", signif(p, 2))),
    label  = paste0("ρ = ", round(rho, 2), "\n", p_text)
  )

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
    hjust   = 0.1, vjust = 1.2,
    size    = 3.5
  ) +
  labs(
    x = "Tumour fraction (%)",
    y = "Feature value"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 4) Save to output directory
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_BM_metrics_spearman_actual_updated2.png"),
  plot     = p_bm_spearman_actual,
  width    = 8,
  height   = 4,
  dpi      = 500
)


# 3) Create plot with actual values and Spearman rho (log–log axes)
p_bm_spearman_actual <- ggplot(plot_df_bm, aes(x = LOD, y = value)) +
  geom_point(size = 2, alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, size = 0.7) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels)) +
  geom_text(
    data    = annot_spear,
    aes(x = 0.007, y = Inf, label = label),
    hjust   = 0, vjust = 1.2,
    size    = 3
  ) +
  scale_x_log10() +
  scale_y_continuous() +
  labs(
    x = "Log tumour fraction (%)",
    y = "Feature value"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 4) Save to output directory
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_BM_metrics_spearman_actual_updated2_loglog.png"),
  plot     = p_bm_spearman_actual,
  width    = 8,
  height   = 4, 
  dpi      = 500
)


 ## Add call features 
# 1) specify which of your features are the model‐probabilities
prob_features <- c(
  "BM_zscore_only_detection_rate_prob",  # BM cVAF Z-score model prob.
  "BM_zscore_only_sites_prob"            # BM sites Z-score model prob.
)

# 2) define your thresholds (one per feature)
thresholds <- tibble(
  feature = c("BM_zscore_only_sites_prob",
              "BM_zscore_only_detection_rate_prob"),
  thr     = c(
    bm_obj$thresholds["BM_zscore_only_sites"],
    bm_obj$thresholds["BM_zscore_only_detection_rate"]
  )
)

# 3) pull out only the prob data, join thresholds, and add a call flag
plot_df_prob <- plot_df %>%
  filter(feature %in% prob_features) %>%
  left_join(thresholds, by = "feature") %>%
  mutate(call = if_else(value > thr, "Positive", "Negative"))

# 4) make the ggplot
p_bm_prob <- ggplot(plot_df_prob, aes(x = LOD, y = value, color = call)) +
  # horizontal line at the threshold, drawn in the right facet
  geom_hline(
    data = thresholds,
    aes(yintercept = thr),
    linetype = "dashed",
    color    = "gray40"
  ) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels)) +
  scale_x_log10() +          # still log the x-axis if you like
  scale_color_manual(
    values = c(Positive = "forestgreen", Negative = "gray60"),
    na.value = "black"
  ) +
  labs(
    x     = "Log tumour fraction (%)",
    y     = "Model probability",
    color = "Call"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 5) preview or save
p_bm_prob

ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_BM_metrics_prob_calls.png"),
  plot     = p_bm_prob,
  width    = 8,
  height   = 4,
  dpi      = 500
)


combined_plot <- p_bm_spearman_actual + p_bm_prob +
  plot_layout(nrow = 1, widths = c(3, 2)) +
  plot_annotation(
    title = "Feature Concordance with Dilution Series",
    theme = theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin  = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(
    # this & theme() still applies to the individual panels
    panel.border   = element_rect(colour = "black", fill = NA),
    strip.text     = element_text(face = "bold")
  )

# draw it
print(combined_plot)

### Figure 4G
# and save
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig4G_LOD_combined.png"),
  plot     = combined_plot,
  width    = 12,        # adjust as needed
  height   = 4,         # one‐line panel
  dpi      = 600
)


#### Now do the overlap across all the features for 4H
# 1) Select just the features you care about, in the order you want them
selected_feats <- c(
  "detect_rate_BM",                       # BM mutation detection rate
  "zscore_BM",                            # BM sites z-score
  "z_score_detection_rate_BM",            # BM cVAF Z-score
  "Mean.Coverage",                        # enhancer coverage
  "Proportion.Short",                     # short-fragment proportion
  "FS",                                   # fragment‐size score
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA" # ichorCNA tumour fraction
)

# your custom labels
custom_labels <- c(
  detect_rate_BM                              = "BM cVAF",
  z_score_detection_rate_BM                   = "BM cVAF Z-score",
  zscore_BM                                   = "BM sites Z-score",
  FS                                          = "Fragment-size score",
  Mean.Coverage                               = "MM regulatory\ncoverage",
  Proportion.Short                            = "Short-fragment\nproportion",
  WGS_Tumor_Fraction_Blood_plasma_cfDNA       = "CNA tumour fraction\n(ichorCNA)"
)

# pick & sort
plot_df2 <- corr_tbl %>%
  filter(feature %in% names(custom_labels)) %>%
  arrange(desc(r)) %>%
  mutate(feature = factor(feature, levels = rev(feature))) 

p_corr_nice <- ggplot(plot_df2, aes(x = r, y = feature)) +
  # grey baseline at zero
  geom_vline(xintercept = 0, colour = "grey80", linetype = "dashed") +
  # segment from zero to rho
  geom_segment(aes(x = 0, xend = r, y = feature, yend = feature),
               colour = "grey70", size = 0.4) +
  # coloured dot
  geom_point(aes(colour = r), size = 3) +
  scale_colour_viridis_c(
    option = "D", end = 0.9,
    name = expression(rho~"(Spearman)"),
    limits = c(-1,1)
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_discrete(
    labels = custom_labels
  ) +
  labs(
    title = "Feature Correlation with Dilution Series",
    x     = expression(rho~"(Spearman)"),
    y     = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title.position = "plot",
    plot.title        = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_text(size = 8),
    axis.title.x      = element_text(size = 8, margin = margin(t = 4)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    legend.key.height  = unit(1, "cm"),
    legend.key.width   = unit(0.3, "cm"),
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8)
  ) 

# print or save
print(p_corr_nice)
ggsave(filename = file.path(OUTPUT_DIR_FIGURES, "Fig4H_feature_corr_lollipop_nice.png"), 
       p_corr_nice, 
       width = 3.55, height = 4, dpi = 600)









# ───────────────────────────────────────────
# 2. Blood Features
# ───────────────────────────────────────────
## Now redo for blood derived muts
blood_features <- c(
  "detect_rate_blood",                          
  "zscore_blood",                               
  "z_score_detection_rate_blood",              
  "Blood_zscore_only_detection_rate_prob",   
  "Blood_zscore_only_sites_prob",
 "Blood_rate_only_prob",                      
  "Blood_plus_fragment_min_prob"
)


# 2) prepare long data
plot_df_blood <- dilution_df %>%
  select(LOD, all_of(blood_features)) %>%
  pivot_longer(-LOD, names_to = "feature", values_to = "value")

# 3) get correlations & annotate
corr_blood <- plot_df_blood %>%
  group_by(feature) %>%
  summarise(
    r   = cor(value, LOD, method = "spearman"),
    r2  = r^2,
    p   = cor.test(value, LOD, method = "spearman")$p.value,
    .groups = "drop"
  )

annot_df <- corr_blood %>%
  mutate(label = sprintf("r² = %.2f", r2))

# 4) nicer facet labels
facet_labels <- c(
  Blood_zscore_only_sites_prob                 = "Blood sites model prob.",
  Blood_zscore_only_detection_rate_prob             = "Blood cVAF model prob.",
  detect_rate_blood                     = "Blood cVAF",
  z_score_detection_rate_blood    = "Blood cVAF Z-score",
  zscore_blood             = "Blood sites Z-score"
)
# 5) clean ggplot - pearson 
p_blood <- ggplot(plot_df_blood, aes(x = LOD, y = value)) +
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
corr_blood_spearman <- plot_df_blood %>%
  group_by(feature) %>%
  summarise(
    rho = cor(value, LOD, method = "spearman"),
    p   = cor.test(value, LOD, method = "spearman")$p.value,
    .groups = "drop"
  )

# 2) Prepare annotation
annot_spear <- corr_blood_spearman %>%
  mutate(
    p_text = if_else(p < 0.01, "p < 0.01", paste0("p = ", signif(p, 2))),
    label  = paste0("ρ = ", round(rho, 2), "\n", p_text)
  )

# 3) Create plot with actual values and Spearman rho
p_blood_spearman_actual <- ggplot(plot_df_blood, aes(x = LOD, y = value)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, size = 0.7) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels)) +
  geom_text(
    data    = annot_spear,
    aes(x = Inf, y = Inf, label = label),
    hjust   = 0.1, vjust = 1.2,
    size    = 3.5
  ) +
  labs(
    x = "Tumour fraction (%)",
    y = "Feature value"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 4) Save to output directory
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_blood_metrics_spearman_actual_updated2.png"),
  plot     = p_blood_spearman_actual,
  width    = 8,
  height   = 4,
  dpi      = 500
)


# 3) Create plot with actual values and Spearman rho (log–log axes)
p_blood_spearman_actual <- ggplot(plot_df_blood %>% filter(feature %in% c("detect_rate_blood", "z_score_detection_rate_blood", "zscore_blood")), aes(x = LOD, y = value)) +
  geom_point(size = 2, alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, size = 0.7) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels)) +
  geom_text(
    data    = annot_spear %>% filter(feature %in% c("detect_rate_blood", "z_score_detection_rate_blood", "zscore_blood")),
    aes(x = 0.007, y = Inf, label = label),
    hjust   = 0, vjust = 1.2,
    size    = 3
  ) +
  scale_x_log10() +
  scale_y_continuous() +
  labs(
    x = "Log tumour fraction (%)",
    y = "Feature value"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 4) Save to output directory
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_blood_metrics_spearman_actual_updated2_loglog.png"),
  plot     = p_blood_spearman_actual,
  width    = 8,
  height   = 4, 
  dpi      = 500
)


## Add call features 
# 1) specify which of your features are the model‐probabilities
prob_features <- c(
  "Blood_zscore_only_detection_rate_prob",  # Blood cVAF Z-score model prob.
  "Blood_zscore_only_sites_prob"            # Blood sites Z-score model prob.
)

# 2) define your thresholds (one per feature)
thresholds <- tibble(
  feature = c("Blood_zscore_only_sites_prob",
              "Blood_zscore_only_detection_rate_prob"),
  thr     = c(
    blood_obj$thresholds["Blood_zscore_only_sites"],
    blood_obj$thresholds["Blood_zscore_only_detection_rate"]
  )
)

# 3) pull out only the prob data, join thresholds, and add a call flag
plot_df_prob <- plot_df_blood %>%
  filter(feature %in% prob_features) %>%
  left_join(thresholds, by = "feature") %>%
  mutate(call = if_else(value > thr, "Positive", "Negative"))

# 4) make the ggplot
p_blood_prob <- ggplot(plot_df_prob, aes(x = LOD, y = value, color = call)) +
  # horizontal line at the threshold, drawn in the right facet
  geom_hline(
    data = thresholds,
    aes(yintercept = thr),
    linetype = "dashed",
    color    = "gray40"
  ) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels)) +
  scale_x_log10() +          # still log the x-axis if you like
  scale_color_manual(
    values = c(Positive = "forestgreen", Negative = "gray60"),
    na.value = "black"
  ) +
  labs(
    x     = "Log tumour fraction (%)",
    y     = "Model probability",
    color = "Call"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 5) preview or save
p_blood_prob

ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_blood_metrics_prob_calls.png"),
  plot     = p_blood_prob,
  width    = 8,
  height   = 4,
  dpi      = 500
)


combined_plot <- p_blood_spearman_actual + p_blood_prob +
  plot_layout(nrow = 1, widths = c(3, 2)) +
  plot_annotation(
    title = "Feature Concordance with Dilution Series",
    theme = theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin  = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(
    # this & theme() still applies to the individual panels
    panel.border   = element_rect(colour = "black", fill = NA),
    strip.text     = element_text(face = "bold")
  )

# draw it
print(combined_plot)

### Figure 5G
# and save
ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig5G_LOD_combined.png"),
  plot     = combined_plot,
  width    = 12,        # adjust as needed
  height   = 4,         # one‐line panel
  dpi      = 600
)

### Add the confirm threshold 
plot_df_prob <- plot_df_blood %>%
  filter(feature %in% prob_features) %>%
  left_join(thresholds, by = "feature") %>%
  mutate(
    call = case_when(
      ## SITES pane: 3 categories
      feature == "Blood_zscore_only_sites_prob" & value >= thr                ~ "Positive",
      feature == "Blood_zscore_only_sites_prob" & value >= 0.457 & value < thr ~ "Confirmatory",
      feature == "Blood_zscore_only_sites_prob"                                   ~ "Negative",
      
      ## Other features: just Positive vs Negative
      value >= thr ~ "Positive",
      TRUE         ~ "Negative"
    ),
    call = factor(call, levels = c("Negative","Confirmatory","Positive"))
  )

p_blood_prob <- ggplot(plot_df_prob, aes(x = LOD, y = value, color = call)) +
  # model-threshold lines for each facet
  geom_hline(
    data = thresholds,
    aes(yintercept = thr),
    linetype = "dashed",
    color    = "gray40"
  ) +
  # confirmatory cutoff (0.457) only in the sites facet
  geom_hline(
    data = data.frame(feature = "Blood_zscore_only_sites_prob",
                      thr     = 0.457),
    aes(yintercept = thr),
    linetype = "dashed",
    color    = "steelblue"
  ) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels)
  ) +
  scale_x_log10() +
  scale_color_manual(
    name   = "Call",
    values = c(Negative     = "gray60",
               Confirmatory = "steelblue",
               Positive     = "forestgreen")
  ) +
  labs(
    x     = "Log tumour fraction (%)",
    y     = "Model probability",
    color = "Call"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black"),
    legend.position = "right"
  )

combined_plot <- p_blood_spearman_actual + p_blood_prob +
  plot_layout(nrow = 1, widths = c(3, 2)) +
  plot_annotation(
    title = "Feature Concordance with Dilution Series",
    theme = theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin  = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(
    # this & theme() still applies to the individual panels
    panel.border   = element_rect(colour = "black", fill = NA),
    strip.text     = element_text(face = "bold")
  )


ggsave(
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig5G_LOD_combined_updated.png"),
  plot     = combined_plot,
  width    = 12,        # adjust as needed
  height   = 4,         # one‐line panel
  dpi      = 600
)



#### Now do the overlap across all the features for 5H
# 1) Select just the features you care about, in the order you want them
selected_feats <- c(
  "detect_rate_blood",                       # BM mutation detection rate
  "zscore_blood",                            # BM sites z-score
  "z_score_detection_rate_blood",            # BM cVAF Z-score
  "Mean.Coverage",                        # enhancer coverage
  "Proportion.Short",                     # short-fragment proportion
  "FS",                                   # fragment‐size score
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA" # ichorCNA tumour fraction
)

# your custom labels
custom_labels <- c(
  detect_rate_blood                              = "Blood cVAF",
  z_score_detection_rate_blood                   = "Blood cVAF Z-score",
  zscore_blood                                   = "Blood sites Z-score",
  FS                                          = "Fragment-size score",
  Mean.Coverage                               = "MM regulatory\ncoverage",
  Proportion.Short                            = "Short-fragment\nproportion",
  WGS_Tumor_Fraction_Blood_plasma_cfDNA       = "CNA tumour fraction\n(ichorCNA)"
)

# pick & sort
plot_df2 <- corr_tbl %>%
  filter(feature %in% names(custom_labels)) %>%
  arrange(desc(r)) %>%
  mutate(feature = factor(feature, levels = rev(feature))) 

p_corr_nice <- ggplot(plot_df2, aes(x = r, y = feature)) +
  # grey baseline at zero
  geom_vline(xintercept = 0, colour = "grey80", linetype = "dashed") +
  # segment from zero to rho
  geom_segment(aes(x = 0, xend = r, y = feature, yend = feature),
               colour = "grey70", size = 0.4) +
  # coloured dot
  geom_point(aes(colour = r), size = 3) +
  scale_colour_viridis_c(
    option = "D", end = 0.9,
    name = expression(rho~"(Spearman)"),
    limits = c(-1,1)
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_discrete(
    labels = custom_labels
  ) +
  labs(
    title = "Feature Correlation with Dilution Series",
    x     = expression(rho~"(Spearman)"),
    y     = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title.position = "plot",
    plot.title        = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_text(size = 8),
    axis.title.x      = element_text(size = 8, margin = margin(t = 4)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    legend.key.height  = unit(1, "cm"),
    legend.key.width   = unit(0.3, "cm"),
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8)
  ) 

# print or save
print(p_corr_nice)
ggsave(filename = file.path(OUTPUT_DIR_FIGURES, "Fig5H_feature_corr_lollipop_nice.png"), 
       p_corr_nice, 
       width = 3.55, height = 4, dpi = 600)

## Rescore based on new LOD 
dilution_df <- dilution_df %>%
  mutate(
    Blood_zscore_only_sites_call_rescored = if_else(
      Blood_zscore_only_sites_prob >= 0.457,
      1L,    # “positive” call
      0L     # “negative” call
    )
  )

# ───────────────────────────────────────────
# 3. Fragmentomics Features
# ───────────────────────────────────────────
## Update this - just template 
fragmentomics_features <- c(
  "FS",                                        
  "Proportion.Short",                          
  "Mean.Coverage",                             
  "WGS_Tumor_Fraction_Blood_plasma_cfDNA",     
  "Fragmentomics_min_prob",                    
  "Fragmentomics_mean_coverage_only_prob"
)


# 2) prepare long data
plot_df_blood <- dilution_df %>%
  select(LOD, all_of(blood_features)) %>%
  pivot_longer(-LOD, names_to = "feature", values_to = "value")

# 3) get correlations & annotate
corr_blood <- plot_df_blood %>%
  group_by(feature) %>%
  summarise(
    r   = cor(value, LOD, method = "spearman"),
    r2  = r^2,
    p   = cor.test(value, LOD, method = "spearman")$p.value,
    .groups = "drop"
  )

annot_df <- corr_blood %>%
  mutate(label = sprintf("r² = %.2f", r2))

# 4) nicer facet labels
facet_labels_fragmentomics <- c(
  FS                                            = "Fragment size\nscore",
  Proportion.Short                              = "Short fragment\nproportion",
  Mean.Coverage                                 = "Mean coverage",
  WGS_Tumor_Fraction_Blood_plasma_cfDNA         = "Tumor fraction",
  Fragmentomics_min_prob                        = "Fragmentomics\nProbability (min)",
  Fragmentomics_mean_coverage_only_prob         = "Fragmentomics\nProbability (mean coverage)"
)

# 5) clean ggplot - pearson 
p_blood <- ggplot(plot_df_blood, aes(x = LOD, y = value)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, size = 0.7) +
  facet_wrap(~ feature,
             scales = "free_y",
             labeller = as_labeller(facet_labels_blood)) +
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
corr_blood_spearman <- plot_df_blood %>%
  group_by(feature) %>%
  summarise(
    rho = cor(value, LOD, method = "spearman"),
    .groups = "drop"
  )

# 2) Prepare annotation
annot_spear <- corr_blood_spearman %>%
  mutate(label = sprintf("ρ = %.2f", rho))

# 3) Create plot with actual values and Spearman rho
p_blood_spearman_actual <- ggplot(plot_df_blood, aes(x = LOD, y = value)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, size = 0.7) +
  facet_wrap(~ feature,
             scales   = "free_y",
             labeller = as_labeller(facet_labels_blood)) +
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
  filename = file.path(OUTPUT_DIR_FIGURES, "Fig_LOD_blood_metrics_spearman_actual_updated.png"),
  plot     = p_blood_spearman_actual,
  width    = 10,
  height   = 5,
  dpi      = 500
)
