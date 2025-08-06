# ──────────────────────────────────────────────────────────────────────────────
# 1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R
#
# Purpose:
#   • Load and process “insert_size_summary.tsv” and “fragment_scores.tsv” for cfWGS samples 
#     (MM baseline/monitoring) and PON (healthy controls).
#   • Harmonize sample names, merge with “combined_clinical_data_updated_April2025.csv” 
#   • Define comparison groups:
#        – Diagnosis/Baseline vs Monitoring vs Healthy (insert size)
#   • For insert‐size (“Proportion.Short”):
#        – Test normality (Shapiro–Wilk), run ANOVA or Kruskal–Wallis
#        – Run all pairwise t‐tests (or Wilcoxon) and save results to CSV.
#   • For fragment‐score (“FS”):
#        – Merge Fragment Score to clinical data
#        – Deduplicate, save “combined_data_fragmentomics_cfWGS.csv”
#   • Load MM‐DARs chromatin activation metrics (MM_DARs_chromatin_activation_data.csv),
#     join into the fragment‐score data frame, and export “Key_fragmentomics_data_updated.csv”
#
# Outputs (to working directory):
#   • Pairwise_t_Test_Results_fragment_scores.csv
#   • Key_fragmentomics_data_updated.csv
#
# Notes:
#   – Plotting will happen later in a separate script.
# ──────────────────────────────────────────────────────────────────────────────

### PREPARE SESSION ################################################################################
library(BoutrosLab.plotting.general)
library(ggplot2)        # used only for tibble/data‐frame handling (no plots)
library(conflicted)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(stringr)
library(rstatix)        # for pairwise_wilcox_test, if needed
library(tidyr)

source('session.functions.R')


### SET PATHS #######################################################################################
input.dir <- '~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data'
pon.dir   <- '~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Normals'
out.dir   <- '~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Results_Fragmentomics/Insert_size'


### STEP 1: PROCESS INSERT‐SIZE (“Proportion.Short”) ###############################################
# – Read all “insert_size_summary.tsv” files (cfWGS + PON)
results.files <- list.files(path = input.dir,
                            pattern = "insert_size_summary.tsv",
                            full.names = TRUE)
pon.files <- rev(sort(
  list.files(path = pon.dir,
             pattern = "insert_size_summary.tsv",
             full.names = TRUE)
))

cfWGS_insert <- results.files %>%
  lapply(read_tsv) %>%
  bind_rows()

pon_insert <- pon.files %>%
  lapply(read_tsv) %>%
  bind_rows()

# Load clinical tables
combined_clinical <- read_csv("combined_clinical_data_updated_April2025.csv")

# Construct a “Merge_ID” for joining
tmp.clin <- combined_clinical %>%
  mutate(
    Merge_ID = if_else(
      str_starts(Bam, "SPORE"),
      str_extract(Bam, "^[^\\.]+"),
      str_replace(Sample_ID, "_Blood_plasma_cfDNA", "-P")
    )
  )

# Standardize cfWGS_insert sample names, then join to clinical keys
cfWGS_insert <- cfWGS_insert %>%
  mutate(
    Sample = str_replace(Sample, "_Blood_plasma_cfDNA", "-P"),
    Sample = if_else(
      !str_starts(Sample, "SPORE"),
      str_replace_all(Sample, "_", "-"),
      Sample
    ),
    Sample = str_replace(Sample, "IMG", "MyP")
  ) %>%
  left_join(tmp.clin, by = c("Sample" = "Merge_ID"))

# Mark PON samples as “Healthy”
pon_insert <- pon_insert %>% mutate(Group = "Healthy")

# Define “Monitoring” vs “Diagnosis/Baseline” groups
group_1 <- c("Post_transplant", "Post_induction", "Maintenance",
             "1.5yr maintenance", "2yr maintenance")
group_2 <- c("Diagnosis", "Baseline")

data_insert <- cfWGS_insert %>%
  mutate(Group = ifelse(timepoint_info %in% group_1,
                        "Monitoring", "Diagnosis/Baseline"))

# Read patient lists
patients_with_either <- read_csv("patients_with_either_cfDNA_at_baseline_and_monitoring.csv")
patients_with_both   <- read_csv("patients_with_cfDNA_at_baseline_and_monitoring.csv")

# Combine cfWGS + PON (dropping one duplicate PON sample if present)
pon_insert <- pon_insert %>%
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample))

# Prepare your pon_insert with the same columns (as NA)
pon_insert2 <- pon_insert %>%
  mutate(
    Patient        = NA_character_,   # same type as data_insert$Patient
    timepoint_info = NA_character_    # same type as data_insert$timepoint_info
  )

# Bind together
all_insert <- bind_rows(
  data_insert %>%
    select(Sample, Patient, timepoint_info, Group, Proportion.Short),
  
  pon_insert2 %>%
    select(Sample, Patient, timepoint_info, Group, Proportion.Short)
)


# TEST NORMALITY (Shapiro–Wilk per group)
shapiro_results <- all_insert %>%
  group_by(Group) %>%
  summarise(p_value = shapiro.test(Proportion.Short)$p.value)

print(shapiro_results)

# DECIDE whether to use ANOVA or Kruskal–Wallis based on shapiro_results
if (all(shapiro_results$p_value > 0.05)) {
  # parametric ANOVA
  test_insert_global <- aov(Proportion.Short ~ Group, data = all_insert)
  print(summary(test_insert_global))
} else {
  # nonparametric Kruskal–Wallis
  test_insert_global <- kruskal.test(Proportion.Short ~ Group, data = all_insert)
  print(test_insert_global)
}

# PAIRWISE comparisons – use t.test if normal, else Wilcoxon
pairwise_t_results <- list()
groups <- unique(all_insert$Group)

for (i in seq_len(length(groups) - 1)) {
  for (j in seq((i + 1), length(groups))) {
    g1 <- groups[i]
    g2 <- groups[j]
    vec1 <- all_insert %>%
      filter(Group == g1) %>%
      pull(Proportion.Short)
    vec2 <- all_insert %>%
      filter(Group == g2) %>%
      pull(Proportion.Short)
    
    if (all(shapiro_results$p_value > 0.05)) {
      ttest <- t.test(vec1, vec2, na.rm = TRUE)
      pval <- ttest$p.value
      method_used <- "t.test"
    } else {
      wt <- wilcox.test(vec1, vec2, na.rm = TRUE)
      pval <- wt$p.value
      method_used <- "wilcox.test"
    }
    
    pairwise_t_results[[paste(g1, "vs", g2)]] <- data.frame(
      Group1 = g1,
      Group2 = g2,
      Mean_Group1 = mean(vec1, na.rm = TRUE),
      Mean_Group2 = mean(vec2, na.rm = TRUE),
      Mean_Difference = mean(vec1, na.rm = TRUE) - mean(vec2, na.rm = TRUE),
      P_Value = pval,
      CI_Lower = if (method_used == "t.test") ttest$conf.int[1] else NA,
      CI_Upper = if (method_used == "t.test") ttest$conf.int[2] else NA,
      Method = method_used,
      stringsAsFactors = FALSE
    )
  }
}

pairwise_t_results_df <- do.call(rbind, pairwise_t_results)

write_csv(
  pairwise_t_results_df,
  file = "Pairwise_t_Test_Results_fragment_scores.csv"
)


### STEP 2: PROCESS FRAGMENT SCORE (“FS”) ############################################################
# – Read all “fragment_scores.tsv” (cfWGS + PON)
fs.files <- list.files(path = input.dir,
                       pattern = "fragment_scores.tsv",
                       full.names = TRUE)
pon.fs.files <- rev(sort(
  list.files(path = pon.dir,
             pattern = "fragment_scores.tsv",
             full.names = TRUE)
))

cfWGS_fs <- fs.files %>%
  lapply(read_tsv) %>%
  bind_rows()

# Drop any “TGL49” rows if they snuck in
if (any(grepl("TGL49", cfWGS_fs$Sample))) {
  cfWGS_fs <- cfWGS_fs[!grepl("TGL49", cfWGS_fs$Sample), ]
}

pon_fs <- read_tsv(pon.fs.files[[1]])

# Standardize sample names & join clinical keys just as above
cfWGS_fs <- cfWGS_fs %>%
  mutate(
    Sample = str_replace(Sample, "_Blood_plasma_cfDNA", "-P"),
    Sample = if_else(
      !str_starts(Sample, "SPORE"),
      str_replace_all(Sample, "_", "-"),
      Sample
    ),
    Sample = str_replace(Sample, "IMG", "MyP")
  ) %>%
  left_join(tmp.clin %>% dplyr::select(Patient, Timepoint, timepoint_info, Merge_ID), by = c("Sample" = "Merge_ID")) %>%
  mutate(Cohort = "cfWGS")

cfWGS_fs <- cfWGS_fs %>% unique()

pon_fs <- pon_fs %>%
  mutate(Cohort = "HBC") %>%
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample))

# Combine them
combined_fs <- bind_rows(cfWGS_fs, pon_fs)

# Join MRD data and drop unnecessary columns
combined_fs <- combined_fs %>%
 # dplyr::select(-Relapsed, -Num_days_to_closest_relapse, -Num_days_to_closest_relapse_absolute, -Num_days_to_closest_relapse_non_absolute) %>%
  left_join(combined_clinical %>% select(-Bam),
            by = c("Patient", "Timepoint", "timepoint_info"))

combined_fs <- combined_fs %>% dplyr::select(-`...1`)
  
# Deduplicate numeric fields (take mean) and non‐numeric (first non‐NA)
combined_fs_dedup <- combined_fs %>%
  group_by(Sample) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    across(where(~ !is.numeric(.x)), ~ first(.x)),
    .groups = "drop"
  )




# STEP 2b: DEFINE HEALTHY CONTROL CUTOFF FOR FRAGMENT SCORE ("FS")

# 1) Pull FS from healthy (PON) samples
fs_healthy <- pon_fs %>%
  filter(Cohort == "HBC") %>%
  pull(FS) %>%
  na.omit()

# 2) Test normality (Shapiro–Wilk)
shp_p <- shapiro.test(fs_healthy)$p.value
message(sprintf("Shapiro–Wilk p‑value for healthy FS = %.3g", shp_p))

# 3) Choose cut‑offs
if (shp_p > 0.05) {
  # Normal‑ish → use mean ± 1.96 SD (≈ 95 % reference range)
  mu  <- mean(fs_healthy)
  sig <- sd(fs_healthy)
  lower_cut <- mu - 1.96 * sig
  upper_cut <- mu + 1.96 * sig
  method_used <- "mean ± 1.96·SD"
} else {
  # Non‑normal → use empirical 2.5th and 97.5th percentiles
  qs <- quantile(fs_healthy, probs = c(0.025, 0.975), na.rm = TRUE)
  lower_cut <- qs[1]; upper_cut <- qs[2]
  method_used <- "2.5th / 97.5th percentiles"
}

cat(sprintf(
  "FS cut‑offs (%s): lower = %.3f, upper = %.3f\n",
  method_used, lower_cut, upper_cut
))

## Now export
# 1) Build a tibble with  cut-off info
fs_cutoffs_tbl <- tibble::tibble(
  Method        = method_used,
  Lower_Cutoff  = lower_cut,
  Upper_Cutoff  = upper_cut
)

# 2) (Optional) print to console
print(fs_cutoffs_tbl)

# 3) Write out to CSV in your output directory
readr::write_csv(
  fs_cutoffs_tbl,
  file.path(outdir, "FS_cutoffs_table_healthy_controls.csv")
)





### STEP 3: MERGE MM-DARs CHROMATIN ACTIVATION METRICS ###############################################
# – Read “MM_DARs_chromatin_activation_data” (generated by 1_7A script)
mm_dars <- readRDS("Results_Fragmentomics/MM_DARs_chromatin_activation_data.rds")

## In case run in two batches (optional)
mm_dars2 <- read.csv("Results_Fragmentomics/MM_Griffin_all_relevant_sites_data_updated_SPORE.csv")

mm_dars2_small <- mm_dars2 %>%
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection),
    Date                     = as.Date(Date)
  ) %>%
  select(
    Sample, Site,
    Mean.Coverage, Midpoint.Coverage, Midpoint.normalized,
    Amplitude, Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude,
    Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
    Bam, Patient, Date_of_sample_collection
  ) %>%
  distinct() %>%
  filter(Site == "MM_DARs_chromatin_activation")

# Keep only the “MM_DARs_chromatin_activation” site rows
mm_dars_small <- mm_dars %>%
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection, format = "%Y-%m-%d"),
    Date                     = as.Date(Date, format = "%Y-%m-%d")
  ) %>%
  select(
    Sample, Site,
    Mean.Coverage, Midpoint.Coverage, Midpoint.normalized,
    Amplitude, Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude,
    Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
    Bam, Patient, Date_of_sample_collection
  ) %>%
  distinct() %>%
  filter(Site == "MM_DARs_chromatin_activation")

# Combine (optional)
mm_dars_combined <- bind_rows(mm_dars_small, mm_dars2_small) %>%
  # normalize your date columns
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection)
  ) %>%
  # pick first occurrence per key
  group_by(Sample, Site, Bam, Date_of_sample_collection) %>%
  slice_head(n = 1) %>%
  ungroup()

# For consistency
mm_dars_small <- mm_dars_combined %>% unique()

# Merge into “combined_fs_dedup” by keys (Sample, Patient, Date_of_sample_collection)
# (First, restore combined_fs_dedup into a data frame we can join onto)
fragdata_for_merge <- combined_fs_dedup

fragdata_for_merge <- fragdata_for_merge %>% dplyr::select(Sample, Patient, Date_of_sample_collection, FS)


merged_data <- fragdata_for_merge %>%
  full_join(mm_dars_small %>% dplyr::select(Sample, Patient, Date_of_sample_collection, Site, Mean.Coverage, Midpoint.Coverage, Midpoint.normalized, Amplitude, Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude, Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude),
            by = c("Sample", "Patient", "Date_of_sample_collection"))


### If minor difference in dates due to different versions of log being used 
merged_data_grouped <- merged_data %>%
  arrange(Patient, Date_of_sample_collection) %>%
  group_by(Patient) %>%
  mutate(
    date_diff = as.numeric(Date_of_sample_collection - lag(Date_of_sample_collection)),
    new_group = is.na(date_diff) | date_diff > 7,
    group_id = cumsum(new_group)
  ) %>%
  group_by(Patient, group_id) %>%
  summarise(
    Sample = first(Sample),  # Pick one sample name (can be adapted to paste0 if needed)
    Date_of_sample_collection = min(Date_of_sample_collection),
    Site = first(Site, na.rm = TRUE),
    Threshold.Coverage = first(Threshold.Coverage, na.rm = TRUE),
    Threshold.Midpoint = first(Threshold.Midpoint, na.rm = TRUE),
    Threshold.Amplitude = first(Threshold.Amplitude, na.rm = TRUE),
    Mean.Coverage = mean(Mean.Coverage, na.rm = TRUE),
    Midpoint.Coverage = mean(Midpoint.Coverage, na.rm = TRUE),
    Midpoint.normalized = mean(Midpoint.normalized, na.rm = TRUE),
    Amplitude = mean(Amplitude, na.rm = TRUE),
    Zscore.Coverage = mean(Zscore.Coverage, na.rm = TRUE),
    Zscore.Midpoint = mean(Zscore.Midpoint, na.rm = TRUE),
    Zscore.Amplitude = mean(Zscore.Amplitude, na.rm = TRUE),
    FS = mean(FS, na.rm = T),
    .groups = "drop"
  )

merged_data_grouped <- merged_data_grouped %>% select(-group_id) %>% unique()

# Identify any missing (should be only healthy controls)
missing_samples <- setdiff(
  fragdata_for_merge$Sample,
  mm_dars_small$Sample
)
cat("Samples missing from MM_DARs data (likely healthy controls):\n")
print(missing_samples)


## Add missing data for VA-07 - changed when 05 timepoint became R in script 1_0 due to relapse at this timepoint update
merged_data_grouped <- merged_data_grouped %>%
  mutate(
    Patient = case_when(
      Sample == "VA-07-05-P" ~ "VA-07",
      TRUE ~ Patient
    ),
    Date_of_sample_collection = case_when(
      Sample == "VA-07-05-P" ~ as.Date("2021-01-01"),
      TRUE ~ Date_of_sample_collection
    )
  )

## Add insert size data 
# 1) Deduplicate all_insert so there's only one Proportion.Short per Sample
all_insert_unique <- all_insert %>%
  distinct(Sample, Proportion.Short)

# 2) Left‐join onto your merged_data_grouped by Sample
merged_with_short <- merged_data_grouped %>%
  left_join(all_insert_unique, by = "Sample")

# 3) Quick check: how many got filled?
merged_with_short %>%
  summarize(
    total       = dplyr::n(),
    has_short   = sum(!is.na(Proportion.Short)),
    missing     = sum(is.na(Proportion.Short))
  )

# Correct naming
merged_with_short <- merged_with_short %>%
  mutate(
    Sample = if_else(
      str_detect(Sample, "^TGL49"),
      str_replace_all(Sample, "-", "_"),
      Sample
    )
  )
## Consolidate 
# 1) The list of metric columns want to merge
metrics <- c(
  "Mean.Coverage",       "Midpoint.Coverage",    "Midpoint.normalized",
  "Amplitude",           "Zscore.Coverage",      "Zscore.Midpoint",
  "Zscore.Amplitude",    "FS",                   "Proportion.Short"
)

# 2) A helper that turns NaN → NA then picks the first non-NA
first_non_na_numeric <- function(x) {
  # replace NaN with NA
  x <- replace(x, is.nan(x), NA_real_)
  # drop NA
  y <- x[!is.na(x)]
  if (length(y)) y[[1]] else NA_real_
}

# 3) Summarise per Sample
merged_final_metrics <- merged_with_short %>%
  group_by(Sample) %>%
  summarise(
    across(
      all_of(metrics),
      first_non_na_numeric
    ),
    .groups = "drop"
  )

# 4) Add columns back 
merged_final_metrics <- merged_final_metrics %>% left_join(merged_with_short %>% select(Patient, Sample, Date_of_sample_collection, Site, Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude)) 

# Remove healthy duplicates
# columns to blank for TGL* samples
threshold_cols <- c(
  "Threshold.Coverage",
  "Threshold.Midpoint",
  "Threshold.Amplitude",
  "Site"
)

merged_final_metrics_clean <- merged_final_metrics %>%
  mutate(
    across(
      all_of(threshold_cols),
      ~ replace(.x, str_detect(Sample, "^TGL"), NA)    # NA here is logical NA
    )
  ) %>%
  distinct()   # drop exact duplicates


# Save “Key_fragmentomics_data.csv” 
write_csv(
  merged_final_metrics_clean,
  file = "Results_Fragmentomics/Key_fragmentomics_data_updated2.csv"
)

saveRDS(
  merged_final_metrics_clean,
  file = "Results_Fragmentomics/Key_fragmentomics_data_updated2.rds"
)



### Now get range of healthy controls 
# Define the specific metrics of interest
selected_metrics <- c(
  "Mean.Coverage",
  "Midpoint.Coverage",
  "Midpoint.normalized",
  "Amplitude"
)

# Compute healthy control ranges for those metrics
healthy_dars <- mm_dars_small %>% 
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample)) %>% 
  filter(grepl("TGL49", Sample))

# Test normality 
# 1) Gather into long form
metric_vals <- healthy_dars %>%
  select(all_of(selected_metrics)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value")

# 2) Run Shapiro–Wilk by metric
normality <- metric_vals %>%
  group_by(metric) %>%
  summarise(
    shapiro_p = shapiro.test(value)$p.value,
    .groups   = "drop"
  )

## Assess cutoff
hc_ranges_selected <- healthy_dars %>%
  select(all_of(selected_metrics)) %>%
  pivot_longer(
    cols     = everything(),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  summarise(
    mean_value      = mean(value, na.rm = TRUE),
    sd_value        = sd(value,   na.rm = TRUE),
    lower_gaussian  = mean_value - 1.96 * sd_value,
    upper_gaussian  = mean_value + 1.96 * sd_value,
    lower_empirical = quantile(value, 0.025, na.rm = TRUE),
    upper_empirical = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# View the table
print(hc_ranges_selected)

# Save to CSV in your output directory
readr::write_csv(
  hc_ranges_selected,
  file.path(outdir, "HC_Ranges_Selected_MM_DARs_Metrics.csv")
)
