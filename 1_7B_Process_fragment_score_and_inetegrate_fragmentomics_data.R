# ──────────────────────────────────────────────────────────────────────────────
# 01_process_fragmentomics_data.R
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
  filter(timepoint_info %in% c(group_1, group_2)) %>%
  mutate(Group = ifelse(timepoint_info %in% group_1,
                        "Monitoring", "Diagnosis/Baseline"))

# Read patient lists
patients_with_either <- read_csv("patients_with_either_cfDNA_at_baseline_and_monitoring.csv")
patients_with_both   <- read_csv("patients_with_cfDNA_at_baseline_and_monitoring.csv")

# Keep only patients that have both baseline & follow‐up
data_insert <- data_insert %>%
  filter(Patient %in% patients_with_both$x)

# Combine cfWGS + PON (dropping one duplicate PON sample if present)
pon_insert <- pon_insert %>%
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample))

all_insert <- bind_rows(
  data_insert %>% select(Sample, Patient, timepoint_info, Group, Proportion.Short),
  pon_insert %>% select(Sample, Patient = NA, timepoint_info = NA, Group, Proportion.Short)
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
  left_join(tmp.clin, by = c("Sample" = "Merge_ID")) %>%
  mutate(Cohort = "cfWGS")

pon_fs <- pon_fs %>%
  mutate(Cohort = "HBC") %>%
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample))

# Combine them
combined_fs <- bind_rows(cfWGS_fs, pon_fs)

# Join MRD data and drop unnecessary columns
combined_fs <- combined_fs %>%
  select(-Relapsed, -Num_days_to_closest_relapse) %>%
  left_join(combined_clinical,
            by = c("Patient", "Timepoint", "timepoint_info"))


# Deduplicate numeric fields (take mean) and non‐numeric (first non‐NA)
combined_fs_dedup <- combined_fs %>%
  group_by(Sample) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    across(where(~ !is.numeric(.x)), ~ first(.x)),
    .groups = "drop"
  )


### STEP 3: MERGE MM-DARs CHROMATIN ACTIVATION METRICS ###############################################
# – Read “MM_DARs_chromatin_activation_data” (generated by 1_7A script)
mm_dars <- read_csv("Results_Fragmentomics/MM_DARs_chromatin_activation_data.csv")

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

# Merge into “combined_fs_dedup” by keys (Sample, Bam, Patient, Date_of_sample_collection)
# (First, restore combined_fs_dedup into a data frame we can join onto)
fragdata_for_merge <- combined_fs_dedup

merged_data <- fragdata_for_merge %>%
  left_join(mm_dars_small,
            by = c("Sample", "Bam", "Patient", "Date_of_sample_collection"))

# Identify any missing (should be only healthy controls)
missing_samples <- setdiff(
  fragdata_for_merge$Sample,
  mm_dars_small$Sample
)
cat("Samples missing from MM_DARs data (likely healthy controls):\n")
print(missing_samples)

# Save “Key_fragmentomics_data.csv” 
write_csv(
  merged_data,
  file = "Key_fragmentomics_data_updated.csv"
)

