# ==============================================================================
# 2_3_Feature_Concordance_And_Mutation_Counts.R
#
# Purpose:
#   1. Load baseline WGS + clinical/feature table.
#   2. Compute FISH ↔ WGS concordance (overall & by ctDNA fraction) for CNAs & SVs.
#   3. Summarise baseline mutation counts in BM and cfDNA by cohort (means, ranges, tests).
#   4. Calculate Spearman correlations between mutation burden and clinical/fragmentomic features.
#   5. Fit a simple multivariable linear model on BM mutation counts.
#   6. Generate publication-ready figures:
#        • Boxplots of mutation counts and cfDNA tumor fraction by cohort
#        • Scatterplots of mutation burden vs. tumor fraction, fragment-size score, albumin
#        • Dumbbell plots of event-level concordance, sensitivity & specificity by ctDNA fraction
#
# Inputs:
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated5.rds
#   - cohort_assignment_table_updated.rds
#   - Jan2025_exported_data/mutation_export_updated.rds
#   - Jan2025_exported_data/All_feature_data_June2025.rds
#   - combined_clinical_data_updated_April2025.csv
#
# Outputs:
#   - Tables (CSV/XLSX) in Output_tables_2025/
#   - Figures in Final Tables and Figures/Baseline_concordance/
#   - R objects (RDS) for downstream use
#
# Required packages:
#   tidyverse, purrr, stringr, writexl, glue, Hmisc, broom,
#   ggpubr, patchwork, viridis, scales
# ==============================================================================

library(tidyverse)   # dplyr, tidyr, readr, etc.
library(purrr)       # for pmap_dfr
library(stringr)     # str_detect, str_to_lower, etc.
library(writexl)     # write_xlsx()
library(glue)        # for building sentence output
library(Hmisc)    # for rcorr()
library(broom)    # for tidy()
library(purrr)    # for map_df()
library(tibble)   # for tibble()

file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated6.rds")



##### PART 1: See concordance to FISH 
## 1.  PARAMETERS  ---------------------------------------------------
## ------------------------------------------------------------------
tf_cut    <- 0.05            # ≥ 0.05 → “high TF”; change if needed
baseline  <- c("Diagnosis","Baseline")   # recognise baseline labels

## ------------------------------------------------------------------
## 2.  STARTING DATA  ------------------------------------------------
## ------------------------------------------------------------------
dat <- file   # <- your tibble

## ------------------------------------------------------------------
## 3.  KEEP BASELINE SAMPLES AND ENSURE NO DUPLICATES ---------------
## ------------------------------------------------------------------
# 1. Subset to only Diagnosis or Baseline timepoints
dat_tb <- dat %>% 
  filter(
    str_to_lower(timepoint_info) %in% str_to_lower(baseline)
  )


# 2. Check for duplicate rows per patient in this subset
dup_patients <- dat_tb %>%
  count(Patient) %>%
  filter(n > 1) %>%
  pull(Patient)

# ensure Date is Date class
dat_tb <- dat_tb %>%
  mutate(Date = as.Date(Date))

# 1) Remove CA-03 timepoint 02
dat_tb2 <- dat_tb %>%
  filter(!(Patient == "CA-03" & timepoint_info == "02"))

# 2) Consolidate the two CA-02 rows
resp_CA02 <- dat_tb2 %>%
  filter(Patient == "CA-02") %>%
  # order so timepoint “01” comes before “02”
  arrange(factor(timepoint_info, levels = c("01","02"))) %>%
  summarise(across(everything(), ~{
    vals <- .
    # first non-NA in order
    first_val <- vals[which(!is.na(vals))[1]]
    # any other non-NA ≠ first_val
    other    <- vals[!is.na(vals) & vals != first_val][1]
    # if first_val is Unknown/Other but other is a “real” value, use other
    if (!is.na(first_val) &&
        first_val %in% c("Unknown","Other") &&
        !is.na(other) &&
        !other %in% c("Unknown","Other")) {
      other
    } else {
      first_val
    }
  }))

# 3) Drop all CA-02 originals
dat_tb3 <- dat_tb2 %>% filter(Patient != "CA-02")

# 4) For SPORE_0009, keep only the 2016-08-24 Baseline row
dat_tb4 <- dat_tb3 %>%
  filter(!(Patient == "SPORE_0009" &
             !(Date == as.Date("2016-08-24") & timepoint_info == "Baseline")))

# 5) Re-bind the collapsed CA-02 row
dat_tb_final <- bind_rows(dat_tb4, resp_CA02) %>%
  arrange(Patient)

dat_base <- dat_tb_final 

## Remove dup
dat_base <- dat_base %>%
  filter(!(Patient == "CA-03" & Timepoint == "02"))

# 6) Quick duplicate check
dups <- dat_base %>%
  count(Patient, timepoint_info) %>%
  filter(n > 1)

if (nrow(dups)) {
  warning("Still multiple Diagnosis/Baseline rows for: ",
          paste(unique(dups$Patient), collapse = ", "))
} else {
  message("All patients now have at most one Diagnosis/Baseline row.")
}



### Filter to the ones interested in 
## Pull from previous export 
cohort_df <- readRDS("cohort_assignment_table_updated.rds")
keep_patients <- cohort_df$Patient

## Keep only interested patients 
dat_base <- dat_base %>% filter(Patient %in% keep_patients)


# -----------------------------------------------------------
# 2.  DEFINE COHORTS  ------------
# -----------------------------------------------------------
dat_base <- dat_base %>%
  left_join(cohort_df, by = "Patient")


dat_base <- dat_base %>%
  mutate(Cohort = case_when(
    Cohort == "Frontline"     ~ "Frontline induction-transplant",
    TRUE                      ~ Cohort
  ))

dat_base$cohort <- dat_base$Cohort ## for consistency 

## Edit low confidence call
dat_base <- dat_base %>%
  mutate(across(
    starts_with("WGS_IGH_"),
    ~ if_else(grepl("ZC-02", Patient), 0L, .)
  ))

## ------------------------------------------------------------------
## 4.  MAPPINGS  -----------------------------------------------------
## ------------------------------------------------------------------
map <- tribble(
  ~fish,      ~wgs_bm,                  ~wgs_cf,                          ~type,
  "DEL_1P",   "WGS_del1p_BM_cells",     "WGS_del1p_Blood_plasma_cfDNA",   "CNA",
  "AMP_1Q",   "WGS_amp1q_BM_cells",     "WGS_amp1q_Blood_plasma_cfDNA",   "CNA",
  "DEL_17P",  "WGS_del17p_BM_cells",    "WGS_del17p_Blood_plasma_cfDNA",  "CNA",
  "T_4_14",   "WGS_IGH_FGFR3_BM_cells", "WGS_IGH_FGFR3_Blood_plasma_cfDNA","Translocation",
  "T_11_14",  "WGS_IGH_CCND1_BM_cells", "WGS_IGH_CCND1_Blood_plasma_cfDNA","Translocation",
  "T_14_16",  "WGS_IGH_MAF_BM_cells",   "WGS_IGH_MAF_Blood_plasma_cfDNA",  "Translocation"
)

## ------------------------------------------------------------------
## 5.  TIDY TO LONG FORMAT  -----------------------------------------
## ------------------------------------------------------------------
long <- map %>% 
  pmap_dfr(function(fish, wgs_bm, wgs_cf, type){
    
    dat_base %>% 
      # grab cohort here
      select(
        Patient, 
        cohort,
        Sample_Code, 
        !!sym(fish), 
        !!sym(wgs_bm), 
        !!sym(wgs_cf),
        WGS_Tumor_Fraction_Blood_plasma_cfDNA
      ) %>% 
      rename(
        fish_call = !!sym(fish),
        wgs_bm    = !!sym(wgs_bm),
        wgs_cf    = !!sym(wgs_cf),
        tf        = WGS_Tumor_Fraction_Blood_plasma_cfDNA
      ) %>% 
      mutate(
        event      = fish,
        type       = type,
        # standardise calls…
        fish_call = case_when(
          str_detect(str_to_lower(fish_call), "pos|^1$|true") ~ 1,
          str_detect(str_to_lower(fish_call), "neg|^0$|false") ~ 0,
          TRUE ~ NA_real_
        ),
        # make WGS calls logical 0/1
        across(
          c(wgs_bm, wgs_cf),
          ~ case_when(
            str_detect(str_to_lower(.), "^1$|true")  ~ 1,
            str_detect(str_to_lower(.), "^0$|false") ~ 0,
            TRUE                                   ~ NA_real_
          )
        ),
        # tumour‐fraction bucket
        tf_group = case_when(
          is.na(tf)      ~ "tf_unknown",
          tf >= tf_cut   ~ "high_tf",
          TRUE           ~ "low_tf"
        )
      ) %>% 
      pivot_longer(
        cols      = c(wgs_bm, wgs_cf),
        names_to  = "wgs_source",
        values_to = "wgs_call"
      ) %>% 
      mutate(
        wgs_source = recode(
          wgs_source,
          wgs_bm = "BM_cells",
          wgs_cf = "cfDNA"
        )
      )
    
  })



## ------------------------------------------------------------------
## 6.  FUNCTION TO SUMMARISE ACCURACY  -------------------------------
## ------------------------------------------------------------------
summarise_concord <- function(df){
  df %>% 
    filter(!is.na(fish_call) & !is.na(wgs_call)) %>%       # both performed
    summarise(
      n          = dplyr::n(),
      tp         = sum(fish_call == 1 & wgs_call == 1),
      tn         = sum(fish_call == 0 & wgs_call == 0),
      fp         = sum(fish_call == 0 & wgs_call == 1),
      fn         = sum(fish_call == 1 & wgs_call == 0),
      concord    = (tp + tn) / n,
      sens       = tp / (tp + fn),
      spec       = tn / (tn + fp)
    )
}

## ------------------------------------------------------------------
## 7.  OVERALL, & BY TF GROUP  --------------------------------------
## ------------------------------------------------------------------

# A)  per‑COHORT × TF‑GROUP  (all combinations)
conc_cohort_tf <- long %>%                     # <‑‑ the “long” object you built
  group_by(event, type, wgs_source,
           cohort,                       # Frontline induction-transplant / pre‑treated
           tf_group)   %>% summarise_concord() %>%                   # high_tf / low_tf / tf_unknow %>% summarise_concord() %>% 
  ungroup()

# B)  per‑COHORT (ignore TF bucket)  → store tf_group == "all"
conc_cohort_overall <- long %>%
  group_by(event, type, wgs_source, cohort) %>% 
  summarise_concord() %>% 
  mutate(tf_group = "all") %>%            # sentinel level
  ungroup()

# C)  per‑TF‑GROUP *irrespective* of cohort (optional – drop if you don’t need it)
conc_tf_overall <- long %>%
  group_by(event, type, wgs_source, tf_group) %>% 
  summarise_concord() %>% 
  mutate(cohort = NA) %>%                 # keep same column set
  ungroup()

# D)  bind them all and order nicely
concordance_tbl <- bind_rows(
  conc_cohort_overall,
  conc_cohort_tf,
  conc_tf_overall
) %>% 
  arrange(type, event, wgs_source, cohort, tf_group)

# have a quick look
print(concordance_tbl, n = 20)

### Get overall concordance for translocations and CNAs 
# Translocation concordance on bone‐marrow WGS vs FISH - used in manuscript
transloc_BM <- long %>% 
  filter(type == "Translocation", wgs_source == "BM_cells") %>% 
  summarise_concord()

# Translocation concordance on cfDNA WGS vs FISH
transloc_cf <- long %>% 
  filter(type == "Translocation", wgs_source == "cfDNA") %>% 
  summarise_concord()

# CNA concordance on bone‐marrow WGS vs FISH
cna_BM <- long %>% 
  filter(type == "CNA", wgs_source == "BM_cells") %>% 
  summarise_concord()

# CNA concordance on cfDNA WGS vs FISH
cna_cf <- long %>% 
  filter(type == "CNA", wgs_source == "cfDNA") %>% 
  summarise_concord()



## ------------------------------------------------------------------
## 8. SAVE RESULTS  -----------------------------------------
## ------------------------------------------------------------------
writexl::write_xlsx(concordance_tbl, "Output_tables_2025/FISH_WGS_concordance_with_cohort_updated4.xlsx")



## Put everything together for manuscript 

# ------------------------------------------------------------------
# helper : grab one concordance number
# ------------------------------------------------------------------
get_conc <- function(df,
                     aberration,          # "Translocation" / "CNA"
                     source,              # "BM_cells" / "cfDNA"
                     tf_grp    = "all",   # "all" / "high_tf" / "low_tf" / "tf_unknown"
                     cohort_val = NA      # "Frontline induction-transplant" / "Non-frontline"
){
  out <- df %>% 
    filter(
      type       == aberration,
      wgs_source == source,
      ( is.na(cohort_val) & is.na(cohort) ) | (cohort == cohort_val),
      tf_group   == tf_grp
    ) %>% 
    pull(concord)
  if(length(out)==0) NA_real_ else out[1]
}

# -----------------------------------------------------------------------------
# A) build the full grid of parameters we want
# -----------------------------------------------------------------------------
cohort_vals  <- sort(unique(na.omit(concordance_tbl$cohort)))
sources      <- unique(concordance_tbl$wgs_source)
aberrations  <- c("Translocation","CNA")
tf_groups    <- c("all","high_tf","low_tf")

param_grid <- expand_grid(
  cohort     = cohort_vals,
  wgs_source = sources,
  aberration = aberrations,
  tf_group   = tf_groups
)

# -----------------------------------------------------------------------------
# B) compute concordance for each combination
# -----------------------------------------------------------------------------
results2 <- param_grid %>%
  rowwise() %>%
  mutate(
    concordance = get_conc(
      concordance_tbl,
      aberration, wgs_source,
      tf_grp    = tf_group,
      cohort_val= cohort
    )
  ) %>%
  ungroup()

print(results2)

# -----------------------------------------------------------------------------
# 4) write everything out
# -----------------------------------------------------------------------------
if(!dir.exists("Output_tables_2025")) {
  dir.create("Output_tables_2025")
}

### This is the overall concordance to FISH
write.csv(
  results2,
  "Output_tables_2025/concordance_by_cohort_and_source_and_tf_to_FISH_updated4.csv",
  row.names = FALSE
)






#### Now see the proportion of cases with evidence of disease, calculated in earlier script
## For reporting in manuscript
evidence_summary <- dat_base %>%
  summarise(
    BM_non_na      = sum(!is.na(WGS_Evidence_of_Disease_BM_cells)),
    BM_positive    = sum(WGS_Evidence_of_Disease_BM_cells == 1, na.rm = TRUE),
    BM_pct         = BM_positive / BM_non_na * 100,
    
    Blood_non_na   = sum(!is.na(WGS_Evidence_of_Disease_Blood_plasma_cfDNA)),
    Blood_positive = sum(WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 1, na.rm = TRUE),
    Blood_pct      = Blood_positive / Blood_non_na * 100
  )







###### PART 2: See mutation overlap based on the specific base change 
mutation_data_total <- readRDS("Jan2025_exported_data/mutation_export_updated_more_info.rds")
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_June2025.rds")
combined_clinical_data_updated <- read.csv("combined_clinical_data_updated_April2025.csv")


# 1) Annotate your mutation table with clinical metadata
mut_feat <- mutation_data_total %>%
  left_join(
    All_feature_data %>%
      select(Sample, Patient, Sample_type, Tumor_Fraction, Timepoint),
    by = "Sample"
  ) %>%
  # restrict to only the two sample types of interest
  filter(Sample_type %in% c("BM_cells","Blood_plasma_cfDNA")) 

# 2) Identify only patient×timepoints that have both BM & cfDNA  
matched_pts <- combined_clinical_data_updated %>%
  filter(Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA")) %>%
  distinct(Patient, Timepoint, Sample_type) %>%
  group_by(Patient, Timepoint) %>%
  filter(n_distinct(Sample_type) == 2) %>%
  ungroup() %>%
  select(Patient, Timepoint) %>% 
  unique()

# 3) Keep only those matched cases but add doulbe negatives
mut_matched <- mut_feat %>%
  semi_join(matched_pts, by = c("Patient","Timepoint")) %>%
  inner_join(matched_pts,  by = c("Patient","Timepoint"))


# filter to in cohort 
mut_matched <- mut_matched %>% left_join(cohort_df)
mut_matched <- mut_matched %>% left_join(combined_clinical_data_updated)
mut_matched <- mut_matched %>% filter(timepoint_info %in% c("Baseline", "Diagnosis")) # get baseline
mut_matched <- mut_matched %>% filter(!is.na(Cohort))

# 4) Build per‐patient×timepoint sets of genes
mut_sets <- mut_matched %>%
  group_by(Patient, Timepoint, Cohort, Sample_type) %>%
  summarise(
    muts = list(unique(Mutation_cDNA)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Sample_type,
    values_from = muts,
    values_fill = list(muts = list(character(0)))  # in case one arm has zero
  )

tf_cut <- 0.05   # 5 %

# 1) grab only the cfDNA rows
cfDNA_tf_all <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  select(Patient, Timepoint, Tumor_Fraction) %>%
  distinct()   # in case you have duplicate

# grab cfDNA TF for every matched Pt×TP
cfDNA_tf_all <- cfDNA_tf_all %>%
  mutate(
    tf_group = case_when(
      is.na(Tumor_Fraction)          ~ "tf_unknown",
      Tumor_Fraction >  tf_cut       ~ "high_tf",
      TRUE                     ~ "low_tf"
    )
  )

# add to your mutation sets
mut_sets <- mut_sets %>%
  left_join(tf_df, by = c("Patient","Timepoint"))


## 5)  per-row concordance 
concordance_row <- mut_sets %>%
  rowwise() %>%
  mutate(
    tp  = length(intersect(BM_cells, Blood_plasma_cfDNA)),
    fn  = length(setdiff(BM_cells, Blood_plasma_cfDNA)),
    fp  = length(setdiff(Blood_plasma_cfDNA, BM_cells)),
    tn  = N_muts - length(union(BM_cells, Blood_plasma_cfDNA)),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    jaccard     = tp / length(union(BM_cells, Blood_plasma_cfDNA))
  ) %>%
  ungroup()

## 6)  global concordance **within each TF bucket** and **overall**
concordance_tf <- concordance_row %>%
  group_by(Cohort) %>%
  mutate(tf_group = replace_na(tf_group, "tf_unknown")) %>%
  group_by(tf_group, Cohort) %>%
  summarise(
    tp  = sum(tp),
    fn  = sum(fn),
    fp  = sum(fp),
    tn  = sum(tn),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    jaccard     = tp / (tp + fn + fp),   # global Jaccard
    .groups = "drop"
  )

# add an “overall” row
concordance_overall <- concordance_row %>%
  group_by(Cohort) %>%
  summarise(
    tp = sum(tp),
    fn = sum(fn),
    fp = sum(fp),
    tn = sum(tn),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    jaccard     = tp / (tp + fn + fp)
  ) %>%
  mutate(tf_group = "all")

mean_jaccard <- concordance_row %>%
  group_by(Cohort) %>%
  summarise(
    avg_jaccard = mean(jaccard, na.rm = TRUE),
    sd_jaccard  = sd(jaccard, na.rm = TRUE),
    n_samples   = dplyr::n()
  )

concordance_global <- bind_rows(concordance_tf, concordance_overall) %>%
  mutate(
    concordance = (tp + tn) / (tp + tn + fp + fn)
  ) %>%
  select(tf_group, tp, fn, fp, tn,
         sensitivity, specificity,
         jaccard, concordance, Cohort)

cat("\nGlobal concordance by tumour-fraction bucket:\n")
print(concordance_global)


## Export
# 1) CSV
write.csv(
  concordance_global,
  file = file.path(outdir, "mutation_concordance_global.csv"),
  row.names = FALSE
)


### This above does not include double negatives 


### Export the samples with mutations for checking 
# Merge to get Patient info
# 1. Create join_id by removing ".bam" from Bam in clinical data
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(join_id = str_remove(Bam, "\\.bam$"))

# 2. Left join on Sample and join_id
merged1 <- mutation_data_total 

# Merge to get cohort info
merged1 <- merged1 %>%
  left_join(cohort_df)

# Filter for Baseline or Diagnosis and patients in cohort
merged1 <- merged1 %>% 
  filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  filter(!is.na(Cohort))

write_csv(merged1, "Mutation_iGV_verification.csv")

all_bam_storage <- read_excel("All_bam_storage_locations.xlsx")

merged1 <- merged1 %>%
  mutate(Bam = paste0(Sample, ".bam"))

merged1 <- merged1 %>% 
  left_join(all_bam_storage, by = c("Bam"))

missing_path <- merged1 %>%
  filter(is.na(Path) & is.na(Location))

# This gives you which alternate name from merged1 joined to which real BAM in storage
all_bam_storage <- all_bam_storage %>%
  mutate(Bam_nosuffix = str_remove_all(Bam, "_WG|_PG"))

# Same for missing_path
missing_path <- missing_path %>%
  mutate(Bam_nosuffix = str_remove_all(Bam, "_WG|_PG"))

joined <- missing_path %>%
  left_join(
    all_bam_storage %>% select(Bam_storage = Bam, Path, Location, Bam_nosuffix, Cohort),
    by = "Bam_nosuffix"
  )

write_csv(joined, "Bam_to_unarchive1.csv")
write_csv(merged1, "Bam_to_unarchive2.csv")


## Now get most similar for still unmatched
unmatched <- joined %>% filter(is.na(Bam_storage))

# Expand unmatched against all storage BAMs, calculate distance, and get closest
library(stringdist)
similar_matches <- unmatched %>%
  rowwise() %>%
  mutate(
    best_match = {
      dists <- stringdist(Bam, all_bam_storage$Bam_noWG)
      i <- which.min(dists)
      all_bam_storage$Bam[i]
    },
    min_dist = {
      dists <- stringdist(Bam, all_bam_storage$Bam_noWG)
      min(dists)
    },
    best_path = {
      dists <- stringdist(Bam, all_bam_storage$Bam_noWG)
      i <- which.min(dists)
      all_bam_storage$Path[i]
    },
    best_location = {
      dists <- stringdist(Bam, all_bam_storage$Bam_noWG)
      i <- which.min(dists)
      all_bam_storage$Location[i]
    }
  ) %>%
  ungroup()




## Now add the location 
# Export Sample and Tumor_Sample_Barcode as CSV
export_df <- filtered[, .(Sample, Tumor_Sample_Barcode)]
fwrite(export_df, "bam_list_and_barcodes.csv")






#### PART 3: Mutation counts and other info 

##### Now get the mutation counts 
# 1) Summary of mutations detected at baseline
baseline_summary <- dat_base %>%
  group_by(cohort) %>%
  summarise(
    n               = dplyr::n(),
    mean_BM         = mean(BM_Mutation_Count,   na.rm = TRUE),
    median_BM       = median(BM_Mutation_Count, na.rm = TRUE),
    sd_BM           = sd(BM_Mutation_Count,     na.rm = TRUE),
    range_BM        = paste0(min(BM_Mutation_Count, na.rm = TRUE),
                             "–",
                             max(BM_Mutation_Count, na.rm = TRUE)),
    mean_Blood      = mean(Blood_Mutation_Count,   na.rm = TRUE),
    median_Blood    = median(Blood_Mutation_Count, na.rm = TRUE),
    sd_Blood        = sd(Blood_Mutation_Count,     na.rm = TRUE),
    range_Blood     = paste0(min(Blood_Mutation_Count, na.rm = TRUE),
                             "–",
                             max(Blood_Mutation_Count, na.rm = TRUE))
  )

print(baseline_summary)

### Add sentences and stats 
# 2) Descriptive sentences––––––––––––––––––––––––––––––––––––––––
sentences <- baseline_summary %>%
  transmute(
    sentence = glue(
      "In the {cohort} cohort (n = {n}), the mean bone marrow mutation count was ",
      "{round(mean_BM,1)} (SD = {round(sd_BM,1)}, range {range_BM}), and the mean ",
      "blood mutation count was {round(mean_Blood,1)} (SD = {round(sd_Blood,1)}, ",
      "range {range_Blood})."
    )
  ) %>%
  pull(sentence)

# Print them to console (you can copy-paste this block into Word)
cat(sentences, sep = "\n\n")


#––– 3) Statistical testing––––––––––––––––––––––––––––––––––––––––––––

# BM mutation counts: Frontline vs Non-frontline
bm_ttest <- t.test(BM_Mutation_Count ~ cohort, data = dat_base)
bm_wilcox <- wilcox.test(BM_Mutation_Count ~ cohort, data = dat_base)

# Blood mutation counts: Frontline vs Non-frontline
blood_ttest <- t.test(Blood_Mutation_Count ~ cohort, data = dat_base)
blood_wilcox <- wilcox.test(Blood_Mutation_Count ~ cohort, data = dat_base)

# Show results
bm_ttest
bm_wilcox
blood_ttest
blood_wilcox

#––– 4) Optional: extract p-values–––––––––––––––––––––––––––––––––––
pvals <- tibble(
  assay = c("BM", "BM (Wilcox)", "Blood", "Blood (Wilcox)"),
  p_value = c(bm_ttest$p.value,
              bm_wilcox$p.value,
              blood_ttest$p.value,
              blood_wilcox$p.value)
)
pvals


#### Asses other correlations with number of mutations detected at baseline 

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 1) Correlation matrix (Spearman) + p‐values across all numeric vars  
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# Select only numeric columns
num_df <- dat_base %>% select(where(is.numeric))

# rcorr returns:
#  • r : correlation matrix
#  • P : p-value matrix
#  • n : matrix of counts used in each pairwise test
rc  <- rcorr(as.matrix(num_df), type = "spearman")

r_mat <- rc$r
p_mat <- rc$P
n_mat <- rc$n

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# Flatten and include the sample count for each pair
flatten_corr <- function(r_mat, p_mat, n_mat) {
  idx <- which(lower.tri(p_mat), arr.ind = TRUE)
  tibble(
    var1     = rownames(p_mat)[idx[,1]],
    var2     = colnames(p_mat)[idx[,2]],
    rho      = r_mat[idx],
    p_val    = p_mat[idx],
    n_pairs  = n_mat[idx]
  )
}

all_corrs <- flatten_corr(r_mat, p_mat, n_mat)

## Export this 
write.csv(all_corrs %>% filter(!is.na(p_adj)), file = "Final Tables and Figures/Suplementary_Table_2_All_Feature_Correlations.csv")
# Adjust p-value for multiple hypothesis test - although not needed since exploratory 
all_corrs <- all_corrs %>%
  mutate(
    p_adj = p.adjust(p_val, method = "BH")
  )

# now significant by raw p-value
sig_raw <- all_corrs %>% filter(p_val < 0.05)

# significant by FDR
sig_fdr <- all_corrs %>% filter(p_adj < 0.05)

# view both
print(sig_raw)   # exploratory list
print(sig_fdr)   # more stringent list




#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 3) Specific Spearman tests of mutation count vs. tumor fraction + FS  
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
specific_tests <- list(
  BM_vs_BM_TF = cor.test(
    dat_base$BM_Mutation_Count,
    dat_base$WGS_Tumor_Fraction_BM_cells,
    method = "spearman",
    use = "complete.obs"
  ),
  Blood_vs_Blood_TF = cor.test(
    dat_base$Blood_Mutation_Count,
    dat_base$WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    method = "spearman",
    use = "complete.obs"
  ),
  BM_vs_FS = cor.test(
    dat_base$BM_Mutation_Count,
    dat_base$FS,
    method = "spearman",
    use = "complete.obs"
  ),
  Blood_vs_FS = cor.test(
    dat_base$Blood_Mutation_Count,
    dat_base$FS,
    method = "spearman",
    use = "complete.obs"
  )
)

# gather into one tidy table
tidy_specific <- map_df(specific_tests, tidy, .id = "comparison")
print(tidy_specific)


### Other tests 
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 1) ISS stage (ordinal) → BM mutation count
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# overall difference
kruskal.test(BM_Mutation_Count ~ ISS_STAGE, data = dat_base)
kruskal.test(Blood_Mutation_Count ~ ISS_STAGE, data = dat_base)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 2) Cytogenetic risk (binary: high vs standard)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
wilcox.test(BM_Mutation_Count ~ Cytogenetic_Risk,
            data = dat_base,
            na.action = na.exclude)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 3) Ig subtype (multilevel categorical)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
kruskal.test(BM_Mutation_Count ~ Subtype,
             data = dat_base)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 4) Spearman correlations for continuous predictors
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
library(Hmisc)
cont_vars <- c("WGS_Tumor_Fraction_BM_cells",
               "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
               "FS",
               "AGE",
               "Plasma_pct",
               "dFLC",
               "LDH")
m <- rcorr(
  as.matrix(dat_base %>% select(BM_Mutation_Count, all_of(cont_vars))),
  type = "spearman"
)
# view rho matrix and p-values
m$r
m$P


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 5) Simple multivariable model
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# log-transform counts + adjust for stage + TF + FS + age
mod <- lm(log1p(BM_Mutation_Count) ~ 
            factor(ISS_STAGE) +
            WGS_Tumor_Fraction_BM_cells +
            FS +
            AGE,
          data = dat_base)
summary(mod)

### Not significant



### Add some plots

# Create a directory for figures
if (!dir.exists("Final Tables and Figures/Baseline_concordance")) dir.create("Final Tables and Figures/Baseline_concordance")

### Updated style 
# palette
cohort_cols <- c(
  `Frontline induction-transplant` = "#3182bd",
  `Non-frontline`    = "#e6550d"
)

format_p <- function(p) {
  if (p < 0.01) {
    "<0.01"
  } else {
    sprintf("%.2f", p)
  }
}

# 1. Figure 2D — Boxplots of baseline BM vs cfDNA mutation counts
plot_df <- dat_base %>%
  select(cohort, BM_Mutation_Count, Blood_Mutation_Count) %>%
  pivot_longer(
    cols      = c(BM_Mutation_Count, Blood_Mutation_Count),
    names_to  = "Assay",
    values_to = "MutCount"
  ) %>%
  mutate(Assay = recode(Assay,
                        BM_Mutation_Count    = "Bone marrow",
                        Blood_Mutation_Count = "cfDNA"))

p1 <- ggplot(plot_df, aes(cohort, MutCount, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, colour = "black", size = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, colour = "black") +
  facet_wrap(~Assay) +
  stat_compare_means(method = "wilcox.test",
                     label    = "p.format",
                     label.y  = max(plot_df$MutCount) * 1.05) +
  scale_fill_manual(values = cohort_cols) +
  labs(
    title    = "Baseline mutation counts by cohort",
   # subtitle = "Primary vs Test, in BM and cfDNA",
    x        = "Cohort",
    y        = "Number of mutations"
  ) +
  scale_x_discrete(
    labels = c(
      "Frontline induction-transplant" = "Training",
      "Non-frontline"                  = "Test"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", size = 12,  hjust = 0.5),
    plot.subtitle = element_text(size = 12),
    strip.text    = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2F_boxplot.png", p1, width = 5, height = 4, dpi = 500)

library(ggpubr)  # for stat_compare_means

p1 <- ggplot(plot_df, aes(cohort, MutCount, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, colour = "black", size = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, colour = "black") +
  facet_wrap(~Assay) +
  # add the Wilcoxon bracket + star
  stat_compare_means(
    comparisons = list(c("Frontline induction-transplant", "Non-frontline")),
    method      = "wilcox.test",
    label       = "p.signif",     # will show * / ** / ***  
    tip.length  = 0.02,           # how far past the box ends the little ticks go
    bracket.size = 0.4            # thickness of the bracket line
  ) +
  
  scale_fill_manual(values = cohort_cols) +
  scale_x_discrete(
    labels = c(
      "Frontline induction-transplant" = "Training",
      "Non-frontline"                  = "Test"
    )
  ) +
  labs(
    title = "Baseline mutation counts by cohort",
    x     = "Cohort",
    y     = "Number of mutations"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", size = 14,  hjust = 0.5),
    strip.text     = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2F_boxplot_with_bracket.png", p1, width = 5, height = 4.25, dpi = 600)



#### Now redo but percent high tumor fraction instead 
### This is 2B
plot_df <- dat_base %>%
  select(cohort, WGS_Tumor_Fraction_Blood_plasma_cfDNA) 

p2 <- ggplot(plot_df, aes(cohort, WGS_Tumor_Fraction_Blood_plasma_cfDNA, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, colour = "black", size = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, colour = "black") +
  # add the Wilcoxon bracket + star
  stat_compare_means(
    comparisons = list(c("Frontline induction-transplant", "Non-frontline")),
    method      = "wilcox.test",
    label       = "p.signif",     # will show * / ** / ***  
    tip.length  = 0.02,           # how far past the box ends the little ticks go
    bracket.size = 0.4            # thickness of the bracket line
  ) +
   scale_y_continuous(
    labels = percent_format(accuracy = 1) #,
  #  limits = c(0, 1)       # if you want the axis to go from 0% to 100%
  ) +
  scale_fill_manual(values = cohort_cols) +
  scale_x_discrete(
    labels = c(
      "Frontline induction-transplant" = "Training",
      "Non-frontline"                  = "Test"
    )
  ) +
  labs(
    title = "cfDNA tumor fraction by cohort",
    x     = "Cohort",
    y     = "cfDNA tumor fraction (ichorCNA)"
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black")+
  theme_classic(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", size = 12,  hjust = 0.5),
    strip.text     = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2B_boxplot_with_bracket_tumor_fraction.png", p2, width = 4, height = 4, dpi = 600)


# 2. Figure 2E — BM vs cfDNA mutation counts scatter
rho_test <- cor.test(dat_base$BM_Mutation_Count,
                     dat_base$Blood_Mutation_Count,
                     method = "spearman")
rho  <- round(rho_test$estimate, 2)
pval_raw <- tf_test$p.value
pval     <- format_p(pval_raw)

p2 <- ggplot(dat_base, aes(BM_Mutation_Count, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho, "\np = ", pval),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation burden: BM vs cfDNA",
    x     = "BM mutation count",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("Figure2B_scatter_BM_vs_cfDNA.png", p2, width = 4, height = 4, dpi = 500)


# 3. Figure 2F — cfDNA mutation count vs ichorCNA tumour fraction
tf_test <- cor.test(dat_base$Blood_Mutation_Count,
                    dat_base$WGS_Tumor_Fraction_Blood_plasma_cfDNA,
                    method = "spearman")
rho_tf <- round(tf_test$estimate, 2)
p_tf   <- signif(tf_test$p.value, 2)
p_tf     <- format_p(p_tf)

p3 <- ggplot(dat_base,
             aes(WGS_Tumor_Fraction_Blood_plasma_cfDNA, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho_tf, "\np = ", p_tf),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation count vs tumour fraction",
    x     = "cfDNA tumour fraction (ichorCNA)",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("Figure2C_scatter_tf.png", p3, width = 4, height = 4, dpi = 500)


# 4. Figure 2G — cfDNA mutation count vs fragment-size score (FS)
fs_test <- cor.test(dat_base$Blood_Mutation_Count,
                    dat_base$FS, method = "spearman")
rho_fs <- round(fs_test$estimate, 2)
p_fs   <- signif(fs_test$p.value, 2)
p_fs   <- format_p(p_fs)

p4 <- ggplot(dat_base,
             aes(FS, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho_fs, "\np = ", p_fs),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation count vs fragment-size score",
    x     = "Fragment-size score (FS)",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("Figure2D_scatter_FS.png", p4, width = 4, height = 4, dpi = 500)


# 5. Figure 2H — cfDNA mutation count vs serum albumin
alb_test <- cor.test(dat_base$Blood_Mutation_Count,
                     dat_base$Albumin, method = "spearman")
rho_alb <- round(alb_test$estimate, 2)
p_alb   <- signif(alb_test$p.value, 2)
p_alb     <- format_p(p_alb)

p5 <- ggplot(dat_base,
             aes(Albumin, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho_alb, "\np = ", p_alb),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation count vs serum albumin",
    x     = "Serum albumin (g/L)",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("Figure2E_scatter_albumin.png", p5, width = 4, height = 4, dpi = 500)


## Combine  
combined <- p1 + plot_spacer() + p2 + p3 + p4 + p5 +
  plot_layout(widths  = c(1.7, .2, 1, 1, 1, 1), 
              nrow = 1, # 
              heights = c(1)) +   # bottom row slightly shorter
  plot_annotation(
    title = "Baseline Concordance and Clinical Correlates",
    theme = theme(
      plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2D_combined.png",
       combined,
       width  = 16,
       height = 4,
       dpi    = 600)



### As facet 
# 1) define which x‐vars go in which panel
panels <- tibble(
  var       = c(
    "BM_Mutation_Count",
    "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
    "FS",
    "Albumin"
  ),
  panel_lab = c(
    "BM mutation count",
    "cfDNA tumour fraction\n(ichorCNA)",
    "Fragment-size score (FS)",
    "Serum albumin (g/L)"
  )
)

# 2) pivot your dat_base into long form
df_long <- dat_base %>%
  pivot_longer(
    cols      = panels$var,
    names_to  = "var",
    values_to = "x"
  ) %>%
  left_join(panels, by="var")

# 3) compute Spearman ρ and p for each panel
corr_df <- df_long %>%
  group_by(panel_lab) %>%
  summarise(
    rho = cor(
      x, Blood_Mutation_Count,
      method = "spearman",
      use    = "complete.obs"   # <- NEW
    ),
    p   = cor.test(
      x, Blood_Mutation_Count,
      method = "spearman",
      use    = "complete.obs"   # <- cor.test() understands it too
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_text = if_else(p < 0.01, "p < 0.01", sprintf("p = %.2f", p)),
    label  = sprintf("ρ = %.2f\n%s", rho, p_text)
  )

# 4) now make the faceted scatter
p_combined <- ggplot(df_long, aes(x = x, y = Blood_Mutation_Count, colour = cohort)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  facet_wrap(~ panel_lab, scales = "free_x", nrow = 1) +
  # add the per‐panel ρ/p annotation
  geom_text(
    data = corr_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = cohort_cols,   # your existing colours
    labels = c(
      "Frontline induction-transplant" = "Training",
      "Non-frontline"                  = "Test"
    ),
    name = "Cohort"
  ) +
  labs(
    title = "Clinical correlates of cfDNA mutation burden",
    x     = NULL,
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 11),
    axis.text       = element_text(size = 9),
    legend.position = "top"
  )

# 5) save
ggsave("Final Tables and Figures/Baseline_concordance/Figure2D_facetted_scatter.png",
       p_combined,
       width  = 10,  # for a 4‐panel row
       height = 4,
       dpi    = 600)






#### Now make a dumbell plot to show concordance between BM and cfDNA
# 1) reshape, keep only frontline + TF strata
event_tf_conc <- concordance_tbl %>%
  filter(cohort == "Frontline induction-transplant",
         tf_group %in% c("high_tf","low_tf")) %>%
  transmute(
    event    = event,
    TF_group = factor(tf_group,
                      levels = c("high_tf","low_tf"),
                      labels = c("High TF","Low TF")),
    sample   = if_else(wgs_source=="BM_cells","BM","cfDNA"),
    conc     = concord * 100      # percent
  )

# 2) extract BM concordance at High TF, to define ordering
bm_high <- event_tf_conc %>%
  filter(TF_group=="High TF", sample=="BM") %>%
  arrange(conc) %>%
  pull(event)

# 3) apply that ordering to the factor
event_tf_conc <- event_tf_conc %>%
  mutate(
    event = factor(event, levels = bm_high)
  )

# Get overall
overall_conc <- concordance_tbl %>%
  filter(
    cohort   == "Frontline induction-transplant",
    tf_group == "all"
  ) %>%
  transmute(
    event  = factor(event, levels = bm_high),    # same ordering
    sample = if_else(wgs_source=="BM_cells","BM","cfDNA"),
    conc   = concord * 100                       # percent
  )

# Overall sensetivity 
sens_overall <- long %>%
  filter(cohort   == "Frontline induction-transplant") %>%
  group_by(event, wgs_source) %>%
  summarise_concord() %>%         # gives sens, etc
  ungroup() %>%
  transmute(
    event      = factor(event, levels = bm_high),
    sample     = if_else(wgs_source=="BM_cells","BM","cfDNA"),
    sens_pct   = sens * 100       # percent
  )


# 4) plot
p_tf <- ggplot(event_tf_conc,
               aes(x = conc, y = event, group = event)) +
  # grey horizontal connector
  geom_line(color="grey80", size=0.6) +
  # two TF‐group points
  geom_point(aes(colour = TF_group), size = 3) +
  # one facet per sample
  facet_wrap(~sample, nrow = 1) +
  scale_colour_viridis_d(
    option = "D", end = 0.8,
    name = "Tumour fraction"
  ) +
  scale_x_continuous(
    limits = c(0,100),
    breaks = seq(0,100, by=25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "SV/CNA concordance with FISH by tumour fraction",
    x     = "Concordance with FISH",
    y     = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    panel.spacing   = unit(1, "lines")
  )


# 5) save
ggsave("Final Tables and Figures/Baseline_concordance/Fig2B_event_concordance_by_TF_updated2.png", p_tf,
       width = 5, height = 4, dpi = 600)


## With sensetivity as red star
tf_plot_df <- bind_rows(
  # high/low TF
  event_tf_conc %>%
    rename(value = conc) %>%
    transmute(event, sample, Measure = TF_group, value),
  # overall sensitivity
  sens_overall %>%
    transmute(event, sample,
              Measure = factor("Overall sensitivity",
                               levels=c("High TF","Low TF","Overall sensitivity")),
              value = sens_pct)
) %>%
  mutate(
    Measure = factor(Measure, 
                     levels=c("High TF","Low TF","Overall sensitivity"))
  )

# 2) the plot
pretty_events <- c(
  T_4_14   = "t(4;14)",
  T_11_14  = "t(11;14)",
  T_14_16  = "t(14;16)",
  AMP_1Q   = "amp(1q)",
  DEL_17P  = "del(17p)",
  DEL_1P   = "del(1p)"
)

p_tf_sens2 <- ggplot(tf_plot_df, aes(x = value, y = event, group = event)) +
  
  # grey connector only for the TF strata
  geom_line(
    data = filter(tf_plot_df, Measure %in% c("High TF","Low TF")),
    aes(x = value, y = event, group = event),
    colour = "grey80", size = 0.6
  ) +
  
  # all three measures as points
  geom_point(aes(colour = Measure, shape = Measure),
             size = 3, stroke = 1) +
  
  # facet by BM vs cfDNA
  facet_wrap(~ sample, nrow = 1) +
  
  # Nice labels
  scale_y_discrete(
    labels = pretty_events
  ) +
  # single legend with both colour + shape
  scale_colour_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = viridis(2, end = 0.8)[1],
      "Low TF"              = viridis(2, end = 0.8)[2],
      "Overall sensitivity" = "black"
    )
  ) +
  scale_shape_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = 16,
      "Low TF"              = 16,
      "Overall sensitivity" = 8
    )
  ) +
  
  # nice breathing room at 0% and 100%
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0,100,25), "%"),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  
  labs(
    title = "SV/CNA concordance with FISH (●) and overall sensitivity",
    x     = "Percent",
    y     = NULL
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    legend.box      = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

# 3) save
ggsave(
  "Final Tables and Figures/Baseline_concordance/Fig2B_event_concordance_with_sensitivity_updated2.png",
  p_tf_sens2, width = 5.5, height = 4, dpi = 600
)


# Calculate sensitivity for CNAs and translocations separately, by WGS source
# Add type_group column and filter for events we care about
long_metrics <- long %>%
  mutate(type_group = case_when(
    grepl("^IGH_", event) ~ "Translocation",
    type == "CNA" ~ "CNA",
    TRUE ~ type
  )) %>%
  filter(!is.na(fish_call) & !is.na(wgs_call)) %>%
  filter(type_group %in% c("CNA", "Translocation"))

# Function to compute metrics for any grouping
compute_metrics <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars)), .add = FALSE) %>%
    summarise(
      tp = sum(wgs_call == 1 & fish_call == 1),
      tn = sum(wgs_call == 0 & fish_call == 0),
      fp = sum(wgs_call == 1 & fish_call == 0),
      fn = sum(wgs_call == 0 & fish_call == 1),
      sensitivity  = tp / (tp + fn),
      specificity  = tn / (tn + fp),
      concordance  = (tp + tn) / (tp + tn + fp + fn),
      ppv          = tp / (tp + fp),
      npv          = tn / (tn + fn),
      .groups = "drop"
    )
}

# Metrics per variant category (CNA vs Translocation)
metrics_by_category <- compute_metrics(long_metrics, c("type_group", "wgs_source"))

# Metrics per feature (event)
metrics_by_feature <- compute_metrics(long_metrics, c("type_group", "event", "wgs_source"))


# Export to CSVs
write_csv(metrics_by_category, "Exported_data_tables_clinical/Supp_table_wgs_fish_metrics_by_category.csv")
write_csv(metrics_by_feature, "Exported_data_tables_clinical/Supp_table_wgs_fish_metrics_by_feature.csv")


## Put together 
metrics_by_category <- metrics_by_category %>%
  mutate(
    event = dplyr::case_when(
      type_group == "CNA" ~ "All_CNAs",
      type_group == "Translocation" ~ "All_Translocations",
      TRUE ~ "All_Other"
    )
  ) %>%
  relocate(event, .after = type_group)

metrics_combined <- bind_rows(
  metrics_by_category,
  metrics_by_feature
) %>%
  arrange(type_group, factor(event, levels = c("All_CNAs","All_Translocations")), wgs_source)

## Round 
metrics_combined <- metrics_combined %>%
  mutate(
    dplyr::across(
      c(sensitivity, specificity, concordance, ppv, npv),
      ~ round(.x, 3)
    )
  )

write_csv(metrics_combined, "Exported_data_tables_clinical/Supp_table_2_WGS_vs_FISH_metrics_combined.csv")



### Now go more into depth for the cfDNA-BM concordance and plot
dir <- "Output_tables_2025"
merged_mut   <- read_rds(file.path(dir, "merged_mut.rds"))
merged_trans <- read_rds(file.path(dir, "merged_trans.rds"))
merged_CNA   <- read_rds(file.path(dir, "merged_CNA.rds"))

tf_cutoff <- 0.05   # 5% ctDNA threshold

# 1) Build a per‐event × TF‐group performance table
## First CNA
perf_by_tf <- merged_CNA %>%
  # only frontline, baseline blood samples
  filter(
    cohort                == "Frontline induction-transplant",
    timepoint_info_blood  == "Baseline"
  ) %>%
  # assign each sample to Low / High TF
  mutate(
    tf_group = case_when(
      is.na(Tumor_Fraction_blood)    ~ NA_character_,
      Tumor_Fraction_blood >= tf_cutoff ~ "High TF",
      TRUE                            ~ "Low TF"
    )
  ) %>%
  
  # pivot the five CNAs into long form
  pivot_longer(
    cols = matches("^(del1p|amp1q|del13q|del17p|hyperdiploid)_(BM|blood)$"),
    names_to      = c("event","source"),
    names_pattern = "(.*)_(BM|blood)$",
    values_to     = "call"
  ) %>%
  # normalize calls to logical
  mutate(
    call   = call == "Yes",
    source = if_else(source=="BM",   "BM", "cfDNA")
  ) %>%
  # spread BM vs blood side by side
  pivot_wider(
    id_cols    = c(Patient, event, tf_group),
    names_from = source,
    values_from= call
  ) %>%
  
  # now summarise per‐event × TF‐group
  group_by(event, tf_group) %>%
  filter(!is.na(cfDNA), !is.na(BM)) %>% # remove when either is NA
  summarise(
    n            = dplyr::n(),                                    # samples
    tp           = sum(cfDNA & BM,   na.rm=TRUE),
    tn           = sum(!cfDNA & !BM, na.rm=TRUE),
    fp           = sum(cfDNA & !BM,  na.rm=TRUE),
    fn           = sum(!cfDNA & BM,  na.rm=TRUE),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  )

# 2) Add the “All”‐TF row for each event
perf_all <- perf_by_tf %>%
  filter(!is.na(tf_group)) %>%         # drop any NA‐TF rows
  group_by(event) %>%
  summarise(
    n            = sum(n),
    tp           = sum(tp),
    tn           = sum(tn),
    fp           = sum(fp),
    fn           = sum(fn),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  ) %>%
  mutate(tf_group = "All")

perf_tf_complete <- bind_rows(perf_by_tf, perf_all) %>%
  mutate(
    tf_group = factor(tf_group, levels = c("Low TF","High TF","All"))
  )


### Now redo for translocations 
merged_trans <- merged_trans %>%
  mutate(
    Patient = str_remove(Patient_Timepoint, "_Baseline$")
  )

perf_by_tf <- merged_trans %>%
  # only frontline, baseline blood samples
  filter(
    cohort                == "Frontline induction-transplant",
    timepoint_info_blood  == "Baseline"
  ) %>%
  # assign each sample to Low / High TF
  mutate(
    tf_group = case_when(
      is.na(Tumor_Fraction_blood)    ~ NA_character_,
      Tumor_Fraction_blood >= tf_cutoff ~ "High TF",
      TRUE                            ~ "Low TF"
    )
  ) %>%
  
  # pivot the five CNAs into long form
  pivot_longer(
    cols = matches("^(IGH_MAF|IGH_MYC|IGH_CCND1|IGH_FGFR3)_(BM|blood)$"),
    names_to      = c("event","source"),
    names_pattern = "(.*)_(BM|blood)$",
    values_to     = "call"
  ) %>%
  # normalize calls to logical
  mutate(
    call   = call == "Yes",
    source = if_else(source=="BM",   "BM", "cfDNA")
  ) %>%
  # spread BM vs blood side by side
  pivot_wider(
    id_cols    = c(Patient, event, tf_group),
    names_from = source,
    values_from= call
  ) %>%
  
  # now summarise per‐event × TF‐group
  group_by(event, tf_group) %>%
  filter(!is.na(cfDNA), !is.na(BM)) %>% # remove when either is NA
  summarise(
    n            = dplyr::n(),                                    # samples
    tp           = sum(cfDNA & BM,   na.rm=TRUE),
    tn           = sum(!cfDNA & !BM, na.rm=TRUE),
    fp           = sum(cfDNA & !BM,  na.rm=TRUE),
    fn           = sum(!cfDNA & BM,  na.rm=TRUE),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  )

# 2) Add the “All”‐TF row for each event
perf_all <- perf_by_tf %>%
  filter(!is.na(tf_group)) %>%         # drop any NA‐TF rows
  group_by(event) %>%
  summarise(
    n            = sum(n),
    tp           = sum(tp),
    tn           = sum(tn),
    fp           = sum(fp),
    fn           = sum(fn),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  ) %>%
  mutate(tf_group = "All")

perf_tf_complete_trans <- bind_rows(perf_by_tf, perf_all) %>%
  mutate(
    tf_group = factor(tf_group, levels = c("Low TF","High TF","All"))
  )

## Now bind the two rows together 
perf_tf_complete <- bind_rows(perf_tf_complete, perf_tf_complete_trans)

## Add the mutations 
concordance_global$event <- "Mutations"
tmp <- concordance_global %>% filter(Cohort == "Frontline")
tmp <- tmp %>% 
  mutate(
    tf_group = case_when(
      tf_group == "high_tf" ~ "High TF",
      tf_group == "low_tf"  ~ "Low TF",
      tf_group == "all"     ~ "All",
      TRUE             ~ tf_group
    )
  )

perf_tf_complete <- bind_rows(perf_tf_complete, tmp)


## Now make plot
# 1) reshape, keep only frontline + TF strata
event_conc <- perf_tf_complete %>%
  filter(tf_group %in% c("High TF","Low TF")) %>%
  transmute(
    event    = event,
    TF_group = factor(tf_group,
                      levels = c("High TF","Low TF"),
                      labels = c("High TF","Low TF")),
    conc     = concordance * 100      # percent
  )

# 2) extract BM concordance at High TF, to define ordering
bm_high <- event_conc %>%
  filter(TF_group=="High TF") %>%
  arrange(conc) %>%
  pull(event)

# 3) apply that ordering to the factor
event_conc <- event_conc %>%
  mutate(
    event = factor(event, levels = bm_high)
  )

# Get overall
overall_conc <- perf_tf_complete %>%
  filter(tf_group %in% c("All")) %>%
  transmute(
    event  = factor(event, levels = bm_high),    # same ordering
    conc     = concordance * 100      # percent
  )

# Overall sensetivity 
sens_overall <- perf_tf_complete %>%
  filter(tf_group %in% c("All")) %>%
  transmute(
    event  = factor(event, levels = bm_high),    # same ordering
    sens_pct   = sensitivity* 100      # percent
  )
  
## With sensetivity as red star
tf_plot_df_BM <- bind_rows(
  # high/low TF
  event_conc %>%
    rename(value = conc) %>%
    transmute(event, Measure = TF_group, value),
  # overall sensitivity
  sens_overall %>%
    transmute(event,
              Measure = factor("Overall sensitivity",
                               levels=c("High TF","Low TF","Overall sensitivity")),
              value = sens_pct)
) %>%
  mutate(
    Measure = factor(Measure, 
                     levels=c("High TF","Low TF","Overall sensitivity"))
  )

# 2) the plot
pretty_events <- c(
  amp1q        = "amp(1q)",
  del13q       = "del(13q)",
  del17p       = "del(17p)",
  del1p        = "del(1p)",
  hyperdiploid = "hyperdiploid",
  IGH_CCND1    = "t(11;14) IGH-CCND1",
  IGH_FGFR3    = "t(4;14) IGH-FGFR3",
  IGH_MAF      = "t(14;16) IGH-MAF",
  IGH_MYC      = "t(8;14) IGH-MYC",
  Mutations      = "Mutations"
  
)

p_tf_sens <- ggplot(tf_plot_df_BM, aes(x = value, y = event, group = event)) +
  
  # grey connector only for the TF strata
  geom_line(
    data = filter(tf_plot_df_BM, Measure %in% c("High TF","Low TF")),
    aes(x = value, y = event, group = event),
    colour = "grey80", size = 0.6
  ) +
  
  # all three measures as points
  geom_point(aes(colour = Measure, shape = Measure),
             size = 3, stroke = 1) +

  # Nice labels
  scale_y_discrete(
    labels = pretty_events
  ) +
  # single legend with both colour + shape
  scale_colour_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = viridis(2, end = 0.8)[1],
      "Low TF"              = viridis(2, end = 0.8)[2],
      "Overall sensitivity" = "black"
    )
  ) +
  scale_shape_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = 16,
      "Low TF"              = 16,
      "Overall sensitivity" = 8
    )
  ) +
  
  # nice breathing room at 0% and 100%
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0,100,25), "%"),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  
  labs(
    title = "BM vs cfDNA SV/CNA concordance (●) and sensitivity",
    x     = "Percent",
    y     = NULL
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    legend.box      = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

# 3) save
ggsave(
  "Final Tables and Figures/Baseline_concordance/Fig2C_event_concordance_between_BM_and_cfDNA_with_sensitivity2_updated.png",
  p_tf_sens, width = 5, height = 4, dpi = 600
)





### Edit so it is 3 panel - one concordance, one sensetivity, one specificity 
# 1.  Reshape:  Concordance, Sensitivity, Specificity  ×  TF-group
perf_long <- perf_tf_complete %>%
  filter(tf_group %in% c("High TF", "Low TF")) %>%         # keep only strata
  pivot_longer(
    cols      = c(concordance, sensitivity, specificity),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Percent  = Value * 100,
    TF_group = factor(tf_group, levels = c("High TF", "Low TF"))
  )


# 2.  Event ordering (by High-TF Concordance)
event_order <- perf_long %>%
  filter(Metric == "concordance", TF_group == "High TF") %>%
  arrange(Percent) %>%
  pull(event)

perf_long <- perf_long %>%
  mutate(event = factor(event, levels = event_order))

# 4.  Plot
p_3panel <- ggplot(perf_long,
                   aes(x = Percent, y = event, group = event)) +
  # grey connector between High / Low TF
  geom_line(aes(group = interaction(event, Metric)),
            colour = "grey80", linewidth = 0.6) +
  # points for the two strata
  geom_point(aes(colour = TF_group, shape = TF_group),
             size = 3, stroke = 0.8) +
  facet_wrap(~ Metric, nrow = 1,
             labeller = labeller(Metric = c(
               concordance = "Concordance",
               sensitivity = "Sensitivity",
               specificity = "Specificity"
             ))) +
  scale_y_discrete(labels = pretty_events) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_colour_manual(
    values = c("High TF" = viridis(2, end = 0.8)[1],
               "Low TF"  = viridis(2, end = 0.8)[2]),
    name   = "Tumour fraction"
  ) +
  scale_shape_manual(
    values = c("High TF" = 16, "Low TF" = 16),
    name   = "Tumour fraction"
  ) +
  labs(
    title = "BM vs cfDNA performance by ctDNA fraction",
    x     = "Percent",
    y     = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    legend.box      = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

ggsave(
  "Final Tables and Figures/Baseline_concordance/Fig2C_BM_cfDNA_conc_sens_spec_byTF_updated2.png",
  p_3panel, width = 5.5, height = 4, dpi = 600
)


# Get overall summary
# --- 0) Define event groups (adjust if your names differ) ---------------------
cna_events           <- c("amp1q","del13q","del17p","del1p","hyperdiploid")
translocation_events <- c("IGH_CCND1","IGH_FGFR3","IGH_MAF","IGH_MYC")
mutation_events      <- c("Mutations")

# --- 1) Per-event metrics  -----------------------------------
perf_event_metrics <- perf_tf_complete %>%
  mutate(
    # ensure counts are integers and define total_n robustly
    total_n = tp + tn + fp + fn,
    prevalence = ifelse(total_n > 0, (tp + fn) / total_n, NA_real_),
    ppv       = ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_),  # precision
    npv       = ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_),
    accuracy  = ifelse(total_n > 0, (tp + tn) / total_n, NA_real_),
    f1_score  = ifelse((ppv + sensitivity) > 0, 2 * (ppv * sensitivity) / (ppv + sensitivity), NA_real_),
    fdr       = ifelse((tp + fp) > 0, fp / (tp + fp), NA_real_),
    forate    = ifelse((tn + fn) > 0, fn / (tn + fn), NA_real_),  # false omission rate
    mcc = ifelse(
      (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0,
      (tp*tn - fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)),
      NA_real_
    ),
    # Jaccard on the confusion matrix
    jaccard = ifelse((tp + fp + fn) > 0, tp / (tp + fp + fn), NA_real_)
  ) %>%
  transmute(
    category = case_when(
      event %in% cna_events ~ "CNA",
      event %in% translocation_events ~ "Translocation",
      event %in% mutation_events ~ "Mutation",
      TRUE ~ "Other"
    ),
    event,
    tf_group,
    n  = total_n,
    tp, tn, fp, fn,
    sensitivity, specificity, accuracy, ppv, npv, f1_score,
    fdr, forate, mcc, jaccard
  )

# --- 2) Category-level metrics (label as All_*) -------------------------------
perf_summary <- perf_tf_complete %>%
  mutate(
    category = case_when(
      event %in% cna_events ~ "CNA",
      event %in% translocation_events ~ "Translocation",
      event %in% mutation_events ~ "Mutation",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  group_by(category, tf_group) %>%
  summarise(
    tp = sum(tp, na.rm = TRUE),
    tn = sum(tn, na.rm = TRUE),
    fp = sum(fp, na.rm = TRUE),
    fn = sum(fn, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  mutate(
    n           = tp + tn + fp + fn,
    sensitivity = ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_),
    specificity = ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_),
    ppv         = ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_),
    npv         = ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_),
    accuracy    = ifelse(n > 0, (tp + tn) / n, NA_real_),
    f1_score    = ifelse((ppv + sensitivity) > 0, 2 * (ppv * sensitivity) / (ppv + sensitivity), NA_real_),
    fdr         = ifelse((tp + fp) > 0, fp / (tp + fp), NA_real_),
    forate      = ifelse((tn + fn) > 0, fn / (tn + fn), NA_real_),
    mcc = ifelse(
      (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0,
      (tp*tn - fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)),
      NA_real_
    ),
    jaccard     = ifelse((tp + fp + fn) > 0, tp / (tp + fp + fn), NA_real_),
    event = dplyr::case_when(
      category == "CNA" ~ "All_CNAs",
      category == "Translocation" ~ "All_Translocations",
      category == "Mutation" ~ "All_Mutations",
      TRUE ~ "All_Other"
    )
  ) %>%
  ungroup() %>%
  select(category, event, tf_group, n, tp, tn, fp, fn,
         sensitivity, specificity, accuracy, ppv, npv, f1_score,
         fdr, forate, mcc, jaccard)

# --- 3) Bind category-level + per-event into one table ------------------------
perf_combined <- bind_rows(perf_summary, perf_event_metrics) %>%
  arrange(match(category, c("Translocation","CNA","Mutation","Other")),
          match(tf_group, c("High TF","Low TF","high_tf","low_tf","BM","Blood","All","NA","NA ")),
          event)

# Round
perf_combined <- perf_combined %>%
  mutate(across(c(sensitivity, specificity, accuracy, ppv, npv, f1_score,
                  fdr, forate, mcc, jaccard),
                ~ round(.x, 3)))
# --- 4) Export ----------------------------------------------------------------
write_csv(perf_combined,
          "Final Tables and Figures/Supplentary_table_BM_vs_cfDNA_performance_by_event_and_category.csv") ## this makes second part of supp table 2
saveRDS(perf_combined,
        "Final Tables and Figures/Supplentary_table_BM_vs_cfDNA_performance_by_event_and_category.rds")




### Export additional important things 
# =====================================================================
# FINAL EXPORTS – put this block right before the script ends
# =====================================================================
outdir <- "Final Tables and Figures/Baseline_concordance/"

## 1. R-objects (RDS) --------------------------------------------------
saveRDS(long,                 file.path(outdir, "FISH_WGS_long_call_table.rds"))
saveRDS(perf_tf_complete,     file.path(outdir, "BM_cfDNA_performance_byTF.rds"))
saveRDS(perf_long,            file.path(outdir, "BM_cfDNA_performance_long_forPlot.rds"))
saveRDS(perf_summary,            file.path(outdir, "BM_cfDNA_performance_summary.rds"))
saveRDS(tf_plot_df,            file.path(outdir, "FISH_BM_cfDNA_performance.rds"))

## 2. Simple CSV / XLSX tables ----------------------------------------
readr::write_csv(baseline_summary,
                 file.path(outdir, "baseline_mutation_summary.csv"))

readr::write_csv(pvals,
                 file.path(outdir, "baseline_mutation_count_pvals.csv"))
readr::write_csv(mean_jaccard,
                 file.path(outdir, "mean_jaccard_baseline.csv"))

## 3. (Optional) one XLSX workbook with the key performance tables -----
writexl::write_xlsx(
  list(
    concordance_all          = concordance_tbl,
    concordance_byTF_cohort  = results2,
    BMcfDNA_perf_byTF        = perf_tf_complete,
    FISH_perf = tf_plot_df
  ),
  path = file.path(outdir, "SV_CNA_performance_summary.xlsx")
)

message("✓ Additional outputs written to ", outdir)

