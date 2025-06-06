# ==============================================================================
# 2_A_Feature_Concordance_And_Mutation_Counts
#
# Purpose:
#   Load the full, assembled feature table and run:
#     • FISH vs WGS concordance analyses
#     • Mutation count summaries by timepoint
#
# Inputs:
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.rds
#
# Required packages:
#   tidyverse, purrr, stringr, writexl, glue
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

file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.rds")


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
      n          = n(),
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


## ------------------------------------------------------------------
## 8. SAVE RESULTS  -----------------------------------------
## ------------------------------------------------------------------
writexl::write_xlsx(concordance_tbl, "Output_tables_2025/FISH_WGS_concordance_with_cohort_updated2.xlsx")



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

write.csv(
  results2,
  "Output_tables_2025/concordance_by_cohort_and_source_and_tf_to_FISH_updated2.csv",
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







##### Now get the mutation counts 
# 1) Summary of mutations detected at baseline
baseline_summary <- dat_base %>%
  group_by(cohort) %>%
  summarise(
    n               = n(),
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
sentences <- summary_by_tp %>%
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