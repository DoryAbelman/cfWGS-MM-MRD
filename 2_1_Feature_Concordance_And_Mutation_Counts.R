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
library(scales)      # percent()


dat <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.rds")

## ------------------------------------------------------------------
## 1.  PARAMETERS  ---------------------------------------------------
## ------------------------------------------------------------------
tf_cut    <- 0.05            # ≥ 0.05 → “high TF”; change if needed
baseline  <- c("Diagnosis","Baseline")   # recognise baseline labels

## ------------------------------------------------------------------
## 2.  STARTING DATA  ------------------------------------------------
## ------------------------------------------------------------------
## The feature table was saved as "joined_clean2" in the assembly script.
## Here we load it directly from the RDS above.
# dat <- joined_clean2   # <- your tibble

## ------------------------------------------------------------------
## 3.  KEEP BASELINE SAMPLES  ---------------------------------------
## ------------------------------------------------------------------
dat_base <- dat %>% 
  filter(
    str_to_lower(timepoint_info) %in% str_to_lower(baseline)
  )

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
      select(Patient, Sample_Code, 
             !!sym(fish), !!sym(wgs_bm), !!sym(wgs_cf),
             WGS_Tumor_Fraction_Blood_plasma_cfDNA) %>% 
      rename(fish_call = !!sym(fish),
             wgs_bm    = !!sym(wgs_bm),
             wgs_cf    = !!sym(wgs_cf),
             tf        = WGS_Tumor_Fraction_Blood_plasma_cfDNA) %>% 
      mutate(
        event      = fish,
        type       = type,        ## ← carry the type column through
        # standardise calls…
        fish_call = case_when(
          str_detect(str_to_lower(fish_call), "pos|^1$|true") ~ 1,
          str_detect(str_to_lower(fish_call), "neg|^0$|false") ~ 0,
          TRUE ~ NA_real_
        ),
        across(c(wgs_bm, wgs_cf), ~ case_when(
          str_detect(str_to_lower(.), "^1$|true") ~ 1,
          str_detect(str_to_lower(.), "^0$|false") ~ 0,
          TRUE ~ NA_real_
        )),
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
        wgs_source = recode(wgs_source,
                            wgs_bm = "BM_cells",
                            wgs_cf = "cfDNA")
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
concordance_tbl <- long %>% 
  group_by(event, type, wgs_source) %>% 
  summarise_concord() %>% 
  bind_rows(
    long %>% 
      group_by(event, type, wgs_source, tf_group) %>% 
      summarise_concord() %>% 
      rename(group = tf_group)
  ) %>% 
  arrange(type, event, wgs_source, group)

## ------------------------------------------------------------------
## 8.  VIEW / SAVE RESULTS  -----------------------------------------
## ------------------------------------------------------------------
concordance_tbl
writexl::write_xlsx(concordance_tbl, "FISH_WGS_concordance.xlsx")


## Make sentence 

# helper to grab concordance
get_conc <- function(df, aberration, source, grp) {
  df %>%
    filter(
      type       == aberration,
      wgs_source == source,
      (is.na(group) & is.na(grp)) | (group == grp)
    ) %>%
    pull(concord) %>% .[1]
}

# overall (group == NA) for BM_cells
overall_trans_BM <- get_conc(concordance_tbl, "Translocation", "BM_cells", NA)
overall_CNA_BM   <- get_conc(concordance_tbl, "CNA"          , "BM_cells", NA)

# high tumor fraction for BM_cells
high_trans_BM <- get_conc(concordance_tbl, "Translocation", "BM_cells", "high_tf")
high_CNA_BM   <- get_conc(concordance_tbl, "CNA"          , "BM_cells", "high_tf")

# build sentences
sent_overall <- glue(
  "Overall concordance in BM_cells was ",
  "{percent(overall_trans_BM, accuracy = 0.1)} for translocations ",
  "and {percent(overall_CNA_BM, accuracy = 0.1)} for CNAs."
)

sent_high_tf <- glue(
  "In high–tumour–fraction samples (cfDNA TF ≥ {tf_cut * 100}%),
  concordance in BM_cells increased to ",
  "{percent(high_trans_BM, accuracy = 0.1)} for translocations ",
  "and {percent(high_CNA_BM, accuracy = 0.1)} for CNAs."
)

# print
cat(sent_overall, "\n", sent_high_tf, "\n")



### Now get for cfDNA
# overall (group == NA) for cfDNA
overall_trans_cf  <- get_conc(concordance_tbl, "Translocation", "cfDNA", NA)
overall_CNA_cf    <- get_conc(concordance_tbl, "CNA"          , "cfDNA", NA)

# high tumor fraction for cfDNA
high_trans_cf     <- get_conc(concordance_tbl, "Translocation", "cfDNA", "high_tf")
high_CNA_cf       <- get_conc(concordance_tbl, "CNA"          , "cfDNA", "high_tf")

# build sentences
sent_overall_cf <- glue(
  "Overall concordance in blood cfDNA was ",
  "{percent(overall_trans_cf, accuracy = 0.1)} for translocations ",
  "and {percent(overall_CNA_cf, accuracy = 0.1)} for CNAs."
)

sent_high_tf_cf <- glue(
  "In high–tumour–fraction samples (cfDNA TF ≥ {tf_cut * 100}%),
  concordance in blood cfDNA increased to ",
  "{percent(high_trans_cf, accuracy = 0.1)} for translocations ",
  "and {percent(high_CNA_cf, accuracy = 0.1)} for CNAs."
)

# print
cat(sent_overall_cf, "\n", sent_high_tf_cf, "\n")

## Export 
write.csv(concordance_tbl, "FISH_WGS_concordance.csv", row.names = FALSE)







##### Now get the mutation counts 
# 1) Summary by each Timepoint (e.g. "01", "02", …) for both BM and Blood mutation counts
summary_by_tp <- joined_consolidated %>%
  group_by(Timepoint) %>%
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

print(summary_by_tp)


# 2) Overall summary restricted to baseline/diagnosis timepoints
baseline_summary <- joined_consolidated %>%
  filter(timepoint_info == "Diagnosis" | timepoint_info == "Baseline") %>%
  summarise(
    n_baseline            = n(),
    mean_BM_baseline      = mean(BM_Mutation_Count,   na.rm = TRUE),
    sd_BM_baseline        = sd(BM_Mutation_Count,     na.rm = TRUE),
    median_BM_baseline    = median(BM_Mutation_Count, na.rm = TRUE),
    range_BM_baseline     = paste0(min(BM_Mutation_Count, na.rm = TRUE),
                                   "–",
                                   max(BM_Mutation_Count, na.rm = TRUE)),
    mean_Blood_baseline   = mean(Blood_Mutation_Count,   na.rm = TRUE),
    sd_Blood_baseline     = sd(Blood_Mutation_Count,     na.rm = TRUE),
    median_Blood_baseline = median(Blood_Mutation_Count, na.rm = TRUE),
    range_Blood_baseline  = paste0(min(Blood_Mutation_Count, na.rm = TRUE),
                                   "–",
                                   max(Blood_Mutation_Count, na.rm = TRUE))
  )

print(baseline_summary)