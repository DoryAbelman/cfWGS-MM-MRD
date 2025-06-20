# =============================================================================
# Script: integrate_features.R
#
# Description:
#   Integrates copy-number (CNA), translocation, tumor-fraction, and mutation
#   features into a unified `All_feature_data` table for downstream analysis.
#   Steps:
#     1. Load and harmonize clinical metadata (sample IDs).
#     2. Read pre-exported CNA and translocation tables.
#     3. Read ichorCNA tumor-fraction data.
#     4. Load mutation calls from combined MAFs (BM and blood) via `read.maf()`.
#     5. Merge CNA + translocation with metadata and tumor fraction.
#     6. Subset MAF objects to `myeloma_genes`, classify VAF and mutation types,
#        and summarize per sample.
#     7. Join mutation summary into the feature table, fill missing values.
#     8. Define two levels of `Evidence_of_Disease` flags (high/low stringency).
#     9. Save `mutation_export` and `All_feature_data` as RDS and TSV.
#
# Inputs:
#   • combined_clinical_data_updated_Feb5_2025.csv
#   • Jan2025_exported_data/cna_data.rds
#   • Jan2025_exported_data/translocation_data.rds
#   • Oct 2024 data/tumor_fraction_cfWGS.txt
#   • combined_maf_temp_blood_Jan2025.maf
#   • combined_maf_temp_bm_Jan2025.maf
#
# Outputs:
#   • Jan2025_exported_data/mutation_export.rds/.txt
#   • Jan2025_exported_data/All_feature_data_Feb2025.rds/.txt
#
# Dependencies:
#   library(dplyr)
#   library(tidyr)
#   library(readr)
#   library(stringr)
#   library(maftools)
#
# Usage:
#   source("integrate_features.R")
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================

# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(maftools)

# Define export directory
export_dir <- "Jan2025_exported_data"

# Load clinical metadata
metada_df_mutation_comparison <- read_csv("combined_clinical_data_updated_April2025.csv") %>%
  mutate(
    Tumor_Sample_Barcode = Bam %>%
      str_remove_all("_PG|_WG") %>%
      str_replace_all("\\.filter.*|\\.ded.*|\\.recalibrate.*", ""),
    Bam_clean_tmp = str_remove(Bam, "\\.bam$")
  )

# Load CNA, translocation, and tumor fraction data
cna_data           <- readRDS(file.path(export_dir, "cna_data.rds"))
translocation_data <- readRDS(file.path(export_dir, "translocation_data.rds"))
tumor_fraction     <- read_tsv("Oct 2024 data/tumor_fraction_cfWGS.txt")

# Load mutation MAF objects
maf_object_blood <- read.maf("combined_maf_temp_blood_Jan2025.maf")
maf_object_bm    <- read.maf("combined_maf_temp_bm_Jan2025.maf")


### Process CNA and translocation data 
#saveRDS(CNA_translocation, "CNA_translocation_original_Feb2025.rds")
CNA_translocation <- 
  full_join(cna_data, translocation_data)

## Now join this with the metadata to ojoin with the and keep only the diagnosis ones
# First, remove the '.bam' suffix from the 'Bam' column in 'metada_df_mutation_comparison'
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Bam_clean_tmp = gsub(".bam$", "", Bam))  # Remove the '.bam' suffix

# Perform the left join
CNA_translocation <- left_join(CNA_translocation, 
                               metada_df_mutation_comparison, 
                               by = c("Sample" = "Bam_clean_tmp"))

## Add the tumor fraction info 
# Keep only the max Tumor_Fraction for each Bam in tumor_fraction
tumor_fraction_max <- tumor_fraction %>%
  group_by(Bam) %>%
  summarise(Tumor_Fraction = max(Tumor_fraction, na.rm = TRUE))  # Ensure to handle NA values


CNA_translocation <- left_join(CNA_translocation, 
                               tumor_fraction_max, 
                               by = c("Sample" = "Bam"))

# Filter rows where Tumor_Sample_Barcode is NA
## Can get to these later if decide to include**
samples_with_na_barcode <- CNA_translocation %>%
  filter(is.na(Tumor_Sample_Barcode)) %>% 
  unique()


##### Extract mutation data for specified genes
myeloma_genes <- c(
  "TP53",    # ~10-15%; high-risk MM
  "KRAS",    # ~20-25%; MAPK/ERK pathway
  "NRAS",    # ~20-25%; MAPK/ERK pathway
  "BRAF",    # ~5-10%; MAPK/ERK pathway
  "FAM46C",  # ~10-15%; RNA stability
  "DIS3",    # ~10-15%; RNA degradation
  "CYLD",    # ~5-10%; NF-κB regulator
  "ATM",     # ~5%; DNA damage repair
  "CCND1",   # ~15-20%; t(11;14), cyclin D1
  "MYC",     # ~15-20%; MYC translocations
  "RB1",     # ~5-10%; cell cycle control
  "TRAF3",   # ~5%; NF-κB regulator
  "IRF4",    # ~5%; plasma cell differentiation
  "FGFR3",   # ~10-15%; t(4;14), receptor tyrosine kinase
  "MMSET",   # ~10-15%; t(4;14), epigenetics
  "BCL2",    # ~15-20%; t(11;14), venetoclax target
  "IKZF1",   # ~5%; transcription regulation
  "IKZF3",   # ~5%; transcription regulation
  "CDKN2C",  # ~5-10%; cell cycle regulation
  "KDM6A",   # ~5%; epigenetics
  "SETD2",   # ~5%; histone modification
  "PTEN",    # ~5%; tumor suppressor
  "XBP1",    # ~5%; plasma cell differentiation
  "MAX",     # ~5%; MYC regulatory partner
  "SP140",   # ~5%; immune dysregulation
  "NFKBIA",  # ~5%; NF-κB inhibitor
  "NFKB2",   # ~5%; NF-κB activator
  "PRDM1",   # ~5%; plasma cell differentiation
  "EGR1",    # ~5%; early growth response
  "LTB"      # <5%; rare but part of NF-κB signaling
)


## Combine mutation_data and mutation_data_blood and export 
maf_subset <- subsetMaf(maf = maf_object_bm, genes = myeloma_genes, includeSyn = FALSE)
maf_subset_blood <- subsetMaf(maf = maf_object_blood, genes = myeloma_genes, includeSyn = FALSE)

temp_bm <- maf_subset@data %>%
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, t_depth, VAF, Bam) %>%
  mutate(Mutation_Type = case_when(
    Variant_Classification %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins") ~ "Truncating",
    Variant_Classification %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins") ~ "Missense",
    Variant_Classification == "Splice_Site" ~ "Splice_Site",
    TRUE ~ "Other"
  )) %>%
  select(-Variant_Classification) %>%
  distinct()

temp_blood <- maf_subset_blood@data %>%
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, t_depth, VAF, Bam) %>%
  mutate(Mutation_Type = case_when(
    Variant_Classification %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins") ~ "Truncating",
    Variant_Classification %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins") ~ "Missense",
    Variant_Classification == "Splice_Site" ~ "Splice_Site",
    TRUE ~ "Other"
  )) %>%
  select(-Variant_Classification) %>%
  distinct()

mutation_export <- bind_rows(temp_bm, temp_blood)

saveRDS(mutation_export, file = file.path(export_dir, "mutation_export.rds"))
write.table(mutation_export, file = file.path(export_dir, "mutation_export.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Step 1: Create a helper table with required mutation information
# Filter mutations with t_depth > 10
filtered_mutations <- mutation_export # %>%
#  filter(t_depth > 10)

# Group by Tumor_Sample_Barcode to summarize mutations for each sample
mutation_summary <- filtered_mutations %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    Mut_identified = ifelse(n() > 0, "Y", "N"),
    Mut_genes = paste(unique(Hugo_Symbol), collapse = ", "),  # List of mutated genes
    Mut_highest_VAF = max(VAF, na.rm = TRUE),                # Highest VAF among mutations
    Mut_type = paste(unique(Mutation_Type), collapse = ", ")  # List of mutation types
  )

# Step 2: Merge the mutation summary back with CNA_translocation
All_feature_data <- CNA_translocation %>%
  left_join(mutation_summary, by = "Tumor_Sample_Barcode")

# Step 3: Fill in missing values for new columns as "N" or NA, as appropriate
All_feature_data <- All_feature_data %>%
  mutate(
    Mut_identified = ifelse(is.na(Mut_identified), "N", Mut_identified),
    Mut_genes = ifelse(is.na(Mut_genes), NA, Mut_genes),
    Mut_highest_VAF = ifelse(is.na(Mut_highest_VAF), NA, Mut_highest_VAF),
    Mut_type = ifelse(is.na(Mut_type), NA, Mut_type)
  )

All_feature_data <- All_feature_data %>% select(-Bam_File)

### Version 1 - higher stringency 
# Add the Evidence_of_Disease column based on the specified conditions
All_feature_data <- All_feature_data %>%
  mutate(
    Evidence_of_Disease = case_when(
      # Condition 1: Tumor_Fraction > 10% and Sample_type is BM_cells
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells" ~ 1,
      
      # Condition 2: Tumor_Fraction > 5% and Sample_type is not BM_cells
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells" ~ 1,
      
      # Condition 3: At least one translocation or mutation with VAF > 0.1
      (IGH_MAF == 1 | IGH_CCND1 == 1 | IGH_MYC == 1 | IGH_FGFR3 == 1 | Mut_highest_VAF > 0.1) ~ 1,
      
      # If none of the above conditions are met, set to 0
      TRUE ~ 0
    )
  )

## Version 2 - less stringency on the cfDNA cases, used in final version
# Add the Evidence_of_Disease column based on the specified conditions
All_feature_data <- All_feature_data %>%
  mutate(
    Evidence_of_Disease = case_when(
      # Condition 0: any of the key translocations OR very high‐VAF (10%) OR cfDNA VAF > 5%
      IGH_MAF == 1 |
        IGH_CCND1 == 1 |
        IGH_MYC == 1 |
        IGH_FGFR3 == 1 |
        Mut_highest_VAF > 0.10 |
        (Sample_type != "BM_cells" & Mut_highest_VAF > 0.05)    ~ 1,
      
      # Primary Condition 1: high Tumor_Fraction in bone‐marrow
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells"     ~ 1,
      
      # Primary Condition 2: tumor fraction > 5% in non‐BM samples
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells"     ~ 1,
      
      # Additional Condition 3: moderate tumor fraction + supportive cytogenetics
      Tumor_Fraction > 0.03 & Tumor_Fraction <= 0.05 & Sample_type != "BM_cells" & (
        hyperdiploid == "TRUE" |
          del1p     == "1" |
          amp1q     == "1" |
          del17p    == "1"
      )                                                     ~ 1,
      
      # else
      TRUE                                                  ~ 0
    )
  )


### See difference to old version (Feb2025 version)
# 1) Define a helper that computes the old Evidence_of_Disease
compute_old_evidence <- function(df) {
  df %>% mutate(
    Evidence_old = case_when(
      IGH_MAF    == 1 | IGH_CCND1 == 1 | IGH_MYC  == 1 | IGH_FGFR3 == 1 | Mut_highest_VAF > 0.10  ~ 1,
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells"                                ~ 1,
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells"                                ~ 1,
      Tumor_Fraction > 0.03 & Tumor_Fraction <= 0.05 & Sample_type != "BM_cells" &
        (hyperdiploid == "TRUE" |
           del1p       == "1"    |
           amp1q       == "1"    |
           del17p      == "1")                                                 ~ 1,
      TRUE                                                                            ~ 0
    )
  )
}

# 2) Define a helper that computes the new Evidence_of_Disease
compute_new_evidence <- function(df) {
  df %>% mutate(
    Evidence_new = case_when(
      IGH_MAF    == 1 | IGH_CCND1 == 1 | IGH_MYC  == 1 | IGH_FGFR3 == 1 |
        Mut_highest_VAF > 0.10 |
        (Sample_type != "BM_cells" & Mut_highest_VAF > 0.05)                         ~ 1,
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.03 & Tumor_Fraction <= 0.05 & Sample_type != "BM_cells" &
        (hyperdiploid == "TRUE" |
           del1p       == "1"    |
           amp1q       == "1"    |
           del17p      == "1")                                                 ~ 1,
      TRUE                                                                        ~ 0
    )
  )
}

# 3) Chain them together and filter for differences
comparison <- All_feature_data %>%
  compute_old_evidence() %>%
  compute_new_evidence() %>%
  filter(Evidence_old != Evidence_new) %>%
  select(
    Sample_ID, Sample_type, Tumor_Fraction, Mut_highest_VAF,
    Evidence_old, Evidence_new
  )

# 4) View
print(comparison)



## Export this 
# Save All_feature_data as an RDS file
saveRDS(All_feature_data, file = file.path(export_dir, "All_feature_data_May2025.rds"))

# Save All_feature_data as a text file with tab-separated values
write.table(All_feature_data, file = file.path(export_dir, "All_feature_data_May2025.txt"), sep = "\t", row.names = TRUE, quote = FALSE)


### Save the CNA_Translocation file 
# Save All_feature_data as an RDS file
saveRDS(CNA_translocation, file = file.path(export_dir, "CNA_translocation_Feb2025.rds"))

# Save All_feature_data as a text file with tab-separated values
write.table(CNA_translocation, file = file.path(export_dir, "CNA_translocation_Feb2025.txt"), sep = "\t", row.names = TRUE, quote = FALSE)