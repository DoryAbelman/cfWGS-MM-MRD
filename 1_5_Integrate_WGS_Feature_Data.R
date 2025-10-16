# =============================================================================
# Script: 1_5_Integrate_WGS_Feature_Data.R
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
#   source("1_5_Integrate_WGS_Feature_Data.R")
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
cna_data           <- readRDS(file.path(export_dir, "cna_data_ichorCNA.rds"))
cna_data_sequenza <- readRDS(file.path(export_dir, "cna_data_from_sequenza_400_updated.rds"))
translocation_data <- readRDS(file.path(export_dir, "translocation_data_cytoband_updated.rds"))
tumor_fraction     <- read_tsv("Oct 2024 data/tumor_fraction_cfWGS.txt")

# Load mutation MAF objects
maf_object_blood <- read.maf("~/OneDrive - University of Toronto/Project data/cfWGS project data/R outputs/combined_maf_temp_blood_Jan2025.maf")
maf_object_bm    <- read.maf("combined_maf_temp_bm_May2025.maf")


### Process CNA and translocation data 
#saveRDS(CNA_translocation, "CNA_translocation_original_Feb2025.rds")

# 1) Rename Sequenza ID to Sample so downstream stays consistent
cna_seq_renamed <- cna_data_sequenza %>% select(-Sample) %>%
  rename(Sample = Bam_clean_tmp) %>%
  # align columns to legacy CNA table if needed
  select(any_of(names(cna_data)))

# 2) Drop legacy rows that are present in Sequenza
overlap_samples <- intersect(cna_data$Sample, cna_seq_renamed$Sample)

cna_legacy_filtered <- cna_data %>%
  filter(!Sample %in% overlap_samples)

# 3) Combine CNA calls (Sequenza takes precedence)
cna_combined <- bind_rows(cna_legacy_filtered, cna_seq_renamed)

# (Optional) In case of any lingering duplicates, keep first occurrence
# cna_combined <- cna_combined %>% distinct(Sample, .keep_all = TRUE)

# 4) Join to translocation data (expects translocation_data$Sample to match Bam_clean_tmp)
CNA_translocation <- full_join(cna_combined, translocation_data, by = "Sample")

# QC
message("✅ Combined CNA rows: ", nrow(cna_combined),
        " | After translocation join: ", nrow(CNA_translocation))

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



### Now redo with the calls at the specific FISH sites 

# 1) Rename Sequenza ID to Sample so downstream stays consistent
FISH_sequenza <- readRDS(file.path(export_dir, "FISH_data_from_sequenza_400_updated.rds"))
FISH_ichor <- readRDS(file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.rds"))

FISH_sequenza <- FISH_sequenza %>% select(-Sample) %>%
  rename(Sample = Bam_clean_tmp) 

# 2) Drop legacy rows that are present in Sequenza
overlap_samples <- intersect(FISH_ichor$Sample, FISH_sequenza$Sample)

FISH_ichor_filtered <- FISH_ichor %>%
  filter(!Sample %in% overlap_samples)

# 3) Combine CNA calls (Sequenza takes precedence)
FISH_CNA_combined <- bind_rows(FISH_ichor_filtered, FISH_sequenza)

## Now join this with the metadata to ojoin with the and keep only the diagnosis ones
# First, remove the '.bam' suffix from the 'Bam' column in 'metada_df_mutation_comparison'
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Bam_clean_tmp = gsub(".bam$", "", Bam))  # Remove the '.bam' suffix

# Perform the left join
FISH_CNA_combined <- left_join(FISH_CNA_combined, 
                               metada_df_mutation_comparison, 
                               by = c("Sample" = "Bam_clean_tmp"))

## Add the tumor fraction info 
# Keep only the max Tumor_Fraction for each Bam in tumor_fraction
tumor_fraction_max <- tumor_fraction %>%
  group_by(Bam) %>%
  summarise(Tumor_Fraction = max(Tumor_fraction, na.rm = TRUE))  # Ensure to handle NA values


FISH_CNA_combined <- left_join(FISH_CNA_combined, 
                               tumor_fraction_max, 
                               by = c("Sample" = "Bam"))

# Now recalculate the lables 
gain_labels <- c("GAIN","AMP","HLAMP")
loss_labels <- c("HOMD","HETD","LOSS","CNLOH")

FISH_CNA_combined <- FISH_CNA_combined %>%
  mutate(
    # normalize case
    across(c(probe_call_amp1q, probe_call_del17p, probe_call_del1p),
           ~ toupper(as.character(.))),
    
    # 1 if gain label for amp1q; else 0 (including NA/NEUT)
    is_altered_at_probe_amp1q  = as.integer(!is.na(probe_call_amp1q)  &
                                              probe_call_amp1q  %in% gain_labels),
    
    # 1 if loss label for 17p and 1p; else 0
    is_altered_at_probe_del17p = as.integer(!is.na(probe_call_del17p) &
                                              probe_call_del17p %in% loss_labels),
    is_altered_at_probe_del1p  = as.integer(!is.na(probe_call_del1p)  &
                                              probe_call_del1p  %in% loss_labels)
  )


## Now export this 
saveRDS(FISH_CNA_combined, file = file.path(export_dir, "CNA_at_FISH_sites_combined.rds"))
write.table(FISH_CNA_combined,
            file = file.path(export_dir, "CNA_at_FISH_sites_combined.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)



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


# Bone‐marrow
temp_bm <- maf_subset@data %>%
  filter(t_depth > 10) %>%                           # only well‐supported calls
  mutate(
    Sample          = sub("\\.bam$", "", Bam),       # drop .bam
    Mutation_cDNA   = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic= paste(Chromosome,
                            Start_Position,
                            Reference_Allele,
                            Tumor_Seq_Allele2,
                            sep = "_"),
    Mutation_Type   = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")       ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")           ~ "Missense",
      Variant_Classification == "Splice_Site"                          ~ "Splice_Site",
      TRUE                                                             ~ "Other"
    )
  ) %>%
  select(
    Tumor_Sample_Barcode, Sample, Hugo_Symbol,
    Mutation_cDNA, Mutation_Genomic,
    Mutation_Type, t_depth, VAF
  ) %>%
  distinct()

# Blood
temp_blood <- maf_subset_blood@data %>%
  filter(t_depth > 10) %>%
  mutate(
    Sample          = sub("\\.bam$", "", Bam),
    Mutation_cDNA   = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic= paste(Chromosome,
                            Start_Position,
                            Reference_Allele,
                            Tumor_Seq_Allele2,
                            sep = "_"),
    Mutation_Type   = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")       ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")           ~ "Missense",
      Variant_Classification == "Splice_Site"                          ~ "Splice_Site",
      TRUE                                                             ~ "Other"
    )
  ) %>%
  select(
    Tumor_Sample_Barcode, Sample, Hugo_Symbol,
    Mutation_cDNA, Mutation_Genomic,
    Mutation_Type, t_depth, VAF
  ) %>%
  distinct()

mutation_export <- bind_rows(temp_bm, temp_blood)

saveRDS(mutation_export, file = file.path(export_dir, "mutation_export_updated2.rds"))
write.table(mutation_export, file = file.path(export_dir, "mutation_export_updated.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## Get more metrics 
temp_qc_blood <- maf_subset_blood@data %>%
  # 1) only well‐supported tumour calls
  filter(t_depth > 10) %>%
  
  # 2) your extra annotation columns
  mutate(
    Sample           = sub("\\.bam$", "", BAM_File),
    Mutation_cDNA    = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic = paste(
      Chromosome, Start_Position,
      Reference_Allele, Tumor_Seq_Allele2,
      sep = "_"
    ),
    Mutation_Type = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")   ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")     ~ "Missense",
      Variant_Classification == "Splice_Site"             ~ "Splice_Site",
      TRUE                                                ~ "Other"
    )
  ) %>%
  
  # 3) pick your key QC columns first, then grab everything else
  select(
    # core sample & variant IDs
    Tumor_Sample_Barcode, Sample, Patient, Timepoint, Sample_type,
    # gene / transcript annotation
    Hugo_Symbol, Entrez_Gene_Id, HGVSc, HGVSp, Transcript_ID, Exon_Number,
    # location & reference info
    Chromosome, Start_Position, End_Position, Strand,
    Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
    dbSNP_RS, Existing_variation, 
    # caller annotations
    Variant_Classification, Variant_Type,
    # depths & counts
    t_depth, t_ref_count, t_alt_count,
    n_depth, n_ref_count, n_alt_count,
    # allele frequencies & quality
    VAF, Score, FILTER, vcf_qual,
    # your new fields
    Mutation_cDNA, Mutation_Genomic, Mutation_Type,
    
    #—and now everything else for downstream QC
    everything()
  ) %>%
  distinct()

## Get more metrics 
temp_qc_bm <- maf_subset@data %>%
  # 1) only well‐supported tumour calls
  filter(t_depth >= 10) %>%
  
  # 2) your extra annotation columns
  mutate(
    Sample           = sub("\\.bam$", "", BAM_File),
    Mutation_cDNA    = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic = paste(
      Chromosome, Start_Position,
      Reference_Allele, Tumor_Seq_Allele2,
      sep = "_"
    ),
    Mutation_Type = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")   ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")     ~ "Missense",
      Variant_Classification == "Splice_Site"             ~ "Splice_Site",
      TRUE                                                ~ "Other"
    )
  ) %>%
  
  # 3) pick your key QC columns first, then grab everything else
  select(
    # core sample & variant IDs
    Tumor_Sample_Barcode, Sample, Patient, Timepoint, Sample_type,
    # gene / transcript annotation
    Hugo_Symbol, Entrez_Gene_Id, HGVSc, HGVSp, Transcript_ID, Exon_Number,
    # location & reference info
    Chromosome, Start_Position, End_Position, Strand,
    Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
    dbSNP_RS, Existing_variation, 
    # caller annotations
    Variant_Classification, Variant_Type,
    # depths & counts
    t_depth, t_ref_count, t_alt_count,
    n_depth, n_ref_count, n_alt_count,
    # allele frequencies & quality
    VAF, Score, FILTER, vcf_qual,
    # your new fields
    Mutation_cDNA, Mutation_Genomic, Mutation_Type,
    
    #—and now everything else for downstream QC
    everything()
  ) %>%
  distinct()

mutation_export2 <- bind_rows(temp_qc_bm, temp_qc_blood)

saveRDS(mutation_export2, file = file.path(export_dir, "mutation_export_updated_more_info2.rds"))
write.table(mutation_export2, file = file.path(export_dir, "mutation_export_updated_more_info.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# Step 1: Create a helper table with required mutation information
# Filter mutations with t_depth > 10 (all meet this criteria anyway)
filtered_mutations <- mutation_export  %>%
  filter(t_depth > 10)

# Group by Tumor_Sample_Barcode to summarize mutations for each sample
mutation_summary <- filtered_mutations %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    Mut_identified = ifelse(dplyr::n() > 0, "Y", "N"),
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
      # Tier 1 – any clear genomic hit
      IGH_MAF == 1 | IGH_CCND1 == 1 | IGH_MYC == 1 | IGH_FGFR3 == 1 |
        Mut_highest_VAF >= 0.10 ~ 1L,
      
      # Tier 2 – cfDNA SNV evidence even if TF low
      Sample_type != "BM_cells" & Mut_highest_VAF >= 0.05 ~ 1L,
      
      # Tier 3 – TF-based
      (Sample_type == "BM_cells" & Tumor_Fraction >= 0.10) |
        (Sample_type != "BM_cells" & Tumor_Fraction >= 0.045) ~ 1L,
      
      # Tier 4 – moderate TF + cytogenetics
      Tumor_Fraction >= 0.03 & Sample_type != "BM_cells" &
        (hyperdiploid == TRUE | del1p == 1 | amp1q == 1 | del17p == 1) ~ 1L,
      
      TRUE ~ 0L
    )
  )

to_logical_bin <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x > 0, FALSE))
  if (is.character(x)) return(replace_na(x %in% c("1","TRUE","T","Yes","Y"), FALSE))
  replace_na(as.logical(x), FALSE)
}

All_feature_data_logical <- All_feature_data %>%
  mutate(
    # 1) Normalize flags used in rules
    del1p        = to_logical_bin(del1p),
    amp1q        = to_logical_bin(amp1q),
    del13q       = to_logical_bin(del13q),
    del17p       = to_logical_bin(del17p),
    hyperdiploid = to_logical_bin(hyperdiploid),
    IGH_CCND1    = to_logical_bin(IGH_CCND1),
    IGH_FGFR3    = to_logical_bin(IGH_FGFR3),
    IGH_MAF      = to_logical_bin(IGH_MAF),
    IGH_MYC      = to_logical_bin(IGH_MYC),
    
    # 2) NA-safe numerics
    Tumor_Fraction  = coalesce(Tumor_Fraction, 0),
    Mut_highest_VAF = coalesce(Mut_highest_VAF, 0),
    
    # 3) Normalize sample type labels (restrict cfDNA tiers to plasma)
    Sample_type = case_when(
      Sample_type %in% c("Blood_plasma_cfDNA","cfDNA","Plasma_cfDNA") ~ "Blood_plasma_cfDNA",
      Sample_type %in% c("BM_cells","Bone_marrow_cells","BM") ~ "BM_cells",
      Sample_type %in% c("Blood_Buffy_coat","Buffy","Buffy_coat") ~ "Blood_Buffy_coat",
      TRUE ~ Sample_type
    ),
    
    Evidence_of_Disease = case_when(
      # Tier 1 – canonical drivers
      IGH_MAF | IGH_CCND1 | IGH_MYC | IGH_FGFR3 | (Mut_highest_VAF >= 0.10) ~ 1L,
      
      # Tier 2 – cfDNA SNV evidence
      Sample_type == "Blood_plasma_cfDNA" & Mut_highest_VAF >= 0.05 ~ 1L,
      
      # Tier 3 – TF-based (matrix-specific)
      (Sample_type == "BM_cells"           & Tumor_Fraction >= 0.10) |
        (Sample_type == "Blood_plasma_cfDNA" & Tumor_Fraction >= 0.05) ~ 1L,
      
      # Tier 4 – moderate cfDNA TF + cytogenetics
      Sample_type == "Blood_plasma_cfDNA" & Tumor_Fraction >= 0.03 &
        (del1p | amp1q | del17p | del13q) ~ 1L,
      
      TRUE ~ 0L
    )
  )

## Set evidence of disease to cases that did show translocations we were just not certain of them
iGV_verified <- read_excel("Jan2025_exported_data/Ig_caller_df_cfWGS_filtered_aggressive2_iGV_check.xlsm")
tmp <- iGV_verified %>% 
  filter(Looks_real > 0.7) %>% 
  select(Bam_clean_tmp) %>% 
  unique()

## If we saw a good transloation just with low read support, set evidence of disease to 1
All_feature_data_logical <- All_feature_data_logical %>%
  mutate(
    Evidence_of_Disease = if_else(
      Sample %in% tmp$Bam_clean_tmp,
      1L,                      # set to integer 1
      Evidence_of_Disease      # otherwise keep original
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

# 2) Define a helper that computes the new Evidence_of_Disease
compute_new_evidence_test <- function(df) {
  df %>% mutate(
    Evidence_new_test = case_when(
      IGH_MAF    == 1 | IGH_CCND1 == 1 | IGH_MYC  == 1 | IGH_FGFR3 == 1 |
        Mut_highest_VAF > 0.10 |
        (Sample_type != "BM_cells" & Mut_highest_VAF > 0.05)                         ~ 1,
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.03 & Tumor_Fraction <= 0.05 & Sample_type != "BM_cells" &
        (hyperdiploid == "TRUE" |
           del1p       == "1"    |
           amp1q       == "1"    |
           del13q == "1"   |
           del17p      == "1")                                                 ~ 1,
      TRUE                                                                        ~ 0
    )
  )
}

# 3) Chain them together and filter for differences
comparison <- All_feature_data %>%
  compute_new_evidence() %>%
  compute_new_evidence_test() %>%
  filter(Evidence_new != Evidence_new_test) %>%
  select(
    Sample_ID, Sample_type, Tumor_Fraction, Mut_highest_VAF,
    Evidence_new, Evidence_new_test
  )

# 4) View
print(comparison)



## Export this 
# Save All_feature_data as an RDS file
saveRDS(All_feature_data_logical, file = file.path(export_dir, "All_feature_data_Sep2025_updated2.rds"))

# Save All_feature_data as a text file with tab-separated values
write.table(All_feature_data_logical, file = file.path(export_dir, "All_feature_data_Sep2025_updated2.txt"), sep = "\t", row.names = TRUE, quote = FALSE)


### Save the CNA_Translocation file 
# Save All_feature_data as an RDS file
saveRDS(CNA_translocation, file = file.path(export_dir, "CNA_translocation_Sep2025_updated2.rds"))
#saveRDS(CNA_translocation, file = file.path(export_dir, "CNA_translocation_June2025.rds"))


# Save All_feature_data as a text file with tab-separated values
write.table(CNA_translocation, file = file.path(export_dir, "CNA_translocation_Sep2025_updated2.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
#write.table(CNA_translocation, file = file.path(export_dir, "CNA_translocation_June2025.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
