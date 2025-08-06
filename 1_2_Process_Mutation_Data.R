# =============================================================================
# Script: 1_2_Process_Mutation_Data.R
#
# Description:
#   This script ingests and harmonizes mutation calls from MAF files for both
#   bone-marrow (BM) and peripheral-blood cfDNA (PB cfDNA), annotates them with
#   clinical metadata, computes variant allele frequencies (VAF), and produces
#   filtered RDS outputs and summary plots.  Key steps:
#     1. Define the list of myeloma-relevant genes (`myeloma_genes`).
#     2. Read all BM‐MAF files from `maf_directory`, parse with `read_tsv()`, bind
#        into `combined_maf`.
#     3. Load clinical metadata from
#        “combined_clinical_data_updated_Feb5_2025.csv” into
#        `metada_df_mutation_comparison` and join on `Tumor_Sample_Barcode`.
#     4. Compute VAF (`t_alt_count/(t_ref_count+t_alt_count)`) and standardize
#        sample IDs.
#     5. Filter BM data to diagnosis/baseline and save `combined_maf_bm_dx`.
#     6. Repeat steps 2–5 for PB cfDNA MAFs into `combined_maf_blood`.
#     7. Export combined objects as RDS:
#        • combined_maf_bm_all_muts.rds
#        • combined_maf_bm_dx.rds
#        • combined_maf_blood_all_muts.rds
#     8. Generate VAF density‐ridge plots (BM and PB), histograms of `t_depth`,
#        and save figures.
#
# Inputs:
#   • maf_directory: path to “*.maf” BM and PB directories
#   • combined_clinical_data_updated_Feb5_2025.csv
#
# Outputs:
#   • RDS: combined_maf_bm_all_muts.rds
#   • RDS: combined_maf_bm_dx.rds
#   • RDS: combined_maf_blood_all_muts.rds
#   • Figures: VAF ridgeplots & t_depth histograms
#
# Dependencies:
#   library(maftools)
#   library(dplyr)
#   library(tidyr)
#   library(readr)
#   library(ggplot2)
#   library(ggridges)
#   library(viridis)
#   library(scales)
#   library(stringr)
#   library(purrr)
#   library(ComplexHeatmap)  # if downstream heatmaps are created
#   library(circlize)        # if downstream heatmaps are created
#
# Usage:
#   source("1_2_Process_Mutation_Data.R")
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================

# Load required libraries
library(maftools)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggridges)
library(viridis)
library(scales)
library(stringr)
library(purrr)


### Define the mutation gene list (used for both BM and blood)
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


#### FIRST Load in BM Data

# Define the directory containing the MAF files
maf_directory <- "~/OneDrive - University of Toronto/Project data/cfWGS project data/MAF files/BM/"

# List all MAF files in the directory
maf_files <- list.files(maf_directory, pattern = "\\.maf$", full.names = TRUE)

# Read each MAF file into a dataframe and correct column types
dfs <- lapply(maf_files, function(file) {
  df <- read_tsv(file, comment = "#", col_types = cols(
    Hugo_Symbol = col_character(),
    Entrez_Gene_Id = col_integer(),
    Center = col_character(),
    NCBI_Build = col_character(),
    Chromosome = col_character(),
    Start_Position = col_integer(),
    End_Position = col_integer(),
    Strand = col_character(),
    Variant_Classification = col_character(),
    Variant_Type = col_character(),
    Reference_Allele = col_character(),
    Tumor_Seq_Allele1 = col_character(),
    Tumor_Seq_Allele2 = col_character(),
    dbSNP_RS = col_character(),
    dbSNP_Val_Status = col_character(),
    Tumor_Sample_Barcode = col_character(),
    Matched_Norm_Sample_Barcode = col_character(),
    Match_Norm_Seq_Allele1 = col_character(),
    Match_Norm_Seq_Allele2 = col_character(),
    Tumor_Validation_Allele1 = col_character(),
    Tumor_Validation_Allele2 = col_character(),
    Match_Norm_Validation_Allele1 = col_character(),
    Match_Norm_Validation_Allele2 = col_character(),
    Verification_Status = col_character(),
    Validation_Status = col_character(),
    Mutation_Status = col_character(),
    Sequencing_Phase = col_character(),
    Sequence_Source = col_character(),
    Validation_Method = col_character(),
    Score = col_double(),
    BAM_File = col_character(),
    Sequencer = col_character(),
    Tumor_Sample_UUID = col_character(),
    Matched_Norm_Sample_UUID = col_character(),
    HGVSc = col_character(),
    HGVSp = col_character(),
    HGVSp_Short = col_character(),
    Transcript_ID = col_character(),
    Exon_Number = col_character(),
    t_depth = col_integer(),
    t_ref_count = col_integer(),
    t_alt_count = col_integer(),
    n_depth = col_integer(),
    n_ref_count = col_integer(),
    n_alt_count = col_integer(),
    all_effects = col_character(),
    Allele = col_character(),
    Gene = col_character(),
    Feature = col_character(),
    Feature_type = col_character(),
    Consequence = col_character(),
    cDNA_position = col_character(),
    CDS_position = col_character(),
    Protein_position = col_character(),
    Amino_acids = col_character(),
    Codons = col_character(),
    Existing_variation = col_character(),
    ALLELE_NUM = col_integer(),
    DISTANCE = col_double(),
    STRAND_VEP = col_character(),
    SYMBOL = col_character(),
    SYMBOL_SOURCE = col_character(),
    HGNC_ID = col_character(),
    BIOTYPE = col_character(),
    CANONICAL = col_character(),
    CCDS = col_character(),
    ENSP = col_character(),
    SWISSPROT = col_character(),
    TREMBL = col_character(),
    UNIPARC = col_character(),
    RefSeq = col_character(),
    SIFT = col_character(),
    PolyPhen = col_character(),
    EXON = col_character(),
    INTRON = col_character(),
    DOMAINS = col_character(),
    AF = col_double(),
    AFR_AF = col_double(),
    AMR_AF = col_double(),
    ASN_AF = col_double(),
    EAS_AF = col_double(),
    EUR_AF = col_double(),
    SAS_AF = col_double(),
    AA_AF = col_double(),
    EA_AF = col_double(),
    CLIN_SIG = col_character(),
    SOMATIC = col_character(),
    PUBMED = col_character(),
    MOTIF_NAME = col_character(),
    MOTIF_POS = col_character(),
    HIGH_INF_POS = col_character(),
    MOTIF_SCORE_CHANGE = col_double(),
    IMPACT = col_character(),
    PICK = col_character(),
    VARIANT_CLASS = col_character(),
    TSL = col_character(),
    HGVS_OFFSET = col_character(),
    PHENO = col_character(),
    MINIMISED = col_character(),
    ExAC_AF = col_double(),
    ExAC_AF_AFR = col_double(),
    ExAC_AF_AMR = col_double(),
    ExAC_AF_EAS = col_double(),
    ExAC_AF_FIN = col_double(),
    ExAC_AF_NFE = col_double(),
    ExAC_AF_OTH = col_double(),
    ExAC_AF_SAS = col_double(),
    GENE_PHENO = col_character(),
    FILTER = col_character(),
    flanking_bps = col_character(),
    vcf_id = col_character(),
    vcf_qual = col_double(),
    ExAC_AF_Adj = col_double(),
    ExAC_AC_AN_Adj = col_character(),
    ExAC_AC_AN = col_character(),
    ExAC_AC_AN_AFR = col_character(),
    ExAC_AC_AN_AMR = col_character(),
    ExAC_AC_AN_EAS = col_character(),
    ExAC_AC_AN_FIN = col_character(),
    ExAC_AC_AN_NFE = col_character(),
    ExAC_AC_AN_OTH = col_character(),
    ExAC_AC_AN_SAS = col_character(),
    ExAC_FILTER = col_character(),
    gnomAD_AF = col_double(),
    gnomAD_AFR_AF = col_double(),
    gnomAD_AMR_AF = col_double(),
    gnomAD_ASJ_AF = col_double(),
    gnomAD_EAS_AF = col_double(),
    gnomAD_FIN_AF = col_double(),
    gnomAD_NFE_AF = col_double(),
    gnomAD_OTH_AF = col_double(),
    gnomAD_SAS_AF = col_double(),
    vcf_pos = col_integer()
  ))
  return(df)
})

# Combine all dataframes into one
combined_maf <- bind_rows(dfs)
rm(dfs)

# Load in the patient info 
metada_df_mutation_comparison <- read_csv("combined_clinical_data_updated_April2025.csv")

# Add a Tumor_Sample_Barcode column to metada_df_mutation_comparison
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Tumor_Sample_Barcode = Bam %>%
           # Remove _PG or _WG
           str_remove_all("_PG|_WG") %>%
           # Remove anything after ".filter", ".ded", or ".recalibrate"
           str_replace_all("\\.filter.*|\\.ded.*|\\.recalibrate.*", ""))

# Add VAF column
combined_maf <- combined_maf %>%
  mutate(VAF = t_alt_count / (t_ref_count + t_alt_count))

# Add the Bam column to combined_maf
combined_maf <- combined_maf %>%
  mutate(Bam = paste0(Tumor_Sample_Barcode, ".filter.deduped.recalibrated.bam"))

# Modify the specific Tumor_Sample_Barcode in combined_maf with error
combined_maf <- combined_maf %>%
  mutate(Tumor_Sample_Barcode = ifelse(Tumor_Sample_Barcode == "TFRIM4_0189_Bm_P_ZC-02", 
                                       paste0(Tumor_Sample_Barcode, "-01-O-DNA"), 
                                       Tumor_Sample_Barcode))

# Join with the metadata dataframe
combined_maf <- left_join(combined_maf %>% select(-Bam), metada_df_mutation_comparison, by = "Tumor_Sample_Barcode")


## First filter to just diagnosis 
combined_maf_bm_dx <- combined_maf %>% filter(timepoint_info %in% c("Diagnosis", "Baseline", "Relapse", "Progression"))
rm(combined_maf)




#### Next load in PB cfDNA Data
# Read each MAF file into a dataframe and correct column types
dfs <- lapply(maf_files, function(file) {
  df <- read_tsv(file, comment = "#", col_types = cols(
    Hugo_Symbol = col_character(),
    Entrez_Gene_Id = col_integer(),
    Center = col_character(),
    NCBI_Build = col_character(),
    Chromosome = col_character(),
    Start_Position = col_integer(),
    End_Position = col_integer(),
    Strand = col_character(),
    Variant_Classification = col_character(),
    Variant_Type = col_character(),
    Reference_Allele = col_character(),
    Tumor_Seq_Allele1 = col_character(),
    Tumor_Seq_Allele2 = col_character(),
    dbSNP_RS = col_character(),
    dbSNP_Val_Status = col_character(),
    Tumor_Sample_Barcode = col_character(),
    Matched_Norm_Sample_Barcode = col_character(),
    Match_Norm_Seq_Allele1 = col_character(),
    Match_Norm_Seq_Allele2 = col_character(),
    Tumor_Validation_Allele1 = col_character(),
    Tumor_Validation_Allele2 = col_character(),
    Match_Norm_Validation_Allele1 = col_character(),
    Match_Norm_Validation_Allele2 = col_character(),
    Verification_Status = col_character(),
    Validation_Status = col_character(),
    Mutation_Status = col_character(),
    Sequencing_Phase = col_character(),
    Sequence_Source = col_character(),
    Validation_Method = col_character(),
    Score = col_double(),
    BAM_File = col_character(),
    Sequencer = col_character(),
    Tumor_Sample_UUID = col_character(),
    Matched_Norm_Sample_UUID = col_character(),
    HGVSc = col_character(),
    HGVSp = col_character(),
    HGVSp_Short = col_character(),
    Transcript_ID = col_character(),
    Exon_Number = col_character(),
    t_depth = col_integer(),
    t_ref_count = col_integer(),
    t_alt_count = col_integer(),
    n_depth = col_integer(),
    n_ref_count = col_integer(),
    n_alt_count = col_integer(),
    all_effects = col_character(),
    Allele = col_character(),
    Gene = col_character(),
    Feature = col_character(),
    Feature_type = col_character(),
    Consequence = col_character(),
    cDNA_position = col_character(),
    CDS_position = col_character(),
    Protein_position = col_character(),
    Amino_acids = col_character(),
    Codons = col_character(),
    Existing_variation = col_character(),
    ALLELE_NUM = col_integer(),
    DISTANCE = col_double(),
    STRAND_VEP = col_character(),
    SYMBOL = col_character(),
    SYMBOL_SOURCE = col_character(),
    HGNC_ID = col_character(),
    BIOTYPE = col_character(),
    CANONICAL = col_character(),
    CCDS = col_character(),
    ENSP = col_character(),
    SWISSPROT = col_character(),
    TREMBL = col_character(),
    UNIPARC = col_character(),
    RefSeq = col_character(),
    SIFT = col_character(),
    PolyPhen = col_character(),
    EXON = col_character(),
    INTRON = col_character(),
    DOMAINS = col_character(),
    AF = col_double(),
    AFR_AF = col_double(),
    AMR_AF = col_double(),
    ASN_AF = col_double(),
    EAS_AF = col_double(),
    EUR_AF = col_double(),
    SAS_AF = col_double(),
    AA_AF = col_double(),
    EA_AF = col_double(),
    CLIN_SIG = col_character(),
    SOMATIC = col_character(),
    PUBMED = col_character(),
    MOTIF_NAME = col_character(),
    MOTIF_POS = col_character(),
    HIGH_INF_POS = col_character(),
    MOTIF_SCORE_CHANGE = col_double(),
    IMPACT = col_character(),
    PICK = col_character(),
    VARIANT_CLASS = col_character(),
    TSL = col_character(),
    HGVS_OFFSET = col_character(),
    PHENO = col_character(),
    MINIMISED = col_character(),
    ExAC_AF = col_double(),
    ExAC_AF_AFR = col_double(),
    ExAC_AF_AMR = col_double(),
    ExAC_AF_EAS = col_double(),
    ExAC_AF_FIN = col_double(),
    ExAC_AF_NFE = col_double(),
    ExAC_AF_OTH = col_double(),
    ExAC_AF_SAS = col_double(),
    GENE_PHENO = col_character(),
    FILTER = col_character(),
    flanking_bps = col_character(),
    vcf_id = col_character(),
    vcf_qual = col_double(),
    ExAC_AF_Adj = col_double(),
    ExAC_AC_AN_Adj = col_character(),
    ExAC_AC_AN = col_character(),
    ExAC_AC_AN_AFR = col_character(),
    ExAC_AC_AN_AMR = col_character(),
    ExAC_AC_AN_EAS = col_character(),
    ExAC_AC_AN_FIN = col_character(),
    ExAC_AC_AN_NFE = col_character(),
    ExAC_AC_AN_OTH = col_character(),
    ExAC_AC_AN_SAS = col_character(),
    ExAC_FILTER = col_character(),
    gnomAD_AF = col_double(),
    gnomAD_AFR_AF = col_double(),
    gnomAD_AMR_AF = col_double(),
    gnomAD_ASJ_AF = col_double(),
    gnomAD_EAS_AF = col_double(),
    gnomAD_FIN_AF = col_double(),
    gnomAD_NFE_AF = col_double(),
    gnomAD_OTH_AF = col_double(),
    gnomAD_SAS_AF = col_double(),
    vcf_pos = col_integer()
  ))
  return(df)
})

# Combine all dataframes into one
combined_maf_blood <- bind_rows(dfs)

rm(dfs)

# Add VAF column
combined_maf_blood <- combined_maf_blood %>%
  mutate(VAF = t_alt_count / (t_ref_count + t_alt_count))


# Join with the metadata dataframe

combined_maf_blood <- left_join(combined_maf_blood, metada_df_mutation_comparison, by = "Tumor_Sample_Barcode")

# Filter rows where timepoint_info is NA and get unique Tumor_Sample_Barcode values
unique_barcodes_na_timepoint <- combined_maf_blood %>%
  dplyr::filter(is.na(timepoint_info)) %>%
  distinct(Tumor_Sample_Barcode)

# Display the unique barcodes
unique_barcodes_na_timepoint

# Save combined_maf_blood as an RDS file
saveRDS(combined_maf_blood, file = "combined_maf_blood_all_muts.rds")

# Save combined_maf_bm_dx as an RDS file
saveRDS(combined_maf_bm_dx, file = "combined_maf_bm_dx.rds")
saveRDS(combined_maf, file = "combined_maf_bm_all_muts.rds")






########### Plot histogram of mutations 
## VAF plot 

# BM samples
vaf_plot <- ggplot(combined_maf_bm_dx %>% filter(!is.na(VAF)) %>%
                     filter(timepoint_info %in% c("Diagnosis", "Baseline")), aes(x = VAF, y = Patient)) +
  geom_density_ridges(scale = 2) +
  scale_fill_viridis_d(alpha = 0.9) +
  theme_classic() +  # Using a minimal theme as an example
  labs(
    title = "VAF Distribution for Each Patient",
    subtitle = "Density ridgeline plots of VAFs of BM cells at diagnosis"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis text
  scale_x_continuous(breaks = scales::breaks_width(0.05), limits = c(0, 1))  # Adjust the x-axis as needed

vaf_plot
ggsave("Vaf_plot_BM_cell_dx_updated_3.png", plot = vaf_plot, width = 15, height = 12, dpi = 500)

## Order 
# Step 1: Compute the mean VAF for each patient
mean_vaf_per_patient <- combined_maf_bm_dx %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  filter(!is.na(VAF)) %>%
  group_by(Patient) %>%
  summarize(mean_vaf = mean(VAF, na.rm = TRUE))

# Step 2: Reorder the patients based on the mean VAF
combined_maf_bm_dx <- combined_maf_bm_dx %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline", "Progression"))   %>%
  filter(!is.na(VAF)) %>%
  mutate(Patient = factor(Patient, levels = mean_vaf_per_patient$Patient[order(mean_vaf_per_patient$mean_vaf)]))

# Step 3: Create the ridgeline plot with reordered patients
vaf_plot <- ggplot(combined_maf_bm_dx, aes(x = VAF, y = Patient)) +
  geom_density_ridges(scale = 2) +
  scale_fill_viridis_d(alpha = 0.9) +
  theme_classic() +
  labs(
    title = "VAF Distribution for Each Patient",
    subtitle = "Density ridgeline plots of VAFs of BM cells at diagnosis, ordered by mean VAF"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = scales::breaks_width(0.05), limits = c(0, 1))

# Display the plot
vaf_plot

# Save the plot
ggsave("Vaf_plot_BM_cell_dx_ordered_updated_4.png", plot = vaf_plot, width = 6, height = 12, dpi = 500)




# Blood samples
combined_maf_blood <- combined_maf_blood %>% dplyr::filter(timepoint_info %in% c("Diagnosis", "Baseline", "Progression", "Relapse")) %>% dplyr::filter(Sample_type == "Blood_plasma_cfDNA")
vaf_plot <- ggplot(combined_maf_blood %>% dplyr::filter(!is.na(VAF)) %>% dplyr::filter(timepoint_info %in% c("Diagnosis", "Baseline")), aes(x = VAF, y = Patient)) +
  geom_density_ridges(scale = 2) +
  scale_fill_viridis_d(alpha = 0.9) +
  theme_classic() +  # Using a minimal theme as an example
  labs(
    title = "VAF Distribution for Each Patient",
    subtitle = "Density ridgeline plots of VAFs of PB cfDNA samples"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis text
  scale_x_continuous(breaks = scales::breaks_width(0.05), limits = c(0, 1))  # Adjust the x-axis as needed

ggsave("Vaf_plot_Blood_cfDNA_dx_updated_4.png", plot = vaf_plot, width = 15, height = 22, dpi = 500)


## Reorder blood VAFs 
# Step 1: Compute the mean VAF for each patient
blood_dx_maf <- combined_maf_blood 

mean_vaf_per_patient <- blood_dx_maf %>%
  filter(!is.na(VAF)) %>%
  group_by(Patient) %>%
  summarize(mean_vaf = mean(VAF, na.rm = TRUE))

# Step 2: Reorder the patients based on the mean VAF
blood_dx_maf <- blood_dx_maf %>%
  filter(!is.na(VAF)) %>%
  mutate(Patient = factor(Patient, levels = mean_vaf_per_patient$Patient[order(mean_vaf_per_patient$mean_vaf)]))

# Step 3: Create the ridgeline plot with reordered patients
vaf_plot <- ggplot(blood_dx_maf %>% dplyr::filter(timepoint_info %in% c("Diagnosis", "Baseline")), aes(x = VAF, y = Patient)) +
  geom_density_ridges(scale = 2) +
  scale_fill_viridis_d(alpha = 0.9) +
  theme_classic() +
  labs(
    title = "VAF Distribution for Each Patient",
    subtitle = "Density ridgeline plots of VAFs of PB cfDNA, ordered by mean VAF"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = scales::breaks_width(0.05), limits = c(0, 0.35))

# Display the plot
vaf_plot

# Save the plot
ggsave("Vaf_plot_blood_dx_ordered_updated_4.png", plot = vaf_plot, width = 6, height = 12, dpi = 500)


# Filter the data to include only t_depth values greater than 1
filtered_data <- combined_maf_bm_dx[combined_maf_bm_dx$t_depth >= 1, ]

# Create the histogram with ggplot2
histogram <- ggplot(combined_maf_bm_dx, aes(x = t_depth)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of t_depth values in BM",
       x = "t_depth",
       y = "Frequency") +
  theme_classic() +
  xlim(1, 200) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggsave("Histogram_BM_muts.png", plot = histogram, width = 6, height = 6, dpi = 500)

# Create the histogram with ggplot2
histogram <- ggplot(combined_maf_blood, aes(x = t_depth)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of t_depth values in Blood cfDNA",
       x = "t_depth",
       y = "Frequency") +
  theme_classic() +
  xlim(1, 200) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggsave("Histogram_blood_muts.png", plot = histogram, width = 6, height = 6, dpi = 500)


### Cleaning up 
rm(combined_maf)
rm(blood_dx_maf)
rm(vaf_plot)


# Write combined_maf to a temporary MAF file 
write.table(as.data.frame(combined_maf_blood), "combined_maf_temp_blood_Jan2025.maf", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(combined_maf_bm_dx), "combined_maf_temp_bm_Jan2025.maf", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(as.data.frame(combined_maf_bm_dx), "combined_maf_temp_bm_May2025.maf", sep = "\t", quote = FALSE, row.names = FALSE)


#### Below here is optional

# Read the MAF file using read.maf
maf_object_blood <- read.maf(maf = "combined_maf_temp_blood_Jan2025.maf")
maf_object_bm <- read.maf(maf = "combined_maf_temp_bm_Jan2025.maf")

#### Transform for heatmaps (optional, redone in heatmap script)
maf_subset <- subsetMaf(maf = maf_object_bm, genes = myeloma_genes, includeSyn = FALSE)
mutation_data <- maf_subset@data %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
  mutate(Mutation_Type = case_when(
    Variant_Classification %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins") ~ "Truncating",
    Variant_Classification %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins") ~ "Missense",
    Variant_Classification == "Splice_Site" ~ "Splice_Site",
    TRUE ~ "Other"
  )) %>%
  select(-Variant_Classification) %>%
  distinct()

mutation_matrix <- mutation_data %>%
  pivot_wider(
    names_from = Hugo_Symbol, 
    values_from = Mutation_Type, 
    values_fill = NA,
    values_fn = list(Mutation_Type = function(x) {
      if ("Truncating" %in% x) {
        return("Truncating")
      } else if ("Missense" %in% x) {
        return("Missense")
      } else {
        return(unique(x)[1])
      }
    })
  )


