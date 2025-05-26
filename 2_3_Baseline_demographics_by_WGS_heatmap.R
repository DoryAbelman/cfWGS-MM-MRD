# figure1_integrated_heatmaps.R
# ------------------------------------------------------------
# Generate Figure 1: Integrated alteration heatmaps
# – Overlaid bone-marrow (upper triangles) vs plasma cfDNA (lower triangles)
# – Samples split by clinical cohort (Newly Diagnosed vs Pre-treated)
# – Sorted by tumour fraction and paired status, annotated with:
#     • cfDNA tumour fraction (points)
#     • BM & blood mutation counts (continuous colour bar)
#     • FISH positivity (★)
# – Rows split into Mutations, CNAs, Translocations
# – Saves individual BM, blood and combined overlay heatmaps + underlying data
# Author: Dory Abelman · Last edit: 2025-05-26
# ------------------------------------------------------------


# Load required libraries
library(maftools)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(purrr)
library(stringr)
library(readr)
library(grid)


### Load data 
# 0. clinical metadata
metada_df_mutation_comparison <- read_csv(
  "combined_clinical_data_updated_Feb5_2025.csv"
) %>%
  # recreate the “clean” BAM key that all your joins rely on
  mutate(
    Tumor_Sample_Barcode = Bam %>%
      str_remove_all("_PG|_WG") %>%
      str_replace_all("\\.filter.*|\\.ded.*|\\.recalibrate.*", ""),
    Bam_clean_tmp = str_remove(Bam, "\\.bam$")
  )

# Load cohort assignments
cohort_df <- read.table("cohort_assignment_table.txt", sep = "\t", header = TRUE)

# Merge into metadata
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  left_join(cohort_df, by = "Patient")


# 1. tumour‐fraction from ichor
tumor_fraction <- read_tsv("Oct 2024 data/tumor_fraction_cfWGS.txt")

# 2. arm‐level CNA matrix (was saved as cna_data.rds in your export_dir)
export_dir <- "Jan2025_exported_data"  
cna_data <- readRDS(file.path(export_dir, "cna_data.rds"))
# your code expects myeloma_CNA_matrix_with_HRD
myeloma_CNA_matrix_with_HRD <- cna_data  

# 3. cfDNA translocation tab (the 4‐call binary table you exported)
translocation_data <- readRDS(
  file.path(export_dir, "translocation_data.rds")
)

# 4. MAF objects for BM and blood
#    these point at the .maf you dumped out in your mutation script
maf_object_bm    <- read.maf(maf = "combined_maf_temp_bm_Jan2025.maf")
maf_object_blood <- read.maf(maf = "combined_maf_temp_blood_Jan2025.maf")


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


### This script first does seperate heatmaps organized in the same order and then a combined heatmap
## The combined heatmap has updated sample groupings and is the figure used in the manuscript

###########################################
### 1. Bone Marrow (BM) Heatmap Section ###
###########################################

# 1a. Process mutation data for BM
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

# 1b. Prepare CNA and translocation data for BM
cna_data <- myeloma_CNA_matrix_with_HRD %>%
  select(Sample, del1p, amp1q, del13q, del17p, hyperdiploid) %>%
  mutate_at(vars(-Sample), as.character)

translocation_data <- translocation_matrix %>%
  select(Bam_File, IGH_MAF = `IGH-MAF`, IGH_CCND1 = `IGH-CCND1`, IGH_MYC = `IGH-MYC`, IGH_FGFR3 = `IGH-FGFR3`) %>%
  mutate(Sample = gsub(".bam$", "", Bam_File))

CNA_translocation <- full_join(cna_data, translocation_data)

# 1c. Merge with metadata and tumor fraction info (metadata must contain a "Patient" column)
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Bam_clean_tmp = gsub(".bam$", "", Bam))

CNA_translocation <- left_join(CNA_translocation, 
                               metada_df_mutation_comparison, 
                               by = c("Sample" = "Bam_clean_tmp"))

tumor_fraction_max <- tumor_fraction %>%
  group_by(Bam) %>%
  summarise(Tumor_Fraction = max(Tumor_fraction, na.rm = TRUE))

CNA_translocation <- left_join(CNA_translocation, 
                               tumor_fraction_max, 
                               by = c("Sample" = "Bam"))

# 1d. Merge mutation with CNA/translocation data for BM
combined_data_heatmap_BM <- mutation_matrix %>%
  full_join(CNA_translocation, by = "Tumor_Sample_Barcode")

# 1e. Filter for Baseline/Progression and BM cells; reformat timepoint labels
combined_data_heatmap_BM <- combined_data_heatmap_BM %>%
  mutate(
    timepoint_info = ifelse(timepoint_info == "Relapse", "Progression", timepoint_info),
    timepoint_info = ifelse(timepoint_info == "Diagnosis", "Baseline", timepoint_info)
  ) %>%
  filter(timepoint_info %in% c("Baseline", "Progression")) %>%
  filter(Sample_type %in% c("BM_cells"))

# 1f. Replace NAs: for mutation columns use "No Mutation" and for CNA/translocation use 0
existing_cols <- intersect(myeloma_genes, colnames(combined_data_heatmap_BM))
combined_data_heatmap_BM <- combined_data_heatmap_BM %>%
  mutate(across(all_of(existing_cols), ~ ifelse(is.na(.), "No Mutation", .)))

cna_trans_cols <- c("del1p", "amp1q", "del13q", "del17p", "hyperdiploid", 
                    "Bam_File", "IGH_MAF", "IGH_CCND1", "IGH_MYC", "IGH_FGFR3")
combined_data_heatmap_BM <- combined_data_heatmap_BM %>%
  mutate(across(all_of(cna_trans_cols), ~ ifelse(is.na(.), 0, .)))

# 1g. Convert CNA and translocation indicators to factors
cna_cols <- c("del1p", "amp1q", "del13q", "del17p")
translocation_cols <- c("IGH_MAF", "IGH_CCND1", "IGH_MYC", "IGH_FGFR3")
combined_data_heatmap_BM <- combined_data_heatmap_BM %>%
  mutate_at(vars(one_of(cna_cols)), ~ ifelse(. == "1", "Yes", "No")) %>%
  mutate_at(vars(one_of(translocation_cols)), ~ ifelse(. == 1, "Yes", "No")) %>% 
  mutate(hyperdiploid = ifelse(hyperdiploid == 1, "Yes", "No")) ## used to be TRUE

# 1h. Remove specific low-confidence duplicate call and ensure uniqueness
combined_data_heatmap_BM <- combined_data_heatmap_BM %>%
  filter(!(Tumor_Sample_Barcode == "TFRIM4_0187_Bm_P_VA-13-01-O-DNA" & IGH_CCND1 == "Yes")) %>%
  unique()

# 1i. ***Order BM samples by Tumor_Fraction (high to low)***
# Order BM samples by descending Tumor_Fraction and then create the composite key
combined_data_heatmap_BM <- combined_data_heatmap_BM %>%
  arrange(desc(Tumor_Fraction)) %>% 
  mutate(Patient_Timepoint = paste0(Patient, "_", timepoint_info))

## Change duplicate SPORE entries 
combined_data_heatmap_BM$Patient_Timepoint[combined_data_heatmap_BM$Tumor_Sample_Barcode == "SPORE_0009_Bm_T_BM066745"] <- "SPORE_0009_Baseline2"
combined_data_heatmap_BM$Patient_Timepoint[combined_data_heatmap_BM$Tumor_Sample_Barcode == "SPORE_0012_Bm_T_BM069319"] <- "SPORE_0012_Progression2"


### Adjust based on blood (run later)
# Identify which BM samples have a matching blood sample
blood_pt <- combined_data_heatmap_blood$Patient_Timepoint

BM_matched <- combined_data_heatmap_BM %>% 
  filter(Patient_Timepoint %in% blood_pt)

BM_unmatched <- combined_data_heatmap_BM %>% 
  filter(!(Patient_Timepoint %in% blood_pt))

# Combine them: matched first (left) then unmatched (right)
combined_data_heatmap_BM <- bind_rows(BM_matched, BM_unmatched)

# Save the composite order
bm_order_composite <- combined_data_heatmap_BM$Patient_Timepoint

# 1j. Prepare the BM heatmap matrix
# Keep the Patient column for ordering; then drop it when forming the matrix
rownames(combined_data_heatmap_BM) <- combined_data_heatmap_BM$Tumor_Sample_Barcode

temp1 <- combined_data_heatmap_BM$Patient_Timepoint

# 2. Remove the unneeded columns but retain the row names
temp2 <- combined_data_heatmap_BM %>%
  select(-c(Bam, Date_of_sample_collection, Sample_type, Timepoint, Study, 
            Sample_ID, Tumor_Sample_Barcode, Sample, timepoint_info, Tumor_Fraction, Patient, Patient_Timepoint))
temp2 <- temp2 %>%
  select(-c(Num_days_to_closest_relapse_absolute, Num_days_to_closest_relapse_non_absolute, Num_days_to_closest_relapse, Relapsed, `...27`, Bam_File))

rownames(temp2) <- temp1

# 3. Transpose the data to make patients as columns, keeping the row names
heatmap_matrix_BM <- as.matrix(temp2)
heatmap_matrix_BM <- t(heatmap_matrix_BM)

# 1k. Create top annotation for BM
top_annotation_BM <- HeatmapAnnotation(
  Timepoint = combined_data_heatmap_BM$timepoint_info,
  Study = combined_data_heatmap_BM$Study,
  col = list(
    Timepoint = c("Baseline" = "#377eb8", "Progression" = "#e41a1c"),
    Study = c("M4" = "#984ea3", "MyC" = "#ff7f00", "SPORE" = "#a65628")
  ),
  Tumor_Fraction = anno_points(
    combined_data_heatmap_BM$Tumor_Fraction,
    gp = gpar(col = "darkgrey"),
    size = unit(2, "mm"),
    axis = TRUE,
    pch = 16,
    ylim = c(0, 1),
    add_lines = TRUE,
    baseline = 0.3,
    baseline_gp = gpar(col = "black", lty = 2)
  ),
  show_annotation_name = TRUE
)

combined_colors <- c(
  "Truncating" = "#d73027",
  "Missense" = "#fc8d59",
  "Splice_Site" = "#fee090",
  "Other" = "#ffffbf",
  "No Mutation" = "#f0f0f0",
  "Yes" = "#1a9641",
  "No" = "#f0f0f0"
)

heatmap_BM <- Heatmap(
  heatmap_matrix_BM,
  name = "Alteration",
  col = combined_colors,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation = top_annotation_BM,
  heatmap_legend_param = list(
    title = "Alteration",
    at = c("Truncating", "Missense", "Splice_Site", "Other", "No Mutation", "Yes", "No"),
    labels = c("Truncating", "Missense", "Splice Site", "Other", "No Mutation", "Yes", "No")
  )  
)

draw(heatmap_BM)
png("heatmap_output_BM_baseline_updated_11.png", width = 11.12, height = 8, units = "in", res = 500)
draw(heatmap_BM)
dev.off()

###########################################
### 2. Blood cfDNA Heatmap Section      ###
###########################################

# 2a. Process mutation data for blood
maf_subset_blood <- subsetMaf(maf = maf_object_blood, genes = myeloma_genes, includeSyn = FALSE)
mutation_data_blood <- maf_subset_blood@data %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
  mutate(Mutation_Type = case_when(
    Variant_Classification %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins") ~ "Truncating",
    Variant_Classification %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins") ~ "Missense",
    Variant_Classification == "Splice_Site" ~ "Splice_Site",
    TRUE ~ "Other"
  )) %>%
  select(-Variant_Classification) %>%
  distinct()

mutation_matrix_blood <- mutation_data_blood %>%
  pivot_wider(
    names_from = Hugo_Symbol, 
    values_from = Mutation_Type, 
    values_fill = NA
  )
# 2b. ***Subset to include only BM mutation genes for alignment***
mutation_matrix_blood <- mutation_matrix_blood %>%
  select(Tumor_Sample_Barcode, all_of(intersect(myeloma_genes, colnames(.))))

# 2c. Merge with CNA/translocation data (same as BM)
combined_data_heatmap_blood <- mutation_matrix_blood %>%
  full_join(CNA_translocation, by = "Tumor_Sample_Barcode")

combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  mutate(
    timepoint_info = ifelse(timepoint_info == "Relapse", "Progression", timepoint_info),
    timepoint_info = ifelse(timepoint_info == "Diagnosis", "Baseline", timepoint_info)
  ) %>%
  filter(timepoint_info %in% c("Baseline", "Progression")) %>%
  filter(Sample_type %in% c("Blood_plasma_cfDNA"))

# 2d. Replace NAs for mutation columns and for CNA/translocation columns
existing_cols_blood <- intersect(myeloma_genes, colnames(combined_data_heatmap_blood))
combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  mutate(across(all_of(existing_cols_blood), ~ ifelse(is.na(.), "No Mutation", .)))

combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  mutate(across(all_of(cna_trans_cols), ~ ifelse(is.na(.), 0, .)))

combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  select(-Bam_File, -Relapsed, -Num_days_to_closest_relapse_absolute, 
         -Num_days_to_closest_relapse_non_absolute, -Num_days_to_closest_relapse, -`...1`)

combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  mutate_at(vars(one_of(cna_cols)), ~ ifelse(. == "1", "Yes", "No")) %>%
  mutate_at(vars(one_of(translocation_cols)), ~ ifelse(. == 1, "Yes", "No")) %>% 
  mutate(hyperdiploid = ifelse(hyperdiploid == 1, "Yes", "No"))

# 2e. Remove low-quality IGH-MAF call from blood
combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  filter(!(Tumor_Sample_Barcode == "TFRIM4_0181_Cf_P_VA-07-05-P-DNA" & IGH_MAF == "Yes"))

# 2f. ***Reorder blood samples using Patient info***
# (Assuming the "Patient" column is available from the join)
# First, extract the blood Patient IDs
# Ensure Patient is a character and create a composite key for blood samples
combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  mutate(Patient = as.character(Patient),
         Patient_Timepoint = paste0(Patient, "_", timepoint_info))

# Split blood samples into those that are common with BM (based on composite key) and extras
common_blood <- combined_data_heatmap_blood %>% 
  filter(Patient_Timepoint %in% bm_order_composite) %>%
  arrange(factor(Patient_Timepoint, levels = unique(bm_order_composite)))

extra_blood <- combined_data_heatmap_blood %>% 
  filter(!(Patient_Timepoint %in% bm_order_composite))

# Define final ordering for blood samples
final_blood_order <- c(common_blood$Tumor_Sample_Barcode, extra_blood$Tumor_Sample_Barcode)

# Reorder blood samples based on the composite key order
combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = final_blood_order)) %>%
  arrange(Tumor_Sample_Barcode)

# 2g. Prepare the blood heatmap matrix (retain Patient for ordering then drop)

## Remove duplicate progression 
combined_data_heatmap_blood <- combined_data_heatmap_blood %>%
  mutate(Patient_Timepoint = ifelse(
    Sample_ID == "SPORE_0006_T4_Blood_plasma_cfDNA" & Patient_Timepoint == "SPORE_0006_Progression",
    "SPORE_0006_Progression2",
    Patient_Timepoint
  ))

rownames(combined_data_heatmap_blood) <- combined_data_heatmap_blood$Patient_Timepoint

temp1 <- combined_data_heatmap_blood$Patient_Timepoint

# 2. Remove the unneeded columns but retain the row names
temp2 <- combined_data_heatmap_blood %>%
  select(-c(Bam, Date_of_sample_collection, Sample_type, Timepoint, Study, 
            Sample_ID, Tumor_Sample_Barcode, Sample, timepoint_info, Tumor_Fraction, Patient, Patient_Timepoint))

rownames(temp2) <- temp1

# 3. Transpose the data to make patients as columns, keeping the row names
heatmap_matrix_blood <- as.matrix(temp2)

#rownames(heatmap_matrix) <- rownames(combined_data_heatmap)  # Keep row names intact
heatmap_matrix_blood <- t(heatmap_matrix_blood)  # Transpose the matrix

### Align the matrix to the BM
#    Initialize everything as "No Mutation"
# Capture the gene (row) names from your BM heatmap matrix
bm_genes <- rownames(heatmap_matrix_BM)

blood_aligned <- matrix(
  "No Mutation",
  nrow = length(bm_genes),
  ncol = ncol(heatmap_matrix_blood),
  dimnames = list(bm_genes, colnames(heatmap_matrix_blood))
)

# 3) For the genes that do exist in the Blood matrix, fill them in
common_genes <- intersect(bm_genes, rownames(heatmap_matrix_blood))
blood_aligned[common_genes, ] <- heatmap_matrix_blood[common_genes, ]

heatmap_matrix_blood <- blood_aligned

# 2h. Create top annotation for blood
top_annotation_blood <- HeatmapAnnotation(
  Timepoint = combined_data_heatmap_blood$timepoint_info,
  Study = combined_data_heatmap_blood$Study,
  col = list(
    Timepoint = c("Baseline" = "#377eb8", "Progression" = "#e41a1c"),
    Study = c("M4" = "#984ea3", "MyC" = "#ff7f00", "SPORE" = "#a65628")
  ),
  Tumor_Fraction = anno_points(
    combined_data_heatmap_blood$Tumor_Fraction,
    gp = gpar(col = "darkgrey"),
    size = unit(2, "mm"),
    axis = TRUE,
    pch = 16,
    ylim = c(0, 1),
    add_lines = TRUE,
    baseline = 0.3,
    baseline_gp = gpar(col = "black", lty = 2)
  ),
  show_annotation_name = TRUE
)

heatmap_blood <- Heatmap(
  heatmap_matrix_blood,
  name = "Alteration",
  col = combined_colors,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation = top_annotation_blood,
  heatmap_legend_param = list(
    title = "Alteration",
    at = c("Truncating", "Missense", "Splice_Site", "Other", "No Mutation", "Yes", "No"),
    labels = c("Truncating", "Missense", "Splice Site", "Other", "No Mutation", "Yes", "No")
  )
)

draw(heatmap_blood)
png("heatmap_output_Blood_baseline_updated_11.png", width = 13.5, height = 8, units = "in", res = 500)
draw(heatmap_blood)
dev.off()


### See difference 
# Rows in blood but not in BM
setdiff(rownames(heatmap_matrix_blood), rownames(heatmap_matrix_BM))

# Rows in BM but not in blood
setdiff(rownames(heatmap_matrix_BM), rownames(heatmap_matrix_blood))


## Export 
# save the filtered matrices for use elsewhere
saveRDS(heatmap_matrix_BM,    file = "heatmap_matrix_BM_March2025.rds")
saveRDS(heatmap_matrix_blood, file = "heatmap_matrix_blood_March2025.rds")

# Save bone marrow combined data
saveRDS(combined_data_heatmap_BM, file = "combined_data_heatmap_BM_March2025.rds")
# Save cfDNA (blood) combined data
saveRDS(combined_data_heatmap_blood, file = "combined_data_heatmap_blood_March2025.rds")


################################################################################
### 3.  OVERLAY HEATMAP (BM ⬆︎ / cfDNA ⬇︎)  ####################################
################################################################################

## ---------- 3a.  make sure both matrices have identical dims -----------------
# rows = features (genes + CNA + translocations)
# cols = Patient_Timepoint
# 0) load your cohort assignment (from the previous export)
cohort_df     <- readRDS("cohort_assignment_table.rds")
keep_patients <- cohort_df$Patient

# helper to pull the patient ID off a column like "SPORE_0009_Baseline"
extract_pid <- function(cols) {
  sapply(strsplit(cols, "_"), function(x) paste(head(x, -1), collapse = "_"))
}

# 1) figure out which columns to keep
bm_cols    <- colnames(heatmap_matrix_BM)
bm_pids    <- extract_pid(bm_cols)
keep_bm    <- bm_cols[
  bm_pids %in% keep_patients &
    grepl("Baseline|Diagnosis", bm_cols) &
    !grepl("_Baseline2$", bm_cols)   # drop any Baseline2
]

blood_cols <- colnames(heatmap_matrix_blood)
blood_pids <- extract_pid(blood_cols)
keep_blood <- blood_cols[blood_pids %in% keep_patients & grepl("Baseline|Diagnosis", blood_cols)]

# 2) unify the final column set
all_cols <- base::union(keep_bm, keep_blood)
pre_order_cols <- all_cols  

# 3) subset your matrices
heatmap_matrix_BM_subset <- heatmap_matrix_BM[, keep_bm, drop = FALSE]
heatmap_matrix_blood_subset <- heatmap_matrix_blood[, keep_blood, drop = FALSE]

## 3a. scaffold with NA
bm_mat    <- matrix(NA_character_, nrow = length(all_rows), ncol = length(all_cols),
                    dimnames = list(all_rows, all_cols))
cfDNA_mat <- bm_mat

bm_mat[rownames(heatmap_matrix_BM_subset),   colnames(heatmap_matrix_BM_subset)]       <-
  heatmap_matrix_BM_subset
cfDNA_mat[rownames(heatmap_matrix_blood_subset), colnames(heatmap_matrix_blood_subset)] <-
  heatmap_matrix_blood_subset

## ---------- 3b.  colour keys --------------------------------------------------
pal_mut <- c(
  "Truncating"  = "#d73027",
  "Missense"    = "#fc8d59",
  "Splice_Site" = "#fee090",
  "Other"       = "#ffffbf",
  "No Mutation" = "#f0f0f0"
)
pal_bin <- c("Yes" = "#1a9641", "No" = "#f0f0f0")
pal_na  <- "#d9d9d9"

value_to_col <- function(v) {
  if (is.na(v)) {
    pal_na
  } else if (v %in% names(pal_mut)) {
    pal_mut[v]
  } else if (v %in% names(pal_bin)) {
    pal_bin[v]
  } else {
    pal_na
  }
}


## ---------- 3c.  column annotation: cohort & tumour fraction -----------------
## Ensure a match 
unmatched <- all_cols[ is.na( match(all_cols,
                                    combined_data_heatmap_blood$Patient_Timepoint) ) ]
print(unmatched)

# then rebuild col_tf:
col_tf <- combined_data_heatmap_BM$Tumor_Fraction[
  match(all_cols, combined_data_heatmap_BM$Patient_Timepoint)
]

col_tf      <- combined_data_heatmap_blood$Tumor_Fraction[match(all_cols,
                                                                combined_data_heatmap_blood$Patient_Timepoint)]  # Blood TF (use NA if none)

# Extract the patient ID from each Patient_Timepoint by dropping only the last “_…” part:
patients <- sapply(
  strsplit(all_cols, "_"), 
  function(x) paste(head(x, -1), collapse = "_")
)

# 2. Look up the NewDx flag in your cohort_df
flags <- cohort_df$NewDx[match(patients, cohort_df$Patient)]

# 3. Build the cohort vector: 1 → "NewDx", 0 → "Pre-treated", NA → NA
col_cohort <- ifelse(
  is.na(flags), 
  NA_character_, 
  ifelse(flags == 1, "NewDx", "Pre-treated")
)

# A. Force the factor to have exactly these two levels (and nothing else):
col_cohort <- factor(col_cohort, levels = c("NewDx", "Pre-treated"))

# B. Double-check what you have:
table(col_cohort, useNA = "ifany")


## ---------- 3d.  custom cell_fun draws the split square ----------------------
# 1) compute which patients really *are* paired
extract_pid <- function(x) sub("_(?!.*_).*", "", x, perl = TRUE)
bm_pids    <- extract_pid(keep_bm)
blood_pids <- extract_pid(keep_blood)
paired_pids <- intersect(bm_pids, blood_pids)

# 2) build a little ordering data.frame
ord_df <- data.frame(
  Sample        = all_cols,
  Cohort        = col_cohort,
  TumourFraction= col_tf,
  Paired        = extract_pid(all_cols) %in% paired_pids,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Cohort)) %>%                            # drop anything outside your two cohorts
  arrange(Cohort, desc(Paired), desc(TumourFraction))  # within each cohort: paired first, by TF

col_order      <- ord_df$Sample
col_cohort_ord <- ord_df$Cohort
col_tf_ord     <- ord_df$TumourFraction

# rename  cohort factor levels
col_cohort_ord <- factor(
  col_cohort_ord,
  levels = c("NewDx", "Pre-treated"),
  labels = c("Newly Diagnosed", "Pre-treated")
)

# 3) re‐order your matrices
bm_mat    <- bm_mat   [, col_order, drop=FALSE]
cfDNA_mat <- cfDNA_mat[, col_order, drop=FALSE]
all_cols  <- col_order

# 4) make a short label vector for the bottom
patient_labels <- extract_pid(all_cols)

# 5) re‐define your colour lookup
pal_mut <- c(
  Truncating  = "#d73027", Missense    = "#fc8d59",
  Splice_Site = "#fee090", Other       = "#ffffbf",
  `No Mutation` = "#f0f0f0"
)
pal_bin <- c(Yes = "#1a9641", No = "#f0f0f0")
pal_na  <- "#d9d9d9"
value_to_col <- function(v) {
  if      (is.na(v))     pal_na
  else if (v %in% names(pal_mut)) pal_mut[v]
  else if (v %in% names(pal_bin)) pal_bin[v]
  else                    pal_na
}

# 6) Set row splits 
# assuming you still have these defined
mut_genes   <- myeloma_genes
cna_features  <- c("del1p","amp1q","del13q","del17p","hyperdiploid")
trans_features <- c("IGH_MAF","IGH_CCND1","IGH_MYC","IGH_FGFR3")

row_group <- factor(
  ifelse(all_rows %in% mut_genes,      "Mutations",
         ifelse(all_rows %in% cna_features,   "CNAs",
                ifelse(all_rows %in% trans_features, "Translocations",
                       NA))),
  levels = c("Mutations","CNAs","Translocations")
)


# 7) build the top‐row annotation
# 7a) Bring in mutation count dataframe
#   (assumes you have a data.frame `mut_counts_df` with columns Patient, MutCount)
# (a) Load your baseline clinical table (if not already loaded)
#     and keep only one Diagnosis/Baseline row per patient.
dat_base <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.rds") %>%
  filter(timepoint_info %in% c("Diagnosis","Baseline")) %>%
  arrange(Patient, timepoint_info) %>%
  group_by(Patient) %>%
  slice(1) %>%    # first row
  ungroup()

# (b) Extract just the patient IDs from your heatmap columns:
patient_ids <- extract_pid(all_cols)

# (c) Pull out the two mutation‐count columns from dat_base
#     (these names must match exactly what's in your data frame)
bm_counts    <- dat_base$BM_Mutation_Count[    match(patient_ids, dat_base$Patient) ]
blood_counts <- dat_base$Blood_Mutation_Count[ match(patient_ids, dat_base$Patient) ]

# 2) set up continuous palette for all muts (you can pick your own colours):
library(circlize)
all_counts <- c(bm_counts, blood_counts)
qnts       <- quantile(all_counts, probs = c(0, 0.5, 1), na.rm = TRUE)

count_pal <- colorRamp2(
  qnts,
  c("#f7fbff",   # light 
    "#6baed6",   # mid-tone blue 
    "#08306b")   # deep blue
)

# 7b) Re-build the top annotation to include a little barplot of mutation counts
top_ha <- HeatmapAnnotation(
  Cohort         = col_cohort_ord,
  `cfDNA Tumour\nFraction` = anno_points(
    col_tf_ord,
    axis = TRUE,            # this one is valid for anno_points
    ylim = c(0, 0.4),
    gp = gpar(col = "darkgrey"),
    pch = 16, size = unit(2, "mm")
  ),
  BM_Muts        = anno_simple(
    bm_counts,
    col    = count_pal,    # continuous colour mapping
    na_col = "grey90"
  ),
  Blood_Muts     = anno_simple(
    blood_counts,
    col    = count_pal,
    na_col = "grey90"
  ),
  col = list(
    Cohort = c(`Newly Diagnosed` = "#3182bd", `Pre-treated` = "#e6550d")
  ),
  annotation_name_side = "left",
  annotation_name_rot  = 0,
  show_annotation_name = TRUE
)



# 8) draw the overlay heatmap

## set up legend
# 1) Build one Legend for the mutation types
lgd_mut <- Legend(
  title     = "Mutation Type",
  at        = c("Truncating", "Missense", "Splice_Site", "Other", "No Mutation"),
  legend_gp = gpar(fill = pal_mut[c("Truncating","Missense","Splice_Site","Other","No Mutation")]),
  direction = "vertical"
)

# 2) Build a second Legend for the binary CNA / translocation calls
lgd_cna <- Legend(
  title     = "CNA / Translocation",
  at        = c("Yes","No","Missing"),
  legend_gp = gpar(fill = c(pal_bin["Yes"], pal_bin["No"], pal_na)),
  direction = "vertical"
)

# continuous‐counts legend:
lgd_count <- Legend(
  title   = "Mutation count",
  col_fun = count_pal,
  at      = qnts,
  labels  = round(qnts,2),
  direction = "vertical"
)


overlay_ht <- Heatmap(
  matrix("", nrow=length(all_rows), ncol=length(all_cols),
         dimnames=list(all_rows, all_cols)),
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  show_row_names   = TRUE,
  show_column_names= TRUE,
  column_labels    = patient_labels,       # <-- show only the patient
  column_split     = col_cohort_ord,       # <-- vertical split by cohort,
  show_heatmap_legend = FALSE, # Since building own
  
  ## split rows 
  row_split        = row_group,
  row_title        = c("Mutations","CNAs","Translocations"),  # optional titles
  row_gap          = unit(c(2, 2), "mm"),                    # space between the splits
  
  row_names_gp     = gpar(fontsize=8),
  column_names_gp  = gpar(fontsize=8, fontface="bold"),
  top_annotation   = top_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # draw BM / cfDNA triangles
    col_bm <- value_to_col(bm_mat[i,j])
    col_cf <- value_to_col(cfDNA_mat[i,j])
    grid.polygon(
      x = unit.c(x - width*0.5, x + width*0.5, x - width*0.5),
      y = unit.c(y + height*0.5, y + height*0.5, y - height*0.5),
      gp = gpar(fill=col_bm, col=NA)
    )
    grid.polygon(
      x = unit.c(x + width*0.5, x + width*0.5, x - width*0.5),
      y = unit.c(y - height*0.5, y + height*0.5, y - height*0.5),
      gp = gpar(fill=col_cf, col=NA)
    )
  },
  heatmap_legend_param = list()
  )


draw(overlay_ht)

draw(
  overlay_ht,
  annotation_legend_side = "right",  # where cohort legend goes
  heatmap_legend_list    = list(lgd_mut, lgd_cna, lgd_count),
)

# save
png("overlay_heatmap_BM_vs_cfDNA_updated3.png", width = 14, height = 8, units = "in", res = 450)
draw(
  overlay_ht,
  annotation_legend_side = "right",  # where cohort legend goes
  heatmap_legend_list    = list(lgd_mut, lgd_cna, lgd_count),
)
dev.off()




#### Now add a star if also found by FISH 
#### Since so many FISH unknown, not sure if good idea or concentration of effort for now

# (you already have this loaded as `dat`; if not, the CSV path:)
dat_baseline <- dat_base %>%
  # only keep diagnosis / baseline rows
  filter(timepoint_info %in% c("Diagnosis","Baseline")) %>%
  
  # pick exactly the FISH columns that map to your heatmap rows:
  #   • T_4_14   → IGH_FGFR3
  #   • T_11_14  → IGH_CCND1
  #   • T_14_16  → IGH_MAF
  #   • DEL_17P  → del17p
  #   • DEL_1P   → del1p
  #   • AMP_1Q   → amp1q
  #   • DEL_13   → del13q
  #   • hyperdiploid → hyperdiploid
  select(
    Patient,
    T_4_14,
    T_11_14,
    T_14_16,
    DEL_17P,
    DEL_1P,
    AMP_1Q,
    DEL_13,
    hyperdiploid
  ) %>%
  
  mutate(across(-Patient, ~ case_when(
    .x == "Positive" ~ TRUE,
    .x == "Negative" ~ FALSE,
    TRUE             ~ NA   # Unknown or anything else becomes NA
  ))) %>%

  # (optional) rename to match your heatmap row‐names exactly:
  rename(
    IGH_FGFR3 = T_4_14,
    IGH_CCND1 = T_11_14,
    IGH_MAF   = T_14_16,
    del17p     = DEL_17P,
    del1p      = DEL_1P,
    amp1q      = AMP_1Q,
    del13q     = DEL_13,
    hyperdiploid = hyperdiploid
  )

# make a matrix of flags, rows = patient, cols = fish assay
fish_flags <- dat_baseline %>%
  column_to_rownames("Patient") %>%
  as.matrix()   

## Now keep only relevant patients 
fish_flags <- fish_flags[rownames(fish_flags) %in% patient_ids, , drop = FALSE]

fish_flags <- fish_flags %>% as.matrix()

### Redraw heatmap
overlay_ht_2 <- Heatmap(
  matrix("", nrow = length(all_rows), ncol = length(all_cols),
         dimnames = list(all_rows, all_cols)),
  cluster_rows        = FALSE,
  cluster_columns     = FALSE,
  show_row_names      = TRUE,
  show_column_names   = TRUE,
  column_labels       = patient_labels,
  column_split        = col_cohort_ord,
  show_heatmap_legend = FALSE,
  
  ## row split, fonts, annotation … (unchanged) ----------------
  row_split           = row_group,
  row_title           = c("Mutations","CNAs","Translocations"),
  row_gap             = unit(c(2, 2), "mm"),
  row_names_gp        = gpar(fontsize = 8),
  column_names_gp     = gpar(fontsize = 8, fontface = "bold"),
  top_annotation      = top_ha,
  
  ## ------------------  cell_fun  -----------------------------
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    ## draw BM / cfDNA halves ---------------------------------
    grid.polygon(
      x = unit.c(x - width*0.5, x + width*0.5, x - width*0.5),
      y = unit.c(y + height*0.5, y + height*0.5, y - height*0.5),
      gp = gpar(fill = value_to_col(bm_mat[i, j]), col = NA)
    )
    grid.polygon(
      x = unit.c(x + width*0.5, x + width*0.5, x - width*0.5),
      y = unit.c(y - height*0.5, y + height*0.5, y - height*0.5),
      gp = gpar(fill = value_to_col(cfDNA_mat[i, j]), col = NA)
    )
    
    ## overlay a star if that FISH call is positive -----------
    pid     <- extract_pid(all_cols[j])
    feature <- all_rows[i]
    
    if (pid %in% rownames(fish_flags) &&
        feature %in% colnames(fish_flags) &&
        isTRUE(fish_flags[pid, feature])) {
      
      grid.points(
        x, y,
        pch  = 8,                  # 8 is the star-glyph
        size = unit(2, "mm"),      
        gp   = gpar(col = "black")
      )
    }
  },
  
  heatmap_legend_param = list()
)

## Add FISH legend 
lgd_fish <- Legend(
  title     = "FISH call",
  labels    = "Positive by FISH",
  type      = "points",
  pch       = 8,                         # the same shape you used in cell_fun
  legend_gp = gpar(col = "black"),   # same colour
  direction = "vertical"
)


# save
png("overlay_heatmap_BM_vs_cfDNA_updated_with_FISH.png", width = 14, height = 8, units = "in", res = 450)
draw(
  overlay_ht_2,
  annotation_legend_side = "right",  # where cohort legend goes
  heatmap_legend_list    = list(lgd_mut, lgd_cna, lgd_count, lgd_fish),
)
dev.off()









########  Now summarise overlap and correlations 
#####
# 1. Build Long-form Mutation Data
#####

# Define gene list for mutations (make sure these match your data)
myeloma_genes <- c(
  "TP53", "KRAS", "NRAS", "BRAF", "FAM46C", "DIS3", "CYLD", "ATM",
  "CCND1", "MYC", "RB1", "TRAF3", "IRF4", "FGFR3", "MMSET", "BCL2",
  "IKZF1", "IKZF3", "CDKN2C", "KDM6A", "SETD2", "PTEN", "XBP1", "MAX",
  "SP140", "NFKBIA", "NFKB2", "PRDM1", "EGR1", "LTB"
)

# Define CNA and translocation column names
cna_cols <- c("del1p", "amp1q", "del13q", "del17p", "hyperdiploid")
transloc_cols <- c("IGH_MAF", "IGH_CCND1", "IGH_MYC", "IGH_FGFR3")

# --- Step 1. Build Full Mutation Lists for BM and cfDNA ----------------------
# For each sample, create a list of mutated genes; if no gene is mutated, an empty vector is returned.

BM_long_all <- combined_data_heatmap_BM %>%
  select(Patient_Timepoint, Tumor_Fraction, timepoint_info, one_of(myeloma_genes)) %>%
  rename(Tumor_Fraction_BM = Tumor_Fraction,
         timepoint_info_BM = timepoint_info) %>%
  pivot_longer(
    cols = one_of(myeloma_genes),
    names_to = "Gene",
    values_to = "Mutation_Type_BM"
  ) %>%
  group_by(Patient_Timepoint, Tumor_Fraction_BM, timepoint_info_BM) %>%
  summarize(
    MutatedGenes_BM = list(unique(Gene[Mutation_Type_BM != "No Mutation"])),
    .groups = "drop"
  )

blood_long_all <- combined_data_heatmap_blood %>%
  select(Patient_Timepoint, Tumor_Fraction, timepoint_info, one_of(myeloma_genes)) %>%
  rename(Tumor_Fraction_blood = Tumor_Fraction,
         timepoint_info_blood = timepoint_info) %>%
  pivot_longer(
    cols = one_of(myeloma_genes),
    names_to = "Gene",
    values_to = "Mutation_Type_blood"
  ) %>%
  group_by(Patient_Timepoint, Tumor_Fraction_blood, timepoint_info_blood) %>%
  summarize(
    MutatedGenes_blood = list(unique(Gene[Mutation_Type_blood != "No Mutation"])),
    .groups = "drop"
  )

# --- Step 2. Merge Mutation Data (Include Negatives) -------------------------
# A full join ensures that even samples with no mutations in one compartment are included.
merged_mut <- full_join(BM_long_all, blood_long_all, by = "Patient_Timepoint", 
                        suffix = c("_BM", "_blood")) %>%
  mutate(
    # Use coalesce so that if one compartment is missing, we still have the timepoint and TF info.
    Tumor_Fraction_BM = coalesce(Tumor_Fraction_BM, NA_real_),
    Tumor_Fraction_blood = coalesce(Tumor_Fraction_blood, NA_real_),
    timepoint = coalesce(timepoint_info_BM, timepoint_info_blood)
  ) %>%
  # Replace NULL lists with empty vectors.
  mutate(
    MutatedGenes_BM = map(MutatedGenes_BM, ~ if (is.null(.x)) character(0) else .x),
    MutatedGenes_blood = map(MutatedGenes_blood, ~ if (is.null(.x)) character(0) else .x)
  )

# --- Step 3. Compute Per-Sample Mutation Metrics ------------------------------
# For each sample, count the number of BM mutations, blood mutations, their intersection, union, and derive an overlap ratio.
# Define overlap = 1 if both compartments are negative.
merged_mut <- merged_mut %>%
  rowwise() %>%
  mutate(
    n_BM = length(MutatedGenes_BM),
    n_blood = length(MutatedGenes_blood),
    intersect_genes = list(base::intersect(MutatedGenes_BM, MutatedGenes_blood)),
    n_intersect = length(intersect_genes),
    union_genes = list(union(MutatedGenes_BM, MutatedGenes_blood)),
    n_union = length(union_genes),
    overlap_ratio = ifelse(n_union == 0, 1, n_intersect / n_union),
    BM_only = list(setdiff(MutatedGenes_BM, MutatedGenes_blood)),
    blood_only = list(setdiff(MutatedGenes_blood, MutatedGenes_BM)),
    n_BM_only = length(BM_only),
    n_blood_only = length(blood_only)
  ) %>%
  ungroup()

# --- Step 4. Define TF Category for Each Matched Sample ----------------------
# TF_category = "High TF" if both BM and blood tumor fractions >5% (i.e. >0.05), otherwise "Low TF".
merged_mut <- merged_mut %>%
  mutate(
    TF_category = case_when(
      is.na(Tumor_Fraction_BM) ~ "BM NA",                # Exclude BM missing
      is.na(Tumor_Fraction_blood) ~ "Blood NA",           # Exclude blood missing
      Tumor_Fraction_BM > 0.05 & Tumor_Fraction_blood > 0.05 ~ "High TF",
      !is.na(Tumor_Fraction_BM) & Tumor_Fraction_blood <= 0.05 ~ "Low TF",
      TRUE ~ "Discordant"
    )
  )


# --- Step 5. Summarize Mutation Concordance by Timepoint and TF Category --------
mutation_summary <- merged_mut %>%
  group_by(timepoint, TF_category) %>%
  summarize(
    n_samples = n(),
    total_BM_mut = sum(n_BM, na.rm = TRUE),
    total_blood_mut = sum(n_blood, na.rm = TRUE),
    total_intersect = sum(n_intersect, na.rm = TRUE),
    avg_overlap = mean(overlap_ratio, na.rm = TRUE),
    total_BM_only = sum(n_BM_only, na.rm = TRUE),
    total_blood_only = sum(n_blood_only, na.rm = TRUE),
    false_positive_rate = ifelse(total_blood_mut == 0, NA, 100 * total_blood_only / total_blood_mut),
    avg_TF_BM = mean(Tumor_Fraction_BM, na.rm = TRUE),
    avg_TF_blood = mean(Tumor_Fraction_blood, na.rm = TRUE)
  ) %>% 
  ungroup()

# Print mutation summary table
cat("### Mutation Summary by Timepoint and Tumor Fraction Category ###\n")
print(mutation_summary)

# --- Step 6. Translocation Concordance ----------------------------------------
# For translocations, we use the columns in each dataset. First, extract translocation calls.
BM_trans <- combined_data_heatmap_BM %>%
  select(Patient_Timepoint, Tumor_Fraction, timepoint_info, one_of(transloc_cols)) %>%
  rename(Tumor_Fraction_BM = Tumor_Fraction,
         timepoint_info_BM = timepoint_info) %>%
  mutate(across(one_of(transloc_cols), ~ ifelse(.=="Yes", "Yes", "No")))

blood_trans <- combined_data_heatmap_blood %>%
  select(Patient_Timepoint, Tumor_Fraction, timepoint_info, one_of(transloc_cols)) %>%
  rename(Tumor_Fraction_blood = Tumor_Fraction,
         timepoint_info_blood = timepoint_info) %>%
  mutate(across(one_of(transloc_cols), ~ ifelse(.=="Yes", "Yes", "No")))

# Merge translocation data
merged_trans <- full_join(BM_trans, blood_trans, by="Patient_Timepoint", 
                          suffix = c("_BM", "_blood")) %>%
  mutate(
    Tumor_Fraction_BM = coalesce(Tumor_Fraction_BM, NA_real_),
    Tumor_Fraction_blood = coalesce(Tumor_Fraction_blood, NA_real_),
    timepoint = coalesce(timepoint_info_BM, timepoint_info_blood),
    TF_category = case_when(
      is.na(Tumor_Fraction_BM) ~ "BM NA",
      is.na(Tumor_Fraction_blood) ~ "Blood NA",
      Tumor_Fraction_BM > 0.05 & Tumor_Fraction_blood > 0.05 ~ "High TF",
      !is.na(Tumor_Fraction_BM) & Tumor_Fraction_blood <= 0.05 ~ "Low TF",
      TRUE ~ "Discordant"
    )
  )
# For each sample, compute the number of positive translocation calls in BM and blood.
# After merging, the BM columns become "IGH_MAF_BM", etc., and blood columns "IGH_MAF_blood", etc.
BM_trans_cols <- paste0(transloc_cols, "_BM")
blood_trans_cols <- paste0(transloc_cols, "_blood")

merged_trans <- merged_trans %>%
  rowwise() %>%
  mutate(
    n_trans_BM = sum(c_across(all_of(BM_trans_cols)) == "Yes", na.rm = TRUE),
    n_trans_blood = sum(c_across(all_of(blood_trans_cols)) == "Yes", na.rm = TRUE),
    n_trans_intersect = sum((c_across(all_of(BM_trans_cols)) == "Yes") &
                              (c_across(all_of(blood_trans_cols)) == "Yes"), na.rm = TRUE),
    n_trans_union = sum((c_across(all_of(BM_trans_cols)) == "Yes") |
                          (c_across(all_of(blood_trans_cols)) == "Yes"), na.rm = TRUE),
    trans_overlap = ifelse(n_trans_union == 0, 1, n_trans_intersect / n_trans_union),
    trans_BM_only = n_trans_BM - n_trans_intersect,
    trans_blood_only = n_trans_blood - n_trans_intersect
  ) %>%
  ungroup()

trans_summary <- merged_trans %>%
  group_by(timepoint, TF_category) %>%
  summarize(
    n_samples = n(),
    total_trans_BM = sum(n_trans_BM, na.rm = TRUE),
    total_trans_blood = sum(n_trans_blood, na.rm = TRUE),
    total_trans_intersect = sum(n_trans_intersect, na.rm = TRUE),
    avg_trans_overlap = mean(trans_overlap, na.rm = TRUE),
    total_trans_BM_only = sum(trans_BM_only, na.rm = TRUE),
    total_trans_blood_only = sum(trans_blood_only, na.rm = TRUE),
    false_positive_rate = ifelse(total_trans_blood == 0, NA, 100 * total_trans_blood_only / total_trans_blood),
    avg_TF_BM = mean(Tumor_Fraction_BM, na.rm = TRUE),
    avg_TF_blood = mean(Tumor_Fraction_blood, na.rm = TRUE)
  ) %>%
  ungroup()

cat("\n### Translocation Summary by Timepoint and Tumor Fraction Category ###\n")
print(trans_summary)

# --- Step 7. CNA Concordance --------------------------------------------------
BM_CNA <- combined_data_heatmap_BM %>%
  select(Patient_Timepoint, Tumor_Fraction, timepoint_info, one_of(cna_cols)) %>%
  rename(Tumor_Fraction_BM = Tumor_Fraction,
         timepoint_info_BM = timepoint_info) %>%
  mutate(across(one_of(cna_cols), ~ ifelse(.=="Yes", "Yes", "No")))
blood_CNA <- combined_data_heatmap_blood %>%
  select(Patient_Timepoint, Tumor_Fraction, timepoint_info, one_of(cna_cols)) %>%
  rename(Tumor_Fraction_blood = Tumor_Fraction,
         timepoint_info_blood = timepoint_info) %>%
  mutate(across(one_of(cna_cols), ~ ifelse(.=="Yes", "Yes", "No")))

merged_CNA <- full_join(BM_CNA, blood_CNA, by = "Patient_Timepoint", 
                        suffix = c("_BM", "_blood")) %>%
  mutate(
    Tumor_Fraction_BM = coalesce(Tumor_Fraction_BM, NA_real_),
    Tumor_Fraction_blood = coalesce(Tumor_Fraction_blood, NA_real_),
    timepoint = coalesce(timepoint_info_BM, timepoint_info_blood),
    TF_category = case_when(
      is.na(Tumor_Fraction_BM) ~ "BM NA",
      is.na(Tumor_Fraction_blood) ~ "Blood NA",
      Tumor_Fraction_BM > 0.05 & Tumor_Fraction_blood > 0.05 ~ "High TF",
      !is.na(Tumor_Fraction_BM) & Tumor_Fraction_blood <= 0.05 ~ "Low TF",
      TRUE ~ "Discordant"
    )
  )
# For each sample, calculate CNA metrics
# Create the suffixed names for BM and blood calls from the merged dataset
BM_CNA_cols <- paste0(cna_cols, "_BM")
blood_CNA_cols <- paste0(cna_cols, "_blood")

# Now calculate CNA metrics using these new column names:
merged_CNA <- merged_CNA %>%
  rowwise() %>%
  mutate(
    n_CNA_BM = sum(c_across(all_of(BM_CNA_cols)) == "Yes", na.rm = TRUE),
    n_CNA_blood = sum(c_across(all_of(blood_CNA_cols)) == "Yes", na.rm = TRUE),
    n_CNA_intersect = sum((c_across(all_of(BM_CNA_cols)) == "Yes") &
                            (c_across(all_of(blood_CNA_cols)) == "Yes"), na.rm = TRUE),
    n_CNA_union = sum((c_across(all_of(BM_CNA_cols)) == "Yes") |
                        (c_across(all_of(blood_CNA_cols)) == "Yes"), na.rm = TRUE),
    CNA_overlap = ifelse(n_CNA_union == 0, 1, n_CNA_intersect / n_CNA_union),
    CNA_BM_only = n_CNA_BM - n_CNA_intersect,
    CNA_blood_only = n_CNA_blood - n_CNA_intersect
  ) %>%
  ungroup()

CNA_summary <- merged_CNA %>%
  group_by(timepoint, TF_category) %>%
  summarize(
    n_samples = n(),
    total_CNA_BM = sum(n_CNA_BM, na.rm = TRUE),
    total_CNA_blood = sum(n_CNA_blood, na.rm = TRUE),
    total_CNA_intersect = sum(n_CNA_intersect, na.rm = TRUE),
    avg_CNA_overlap = mean(CNA_overlap, na.rm = TRUE),
    total_CNA_BM_only = sum(CNA_BM_only, na.rm = TRUE),
    total_CNA_blood_only = sum(CNA_blood_only, na.rm = TRUE),
    false_positive_rate = ifelse(total_CNA_blood == 0, NA, 100 * total_CNA_blood_only / total_CNA_blood),
    avg_TF_BM = mean(Tumor_Fraction_BM, na.rm = TRUE),
    avg_TF_blood = mean(Tumor_Fraction_blood, na.rm = TRUE)
  ) %>% ungroup()

cat("\n### CNA Summary by Timepoint and Tumor Fraction Category ###\n")
print(CNA_summary)

# --- Step 8. Generate Summary Paragraphs -------------------------------------
# Example paragraphs (you can further modify the text based on your desired style)

# Mutation paragraph example:
mutation_paragraph <- with(mutation_summary, {
  # For each timepoint and TF category, report:
  paste0("For BM mutations at ", timepoint, " with ", TF_category, " samples (n=", n_samples, 
         "): Out of a total of ", total_BM_mut, " mutations identified in BM, ", total_intersect, 
         " were also detected in cfDNA (concordance = ", round(avg_overlap*100,1), 
         "%). There were ", total_blood_only, " mutations found exclusively in cfDNA, representing a false positive rate of ", 
         round(false_positive_rate,1), "%. The average tumor fraction was ", round(avg_TF_BM*100,1),
         "% in BM and ", round(avg_TF_blood*100,1), "% in cfDNA."
  )
})

cat("\n--- Mutation Summary Paragraphs ---\n")
cat(paste(mutation_paragraph, collapse = "\n"), "\n\n")

# Translocation paragraph example:
trans_paragraph <- with(trans_summary, {
  paste0("For translocations at ", timepoint, " with ", TF_category, " samples (n=", n_samples, 
         "): Out of ", total_trans_BM, " translocation calls in BM, ", total_trans_intersect, 
         " were concordantly detected in cfDNA (concordance = ", round(avg_trans_overlap*100,1), 
         "%). There were ", total_trans_blood_only, " cfDNA-only calls (false positive rate = ", 
         round(false_positive_rate,1), "%). The average tumor fractions were ", round(avg_TF_BM*100,1),
         "% (BM) and ", round(avg_TF_blood*100,1), "% (cfDNA)."
  )
})

cat("--- Translocation Summary Paragraphs ---\n")
cat(paste(trans_paragraph, collapse = "\n"), "\n\n")

# CNA paragraph example:
CNA_paragraph <- with(CNA_summary, {
  paste0("For CNAs at ", timepoint, " with ", TF_category, " samples (n=", n_samples, 
         "): Out of ", total_CNA_BM, " CNA calls in BM, ", total_CNA_intersect, 
         " were also observed in cfDNA (concordance = ", round(avg_CNA_overlap*100,1), 
         "%). There were ", total_CNA_blood_only, " cfDNA-only calls (false positive rate = ", 
         round(false_positive_rate,1), "%). Average tumor fractions were ", round(avg_TF_BM*100,1),
         "% in BM and ", round(avg_TF_blood*100,1), "% in cfDNA."
  )
})

cat("--- CNA Summary Paragraphs ---\n")
cat(paste(CNA_paragraph, collapse = "\n"), "\n")


# --- Global Percentages for Baseline Samples with High Tumor Fraction ---
# Filter baseline samples from BM and cfDNA datasets
baseline_BM <- combined_data_heatmap_BM %>% filter(timepoint_info == "Baseline")
baseline_cfDNA <- combined_data_heatmap_blood %>% filter(timepoint_info == "Baseline")

# Calculate percentage of baseline BM samples with tumor fraction > 5%
percent_BM_highTF <- sum(baseline_BM$Tumor_Fraction > 0.05, na.rm = TRUE) / nrow(baseline_BM) * 100

# Calculate percentage of baseline cfDNA samples with tumor fraction > 5%
percent_cfDNA_highTF <- sum(baseline_cfDNA$Tumor_Fraction > 0.05, na.rm = TRUE) / nrow(baseline_cfDNA) * 100

cat("\nGlobally, at Baseline, ", round(percent_BM_highTF,1), 
    "% of BM samples and ", round(percent_cfDNA_highTF,1), "% of cfDNA samples have high tumor fractions (>5%).\n")

# --- Compute Median Tumor Fraction and Range at Baseline ---
median_TF_BM <- median(baseline_BM$Tumor_Fraction, na.rm = TRUE) * 100
range_TF_BM <- range(baseline_BM$Tumor_Fraction, na.rm = TRUE) * 100

median_TF_cfDNA <- median(baseline_cfDNA$Tumor_Fraction, na.rm = TRUE) * 100
range_TF_cfDNA <- range(baseline_cfDNA$Tumor_Fraction, na.rm = TRUE) * 100

# --- Print the Summary ---
cat("\nAt baseline, the median tumor fraction estimated by ichorCNA was ", 
    round(median_TF_BM,1), "% (range ", round(range_TF_BM[1],1), "-", round(range_TF_BM[2],1),
    "%) in BM and ", round(median_TF_cfDNA,1), "% (range ", round(range_TF_cfDNA[1],1), "-",
    round(range_TF_cfDNA[2],1), "%) in cfDNA.\n")


# Filter for patients with High Tumor Fraction (>5%) in BM and cfDNA
valid_genes_BM <- intersect(myeloma_genes, colnames(combined_data_heatmap_BM))
valid_genes_cfDNA <- intersect(myeloma_genes, colnames(combined_data_heatmap_blood))

cat("Valid genes in BM dataset:", paste(valid_genes_BM, collapse = ", "), "\n")
cat("Valid genes in cfDNA dataset:", paste(valid_genes_cfDNA, collapse = ", "), "\n")


# Filter for patients with High Tumor Fraction (>5%) in BM and cfDNA
high_TF_BM <- combined_data_heatmap_BM %>% filter(Tumor_Fraction > 0.05)
high_TF_cfDNA <- combined_data_heatmap_blood %>% filter(Tumor_Fraction > 0.05)

# Compute percentage of patients with at least one MM-associated mutation in BM
BM_high_TF_with_mut <- high_TF_BM %>%
  mutate(has_mutation = rowSums(select(., all_of(valid_genes_BM)) != "No Mutation") > 0) %>%
  summarise(percent_with_mut = mean(has_mutation, na.rm = TRUE) * 100) %>%
  pull(percent_with_mut)

# Compute percentage of patients with at least one MM-associated mutation in cfDNA
cfDNA_high_TF_with_mut <- high_TF_cfDNA %>%
  mutate(has_mutation = rowSums(select(., all_of(valid_genes_cfDNA)) != "No Mutation") > 0) %>%
  summarise(percent_with_mut = mean(has_mutation, na.rm = TRUE) * 100) %>%
  pull(percent_with_mut)

# Print the results
cat("In patients with tumor fractions of >5%,",
    round(BM_high_TF_with_mut, 1), "% demonstrated MM-associated mutations in BM, and",
    round(cfDNA_high_TF_with_mut, 1), "% in cfDNA.\n")


