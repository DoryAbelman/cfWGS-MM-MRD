# =============================================================================
# Script: 1_4_Process_CNA_Data.R
#
# Description:
#   This script imports, harmonizes, and summarizes copy-number segment (SEG)
#   calls from ichorCNA outputs, mapping them to chromosomal arms of interest
#   and generating a binary matrix of arm-level alterations for myeloma analyses.
#   Steps include:
#     1. Recursively list and read all “*.seg” files in seg_dir
#        – Extracts chr, start, end, and the sample’s “Corrected_Call” column
#        – Renames each call column to the sample name (basename without .seg)
#     2. Merges every SEG table by (chr, start, end) via full joins into
#        `combined_seg_data`
#     3. Loads cytoband definitions (cb_frame) and derives an “arm” label
#        (e.g. “1p”, “1q”, “13q”, “17p”) by overlapping each segment
#     4. Saves the wide format as:
#        – CSV “Oct_2024_combined_corrected_calls.csv”
#        – RDS “Oct_2024_combined_corrected_calls.rds”
#     5. Pivots to long format and uppercases calls for consistency
#     6. Defines a list of myeloma-relevant arms (del1p, amp1q,
#        del13q, del17p) and for each:
#        – Flags segments as altered if “Corrected_Call” matches gain/loss labels
#        – Computes the proportion of altered bins per sample
#        – Sets a binary indicator (1 if >33% of that arm is altered)
#     7. Adds hyperdiploidy status by checking full gains (>80%) across
#        chromosomes 3,5,7,9,11,15,19,21
#     8. Outputs `myeloma_CNA_matrix_with_HRD`, a Sample × feature binary table
#
# Inputs:
#   • seg_dir          = "Oct 2024 data/Ichor_CNA"          # ichorCNA *.seg files
#   • cb_frame         = data frame of hg38 cytobands
#
# Outputs:
#   • RDS: Oct_2024_combined_corrected_calls.rds
#   • CSV: Oct_2024_combined_corrected_calls.csv
#   • R object: myeloma_CNA_matrix_with_HRD (binary arm-level matrix)
#
# Dependencies:
#   library(tidyverse); library(purrr)
#   library(GenomicRanges)
#
# Usage:
#   source("1_4_Process_CNA_Data.R")
#   # creates combined_seg_data and myeloma_CNA_matrix_with_HRD in working dir
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================


# Load Libraries:


source("setup_packages.R")
source("config.R")
source("helpers.R")
## Load clinical info
# Load in the patient info 
metada_df_mutation_comparison <- read_csv("combined_clinical_data_updated_Feb5_2025.csv")

# Add a Tumor_Sample_Barcode column to metada_df_mutation_comparison
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Tumor_Sample_Barcode = Bam %>%
           # Remove _PG or _WG
           str_remove_all("_PG|_WG") %>%
           # Remove anything after ".filter", ".ded", or ".recalibrate"
           str_replace_all("\\.filter.*|\\.ded.*|\\.recalibrate.*", ""))

metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Bam_clean_tmp = gsub(".bam$", "", Bam))  # Remove the '.bam' suffix



#### Now get the ichorCNA data loaded 
# Set the directory containing the SEG files
seg_dir <- "Oct 2024 data/Ichor_CNA"

# Get a list of all SEG files in the directory
seg_files <- list.files(seg_dir, pattern = "*.seg", full.names = TRUE)

# Initialize an empty list to store data for each file
seg_data_list <- list()

# Loop over each SEG file and process it
for (file in seg_files) {
  
  # Read the SEG file into a dataframe
  seg_data <- read.delim(file, header = TRUE)
  
  # Extract the filename without the .seg extension
  sample_name <- gsub(".cna.seg$", "", basename(file))
  
  # Extract the relevant columns: chr, start, end, and the exact Corrected_Call column
  seg_corrected_call <- seg_data %>%
    dplyr::select(chr, start, end, Corrected_Call = matches("Corrected_Call$"))
  
  # Rename the Corrected_Call column to include the sample name
  colnames(seg_corrected_call)[4] <- paste0(sample_name)
  
  # Add the processed data to the list
  seg_data_list[[file]] <- seg_corrected_call
}


# Combine all the data into one dataframe by chr, start, and end
combined_seg_data <- purrr::reduce(seg_data_list, full_join, by = c("chr", "start", "end"))

## Add the arm info to the directory 
cb_frame_arm <- cb_frame %>%
  mutate(arm = paste0(gsub("chr", "", chr), substr(band, 1, 1))) %>%
  dplyr::select(chr, start, end, arm)  # Select only chr, start, end, and

cb_frame_arm <- cb_frame_arm %>%
  mutate(chr = gsub("chr", "", chr))  # edit to match other dataframe

# Function to determine if two ranges overlap
check_overlap <- function(chr, start1, end1, cb_frame_arm) {
  cb_matches <- cb_frame_arm %>%
    filter(chr == !!chr & start <= !!end1 & end >= !!start1)
  
  if (nrow(cb_matches) > 0) {
    return(cb_matches$arm[1])  # Return the first matching arm if there's an overlap
  } else {
    return(NA)  # Return NA if no overlap is found
  }
}

# Apply the function to add the arm column based on overlap
combined_seg_data <- combined_seg_data %>%
  rowwise() %>%
  mutate(arm = check_overlap(chr, start, end, cb_frame_arm)) %>%
  ungroup() %>%
  relocate(arm, .after = end)  # Move the arm column to the 4th position (after 'end')

# Save the combined dataframe as a CSV file for later use
write.csv(combined_seg_data, file = "Oct_2024_combined_corrected_calls.csv", row.names = FALSE)

# Save as RDS for efficient reloading later
saveRDS(combined_seg_data, file = "Oct_2024_combined_corrected_calls.rds")

## Transform data 
# Pivot the data to long format

combined_seg_data <- readRDS("Oct_2024_combined_corrected_calls.rds")

# Pivot the data (assuming combined_seg_data is already defined)
long_data <- combined_seg_data %>%
  pivot_longer(cols = -(1:4), names_to = "Sample", values_to = "Value") %>%
  mutate(Value = toupper(Value))


# Define the list of chromosomal arms with corrected chr values (no "chr" prefix)
myeloma_cna_arms <- list(
  del1p = list(chr = "1", arm = "1p"),    # Deletion in 1p
  amp1q = list(chr = "1", arm = "1q"),      # Amplification in 1q
  del13q = list(chr = "13", arm = "13q"),   # Deletion in 13q
  del17p = list(chr = "17", arm = "17p")     # Deletion in 17p
)

# Flatten the list into a data frame
myeloma_arms_df <- data.frame(
  Arm = names(myeloma_cna_arms),
  chr = sapply(myeloma_cna_arms, function(x) x$chr),
  arm = sapply(myeloma_cna_arms, function(x) x$arm),
  stringsAsFactors = FALSE
)

# Define gain and loss labels
gain_labels <- c("GAIN", "AMP", "HLAMP")
loss_labels <- c("HOMD", "HETD")

# Get the list of samples
samples <- unique(long_data$Sample)

# Initialize the results data frame
results <- data.frame(Sample = samples, stringsAsFactors = FALSE)

# Loop over each chromosomal arm
for (i in 1:nrow(myeloma_arms_df)) {
  arm_name <- myeloma_arms_df$Arm[i]
  chr_value <- myeloma_arms_df$chr[i]
  arm_value <- myeloma_arms_df$arm[i]
  
  # Filter the data for the specific chromosome and arm
  data_arm <- long_data %>%
    filter(chr == chr_value, arm == arm_value)
  
  # Debug: print dimensions and unique values for one arm (optional)
  # print(dim(data_arm))
  # print(unique(data_arm$Value))
  
  # Determine if it's an amplification or deletion
  if (startsWith(arm_name, "amp")) {
    data_arm <- data_arm %>%
      mutate(IsAltered = ifelse(Value %in% gain_labels, 1, 0))
  } else if (startsWith(arm_name, "del")) {
    data_arm <- data_arm %>%
      mutate(IsAltered = ifelse(Value %in% loss_labels, 1, 0))
  } else {
    next
  }
  
  # Calculate the proportion of altered segments per sample
  sample_proportions <- data_arm %>%
    group_by(Sample) %>%
    summarise(ProportionAltered = mean(IsAltered))
  
  # Ensure all samples are included (set missing samples to 0)
  missing_samples <- setdiff(samples, sample_proportions$Sample)
  if (length(missing_samples) > 0) {
    sample_proportions <- bind_rows(
      sample_proportions,
      data.frame(Sample = missing_samples, ProportionAltered = 0, stringsAsFactors = FALSE)
    )
  }
  
  # Merge with results and assign 1 if proportion > 1/3
  results <- left_join(results, sample_proportions, by = "Sample")
  results[[arm_name]] <- ifelse(results$ProportionAltered > 1/3, 1, 0)
  results <- results %>% select(-ProportionAltered)
}

# View results to check the output
print(results)

myeloma_CNA_matrix_with_HRD <- results



## Add hyperdiploid info:
## Edit based on Suzanne new criteria 
## Define hyperdiploid chromosomes (full gains required for all eight)
myeloma_cna_HRD <- list(
  hyperdiploid_chr3  = list(chr = "3"),
  hyperdiploid_chr5  = list(chr = "5"),
  hyperdiploid_chr7  = list(chr = "7"),
  hyperdiploid_chr9  = list(chr = "9"),
  hyperdiploid_chr11 = list(chr = "11"),
  hyperdiploid_chr15 = list(chr = "15"),
  hyperdiploid_chr19 = list(chr = "19"),
  hyperdiploid_chr21 = list(chr = "21")
)

# Flatten the hyperdiploid list into a data frame
hyperdiploid_arms_df <- data.frame(
  Chr = names(myeloma_cna_HRD),
  chr = sapply(myeloma_cna_HRD, function(x) x$chr),
  stringsAsFactors = FALSE
)

# Loop over hyperdiploid chromosomes and assess full gains (threshold > 90%)
for (i in 1:nrow(hyperdiploid_arms_df)) {
  chr_value <- hyperdiploid_arms_df$chr[i]
  
  # Filter the data for the specific chromosome
  data_chr <- long_data %>%
    filter(chr == chr_value)
  
  # Mark segments as gained if they match the gain labels
  data_chr <- data_chr %>%
    mutate(IsGain = ifelse(Value %in% gain_labels, 1, 0))
  
  # Calculate the proportion of segments gained per sample for this chromosome
  sample_gains <- data_chr %>%
    group_by(Sample) %>%
    summarise(ProportionGained = mean(IsGain))
  
  # Ensure all samples are included (assign 0 if no data)
  missing_samples <- setdiff(samples, sample_gains$Sample)
  if (length(missing_samples) > 0) {
    sample_gains <- bind_rows(
      sample_gains,
      data.frame(Sample = missing_samples, ProportionGained = 0)
    )
  }
  
  # Merge with the results table and assign gain indicator if >80% of segments are gained
  results <- left_join(results, sample_gains, by = "Sample")
  results[[hyperdiploid_arms_df$Chr[i]]] <- ifelse(results$ProportionGained > 0.8, 1, 0)
  results <- results %>% select(-ProportionGained)
}

# Mark a sample as hyperdiploid only if all eight hyperdiploid chromosomes have full gains
results <- results %>%
  mutate(
    hyperdiploid = if_else(
      rowSums(select(., starts_with("hyperdiploid_"))) == 8,
      1, 0
    )
  )

myeloma_CNA_matrix_with_HRD <- results


# Prepare CNA data
# We'll select relevant columns and set Sample as row names
cna_data <- myeloma_CNA_matrix_with_HRD %>%
  dplyr::select(Sample, del1p, amp1q, del13q, del17p, hyperdiploid) %>%
  mutate_at(vars(-Sample), as.character)  # Convert numeric to character for consistency



## Export the important tables 
# Define the directory for export
export_dir <- "Jan2025_exported_data"

# Check if the directory exists, if not, create it
if (!dir.exists(export_dir)) {
  dir.create(export_dir)
}

# Exporting dataframes as RDS files
saveRDS(cna_data, file = file.path(export_dir, "cna_data.rds"))

# Exporting dataframes as text files
write.table(cna_data, file = file.path(export_dir, "cna_data.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

