# =============================================================================
# Script: process_cna_seg_files.R
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
#   source("process_cna_seg_files.R")
#   # creates combined_seg_data and myeloma_CNA_matrix_with_HRD in working dir
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================


# Load Libraries:
library(tidyverse)       # dplyr, purrr, tidyr, readr, etc.
library(purrr)           # for reduce()
library(tidyr)           # for pivot_longer()
library(GenomicRanges)
library(IRanges)
library(S4Vectors)


## Load clinical info
# Load in the patient info 
metada_df_mutation_comparison <- read_csv("combined_clinical_data_updated_April2025.csv")

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



## Now get FISH concordance with the specific probes
# =====================================================================
# Purpose: Recompute per-sample CNA calls at FISH probe cytobands
#          using ichorCNA segments, independent of other scripts.
# Output:  probe_calls_bin_cytoband (in memory) + optional exports
# =====================================================================

fish_probe_xlsx    <- "Clinical data/FISH probe locations.xlsx" # must contain Chromosome, Target, Location in hg38
# EITHER provide cb_frame as an RDS with columns chr,start,end,band
cb_frame_rds       <- NULL  # e.g., "cb_frame_hg38.rds"  # set to a path OR leave NULL
# OR provide a UCSC cytoband txt (no header), which we'll parse
cytoband_txt       <- "cytoband.txt"  # fallback if cb_frame_rds is NULL



## ---- Tunables & label maps ----
pad_bp  <- 150000L  # ±150 kb padding around band windows
min_bp  <- 1000L    # require >= 1 kb overlap to count
max_nearest_bp  <- 10e6L     # allow nearest-segment fallback within 10 Mb 

# FISH feature mapping by arm label in fish table
feature_map <- c("1p" = "del1p", "1q" = "amp1q", "17p" = "del17p", "13q" = "del13q")
# Direction map for binning
dir_map     <- c(amp1q = "gain", del1p = "loss", del13q = "loss", del17p = "loss")

# Call label sets
gain_labels <- c("GAIN", "AMP", "HLAMP")
loss_labels <- c("HOMD", "HETD")

if (!file.exists(fish_probe_xlsx)) stop("Missing: ", fish_probe_xlsx)
fish_probe_locations <- readxl::read_excel(fish_probe_xlsx)

# Cytobands: prefer cb_frame RDS if provided; otherwise parse UCSC txt
if (!is.null(cb_frame_rds)) {
  if (!file.exists(cb_frame_rds)) stop("Missing: ", cb_frame_rds)
  cb_frame <- readRDS(cb_frame_rds)
  cyto <- cb_frame %>%
    transmute(
      chrom = gsub("^chr","", as.character(chr)),
      start = as.integer(start),
      end   = as.integer(end),
      band  = as.character(band)
    ) %>% arrange(chrom, start)
} else {
  if (!file.exists(cytoband_txt)) stop("Missing: ", cytoband_txt)
  cyto <- readr::read_tsv(
    cytoband_txt,
    col_names = c("chrom", "start", "end", "band", "gieStain"),
    col_types = cols(
      chrom = col_character(),
      start = col_double(),
      end   = col_double(),
      band  = col_character(),
      gieStain = col_character()
    ),
    progress = FALSE
  ) %>%
    mutate(
      chrom = gsub("^chr","", chrom),
      start = as.integer(start),
      end   = as.integer(end)
    ) %>%
    arrange(chrom, start)
}
stopifnot(all(c("chrom","start","end","band") %in% names(cyto)))

## ---- Build GRanges from ichorCNA segments (wide) ----
# Expect first 4 columns: chr, start, end, arm (your ichor pipeline wrote this)
need_cols <- c("chr","start","end","arm")
if (!all(need_cols %in% names(combined_seg_data))) {
  stop("combined_seg_data must have columns: ", paste(need_cols, collapse = ", "))
}

sample_cols_seg <- setdiff(names(combined_seg_data), need_cols)
if (length(sample_cols_seg) == 0L) stop("No sample columns found in combined_seg_data")

# GRanges of all merged segments
seg_gr <- GRanges(
  seqnames = as.character(combined_seg_data$chr),
  ranges   = IRanges(as.integer(combined_seg_data$start),
                     as.integer(combined_seg_data$end))
)

# Attach calls (one column per sample) to mcols; uppercase for consistency
calls_df <- combined_seg_data[, sample_cols_seg] %>%
  mutate(across(everything(), ~ toupper(as.character(.))))
if (anyDuplicated(names(calls_df))) {
  names(calls_df) <- make.unique(names(calls_df), sep = "_dup")
}
mcols(seg_gr) <- S4Vectors::DataFrame(as.list(calls_df), check.names = FALSE)
samples <- colnames(mcols(seg_gr))

## ---- Helpers: band parsing & expansion ----
.band_span <- function(cyto_chr, arm, band_tag) {
  tag <- paste0(arm, band_tag)  # e.g., "q21", "p13.1"
  rows <- which(startsWith(cyto_chr$band, tag))
  if (!length(rows)) return(NULL)
  tibble(start = min(cyto_chr$start[rows]),
         end   = max(cyto_chr$end[rows]))
}

.parse_target <- function(s) {
  s <- gsub("\\s", "", s)
  m <- regexec("^([0-9XYM]+)([pq])([0-9]+(?:\\.[0-9]+)?)(?:-([pq])?([0-9]+(?:\\.[0-9]+)?))?$", s)
  g <- regmatches(s, m)[[1]]
  if (!length(g)) return(NULL)
  chr   <- g[2]
  arm1  <- g[3]; band1 <- g[4]
  arm2  <- g[5]; band2 <- g[6]
  if (is.na(arm2) || identical(arm2, "")) arm2 <- arm1
  if (is.na(band2) || identical(band2, "")) {
    list(chr = chr, arm1 = arm1, band1 = band1, arm2 = arm1, band2 = band1, is_range = FALSE)
  } else {
    list(chr = chr, arm1 = arm1, band1 = band1, arm2 = arm2, band2 = band2, is_range = TRUE)
  }
}

## ---- Turn FISH Targets into genomic windows via cytobands ----
stopifnot(all(c("Chromosome","Target") %in% names(fish_probe_locations)))

probe_windows_cytoband <- fish_probe_locations %>%
  transmute(
    feature = dplyr::recode(Chromosome, !!!feature_map, .default = NA_character_),
    Target  = Target
  ) %>%
  rowwise() %>%
  mutate(.pt = list(.parse_target(Target))) %>%
  ungroup() %>%
  filter(!purrr::map_lgl(.pt, is.null), !is.na(feature)) %>%
  mutate(
    chr  = purrr::map_chr(.pt, ~ .x$chr),
    arm1 = purrr::map_chr(.pt, ~ .x$arm1),
    b1   = purrr::map_chr(.pt, ~ .x$band1),
    arm2 = purrr::map_chr(.pt, ~ .x$arm2),
    b2   = purrr::map_chr(.pt, ~ .x$band2)
  ) %>%
  select(-.pt) %>%
  group_by(feature, chr, arm1, b1, arm2, b2) %>%
  reframe({
    cy_chr <- dplyr::filter(cyto, chrom == chr) %>% arrange(start)
    sspan  <- .band_span(cy_chr, arm1, b1)
    espan  <- .band_span(cy_chr, arm2, b2)
    if (is.null(sspan) || is.null(espan)) tibble(start = NA_integer_, end = NA_integer_) else
      tibble(start = min(sspan$start, espan$start),
             end   = max(sspan$end,   espan$end))
  }) %>%
  ungroup() %>%
  filter(!is.na(start), !is.na(end)) %>%
  mutate(
    start_pad = pmax(1L, start - pad_bp),
    end_pad   = end + pad_bp
  ) %>%
  distinct(feature, chr, start_pad, end_pad, .keep_all = TRUE)

if (!nrow(probe_windows_cytoband)) stop("No valid probe windows built from Target field")

## ---- Overlap probe windows with segments ----
probe_gr_cyto <- GRanges(
  seqnames = probe_windows_cytoband$chr,
  ranges   = IRanges(probe_windows_cytoband$start_pad, probe_windows_cytoband$end_pad),
  feature  = probe_windows_cytoband$feature
)

hits_cyto <- findOverlaps(probe_gr_cyto, seg_gr, ignore.strand = TRUE)

ovl_bp_cyto <- if (length(hits_cyto)) width(pintersect(
  ranges(probe_gr_cyto)[queryHits(hits_cyto)],
  ranges(seg_gr)[subjectHits(hits_cyto)]
)) else integer(0)

hits_df_cyto <- tibble::tibble(
  probe_idx = as.integer(queryHits(hits_cyto)),
  seg_idx   = as.integer(subjectHits(hits_cyto)),
  ovl       = as.integer(ovl_bp_cyto),
  src       = "overlap"
) %>%
  dplyr::filter(ovl >= min_bp)

# --- Fallback: nearest segment per probe when no overlaps made it, since ichorCNA skips problematic regions of the genome
# distanceToNearest returns the *closest* segment (ties broken arbitrarily), with distance in bp.
nearest_hits <- GenomicRanges::distanceToNearest(probe_gr_cyto, seg_gr, ignore.strand = TRUE)

nearest_df <- tibble::tibble(
  probe_idx = as.integer(S4Vectors::queryHits(nearest_hits)),
  seg_idx   = as.integer(S4Vectors::subjectHits(nearest_hits)),
  dist_bp   = as.integer(S4Vectors::mcols(nearest_hits)$distance)
) %>%
  # keep only probes that have no overlapping rows already
  dplyr::anti_join(hits_df_cyto %>% dplyr::select(probe_idx) %>% dplyr::distinct(),
                   by = "probe_idx") %>%
  # enforce a maximum distance cap to avoid silly jumps
  dplyr::filter(dist_bp <= max_nearest_bp) %>%
  # encode "overlap score" as negative distance so sorting still prefers overlaps
  dplyr::transmute(
    probe_idx,
    seg_idx,
    ovl = -dist_bp,   # negative = farther; less negative = closer
    src = "nearest"
  )

# Combine, order by probe then "best support" (overlap first, then nearest by distance)
hits_df_cyto <- dplyr::bind_rows(hits_df_cyto, nearest_df) %>%
  dplyr::arrange(probe_idx, dplyr::desc(ovl))

## ---- Choose per-sample call by severity first, then overlap ----
probe_ids_cyto <- seq_along(probe_gr_cyto)
final_calls_mat_cyto <- matrix(
  NA_character_, nrow = length(probe_ids_cyto), ncol = length(samples),
  dimnames = list(NULL, samples)
)

gain_priority <- c("HLAMP","AMP","GAIN")
loss_priority <- c("HOMD","HETD","LOSS","CNLOH")

feature_vec_cyto <- as.character(mcols(probe_gr_cyto)$feature)

for (p in probe_ids_cyto) {
  seg_rows <- hits_df_cyto %>% filter(probe_idx == p)
  segs     <- seg_rows$seg_idx
  if (!length(segs)) next
  
  seg_calls <- as.matrix(mcols(seg_gr)[segs, samples, drop = FALSE])
  ovl_here  <- seg_rows$ovl
  
  feat <- feature_vec_cyto[p]
  dir  <- unname(dir_map[feat])   # "gain" or "loss"
  pri  <- if (identical(dir, "gain")) gain_priority else loss_priority
  
  lab <- toupper(seg_calls)  # same shape as seg_calls
  sev_rank <- array(
    match(lab, pri, nomatch = length(pri) + 1L),  # unknown => worst rank
    dim       = dim(lab),
    dimnames  = dimnames(lab)
  )
  
  for (j in seq_along(samples)) {
    calls_j <- seg_calls[, j]
    valid   <- which(!is.na(calls_j) & nzchar(calls_j))
    if (!length(valid)) next
    
    rnk  <- sev_rank[valid, j]
    best <- which(rnk == min(rnk))
    if (length(best) > 1L) best <- best[which.max(ovl_here[valid][best])]
    final_calls_mat_cyto[p, j] <- calls_j[valid][best]
  }
}

## ---- Tidy & bin by expected direction ----
probe_meta_cyto <- tibble(
  probe_idx = probe_ids_cyto,
  feature   = as.character(mcols(probe_gr_cyto)$feature)
)

probe_calls_long_cyto <- as_tibble(final_calls_mat_cyto, .name_repair = "minimal") %>%
  mutate(probe_idx = probe_ids_cyto) %>%
  pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "probe_call") %>%
  left_join(probe_meta_cyto, by = "probe_idx") %>%
  select(Sample, feature, probe_call)

probe_calls_bin_cytoband <- probe_calls_long_cyto %>%
  mutate(
    direction = unname(dir_map[feature]),
    is_altered_at_probe = if_else(
      direction == "gain",
      as.integer(probe_call %in% gain_labels),
      as.integer(probe_call %in% loss_labels)
    )
  ) %>%
  select(-direction) %>%
  pivot_wider(
    names_from  = feature,
    values_from = c(probe_call, is_altered_at_probe),
    names_sep   = "_"
  )

# View head
print(head(probe_calls_bin_cytoband, 10))

## ---- Optional exports ----
export_dir <- "Jan2025_exported_data"
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
saveRDS(probe_calls_bin_cytoband, file = file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.rds"))
write.table(probe_calls_bin_cytoband,
            file = file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


### Now continue with regular data transformations

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
  
  # Merge with the results table and assign gain indicator if >65% of segments are gained
  results <- left_join(results, sample_gains, by = "Sample")
  results[[hyperdiploid_arms_df$Chr[i]]] <- ifelse(results$ProportionGained > 0.65, 1, 0)
  results <- results %>% select(-ProportionGained)
}

# Mark a sample as hyperdiploid only if all eight hyperdiploid chromosomes have full gains
results <- results %>%
  mutate(
    hyperdiploid = if_else(
      rowSums(select(., starts_with("hyperdiploid_"))) >= 5,
      1, 0
    )
  )

myeloma_CNA_matrix_with_HRD <- results


# Prepare CNA data
# We'll select relevant columns and set Sample as row names
cna_data <- myeloma_CNA_matrix_with_HRD %>%
  dplyr::select(Sample, del1p, amp1q, del13q, del17p, hyperdiploid) %>%
  mutate_at(vars(-Sample), as.character)  # Convert numeric to character for consistency

cna_data_backup <- cna_data

## Export the important tables 
# Define the directory for export
export_dir <- "Jan2025_exported_data"

# Check if the directory exists, if not, create it
if (!dir.exists(export_dir)) {
  dir.create(export_dir)
}

# Exporting dataframes as RDS files
saveRDS(cna_data, file = file.path(export_dir, "cna_data_ichorCNA.rds"))

# Exporting dataframes as text files
write.table(cna_data, file = file.path(export_dir, "cna_data_ichorCNA.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cna_data_backup, file = file.path(export_dir, "cna_data_backup.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

