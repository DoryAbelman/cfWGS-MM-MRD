# =============================================================================
# Script: process_sequenza_cna_seg_files.R
#
# Description:
#   End-to-end arm-level CNA caller for Sequenza segments with **ploidy-aware**
#   labels and a realistic **hyperdiploidy** definition for myeloma.
#
#   What it does:
#     1) Reads all Sequenza “*_segments.txt(.gz)” in `seg_dir`
#        – Requires columns: chromosome, start.pos, end.pos, CNt, A, B
#     2) Estimates per-sample baseline ploidy (rounded) from autosomes using a
#        length-weighted mode of integer CNt; saves a ploidy table
#     3) Converts each segment to a **baseline-aware categorical call**
#        (LOSS / GAIN / HLAMP; plus HOMD/HETD; optional CNLOH at baseline with B=0)
#     4) Merges all samples on (chr, start, end) into `combined_seg_data`
#     5) Maps segments to chromosomal **arms** (1p/1q, …) via cytobands using
#        GRanges and assigns the arm with **maximum bp overlap**
#     6) Precomputes **fixed denominators**:
#        – `arm_lengths` for each p/q arm (bp)
#        – `chr_lengths` for whole chromosomes (bp)
#     7) Computes **length-weighted altered fraction per arm** using fixed arm
#        lengths as denominators, then binarizes:
#        – del1p, amp1q, del13q, del17p (default threshold: > 1/3 of arm length)
#     8) Calls **hyperdiploidy** using relaxed myeloma-appropriate criteria:
#        – A chromosome is considered gained if > `gain_thresh_chr` (default 0.65)
#          of its length is GAIN/HLAMP
#        – A sample is hyperdiploid if ≥ `k_trisomies` (default 5) of
#          {3,5,7,9,11,15,19,21} are gained
#     9) Exports:
#        – Combined per-segment calls: CSV/RDS
#        – Per-sample ploidy estimates: CSV
#        – Final Sample × feature matrix (`cna_data`) with
#          {del1p, amp1q, del13q, del17p, hyperdiploid}: RDS/TXT
#    10) Prints cohort-level frequencies for quick QC
#
# Inputs:
#   • seg_dir   = directory with Sequenza segment files (e.g., “…/All_Segments_400/”)
#   • cb_frame  = hg19/hg38 cytobands data.frame with columns: chr, start, end, band
#                 (chr names like “chr1”, “chr2”, …; script harmonizes to 1..22,X,Y)
#
# Outputs (example file names shown for γ=400 run, adjust to your run/date):
#   • CSV/RDS:  "Sep_2025_combined_sequenza_calls_400.csv/.rds"
#   • CSV:      "Sep_2025_sequenza_ploidy_estimates.csv"
#   • RDS/TXT:  "Jan2025_exported_data/cna_data_from_sequenza_400.rds/.txt"
#   • Console:  one-row tibble with cohort proportions for each feature
#
# Key parameters (tunable):
#   • Arm binary threshold:        > 1/3 of arm length (per-arm override optional)
#   • Hyperdiploid per-chr gain:   gain_thresh_chr = 0.65 (0.60–0.70 typical)
#   • Hyperdiploid sample rule:    k_trisomies = 5 (≥ 4–5 common in MM)
#
# Assumptions & notes:
#   • Calls are **baseline-aware** per sample (neutral ≈ round(ploidy)).
#     Use `gain_labels = {GAIN, HLAMP}`; `loss_labels = {LOSS, HETD, HOMD}`.
#   • Denominators are **fixed** (full arm/chrom lengths) to avoid NA inflation.
#   • Genome build of `cb_frame` must match BAM/Sequenza build.
#   • Recommended Sequenza segmentation for MM: γ ≈ 300–500; this run uses γ=400.
#
# Dependencies:
#   library(tidyverse)  # dplyr, tidyr, readr, purrr
#   library(GenomicRanges)
#
# Usage:
#   source("process_sequenza_cna_seg_files.R")
#   # Produces: combined_seg_data, sample_ploidy, results, cna_data, and exports
#
# Author: Dory Abelman
# Date:   2025-10-06
# =============================================================================



# Load Libraries:
library(tidyverse)       # dplyr, purrr, tidyr, readr, etc.
library(purrr)           # for reduce()
library(tidyr)           # for pivot_longer()
library(GenomicRanges)


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


seg_dir <- "Oct 2024 data/Sequenza/All_Segments_400/"

# ---- helper: derive ichor-like categorical calls from Sequenza fields
call_from_CNt_AB <- function(CNt, A, B) {
  case_when(
    is.na(CNt)               ~ NA_character_,
    CNt == 0                 ~ "HOMD",
    CNt == 1                 ~ "HETD",
    CNt == 2 & !is.na(B) & B == 0 ~ "CNLOH",        # not used downstream
    CNt == 3                 ~ "GAIN",
    CNt >= 4 & CNt < 6       ~ "AMP",
    CNt >= 6                 ~ "HLAMP",
    TRUE                     ~ "NEUT"
  )
}

# per-sample baseline ploidy-aware caller
call_from_CNt_baseline <- function(CNt, A, B, baseline) {
  case_when(
    is.na(CNt)               ~ NA_character_,
    CNt == 0                 ~ "HOMD",
    CNt == 1                 ~ "HETD",
    CNt >= baseline + 3      ~ "HLAMP",
    CNt >  baseline          ~ "GAIN",
    CNt <  baseline          ~ "LOSS",
    # optional: CNLOH when total equals baseline but B==0
    CNt == baseline & !is.na(B) & B == 0 ~ "CNLOH",
    TRUE                      ~ "NEUT"
  )
}

# length-weighted mode for discrete values
wmode <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  t <- tapply(w[ok], x[ok], sum, na.rm = TRUE)
  as.numeric(names(t)[which.max(t)])
}

# robust weighted mode on integers with dominance score
robust_ploidy <- function(CNt, start, end, chr, min_seg_bp = 5e6) {
  # filter autosomes + long enough segments
  ok <- chr %in% as.character(1:22)
  L  <- pmax(end - start + 1, 1)
  CNt_round <- round(CNt)
  
  ok <- ok & !is.na(CNt_round) & (L >= min_seg_bp)
  if (!any(ok)) return(list(ploidy_est = NA_real_, dominance = NA_real_))
  
  # length-weighted counts per integer CN
  tab <- tapply(L[ok], CNt_round[ok], sum, na.rm = TRUE)
  tab <- tab[order(as.numeric(names(tab)))]
  # winner and dominance (mass of winner / total mass)
  winner <- as.numeric(names(tab)[which.max(tab)])
  dom    <- as.numeric(max(tab) / sum(tab))
  
  list(ploidy_est = winner, dominance = dom)
}


# ---- read Sequenza *_segments.txt (optionally .gz), rename cols, make calls
seg_files <- list.files(seg_dir, pattern = "_segments\\.txt(\\.gz)?$", full.names = TRUE)

seg_data_list <- vector("list", length(seg_files))
names(seg_data_list) <- basename(seg_files)

ploidy_list <- vector("list", length(seg_files))   


for (i in seq_along(seg_files)) {
  f <- seg_files[i]
  df <- readr::read_tsv(f, show_col_types = FALSE)
  
  # required columns in Sequenza segments
  req <- c("chromosome", "start.pos", "end.pos", "CNt", "A", "B")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", basename(f), ": ", paste(miss, collapse = ", "))
  
  sample_name <- sub("_segments\\.txt(\\.gz)?$", "", basename(f))
  
  # ----- estimate baseline (rounded ploidy) from autosomes, length-weighted mode of rounded CNt
  df_auto <- df %>%
    transmute(
      chr   = stringr::str_remove(chromosome, "^chr"),
      start = as.numeric(start.pos),
      end   = as.numeric(end.pos),
      CNt   = as.numeric(CNt)
    ) %>%
    filter(chr %in% as.character(1:22)) %>%
    mutate(seg_len = pmax(end - start + 1, 1),
           CNt_round = round(CNt))
  
  est_ploidy <- wmode(df_auto$CNt_round, df_auto$seg_len)  # assumes wmode() defined earlier
  baseline   <- ifelse(is.na(est_ploidy), NA_real_, pmax(1, round(est_ploidy)))
  
  # record ploidy for reporting
  ploidy_list[[i]] <- tibble(
    Sample     = sample_name,
    ploidy_est = est_ploidy,
    baseline   = baseline
  )
  
  # ----- build per-segment categorical calls (baseline-aware; fallback to original if baseline is NA)
  seg_core <- df %>%
    transmute(
      chr   = stringr::str_remove(chromosome, "^chr"),
      start = as.numeric(start.pos),
      end   = as.numeric(end.pos),
      CNt   = as.numeric(CNt),
      A     = as.numeric(A),
      B     = as.numeric(B)
    )
  
  if (is.na(baseline)) {
    # fallback: your original diploid-centric rule if baseline couldn't be estimated
    seg_calls <- seg_core %>%
      mutate(Call = call_from_CNt_AB(CNt, A, B))
  } else {
    # baseline-aware labels
    seg_calls <- seg_core %>%
      mutate(Call = case_when(
        is.na(CNt)                 ~ NA_character_,
        CNt == 0                   ~ "HOMD",
        CNt == 1                   ~ "HETD",
        CNt >= baseline + 3        ~ "HLAMP",
        CNt >  baseline            ~ "GAIN",
        CNt <  baseline            ~ "LOSS",
        CNt == baseline & !is.na(B) & B == 0 ~ "CNLOH",  # optional: keep CNLOH tag
        TRUE                       ~ "NEUT"
      ))
  }
  
  seg_df <- seg_calls %>% select(chr, start, end, Call)
  
  # stash with sample-specific column name
  colnames(seg_df)[4] <- sample_name
  seg_data_list[[i]] <- seg_df
}

# ---- combine and ploidy table (unchanged from your flow)
combined_seg_data <- purrr::reduce(seg_data_list, full_join, by = c("chr","start","end"))

sample_ploidy <- bind_rows(ploidy_list) %>%
  arrange(Sample)


print(sample_ploidy, n = 50)
write.csv(sample_ploidy, "Sep_2025_sequenza_ploidy_estimates.csv", row.names = FALSE)


# Harmonize chr labels in segments
combined_seg_data <- combined_seg_data %>%
  mutate(chr = toupper(as.character(chr))) %>%
  mutate(chr = gsub("^CHR", "", chr)) %>%
  mutate(chr = dplyr::recode(chr, `23`="X", `24`="Y"))

valid_chr <- c(as.character(1:22), "X", "Y")
combined_seg_data <- combined_seg_data %>% filter(chr %in% valid_chr)

# Ensure cytoband column is named 'band' (some tables use 'name')
if (!"band" %in% names(cb_frame) && "name" %in% names(cb_frame)) {
  cb_frame <- cb_frame %>% rename(band = name)
}
stopifnot(all(c("chr","start","end","band") %in% names(cb_frame)))

# Collapse cytobands to one p- and one q-range per chromosome
cb_arms <- cb_frame %>%
  mutate(chr = gsub("^chr", "", chr)) %>%
  filter(chr %in% valid_chr) %>%
  mutate(arm_letter = substr(band, 1, 1)) %>%
  filter(arm_letter %in% c("p","q")) %>%
  group_by(chr, arm_letter) %>%
  summarise(start = min(start), end = max(end), .groups = "drop") %>%
  mutate(arm = paste0(chr, arm_letter))

# Arm lengths (bp)
arm_lengths <- cb_arms %>%
  mutate(arm_len = end - start + 1) %>%
  rename(arm_start = start, arm_end = end) %>%
  select(chr, arm, arm_start, arm_end, arm_len)

# Chromosome lengths (p+q)
chr_lengths <- cb_arms %>%
  group_by(chr) %>%
  summarise(chr_len = max(end) - min(start) + 1, .groups = "drop")

# Build GRanges
arms_gr <- GRanges(
  seqnames = cb_arms$chr,
  ranges   = IRanges(cb_arms$start, cb_arms$end),
  arm      = cb_arms$arm
)
seg_gr <- GRanges(
  seqnames = combined_seg_data$chr,
  ranges   = IRanges(combined_seg_data$start, combined_seg_data$end)
)

# Overlap and assign arm with maximum basepair overlap
hits <- findOverlaps(seg_gr, arms_gr, ignore.strand = TRUE)
if (length(hits) == 0L) stop("No overlaps between segments and cytoband arms. Check genome build / chr labels.")

ov_w <- width(pintersect(ranges(seg_gr)[queryHits(hits)],
                         ranges(arms_gr)[subjectHits(hits)]))

best_by_row <- tapply(
  X     = seq_along(ov_w),
  INDEX = queryHits(hits),
  FUN   = function(idx) {
    arm_ids <- subjectHits(hits)[idx]
    mcols(arms_gr)$arm[ arm_ids[ which.max(ov_w[idx]) ] ]
  }
)

arm_vec <- rep(NA_character_, nrow(combined_seg_data))
arm_vec[as.integer(names(best_by_row))] <- unname(unlist(best_by_row))

combined_seg_data <- combined_seg_data %>%
  mutate(arm = arm_vec) %>%
  relocate(arm, .after = end)

combined_seg_data <- combined_seg_data %>%
  left_join(arm_lengths, by = c("chr","arm")) %>%
  mutate(
    seg_len_in_arm = pmax(pmin(end, arm_end) - pmax(start, arm_start) + 1, 0),
    seg_len_in_arm = ifelse(is.na(seg_len_in_arm), 0, seg_len_in_arm)
  )

# (Optional) quick sanity checks
# combined_seg_data %>% summarize(prop_assigned = mean(!is.na(arm)))
# combined_seg_data %>% filter(is.na(arm)) %>% count(chr, sort = TRUE)

# ---- save combined (optional)
write.csv(combined_seg_data, "Sep_2025_combined_sequenza_calls_400.csv", row.names = FALSE)
saveRDS(combined_seg_data,  "Sep_2025_combined_sequenza_calls_400.rds")

# ---- long format with segment length (for length-weighting)
sample_cols <- setdiff(
  names(combined_seg_data),
  c("chr","start","end","arm","arm_start","arm_end","arm_len","seg_len_in_arm")
)

long_data <- combined_seg_data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Value") %>%
  mutate(
    Value   = toupper(as.character(Value)),
    seg_len = seg_len_in_arm
  )

# ---- arms of interest
myeloma_arms_df <- tibble::tribble(
  ~Arm,     ~chr,  ~arm,
  "del1p",  "1",   "1p",
  "amp1q",  "1",   "1q",
  "del13q", "13",  "13q",
  "del17p", "17",  "17p"
)

gain_labels <- c("GAIN", "AMP", "HLAMP")
loss_labels <- c("HOMD", "HETD", "LOSS")

samples <- unique(long_data$Sample)
results <- tibble(Sample = samples)

# ---- length-weighted % altered per arm
for (i in seq_len(nrow(myeloma_arms_df))) {
  arm_name  <- myeloma_arms_df$Arm[i]
  chr_value <- myeloma_arms_df$chr[i]
  arm_value <- myeloma_arms_df$arm[i]
  
  denom <- arm_lengths %>%
    filter(chr == chr_value, arm == arm_value) %>%
    pull(arm_len) %>% unique()
  if (length(denom) != 1) stop("Arm length lookup failed for ", arm_value)
  
  data_arm <- long_data %>%
    filter(chr == chr_value, arm == arm_value, !is.na(Value)) %>%  # key!
    mutate(IsAltered = case_when(
      startsWith(arm_name, "amp") ~ as.integer(Value %in% gain_labels),
      startsWith(arm_name, "del") ~ as.integer(Value %in% loss_labels),
      TRUE ~ 0L
    ))
  
  sample_props <- data_arm %>%
    group_by(Sample) %>%
    summarise(PropAltered = sum(IsAltered * seg_len, na.rm = TRUE) / denom,
              .groups = "drop")
  
  sample_props <- full_join(tibble(Sample = samples), sample_props, by = "Sample") %>%
    mutate(PropAltered = replace_na(PropAltered, 0))
  
  results <- results %>%
    left_join(sample_props, by = "Sample") %>%
    mutate(!!arm_name := as.integer(PropAltered > (1/3))) %>%
    select(-PropAltered)
}

# ---- Hyperdiploidy: length-weighted full-gain on 3,5,7,9,11,15,19,21
hyperdiploid_chrs <- c("3","5","7","9","11","15","19","21")

# Adjustable for caller 
gain_thresh_chr <- 0.65   # was 0.80. 0.60–0.70 is typical for "whole chr gain"
k_trisomies     <- 5      # call hyperdiploid if at least k of 8 are gained

# compute PropGained per chr
gain_prop_tbl <- purrr::map_dfr(hyperdiploid_chrs, function(chr_value){
  # fixed denominator (whole chromosome length) to avoid NA-driven inflation
  denom_chr <- chr_lengths %>%
    filter(chr == chr_value) %>% pull(chr_len)
  if (length(denom_chr) != 1) stop("Chromosome length lookup failed for chr", chr_value)
  
  long_data %>%
    dplyr::filter(chr == chr_value, !is.na(Value)) %>%
    mutate(IsGain = as.integer(Value %in% gain_labels)) %>%
    group_by(Sample) %>%
    summarise(PropGained = sum(IsGain * seg_len, na.rm = TRUE) / denom_chr,
              .groups = "drop") %>%
    mutate(chr = chr_value)
}) %>%
  tidyr::pivot_wider(names_from = chr, values_from = PropGained,
                     names_prefix = "chr") %>%
  # ensure all samples are present
  right_join(tibble(Sample = unique(long_data$Sample)), by = "Sample") %>%
  mutate(across(starts_with("chr"), ~ tidyr::replace_na(.x, 0)))

# binary per-chr gains with relaxed threshold
hd_bin <- gain_prop_tbl %>%
  mutate(across(starts_with("chr"),
                ~ as.integer(.x > gain_thresh_chr),
                .names = "hyperdiploid_chr{substr(.col,4,6)}"))

# count trisomies and call hyperdiploid
hd_count <- hd_bin %>%
  mutate(hd_count = rowSums(select(., starts_with("hyperdiploid_chr")), na.rm = TRUE),
         hyperdiploid = as.integer(hd_count >= k_trisomies)) %>%
  select(Sample, starts_with("hyperdiploid_chr"), hd_count, hyperdiploid)

# merge into results (drop older strict cols if present)
results <- results %>%
  select(-dplyr::any_of(c("hyperdiploid")),
         -dplyr::starts_with("hyperdiploid_chr")) %>%
  left_join(hd_count, by = "Sample")

## Old way - undercalled with Sequenza
# for (chr_value in hyperdiploid_chrs) {
#   colname <- paste0("hyperdiploid_chr", chr_value)
#   
#   denom_chr <- chr_lengths %>%
#     filter(chr == chr_value) %>% pull(chr_len)
#   if (length(denom_chr) != 1) stop("Chromosome length lookup failed for chr", chr_value)
#   
#   data_chr <- long_data %>%
#     filter(chr == chr_value, !is.na(Value)) %>%          # key!
#     mutate(IsGain = as.integer(Value %in% gain_labels))
#   
#   sample_gain <- data_chr %>%
#     group_by(Sample) %>%
#     summarise(PropGained = sum(IsGain * seg_len, na.rm = TRUE) / denom_chr,
#               .groups = "drop")
#   
#   sample_gain <- full_join(tibble(Sample = samples), sample_gain, by = "Sample") %>%
#     mutate(PropGained = replace_na(PropGained, 0))
#   
#   results <- results %>%
#     left_join(sample_gain, by = "Sample") %>%
#     mutate(!!colname := as.integer(PropGained > 0.80)) %>%
#     select(-PropGained)
# # }
# 
# results <- results %>%
#   mutate(hyperdiploid = as.integer(rowSums(select(., starts_with("hyperdiploid_chr"))) == 8))

myeloma_CNA_matrix_with_HRD <- results

# ---- export final tables (same as before)
cna_data <- myeloma_CNA_matrix_with_HRD %>%
  select(Sample, del1p, amp1q, del13q, del17p, hyperdiploid) %>%
  mutate(across(-Sample, as.character))

export_dir <- "Jan2025_exported_data"
if (!dir.exists(export_dir)) dir.create(export_dir)


## Check if what have is expected
# Convert all alteration columns to numeric
cna_data_summary <- cna_data %>%
  mutate(across(-Sample, as.numeric)) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    prop_del1p       = mean(del1p == 1, na.rm = TRUE),
    prop_amp1q       = mean(amp1q == 1, na.rm = TRUE),
    prop_del13q      = mean(del13q == 1, na.rm = TRUE),
    prop_del17p      = mean(del17p == 1, na.rm = TRUE),
    prop_hyperdiploid= mean(hyperdiploid == 1, na.rm = TRUE)
  )

print(cna_data_summary)


### Edit the same so in same format as other tools expect 
# Add a Tumor_Sample_Barcode column to metada_df_mutation_comparison
# Clean up Sample names in CNA data to match metadata format
cna_data_cleaned <- cna_data %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

cna_data_cleaned <- cna_data_cleaned %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join to metadata by Sample ↔ Tumor_Sample_Barcode
cna_data_merged <- cna_data_cleaned %>%
  left_join(
    metada_df_mutation_comparison %>%
      select(Tumor_Sample_Barcode, Bam_clean_tmp),
    by = c("Sample" = "Tumor_Sample_Barcode")
  )

# Check results
message("✅ Merged CNA data with metadata: added Bam_clean_tmp column.")


## ---- Export Results ----
## The number in filenames corresponds to the gamma parameter used by Sequenza (γ=400)

# Primary CNA matrix
saveRDS(cna_data_merged, file = file.path(export_dir, "cna_data_from_sequenza_400.rds"))
write.table(cna_data_merged,
            file = file.path(export_dir, "cna_data_from_sequenza_400.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Summary statistics (proportion of samples per alteration)
saveRDS(cna_data_summary, file = file.path(export_dir, "cna_data_summary_400.rds"))
write.table(cna_data_summary,
            file = file.path(export_dir, "cna_data_summary_400.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Console confirmation
message("✅ Export complete:")