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
library(readxl)

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

export_dir <- "Jan2025_exported_data"
if (!dir.exists(export_dir)) dir.create(export_dir)

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


# ---- read Sequenza *_segments.txt (optionally .gz), rename cols, make calls
seg_files <- list.files(seg_dir, pattern = "_segments\\.txt(\\.gz)?$", full.names = TRUE)

seg_data_list <- vector("list", length(seg_files))
names(seg_data_list) <- basename(seg_files)

ploidy_list <- vector("list", length(seg_files))   


### Get the ploidy from confints
### Lastly, get the purity and ploidy estimates
confints_dir <- "Oct 2024 data/Sequenza/All_confints_400/"

# Discover files that end with _confints_CP.txt
confints <- list.files(
  confints_dir,
  pattern = "_confints_CP\\.txt$",
  full.names = TRUE
)

if (length(confints) == 0L) {
  stop("No *_confints_CP.txt files found in: ", confints_dir)
}

# Safe reader: grab the 2nd row if present; return NA otherwise
read_pp_safe <- function(path) {
  df <- suppressMessages(read_tsv(path, show_col_types = FALSE, progress = FALSE))
  if (nrow(df) < 2L) {
    warning("File has <2 rows, skipping values: ", basename(path))
    return(tibble(
      Sample = NA_character_,
      Purity = NA_real_,
      Ploidy = NA_real_,
      File   = path
    ))
  }
  
  # Some sequenza outputs use 'cellularity' and 'ploidy.estimate'
  # Normalize possible name variants just in case.
  nm <- names(df)
  nm <- nm |>
    str_replace("^cellularity$", "cellularity") |>
    str_replace("^ploidy\\.estimate$", "ploidy.estimate")
  names(df) <- nm
  
  row2 <- df |> slice(2)
  
  # Extract clean sample ID from filename
  sample_id <- basename(path) |> str_remove("_confints_CP\\.txt$")
  
  tibble(
    Sample = sample_id,
    Purity = suppressWarnings(as.numeric(row2$cellularity)),
    Ploidy = suppressWarnings(as.numeric(row2$`ploidy.estimate`)),
    File   = path
  )
}

pp_table <- map_dfr(confints, read_pp_safe) |>
  # Keep only labeled columns for downstream use
  select(Sample, Purity, Ploidy)

# Basic sanity checks + labeling polish
pp_table <- pp_table |>
  mutate(
    Purity = round(Purity, 4),
    Ploidy = round(Ploidy, 3)
  ) |>
  arrange(Sample)

write.table(pp_table,
            paste0(export_dir,"/sequenza_purity_ploidy_estimates.txt"),sep="\t", 
            row.names = FALSE, quote = FALSE)



ploidy_map <- setNames(pp_table$Ploidy, pp_table$Sample)

ploidy_list   <- vector("list", length(seg_files))
seg_data_list <- vector("list", length(seg_files))

for (i in seq_along(seg_files)) {
  f  <- seg_files[i]
  df <- readr::read_tsv(f, show_col_types = FALSE)
  
  # required columns in Sequenza segments
  req <- c("chromosome", "start.pos", "end.pos", "CNt", "A", "B")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", basename(f), ": ", paste(miss, collapse = ", "))
  
  # sample id must match the id used in pp_table$Sample
  sample_name <- sub("_segments\\.txt(\\.gz)?$", "", basename(f))
  
  # ---- core numeric frame
  seg_core <- df %>%
    dplyr::transmute(
      chr   = stringr::str_remove(chromosome, "^chr"),
      start = suppressWarnings(as.numeric(start.pos)),
      end   = suppressWarnings(as.numeric(end.pos)),
      CNt   = suppressWarnings(as.numeric(CNt)),
      A     = suppressWarnings(as.numeric(A)),
      B     = suppressWarnings(as.numeric(B))
    ) %>%
    dplyr::mutate(seg_len = pmax(end - start + 1, 1))
  
  # ---- try to use confints ploidy
  ploidy_conf <- unname(ploidy_map[sample_name])
  
  # If confints missing, fall back to autosome-weighted mode of rounded CNt
  if (!is.finite(ploidy_conf)) {
    df_auto <- seg_core %>%
      dplyr::filter(chr %in% as.character(1:22)) %>%
      dplyr::mutate(CNt_round = round(CNt))
    est_ploidy <- wmode(df_auto$CNt_round, df_auto$seg_len)  # assumes wmode() defined upstream
    baseline_int <- ifelse(is.na(est_ploidy), NA_real_, pmax(1, round(est_ploidy)))
    source_tag   <- "segments_fallback"
  } else {
    # Use Sequenza confints ploidy; snap to nearest integer for labeling
    baseline_int <- pmax(1, round(ploidy_conf))
    source_tag   <- "confints"
  }
  
  # record ploidy sources for QC
  ploidy_list[[i]] <- tibble::tibble(
    Sample        = sample_name,
    ploidy_conf   = ifelse(is.finite(ploidy_conf), ploidy_conf, NA_real_),
    baseline_int  = baseline_int,
    baseline_src  = source_tag
  )
  
  # ---- build per-segment categorical calls (baseline-aware; fallback to diploid if NA)
  if (is.na(baseline_int)) {
    seg_calls <- seg_core %>%
      dplyr::mutate(Call = call_from_CNt_AB(CNt, A, B))  # your original backup rule
  } else {
    seg_calls <- seg_core %>%
      dplyr::mutate(
        Call = dplyr::case_when(
          !is.finite(CNt)                    ~ NA_character_,
          CNt == 0                           ~ "HOMD",
          CNt == 1                           ~ "HETD",
          CNt >= baseline_int + 3            ~ "HLAMP",
          CNt >  baseline_int + 1            ~ "AMP",
          CNt >  baseline_int                ~ "GAIN",
          CNt <  baseline_int                ~ "LOSS",
          CNt == baseline_int & !is.na(B) & B == 0 ~ "CNLOH",
          CNt == baseline_int                ~ "NEUT",
          TRUE                               ~ NA_character_
        )
      )
  }
  
  seg_df <- seg_calls %>% dplyr::select(chr, start, end, Call)
  colnames(seg_df)[4] <- sample_name
  seg_data_list[[i]] <- seg_df
}

combined_seg_data <- purrr::reduce(seg_data_list, dplyr::full_join, by = c("chr","start","end")) %>%
  dplyr::arrange(factor(chr, levels = c(as.character(1:22), "X", "Y")), start)

sample_ploidy <- dplyr::bind_rows(ploidy_list) %>%
  dplyr::arrange(Sample)


print(sample_ploidy, n = 50)
write.csv(sample_ploidy, "Sep_2025_sequenza_ploidy_estimates_updated.csv", row.names = FALSE)


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



### Get FISH probe overlap
# Tunables
pad_bp  <- 150000L    # ±150 kb padding around each vendor span
min_bp  <- 1000L     # require at least 1 kb overlap 

# 1) Load probe table (if not already loaded)
if (!exists("fish_probe_locations")) {
  fish_probe_locations <- read_excel("Clinical data/FISH probe locations.xlsx")
}
stopifnot(all(c("Chromosome","Location in hg38") %in% names(fish_probe_locations)))

# 2) Map Chromosome -> feature and parse coordinates
feature_map <- c("1p" = "del1p", "1q" = "amp1q", "17p" = "del17p", "13q" = "del13q")

probe_windows <- fish_probe_locations %>%
  transmute(
    feature = dplyr::recode(Chromosome, !!!feature_map, .default = NA_character_),
    loc     = gsub(",", "", `Location in hg38`)
  ) %>%
  tidyr::extract(loc, into = c("chr","start","end"),
                 regex = "^chr([^:]+):(\\d+)-(\\d+)$", remove = TRUE) %>%
  mutate(
    chr   = gsub("^chr","", chr),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  filter(!is.na(feature), !is.na(chr), !is.na(start), !is.na(end)) %>%
  mutate(
    start_pad = pmax(1L, start - pad_bp),
    end_pad   = end + pad_bp
  ) %>%
  distinct(feature, chr, start_pad, end_pad, .keep_all = TRUE)

# 3) Build GRanges: segments and padded probe windows
sample_cols_seg <- setdiff(
  names(combined_seg_data),
  c("chr","start","end","arm","arm_start","arm_end","arm_len","seg_len_in_arm")
)

seg_gr <- GRanges(
  seqnames = as.character(combined_seg_data$chr),
  ranges   = IRanges(as.integer(combined_seg_data$start),
                     as.integer(combined_seg_data$end))
)

# Ensure calls are uppercase characters
calls_df <- combined_seg_data[, sample_cols_seg] |>
  dplyr::mutate(across(everything(), ~ toupper(as.character(.))))

# If any duplicate sample names, make them unique (prevents subsetting errors)
if (anyDuplicated(names(calls_df))) {
  names(calls_df) <- make.unique(names(calls_df), sep = "_dup")
}

# Correct: one column per sample, no name mangling
mcols(seg_gr) <- S4Vectors::DataFrame(as.list(calls_df), check.names = FALSE)

# Use the actual names from mcols downstream
samples <- colnames(mcols(seg_gr))

probe_gr <- GRanges(
  seqnames = probe_windows$chr,
  ranges   = IRanges(probe_windows$start_pad, probe_windows$end_pad),
  feature  = probe_windows$feature
)

# 4) Overlap probes with segments; compute bp overlap
hits <- findOverlaps(probe_gr, seg_gr, ignore.strand = TRUE)
if (length(hits) == 0L) {
  warning("No overlaps between padded probe windows and segments. Check genome build / chr labels.")
}

ovl_bp <- width(pintersect(ranges(probe_gr)[queryHits(hits)],
                           ranges(seg_gr)[subjectHits(hits)]))

hits_df <- tibble::tibble(
  probe_idx = as.integer(queryHits(hits)),
  seg_idx   = as.integer(subjectHits(hits)),
  ovl       = as.integer(ovl_bp)
) %>%
  dplyr::filter(ovl >= min_bp) %>%                # apply minimum overlap filter
  dplyr::arrange(probe_idx, dplyr::desc(ovl))


# 5) For each probe and each sample, take the first non-NA call among overlapping
#    segments ordered by decreasing overlap.
probe_ids <- seq_along(probe_gr)

final_calls_mat <- matrix(NA_character_, nrow = length(probe_ids), ncol = length(samples),
                          dimnames = list(NULL, samples))

for (p in probe_ids) {
  segs <- hits_df %>% dplyr::filter(probe_idx == p) %>% dplyr::pull(seg_idx)
  if (length(segs) == 0L) next
  # Subset the mcols matrix for these segments
  seg_calls <- as.matrix(mcols(seg_gr)[segs, samples, drop = FALSE])
  # For each sample, pick first non-NA down the rows
  for (j in seq_along(samples)) {
    col_vals <- seg_calls[, j]
    nn <- which(!is.na(col_vals) & nzchar(col_vals))
    if (length(nn)) final_calls_mat[p, j] <- col_vals[nn[1]]
  }
}

# 6) Tidy to long and bin by expected direction
probe_meta <- tibble::tibble(
  probe_idx = probe_ids,
  feature   = as.character(mcols(probe_gr)$feature)
)

probe_calls_long <- as_tibble(final_calls_mat, .name_repair = "minimal") %>%
  mutate(probe_idx = probe_ids) %>%
  tidyr::pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "probe_call") %>%
  left_join(probe_meta, by = "probe_idx") %>%
  select(Sample, feature, probe_call)

gain_labels <- c("GAIN","AMP","HLAMP")
loss_labels <- c("HOMD","HETD","LOSS","CNLOH")
dir_map     <- c(amp1q = "gain", del1p = "loss", del13q = "loss", del17p = "loss")

probe_calls_bin <- probe_calls_long %>%
  mutate(direction = unname(dir_map[feature]),
         is_altered_at_probe = dplyr::if_else(
           direction == "gain",
           as.integer(probe_call %in% gain_labels),
           as.integer(probe_call %in% loss_labels)
         )) %>%
  select(-direction) %>%
  tidyr::pivot_wider(names_from = feature,
                     values_from = c(probe_call, is_altered_at_probe),
                     names_sep = "_")
# =======================================================================

### Now redo, but using the cytoband instead 
# ---- Cytoband-based concordance ----

# Expected cytoband file format (UCSC): chrom  start  end  band  gieStain
# Example row:                          chr1   0      2300000  p36.33  gneg
cytoband_file <- "cytoband.txt"

cyto <- readr::read_tsv(
  cytoband_file,
  col_names = c("chrom", "start", "end", "band", "gieStain"),
  col_types = readr::cols(
    chrom = readr::col_character(),
    start = readr::col_double(),
    end = readr::col_double(),
    band = readr::col_character(),
    gieStain = readr::col_character()
  ),
  progress = FALSE
)

# Clean up chromosome column (remove "chr" prefix if present)
cyto <- cyto %>%
  mutate(
    chrom = gsub("^chr", "", chrom),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  arrange(chrom, start)


# Be forgiving about column names
nm <- names(cyto)
stopifnot(
  any(grepl("^chrom", nm, ignore.case = TRUE)),
  any(grepl("^start$", nm, ignore.case = TRUE)),
  any(grepl("^end$",   nm, ignore.case = TRUE)),
  any(grepl("^band",   nm, ignore.case = TRUE))
)

# Normalize column names we need
names(cyto)[grepl("^chrom", names(cyto), ignore.case = TRUE)] <- "chrom"
names(cyto)[grepl("^start$", names(cyto), ignore.case = TRUE)] <- "start"
names(cyto)[grepl("^end$",   names(cyto), ignore.case = TRUE)] <- "end"
names(cyto)[grepl("^band",   names(cyto), ignore.case = TRUE)]  <- "band"

cyto <- cyto %>%
  mutate(
    chrom = gsub("^chr", "", chrom),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  arrange(chrom, start)

# Helper: expand a band tag like "q21" or "q21.3" to the span covering all matching rows
.band_span <- function(cyto_chr, arm, band_tag) {
  tag <- paste0(arm, band_tag)            # e.g., "q21"
  rows <- which(startsWith(cyto_chr$band, tag))
  if (!length(rows)) return(NULL)
  tibble::tibble(
    start = min(cyto_chr$start[rows]),
    end   = max(cyto_chr$end[rows])
  )
}

# Parse Target like:
#   "17p13.1"            -> chr=17, arm1=p, band1=13.1 (single)
#   "1q21-q22"           -> chr=1,  arm1=q, band1=21 ; arm2=q, band2=22
#   "1q21.1-q22.2"       -> range with decimals
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

# Map each probe row to a cytoband-derived window
probe_windows_cytoband <- fish_probe_locations %>%
  transmute(
    feature = dplyr::recode(Chromosome, !!!feature_map, .default = NA_character_),
    Target  = Target
  ) %>%
  rowwise() %>%
  mutate(
    .pt = list(.parse_target(Target))
  ) %>%
  ungroup() %>%
  filter(!purrr::map_lgl(.pt, is.null)) %>%
  mutate(
    chr  = purrr::map_chr(.pt, ~ .x$chr),
    arm1 = purrr::map_chr(.pt, ~ .x$arm1),
    b1   = purrr::map_chr(.pt, ~ .x$band1),
    arm2 = purrr::map_chr(.pt, ~ .x$arm2),
    b2   = purrr::map_chr(.pt, ~ .x$band2)
  ) %>%
  select(-.pt) %>%
  group_by(feature, chr, arm1, b1, arm2, b2) %>%
  # Resolve to hg38 span using cytobands
  reframe({
    cy_chr <- dplyr::filter(cyto, chrom == chr) %>% arrange(start)
    sspan  <- .band_span(cy_chr, arm1, b1)
    espan  <- .band_span(cy_chr, arm2, b2)
    if (is.null(sspan) || is.null(espan)) {
      tibble::tibble(start = NA_integer_, end = NA_integer_)
    } else {
      tibble::tibble(
        start = min(sspan$start, espan$start),
        end   = max(sspan$end,   espan$end)
      )
    }
  }) %>%
  ungroup() %>%
  filter(!is.na(start), !is.na(end)) %>%
  mutate(
    start_pad = pmax(1L, start - pad_bp),
    end_pad   = end + pad_bp
  ) %>%
  distinct(feature, chr, start_pad, end_pad, .keep_all = TRUE)

# Quick sanity check of the tricky one (1q21-q22): should span all q21.* through q22.* on chr1
# print(dplyr::filter(probe_windows_cytoband, grepl("^amp1q$", feature)))

# Build GRanges for cytoband windows
probe_gr_cyto <- GRanges(
  seqnames = probe_windows_cytoband$chr,
  ranges   = IRanges(probe_windows_cytoband$start_pad, probe_windows_cytoband$end_pad),
  feature  = probe_windows_cytoband$feature
)

# Overlaps and bp widths
hits_cyto <- findOverlaps(probe_gr_cyto, seg_gr, ignore.strand = TRUE)
if (length(hits_cyto) == 0L) {
  warning("No overlaps between cytoband windows and segments. Check build/labels.")
}

ovl_bp_cyto <- width(pintersect(
  ranges(probe_gr_cyto)[queryHits(hits_cyto)],
  ranges(seg_gr)[subjectHits(hits_cyto)]
))

hits_df_cyto <- tibble::tibble(
  probe_idx = as.integer(queryHits(hits_cyto)),
  seg_idx   = as.integer(subjectHits(hits_cyto)),
  ovl       = as.integer(ovl_bp_cyto)
) %>%
  dplyr::filter(ovl >= min_bp) %>%
  dplyr::arrange(probe_idx, dplyr::desc(ovl))


# Senestivity/overlap chooser
probe_ids_cyto <- seq_along(probe_gr_cyto)
final_calls_mat_cyto <- matrix(
  NA_character_, nrow = length(probe_ids_cyto), ncol = length(samples),
  dimnames = list(NULL, samples)
)

# Severity ladders (edit if your labels differ)
gain_priority <- c("HLAMP","AMP","GAIN")
loss_priority <- c("HOMD","HETD","LOSS","CNLOH")

# Feature per cytoband probe
feature_vec_cyto <- as.character(mcols(probe_gr_cyto)$feature)

for (p in probe_ids_cyto) {
  seg_rows <- hits_df_cyto %>% dplyr::filter(probe_idx == p)
  segs     <- seg_rows$seg_idx
  if (!length(segs)) next
  
  # Segment calls for just these overlaps (rows = segments, cols = samples)
  seg_calls <- as.matrix(mcols(seg_gr)[segs, samples, drop = FALSE])
  
  # Overlap widths aligned to 'segs'
  ovl_here <- seg_rows$ovl
  
  # Direction-specific priority for this probe
  feat <- feature_vec_cyto[p]
  dir  <- unname(dir_map[feat])   # "gain" or "loss" from your existing dir_map
  pri  <- if (identical(dir, "gain")) gain_priority else loss_priority
  
  # Precompute severity rank per cell (lower is more severe)
  sev_rank <- apply(seg_calls, 2, function(col) {
    match(toupper(col), pri, nomatch = length(pri) + 1L)
  })
  
  # Pick, per sample: most severe; tie-break by larger bp overlap
  for (j in seq_along(samples)) {
    calls_j <- seg_calls[, j]
    valid   <- which(!is.na(calls_j) & nzchar(calls_j))
    if (!length(valid)) next
    
    rnk  <- sev_rank[valid, j]
    best <- which(rnk == min(rnk))
    if (length(best) > 1L) {
      best <- best[which.max(ovl_here[valid][best])]
    }
    final_calls_mat_cyto[p, j] <- calls_j[valid][best]
  }
}

# Tidy long and bin calls by expected direction, mirroring your probe_calls_bin
probe_meta_cyto <- tibble::tibble(
  probe_idx = probe_ids_cyto,
  feature   = as.character(mcols(probe_gr_cyto)$feature)
)

probe_calls_long_cyto <- as_tibble(final_calls_mat_cyto, .name_repair = "minimal") %>%
  mutate(probe_idx = probe_ids_cyto) %>%
  tidyr::pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "probe_call") %>%
  dplyr::left_join(probe_meta_cyto, by = "probe_idx") %>%
  dplyr::select(Sample, feature, probe_call)

probe_calls_bin_cytoband <- probe_calls_long_cyto %>%
  mutate(
    direction = unname(dir_map[feature]),
    is_altered_at_probe = dplyr::if_else(
      direction == "gain",
      as.integer(probe_call %in% gain_labels),
      as.integer(probe_call %in% loss_labels)
    )
  ) %>%
  select(-direction) %>%
  tidyr::pivot_wider(
    names_from  = feature,
    values_from = c(probe_call, is_altered_at_probe),
    names_sep   = "_"
  )





# (Optional) quick sanity checks
# combined_seg_data %>% summarize(prop_assigned = mean(!is.na(arm)))
# combined_seg_data %>% filter(is.na(arm)) %>% count(chr, sort = TRUE)

# ---- save combined (optional)
write.csv(combined_seg_data, "Sep_2025_combined_sequenza_calls_400_updated.csv", row.names = FALSE)
saveRDS(combined_seg_data,  "Sep_2025_combined_sequenza_calls_400_updated.rds")

write.csv(probe_calls_bin, "Sep_2025_FISH_probe_calls_400_updated.csv", row.names = FALSE)
saveRDS(probe_calls_bin,  "Sep_2025_FISH_probe_calls_400_updated.rds")

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
loss_labels <- c("HOMD", "HETD", "LOSS", "CNLOH")

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
k_trisomies     <- 5      # call hyperdiploid if at least k of 5 are gained

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
saveRDS(cna_data_merged, file = file.path(export_dir, "cna_data_from_sequenza_400_updated.rds"))
write.table(cna_data_merged,
            file = file.path(export_dir, "cna_data_from_sequenza_400_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Summary statistics (proportion of samples per alteration)
saveRDS(cna_data_summary, file = file.path(export_dir, "cna_data_summary_400_updated.rds"))
write.table(cna_data_summary,
            file = file.path(export_dir, "cna_data_summary_400_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Console confirmation
message("✅ Export complete:")


### Now do the same for FISH calls 
FISH_data_cleaned <- probe_calls_bin %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

FISH_data_cleaned <- FISH_data_cleaned %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join to metadata by Sample ↔ Tumor_Sample_Barcode
FISH_data_cleaned <- FISH_data_cleaned %>%
  left_join(
    metada_df_mutation_comparison %>%
      select(Tumor_Sample_Barcode, Bam_clean_tmp, Patient),
    by = c("Sample" = "Tumor_Sample_Barcode")
  )


# Export
saveRDS(FISH_data_cleaned, file = file.path(export_dir, "FISH_data_from_sequenza_400_updated.rds"))
write.table(FISH_data_cleaned,
            file = file.path(export_dir, "FISH_data_from_sequenza_400_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


## Try again but by cytoband 
### Now do the same for FISH calls 
FISH_data_cleaned <- probe_calls_bin_cytoband %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

FISH_data_cleaned <- FISH_data_cleaned %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join to metadata by Sample ↔ Tumor_Sample_Barcode
FISH_data_cleaned <- FISH_data_cleaned %>%
  left_join(
    metada_df_mutation_comparison %>%
      select(Tumor_Sample_Barcode, Bam_clean_tmp, Patient),
    by = c("Sample" = "Tumor_Sample_Barcode")
  )


# Export
saveRDS(FISH_data_cleaned, file = file.path(export_dir, "FISH_data_from_sequenza_400_by_cytoband_updated.rds"))
write.table(FISH_data_cleaned,
            file = file.path(export_dir, "FISH_data_from_sequenza_400_by_cytoband_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


## Lastly for ploidy 

### Now do the same for FISH calls 
sample_ploidy <- sample_ploidy %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

sample_ploidy <- sample_ploidy %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join to metadata by Sample ↔ Tumor_Sample_Barcode
sample_ploidy <- sample_ploidy %>%
  left_join(
    metada_df_mutation_comparison %>%
      select(Tumor_Sample_Barcode, Bam_clean_tmp, Patient, Timepoint),
    by = c("Sample" = "Tumor_Sample_Barcode")
  )


## Export 
saveRDS(sample_ploidy, file = file.path(export_dir, "Sample_ploidy_from_sequenza_400.rds"))
write.table(sample_ploidy,
            file = file.path(export_dir, "Sample_ploidy_from_sequenza.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)




## Quick check
fish_summary2 <- filled_df %>% ## this is generated in later script 
  dplyr::select(Patient, DEL_17P, DEL_1P, AMP_1Q) %>%
  dplyr::left_join(cohort_df, by = "Patient") %>%
  dplyr::filter(!is.na(Cohort)) %>%
  dplyr::filter(
    DEL_17P == "Positive" |
      DEL_1P  == "Positive" |
      AMP_1Q  == "Positive"
  )

fish_summary2 <- fish_summary %>% 
  left_join(FISH_data_cleaned)
