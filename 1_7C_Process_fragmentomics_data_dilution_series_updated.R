# ──────────────────────────────────────────────────────────────────────────────
# process_fragmentomics_dilution_series.R
#
# Purpose  ────────────────────────────────────────────────────────────────────
#   1. Load *dilution-series* nucleosome-accessibility distance files (cfWGS) 
#      and matching PON (healthy) files.
#   2. Harmonize sample IDs, join clinical keys, compute per-site:
#        • Mean.Coverage
#        • Midpoint.Coverage  (at Positions -30, -15, 0, +15, +30)
#        • Midpoint.normalized
#        • Amplitude   (max of the periodogram)
#   3. For each site, calculate:
#        • Z-scores vs. healthy (PON) for Coverage, Midpoint, Amplitude
#        • Fold-changes (tumour vs. healthy) for the three metrics
#        • ±10% FDR thresholds based on healthy z-scores
#        • Binary “Threshold” flags per sample for each metric
#   4. Extract the “MM_DARs_chromatin_activation” site rows, re-attach  
#      clinical metadata (Bam, Patient, Date_of_sample_collection).
#   5. Load insert_size_summary.tsv and fragment_scores.tsv (dilution cohort), 
#      harmonize sample IDs, and merge with the MM_DARs data.
#   6. Optionally attach any “dilution_series_metadata.csv” if present.
#   7. Write out a single CSV:
#         Results_Fragmentomics/Dilution_series/key_fragmentomics_info_dilution_series.csv
#
# Author:  Dory Abelman
# Updated: 2025-06-06
# ──────────────────────────────────────────────────────────────────────────────

### 0.  PACKAGES & HELPER FUNCTIONS  ############################################
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(data.table)
  library(GeneCycle)   # for periodogram()
})

# (If you have any custom functions in session.functions.R, load them here:)
if (file.exists("session.functions.R")) {
  source("session.functions.R")
}

# Compute a z-score of x versus the vector y
calculate.zscore <- function(x, y) {
  (x - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
}

# Helper to parse Fold-Change + p-value from t.test (group2 vs group1)
get.foldchange.p <- function(i, group1.idx, group2.idx) {
  test <- t.test(i[group2.idx], i[group1.idx], na.rm = TRUE)
  fc   <- mean(i[group2.idx], na.rm = TRUE) / mean(i[group1.idx], na.rm = TRUE)
  return(c(fc = fc, p = test$p.value))
}

# Clean sample names exactly as in your main scripts:
#  • Replace “_Blood_plasma_cfDNA” → “-P”
#  • Swap “_” → “-” for non-SPORE
#  • Swap “IMG” → “MyP”
clean_sample <- function(x) {
  x %>%
    str_replace("_Blood_plasma_cfDNA$", "-P") %>%
    ifelse(!str_starts(x, "SPORE"),
           str_replace_all(., "_", "-"),
           .) %>%
    str_replace("IMG", "MyP")
}


### 1.  DEFINE PATHS  ###########################################################
# Point these to your *dilution-series* folders:

# 1a) Nucleosome-accessibility input (cfWGS for dilution series)
nuc_input.dir <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series"

# 1b) Nucleosome-accessibility PON folder (healthy controls)
nuc_pon.dir   <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Normals"

# 1c) Insert-size + FS input (dilution series)
ins_fs.dir    <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series"

# 1d) Combined clinical CSV (used for joining Sample → Bam, Patient, Date…)
clinical.csv  <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Dilution_series/Metadata_dilution_series.csv"

# 1e) If you have additional dilution series metadata file (e.g. dilution percentages), put it here:
# meta.csv      <- file.path(ins_fs.dir, "Metadata_dilution_series.csv") ## Not used

# 1f) Output folder for dilution series results
out.dir       <- "~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Results_Fragmentomics/Dilution_series"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)


### 2.  LOAD CLINICAL TABLE   ##################################
combined_clinical <- read_csv(clinical.csv, show_col_types = FALSE)


### 3.  READ & PROCESS NUCLEOSOME-ACCESSIBILITY (cfWGS + PON)  ##################

# 3a) List all cfWGS dilution files
results.files <- list.files(path = nuc_input.dir,
                            pattern = "nucleosome_accessibility_distances.tsv$",
                            full.names = TRUE)

# 3b) List PON (healthy) files; we only need the first one if they’re identical format
pon.files <- rev(sort(
  list.files(path = nuc_pon.dir,
             pattern = "nucleosome_accessibility_distances.tsv$",
             full.names = TRUE)
))

# 3c) Read in all cfWGS dilution series data, bind rows
cfWGS.data <- results.files %>%
  lapply(read_tsv, show_col_types = FALSE) %>%
  bind_rows()

# 3d) If any “TGL49” columns snuck in as columns, drop them (same as main script)
if (any(grepl("TGL49", colnames(cfWGS.data)))) {
  cfWGS.data <- cfWGS.data[, !grepl("TGL49", colnames(cfWGS.data))]
}

# 3e) Read in PON once
pon.data <- read_tsv(pon.files[[1]], show_col_types = FALSE)

# 3f) Harmonize Sample names in cfWGS.data exactly as before,
#     then left_join(tmp.clin) to grab Bam/Patient/Date/timepoint_info
cfWGS.data <- cfWGS.data %>%
  mutate(
    Sample = clean_sample(Sample),
    Sample = if_else(! str_starts(Sample, "SPORE"),
                     Sample,      # keep SPORE as-is
                     Sample),
    Sample = str_replace(Sample, "IMG", "MyP")
  ) %>%
  left_join(combined_clinical, by = c("Sample" = "Merge")) %>%
  mutate(Cohort = "cfWGS")

# 3g) Mark PON samples as “HBC” (healthy bank controls)
pon.data <- pon.data %>%
  mutate(
    Sample = clean_sample(Sample),
    Cohort = "HBC"
  )

# 3h) Build sample lists
all.tumour.samples <- unique(cfWGS.data$Sample)
all.pon.samples    <- setdiff(unique(pon.data$Sample), "TGL49_0267_Pb_U_PE_428") # remove duplicate longitudinal
all.samples        <- c(all.tumour.samples, all.pon.samples)

# 3i) Grab all sites present in the cfWGS (dilution) data
all.sites <- unique(cfWGS.data$site_name)


### 4.  INITIALIZE STORAGE FOR METRICS, SCORES, STATS  #########################
metrics.per.site <- list()
scores.per.site  <- list()

### BUILD SAMPLE LISTS FOR GROUP COMPARISONS #########################################################
tumour.samples     <- all.tumour.samples
pon.samples        <- all.pon.samples

all.samples <- c(tumour.samples, pon.samples)
all.sites   <- unique(cfWGS.data$site_name)


### UTILITY FUNCTIONS ##############################################################################
calculate.zscore <- function(x, y) {
  (x - mean(y)) / sd(y)
}

get.ttest.p.and.foldchange <- function(i, group1, group2) {
  p  <- t.test(i[group2], i[group1])$p.value
  fc <- mean(i[group2]) / mean(i[group1])
  return(c(fc, p))
}


### KEEP ONLY COLUMNS THAT APPEAR IN BOTH cfWGS & PON ##############################################
common_columns <- intersect(colnames(cfWGS.data), colnames(pon.data))
cfWGS.data      <- cfWGS.data[, common_columns, drop = FALSE]
pon.data        <- pon.data[, common_columns, drop = FALSE]


### 5.  LOOP OVER EACH SITE: COMPUTE METRICS + Z-SCORES + FC + THRESHOLDS  #####
### LOOP OVER EACH SITE: COMPUTE METRICS + Z-SCORES + T-TESTS #######################################
metrics.per.site <- list()
scores.per.site  <- list()

stats.data <- data.frame(
  Site            = unique(cfWGS.data[, c("site_name", "site_type")]     )$site_name,
  Type            = unique(cfWGS.data[, c("site_name", "site_type")]     )$site_type,
  Coverage.fc     = NA,  Midpoint.fc     = NA,  Amplitude.fc     = NA,
  Coverage.p      = NA,  Midpoint.p      = NA,  Amplitude.p      = NA,
  stringsAsFactors = FALSE
)

for (site in all.sites) {
  # combine cfWGS + PON rows for this site
  site.data <- rbind(
    cfWGS.data %>% filter(site_name == site),
    pon.data   %>% filter(site_name == site)
  )
  site.data <- as.data.table(site.data) %>% unique()
  
  # collapse duplicates and convert to data.table
  site.data <- site.data %>%
    group_by(Sample, site_name, site_type, Position, Cohort) %>%
    summarise(
      Coverage = mean(Coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.table()
  
  # cast to wide format and convert to data.frame
  gc.distances <- data.table::dcast(
    site.data,
    Position ~ Sample,
    value.var = "Coverage"
  ) %>% as.data.frame()
  colnames(gc.distances) <- gsub("^Coverage\\.", "", colnames(gc.distances))
  
  # keep only valid samples in both places
  valid.samples <- intersect(colnames(gc.distances), all.samples)
  all.samples <- valid.samples
  
  # store this distance matrix
  metrics.per.site[[site]] <- data.frame(
    Sample            = valid.samples,
    Site              = rep(site, length(valid.samples)),
    Mean.Coverage     = colMeans(gc.distances[, valid.samples, drop = FALSE]),
    Midpoint.Coverage = colMeans(
      gc.distances[gc.distances$Position %in% c(-30, -15, 0, 15, 30), valid.samples, drop = FALSE]
    ),
    Midpoint.normalized = NA,
    Amplitude         = apply(GeneCycle::periodogram(gc.distances[, valid.samples, drop = FALSE])[["spec"]], 2, max)
  )
  
  metrics.per.site[[site]]$Midpoint.normalized <-
    (metrics.per.site[[site]]$Midpoint.Coverage - metrics.per.site[[site]]$Mean.Coverage) + 1
  
  # z-scores vs. PON
  score.data <- data.frame(
    Sample        = valid.samples,
    Site          = rep(site, length(valid.samples)),
    Zscore.Coverage =
      sapply(metrics.per.site[[site]]$Mean.Coverage,
             calculate.zscore,
             y = metrics.per.site[[site]]$Mean.Coverage[
               metrics.per.site[[site]]$Sample %in% pon.samples
             ]
      ),
    Zscore.Midpoint =
      sapply(metrics.per.site[[site]]$Midpoint.normalized,
             calculate.zscore,
             y = metrics.per.site[[site]]$Midpoint.normalized[
               metrics.per.site[[site]]$Sample %in% pon.samples
             ]
      ),
    Zscore.Amplitude =
      sapply(metrics.per.site[[site]]$Amplitude,
             calculate.zscore,
             y = metrics.per.site[[site]]$Amplitude[
               metrics.per.site[[site]]$Sample %in% pon.samples
             ]
      )
  )
  scores.per.site[[site]] <- score.data
  
  # group comparisons (fold-change & p-value)
  mdata <- metrics.per.site[[site]]
  
  # tumour vs. healthy
  idx_pon   <- which(mdata$Sample %in% pon.samples)
  idx_tumor <- which(mdata$Sample %in% tumour.samples)
  stats.data[stats.data$Site == site, c("Coverage.fc", "Coverage.p")] <-
    get.ttest.p.and.foldchange(mdata$Mean.Coverage, idx_pon, idx_tumor)
  stats.data[stats.data$Site == site, c("Midpoint.fc", "Midpoint.p")] <-
    get.ttest.p.and.foldchange(mdata$Midpoint.normalized, idx_pon, idx_tumor)
  stats.data[stats.data$Site == site, c("Amplitude.fc", "Amplitude.p")] <-
    get.ttest.p.and.foldchange(mdata$Amplitude, idx_pon, idx_tumor)
  }

### 6.  ADJUST P-VALUES (FDR) #########################
### ADJUST P-VALUES (FDR) ############################################################################
fdr.data <- data.frame(
  Cohort.Cov.fdr   = p.adjust(stats.data$Coverage.p,    "fdr"),
  Cohort.Mid.fdr   = p.adjust(stats.data$Midpoint.p,    "fdr"),
  Cohort.Amp.fdr   = p.adjust(stats.data$Amplitude.p,   "fdr")
)


### 7.  COMBINE METRICS + SCORES + APPLY CLASSIFICATION FLAGS  #################
results.data <- merge(
  do.call(rbind, metrics.per.site),
  do.call(rbind, scores.per.site),
  all = TRUE
)

# Derive HBC-based thresholds (±10% of max/min z-score among healthy samples)
thresholds.up <- aggregate(
  results.data[grep("TGL49", results.data$Sample), grep("Zscore", colnames(results.data))],
  by = list(Site = results.data[grep("TGL49", results.data$Sample), ]$Site),
  FUN = function(i) { max(i) * 1.10 }
)
thresholds.down <- aggregate(
  results.data[grep("TGL49", results.data$Sample), grep("Zscore", colnames(results.data))],
  by = list(Site = results.data[grep("TGL49", results.data$Sample), ]$Site),
  FUN = function(i) { min(i) * 1.10 }
)

thresholds.up$Site   <- factor(thresholds.up$Site, levels = stats.data$Site)
thresholds.down$Site <- factor(thresholds.down$Site, levels = stats.data$Site)
thresholds.up         <- thresholds.up[order(thresholds.up$Site), ]
thresholds.down       <- thresholds.down[order(thresholds.down$Site), ]

stats.data$Coverage.threshold  <- thresholds.up$Zscore.Coverage
stats.data$Midpoint.threshold  <- thresholds.up$Zscore.Midpoint
stats.data$Amplitude.threshold <- thresholds.up$Zscore.Amplitude

# For sites where fold-change<0, use “down” threshold instead:
neg.fc <- which(stats.data$Coverage.fc < 0)
if (length(neg.fc)) {
  stats.data$Coverage.threshold[neg.fc] <-
    thresholds.down$Zscore.Coverage[neg.fc]
}
neg.fc2 <- which(stats.data$Midpoint.fc < 0)
if (length(neg.fc2)) {
  stats.data$Midpoint.threshold[neg.fc2] <-
    thresholds.down$Zscore.Midpoint[neg.fc2]
}
neg.fc3 <- which(stats.data$Amplitude.fc < 0)
if (length(neg.fc3)) {
  stats.data$Amplitude.threshold[neg.fc3] <-
    thresholds.down$Zscore.Amplitude[neg.fc3]
}

# apply thresholds to classify each sample’s z-score
results.data <- results.data %>%
  group_by(Site) %>%
  mutate(
    Threshold.Coverage = if (sign(stats.data$Coverage.fc[stats.data$Site == Site]) > 0) {
      Zscore.Coverage > stats.data$Coverage.threshold[stats.data$Site == Site]
    } else {
      Zscore.Coverage < stats.data$Coverage.threshold[stats.data$Site == Site]
    },
    Threshold.Midpoint = if (sign(stats.data$Midpoint.fc[stats.data$Site == Site]) > 0) {
      Zscore.Midpoint > stats.data$Midpoint.threshold[stats.data$Site == Site]
    } else {
      Zscore.Midpoint < stats.data$Midpoint.threshold[stats.data$Site == Site]
    },
    Threshold.Amplitude = if (sign(stats.data$Amplitude.fc[stats.data$Site == Site]) > 0) {
      Zscore.Amplitude > stats.data$Amplitude.threshold[stats.data$Site == Site]
    } else {
      Zscore.Amplitude < stats.data$Amplitude.threshold[stats.data$Site == Site]
    }
  ) %>%
  ungroup()

# Save out per-site metrics + stats
dir.create(file.path(out.dir, "Griffin"), recursive = TRUE, showWarnings = FALSE)

write_tsv(results.data,
          file.path(out.dir, "Griffin/griffin_per_site_metrics_dilution.tsv"))

write_tsv(stats.data,
          file.path(out.dir, "Griffin/griffin_per_site_stats_dilution.tsv"))


### 8.  EXTRACT “MM_DARs_chromatin_activation” SITE + ATTACH CLINICAL KEYS #######
mm_dars_small <- results.data %>%
  filter(Site == "MM_DARs_chromatin_activation") %>%
  # re-attach Bam, Patient, Date_of_sample_collection from clinical.csv
  left_join(
    combined_clinical %>% 
      dplyr::select(Merge, Bam, Patient, Sample_ID, LOD),
    by = c("Sample" = "Merge")
  ) %>%
  dplyr::select(
    Sample, Site,
    Mean.Coverage, Midpoint.Coverage, Midpoint.normalized, Amplitude,
    Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude,
    Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
    Bam, Patient, Sample_ID, LOD
  ) %>%
  distinct()


### 9.  READ & MERGE INSERT-SIZE + FRAGMENT-SCORE (DILUTION) ###################
# 9a) INSERT-SIZE (“Proportion.Short”)
ins.files <- list.files(path = ins_fs.dir,
                        pattern = "insert_size_summary.tsv$",
                        full.names = TRUE)

insert_df <- ins.files %>%
  lapply(read_tsv, show_col_types = FALSE) %>%
  bind_rows() %>%
  mutate(Sample = clean_sample(Sample)) %>%
  dplyr::select(Sample, Proportion.Short)

# 9b) FRAGMENT SCORE (“FS”)
fs.files <- list.files(path = ins_fs.dir,
                       pattern = "fragment_scores.tsv$",
                       full.names = TRUE)

fs_df <- fs.files %>%
  lapply(read_tsv, show_col_types = FALSE) %>%
  bind_rows() %>%
  mutate(Sample = clean_sample(Sample)) %>%
  dplyr::select(Sample, FS)


### 10. OPTIONAL: LOAD ADDITIONAL DILUTION METADATA  ######################################
if (file.exists(meta.csv)) {
  meta_df <- read_csv(meta.csv, show_col_types = FALSE) %>%
    mutate(Sample = clean_sample(Sample))
} else {
  meta_df <- tibble(Sample = character(0))
}


### 11. MERGE EVERYTHING INTO ONE “KEY” TABLE  ##################################
key_frag_dilution <- mm_dars_small %>%
  dplyr::select(-Site) %>%
  # now join insert-size and FS
  left_join(insert_df, by = "Sample") %>%
  left_join(fs_df,     by = "Sample") %>%
  # If any samples miss DARs entries (e.g. PON-only), they'll still appear with NA for DARs columns
  arrange(Sample)


## Remove healthy control info since already have in other table 
key_frag_dilution <- key_frag_dilution %>% filter(!is.na(Bam))

### 12. WRITE OUT FINAL CSV & RDS  #############################################
write_csv(
  key_frag_dilution,
  file = file.path(out.dir, "key_fragmentomics_info_dilution_series.csv")
)

write_rds(
  key_frag_dilution,
  file = file.path(out.dir, "key_fragmentomics_info_dilution_series.rds")
)


message("✔️  Done! Wrote: ", 
        file.path(out.dir, "key_fragmentomics_info_dilution_series.csv"))
