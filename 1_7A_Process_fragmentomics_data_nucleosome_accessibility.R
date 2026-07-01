# ──────────────────────────────────────────────────────────────────────────────
# 1_7A_Process_fragmentomics_data_nucleosome_accessibility.R
#
# Purpose:
#   • Load nucleosome‐accessibility distance data for cfWGS plasma samples (MM patients)
#     and PON (healthy controls).
#   • Harmonize sample IDs against “combined_clinical_data_updated_April2025.csv”.
#   • Compute per‐site metrics (mean coverage, midpoint coverage, amplitude).
#   • Calculate z‐scores vs. healthy (PON) for each metric.
#   • Run group‐wise t‐tests:
#        – Cohort (tumour vs. healthy)
#        – Baseline vs. healthy
#        – Maintenance vs. healthy
#        – Baseline vs. Maintenance
#        – M4/IMG vs. SPORE
#   • Adjust p‐values (FDR) and derive ±10% z‐score thresholds.
#   • Save:
#        – “Results_Fragmentomics/Griffin/griffin_per_site_metrics.tsv”
#        – “Results_Fragmentomics/Griffin/griffin_per_site_stats.tsv”
#   • Build “MM_DARs_chromatin_activation_data.csv” and save.
#
# Assumptions & file paths:
#   • input.dir <- "Fragmentomics_data"
#   • pon.dir   <- "Fragmentomics_data/Normals"
#   • out.dir   <- "Results_Fragmentomics"
#   • Clinical CSVs in working dir:
#        – combined_clinical_data_updated_April2025.csv
#        – Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.csv
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_7A_Process_fragmentomics_data_nucleosome_accessibility.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script processes nucleosome-accessibility
#     features for fragmentomics model training, longitudinal summaries, and
#     dilution-series feature processing.
#
# Author: Dory Abelman
# Last updated: 2025-05-23
# ──────────────────────────────────────────────────────────────────────────────
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#


### PREPARE SESSION ################################################################################
library(BoutrosLab.plotting.general)
library(GeneCycle)
library(pracma)
library(dplyr)
library(readr)
library(stringr)
library(data.table)

session_functions_path <- c(
  "session.functions.R",
  file.path("..", "..", "session.functions.R")
)
session_functions_path <- session_functions_path[file.exists(session_functions_path)][1]
if (is.na(session_functions_path)) {
  stop(
    "Could not find session.functions.R. Run from the project root or stage the helper at the project root.",
    call. = FALSE
  )
}
source(session_functions_path)
rm(session_functions_path)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

input.dir <- "Fragmentomics_data"
pon.dir   <- file.path("Fragmentomics_data", "Normals")
out.dir   <- "Results_Fragmentomics"

if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

### Load clinical tables ############################################################################
combined_clinical_data_updated <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
)

### READ NUCLEOSOME-ACCESSIBILITY DISTANCE FILES ####################################################
# GRIFFIN outputs one "nucleosome_accessibility_distances.tsv" per sample.
# Each row is a genomic Position (bp offset from site centre, e.g. -500 to +500),
# and the value column "Coverage" is the normalised cfDNA read depth at that offset.
# Together these rows form the nucleosome-accessibility coverage profile around
# each regulatory site; a dip at position 0 with flanking peaks indicates an
# accessible (nucleosome-depleted) region, as seen in myeloma-specific DARs.
results.files <- list.files(path = input.dir,
                            pattern = "nucleosome_accessibility_distances.tsv",
                            full.names = TRUE)
results.files <- unique(c(
  results.files,
  spring2026_revision_files(
    "Fragmentomics_Pipeeline_Suite_all_outputs",
    "^2026-06-25_cfWGS_MM_fragmentomics_Revisions_Spring2026_nucleosome_accessibility_distances[.]tsv$"
  )
))
pon.files <- rev(sort(
  list.files(path = pon.dir,
             pattern = "nucleosome_accessibility_distances.tsv",
             full.names = TRUE)
))

## Read cfWGS data
cfWGS.data <- results.files %>%
  lapply(read_tsv) %>%
  bind_rows()

# If TGL49 columns sneaked in as columns, drop them
if (any(grepl("TGL49", colnames(cfWGS.data)))) {
  cfWGS.data <- cfWGS.data[, !grepl("TGL49", colnames(cfWGS.data))]
}

# Read PON once
pon.data <- read_tsv(pon.files[[1]])


### FORMAT SAMPLE NAMES & JOIN CLINICAL KEYS ########################################################
tmp <- combined_clinical_data_updated %>%
  mutate(
    Merge_ID = if_else(
      str_starts(Bam, "SPORE"),
      str_extract(Bam, "^[^\\.]+"),
      str_replace(Sample_ID, "_Blood_plasma_cfDNA", "-P")
    ),
    Cohort = dplyr::coalesce(Study, Study.x, Study.y)
  )
tmp <- bind_rows(
  tmp,
  tmp %>% mutate(Merge_ID = str_replace(Merge_ID, "^IMG", "MyP"))
) %>%
  distinct(Merge_ID, Bam, Patient, Sample_ID, .keep_all = TRUE)

sample_metadata <- tmp %>%
  dplyr::select(Sample = Merge_ID, Bam, Patient, Date_of_sample_collection) %>%
  dplyr::distinct()

cfWGS.data <- cfWGS.data %>%
  mutate(
    Sample = str_replace(Sample, "_Blood_plasma_cfDNA", "-P"),
    Sample = if_else(
      !str_starts(Sample, "SPORE"),
      str_replace_all(Sample, "_", "-"),
      Sample
    ),
    Sample = str_replace(Sample, "IMG", "MyP")
  ) %>%
  left_join(tmp, by = c("Sample" = "Merge_ID"))

pon.data <- pon.data %>%
  mutate(Cohort = "HBC")


### BUILD SAMPLE LISTS FOR GROUP COMPARISONS #########################################################
tumour.samples     <- unique(cfWGS.data$Sample)
baseline.samples   <- unique(cfWGS.data %>% filter(timepoint_info %in% c("Baseline", "Diagnosis"))   %>% pull(Sample))
maintenance.samples<- unique(cfWGS.data %>% filter(timepoint_info %in% c("Post_induction", "Post_transplant",
                                                                         "Maintenance", "1.5yr maintenance",
                                                                         "2yr maintenance"))         %>% pull(Sample))
M4_IMG.samples     <- tumour.samples[!grepl("SPORE", tumour.samples)]
SPORE.samples      <- tumour.samples[grepl("SPORE", tumour.samples)]
pon.samples        <- setdiff(unique(pon.data$Sample), "TGL49_0267_Pb_U_PE_428")

all.samples <- c(tumour.samples, pon.samples)
all.sites   <- unique(cfWGS.data$site_name)


### UTILITY FUNCTIONS ##############################################################################
# calculate.zscore: standard z-score of a single value x relative to reference
# vector y (the PON/healthy-control distribution for this site & metric).
# Formula: z = (x - mean_PON) / sd_PON
# A positive z means the sample is MORE accessible than healthy controls;
# for tumour-enriched DARs this should be elevated in MM plasma samples.
calculate.zscore <- function(x, y) {
  (x - mean(y)) / sd(y)
}

get.ttest.p.and.foldchange <- function(i, group1, group2) {
  x <- i[group2]
  y <- i[group1]
  p  <- safe_ttest_p(x, y)
  fc <- mean(x, na.rm = TRUE) / mean(y, na.rm = TRUE)
  return(c(fc, p))
}

safe_ttest_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  if (stats::sd(x) == 0 && stats::sd(y) == 0) return(NA_real_)
  stats::t.test(x, y)$p.value
}

safe_min_p <- function(...) {
  p <- c(...)
  if (all(is.na(p))) return(NA_real_)
  min(p, na.rm = TRUE)
}


### KEEP ONLY COLUMNS THAT APPEAR IN BOTH cfWGS & PON ##############################################
common_columns <- intersect(colnames(cfWGS.data), colnames(pon.data))
cfWGS.data      <- cfWGS.data[, common_columns, drop = FALSE]
pon.data        <- pon.data[, common_columns, drop = FALSE]


### LOOP OVER EACH SITE: COMPUTE METRICS + Z-SCORES + T-TESTS #######################################
metrics.per.site <- list()
scores.per.site  <- list()

stats.data <- data.frame(
  Site            = unique(cfWGS.data[, c("site_name", "site_type")]     )$site_name,
  Type            = unique(cfWGS.data[, c("site_name", "site_type")]     )$site_type,
  Coverage.fc     = NA,  Midpoint.fc     = NA,  Amplitude.fc     = NA,
  Coverage.p      = NA,  Midpoint.p      = NA,  Amplitude.p      = NA,
  PR.coverage.p   = NA,  PR.midpoint.p   = NA,  PR.amplitude.p   = NA,
  PS.coverage.p   = NA,  PS.midpoint.p   = NA,  PS.amplitude.p   = NA,
  response.p      = NA,
  trial.p         = NA,
  stringsAsFactors = FALSE
)

for (site in all.sites) {
  # combine cfWGS + PON rows for this site
  site.data <- rbind(
    cfWGS.data %>% filter(site_name == site),
    pon.data   %>% filter(site_name == site)
  )
  site.data <- as.data.table(site.data) %>% unique()
  
  # collapse duplicates: mean(Coverage) per Sample
  site.data <- site.data %>%
    group_by(Sample, site_name, site_type, Position, Cohort) %>%
    summarise(Coverage = mean(Coverage, na.rm = TRUE), .groups = "drop")
  site.data <- as.data.table(site.data)
  
  # wide format: Position × Sample
  gc.distances <- data.table::dcast(site.data,
                                    Position ~ Sample,
                                    value.var = "Coverage")
  gc.distances <- as.data.frame(gc.distances)
  colnames(gc.distances) <- gsub("Coverage\\.", "", colnames(gc.distances))
  
  # keep only valid samples in both places
  valid.samples <- intersect(colnames(gc.distances), all.samples)
  all.samples <- valid.samples
  
  # store this distance matrix
  metrics.per.site[[site]] <- data.frame(
    Sample            = valid.samples,
    Site              = rep(site, length(valid.samples)),
    # Mean.Coverage: average cfDNA coverage across the entire ±500 bp window
    # around the regulatory site centre - a global accessibility signal.
    Mean.Coverage     = colMeans(gc.distances[, valid.samples, drop = FALSE]),
    # Midpoint.Coverage: average coverage restricted to the five positions
    # nearest the site centre (±30 bp). This captures the nucleosome-free
    # region signal most directly relevant to chromatin accessibility.
    Midpoint.Coverage = colMeans(
      gc.distances[gc.distances$Position %in% c(-30, -15, 0, 15, 30), valid.samples, drop = FALSE]
    ),
    Midpoint.normalized = NA,
    # Amplitude: peak power from a periodogram of the coverage profile.
    # A high amplitude reflects the periodic nucleosome phasing pattern
    # (typically ~200 bp repeat) flanking an accessible site.
    Amplitude         = apply(GeneCycle::periodogram(gc.distances[, valid.samples, drop = FALSE])[["spec"]], 2, max)
  )
  
  # Midpoint.normalized: centre-relative accessibility score.
  # Subtracts global mean coverage so that values > 1 indicate the site centre
  # is MORE covered (more accessible) than the average flanking background.
  # Adding 1 keeps the scale positive and interpretable as a ratio around unity.
  metrics.per.site[[site]]$Midpoint.normalized <-
    (metrics.per.site[[site]]$Midpoint.Coverage - metrics.per.site[[site]]$Mean.Coverage) + 1
  
  # Z-scores for each metric: each sample's metric value is standardised
  # relative to the PON (healthy-control) distribution for this site.
  # Reference vector y is restricted to pon.samples, so tumour samples
  # are scored against the healthy baseline (not against themselves).
  # Positive z → more open chromatin than healthy; negative → less open.
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
  # Fold-change direction: group2 / group1 = tumour / healthy.
  # FC > 1 means the metric is elevated in tumour (open chromatin in MM);
  # FC < 1 means it is depleted relative to healthy controls.
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
  
  # baseline vs. healthy
  idx_baseline <- which(mdata$Sample %in% baseline.samples)
  stats.data[stats.data$Site == site, "PR.coverage.p"]   <-
    safe_ttest_p(mdata$Mean.Coverage[idx_baseline], mdata$Mean.Coverage[idx_pon])
  stats.data[stats.data$Site == site, "PR.midpoint.p"]   <-
    safe_ttest_p(mdata$Midpoint.normalized[idx_baseline],
                 mdata$Midpoint.normalized[idx_pon])
  stats.data[stats.data$Site == site, "PR.amplitude.p"]  <-
    safe_ttest_p(mdata$Amplitude[idx_baseline], mdata$Amplitude[idx_pon])
  
  # maintenance vs. healthy
  idx_maint <- which(mdata$Sample %in% maintenance.samples)
  stats.data[stats.data$Site == site, "PS.coverage.p"]   <-
    safe_ttest_p(mdata$Mean.Coverage[idx_maint], mdata$Mean.Coverage[idx_pon])
  stats.data[stats.data$Site == site, "PS.midpoint.p"]   <-
    safe_ttest_p(mdata$Midpoint.normalized[idx_maint],
                 mdata$Midpoint.normalized[idx_pon])
  stats.data[stats.data$Site == site, "PS.amplitude.p"]  <-
    safe_ttest_p(mdata$Amplitude[idx_maint], mdata$Amplitude[idx_pon])
  
  # baseline vs. maintenance → “response”
  p1 <- safe_ttest_p(mdata$Mean.Coverage[idx_baseline], mdata$Mean.Coverage[idx_maint])
  p2 <- safe_ttest_p(mdata$Midpoint.normalized[idx_baseline],
                     mdata$Midpoint.normalized[idx_maint])
  p3 <- safe_ttest_p(mdata$Amplitude[idx_baseline], mdata$Amplitude[idx_maint])
  stats.data[stats.data$Site == site, "response.p"] <- safe_min_p(p1, p2, p3)
  
  # M4_IMG vs. SPORE → “trial”
  idx_M4  <- which(mdata$Sample %in% M4_IMG.samples)
  idx_SPO <- which(mdata$Sample %in% SPORE.samples)
  p4 <- safe_ttest_p(mdata$Mean.Coverage[idx_M4], mdata$Mean.Coverage[idx_SPO])
  p5 <- safe_ttest_p(mdata$Midpoint.normalized[idx_M4],
                     mdata$Midpoint.normalized[idx_SPO])
  p6 <- safe_ttest_p(mdata$Amplitude[idx_M4], mdata$Amplitude[idx_SPO])
  stats.data[stats.data$Site == site, "trial.p"] <- safe_min_p(p4, p5, p6)
}


### ADJUST P-VALUES (FDR) ############################################################################
# Benjamini-Hochberg FDR correction applied across all genomic sites
# for each comparison type (Cohort, Baseline/PR, Maintenance/PS, Response,
# Trial). Corrects for the number of regulatory sites tested in parallel.
fdr.data <- data.frame(
  Cohort.Cov.fdr   = p.adjust(stats.data$Coverage.p,    "fdr"),
  Cohort.Mid.fdr   = p.adjust(stats.data$Midpoint.p,    "fdr"),
  Cohort.Amp.fdr   = p.adjust(stats.data$Amplitude.p,   "fdr"),
  PR.Cov.fdr       = p.adjust(stats.data$PR.coverage.p, "fdr"),
  PR.Mid.fdr       = p.adjust(stats.data$PR.midpoint.p, "fdr"),
  PR.Amp.fdr       = p.adjust(stats.data$PR.amplitude.p,"fdr"),
  PS.Cov.fdr       = p.adjust(stats.data$PS.coverage.p, "fdr"),
  PS.Mid.fdr       = p.adjust(stats.data$PS.midpoint.p, "fdr"),
  PS.Amp.fdr       = p.adjust(stats.data$PS.amplitude.p,"fdr"),
  Response.fdr     = p.adjust(stats.data$response.p,    "fdr"),
  Trial.fdr        = p.adjust(stats.data$trial.p,       "fdr")
)


### COMBINE METRICS + Z-SCORES & SAVE TSVs ##########################################################
results.data <- merge(
  do.call(rbind, metrics.per.site),
  do.call(rbind, scores.per.site),
  all = TRUE
)

# Derive HBC-based thresholds (±10% of max/min z-score among healthy samples)
# Strategy: take the maximum (or minimum) z-score observed across all
# healthy controls (TGL49 samples), then inflate by 10 % as a buffer.
# This defines the outer boundary of "normal" accessibility variation;
# tumour samples exceeding these bounds are flagged as aberrant.
# The 10 % buffer reduces false positives from healthy-control outliers.
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

# Directionality correction: if the tumour fold-change is NEGATIVE (metric
# is lower in tumour than in healthy controls), the relevant extreme is the
# minimum healthy z-score (thresholds.down) rather than the maximum.
# This ensures the binary flag always tests in the direction of the
# observed tumour-vs-healthy difference for each site.
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

# apply thresholds to classify each sample's z-score. A per-site lookup avoids
# ambiguous grouped-vector recycling and records NA when a comparison was
# underpowered or a threshold could not be estimated.
threshold_lookup <- stats.data %>%
  dplyr::select(
    Site,
    Coverage.fc, Coverage.threshold,
    Midpoint.fc, Midpoint.threshold,
    Amplitude.fc, Amplitude.threshold
  )

results.data <- results.data %>%
  dplyr::left_join(threshold_lookup, by = "Site") %>%
  dplyr::mutate(
    Threshold.Coverage = dplyr::case_when(
      is.na(Coverage.fc) | is.na(Coverage.threshold) ~ NA,
      Coverage.fc > 0 ~ Zscore.Coverage > Coverage.threshold,
      TRUE ~ Zscore.Coverage < Coverage.threshold
    ),
    Threshold.Midpoint = dplyr::case_when(
      is.na(Midpoint.fc) | is.na(Midpoint.threshold) ~ NA,
      Midpoint.fc > 0 ~ Zscore.Midpoint > Midpoint.threshold,
      TRUE ~ Zscore.Midpoint < Midpoint.threshold
    ),
    Threshold.Amplitude = dplyr::case_when(
      is.na(Amplitude.fc) | is.na(Amplitude.threshold) ~ NA,
      Amplitude.fc > 0 ~ Zscore.Amplitude > Amplitude.threshold,
      TRUE ~ Zscore.Amplitude < Amplitude.threshold
    )
  ) %>%
  dplyr::select(
    -Coverage.fc, -Coverage.threshold,
    -Midpoint.fc, -Midpoint.threshold,
    -Amplitude.fc, -Amplitude.threshold
  ) %>%
  dplyr::left_join(sample_metadata, by = "Sample")

# Save out per-site metrics + stats
dir.create(file.path(out.dir, "Griffin"), recursive = TRUE, showWarnings = FALSE)
write_tsv(results.data,
          file.path(out.dir, "Griffin/griffin_per_site_metrics.tsv"))

write_tsv(stats.data,
          file.path(out.dir, "Griffin/griffin_per_site_stats.tsv"))


### EXPORT “MM_DARs_chromatin_activation_data.csv” ################################################
# MM_DARs_chromatin_activation: a curated set of Differentially Accessible
# Regions (DARs) identified from ATAC-seq in multiple myeloma vs. normal
# plasma cells. Coverage at these sites in cfDNA reflects MM chromatin
# activity and is the primary GRIFFIN-based MRD feature used in this study.
# This export feeds directly into scripts 1_7B and 2_0.
MM_DARs_data <- results.data %>%
  filter(Site == "MM_DARs_chromatin_activation") %>%
  select(
    Sample, Site,
    Mean.Coverage, Midpoint.Coverage, Midpoint.normalized,
    Amplitude, Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude,
    Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
    Bam, Patient, Date_of_sample_collection
  ) %>%
  distinct()

write_csv(MM_DARs_data,
          file.path(out.dir, "MM_DARs_chromatin_activation_data.csv"))

# End of processing script
