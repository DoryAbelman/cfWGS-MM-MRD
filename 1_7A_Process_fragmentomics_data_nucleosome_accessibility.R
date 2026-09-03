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
#   • Adjust p‐values (FDR) and place thresholds 10% beyond the observed
#     healthy-control z-score extremes.
#   • Save:
#        – “Results_Fragmentomics/Griffin/griffin_per_site_metrics.tsv”
#        – “Results_Fragmentomics/Griffin/griffin_per_site_stats.tsv”
#        – “Results_Fragmentomics/fragmentomics_platform_reference_thresholds.csv”
#        – “Results_Fragmentomics/fragmentomics_threshold_direction_sensitivity_audit.csv”
#   • Build and save “MM_DARs_chromatin_activation_data.csv/.rds”.
#
# Assumptions & file paths:
#   • input.dir <- "Fragmentomics_data"
#   • pon.dir   <- "Fragmentomics_data/Normals"
#   • out.dir   <- "Results_Fragmentomics"
#   • Clinical CSVs in working dir:
#        – combined_clinical_data_updated_April2025.csv
#   • Optional Spring 2026 patient and matched-control GRIFFIN exports are
#     resolved by helpers.R when present.
#   • session.functions.R and Scripts_2025/Final_Scripts/helpers.R
#
# Analysis units and statistical scope:
#   Metric tables contain one sample-site row. Site-level t-tests treat sample
#   rows as independent, so patients with repeated longitudinal specimens can
#   contribute more than once; these tests are exploratory and are not
#   patient-clustered or paired. `response.p` and `trial.p` are the unadjusted
#   minimum of three metric-specific p-values before across-site BH correction,
#   so their reported FDR values do not control the within-site metric search.
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
# Reproducibility boundary:
#   The historical 26-control reference and XPlus 19-control reference are
#   intentionally frozen and platform-specific. They include the material types
#   selected by the helper filters, not necessarily plasma cfDNA only.
#
# Threshold direction:
#   Fold change is a tumour/healthy ratio. Values > 1 therefore use the upper
#   healthy-control boundary and values <= 1 use the lower boundary. The prior
#   >= 0 / > 0 implementation is retained only in the sensitivity audit.
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
historical_results_files <- list.files(
  path = input.dir,
  pattern = "nucleosome_accessibility_distances.tsv",
  full.names = TRUE
)
revision_results_files <- spring2026_revision_files(
  "Fragmentomics_Pipeeline_Suite_all_outputs",
  "^2026-06-25_cfWGS_MM_fragmentomics_Revisions_Spring2026_nucleosome_accessibility_distances[.]tsv$"
)
pon.files <- rev(sort(
  list.files(path = pon.dir,
             pattern = "nucleosome_accessibility_distances.tsv",
             full.names = TRUE)
))
xplus_pon_files <- spring2026_revision_files(
  "Fragmentomics_Pipeeline_Suite_all_outputs",
  "^2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_nucleosome_accessibility_distances[.]tsv$"
)

if (!length(pon.files) || !length(xplus_pon_files)) {
  stop("Both historical and XPlus nucleosome healthy-control inputs are required.", call. = FALSE)
}

read_nucleosome_batch <- function(paths, platform, role) {
  if (!length(paths)) return(tibble())
  bind_rows(lapply(paths, function(path) {
    read_tsv(path, show_col_types = FALSE) %>%
      mutate(
        fragmentomics_sequencing_platform = platform,
        fragmentomics_sample_role = role,
        fragmentomics_source_file = basename(path)
      )
  }))
}

## Read cfWGS data
cfWGS.data <- bind_rows(
  read_nucleosome_batch(historical_results_files, "NovaSeq 6000", "patient"),
  read_nucleosome_batch(revision_results_files, "NovaSeq XPlus", "patient_revision")
) %>%
  filter(
    !(.data$fragmentomics_sample_role == "patient_revision" &
        is_spring2026_revision_primary_analysis_excluded(.data$Sample))
  )

# If TGL49 columns sneaked in as columns, drop them
if (any(grepl("TGL49", colnames(cfWGS.data)))) {
  cfWGS.data <- cfWGS.data[, !grepl("TGL49", colnames(cfWGS.data))]
}

# Preserve the original 26-library NovaSeq 6000 reference exactly as used by
# the locked historical analysis. The XPlus reference is restricted to the 19
# identity-matched controls available on both platforms. The matched historical
# subset is used later only for platform-harmonization parameter estimation;
# it must not replace the original historical scoring reference.
classify_fragmentomics_material <- function(sample) {
  case_when(
    str_detect(sample, "(^|[-_])Cf([-_]|$)|_Cf_|-Cf-") ~ "cfDNA",
    str_detect(sample, "(^|[-_])Pb([-_]|$)|_Pb_|-Pb-") ~ "PB_cellular_or_buffy",
    str_detect(sample, "(^|[-_])Ct([-_]|$)|_Ct_|-Ct-") ~ "cellular_or_tissue",
    TRUE ~ "unknown"
  )
}

historical_pon_data <- read_nucleosome_batch(
  pon.files[[1]], "NovaSeq 6000", "healthy_control"
) %>%
  mutate(fragmentomics_material = classify_fragmentomics_material(.data$Sample)) %>%
  filter_fragmentomics_original_historical_controls()

xplus_pon_data <- read_nucleosome_batch(
  xplus_pon_files[[1]], "NovaSeq XPlus", "healthy_control"
) %>%
  mutate(fragmentomics_material = classify_fragmentomics_material(.data$Sample)) %>%
  filter_fragmentomics_shared_healthy_controls()

pon.data <- bind_rows(historical_pon_data, xplus_pon_data)


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
spring2026_pwgval_metadata <- load_spring2026_pwgval_dilution_metadata(required = FALSE)
if (!is.null(spring2026_pwgval_metadata)) {
  tmp <- bind_rows(
    tmp,
    spring2026_pwgval_metadata %>%
      transmute(
        Merge_ID = paste0(.data$Group_ID, "-P"),
        Bam = .data$Bam,
        Patient = .data$Patient,
        Sample_ID = .data$Sample_ID,
        Date_of_sample_collection = as.Date(NA),
        Timepoint = as.character(.data$Timepoint),
        timepoint_info = .data$timepoint_info,
        Cohort = "PWGVAL_M4CHIP"
      )
  )
}
tmp <- bind_rows(
  tmp,
  tmp %>% mutate(Merge_ID = str_replace(Merge_ID, "^IMG", "MyP"))
) %>%
  distinct(Merge_ID, Bam, Patient, Sample_ID, .keep_all = TRUE)

# One fragmentomics sample must map to one clinical row. Retaining multiple BAM
# metadata rows here multiplies every per-position GRIFFIN record before the site
# loop and can silently average technical/platform replicates.
tmp_join <- tmp %>%
  arrange(.data$Merge_ID, desc(!is.na(.data$Date_of_sample_collection))) %>%
  distinct(.data$Merge_ID, .keep_all = TRUE)

sample_metadata <- tmp_join %>%
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
  left_join(tmp_join, by = c("Sample" = "Merge_ID"), relationship = "many-to-one")

pon.data <- pon.data %>% mutate(Cohort = "HBC")


### BUILD SAMPLE LISTS FOR GROUP COMPARISONS #########################################################
tumour.samples     <- unique(cfWGS.data$Sample)
baseline.samples   <- unique(cfWGS.data %>% filter(timepoint_info %in% c("Baseline", "Diagnosis"))   %>% pull(Sample))
maintenance.samples<- unique(cfWGS.data %>% filter(timepoint_info %in% c("Post_induction", "Post_transplant",
                                                                         "Maintenance", "1.5yr maintenance",
                                                                         "2yr maintenance"))         %>% pull(Sample))
M4_IMG.samples     <- tumour.samples[!grepl("SPORE", tumour.samples)]
SPORE.samples      <- tumour.samples[grepl("SPORE", tumour.samples)]
pon.samples        <- unique(pon.data$Sample)

sample_reference_map <- bind_rows(
  cfWGS.data %>%
    distinct(Sample, fragmentomics_sequencing_platform, fragmentomics_sample_role),
  pon.data %>%
    distinct(Sample, fragmentomics_sequencing_platform, fragmentomics_sample_role)
) %>%
  distinct(Sample, .keep_all = TRUE)
sample_platform <- setNames(
  sample_reference_map$fragmentomics_sequencing_platform,
  sample_reference_map$Sample
)
sample_role <- setNames(sample_reference_map$fragmentomics_sample_role, sample_reference_map$Sample)

all.samples <- c(tumour.samples, pon.samples)
all.sites   <- unique(cfWGS.data$site_name)


### UTILITY FUNCTIONS ##############################################################################
# calculate.zscore: standard z-score of a single value x relative to reference
# vector y (the PON/healthy-control distribution for this site & metric).
# Formula: z = (x - mean_PON) / sd_PON
# A positive z means this particular metric is higher than the platform-matched
# healthy-control mean. Accessibility interpretation depends on the metric and
# the site-specific tumour-versus-control direction calculated below.
calculate.zscore <- function(x, y) {
  y <- y[is.finite(y)]
  if (length(y) < 2L || !is.finite(stats::sd(y)) || stats::sd(y) == 0) return(NA_real_)
  (x - mean(y)) / sd(y)
}

calculate_platform_zscores <- function(values, samples, platforms, roles) {
  vapply(seq_along(values), function(i) {
    reference <- values[
      roles == "healthy_control" & platforms == platforms[[i]]
    ]
    calculate.zscore(values[[i]], reference)
  }, numeric(1))
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
    fragmentomics_sequencing_platform = unname(sample_platform[valid.samples]),
    fragmentomics_sample_role = unname(sample_role[valid.samples]),
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
  # is more covered than the average window background. Adding 1 centers the
  # difference at unity; despite that scale, this is not a ratio.
  metrics.per.site[[site]]$Midpoint.normalized <-
    (metrics.per.site[[site]]$Midpoint.Coverage - metrics.per.site[[site]]$Mean.Coverage) + 1
  
  # Z-scores for each metric: each sample's metric value is standardised
  # relative to the PON (healthy-control) distribution for this site.
  # Reference vector y is restricted to pon.samples, so tumour samples
  # are scored against the healthy baseline (not against themselves).
  # Positive z means a higher metric value than healthy; it does not by itself
  # establish greater biological accessibility for every metric and site.
  metric_platform <- metrics.per.site[[site]]$fragmentomics_sequencing_platform
  metric_role <- metrics.per.site[[site]]$fragmentomics_sample_role
  score.data <- data.frame(
    Sample        = valid.samples,
    Site          = rep(site, length(valid.samples)),
    Zscore.Coverage = calculate_platform_zscores(
      metrics.per.site[[site]]$Mean.Coverage,
      valid.samples,
      metric_platform,
      metric_role
    ),
    Zscore.Midpoint = calculate_platform_zscores(
      metrics.per.site[[site]]$Midpoint.normalized,
      valid.samples,
      metric_platform,
      metric_role
    ),
    Zscore.Amplitude = calculate_platform_zscores(
      metrics.per.site[[site]]$Amplitude,
      valid.samples,
      metric_platform,
      metric_role
    )
  )
  scores.per.site[[site]] <- score.data
  
  # group comparisons (fold-change & p-value)
  # Fold-change direction: group2 / group1 = tumour / healthy.
  # FC > 1 means the metric is elevated in tumour; FC < 1 means it is depleted.
  # The accessibility meaning of that direction is metric- and site-dependent.
  mdata <- metrics.per.site[[site]]
  # tumour vs. healthy
  # Preserve the original feature-discovery comparison: only historical tumour
  # samples and historical healthy controls determine directionality/FDR. The
  # XPlus revision set remains an external test-cohort expansion.
  idx_pon <- which(
    mdata$fragmentomics_sample_role == "healthy_control" &
      mdata$fragmentomics_sequencing_platform == "NovaSeq 6000"
  )
  idx_tumor <- which(
    mdata$fragmentomics_sample_role == "patient" &
      mdata$fragmentomics_sequencing_platform == "NovaSeq 6000"
  )
  stats.data[stats.data$Site == site, c("Coverage.fc", "Coverage.p")] <-
    get.ttest.p.and.foldchange(mdata$Mean.Coverage, idx_pon, idx_tumor)
  stats.data[stats.data$Site == site, c("Midpoint.fc", "Midpoint.p")] <-
    get.ttest.p.and.foldchange(mdata$Midpoint.normalized, idx_pon, idx_tumor)
  stats.data[stats.data$Site == site, c("Amplitude.fc", "Amplitude.p")] <-
    get.ttest.p.and.foldchange(mdata$Amplitude, idx_pon, idx_tumor)
  
  # baseline vs. healthy
  idx_baseline <- which(
    mdata$Sample %in% baseline.samples &
      mdata$fragmentomics_sample_role == "patient" &
      mdata$fragmentomics_sequencing_platform == "NovaSeq 6000"
  )
  stats.data[stats.data$Site == site, "PR.coverage.p"]   <-
    safe_ttest_p(mdata$Mean.Coverage[idx_baseline], mdata$Mean.Coverage[idx_pon])
  stats.data[stats.data$Site == site, "PR.midpoint.p"]   <-
    safe_ttest_p(mdata$Midpoint.normalized[idx_baseline],
                 mdata$Midpoint.normalized[idx_pon])
  stats.data[stats.data$Site == site, "PR.amplitude.p"]  <-
    safe_ttest_p(mdata$Amplitude[idx_baseline], mdata$Amplitude[idx_pon])
  
  # maintenance vs. healthy
  idx_maint <- which(
    mdata$Sample %in% maintenance.samples &
      mdata$fragmentomics_sample_role == "patient" &
      mdata$fragmentomics_sequencing_platform == "NovaSeq 6000"
  )
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
  idx_M4 <- which(
    mdata$Sample %in% M4_IMG.samples &
      mdata$fragmentomics_sample_role == "patient" &
      mdata$fragmentomics_sequencing_platform == "NovaSeq 6000"
  )
  idx_SPO <- which(
    mdata$Sample %in% SPORE.samples &
      mdata$fragmentomics_sample_role == "patient" &
      mdata$fragmentomics_sequencing_platform == "NovaSeq 6000"
  )
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

# Derive platform-specific HBC thresholds. Each healthy control is first scored
# against controls from its own platform; the outer healthy z-score boundary is
# therefore estimated independently for NovaSeq 6000 and NovaSeq XPlus.
platform_thresholds <- results.data %>%
  filter(.data$fragmentomics_sample_role == "healthy_control") %>%
  group_by(.data$Site, .data$fragmentomics_sequencing_platform) %>%
  summarise(
    Coverage.threshold.up = max(.data$Zscore.Coverage, na.rm = TRUE) * 1.10,
    Coverage.threshold.down = min(.data$Zscore.Coverage, na.rm = TRUE) * 1.10,
    Midpoint.threshold.up = max(.data$Zscore.Midpoint, na.rm = TRUE) * 1.10,
    Midpoint.threshold.down = min(.data$Zscore.Midpoint, na.rm = TRUE) * 1.10,
    Amplitude.threshold.up = max(.data$Zscore.Amplitude, na.rm = TRUE) * 1.10,
    Amplitude.threshold.down = min(.data$Zscore.Amplitude, na.rm = TRUE) * 1.10,
    n_healthy_reference = n_distinct(.data$Sample),
    .groups = "drop"
  ) %>%
  left_join(
    stats.data %>%
      select(.data$Site, .data$Coverage.fc, .data$Midpoint.fc, .data$Amplitude.fc),
    by = "Site"
  ) %>%
  mutate(
    # Fold change is a ratio: > 1 indicates elevation in tumour and therefore
    # uses the upper reference boundary; <= 1 uses the lower boundary.
    Coverage.threshold = if_else(.data$Coverage.fc > 1, .data$Coverage.threshold.up, .data$Coverage.threshold.down),
    Midpoint.threshold = if_else(.data$Midpoint.fc > 1, .data$Midpoint.threshold.up, .data$Midpoint.threshold.down),
    Amplitude.threshold = if_else(.data$Amplitude.fc > 1, .data$Amplitude.threshold.up, .data$Amplitude.threshold.down)
  )

historical_thresholds <- platform_thresholds %>%
  filter(.data$fragmentomics_sequencing_platform == "NovaSeq 6000") %>%
  select(.data$Site, .data$Coverage.threshold, .data$Midpoint.threshold, .data$Amplitude.threshold)
stats.data <- stats.data %>%
  left_join(historical_thresholds, by = "Site")

write_csv(
  platform_thresholds,
  file.path(out.dir, "fragmentomics_platform_reference_thresholds.csv")
)

results.data <- results.data %>%
  dplyr::left_join(
    platform_thresholds %>%
      dplyr::select(
        .data$Site, .data$fragmentomics_sequencing_platform,
        .data$Coverage.fc, .data$Coverage.threshold,
        .data$Midpoint.fc, .data$Midpoint.threshold,
        .data$Amplitude.fc, .data$Amplitude.threshold,
        .data$n_healthy_reference
      ),
    by = c("Site", "fragmentomics_sequencing_platform")
  ) %>%
  dplyr::mutate(
    Threshold.Coverage = dplyr::case_when(
      is.na(Coverage.fc) | is.na(Coverage.threshold) ~ NA,
      Coverage.fc > 1 ~ Zscore.Coverage > Coverage.threshold,
      TRUE ~ Zscore.Coverage < Coverage.threshold
    ),
    Threshold.Midpoint = dplyr::case_when(
      is.na(Midpoint.fc) | is.na(Midpoint.threshold) ~ NA,
      Midpoint.fc > 1 ~ Zscore.Midpoint > Midpoint.threshold,
      TRUE ~ Zscore.Midpoint < Midpoint.threshold
    ),
    Threshold.Amplitude = dplyr::case_when(
      is.na(Amplitude.fc) | is.na(Amplitude.threshold) ~ NA,
      Amplitude.fc > 1 ~ Zscore.Amplitude > Amplitude.threshold,
      TRUE ~ Zscore.Amplitude < Amplitude.threshold
    )
  ) %>%
  dplyr::select(
    -Coverage.fc, -Coverage.threshold,
    -Midpoint.fc, -Midpoint.threshold,
    -Amplitude.fc, -Amplitude.threshold
  ) %>%
  dplyr::left_join(sample_metadata, by = "Sample")

# Audit the corrected active rule against the historical direction rule.
# `Threshold.*` above uses the ratio-appropriate split at 1. The legacy audit
# columns reproduce the prior > 0 upper-tail implementation without feeding it
# back into the active table. Only rows whose binary call changes are exported.
threshold_direction_sensitivity_audit <- results.data %>%
  left_join(
    platform_thresholds %>%
      select(
        Site, fragmentomics_sequencing_platform,
        Coverage.fc,
        Coverage.threshold.up, Coverage.threshold.down,
        Midpoint.fc,
        Midpoint.threshold.up, Midpoint.threshold.down,
        Amplitude.fc,
        Amplitude.threshold.up, Amplitude.threshold.down
      ),
    by = c("Site", "fragmentomics_sequencing_platform")
  ) %>%
  mutate(
    Threshold.Coverage.legacy = case_when(
      is.na(Coverage.fc) | is.na(Zscore.Coverage) ~ NA,
      Coverage.fc > 0 ~ Zscore.Coverage > Coverage.threshold.up,
      TRUE ~ Zscore.Coverage < Coverage.threshold.down
    ),
    Threshold.Midpoint.legacy = case_when(
      is.na(Midpoint.fc) | is.na(Zscore.Midpoint) ~ NA,
      Midpoint.fc > 0 ~ Zscore.Midpoint > Midpoint.threshold.up,
      TRUE ~ Zscore.Midpoint < Midpoint.threshold.down
    ),
    Threshold.Amplitude.legacy = case_when(
      is.na(Amplitude.fc) | is.na(Zscore.Amplitude) ~ NA,
      Amplitude.fc > 0 ~ Zscore.Amplitude > Amplitude.threshold.up,
      TRUE ~ Zscore.Amplitude < Amplitude.threshold.down
    ),
    Coverage.call_changed = Threshold.Coverage != Threshold.Coverage.legacy,
    Midpoint.call_changed = Threshold.Midpoint != Threshold.Midpoint.legacy,
    Amplitude.call_changed = Threshold.Amplitude != Threshold.Amplitude.legacy
  ) %>%
  filter(
    coalesce(Coverage.call_changed, FALSE) |
      coalesce(Midpoint.call_changed, FALSE) |
      coalesce(Amplitude.call_changed, FALSE)
  ) %>%
  select(
    Sample, Patient, Site,
    fragmentomics_sequencing_platform, fragmentomics_sample_role,
    starts_with("Zscore."),
    Coverage.fc, Coverage.threshold.up, Coverage.threshold.down,
    Threshold.Coverage.legacy, Threshold.Coverage,
    Coverage.call_changed,
    Midpoint.fc, Midpoint.threshold.up, Midpoint.threshold.down,
    Threshold.Midpoint.legacy, Threshold.Midpoint,
    Midpoint.call_changed,
    Amplitude.fc, Amplitude.threshold.up, Amplitude.threshold.down,
    Threshold.Amplitude.legacy, Threshold.Amplitude,
    Amplitude.call_changed
  )

write_csv(
  threshold_direction_sensitivity_audit,
  file.path(out.dir, "fragmentomics_threshold_direction_sensitivity_audit.csv")
)

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
    fragmentomics_sequencing_platform, fragmentomics_sample_role,
    n_healthy_reference, Bam, Patient, Date_of_sample_collection
  ) %>%
  distinct()

write_csv(MM_DARs_data,
          file.path(out.dir, "MM_DARs_chromatin_activation_data.csv"))

saveRDS(MM_DARs_data,
        file.path(out.dir, "MM_DARs_chromatin_activation_data.rds"))

# End of processing script
