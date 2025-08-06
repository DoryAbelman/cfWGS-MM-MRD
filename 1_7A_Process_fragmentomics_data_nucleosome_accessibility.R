source("setup_packages.R")
source("config.R")
source("helpers.R")

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
#   • input.dir ← "~/…/Fragmentomics_data"
#   • pon.dir   ← "~/…/Fragmentomics_data/Normals"
#   • out.dir   ← "~/…/Results_Fragmentomics/"
#   • Clinical CSVs in working dir:
#        – combined_clinical_data_updated_April2025.csv
#        – Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated3.csv
#
# Author: Dory Abelman
# Last updated: 2025-05-23
# ──────────────────────────────────────────────────────────────────────────────


### PREPARE SESSION ################################################################################

source('session.functions.R')

input.dir <- '~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data'
pon.dir   <- '~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Fragmentomics_data/Normals'
out.dir   <- '~/Documents/Thesis_work/R/M4/Projects/High_risk_MM_baselinbe_relapse_marrow/Results_Fragmentomics'

if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

### Load clinical tables ############################################################################
combined_clinical_data_updated <- read_csv("combined_clinical_data_updated_April2025.csv")

### READ NUCLEOSOME-ACCESSIBILITY DISTANCE FILES ####################################################
results.files <- list.files(path = input.dir,
                            pattern = "nucleosome_accessibility_distances.tsv",
                            full.names = TRUE)
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
    )
  )

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
  
  # wide format: Position × Sample
  gc.distances <- dcast(site.data,
                        Position ~ Sample,
                        value.var = "Coverage")
  colnames(gc.distances) <- gsub("Coverage\\.", "", colnames(gc.distances))
  
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
  
  # baseline vs. healthy
  idx_baseline <- which(mdata$Sample %in% baseline.samples)
  stats.data[stats.data$Site == site, "PR.coverage.p"]   <-
    t.test(mdata$Mean.Coverage[idx_baseline], mdata$Mean.Coverage[idx_pon])$p.value
  stats.data[stats.data$Site == site, "PR.midpoint.p"]   <-
    t.test(mdata$Midpoint.normalized[idx_baseline],
           mdata$Midpoint.normalized[idx_pon])$p.value
  stats.data[stats.data$Site == site, "PR.amplitude.p"]  <-
    t.test(mdata$Amplitude[idx_baseline], mdata$Amplitude[idx_pon])$p.value
  
  # maintenance vs. healthy
  idx_maint <- which(mdata$Sample %in% maintenance.samples)
  stats.data[stats.data$Site == site, "PS.coverage.p"]   <-
    t.test(mdata$Mean.Coverage[idx_maint], mdata$Mean.Coverage[idx_pon])$p.value
  stats.data[stats.data$Site == site, "PS.midpoint.p"]   <-
    t.test(mdata$Midpoint.normalized[idx_maint],
           mdata$Midpoint.normalized[idx_pon])$p.value
  stats.data[stats.data$Site == site, "PS.amplitude.p"]  <-
    t.test(mdata$Amplitude[idx_maint], mdata$Amplitude[idx_pon])$p.value
  
  # baseline vs. maintenance → “response”
  p1 <- t.test(mdata$Mean.Coverage[idx_baseline], mdata$Mean.Coverage[idx_maint])$p.value
  p2 <- t.test(mdata$Midpoint.normalized[idx_baseline],
               mdata$Midpoint.normalized[idx_maint])$p.value
  p3 <- t.test(mdata$Amplitude[idx_baseline], mdata$Amplitude[idx_maint])$p.value
  stats.data[stats.data$Site == site, "response.p"] <- min(c(p1, p2, p3))
  
  # M4_IMG vs. SPORE → “trial”
  idx_M4  <- which(mdata$Sample %in% M4_IMG.samples)
  idx_SPO <- which(mdata$Sample %in% SPORE.samples)
  p4 <- t.test(mdata$Mean.Coverage[idx_M4], mdata$Mean.Coverage[idx_SPO])$p.value
  p5 <- t.test(mdata$Midpoint.normalized[idx_M4],
               mdata$Midpoint.normalized[idx_SPO])$p.value
  p6 <- t.test(mdata$Amplitude[idx_M4], mdata$Amplitude[idx_SPO])$p.value
  stats.data[stats.data$Site == site, "trial.p"] <- min(c(p4, p5, p6))
}


### ADJUST P-VALUES (FDR) ############################################################################
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
write_tsv(results.data,
          file.path(out.dir, "Griffin/griffin_per_site_metrics.tsv"))

write_tsv(stats.data,
          file.path(out.dir, "Griffin/griffin_per_site_stats.tsv"))


### EXPORT “MM_DARs_chromatin_activation_data.csv” ################################################
# Filter only the “MM_DARs_chromatin_activation” site from results.data
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
