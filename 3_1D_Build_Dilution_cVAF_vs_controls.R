# ==============================================================================
# Build a shareable cVAF-only dilution-series comparison panel
#
# Goal:
#   Plot all four controlled dilution series using raw cumulative VAF (cVAF),
#   followed by platform-matched healthy-control measurements and unmatched
#   plasma negative-control pairings from the revision all-by-all MRDetect run.
#
# Scientific interpretation:
#   - Dilution points are patient-specific trajectories and technical replicates.
#   - Healthy controls are evaluated against the four dilution-patient panels.
#   - Unmatched plasma points are negative BAM-panel pairings, not biological
#     longitudinal observations. The XPlus component is pooled across the 17
#     baseline BM mutation panels included in the revision all-by-all run.
#   - Controls are excluded from dilution-series correlation calculations.
#
# Inputs:
#   Output_tables_2025/Source_Data_Extended_Data/
#     SourceData_Fig4G_BM_patient_series_points_zero_final_tf.csv
#     SourceData_ED5D_BM_healthy_control_background_after_zero.csv
#   MRDetect_output_winter_2025/Processed_R_outputs/
#     cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.rds
#
# Outputs:
#   Final Tables and Figures/
#     Dilution_series_cVAF_vs_healthy_and_unmatched_plasma.png
#     Dilution_series_cVAF_vs_healthy_and_unmatched_plasma.pdf
#     Dilution_series_cVAF_vs_healthy_and_unmatched_plasma_pseudolog.png/.pdf
#     Dilution_series_cVAF_vs_healthy_and_unmatched_plasma_pseudolog_low_LOD_zoom.png/.pdf
#     Dilution_series_BM_cVAF_and_sites_zscores_vs_controls.png/.pdf
#     Dilution_series_BM_cVAF_and_sites_zscores_vs_controls_pseudolog.png/.pdf
#   Output_tables_2025/Source_Data_Extended_Data/
#     SourceData_Dilution_cVAF_vs_healthy_and_unmatched_plasma.csv
#     SourceData_Dilution_cVAF_vs_healthy_and_unmatched_plasma_summary.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

set.seed(20260802)

source_dir <- file.path("Output_tables_2025", "Source_Data_Extended_Data")
figure_dir <- "Final Tables and Figures"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dilution_path <- file.path(
  source_dir,
  "SourceData_Fig4G_BM_patient_series_points_zero_final_tf.csv"
)
healthy_path <- file.path(
  source_dir,
  "SourceData_ED5D_BM_healthy_control_background_after_zero.csv"
)
all_mrdetect_path <- file.path(
  "MRDetect_output_winter_2025", "Processed_R_outputs",
  "cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.rds"
)
high_quality_panel_path <- file.path(
  "Output_tables_2025", "high_quality_patients_list_for_baseline_mut_calling2.csv"
)
sample_manifest_path <- file.path(
  "Output_tables_2025", "clinical_support", "sample_scoring_status_manifest.csv"
)

required_files <- c(
  dilution_path, healthy_path, all_mrdetect_path,
  high_quality_panel_path, sample_manifest_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

assert_columns <- function(dat, required, label) {
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

# Positive dilution-series points. The existing source has the three XPlus 0%
# controls expanded once per replicate line, which is intentional for drawing
# both technical-replicate trajectories.
dilution_raw <- read_csv(dilution_path, show_col_types = FALSE)
assert_columns(
  dilution_raw,
  c("Patient", "Sample_ID", "Bam", "LOD_zero_final", "replicate_id",
    "line_group", "feature", "value"),
  "Dilution source"
)

dilution <- dilution_raw %>%
  filter(.data$feature == "detect_rate_BM", !is.na(.data$value)) %>%
  transmute(
    record_type = "Dilution series",
    patient = .data$Patient,
    query_bam = basename(.data$Bam),
    mutation_panel = .data$Patient,
    platform = if_else(.data$Patient == "M4-28", "NovaSeq 6000", "NovaSeq XPlus"),
    tumor_fraction_percent = .data$LOD_zero_final,
    x_position = if_else(.data$LOD_zero_final > 0, -log10(.data$LOD_zero_final), 5),
    cvaf = .data$value,
    replicate_id = as.character(.data$replicate_id),
    line_group = .data$line_group,
    comparator_scope = "Patient-specific dilution trajectory"
  )

expected_patients <- c("M4-28", "M4-34", "M4-37", "M4-38")
if (!setequal(unique(dilution$patient), expected_patients)) {
  stop("Dilution source does not contain exactly the expected four de-identified patients.", call. = FALSE)
}
if (nrow(distinct(dilution, .data$patient, .data$query_bam)) != 48L) {
  stop("Expected 48 unique dilution patient-BAM combinations.", call. = FALSE)
}

# Healthy-control measurements are already panel- and platform-matched by the
# upstream dilution workflow (26 historical measurements for M4-28 and 22
# XPlus measurements for each of M4-34/M4-37/M4-38).
healthy_raw <- read_csv(healthy_path, show_col_types = FALSE)
assert_columns(
  healthy_raw,
  c("Patient", "control_bam", "VCF_clean", "healthy_reference_tier", "feature", "value"),
  "Healthy-control source"
)

healthy <- healthy_raw %>%
  filter(.data$feature == "detect_rate_BM", !is.na(.data$value)) %>%
  transmute(
    record_type = "Healthy controls",
    patient = .data$Patient,
    query_bam = .data$control_bam,
    mutation_panel = .data$VCF_clean,
    platform = if_else(
      .data$healthy_reference_tier == "XPLUS_CHARM_healthy",
      "NovaSeq XPlus", "NovaSeq 6000"
    ),
    tumor_fraction_percent = NA_real_,
    x_position = 6,
    cvaf = .data$value,
    replicate_id = NA_character_,
    line_group = NA_character_,
    comparator_scope = "Platform- and dilution-panel-matched healthy BAM-panel pairing"
  )

if (nrow(healthy) != 92L || n_distinct(healthy$query_bam) != 48L) {
  stop("Expected 92 healthy BAM-panel measurements from 48 unique control BAMs.", call. = FALSE)
}

# Revision XPlus all-by-all negative controls. These are unrelated plasma BAM
# and baseline BM mutation-panel pairs. Require both (1) a mutation-source
# patient in the canonical high-quality baseline helper and (2) a query plasma
# row with an evaluable blood cfWGS score in the current scoring manifest.
# Restricting to available XPlus healthy normalization prevents legacy/revision
# platform mixing in this component.
all_mrdetect <- readRDS(all_mrdetect_path)
assert_columns(
  all_mrdetect,
  c("BAM", "VCF_clean", "Patient", "Patient_Bam", "plotting_type",
    "revision_new_xplus_bam_Bam", "Sample_type_Bam", "Mut_source",
    "Filter_source", "timepoint_info", "Study_Bam", "detection_rate",
    "healthy_reference_requested", "healthy_reference_status"),
  "All-patient MRDetect object"
)

high_quality_panels <- read_csv(high_quality_panel_path, show_col_types = FALSE)
assert_columns(high_quality_panels, "Patient", "High-quality panel helper")
eligible_panel_patients <- unique(high_quality_panels$Patient)

sample_manifest <- read_csv(sample_manifest_path, show_col_types = FALSE)
assert_columns(
  sample_manifest,
  c("Sample_Code", "has_blood_cfWGS_score"),
  "Sample scoring manifest"
)
evaluable_plasma_samples <- sample_manifest %>%
  filter(.data$has_blood_cfWGS_score %in% TRUE) %>%
  distinct(.data$Sample_Code)

unmatched <- all_mrdetect %>%
  filter(
    .data$plotting_type == "Unmatched_plasma",
    .data$revision_new_xplus_bam_Bam %in% TRUE,
    .data$Sample_type_Bam == "Blood_plasma_cfDNA",
    .data$Mut_source == "BM_cells",
    .data$Filter_source == "STR_encode",
    .data$timepoint_info == "Baseline",
    .data$Study_Bam == "IMMAGINE_revision_OICR",
    .data$healthy_reference_requested == "XPLUS_CHARM_healthy",
    .data$healthy_reference_status == "available",
    .data$Patient %in% .env$eligible_panel_patients,
    !is.na(.data$detection_rate)
  ) %>%
  distinct(.data$BAM, .data$VCF_clean, .keep_all = TRUE) %>%
  mutate(query_sample_code = str_remove(.data$Sample_ID_Bam, "-P$")) %>%
  semi_join(
    evaluable_plasma_samples,
    by = c("query_sample_code" = "Sample_Code")
  ) %>%
  transmute(
    record_type = "Unmatched plasma",
    patient = NA_character_,
    query_bam = basename(.data$BAM),
    mutation_panel = .data$VCF_clean,
    platform = "NovaSeq XPlus",
    tumor_fraction_percent = NA_real_,
    x_position = 7,
    cvaf = .data$detection_rate,
    replicate_id = NA_character_,
    line_group = NA_character_,
    comparator_scope = paste(
      "Revision XPlus all-by-all unrelated plasma BAM-baseline BM panel pairing;",
      "high-quality mutation-source patient and evaluable query plasma required"
    )
  )

if (nrow(unmatched) != 725L || n_distinct(unmatched$query_bam) != 54L ||
    n_distinct(unmatched$mutation_panel) != 14L) {
  stop(
    "Unexpected evaluable revision unmatched-control denominator: expected 725 pairings, 54 plasma BAMs, and 14 panels.",
    call. = FALSE
  )
}

source_data <- bind_rows(dilution, healthy, unmatched) %>%
  mutate(record_type = factor(
    .data$record_type,
    levels = c("Dilution series", "Healthy controls", "Unmatched plasma")
  ))

write_csv(
  source_data,
  file.path(source_dir, "SourceData_Dilution_cVAF_vs_healthy_and_unmatched_plasma.csv")
)

summary_data <- source_data %>%
  group_by(.data$record_type, .data$platform) %>%
  summarise(
    n_measurements = n(),
    n_query_bams = n_distinct(.data$query_bam),
    n_mutation_panels = n_distinct(.data$mutation_panel),
    median_cvaf = median(.data$cvaf, na.rm = TRUE),
    q1_cvaf = quantile(.data$cvaf, 0.25, na.rm = TRUE),
    q3_cvaf = quantile(.data$cvaf, 0.75, na.rm = TRUE),
    min_cvaf = min(.data$cvaf, na.rm = TRUE),
    max_cvaf = max(.data$cvaf, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(
  summary_data,
  file.path(source_dir, "SourceData_Dilution_cVAF_vs_healthy_and_unmatched_plasma_summary.csv")
)

patient_colors <- c(
  "M4-28" = "#F8766D",
  "M4-34" = "#7CAE00",
  "M4-37" = "#00BFC4",
  "M4-38" = "#C77CFF"
)

control_plot_data <- bind_rows(healthy, unmatched) %>%
  mutate(
    control_group = factor(
      as.character(.data$record_type),
      levels = c("Healthy controls", "Unmatched plasma")
    )
  )

# Distributions summarize all measurements while the translucent points retain
# the observed spread. Jitter is deterministic because the seed is fixed.
p <- ggplot() +
  geom_line(
    data = dilution,
    aes(x = .data$x_position, y = .data$cvaf,
        group = .data$line_group, color = .data$patient),
    linewidth = 0.65,
    alpha = 0.75
  ) +
  geom_point(
    data = dilution,
    aes(x = .data$x_position, y = .data$cvaf, color = .data$patient),
    size = 2.1,
    alpha = 0.9
  ) +
  geom_violin(
    data = control_plot_data,
    aes(x = .data$x_position, y = .data$cvaf, fill = .data$control_group),
    width = 0.68,
    scale = "width",
    trim = TRUE,
    alpha = 0.16,
    color = NA
  ) +
  geom_boxplot(
    data = control_plot_data,
    aes(x = .data$x_position, y = .data$cvaf, color = .data$control_group),
    width = 0.25,
    outlier.shape = NA,
    linewidth = 0.55,
    alpha = 0.85
  ) +
  geom_point(
    data = control_plot_data,
    aes(x = .data$x_position, y = .data$cvaf, color = .data$control_group),
    position = position_jitter(width = 0.22, height = 0, seed = 20260802),
    size = 0.85,
    alpha = 0.22
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      patient_colors,
      "Healthy controls" = "#6A3D9A",
      "Unmatched plasma" = "#111111"
    ),
    breaks = c(names(patient_colors), "Healthy controls", "Unmatched plasma"),
    labels = c(
      "M4-28 (n=5,291 mutations)",
      "M4-34 (n=3,439 mutations)",
      "M4-37 (n=3,654 mutations)",
      "M4-38 (n=3,070 mutations)",
      "Healthy controls",
      "Unmatched plasma"
    )
  ) +
  scale_fill_manual(
    values = c("Healthy controls" = "#6A3D9A", "Unmatched plasma" = "#111111"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = -1:7,
    labels = c(
      "10%", "1%", "0.1%", "0.01%", "0.001%", "0.0001%", "0%",
      "Healthy\ncontrols", "Unmatched\nplasma"
    ),
    limits = c(-1.18, 7.35),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(ylim = c(0, max(dilution$cvaf, control_plot_data$cvaf, na.rm = TRUE) * 1.06)) +
  labs(
    title = "Cumulative VAF across controlled dilution series and negative controls",
    subtitle = paste0(
      "Four patient dilution series; healthy controls are platform matched; ",
      "unmatched plasma represents 725 evaluable revision XPlus all-by-all pairings"
    ),
    x = "Tumor fraction in dilution (%) and negative-control groups",
    y = "Cumulative VAF (cVAF)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5, color = "#333333"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.key.width = grid::unit(1.2, "lines")
  )

png_path <- file.path(figure_dir, "Dilution_series_cVAF_vs_healthy_and_unmatched_plasma.png")
pdf_path <- file.path(figure_dir, "Dilution_series_cVAF_vs_healthy_and_unmatched_plasma.pdf")
ggsave(png_path, p, width = 12.5, height = 6.4, units = "in", dpi = 600, bg = "white")
ggsave(pdf_path, p, width = 12.5, height = 6.4, units = "in", device = cairo_pdf)

# Zero-preserving pseudo-log views. A conventional log10 scale cannot display
# true zeros; pseudo_log_trans is linear in a narrow interval around zero and
# log-like above it. Sigma is anchored to the smallest observed positive cVAF
# (approximately 2e-5), so no pseudocount is added to any measurement.
pseudolog_sigma <- 2e-5
pseudolog_breaks <- c(0, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2)
pseudolog_labels <- c(
  "0", "2e-5", "5e-5", "1e-4", "2e-4", "5e-4",
  "1e-3", "2e-3", "5e-3", "1e-2", "2e-2"
)

p_pseudolog <- p +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10, sigma = pseudolog_sigma),
    breaks = pseudolog_breaks,
    labels = pseudolog_labels
  ) +
  labs(
    subtitle = paste0(
      "Four patient dilution series; 725 evaluable revision XPlus all-by-all pairings; ",
      "pseudo-log cVAF axis preserves zero"
    )
  )

pseudolog_png_path <- file.path(
  figure_dir,
  "Dilution_series_cVAF_vs_healthy_and_unmatched_plasma_pseudolog.png"
)
pseudolog_pdf_path <- file.path(
  figure_dir,
  "Dilution_series_cVAF_vs_healthy_and_unmatched_plasma_pseudolog.pdf"
)
ggsave(
  pseudolog_png_path, p_pseudolog,
  width = 12.5, height = 6.4, units = "in", dpi = 600, bg = "white"
)
ggsave(
  pseudolog_pdf_path, p_pseudolog,
  width = 12.5, height = 6.4, units = "in", device = cairo_pdf
)

# Focus the x-axis on the low-tumor-fraction region: 0.1%, 0.01%, 0.001%,
# 0.0001%, the measured 0% controls, and the two external comparator groups.
p_pseudolog_low_lod <- p_pseudolog +
  coord_cartesian(xlim = c(0.78, 7.35), clip = "on") +
  labs(
    title = "Cumulative VAF near the lowest tested dilution levels",
    x = "Low tumor fraction in dilution (%) and negative-control groups"
  )

pseudolog_low_lod_png_path <- file.path(
  figure_dir,
  "Dilution_series_cVAF_vs_healthy_and_unmatched_plasma_pseudolog_low_LOD_zoom.png"
)
pseudolog_low_lod_pdf_path <- file.path(
  figure_dir,
  "Dilution_series_cVAF_vs_healthy_and_unmatched_plasma_pseudolog_low_LOD_zoom.pdf"
)
ggsave(
  pseudolog_low_lod_png_path, p_pseudolog_low_lod,
  width = 11.5, height = 6.4, units = "in", dpi = 600, bg = "white"
)
ggsave(
  pseudolog_low_lod_pdf_path, p_pseudolog_low_lod,
  width = 11.5, height = 6.4, units = "in", device = cairo_pdf
)

# BM-derived cVAF and mutant-sites z-score companion figure. Field mapping is
# inherited from the canonical dilution workflow:
#   z_score_detection_rate_BM = cVAF z-score
#   zscore_BM                 = proportion of sites detected z-score
zscore_labels <- c(
  z_score_detection_rate_BM = "BM-derived cVAF z-score",
  zscore_BM = "BM-derived proportion of sites detected z-score"
)

zscore_dilution <- dilution_raw %>%
  filter(.data$feature %in% names(zscore_labels), !is.na(.data$value)) %>%
  transmute(
    record_type = "Dilution series",
    feature = .data$feature,
    feature_label = unname(zscore_labels[.data$feature]),
    patient = .data$Patient,
    query_bam = basename(.data$Bam),
    mutation_panel = .data$Patient,
    platform = if_else(.data$Patient == "M4-28", "NovaSeq 6000", "NovaSeq XPlus"),
    tumor_fraction_percent = .data$LOD_zero_final,
    x_position = if_else(.data$LOD_zero_final > 0, -log10(.data$LOD_zero_final), 5),
    zscore = .data$value,
    line_group = .data$line_group
  )

zscore_healthy <- healthy_raw %>%
  filter(.data$feature %in% names(zscore_labels), !is.na(.data$value)) %>%
  transmute(
    record_type = "Healthy controls",
    feature = .data$feature,
    feature_label = unname(zscore_labels[.data$feature]),
    patient = .data$Patient,
    query_bam = .data$control_bam,
    mutation_panel = .data$VCF_clean,
    platform = if_else(
      .data$healthy_reference_tier == "XPLUS_CHARM_healthy",
      "NovaSeq XPlus", "NovaSeq 6000"
    ),
    tumor_fraction_percent = NA_real_,
    x_position = 6,
    zscore = .data$value,
    line_group = NA_character_
  )

zscore_unmatched <- all_mrdetect %>%
  filter(
    .data$plotting_type == "Unmatched_plasma",
    .data$revision_new_xplus_bam_Bam %in% TRUE,
    .data$Sample_type_Bam == "Blood_plasma_cfDNA",
    .data$Mut_source == "BM_cells",
    .data$Filter_source == "STR_encode",
    .data$timepoint_info == "Baseline",
    .data$Study_Bam == "IMMAGINE_revision_OICR",
    .data$healthy_reference_requested == "XPLUS_CHARM_healthy",
    .data$healthy_reference_status == "available",
    .data$Patient %in% .env$eligible_panel_patients
  ) %>%
  distinct(.data$BAM, .data$VCF_clean, .keep_all = TRUE) %>%
  mutate(query_sample_code = str_remove(.data$Sample_ID_Bam, "-P$")) %>%
  semi_join(evaluable_plasma_samples, by = c("query_sample_code" = "Sample_Code")) %>%
  dplyr::select(
    "BAM",
    "VCF_clean",
    cVAF_zscore = "detection_rate_zscore_reads_checked_charm",
    sites_zscore = "sites_rate_zscore_charm"
  ) %>%
  pivot_longer(
    cols = c("cVAF_zscore", "sites_zscore"),
    names_to = "zscore_type",
    values_to = "zscore"
  ) %>%
  filter(!is.na(.data$zscore)) %>%
  mutate(
    feature = recode(
      .data$zscore_type,
      cVAF_zscore = "z_score_detection_rate_BM",
      sites_zscore = "zscore_BM"
    )
  ) %>%
  transmute(
    record_type = "Unmatched plasma",
    feature = .data$feature,
    feature_label = unname(zscore_labels[.data$feature]),
    patient = NA_character_,
    query_bam = basename(.data$BAM),
    mutation_panel = .data$VCF_clean,
    platform = "NovaSeq XPlus",
    tumor_fraction_percent = NA_real_,
    x_position = 7,
    zscore = .data$zscore,
    line_group = NA_character_
  )

zscore_unmatched_counts <- zscore_unmatched %>% count(.data$feature)
if (nrow(zscore_unmatched_counts) != 2L || any(zscore_unmatched_counts$n != 725L)) {
  stop("Expected 725 evaluable unmatched values for each BM z-score feature.", call. = FALSE)
}

zscore_source <- bind_rows(zscore_dilution, zscore_healthy, zscore_unmatched) %>%
  mutate(
    record_type = factor(
      .data$record_type,
      levels = c("Dilution series", "Healthy controls", "Unmatched plasma")
    ),
    feature_label = factor(.data$feature_label, levels = unname(zscore_labels))
  )
write_csv(
  zscore_source,
  file.path(
    source_dir,
    "SourceData_Dilution_BM_cVAF_and_sites_zscores_vs_healthy_and_unmatched_plasma.csv"
  )
)

zscore_controls <- zscore_source %>%
  filter(.data$record_type != "Dilution series") %>%
  mutate(control_group = factor(
    as.character(.data$record_type),
    levels = c("Healthy controls", "Unmatched plasma")
  ))
zscore_lines <- zscore_source %>% filter(.data$record_type == "Dilution series")

p_zscores <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "#666666", linewidth = 0.45) +
  geom_line(
    data = zscore_lines,
    aes(x = .data$x_position, y = .data$zscore,
        group = .data$line_group, color = .data$patient),
    linewidth = 0.65, alpha = 0.75
  ) +
  geom_point(
    data = zscore_lines,
    aes(x = .data$x_position, y = .data$zscore, color = .data$patient),
    size = 1.8, alpha = 0.9
  ) +
  geom_violin(
    data = zscore_controls,
    aes(x = .data$x_position, y = .data$zscore, fill = .data$control_group),
    width = 0.7, scale = "width", trim = TRUE, alpha = 0.16, color = NA
  ) +
  geom_boxplot(
    data = zscore_controls,
    aes(x = .data$x_position, y = .data$zscore, color = .data$control_group),
    width = 0.25, outlier.shape = NA, linewidth = 0.5, alpha = 0.85
  ) +
  geom_point(
    data = zscore_controls,
    aes(x = .data$x_position, y = .data$zscore, color = .data$control_group),
    position = position_jitter(width = 0.22, height = 0, seed = 20260802),
    size = 0.65, alpha = 0.18
  ) +
  facet_wrap(vars(.data$feature_label), nrow = 1, scales = "free_y") +
  scale_color_manual(
    name = NULL,
    values = c(
      patient_colors,
      "Healthy controls" = "#6A3D9A",
      "Unmatched plasma" = "#111111"
    ),
    breaks = c(names(patient_colors), "Healthy controls", "Unmatched plasma"),
    labels = c(
      "M4-28 (n=5,291 mutations)",
      "M4-34 (n=3,439 mutations)",
      "M4-37 (n=3,654 mutations)",
      "M4-38 (n=3,070 mutations)",
      "Healthy controls", "Unmatched plasma"
    )
  ) +
  scale_fill_manual(
    values = c("Healthy controls" = "#6A3D9A", "Unmatched plasma" = "#111111"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = -1:7,
    labels = c(
      "10%", "1%", "0.1%", "0.01%", "0.001%", "0.0001%", "0%",
      "Healthy\ncontrols", "Unmatched\nplasma"
    ),
    limits = c(-1.18, 7.35),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "BM-derived MRDetect z-scores across controlled dilution series",
    subtitle = paste0(
      "Platform-matched healthy controls and 725 evaluable revision XPlus ",
      "all-by-all unmatched pairings per feature"
    ),
    x = "Tumor fraction in dilution (%) and negative-control groups",
    y = "MRDetect z-score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5, color = "#333333"),
    strip.background = element_rect(fill = "#F2F2F2", color = "#333333"),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

zscore_png_path <- file.path(figure_dir, "Dilution_series_BM_cVAF_and_sites_zscores_vs_controls.png")
zscore_pdf_path <- file.path(figure_dir, "Dilution_series_BM_cVAF_and_sites_zscores_vs_controls.pdf")
ggsave(zscore_png_path, p_zscores, width = 14.5, height = 6.2, units = "in", dpi = 600, bg = "white")
ggsave(zscore_pdf_path, p_zscores, width = 14.5, height = 6.2, units = "in", device = cairo_pdf)

# Symmetric pseudo-log preserves negative z-scores and zero while expanding the
# control-scale region near zero.
p_zscores_pseudolog <- p_zscores +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(-2, -1, 0, 1, 2, 5, 10, 20, 50, 100, 200, 500),
    labels = c("-2", "-1", "0", "1", "2", "5", "10", "20", "50", "100", "200", "500")
  ) +
  labs(subtitle = paste0(
    "Platform-matched healthy controls and 725 evaluable unmatched pairings per feature; ",
    "symmetric pseudo-log axis preserves negative values and zero"
  ))

zscore_pseudolog_png_path <- file.path(
  figure_dir, "Dilution_series_BM_cVAF_and_sites_zscores_vs_controls_pseudolog.png"
)
zscore_pseudolog_pdf_path <- file.path(
  figure_dir, "Dilution_series_BM_cVAF_and_sites_zscores_vs_controls_pseudolog.pdf"
)
ggsave(
  zscore_pseudolog_png_path, p_zscores_pseudolog,
  width = 14.5, height = 6.2, units = "in", dpi = 600, bg = "white"
)
ggsave(
  zscore_pseudolog_pdf_path, p_zscores_pseudolog,
  width = 14.5, height = 6.2, units = "in", device = cairo_pdf
)

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Saved: ", pseudolog_png_path)
message("Saved: ", pseudolog_pdf_path)
message("Saved: ", pseudolog_low_lod_png_path)
message("Saved: ", pseudolog_low_lod_pdf_path)
message("Saved: ", zscore_png_path)
message("Saved: ", zscore_pdf_path)
message("Saved: ", zscore_pseudolog_png_path)
message("Saved: ", zscore_pseudolog_pdf_path)
message("Source rows: ", nrow(source_data),
        " (dilution=", nrow(dilution),
        ", healthy=", nrow(healthy),
        ", unmatched=", nrow(unmatched), ")")
