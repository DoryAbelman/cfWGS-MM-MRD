#!/usr/bin/env Rscript

# 1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R
# Build Extended Data Figure 3F-G: fragmentomics platform shift and
# cross-validated harmonization in identity-matched healthy-control cfDNA.
#
# Panel F shows the paired raw NovaSeq 6000 versus NovaSeq XPlus values for the
# three fragmentomics inputs used by the frozen models. Panel G validates the
# location/scale correction without reusing the held-out control: for each
# identity and feature, harmonization parameters are estimated from the other
# 18 matched identities, and the held-out XPlus value is mapped to the NovaSeq
# 6000 scale. This leave-one-control-out display is validation only; production
# scoring continues to use the prespecified all-19-control transform.
#
# Inputs:
#   Results_Fragmentomics/CHARM_Xplus_HC_comparison/
#     sample_level_feature_values.csv
#   Results_Fragmentomics/MM_DARs_chromatin_activation_data.rds
# Outputs:
#   ED3F/ED3G figures, source data, statistics, manifests, and captions in the
#   generated-component and manuscript-object trees.
# Unit of inference:
#   One paired healthy-control identity (n = 19); three features are tested and
#   their P values are adjusted together by Benjamini-Hochberg.
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(helpers_path)) helpers_path <- "helpers.R"
source(helpers_path)
rm(helpers_path)

comparison_path <- file.path(
  "Results_Fragmentomics", "CHARM_Xplus_HC_comparison",
  "sample_level_feature_values.csv"
)
coverage_path <- file.path(
  "Results_Fragmentomics", "MM_DARs_chromatin_activation_data.rds"
)
component_root <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects", "generated",
  "figure_components", "Extended_Data_Figure_3"
)
manuscript_root <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects",
  "02_extended_data_figures", "Extended_Data_Figure_3"
)
panel_f_dir <- file.path(component_root, "panel_F")
panel_g_dir <- file.path(component_root, "panel_G")
manuscript_f_dir <- file.path(manuscript_root, "Extended_Data_Figure_3F")
manuscript_g_dir <- file.path(manuscript_root, "Extended_Data_Figure_3G")
combined_dir <- file.path(
  "Results_Fragmentomics", "Extended_Data_Figure_3F_G"
)

required_paths <- c(comparison_path, coverage_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop(
    "Missing fragmentomics calibration input(s): ",
    paste(missing_paths, collapse = ", "),
    call. = FALSE
  )
}

dir.create(panel_f_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_g_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_f_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_g_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

platform_levels <- c("NovaSeq 6000", "NovaSeq XPlus")
feature_levels <- c(
  "Fragment score",
  "Short-fragment proportion",
  "Mean regulatory coverage"
)
feature_labels <- c(
  FS = "Fragment score",
  Proportion.Short = "Short-fragment proportion",
  Mean.Coverage = "Mean regulatory coverage"
)

coverage <- readRDS(coverage_path) %>%
  filter(.data$fragmentomics_sample_role == "healthy_control") %>%
  filter_fragmentomics_shared_healthy_controls() %>%
  select(
    Sample,
    fragmentomics_sequencing_platform,
    Mean.Coverage
  )

controls <- read_csv(comparison_path, show_col_types = FALSE) %>%
  mutate(
    fragmentomics_sequencing_platform = recode(
      .data$Platform,
      Historical_NovaSeq6000 = "NovaSeq 6000",
      Spring2026_Xplus = "NovaSeq XPlus"
    )
  ) %>%
  filter_fragmentomics_shared_healthy_controls() %>%
  select(
    fragmentomics_healthy_control_id,
    Sample,
    fragmentomics_sequencing_platform,
    FS,
    Proportion.Short
  ) %>%
  left_join(
    coverage,
    by = c("Sample", "fragmentomics_sequencing_platform")
  )

control_audit <- controls %>%
  count(
    .data$fragmentomics_healthy_control_id,
    .data$fragmentomics_sequencing_platform,
    name = "n_libraries"
  )
if (nrow(controls) != 38L ||
    n_distinct(controls$fragmentomics_healthy_control_id) != 19L ||
    anyNA(controls[c("FS", "Proportion.Short", "Mean.Coverage")]) ||
    any(control_audit$n_libraries != 1L) ||
    any(count(control_audit, .data$fragmentomics_healthy_control_id)$n != 2L)) {
  stop(
    "ED3F-G requires 19 complete one-to-one healthy-control identities on both platforms.",
    call. = FALSE
  )
}

wide <- controls %>%
  pivot_longer(
    c(FS, Proportion.Short, Mean.Coverage),
    names_to = "feature",
    values_to = "value"
  ) %>%
  mutate(
    feature_label = factor(
      unname(feature_labels[.data$feature]),
      levels = feature_levels
    )
  ) %>%
  select(
    control_id = fragmentomics_healthy_control_id,
    feature,
    feature_label,
    fragmentomics_sequencing_platform,
    value
  ) %>%
  pivot_wider(
    names_from = fragmentomics_sequencing_platform,
    values_from = value
  )

if (nrow(wide) != 19L * 3L ||
    anyNA(wide[platform_levels]) ||
    any(!is.finite(as.matrix(wide[platform_levels])))) {
  stop("Incomplete paired fragmentomics feature matrix.", call. = FALSE)
}

# Panel F: raw paired platform shift ----------------------------------------
panel_f_source <- wide %>%
  pivot_longer(
    all_of(platform_levels),
    names_to = "platform",
    values_to = "raw_value"
  ) %>%
  mutate(platform = factor(.data$platform, levels = platform_levels))

panel_f_stats <- wide %>%
  group_by(.data$feature, .data$feature_label) %>%
  summarise(
    n_pairs = n(),
    historical_sd = sd(.data$`NovaSeq 6000`),
    median_standardized_shift = median(
      (.data$`NovaSeq XPlus` - .data$`NovaSeq 6000`) / .data$historical_sd
    ),
    paired_wilcoxon_p = wilcox.test(
      .data$`NovaSeq XPlus`, .data$`NovaSeq 6000`,
      paired = TRUE, exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    paired_wilcoxon_fdr = p.adjust(.data$paired_wilcoxon_p, method = "BH"),
    annotation = sprintf(
      "Median Δ = %+.2f SD; FDR = %.2g",
      .data$median_standardized_shift,
      .data$paired_wilcoxon_fdr
    )
  )

# Panel G: leave-one-control-out harmonization ------------------------------
loo_rows <- lapply(seq_len(nrow(wide)), function(i) {
  held_out <- wide[i, , drop = FALSE]
  training <- wide %>%
    filter(
      .data$feature == held_out$feature[[1]],
      .data$control_id != held_out$control_id[[1]]
    )
  if (nrow(training) != 18L) {
    stop("Each leave-one-out transform must be estimated from 18 identities.", call. = FALSE)
  }
  mean_6000 <- mean(training$`NovaSeq 6000`)
  sd_6000 <- sd(training$`NovaSeq 6000`)
  mean_xplus <- mean(training$`NovaSeq XPlus`)
  sd_xplus <- sd(training$`NovaSeq XPlus`)
  if (any(!is.finite(c(mean_6000, sd_6000, mean_xplus, sd_xplus))) ||
      sd_6000 <= 0 || sd_xplus <= 0) {
    stop("Invalid leave-one-out harmonization parameters.", call. = FALSE)
  }
  corrected_xplus <- (
    (held_out$`NovaSeq XPlus`[[1]] - mean_xplus) / sd_xplus
  ) * sd_6000 + mean_6000
  held_out %>%
    transmute(
      .data$control_id,
      .data$feature,
      .data$feature_label,
      historical_value = .data$`NovaSeq 6000`,
      xplus_raw_value = .data$`NovaSeq XPlus`,
      xplus_loo_corrected_value = corrected_xplus,
      training_n = nrow(training),
      training_mean_6000 = mean_6000,
      training_sd_6000 = sd_6000,
      training_mean_xplus = mean_xplus,
      training_sd_xplus = sd_xplus,
      difference_before_sd =
        (.data$`NovaSeq XPlus` - .data$`NovaSeq 6000`) / sd_6000,
      difference_after_sd =
        (corrected_xplus - .data$`NovaSeq 6000`) / sd_6000
    )
})
loo_wide <- bind_rows(loo_rows)

if (nrow(loo_wide) != 19L * 3L ||
    any(loo_wide$training_n != 18L) ||
    any(!is.finite(as.matrix(
      loo_wide[c("difference_before_sd", "difference_after_sd")]
    )))) {
  stop("Incomplete leave-one-control-out fragmentomics audit.", call. = FALSE)
}

panel_g_source <- loo_wide %>%
  select(
    control_id,
    feature,
    feature_label,
    historical_value,
    xplus_raw_value,
    xplus_loo_corrected_value,
    training_n,
    training_mean_6000,
    training_sd_6000,
    training_mean_xplus,
    training_sd_xplus,
    `Before correction` = difference_before_sd,
    `Leave-one-out corrected` = difference_after_sd
  ) %>%
  pivot_longer(
    c(`Before correction`, `Leave-one-out corrected`),
    names_to = "correction_state",
    values_to = "standardized_paired_difference"
  ) %>%
  mutate(
    correction_state = factor(
      .data$correction_state,
      levels = c("Before correction", "Leave-one-out corrected")
    )
  )

panel_g_stats <- loo_wide %>%
  group_by(.data$feature, .data$feature_label) %>%
  summarise(
    n_pairs = n(),
    median_abs_difference_before_sd = median(abs(.data$difference_before_sd)),
    median_abs_difference_after_sd = median(abs(.data$difference_after_sd)),
    abs_difference_wilcoxon_p = wilcox.test(
      abs(.data$difference_before_sd),
      abs(.data$difference_after_sd),
      paired = TRUE, exact = FALSE
    )$p.value,
    corrected_shift_wilcoxon_p = wilcox.test(
      .data$difference_after_sd,
      mu = 0, exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    abs_difference_wilcoxon_fdr = p.adjust(
      .data$abs_difference_wilcoxon_p, method = "BH"
    ),
    corrected_shift_wilcoxon_fdr = p.adjust(
      .data$corrected_shift_wilcoxon_p, method = "BH"
    ),
    annotation = sprintf(
      "Median |\u0394|: %.2f \u2192 %.2f SD\nCorrected shift FDR = %.2g",
      .data$median_abs_difference_before_sd,
      .data$median_abs_difference_after_sd,
      .data$corrected_shift_wilcoxon_fdr
    )
  )

# Nature-style plotting ------------------------------------------------------
platform_colors <- c(
  "NovaSeq 6000" = "#222222",
  "NovaSeq XPlus" = "#FF3030"
)
correction_colors <- c(
  "Before correction" = "#FF3030",
  "Leave-one-out corrected" = "#222222"
)

base_theme <- theme_classic(base_size = 11.5, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text = element_text(
      face = "plain", size = 10.8,
      margin = margin(2.5, 2, 2.5, 2)
    ),
    plot.title = element_text(
      face = "bold", size = 13.5, hjust = 0.5,
      margin = margin(b = 5)
    ),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 9.8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3),
    axis.line = element_line(linewidth = 0.35),
    panel.spacing.x = unit(7, "pt"),
    plot.margin = margin(8, 5, 4, 20)
  )

panel_f <- ggplot(
  panel_f_source,
  aes(
    x = .data$platform,
    y = .data$raw_value,
    group = .data$control_id
  )
) +
  geom_line(colour = "grey68", linewidth = 0.32, alpha = 0.68) +
  geom_point(
    aes(fill = .data$platform),
    shape = 21, size = 1.45, stroke = 0.22,
    colour = "white", alpha = 0.95
  ) +
  geom_text(
    data = panel_f_stats,
    aes(x = 1.5, y = Inf, label = .data$annotation),
    inherit.aes = FALSE,
    vjust = 1.15,
    size = 3.25,
    family = "Arial"
  ) +
  facet_wrap(~feature_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = platform_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.20))) +
  labs(
    title = "Fragmentomics healthy-control platform shift",
    x = NULL,
    y = "Raw feature value"
  ) +
  base_theme

panel_g <- ggplot(
  panel_g_source,
  aes(
    x = .data$correction_state,
    y = .data$standardized_paired_difference,
    group = .data$control_id
  )
) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.34) +
  geom_line(colour = "grey70", linewidth = 0.28, alpha = 0.55) +
  geom_boxplot(
    aes(group = .data$correction_state),
    width = 0.30, outlier.shape = NA,
    fill = "white", colour = "black", linewidth = 0.38
  ) +
  geom_point(
    aes(fill = .data$correction_state),
    position = position_jitter(width = 0.075, height = 0, seed = 20260730),
    shape = 21, size = 1.05, stroke = 0.18,
    alpha = 0.72, colour = "grey20"
  ) +
  geom_text(
    data = panel_g_stats,
    aes(x = 1, y = Inf, label = .data$annotation),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.10,
    size = 3.10,
    family = "Arial"
  ) +
  facet_wrap(~feature_label, nrow = 1) +
  scale_fill_manual(values = correction_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.24))) +
  labs(
    title = "Leave-one-control-out fragmentomics calibration",
    x = NULL,
    y = "XPlus \u2212 6000 (6000 SD)"
  ) +
  base_theme +
  theme(axis.text.x = element_text(size = 9.2, lineheight = 0.90))

combined <- panel_f / panel_g +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(8, 5, 4, 20))

# Compact half-page variants for the complete ED3 assembly. These use the same
# data and statistics but shorter labels so text remains legible when F and G
# are placed side by side at final page width.
compact_feature_labels <- c(
  "Fragment score" = "Fragment score",
  "Short-fragment proportion" = "Short fragments",
  "Mean regulatory coverage" = "Mean coverage"
)
panel_f_stats_compact <- panel_f_stats %>%
  mutate(annotation = sprintf(
    "\u0394 %+.2f SD\nFDR %.2g",
    .data$median_standardized_shift,
    .data$paired_wilcoxon_fdr
  ))
panel_g_stats_compact <- panel_g_stats %>%
  mutate(annotation = sprintf(
    "|\u0394| %.2f\u2192%.2f\nShift FDR %.2g",
    .data$median_abs_difference_before_sd,
    .data$median_abs_difference_after_sd,
    .data$corrected_shift_wilcoxon_fdr
  ))

compact_theme <- base_theme +
  theme(
    strip.text = element_text(size = 10.8, margin = margin(2, 1, 2, 1)),
    plot.title = element_text(size = 14.5),
    axis.title = element_text(size = 11.2),
    axis.text = element_text(size = 10.2),
    panel.spacing.x = unit(4, "pt"),
    plot.margin = margin(5, 3, 3, 10)
  )

panel_f_compact <- ggplot(
  panel_f_source,
  aes(x = .data$platform, y = .data$raw_value, group = .data$control_id)
) +
  geom_line(colour = "grey68", linewidth = 0.30, alpha = 0.68) +
  geom_point(
    aes(fill = .data$platform),
    shape = 21, size = 1.15, stroke = 0.18,
    colour = "white", alpha = 0.95
  ) +
  geom_text(
    data = panel_f_stats_compact,
    aes(x = 1.5, y = Inf, label = .data$annotation),
    inherit.aes = FALSE,
    vjust = 1.08,
    size = 3.05,
    family = "Arial"
  ) +
  facet_wrap(
    ~feature_label, scales = "free_y", nrow = 1,
    labeller = as_labeller(compact_feature_labels)
  ) +
  scale_fill_manual(values = platform_colors, guide = "none") +
  scale_x_discrete(labels = c("NovaSeq 6000" = "6000", "NovaSeq XPlus" = "XPlus")) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.23))) +
  labs(
    title = "Fragmentomics platform shift",
    x = NULL,
    y = "Raw value"
  ) +
  compact_theme

panel_g_compact <- ggplot(
  panel_g_source,
  aes(
    x = .data$correction_state,
    y = .data$standardized_paired_difference,
    group = .data$control_id
  )
) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.32) +
  geom_line(colour = "grey70", linewidth = 0.26, alpha = 0.55) +
  geom_boxplot(
    aes(group = .data$correction_state),
    width = 0.30, outlier.shape = NA,
    fill = "white", colour = "black", linewidth = 0.35
  ) +
  geom_point(
    aes(fill = .data$correction_state),
    position = position_jitter(width = 0.065, height = 0, seed = 20260730),
    shape = 21, size = 0.85, stroke = 0.14,
    alpha = 0.70, colour = "grey20"
  ) +
  geom_text(
    data = panel_g_stats_compact,
    aes(x = 1, y = Inf, label = .data$annotation),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.08,
    size = 2.95,
    family = "Arial"
  ) +
  facet_wrap(
    ~feature_label, nrow = 1,
    labeller = as_labeller(compact_feature_labels)
  ) +
  scale_fill_manual(values = correction_colors, guide = "none") +
  scale_x_discrete(labels = c(
    "Before correction" = "Raw",
    "Leave-one-out corrected" = "LOO"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.25))) +
  labs(
    title = "LOO fragmentomics calibration",
    x = NULL,
    y = "XPlus \u2212 6000 (SD)"
  ) +
  compact_theme

save_plot <- function(plot, path_stub, width, height) {
  ggsave(
    paste0(path_stub, ".png"), plot,
    width = width, height = height, dpi = 600, bg = "white",
    device = ragg::agg_png
  )
  ggsave(
    paste0(path_stub, ".pdf"), plot,
    width = width, height = height, bg = "white", device = cairo_pdf
  )
  ggsave(
    paste0(path_stub, ".tiff"), plot,
    width = width, height = height, dpi = 600, bg = "white",
    device = ragg::agg_tiff, compression = "lzw"
  )
}

panel_f_stub <- file.path(
  panel_f_dir,
  "Extended_Data_Figure_3F_fragmentomics_raw_platform_shift"
)
panel_g_stub <- file.path(
  panel_g_dir,
  "Extended_Data_Figure_3G_fragmentomics_leave_one_out_calibration"
)
manuscript_f_stub <- file.path(
  manuscript_f_dir,
  "Extended_Data_Figure_3F_fragmentomics_raw_platform_shift"
)
manuscript_g_stub <- file.path(
  manuscript_g_dir,
  "Extended_Data_Figure_3G_fragmentomics_leave_one_out_calibration"
)
combined_stub <- file.path(
  combined_dir,
  "Extended_Data_Figure_3F_G_fragmentomics_platform_calibration"
)
manuscript_combined_stub <- file.path(
  manuscript_root,
  "Extended_Data_Figure_3F_G_fragmentomics_platform_calibration"
)

save_plot(panel_f, panel_f_stub, width = 10.6, height = 4.0)
save_plot(panel_g, panel_g_stub, width = 10.6, height = 4.0)
save_plot(panel_f, manuscript_f_stub, width = 10.6, height = 4.0)
save_plot(panel_g, manuscript_g_stub, width = 10.6, height = 4.0)

ggsave(
  file.path(manuscript_f_dir, "ED3F.png"),
  panel_f_compact,
  width = 5.3, height = 2.2, dpi = 600, bg = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(manuscript_g_dir, "ED3G.png"),
  panel_g_compact,
  width = 5.3, height = 2.2, dpi = 600, bg = "white",
  device = ragg::agg_png
)
# Half-page compact rasters used only by the full ED3 assembly. Exporting at
# the intended half-width aspect preserves readable type after placement.
ggsave(
  file.path(panel_f_dir, "Extended_Data_Figure_3F_compact_for_assembly.png"),
  panel_f_compact, width = 5.3, height = 3.1, dpi = 600, bg = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(panel_g_dir, "Extended_Data_Figure_3G_compact_for_assembly.png"),
  panel_g_compact, width = 5.3, height = 3.1, dpi = 600, bg = "white",
  device = ragg::agg_png
)
save_plot(combined, combined_stub, width = 10.6, height = 8.1)
save_plot(combined, manuscript_combined_stub, width = 10.6, height = 8.1)

# Exact 180-mm, 300-dpi production raster for manuscript placement.
combined_png_image <- png::readPNG(paste0(combined_stub, ".png"))
nature_width_mm <- 180
nature_height_mm <- nature_width_mm *
  dim(combined_png_image)[1] / dim(combined_png_image)[2]
ragg::agg_tiff(
  paste0(combined_stub, "_Nature_180mm_300dpi.tiff"),
  width = nature_width_mm,
  height = nature_height_mm,
  units = "mm",
  res = 300,
  compression = "lzw",
  bg = "white"
)
grid::grid.newpage()
grid::grid.raster(combined_png_image, interpolate = TRUE)
grDevices::dev.off()
if (!file.copy(
  paste0(combined_stub, "_Nature_180mm_300dpi.tiff"),
  paste0(manuscript_combined_stub, "_Nature_180mm_300dpi.tiff"),
  overwrite = TRUE
)) {
  stop("Failed to copy the Nature-resolution ED3F-G TIFF.", call. = FALSE)
}
write_csv(panel_f_source, file.path(panel_f_dir, "Extended_Data_Figure_3F_source_data.csv"))
write_csv(panel_f_stats, file.path(panel_f_dir, "Extended_Data_Figure_3F_statistics.csv"))
write_csv(panel_g_source, file.path(panel_g_dir, "Extended_Data_Figure_3G_source_data.csv"))
write_csv(panel_g_stats, file.path(panel_g_dir, "Extended_Data_Figure_3G_statistics.csv"))
write_csv(panel_f_source, file.path(manuscript_f_dir, "Extended_Data_Figure_3F_source_data.csv"))
write_csv(panel_f_stats, file.path(manuscript_f_dir, "Extended_Data_Figure_3F_statistics.csv"))
write_csv(panel_g_source, file.path(manuscript_g_dir, "Extended_Data_Figure_3G_source_data.csv"))
write_csv(panel_g_stats, file.path(manuscript_g_dir, "Extended_Data_Figure_3G_statistics.csv"))

manifest_paths <- c(
  file.path(component_root, "Extended_Data_Figure_3_panel_manifest.csv"),
  file.path(manuscript_root, "Extended_Data_Figure_3_panel_manifest.csv")
)
if (any(!file.exists(manifest_paths))) {
  stop("Run 1_8D before 1_8E so the ED3A-E manifest exists.", call. = FALSE)
}

component_manifest <- read_csv(manifest_paths[[1]], show_col_types = FALSE) %>%
  filter(!.data$panel %in% c("F", "G"))
if (!identical(component_manifest$panel, LETTERS[1:5])) {
  stop("Expected the existing ED3 manifest to contain panels A-E.", call. = FALSE)
}
component_additions <- tibble::tribble(
  ~figure, ~panel, ~status, ~description, ~generator_script, ~primary_component_png, ~source_data, ~notes,
  "Extended Data Figure 3", "F", "generated", "Paired raw fragmentomics healthy-control platform shift", "1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R", "panel_F/Extended_Data_Figure_3F_fragmentomics_raw_platform_shift.png", "panel_F/Extended_Data_Figure_3F_source_data.csv", "Nineteen paired healthy-control identities; three frozen-model fragmentomics inputs.",
  "Extended Data Figure 3", "G", "generated", "Leave-one-control-out fragmentomics platform calibration", "1_8E_Build_ED3FG_Fragmentomics_Platform_Calibration.R", "panel_G/Extended_Data_Figure_3G_fragmentomics_leave_one_out_calibration.png", "panel_G/Extended_Data_Figure_3G_source_data.csv", "Each held-out identity is corrected using parameters estimated from the other 18 identities."
)
write_csv(
  bind_rows(component_manifest, component_additions),
  manifest_paths[[1]],
  na = ""
)

manuscript_manifest <- read_csv(manifest_paths[[2]], show_col_types = FALSE) %>%
  filter(!.data$panel %in% c("F", "G"))
if (!identical(manuscript_manifest$panel, LETTERS[1:5])) {
  stop("Expected the Extended Data Figure 3 manifest to contain panels A-E.", call. = FALSE)
}
manuscript_additions <- component_additions %>%
  mutate(
    primary_component_png = c(
      "Extended_Data_Figure_3F/Extended_Data_Figure_3F_fragmentomics_raw_platform_shift.png",
      "Extended_Data_Figure_3G/Extended_Data_Figure_3G_fragmentomics_leave_one_out_calibration.png"
    ),
    source_data = c(
      "Extended_Data_Figure_3F/Extended_Data_Figure_3F_source_data.csv",
      "Extended_Data_Figure_3G/Extended_Data_Figure_3G_source_data.csv"
    )
  )
write_csv(
  bind_rows(manuscript_manifest, manuscript_additions),
  manifest_paths[[2]],
  na = ""
)

caption_lines <- c(
  "Extended Data Fig. 3f,g | Cross-platform harmonization of fragmentomics features in healthy-control cfDNA.",
  "f, Paired raw fragmentomics values from 19 identity-matched healthy controls sequenced on NovaSeq 6000 and NovaSeq XPlus. The three displayed features are inputs to the frozen fragmentomics-containing models. Lines connect the same control identity across platforms. Standardized median shifts use the NovaSeq 6000 standard deviation; FDR values are Benjamini-Hochberg-adjusted two-sided paired Wilcoxon P values across the three features.",
  "g, Leave-one-control-out validation of the XPlus-to-NovaSeq 6000 location/scale correction. For each held-out identity and feature, transformation parameters were estimated using the other 18 matched controls. Points show the held-out XPlus-minus-NovaSeq 6000 paired difference in units of the training-set NovaSeq 6000 standard deviation, before and after correction; lines connect the same identity, boxes show medians and interquartile ranges, and whiskers extend to 1.5 times the interquartile range. Corrected-shift FDR values test the corrected paired differences against zero across the three features. Production scoring uses the prespecified transform estimated from all 19 matched controls; leave-one-out estimates are used only for this validation."
)
writeLines(caption_lines, file.path(combined_dir, "Extended_Data_Figure_3F_G_caption.txt"))
writeLines(caption_lines, file.path(manuscript_root, "Extended_Data_Figure_3F_G_caption.txt"))

message("ED3F component: ", paste0(manuscript_f_stub, ".png"))
message("ED3G component: ", paste0(manuscript_g_stub, ".png"))
message("Combined ED3F-G: ", paste0(manuscript_combined_stub, ".png"))
