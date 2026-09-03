#!/usr/bin/env Rscript
# =============================================================================
# 1_7E_Harmonize_Fragmentomics_Sequencing_Platforms.R
#
# Purpose
#   Visualize the NovaSeq 6000 versus NovaSeq XPlus healthy-control shift in
#   three fragmentomics features and the prespecified location/scale mapping
#   applied before frozen NovaSeq 6000 model scoring.
#
# Inputs
#   - Identity-matched healthy-control FS and short-fragment values.
#   - Healthy-control regulatory coverage values.
#   - Versioned platform mean/SD reference parameters.
#
# Analysis unit and method
#   One row per matched healthy-control identity, platform, and feature. Sample
#   names are platform-specific library labels; `fragmentomics_healthy_control_id`
#   is the cross-platform identity key. The XPlus value is standardized on its
#   platform and mapped to the historical NovaSeq 6000 mean/SD. This is
#   calibration QC; it does not retrain a model. Because these same 19 controls
#   estimate and display the mapping, the after-correction panel is an in-sample
#   calibration diagnostic rather than independent validation.
#
# Outputs
#   Source data, platform summaries, and before/after QC figures under
#   Results_Fragmentomics/Healthy_control_platform_harmonization/.
#
# Script roadmap
#   1. Load and validate the matched 19-control input set.
#   2. Join complete reference parameters for all three features/platforms.
#   3. Apply the fixed location/scale transformation.
#   4. Summarize and export auditable source data.
#   5. Plot values before and after harmonization.
#
# Failure conditions
#   Stop on missing inputs, incomplete regulatory coverage, a control count
#   other than 19 per platform, or incomplete calibration parameters.
#
# Usage
#   Rscript Scripts_2025/Final_Scripts/1_7E_Harmonize_Fragmentomics_Sequencing_Platforms.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# 1. Resolve helpers, inputs, and output directory
# -----------------------------------------------------------------------------

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) .helpers_path <- "helpers.R"
source(.helpers_path)
rm(.helpers_path)

comparison_path <- file.path(
  "Results_Fragmentomics", "CHARM_Xplus_HC_comparison",
  "sample_level_feature_values.csv"
)
parameter_path <- file.path(
  "Results_Fragmentomics",
  "fragmentomics_platform_harmonization_reference_parameters.csv"
)
coverage_path <- file.path(
  "Results_Fragmentomics",
  "MM_DARs_chromatin_activation_data.rds"
)
out_dir <- file.path(
  "Results_Fragmentomics",
  "Healthy_control_platform_harmonization"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(comparison_path, parameter_path, coverage_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing platform-harmonization input(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

# -----------------------------------------------------------------------------
# 2. Construct and validate the identity-matched healthy-control analysis frame
# -----------------------------------------------------------------------------

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
    .data$fragmentomics_healthy_control_id,
    .data$Sample,
    .data$fragmentomics_sequencing_platform,
    .data$FS,
    .data$Proportion.Short
  ) %>%
  left_join(
    readRDS(coverage_path) %>%
      filter(.data$fragmentomics_sample_role == "healthy_control") %>%
      select(.data$Sample, .data$fragmentomics_sequencing_platform, .data$Mean.Coverage),
    by = c("Sample", "fragmentomics_sequencing_platform"),
    relationship = "one-to-one"
  )

control_key_counts <- controls %>%
  count(
    .data$fragmentomics_healthy_control_id,
    .data$fragmentomics_sequencing_platform,
    name = "n"
  )
if (anyNA(controls[c("FS", "Proportion.Short", "Mean.Coverage")]) ||
    any(count(controls, .data$fragmentomics_sequencing_platform)$n != 19L) ||
    any(control_key_counts$n != 1L)) {
  stop("Platform harmonization requires 19 complete identity-matched controls per platform.", call. = FALSE)
}

parameters <- read_csv(parameter_path, show_col_types = FALSE)
required_parameter_rows <- tidyr::expand_grid(
  feature = c("FS", "Proportion.Short", "Mean.Coverage"),
  fragmentomics_sequencing_platform = c("NovaSeq 6000", "NovaSeq XPlus")
)
parameter_key_counts <- parameters %>%
  semi_join(required_parameter_rows, by = names(required_parameter_rows)) %>%
  count(.data$feature, .data$fragmentomics_sequencing_platform, name = "n_rows")
invalid_parameters <- parameters %>%
  semi_join(required_parameter_rows, by = names(required_parameter_rows)) %>%
  filter(
    !is.finite(.data$mean) |
      !is.finite(.data$sd) | .data$sd <= 0 |
      !is.finite(.data$n) | .data$n != 19
  )
if (nrow(anti_join(required_parameter_rows, parameters, by = names(required_parameter_rows))) ||
    any(parameter_key_counts$n_rows != 1L) ||
    nrow(invalid_parameters)) {
  stop(
    "Harmonization parameter table must contain one complete, finite 19-control row ",
    "with positive SD for every feature/platform combination.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 3. Apply the fixed platform location/scale harmonization
# -----------------------------------------------------------------------------

long <- controls %>%
  pivot_longer(
    c(.data$FS, .data$Proportion.Short, .data$Mean.Coverage),
    names_to = "feature",
    values_to = "value_before"
  ) %>%
  left_join(
    parameters %>% select(.data$feature, .data$fragmentomics_sequencing_platform, platform_mean = .data$mean, platform_sd = .data$sd),
    by = c("feature", "fragmentomics_sequencing_platform")
  ) %>%
  left_join(
    parameters %>%
      filter(.data$fragmentomics_sequencing_platform == "NovaSeq 6000") %>%
      select(.data$feature, historical_mean = .data$mean, historical_sd = .data$sd),
    by = "feature"
  ) %>%
  mutate(
    value_after = if_else(
      .data$fragmentomics_sequencing_platform == "NovaSeq XPlus",
      ((.data$value_before - .data$platform_mean) / .data$platform_sd) * .data$historical_sd + .data$historical_mean,
      .data$value_before
    ),
    feature = recode(
      .data$feature,
      FS = "Fragment score",
      Proportion.Short = "Proportion short",
      Mean.Coverage = "MM regulatory coverage"
    ),
    fragmentomics_sequencing_platform = factor(
      .data$fragmentomics_sequencing_platform,
      levels = c("NovaSeq 6000", "NovaSeq XPlus")
    )
  )

plot_data <- long %>%
  select(
    .data$fragmentomics_healthy_control_id, .data$Sample,
    .data$fragmentomics_sequencing_platform, .data$feature,
    `Before correction` = .data$value_before,
    `After correction` = .data$value_after
  ) %>%
  pivot_longer(
    c(`Before correction`, `After correction`),
    names_to = "correction_state",
    values_to = "value"
  ) %>%
  mutate(
    correction_state = factor(
      .data$correction_state,
      levels = c("Before correction", "After correction")
    )
  )

summary_table <- plot_data %>%
  group_by(.data$feature, .data$correction_state, .data$fragmentomics_sequencing_platform) %>%
  summarise(
    n = sum(is.finite(.data$value)),
    mean = mean(.data$value, na.rm = TRUE),
    sd = sd(.data$value, na.rm = TRUE),
    median = median(.data$value, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(plot_data, file.path(out_dir, "fragmentomics_platform_harmonization_source_data.csv"))
write_csv(summary_table, file.path(out_dir, "fragmentomics_platform_harmonization_summary.csv"))

palette <- c("NovaSeq 6000" = "#3B6FB6", "NovaSeq XPlus" = "#D55E00")
harmonization_plot <- ggplot(
  plot_data,
  aes(x = .data$fragmentomics_sequencing_platform, y = .data$value,
      fill = .data$fragmentomics_sequencing_platform)
) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.62, linewidth = 0.35) +
  geom_jitter(width = 0.10, height = 0, size = 1.65, alpha = 0.9, shape = 21, stroke = 0.25) +
  facet_grid(.data$feature ~ .data$correction_state, scales = "free_y") +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(
    title = "Fragmentomics platform shift and correction",
    subtitle = "Matched healthy controls (n = 19/platform): XPlus mapped to the NovaSeq 6000 scale",
    x = NULL,
    y = "Feature value",
    caption = "Nineteen identity-matched healthy controls per platform; boxes show median and interquartile range.",
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "Fragmentomics_healthy_control_platform_harmonization.pdf"),
  harmonization_plot, width = 7.2, height = 7.0, units = "in", device = cairo_pdf
)
ggsave(
  file.path(out_dir, "Fragmentomics_healthy_control_platform_harmonization.png"),
  harmonization_plot, width = 7.2, height = 7.0, units = "in", dpi = 400, bg = "white"
)

message("Platform-harmonization figure and source tables written to: ", out_dir)
