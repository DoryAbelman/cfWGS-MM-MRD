#!/usr/bin/env Rscript
# Reviewer-only QC figure for NovaSeq 6000 versus NovaSeq XPlus fragmentomics.
# This script visualizes the observed healthy-control platform shift and the
# location/scale correction used before applying the frozen NovaSeq 6000 model.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

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
  "Reviewer_platform_harmonization"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(comparison_path, parameter_path, coverage_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing reviewer-figure input(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

controls <- read_csv(comparison_path, show_col_types = FALSE) %>%
  mutate(
    fragmentomics_sequencing_platform = recode(
      .data$Platform,
      Historical_NovaSeq6000 = "NovaSeq 6000",
      Spring2026_Xplus = "NovaSeq XPlus"
    )
  ) %>%
  filter_fragmentomics_shared_healthy_controls() %>%
  select(.data$Sample, .data$fragmentomics_sequencing_platform, .data$FS, .data$Proportion.Short) %>%
  left_join(
    readRDS(coverage_path) %>%
      filter(.data$fragmentomics_sample_role == "healthy_control") %>%
      select(.data$Sample, .data$fragmentomics_sequencing_platform, .data$Mean.Coverage),
    by = c("Sample", "fragmentomics_sequencing_platform")
  )

if (anyNA(controls$Mean.Coverage) ||
    any(count(controls, .data$fragmentomics_sequencing_platform)$n != 19L)) {
  stop("Reviewer figure requires 19 complete identity-matched controls per platform.", call. = FALSE)
}

parameters <- read_csv(parameter_path, show_col_types = FALSE)
required_parameter_rows <- tidyr::expand_grid(
  feature = c("FS", "Proportion.Short", "Mean.Coverage"),
  fragmentomics_sequencing_platform = c("NovaSeq 6000", "NovaSeq XPlus")
)
if (nrow(anti_join(required_parameter_rows, parameters, by = names(required_parameter_rows)))) {
  stop("Harmonization parameter table is incomplete for reviewer figure.", call. = FALSE)
}

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
    .data$Sample, .data$fragmentomics_sequencing_platform, .data$feature,
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
write_csv(plot_data, file.path(out_dir, "reviewer_platform_harmonization_source_data.csv"))
write_csv(summary_table, file.path(out_dir, "reviewer_platform_harmonization_summary.csv"))

palette <- c("NovaSeq 6000" = "#3B6FB6", "NovaSeq XPlus" = "#D55E00")
reviewer_plot <- ggplot(
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
  file.path(out_dir, "Reviewer_Figure_fragmentomics_platform_shift_and_correction.pdf"),
  reviewer_plot, width = 7.2, height = 7.0, units = "in", device = cairo_pdf
)
ggsave(
  file.path(out_dir, "Reviewer_Figure_fragmentomics_platform_shift_and_correction.png"),
  reviewer_plot, width = 7.2, height = 7.0, units = "in", dpi = 400, bg = "white"
)

message("Reviewer platform-harmonization figure and source tables written to: ", out_dir)
