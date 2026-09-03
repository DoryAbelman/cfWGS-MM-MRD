#!/usr/bin/env Rscript

# =============================================================================
# 6_14_Build_Corrected_Grouped_CV_Panels.R
#
# Build replacement Figure 3A and Figure 4A components from the fully
# patient-grouped repeated nested-CV summaries. These panels report
# threshold-free discrimination (AUC) with patient-clustered bootstrap
# intervals. They do not overwrite the historical fold-wise
# sensitivity/specificity panels.
#
# Script roadmap
#   1. Read the locked grouped-CV summaries for the headline comparisons.
#   2. Require exactly four models and the locked resampling design.
#   3. Export the exact plotted values and input checksums.
#   4. Render Figure 3A and Figure 4A in wide or compact layout.
#   5. Hash outputs and write a versioned RUN_COMPLETE marker.
#
# Why plotting is separate from fitting
#   Figure construction should be deterministic from a compact summary table.
#   This prevents a cosmetic layout change from refitting a model or changing
#   the underlying validation estimates.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ----------------------------------------------------------------------------
# 1. Resolve locked summaries and a new versioned output directory
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) && length(args) > hit[1]) args[hit[1] + 1] else default
}

RUN_ID <- get_arg("--run-id", format(Sys.Date(), "%Y-%m-%d_v1"))
PANEL_LAYOUT <- get_arg("--layout", "wide")
if (!PANEL_LAYOUT %in% c("wide", "compact")) {
  stop("--layout must be either 'wide' or 'compact'.", call. = FALSE)
}
HEADLINE_RUN <- get_arg(
  "--headline-run",
  file.path("Output_tables_2025", "patient_grouped_repeated_nested_cv",
            "2026-07-31_v1", "primary_performance_summary.csv")
)
BLOOD_SITES_RUN <- get_arg(
  "--blood-sites-run",
  file.path("Output_tables_2025", "patient_grouped_repeated_nested_cv",
            "2026-07-31_blood_sites_v1", "primary_performance_summary.csv")
)
OUT_DIR <- file.path(
  get_arg("--output-root",
          file.path("Output_figures_2025", "patient_grouped_repeated_nested_cv")),
  RUN_ID
)

if (dir.exists(OUT_DIR) || file.exists(OUT_DIR)) {
  stop("Refusing to overwrite existing panel directory: ", OUT_DIR, call. = FALSE)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

inputs <- c(HEADLINE_RUN, BLOOD_SITES_RUN)
if (any(!file.exists(inputs))) {
  stop("Missing grouped-CV summary input(s): ",
       paste(inputs[!file.exists(inputs)], collapse = ", "), call. = FALSE)
}

# ----------------------------------------------------------------------------
# 2. Validate the exact models and resampling design before plotting
# ----------------------------------------------------------------------------
raw <- bind_rows(
  read_csv(HEADLINE_RUN, show_col_types = FALSE),
  read_csv(BLOOD_SITES_RUN, show_col_types = FALSE)
)
required_models <- c("BM_Sites", "BM_cVAF", "Blood_Combined", "Blood_Sites")
if (!setequal(raw$model, required_models)) {
  stop("Grouped-CV summaries do not contain exactly the four expected models.",
       call. = FALSE)
}
if (any(raw$outer_repeats != 50L) || any(raw$inner_repeats != 5L) ||
    any(raw$bootstrap_reps != 2000L)) {
  stop("Grouped-CV run parameters differ from the locked panel specification.",
       call. = FALSE)
}

# ----------------------------------------------------------------------------
# 3. Create canonical panel source data and hash both inputs
# ----------------------------------------------------------------------------
source_data <- raw %>%
  transmute(
    model,
    assay_family = if_else(grepl("^BM_", model), "BM-informed", "Baseline-plasma-informed"),
    display_label = recode(
      model,
      BM_Sites = "Mutant-sites z-score",
      BM_cVAF = "cVAF z-score",
      Blood_Combined = "Mutation + fragmentomics",
      Blood_Sites = "Mutant-sites z-score"
    ),
    auc = auc_estimate,
    auc_ci_lower = auc_cluster_boot_q025,
    auc_ci_upper = auc_cluster_boot_q975,
    n_samples,
    n_patients,
    n_positive_samples,
    n_negative_samples,
    outer_repeats,
    outer_folds,
    inner_repeats,
    inner_folds,
    bootstrap_reps,
    estimand
  ) %>%
  arrange(assay_family, desc(auc))

write_csv(source_data, file.path(OUT_DIR, "corrected_grouped_cv_auc_panel_source_data.csv"))

input_manifest <- tibble(
  role = c("headline_grouped_cv_summary", "blood_sites_grouped_cv_summary"),
  path = inputs,
  size_bytes = as.numeric(file.info(inputs)$size),
  modified_time = format(file.info(inputs)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  md5 = unname(tools::md5sum(inputs))
)
write_csv(input_manifest, file.path(OUT_DIR, "input_manifest.csv"))

# ----------------------------------------------------------------------------
# 4. Render panels from the validated source-data table
# ----------------------------------------------------------------------------
# Layout affects presentation only. Both forms show identical AUCs,
# patient-clustered intervals, and sample/patient denominators.
build_panel <- function(data, assay_subtitle, filename_stem, colours) {
  if (identical(PANEL_LAYOUT, "compact")) {
    data <- data %>%
      mutate(
        axis_label = sprintf("%s\n%d samples; %d patients", display_label, n_samples, n_patients),
        axis_label = factor(axis_label, levels = rev(axis_label)),
        estimate_label = sprintf("%.2f [%.2f-%.2f]", auc, auc_ci_lower, auc_ci_upper)
      )

    plot <- ggplot(data, aes(x = auc, y = axis_label, colour = display_label)) +
      geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.6,
                 colour = "grey65") +
      geom_errorbar(aes(xmin = auc_ci_lower, xmax = auc_ci_upper),
                    orientation = "y", width = 0.12, linewidth = 1.15) +
      geom_point(size = 4.4) +
      geom_text(aes(x = 0.51, label = estimate_label), colour = "black",
                hjust = 0, nudge_y = 0.22, size = 3.6, fontface = "bold") +
      scale_x_continuous(
        limits = c(0.45, 1.01), breaks = c(0.5, 0.75, 1),
        expand = expansion(mult = c(0, 0.01))
      ) +
      scale_colour_manual(values = colours, guide = "none") +
      labs(
        title = "Patient-grouped nested-CV AUC",
        subtitle = assay_subtitle,
        x = "Cross-validated AUC",
        y = NULL,
        caption = "50 repeated 5-fold splits; outer and inner folds grouped by patient."
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10.5, colour = "grey25"),
        plot.caption = element_text(size = 7.5, hjust = 0, colour = "grey30"),
        axis.text.y = element_text(size = 9.5, face = "bold", lineheight = 1.0),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10.5),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.margin = margin(8, 8, 8, 8)
      )

    ggsave(file.path(OUT_DIR, paste0(filename_stem, ".png")), plot,
           width = 5.0, height = 4.2, units = "in", dpi = 400, bg = "white")
    ggsave(file.path(OUT_DIR, paste0(filename_stem, ".pdf")), plot,
           width = 5.0, height = 4.2, units = "in", device = cairo_pdf)
    return(invisible(NULL))
  }

  data <- data %>%
    mutate(
      axis_label = sprintf(
        "%s\nAUC %.2f (95%% CI %.2f-%.2f); %d samples, %d patients",
        display_label,
        auc, auc_ci_lower, auc_ci_upper, n_samples, n_patients
      ),
      axis_label = factor(axis_label, levels = rev(axis_label))
    )

  plot <- ggplot(data, aes(x = auc, y = axis_label, colour = display_label)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.7,
               colour = "grey65") +
    geom_errorbar(aes(xmin = auc_ci_lower, xmax = auc_ci_upper),
                  orientation = "y", width = 0.13, linewidth = 1.1) +
    geom_point(size = 4.2) +
    scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_colour_manual(values = colours, guide = "none") +
    labs(
      title = "Patient-grouped nested-CV AUC",
      subtitle = paste0(
        assay_subtitle,
        " | 50 repeated 5-fold splits | patient-clustered 95% CI"
      ),
      x = "Cross-validated AUC",
      y = NULL,
      caption = paste0(
        "Outer and inner folds grouped by patient.\n",
        "Feature specifications fixed before resampling."
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, colour = "grey25"),
      plot.caption = element_text(size = 9, hjust = 0, colour = "grey30"),
      axis.text.y = element_text(size = 11, face = "bold", lineheight = 1.05),
      axis.title.x = element_text(size = 13),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = margin(12, 16, 12, 12)
    )

  ggsave(file.path(OUT_DIR, paste0(filename_stem, ".png")), plot,
         width = 8.5, height = 4.8, units = "in", dpi = 400, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(filename_stem, ".pdf")), plot,
         width = 8.5, height = 4.8, units = "in", device = cairo_pdf)
}

build_panel(
  source_data %>% filter(assay_family == "BM-informed"),
  "BM-informed models",
  "F3A_patient_grouped_nested_cv_auc",
  c("Mutant-sites z-score" = "#D55E00", "cVAF z-score" = "#0072B2")
)

build_panel(
  source_data %>% filter(assay_family == "Baseline-plasma-informed"),
  "Baseline-plasma-informed models",
  "F4A_patient_grouped_nested_cv_auc",
  c("Mutation + fragmentomics" = "#56B4E9", "Mutant-sites z-score" = "#009E73")
)

# ----------------------------------------------------------------------------
# 5. Hash outputs and mark the versioned panel run complete
# ----------------------------------------------------------------------------
output_paths <- list.files(OUT_DIR, full.names = TRUE)
output_paths <- output_paths[basename(output_paths) != "output_manifest.csv"]
write_csv(
  tibble(
    path = output_paths,
    size_bytes = as.numeric(file.info(output_paths)$size),
    md5 = unname(tools::md5sum(output_paths))
  ),
  file.path(OUT_DIR, "output_manifest.csv")
)
writeLines(
  c(
    "Corrected patient-grouped CV panels completed successfully.",
    paste("Run ID:", RUN_ID),
    paste("Output directory:", OUT_DIR),
    "Historical Figure 3A and Figure 4A files were not overwritten."
  ),
  file.path(OUT_DIR, "RUN_COMPLETE")
)

message("Corrected grouped-CV panels written to: ", OUT_DIR)
