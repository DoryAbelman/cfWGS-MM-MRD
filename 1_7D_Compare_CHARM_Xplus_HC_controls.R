#!/usr/bin/env Rscript
# =============================================================================
# 1_7D_Compare_CHARM_Xplus_HC_controls.R
#
# Purpose:
#   Compare the historical CHARM/TGL49 healthy-control fragmentomics reference
#   set against the Spring 2026 CHARM_Xplus_HC healthy controls. The matched
#   scalar-feature table written here is used by 1_7E and 1_8E to display and
#   validate the platform harmonization used for the Spring 2026 cohort.
#
# Rationale:
#   The Spring 2026 test-cohort samples were sequenced on NovaSeq X Plus,
#   whereas the historical healthy controls used by the main pipeline were run
#   on NovaSeq 6000-era data. This script quantifies whether the new X Plus
#   healthy controls differ materially from the historical controls before any
#   future decision to use them as the reference distribution for new cases.
#
# Outputs:
#   Results_Fragmentomics/CHARM_Xplus_HC_comparison/
#     - insert_size_group_comparison.csv
#     - fragment_score_group_comparison.csv
#     - nucleosome_site_metric_group_comparison.csv
#     - fragmentomics_domain_feature_shift_tests.csv
#     - fragmentomics_domain_pca_scores.csv
#     - fragmentomics_domain_pca_variance.csv
#     - fragmentomics_reference_likeness_distances.csv
#     - fragmentomics_reference_likeness_summary.csv
#     - fragmentomics_input_domain_availability.csv
#     - sample_level_feature_values.csv
#     - xplus_hc_comparison_boxplots.pdf
#     - xplus_hc_nucleosome_site_shift.pdf
#     - xplus_hc_domain_pca_facets.pdf
#     - xplus_hc_reference_likeness.pdf
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_7D_Compare_CHARM_Xplus_HC_controls.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) .helpers_path <- "helpers.R"
source(.helpers_path)
rm(.helpers_path)

out_dir <- file.path("Results_Fragmentomics", "CHARM_Xplus_HC_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

historical_dir <- file.path("Fragmentomics_data", "Normals")
xplus_dir <- file.path("Data_Spring_2026_Revisions", "Fragmentomics_Pipeeline_Suite_all_outputs")
# Historical inputs are the existing CHARM/TGL49 healthy-control reference files.
# X Plus inputs are Spring 2026 support-analysis outputs. The script compares
# distributions only; it does not add X Plus controls to the production
# reference set used by 1_7A/1_7B.

historical_insert_path <- file.path(historical_dir, "2023-12-04_TGL49_HBC_insert_size_summary.tsv")
historical_insert_counts_path <- file.path(historical_dir, "2023-12-04_TGL49_HBC_insert_size_counts.tsv")
historical_fs_path <- file.path(historical_dir, "2023-12-04_TGL49_HBC_fragment_scores.tsv")
historical_breakpoint_path <- file.path(historical_dir, "2023-12-04_TGL49_HBC_breakpoint_frequencies.tsv")
historical_endmotif_3prime_path <- file.path(
  historical_dir,
  "2024-08-14_TGL49_HBCs_3prime_endmotif_frequencies.tsv"
)
historical_endmotif_5prime_path <- file.path(
  historical_dir,
  "2024-08-14_TGL49_HBCs_5prime_endmotif_frequencies.tsv"
)
historical_per5mb_path <- file.path(historical_dir, "2023-12-04_TGL49_HBC_per5Mb_fragment_ratios.tsv")
historical_nuc_path <- file.path(
  historical_dir,
  "2024-09-07_cfWGS_fragmentomics_Healthy_control_MM_signature_nucleosome_accessibility_distances.tsv"
)

xplus_insert_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_insert_size_summary.tsv"
)
xplus_insert_counts_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_insert_size_counts.tsv"
)
xplus_fs_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_fragment_scores.tsv"
)
xplus_breakpoint_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_breakpoint_frequencies.tsv"
)
xplus_dinucleotide_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_dinucleotide_frequencies.tsv"
)
xplus_endmotif_5prime_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_5prime_endmotif_frequencies.tsv"
)
xplus_endmotif_3prime_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_3prime_endmotif_frequencies.tsv"
)
xplus_per5mb_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_per5Mb_fragment_ratios.tsv"
)
xplus_nuc_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_nucleosome_accessibility_distances.tsv"
)
xplus_nucleosome_peak_path <- file.path(
  xplus_dir,
  "2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_nucleosome_peak_distances.tsv"
)

required_paths <- c(
  historical_insert_path,
  historical_insert_counts_path,
  historical_fs_path,
  historical_breakpoint_path,
  historical_endmotif_3prime_path,
  historical_endmotif_5prime_path,
  historical_per5mb_path,
  historical_nuc_path,
  xplus_insert_path,
  xplus_insert_counts_path,
  xplus_fs_path,
  xplus_breakpoint_path,
  xplus_endmotif_5prime_path,
  xplus_endmotif_3prime_path,
  xplus_per5mb_path,
  xplus_nuc_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "Missing required CHARM healthy-control comparison input(s):\n",
    paste(missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

domain_availability <- tibble::tribble(
  # Record which fragmentomics domains can be compared directly. Some Spring
  # 2026 outputs have no historical counterpart, so they are inventoried but not
  # included in formal historical-vs-X Plus statistics.
  ~Domain, ~Historical_input, ~Xplus_input, ~Comparison_status,
  "insert_summary", historical_insert_path, xplus_insert_path, "compared",
  "fragment_score", historical_fs_path, xplus_fs_path, "compared",
  "insert_profile", historical_insert_counts_path, xplus_insert_counts_path, "compared",
  "breakpoint_frequency", historical_breakpoint_path, xplus_breakpoint_path, "compared",
  "per5mb_fragment_ratio", historical_per5mb_path, xplus_per5mb_path, "compared",
  "nucleosome_accessibility", historical_nuc_path, xplus_nuc_path, "compared",
  "end_motif_5prime", historical_endmotif_5prime_path, xplus_endmotif_5prime_path,
  "compared_strand_specific",
  "end_motif_3prime", historical_endmotif_3prime_path, xplus_endmotif_3prime_path,
  "compared_strand_specific",
  "dinucleotide_frequency", NA_character_, xplus_dinucleotide_path, "xplus_only_not_compared",
  "nucleosome_peak_distance", NA_character_, xplus_nucleosome_peak_path, "xplus_only_not_compared"
) %>%
  mutate(
    Historical_exists = if_else(is.na(Historical_input), FALSE, file.exists(Historical_input)),
    Xplus_exists = if_else(is.na(Xplus_input), FALSE, file.exists(Xplus_input))
  )

readr::write_csv(
  domain_availability,
  file.path(out_dir, "fragmentomics_input_domain_availability.csv")
)

theme_fragmentomics_qc <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      plot.title = element_text(face = "bold", hjust = 0),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.position = "top"
    )
}

classify_material <- function(sample) {
  # Material class is retained for descriptive provenance. Platform contrasts
  # use the same 19 control identities on both instruments, so material is
  # matched by identity rather than used as an asymmetric exclusion rule.
  case_when(
    str_detect(sample, "(^|[-_])Cf([-_]|$)|_Cf_|-Cf-") ~ "cfDNA",
    str_detect(sample, "(^|[-_])Pb([-_]|$)|_Pb_|-Pb-") ~ "PB_cellular_or_buffy",
    str_detect(sample, "(^|[-_])Ct([-_]|$)|_Ct_|-Ct-") ~ "cellular_or_tissue",
    TRUE ~ "unknown"
  )
}

filter_included_materials <- function(data, context) {
  # Apply the shared-identity gate after each reader standardizes sample names.
  required <- c("Sample", "Platform")
  missing_required <- setdiff(required, names(data))
  if (length(missing_required) > 0) {
    stop(
      "Cannot apply shared healthy-control identity filter to ",
      context,
      "; missing column(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  filtered <- filter_fragmentomics_shared_healthy_controls(data)
  message(context, ": retained 19 matched healthy-control identities per platform.")
  filtered
}

filter_matrix_object_to_included_materials <- function(obj, context) {
  # Matrix readers return both a numeric matrix and sample metadata. This helper
  # keeps those synchronized after matched-identity filtering so PCA/distance analyses
  # do not compare mismatched rows.
  mat <- obj$mat[[1]]
  sample_info <- obj[["sample_info"]]
  while (!is.data.frame(sample_info) && is.list(sample_info) && length(sample_info) == 1L) {
    sample_info <- sample_info[[1]]
  }
  if (is.data.frame(sample_info) && "sample_info" %in% names(sample_info)) {
    sample_info <- bind_rows(sample_info$sample_info)
  }
  if (is.data.frame(sample_info) && !"Sample" %in% names(sample_info) && nrow(sample_info) == nrow(mat)) {
    sample_info <- sample_info %>%
      mutate(Sample = rownames(mat), .before = 1)
  }
  if (!is.data.frame(sample_info) || !"Sample" %in% names(sample_info)) {
    stop(
      "Matrix object is missing sample metadata for ",
      context,
      ". Available metadata columns: ",
      paste(names(sample_info), collapse = ", "),
      call. = FALSE
    )
  }

  sample_info <- filter_included_materials(sample_info, context)
  keep_samples <- intersect(sample_info$Sample, rownames(mat))
  sample_info <- sample_info %>%
    filter(Sample %in% keep_samples)
  mat <- mat[keep_samples, , drop = FALSE]

  tibble(
    mat = list(mat),
    sample_info = list(sample_info)
  )
}

read_insert <- function(path, platform) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      Sample = as.character(Sample),
      Platform = platform,
      Material = classify_material(Sample),
      Proportion.Short = as.numeric(Proportion.Short),
      Proportion.Long = as.numeric(Proportion.Long),
      Ratio.Short.Normal = as.numeric(Ratio.Short.Normal),
      MEDIAN_INSERT_SIZE = as.numeric(MEDIAN_INSERT_SIZE),
      MEAN_INSERT_SIZE = as.numeric(MEAN_INSERT_SIZE),
      READ_PAIRS = as.numeric(READ_PAIRS)
    )
}

read_fs <- function(path, platform) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      Sample = as.character(Sample),
      Platform = platform,
      Material = classify_material(Sample),
      FS = as.numeric(FS)
    )
}

safe_wilcox <- function(x, y) {
  # Use a non-parametric Wilcoxon test because these are small healthy-control
  # groups and feature distributions are not assumed Gaussian. These legacy
  # comparisons are unpaired even though the retained controls are matched by
  # identity; interpret them as exploratory platform summaries. Return NA instead
  # of forcing a p-value when either group has fewer than two finite values.
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
}

impute_feature_median_matrix <- function(mat) {
  # PCA and distance metrics require a complete numeric matrix. Missing or
  # non-finite feature values are filled with the within-feature median so an
  # individual sample is not removed across an entire domain; all-missing
  # features fall back to zero and later zero-variance columns are removed.
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(mat)
  out <- mat
  rn <- rownames(out)
  cn <- colnames(out)
  for (j in seq_len(ncol(out))) {
    col_vals <- suppressWarnings(as.numeric(out[, j]))
    fill_value <- median(col_vals[is.finite(col_vals)], na.rm = TRUE)
    if (!is.finite(fill_value)) fill_value <- 0
    col_vals[!is.finite(col_vals)] <- fill_value
    out[, j] <- col_vals
  }
  mode(out) <- "numeric"
  rownames(out) <- rn
  colnames(out) <- cn
  out
}

clean_matrix_sample_name <- function(x) {
  x %>%
    str_remove("_5prime$") %>%
    str_remove("_3prime$")
}

standardize_feature_matrix <- function(mat) {
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 1) return(mat)
  mat <- impute_feature_median_matrix(mat)
  keep <- apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) == 0) return(mat)
  scaled <- scale(mat)
  scaled[!is.finite(scaled)] <- 0
  scaled
}

read_insert_profile_matrix <- function(path, platform) {
  # Convert insert-size counts to per-sample proportions over 20-400 bp. This
  # focuses the profile on the cfDNA fragment-size range used by the project and
  # avoids total-depth differences driving the comparison.
  count_df <- readr::read_tsv(path, show_col_types = FALSE)
  if (!"insert_size" %in% names(count_df)) {
    stop("Insert-count matrix is missing insert_size: ", path, call. = FALSE)
  }
  profile_wide <- count_df %>%
    pivot_longer(cols = -insert_size, names_to = "Sample", values_to = "count") %>%
    mutate(
      Sample = clean_matrix_sample_name(Sample),
      insert_size = suppressWarnings(as.numeric(insert_size)),
      count = suppressWarnings(as.numeric(count))
    ) %>%
    filter(is.finite(insert_size), is.finite(count), insert_size >= 20, insert_size <= 400) %>%
    group_by(Sample) %>%
    mutate(total_count = sum(count, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(prop = if_else(total_count > 0, count / total_count, NA_real_)) %>%
    select(Sample, insert_size, prop) %>%
    mutate(feature_id = paste0("insert_", insert_size)) %>%
    select(Sample, feature_id, prop) %>%
    pivot_wider(names_from = feature_id, values_from = prop, values_fill = 0)

  mat <- profile_wide %>% select(-Sample) %>% as.matrix()
  mode(mat) <- "numeric"
  rownames(mat) <- profile_wide$Sample
  sample_names <- rownames(mat)
  tibble(
    mat = list(mat),
    sample_info = list(tibble(
      Sample = sample_names,
      Platform = platform,
      Material = classify_material(sample_names)
    ))
  )
}

read_feature_wide_matrix <- function(path, platform, id_cols, feature_builder) {
  # Feature-wide pipeline outputs have one row per genomic/motif/bin feature and
  # one column per sample. The feature_builder creates stable IDs so historical
  # and X Plus matrices can be intersected by exact feature identity.
  raw_df <- readr::read_tsv(path, show_col_types = FALSE)
  missing_id <- setdiff(id_cols, names(raw_df))
  if (length(missing_id) > 0) {
    stop("Missing required feature ID columns in ", path, ": ", paste(missing_id, collapse = ", "), call. = FALSE)
  }
  feature_df <- raw_df %>%
    mutate(feature_id = feature_builder(.)) %>%
    select(-any_of(id_cols)) %>%
    relocate(feature_id)
  sample_cols <- setdiff(names(feature_df), "feature_id")
  mat <- feature_df %>%
    select(all_of(sample_cols)) %>%
    as.matrix()
  rownames(mat) <- feature_df$feature_id
  mat <- t(mat)
  sample_names <- clean_matrix_sample_name(rownames(mat))
  mode(mat) <- "numeric"
  rownames(mat) <- sample_names
  tibble(
    mat = list(impute_feature_median_matrix(mat)),
    sample_info = list(tibble(
      Sample = sample_names,
      Platform = platform,
      Material = classify_material(sample_names)
    ))
  )
}

read_long_feature_matrix <- function(path, platform, sample_col, feature_cols, value_col) {
  # Long-format outputs are collapsed to one value per Sample x feature_id before
  # matrix conversion. Duplicate rows are averaged to handle any repeated
  # nucleotide/position records without dropping the whole feature.
  raw_df <- readr::read_tsv(path, show_col_types = FALSE)
  required <- c(sample_col, feature_cols, value_col)
  missing_required <- setdiff(required, names(raw_df))
  if (length(missing_required) > 0) {
    stop("Missing required long-format columns in ", path, ": ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  wide_df <- raw_df %>%
    transmute(
      Sample = clean_matrix_sample_name(as.character(.data[[sample_col]])),
      feature_id = do.call(paste, c(across(all_of(feature_cols)), sep = "|")),
      value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    group_by(Sample, feature_id) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = feature_id, values_from = value)
  mat <- wide_df %>% select(-Sample) %>% as.matrix()
  mode(mat) <- "numeric"
  rownames(mat) <- wide_df$Sample
  sample_names <- rownames(mat)
  tibble(
    mat = list(impute_feature_median_matrix(mat)),
    sample_info = list(tibble(
      Sample = sample_names,
      Platform = platform,
      Material = classify_material(sample_names)
    ))
  )
}

combine_domain_matrices <- function(historical_obj, xplus_obj, domain) {
  # Compare only common features between historical and X Plus outputs. This
  # prevents platform-specific file schema differences from creating artificial
  # missingness or giving one group extra dimensions in PCA/distance summaries.
  h_mat <- historical_obj$mat[[1]]
  x_mat <- xplus_obj$mat[[1]]
  common_features <- intersect(colnames(h_mat), colnames(x_mat))
  if (length(common_features) < 2) {
    warning("Skipping domain with fewer than two common features: ", domain, call. = FALSE)
    return(NULL)
  }
  mat <- rbind(
    h_mat[, common_features, drop = FALSE],
    x_mat[, common_features, drop = FALSE]
  )

  h_info <- historical_obj$sample_info
  x_info <- xplus_obj$sample_info
  if (!is.data.frame(h_info) && length(h_info) == 1L && is.data.frame(h_info[[1]])) h_info <- h_info[[1]]
  if (!is.data.frame(x_info) && length(x_info) == 1L && is.data.frame(x_info[[1]])) x_info <- x_info[[1]]
  if (is.data.frame(h_info) && "sample_info" %in% names(h_info)) h_info <- bind_rows(h_info$sample_info)
  if (is.data.frame(x_info) && "sample_info" %in% names(x_info)) x_info <- bind_rows(x_info$sample_info)
  if (!"Sample" %in% names(h_info)) {
    h_info <- tibble(
      Sample = rownames(h_mat),
      Platform = "Historical_NovaSeq6000",
      Material = classify_material(rownames(h_mat))
    )
  }
  if (!"Sample" %in% names(x_info)) {
    x_info <- tibble(
      Sample = rownames(x_mat),
      Platform = "Spring2026_Xplus",
      Material = classify_material(rownames(x_mat))
    )
  }

  sample_info <- bind_rows(h_info, x_info) %>%
    mutate(Domain = domain)
  list(
    domain = domain,
    mat = impute_feature_median_matrix(mat),
    sample_info = sample_info,
    n_features = length(common_features)
  )
}

get_domain_sample_info <- function(domain_obj) {
  sample_info <- domain_obj$sample_info
  if (!is.data.frame(sample_info) && length(sample_info) == 1L && is.data.frame(sample_info[[1]])) {
    sample_info <- sample_info[[1]]
  }
  if (is.data.frame(sample_info) && "sample_info" %in% names(sample_info)) {
    sample_info <- bind_rows(sample_info$sample_info)
  }
  if (!"Sample" %in% names(sample_info)) {
    stop(
      "Domain sample metadata is missing Sample for domain: ",
      domain_obj$domain,
      ". Available columns: ",
      paste(names(sample_info), collapse = ", "),
      call. = FALSE
    )
  }
  sample_info
}

compute_domain_pca <- function(domain_obj) {
  # PCA is descriptive QC, not a classifier. Features are already standardized by
  # standardize_feature_matrix(), so prcomp is run without additional centering
  # or scaling to avoid double-transforming the matrix.
  mat <- standardize_feature_matrix(domain_obj$mat)
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
  sample_info <- get_domain_sample_info(domain_obj)
  pca <- stats::prcomp(mat, center = FALSE, scale. = FALSE)
  variance <- (pca$sdev^2) / sum(pca$sdev^2)
  scores <- as_tibble(pca$x[, seq_len(min(5, ncol(pca$x))), drop = FALSE], rownames = "Sample") %>%
    left_join(sample_info, by = "Sample") %>%
    mutate(
      Domain = domain_obj$domain,
      n_features = domain_obj$n_features
    )
  variance_df <- tibble(
    Domain = domain_obj$domain,
    PC = paste0("PC", seq_along(variance)),
    variance_explained = variance,
    n_features = domain_obj$n_features
  )
  list(scores = scores, variance = variance_df, pca = pca, scaled_mat = mat)
}

compute_reference_likeness <- function(domain_obj) {
  # "Reference likeness" is the Euclidean distance from each sample to the
  # historical-control centroid after feature scaling. It is a support metric for
  # platform shift review, not an MRD/tumor-fraction decision rule.
  mat <- standardize_feature_matrix(domain_obj$mat)
  info <- get_domain_sample_info(domain_obj) %>% distinct(Sample, Platform, Material, Domain)
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) return(tibble())
  historical_samples <- info %>%
    filter(Platform == "Historical_NovaSeq6000") %>%
    pull(Sample)
  historical_samples <- intersect(historical_samples, rownames(mat))
  if (length(historical_samples) < 2) return(tibble())
  historical_centroid <- colMeans(mat[historical_samples, , drop = FALSE], na.rm = TRUE)
  distance_to_historical_centroid <- apply(mat, 1, function(x) sqrt(sum((x - historical_centroid)^2, na.rm = TRUE)))
  historical_mat <- mat[historical_samples, , drop = FALSE]
  nearest_historical_distance <- apply(mat, 1, function(x) {
    min(apply(historical_mat, 1, function(y) sqrt(sum((x - y)^2, na.rm = TRUE))), na.rm = TRUE)
  })
  tibble(
    Sample = rownames(mat),
    Domain = domain_obj$domain,
    distance_to_historical_centroid = as.numeric(distance_to_historical_centroid),
    nearest_historical_distance = as.numeric(nearest_historical_distance),
    n_features = domain_obj$n_features
  ) %>%
    left_join(info, by = c("Sample", "Domain"))
}

compute_domain_feature_tests <- function(domain_obj) {
  # Per-feature tests are exploratory domain diagnostics. FDR is applied within
  # each domain so reviewers can identify features most responsible for a shift,
  # but these p-values do not update the production fragmentomics model.
  mat <- impute_feature_median_matrix(domain_obj$mat)
  info <- get_domain_sample_info(domain_obj) %>%
    distinct(Sample, Platform, Material, Domain)
  feature_long <- as_tibble(mat, rownames = "Sample") %>%
    pivot_longer(-Sample, names_to = "feature_id", values_to = "value") %>%
    left_join(info, by = "Sample")

  feature_long %>%
    group_by(Domain, feature_id) %>%
    summarise(
      n_historical = sum(Platform == "Historical_NovaSeq6000" & is.finite(value)),
      n_xplus = sum(Platform == "Spring2026_Xplus" & is.finite(value)),
      median_historical = median(value[Platform == "Historical_NovaSeq6000"], na.rm = TRUE),
      median_xplus = median(value[Platform == "Spring2026_Xplus"], na.rm = TRUE),
      median_difference_xplus_minus_historical = median_xplus - median_historical,
      wilcox_p = safe_wilcox(
        value[Platform == "Historical_NovaSeq6000"],
        value[Platform == "Spring2026_Xplus"]
      ),
      .groups = "drop"
    ) %>%
    group_by(Domain) %>%
    mutate(wilcox_fdr = p.adjust(wilcox_p, method = "fdr")) %>%
    ungroup() %>%
    arrange(Domain, wilcox_fdr, desc(abs(median_difference_xplus_minus_historical)))
}

summarise_group_difference <- function(data, feature_name, strata_cols = character()) {
  # Shared summary helper for scalar features. It reports medians/IQRs and a
  # Wilcoxon p-value for Spring2026_Xplus versus Historical_NovaSeq6000, with
  # optional strata such as nucleosome site/type.
  group_cols <- c(strata_cols, "Platform")

  summary_df <- data %>%
    filter(!is.na(.data[[feature_name]])) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n = dplyr::n(),
      median = median(.data[[feature_name]], na.rm = TRUE),
      iqr = IQR(.data[[feature_name]], na.rm = TRUE),
      mean = mean(.data[[feature_name]], na.rm = TRUE),
      sd = sd(.data[[feature_name]], na.rm = TRUE),
      .groups = "drop"
    )

  comparison_df <- data %>%
    filter(!is.na(.data[[feature_name]])) %>%
    group_by(across(all_of(strata_cols))) %>%
    summarise(
      n_historical = sum(Platform == "Historical_NovaSeq6000"),
      n_xplus = sum(Platform == "Spring2026_Xplus"),
      median_historical = median(.data[[feature_name]][Platform == "Historical_NovaSeq6000"], na.rm = TRUE),
      median_xplus = median(.data[[feature_name]][Platform == "Spring2026_Xplus"], na.rm = TRUE),
      median_difference_xplus_minus_historical = median_xplus - median_historical,
      wilcox_p = safe_wilcox(
        .data[[feature_name]][Platform == "Historical_NovaSeq6000"],
        .data[[feature_name]][Platform == "Spring2026_Xplus"]
      ),
      .groups = "drop"
    )

  list(summary = summary_df, comparison = comparison_df)
}

insert_df <- bind_rows(
  read_insert(historical_insert_path, "Historical_NovaSeq6000"),
  read_insert(xplus_insert_path, "Spring2026_Xplus")
) %>%
  filter_included_materials("insert-size summary")

fs_df <- bind_rows(
  read_fs(historical_fs_path, "Historical_NovaSeq6000"),
  read_fs(xplus_fs_path, "Spring2026_Xplus")
) %>%
  filter_included_materials("fragment score")

sample_level_features <- full_join(
  insert_df,
  fs_df,
  by = c("Sample", "Platform", "Material")
)

readr::write_csv(
  sample_level_features,
  file.path(out_dir, "sample_level_feature_values.csv")
)

insert_features <- c(
  "Proportion.Short",
  "Proportion.Long",
  "Ratio.Short.Normal",
  "MEDIAN_INSERT_SIZE",
  "MEAN_INSERT_SIZE"
)

insert_comparison <- bind_rows(lapply(insert_features, function(feature) {
  all_samples <- summarise_group_difference(insert_df, feature)$comparison %>%
    mutate(Feature = feature, Stratum = "cfDNA_only")
  all_samples
}))

fs_comparison <- summarise_group_difference(fs_df, "FS")$comparison %>%
  mutate(Feature = "FS", Stratum = "cfDNA_only")

readr::write_csv(insert_comparison, file.path(out_dir, "insert_size_group_comparison.csv"))
readr::write_csv(fs_comparison, file.path(out_dir, "fragment_score_group_comparison.csv"))

compute_nucleosome_metrics <- function(path, platform) {
  # Reduce per-position nucleosome accessibility traces to sample/site metrics.
  # Midpoint.normalized centers the midpoint coverage relative to mean local
  # coverage so between-sample depth shifts do not dominate the comparison.
  nuc <- readr::read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      Sample = as.character(Sample),
      site_name = as.character(site_name),
      site_type = as.character(site_type),
      Position = as.numeric(Position),
      Coverage = as.numeric(Coverage)
    )

  metric_df <- nuc %>%
    group_by(Sample, site_name, site_type) %>%
    summarise(
      Mean.Coverage = mean(Coverage, na.rm = TRUE),
      Midpoint.Coverage = mean(Coverage[Position %in% c(-30, -15, 0, 15, 30)], na.rm = TRUE),
      Midpoint.normalized = Midpoint.Coverage - Mean.Coverage + 1,
      .groups = "drop"
    ) %>%
    mutate(
      Platform = platform,
      Material = classify_material(Sample)
    )

  if (requireNamespace("GeneCycle", quietly = TRUE)) {
    amplitude_df <- nuc %>%
      group_by(Sample, site_name, site_type) %>%
      summarise(
        Amplitude = max(GeneCycle::periodogram(Coverage)[["spec"]], na.rm = TRUE),
        .groups = "drop"
      )
    metric_df <- left_join(metric_df, amplitude_df, by = c("Sample", "site_name", "site_type"))
  } else {
    metric_df$Amplitude <- NA_real_
    warning("Package GeneCycle is unavailable; nucleosome Amplitude was not calculated.", call. = FALSE)
  }

  metric_df
}

nuc_metrics <- bind_rows(
  compute_nucleosome_metrics(historical_nuc_path, "Historical_NovaSeq6000"),
  compute_nucleosome_metrics(xplus_nuc_path, "Spring2026_Xplus")
) %>%
  filter_included_materials("nucleosome accessibility metrics")

nuc_metric_features <- c("Mean.Coverage", "Midpoint.Coverage", "Midpoint.normalized", "Amplitude")
nuc_comparison <- bind_rows(lapply(nuc_metric_features, function(feature) {
  summarise_group_difference(nuc_metrics, feature, c("site_name", "site_type"))$comparison %>%
    mutate(Feature = feature)
}))

nuc_comparison <- nuc_comparison %>%
  group_by(Feature) %>%
  mutate(wilcox_fdr = p.adjust(wilcox_p, method = "fdr")) %>%
  ungroup() %>%
  arrange(Feature, wilcox_fdr, desc(abs(median_difference_xplus_minus_historical)))

readr::write_csv(
  nuc_comparison,
  file.path(out_dir, "nucleosome_site_metric_group_comparison.csv")
)

summary_feature_cols <- c(
  "Proportion.Short",
  "Proportion.Long",
  "Ratio.Short.Normal",
  "MEDIAN_INSERT_SIZE",
  "MEAN_INSERT_SIZE",
  "FS"
)
summary_feature_mat <- sample_level_features %>%
  select(Sample, all_of(summary_feature_cols)) %>%
  as.data.frame()
summary_feature_samples <- summary_feature_mat$Sample
summary_feature_mat <- summary_feature_mat[, summary_feature_cols, drop = FALSE] %>%
  as.matrix()
mode(summary_feature_mat) <- "numeric"
rownames(summary_feature_mat) <- summary_feature_samples

summary_domain <- list(
  # The summary domain combines the scalar fragment-size and fragment-score
  # features used in earlier QC plots, giving a compact PCA/distance check
  # alongside the higher-dimensional domains below.
  domain = "summary_insert_size_and_fragment_score",
  mat = impute_feature_median_matrix(summary_feature_mat),
  sample_info = sample_level_features %>%
    distinct(Sample, Platform, Material) %>%
    mutate(Domain = "summary_insert_size_and_fragment_score"),
  n_features = length(summary_feature_cols)
)

historical_insert_profile <- read_insert_profile_matrix(
  historical_insert_counts_path,
  "Historical_NovaSeq6000"
) %>%
  filter_matrix_object_to_included_materials("historical insert-size profile")
xplus_insert_profile <- read_insert_profile_matrix(
  xplus_insert_counts_path,
  "Spring2026_Xplus"
) %>%
  filter_matrix_object_to_included_materials("X Plus insert-size profile")

historical_breakpoint <- read_long_feature_matrix(
  historical_breakpoint_path,
  "Historical_NovaSeq6000",
  sample_col = "Sample",
  feature_cols = c("nucleotide", "Position"),
  value_col = "Frequency"
) %>%
  filter_matrix_object_to_included_materials("historical breakpoint frequency")
xplus_breakpoint <- read_long_feature_matrix(
  xplus_breakpoint_path,
  "Spring2026_Xplus",
  sample_col = "Sample",
  feature_cols = c("nucleotide", "Position"),
  value_col = "Frequency"
) %>%
  filter_matrix_object_to_included_materials("X Plus breakpoint frequency")

historical_per5mb <- read_feature_wide_matrix(
  historical_per5mb_path,
  "Historical_NovaSeq6000",
  id_cols = c("bin", "seqnames", "arm", "start", "end"),
  feature_builder = function(df) paste0("bin_", df$bin, "|", df$seqnames, ":", df$start, "-", df$end)
) %>%
  filter_matrix_object_to_included_materials("historical per-5Mb fragment ratio")
xplus_per5mb <- read_feature_wide_matrix(
  xplus_per5mb_path,
  "Spring2026_Xplus",
  id_cols = c("bin", "seqnames", "arm", "start", "end"),
  feature_builder = function(df) paste0("bin_", df$bin, "|", df$seqnames, ":", df$start, "-", df$end)
) %>%
  filter_matrix_object_to_included_materials("X Plus per-5Mb fragment ratio")

historical_endmotif_3prime <- read_feature_wide_matrix(
  historical_endmotif_3prime_path,
  "Historical_NovaSeq6000",
  id_cols = "motif",
  feature_builder = function(df) paste0("motif_", df$motif)
) %>%
  filter_matrix_object_to_included_materials("historical 3-prime end motif")
historical_endmotif_5prime <- read_feature_wide_matrix(
  historical_endmotif_5prime_path,
  "Historical_NovaSeq6000",
  id_cols = "motif",
  feature_builder = function(df) paste0("motif_", df$motif)
) %>%
  filter_matrix_object_to_included_materials("historical 5-prime end motif")
xplus_endmotif_5prime <- read_feature_wide_matrix(
  xplus_endmotif_5prime_path,
  "Spring2026_Xplus",
  id_cols = "motif",
  feature_builder = function(df) paste0("motif_", df$motif)
) %>%
  filter_matrix_object_to_included_materials("X Plus 5-prime end motif")
xplus_endmotif_3prime <- read_feature_wide_matrix(
  xplus_endmotif_3prime_path,
  "Spring2026_Xplus",
  id_cols = "motif",
  feature_builder = function(df) paste0("motif_", df$motif)
) %>%
  filter_matrix_object_to_included_materials("X Plus 3-prime end motif")

nuc_metric_mat <- nuc_metrics %>%
  transmute(
    Sample,
    feature_id = paste(site_name, site_type, "Midpoint.normalized", sep = "|"),
    value = Midpoint.normalized
  ) %>%
  pivot_wider(names_from = feature_id, values_from = value) %>%
  as.data.frame()
nuc_metric_samples <- nuc_metric_mat$Sample
nuc_metric_mat <- nuc_metric_mat[, setdiff(names(nuc_metric_mat), "Sample"), drop = FALSE] %>%
  as.matrix()
mode(nuc_metric_mat) <- "numeric"
rownames(nuc_metric_mat) <- nuc_metric_samples
nuc_domain <- list(
  domain = "nucleosome_midpoint_normalized",
  mat = impute_feature_median_matrix(nuc_metric_mat),
  sample_info = nuc_metrics %>%
    distinct(Sample, Platform, Material) %>%
    mutate(Domain = "nucleosome_midpoint_normalized"),
  n_features = ncol(nuc_metric_mat)
)

domain_objects <- list(
  # Each domain object contains a matrix and matching sample metadata. Objects
  # with too few shared features are dropped by combine_domain_matrices() and
  # excluded from downstream PCA/reference-likeness summaries.
  summary_domain,
  combine_domain_matrices(historical_insert_profile, xplus_insert_profile, "insert_size_profile_20_400bp"),
  combine_domain_matrices(historical_breakpoint, xplus_breakpoint, "breakpoint_frequency"),
  combine_domain_matrices(historical_per5mb, xplus_per5mb, "per5mb_fragment_ratio"),
  nuc_domain,
  combine_domain_matrices(
    historical_endmotif_5prime,
    xplus_endmotif_5prime,
    "end_motif_5prime"
  ),
  combine_domain_matrices(
    historical_endmotif_3prime,
    xplus_endmotif_3prime,
    "end_motif_3prime"
  )
) %>%
  Filter(Negate(is.null), .)

domain_pca <- lapply(domain_objects, compute_domain_pca) %>%
  Filter(Negate(is.null), .)

pca_scores <- bind_rows(lapply(domain_pca, `[[`, "scores"))
pca_variance <- bind_rows(lapply(domain_pca, `[[`, "variance"))
reference_likeness <- bind_rows(lapply(domain_objects, compute_reference_likeness))
domain_feature_tests <- bind_rows(lapply(domain_objects, compute_domain_feature_tests))

readr::write_csv(
  pca_scores,
  file.path(out_dir, "fragmentomics_domain_pca_scores.csv")
)
readr::write_csv(
  pca_variance,
  file.path(out_dir, "fragmentomics_domain_pca_variance.csv")
)
readr::write_csv(
  reference_likeness,
  file.path(out_dir, "fragmentomics_reference_likeness_distances.csv")
)
readr::write_csv(
  domain_feature_tests,
  file.path(out_dir, "fragmentomics_domain_feature_shift_tests.csv")
)

plot_features <- sample_level_features %>%
  # Plot only scalar sample-level features here; high-dimensional domain shifts
  # are summarized separately by PCA, reference-likeness, and per-feature CSVs.
  select(Sample, Platform, Material, Proportion.Short, Ratio.Short.Normal, MEDIAN_INSERT_SIZE, FS) %>%
  pivot_longer(
    cols = c(Proportion.Short, Ratio.Short.Normal, MEDIAN_INSERT_SIZE, FS),
    names_to = "Feature",
    values_to = "Value"
  )

p_box <- ggplot(plot_features, aes(x = Platform, y = Value, colour = Platform)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.3, alpha = 0.8) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("Historical_NovaSeq6000" = "#4477AA", "Spring2026_Xplus" = "#CC6677")) +
  labs(
    title = "Historical vs X Plus CHARM healthy controls",
    x = NULL,
    y = "Feature value"
  ) +
  theme_fragmentomics_qc() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  filename = file.path(out_dir, "xplus_hc_comparison_boxplots.pdf"),
  plot = p_box,
  width = 7,
  height = 5.5,
  units = "in"
)

top_nuc_shift <- nuc_comparison %>%
  filter(Feature == "Midpoint.normalized") %>%
  slice_max(order_by = abs(median_difference_xplus_minus_historical), n = 20, with_ties = FALSE) %>%
  mutate(site_name = reorder(site_name, median_difference_xplus_minus_historical))

p_nuc <- ggplot(
  top_nuc_shift,
  aes(x = median_difference_xplus_minus_historical, y = site_name)
) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(size = -log10(wilcox_p), colour = wilcox_fdr < 0.05), alpha = 0.85) +
  scale_colour_manual(values = c("FALSE" = "grey45", "TRUE" = "#CC6677"), na.value = "grey70") +
  labs(
    title = "Largest X Plus shifts in healthy-control nucleosome metrics",
    x = "Median difference in midpoint-normalized coverage\n(X Plus - historical)",
    y = NULL,
    colour = "FDR < 0.05",
    size = "-log10 p"
  ) +
  theme_fragmentomics_qc()

ggsave(
  filename = file.path(out_dir, "xplus_hc_nucleosome_site_shift.pdf"),
  plot = p_nuc,
  width = 7,
  height = 5.5,
  units = "in"
)

domain_labels <- pca_scores %>%
  distinct(Domain, n_features) %>%
  left_join(
    pca_variance %>%
      filter(PC == "PC1") %>%
      transmute(Domain, pc1_var = variance_explained),
    by = "Domain"
  ) %>%
  left_join(
    pca_variance %>%
      filter(PC == "PC2") %>%
      transmute(Domain, pc2_var = variance_explained),
    by = "Domain"
  ) %>%
  mutate(
    facet_label = paste0(
      Domain,
      "\n", n_features, " features; PC1 ",
      round(pc1_var * 100, 1),
      "%, PC2 ",
      round(pc2_var * 100, 1),
      "%"
    )
  )

pca_plot_df <- pca_scores %>%
  left_join(domain_labels %>% select(Domain, facet_label), by = "Domain") %>%
  mutate(
    Platform = factor(Platform, levels = c("Historical_NovaSeq6000", "Spring2026_Xplus"))
  )

p_pca <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, colour = Platform)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_point(size = 1.9, alpha = 0.82) +
  facet_wrap(~ facet_label, scales = "free", ncol = 2) +
  scale_colour_manual(values = c("Historical_NovaSeq6000" = "#4477AA", "Spring2026_Xplus" = "#CC6677")) +
  labs(
    title = "PCA of healthy-control fragmentomics feature domains",
    x = "PC1",
    y = "PC2",
    colour = "Healthy-control set"
  ) +
  theme_fragmentomics_qc(base_size = 7.5) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 6.8)
  )

ggsave(
  filename = file.path(out_dir, "xplus_hc_domain_pca_facets.pdf"),
  plot = p_pca,
  width = 8.5,
  height = 9.5,
  units = "in"
)

reference_summary <- reference_likeness %>%
  group_by(Domain, Platform) %>%
  summarise(
    n = sum(is.finite(distance_to_historical_centroid)),
    median_distance_to_historical_centroid = median(distance_to_historical_centroid, na.rm = TRUE),
    iqr_distance_to_historical_centroid = IQR(distance_to_historical_centroid, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(
  reference_summary,
  file.path(out_dir, "fragmentomics_reference_likeness_summary.csv")
)

p_ref <- reference_likeness %>%
  mutate(
    Platform = factor(Platform, levels = c("Historical_NovaSeq6000", "Spring2026_Xplus")),
    Domain = factor(Domain, levels = unique(domain_labels$Domain))
  ) %>%
  ggplot(aes(x = Platform, y = distance_to_historical_centroid, colour = Platform)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    size = 1.5,
    alpha = 0.8
  ) +
  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("Historical_NovaSeq6000" = "#4477AA", "Spring2026_Xplus" = "#CC6677")) +
  labs(
    title = "Distance to historical healthy-control centroid",
    x = NULL,
    y = "Euclidean distance after feature scaling",
    colour = "Healthy-control set"
  ) +
  theme_fragmentomics_qc(base_size = 7.5) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 6.8)
  )

ggsave(
  filename = file.path(out_dir, "xplus_hc_reference_likeness.pdf"),
  plot = p_ref,
  width = 8.5,
  height = 9.5,
  units = "in"
)

message("Wrote CHARM X Plus healthy-control comparison outputs to: ", out_dir)
