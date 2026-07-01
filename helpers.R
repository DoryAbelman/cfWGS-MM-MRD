# =============================================================================
# helpers.R
# Project:  cfWGS MRD Detection in Multiple Myeloma
# Author:   Dory Abelman
#
# Purpose:
#   Shared utility functions used across the numbered analysis scripts.
#   Source this file at the start of any script that calls these functions.
#
# Functions defined here:
#   clean_sample_id(x)  - strips whitespace and converts sample IDs to
#                         uppercase for consistent joining across data sources.
#
# Usage:
#   source("helpers.R")  # called at the top of each numbered script
#
# Manuscript outputs created/updated:
#   - None directly. This support file centralizes reusable helper functions
#     for consistent sample-ID handling.
# =============================================================================

# Collection of helper functions used across scripts

clean_sample_id <- function(x) {
  # remove spaces and convert to uppercase
  toupper(gsub("\\s+", "", x))
}

if (!exists("%>%", mode = "function")) {
  `%>%` <- magrittr::`%>%`
}

project_file <- function(...) {
  file.path(...)
}

spring2026_revision_data_dir <- function() {
  Sys.getenv(
    "CFWGS_SPRING2026_REVISION_DATA_DIR",
    unset = "Data_Spring_2026_Revisions"
  )
}

spring2026_revision_metadata_path <- function() {
  Sys.getenv(
    "CFWGS_SPRING2026_REVISION_METADATA",
    unset = file.path(
      "New OICR Submissions",
      "derived_metadata",
      "oicr_revision_repo_style_metadata.csv"
    )
  )
}

spring2026_revision_timepoint_override_path <- function() {
  Sys.getenv(
    "CFWGS_SPRING2026_TIMEPOINT_OVERRIDES",
    unset = file.path(
      "New OICR Submissions",
      "derived_metadata",
      "oicr_revision_timepoint_overrides.csv"
    )
  )
}

spring2026_revision_primary_analysis_exclusion_path <- function() {
  Sys.getenv(
    "CFWGS_SPRING2026_PRIMARY_ANALYSIS_EXCLUSIONS",
    unset = file.path(
      "New OICR Submissions",
      "derived_metadata",
      "oicr_revision_primary_analysis_exclusions.csv"
    )
  )
}

spring2026_revision_enabled <- function() {
  identical(Sys.getenv("CFWGS_USE_SPRING2026_REVISION", unset = "1"), "1")
}

require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

xplus_charm_healthy_bams <- function() {
  c(
    "TGL49_0197_Cf_R_WG_HBC-01-0002-T2-P-DNA_NovaX.filter.deduped.recalibrated.0.6.bam",
    "TGL49_0199_Cf_R_WG_HBC-01-0004-T2-P-DNA_NovaX.filter.deduped.recalibrated.0.3.bam",
    "TGL49_0255_Ct_n_WG.filter.deduped.recalibrated.0.6.bam",
    "TGL49_0256_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0257_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0259_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0260_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0261_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0262_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0263_Ct_n_WG.filter.deduped.recalibrated.bam",
    "TGL49_0264_Pb_U_WG_HBC-01-0018-T0-P-DNA-1_NovaX.filter.deduped.recalibrated.0.6.bam",
    "TGL49_0265_Pb_U_WG_HBC-01-0019-T0-P-DNA-1_NovaX.filter.deduped.recalibrated.0.7.bam",
    "TGL49_0266_Pb_U_WG_HBC-01-0020-T0-P-DNA-1_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0267_Cf_R_WG_HBC-01-0008-T1_NovaX.filter.deduped.recalibrated.0.55.bam",
    "TGL49_0267_Pb_U_WG_HBC-01-0021-T0-P-DNA-1.filter.deduped.recalibrated.0.7.bam",
    "TGL49_0268_Pb_U_WG_HBC-01-0022-T0-P-DNA-1.filter.deduped.recalibrated.bam",
    "TGL49_0302_Cf_R_WG_HBC-01-0023-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0303_Cf_R_WG_HBC-01-0024-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0304_Cf_R_WG_HBC-01-0025-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0305_Cf_R_WG_HBC-01-0026-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0306_Cf_R_WG_HBC-01-0027-T0-P-DNA_NovaX.filter.deduped.recalibrated.bam",
    "TGL49_0307_Cf_R_WG_HBC-01-0028-T0-P-DNA.filter.deduped.recalibrated.bam"
  )
}

is_xplus_charm_healthy_bam <- function(x) {
  basename(as.character(x)) %in% xplus_charm_healthy_bams()
}

load_spring2026_revision_metadata <- function(required = FALSE) {
  path <- spring2026_revision_metadata_path()
  if (!spring2026_revision_enabled()) return(NULL)
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing Spring 2026 revision metadata: ", path, call. = FALSE)
    return(NULL)
  }
  metadata <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  require_columns(
    metadata,
    c("Bam", "Patient", "Date_of_sample_collection", "Sample_type", "Timepoint",
      "Study", "Sample_ID", "timepoint_info", "integration_ready_for_combined_clinical_data"),
    "Spring 2026 revision metadata"
  )
  not_ready <- metadata$integration_ready_for_combined_clinical_data %in% FALSE
  if (any(not_ready, na.rm = TRUE)) {
    stop(
      "Spring 2026 revision metadata has rows not ready for combined clinical integration: ",
      paste(metadata$Sample_ID[not_ready], collapse = ", "),
      call. = FALSE
    )
  }
  override_path <- spring2026_revision_timepoint_override_path()
  if (file.exists(override_path)) {
    overrides <- utils::read.csv(override_path, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(
      overrides,
      c("Sample_ID", "timepoint_info"),
      "Spring 2026 revision timepoint override table"
    )
    duplicate_override_sample_id <- overrides$Sample_ID[
      !is.na(overrides$Sample_ID) & duplicated(overrides$Sample_ID)
    ]
    if (length(duplicate_override_sample_id)) {
      stop(
        "Spring 2026 revision timepoint override table has duplicate Sample_ID values: ",
        paste(unique(duplicate_override_sample_id), collapse = ", "),
        call. = FALSE
      )
    }
    unmatched_override_sample_id <- setdiff(overrides$Sample_ID, metadata$Sample_ID)
    if (length(unmatched_override_sample_id)) {
      stop(
        "Spring 2026 revision timepoint override table has Sample_ID values absent from revision metadata: ",
        paste(unmatched_override_sample_id, collapse = ", "),
        call. = FALSE
      )
    }
    idx <- match(overrides$Sample_ID, metadata$Sample_ID)
    metadata$timepoint_info_original_before_override <- metadata$timepoint_info
    metadata$timepoint_info[idx] <- overrides$timepoint_info
    metadata$timepoint_info_override_reason <- NA_character_
    if ("override_reason" %in% names(overrides)) {
      metadata$timepoint_info_override_reason[idx] <- overrides$override_reason
    }
  }
  exclusion_path <- spring2026_revision_primary_analysis_exclusion_path()
  if (file.exists(exclusion_path)) {
    exclusions <- utils::read.csv(exclusion_path, check.names = FALSE, stringsAsFactors = FALSE)
    require_columns(
      exclusions,
      c("Sample_ID", "exclusion_reason"),
      "Spring 2026 primary-analysis exclusion table"
    )
    duplicate_exclusion_sample_id <- exclusions$Sample_ID[
      !is.na(exclusions$Sample_ID) & duplicated(exclusions$Sample_ID)
    ]
    if (length(duplicate_exclusion_sample_id)) {
      stop(
        "Spring 2026 primary-analysis exclusion table has duplicate Sample_ID values: ",
        paste(unique(duplicate_exclusion_sample_id), collapse = ", "),
        call. = FALSE
      )
    }
    unmatched_exclusion_sample_id <- setdiff(exclusions$Sample_ID, metadata$Sample_ID)
    if (length(unmatched_exclusion_sample_id)) {
      stop(
        "Spring 2026 primary-analysis exclusion table has Sample_ID values absent from revision metadata: ",
        paste(unmatched_exclusion_sample_id, collapse = ", "),
        call. = FALSE
      )
    }
    metadata <- metadata[
      !(metadata$Sample_ID %in% exclusions$Sample_ID),
      ,
      drop = FALSE
    ]
  }
  metadata
}

coerce_revision_column_like_current <- function(current_col, revision_col, column_name) {
  if (grepl("date", column_name, ignore.case = TRUE)) {
    return(parse_date_safely(revision_col))
  }
  if (inherits(current_col, "Date")) {
    return(parse_date_safely(revision_col))
  }
  if (inherits(current_col, "POSIXct") || inherits(current_col, "POSIXt")) {
    return(as.POSIXct(revision_col, tz = attr(current_col, "tzone", exact = TRUE) %||% "UTC"))
  }
  if (is.character(current_col)) {
    return(as.character(revision_col))
  }
  if (is.numeric(current_col) && !is.logical(current_col)) {
    return(suppressWarnings(as.numeric(revision_col)))
  }
  if (is.integer(current_col)) {
    return(suppressWarnings(as.integer(revision_col)))
  }
  if (is.logical(current_col)) {
    if (all(is.na(revision_col) | revision_col %in% c(TRUE, FALSE, "TRUE", "FALSE", "true", "false", "1", "0"))) {
      return(as.logical(revision_col))
    }
    return(as.character(revision_col))
  }
  revision_col
}

parse_date_safely <- function(x) {
  if (inherits(x, "Date")) return(x)
  x <- as.character(x)
  x[!nzchar(x)] <- NA_character_
  out <- suppressWarnings(tryCatch(
    as.Date(x),
    error = function(e) rep(as.Date(NA), length(x))
  ))
  unresolved <- is.na(out) & !is.na(x)
  if (any(unresolved)) {
    formats <- c("%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d", "%m-%d-%Y", "%d-%m-%Y")
    for (fmt in formats) {
      idx <- unresolved
      parsed <- suppressWarnings(as.Date(x[idx], format = fmt))
      idx_positions <- which(idx)
      out[idx_positions[!is.na(parsed)]] <- parsed[!is.na(parsed)]
      unresolved <- is.na(out) & !is.na(x)
      if (!any(unresolved)) break
    }
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

spring2026_normalize_patient_alias <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[is.na(x)] <- NA_character_
  # MyP and IMG identifiers can name the same patient/sample in different
  # input systems. Normalize only patient-level overlap keys; do not rewrite
  # displayed Sample_ID/Bam values.
  y <- sub("^(MYP|IMG)[-_]?0*([0-9]+)$", "IMG-\\2", y)
  y
}

read_combined_clinical_metadata_with_revision <- function(
    path = "combined_clinical_data_updated_April2025.csv",
    include_revision_extra = FALSE) {
  if (!file.exists(path)) stop("Missing combined clinical metadata: ", path, call. = FALSE)
  current <- readr::read_csv(path, show_col_types = FALSE)
  revision <- load_spring2026_revision_metadata(required = FALSE)
  if (is.null(revision)) return(current)

  current_original <- current
  revision_original <- revision

  shared_cols <- if (isTRUE(include_revision_extra)) union(names(current), names(revision)) else names(current)
  extra_cols <- setdiff(names(revision), shared_cols)
  for (col in setdiff(shared_cols, names(current))) current[[col]] <- NA
  for (col in setdiff(shared_cols, names(revision))) revision[[col]] <- NA
  for (col in intersect(shared_cols, intersect(names(current), names(revision)))) {
    revision[[col]] <- coerce_revision_column_like_current(current[[col]], revision[[col]], col)
    if (grepl("date", col, ignore.case = TRUE) || inherits(revision[[col]], "Date")) {
      current[[col]] <- parse_date_safely(current[[col]])
    }
    if (is.character(revision[[col]]) && is.logical(current[[col]])) {
      current[[col]] <- as.character(current[[col]])
    }
  }
  revision_core <- revision[, shared_cols, drop = FALSE]
  current_core <- current[, shared_cols, drop = FALSE]

  duplicate_revision_sample_id <- revision_original$Sample_ID[
    !is.na(revision_original$Sample_ID) & duplicated(revision_original$Sample_ID)
  ]
  duplicate_revision_bam <- revision_original$Bam[
    !is.na(revision_original$Bam) & duplicated(revision_original$Bam)
  ]
  current_original <- current_original %>%
    dplyr::mutate(Date_of_sample_collection = parse_date_safely(.data$Date_of_sample_collection))
  revision_original <- revision_original %>%
    dplyr::mutate(Date_of_sample_collection = parse_date_safely(.data$Date_of_sample_collection))

  cross_sample_id <- intersect(
    unique(stats::na.omit(as.character(current_original$Sample_ID))),
    unique(stats::na.omit(as.character(revision_original$Sample_ID)))
  )
  cross_bam <- intersect(
    unique(stats::na.omit(as.character(current_original$Bam))),
    unique(stats::na.omit(as.character(revision_original$Bam)))
  )
  current_alias_keys <- current_original %>%
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Date_of_sample_collection), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      alias_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Date_of_sample_collection,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$alias_key) %>%
    unique()
  revision_alias_keys <- revision_original %>%
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Date_of_sample_collection), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      alias_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Date_of_sample_collection,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$alias_key) %>%
    unique()
  cross_alias_key <- intersect(current_alias_keys, revision_alias_keys)

  current_timepoint_type_keys <- current_original %>%
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Timepoint), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      timepoint_type_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Timepoint,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$timepoint_type_key) %>%
    unique()
  revision_timepoint_type_keys <- revision_original %>%
    dplyr::filter(!is.na(.data$Patient), !is.na(.data$Timepoint), !is.na(.data$Sample_type)) %>%
    dplyr::transmute(
      timepoint_type_key = paste(
        spring2026_normalize_patient_alias(.data$Patient),
        .data$Timepoint,
        .data$Sample_type,
        sep = "|"
      )
    ) %>%
    dplyr::pull(.data$timepoint_type_key) %>%
    unique()
  cross_timepoint_type_key <- intersect(current_timepoint_type_keys, revision_timepoint_type_keys)

  if (length(duplicate_revision_sample_id) || length(duplicate_revision_bam)) {
    stop(
      "Spring 2026 revision metadata has internal duplicate identifiers that need adjudication before automatic append.",
      "\nDuplicate revision Sample_ID: ", paste(unique(duplicate_revision_sample_id), collapse = ", "),
      "\nDuplicate revision Bam: ", paste(unique(duplicate_revision_bam), collapse = ", "),
      call. = FALSE
    )
  }

  if (length(cross_sample_id) || length(cross_bam) || length(cross_alias_key) || length(cross_timepoint_type_key)) {
    audit_dir <- "Output_tables_2025"
    if (!dir.exists(audit_dir)) dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
    overlap_audit <- dplyr::bind_rows(
      current_original %>%
        dplyr::mutate(
          spring2026_alias_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Date_of_sample_collection,
            .data$Sample_type,
            sep = "|"
          ),
          spring2026_timepoint_type_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Timepoint,
            .data$Sample_type,
            sep = "|"
          )
        ) %>%
        dplyr::filter(
          Sample_ID %in% cross_sample_id |
            Bam %in% cross_bam |
            spring2026_alias_key %in% cross_alias_key |
            spring2026_timepoint_type_key %in% cross_timepoint_type_key
        ) %>%
        dplyr::mutate(metadata_source_for_overlap_audit = "current") %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character)),
      revision_original %>%
        dplyr::mutate(
          spring2026_alias_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Date_of_sample_collection,
            .data$Sample_type,
            sep = "|"
          ),
          spring2026_timepoint_type_key = paste(
            spring2026_normalize_patient_alias(.data$Patient),
            .data$Timepoint,
            .data$Sample_type,
            sep = "|"
          )
        ) %>%
        dplyr::filter(
          Sample_ID %in% cross_sample_id |
            Bam %in% cross_bam |
            spring2026_alias_key %in% cross_alias_key |
            spring2026_timepoint_type_key %in% cross_timepoint_type_key
        ) %>%
        dplyr::mutate(metadata_source_for_overlap_audit = "spring2026_revision") %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    )
    readr::write_csv(
      overlap_audit,
      file.path(audit_dir, "spring2026_current_revision_identifier_overlap_audit.csv")
    )

    current_core <- current_core %>%
      dplyr::mutate(
        spring2026_alias_key = paste(
          spring2026_normalize_patient_alias(.data$Patient),
          .data$Date_of_sample_collection,
          .data$Sample_type,
          sep = "|"
        ),
        spring2026_timepoint_type_key = paste(
          spring2026_normalize_patient_alias(.data$Patient),
          .data$Timepoint,
          .data$Sample_type,
          sep = "|"
        )
      ) %>%
      dplyr::filter(
        !(Sample_ID %in% cross_sample_id),
        !(Bam %in% cross_bam),
        !(spring2026_alias_key %in% cross_alias_key),
        !(spring2026_timepoint_type_key %in% cross_timepoint_type_key)
      ) %>%
      dplyr::select(-spring2026_alias_key, -spring2026_timepoint_type_key)
  }

  combined <- dplyr::bind_rows(current_core, revision_core)

  if (length(extra_cols)) {
    attr(combined, "spring2026_revision_extra_columns") <- extra_cols
  }
  combined
}

spring2026_revision_maf_files <- function(sample_types) {
  metadata <- load_spring2026_revision_metadata(required = FALSE)
  if (is.null(metadata)) return(character())
  maf_dir <- file.path(spring2026_revision_data_dir(), "MuTect2_All_Mafs")
  if (!dir.exists(maf_dir)) return(character())

  wanted_pairs <- unique(metadata$mutect2_pair_id[
    metadata$Sample_type %in% sample_types &
      !is.na(metadata$mutect2_pair_id) &
      nzchar(metadata$mutect2_pair_id)
  ])
  files <- list.files(maf_dir, pattern = "[.]maf$", full.names = TRUE)
  if (!length(files) || !length(wanted_pairs)) return(character())
  keep <- vapply(files, function(path) {
    any(startsWith(basename(path), paste0(wanted_pairs, ".")))
  }, logical(1))
  sort(files[keep])
}

spring2026_revision_files <- function(subdir, pattern) {
  root <- file.path(spring2026_revision_data_dir(), subdir)
  if (!dir.exists(root)) return(character())
  sort(list.files(root, pattern = pattern, full.names = TRUE, recursive = FALSE))
}

spring2026_pwgval_dilution_metadata_path <- function() {
  file.path(
    spring2026_revision_data_dir(),
    "derived_metadata",
    "pwgval_dilution_metadata.csv"
  )
}

spring2026_pwgval_dilution_tumor_fraction_path <- function() {
  file.path(
    spring2026_revision_data_dir(),
    "derived_metadata",
    "pwgval_dilution_tumor_fraction.tsv"
  )
}

spring2026_pwgval_dilution_workbook_path <- function() {
  Sys.getenv(
    "CFWGS_SPRING2026_PWGVAL_DILUTION_WORKBOOK",
    unset = file.path(
      "Preparing shell scripts for cluster",
      "Dilution amounts for M4CHIP libs.xlsx"
    )
  )
}

build_spring2026_pwgval_dilution_metadata_from_workbook <- function(path = spring2026_pwgval_dilution_workbook_path()) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop(
      "Package readxl is required to read Spring 2026 PWGVAL dilution workbook: ",
      path,
      call. = FALSE
    )
  }
  workbook <- readxl::read_excel(path, sheet = 1)
  names(workbook) <- gsub("[^A-Za-z0-9]+", "_", names(workbook))
  names(workbook) <- gsub("_$", "", names(workbook))
  require_columns(
    workbook,
    c("Name", "Alias", "Group_ID", "TF_Dilution", "Tissue_Origin", "Tissue_Type",
      "Library_Design", "Index_es", "i7_Index", "i5_Index", "External_Identifier"),
    "Spring 2026 PWGVAL dilution workbook"
  )

  group_id <- as.character(workbook$Group_ID)
  tf_dilution <- suppressWarnings(as.numeric(workbook$TF_Dilution))
  replicate <- sub("^.*-", "", group_id)
  patient <- sub("^([A-Za-z]+-[0-9]+).*$", "\\1", group_id)
  m4chip_id <- sub("^(M4CHIP_[0-9]+).*$", "\\1", as.character(workbook$Alias))
  bam <- paste0(
    m4chip_id,
    "_",
    as.character(workbook$Tissue_Origin),
    "_",
    as.character(workbook$Tissue_Type),
    "_",
    as.character(workbook$Library_Design),
    "_",
    group_id,
    ".filter.deduped.recalibrated.bam"
  )

  workbook %>%
    dplyr::mutate(
      M4CHIP_ID = m4chip_id,
      Merge = paste0(group_id, "-P"),
      Patient = patient,
      replicate = replicate,
      Bam = bam,
      BAM = bam,
      LOD = tf_dilution,
      Sample_type = "Blood_plasma_cfDNA",
      Sample_type_Bam = "Blood_plasma_cfDNA",
      Timepoint = as.character(tf_dilution),
      timepoint_info = paste0("Dilution: ", as.character(tf_dilution)),
      Sample_ID = paste0(patient, "_PWGVAL_LOD_", as.character(tf_dilution), "_rep", replicate),
      dilution_series = "PWGVAL_M4CHIP",
      metadata_source = path
    ) %>%
    dplyr::rename(
      `TF_Dilution_pct` = "TF_Dilution",
      `Tissue Origin` = "Tissue_Origin",
      `Tissue Type` = "Tissue_Type",
      `Library Design` = "Library_Design",
      `Index(es)` = "Index_es",
      `i7 Index` = "i7_Index",
      `i5 Index` = "i5_Index",
      `External_Identifier` = "External_Identifier"
    )
}

load_spring2026_pwgval_dilution_metadata <- function(required = FALSE) {
  path <- spring2026_pwgval_dilution_metadata_path()
  if (!spring2026_revision_enabled()) return(NULL)
  if (!file.exists(path)) {
    metadata <- build_spring2026_pwgval_dilution_metadata_from_workbook()
    if (is.null(metadata)) {
      if (isTRUE(required)) {
        stop(
          "Missing Spring 2026 PWGVAL dilution metadata and workbook fallback: ",
          path,
          " / ",
          spring2026_pwgval_dilution_workbook_path(),
          call. = FALSE
        )
      }
      return(NULL)
    }
  } else {
    metadata <- readr::read_csv(path, show_col_types = FALSE)
  }
  require_columns(
    metadata,
    c("BAM", "Bam", "LOD", "Patient", "Sample_type", "Timepoint", "timepoint_info", "Sample_ID", "Merge"),
    "Spring 2026 PWGVAL dilution metadata"
  )
  metadata
}

build_spring2026_pwgval_dilution_tumor_fraction <- function() {
  metadata <- load_spring2026_pwgval_dilution_metadata(required = FALSE)
  params_path <- file.path(
    spring2026_revision_data_dir(),
    "Ichor_CNA_combined_params_summary_tumor_fraction_sex.tsv"
  )
  if (is.null(metadata) || !file.exists(params_path)) return(NULL)
  params <- readr::read_tsv(params_path, show_col_types = FALSE) %>%
    dplyr::filter(grepl("^M4CHIP_", .data$file)) %>%
    dplyr::mutate(Bam = sub("[.]params[.]txt$", "", .data$file))
  metadata %>%
    dplyr::mutate(Bam = sub("[.]bam$", "", .data$Bam)) %>%
    dplyr::inner_join(params, by = "Bam") %>%
    dplyr::transmute(
      Bam,
      Tumor_fraction = .data$tumor_fraction,
      Ploidy = .data$ploidy,
      Metadata_TF_Dilution_pct = .data$LOD,
      Patient = .data$Patient,
      Merge = .data$Merge,
      Group_ID = .data$Group_ID
    )
}

read_dilution_metadata_with_spring2026 <- function(path) {
  if (!file.exists(path)) stop("Missing dilution metadata: ", path, call. = FALSE)
  current <- readr::read_csv(path, show_col_types = FALSE)
  revision <- load_spring2026_pwgval_dilution_metadata(required = FALSE)
  if (is.null(revision)) return(current)

  shared_cols <- union(names(current), names(revision))
  for (col in setdiff(shared_cols, names(current))) current[[col]] <- NA
  for (col in setdiff(shared_cols, names(revision))) revision[[col]] <- NA
  current <- current[, shared_cols, drop = FALSE]
  revision <- revision[, shared_cols, drop = FALSE]

  combined <- dplyr::bind_rows(current, revision)
  if (!"Patient_Bam" %in% names(combined)) combined$Patient_Bam <- NA_character_
  if (!"Sample_type_Bam" %in% names(combined)) combined$Sample_type_Bam <- NA_character_
  if (!"Patient" %in% names(combined)) combined$Patient <- NA_character_
  if (!"Sample_type" %in% names(combined)) combined$Sample_type <- NA_character_

  combined %>%
    dplyr::mutate(
      Patient_Bam = dplyr::coalesce(.data$Patient_Bam, .data$Patient),
      Sample_type_Bam = dplyr::coalesce(.data$Sample_type_Bam, .data$Sample_type)
    ) %>%
    dplyr::distinct(BAM, Bam, Merge, .keep_all = TRUE)
}

read_dilution_tumor_fraction_with_spring2026 <- function(path) {
  if (!file.exists(path)) stop("Missing dilution tumor-fraction table: ", path, call. = FALSE)
  current <- readr::read_tsv(path, show_col_types = FALSE)
  revision_path <- spring2026_pwgval_dilution_tumor_fraction_path()
  if (!spring2026_revision_enabled()) return(current)
  revision <- if (file.exists(revision_path)) {
    readr::read_tsv(revision_path, show_col_types = FALSE)
  } else {
    build_spring2026_pwgval_dilution_tumor_fraction()
  }
  if (is.null(revision)) return(current)

  shared_cols <- union(names(current), names(revision))
  for (col in setdiff(shared_cols, names(current))) current[[col]] <- NA
  for (col in setdiff(shared_cols, names(revision))) revision[[col]] <- NA
  current <- current[, shared_cols, drop = FALSE]
  revision <- revision[, shared_cols, drop = FALSE]

  dplyr::bind_rows(current, revision) %>%
    dplyr::distinct(Bam, .keep_all = TRUE)
}
